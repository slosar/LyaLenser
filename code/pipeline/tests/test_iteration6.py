"""Regression checks for the review-5 fixes (iteration 6)."""
import numpy as np
import pytest
from config import Config, SightlineSet


def _two_forests(theta_deg, npix=60, chi0=4300., step=.55, chi_ref=3955.7):
    """Two full-weight forests at angular separation theta (degrees) on the equator."""
    chi=chi0+step*np.arange(npix)
    ra=np.array([180., 180.+theta_deg]); dec=np.zeros(2)
    return SightlineSet(np.arange(2),ra,dec,np.full(2,2.9),np.array([0,npix,2*npix]),np.tile(chi,2),
                        np.zeros(2*npix,np.float32),np.ones(2*npix,np.float32),np.zeros(2*npix,np.int8),{'chi_ref':chi_ref})


def _table():
    from xi_model import XiTable
    g=np.arange(0,40.25,.25); x=np.exp(-np.hypot(g[:,None],g[None,:])/6)
    return XiTable(g,g,x,np.gradient(x,.25,axis=0))


def test_response_prediction_applies_the_production_selection():
    """review 5, finding 9: the independent prediction must use cfg.r_perp_min/r_perp_max like the pair kernel."""
    from pairs import find_pairs, accumulate
    from response import prediction_catalogue
    t=_table(); sl=_two_forests(np.rad2deg(3.0/4316.))   # r_perp ~ 2.99-3.02 Mpc/h along the forests
    pairs=find_pairs(sl,30/sl.chi.min()); mod=np.full(len(sl.chi),.3)
    def score(rmin,rmax=30.):
        cfg=Config(chi_ref=3955.7,r_perp_min=rmin,r_perp_max=rmax)
        cat=accumulate(sl,pairs,t,cfg)
        return prediction_catalogue(sl,cat,mod,t,t,cfg).accum[:,0,:].sum()
    full=score(0.); assert full!=0
    assert score(3.2)==0                         # every pixel pair lies below the lower cut
    assert score(0.,2.8)==0                      # ... and above the upper cut
    partial=score(3.0)
    assert partial!=0 and abs(partial)<abs(full) # a straddling selection keeps a strict subset


def test_pair_kernel_expectation_mode_reproduces_mean_field():
    """delta delta -> xi(true separation): with true = observed positions the score equals the mean field."""
    from pairs import find_pairs, accumulate
    t=_table(); sl=_two_forests(.08); cfg=Config(chi_ref=3955.7,r_perp_min=3.)
    pairs=find_pairs(sl,30/sl.chi.min()); pos=np.column_stack((sl.ra,sl.dec))
    same=accumulate(sl,pairs,t,cfg,true_positions=pos)
    assert same.attrs['expectation'] and np.allclose(same.accum[:,0:3],same.accum[:,8:11],rtol=1e-6,atol=0)
    # A true separation larger than the observed one lowers xi(true) below xi(observed); the kernel G = chi xi' is
    # negative for this decaying table, so the score becomes less negative than the mean field: |q| < |mf|.
    wider=pos.copy(); wider[1,0]+=.01
    moved=accumulate(sl,pairs,t,cfg,true_positions=wider)
    q=moved.accum[:,0,:].sum(); mf=moved.accum[:,8,:].sum()
    assert mf<0 and q>mf and abs(q)<abs(mf)


def test_injection_slopes_by_amplitude():
    from inject import paired_slopes
    x=np.array([-2,-1,-.5,0,.5,1,2]); y=.3+.9*x+.05*x**2   # even parts (0.3, quadratic) drop out; odd slope 0.9
    slope,per=paired_slopes(x,y)
    assert slope==pytest.approx(.9) and set(per)=={.5,1.,2.} and all(v==pytest.approx(.9) for v in per.values())


def test_band_regression_masks_exactly_once():
    """review 5, finding 5: the sampled coefficient was computed on an already masked map. band_regression takes
    unmasked inputs; a soft mask applied twice biases the coefficient of an otherwise perfect template."""
    from template_audit import band_regression
    rng=np.random.default_rng(3); n=128; side=np.deg2rad(20.)
    L=2*np.pi*np.fft.fftfreq(n,d=side/n); LL=np.hypot(L[:,None],L[None,:])
    f=np.fft.fft2(rng.normal(size=(n,n))); f[(LL<40)|(LL>300)]=0; truth=np.fft.ifft2(f).real
    mask=np.clip(np.linspace(0,1.,n)[:,None]*np.ones((1,n)),0,1)      # soft apodisation
    assert band_regression(truth,truth,mask,side)==pytest.approx(1.,abs=1e-9)
    assert abs(band_regression(truth*mask,truth,mask,side)-1)>.05


def test_patch_side_is_the_generator_side():
    from types import SimpleNamespace
    from mock import patch_side_rad, grid_geometry
    cfg=Config(scale=1.); geo=grid_geometry(cfg); nx=geo['shape'][0]
    mock=SimpleNamespace(maps={'kappa_lya':np.zeros((nx,nx))},attrs={'dx':geo['dx']},sightlines=SimpleNamespace(attrs={'chi_ref':cfg.chi_ref}))
    side=patch_side_rad(mock)
    assert side==pytest.approx(nx*geo['dx']/cfg.chi_ref) and side>=np.deg2rad(20.) and side<np.deg2rad(20.)*(1+2./nx)


def test_rows_carry_the_required_flag():
    from run_mock_validation import record
    rows=[]
    record(rows,1,'a','|mean| <= 2 SEM',{},True); record(rows,2,'b','report only',{},False); record(rows,3,'c','x',{},True,required=False)
    assert [r['required'] for r in rows]==[True,False,False]


def test_fingerprint_records_the_campaign_configuration_and_nested_sources():
    import campaign4 as c
    fp=c.fingerprint(1.)
    assert fp['config']['r_perp_min']==3. and fp['config']['fit_rperp_min']==3. and fp['iteration']==6
    assert any(s.startswith('code/pipeline/tests/') for s in fp['sources']) and 'code/pipeline/campaign4.py' in fp['sources']
    assert fp['seeds']['sparse']==[1000,1399] and fp['deprojected_bound']==.5
