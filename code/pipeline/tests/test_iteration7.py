"""Iteration 7: low-redshift tracers (lowz.py), generator branch, bias fit, combination."""
import numpy as np
import pytest
from config import Config


def test_stream_names_prefix_unchanged():
    """The three iteration-7 streams are appended; the frozen first ten (and their draws) are unchanged."""
    from random_streams import STREAM_NAMES, seed_streams
    assert STREAM_NAMES[:10]==('field','outside_lya','outside_cmb','cmb_noise','quasar_sampling','sightline_selection',
                               'forest_noise','random_catalogue','catalogue_split','random_templates')
    assert STREAM_NAMES[10:]==('lowz_fields','lowz_sampling','lowz_randoms')
    old=np.random.default_rng(np.random.SeedSequence(5).spawn(10)[0]).normal(size=4)
    assert np.allclose(seed_streams(5)['field'].normal(size=4),old)


def test_bias_band_is_kmax_times_chi():
    from lowz import bias_band, KMAX_BIAS
    from cosmo import chi
    lo,hi=bias_band(.4,.6)
    assert lo==40. and abs(hi-(KMAX_BIAS*float(chi(.5))-.5))<1e-6
    assert bias_band(1.1,1.6)[1]>bias_band(.4,.6)[1]


def test_gaussian_maps_reproduce_the_covariance():
    from lowz import gaussian_maps
    from mock import flat_sky_power
    nx=256; side=.35; L=np.linspace(2,3000,300); cov=np.zeros((len(L),2,2))
    cov[:,0,0]=1e-8/(1+L/100); cov[:,1,1]=4e-8/(1+L/100); cov[:,0,1]=cov[:,1,0]=.6*np.sqrt(cov[:,0,0]*cov[:,1,1])
    rng=np.random.default_rng(3); acc=np.zeros(3)
    for _ in range(6):
        a,b=gaussian_maps(nx,nx,side,L,cov,rng)
        acc+=[flat_sky_power(a,a,side,40,300)[0],flat_sky_power(a,b,side,40,300)[0],flat_sky_power(b,b,side,40,300)[0]]
    acc/=6
    # expectation: the mode-weighted band average of the input covariance on the map's own Fourier grid
    lx=2*np.pi*np.fft.fftfreq(nx,side/nx); ly=2*np.pi*np.fft.rfftfreq(nx,side/nx); ell=np.hypot(lx[:,None],ly[None,:]); use=(ell>40)&(ell<300)
    expect=[np.mean(np.interp(ell[use],L,cov[:,i,j])) for i,j in ((0,0),(0,1),(1,1))]
    assert np.allclose(acc/expect,1,atol=.08)


def test_slices_plus_rest_equal_the_full_foreground():
    """Sum of the slice kappa_lya spectra and the uncovered rest equals the Limber integral over the whole foreground."""
    from lowz import slices_of, slice_spectra, rest_spectra, front_ranges, TRACERS
    from three_tracer import limber
    from cross_spectrum import kernel
    from cosmo import chi
    cref=float(chi(2.4)); cbox=(float(chi(2.1))-300.,float(chi(3.))+300.); sl=slices_of(TRACERS)
    L=np.array([40.,100.,300.]); tot=np.zeros(3)
    for s in sl:
        Ls,C=slice_spectra(s['zmin'],s['zmax'],cref,400); tot+=np.interp(L,Ls,C[:,1,1])
    Lr,ll,_,_=rest_spectra(front_ranges(sl,cbox[0]),cbox[1],cref,400); tot+=np.interp(L,Lr,ll)
    wl=lambda c: kernel(c,cref); full=limber(L,wl,wl,1.,cbox[0],nchi=2000,to_recombination=False)
    assert np.allclose(tot/full,1,atol=.01)


def test_lognormal_sampler_mean_density_and_correlation():
    from lowz import sample_lognormal_2d
    rng=np.random.default_rng(1); nx=128; side=.2; d=rng.normal(scale=.1,size=(nx,nx))
    g=sample_lognormal_2d(d,2.,500.,side,rng)
    expected=500.*np.rad2deg(side)**2
    assert abs(len(g['x'])-expected)<4*np.sqrt(expected)
    counts=np.zeros((nx,nx)); np.add.at(counts,(g['ix'],g['iy']),1)
    assert np.corrcoef(counts.ravel(),d.ravel())[0,1]>.3
    assert np.all(abs(g['x'])<=side/2+1e-9)


def test_fit_bias_recovers_b_on_synthetic_annuli():
    from lowz import fit_bias
    L=np.arange(20,700,40.); T=1e-9/(1+L/80); shot=2e-10; nmodes=np.full(len(L),400.)
    rng=np.random.default_rng(0); b=2.1
    meas=b*b*T+shot; meas=meas*(1+rng.normal(scale=np.sqrt(2/(nmodes/2))))
    fit=fit_bias(meas,shot,T,L,nmodes,band=(40,350))
    assert abs(fit['b']-b)<3*fit['sigma_b'] and fit['sigma_b']<.1 and fit['annuli_used']==8


def test_optimal_combination_of_correlated_slices():
    from lowz import optimal_combination
    rng=np.random.default_rng(2); k=4; nr=40
    common=rng.normal(size=nr); jk=np.array([1+.05*common+.03*rng.normal(size=nr) for _ in range(k)])
    A=jk.mean(axis=1)
    out=optimal_combination(A,jk)
    # the optimal error never exceeds the best single slice (up to the Hartlap correction of the inverse)
    assert np.isfinite(out['A']) and out['error']<=np.sqrt(np.diag(np.array(out['covariance']))/out['hartlap']).min()*(1+1e-9)
    assert abs(out['A']-1)<.5 and 0<out['hartlap']<1 and out['naive_error']>0


def test_matched_template_source_default_is_the_cmb():
    from templates import matched_template_flat
    from cross_spectrum import Z_CMB
    from cosmo import chi
    rng=np.random.default_rng(4); n=64; pix=np.deg2rad(20/n)
    q={'ra':180+rng.uniform(-9,9,3000)/np.cos(np.deg2rad(30)),'dec':30+rng.uniform(-9,9,3000),'z':rng.uniform(2.,3.,3000)}
    r={'ra':180+rng.uniform(-9,9,30000)/np.cos(np.deg2rad(30)),'dec':30+rng.uniform(-9,9,30000),'z':rng.uniform(2.,3.,30000)}
    m1,_,_=matched_template_flat(q,r,lambda z: np.full_like(np.asarray(z,float),3.5),(n,n),pix)
    m2,_,_=matched_template_flat(q,r,lambda z: np.full_like(np.asarray(z,float),3.5),(n,n),pix,source_chi=float(chi(Z_CMB)))
    m3,_,_=matched_template_flat(q,r,lambda z: np.full_like(np.asarray(z,float),3.5),(n,n),pix,source_chi=float(chi(2.4)))
    assert np.array_equal(m1,m2) and not np.array_equal(m1,m3)


def test_mock_with_lowz_tracers_end_to_end(tmp_path):
    """Generator branch, save/load, templates: kappa_lya_rest is the sum of the slice and rest maps; every template
    set has the seven components; the fixed templates of another map sample at this mock's sightlines."""
    from mock import generate_mock, save_mock
    from run_mock_validation import load_mock
    from lowz import lowz_bundles, tracer_maps, tracer_summary
    cfg=Config(scale=.08,r_perp_min=3.,fit_rperp_min=3.)
    m=generate_mock(cfg,5,A_true=1,response=True,completeness=True,magnification=True,cmb_noise=True,
                    disjoint_selection=True,variant_A_values=[0,1],lowz=True)
    nsl=len(m.attrs['lowz']['slices'])
    total=sum(m.maps[f'lowz_kappa_lya_{i}'] for i in range(nsl))+m.maps['lowz_rest_lya']
    assert np.allclose(total,m.maps['kappa_lya_rest'],atol=1e-6)
    assert len(m.lowz_catalogue['ra'])>0 and set(np.unique(m.lowz_catalogue['slice']))=={*range(nsl)}
    assert np.all(m.lowz_catalogue['z']<1.8)
    save_mock(m,tmp_path/'m.h5'); m2=load_mock(tmp_path/'m.h5')
    assert isinstance(m2.attrs['lowz'],dict) and len(m2.lowz_catalogue['ra'])==len(m.lowz_catalogue['ra'])
    tm=tracer_maps(m2,cfg); b=lowz_bundles(m2,cfg,fixed=tm)
    assert set(b['templates'])=={*(f'slice{i}' for i in range(nsl)),'combined','truth','fixed_combined',*(f'fixed_slice{i}' for i in range(nsl))}
    for name,ts in b['templates'].items():
        assert len(ts)==7 and sum(t.kind=='curl' for t in ts)==3 and sum(t.kind=='junk' for t in ts)==1
        assert ts[0].alpha.shape==(m2.sightlines.nq,2)
    summary=tracer_summary(tm)
    assert all('bias_source' in t for t in summary['tracers'].values())
