import json
from pathlib import Path
import numpy as np
import pytest
from random_streams import STREAM_NAMES, seed_streams
from template_audit import overlap_matrix, count_pixel_map, regression, continuous_projections
from campaign4 import A_GRID, completion, finish, freeze, fingerprint, seed_phase
from validation_stats import slope_statistics


def test_stream_prefixes_and_reproducibility():
    a=seed_streams(104); b=seed_streams(104)
    draws=[a[name].normal(size=64) for name in STREAM_NAMES]
    assert all(np.array_equal(x,b[name].normal(size=64)) for name,x in zip(STREAM_NAMES,draws))
    assert all(not np.array_equal(x,y) for i,x in enumerate(draws) for y in draws[i+1:])


def test_centered_cells_and_exact_count_pixel_operator():
    for nin,nout in ((17,12),(16,16),(691,128)):
        w=overlap_matrix(nin,nout)
        assert np.allclose(w.sum(axis=1),1)
        assert np.allclose(w.sum(axis=0),nout/nin)
    n=64; signal=np.cos(2*np.pi*10*np.arange(n)/n)[:,None]*np.ones((1,n))
    truth=count_pixel_map(signal,32)
    assert regression(truth,truth,np.ones_like(truth))==pytest.approx(1)
    # The historical offset changes phase; its cross attenuation depends on rebinning.
    assert np.linalg.norm(count_pixel_map(signal,32,.5)-truth)/np.linalg.norm(truth)>.1


def test_continuous_template_range_and_linear_bias():
    n=32; z=np.arange(160)*.5
    a=.001*np.cos(2*np.pi*np.arange(n)/n)
    density=np.broadcast_to(a[:,None,None],(n,n,len(z))).copy()
    maps=continuous_projections(density,z,np.ones(len(z)),[10.,60.],0.,3.5)
    for margin in (0,150,300):
        truth=maps[f'kappa_range_{margin}']; continuous=maps[f'intensity_normalized_range_{margin}']
        assert regression(continuous,truth,np.ones_like(truth))==pytest.approx(1,abs=1e-5)
    assert np.allclose(maps['kappa_range_0'],density[:,:,0]*50)


def test_physical_grid_absolute_slope_and_offset():
    assert A_GRID==(0,.5,1,2)
    values=7+np.arange(10)[:,None]*.01+np.array(A_GRID)
    s=slope_statistics(values)
    assert s['mean']==pytest.approx(1)
    assert min(s['intercepts'])>=7


def test_completion_idempotency_and_tamper_rejection(tmp_path):
    prov={'source':'a'}; (tmp_path/'fit.h5').write_bytes(b'fit')
    result={'seed':0}; finish(tmp_path,prov,result)
    marker=tmp_path/'complete.json'; before=marker.stat().st_mtime_ns
    assert completion(tmp_path,prov)==result
    assert marker.stat().st_mtime_ns==before
    with pytest.raises(RuntimeError,match='provenance'): completion(tmp_path,{'source':'b'})
    (tmp_path/'fit.h5').write_bytes(b'changed')
    with pytest.raises(RuntimeError,match='artifact'): completion(tmp_path,prov)


def test_freeze_refuses_missing_development_seeds(tmp_path):
    with pytest.raises(RuntimeError,match='development seed 100'): freeze(tmp_path,.25)


def test_seed_phase_skips_before_heavy_work(tmp_path,monkeypatch):
    import campaign4 as c
    prov={'frozen':'test'}
    monkeypatch.setattr(c,'frozen',lambda root,scale:(None,prov))
    directory=tmp_path/'sparse/0'; directory.mkdir(parents=True)
    finish(directory,prov|{'variant':'sparse','seed':0},{'seed':0})
    monkeypatch.setattr(c.v,'process_seed',lambda *a:pytest.fail('completed seed regenerated'))
    assert seed_phase(tmp_path,0,'sparse',.25)=={'seed':0}


def test_disjoint_baseline_and_unchanged_randoms():
    # Small physical generation checks the actual marking/selection path.
    from config import Config
    from mock import generate_mock
    cfg=Config(scale=.03)
    a=generate_mock(cfg,104,disjoint_selection=True,pixel_noise_power=0,response=False)
    b=generate_mock(cfg,104,disjoint_selection=False,pixel_noise_power=0,response=False)
    assert np.intersect1d(a.quasars['qid'],a.sightlines.qid).size==0
    assert np.array_equal(a.sightlines.qid,b.sightlines.qid)
    assert all(np.array_equal(a.randoms[k],b.randoms[k]) for k in a.randoms)
    assert np.all(np.isin(a.sightlines.qid,b.quasars['qid']))
    assert a.attrs['quasar_radial_smoothing_mpc']==8.0
    # The smoothed lognormal input must have b_q sigma_s well below the raw-cell value.
    assert 3.5*np.sqrt(a.attrs['delta_g_variance'])<1.6


def test_band_regression_ignores_scales_outside_the_science_band():
    from template_audit import band_regression,regression
    rng=np.random.default_rng(11); n=128; side=np.deg2rad(20)
    L=2*np.pi*np.fft.fftfreq(n,d=side/n); LL=np.hypot(L[:,None],L[None,:])
    def field(lo,hi):
        f=np.fft.fft2(rng.normal(size=(n,n))); f[(LL<lo)|(LL>hi)]=0; return np.fft.ifft2(f).real
    low=field(20,300); high=field(400,1500); truth=low+high
    # A tracer that follows the truth with 0.9 in the science band but 1.5 on small scales.
    template=.9*low+1.5*high
    mask=np.ones((n,n)); mask[:8]=0
    assert band_regression(template,truth,mask,side)==pytest.approx(.9,abs=.02)
    assert regression(template*mask,truth*mask,mask)>1.2


def test_radially_smoothed_lognormal_keeps_transverse_response():
    from template_audit import continuous_projections,regression
    from scipy.ndimage import gaussian_filter1d
    rng=np.random.default_rng(7); n=48; nz=400; z=np.arange(nz)*.5
    # Transverse-only modes: radial smoothing must leave the projected truth/template relation exact.
    a=.02*np.cos(2*np.pi*3*np.arange(n)/n)[:,None]+.02*np.sin(2*np.pi*5*np.arange(n)/n)[None,:]
    density=np.broadcast_to(a[:,:,None],(n,n,nz)).copy()+.3*rng.normal(size=(n,n,nz))
    tracer=gaussian_filter1d(density,16,axis=2,mode='wrap')
    var=float(np.var(tracer)); maps=continuous_projections(density,z,np.ones(nz),[20.,180.],.5*3.5**2*var,3.5,tracer=tracer)
    truth=maps['kappa_range_0']; template=maps['intensity_normalized_range_0']
    assert regression(template,truth,np.ones_like(truth))==pytest.approx(1,abs=.03)


def test_discrete_grid_derivative_and_pixel_window_at_smoke_scale():
    from config import Config
    from grid_covariance import xi_from_mock_grid,derivative_comparison
    # Same transverse shape as scale .25; shorter periodic radial domain for unit cost.
    cfg=Config(scale=.25,xi_max=32,xi_step=.5)
    xi=xi_from_mock_grid(cfg,(173,173,192),nangle=32)
    inner=(xi.r_perp>5)&(xi.r_perp<30)
    fd=np.gradient(xi.xi,.5,axis=0)
    assert np.linalg.norm((fd-xi.xi_rp)[inner])/np.linalg.norm(xi.xi_rp[inner])<.03
    assert abs(derivative_comparison(xi,xi)['residual'])<1e-12
    assert xi.meta['empirical_factor']==1


def test_sampler_uses_node_centered_cells():
    from mock import sample_lognormal_quasars
    # Conditional offset must be centered on ix/iy, with periodic wrap.
    rng=np.random.default_rng(4)
    q=sample_lognormal_quasars(np.zeros((16,16,4)),2.,.5,100.,100.,1000.,rng)
    off=((q['x']/2+8-q['ix']+.5)%16)-.5
    assert abs(off.mean())<.005
    assert off.min()>=-.5 and off.max()<.5


def test_interrupted_phase_rejects_changed_provenance(tmp_path):
    from campaign4 import begin
    begin(tmp_path,{'source':'a'})
    (tmp_path/'partial.h5').write_bytes(b'partial')
    with pytest.raises(RuntimeError,match='interrupted'): begin(tmp_path,{'source':'b'})
    assert (tmp_path/'partial.h5').exists()
    begin(tmp_path,{'source':'a'})
    assert not (tmp_path/'partial.h5').exists()


def test_continuous_catalogue_includes_same_radial_selection():
    from template_audit import continuous_catalogue_projection
    from cross_spectrum import kernel,Z_CMB
    from cosmo import chi
    z=np.arange(3100.,3800.5,.5)
    a=.001*np.cos(2*np.pi*np.arange(8)/8)
    density=np.broadcast_to(a[:,None,None],(8,8,len(z))).copy()
    predicted=continuous_catalogue_projection(density,np.zeros_like(density),z,2.,.5,3400.,
                                              [3400.,3500.],0.,np.ones((8,8)),None)
    truth=continuous_projections(density,z,kernel(z,float(chi(Z_CMB))),[3400.,3500.],0.,3.5)
    for margin in (0,150,300):
        x=truth[f'kappa_range_{margin}']; y=predicted[f'continuous_catalogue_range_{margin}']
        # Forty-bin radial normalization has O(dchi/bin_width) edge quadrature.
        assert regression(y,x,np.ones_like(x))==pytest.approx(1,abs=2e-4)


def test_measured_xi_respects_even_separation_boundary():
    from config import Config
    from xi_model import xi_from_counts
    r=np.arange(40)+.5
    values=np.exp(-(r[:,None]**2+r[None,:]**2)/100)
    xi=xi_from_counts(values,np.ones_like(values),Config(xi_smoothing=.76))
    assert np.max(abs(xi.xi_rp[0]))<1e-12
    expected=-2*xi.r_perp[:,None]/100*np.exp(-(xi.r_perp[:,None]**2+xi.r_par[None,:]**2)/100)
    use=(xi.r_perp>5)&(xi.r_perp<30)
    assert np.linalg.norm((xi.xi_rp-expected)[use])/np.linalg.norm(expected[use])<.03
