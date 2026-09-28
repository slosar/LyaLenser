"""The source-distance derivative of the templates: kernel derivative, derivative map of a tracer catalogue and
the band templates built from it (iteration 15)."""
import numpy as np
import healpy as hp
from lyalenser.lensing import kernel, kernel_dsource
from lyalenser.cosmo import chi as chi_of_z, z_of_chi


def test_kernel_derivative_is_the_source_distance_derivative():
    cs=float(chi_of_z(2.35)); chi=np.linspace(300.,3500.,50); eps=1.0
    fd=(kernel(chi,cs+eps)-kernel(chi,cs-eps))/(2*eps)
    np.testing.assert_allclose(kernel_dsource(chi,cs),fd,rtol=1e-6)
    assert kernel_dsource(np.array([cs+10.]),cs)[0]==0.


def _toy_catalogue(n=4000,seed=1):
    rng=np.random.default_rng(seed)
    ra=rng.uniform(20.,40.,n); dec=rng.uniform(-5.,5.,n); z=rng.uniform(.4,.6,n)
    return {"ra":ra,"dec":dec,"z":z}


def test_matched_template_derivative_map_is_finite_difference():
    """The dkappa map of matched_template equals the derivative of its kappa map with respect to the source
    distance: same objects, randoms, completeness and mask, dW/dchi_s in place of W."""
    from lyalenser.templates import matched_template
    from lyalenser.config import Config
    nside=32; cfg=Config(nside_alpha=nside,lmax_alpha=64); cs=float(chi_of_z(2.35)); eps=2.0
    data=_toy_catalogue(3000,1); rnd=_toy_catalogue(30000,2)
    pix=hp.ang2pix(nside,rnd["ra"],rnd["dec"],lonlat=True); fp=np.bincount(pix,minlength=hp.nside2npix(nside))>0
    unit=lambda z: np.ones_like(np.asarray(z,float))
    kw=dict(footprint_mask=fp,nside=nside,lmax=64,radial_bins=5)
    _,mask,meta=matched_template(data,rnd,unit,cfg,source_chi=cs,**kw)
    _,_,up=matched_template(data,rnd,unit,cfg,source_chi=cs+eps,**kw)
    _,_,dn=matched_template(data,rnd,unit,cfg,source_chi=cs-eps,**kw)
    fd=(up["kappa_map"]-dn["kappa_map"])/(2*eps)
    np.testing.assert_allclose(meta["dkappa_map"][mask],fd[mask],rtol=1e-5,atol=1e-12*np.abs(fd[mask]).max())
    assert np.abs(meta["dkappa_map"][mask]).max()>0


def test_band_templates_carry_the_derivative_of_the_same_filter():
    from lyalenser.templates import sphere_band_templates, phi_from_kappa, alpha_at, cosine_band, curl
    lmax=48; nside=32; rng=np.random.default_rng(3); size=hp.Alm.getsize(lmax)
    kalm=(rng.normal(size=size)+1j*rng.normal(size=size))*1e-5; dalm=(rng.normal(size=size)+1j*rng.normal(size=size))*1e-8
    ra=np.array([30.,31.,35.]); dec=np.array([0.,2.,-1.])
    ts,filt=sphere_band_templates(kalm,ra,dec,nside=nside,science_bands=((10,20),(20,30)),dkappa_alm=dalm)
    pos=np.column_stack((ra,dec))
    for t in ts:
        assert t.dalpha is not None and t.dalpha.shape==t.alpha.shape
        name=t.name.replace("_curl","") if t.kind=="curl" else t.name
        expect=alpha_at(pos,phi_from_kappa(hp.almxfl(dalm,filt[name]),lmax),nside)
        if t.kind=="curl": expect=curl(expect)
        np.testing.assert_allclose(t.dalpha,expect,rtol=1e-6,atol=1e-14)
    plain,_=sphere_band_templates(kalm,ra,dec,nside=nside,science_bands=((10,20),(20,30)))
    assert all(t.dalpha is None for t in plain)
    assert all(np.array_equal(a.alpha,b.alpha) for a,b in zip(ts,plain))
