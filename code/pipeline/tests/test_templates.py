import numpy as np
import healpy as hp
from templates import phi_from_kappa,alpha_at,flat_sky_band_filters


def test_phi_and_alm2map_der1_convention():
    nside=32;lmax=8
    kalm=np.zeros(hp.Alm.getsize(lmax),complex); kalm[hp.Alm.getidx(lmax,3,1)]=2e-4+1e-4j
    phi=phi_from_kappa(kalm,lmax)
    ell,m=hp.Alm.getlm(lmax)
    assert np.all(phi[ell<2]==0)
    pos=np.array([[40.,25.],[121.,-20.],[250.,52.]])
    got=alpha_at(pos,phi,nside)
    _,dt,dp=hp.alm2map_der1(phi,nside)
    th=np.deg2rad(90-pos[:,1]); ph=np.deg2rad(pos[:,0])
    expected=np.column_stack((hp.get_interp_val(dp,th,ph),-hp.get_interp_val(dt,th,ph)))
    assert np.allclose(got,expected,rtol=1e-6,atol=1e-10)


def test_low_l_divergence_equals_minus_two_kappa():
    nside=32;lmax=6
    rng=np.random.default_rng(2); kalm=(rng.normal(size=hp.Alm.getsize(lmax))+1j*rng.normal(size=hp.Alm.getsize(lmax)))*1e-6
    ell,m=hp.Alm.getlm(lmax); kalm[ell<2]=0
    phi=phi_from_kappa(kalm,lmax)
    # Spectral divergence of grad(phi) is laplacian(phi) = -2 kappa.
    lap=hp.almxfl(phi,-np.arange(lmax+1)*(np.arange(lmax+1)+1))
    assert np.allclose(lap,-2*kalm,rtol=1e-13,atol=1e-18)
    # Explicit central differences at an equatorial point.  At theta=pi/2,
    # div(alpha)=d(alpha_north)/d(dec)+d(alpha_east)/d(RA).
    eps=2e-4; p=np.array([[73.,0.],[73.+np.rad2deg(eps),0.],[73.-np.rad2deg(eps),0.],
                          [73.,np.rad2deg(eps)],[73.,-np.rad2deg(eps)]])
    aa=alpha_at(p,phi,nside)
    div=(aa[3,1]-aa[4,1])/(2*eps)+(aa[1,0]-aa[2,0])/(2*eps)
    km=hp.get_interp_val(hp.alm2map(kalm,nside),np.pi/2,np.deg2rad(73.))
    assert np.isclose(div,-2*km,rtol=.08,atol=2e-7)


def test_flat_sky_cosine_sign_amplitude():
    phi0=2e-6; L=80.; x=np.linspace(-.02,.02,1001); phi=phi0*np.cos(L*x)
    numeric=np.gradient(phi,x)
    truth=-phi0*L*np.sin(L*x)
    assert np.max(np.abs(numeric[2:-2]-truth[2:-2]))<3e-9


def test_flat_sky_band_filters_are_orthogonal():
    ell,filters=flat_sky_band_filters((128,128),np.deg2rad(6/128),taper=8)
    names=["L40_100","L100_200","L200_300","junk"]
    for i,a in enumerate(names):
        assert np.sum(filters[a]**2)>0
        for b in names[i+1:]:
            assert np.sum(filters[a]*filters[b])==0
    covered=(ell<40)|(ell>300)
    assert np.all(filters["junk"][covered]==1)
