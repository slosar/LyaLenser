"""Gradient/curl, Wiener, Gaussian and matched-density templates."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import sys
import numpy as np
import healpy as hp

HERE=Path(__file__).resolve().parent; CODE=HERE.parent
if str(CODE) not in sys.path: sys.path.insert(0,str(CODE))
from cosmo import chi as chi_of_z
from cross_spectrum import kernel, Z_CMB


@dataclass
class Template:
    alpha: np.ndarray
    name: str
    kind: str
    phi_lm: np.ndarray | None = None
    Lmin: int = 2
    Lmax: int = 0
    filter: str = "none"
    source: str = ""
    attrs: dict = field(default_factory=dict)

    def __post_init__(self):
        self.alpha=np.asarray(self.alpha,np.float32)
        if self.alpha.ndim != 2 or self.alpha.shape[1] != 2:
            raise ValueError("alpha must be [Nq,2] (east,north)")
        if self.kind not in {"signal","response","curl","random","injection","truth","junk"}:
            raise ValueError("invalid template kind")


def phi_from_kappa(kappa_alm, lmax=None):
    a=np.asarray(kappa_alm).copy()
    lm=hp.Alm.getlmax(len(a)) if lmax is None else lmax
    ell,_=hp.Alm.getlm(lm)
    factor=np.zeros_like(ell,dtype=float); m=ell>=2
    factor[m]=2.0/(ell[m]*(ell[m]+1.0))
    a*=factor
    a[ell<2]=0
    return a


def alpha_at(positions, phi_lm, nside):
    """Deflection (east,north); healpy's third derivative is already /sin(theta)."""
    pos=np.asarray(positions)
    if hasattr(positions,"ra"):
        ra=np.asarray(positions.ra); dec=np.asarray(positions.dec)
    else:
        ra=pos[:,0]; dec=pos[:,1]
    _,dtheta,dphi_over_sin=hp.alm2map_der1(phi_lm,nside)
    theta=np.deg2rad(90.0-dec); phi=np.deg2rad(ra)% (2*np.pi)
    east=hp.get_interp_val(dphi_over_sin,theta,phi)
    north=-hp.get_interp_val(dtheta,theta,phi)
    return np.column_stack((east,north)).astype(np.float32)


def cosine_band(ell,Lmin,Lmax,taper=10):
    ell=np.asarray(ell,float); w=np.zeros_like(ell)
    if Lmax <= Lmin: return w
    lo=min(float(taper),.5*(Lmax-Lmin)); hi=lo
    core=(ell>=Lmin+lo)&(ell<=Lmax-hi); w[core]=1
    if lo>0:
        m=(ell>=Lmin)&(ell<Lmin+lo); w[m]=.5*(1-np.cos(np.pi*(ell[m]-Lmin)/lo))
    if hi>0:
        m=(ell>Lmax-hi)&(ell<=Lmax); w[m]=.5*(1-np.cos(np.pi*(Lmax-ell[m])/hi))
    return w


def wiener_filter(X_alm,S_L=None,C_XX_L=None,Lmin=2,Lmax=None,taper=10,
                  spectra_args=None):
    lmax=hp.Alm.getlmax(len(X_alm)); ell,_=hp.Alm.getlm(lmax)
    if S_L is None:
        from three_tracer import spectra
        kw={} if spectra_args is None else dict(spectra_args)
        s=spectra(L_values=np.arange(2,lmax+1),**kw)
        S_L=np.r_[np.zeros(2),s["klkc"]-s["skl"]]
    S_L=np.asarray(S_L,float)
    if C_XX_L is None:
        C_XX_L=hp.alm2cl(X_alm)
    C_XX_L=np.asarray(C_XX_L,float)
    h=np.divide(S_L[:lmax+1],C_XX_L[:lmax+1],out=np.zeros(lmax+1),where=C_XX_L[:lmax+1]>0)
    top=lmax if Lmax is None else Lmax
    h*=cosine_band(np.arange(lmax+1),Lmin,top,taper)
    return hp.almxfl(X_alm,h),h


def curl(alpha):
    a=np.asarray(alpha)
    return np.column_stack((-a[:,1],a[:,0])).astype(np.float32)


def _catalog_columns(cat):
    if isinstance(cat,dict): return (np.asarray(cat[k]) for k in ("ra","dec","z"))
    return np.asarray(cat.ra),np.asarray(cat.dec),np.asarray(cat.zq)


def matched_template(quasars,randoms,b_q_of_z,cfg,nmin_rand=1):
    """Kernel-matched quasar overdensity and its mask/white shot-noise model."""
    qra,qdec,qz=_catalog_columns(quasars); rra,rdec,rz=_catalog_columns(randoms)
    nside=cfg.nside_alpha; npix=hp.nside2npix(nside); area=4*np.pi/npix
    qpix=hp.ang2pix(nside,qra,qdec,lonlat=True); rpix=hp.ang2pix(nside,rra,rdec,lonlat=True)
    cq=np.bincount(qpix,minlength=npix); cr=np.bincount(rpix,minlength=npix)
    occupied=cr>0; meanr=cr[occupied].mean() if occupied.any() else 1
    comp=cr/meanr
    comp=hp.smoothing(comp,fwhm=np.deg2rad(1),verbose=False)
    c1,c2=(float(chi_of_z(z)) for z in cfg.slabs[0]); width=c2-c1
    omega=np.count_nonzero(comp>=.5)*area
    # Histogram estimate nbar=dN/(dchi dOmega), evaluated per object.
    qc=chi_of_z(qz); edges=np.linspace(c1,c2,41); hist,_=np.histogram(qc,edges)
    ib=np.clip(np.searchsorted(edges,qc,side="right")-1,0,len(hist)-1)
    nbar=np.maximum(hist[ib]/(np.diff(edges)[ib]*max(omega,area)),1e-30)
    wcmb=kernel(qc,float(chi_of_z(Z_CMB)))
    uq=wcmb/(np.asarray(b_q_of_z(qz))*nbar*area)
    qmap=np.bincount(qpix,weights=uq,minlength=npix)
    # Random radial weights use the q nbar estimate and are normalised by Nq/Nr.
    rc=chi_of_z(rz); ir=np.clip(np.searchsorted(edges,rc,side="right")-1,0,len(hist)-1)
    nr=np.maximum(hist[ir]/(np.diff(edges)[ir]*max(omega,area)),1e-30)
    ur=kernel(rc,float(chi_of_z(Z_CMB)))/(np.asarray(b_q_of_z(rz))*nr*area)
    rmap=np.bincount(rpix,weights=ur,minlength=npix)
    ratio=len(qra)/max(len(rra),1)
    mask=(comp>=.5)&(cr>=nmin_rand)
    kmap=np.zeros(npix); kmap[mask]=(qmap[mask]-ratio*rmap[mask])/comp[mask]
    alm=hp.map2alm(kmap,lmax=cfg.lmax_alpha,iter=0)
    shot=(1+ratio)*width*np.mean(kernel(np.linspace(c1,c2,200),float(chi_of_z(Z_CMB)))**2)/max(len(qra)/max(omega,area),1e-30)
    return alm,mask,{"shot_s":float(shot),"data_random_ratio":ratio,"fsky":float(mask.mean()),"completeness":comp}


def gaussian_realisations(Cls,nside,lmax,seed,n=1):
    """Generate correlated scalar Gaussian fields; returns [n,nfield,npix]."""
    cov=np.asarray(Cls,float)
    if cov.ndim==2: cov=cov[None,:,:]
    nf=cov.shape[1]; rng=np.random.default_rng(seed); out=[]
    # synalm's old global RNG API is isolated/restored for reproducibility.
    state=np.random.get_state()
    try:
        for _ in range(n):
            alms=[np.zeros(hp.Alm.getsize(lmax),complex) for __ in range(nf)]
            ell,m=hp.Alm.getlm(lmax)
            for L in range(lmax+1):
                C=cov[min(L,len(cov)-1)]; vals,vec=np.linalg.eigh((C+C.T)/2)
                root=vec@np.diag(np.sqrt(np.maximum(vals,0)))
                ids=np.where(ell==L)[0]
                for idx in ids:
                    z=rng.normal(size=nf) if m[idx]==0 else (rng.normal(size=nf)+1j*rng.normal(size=nf))/np.sqrt(2)
                    v=root@z
                    for f in range(nf): alms[f][idx]=v[f]
            out.append(np.asarray([hp.alm2map(a,nside,verbose=False) for a in alms]))
    finally: np.random.set_state(state)
    return np.asarray(out)

