"""Deflection templates on the sphere: phi from kappa, alpha = grad phi at the sightline positions, the cosine-tapered
multipole bands, the curl partner, and the kernel-weighted density map of a tracer catalogue (`matched_template`)."""
from __future__ import annotations

from dataclasses import dataclass, field
import numpy as np
import healpy as hp

from lyalenser.cosmo import chi as chi_of_z
from lyalenser.lensing import kernel, Z_CMB

# Science bands of the deflection template (iteration 9: extended from 300 to 500).  The tracer auto-spectra
# follow linear theory to ell ~ 600 (report/lowz.tex, tracer spectra figure), so the 300-400 and 400-500 bands carry signal; each
# band is fitted independently and collapsed to one amplitude with the curl and junk components marginalised, so
# a band that turns out to be noise costs nothing but its own error.
SCIENCE_BANDS = ((40, 100), (100, 200), (200, 300), (300, 400), (400, 500))


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


def curl(alpha):
    a=np.asarray(alpha)
    return np.column_stack((-a[:,1],a[:,0])).astype(np.float32)


def sphere_band_templates(kappa_alm,ra,dec,nside=1024,science_bands=SCIENCE_BANDS,taper=10,source="alm",
                          transfer=None,lmax=None):
    """Spherical analogue of flat_sky_band_templates: science, curl-partner and junk deflection templates at (ra, dec)
    from a convergence alm (already Wiener-filtered if ``transfer`` is None; otherwise ``transfer[ell]`` is applied).
    The junk template is the complement of the science windows up to the alm's lmax."""
    alm=np.asarray(kappa_alm,complex); lm=hp.Alm.getlmax(len(alm)) if lmax is None else int(lmax)
    ell=np.arange(lm+1); h=np.ones(lm+1) if transfer is None else np.asarray(transfer,float)[:lm+1]
    filters={f"L{lo}_{hi}":cosine_band(ell,lo,hi,taper) for lo,hi in science_bands}
    filters["junk"]=1.-sum(filters.values())
    pos=np.column_stack((np.asarray(ra,float),np.asarray(dec,float)))
    signals=[]; curls=[]
    for lo,hi in science_bands:
        name=f"L{lo}_{hi}"; phi=phi_from_kappa(hp.almxfl(alm,filters[name]*h),lm); alpha=alpha_at(pos,phi,nside)
        signals.append(Template(alpha,name,"signal",phi_lm=phi,Lmin=lo,Lmax=hi,filter=f"cosine taper {taper}",source=source))
        curls.append(Template(curl(alpha),name+"_curl","curl",phi_lm=phi,Lmin=lo,Lmax=hi,filter=f"90-degree rotation; cosine taper {taper}",source=source))
    phi=phi_from_kappa(hp.almxfl(alm,filters["junk"]*h),lm); junk=alpha_at(pos,phi,nside)
    return signals+curls+[Template(junk,"junk","junk",phi_lm=phi,Lmin=0,Lmax=lm,filter="complement of science windows, including taper wings",source=source)],filters


def _catalog_columns(cat):
    if isinstance(cat,dict): return (np.asarray(cat[k]) for k in ("ra","dec","z"))
    return np.asarray(cat.ra),np.asarray(cat.dec),np.asarray(cat.zq)


def matched_template(quasars,randoms,b_q_of_z,cfg,nmin_rand=1,footprint_mask=None,source_chi=None,radial_bins=40,
                     data_weights=None,random_weights=None,nside=None,lmax=None,comp_fwhm_deg=1.):
    """Kernel-matched (CMB source plane, default) or forest-source (``source_chi``) tracer overdensity on the sphere,
    its mask, white shot-noise model and completeness. Optional catalogue weights (DESI WEIGHT columns) and the
    number of radial bins of the nbar(chi) estimate (1 = uniform over the range: thin low-z slices)."""
    qra,qdec,qz=_catalog_columns(quasars); rra,rdec,rz=_catalog_columns(randoms)
    nside=cfg.nside_alpha if nside is None else int(nside); lmax=cfg.lmax_alpha if lmax is None else int(lmax)
    npix=hp.nside2npix(nside); area=4*np.pi/npix
    wd=np.ones(len(qra)) if data_weights is None else np.asarray(data_weights,float)
    wr=np.ones(len(rra)) if random_weights is None else np.asarray(random_weights,float)
    cs=float(chi_of_z(Z_CMB)) if source_chi is None else float(source_chi)
    qpix=hp.ang2pix(nside,qra,qdec,lonlat=True); rpix=hp.ang2pix(nside,rra,rdec,lonlat=True)
    cq=np.bincount(qpix,minlength=npix); cr=np.bincount(rpix,weights=wr,minlength=npix)
    if footprint_mask is None:
        raise ValueError("spherical matched template requires an independently defined footprint_mask")
    footprint_mask=np.asarray(footprint_mask,bool)
    meanr=cr[footprint_mask].mean()
    smooth=hp.smoothing(cr.astype(float),fwhm=np.deg2rad(comp_fwhm_deg))
    support=hp.smoothing(footprint_mask.astype(float),fwhm=np.deg2rad(comp_fwhm_deg))
    comp=np.divide(smooth,meanr*support,out=np.zeros(npix),where=support>1e-6)
    comp/=max(comp[footprint_mask].mean(),1e-30)
    qc=chi_of_z(qz); rc=chi_of_z(rz)
    c1=min(float(np.min(qc)),float(np.min(rc))); c2=max(float(np.max(qc)),float(np.max(rc)))
    omega=np.count_nonzero(footprint_mask)*area
    # Weighted histogram estimate nbar=dN/(dchi dOmega), evaluated per object.
    edges=np.linspace(c1,c2,int(radial_bins)+1); hist,_=np.histogram(qc,edges,weights=wd)
    ib=np.clip(np.searchsorted(edges,qc,side="right")-1,0,len(hist)-1)
    nbar=np.maximum(hist[ib]/(np.diff(edges)[ib]*max(omega,area)),1e-30)
    uq=kernel(qc,cs)/(np.asarray(b_q_of_z(qz))*nbar*area)
    qmap=np.bincount(qpix,weights=uq*wd,minlength=npix)
    ir=np.clip(np.searchsorted(edges,rc,side="right")-1,0,len(hist)-1)
    nr=np.maximum(hist[ir]/(np.diff(edges)[ir]*max(omega,area)),1e-30)
    ur=kernel(rc,cs)/(np.asarray(b_q_of_z(rz))*nr*area)
    rmap=np.bincount(rpix,weights=ur*wr,minlength=npix)
    ratio=float(wd.sum())/max(float(wr.sum()),1e-30)
    mask=footprint_mask&(comp>=.5)&(np.bincount(rpix,minlength=npix)>=nmin_rand)
    kmap=np.zeros(npix); kmap[mask]=(qmap[mask]-ratio*rmap[mask])/comp[mask]
    alm=hp.map2alm(kmap,lmax=lmax,iter=0)
    # White shot noise of the weighted map: N = Var(pixel) x Omega_pix = (sum_obj (u w)^2 / N_pix) x Omega_pix
    # = sum (u w)^2 x Omega_pix^2 / Omega_mask, objects plus randoms (scaled by the data/random ratio).
    shot=float((np.sum((uq*wd)**2)+ratio**2*np.sum((ur*wr)**2))*area*area/max(omega,area))
    return alm,mask,{"shot_s":shot,"shot_s_uniform_model":matched_shot_noise(edges,hist/np.diff(edges)/omega,b_q_of_z,ratio,source_chi=cs),
                     "data_random_ratio":ratio,"fsky":float(mask.mean()),"completeness":comp,"kappa_map":kmap,"nside":nside,"lmax":lmax,
                     "radial_edges":edges,"nbar_chi":hist/np.diff(edges)/omega,"omega_sr":omega}


