"""Gradient/curl, Wiener, Gaussian and matched-density templates."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import sys
import numpy as np
import healpy as hp
from scipy.ndimage import gaussian_filter

HERE=Path(__file__).resolve().parent; CODE=HERE.parent
if str(CODE) not in sys.path: sys.path.insert(0,str(CODE))
from cosmo import chi as chi_of_z
from cross_spectrum import kernel, Z_CMB

# Science bands of the deflection template (iteration 9: extended from 300 to 500).  The tracer auto-spectra
# follow linear theory to ell ~ 600 (report/lowz Figure 5), so the 300-400 and 400-500 bands carry signal; each
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


def flat_sky_band_filters(shape,pixel_size_rad,
                          science_bands=SCIENCE_BANDS,taper=10):
    """Disjoint science filters plus the mandatory outside-band junk filter."""
    nx,ny=shape
    lx=2*np.pi*np.fft.fftfreq(nx,pixel_size_rad)
    ly=2*np.pi*np.fft.rfftfreq(ny,pixel_size_rad)
    ell=np.hypot(lx[:,None],ly[None,:])
    filters={f"L{lo}_{hi}":cosine_band(ell,lo,hi,taper) for lo,hi in science_bands}
    lo=min(x[0] for x in science_bands); hi=max(x[1] for x in science_bands)
    filters["junk"]=1.0-sum(filters.values())
    return ell,filters


def flat_sky_wiener_transfer(shape,pixel_size_rad,L,S_L,C_XX_L):
    """Evaluate h_L=S_L/C_XX_L on the Fourier grid of a flat patch."""
    nx,ny=shape
    lx=2*np.pi*np.fft.fftfreq(nx,pixel_size_rad)
    ly=2*np.pi*np.fft.rfftfreq(ny,pixel_size_rad)
    ell=np.hypot(lx[:,None],ly[None,:])
    s=np.interp(ell,np.asarray(L,float),np.asarray(S_L,float),left=0,right=0)
    c=np.interp(ell,np.asarray(L,float),np.asarray(C_XX_L,float),left=np.inf,right=np.inf)
    return np.divide(s,c,out=np.zeros_like(s),where=c>0)


def _flat_alpha_grid(kappa,pixel_size_rad,filt=None):
    nx,ny=kappa.shape; fk=np.fft.rfft2(kappa)
    if filt is not None: fk=fk*np.asarray(filt)
    lx=2*np.pi*np.fft.fftfreq(nx,pixel_size_rad)
    ly=2*np.pi*np.fft.rfftfreq(ny,pixel_size_rad)
    l2=lx[:,None]**2+ly[None,:]**2
    phi=np.divide(2*fk,l2,out=np.zeros_like(fk),where=l2>0)
    east=np.fft.irfft2(1j*lx[:,None]*phi,s=(nx,ny)).real
    north=np.fft.irfft2(1j*ly[None,:]*phi,s=(nx,ny)).real
    return east,north


def _flat_interp(grid,xrad,yrad,pixel_size_rad):
    nx,ny=grid.shape
    ix=(np.asarray(xrad)/pixel_size_rad+nx/2)%nx
    iy=(np.asarray(yrad)/pixel_size_rad+ny/2)%ny
    i0=np.floor(ix).astype(int); j0=np.floor(iy).astype(int)
    tx=ix-i0; ty=iy-j0; i1=(i0+1)%nx; j1=(j0+1)%ny
    return ((1-tx)*(1-ty)*grid[i0,j0]+tx*(1-ty)*grid[i1,j0]+
            (1-tx)*ty*grid[i0,j1]+tx*ty*grid[i1,j1])


def flat_sky_band_templates(kappa,ra,dec,pixel_size_rad,center=(180.,30.),
                            science_bands=SCIENCE_BANDS,taper=10,
                            source="flat-sky map",transfer=None):
    """Real science, curl-partner and junk templates sampled at a catalogue."""
    ell,filters=flat_sky_band_filters(np.shape(kappa),pixel_size_rad,science_bands,taper)
    x=np.deg2rad(np.asarray(ra)-center[0])*np.cos(np.deg2rad(center[1]))
    y=np.deg2rad(np.asarray(dec)-center[1])
    signals=[]; curls=[]
    for lo,hi in science_bands:
        name=f"L{lo}_{hi}"; filt=filters[name] if transfer is None else filters[name]*np.asarray(transfer)
        ae,an=_flat_alpha_grid(kappa,pixel_size_rad,filt)
        alpha=np.column_stack((_flat_interp(ae,x,y,pixel_size_rad),
                               _flat_interp(an,x,y,pixel_size_rad))).astype(np.float32)
        signals.append(Template(alpha,name,"signal",phi_lm=flat_potential(kappa,pixel_size_rad,filt),Lmin=lo,Lmax=hi,
                                filter=f"cosine taper {taper}",source=source))
        curls.append(Template(curl(alpha),name+"_curl","curl",phi_lm=flat_potential(kappa,pixel_size_rad,filt),Lmin=lo,Lmax=hi,
                              filter=f"90-degree rotation; cosine taper {taper}",source=source))
    jf=filters["junk"] if transfer is None else filters["junk"]*np.asarray(transfer)
    ae,an=_flat_alpha_grid(kappa,pixel_size_rad,jf)
    junk=np.column_stack((_flat_interp(ae,x,y,pixel_size_rad),
                          _flat_interp(an,x,y,pixel_size_rad))).astype(np.float32)
    return signals+curls+[Template(junk,"junk","junk",phi_lm=flat_potential(kappa,pixel_size_rad,jf),Lmin=0,Lmax=int(np.ceil(ell.max())),
                                    filter="complement of science windows, including taper wings",source=source)],filters


def matched_template_flat(quasars,randoms,b_q_of_z,shape,pixel_size_rad,
                          center=(180.,30.),nmin_rand=1,radial_range=None,source_chi=None,radial_bins=40):
    """Flat-sky implementation of the matched catalogue template. ``source_chi`` selects the lensing kernel's
    source plane: the CMB (default, the kernel-matched quasar template) or the forest (chi_ref; iteration 7
    low-redshift tracers, whose weighted map estimates the tracer slice's contribution to kappa_lya)."""
    cs=float(chi_of_z(Z_CMB)) if source_chi is None else float(source_chi)
    qra,qdec,qz=_catalog_columns(quasars); rra,rdec,rz=_catalog_columns(randoms)
    nx,ny=shape; omega_pix=pixel_size_rad**2; footprint=nx*ny*omega_pix
    def indices(ra,dec):
        xr=np.deg2rad(ra-center[0])*np.cos(np.deg2rad(center[1])); yr=np.deg2rad(dec-center[1])
        # Pixel i is centred at (i - n/2) pixels, the convention of _flat_interp / flat_sky_band_templates and of the
        # Fourier-resampled maps; nearest-pixel binning is floor(x + 1/2). floor(x) (used until 2026-09-15) shifted
        # every catalogue map by half a pixel, suppressing its cross-correlation with the field by J0(l pix / 2)
        # (~3 % averaged over 40 <= L <= 300 at 128 pixels on 20 degrees).
        ix=np.floor(xr/pixel_size_rad+nx/2+.5).astype(int); iy=np.floor(yr/pixel_size_rad+ny/2+.5).astype(int)
        keep=(ix>=0)&(ix<nx)&(iy>=0)&(iy<ny)
        return ix,iy,keep
    qi,qj,qkeep=indices(qra,qdec); ri,rj,rkeep=indices(rra,rdec)
    qc=np.asarray(chi_of_z(qz)); rc=np.asarray(chi_of_z(rz))
    c1=min(float(qc.min()),float(rc.min())); c2=max(float(qc.max()),float(rc.max()))
    if radial_range is not None: c1,c2=map(float,radial_range)
    # radial_bins=1: uniform nbar over the range (thin low-z slices, iteration 7; the data path uses n(z) files).
    edges=np.linspace(c1,c2,int(radial_bins)+1); hist,_=np.histogram(qc[qkeep],edges)
    ib=np.clip(np.searchsorted(edges,qc,side="right")-1,0,len(hist)-1)
    ir=np.clip(np.searchsorted(edges,rc,side="right")-1,0,len(hist)-1)
    nbarq=np.maximum(hist[ib]/(np.diff(edges)[ib]*footprint),1e-30)
    nbarr=np.maximum(hist[ir]/(np.diff(edges)[ir]*footprint),1e-30)
    uq=kernel(qc,cs)/(np.asarray(b_q_of_z(qz))*nbarq*omega_pix)
    ur=kernel(rc,cs)/(np.asarray(b_q_of_z(rz))*nbarr*omega_pix)
    qmap=np.zeros(shape); rmap=np.zeros(shape); counts=np.zeros(shape)
    np.add.at(qmap,(qi[qkeep],qj[qkeep]),uq[qkeep])
    np.add.at(rmap,(ri[rkeep],rj[rkeep]),ur[rkeep])
    np.add.at(counts,(ri[rkeep],rj[rkeep]),1)
    sigma=np.deg2rad(1)/(2.355*pixel_size_rad)
    comp=gaussian_filter(counts,sigma=sigma,mode="wrap")
    comp/=max(comp.mean(),1e-30)
    ratio=qkeep.sum()/max(rkeep.sum(),1); mask=(comp>=.5)&(counts>=nmin_rand)
    out=np.zeros(shape); out[mask]=(qmap[mask]-ratio*rmap[mask])/comp[mask]
    return out.astype(np.float32),mask,{"data_random_ratio":float(ratio),"completeness":comp,
                                       "template_chi_range":[c1,c2],"fsky":float(mask.mean()),
                                       "shot_s":matched_shot_noise(edges,hist/np.diff(edges)/footprint,b_q_of_z,ratio,source_chi=cs),
                                       "radial_edges":edges,"nbar_chi":hist/np.diff(edges)/footprint}


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


def flat_potential(kappa,pixel_size_rad,transfer=None):
    nx,ny=kappa.shape
    lx=2*np.pi*np.fft.fftfreq(nx,pixel_size_rad)
    ly=2*np.pi*np.fft.rfftfreq(ny,pixel_size_rad)
    l2=lx[:,None]**2+ly[None,:]**2
    f=np.fft.rfft2(kappa)
    if transfer is not None: f*=transfer
    return np.divide(2*f,l2,out=np.zeros_like(f),where=l2>0)


def matched_shot_noise(edges,nbar_chi,bias,ratio,source_chi=None):
    """(1+r) integral W^2/(b^2 nbar_3D chi^2) dchi; nbar_chi=n3D chi^2. ``source_chi``: kernel source plane (CMB default)."""
    from cosmo import z_of_chi
    cs=float(chi_of_z(Z_CMB)) if source_chi is None else float(source_chi)
    # Gauss-Legendre integrates each radial histogram cell without boundary ambiguity.
    x,w=np.polynomial.legendre.leggauss(8)
    c=.5*(edges[1:]+edges[:-1])[:,None]+.5*np.diff(edges)[:,None]*x
    b=np.asarray(bias(z_of_chi(c.ravel()).reshape(c.shape)))
    integrand=kernel(c.ravel(),cs).reshape(c.shape)**2/(b*b*np.maximum(np.asarray(nbar_chi)[:,None],1e-30))
    return float((1+ratio)*np.sum(.5*np.diff(edges)*np.sum(integrand*w,axis=1)))


def map_spectrum(map_a,pixel_size_rad,edges,mask=None,map_b=None):
    """Binned pseudo spectrum of maps already passed through the common operator."""
    a=np.asarray(map_a); b=a if map_b is None else np.asarray(map_b)
    f=np.fft.rfft2(a); g=np.fft.rfft2(b)
    ell,_=flat_sky_band_filters(a.shape,pixel_size_rad)
    multiplicity=np.full(ell.shape,2.); multiplicity[:,0]=1
    if a.shape[1]%2==0: multiplicity[:,-1]=1
    fsky=1. if mask is None else np.mean(np.asarray(mask)**2)
    power=(f*g.conj()).real*pixel_size_rad**2/a.size/max(fsky,1e-30)
    count=np.histogram(ell,bins=edges,weights=multiplicity)[0]
    num=np.histogram(ell,bins=edges,weights=power*multiplicity)[0]
    return np.divide(num,count,out=np.zeros_like(num),where=count>0)
