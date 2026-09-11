"""Fixed-observed-geometry flat-sky Stage-A mocks.

The implementation uses the specified 2 Mpc/h transverse and 0.5 Mpc/h
radial FFT grid.  ``scale`` shrinks the angular patch only, making validation
and unit-test runs inexpensive while retaining the full forest depth.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
import argparse,gc,json,sys,time
import numpy as np
from scipy.fft import rfftn,irfftn,fftfreq,rfftfreq
from scipy.ndimage import gaussian_filter, distance_transform_edt
from numba import njit, prange
import h5py

HERE=Path(__file__).resolve().parent; CODE=HERE.parent
if str(CODE) not in sys.path: sys.path.insert(0,str(CODE))
from cosmo import chi as chi_of_z,z_of_chi,linear_pk_interp
from cross_spectrum import kernel,Z_CMB
from forest_power import FID
try:
    from .config import Config,SightlineSet,kernel_product_g1
except ImportError:
    from config import Config,SightlineSet,kernel_product_g1


@dataclass
class MockResult:
    sightlines: SightlineSet
    alpha_lya: np.ndarray
    truth: dict
    quasars: dict
    randoms: dict
    maps: dict
    attrs: dict


def _fft_fields(nx,ny,nz,dx,dz,seed,workers=4):
    """Matter, forest, long mode and radial RSD displacement.

    A unit-variance white field filtered by sqrt(P / cell_volume) has the
    continuum FFT normalization int d^3k P/(2pi)^3.  In particular, delta_m
    is never rescaled by its realization standard deviation: it is the
    physical 2 Mpc/h density requested by the mock specification.
    """
    def progress(label):
        import resource
        print(f"mock FFT {label}: peak_RSS={resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2:.2f} GB",flush=True)
    rng=np.random.default_rng(seed)
    white=rng.normal(size=(nx,ny,nz)).astype(np.float32)
    progress("white")
    fk=rfftn(white,workers=workers,overwrite_x=True)
    del white
    progress("forward")
    kx=2*np.pi*fftfreq(nx,dx); ky=2*np.pi*fftfreq(ny,dx); kz=2*np.pi*rfftfreq(nz,dz)
    kt=np.sqrt(kx[:,None]**2+ky[None,:]**2)
    pk=linear_pk_interp(zmax=4,kmax=50,nonlinear=False)
    kg=np.geomspace(1e-4,49,3000); pg=pk.P(2.4,kg)
    cell_volume=dx*dx*dz
    for iz,kzi in enumerate(kz):
        kk=np.sqrt(kt*kt+kzi*kzi); pl=np.interp(np.clip(kk,kg[0],kg[-1]),kg,pg)
        filt=np.sqrt(pl/cell_volume); filt[0,0]=0 if iz==0 else filt[0,0]
        fk[:,:,iz]*=filt
    # Derive each secondary field sequentially from the one matter spectrum;
    # avoiding three simultaneous complex copies saves ~8 GB at scale=1.
    ff=fk.copy()
    progress("forest spectrum copy")
    for iz,kzi in enumerate(kz):
        kk=np.sqrt(kt*kt+kzi*kzi); mu2=kzi*kzi/np.maximum(kk*kk,1e-30)
        ff[:,:,iz]*=FID["b_F"]*(1+FID["beta_F"]*mu2)*np.exp(-.5*(kk/FID["kp"])**2)
    df=irfftn(ff,s=(nx,ny,nz),workers=workers,overwrite_x=True).astype(np.float32,copy=False)
    del ff
    progress("forest inverse")
    dm=irfftn(fk,s=(nx,ny,nz),workers=workers,overwrite_x=False).astype(np.float32,copy=False)
    progress("matter inverse")
    # Solve on the ORIGINAL Fourier grid: no density decimation or aliasing.
    for iz,kzi in enumerate(kz):
        k2=kt*kt+kzi*kzi
        fk[:,:,iz]*=1j*.97*kzi*np.divide(1.,k2,out=np.zeros_like(k2),where=k2>0)
    rsd=irfftn(fk,s=(nx,ny,nz),workers=workers,overwrite_x=True).astype(np.float32,copy=False)
    del fk
    stride=1
    progress("full RSD inverse")
    long=gaussian_filter(dm,sigma=(5.0,5.0,20.0),mode="wrap",output=np.float32)
    progress("long smoothing")
    return dm,df,long,rsd,stride


def completeness_pattern(nx,ny):
    """Smooth, positive angular selection pattern with a 30 percent depth."""
    x=2*np.pi*(np.arange(nx)+.5)/nx
    y=2*np.pi*(np.arange(ny)+.5)/ny
    return (1.0-.3*np.sin(x[:,None])**2*np.sin(y[None,:])**2).astype(np.float32)


def sample_lognormal_quasars(delta_g,dx,dz,chi0,area_deg2,nbar_deg2,
                              rng,b_q=3.5,angular_weight=None,
                              magnification_map=None,magnification=0.5):
    """Poisson-sample the physical 3D lognormal quasar intensity.

    ``nbar_deg2`` is the mean surface density over this radial volume.  The
    hierarchical radial/transverse draw avoids materializing a second copy of
    a potentially multi-billion-cell field.
    """
    nx,ny,nz=delta_g.shape
    # Streaming moments avoid a second full-volume float64 subtraction array.
    total=total2=0.
    for iz in range(0,nz,16):
        block=np.asarray(delta_g[:,:,iz:iz+16],dtype=np.float64)
        total+=float(block.sum()); total2+=float(np.sum(block*block))
    del block
    variance=max(0.,total2/delta_g.size-(total/delta_g.size)**2)
    correction=.5*b_q*b_q*variance
    aw=np.ones((nx,ny),np.float32) if angular_weight is None else np.asarray(angular_weight,np.float32)
    if magnification_map is not None:
        aw=aw*np.clip(1.0+magnification*np.asarray(magnification_map,np.float32),0.05,None)
    zsum=np.empty(nz,np.float64)
    for iz in range(nz):
        w=np.exp(np.clip(b_q*delta_g[:,:,iz]-correction,-30,30))*aw
        zsum[iz]=np.sum(w,dtype=np.float64)
    intensity=area_deg2*nbar_deg2*zsum.sum()/delta_g.size
    nobj=int(rng.poisson(intensity))
    zdraw=rng.choice(nz,size=nobj,p=zsum/zsum.sum())
    ix=np.empty(nobj,np.int32); iy=np.empty(nobj,np.int32)
    for iz in np.unique(zdraw):
        take=np.flatnonzero(zdraw==iz)
        w=np.exp(np.clip(b_q*delta_g[:,:,iz]-correction,-30,30))*aw
        cell=rng.choice(nx*ny,size=len(take),replace=True,p=(w/w.sum()).ravel())
        ix[take]=cell//ny; iy[take]=cell%ny
    x=(ix+rng.random(nobj)-nx/2)*dx
    y=(iy+rng.random(nobj)-ny/2)*dx
    chi=chi0+(zdraw+rng.random(nobj))*dz
    return {"x":x,"y":y,"chi":chi,"ix":ix,"iy":iy,"iz":zdraw,
            "density_variance":variance,"lognormal_correction":correction,"integrated_intensity":float(intensity)}


def _sample_angular_selection(n,rng,side,weight):
    """Uniform radial/angular randoms accepted by a smooth completeness."""
    nx,ny=weight.shape; xs=[]; ys=[]
    while sum(map(len,xs)) < n:
        batch=max(1024,int(1.5*(n-sum(map(len,xs)))))
        x=rng.uniform(-side/2,side/2,batch); y=rng.uniform(-side/2,side/2,batch)
        ix=np.clip(((x/side+.5)*nx).astype(int),0,nx-1)
        iy=np.clip(((y/side+.5)*ny).astype(int),0,ny-1)
        keep=rng.random(batch)<weight[ix,iy]
        xs.append(x[keep]); ys.append(y[keep])
    return np.concatenate(xs)[:n],np.concatenate(ys)[:n]


@lru_cache(maxsize=4)
def _outside_limber_grid(cmin,cmax,cref,lmax):
    """Three restricted Limber spectra for convergence outside the box."""
    from three_tracer import limber
    ccmb=float(chi_of_z(Z_CMB))
    L=np.unique(np.r_[2.,np.linspace(2,max(2,lmax),640)])
    wlya=lambda c: kernel(c,cref)
    wcmb=lambda c: kernel(c,ccmb)
    kw=dict(nchi=500,to_recombination=False)
    ll=limber(L,wlya,wlya,1.0,cmin,**kw)
    lc=limber(L,wlya,wcmb,1.0,cmin,**kw)
    cc=limber(L,wcmb,wcmb,1.0,cmin,**kw)
    if cmax < ccmb:
        cc+=limber(L,wcmb,wcmb,cmax,ccmb,nchi=500,to_recombination=True)
    return L,ll,lc,cc


def _correlated_outside_pair(nx,ny,side_angle,cmin,cmax,cref,seed):
    """Flat-sky Gaussian (kappa_lya_rest, kappa_CMB_rest) from Limber C_L."""
    lx=2*np.pi*np.fft.fftfreq(nx,side_angle/nx)
    ly=2*np.pi*np.fft.rfftfreq(ny,side_angle/ny)
    ell=np.hypot(lx[:,None],ly[None,:]); lmax=int(np.ceil(ell.max()))
    L,c11g,c12g,c22g=_outside_limber_grid(float(cmin),float(cmax),float(cref),lmax)
    c11=np.interp(ell,L,c11g,left=0,right=c11g[-1])
    c12=np.interp(ell,L,c12g,left=0,right=c12g[-1])
    c22=np.interp(ell,L,c22g,left=0,right=c22g[-1])
    rng=np.random.default_rng(seed)
    z1=np.fft.rfft2(rng.normal(size=(nx,ny)).astype(np.float32))
    z2=np.fft.rfft2(rng.normal(size=(nx,ny)).astype(np.float32))
    pixarea=(side_angle/nx)*(side_angle/ny)
    a=np.sqrt(np.maximum(c11,0)/pixarea)
    b=np.divide(c12,np.sqrt(np.maximum(c11,0)*pixarea),out=np.zeros_like(c12),where=c11>0)
    cres=np.maximum(c22-np.divide(c12*c12,c11,out=np.zeros_like(c12),where=c11>0),0)
    f1=z1*a; f2=z1*b+z2*np.sqrt(cres/pixarea)
    f1[0,0]=0; f2[0,0]=0
    return (np.fft.irfft2(f1,s=(nx,ny)).astype(np.float32),
            np.fft.irfft2(f2,s=(nx,ny)).astype(np.float32),
            {"L":L,"klkl_rest":c11g,"klkc_rest":c12g,"kckc_rest":c22g})


def _act_mask_cutout(nx,ny,scale):
    """Actual ACT DR6 mask cutout at (180,10), rotated onto the mock patch."""
    import fitsio,healpy as hp
    path=str(__import__("pathlib").Path(__import__("os").environ.get("LYALENSER_DATA","/data/LyaLenser"))/"raw/act/baseline/mask_act_dr6_lensing_v1_healpix_nside_4096_baseline.fits")
    off=10.0*scale
    ra=180+np.linspace(-off,off,nx,endpoint=False)+off/nx
    dec=10+np.linspace(-off,off,ny,endpoint=False)+off/ny
    rr,dd=np.meshgrid(ra,dec,indexing="ij")
    pix=hp.ang2pix(4096,rr.ravel(),dd.ravel(),lonlat=True)
    rows=np.unique(pix//1024)
    values=fitsio.FITS(path)[1].read(rows=rows,columns=["T"])["T"]
    out=values[np.searchsorted(rows,pix//1024),pix%1024].reshape(nx,ny)
    return np.asarray(out,np.float32),path


def _deflection(kappa,dx,chi_ref):
    nx,ny=kappa.shape; f=np.fft.rfft2(kappa)
    lx=2*np.pi*np.fft.fftfreq(nx,dx)*chi_ref
    ly=2*np.pi*np.fft.rfftfreq(ny,dx)*chi_ref
    L2=lx[:,None]**2+ly[None,:]**2
    phi=np.divide(2*f,L2,out=np.zeros_like(f),where=L2>0)
    ae=np.fft.irfft2(1j*lx[:,None]*phi,s=kappa.shape).real
    an=np.fft.irfft2(1j*ly[None,:]*phi,s=kappa.shape).real
    return ae.astype(np.float32),an.astype(np.float32)


def _interp2(grid,x,y,dx):
    nx,ny=grid.shape
    ix=(x/dx+nx/2)%nx; iy=(y/dx+ny/2)%ny
    i0=np.floor(ix).astype(int); j0=np.floor(iy).astype(int); tx=ix-i0; ty=iy-j0
    i1=(i0+1)%nx; j1=(j0+1)%ny
    return ((1-tx)*(1-ty)*grid[i0,j0]+tx*(1-ty)*grid[i1,j0]+(1-tx)*ty*grid[i0,j1]+tx*ty*grid[i1,j1])


def alpha_from_map(mock,map_name):
    """Flat-sky deflection of one saved convergence map at the sightlines."""
    dx=float(mock.attrs["dx"]); cref=float(mock.sightlines.attrs["chi_ref"])
    ae,an=_deflection(mock.maps[map_name],dx,cref)
    x=np.deg2rad(mock.sightlines.ra-180)*np.cos(np.deg2rad(30))*cref
    y=np.deg2rad(mock.sightlines.dec-30)*cref
    return np.column_stack((_interp2(ae,x,y,dx),_interp2(an,x,y,dx))).astype(np.float32)


def flat_sky_power(a,b,side_angle,Lmin=40,Lmax=300):
    """Mean flat-sky cross spectrum over an annular multipole band."""
    nx,ny=a.shape; pixarea=(side_angle/nx)*(side_angle/ny); area=side_angle**2
    fa=np.fft.rfft2(np.asarray(a))*pixarea; fb=np.fft.rfft2(np.asarray(b))*pixarea
    lx=2*np.pi*np.fft.fftfreq(nx,side_angle/nx); ly=2*np.pi*np.fft.rfftfreq(ny,side_angle/ny)
    ell=np.hypot(lx[:,None],ly[None,:]); use=(ell>Lmin)&(ell<Lmax)
    return float(np.mean((fa*np.conj(fb)).real[use]/area)),float(np.mean(ell[use])),int(use.sum())


@lru_cache(maxsize=8)
def _spectrum_theory(z1,z2,margin,Lmin,Lmax):
    from three_tracer import spectra
    cf=np.array([float(chi_of_z(z1)),float(chi_of_z(z2))]); Lgrid=np.linspace(Lmin,Lmax,261)
    theory=spectra(z1=z1,z2=z2,zs=2.4,template_z1=float(z_of_chi(cf[0]-margin)),
                   template_z2=float(z_of_chi(cf[1]+margin)),L_values=Lgrid)
    return Lgrid,theory


def mock_spectrum_check(mock,Lmin=40,Lmax=300):
    """Compare the three signal-map spectra with three_tracer.spectra."""
    side_angle=np.deg2rad(20*float(mock.attrs["scale"])); cbox=np.asarray(mock.attrs.get("template_chi",[]))
    kl=mock.maps["kappa_lya"]; kc=mock.maps["kappa_CMB_signal"]
    ll,Leff,nmode=flat_sky_power(kl,kl,side_angle,Lmin,Lmax)
    lc,_,_=flat_sky_power(kl,kc,side_angle,Lmin,Lmax)
    cc,_,_=flat_sky_power(kc,kc,side_angle,Lmin,Lmax)
    z1,z2=mock.sightlines.attrs["zmin"],mock.sightlines.attrs["zmax"]
    margin=float(mock.attrs["template_margin"])
    Lgrid,theory=_spectrum_theory(float(z1),float(z2),margin,float(Lmin),float(Lmax))
    nx,ny=kl.shape
    lx=2*np.pi*np.fft.fftfreq(nx,side_angle/nx); ly=2*np.pi*np.fft.rfftfreq(ny,side_angle/ny)
    ell=np.hypot(lx[:,None],ly[None,:]); modes=ell[(ell>Lmin)&(ell<Lmax)]
    expected=np.array([np.mean(np.interp(modes,Lgrid,theory[k])) for k in ("klkl","klkc","kckc")])
    measured=np.array([ll,lc,cc])
    return {"L_effective":Leff,"n_modes":nmode,"measured":measured,
            "expected":expected,"ratio":measured/expected}


def sightlines_for_variant(mock,A_true,response):
    """Return a sightline view for a precomputed fixed-geometry variant."""
    label=("%g"%float(A_true)).replace("-","m").replace(".","p")
    key=f"delta_A{label}_response{int(bool(response))}"
    if key not in mock.truth: raise KeyError(f"mock does not contain {key}")
    s=mock.sightlines
    return SightlineSet(s.qid,s.ra,s.dec,s.zq,s.pix_start,s.chi,mock.truth[key],s.w,s.slab,
                        {**s.attrs,"A_true":float(A_true),"response":bool(response)})


@njit(parallel=True,cache=True)
def _interp3_chunked(grid,points,chunk=1_000_000):
    """Periodic trilinear interpolation, including the last-to-first cell."""
    out=np.empty(len(points),np.float32)
    nx,ny,nz=grid.shape
    for p in prange(len(points)):
        x,y,z=points[p,0]%nx,points[p,1]%ny,points[p,2]%nz
        i,j,k=int(x),int(y),int(z); tx,ty,tz=x-i,y-j,z-k
        value=0.
        for di in range(2):
            for dj in range(2):
                for dk in range(2):
                    value+=grid[(i+di)%nx,(j+dj)%ny,(k+dk)%nz]*(tx if di else 1-tx)*(ty if dj else 1-ty)*(tz if dk else 1-tz)
        out[p]=value
    return out


def project_continuum(values,chi,weights):
    """Weighted joint intercept/slope projection along the final axis."""
    values=np.asarray(values,float); weights=np.broadcast_to(weights,values.shape)
    chi=np.broadcast_to(chi,values.shape)
    norm=weights.sum(axis=-1,keepdims=True)
    if np.any(norm<=0): raise ValueError("empty continuum weights")
    offset=chi-chi[..., :1]
    x=offset-(weights*offset).sum(axis=-1,keepdims=True)/norm
    mean=(weights*values).sum(axis=-1,keepdims=True)/norm
    denom=(weights*x*x).sum(axis=-1,keepdims=True)
    slope=np.divide((weights*x*(values-mean)).sum(axis=-1,keepdims=True),denom,out=np.zeros_like(mean),where=denom>0)
    return values-mean-slope*x


def sky_rays(ra,dec):
    """Unit rays in the Cartesian basis east,north,LOS at (180,30).

    Their local metric is dtheta^2=cos(dec_mid)^2 dra^2+ddec^2,
    identical to pair_geometry to the specified small-pair accuracy.
    """
    r=np.deg2rad(np.asarray(ra)-180); d=np.deg2rad(dec); d0=np.deg2rad(30.)
    return np.stack((np.cos(d)*np.sin(r),np.sin(d)*np.cos(d0)-np.cos(d)*np.cos(r)*np.sin(d0),
                     np.sin(d)*np.sin(d0)+np.cos(d)*np.cos(r)*np.cos(d0)),axis=-1)


def ray_points(ra,dec,chi,dx,dz,shape,cmin):
    v=sky_rays(ra,dec)*np.asarray(chi)[...,None]
    return np.stack((v[...,0]/dx+shape[0]/2,v[...,1]/dx+shape[1]/2,(v[...,2]-cmin)/dz),axis=-1)


def rsd_displacement(delta,dx,dz,f=.97):
    """Reference full Fourier solution used in convergence tests."""
    nx,ny,nz=delta.shape; ft=rfftn(delta)
    kx=2*np.pi*fftfreq(nx,dx); ky=2*np.pi*fftfreq(ny,dx); kz=2*np.pi*rfftfreq(nz,dz)
    k2=kx[:,None,None]**2+ky[None,:,None]**2+kz[None,None,:]**2
    ft*=1j*f*np.divide(kz[None,None,:],k2,out=np.zeros_like(k2),where=k2>0)
    return irfftn(ft,s=delta.shape)


def project_lightcone_fields(dm,long,chis,dx,dz,cref,wc,wl,trap,cforest):
    """Joint density-intensity and convergence operator, exposed for regression tests."""
    nx,ny,nz=dm.shape; cbox=(float(chis[0]),float(chis[-1]))
    # All projected fields and catalogue intensities follow exactly these rays.
    ax=(np.arange(nx)-nx/2)*dx/cref; ay=(np.arange(ny)-ny/2)*dx/cref
    rr,dd=np.meshgrid(180+np.rad2deg(ax)/np.cos(np.deg2rad(30)),30+np.rad2deg(ay),indexing="ij")
    raydm=np.empty_like(dm); kslab=np.zeros((nx,ny),np.float32); klya_box=np.zeros_like(kslab)
    delta_L_map=np.zeros_like(kslab); forest_slices=0
    rays=sky_rays(rr.ravel(),dd.ravel())
    # Radial batches preserve cache locality in the Cartesian cube and avoid
    # thousands of launches of a parallel interpolation kernel.
    for iz0 in range(0,nz,16):
        iz1=min(nz,iz0+16); cs=chis[iz0:iz1]
        v=rays[:,None,:]*cs[None,:,None]
        v[...,0]=v[...,0]/dx+nx/2; v[...,1]=v[...,1]/dx+ny/2; v[...,2]=(v[...,2]-cbox[0])/dz
        points=v.reshape(-1,3)
        planes=_interp3_chunked(dm,points).reshape(nx,ny,len(cs))
        raydm[:,:,iz0:iz1]=planes
        kslab+=planes@np.asarray(wc[iz0:iz1]*trap[iz0:iz1],np.float32)
        klya_box+=planes@np.asarray(wl[iz0:iz1]*trap[iz0:iz1],np.float32)
        inside=(cs>=cforest[0])&(cs<=cforest[1])
        if inside.any():
            longplanes=_interp3_chunked(long,points).reshape(nx,ny,len(cs))
            delta_L_map+=longplanes@inside.astype(np.float32); forest_slices+=int(inside.sum())
        del points,v,planes
    delta_L_map/=forest_slices
    return raydm,kslab,klya_box,delta_L_map


def generate_mock(cfg=None,seed=0,scale=None,A_true=1.0,g_on=True,response=True,
                  magnification=False,completeness=False,real_mask=False,
                  n_los=None,pixel_noise_power=None,template_margin=150.0,
                  disjoint_selection=False,cmb_noise=False,variant_A_values=None,box_margin=300.0):
    cfg=Config() if cfg is None else cfg
    scale=cfg.scale if scale is None else scale; n_los=cfg.n_los if n_los is None else n_los
    pn=cfg.pixel_noise_power if pixel_noise_power is None else pixel_noise_power
    rng=np.random.default_rng(seed); cref=cfg.chi_ref
    forest_z=(min(b[0] for b in cfg.slabs),max(b[1] for b in cfg.slabs))
    cforest=np.array([float(chi_of_z(z)) for z in forest_z])
    cbox=np.array([cforest[0]-box_margin,cforest[1]+box_margin]); dx=2.0; dz=.5
    side=cref*np.deg2rad(20*scale)
    nx=max(16,int(np.ceil(side/dx))); ny=nx; nz=int(np.ceil((cbox[1]-cbox[0])/dz))+1
    stage_start=time.perf_counter(); timings={}
    dm,df,long,rsd,rsd_stride=_fft_fields(nx,ny,nz,dx,dz,seed)
    timings["fft_s"]=time.perf_counter()-stage_start; stage_start=time.perf_counter()
    volume_field_bytes=dm.nbytes
    chis=cbox[0]+np.arange(nz)*dz
    wc=kernel(chis,float(chi_of_z(Z_CMB)))
    wl=kernel(chis,cref)*(chis<cref)
    trap=np.full(nz,dz,np.float32); trap[[0,-1]]*=.5
    raydm,kslab,klya_box,delta_L_map=project_lightcone_fields(dm,long,chis,dx,dz,cref,wc,wl,trap,cforest)
    del dm
    timings["lightcone_projection_s"]=time.perf_counter()-stage_start; stage_start=time.perf_counter()
    rest_lya,rest_cmb,rest_cls=_correlated_outside_pair(nx,ny,side/cref,cbox[0],cbox[1],cref,seed+71003)
    klya=klya_box+rest_lya; kcmb_signal=kslab+rest_cmb; kcmb=kcmb_signal.copy()
    cmb_noise_level=0.0
    if cmb_noise:
        numbers=CODE.parent/"report"/"numbers3.json"
        cmb_noise_level=float(json.loads(numbers.read_text())["Nc"]["ACT_white"])
        pixarea=(side/cref/nx)*(side/cref/ny)
        kcmb+=rng.normal(scale=np.sqrt(cmb_noise_level/pixarea),size=(nx,ny)).astype(np.float32)
    kcmb_unmasked=kcmb.copy()
    mask=np.ones((nx,ny),np.float32); mask_path=""
    if real_mask:
        mask,mask_path=_act_mask_cutout(nx,ny,scale)
        kcmb=kcmb*mask
    ae,an=_deflection(klya,dx,cref)
    area=(20*scale)**2
    comp=completeness_pattern(nx,ny) if completeness else np.ones((nx,ny),np.float32)
    radial_ratio=(cbox[1]-cbox[0])/(cforest[1]-cforest[0])
    nq_target=max(20,int(round(n_los*area)))
    base_density=max(25.0,1.15*n_los)
    if area<4: base_density=max(base_density,1.5*nq_target/area)
    template_density=base_density*radial_ratio
    qc=sample_lognormal_quasars(raydm,dx,dz,cbox[0],area,template_density,rng,b_q=3.5,
                                angular_weight=comp,
                                magnification_map=(klya if magnification else None))
    # Linear-theory radial redshift-space displacement, v_parallel/(aH).
    qra=180+np.rad2deg(qc["x"]/cref)/np.cos(np.deg2rad(30)); qdec=30+np.rad2deg(qc["y"]/cref)
    qchi=qc["chi"]+_interp3_chunked(rsd,ray_points(qra,qdec,qc["chi"],dx,dz,rsd.shape,cbox[0]))
    del rsd
    in_template=(qchi>=cbox[0])&(qchi<=cbox[1])
    for key in ("x","y","chi","ix","iy","iz"): qc[key]=np.asarray(qc[key])[in_template]
    qchi=qchi[in_template]
    eligible=(qchi>=cforest[0])&(qchi<=cforest[1])
    if real_mask:
        # A zero pad makes the boundary of the extracted ACT cutout an edge.
        dist=distance_transform_edt(np.pad(mask>.5,1))[1:-1,1:-1]*(20*scale/nx)
        edge_cut=min(2.0,2.0*scale)
        eligible &= dist[qc["ix"],qc["iy"]] >= edge_cut
    else:
        edge_cut=0.0
    ids=np.flatnonzero(eligible)
    if len(ids)<nq_target and not real_mask:
        raise RuntimeError(f"only {len(ids)} eligible clustered quasars for {nq_target} sightlines")
    if real_mask: nq_target=min(nq_target,len(ids))
    if nq_target<20: raise RuntimeError("real-mask edge rejection leaves fewer than 20 sightlines")
    sight_ids=np.sort(rng.choice(ids,nq_target,replace=False)); nq=len(sight_ids)
    del raydm; gc.collect()
    x=qc["x"][sight_ids]; y=qc["y"][sight_ids]
    alpha=np.column_stack((_interp2(ae,x,y,dx),_interp2(an,x,y,dx))).astype(np.float32)
    ra=180+np.rad2deg(x/cref)/np.cos(np.deg2rad(30)); dec=30+np.rad2deg(y/cref)
    timings["maps_and_catalogue_s"]=time.perf_counter()-stage_start; stage_start=time.perf_counter()
    cpix=np.arange(cforest[0],cforest[1]+.25,.55,dtype=np.float32); npc=len(cpix)
    allchi=np.tile(cpix,nq); xs=np.repeat(x,npc); ys=np.repeat(y,npc)
    aa=np.repeat(alpha[:,0],npc); bb=np.repeat(alpha[:,1],npc)
    # Log-normal per-skewer P_N; weighted continuum mean+slope removal.
    pnsk=pn*np.exp(2*rng.normal(size=nq)) if pn>0 else np.zeros(nq)
    sig=np.sqrt(np.repeat(pnsk,npc)/.55); noise=(rng.normal(size=len(sig))*sig).astype(np.float32)
    varf=float(np.var(df)); w=(1/(sig*sig+varf)).astype(np.float32)
    def continuum(values):
        return project_continuum(np.asarray(values).reshape(nq,npc),cpix,w.reshape(nq,npc)).astype(np.float32).ravel()
    requested=sorted(set([float(A_true)]+([] if variant_A_values is None else [float(a) for a in variant_A_values])))
    # Pixel modulation is evaluated after interpolation, BEFORE continuum.
    # This makes the actual operator exactly P diag(1+R delta_L(pixel)/2).
    variants={}
    g=1+cfg.g1*(allchi-cref) if g_on else np.ones_like(allchi)
    for Aval in requested:
        # Mock lensing samples delta_F(theta_obs + alpha), hence +chi*alpha in x/y.
        # The transverse FFT coordinates are comoving: an angular sightline
        # at x_ref/chi_ref intersects a radial slice at x_ref*chi/chi_ref.
        ral=np.repeat(ra,npc)+np.rad2deg(Aval*aa*g/np.cos(np.deg2rad(np.repeat(dec,npc))))
        decl=np.repeat(dec,npc)+np.rad2deg(Aval*bb*g)
        points=ray_points(ral,decl,allchi,dx,dz,df.shape,cbox[0])
        base_values=_interp3_chunked(df,points)
        modulation=_interp3_chunked(long,points)
        variants[(Aval,False)]=continuum(base_values+noise)
        variants[(Aval,True)]=continuum(base_values*(1+.5*cfg.response_delta*modulation)+noise)
        if Aval==0: pixel_long=modulation.copy()
        del points,base_values
    if 0. not in requested:
        pixel_long=_interp3_chunked(long,ray_points(np.repeat(ra,npc),np.repeat(dec,npc),allchi,dx,dz,df.shape,cbox[0]))
    del df,long
    delta=variants[(float(A_true),bool(response))]
    starts=np.arange(nq+1,dtype=np.int64)*npc
    pixel_slabs=np.full(len(allchi),-1,np.int8)
    for index,bounds in enumerate(cfg.slabs):
        lo,hi=(float(chi_of_z(z)) for z in bounds)
        pixel_slabs[(allchi>=lo)&(allchi<hi)]=index
    sight=SightlineSet(sight_ids,ra,dec,np.full(nq,z_of_chi(cforest[1]),np.float32),starts,
                       allchi,delta,w,pixel_slabs,
                       {"zmin":forest_z[0],"zmax":forest_z[1],"description":"Stage A flat-sky FFT mock",
                        "chi_ref":cref,"A_true":A_true,"g_on":g_on,"response":response})
    keep_q=np.ones(len(qchi),bool)  # retain the largest catalogue for nested-margin selections
    if disjoint_selection: keep_q[sight_ids]=False
    tqx=qc["x"][keep_q]; tqy=qc["y"][keep_q]; tqchi=qchi[keep_q]
    qcat={"ra":180+np.rad2deg(tqx/cref)/np.cos(np.deg2rad(30)),
          "dec":30+np.rad2deg(tqy/cref),"z":z_of_chi(tqchi),"chi":tqchi,"qid":np.flatnonzero(keep_q)}
    nr=20*len(tqx); rw=comp/float(comp.max())
    rx,ry=_sample_angular_selection(nr,rng,side,rw)
    rchi=rng.uniform(cbox[0],cbox[1],nr)
    rcat={"ra":180+np.rad2deg(rx/cref)/np.cos(np.deg2rad(30)),
          "dec":30+np.rad2deg(ry/cref),"z":z_of_chi(rchi),"chi":rchi}
    maps={"kappa_lya":klya,"kappa_lya_box":klya_box,"kappa_lya_rest":rest_lya,
          "kappa_slab":kslab,"kappa_rest":rest_cmb,"kappa_CMB_signal":kcmb_signal,"kappa_CMB":kcmb,
          "kappa_CMB_unmasked":kcmb_unmasked,"alpha_east":ae,"alpha_north":an,"completeness":comp,"act_mask":mask}
    maps["delta_L"]=delta_L_map
    truth={"alpha_lya":alpha,"delta_L_pixel":pixel_long,"quasar_sightline_ids":sight_ids,
           "outside_L":rest_cls["L"],"outside_klkl":rest_cls["klkl_rest"],
           "outside_klkc":rest_cls["klkc_rest"],"outside_kckc":rest_cls["kckc_rest"]}
    if variant_A_values is not None:
        for (Aval,respflag),values in variants.items():
            label=("%g"%Aval).replace("-","m").replace(".","p")
            truth[f"delta_A{label}_response{int(respflag)}"]=values
    timings["skewers_s"]=time.perf_counter()-stage_start
    g1_value,g1_mean_chi=kernel_product_g1()
    attrs={"seed":seed,"scale":scale,"timings":timings,"grid":[nx,ny,nz],"dx":dx,"dz":dz,"single_growth_z":2.4,
           "magnification":magnification,"completeness":completeness,"real_mask":real_mask,
           "template_margin":template_margin,"box_margin":box_margin,"template_density_deg2":int(np.sum((tqchi>=cforest[0]-template_margin)&(tqchi<=cforest[1]+template_margin)))/area,
           "integrated_quasar_intensity":qc["integrated_intensity"],
           "quasar_bias":3.5,"delta_g_variance":qc["density_variance"],"rsd_f":.97,
           "rsd_grid_stride":rsd_stride,
           "cmb_noise":cmb_noise,"cmb_noise_level":cmb_noise_level,
           "variant_A_values":requested,
           "g1":g1_value,"g1_mean_lens_chi":g1_mean_chi,"act_mask_source":mask_path,
           "mask_edge_rejection_deg":edge_cut,
           "lensing_sampling_sign":"theta_obs + alpha",
           "peak_array_bytes_estimate":int(4*volume_field_bytes+
                                            0+
                                            sum(v.nbytes for v in variants.values())+8*nx*ny)}
    gc.collect()
    return MockResult(sight,alpha,truth,qcat,rcat,maps,attrs)


def save_mock(mock,path):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    with h5py.File(path,"w") as f:
        s=f.create_group("sightlines")
        for k in ("qid","ra","dec","zq","pix_start","chi","delta","w","slab"): s[k]=getattr(mock.sightlines,k)
        for k,v in mock.sightlines.attrs.items(): s.attrs[k]=v
        for name,d in (("truth",mock.truth),("quasars",mock.quasars),("randoms",mock.randoms),("maps",mock.maps)):
            g=f.create_group(name)
            for k,v in d.items(): g[k]=v
        for k,v in mock.attrs.items(): f.attrs[k]=json.dumps(v) if isinstance(v,(list,dict)) else v


def load_sightlines(path):
    with h5py.File(path,"r") as f:
        g=f["sightlines"]; vals=[g[k][()] for k in ("qid","ra","dec","zq","pix_start","chi","delta","w","slab")]
        attrs=dict(g.attrs)
    return SightlineSet(*vals,attrs)


if __name__=="__main__":
    ap=argparse.ArgumentParser(); ap.add_argument("--seed",type=int,default=0); ap.add_argument("--scale",type=float,default=.08)
    ap.add_argument("--output",default="/tmp/LyaLenser/mocks/mock.h5"); args=ap.parse_args()
    m=generate_mock(seed=args.seed,scale=args.scale); save_mock(m,args.output); print(args.output,m.attrs)
