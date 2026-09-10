"""Fixed-observed-geometry flat-sky Stage-A mocks.

The implementation uses the specified 2 Mpc/h transverse and 0.5 Mpc/h
radial FFT grid.  ``scale`` shrinks the angular patch only, making validation
and unit-test runs inexpensive while retaining the full forest depth.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import argparse,json,sys
import numpy as np
from scipy.fft import rfftn,irfftn,fftfreq,rfftfreq
from scipy.ndimage import gaussian_filter
from scipy.interpolate import RegularGridInterpolator
import h5py

HERE=Path(__file__).resolve().parent; CODE=HERE.parent
if str(CODE) not in sys.path: sys.path.insert(0,str(CODE))
from cosmo import chi as chi_of_z,z_of_chi,linear_pk_interp
from cross_spectrum import kernel,Z_CMB
from forest_power import FID
try:
    from .config import Config,SightlineSet
except ImportError:
    from config import Config,SightlineSet


@dataclass
class MockResult:
    sightlines: SightlineSet
    alpha_lya: np.ndarray
    truth: dict
    quasars: dict
    randoms: dict
    maps: dict
    attrs: dict


def _fft_fields(nx,ny,nz,dx,dz,seed):
    rng=np.random.default_rng(seed)
    white=rng.normal(size=(nx,ny,nz)).astype(np.float32)
    fk=rfftn(white,workers=24,overwrite_x=True)
    kx=2*np.pi*fftfreq(nx,dx); ky=2*np.pi*fftfreq(ny,dx); kz=2*np.pi*rfftfreq(nz,dz)
    kt=np.sqrt(kx[:,None]**2+ky[None,:]**2)
    pk=linear_pk_interp(zmax=4,kmax=50,nonlinear=False)
    kg=np.geomspace(1e-4,49,3000); pg=pk.P(2.4,kg)
    ff=fk.copy()
    for iz,kzi in enumerate(kz):
        kk=np.sqrt(kt*kt+kzi*kzi); pl=np.interp(np.clip(kk,kg[0],kg[-1]),kg,pg)
        filt=np.sqrt(pl); filt[0,0]=0 if iz==0 else filt[0,0]
        fk[:,:,iz]*=filt
        mu2=kzi*kzi/np.maximum(kk*kk,1e-30)
        ff[:,:,iz]*=filt*FID["b_F"]*(1+FID["beta_F"]*mu2)*np.exp(-.5*(kk/FID["kp"])**2)
    dm=irfftn(fk,s=(nx,ny,nz),workers=24).astype(np.float32)
    df=irfftn(ff,s=(nx,ny,nz),workers=24).astype(np.float32)
    dm/=max(float(dm.std()),1e-12)
    # Finite-box/discrete-grid normalisation chosen to match the continuum
    # Kaiser ForestPower variance after the per-skewer continuum projection.
    df*=0.255/max(float(df.std()),1e-12)
    return dm,df


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


def generate_mock(cfg=None,seed=0,scale=None,A_true=1.0,g_on=True,response=True,
                  magnification=False,completeness=False,real_mask=False,
                  n_los=None,pixel_noise_power=None):
    cfg=Config() if cfg is None else cfg
    scale=cfg.scale if scale is None else scale; n_los=cfg.n_los if n_los is None else n_los
    pn=cfg.pixel_noise_power if pixel_noise_power is None else pixel_noise_power
    rng=np.random.default_rng(seed); cref=cfg.chi_ref
    cforest=np.array([float(chi_of_z(z)) for z in cfg.slabs[0]])
    cbox=np.array([cforest[0]-150,cforest[1]+150]); dx=2.0; dz=.5
    side=cref*np.deg2rad(20*scale)
    nx=max(16,int(np.ceil(side/dx))); ny=nx; nz=int(np.ceil((cbox[1]-cbox[0])/dz))+1
    dm,df=_fft_fields(nx,ny,nz,dx,dz,seed)
    long=gaussian_filter(dm,sigma=(2.5,2.5,20.0),mode="wrap")
    if response: df*=1+.5*cfg.response_delta*long
    chis=cbox[0]+np.arange(nz)*dz
    wc=kernel(chis,float(chi_of_z(Z_CMB)))
    wl=kernel(chis,cref)*(chis<cref)
    kslab=np.trapz(dm*wc[None,None,:],chis,axis=2).astype(np.float32)
    klya=np.trapz(dm*wl[None,None,:],chis,axis=2).astype(np.float32)
    # Correlated outside-box modes: a low-pass pair, independent plus shared.
    rest=gaussian_filter(rng.normal(size=(nx,ny)),3,mode="wrap").astype(np.float32)
    rest*=max(float(klya.std()),1e-7)/max(float(rest.std()),1e-12)
    klya=klya+0.45*rest; kcmb=kslab+0.65*rest
    if real_mask:
        xx=(np.arange(nx)-nx/2)[:,None]/(nx/2); yy=(np.arange(ny)-ny/2)[None,:]/(ny/2)
        mask=(xx*xx+(.8*yy)**2<.85**2); kcmb=kcmb*mask
    ae,an=_deflection(klya,dx,cref)
    area=(20*scale)**2; nq=max(20,int(round(n_los*area)))
    # Angular positions are Poisson samples of a projected lognormal density.
    proj=dm.mean(axis=2); prob=np.exp(3.5*proj/np.std(proj)-.5*3.5**2); prob=prob.ravel(); prob/=prob.sum()
    cell=rng.choice(nx*ny,size=nq,replace=True,p=prob)
    ix=cell//ny; iy=cell%ny
    x=(ix+rng.random(nq)-nx/2)*dx; y=(iy+rng.random(nq)-ny/2)*dx
    alpha=np.column_stack((_interp2(ae,x,y,dx),_interp2(an,x,y,dx))).astype(np.float32)
    ra=180+np.rad2deg(x/cref)/np.cos(np.deg2rad(30)); dec=30+np.rad2deg(y/cref)
    cpix=np.arange(cforest[0],cforest[1]+.25,.55,dtype=np.float32); npc=len(cpix)
    allchi=np.tile(cpix,nq); xs=np.repeat(x,npc); ys=np.repeat(y,npc)
    aa=np.repeat(alpha[:,0],npc); bb=np.repeat(alpha[:,1],npc)
    if g_on: g=1+cfg.g1*(allchi-cref)
    else: g=np.ones_like(allchi)
    # Mock lensing samples delta_F(theta_obs + alpha), hence +chi*alpha in x/y.
    xl=xs+A_true*allchi*aa*g; yl=ys+A_true*allchi*bb*g
    gx=((xl/dx+nx/2)%nx)*dx; gy=((yl/dx+ny/2)%ny)*dx; gz=allchi-cbox[0]
    interp=RegularGridInterpolator((np.arange(nx)*dx,np.arange(ny)*dx,np.arange(nz)*dz),df,
                                   bounds_error=False,fill_value=None)
    delta=interp(np.column_stack((gx,gy,gz))).astype(np.float32)
    # Log-normal per-skewer P_N; weighted continuum mean+slope removal.
    pnsk=pn*np.exp(2*rng.normal(size=nq)) if pn>0 else np.zeros(nq)
    sig=np.sqrt(np.repeat(pnsk,npc)/.55); delta+=rng.normal(size=len(delta))*sig
    varf=float(np.var(df)); w=(1/(sig*sig+varf)).astype(np.float32)
    for i in range(nq):
        sl=slice(i*npc,(i+1)*npc); cc=cpix-cref; ww=w[sl].astype(float); yy=delta[sl].astype(float)
        X=np.column_stack((np.ones(npc),cc)); mat=X.T@(ww[:,None]*X)
        coef=np.linalg.solve(mat,X.T@(ww*yy)); delta[sl]-=(X@coef).astype(np.float32)
    starts=np.arange(nq+1,dtype=np.int64)*npc
    sight=SightlineSet(np.arange(nq),ra,dec,np.full(nq,z_of_chi(cforest[1]),np.float32),starts,
                       allchi,delta,w,np.zeros(len(delta),np.int8),
                       {"zmin":cfg.slabs[0][0],"zmax":cfg.slabs[0][1],"description":"Stage A flat-sky FFT mock",
                        "chi_ref":cref,"A_true":A_true,"g_on":g_on,"response":response})
    nt=max(20,int(round(25*area))); tqx=rng.uniform(-side/2,side/2,nt); tqy=rng.uniform(-side/2,side/2,nt)
    tqz=rng.uniform(cfg.slabs[0][0],cfg.slabs[0][1],nt)
    qcat={"ra":180+np.rad2deg(tqx/cref)/np.cos(np.deg2rad(30)),"dec":30+np.rad2deg(tqy/cref),"z":tqz}
    nr=20*nt; rcat={"ra":rng.uniform(ra.min(),ra.max(),nr),"dec":rng.uniform(dec.min(),dec.max(),nr),
                    "z":rng.uniform(cfg.slabs[0][0],cfg.slabs[0][1],nr)}
    maps={"kappa_lya":klya,"kappa_slab":kslab,"kappa_CMB":kcmb,"alpha_east":ae,"alpha_north":an}
    truth={"alpha_lya":alpha,"delta_L":np.asarray([_interp2(long[:,:,nz//2],xx,yy,dx) for xx,yy in zip(x,y)])}
    attrs={"seed":seed,"scale":scale,"grid":[nx,ny,nz],"dx":dx,"dz":dz,"single_growth_z":2.4,
           "magnification":magnification,"completeness":completeness,"real_mask":real_mask,
           "lensing_sampling_sign":"theta_obs + alpha","peak_array_bytes_estimate":int((dm.nbytes+df.nbytes)*4)}
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
