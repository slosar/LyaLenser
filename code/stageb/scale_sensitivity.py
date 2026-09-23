"""Scale sensitivity of the lensing estimator on DR1 (paper Section IV.D): the Fisher-information density of the
amplitude on the (r_perp, r_par) plane, F(cell) = sum_{pairs in cell} w_p w_q G^2 d^2 with G = chi_mid dxi/dr_perp
the response kernel at the pair's mean distance (layered table) and d the pair scalar of the combined science-band
deflection template (sum of the fiducial bands), i.e. the contribution of each 1 Mpc/h cell to the response matrix
element of the science amplitude. The same for the quasar-forest pairs (w_p, signed r_par). The mock version of
this (code/pipeline/signal_profile.py) used a flat table and truth templates; this one uses the production pair
catalogues, weights, cuts and kernels. Writes report/stageb/scale_sensitivity_<tag>.json and
Paper/figures/scale_sensitivity.pdf (+ report/lowz/figures).
Usage: python scale_sensitivity.py --auto $LYALENSER_DATA/stageb/dr1_lowz_v7d --cross $LYALENSER_DATA/stageb/dr1_qso_v1d --lowz $LYALENSER_DATA/lowz_v4 --bands 40 200 400 600 800 1000
"""
from __future__ import annotations
import argparse, json, sys, time
from pathlib import Path
import numpy as np
import h5py, healpy as hp
from numba import njit, prange
HERE=Path(__file__).resolve().parent; CODE=HERE.parent; ROOT=CODE.parent
for p in (CODE,CODE/'pipeline',HERE):
    if str(p) not in sys.path: sys.path.insert(0,str(p))
from paths import DATA
from mock import load_sightlines
from campaign4 import campaign_config, read_xi
from templates import SCIENCE_BANDS, sphere_band_templates
from pairs import _interp_layer, PairCatalogue
from amplitude import pair_scalars
from qso_io import read_quasars
from xi_cross import Positions
from cosmo import chi as chi_of_z


def load_cat(path,group='all'):
    with h5py.File(path) as f:
        g=f[group]; return PairCatalogue(g['a'][()],g['b'][()],g['thx'][()],g['thy'][()],g['theta'][()],g['accum'][()].astype(np.float64),g['npair'][()],{k:g.attrs[k] for k in g.attrs})


@njit(parallel=True)
def _ff(pix_start,chi,weight,slab,region,pa,pb,theta,dsq,rp_grid,rz_grid,xi,xirp,rpmax,rzmax,rpmin,chi0,dchi,nchi,nb,nblock):
    n=pa.size; ff=np.zeros((nblock,nb,nb)); cnt=np.zeros((nblock,nb,nb)); bs=(n+nblock-1)//nblock
    rp0=rp_grid[0]; rz0=rz_grid[0]; drp=rp_grid[1]-rp_grid[0]; drz=rz_grid[1]-rz_grid[0]; nrp=rp_grid.size; nrz=rz_grid.size
    for blk in prange(nblock):
        for ip in range(blk*bs,min((blk+1)*bs,n)):
            a=pa[ip]; b=pb[ip]; qb=pix_start[b]; qend=pix_start[b+1]; d2=dsq[ip]
            for p in range(pix_start[a],pix_start[a+1]):
                cp=np.float64(chi[p])
                while qb<qend and np.float64(chi[qb])<cp-rzmax: qb+=1
                q=qb
                while q<qend and np.float64(chi[q])<=cp+rzmax:
                    cq=np.float64(chi[q]); rz=abs(cp-cq); cm=.5*(cp+cq); rp=cm*np.float64(theta[ip])
                    sel=slab[p]>=0 and slab[q]>=0 and not (region[p]>0 and region[q]>0)
                    if sel and rp<=rpmax and rp>=rpmin and rz<rzmax:
                        xv,xg=_interp_layer(rp,rz,cm,rp0,drp,nrp,rz0,drz,nrz,chi0,dchi,nchi,xi,xirp)
                        G=cm*xg; ww=np.float64(weight[p])*np.float64(weight[q]); i=min(int(rp),nb-1); j=min(int(rz),nb-1)
                        ff[blk,i,j]+=ww*G*G*d2; cnt[blk,i,j]+=1.
                    q+=1
    return ff.sum(axis=0),cnt.sum(axis=0)


@njit(parallel=True)
def _qf(pix_start,chi,weight,slab,chiq,pa,pb,theta,nsl,dsq,rp_grid,rz_grid,xi,xirp,rpmax,rzmax,rpmin,chi0,dchi,nchi,nb,nblock):
    n=pa.size; ff=np.zeros((nblock,nb,2*nb)); cnt=np.zeros((nblock,nb,2*nb)); bs=(n+nblock-1)//nblock
    rp0=rp_grid[0]; rz0=rz_grid[0]; drp=rp_grid[1]-rp_grid[0]; drz=rz_grid[1]-rz_grid[0]; nrp=rp_grid.size; nrz=rz_grid.size
    for blk in prange(nblock):
        for ip in range(blk*bs,min((blk+1)*bs,n)):
            cq=np.float64(chiq[pa[ip]-nsl]); s=pb[ip]; d2=dsq[ip]
            for p in range(pix_start[s],pix_start[s+1]):
                cp=np.float64(chi[p]); rz=cp-cq
                if rz<-rzmax or rz>rzmax: continue
                cm=.5*(cp+cq); rp=cm*np.float64(theta[ip])
                if slab[p]>=0 and rp<=rpmax and rp>=rpmin:
                    xv,xg=_interp_layer(rp,rz,cm,rp0,drp,nrp,rz0,drz,nrz,chi0,dchi,nchi,xi,xirp)
                    G=cm*xg; ww=np.float64(weight[p]); i=min(int(rp),nb-1); j=min(int(rz+nb),2*nb-1)
                    ff[blk,i,j]+=ww*G*G*d2; cnt[blk,i,j]+=1.
    return ff.sum(axis=0),cnt.sum(axis=0)


def science_dsq(cat,alm,ra,dec,nside,bands):
    ts,_=sphere_band_templates(alm,ra,dec,nside=nside,science_bands=bands,source='combined')
    sci=[t for t in ts if getattr(t,'kind','')=='signal']; d,_,_=pair_scalars(cat,sci); return d.sum(axis=0)**2


def summarise(ff,rp_edges,rz_edges):
    tot=ff.sum(); frp=ff.sum(axis=1)/tot; frz=ff.sum(axis=0)/tot; cum=np.cumsum(frp)
    return {'information_2d':ff.tolist(),'fraction_vs_rperp':frp.tolist(),'fraction_vs_rpar':frz.tolist(),
            'rperp_at_25_50_75_percent':[float(np.searchsorted(cum,x)+1) for x in (.25,.5,.75)],
            'fraction_rpar_below_5':float(ff[:,(np.abs(rz_edges[:-1]+.5)<5)].sum()/tot),'fraction_rpar_below_10':float(ff[:,(np.abs(rz_edges[:-1]+.5)<10)].sum()/tot),
            'fraction_rperp_below_10':float(frp[:10].sum()),'fraction_rperp_below_20':float(frp[:20].sum())}


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--auto',type=Path,default=DATA/'stageb/dr1_lowz_v7d'); ap.add_argument('--cross',type=Path,default=DATA/'stageb/dr1_qso_v1d')
    ap.add_argument('--lowz',type=Path,default=DATA/'lowz_v4'); ap.add_argument('--bands',type=float,nargs='+',default=None); ap.add_argument('--nside-alpha',type=int,default=2048); ap.add_argument('--tag',default='v4')
    a=ap.parse_args(); t0=time.perf_counter(); BANDS=tuple((int(a.bands[i]),int(a.bands[i+1])) for i in range(len(a.bands)-1)) if a.bands else SCIENCE_BANDS
    alm=hp.read_alm(str(a.lowz/'kappa_combined_alm.fits')); nb=30; out={'bands':[list(b) for b in BANDS],'cell_mpc':1.0}
    sl=load_sightlines(a.auto/'sightlines.h5'); alog=json.loads((a.auto/'dr1_lowz.json').read_text()); cfg=campaign_config(1.).copy(chi_ref=float(alog['chi_ref']))
    region=sl.region if getattr(sl,'region',None) is not None else np.zeros(len(sl.chi),np.int8)
    # ---- forest x forest
    cat=load_cat(a.auto/'catalogue.h5'); print(f'[ff] {len(cat.a)} pairs ({time.perf_counter()-t0:.0f} s)',flush=True)
    dsq=science_dsq(cat,alm,sl.ra,sl.dec,a.nside_alpha,BANDS); print(f'[ff] template evaluated ({time.perf_counter()-t0:.0f} s)',flush=True)
    T=read_xi(a.auto/'xi.h5','xi'); chi0,dchi,nchi=T.layers()
    ff,cnt=_ff(sl.pix_start,sl.chi,sl.w,sl.slab,region,cat.a,cat.b,cat.theta,dsq,T.r_perp,T.r_par,T.xi.ravel(),T.xi_rp.ravel(),cfg.r_perp_max,cfg.r_par_max,float(getattr(cfg,'r_perp_min',0.)),chi0,dchi,nchi,nb,1024)
    out['ff']=summarise(ff,np.arange(nb+1),np.arange(nb+1)); out['ff']['pixel_pairs']=float(cnt.sum()); print(f"[ff] done: r_perp at 25/50/75 % {out['ff']['rperp_at_25_50_75_percent']}, |r_par|<5: {out['ff']['fraction_rpar_below_5']:.2f} ({time.perf_counter()-t0:.0f} s)",flush=True)
    # ---- quasar x forest
    clog=json.loads((a.cross/'dr1_qso.json').read_text()); qso=read_quasars(*clog['quasar_z'],verbose=False); pos=Positions(sl,qso)
    catq=load_cat(a.cross/'catalogue.h5'); print(f'[qf] {len(catq.a)} pairs ({time.perf_counter()-t0:.0f} s)',flush=True)
    dsq=science_dsq(catq,alm,pos.ra,pos.dec,a.nside_alpha,BANDS); Tq=read_xi(a.cross/'xi_qf.h5','xi'); chi0,dchi,nchi=Tq.layers()
    fq,cq=_qf(sl.pix_start,sl.chi,sl.w,sl.slab,np.asarray(qso.chi,np.float32),catq.a,catq.b,catq.theta,sl.nq,dsq,Tq.r_perp,Tq.r_par,Tq.xi.ravel(),Tq.xi_rp.ravel(),cfg.r_perp_max,cfg.r_par_max,float(getattr(cfg,'r_perp_min',0.)),chi0,dchi,nchi,nb,1024)
    out['qf']=summarise(fq,np.arange(nb+1),np.arange(-nb,nb+1)); out['qf']['pixel_pairs']=float(cq.sum()); print(f"[qf] done: r_perp at 25/50/75 % {out['qf']['rperp_at_25_50_75_percent']}, |r_par|<5: {out['qf']['fraction_rpar_below_5']:.2f} ({time.perf_counter()-t0:.0f} s)",flush=True)
    (ROOT/'report/stageb'/f'scale_sensitivity_{a.tag}.json').write_text(json.dumps(out,indent=1)+'\n')
    plot(out)


def plot(out):
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':9,'axes.grid':False,'legend.frameon':False})
    ff=np.array(out['ff']['information_2d']); fq=np.array(out['qf']['information_2d']); ff/=ff.sum(); fq/=fq.sum()
    fig,ax=plt.subplots(1,2,figsize=(8.6,3.6))
    im=ax[0].imshow(1e3*ff.T,origin='lower',extent=(0,30,0,30),cmap='viridis'); ax[0].set(xlabel=r'$r_\perp$ ($h^{-1}$Mpc)',ylabel=r'$r_\parallel$ ($h^{-1}$Mpc)',title='forest $\\times$ forest'); fig.colorbar(im,ax=ax[0],shrink=.85,label=r'$10^3\times$ fraction of the information per cell')
    im=ax[1].imshow(1e3*fq.T,origin='lower',extent=(0,30,-30,30),aspect='auto',cmap='viridis'); ax[1].set(xlabel=r'$r_\perp$ ($h^{-1}$Mpc)',ylabel=r'$r_\parallel$ (pixel minus quasar, $h^{-1}$Mpc)',title='quasar $\\times$ forest'); fig.colorbar(im,ax=ax[1],shrink=.85,label=r'$10^3\times$ fraction of the information per cell')
    fig.tight_layout()
    for o in (ROOT/'Paper/figures',ROOT/'report/lowz/figures'):
        if o.exists(): fig.savefig(o/'scale_sensitivity.pdf')
    plt.close(fig); print('wrote scale_sensitivity.pdf')


if __name__=='__main__': main()
