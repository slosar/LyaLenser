"""Paper figures of the correlation-function fits that set the response kernels (Section IV): for the forest
auto-correlation and the quasar-forest cross-correlation, one figure each with (left) the measured cells against
the fitted model as a function of r_perp in a few r_par strips, (middle) the measured correlation on the
(r_perp, r_par) plane, (right) the residual in units of the per-cell error. Cells are collapsed over the redshift
bins with the pair weights and the model is the same collapse of the redshift-evolving table (the layer at each
cell's mean distance). Reads $LYALENSER_DATA/stageb/dr1_lowz_v7d/xi.h5 (with the cached z-resolved cells) and
dr1_qso_v1d/xi_qf.h5; writes Paper/figures/xi_ff_fit.pdf and xi_qf_fit.pdf (and copies to report/figures)."""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np
import h5py
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator
ROOT=Path(__file__).resolve().parents[1]
import sys; from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # repository root: `lyalenser` imports without installation
from lyalenser.paths import DATA
plt.rcParams.update({'font.size':9,'axes.grid':True,'grid.alpha':.25,'legend.frameon':False})
C=['#0072B2','#E69F00','#009E73']
OUTS=[ROOT/'Paper/figures',ROOT/'report/figures']


def layered_collapse(f,group,num,den,chisum,nr,nz2,rz_c):
    """den-weighted collapse over the z bins of the layered table evaluated at each cell's mean distance."""
    g=f[group]; arr=g['xi'][()]; nodes=g['chi_nodes'][()]; rp=g['r_perp'][()]; rz=g['r_par'][()]
    RP,RZ=np.meshgrid(np.arange(nr)+.5,rz_c,indexing='ij'); out=np.zeros((num.shape[0],nr,nz2))
    for k in range(num.shape[0]):
        d=den[k]; cm=np.divide(chisum[k],d,out=np.full_like(d,np.nan),where=d>0); cm=np.where(np.isfinite(cm),cm,np.nanmean(cm))
        t=np.clip((cm-nodes[0])/(nodes[1]-nodes[0]),0,len(nodes)-1); i=np.minimum(t.astype(int),len(nodes)-2); fr=t-i
        for j in np.unique(np.r_[i.ravel(),i.ravel()+1]):
            L=RegularGridInterpolator((rp,rz),arr[j],bounds_error=False,fill_value=np.nan)((RP,RZ))
            out[k]+=(np.where(i==j,1-fr,0)+np.where(i+1==j,fr,0))*L
    W=den.sum(axis=0); return np.divide((out*den).sum(axis=0),W,out=np.full((nr,nz2),np.nan),where=W>0)


def panel(fig,axes,rp_c,rz_c,raw,err,model,strips,ok,ylab,title,scale_r2,label_model):
    for (lo,hi),col in zip(strips,C):
        m=(rz_c>=lo)&(rz_c<hi); y=raw[:,m].mean(axis=1); e=np.sqrt((err[:,m]**2).sum(axis=1))/m.sum(); mc=np.nanmean(model[:,m],axis=1)
        sc=rp_c**2 if scale_r2 else 1.
        axes[0].errorbar(rp_c[ok],(sc*y)[ok],(sc*e)[ok],fmt='o',ms=3.5,color=col,lw=1,capsize=0,label=f'$r_\\parallel\\in[{lo},{hi})$')
        axes[0].plot(rp_c[ok],(sc*mc)[ok],'-',color=col,lw=1.6)
    axes[0].plot([],[],'-',color='0.3',lw=1.6,label=label_model); axes[0].legend(fontsize=7,ncol=2,loc='best')
    axes[0].set(xlabel=r'$r_\perp$ ($h^{-1}$Mpc)',ylabel=ylab,title=title,xlim=(0,31))
    ext=(0,rp_c[-1]+.5,rz_c[0]-.5,rz_c[-1]+.5); shown=np.where(ok[:,None],raw,np.nan)
    vmax=np.nanpercentile(np.abs(shown),98)
    im=axes[1].imshow(shown.T,origin='lower',extent=ext,aspect='auto',cmap='RdBu_r',vmin=-vmax,vmax=vmax); axes[1].set(xlabel=r'$r_\perp$ ($h^{-1}$Mpc)',ylabel=r'$r_\parallel$ ($h^{-1}$Mpc)',title='measured'); axes[1].grid(False)
    fig.colorbar(im,ax=axes[1],shrink=.85)
    res=np.divide(raw-model,err,out=np.full_like(raw,np.nan),where=np.isfinite(err)&np.isfinite(model)); res[~ok]=np.nan
    im=axes[2].imshow(res.T,origin='lower',extent=ext,aspect='auto',cmap='RdBu_r',vmin=-4,vmax=4); axes[2].set(xlabel=r'$r_\perp$ ($h^{-1}$Mpc)',ylabel=r'$r_\parallel$ ($h^{-1}$Mpc)',title='(measured $-$ model) / error'); axes[2].grid(False)
    fig.colorbar(im,ax=axes[2],shrink=.85)


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--auto',type=Path,default=DATA/'stageb/dr1_lowz_v7d'); ap.add_argument('--cross',type=Path,default=DATA/'stageb/dr1_qso_v1d'); a=ap.parse_args()
    # ---- forest x forest: coarse 1 Mpc/h cells [40, 40] over r_perp, r_par >= 0; fit range 3 <= r_perp < 30, r_par < 30
    with h5py.File(a.auto/'xi.h5') as f:
        num=f['cells/num'][()]; den=f['cells/den'][()]; cs=f['cells/chisum'][()]; fit=json.loads(f['xi'].attrs['meta'])['fit']
        nr,nz2=num.shape[1],num.shape[2]; rz_c=np.arange(nz2)+.5; rp_c=np.arange(nr)+.5
        model=layered_collapse(f,'xi',num,den,cs,nr,nz2,rz_c)
    W=den.sum(axis=0); raw=np.divide(num.sum(axis=0),W,out=np.zeros_like(W),where=W>0); err=np.divide(1.,np.sqrt(W),out=np.full_like(W,np.inf),where=W>0)
    ok=(rp_c>=3)&(rp_c<30); keep=slice(0,30); 
    fig,axes=plt.subplots(1,3,figsize=(12.5,3.7),gridspec_kw={'width_ratios':[1.3,1,1]})
    panel(fig,axes,rp_c[keep],rz_c[keep],raw[keep,keep],err[keep,keep],model[keep,keep],((0,2),(4,6),(10,12)),ok[keep],r'$r_\perp^2\,\xi_F$',
          f"forest $\\times$ forest ($\\chi^2={fit['chi2']:.0f}$ for {fit['cells']} cells)",True,'fitted model (Eq. xifit)')
    fig.tight_layout()
    for o in OUTS: fig.savefig(o/'xi_ff_fit.pdf')
    plt.close(fig); print('wrote xi_ff_fit.pdf')
    # ---- quasar x forest: cells [nz, 40, 80] over r_perp and SIGNED r_par (pixel minus quasar); fit range 3 <= r_perp < 30, |r_par| < 30
    log=json.loads((a.cross/'dr1_qso.json').read_text()); par=log['xi_fit']
    with h5py.File(a.cross/'xi_qf.h5') as f:
        num=f['cells/num'][()]; den=f['cells/den'][()]; use=f['cells/use'][()]; mod=f['cells/model_cells'][()]
    nz,nr,nz2=num.shape; M=np.full(num.shape,np.nan); M[use]=mod; rz_c=np.arange(nz2)-nz2/2+.5; rp_c=np.arange(nr)+.5
    W=den.sum(axis=0); raw=np.divide(num.sum(axis=0),W,out=np.zeros_like(W),where=W>0); err=np.divide(1.,np.sqrt(W),out=np.full_like(W,np.inf),where=W>0)
    Mw=np.where(use,M,0.); model=np.divide((Mw*den).sum(axis=0),(den*use).sum(axis=0),out=np.full_like(W,np.nan),where=(den*use).sum(axis=0)>0)
    okp=(rp_c>=3)&(rp_c<30); okz=np.abs(rz_c)<30; kp=okp|(rp_c<30); kz=okz
    fig,axes=plt.subplots(1,3,figsize=(12.5,3.7),gridspec_kw={'width_ratios':[1.3,1,1]})
    panel(fig,axes,rp_c[kp],rz_c[kz],raw[kp][:,kz],err[kp][:,kz],model[kp][:,kz],((-1,1),(4,6),(10,12)),okp[kp],r'$\xi_{qF}$',
          f"quasar $\\times$ forest ($\\chi^2={par['chi2']:.0f}$ for {par['cells']} cells)",False,'fitted model (Eq. xiqf)')
    fig.tight_layout()
    for o in OUTS: fig.savefig(o/'xi_qf_fit.pdf')
    plt.close(fig); print('wrote xi_qf_fit.pdf')


if __name__=='__main__': main()
