"""QA plots of the quasar x forest correlation fit (iteration 12): the measured cells against the model in r_par
profiles at several r_perp (offset and smoothing visible), the residual map, the per-redshift-bin amplitude and
the fitted b_q(z) against the DESI relation, and the kernel. Reads --run/xi_qf.h5 and dr1_qso.json; writes
report/lowz/figures/qso_xi_{profiles,residuals,evolution}.pdf.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np
import h5py
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent; CODE=HERE.parent; ROOT=CODE.parent
for p in (CODE,CODE/'pipeline',HERE):
    if str(p) not in sys.path: sys.path.insert(0,str(p))
from paths import DATA
from cosmo import z_of_chi, chi as chi_of_z
OUT=ROOT/'report/lowz/figures'
plt.rcParams.update({'font.size':9,'axes.grid':True,'grid.alpha':.25,'legend.frameon':False})


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--run',type=Path,default=DATA/'stageb/dr1_qso_v1'); a=ap.parse_args()
    log=json.loads((a.run/'dr1_qso.json').read_text()); par=log['xi_fit']
    with h5py.File(a.run/'xi_qf.h5') as f:
        num=f['cells/num'][()]; den=f['cells/den'][()]; cs=f['cells/chisum'][()]; use=f['cells/use'][()]; model=f['cells/model_cells'][()]; zed=f['cells/z_edges'][()]
        xi=f['xi/xi'][()]; xirp=f['xi/xi_rp'][()]; rp=f['xi/r_perp'][()]; rz=f['xi/r_par'][()]; nodes=f['xi/chi_nodes'][()]
    nz,nr,nz2=num.shape; raw=np.divide(num,den,out=np.zeros_like(num),where=den>0); err=np.divide(1.,np.sqrt(den),out=np.full_like(den,np.inf),where=den>0)
    M=np.zeros_like(raw); M[use]=model; rz_c=np.arange(nz2)-nz2/2+.5; rp_c=np.arange(nr)+.5
    # collapse over z (weighted)
    W=den.sum(axis=0); R=np.divide(num.sum(axis=0),W,out=np.zeros_like(W),where=W>0); Mz=np.divide((M*den).sum(axis=0),W,out=np.zeros_like(W),where=W>0); E=np.divide(1.,np.sqrt(W),out=np.full_like(W,np.inf),where=W>0)
    # 1. r_par profiles at several r_perp
    fig,axes=plt.subplots(2,3,figsize=(11,6),sharex=True); sel=[3,5,8,12,18,25]
    for ax,i in zip(axes.ravel(),sel):
        ax.errorbar(rz_c,R[i],E[i],fmt='o',ms=3,color='0.3',label='measured' if i==sel[0] else None); ax.plot(rz_c,Mz[i],color='#D55E00',lw=1.5,label='model' if i==sel[0] else None)
        ax.axvline(par['dr_par'],color='#0072B2',ls='--',lw=1,label=f"$\\Delta r_\\parallel={par['dr_par']:.2f}$" if i==sel[0] else None); ax.axhline(0,color='k',lw=.5)
        ax.set_title(f'$r_\\perp\\in[{i},{i+1})$ Mpc/h',fontsize=9)
    for ax in axes[1]: ax.set_xlabel(r'$r_\parallel$ (pixel minus quasar, $h^{-1}$Mpc)')
    for ax in axes[:,0]: ax.set_ylabel(r'$\xi_{qF}$')
    axes[0,0].legend(fontsize=8); fig.suptitle(f"quasar x forest correlation, cells collapsed over redshift; fit: $b_q={par['b_q']:.2f}$, $\\gamma_q={par['gamma_q']:.2f}$, $\\sigma_\\parallel={par['sigma_par']:.2f}$, $\\chi^2={par['chi2']:.0f}/{par['cells']}$",fontsize=9)
    fig.tight_layout(); fig.savefig(OUT/'qso_xi_profiles.pdf'); plt.close(fig)
    # 2. residual maps: data, model, (data-model)/err, per z bin chi2
    fig,axes=plt.subplots(1,3,figsize=(12,3.6)); ext=[rz_c[0]-.5,rz_c[-1]+.5,0,nr]
    vmax=np.nanpercentile(np.abs(R[use.any(axis=0)]),98)
    for ax,arr,t in zip(axes[:2],(R,Mz),('measured','model')): im=ax.imshow(arr,origin='lower',extent=ext,aspect='auto',cmap='RdBu_r',vmin=-vmax,vmax=vmax); ax.set_title(t); ax.set_xlabel(r'$r_\parallel$'); ax.set_ylabel(r'$r_\perp$')
    fig.colorbar(im,ax=axes[:2],shrink=.8,label=r'$\xi_{qF}$')
    res=np.divide(R-Mz,E,out=np.zeros_like(R),where=np.isfinite(E)); res[~use.any(axis=0)]=np.nan
    im=axes[2].imshow(res,origin='lower',extent=ext,aspect='auto',cmap='RdBu_r',vmin=-4,vmax=4); axes[2].set_title('(measured - model) / error, fit range'); axes[2].set_xlabel(r'$r_\parallel$'); fig.colorbar(im,ax=axes[2],shrink=.8)
    fig.savefig(OUT/'qso_xi_residuals.pdf',bbox_inches='tight'); plt.close(fig)
    # 3. evolution: per z bin chi2 and amplitude ratio (measured / model summed over the fit range), b_q(z) vs DESI relation
    fig,(ax1,ax2)=plt.subplots(1,2,figsize=(9,3.4))
    zm=[b['mean_z'] for b in par['per_z']]; chi=[b['chi2']/max(b['cells'],1) for b in par['per_z']]
    amp=[float(np.sum((num[k]*M[k])[use[k]])/np.sum((M[k]**2*den[k])[use[k]])) for k in range(nz)]
    ax1.plot(zm,amp,'o-',color='#D55E00',label='measured / model amplitude per z bin'); ax1.axhline(1,color='k',lw=.8); ax1.set(xlabel='pair mean redshift',ylabel='amplitude ratio'); ax1.legend(fontsize=8)
    ax1b=ax1.twinx(); ax1b.plot(zm,chi,'s--',color='0.4',label=r'$\chi^2$ per cell'); ax1b.set_ylabel(r'$\chi^2$ per cell'); ax1b.legend(fontsize=8,loc='lower right')
    zz=np.linspace(1.9,3.1,50); bq=par['b_q']*((1+zz)/(1+par['z_ref']))**par['gamma_q']; desi=0.278*((1+zz)**2-6.565)+2.393
    ax2.plot(zz,bq,color='#D55E00',lw=2,label=f"fit: $b_q={par['b_q']:.2f}\\,x^{{{par['gamma_q']:.2f}}}$"); ax2.plot(zz,desi,color='#0072B2',ls='--',label='DESI QSO relation $0.278[(1+z)^2-6.565]+2.393$')
    ax2.set(xlabel='$z$',ylabel='$b_q$'); ax2.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(OUT/'qso_xi_evolution.pdf'); plt.close(fig)
    print('wrote',OUT/'qso_xi_profiles.pdf',OUT/'qso_xi_residuals.pdf',OUT/'qso_xi_evolution.pdf')


if __name__=='__main__': main()
