"""Figure: the combined template against the ACT and Planck CMB lensing maps. Left: the convergence template x
kappa_CMB (spin-0 bandpowers); middle: the deflection template (science window 40 <= L <= 1000, E-mode of the
spin-1 field the estimator consumes) x kappa_CMB; right: the ratio of the deflection cross-spectrum to its prediction.
Reads report/stageb/deflection_cmb_check_<tag>.json (default v4); writes report/lowz/figures/template_cmb_cross.pdf
and Paper/figures/template_cmb_cross.pdf when that directory exists."""
import argparse, json
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2]
ap=argparse.ArgumentParser(); ap.add_argument('--tag',default='v4'); a=ap.parse_args()
d=json.load(open(ROOT/'report/stageb'/f'deflection_cmb_check_{a.tag}.json'))
C={'ACT':'#D55E00','Planck':'#0072B2'}; OFF={'ACT':-5,'Planck':5}
plt.rcParams.update({'font.size':9,'axes.grid':True,'grid.alpha':.25,'legend.frameon':False})
fig,ax=plt.subplots(1,3,figsize=(12,3.6))
for s in ('ACT','Planck'):
    c=d['surveys'][s]['templates']['combined']; ls='-' if s=='ACT' else '--'
    # left: convergence template x kappa_CMB
    k=c['kappa_spin0']; L=np.array(k['L']); cx=np.array(k['cross']); e=np.array(k['cross_err']); th=np.array(k['prediction']); m=L>=40
    ax[0].errorbar(L[m]+OFF[s],(L*cx)[m],(L*e)[m],fmt='o',ms=4,color=C[s],capsize=2,label=f"{s}: $A={k['A']:.2f}\\pm{k['sigma_A']:.2f}$")
    ax[0].plot(L[m],(L*th)[m],color=C[s],lw=1,ls=ls,alpha=.8)
    # middle: deflection E-mode x kappa_CMB over the science window
    b=c['bands']['science']; L=np.array(b['L']); cx=np.array(b['E_cross']); e=np.array(b['E_err']); th=np.array(b['prediction']); used=np.array(b['used'],bool); m=L>=40
    ax[1].errorbar(L[m]+OFF[s],(L**2*cx)[m],(L**2*e)[m],fmt='o',ms=4,color=C[s],capsize=2,label=f"{s}: $A_L={b['A']:.2f}\\pm{b['sigma_A']:.2f}$")
    ax[1].plot(L[m],(L**2*th)[m],color=C[s],lw=1,ls=ls,alpha=.8)
    # right: ratio for the deflection (filled: bandpowers in the fit; open: outside it)
    r=cx/th; re=e/np.abs(th)
    for sel,mfc in ((m&used,C[s]),(m&~used,'white')):
        if sel.any(): ax[2].errorbar(L[sel]+OFF[s],r[sel],re[sel],fmt='o',ms=4,color=C[s],mfc=mfc,capsize=2,label=(f"{s}" if mfc!='white' else None))
ax[0].set(xlim=(20,1000),xlabel=r'$\ell$',ylabel=r'$\ell\,C_\ell^{T\kappa_{\rm CMB}}$',title='convergence template $\\times$ CMB convergence'); ax[0].legend(fontsize=8)
ax[1].set(xlim=(20,1000),xlabel=r'$\ell$',ylabel=r'$\ell^2\,C_\ell^{E\kappa_{\rm CMB}}$',title='deflection template (E, $40\\leq \\ell\\leq1000$) $\\times$ CMB convergence'); ax[1].legend(fontsize=8)
ax[2].axhline(1,color='k',lw=1); ax[2].set(xlim=(20,1000),ylim=(0,2),xlabel=r'$\ell$',ylabel='measured / predicted',title='deflection: ratio to the template prediction'); ax[2].legend(fontsize=8)
for x in ax: x.axvspan(0,40,color='0.9')
fig.suptitle('Combined template $\\times$ CMB lensing (NaMaster bandpowers, Gaussian errors; lines: the prediction for the template through each mask pair)',fontsize=9)
fig.tight_layout()
for out in (ROOT/'report/lowz/figures/template_cmb_cross.pdf',ROOT/'Paper/figures/template_cmb_cross.pdf'):
    if out.parent.exists(): fig.savefig(out); print('wrote',out)
