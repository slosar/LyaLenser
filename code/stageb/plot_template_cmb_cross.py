"""Figure: the combined convergence template x ACT and x Planck bandpowers on one plot, with the prediction
(sum over slices of w_eff C^{l c} pw) and the same-sky comparison. Reads report/stageb/cmb_map_checks.json;
writes report/lowz/figures/template_cmb_cross.pdf."""
import json, sys
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2]
d=json.load(open(ROOT/'report/stageb/cmb_map_checks.json')); L=np.array(d['L'])
C={'ACT':'#D55E00','Planck':'#0072B2'}
plt.rcParams.update({'font.size':9,'axes.grid':True,'grid.alpha':.25,'legend.frameon':False})
fig,ax=plt.subplots(1,3,figsize=(12,3.6))
for panel,key,title in ((0,'full_overlap','full overlap with each map'),(1,'same_sky',f"same sky (fsky {d['same_sky']['fsky']:.3f})")):
    for s in ('ACT','Planck'):
        r=d[key][s]; cx=np.array(r['cross']); e=np.array(r['err']); th=np.array(r['theory']); off=-4 if s=='ACT' else 4; m=L>=40
        ax[panel].errorbar(L[m]+off,(L*cx)[m],(L*e)[m],fmt='o',ms=4,color=C[s],capsize=2,label=f"{s}: $A_L={r['A']:.2f}\\pm{r['sigma_A']:.2f}$")
        ax[panel].plot(L[m],(L*th)[m],color=C[s],lw=1,ls='--' if s=='Planck' else '-',alpha=.8)
    ax[panel].axvspan(0,40,color='0.9'); ax[panel].set(xlim=(20,1000),xlabel=r'$L$',ylabel=r'$L\,C_L^{T\kappa_{\rm CMB}}$',title=title); ax[panel].legend(fontsize=8)
# ratio panel
for s in ('ACT','Planck'):
    for key,mk,ls in (('full_overlap','o','-'),('same_sky','s','--')):
        r=d[key][s]; cx=np.array(r['cross']); e=np.array(r['err']); th=np.array(r['theory']); m=(L>=40)&(L<=600)
        ax[2].errorbar(L[m]+(-6 if s=='ACT' else 6)+(0 if key=='full_overlap' else 3),(cx/th)[m],(e/np.abs(th))[m],fmt=mk,ms=4,color=C[s],mfc=(C[s] if key=='full_overlap' else 'white'),capsize=2,label=f"{s}, {key.replace('_',' ')}")
ax[2].axhline(1,color='k',lw=1); ax[2].set(xlim=(20,620),ylim=(0,1.6),xlabel=r'$L$',ylabel='measured / predicted',title='ratio to the template prediction'); ax[2].legend(fontsize=7,ncol=2)
fig.suptitle('Combined convergence template $\\times$ CMB convergence (NaMaster bandpowers, Gaussian errors); lines: prediction $\\sum_s w_{{\\rm eff},s}C_L^{(\\ell_s,c)}p_L$ through each mask',fontsize=9)
fig.tight_layout(); out=ROOT/'report/lowz/figures/template_cmb_cross.pdf'; fig.savefig(out); print('wrote',out)
