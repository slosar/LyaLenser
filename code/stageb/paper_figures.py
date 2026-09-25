"""Paper figures from the committed JSON products (iteration 14).
  redshift_split.pdf : auto, cross and combined amplitudes per range of the pair mean redshift (fiducial).
  robustness.pdf     : one horizontal 1-sigma bar per analysis variant, lines at A = 0 and 1.
Usage: python paper_figures.py [--split] [--robustness]
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2]; R=ROOT/'report/stageb'; OUT=[ROOT/'Paper/figures',ROOT/'report/lowz/figures']
plt.rcParams.update({'font.size':9,'axes.grid':True,'grid.alpha':.25,'legend.frameon':False})
C={'auto':'#0072B2','cross':'#D55E00','comb':'k'}


def save(fig,name):
    for o in OUT: o.mkdir(parents=True,exist_ok=True); fig.savefig(o/name,bbox_inches='tight')
    plt.close(fig); print('wrote',name)


def redshift_split(comb_json='auto_cross_combination_v4.json'):
    c=json.load(open(R/comb_json)); rows=[(k.replace('subslab_','').replace('_','--'),v) for k,v in c['sub_slabs'].items()]
    rows.sort(key=lambda kv: float(kv[0].split('--')[0])); x=np.arange(len(rows)); fig,ax=plt.subplots(figsize=(5,3.2))
    for off,key,lab,col in ((-.22,'auto','forest $\\times$ forest',C['auto']),(0.,'cross','quasar $\\times$ forest',C['cross']),(.22,'comb','combined',C['comb'])):
        A=[v[key]['A'] if key!='comb' else v['A'] for _,v in rows]; e=[v['errors'][0] if key=='auto' else v['errors'][1] if key=='cross' else v['error'] for _,v in rows]
        ax.errorbar(x+off,A,e,fmt='o' if key!='comb' else 's',color=col,ms=5,capsize=3,label=lab)
    ax.axhline(1,color='0.4',ls='--',lw=1); ax.axhline(0,color='0.4',lw=.8); ax.set_xticks(x); ax.set_xticklabels([k for k,_ in rows]); ax.set(xlabel='pair mean redshift',ylabel='$A_L$'); ax.legend(fontsize=8,loc='upper right')
    g=c['combined']; ax.set_title(f"combined over all redshifts: $A_L={g['A']:.2f}\\pm{g['error']:.2f}$",fontsize=9); save(fig,'redshift_split.pdf')


def slice_split(joint_json='joint_fit_v4.json'):
    """Per tracer slice: the forest x forest, quasar x forest and combined amplitudes of the joint fit (jackknife errors), plus all slices combined."""
    j=json.load(open(R/joint_json)); sl=j['slices']; labels=[s.replace('slice_','').replace('_','--') for s in sl]+['all slices']
    x=np.arange(len(labels)); fig,ax=plt.subplots(figsize=(5.4,3.2))
    for off,key,lab,col in ((-.22,'auto','forest $\\times$ forest',C['auto']),(0.,'cross','quasar $\\times$ forest',C['cross']),(.22,'comb','combined',C['comb'])):
        if key=='comb':
            A=[j['combination']['per_slice'][s]['A'] for s in sl]+[j['combination']['global']['A']]; e=[j['combination']['per_slice'][s]['error'] for s in sl]+[j['combination']['global']['error']]
        else:
            st=j['statistics'][key]; A=list(st['per_slice']['A'])+[st['global']['A'][0]]; e=list(st['per_slice']['jk_error'])+[st['global']['jk_error'][0]]
        ax.errorbar(x+off,A,e,fmt='o' if key!='comb' else 's',color=col,ms=5,capsize=3,label=lab)
    ax.axvline(len(sl)-.5,color='0.6',lw=.8); ax.axhline(1,color='0.4',ls='--',lw=1); ax.axhline(0,color='0.4',lw=.8); ax.set_xticks(x); ax.set_xticklabels(labels,fontsize=8); ax.set(xlabel='tracer slice (redshift)',ylabel='$A_L$'); ax.legend(fontsize=8,loc='lower left')
    save(fig,'slice_split.pdf')


def band_split(joint_json='joint_fit_v4.json'):
    """Per science band: the forest x forest, quasar x forest and combined amplitudes of the joint fit (jackknife errors), plus all bands combined."""
    j=json.load(open(R/joint_json)); bands=j['statistics']['auto']['per_band']['groups']; labels=[b.replace('L','').replace('_','--') for b in bands]+['all bands']
    x=np.arange(len(labels)); fig,ax=plt.subplots(figsize=(5.4,3.2))
    for off,key,lab,col in ((-.22,'auto','forest $\\times$ forest',C['auto']),(0.,'cross','quasar $\\times$ forest',C['cross']),(.22,'comb','combined',C['comb'])):
        if key=='comb':
            A=[j['combination']['per_band'][b]['A'] for b in bands]+[j['combination']['global']['A']]; e=[j['combination']['per_band'][b]['error'] for b in bands]+[j['combination']['global']['error']]
        else:
            st=j['statistics'][key]; A=list(st['per_band']['A'])+[st['global']['A'][0]]; e=list(st['per_band']['jk_error'])+[st['global']['jk_error'][0]]
        ax.errorbar(x+off,A,e,fmt='o' if key!='comb' else 's',color=col,ms=5,capsize=3,label=lab)
    ax.axvline(len(bands)-.5,color='0.6',lw=.8); ax.axhline(1,color='0.4',ls='--',lw=1); ax.axhline(0,color='0.4',lw=.8); ax.set_xticks(x); ax.set_xticklabels(labels,fontsize=8); ax.set(xlabel=r'template multipole band $\ell$',ylabel='$A_L$',ylim=(-6,9)); ax.legend(fontsize=8,loc='upper left')
    save(fig,'band_split.pdf')


def robustness(rows):
    """rows: list of (label, A, err, kind) top to bottom."""
    fig,ax=plt.subplots(figsize=(5.6,0.48*len(rows)+1.2)); y=np.arange(len(rows))[::-1]
    for yi,(lab,A,e,kind) in zip(y,rows):
        ax.errorbar(A,yi,xerr=e,fmt='o',color=C.get(kind,'k'),ms=5,capsize=3,lw=1.5 if kind=='comb' else 1)
    ax.axvline(0,color='0.4',lw=.8); ax.axvline(1,color='0.4',ls='--',lw=1); ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows],fontsize=13); ax.tick_params(axis='x',labelsize=13); ax.set_xlabel('$A_L$',fontsize=18); ax.set_xlim(-1.2,2.2); ax.grid(axis='y',alpha=0)
    save(fig,'robustness.pdf')


def collect_robustness():
    j=json.load(open(R/'joint_fit_v4.json')); c=j['combination']; g=c['global']
    rows=[('fiducial (joint $R$, $40\\leq L\\leq 1000$, $r_\\perp\\leq 30$)',g['A'],g['error'],'comb'),
          ('forest $\\times$ forest only',g['auto']['A'] if 'auto' in g else j['statistics']['auto']['global']['A'][0],g['errors'][0],'auto'),
          ('quasar $\\times$ forest only',j['statistics']['cross']['global']['A'][0],g['errors'][1],'cross')]
    for lab,f in (('$L_{\\max}=1300$','auto_cross_combination_l1300v4.json'),('$L_{\\max}=500$','auto_cross_combination_l500v4.json'),('$r_\\perp\\leq 40\\,h^{-1}$Mpc','auto_cross_combination_rp40.json'),('$r_\\perp\\leq 20\\,h^{-1}$Mpc','auto_cross_combination_rp20.json')):
        if (R/f).exists(): v=json.load(open(R/f))['combined']; rows.append((lab,v['A'],v['error'],'comb'))
        else: print('missing',f)
    for lab,f in (('spline correction (30 / 42 coefficients)','auto_cross_combination_kmed.json'),('spline correction (42 / 63 coefficients)','auto_cross_combination_kfine.json')):
        if (R/f).exists(): v=json.load(open(R/f))['combined']; rows.append((lab,v['A'],v['error'],'comb'))
        else: print('missing',f)
    b=c['block_diagonal_global']; rows.append(('block-diagonal $R$ (per-slice fits)',b['A'],b['error'],'comb'))
    return rows


if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--split',action='store_true'); ap.add_argument('--slices',action='store_true'); ap.add_argument('--bands',action='store_true'); ap.add_argument('--robustness',action='store_true'); a=ap.parse_args()
    if a.split: redshift_split()
    if a.slices: slice_split()
    if a.bands: band_split()
    if a.robustness: robustness(collect_robustness())
