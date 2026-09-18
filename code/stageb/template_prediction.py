"""Prediction of <T kappa'> for the Wiener-combined tracer template on a given sky weight.

The slice template is built per coverage class c (subset of tracers covering a pixel): T_s = sum_c 1_c sum_{k in c}
W^c_k * m_k, with the Wiener weights W^c solving [S pw^2 + N]_{cc} W^c = S pw on the tracers of the class. A
pseudo-spectrum of T against a field kappa' on the weight M (the product of the two NaMaster masks) is
sum_c <1_c M> sum_k W^c_k C^{(l_s, c)}, so the effective weight that predicts it is the class-fraction average
INSIDE M, not over the whole union footprint (which `summary.json`'s `effective_weight` uses). The difference is
large where the classes are unevenly distributed on the sky: the northern DESI imaging region is mostly BOSS-only
coverage with low weights, the ACT overlap is DESI-rich.
"""
from __future__ import annotations
import json, sys
from pathlib import Path
from itertools import combinations
import numpy as np
import healpy as hp
HERE=Path(__file__).resolve().parent; CODE=HERE.parent
for p in (CODE,CODE/'pipeline',HERE):
    if str(p) not in sys.path: sys.path.insert(0,str(p))
from lowz import slice_spectra


class TemplatePrediction:
    def __init__(self,lowz_dir,cref,lmax,nside):
        self.dir=Path(lowz_dir); self.summary=json.loads((self.dir/'summary.json').read_text()); self.lmax=int(lmax); self.nside=int(nside)
        self.cref=float(self.summary.get('chi_ref',cref)); self.pw=hp.pixwin(nside,lmax=lmax); ell=np.arange(lmax+1); self.slices=[]
        for s in self.summary['slices']:
            labels=list(s['tracers']); masks=[hp.read_map(str(self.dir/f'unitbias_{lab}_nside{nside}.fits'),field=1)>0.5 for lab in labels]
            Lth,C=slice_spectra(s['zmin'],s['zmax'],self.cref,lmax); S=np.interp(ell,Lth,C[:,1,1]); Clc=np.interp(ell,Lth,C[:,1,2])
            # iteration-7 templates (lowz_split) had a diagonal noise matrix and ONE class, the intersection of the masks
            legacy='noise_matrix' not in s
            N=np.diag(np.asarray(s['shot_s'],float)) if legacy else np.asarray(s['noise_matrix'],float); classes=[]
            for r in range(1,len(masks)+1):
                for sub in combinations(range(len(masks)),r):
                    if legacy and len(sub)<len(masks): continue
                    m=np.ones(len(masks[0]),bool)
                    for i in range(len(masks)): m&=masks[i] if (i in sub or legacy) else ~masks[i]
                    if not m.any(): continue
                    idx=list(sub); Wt=np.zeros(lmax+1)
                    for l in range(2,lmax+1):
                        Cm=S[l]*self.pw[l]**2*np.ones((len(idx),len(idx)))+N[np.ix_(idx,idx)]
                        Wt[l]=np.linalg.solve(Cm,np.full(len(idx),S[l]*self.pw[l])).sum()
                    classes.append((m,Wt))
            union=np.zeros(len(masks[0]),bool)
            for m in masks: union|=m
            if legacy: union=np.prod(masks,axis=0).astype(bool)
            self.slices.append({'name':f"slice_{s['zmin']:g}_{s['zmax']:g}",'Clc':Clc,'classes':classes,'union':union})

    def effective_weight(self,M,name):
        """w_eff(ell) of one slice on the sky weight M (float map, the PRODUCT of the two NaMaster masks). The class
        fractions are normalised by the area of M itself: NaMaster attributes the cross-spectrum to the whole field
        mask, so where a slice's tracers cover only part of it the slice contributes proportionally less (the
        high-z slices are DESI-only while the combined footprint includes BOSS-only sky)."""
        s=next(x for x in self.slices if x['name']==name); M=np.asarray(M,float); tot=float(np.sum(M))
        w=np.zeros(self.lmax+1)
        for m,Wt in s['classes']: w+=float(np.sum(M*m))/max(tot,1e-30)*Wt
        return w

    def cross(self,M,name='combined'):
        """<T kappa_CMB>_ell for the template ``name`` on the sky weight M (pixel window included)."""
        names=[s['name'] for s in self.slices] if name=='combined' else [name]
        return sum(self.effective_weight(M,n)*next(s for s in self.slices if s['name']==n)['Clc'] for n in names)*self.pw

    def union_weight(self,name):
        s=next(x for x in self.slices if x['name']==name); return self.effective_weight(s['union'].astype(float),name)
