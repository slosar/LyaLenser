"""Optimal combination of the forest x forest and quasar x forest lensing amplitudes (iteration 12), overall, per
template band and per redshift sub-slab, with the joint jackknife covariance (both statistics use the nside-8
pair-midpoint regions; the union of the region ids is used, absent regions contribute zero). Reports the relative
signal-to-noise (inverse-variance weights). Writes results/auto_cross_combination.json.
Usage: python combine_auto_cross.py --auto $LYALENSER_DATA/stageb/dr1_lowz_v7 --cross $LYALENSER_DATA/stageb/dr1_qso_v1
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
import numpy as np
import h5py
import sys; from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # repository root: `lyalenser` imports without installation
from lyalenser.paths import DATA
from lyalenser.amplitude import common_science
from lyalenser.amplitude import AmplitudeResult


def load_fit(path,group):
    with h5py.File(path) as f:
        g=f[group]; d={k:g[k][()] for k in ('q','F','mf','A','sigma_F','jk_samples','jk_cov','regions','partial_q','partial_F','partial_mf')}
        names=json.loads(g.attrs['names']); attrs=json.loads(g.attrs['config'])
    return AmplitudeResult(names,d['q'],d['F'],d['mf'],d['A'],d['sigma_F'],d['jk_samples'],d['jk_cov'],d['regions'],d['partial_q'],d['partial_F'],d['partial_mf'],attrs)


def science_jk(r):
    """Common science amplitude, its jackknife samples (per region), the Fisher scale and the region ids."""
    s=common_science(r); return float(s['A']),np.asarray(s['jk'],float),float(s['sigma_F']),np.asarray(r.regions)


def band_jk(r,k):
    return float(r.A[k]),np.asarray(r.jk_samples[:,k],float),float(r.sigma_F[k]),np.asarray(r.regions)


def combine(A1,jk1,reg1,A2,jk2,reg2):
    """Joint jackknife covariance on the union of regions: a region missing from one statistic leaves that
    statistic at its full value in that pseudo-sample."""
    regs=np.union1d(reg1,reg2); n=len(regs)
    s1=np.full(n,A1); s2=np.full(n,A2); s1[np.searchsorted(regs,reg1)]=jk1; s2[np.searchsorted(regs,reg2)]=jk2
    S=np.column_stack((s1,s2)); d=S-S.mean(axis=0); C=(n-1)/n*d.T@d
    Ci=np.linalg.pinv(C)*(n-2-2)/(n-1); one=np.ones(2); var=1/float(one@Ci@one); Ac=float(one@Ci@np.array([A1,A2])*var)
    w=(Ci@one)*var; rho=C[0,1]/np.sqrt(C[0,0]*C[1,1])
    return {'A':Ac,'error':float(np.sqrt(var)),'weights':w.tolist(),'correlation':float(rho),'errors':[float(np.sqrt(C[0,0])),float(np.sqrt(C[1,1]))],
            'snr_fraction':[float(w[0]),float(w[1])],'regions':int(n)}


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--auto',type=Path,default=DATA/'stageb/dr1_lowz_v7'); ap.add_argument('--cross',type=Path,default=DATA/'stageb/dr1_qso_v1')
    ap.add_argument('--out',type=Path,default=ROOT/'results/auto_cross_combination.json'); a=ap.parse_args()
    out={'auto':str(a.auto),'cross':str(a.cross)}
    ra=load_fit(a.auto/'fits.h5','combined'); rc=load_fit(a.cross/'fits.h5','combined')
    A1,j1,f1,g1=science_jk(ra); A2,j2,f2,g2=science_jk(rc); out['combined']=dict(combine(A1,j1,g1,A2,j2,g2),auto={'A':A1,'sigma_F':f1},cross={'A':A2,'sigma_F':f2})
    c=out['combined']; print(f"combined template: auto {A1:.3f} +- {c['errors'][0]:.3f} (F {f1:.3f}), cross {A2:.3f} +- {c['errors'][1]:.3f} (F {f2:.3f}), corr {c['correlation']:.2f} -> {c['A']:.3f} +- {c['error']:.3f}; weights {np.round(c['weights'],2)}")
    nb=int(ra.attrs.get('n_science',5)); out['bands']={}
    for k in range(nb):
        A1,j1,f1,g1=band_jk(ra,k); A2,j2,f2,g2=band_jk(rc,k); out['bands'][ra.names[k]]=dict(combine(A1,j1,g1,A2,j2,g2),auto={'A':A1,'sigma_F':f1},cross={'A':A2,'sigma_F':f2})
        b=out['bands'][ra.names[k]]; print(f"  band {ra.names[k]:9s}: auto {A1:6.2f} +- {b['errors'][0]:.2f}, cross {A2:6.2f} +- {b['errors'][1]:.2f} -> {b['A']:6.2f} +- {b['error']:.2f}; weights {np.round(b['weights'],2)}")
    # slices
    out['slices']={}
    with h5py.File(a.auto/'fits.h5') as f: groups=[k for k in f if k.startswith('slice_')]
    for gname in groups:
        try: r1=load_fit(a.auto/'fits.h5',gname); r2=load_fit(a.cross/'fits.h5',gname)
        except KeyError: continue
        A1,j1,f1,g1=science_jk(r1); A2,j2,f2,g2=science_jk(r2); out['slices'][gname]=dict(combine(A1,j1,g1,A2,j2,g2),auto={'A':A1,'sigma_F':f1},cross={'A':A2,'sigma_F':f2})
        b=out['slices'][gname]; print(f"  {gname:14s}: auto {A1:6.2f} +- {b['errors'][0]:.2f}, cross {A2:6.2f} +- {b['errors'][1]:.2f} -> {b['A']:6.2f} +- {b['error']:.2f}")
    # sub-slabs
    out['sub_slabs']={}
    with h5py.File(a.auto/'fits.h5') as f: groups=[k for k in f if k.startswith('subslab_')]
    for gname in groups:
        try: r1=load_fit(a.auto/'fits.h5',gname); r2=load_fit(a.cross/'fits.h5',gname)
        except KeyError: continue
        A1,j1,f1,g1=science_jk(r1); A2,j2,f2,g2=science_jk(r2); out['sub_slabs'][gname]=dict(combine(A1,j1,g1,A2,j2,g2),auto={'A':A1,'sigma_F':f1},cross={'A':A2,'sigma_F':f2})
        b=out['sub_slabs'][gname]; print(f"  {gname:18s}: auto {A1:6.2f} +- {b['errors'][0]:.2f}, cross {A2:6.2f} +- {b['errors'][1]:.2f} -> {b['A']:6.2f} +- {b['error']:.2f}; weights {np.round(b['weights'],2)}")
    a.out.write_text(json.dumps(out,indent=1)+'\n'); print('saved',a.out)


if __name__=='__main__': main()
