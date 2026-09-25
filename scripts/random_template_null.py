"""Repeat the combined-template error/null check on cached production pairs.

The seed sequence and map generation match the original DR1 drivers, so the
first 20 draws reproduce their check. Draws are independent, not +/- pairs.
The 11-component fit is a diagnostic, not the full 55-component slice fit.
Use --start/--stop for independent shards; each draw is checkpointed to JSON.
"""
from __future__ import annotations
import argparse
import json
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import healpy as hp
import h5py
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from lyalenser.paths import DATA
from lyalenser.pairs import PairCatalogue, pair_midpoint_regions
from lyalenser.templates import sphere_band_templates
from lyalenser.joint_fit import build_joint, JointFit
from lyalenser.qso_io import read_quasars
from lyalenser.xi_cross import Positions


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--statistic',choices=['auto','cross'],required=True)
    ap.add_argument('--auto',type=Path,default=DATA/'stageb/dr1_lowz_v7d')
    ap.add_argument('--cross',type=Path,default=DATA/'stageb/dr1_qso_v1d')
    ap.add_argument('--lowz',type=Path,default=DATA/'lowz_v4')
    ap.add_argument('--seed',type=int,default=2026)
    ap.add_argument('--start',type=int,default=0)
    ap.add_argument('--stop',type=int,default=100)
    ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args(); t0=time.perf_counter()
    run=a.auto if a.statistic=='auto' else a.cross
    log=json.loads((run/('dr1_lowz.json' if a.statistic=='auto' else 'dr1_qso.json')).read_text())
    bands=tuple(map(tuple,log['science_bands']))
    with h5py.File(a.auto/'sightlines.h5') as f:
        pos=SimpleNamespace(ra=f['sightlines/ra'][()],dec=f['sightlines/dec'][()])
    if a.statistic=='cross': pos=Positions(pos,read_quasars(*log['quasar_z']))
    with h5py.File(run/'catalogue.h5') as f:
        g=f['all']; cat=PairCatalogue(*(g[k][()] for k in ('a','b','thx','thy','theta','accum','npair')),dict(g.attrs))
    cat.accum=cat.accum.astype(np.float64)
    regions=pair_midpoint_regions(cat,pos,8)
    alm=hp.read_alm(str(a.lowz/'kappa_combined_alm.fits')); lmax=hp.Alm.getlmax(len(alm))
    mask=hp.read_map(str(a.lowz/'mask_combined_nside512.fits'))
    cl=hp.alm2cl(alm)/float(np.mean(mask**2))
    seeds=np.random.default_rng(a.seed).integers(2**31,size=a.stop)
    out={'statistic':a.statistic,'seed':a.seed,'start':a.start,'stop':a.stop,
         'fit':'combined template: 5 science bands, 5 curls, junk; common science amplitude',
         'bands':bands,'g1':log['config']['g1'],'nside_alpha':2048,'draws':[]}
    if a.out.exists():
        previous=json.loads(a.out.read_text())
        for key in ('statistic','seed','start','stop'):
            if previous[key]!=out[key]: raise ValueError(f'incompatible checkpoint: {key}')
        out['draws']=previous['draws']
    done={d['index'] for d in out['draws']}
    for i in range(a.start,a.stop):
        if i in done: continue
        np.random.seed(int(seeds[i]))
        sky=hp.synfast(cl,512,lmax=lmax)*mask
        ts,_=sphere_band_templates(hp.map2alm(sky,lmax=lmax,iter=0),pos.ra,pos.dec,
                                   nside=2048,science_bands=bands,source=f'random {i}')
        jf=build_joint(cat,{'combined':ts},log['config']['g1'],regions)
        collapsed=JointFit(jf.names,jf.kinds,['all' if g is not None else None for g in jf.groups],
                           jf.regvals,jf.pq,jf.pF,jf.pmf)
        fit=collapsed.fit(['all'])
        d={'index':i,'map_seed':int(seeds[i]),'A':fit['A'][0],'jk_error':fit['jk_error'][0],
           'sigma_F':fit['sigma_F'][0]}
        # Flipping every component must reverse A and leave its uncertainty fixed.
        # This is a sign-symmetry diagnostic, not another independent realisation.
        if i==a.start:
            neg=JointFit(collapsed.names,collapsed.kinds,collapsed.groups,collapsed.regvals,
                         -collapsed.pq,collapsed.pF,-collapsed.pmf).fit(['all'])
            d['sign_reversal_A_sum']=fit['A'][0]+neg['A'][0]
            d['sign_reversal_error_difference']=fit['jk_error'][0]-neg['jk_error'][0]
        out['draws'].append(d); out['wall_s']=time.perf_counter()-t0
        a.out.parent.mkdir(parents=True,exist_ok=True)
        a.out.write_text(json.dumps(out,indent=1)+'\n')
        print(f'{a.statistic} {i}: {d["A"]:.5f} +/- {d["jk_error"]:.5f}; {out["wall_s"]:.0f} s',flush=True)


if __name__=='__main__': main()
