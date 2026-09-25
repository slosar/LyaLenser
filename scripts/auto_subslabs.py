"""Redshift split of the forest auto-correlation lensing amplitude (iteration 12 companion): re-accumulate the
iteration-11 pair catalogue in sub-slabs of the pair mean redshift with the fitted layered table, and fit the
combined template per sub-slab. Writes --out/auto_subslabs.json (and the fits into --run/fits.h5).
Usage: python auto_subslabs.py --run $LYALENSER_DATA/stageb/dr1_lowz_v7 --lowz $LYALENSER_DATA/lowz_v3 --edges 1.96 2.25 2.55 3.0
"""
from __future__ import annotations
import argparse, json, sys, time
from pathlib import Path
import numpy as np
import healpy as hp
import sys; from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # repository root: `lyalenser` imports without installation
from lyalenser.paths import DATA
from lyalenser.config import production_config
from lyalenser.tables import read_xi
from lyalenser.cosmo import chi as chi_of_z
from lyalenser.pairs import find_pairs, accumulate, pair_midpoint_regions
from lyalenser.templates import sphere_band_templates
from lyalenser.amplitude import amplitude, curl_amplitude
from lyalenser.amplitude import common_science
from lyalenser.desi_io import load_sightlines


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--run',type=Path,default=DATA/'stageb/dr1_lowz_v7'); ap.add_argument('--lowz',type=Path,default=DATA/'lowz_v3')
    ap.add_argument('--edges',type=float,nargs='+',default=[1.96,2.25,2.55,3.0]); ap.add_argument('--nside-jk',type=int,default=8); ap.add_argument('--nside-alpha',type=int,default=1024)
    ap.add_argument('--bands',type=float,nargs='+',default=None); ap.add_argument('--fits',type=Path,default=None,help='fits file (default --run/fits.h5)'); ap.add_argument('--out',type=Path,default=None)
    a=ap.parse_args(); t0=time.perf_counter(); log=json.loads((a.run/'dr1_lowz.json').read_text())
    cfg=production_config().copy(xi_correction='spline',xi_z_evolution=True,slabs=tuple((a.edges[i],a.edges[i+1]) for i in range(len(a.edges)-1)),chi_ref=float(log['chi_ref']))
    sl=load_sightlines(a.run/'sightlines.h5'); tab=read_xi(a.run/'xi.h5','xi')
    pairs=find_pairs(sl,cfg.r_perp_max/float(sl.chi.min())); print(f'{sl.nq} forests, {len(pairs[0])} pairs ({time.perf_counter()-t0:.0f} s)',flush=True)
    from lyalenser.templates import SCIENCE_BANDS
    BANDS=tuple((int(a.bands[i]),int(a.bands[i+1])) for i in range(len(a.bands)-1)) if a.bands else SCIENCE_BANDS
    fits_path=a.fits or (a.run/'fits.h5'); out_path=a.out or (a.run/'auto_subslabs.json')
    alm=hp.read_alm(str(a.lowz/'kappa_combined_alm.fits')); templates,_=sphere_band_templates(alm,sl.ra,sl.dec,nside=a.nside_alpha,science_bands=BANDS,source='combined')
    out={'edges':a.edges,'science_bands':[list(b) for b in BANDS],'sub_slabs':{}}
    for k in range(len(a.edges)-1):
        c=cfg.copy(slab_index=k); cat=accumulate(sl,pairs,tab,c); reg=pair_midpoint_regions(cat,sl,a.nside_jk)
        r=amplitude(cat,templates,c.g1,reg); s=common_science(r); cu,cue=curl_amplitude(r); r.save(fits_path,f'subslab_{a.edges[k]:g}_{a.edges[k+1]:g}')
        out['sub_slabs'][f'{a.edges[k]:g}-{a.edges[k+1]:g}']={'A':s['A'],'jk_error':s['jk_error'],'sigma_F':s['sigma_F'],'curl':cu,'curl_jk_error':cue,'bands':r.A.tolist(),'band_errors':r.jk_error.tolist(),'pairs':int(len(cat.a))}
        print(f"sub-slab {a.edges[k]:g}-{a.edges[k+1]:g}: A = {s['A']:.3f} +- {s['jk_error']:.3f} (sigma_F {s['sigma_F']:.3f}) ({time.perf_counter()-t0:.0f} s)",flush=True)
        out_path.write_text(json.dumps(out,indent=1)+'\n')
    print('done')


if __name__=='__main__': main()
