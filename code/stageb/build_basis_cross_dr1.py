"""Basis tables for the quasar x forest cross-correlation: the three Hankel basis spectra of the forest basis file
(`raw/mu0, mu2, mu4`, DESI pixel windows included) pushed through the ONE-SIDED continuum projection averaged over
real DR1 forests (`xi_cross.project_cross_sample`), on the signed r_par grid [-xi_max, xi_max].
Usage: python build_basis_cross_dr1.py --basis $LYALENSER_DATA/stageb/basis_dr1_ab_z196.h5 --out .../basis_cross_dr1_z196.h5 [--forests 400]
"""
from __future__ import annotations
import argparse, sys, time
from pathlib import Path
import numpy as np
import h5py
HERE=Path(__file__).resolve().parent; CODE=HERE.parent
for p in (CODE,CODE/'pipeline',HERE):
    if str(p) not in sys.path: sys.path.insert(0,str(p))
from paths import DATA
from campaign4 import campaign_config, read_xi
from xi_fit import BASIS
from xi_cross import project_cross_sample, refine_cross_grid
from build_basis_dr1 import forest_geometries


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--basis',type=Path,default=DATA/'stageb/basis_dr1_ab_z196.h5'); ap.add_argument('--out',type=Path,default=DATA/'stageb/basis_cross_dr1_z196.h5')
    ap.add_argument('--forests',type=int,default=400); ap.add_argument('--zmin',type=float,default=1.96); ap.add_argument('--zmax',type=float,default=3.0)
    ap.add_argument('--rp-step',type=float,default=0.5); ap.add_argument('--q-step',type=float,default=4.0); ap.add_argument('--modulus',type=int,default=20)
    a=ap.parse_args(); cfg=campaign_config(1.); t0=time.perf_counter()
    raw={k:read_xi(a.basis,f'raw/{k}') for k in BASIS}
    forests=forest_geometries(a.zmin,a.zmax,seed=7,forest_regions=('lya','lyb'),modulus=a.modulus,verbose=False)
    print(f'{len(forests)} sightlines for the projection ({time.perf_counter()-t0:.0f} s)',flush=True)
    rp_grid=np.arange(0,cfg.xi_max+.5*a.rp_step,a.rp_step)
    proj=project_cross_sample(raw,forests,cfg,n_forests=a.forests,seed=11,rp_grid=rp_grid,q_step=a.q_step,verbose=True)
    proj=refine_cross_grid(proj,cfg)
    if a.out.exists(): a.out.unlink()
    for k in BASIS: raw[k].save(a.out,f'raw/{k}'); proj[k].save(a.out,f'projected/{k}')
    with h5py.File(a.out,'a') as f:
        f.attrs.update(provider='hankel',projection='one-sided over real DR1 forests',n_forests=a.forests,q_step=a.q_step,source_basis=str(a.basis),signed_rpar=True)
    print(f'wrote {a.out} ({time.perf_counter()-t0:.0f} s)')


if __name__=='__main__': main()
