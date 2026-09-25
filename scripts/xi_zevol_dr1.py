"""Redshift evolution of the DR1 forest correlation (Stage B iteration 10): the z-resolved cells of a saved
sightline set, the evolving fit (`xi_zevol`), and the A x A / A x B comparison at fixed redshift.

The question this answers: is the 0.84 amplitude ratio of the A x B to the A x A correlation (iteration 9) the
redshift evolution of the forest correlation (the A x B pairs sit at lower redshift) or a property of the
region-B deltas? Fits reported: the joint evolving model (spline correction) on all pairs; the base-only
evolving model; the same on A x A and A x B separately with the exponents fixed at the joint values (the
amplitude ratio at fixed z); and the non-evolving per-type fits of iteration 9 for reference.

Usage: python xi_zevol_dr1.py --run $LYALENSER_DATA/stageb/dr1_lowz_v4 --basis $LYALENSER_DATA/stageb/basis_dr1_ab.h5
       [--out results/xi_zevol_dr1.json] [--z-edges 2.1 2.2 2.3 2.4 2.55 2.75 3.0]
"""
from __future__ import annotations
import argparse, json, sys, time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
import numpy as np
import sys; from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # repository root: `lyalenser` imports without installation
from lyalenser.paths import DATA
from lyalenser.config import production_config
from lyalenser.xi_model import xi_from_data
from lyalenser.xi_zevol import evolving_table_for, PARAMS
from lyalenser.desi_io import load_sightlines
from dry_run_lowz import load_basis
from lyalenser.tables import table_for


def summarize(ft,keys=('b_F2','beta_F','gamma_b','gamma_beta','gamma_S','amplitude_exponent_total','chi2','cells','dof','n_correction','free','per_z')):
    return {k:ft.params[k] for k in keys if k in ft.params}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--run',type=Path,default=DATA/'stageb/dr1_lowz_v4'); ap.add_argument('--basis',type=Path,default=DATA/'stageb/basis_dr1_ab.h5')
    ap.add_argument('--out',type=Path,default=ROOT/'results/xi_zevol_dr1.json')
    ap.add_argument('--z-edges',type=float,nargs='+',default=None); ap.add_argument('--z-ref',type=float,default=2.4)
    ap.add_argument('--xi-ridge',type=float,default=1e-2); ap.add_argument('--save-table',type=Path,default=None)
    a=ap.parse_args(); t0=time.perf_counter()
    cfg=production_config().copy(xi_correction='spline',xi_correction_ridge=a.xi_ridge,xi_z_evolution=True,xi_z_ref=a.z_ref)
    if a.z_edges: cfg=cfg.copy(xi_z_edges=tuple(a.z_edges))
    sl=load_sightlines(a.run/'sightlines.h5'); basis=load_basis(a.basis)
    print(f'{sl.nq} sightlines, {len(sl.chi)} pixels ({time.perf_counter()-t0:.0f} s)',flush=True)
    measured=xi_from_data(sl,cfg)
    print(f"measured; pair weight by z: {measured.meta['pair_weight_by_z']} ({time.perf_counter()-t0:.0f} s)",flush=True)
    out={'run':str(a.run),'basis':str(a.basis),'z_edges':list(cfg.xi_z_edges),'z_ref':a.z_ref,
         'pair_weight_by_z':measured.meta['pair_weight_by_z'],'pair_weight_by_type':measured.meta['pair_weight_by_type']}
    cz=measured.counts_z
    # joint fits on all pairs
    joint=evolving_table_for(cz['all'],cfg,basis); out['joint_spline']=summarize(joint)
    print(f"joint (spline): {json.dumps({k:round(joint.params[k],4) for k in PARAMS})} chi2 {joint.params['chi2']:.1f}/{joint.params['cells']} ({time.perf_counter()-t0:.0f} s)",flush=True)
    base=evolving_table_for(cz['all'],cfg.copy(xi_correction='none'),basis); out['joint_base']=summarize(base)
    print(f"joint (base):   {json.dumps({k:round(base.params[k],4) for k in PARAMS if k!='gamma_S'})} chi2 {base.params['chi2']:.1f}",flush=True)
    # non-evolving reference on the collapsed cells (iteration 8/9 table) and the per-type non-evolving fits
    flat=table_for(sl,cfg.copy(xi_z_evolution=False),basis,counts=measured.counts); out['flat_spline']={k:flat.params[k] for k in ('b_F2','beta_F','chi2','cells')}
    print(f"flat (spline):  b_F^2 {flat.params['b_F2']:.4f} beta {flat.params['beta_F']:.3f} chi2 {flat.params['chi2']:.1f}",flush=True)
    out['per_type']={}
    fixed={k:base.params[k] for k in ('gamma_b','gamma_beta')}
    fixed_s={k:joint.params[k] for k in ('gamma_b','gamma_beta','gamma_S')}
    for key in ('AA','AB'):
        n_t,d_t,c_t=cz[key]
        if d_t.sum()<=0: continue
        e={}
        f0=table_for(sl,cfg.copy(xi_correction='none',xi_z_evolution=False),basis,counts=(n_t.sum(axis=0),d_t.sum(axis=0)))
        e['flat_base']={k:f0.params[k] for k in ('b_F2','beta_F','chi2','cells')}
        f1=evolving_table_for((n_t,d_t,c_t),cfg.copy(xi_correction='none'),basis,fixed=fixed); e['evolving_base_fixed_exponents']=summarize(f1)
        f2=evolving_table_for((n_t,d_t,c_t),cfg.copy(xi_correction='none'),basis); e['evolving_base_free']=summarize(f2)
        f3=evolving_table_for((n_t,d_t,c_t),cfg,basis,fixed=fixed_s); e['evolving_spline_fixed_exponents']=summarize(f3)
        out['per_type'][key]=e
        print(f"{key}: flat b_F^2 {f0.params['b_F2']:.4f}; evolving (exponents fixed) b_F^2 {f1.params['b_F2']:.4f} beta {f1.params['beta_F']:.3f}; "
              f"free exponents gamma_b {f2.params['gamma_b']:.2f} gamma_beta {f2.params['gamma_beta']:.2f}; spline-fixed b_F^2 {f3.params['b_F2']:.4f} ({time.perf_counter()-t0:.0f} s)",flush=True)
    if {'AA','AB'}<=set(out['per_type']):
        pt=out['per_type']
        out['AB_over_AA']={'flat':pt['AB']['flat_base']['b_F2']/pt['AA']['flat_base']['b_F2'],
                           'evolving_base_fixed_exponents':pt['AB']['evolving_base_fixed_exponents']['b_F2']/pt['AA']['evolving_base_fixed_exponents']['b_F2'],
                           'evolving_spline_fixed_exponents':pt['AB']['evolving_spline_fixed_exponents']['b_F2']/pt['AA']['evolving_spline_fixed_exponents']['b_F2']}
        print('A x B over A x A amplitude:',json.dumps({k:round(x,4) for k,x in out['AB_over_AA'].items()}),flush=True)
    if a.save_table:
        joint.table.counts=measured.counts; joint.save(a.save_table)
    out['wall_s']=time.perf_counter()-t0
    a.out.parent.mkdir(parents=True,exist_ok=True); a.out.write_text(json.dumps(out,indent=1,default=float)+'\n'); print('saved',a.out)


if __name__=='__main__': main()
