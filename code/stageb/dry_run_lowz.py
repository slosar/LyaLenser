"""Stage B dry run: the pair-template estimator on real DR1 forests against the low-redshift-tracer templates,
on one sky disc. Not a measurement: it checks that every piece runs on the real geometry and reports the
amplitudes with their jackknife errors, the curl, and the injection expectation (bookkeeping) on the real
sightlines. Usage: python dry_run_lowz.py --disc 190 30 12 [--out DIR]
"""
from __future__ import annotations
import argparse, json, sys, time
from pathlib import Path
import numpy as np
import healpy as hp
HERE=Path(__file__).resolve().parent; CODE=HERE.parent
for p in (CODE,CODE/'pipeline',HERE):
    if str(p) not in sys.path: sys.path.insert(0,str(p))
from paths import DATA
from campaign4 import campaign_config, read_xi
from xi_fit import BASIS, fit_model_table
from xi_model import xi_from_data
from pairs import find_pairs
from templates import sphere_band_templates
from amplitude import amplitude
from inject import injection_test
import run_mock_validation as v
from desi_io import read_deltas, save_sightlines


def load_basis(path):
    import h5py
    out={'raw':{k:read_xi(path,f'raw/{k}') for k in BASIS},'projected':{k:read_xi(path,f'projected/{k}') for k in BASIS}}
    with h5py.File(path) as f:
        out['coarse']={k:f[f'coarse/{k}'][()] for k in BASIS}
        if 'cpix' in f: out['cpix']=f['cpix'][()]        # iteration-5 bases only; v2 averages over forest pairs
        out['attrs']=dict(f.attrs)
    return out


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--disc',type=float,nargs=3,default=(190.,30.,12.))
    ap.add_argument('--out',type=Path,default=DATA/'stageb/dry_run'); ap.add_argument('--nside',type=int,default=512)
    ap.add_argument('--lowz',type=Path,default=DATA/'lowz'); ap.add_argument('--basis',type=Path,default=DATA/'stageb/basis_hankel_dr1.h5')
    a=ap.parse_args(); a.out.mkdir(parents=True,exist_ok=True); cfg=campaign_config(1.); t0=time.perf_counter(); log={}
    sl=read_deltas(2.1,3.0,region={'disc':tuple(a.disc)},cfg=cfg); save_sightlines(sl,a.out/'sightlines.h5')
    log['forests']=int(sl.nq); log['pixels']=int(len(sl.chi)); log['median_pixels_per_forest']=float(np.median(np.diff(sl.pix_start)))
    area=np.pi*a.disc[2]**2; log['area_deg2']=area; log['forests_per_deg2']=sl.nq/area
    # Measured xi and the fitted table on the DESI basis (Hankel provider, 0.8 A pixel, resolution).
    basis=load_basis(a.basis); num,den=xi_from_data(sl,cfg).counts
    ft=fit_model_table(num,den,basis['projected'],cfg,basis['coarse']); ft.save(a.out/'xi.h5'); log['xi_fit']=ft.params
    print('xi fit',ft.params['b_F2'],ft.params['beta_F'],f'{time.perf_counter()-t0:.0f} s',flush=True)
    pairs=find_pairs(sl,cfg.r_perp_max/float(sl.chi.min())); cat=v.cat_for(sl,ft.table,cfg,pairs); cat.save(a.out/'catalogue.h5','all')
    log['sightline_pairs']=int(len(cat.a)); print('pairs',len(cat.a),f'{time.perf_counter()-t0:.0f} s',flush=True)
    # Templates from the low-z alm (combined and per slice).
    summary=json.loads((a.lowz/'summary.json').read_text()); names={'combined':a.lowz/'kappa_combined_alm.fits'}
    for s in summary['slices']: names[f"slice_{s['zmin']:g}_{s['zmax']:g}"]=a.lowz/f"kappa_slice_{s['zmin']:g}_{s['zmax']:g}_alm.fits"
    bundles={}
    for name,path in names.items():
        alm=hp.read_alm(str(path)); bundles[name],_=sphere_band_templates(alm,sl.ra,sl.dec,nside=a.nside,source=name)
    reg,nside_jk,nreg=v.midpoint_regions(cat,sl); log['jackknife']={'nside':nside_jk,'regions':nreg}
    log['fits']={}
    for name,b in bundles.items():
        r=amplitude(cat,b,cfg.g1,reg); s=v.common_science(r)
        log['fits'][name]={'A':s['A'],'jk_error':s['jk_error'],'sigma_F':s['sigma_F'],'curl':float(np.mean(r.A[3:6])),'curl_jk_error':float(np.sqrt(np.mean(r.jk_error[3:6]**2))),'bands':r.A.tolist(),'band_errors':r.jk_error.tolist()}
        r.save(a.out/'fits.h5',name)
        print(f"{name:18s} A = {s['A']:8.3f} +- {s['jk_error']:.3f} (sigma_F {s['sigma_F']:.3f}); curl {np.mean(r.A[3:6]):8.3f}",flush=True)
    # Injection expectation with the combined template's own deflection: bookkeeping on the real geometry.
    alpha_inj=sum(t.alpha for t in bundles['combined'][:3])
    exp=injection_test(sl,ft.table,alpha_inj,[-.5,-.25,.25,.5],cfg,templates=bundles['combined'],expectation=True)
    log['injection_expectation']={'paired_slopes_by_amplitude':exp['paired_slopes_by_amplitude'],'paired_slope':exp['paired_slope']}
    print('injection expectation slopes',exp['paired_slopes_by_amplitude'],flush=True)
    log['wall_s']=time.perf_counter()-t0
    (a.out/'dry_run.json').write_text(json.dumps(log,indent=1,default=float)+'\n'); print('done',a.out,f'{log["wall_s"]:.0f} s')


if __name__=='__main__': main()
