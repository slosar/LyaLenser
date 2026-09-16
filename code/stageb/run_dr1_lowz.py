"""Stage B, single slab: DR1 forests x low-redshift-tracer templates over the whole footprint.

Steps: all DR1 forests (2.1 <= z <= 3.0, picca inverse-variance weights) -> fitted xi table (DESI-pixel Hankel
basis) -> pair catalogue (3 <= r_perp <= 30, r_par < 30 Mpc/h) -> deflection templates from the per-slice and
combined Wiener-filtered alm (lowz_catalogues.py) -> amplitudes with the seven-component basis (three science
bands, curl partners, junk) and a HEALPix jackknife -> per-slice and combined A with jackknife errors, the
jackknife-covariance combination of the slices, the curl null, the injection expectation (bookkeeping) with the
combined template's own deflection, and a random-template null: Gaussian realisations of the combined map's
measured spectrum through the tracer mask, fitted exactly like the data (mean and scatter vs the jackknife error).
Everything is written to --out (HDF5 products, dr1_lowz.json, dr1_lowz.md). No interpretation here.
Usage: python run_dr1_lowz.py --out DIR [--lowz DIR] [--randoms 40] [--nside-jk 8]
"""
from __future__ import annotations
import argparse, json, sys, time, resource
from pathlib import Path
import numpy as np
import healpy as hp
HERE=Path(__file__).resolve().parent; CODE=HERE.parent
for p in (CODE,CODE/'pipeline',HERE):
    if str(p) not in sys.path: sys.path.insert(0,str(p))
from paths import DATA
from campaign4 import campaign_config
from xi_model import xi_from_data
from pairs import find_pairs, pair_midpoint_regions
from templates import sphere_band_templates
from amplitude import amplitude, curl_amplitude, n_science
from inject import injection_test
from lowz import optimal_combination
import run_mock_validation as v
from desi_io import read_deltas, save_sightlines
from dry_run_lowz import load_basis


def fit(cat,templates,cfg,reg):
    r=amplitude(cat,templates,cfg.g1,reg); s=v.common_science(r); cu,cue=curl_amplitude(r)
    return r,{'A':s['A'],'jk_error':s['jk_error'],'sigma_F':s['sigma_F'],'curl':cu,
              'curl_jk_error':cue,'bands':r.A.tolist(),'band_errors':r.jk_error.tolist(),'jk':np.asarray(s['jk'])}


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--out',type=Path,default=DATA/'stageb/dr1_lowz')
    ap.add_argument('--lowz',type=Path,default=DATA/'lowz_split'); ap.add_argument('--basis',type=Path,default=DATA/'stageb/basis_hankel_dr1.h5')
    ap.add_argument('--nside',type=int,default=512,help='nside of the tracer maps and of the mask files')
    ap.add_argument('--nside-alpha',type=int,default=1024,
                    help='nside at which the deflection is evaluated; must resolve the top science band')
    ap.add_argument('--nside-jk',type=int,default=8); ap.add_argument('--randoms',type=int,default=40)
    ap.add_argument('--xi-correction',choices=('none','spline'),default='spline',
                    help="'none' = iteration-5 two-parameter Kaiser fit; 'spline' = iteration-8 corrected table")
    ap.add_argument('--xi-ridge',type=float,default=1e-2)
    ap.add_argument('--regions',nargs='+',default=['lya'],choices=['lya','lyb'],
                    help="delta regions to use; 'lya lyb' extends every sightline with its Lyb-region segment "
                         "(A x A and A x B pixel pairs; B x B is dropped)")
    ap.add_argument('--region',type=float,nargs=3,default=None,metavar=('RA','DEC','RADIUS')); ap.add_argument('--seed',type=int,default=2026)
    a=ap.parse_args(); a.out.mkdir(parents=True,exist_ok=True)
    cfg=campaign_config(1.).copy(xi_correction=a.xi_correction,xi_correction_ridge=a.xi_ridge); t0=time.perf_counter(); log={'config':{k:(str(v) if isinstance(v,Path) else v) for k,v in vars(cfg).items()}}
    def stamp(msg): print(f'[{time.perf_counter()-t0:6.0f} s, {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2:.1f} GB] {msg}',flush=True)
    region={'disc':tuple(a.region)} if a.region else None
    sl=read_deltas(2.1,3.0,region=region,cfg=cfg,forest_regions=tuple(a.regions)); save_sightlines(sl,a.out/'sightlines.h5')
    log['forests']=int(sl.nq); log['pixels']=int(len(sl.chi)); log['median_pixels_per_forest']=float(np.median(np.diff(sl.pix_start)))
    log['forest_regions']=list(a.regions); log['region_B_pixels']=int((sl.region>0).sum())
    pix64=hp.ang2pix(64,sl.ra,sl.dec,lonlat=True); log['area_deg2_nside64']=float(len(np.unique(pix64))*hp.nside2pixarea(64,degrees=True))
    stamp(f"{sl.nq} forests, {len(sl.chi)} pixels, area {log['area_deg2_nside64']:.0f} deg^2")
    basis=load_basis(a.basis); measured=xi_from_data(sl,cfg); counts=measured.counts
    log['pair_weight_by_type']=measured.meta.get('pair_weight_by_type',{})
    # the single fitted table is what the estimator applies to every pair; the per-type fits say how different
    # the A x A and A x B correlations actually are (the amplitude ratio is the bias of using one table)
    by_type=getattr(measured,'counts_by_type',None)
    if by_type is not None:
        log['xi_fit_by_pair_type']={}
        for key,(n_t,d_t) in by_type.items():
            if d_t.sum()<=0: continue
            ft_t=v.table_for(sl,cfg.copy(xi_correction='none'),basis,counts=(n_t,d_t))
            log['xi_fit_by_pair_type'][key]={k:ft_t.params[k] for k in ('b_F2','beta_F','chi2','cells')}
        if {'AA','AB'}<=set(log['xi_fit_by_pair_type']):
            r=log['xi_fit_by_pair_type']['AB']['b_F2']/log['xi_fit_by_pair_type']['AA']['b_F2']
            log['xi_fit_by_pair_type']['AB_over_AA_amplitude']=float(r)
            stamp(f"amplitude ratio of the A x B to the A x A correlation: {r:.4f}")
    ft=v.table_for(sl,cfg,basis,counts=counts); ft.table.counts=counts
    ft.save(a.out/'xi.h5'); log['xi_fit']=ft.params
    stamp(f"xi fit ({a.xi_correction}) b_F^2 {ft.params['b_F2']:.4f} beta_F {ft.params['beta_F']:.3f} "
          f"chi2 {ft.params['chi2']:.0f} / {ft.params['cells']} cells")
    # the two-parameter fit of the same counts and basis, saved for the residual comparison of the report
    plain=v.table_for(sl,cfg.copy(xi_correction='none'),basis,counts=counts)
    plain.table.counts=counts; plain.save(a.out/'xi.h5','xi_uncorrected'); log['xi_fit_uncorrected']=plain.params
    stamp(f"reference two-parameter fit chi2 {plain.params['chi2']:.0f} / {plain.params['cells']} cells")
    pairs=find_pairs(sl,cfg.r_perp_max/float(sl.chi.min())); cat=v.cat_for(sl,ft.table,cfg,pairs); cat.save(a.out/'catalogue.h5','all')
    log['sightline_pairs']=int(len(cat.a)); stamp(f'{len(cat.a)} sightline pairs')
    reg=pair_midpoint_regions(cat,sl,a.nside_jk); log['jackknife']={'nside':a.nside_jk,'regions':int(len(np.unique(reg)))}
    summary=json.loads((a.lowz/'summary.json').read_text()); names={'combined':a.lowz/'kappa_combined_alm.fits'}
    for s in summary['slices']: names[f"slice_{s['zmin']:g}_{s['zmax']:g}"]=a.lowz/f"kappa_slice_{s['zmin']:g}_{s['zmax']:g}_alm.fits"
    log['fits']={}; jks=[]; slice_names=[]
    for name,path in names.items():
        alm=hp.read_alm(str(path)); b,_=sphere_band_templates(alm,sl.ra,sl.dec,nside=a.nside_alpha,source=name)
        r,s=fit(cat,b,cfg,reg); r.save(a.out/'fits.h5',name); jk=s.pop('jk'); log['fits'][name]=s
        if name!='combined': jks.append(jk); slice_names.append(name)
        else: combined_templates=b
        stamp(f"{name:18s} A = {s['A']:8.3f} +- {s['jk_error']:.3f} (sigma_F {s['sigma_F']:.3f}); curl {s['curl']:8.3f} +- {s['curl_jk_error']:.3f}")
    log['joint']=optimal_combination([log['fits'][n]['A'] for n in slice_names],np.asarray(jks)); log['joint']['slices']=slice_names
    stamp(f"jackknife combination of the slices A = {log['joint']['A']:.3f} +- {log['joint']['error']:.3f}")
    # Injection expectation with the combined template's science deflection (bookkeeping on the real geometry).
    alpha_inj=sum(t.alpha for t in combined_templates if getattr(t,'kind','')=='signal')
    exp=injection_test(sl,ft.table,alpha_inj,[-.5,-.25,.25,.5],cfg,templates=combined_templates,expectation=True)
    log['injection_expectation']={'paired_slopes_by_amplitude':exp['paired_slopes_by_amplitude'],'paired_slope':exp['paired_slope']}
    stamp(f"injection expectation slopes {exp['paired_slopes_by_amplitude']}")
    # Random-template null: Gaussian realisations of the combined map's spectrum on the tracer mask.
    comb_alm=hp.read_alm(str(names['combined'])); lmax=hp.Alm.getlmax(len(comb_alm))
    mask=hp.read_map(str(a.lowz/f'mask_combined_nside{a.nside}.fits')); cl=hp.alm2cl(comb_alm)/max(float(np.mean(mask**2)),1e-30)
    rng=np.random.default_rng(a.seed); rand=[]; state=np.random.get_state()
    try:
        for i in range(a.randoms):
            np.random.seed(int(rng.integers(2**31))); m=hp.synfast(cl,a.nside,lmax=lmax,verbose=False)*mask
            b,_=sphere_band_templates(hp.map2alm(m,lmax=lmax,iter=0),sl.ra,sl.dec,nside=a.nside_alpha,source=f'random {i}')
            _,s=fit(cat,b,cfg,reg); s.pop('jk'); rand.append(s)
            if i%10==0: stamp(f'random {i}: A = {s["A"]:.3f} +- {s["jk_error"]:.3f}')
    finally: np.random.set_state(state)
    A=np.array([r['A'] for r in rand]); log['random_templates']={'n':len(A),'mean':float(A.mean()),'sem':float(A.std(ddof=1)/np.sqrt(len(A))),'scatter':float(A.std(ddof=1)),
                                                                    'rms_jk_error':float(np.sqrt(np.mean([r['jk_error']**2 for r in rand]))),'amplitudes':A.tolist()}
    log['wall_s']=time.perf_counter()-t0; log['peak_gb']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2
    (a.out/'dr1_lowz.json').write_text(json.dumps(log,indent=1,default=lambda x: str(x) if isinstance(x,Path) else float(x))+'\n')
    c=log['fits']['combined']; rt=log['random_templates']
    lines=['# DR1 Lya forest lensing x low-redshift tracers (single slab 2.1 < z < 3.0)','',
           f"{log['forests']} forests, {log['sightline_pairs']} sightline pairs, {log['area_deg2_nside64']:.0f} deg^2; xi fit b_F^2 = {log['xi_fit']['b_F2']:.4f}, beta_F = {log['xi_fit']['beta_F']:.3f}; jackknife nside {a.nside_jk} ({log['jackknife']['regions']} regions).",'',
           '| Template | A | jackknife error | sigma_F | curl |','|---|---|---|---|---|']
    for n,s in log['fits'].items(): lines.append(f"| {n} | {s['A']:.3f} | {s['jk_error']:.3f} | {s['sigma_F']:.3f} | {s['curl']:.3f} +- {s['curl_jk_error']:.3f} |")
    lines+=['',f"Jackknife-covariance combination of the slices: A = {log['joint']['A']:.3f} +- {log['joint']['error']:.3f}.",
            f"Random-template null ({rt['n']} Gaussian realisations of the combined map's spectrum): mean {rt['mean']:.3f} +- {rt['sem']:.3f}, scatter {rt['scatter']:.3f} against the RMS jackknife error {rt['rms_jk_error']:.3f}.",
            f"Injection expectation (bookkeeping) odd slopes: {log['injection_expectation']['paired_slopes_by_amplitude']}.",
            f"Wall {log['wall_s']/60:.0f} min, peak {log['peak_gb']:.1f} GB."]
    (a.out/'dr1_lowz.md').write_text('\n'.join(lines)+'\n'); stamp('done')


if __name__=='__main__': main()
