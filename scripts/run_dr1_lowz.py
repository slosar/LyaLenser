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
Refit (iteration 15): --refit-from RUN reuses RUN's sightlines.h5, xi.h5 and catalogue.h5 (hard-linked into --out)
and redoes only the template fits, the injection expectation and the random-template null with the templates of
--lowz; the other options must match RUN's.
"""
from __future__ import annotations
import argparse, json, sys, time, resource
from pathlib import Path
import numpy as np
import healpy as hp
import sys; from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # repository root: `lyalenser` imports without installation
from lyalenser.paths import DATA
from lyalenser.config import production_config
from lyalenser.xi_model import xi_from_data
from lyalenser.pairs import find_pairs, pair_midpoint_regions
from lyalenser.templates import sphere_band_templates, load_templates, derivative_ratio, with_scalar_derivative
from lyalenser.tables import read_xi
from lyalenser.pairs import PairCatalogue
import os
from lyalenser.amplitude import amplitude, curl_amplitude, n_science
from lyalenser.inject import injection_test
from lyalenser.lowz import optimal_combination
from lyalenser.amplitude import common_science
from lyalenser.pairs import cat_for
from lyalenser.tables import table_for
from lyalenser.desi_io import read_deltas, save_sightlines, load_sightlines
from dry_run_lowz import load_basis


def fit(cat,templates,cfg,reg):
    r=amplitude(cat,templates,cfg.g1,reg); s=common_science(r); cu,cue=curl_amplitude(r)
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
    ap.add_argument('--xi-knots',choices=('bicubic','medium','fine'),default='bicubic',help='knot set of the spline correction (fiducial bicubic: 16 + 6 coefficients; medium 30 + 6; fine 42 + 6)')
    ap.add_argument('--z-evolution',dest='zevol',action='store_true',default=True,
                    help='fit the redshift-evolving table (iteration 10, default): power laws in (1+z) for the bias, beta and the correction')
    ap.add_argument('--no-z-evolution',dest='zevol',action='store_false')
    ap.add_argument('--bands',type=float,nargs='+',default=None,help='edges of the science bands in L, e.g. 40 200 400 600 800 1000 (default: templates.SCIENCE_BANDS)')
    ap.add_argument('--rperp-max',type=float,default=30.,help='transverse separation cut of the pair estimator and of the correlation fit (fiducial 30; robustness rows 20 and 40)')
    ap.add_argument('--regions',nargs='+',default=['lya'],choices=['lya','lyb'],
                    help="delta regions to use; 'lya lyb' extends every sightline with its Lyb-region segment "
                         "(A x A and A x B pixel pairs; B x B is dropped)")
    ap.add_argument('--region',type=float,nargs=3,default=None,metavar=('RA','DEC','RADIUS')); ap.add_argument('--seed',type=int,default=2026)
    ap.add_argument('--zmin',type=float,default=2.1,help='forest slab (iteration 11: 1.96, the blue end of the DR1 grid)'); ap.add_argument('--zmax',type=float,default=3.0)
    ap.add_argument('--zeff',type=float,default=None,help='source-plane redshift of the templates (chi_ref); must match the --zref the templates were built with. '
                                                          'Default: the z of chi_ref in Config (2.4). The weighted mean pixel redshift of the sample is always reported.')
    ap.add_argument('--refit-from',type=Path,default=None,help='reuse the sightlines, the fitted table and the pair catalogue of this run; redo the fits only')
    a=ap.parse_args(); a.out.mkdir(parents=True,exist_ok=True)
    from lyalenser.cosmo import chi as chi_of_z, z_of_chi
    cfg=production_config().copy(xi_correction=a.xi_correction,xi_correction_ridge=a.xi_ridge,xi_z_evolution=bool(a.zevol),slabs=((a.zmin,a.zmax),),
                                 xi_z_edges=tuple(sorted({a.zmin,a.zmax}|{z for z in (2.1,2.2,2.3,2.4,2.55,2.75) if a.zmin+0.05<z<a.zmax-0.05})))
    if a.zeff is not None: cfg=cfg.copy(chi_ref=float(chi_of_z(a.zeff)))
    if a.rperp_max!=30.: cfg=cfg.copy(r_perp_max=float(a.rperp_max))
    if a.xi_knots!='bicubic': cfg=cfg.copy(xi_knots=a.xi_knots)
    from lyalenser.templates import SCIENCE_BANDS
    BANDS=tuple((int(a.bands[i]),int(a.bands[i+1])) for i in range(len(a.bands)-1)) if a.bands else SCIENCE_BANDS
    t0=time.perf_counter(); log={'config':{k:(str(v) if isinstance(v,Path) else v) for k,v in vars(cfg).items()},'zmin':a.zmin,'zmax':a.zmax,'z_source_plane':a.zeff,'science_bands':[list(b) for b in BANDS]}
    def stamp(msg): print(f'[{time.perf_counter()-t0:6.0f} s, {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2:.1f} GB] {msg}',flush=True)
    if a.refit_from is not None:
        prev=json.loads((a.refit_from/'dr1_lowz.json').read_text())
        for k,v in (('zmin',a.zmin),('zmax',a.zmax)):
            if abs(float(prev[k])-v)>1e-9: raise SystemExit(f'--refit-from: {k} differs ({prev[k]} vs {v})')
        for k in ('r_perp_max','xi_knots','xi_correction'):
            if prev['config'].get(k)!=getattr(cfg,k): raise SystemExit(f'--refit-from: {k} differs ({prev["config"].get(k)} vs {getattr(cfg,k)})')
        cfg=cfg.copy(chi_ref=float(prev['chi_ref']))
        for f in ('sightlines.h5','xi.h5','catalogue.h5'):
            dst=a.out/f
            if dst.exists() or dst.is_symlink(): dst.unlink()
            try: os.link(a.refit_from/f,dst)
            except OSError: os.symlink((a.refit_from/f).resolve(),dst)
        log={**{k:v for k,v in prev.items() if k not in ('fits','joint','injection_expectation','random_templates','wall_s','peak_gb')},
             'config':{k:(str(v) if isinstance(v,Path) else v) for k,v in vars(cfg).items()},'science_bands':[list(b) for b in BANDS],'refit_from':str(a.refit_from)}
        sl=load_sightlines(a.out/'sightlines.h5'); tab=read_xi(a.out/'xi.h5','xi'); cat=PairCatalogue.load(a.out/'catalogue.h5')
        stamp(f"refit from {a.refit_from}: {sl.nq} forests, {len(cat.a)} sightline pairs; templates {a.lowz}")
        return fits_and_nulls(a,cfg,BANDS,log,sl,tab,cat,stamp,t0)
    region={'disc':tuple(a.region)} if a.region else None
    sl=read_deltas(a.zmin,a.zmax,region=region,cfg=cfg,forest_regions=tuple(a.regions)); save_sightlines(sl,a.out/'sightlines.h5')
    log['forests']=int(sl.nq); log['pixels']=int(len(sl.chi)); log['median_pixels_per_forest']=float(np.median(np.diff(sl.pix_start)))
    zpix=z_of_chi(sl.chi.astype(float)); wpix=sl.w.astype(float)
    log['z_eff_weighted']=float(np.sum(wpix*zpix)/np.sum(wpix)); log['z_mean_unweighted']=float(zpix.mean()); log['chi_ref']=float(cfg.chi_ref); log['z_ref']=float(z_of_chi(cfg.chi_ref))
    stamp(f"weighted mean pixel redshift {log['z_eff_weighted']:.4f} (unweighted {log['z_mean_unweighted']:.4f}); templates' source plane z = {log['z_ref']:.4f}")
    log['forest_regions']=list(a.regions); log['region_B_pixels']=int((sl.region>0).sum())
    pix64=hp.ang2pix(64,sl.ra,sl.dec,lonlat=True); log['area_deg2_nside64']=float(len(np.unique(pix64))*hp.nside2pixarea(64,degrees=True))
    stamp(f"{sl.nq} forests, {len(sl.chi)} pixels, area {log['area_deg2_nside64']:.0f} deg^2")
    basis=load_basis(a.basis); measured=xi_from_data(sl,cfg); counts=measured.counts
    log['pair_weight_by_type']=measured.meta.get('pair_weight_by_type',{}); log['pair_weight_by_z']=measured.meta.get('pair_weight_by_z')
    if a.zevol:
        # Iteration 10: one evolving table for every pair; the per-type fits with the exponents held at the joint
        # values give the A x A / A x B amplitude ratio AT FIXED REDSHIFT (iteration 9's 0.84 was at the pair
        # types' own mean redshifts). Full diagnostics: scripts/xi_zevol_dr1.py.
        from lyalenser.xi_zevol import evolving_table_for, PARAMS
        cz=measured.counts_z
        ft=evolving_table_for(cz['all'],cfg,basis); ft.table.counts=counts; ft.save(a.out/'xi.h5'); log['xi_fit']=ft.params
        stamp(f"xi fit (evolving, {a.xi_correction}) "+' '.join(f"{k} {ft.params[k]:.4f}" for k in PARAMS)+
              f" chi2 {ft.params['chi2']:.0f} / {ft.params['cells']} cells in {len(cfg.xi_z_edges)-1} z bins")
        plain=evolving_table_for(cz['all'],cfg.copy(xi_correction='none'),basis); plain.table.counts=counts
        plain.save(a.out/'xi.h5','xi_uncorrected'); log['xi_fit_uncorrected']=plain.params
        stamp(f"reference evolving base fit chi2 {plain.params['chi2']:.0f} / {plain.params['cells']} cells")
        # base-only fits for the ratio: with the spline correction the base amplitude is degenerate with the
        # correction (both evolve alike), so the per-type b_F^2 would not be comparable
        fixed={k:plain.params[k] for k in ('gamma_b','gamma_beta')}; log['xi_fit_by_pair_type']={}
        for key in ('AA','AB'):
            n_t,d_t,c_t=cz[key]
            if d_t.sum()<=0: continue
            ft_t=evolving_table_for((n_t,d_t,c_t),cfg.copy(xi_correction='none'),basis,fixed=fixed)
            log['xi_fit_by_pair_type'][key]={k:ft_t.params[k] for k in ('b_F2','beta_F','chi2','cells','fixed')}
        if {'AA','AB'}<=set(log['xi_fit_by_pair_type']):
            r=log['xi_fit_by_pair_type']['AB']['b_F2']/log['xi_fit_by_pair_type']['AA']['b_F2']
            log['xi_fit_by_pair_type']['AB_over_AA_amplitude_at_fixed_z']=float(r)
            stamp(f"amplitude ratio of the A x B to the A x A correlation at fixed redshift: {r:.4f}")
    else:
        # the single fitted table is what the estimator applies to every pair; the per-type fits say how different
        # the A x A and A x B correlations actually are (the amplitude ratio is the bias of using one table)
        by_type=getattr(measured,'counts_by_type',None)
        if by_type is not None:
            log['xi_fit_by_pair_type']={}
            for key,(n_t,d_t) in by_type.items():
                if d_t.sum()<=0: continue
                ft_t=table_for(sl,cfg.copy(xi_correction='none'),basis,counts=(n_t,d_t))
                log['xi_fit_by_pair_type'][key]={k:ft_t.params[k] for k in ('b_F2','beta_F','chi2','cells')}
            if {'AA','AB'}<=set(log['xi_fit_by_pair_type']):
                r=log['xi_fit_by_pair_type']['AB']['b_F2']/log['xi_fit_by_pair_type']['AA']['b_F2']
                log['xi_fit_by_pair_type']['AB_over_AA_amplitude']=float(r)
                stamp(f"amplitude ratio of the A x B to the A x A correlation: {r:.4f}")
        ft=table_for(sl,cfg,basis,counts=counts); ft.table.counts=counts
        ft.save(a.out/'xi.h5'); log['xi_fit']=ft.params
        stamp(f"xi fit ({a.xi_correction}) b_F^2 {ft.params['b_F2']:.4f} beta_F {ft.params['beta_F']:.3f} "
              f"chi2 {ft.params['chi2']:.0f} / {ft.params['cells']} cells")
        # the two-parameter fit of the same counts and basis, saved for the residual comparison of the report
        plain=table_for(sl,cfg.copy(xi_correction='none'),basis,counts=counts)
        plain.table.counts=counts; plain.save(a.out/'xi.h5','xi_uncorrected'); log['xi_fit_uncorrected']=plain.params
        stamp(f"reference two-parameter fit chi2 {plain.params['chi2']:.0f} / {plain.params['cells']} cells")
    pairs=find_pairs(sl,cfg.r_perp_max/float(sl.chi.min())); cat=cat_for(sl,ft.table,cfg,pairs); cat.save(a.out/'catalogue.h5','all')
    log['sightline_pairs']=int(len(cat.a)); stamp(f'{len(cat.a)} sightline pairs')
    return fits_and_nulls(a,cfg,BANDS,log,sl,ft.table,cat,stamp,t0)


def fits_and_nulls(a,cfg,BANDS,log,sl,tab,cat,stamp,t0):
    """Step 5: the template fits (per slice and combined, with the derivative maps), the jackknife combination of the
    slices, the injection expectation and the random-template null; the summary files."""
    reg=pair_midpoint_regions(cat,sl,a.nside_jk); log['jackknife']={'nside':a.nside_jk,'regions':int(len(np.unique(reg)))}
    summary=json.loads((a.lowz/'summary.json').read_text()); names=['combined']+[f"slice_{s['zmin']:g}_{s['zmax']:g}" for s in summary['slices']]
    log['fits']={}; jks=[]; slice_names=[]; log['derivative_ratio']={}
    for name in names:
        b,_=load_templates(a.lowz,name,sl.ra,sl.dec,nside=a.nside_alpha,science_bands=BANDS)
        r,s=fit(cat,b,cfg,reg); r.save(a.out/'fits.h5',name); jk=s.pop('jk'); log['fits'][name]=s; log['derivative_ratio'][name]=derivative_ratio(b)
        if name!='combined': jks.append(jk); slice_names.append(name)
        else: combined_templates=b
        stamp(f"{name:18s} A = {s['A']:8.3f} +- {s['jk_error']:.3f} (sigma_F {s['sigma_F']:.3f}); curl {s['curl']:8.3f} +- {s['curl_jk_error']:.3f}; derivative ratio {log['derivative_ratio'][name]:.2e}")
    log['joint']=optimal_combination([log['fits'][n]['A'] for n in slice_names],np.asarray(jks)); log['joint']['slices']=slice_names
    stamp(f"jackknife combination of the slices A = {log['joint']['A']:.3f} +- {log['joint']['error']:.3f}")
    # Injection expectation with the combined template's science deflection (bookkeeping on the real geometry).
    alpha_inj=sum(t.alpha for t in combined_templates if getattr(t,'kind','')=='signal')
    exp=injection_test(sl,tab,alpha_inj,[-.5,-.25,.25,.5],cfg,templates=combined_templates,expectation=True)
    log['injection_expectation']={k:exp[k] for k in ('paired_slopes_by_amplitude','paired_slope','paired_slope_jk_error')}
    log['injection_expectation']['region_labels_preserved']=True
    stamp(f"injection expectation slopes {exp['paired_slopes_by_amplitude']}")
    # Random-template null: Gaussian realisations of the combined map's spectrum on the tracer mask, with the
    # combined template's effective derivative ratio as their (scalar) derivative map.
    comb_alm=hp.read_alm(str(a.lowz/'kappa_combined_alm.fits')); lmax=hp.Alm.getlmax(len(comb_alm)); g_eff=log['derivative_ratio']['combined']
    mask=hp.read_map(str(a.lowz/f'mask_combined_nside{a.nside}.fits')); cl=hp.alm2cl(comb_alm)/max(float(np.mean(mask**2)),1e-30)
    rng=np.random.default_rng(a.seed); rand=[]; state=np.random.get_state()
    try:
        for i in range(a.randoms):
            np.random.seed(int(rng.integers(2**31))); m=hp.synfast(cl,a.nside,lmax=lmax,verbose=False)*mask
            b,_=sphere_band_templates(hp.map2alm(m,lmax=lmax,iter=0),sl.ra,sl.dec,nside=a.nside_alpha,science_bands=BANDS,source=f'random {i}')
            _,s=fit(cat,with_scalar_derivative(b,g_eff),cfg,reg); s.pop('jk'); rand.append(s)
            if i%10==0: stamp(f'random {i}: A = {s["A"]:.3f} +- {s["jk_error"]:.3f}')
    finally: np.random.set_state(state)
    if rand:
        A=np.array([r['A'] for r in rand]); log['random_templates']={'n':len(A),'mean':float(A.mean()),'sem':float(A.std(ddof=1)/np.sqrt(len(A))) if len(A)>1 else float('nan'),'scatter':float(A.std(ddof=1)) if len(A)>1 else float('nan'),
                                                                        'rms_jk_error':float(np.sqrt(np.mean([r['jk_error']**2 for r in rand]))),'amplitudes':A.tolist(),'derivative_ratio':g_eff}
    else: log['random_templates']={'n':0,'mean':float('nan'),'sem':float('nan'),'scatter':float('nan'),'rms_jk_error':float('nan'),'amplitudes':[]}
    log['wall_s']=time.perf_counter()-t0; log['peak_gb']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2
    (a.out/'dr1_lowz.json').write_text(json.dumps(log,indent=1,default=lambda x: str(x) if isinstance(x,Path) else float(x))+'\n')
    c=log['fits']['combined']; rt=log['random_templates']
    lines=['# DR1 Lya forest lensing x low-redshift tracers (single slab 2.1 < z < 3.0)','',
           f"{log['forests']} forests, {log['sightline_pairs']} sightline pairs, {log['area_deg2_nside64']:.0f} deg^2; xi fit b_F^2 = {log['xi_fit']['b_F2']:.4f}, beta_F = {log['xi_fit']['beta_F']:.3f}"
           +(f", gamma_b = {log['xi_fit']['gamma_b']:.2f}, gamma_beta = {log['xi_fit']['gamma_beta']:.2f}, gamma_S = {log['xi_fit']['gamma_S']:.2f} (z_ref {log['xi_fit']['z_ref']})" if 'gamma_b' in log['xi_fit'] else '')
           +f"; jackknife nside {a.nside_jk} ({log['jackknife']['regions']} regions).",'',
           '| Template | A | jackknife error | sigma_F | curl |','|---|---|---|---|---|']
    for n,s in log['fits'].items(): lines.append(f"| {n} | {s['A']:.3f} | {s['jk_error']:.3f} | {s['sigma_F']:.3f} | {s['curl']:.3f} +- {s['curl_jk_error']:.3f} |")
    lines+=['',f"Jackknife-covariance combination of the slices: A = {log['joint']['A']:.3f} +- {log['joint']['error']:.3f}.",
            f"Random-template null ({rt['n']} Gaussian realisations of the combined map's spectrum): mean {rt['mean']:.3f} +- {rt['sem']:.3f}, scatter {rt['scatter']:.3f} against the RMS jackknife error {rt['rms_jk_error']:.3f}.",
            f"Injection expectation (bookkeeping) odd slopes: {log['injection_expectation']['paired_slopes_by_amplitude']}.",
            f"Wall {log['wall_s']/60:.0f} min, peak {log['peak_gb']:.1f} GB."]
    lines.insert(2,f"Source-distance treatment: derivative maps of the templates ({a.lowz}); effective ratios "+', '.join(f"{n} {g:.2e}" for n,g in log['derivative_ratio'].items())+" per Mpc/h.")
    (a.out/'dr1_lowz.md').write_text('\n'.join(lines)+'\n'); stamp('done')


if __name__=='__main__': main()
