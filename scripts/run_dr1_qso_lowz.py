"""Stage B iteration 12: DR1 quasar x forest cross-correlation lensing, single slab, low-redshift templates.

Steps: forests (the iteration-11 sightline set, 1.96 <= z <= 3.0, both regions) + the DR1 quasar catalogue
(1.9 < z < 3.1) -> quasar-pixel cells in signed (r_perp, r_par) by pair redshift -> evolving cross model fit with
the forest side fixed from the auto fit (`xi_fit_uncorrected` of the auto run), free quasar bias and its
evolution, redshift offset and smoothing, spline correction -> quasar-sightline pair catalogue and accumulation
-> amplitudes against the templates evaluated at the concatenated (sightline, quasar) positions (combined, per
slice, per band), jackknife on the pair midpoints, curl, random templates, optional redshift sub-slabs.
Writes --out/{xi_qf.h5, catalogue.h5, fits.h5, dr1_qso.json, dr1_qso.md} and the QA data for the plots.
Usage: python run_dr1_qso_lowz.py --out DIR --auto-run $LYALENSER_DATA/stageb/dr1_lowz_v7 --lowz $LYALENSER_DATA/lowz_v3 --basis .../basis_cross_dr1_z196.h5
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
from lyalenser.tables import read_xi
from lyalenser.cosmo import chi as chi_of_z, z_of_chi
from lyalenser.xi_fit import BASIS
from lyalenser.xi_cross import find_cross_pairs, xi_qf_from_data, fit_evolving_cross, cross_correction_for, accumulate_cross, Positions
from lyalenser.pairs import pair_midpoint_regions
from lyalenser.templates import sphere_band_templates, load_templates, derivative_ratio, with_scalar_derivative
from lyalenser.pairs import PairCatalogue
import os
from lyalenser.amplitude import amplitude, curl_amplitude
from lyalenser.lowz import optimal_combination
from lyalenser.amplitude import common_science
from lyalenser.desi_io import load_sightlines
from lyalenser.desi_io import read_deltas, in_region
from lyalenser.qso_io import read_quasars, QuasarSet


def fit(cat,templates,cfg,reg):
    r=amplitude(cat,templates,cfg.g1,reg); s=common_science(r); cu,cue=curl_amplitude(r)
    return r,{'A':s['A'],'jk_error':s['jk_error'],'sigma_F':s['sigma_F'],'curl':cu,'curl_jk_error':cue,'bands':r.A.tolist(),'band_errors':r.jk_error.tolist(),'jk':np.asarray(s['jk'])}


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--out',type=Path,default=DATA/'stageb/dr1_qso_v1')
    ap.add_argument('--auto-run',type=Path,default=DATA/'stageb/dr1_lowz_v7',help='auto-correlation run: sightlines.h5 and the forest fit parameters')
    ap.add_argument('--lowz',type=Path,default=DATA/'lowz_v3'); ap.add_argument('--basis',type=Path,default=DATA/'stageb/basis_cross_dr1_z196.h5')
    ap.add_argument('--nside',type=int,default=512); ap.add_argument('--nside-alpha',type=int,default=1024); ap.add_argument('--nside-jk',type=int,default=8)
    ap.add_argument('--randoms',type=int,default=40); ap.add_argument('--xi-ridge',type=float,default=1e-2); ap.add_argument('--seed',type=int,default=2026)
    ap.add_argument('--zmin',type=float,default=1.96); ap.add_argument('--zmax',type=float,default=3.0); ap.add_argument('--zeff',type=float,default=2.3476)
    ap.add_argument('--qzmin',type=float,default=1.9); ap.add_argument('--qzmax',type=float,default=3.1)
    ap.add_argument('--region',type=float,nargs=3,default=None,metavar=('RA','DEC','RADIUS'),help='disc test: read the deltas in the disc instead of the saved sightlines')
    ap.add_argument('--sub-slabs',type=float,nargs='*',default=[1.96,2.25,2.55,3.0],help='edges of the redshift sub-slabs (pair mean redshift) for the split')
    ap.add_argument('--no-spline',action='store_true')
    ap.add_argument('--xi-knots',choices=('bicubic','medium','fine'),default='bicubic',help='knot set of the cross spline correction (fiducial bicubic: 20 coefficients; medium 42; fine 63)')
    ap.add_argument('--spline-fixed-base',action='store_true',help='fit the base model first and hold b_q, gamma_q, dr_par, sigma_par at those values in the spline fit (the correction is then purely residual)')
    ap.add_argument('--bands',type=float,nargs='+',default=None,help='edges of the science bands in L, e.g. 40 200 400 600 800 1000 (default: templates.SCIENCE_BANDS)')
    ap.add_argument('--rperp-max',type=float,default=30.,help='transverse separation cut of the pair estimator and of the correlation fit (fiducial 30; robustness rows 20 and 40)')
    ap.add_argument('--refit-from',type=Path,default=None,help='reuse the fitted cross table and the pair catalogue of this run (hard-linked into --out); redo the fits only (iteration 15)')
    a=ap.parse_args(); a.out.mkdir(parents=True,exist_ok=True); t0=time.perf_counter()
    for f in ('xi_qf.h5','catalogue.h5','fits.h5'):          # products of an earlier (possibly interrupted) run
        if a.refit_from is not None and a.out.resolve()==a.refit_from.resolve(): raise SystemExit('--refit-from must differ from --out')
        if (a.out/f).exists() or (a.out/f).is_symlink(): (a.out/f).unlink()
    cfg=production_config().copy(xi_correction='none' if a.no_spline else 'spline',xi_correction_ridge=a.xi_ridge,xi_z_evolution=True,slabs=((a.zmin,a.zmax),),
                                 xi_z_edges=tuple(sorted({a.zmin,a.zmax}|{z for z in (2.1,2.2,2.3,2.4,2.55,2.75) if a.zmin+0.05<z<a.zmax-0.05})),chi_ref=float(chi_of_z(a.zeff)))
    from lyalenser.templates import SCIENCE_BANDS
    BANDS=tuple((int(a.bands[i]),int(a.bands[i+1])) for i in range(len(a.bands)-1)) if a.bands else SCIENCE_BANDS
    if a.rperp_max!=30.: cfg=cfg.copy(r_perp_max=float(a.rperp_max))
    if a.xi_knots!='bicubic': cfg=cfg.copy(xi_knots=a.xi_knots)
    log={'config':{k:(str(x) if isinstance(x,Path) else x) for k,x in vars(cfg).items()},'zmin':a.zmin,'zmax':a.zmax,'z_source_plane':a.zeff,'quasar_z':[a.qzmin,a.qzmax],'science_bands':[list(b) for b in BANDS]}
    def stamp(msg): print(f'[{time.perf_counter()-t0:6.0f} s, {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2:.1f} GB] {msg}',flush=True)
    if a.refit_from is not None:
        prev=json.loads((a.refit_from/'dr1_qso.json').read_text())
        for k in ('r_perp_max','xi_knots','xi_correction'):
            if prev['config'].get(k)!=getattr(cfg,k): raise SystemExit(f'--refit-from: {k} differs ({prev["config"].get(k)} vs {getattr(cfg,k)})')
        if abs(float(prev['chi_ref'] if 'chi_ref' in prev else prev['config']['chi_ref'])-cfg.chi_ref)>1e-6: raise SystemExit('--refit-from: chi_ref differs')
        for f in ('xi_qf.h5','catalogue.h5'):
            try: os.link(a.refit_from/f,a.out/f)
            except OSError: os.symlink((a.refit_from/f).resolve(),a.out/f)
        log={**{k:v for k,v in prev.items() if k not in ('fits','joint','sub_slabs','random_templates','wall_s','peak_gb')},'config':log['config'],'science_bands':log['science_bands'],'refit_from':str(a.refit_from)}
        sl=load_sightlines(a.auto_run/'sightlines.h5'); qso=read_quasars(a.qzmin,a.qzmax)
        tab=read_xi(a.out/'xi_qf.h5','xi'); par=log['xi_fit']; cat=PairCatalogue.load(a.out/'catalogue.h5')
        pairs=find_cross_pairs(sl,qso,cfg.r_perp_max/max(1.,float(min(sl.chi.min(),qso.chi.min())))) if (a.sub_slabs and len(a.sub_slabs)>2) else None
        stamp(f"refit from {a.refit_from}: {sl.nq} forests, {qso.n} quasars, {len(cat.a)} quasar-sightline pairs; templates {a.lowz}")
        return fits_and_nulls(a,cfg,BANDS,log,sl,qso,tab,par,cat,pairs,stamp,t0)
    # ---- forests and quasars
    if a.region:
        sl=read_deltas(a.zmin,a.zmax,region={'disc':tuple(a.region)},cfg=cfg,forest_regions=('lya','lyb'))
    else: sl=load_sightlines(a.auto_run/'sightlines.h5')
    qso=read_quasars(a.qzmin,a.qzmax)
    if a.region:
        ra0,dec0,rad=a.region; keep=in_region(qso.ra,qso.dec,{'disc':(ra0,dec0,rad+1.)}); qso=QuasarSet(qso.qid[keep],qso.ra[keep],qso.dec[keep],qso.z[keep],qso.chi[keep],qso.attrs)
    log['forests']=int(sl.nq); log['pixels']=int(len(sl.chi)); log['quasars']=int(qso.n); log['quasars_with_forest']=int(np.isin(qso.qid,sl.qid).sum())
    stamp(f'{sl.nq} forests, {len(sl.chi)} pixels; {qso.n} quasars ({log["quasars_with_forest"]} with a forest in the sample)')
    # ---- forest parameters from the auto fit (base model: b_F^2, beta_F and their power laws)
    auto=json.loads((a.auto_run/'dr1_lowz.json').read_text()); u=auto['xi_fit_uncorrected']
    forest={k:float(u[k]) for k in ('b_F2','beta_F','gamma_b','gamma_beta')}; log['forest_fixed']=forest
    # ---- measured cells and the fit
    basis={'projected':{k:read_xi(a.basis,f'projected/{k}') for k in BASIS}}
    theta_max=cfg.xi_max/max(1.,float(min(sl.chi.min(),qso.chi.min()))); xpairs=find_cross_pairs(sl,qso,theta_max)
    meas=xi_qf_from_data(sl,qso,cfg,xpairs); log['xi_pairs_quasar_sightline']=meas['n_pairs']; log['pair_weight_by_z']=meas['pair_weight_by_z']; log['pair_weight_by_region']=meas['pair_weight_by_region']
    stamp(f"{meas['n_pairs']} quasar-sightline pairs within {cfg.xi_max} Mpc/h; cells measured")
    base,parb,_=fit_evolving_cross(meas['all'],basis['projected'],cfg,forest,None)
    base.save(a.out/'xi_qf.h5','xi_uncorrected'); log['xi_fit_uncorrected']=parb
    stamp(f"cross base fit: b_q {parb['b_q']:.3f} gamma_q {parb['gamma_q']:.2f} dr_par {parb['dr_par']:.2f} sigma_par {parb['sigma_par']:.2f} chi2 {parb['chi2']:.0f}")
    corr=None if a.no_spline else cross_correction_for(basis['projected'],forest,cfg)
    fixed={k:parb[k] for k in ('b_q','gamma_q','dr_par','sigma_par')} if (a.spline_fixed_base and not a.no_spline) else None
    tab,par,dbg=fit_evolving_cross(meas['all'],basis['projected'],cfg,forest,corr,ridge=(0. if a.no_spline else a.xi_ridge),fixed=fixed)
    tab.save(a.out/'xi_qf.h5'); log['xi_fit']=par; log['spline_fixed_base']=bool(fixed)
    stamp(f"cross fit{' (base fixed)' if fixed else ''}: b_q {par['b_q']:.3f} gamma_q {par['gamma_q']:.2f} dr_par {par['dr_par']:.2f} sigma_par {par['sigma_par']:.2f} gamma_S {par['gamma_S']:.2f} chi2 {par['chi2']:.0f} / {par['cells']} cells; beta_q(z_ref) {par['beta_q_at_zref']:.3f}")
    # QA data
    import h5py
    with h5py.File(a.out/'xi_qf.h5','a') as f:
        g=f.require_group('cells'); [g.__delitem__(k) for k in list(g)]
        g['num']=meas['all'][0]; g['den']=meas['all'][1]; g['chisum']=meas['all'][2]; g['model_cells']=dbg['model_cells']; g['use']=dbg['use']; g['z_edges']=meas['z_edges']
        g['numA']=meas['A'][0]; g['denA']=meas['A'][1]; g['numB']=meas['B'][0]; g['denB']=meas['B'][1]
    # ---- pairs for the estimator (r_perp <= 30) and the accumulation
    pairs=find_cross_pairs(sl,qso,cfg.r_perp_max/max(1.,float(min(sl.chi.min(),qso.chi.min()))))
    cat=accumulate_cross(sl,qso,pairs,tab,cfg); cat.save(a.out/'catalogue.h5','all'); log['sightline_pairs']=int(len(cat.a)); stamp(f'{len(cat.a)} quasar-sightline pairs accumulated')
    return fits_and_nulls(a,cfg,BANDS,log,sl,qso,tab,par,cat,pairs,stamp,t0)


def fits_and_nulls(a,cfg,BANDS,log,sl,qso,tab,par,cat,pairs,stamp,t0):
    """The template fits (per slice and combined, with the derivative maps), the jackknife combination of the slices,
    the sub-slab split and the random-template null; the summary files."""
    pos=Positions(sl,qso); reg=pair_midpoint_regions(cat,pos,a.nside_jk); log['jackknife']={'nside':a.nside_jk,'regions':int(len(np.unique(reg)))}
    summary=json.loads((a.lowz/'summary.json').read_text()); names=['combined']+[f"slice_{s['zmin']:g}_{s['zmax']:g}" for s in summary['slices']]
    log['fits']={}; jks=[]; slice_names=[]; templates={}; log['derivative_ratio']={}
    for name in names:
        b,_=load_templates(a.lowz,name,pos.ra,pos.dec,nside=a.nside_alpha,science_bands=BANDS); templates[name]=b
        r,s=fit(cat,b,cfg,reg); r.save(a.out/'fits.h5',name); jk=s.pop('jk'); log['fits'][name]=s; log['derivative_ratio'][name]=derivative_ratio(b)
        if name!='combined': jks.append(jk); slice_names.append(name)
        stamp(f"{name:18s} A = {s['A']:8.3f} +- {s['jk_error']:.3f} (sigma_F {s['sigma_F']:.3f}); curl {s['curl']:8.3f} +- {s['curl_jk_error']:.3f}; derivative ratio {log['derivative_ratio'][name]:.2e}")
    log['joint']=optimal_combination([log['fits'][n]['A'] for n in slice_names],np.asarray(jks)); log['joint']['slices']=slice_names
    stamp(f"jackknife combination of the slices A = {log['joint']['A']:.3f} +- {log['joint']['error']:.3f}")
    # ---- redshift sub-slabs (pair mean redshift), combined template
    log['sub_slabs']={}
    if a.sub_slabs and len(a.sub_slabs)>2:
        edges=list(a.sub_slabs)
        for k in range(len(edges)-1):
            c2=cfg.copy(slabs=tuple((edges[i],edges[i+1]) for i in range(len(edges)-1)),slab_index=k)
            catk=accumulate_cross(sl,qso,pairs,tab,c2); regk=pair_midpoint_regions(catk,pos,a.nside_jk)
            rk,sk=fit(catk,templates['combined'],c2,regk); rk.save(a.out/'fits.h5',f'subslab_{edges[k]:g}_{edges[k+1]:g}'); sk.pop('jk')
            log['sub_slabs'][f'{edges[k]:g}-{edges[k+1]:g}']=sk; stamp(f"sub-slab {edges[k]:g}-{edges[k+1]:g}: A = {sk['A']:.3f} +- {sk['jk_error']:.3f} (sigma_F {sk['sigma_F']:.3f})")
    # ---- random-template null
    comb_alm=hp.read_alm(str(a.lowz/'kappa_combined_alm.fits')); lmax=hp.Alm.getlmax(len(comb_alm)); mask=hp.read_map(str(a.lowz/f'mask_combined_nside{a.nside}.fits')); cl=hp.alm2cl(comb_alm)/max(float(np.mean(mask**2)),1e-30)
    g_eff=log['derivative_ratio']['combined']; rng=np.random.default_rng(a.seed); rand=[]; state=np.random.get_state()
    try:
        for i in range(a.randoms):
            np.random.seed(int(rng.integers(2**31))); m=hp.synfast(cl,a.nside,lmax=lmax,verbose=False)*mask
            b,_=sphere_band_templates(hp.map2alm(m,lmax=lmax,iter=0),pos.ra,pos.dec,nside=a.nside_alpha,science_bands=BANDS,source=f'random {i}'); _,s=fit(cat,with_scalar_derivative(b,g_eff),cfg,reg); s.pop('jk'); rand.append(s)
            if i%10==0: stamp(f'random {i}: A = {s["A"]:.3f} +- {s["jk_error"]:.3f}')
    finally: np.random.set_state(state)
    if rand:
        A=np.array([r['A'] for r in rand]); log['random_templates']={'n':len(A),'mean':float(A.mean()),'sem':float(A.std(ddof=1)/np.sqrt(len(A))),'scatter':float(A.std(ddof=1)),'rms_jk_error':float(np.sqrt(np.mean([r['jk_error']**2 for r in rand]))),'amplitudes':A.tolist(),'derivative_ratio':g_eff}
    log['wall_s']=time.perf_counter()-t0; log['peak_gb']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2
    (a.out/'dr1_qso.json').write_text(json.dumps(log,indent=1,default=lambda x: str(x) if isinstance(x,Path) else float(x))+'\n')
    c=log['fits']['combined']; lines=['# DR1 quasar x Lya forest lensing x low-redshift tracers','',
        f"{log['forests']} forests, {log['quasars']} quasars, {log['sightline_pairs']} quasar-sightline pairs; cross fit b_q = {par['b_q']:.3f} (gamma_q {par['gamma_q']:.2f}), dr_par = {par['dr_par']:.2f}, sigma_par = {par['sigma_par']:.2f} Mpc/h, chi2 {par['chi2']:.0f}/{par['cells']}.",'',
        '| Template | A | jackknife error | sigma_F | curl |','|---|---|---|---|---|']
    for n,s in log['fits'].items(): lines.append(f"| {n} | {s['A']:.3f} | {s['jk_error']:.3f} | {s['sigma_F']:.3f} | {s['curl']:.3f} +- {s['curl_jk_error']:.3f} |")
    for n,s in log['sub_slabs'].items(): lines.append(f"| sub-slab {n} | {s['A']:.3f} | {s['jk_error']:.3f} | {s['sigma_F']:.3f} | {s['curl']:.3f} +- {s['curl_jk_error']:.3f} |")
    if rand: rt=log['random_templates']; lines+=['',f"Random-template null ({rt['n']}): mean {rt['mean']:.3f} +- {rt['sem']:.3f}, scatter {rt['scatter']:.3f} vs RMS jackknife {rt['rms_jk_error']:.3f}."]
    lines.insert(2,f"Source-distance treatment: derivative maps of the templates ({a.lowz}); effective ratios "+', '.join(f"{n} {g:.2e}" for n,g in log['derivative_ratio'].items())+" per Mpc/h.")
    lines.append(f"Wall {log['wall_s']/60:.0f} min, peak {log['peak_gb']:.1f} GB."); (a.out/'dr1_qso.md').write_text('\n'.join(lines)+'\n'); stamp('done')


if __name__=='__main__': main()
