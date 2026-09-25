"""Seed-parallel, provenance-checked Stage A iteration-4 campaign.

A completion JSON is the commit marker. Every phase owns a separate directory;
collect never generates a mock. Source or artifact changes are hard errors.
"""
from pathlib import Path
import argparse, dataclasses, hashlib, json, os, resource, shutil, time
import numpy as np
from numba import set_num_threads
import run_mock_validation as v
from config import Config
from mock import generate_mock, save_mock, sightlines_for_variant
from xi_model import XiTable, xi_from_data, xi_from_counts, xi_from_model
from grid_covariance import xi_from_mock_grid, sample_covariance
from validation_stats import absolute_statistics, slope_statistics
from random_streams import STREAM_NAMES, seed_streams
from xi_fit import basis_tables, project_fine, coarse_bin, BASIS
from mock import grid_geometry, patch_side_rad
from paths import MOCKS, ACT_MASK

ITERATION=6
A_GRID=(0,.5,1,2)
# Iteration 6 (GATES v6): previously unexamined realisations. Seeds 0-59, 100-104 and 200-209 of iterations 4-5 are
# development material and are refused by the seed phase.
DEV_SEEDS=tuple(range(3000,3005))
DENSE_SEEDS=tuple(range(2000,2020))
SPARSE_SEEDS=tuple(range(1000,1400))
FULL_SEEDS=tuple(range(1000,1040))      # carry the A grid {0, .5, 1, 2}; the rest carry A in {0, 1}
EXTRAS_SEED=SPARSE_SEEDS[0]             # margin, 100-random-template, injection, flags and benchmark diagnostics
SMOKE_SPARSE=SPARSE_SEEDS[:2]; SMOKE_DENSE=DENSE_SEEDS[:2]
INJECTION_AMPLITUDES=(-2,-1,-.5,-.25,0,.25,.5,1,2)
CONTROL_NAMES=('numerical','injection','flags','random','benchmark')

def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''): h.update(block)
    return h.hexdigest()

def fingerprint(scale):
    # Every python source under code/ (recursive; review 5, finding 7) and the CAMPAIGN configuration, i.e. the
    # frozen numerical choices actually used, not the Config defaults.
    sources=[p for p in v.CODE.rglob('*.py') if '__pycache__' not in p.parts]
    # Path-valued fields (data_root, report_root) are machine-specific and must not enter the fingerprint,
    # otherwise products made on different sites can never be merged; numerically integrated floats (g1) differ
    # between CPUs in the last bits, so floats are rounded to 12 significant digits.
    config={k:(float('%.12g'%val) if isinstance(val,float) else val)
            for k,val in dataclasses.asdict(campaign_config(scale)).items() if not isinstance(val,Path)}
    return {'iteration':ITERATION,'sources':{str(p.relative_to(v.ROOT)):digest(p) for p in sorted(sources)},
            'gates':digest(v.ROOT/'GATES.md'),'scale':float(scale),
            'config':json.loads(json.dumps(config,default=v.serializable)),
            'A_grid':list(A_GRID),'streams':list(STREAM_NAMES),'seeds':{'dev':list(DEV_SEEDS),'dense':list(DENSE_SEEDS),
            'sparse':[SPARSE_SEEDS[0],SPARSE_SEEDS[-1]],'full':[FULL_SEEDS[0],FULL_SEEDS[-1]],'extras':EXTRAS_SEED},
            'injection_amplitudes':list(INJECTION_AMPLITUDES),'deprojected_bound':v.DEPROJECTED_BOUND,
            'numbers3':digest(v.ROOT/'report/numbers3.json'),
            'ACT_mask':digest(ACT_MASK) if ACT_MASK.exists() else None,
            'versions':{name:__import__('importlib.metadata',fromlist=['version']).version(name)
                        for name in ('numpy','scipy','numba','healpy','h5py','camb','fitsio')}}

def completion(directory,expected):
    path=Path(directory)/'complete.json'
    if not path.exists(): return None
    d=json.loads(path.read_text())
    if d['provenance']!=expected: raise RuntimeError(f'provenance mismatch: {path}')
    for name,h in d['artifacts'].items():
        p=Path(directory)/name
        if not p.exists() or digest(p)!=h: raise RuntimeError(f'missing/changed artifact: {p}')
    return d['result']

def begin(directory,provenance):
    """Restart an interrupted phase only under identical inputs, in its own directory."""
    attempt=directory/'attempt.json'
    if attempt.exists() and json.loads(attempt.read_text())!=provenance:
        raise RuntimeError(f'provenance mismatch in interrupted phase: {directory}')
    # Only incomplete products owned by this phase can be removed.
    for p in directory.iterdir():
        if p.is_file() and p.name!='attempt.json': p.unlink()
    v.dump(attempt,provenance)


def finish(directory,provenance,result):
    files={str(p.relative_to(directory)):digest(p) for p in sorted(directory.rglob('*'))
           if p.is_file() and p.name not in ('complete.json',) and not p.name.endswith('.tmp')}
    v.dump(directory/'complete.json',{'provenance':provenance,'artifacts':files,'result':result,
                                           'threads':int(os.environ.get('NUMBA_NUM_THREADS','4')),
                                           'completed_UTC':__import__('datetime').datetime.now(__import__('datetime').timezone.utc).isoformat()})
    return result

def read_xi(path,group='xi'):
    import h5py
    with h5py.File(path) as f:
        g=f[group]
        return XiTable(*(g[k][()] for k in ('r_perp','r_par','xi','xi_rp')),json.loads(g.attrs['meta']),
                       chi_nodes=(g['chi_nodes'][()] if 'chi_nodes' in g else None))

def campaign_config(scale):
    """The frozen numerical choices of iterations 5-6: kernel r_perp in [3, 30], fit range from 3, Kaiser model."""
    return Config(scale=scale,r_perp_min=3.,fit_rperp_min=3.)

def basis_phase(root,scale):
    """Iteration 5: basis spectra of the forest model on the campaign's grid, projected through the continuum
    operator on the campaign's radial pixel grid, plus their 1 Mpc/h binning. One job, reused by every phase."""
    directory=root/'basis'; directory.mkdir(parents=True,exist_ok=True)
    prov=fingerprint(scale); prior=completion(directory,prov)
    if prior is not None: return prior
    begin(directory,prov); start=time.perf_counter(); cfg=campaign_config(scale); geo=grid_geometry(cfg)
    if cfg.forest_model=='kaiser': raw=basis_tables(cfg,'grid',grid_shape=geo['shape'])
    else: raw=basis_tables(cfg,'hankel',nk=cfg.analytic_nk,model=cfg.forest_model,los_pixel=cfg.los_pixel,los_resolution=cfg.los_resolution)
    proj={k:project_fine(raw[k],geo['cpix'],cfg) for k in BASIS}
    coarse={k:coarse_bin(proj[k],cfg) for k in BASIS}
    import h5py
    path=directory/'basis.h5'
    for k in BASIS:
        raw[k].save(path,f'raw/{k}'); proj[k].save(path,f'projected/{k}')
    with h5py.File(path,'a') as f:
        for k in BASIS: f[f'coarse/{k}']=coarse[k]
        f['cpix']=geo['cpix']
    return finish(directory,prov,{'grid':list(geo['shape']),'pixels':int(len(geo['cpix'])),'basis':list(BASIS),
                                  'provider':'grid' if cfg.forest_model=='kaiser' else 'hankel','wall_s':time.perf_counter()-start})

def campaign_basis(root,scale):
    prov=fingerprint(scale)
    if completion(root/'basis',prov) is None: raise RuntimeError('basis phase must complete first')
    import h5py
    path=root/'basis/basis.h5'
    out={'raw':{k:read_xi(path,f'raw/{k}') for k in BASIS},'projected':{k:read_xi(path,f'projected/{k}') for k in BASIS}}
    with h5py.File(path) as f: out['coarse']={k:f[f'coarse/{k}'][()] for k in BASIS}; out['cpix']=f['cpix'][()]
    return out

def frozen(root,scale):
    """Provenance of every downstream product: the fingerprint plus the digests of the freeze AND basis markers
    (the basis marker hashes basis.h5, so a regenerated basis invalidates every seed and control)."""
    current=fingerprint(scale)
    result=completion(root/'freeze',current)
    if result is None: raise RuntimeError('freeze phase must complete first')
    if completion(root/'basis',current) is None: raise RuntimeError('basis phase must complete first')
    basis_marker=digest(root/'basis'/'complete.json')
    if result.get('basis_marker') not in (None,basis_marker): raise RuntimeError('basis changed after the freeze')
    for seed in DEV_SEEDS:
        if completion(root/f'dev/{seed}',current) is None:
            raise RuntimeError(f'missing development seed {seed}')
        if digest(root/f'dev/{seed}/complete.json')!=result['development_markers'][str(seed)]:
            raise RuntimeError(f'changed development seed {seed}')
    return campaign_config(scale),current|{'freeze':digest(root/'freeze'/'complete.json'),'basis':basis_marker}

def dev_seed(root,seed,scale):
    if seed not in DEV_SEEDS: raise ValueError(f'development seeds must be {DEV_SEEDS[0]}--{DEV_SEEDS[-1]}')
    directory=root/f'dev/{seed}'; directory.mkdir(parents=True,exist_ok=True)
    prov=fingerprint(scale); prior=completion(directory,prov)
    if prior is not None: return prior
    if not ACT_MASK.exists(): raise FileNotFoundError(ACT_MASK)
    basis=campaign_basis(root,scale)
    begin(directory,prov)
    start=time.perf_counter(); cfg=campaign_config(scale)
    m=generate_mock(cfg,seed,A_true=1,response=False,real_mask=True,completeness=True,
                    magnification=True,cmb_noise=True,disjoint_selection=True,variant_A_values=A_GRID)
    path=directory/'mock.h5'; save_mock(m,path)
    b=v.make_bundles(m); v.save_bundles(directory/'fits.h5',b)
    pairs=v.find_pairs(m.sightlines,cfg.r_perp_max/m.sightlines.chi.min())
    values=[]; params=[]
    for A in A_GRID:
        sl=sightlines_for_variant(m,A,False); ft=v.table_for(sl,cfg,basis); ft.save(directory/'fits.h5',f'xi/{A}')
        fit=v.fit_save(v.cat_for(sl,ft.table,cfg,pairs),b['truth'],cfg,sl,directory/'fits.h5',f'A{A}')
        values.append(fit['A']); params.append(ft.params)
    from template_audit import audit_mock
    audit=audit_mock(m,directory/'fits.h5',v.make_bundles)
    if seed==DEV_SEEDS[0]:
        raw,cube=xi_from_mock_grid(cfg,m.attrs['grid'],return_cube=True)
        raw.save(directory/'raw.h5')
        _,gd=sample_covariance(cube,raw.r_perp,2.,.5,256)
        take=(raw.r_perp>=10)&(raw.r_perp<=30)
        quadrature=float(np.linalg.norm(raw.xi_rp[take]-gd[take])/np.linalg.norm(gd[take]))
        import h5py
        with h5py.File(directory/'raw.h5','a') as f: f['grid_covariance_cube']=cube
    else: quadrature=None
    return finish(directory,prov,{'seed':seed,'A_values':values,'xi_fit':params,
                                  'template_coefficients':audit,'quadrature':quadrature,'wall_s':time.perf_counter()-start})

def freeze(root,scale):
    prov=fingerprint(scale); directory=root/'freeze'; directory.mkdir(exist_ok=True)
    prior=completion(directory,prov)
    if prior is not None: return prior
    if completion(root/'basis',prov) is None: raise RuntimeError('basis phase must complete first')
    dev=[]
    for seed in DEV_SEEDS:
        d=completion(root/f'dev/{seed}',prov)
        if d is None: raise RuntimeError(f'missing development seed {seed}')
        dev.append(d)
    begin(directory,prov); cfg=campaign_config(scale)
    stats=slope_statistics([d['A_values'] for d in dev],A_GRID)
    if dev[0]['quadrature']>=.01: raise RuntimeError('raw covariance angular quadrature fails')
    result={'seeds':list(DEV_SEEDS),'development_slope':stats,'results':dev,
            'frozen_choices':{'r_perp_min':cfg.r_perp_min,'fit_rperp_min':cfg.fit_rperp_min,'forest_model':cfg.forest_model,
                              'note':'no numerical choice is made from the development seeds in iteration 5; they only check the chain'},
            'development_fit_parameters':[[p for p in d['xi_fit']] for d in dev],
            'development_markers':{str(s):digest(root/f'dev/{s}/complete.json') for s in DEV_SEEDS},
            'basis_marker':digest(root/'basis'/'complete.json')}
    shutil.copyfile(root/f'dev/{DEV_SEEDS[0]}/raw.h5',directory/'raw.h5')
    v.dump(directory/'development.json',result); v.dump(directory/'frozen.json',prov)
    shutil.copyfile(v.ROOT/'GATES.md',directory/'GATES.md')
    return finish(directory,prov,result)

def generator_table(basis):
    """The generator's own forest model (b_F, beta_F of forest_power.FID) on the projected basis: the known-model
    reference kernel against which the fitted table is judged."""
    from forest_power import FID
    b2,beta=FID['b_F']**2,FID['beta_F']; coef={'mu0':b2,'mu2':2*b2*beta,'mu4':b2*beta*beta}
    p=basis['projected']; t=p['mu0']
    return XiTable(t.r_perp,t.r_par,sum(coef[k]*p[k].xi for k in BASIS),sum(coef[k]*p[k].xi_rp for k in BASIS),
                   {'provider':'generator_model_projected','b_F2':b2,'beta_F':beta}),{'b_F2':b2,'beta_F':beta}

def coarse_coefficient(num,den,table,cfg):
    """Response-integral coefficient of the measured coarse table on the fitted table (r_perp^3 weight)."""
    raw=np.divide(num,den,out=np.zeros_like(num),where=den>0); fc=coarse_bin(table,cfg)
    n=raw.shape[0]; rp=np.arange(n)+.5; use=(rp[:,None]>=cfg.fit_rperp_min)&(rp[:,None]<cfg.r_perp_max)&((np.arange(n)+.5)[None,:]<cfg.r_par_max)
    w=rp[:,None]**3*use
    return {'coefficient':float(np.sum(w*raw*fc)/np.sum(w*fc*fc)),'relative_norm':float(np.sqrt(np.sum(w*(raw-fc)**2)/np.sum(w*fc*fc)))}

def dense_seed(directory,seed,cfg,basis):
    from scipy.stats import linregress
    start=time.perf_counter(); result={'seed':seed,'variant':'dense','g':{}}
    gen,gen_params=generator_table(basis); result['generator_model']=gen_params
    for gon in (True,False):
        m=generate_mock(cfg,seed,A_true=1,g_on=gon,response=False,n_los=100,pixel_noise_power=0,
                        disjoint_selection=True,variant_A_values=A_GRID)
        path=directory/f'g{int(gon)}.h5'; save_mock(m,path)
        sl=m.sightlines; c=cfg.copy(g1=cfg.g1 if gon else 0.)
        # Generator pixel size (dx / chi_ref), identical to the deflection the mock was lensed with (review 5, 2).
        tbasis,_=v.flat_sky_band_templates(m.maps['kappa_lya'],sl.ra,sl.dec,
                                        patch_side_rad(m)/m.maps['kappa_lya'].shape[0])
        pairs=v.find_pairs(sl,cfg.r_perp_max/sl.chi.min()); fits=[]; analytic=[]; omitted=[]; comparisons={}; fitparams=[]
        gen.save(path,'xi_generator_model')
        for A in A_GRID:
            sample=sightlines_for_variant(m,A,False); ft=v.table_for(sample,c,basis); ft.save(path,f'xi/{A}'); xi=ft.table; fitparams.append(ft.params)
            cat=v.cat_for(sample,xi,c,pairs); cat.save(path,f'catalogue/{A}')
            fits.append(v.fit_save(cat,tbasis,c,sample,path,f'fits/{A}'))
            analytic.append(v.fit_save(v.cat_for(sample,gen,c,pairs),tbasis,c,sample,path,f'generator/{A}')['A'])
            if gon:
                omitted.append(v.fit_save(cat,tbasis,c.copy(g1=0),sample,path,f'omitted/{A}')['A'])
            if A==0:
                comparisons={'measured_vs_fitted':coarse_coefficient(*xi.counts,xi,c),'measured_vs_generator':coarse_coefficient(*xi.counts,gen,c)}
            if A==1:
                conv=0.
                if gon:
                    from amplitude import pair_scalars
                    from pairs import PairCatalogue
                    d,s,_=pair_scalars(cat,tbasis); dt=d[:3].sum(axis=0)+d[-1]; st=s[:3].sum(axis=0)+s[-1]
                    acc=cat.accum.copy(); acc[:,0]=(acc[:,3]+cfg.g1*acc[:,4])*dt[:,None]+.5*cfg.g1*acc[:,6]*st[:,None]
                    acc[:,1:3]=0; acc[:,8:]=0
                    predcat=PairCatalogue(cat.a,cat.b,cat.thx,cat.thy,cat.theta,acc,cat.npair)
                    moment_prediction=v.fit_save(predcat,tbasis,c.copy(g1=0),sample,path,'moment_prediction')['A']
        aa=[f['A'] for f in fits]
        result['g'][str(gon)]={'fits':fits,'A':aa,'slope':linregress(A_GRID,aa).slope,
                              'generator_A':analytic,'xi_comparison':comparisons,'xi_fit':fitparams,
                              'omitted_A':omitted,'moment_prediction':moment_prediction if gon else None}
        if gon:
            from template_audit import audit_mock
            result['template_coefficients']=audit_mock(m,path,v.make_bundles)
    result.update(wall_s=time.perf_counter()-start,peak_memory_gb=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2)
    return result

def seed_phase(root,seed,variant,scale):
    if variant=='sparse' and seed not in SPARSE_SEEDS: raise ValueError(f'sparse seeds must be {SPARSE_SEEDS[0]}--{SPARSE_SEEDS[-1]}')
    if variant=='dense' and seed not in DENSE_SEEDS: raise ValueError(f'dense seeds must be {DENSE_SEEDS[0]}--{DENSE_SEEDS[-1]}')
    cfg,prov=frozen(root,scale)
    directory=root/f'{variant}/{seed}'; directory.mkdir(parents=True,exist_ok=True)
    prov=prov|{'variant':variant,'seed':seed}; prior=completion(directory,prov)
    if prior is not None: return prior
    basis=campaign_basis(root,scale)
    begin(directory,prov)
    if variant=='dense': result=dense_seed(directory,seed,cfg,basis)
    else:
        raw=read_xi(root/'freeze/raw.h5'); fixed=v.make_bundles(v.load_mock(root/f'dev/{DEV_SEEDS[0]}/mock.h5'))
        # 'full' seeds carry the A grid; every seed carries A in {0, 1} with the response on (GATES v6).
        role='full' if seed in FULL_SEEDS else 'core'
        result=v.process_seed(seed,cfg,directory,role,raw,fixed,basis,extras=(seed==EXTRAS_SEED))
    return finish(directory,prov,result)

def control(root,name,scale):
    cfg,prov=frozen(root,scale)
    if name not in CONTROL_NAMES: raise ValueError(name)
    directory=root/f'controls/{name}'; directory.mkdir(parents=True,exist_ok=True)
    prov=prov|{'control':name}; prior=completion(directory,prov)
    if prior is not None: return prior
    basis=campaign_basis(root,scale)
    begin(directory,prov)
    start=time.perf_counter(); rows=[]; info={}
    def row(name,tol,metrics,ok): v.record(rows,0,name,tol,metrics,ok)
    if name=='numerical':
        model=xi_from_model(v.ForestPower(model='kaiser'),cfg,nk=cfg.analytic_nk)
        doubled=xi_from_model(v.ForestPower(model='kaiser'),cfg.copy(xi_max=30),nk=2*cfg.analytic_nk)
        use=(model.r_perp>=10)&(model.r_perp<=30); ref=doubled.xi_rp[doubled.r_perp>=10]
        change=np.linalg.norm(model.xi_rp[use][:,model.r_par<=30]-ref)/np.linalg.norm(ref)
        model.save(directory/'tables.h5','analytic'); doubled.save(directory/'tables.h5','doubled')
        row('analytic derivative convergence','<1%',{'relative_norm_change':change},change<.01)
        dev=json.loads((root/'freeze/development.json').read_text())
        quad=dev['results'][0]['quadrature']; row('discrete-grid covariance quadrature','<1%',{'relative_norm_change':quad},quad<.01)
        import subprocess,sys
        t=time.perf_counter(); run=subprocess.run([sys.executable,'-m','pytest','-q','tests'],cwd=v.HERE,text=True,capture_output=True)
        (directory/'pytest.log').write_text(run.stdout+'\n'+run.stderr)
        row('unit and regression tests','all pass',{'exit_code':run.returncode,'wall_s':time.perf_counter()-t,'output':run.stdout.strip().splitlines()[-1:]},run.returncode==0)
    else:
        s0=EXTRAS_SEED; seedprov=frozen(root,scale)[1]|{'variant':'sparse','seed':s0}
        d=completion(root/f'sparse/{s0}',seedprov)
        if d is None: raise RuntimeError(f'sparse seed {s0} must complete before controls')
        base=v.load_mock(root/f'sparse/{s0}/seed{s0:03d}.h5'); sl=sightlines_for_variant(base,1,True)
        xi=v.table_for(sl,cfg,basis).table; b=v.make_bundles(base); path=directory/'products.h5'
        if name=='injection':
            # GATES v6. Required: the noise-free expectation of the coordinate injection (delta delta -> xi at the true
            # separation; production selection; exact deflection against the band basis) has odd slope 1 +- 0.05 at
            # the smallest amplitude. Report: its convergence with amplitude, the measured (noisy, one-realisation)
            # odd slopes and their difference from the expectation.
            amps=list(INJECTION_AMPLITUDES)
            exp=v.injection_test(sl,xi,base.alpha_lya,amps,cfg,templates=b['truth'],output=path,expectation=True)
            inj=v.injection_test(sl,xi,base.alpha_lya,amps,cfg,templates=b['truth'],output=path)
            curlerr=float(np.sqrt(np.mean(inj['errors'][:,3:6]**2))); info['injection']=inj; info['injection_expectation']=exp
            small=min(exp['paired_slopes_by_amplitude']); e_small=exp['paired_slopes_by_amplitude'][small]
            row('injection bookkeeping (noise-free expectation, production selection)',f'odd slope at |A| = {small:g} within 1 ± 0.05',
                {'expected_slope_small_amplitude':e_small,'expected_slopes_by_amplitude':exp['paired_slopes_by_amplitude'],'amplitude':small},
                abs(e_small-1)<=.05)
            row('injection convergence and measured slopes','report: expected slope vs |A|; measured odd slopes; measured minus expected; curl',
                {'expected_by_amplitude':exp['paired_slopes_by_amplitude'],'measured_by_amplitude':inj['paired_slopes_by_amplitude'],
                 'measured_combined':inj['paired_slope'],'expected_combined':exp['paired_slope'],
                 'measured_minus_expected_by_amplitude':{k:inj['paired_slopes_by_amplitude'][k]-exp['paired_slopes_by_amplitude'][k] for k in exp['paired_slopes_by_amplitude']},
                 'curl_slope':inj['curl_slope'],'curl_error':curlerr},np.isfinite(inj['paired_slope']))
        elif name=='flags':
            baseline=d['fits']['A1_R1']['deprojected']['A']; info['flags']={}
            for flag in ('magnification','completeness','real_mask'):
                flags=dict(magnification=True,completeness=True,real_mask=True); flags[flag]=False
                m=generate_mock(cfg,s0,A_true=1,response=True,cmb_noise=True,disjoint_selection=True,**flags)
                p=directory/f'{flag}.h5'; save_mock(m,p); bb=v.make_bundles(m); v.save_bundles(p,bb)
                fx=v.table_for(m.sightlines,cfg,basis); fx.save(p); xx=fx.table
                fit=v.fit_save(v.cat_for(m.sightlines,xx,cfg),bb['deprojected'],cfg,m.sightlines,p,'fit')
                info['flags'][flag]={'absolute_A':fit['A'],'baseline_A':baseline,'diagnostic_shift':fit['A']-baseline}
            row('flag and same-realization margin diagnostics','report: persist all (diagnostic)',{'flags':info['flags'],'margins':d['margin']},True)
        elif name=='random':
            rr=d['random']; st=absolute_statistics([r['A'] for r in rr]); sd=np.std([r['A'] for r in rr],ddof=1)
            st.update(scatter_over_rms_jk=sd/np.sqrt(np.mean([r['jk_error']**2 for r in rr])),scatter_over_rms_sigma_F=sd/np.sqrt(np.mean([r['sigma_F']**2 for r in rr])))
            row('100 random-template diagnostic','report: finite statistics; no numerical tolerance',st,len(rr)==100 and all(np.isfinite(st[k]) for k in ('mean','sem','scatter_over_rms_jk','scatter_over_rms_sigma_F')))
        elif name=='benchmark':
            bench=v.benchmark(sl,.1,xi,cfg,directory/'benchmark.json'); peak=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2
            row('benchmark and memory','report: positive throughput; requested threads; peak <40 GB',{'benchmark':bench,'peak_gib':peak},bench['threads']==int(os.environ['NUMBA_NUM_THREADS']) and bench['pixel_pairs_per_s']>0 and peak*1024**3<40e9)
    return finish(directory,prov,info|{'rows':rows,'wall_s':time.perf_counter()-start})


def collect(root,scale,output):
    cfg,prov=frozen(root,scale); smoke=scale<1
    # GATES v6: one fixed ensemble (400 sparse, 20 dense), no data-dependent extension or stopping.
    sparse_ids=list(SMOKE_SPARSE) if smoke else list(SPARSE_SEEDS)
    dense_ids=list(SMOKE_DENSE) if smoke else list(DENSE_SEEDS)
    def read_seeds(ids,variant):
        ds=[]
        for seed in ids:
            d=completion(root/f'{variant}/{seed}',prov|{'variant':variant,'seed':seed})
            if d is None: raise RuntimeError(f'missing {variant} seed {seed}')
            ds.append(d)
        return ds
    sparse=read_seeds(sparse_ids,'sparse')
    decision={'smoke':smoke,'extend':False,'fixed_N':{'sparse':len(sparse_ids),'dense':len(dense_ids)},
              'reason':'GATES v6: fixed ensemble, no extension'}
    dense=read_seeds(dense_ids,'dense'); extras=[]
    for name in CONTROL_NAMES:
        d=completion(root/f'controls/{name}',prov|{'control':name})
        if d is None: raise RuntimeError(f'missing control {name}')
        extras.extend(d['rows'])
    rows=v.build_rows(sparse,{'rows':extras})
    def row(name,tol,st,ok): v.record(rows,0,name,tol,st,ok)
    # Sparse normalization is consistency; dense ensemble supplies precision.
    for r in rows:
        if r['acceptance'].endswith('normalization slope'):
            st=r['measured']; r['tolerance']='|slope-1| <= 2 SEM'; r['pass']=abs(st['residual'])<=2*st['sem']
    template_ok=True
    for margin in (0,150,300):
        for kind in ('continuous','sampled'):
            # Gate in the science band 40 <= L <= 300; the full-resolution pixel regression is a diagnostic.
            st=absolute_statistics([d['template_coefficients'][str(margin)][kind+'_band'] for d in sparse],1)
            ok=abs(st['residual'])<= (.03 if kind=='continuous' else 2*st['sem'])
            row(f'template coefficient {kind} (40<=L<=300), margin {margin}','1 ± .03' if kind=='continuous' else '1 within 2 SEM',st,ok)
            template_ok &= ok
            full=absolute_statistics([d['template_coefficients'][str(margin)][kind] for d in sparse],1)
            row(f'template coefficient {kind} (all pixel scales, diagnostic), margin {margin}','report only',full,np.isfinite(full['mean']))
    # Fitted forest-model parameters of every sparse sample against the generator's values (diagnostic).
    from forest_power import FID
    for key,truth in (('b_F2',FID['b_F']**2),('beta_F',FID['beta_F'])):
        vals=[d['xi_fit']['A0_R0'][key] for d in sparse]; st=absolute_statistics(vals,truth)
        row(f'fitted {key} (sparse A=0 samples) versus generator {truth:.4g}','report only',st,np.isfinite(st['mean']))
    for gon in ('True','False'):
        dd=[d['g'][gon] for d in dense]
        slopes=slope_statistics([d['A'] for d in dd],A_GRID)
        row(f'dense physical g={gon}',f'slope 1 ± .05; >={len(DENSE_SEEDS)} scale-1 seeds',slopes,not smoke and abs(slopes['residual'])<=.05)
        # Independent absolute ensemble errors, never paired-error cancellation.
        a0=absolute_statistics([d['A'][0] for d in dd]); a1=absolute_statistics([d['A'][2] for d in dd],1)
        check={'A0_absolute':a0,'A1_absolute':a1,'difference_of_absolute_means':a1['mean']-a0['mean'],
               'sem_from_absolute_ensembles':float(np.hypot(a0['sem'],a1['sem']))}
        row(f'dense absolute endpoints g={gon}','A0 within 2 SEM of 0 and A1 within 2 SEM of 1',check,
            not smoke and abs(a0['mean'])<=2*a0['sem'] and abs(a1['residual'])<=2*a1['sem'])
        comp=absolute_statistics([d['xi_comparison']['measured_vs_fitted']['coefficient'] for d in dd],1)
        row(f'measured versus fitted table g={gon}','report only (coarse, r_perp^3 weight, fit range)',comp,np.isfinite(comp['mean']))
        for key,truth in (('b_F2',FID['b_F']**2),('beta_F',FID['beta_F'])):
            st=absolute_statistics([d['xi_fit'][0][key] for d in dd],truth)
            row(f'dense fitted {key} versus generator {truth:.4g}, g={gon}','report only',st,np.isfinite(st['mean']))
        analytical=slope_statistics([d['generator_A'] for d in dd],A_GRID)
        ratio=slopes['mean']/analytical['mean']
        row(f'fitted-model versus generator-model normalisation g={gon}','slope ratio within 5%',{'fitted':slopes,'generator':analytical,'ratio':ratio},abs(ratio-1)<=.05)
        if gon=='True':
            omitted=slope_statistics([d['omitted_A'] for d in dd],A_GRID)
            expected=float(np.mean([d['moment_prediction'] for d in dd])); ratio=omitted['mean']/slopes['mean']
            row('first moments','omitted/correct slope ratio agrees with prediction within .05',{'observed_ratio':ratio,'predicted_ratio':expected},abs(ratio-expected)<=.05)
    for r in rows:
        if 'deprojected' in r['acceptance'] and not template_ok:
            r['pass']=False; r['detail']+=' Blocked by prerequisite map-level template gate.'
        r.setdefault('required',not str(r['tolerance']).lower().startswith('report'))
    required=[r for r in rows if r['required']]
    result={'iteration':ITERATION,'smoke':smoke,'scale':scale,'acceptance':rows,'provenance':prov,
            'stopping':decision,'sparse':sparse,'dense':dense,
            'required_passed':sum(r['pass'] for r in required),'required_total':len(required),
            'Stage_B_allowed':not smoke and all(r['pass'] for r in required)}
    v.dump(root/'collection.json',result); render(result,output)
    return result


def render(result,output):
    output=Path(output); output.mkdir(parents=True,exist_ok=True)
    v.dump(output/'mock_validation.json',result)
    label='SMOKE — not acceptance' if result['smoke'] else 'frozen campaign'
    req=result.get('required_total',len(result['acceptance'])); ok=result.get('required_passed',sum(r['pass'] for r in result['acceptance']))
    lines=[f'# Stage A iteration {ITERATION}: {label}','',
           f"Scale {result['scale']}; {len(result['sparse'])} sparse and {len(result['dense'])} dense seeds; {ok}/{req} required gates pass. "
           +("**Stage B remains blocked.**" if not result['Stage_B_allowed'] else '**Stage A accepted.**'),'',
           'All recovery/null statistics are absolute F^-1(q-mf). SEMs and 95% residual bounds use the seed ensemble (every seed carries A_true = 0 and 1, so the null and recovery rows share realisations). Shared/disjoint shifts and cell-offset comparisons are diagnostics. Smoke precision cannot certify scale-1 gates.','',
           '| Gate | Measurement | Frozen tolerance | Type | Result |','|---|---|---|---|---|']
    for r in result['acceptance']:
        s=r['measured']
        metric=f"{s['mean']:.6g} ± {s['sem']:.4g} SEM; bound95 {s['bound95']:.4g}; N={s['n']}" if 'mean' in s else json.dumps(s,default=v.serializable)
        kind='required' if r.get('required',True) else 'report only'
        lines.append(f"| {r['acceptance']} | {metric.replace('|','&#124;')} | {r['tolerance'].replace('|','&#124;')} | {kind} | {'PASS' if r['pass'] else 'FAIL'} |")
    lines+=['','## Template coefficients versus margin','','Band columns are the gated 40<=L<=300 cross/auto coefficients; the other columns are full-resolution pixel regressions (diagnostic).','',
            '| Seed | Margin | Continuous band | Sampled band | Continuous | Sampled | Nominal (no radial normalisation) | Old half-cell offset |','|---|---|---|---|---|---|---|---|']
    for d in result['sparse']:
        for margin,c in d['template_coefficients'].items():
            lines.append(f"| {d['seed']} | {margin} | {c['continuous_band']:.6g} | {c['sampled_band']:.6g} | {c['continuous']:.6g} | {c['sampled']:.6g} | {c['nominal']:.6g} | {c['old_offset']:.6g} |")
    lines+=['','## Shared-selection diagnostic','','| Seed | Absolute disjoint A(1), response on | Absolute shared A(1), response on | Shared minus disjoint |','|---|---|---|---|']
    for d in result['sparse']:
        a=d['fits']['A1_R1']['deprojected']['A']; b=d['shared']['A1_R1']['deprojected']['A']
        lines.append(f"| {d['seed']} | {a:.6g} | {b:.6g} | {b-a:.6g} |")
    lines+=['','## Reconstruction and limits','',
            'Per-seed HDF5 retains maps over each template range, continuum intensity, masks, fitted templates, xi tables, catalogues and q/F/mf partial sums. Completion JSON files contain source/GATES/config provenance and artifact hashes. `collect` validates these before rebuilding this report; it performs no simulation.',
            '', 'See [iteration-4 notes](../code/pipeline/NOTES.md) for the grid-cell diagnosis, derivative comparison approximations and remaining scale-1 checks. Full absolute points, SEMs, predictions and stopping decision are in [JSON](mock_validation.json).']
    (output/'mock_validation.md').write_text('\n'.join(lines)+'\n')
    figs=output/'figures'; figs.mkdir(exist_ok=True)
    fig,ax=v.plt.subplots(figsize=(6,4))
    # Only the A-grid seeds carry A_true in {0, .5, 1, 2}; the grid figure uses those.
    recovery=[d for d in result['sparse'] if d['role'] in v.GRID_ROLES]
    for variant,ds in [('sparse',recovery),('dense',result['dense'])]:
        y=np.array([[d['fits'][f'A{A:g}_R0']['truth']['A'] for A in A_GRID] if variant=='sparse' else d['g']['True']['A'] for d in ds])
        ax.errorbar(A_GRID,y.mean(axis=0),yerr=y.std(axis=0,ddof=1)/np.sqrt(len(y)),marker='o',label=variant)
    ax.plot([0,2],[0,2],'k--'); ax.set(xlabel='Physical A true',ylabel='Absolute fitted A ± SEM',title=label); ax.legend()
    fig.tight_layout(); fig.savefig(figs/'mock_iteration4_normalisation.pdf'); v.plt.close(fig)
    fig,ax=v.plt.subplots(figsize=(6,4))
    for kind in ('continuous_band','sampled_band','continuous','sampled','old_offset'):
        y=np.array([[d['template_coefficients'][str(m)][kind] for m in (0,150,300)] for d in result['sparse']])
        ax.errorbar([0,150,300],y.mean(axis=0),yerr=y.std(axis=0,ddof=1)/np.sqrt(len(y)),marker='o',label=kind)
    ax.axhspan(.97,1.03,alpha=.15,color='grey'); ax.legend(fontsize=8); ax.set(xlabel='Template margin (Mpc/h)',ylabel='Equal-range map coefficient',title=label)
    fig.tight_layout(); fig.savefig(figs/'mock_iteration4_template.pdf'); v.plt.close(fig)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--phase',required=True,choices=['basis','dev-seed','freeze','seed','control','collect','rebuild'])
    ap.add_argument('--seed',type=int); ap.add_argument('--variant',choices=['sparse','dense'],default='sparse')
    ap.add_argument('--name',choices=CONTROL_NAMES); ap.add_argument('--scale',type=float,default=1.)
    ap.add_argument('--mock-root',type=Path,default=MOCKS/f'iteration{ITERATION}'); ap.add_argument('--output',type=Path,default=v.ROOT/'report')
    args=ap.parse_args(); set_num_threads(int(os.environ.get('NUMBA_NUM_THREADS','4')))
    if not 0<args.scale<=1: ap.error('scale must be in (0,1]')
    args.mock_root.mkdir(parents=True,exist_ok=True)
    if args.phase=='basis': result=basis_phase(args.mock_root,args.scale)
    elif args.phase=='dev-seed': result=dev_seed(args.mock_root,args.seed,args.scale)
    elif args.phase=='freeze': result=freeze(args.mock_root,args.scale)
    elif args.phase=='seed': result=seed_phase(args.mock_root,args.seed,args.variant,args.scale)
    elif args.phase=='control': result=control(args.mock_root,args.name,args.scale)
    else: result=collect(args.mock_root,args.scale,args.output)
    print('COMPLETE',args.phase,args.seed,args.variant,flush=True)
