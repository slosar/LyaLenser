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
from grid_covariance import xi_from_mock_grid, sample_covariance, projected_grid_table, derivative_comparison
from validation_stats import absolute_statistics, slope_statistics
from random_streams import STREAM_NAMES, seed_streams
from paths import MOCKS, ACT_MASK

A_GRID=(0,.5,1,2)
DEV_SEEDS=tuple(range(100,105))
DENSE_SEEDS=tuple(range(200,210))
WIDTHS=(.76,.9,1.1,1.3,1.5,1.7,2.)
CONTROL_NAMES=('numerical','injection','flags','random','benchmark')

def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''): h.update(block)
    return h.hexdigest()

def fingerprint(scale):
    sources=list(v.HERE.glob('*.py'))+list((v.HERE/'tests').glob('*.py'))+list(v.CODE.glob('*.py'))
    return {'sources':{str(p.relative_to(v.ROOT)):digest(p) for p in sorted(sources)},
            'gates':digest(v.ROOT/'GATES.md'),'scale':float(scale),
            'config':json.loads(json.dumps(dataclasses.asdict(Config(scale=scale)),default=v.serializable)),
            'A_grid':list(A_GRID),'streams':list(STREAM_NAMES),
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
        return XiTable(*(g[k][()] for k in ('r_perp','r_par','xi','xi_rp')),json.loads(g.attrs['meta']))

def frozen(root,scale):
    current=fingerprint(scale)
    result=completion(root/'freeze',current)
    if result is None: raise RuntimeError('freeze phase must complete first')
    for seed in DEV_SEEDS:
        if completion(root/f'dev/{seed}',current) is None:
            raise RuntimeError(f'missing development seed {seed}')
        if digest(root/f'dev/{seed}/complete.json')!=result['development_markers'][str(seed)]:
            raise RuntimeError(f'changed development seed {seed}')
    return Config(scale=scale,xi_smoothing=result['xi_smoothing']),current|{'freeze':digest(root/'freeze'/'complete.json')}

def dev_seed(root,seed,scale):
    if seed not in DEV_SEEDS: raise ValueError('development seeds must be 100--104')
    directory=root/f'dev/{seed}'; directory.mkdir(parents=True,exist_ok=True)
    prov=fingerprint(scale); prior=completion(directory,prov)
    if prior is not None: return prior
    if not ACT_MASK.exists(): raise FileNotFoundError(ACT_MASK)
    begin(directory,prov)
    start=time.perf_counter(); cfg=Config(scale=scale)
    m=generate_mock(cfg,seed,A_true=1,response=False,real_mask=True,completeness=True,
                    magnification=True,cmb_noise=True,disjoint_selection=True,variant_A_values=A_GRID)
    path=directory/'mock.h5'; save_mock(m,path)
    b=v.make_bundles(m); v.save_bundles(directory/'fits.h5',b)
    pairs=v.find_pairs(m.sightlines,cfg.r_perp_max/m.sightlines.chi.min())
    values={str(w):[] for w in WIDTHS}; convergence={}
    for A in A_GRID:
        sl=sightlines_for_variant(m,A,False); counts=xi_from_data(sl,cfg)
        counts.save(directory/'fits.h5',f'counts/{A}')
        for width in WIDTHS:
            c=cfg.copy(xi_smoothing=width); xi=xi_from_counts(*counts.counts,c)
            fit=v.fit_save(v.cat_for(sl,xi,c,pairs),b['truth'],c,sl,directory/'fits.h5',f'width/{width}/A{A}')
            values[str(width)].append(fit['A'])
            if A==1:
                fine=xi_from_counts(*counts.counts,c.copy(xi_step=.125))
                af=v.fit_bundle(v.cat_for(sl,fine,c,pairs),b['truth'],c.g1,sl)['A']
                convergence[str(width)]=abs(af-fit['A'])/max(abs(af),1e-30)
    from template_audit import audit_mock
    audit=audit_mock(m,directory/'fits.h5',v.make_bundles)
    if seed==100:
        raw,cube=xi_from_mock_grid(cfg,m.attrs['grid'],return_cube=True)
        raw.save(directory/'raw.h5')
        _,gd=sample_covariance(cube,raw.r_perp,2.,.5,256)
        take=(raw.r_perp>=10)&(raw.r_perp<=30)
        quadrature=float(np.linalg.norm(raw.xi_rp[take]-gd[take])/np.linalg.norm(gd[take]))
        import h5py
        with h5py.File(directory/'raw.h5','a') as f: f['grid_covariance_cube']=cube
    else: quadrature=None
    return finish(directory,prov,{'seed':seed,'widths':values,'convergence':convergence,
                                  'template_coefficients':audit,'quadrature':quadrature,'wall_s':time.perf_counter()-start})

def freeze(root,scale):
    prov=fingerprint(scale); directory=root/'freeze'; directory.mkdir(exist_ok=True)
    prior=completion(directory,prov)
    if prior is not None: return prior
    dev=[]
    for seed in DEV_SEEDS:
        d=completion(root/f'dev/{seed}',prov)
        if d is None: raise RuntimeError(f'missing development seed {seed}')
        dev.append(d)
    begin(directory,prov)
    stats={str(w):slope_statistics([d['widths'][str(w)] for d in dev],A_GRID) for w in WIDTHS}
    chosen=min(stats,key=lambda w:abs(stats[w]['mean']-1))
    if dev[0]['quadrature']>=.01: raise RuntimeError('raw covariance angular quadrature fails')
    result={'seeds':list(DEV_SEEDS),'xi_smoothing':float(chosen),'slopes':stats,'results':dev,
            'choice_rule':'minimum absolute development mean-slope residual; ties in declared width order',
            'development_markers':{str(s):digest(root/f'dev/{s}/complete.json') for s in DEV_SEEDS}}
    shutil.copyfile(root/'dev/100/raw.h5',directory/'raw.h5')
    v.dump(directory/'development.json',result); v.dump(directory/'frozen.json',prov)
    shutil.copyfile(v.ROOT/'GATES.md',directory/'GATES.md')
    return finish(directory,prov,result)

def dense_seed(directory,seed,cfg):
    from scipy.stats import linregress
    start=time.perf_counter(); result={'seed':seed,'variant':'dense','g':{}}
    raw=read_xi(directory.parents[1]/'freeze/raw.h5')
    model=xi_from_model(v.ForestPower(model='kaiser'),cfg,nk=cfg.analytic_nk)
    for gon in (True,False):
        m=generate_mock(cfg,seed,A_true=1,g_on=gon,response=False,n_los=100,pixel_noise_power=0,
                        disjoint_selection=True,variant_A_values=A_GRID)
        path=directory/f'g{int(gon)}.h5'; save_mock(m,path)
        sl=m.sightlines; c=cfg.copy(g1=cfg.g1 if gon else 0.)
        basis,_=v.flat_sky_band_templates(m.maps['kappa_lya'],sl.ra,sl.dec,
                                        np.deg2rad(20*cfg.scale)/m.maps['kappa_lya'].shape[0])
        pairs=v.find_pairs(sl,cfg.r_perp_max/sl.chi.min()); fits=[]; analytic=[]; omitted=[]; comparisons={}
        predicted=projected_grid_table(raw,sl,cfg); predicted.save(path,'xi_grid_projected')
        raw.save(path,'xi_grid_raw'); model.save(path,'xi_continuum_model')
        for A in A_GRID:
            sample=sightlines_for_variant(m,A,False); xi=xi_from_data(sample,cfg); xi.save(path,f'xi/{A}')
            cat=v.cat_for(sample,xi,c,pairs); cat.save(path,f'catalogue/{A}')
            fits.append(v.fit_save(cat,basis,c,sample,path,f'fits/{A}'))
            analytic.append(v.fit_save(v.cat_for(sample,model,c,pairs),basis,c,sample,path,f'analytic/{A}')['A'])
            if gon:
                omitted.append(v.fit_save(cat,basis,c.copy(g1=0),sample,path,f'omitted/{A}')['A'])
            if A==0:
                comparisons={'measured_vs_projected_grid':derivative_comparison(xi,predicted),
                             'measured_vs_raw_grid':derivative_comparison(xi,raw),
                             'raw_grid_vs_continuum_model':derivative_comparison(raw,model)}
            if A==1:
                fine=xi_from_data(sample,c.copy(xi_step=.125)); fine.save(path,'xi_fine')
                af=v.fit_save(v.cat_for(sample,fine,c,pairs),basis,c,sample,path,'fit_fine')['A']
                conv=abs(af-fits[-1]['A'])/max(abs(af),1e-30)
                if gon:
                    from amplitude import pair_scalars
                    from pairs import PairCatalogue
                    d,s,_=pair_scalars(cat,basis); dt=d[:3].sum(axis=0)+d[-1]; st=s[:3].sum(axis=0)+s[-1]
                    acc=cat.accum.copy(); acc[:,0]=(acc[:,3]+cfg.g1*acc[:,4])*dt[:,None]+.5*cfg.g1*acc[:,6]*st[:,None]
                    acc[:,1:3]=0; acc[:,8:]=0
                    predcat=PairCatalogue(cat.a,cat.b,cat.thx,cat.thy,cat.theta,acc,cat.npair)
                    moment_prediction=v.fit_save(predcat,basis,c.copy(g1=0),sample,path,'moment_prediction')['A']
        aa=[f['A'] for f in fits]
        result['g'][str(gon)]={'fits':fits,'A':aa,'slope':linregress(A_GRID,aa).slope,
                              'analytic_A':analytic,'xi_comparison':comparisons,'interpolation_change':conv,
                              'omitted_A':omitted,'moment_prediction':moment_prediction if gon else None}
        if gon:
            from template_audit import audit_mock
            result['template_coefficients']=audit_mock(m,path,v.make_bundles)
    result.update(wall_s=time.perf_counter()-start,peak_memory_gb=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2)
    return result

def seed_phase(root,seed,variant,scale):
    cfg,prov=frozen(root,scale)
    if variant=='sparse' and seed not in range(60): raise ValueError('sparse seeds must be 0--59')
    if variant=='dense' and seed not in DENSE_SEEDS: raise ValueError('dense seeds must be 200--209')
    directory=root/f'{variant}/{seed}'; directory.mkdir(parents=True,exist_ok=True)
    prov=prov|{'variant':variant,'seed':seed}; prior=completion(directory,prov)
    if prior is not None: return prior
    begin(directory,prov)
    if variant=='dense': result=dense_seed(directory,seed,cfg)
    else:
        raw=read_xi(root/'freeze/raw.h5'); fixed=v.make_bundles(v.load_mock(root/'dev/100/mock.h5'))
        # Smoke seeds exercise both roles. Full campaign retains the prospective split.
        role='extension' if scale<1 or seed>=40 else ('recovery' if seed<20 else 'null')
        result=v.process_seed(seed,cfg,directory,role,raw,fixed)
    return finish(directory,prov,result)

def control(root,name,scale):
    cfg,prov=frozen(root,scale)
    if name not in CONTROL_NAMES: raise ValueError(name)
    directory=root/f'controls/{name}'; directory.mkdir(parents=True,exist_ok=True)
    prov=prov|{'control':name}; prior=completion(directory,prov)
    if prior is not None: return prior
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
        seedprov=frozen(root,scale)[1]|{'variant':'sparse','seed':0}
        d=completion(root/'sparse/0',seedprov)
        if d is None: raise RuntimeError('sparse seed 0 must complete before controls')
        base=v.load_mock(root/'sparse/0/seed000.h5'); sl=sightlines_for_variant(base,1,True)
        xi=xi_from_data(sl,cfg); b=v.make_bundles(base); path=directory/'products.h5'
        if name=='injection':
            inj=v.injection_test(sl,xi,base.alpha_lya,[-2,-1,0,1,2],cfg,templates=b['truth'],output=path)
            curlerr=float(np.sqrt(np.mean(inj['errors'][:,3:6]**2))); info['injection']=inj
            row('injection bookkeeping','science slope within .05; curl within JK',{'slope':inj['paired_slope'],'curl_slope':inj['curl_slope'],'curl_error':curlerr},abs(inj['paired_slope']-1)<=.05 and abs(inj['curl_slope'])<=curlerr)
        elif name=='flags':
            baseline=d['fits']['A1_R1']['deprojected']['A']; info['flags']={}
            for flag in ('magnification','completeness','real_mask'):
                flags=dict(magnification=True,completeness=True,real_mask=True); flags[flag]=False
                m=generate_mock(cfg,0,A_true=1,response=True,cmb_noise=True,disjoint_selection=True,**flags)
                p=directory/f'{flag}.h5'; save_mock(m,p); bb=v.make_bundles(m); v.save_bundles(p,bb)
                xx=xi_from_data(m.sightlines,cfg); xx.save(p)
                fit=v.fit_save(v.cat_for(m.sightlines,xx,cfg),bb['deprojected'],cfg,m.sightlines,p,'fit')
                info['flags'][flag]={'absolute_A':fit['A'],'baseline_A':baseline,'diagnostic_shift':fit['A']-baseline}
            row('flag and same-realization margin diagnostics','persist all (diagnostic)',{'flags':info['flags'],'margins':d['margin']},True)
        elif name=='random':
            rr=d['random']; st=absolute_statistics([r['A'] for r in rr]); sd=np.std([r['A'] for r in rr],ddof=1)
            st.update(scatter_over_rms_jk=sd/np.sqrt(np.mean([r['jk_error']**2 for r in rr])),scatter_over_rms_sigma_F=sd/np.sqrt(np.mean([r['sigma_F']**2 for r in rr])))
            row('100 random-template diagnostic','finite statistics; no numerical tolerance',st,len(rr)==100 and all(np.isfinite(st[k]) for k in ('mean','sem','scatter_over_rms_jk','scatter_over_rms_sigma_F')))
        elif name=='benchmark':
            bench=v.benchmark(sl,.1,xi,cfg,directory/'benchmark.json'); peak=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2
            row('benchmark and memory','positive throughput; requested threads; peak <40 GB',{'benchmark':bench,'peak_gib':peak},bench['threads']==int(os.environ['NUMBA_NUM_THREADS']) and bench['pixel_pairs_per_s']>0 and peak*1024**3<40e9)
    return finish(directory,prov,info|{'rows':rows,'wall_s':time.perf_counter()-start})


def collect(root,scale,output):
    cfg,prov=frozen(root,scale); smoke=scale<1
    sparse_ids=list(range(2)) if smoke else list(range(40))
    dense_ids=[200,201] if smoke else list(DENSE_SEEDS)
    def read_seeds(ids,variant):
        ds=[]
        for seed in ids:
            d=completion(root/f'{variant}/{seed}',prov|{'variant':variant,'seed':seed})
            if d is None: raise RuntimeError(f'missing {variant} seed {seed}')
            ds.append(d)
        return ds
    sparse=read_seeds(sparse_ids,'sparse')
    if smoke:
        decision={'smoke':True,'extend':False,'reason':'smoke tests phase chaining only; no acceptance stopping decision'}
    else:
        rec=absolute_statistics([d['fits']['A1_R1']['deprojected']['A'] for d in sparse[:20]],1)
        null=absolute_statistics([d['fits']['A0_R1']['deprojected']['A'] for d in sparse[20:]])
        extend=rec['bound95']>.3 or null['bound95']>.3
        decision={'recovery':rec,'null':null,'extend':extend,'seeds':list(range(40,60)) if extend else []}
        v.dump(root/'extension_decision.json',decision)
        if extend: sparse+=read_seeds(list(range(40,60)),'sparse')
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
            ok=abs(st['residual'])<= (.03 if kind=='continuous' else st['sem'])
            row(f'template coefficient {kind} (40<=L<=300), margin {margin}','1 ± .03' if kind=='continuous' else '1 within SEM',st,ok)
            template_ok &= ok
            full=absolute_statistics([d['template_coefficients'][str(margin)][kind] for d in sparse],1)
            row(f'template coefficient {kind} (all pixel scales, diagnostic), margin {margin}','report only',full,np.isfinite(full['mean']))
    for gon in ('True','False'):
        dd=[d['g'][gon] for d in dense]
        slopes=slope_statistics([d['A'] for d in dd],A_GRID)
        row(f'dense physical g={gon}','slope 1 ± .03; >=10 scale-1 seeds',slopes,not smoke and abs(slopes['residual'])<=.03)
        # Independent absolute ensemble errors, never paired-error cancellation.
        a0=absolute_statistics([d['A'][0] for d in dd]); a1=absolute_statistics([d['A'][2] for d in dd],1)
        check={'A0_absolute':a0,'A1_absolute':a1,'difference_of_absolute_means':a1['mean']-a0['mean'],
               'sem_from_absolute_ensembles':float(np.hypot(a0['sem'],a1['sem']))}
        row(f'dense absolute endpoints g={gon}','A0=0 and A1=1 within .03; difference of absolute means 1 ± .03',check,
            not smoke and abs(a0['mean'])<=.03 and abs(a1['residual'])<=.03 and abs(check['difference_of_absolute_means']-1)<=.03)
        comp=absolute_statistics([d['xi_comparison']['measured_vs_projected_grid']['coefficient'] for d in dd],1)
        row(f'measured versus grid-predicted derivative g={gon}','response integral residual <3%',comp,abs(comp['residual'])<.03)
        changes=[d['interpolation_change'] for d in dd]
        row(f'xi interpolation g={gon}','all absolute amplitude changes <.005',{'changes':changes},max(changes)<.005)
        analytical=slope_statistics([d['analytic_A'] for d in dd],A_GRID)
        ratio=slopes['mean']/analytical['mean']
        row(f'data versus analytic baseline g={gon}','slope ratio within 3%',{'data':slopes,'analytic':analytical,'ratio':ratio},abs(ratio-1)<=.03)
        if gon=='True':
            omitted=slope_statistics([d['omitted_A'] for d in dd],A_GRID)
            expected=float(np.mean([d['moment_prediction'] for d in dd])); ratio=omitted['mean']/slopes['mean']
            row('first moments','omitted/correct slope ratio agrees with prediction within .05',{'observed_ratio':ratio,'predicted_ratio':expected},abs(ratio-expected)<=.05)
    for r in rows:
        if 'deprojected' in r['acceptance'] and not template_ok:
            r['pass']=False; r['detail']+=' Blocked by prerequisite map-level template gate.'
    result={'iteration':4,'smoke':smoke,'scale':scale,'acceptance':rows,'provenance':prov,
            'stopping':decision,'sparse':sparse,'dense':dense,'Stage_B_allowed':not smoke and all(r['pass'] for r in rows)}
    v.dump(root/'collection.json',result); render(result,output)
    return result


def render(result,output):
    output=Path(output); output.mkdir(parents=True,exist_ok=True)
    v.dump(output/'mock_validation.json',result)
    label='SMOKE — not acceptance' if result['smoke'] else 'frozen campaign'
    lines=[f'# Stage A iteration 4: {label}','',
           f"Scale {result['scale']}; {len(result['sparse'])} sparse and {len(result['dense'])} dense seeds. **Stage B remains blocked.**" if not result['Stage_B_allowed'] else '**Stage A accepted.**','',
           'All recovery/null statistics are absolute F^-1(q-mf). SEMs and 95% residual bounds use independent seed ensembles. Shared/disjoint shifts and cell-offset comparisons are diagnostics. Smoke precision cannot certify scale-1 gates.','',
           '| Gate | Measurement | Frozen tolerance | Result |','|---|---|---|---|']
    for r in result['acceptance']:
        s=r['measured']
        metric=f"{s['mean']:.6g} ± {s['sem']:.4g} SEM; bound95 {s['bound95']:.4g}; N={s['n']}" if 'mean' in s else json.dumps(s,default=v.serializable)
        lines.append(f"| {r['acceptance']} | {metric.replace('|','&#124;')} | {r['tolerance'].replace('|','&#124;')} | {'PASS' if r['pass'] else 'FAIL'} |")
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
    # Null-role seeds carry only A_true in {0,1}; the grid figure uses recovery-role seeds.
    recovery=[d for d in result['sparse'] if d['role'] in ('recovery','extension')]
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
    ap.add_argument('--phase',required=True,choices=['dev-seed','freeze','seed','control','collect','rebuild'])
    ap.add_argument('--seed',type=int); ap.add_argument('--variant',choices=['sparse','dense'],default='sparse')
    ap.add_argument('--name',choices=CONTROL_NAMES); ap.add_argument('--scale',type=float,default=1.)
    ap.add_argument('--mock-root',type=Path,default=MOCKS/'iteration4'); ap.add_argument('--output',type=Path,default=v.ROOT/'report')
    args=ap.parse_args(); set_num_threads(int(os.environ.get('NUMBA_NUM_THREADS','4')))
    if not 0<args.scale<=1: ap.error('scale must be in (0,1]')
    args.mock_root.mkdir(parents=True,exist_ok=True)
    if args.phase=='dev-seed': result=dev_seed(args.mock_root,args.seed,args.scale)
    elif args.phase=='freeze': result=freeze(args.mock_root,args.scale)
    elif args.phase=='seed': result=seed_phase(args.mock_root,args.seed,args.variant,args.scale)
    elif args.phase=='control': result=control(args.mock_root,args.name,args.scale)
    else: result=collect(args.mock_root,args.scale,args.output)
    print('COMPLETE',args.phase,args.seed,args.variant,flush=True)
