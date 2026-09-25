"""Stage A iteration 7 campaign: low-redshift-tracer cross-correlation (GATES v7).

Phases (each idempotent, provenance-checked, own directory): basis -> dev-seed x5 -> freeze -> seed x N ->
control {numerical, injection, benchmark} -> collect. Reuses the completion/attempt machinery and the basis
tables of campaign4; its own fingerprint (iteration 7, seed ranges, tracer table, bound).
"""
from pathlib import Path
import argparse, dataclasses, json, os, resource, shutil, time
import numpy as np
from numba import set_num_threads
import run_mock_validation as v
import run_lowz_validation as lz
from campaign4 import digest, completion, begin, finish, read_xi, campaign_config, basis_phase as _basis4
from mock import generate_mock, save_mock, sightlines_for_variant
from validation_stats import absolute_statistics, slope_statistics
from random_streams import STREAM_NAMES
from xi_fit import BASIS
from lowz import TRACERS, lowz_bundles, tracer_maps
from paths import MOCKS, ACT_MASK

ITERATION=7
A_GRID=(0,.5,1,2)
DEV_SEEDS=tuple(range(5000,5005))
SPARSE_SEEDS=tuple(range(4000,4400))
FULL_SEEDS=tuple(range(4000,4040))
EXTRAS_SEED=SPARSE_SEEDS[0]
SMOKE_SPARSE=SPARSE_SEEDS[:2]
INJECTION_AMPLITUDES=(-2,-1,-.5,-.25,0,.25,.5,1,2)
CONTROL_NAMES=('numerical','injection','benchmark')


def fingerprint(scale):
    sources=[p for p in v.CODE.rglob('*.py') if '__pycache__' not in p.parts]
    config={k:(float('%.12g'%val) if isinstance(val,float) else val)
            for k,val in dataclasses.asdict(campaign_config(scale)).items() if not isinstance(val,Path)}
    return {'iteration':ITERATION,'sources':{str(p.relative_to(v.ROOT)):digest(p) for p in sorted(sources)},
            'gates':digest(v.ROOT/'GATES.md'),'scale':float(scale),
            'config':json.loads(json.dumps(config,default=v.serializable)),
            'A_grid':list(A_GRID),'streams':list(STREAM_NAMES),
            'seeds':{'dev':list(DEV_SEEDS),'sparse':[SPARSE_SEEDS[0],SPARSE_SEEDS[-1]],'full':[FULL_SEEDS[0],FULL_SEEDS[-1]],'extras':EXTRAS_SEED},
            'tracers':[dataclasses.asdict(t) for t in TRACERS],'combined_bound':lz.COMBINED_BOUND,
            'injection_amplitudes':list(INJECTION_AMPLITUDES),
            'numbers3':digest(v.ROOT/'report/numbers3.json'),
            'ACT_mask':digest(ACT_MASK) if ACT_MASK.exists() else None,
            'versions':{name:__import__('importlib.metadata',fromlist=['version']).version(name)
                        for name in ('numpy','scipy','numba','healpy','h5py','camb','fitsio')}}


def basis_phase(root,scale):
    """Same basis tables as iteration 5-6 (forest model, continuum projection, coarse binning), own provenance."""
    directory=root/'basis'; directory.mkdir(parents=True,exist_ok=True)
    prov=fingerprint(scale); prior=completion(directory,prov)
    if prior is not None: return prior
    from xi_fit import basis_tables, project_fine, coarse_bin
    from mock import grid_geometry
    begin(directory,prov); start=time.perf_counter(); cfg=campaign_config(scale); geo=grid_geometry(cfg)
    raw=basis_tables(cfg,'grid',grid_shape=geo['shape'])
    proj={k:project_fine(raw[k],geo['cpix'],cfg) for k in BASIS}; coarse={k:coarse_bin(proj[k],cfg) for k in BASIS}
    import h5py
    path=directory/'basis.h5'
    for k in BASIS: raw[k].save(path,f'raw/{k}'); proj[k].save(path,f'projected/{k}')
    with h5py.File(path,'a') as f:
        for k in BASIS: f[f'coarse/{k}']=coarse[k]
        f['cpix']=geo['cpix']
    return finish(directory,prov,{'grid':list(geo['shape']),'pixels':int(len(geo['cpix'])),'basis':list(BASIS),'wall_s':time.perf_counter()-start})


def campaign_basis(root,scale):
    if completion(root/'basis',fingerprint(scale)) is None: raise RuntimeError('basis phase must complete first')
    import h5py
    path=root/'basis/basis.h5'
    out={'raw':{k:read_xi(path,f'raw/{k}') for k in BASIS},'projected':{k:read_xi(path,f'projected/{k}') for k in BASIS}}
    with h5py.File(path) as f: out['coarse']={k:f[f'coarse/{k}'][()] for k in BASIS}; out['cpix']=f['cpix'][()]
    return out


def dev_seed(root,seed,scale):
    """Chain check on a development seed: mock with low-z tracers, own templates, truth and combined fits on the
    A grid. No numerical choice is made from it."""
    if seed not in DEV_SEEDS: raise ValueError(f'development seeds must be {DEV_SEEDS[0]}--{DEV_SEEDS[-1]}')
    directory=root/f'dev/{seed}'; directory.mkdir(parents=True,exist_ok=True)
    prov=fingerprint(scale); prior=completion(directory,prov)
    if prior is not None: return prior
    if not ACT_MASK.exists(): raise FileNotFoundError(ACT_MASK)
    basis=campaign_basis(root,scale); begin(directory,prov); start=time.perf_counter(); cfg=campaign_config(scale)
    m=generate_mock(cfg,seed,A_true=1,response=True,real_mask=True,completeness=True,magnification=True,cmb_noise=True,
                    disjoint_selection=True,variant_A_values=A_GRID,lowz=True)
    save_mock(m,directory/'mock.h5'); b=lowz_bundles(m,cfg); lz.save_tracer_maps(directory/'fits.h5',b['maps'])
    pairs=v.find_pairs(m.sightlines,cfg.r_perp_max/m.sightlines.chi.min()); values={'truth':[],'combined':[]}; params=[]
    for A in A_GRID:
        sl=sightlines_for_variant(m,A,True); ft=v.table_for(sl,cfg,basis); ft.save(directory/'fits.h5',f'xi/{A}')
        cat=v.cat_for(sl,ft.table,cfg,pairs)
        for name in values: values[name].append(v.fit_save(cat,b['templates'][name],cfg,sl,directory/'fits.h5',f'{name}/A{A}')['A'])
        params.append(ft.params)
    return finish(directory,prov,{'seed':seed,'A_values':values,'xi_fit':params,'templates':lz.tracer_summary(b['maps']),
                                  'wall_s':time.perf_counter()-start})


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
    result={'seeds':list(DEV_SEEDS),'development_slopes':{n:slope_statistics([d['A_values'][n] for d in dev],A_GRID) for n in ('truth','combined')},
            'results':dev,'frozen_choices':{'r_perp_min':cfg.r_perp_min,'fit_rperp_min':cfg.fit_rperp_min,'forest_model':cfg.forest_model,
                                            'tracers':[dataclasses.asdict(t) for t in TRACERS],'combined_bound':lz.COMBINED_BOUND,
                                            'note':'no numerical choice is made from the development seeds; they check the chain'},
            'development_markers':{str(s):digest(root/f'dev/{s}/complete.json') for s in DEV_SEEDS},
            'basis_marker':digest(root/'basis'/'complete.json')}
    v.dump(directory/'development.json',result); v.dump(directory/'frozen.json',prov)
    shutil.copyfile(v.ROOT/'GATES.md',directory/'GATES.md')
    return finish(directory,prov,result)


def frozen(root,scale):
    current=fingerprint(scale); result=completion(root/'freeze',current)
    if result is None: raise RuntimeError('freeze phase must complete first')
    if completion(root/'basis',current) is None: raise RuntimeError('basis phase must complete first')
    basis_marker=digest(root/'basis'/'complete.json')
    if result.get('basis_marker')!=basis_marker: raise RuntimeError('basis changed after the freeze')
    for seed in DEV_SEEDS:
        if completion(root/f'dev/{seed}',current) is None: raise RuntimeError(f'missing development seed {seed}')
        if digest(root/f'dev/{seed}/complete.json')!=result['development_markers'][str(seed)]: raise RuntimeError(f'changed development seed {seed}')
    return campaign_config(scale),current|{'freeze':digest(root/'freeze'/'complete.json'),'basis':basis_marker}


def fixed_maps(root,cfg):
    """The template maps of the first development seed: the other-realisation (fixed) templates of every seed."""
    return tracer_maps(v.load_mock(root/f'dev/{DEV_SEEDS[0]}/mock.h5'),cfg)


def seed_phase(root,seed,scale):
    if seed not in SPARSE_SEEDS: raise ValueError(f'sparse seeds must be {SPARSE_SEEDS[0]}--{SPARSE_SEEDS[-1]}')
    cfg,prov=frozen(root,scale)
    directory=root/f'sparse/{seed}'; directory.mkdir(parents=True,exist_ok=True)
    prov=prov|{'variant':'sparse','seed':seed}; prior=completion(directory,prov)
    if prior is not None: return prior
    basis=campaign_basis(root,scale); begin(directory,prov)
    role='full' if seed in FULL_SEEDS else 'core'
    result=lz.process_seed(seed,cfg,directory,role,basis,fixed_maps(root,cfg),extras=(seed==EXTRAS_SEED))
    return finish(directory,prov,result)


def control(root,name,scale):
    cfg,prov=frozen(root,scale)
    if name not in CONTROL_NAMES: raise ValueError(name)
    directory=root/f'controls/{name}'; directory.mkdir(parents=True,exist_ok=True)
    prov=prov|{'control':name}; prior=completion(directory,prov)
    if prior is not None: return prior
    basis=campaign_basis(root,scale); begin(directory,prov)
    start=time.perf_counter(); rows=[]; info={}
    def row(name,tol,metrics,ok): v.record(rows,0,name,tol,metrics,ok)
    if name=='numerical':
        import subprocess,sys
        t=time.perf_counter(); run=subprocess.run([sys.executable,'-m','pytest','-q','tests'],cwd=v.HERE,text=True,capture_output=True)
        (directory/'pytest.log').write_text(run.stdout+'\n'+run.stderr)
        row('unit and regression tests','all pass',{'exit_code':run.returncode,'wall_s':time.perf_counter()-t,'output':run.stdout.strip().splitlines()[-1:]},run.returncode==0)
    else:
        s0=EXTRAS_SEED; d=completion(root/f'sparse/{s0}',frozen(root,scale)[1]|{'variant':'sparse','seed':s0})
        if d is None: raise RuntimeError(f'sparse seed {s0} must complete before controls')
        base=v.load_mock(root/f'sparse/{s0}/seed{s0:03d}.h5'); sl=sightlines_for_variant(base,1,True)
        xi=v.table_for(sl,cfg,basis).table; b=lowz_bundles(base,cfg); path=directory/'products.h5'
        if name=='injection':
            amps=list(INJECTION_AMPLITUDES)
            exp=v.injection_test(sl,xi,base.alpha_lya,amps,cfg,templates=b['templates']['truth'],output=path,expectation=True)
            inj=v.injection_test(sl,xi,base.alpha_lya,amps,cfg,templates=b['templates']['truth'],output=path)
            info['injection']=inj; info['injection_expectation']=exp
            small=min(exp['paired_slopes_by_amplitude']); e_small=exp['paired_slopes_by_amplitude'][small]
            row('injection bookkeeping (noise-free expectation, production selection)',f'odd slope at |A| = {small:g} within 1 ± 0.05',
                {'expected_slope_small_amplitude':e_small,'expected_slopes_by_amplitude':exp['paired_slopes_by_amplitude'],'amplitude':small},abs(e_small-1)<=.05)
            row('injection convergence and measured slopes','report: expected slope vs |A|; measured odd slopes; curl',
                {'expected_by_amplitude':exp['paired_slopes_by_amplitude'],'measured_by_amplitude':inj['paired_slopes_by_amplitude'],
                 'measured_combined':inj['paired_slope'],'expected_combined':exp['paired_slope'],'curl_slope':inj['curl_slope']},np.isfinite(inj['paired_slope']))
        elif name=='benchmark':
            bench=v.benchmark(sl,.1,xi,cfg,directory/'benchmark.json'); peak=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2
            row('benchmark and memory','report: positive throughput; requested threads; peak <40 GB',{'benchmark':bench,'peak_gib':peak},
                bench['threads']==int(os.environ['NUMBA_NUM_THREADS']) and bench['pixel_pairs_per_s']>0 and peak*1024**3<40e9)
    return finish(directory,prov,info|{'rows':rows,'wall_s':time.perf_counter()-start})


def collect(root,scale,output):
    cfg,prov=frozen(root,scale); smoke=scale<1
    ids=list(SMOKE_SPARSE) if smoke else list(SPARSE_SEEDS); sparse=[]
    for seed in ids:
        d=completion(root/f'sparse/{seed}',prov|{'variant':'sparse','seed':seed})
        if d is None: raise RuntimeError(f'missing sparse seed {seed}')
        sparse.append(d)
    extras=[]
    for name in CONTROL_NAMES:
        d=completion(root/f'controls/{name}',prov|{'control':name})
        if d is None: raise RuntimeError(f'missing control {name}')
        extras.extend(d['rows'])
    rows=lz.build_rows(sparse,{'rows':extras})
    for r in rows: r.setdefault('required',not str(r['tolerance']).lower().startswith('report'))
    required=[r for r in rows if r['required']]
    result={'iteration':ITERATION,'smoke':smoke,'scale':scale,'acceptance':rows,'provenance':prov,
            'fixed_N':len(ids),'sparse':sparse,'required_passed':sum(r['pass'] for r in required),'required_total':len(required),
            'Stage_B_allowed':not smoke and all(r['pass'] for r in required)}
    v.dump(root/'collection.json',result); render(result,output)
    return result


def render(result,output):
    output=Path(output); output.mkdir(parents=True,exist_ok=True)
    v.dump(output/'lowz_validation.json',result)
    label='SMOKE — not acceptance' if result['smoke'] else 'frozen campaign'
    lines=[f'# Stage A iteration {ITERATION} (low-redshift tracers): {label}','',
           f"Scale {result['scale']}; {len(result['sparse'])} seeds; {result['required_passed']}/{result['required_total']} required gates pass. "
           +("**Stage B remains blocked.**" if not result['Stage_B_allowed'] else '**Stage A accepted (low-redshift path).**'),'',
           'Absolute statistics F^-1(q-mf); every seed carries A_true = 0 and 1 with the response on (null and recovery rows share realisations). Templates: Wiener-combined lognormal tracers of the same realisation per redshift slice, summed over slices. Smoke precision cannot certify scale-1 gates.','',
           '| Gate | Measurement | Frozen tolerance | Type | Result |','|---|---|---|---|---|']
    for r in result['acceptance']:
        s=r['measured']
        metric=f"{s['mean']:.6g} ± {s['sem']:.4g} SEM; bound95 {s['bound95']:.4g}; N={s['n']}" if 'mean' in s else json.dumps(s,default=v.serializable)[:300]
        kind='required' if r.get('required',True) else 'report only'
        lines.append(f"| {r['acceptance']} | {metric.replace('|','&#124;')} | {r['tolerance'].replace('|','&#124;')} | {kind} | {'PASS' if r['pass'] else 'FAIL'} |")
    lines+=['','## Per-seed combined amplitudes','','| Seed | A(0) combined | A(1) combined | A(1) truth | jackknife combination A(1) ± err | slices A(1) |','|---|---|---|---|---|---|']
    for d in result['sparse']:
        sn=sorted(k for k in d['fits']['A1_R1'] if k.startswith('slice'))
        lines.append(f"| {d['seed']} | {d['fits']['A0_R1']['combined']['A']:.3f} | {d['fits']['A1_R1']['combined']['A']:.3f} | {d['fits']['A1_R1']['truth']['A']:.3f} | {d['joint']['A1_R1']['A']:.3f} ± {d['joint']['A1_R1']['error']:.3f} | "+', '.join(f"{d['fits']['A1_R1'][n]['A']:.2f}" for n in sn)+' |')
    lines+=['','## Tracer bias fits (first seed)','','| Tracer | b true | b fit ± sigma | objects |','|---|---|---|---|']
    for lab,t in sorted(result['sparse'][0]['templates']['tracers'].items()):
        lines.append(f"| {lab} | {t['bias_true']:.2f} | {t['bias_fit']['b']:.3f} ± {t['bias_fit']['sigma_b']:.3f} | {t['n_objects']} |")
    (output/'lowz_validation.md').write_text('\n'.join(lines)+'\n')
    figs=output/'figures'; figs.mkdir(exist_ok=True)
    rec=[d for d in result['sparse'] if d['role'] in lz.GRID_ROLES]
    if rec:
        fig,ax=v.plt.subplots(figsize=(6,4))
        for name in ('combined','truth'):
            y=np.array([[d['fits'][f'A{A:g}_R1'][name]['A'] for A in A_GRID] for d in rec])
            ax.errorbar(A_GRID,y.mean(axis=0),yerr=y.std(axis=0,ddof=1)/np.sqrt(len(y)),marker='o',label=name)
        ax.plot([0,2],[0,2],'k--'); ax.set(xlabel='A true',ylabel='fitted A ± SEM',title=label); ax.legend()
        fig.tight_layout(); fig.savefig(figs/'lowz_iteration7_normalisation.pdf'); v.plt.close(fig)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--phase',required=True,choices=['basis','dev-seed','freeze','seed','control','collect'])
    ap.add_argument('--seed',type=int); ap.add_argument('--variant',choices=['sparse'],default='sparse')
    ap.add_argument('--name',choices=CONTROL_NAMES); ap.add_argument('--scale',type=float,default=1.)
    ap.add_argument('--mock-root',type=Path,default=MOCKS/f'iteration{ITERATION}'); ap.add_argument('--output',type=Path,default=v.ROOT/'report')
    args=ap.parse_args(); set_num_threads(int(os.environ.get('NUMBA_NUM_THREADS','4')))
    if not 0<args.scale<=1: ap.error('scale must be in (0,1]')
    args.mock_root.mkdir(parents=True,exist_ok=True)
    if args.phase=='basis': result=basis_phase(args.mock_root,args.scale)
    elif args.phase=='dev-seed': result=dev_seed(args.mock_root,args.seed,args.scale)
    elif args.phase=='freeze': result=freeze(args.mock_root,args.scale)
    elif args.phase=='seed': result=seed_phase(args.mock_root,args.seed,args.scale)
    elif args.phase=='control': result=control(args.mock_root,args.name,args.scale)
    else: result=collect(args.mock_root,args.scale,args.output)
    print('COMPLETE',args.phase,args.seed,flush=True)


if __name__=='__main__': main()
