"""Stage A iteration 7: mock acceptance of the low-redshift-tracer cross-correlation (GATES v7).

Per seed: a mock with the lognormal low-z tracers (lowz.py), the forest lensed by kappa_lya whose foreground is
the sum of the slice contributions; templates from the seed's own tracers (bias from their auto-spectra, Wiener
combination per slice, sum over slices), the truth template, and the templates of another realisation (fixed);
fits at A_true in {0, 1} (A grid on the 'full' seeds) with the response on (the response term does not correlate
with a low-z template; the response-off pair is a diagnostic on the full seeds). Rows: bias recovery per tracer,
per-slice and combined nulls and recoveries, the jackknife combination of the slice amplitudes against the
combined-template amplitude, normalisation slope, covariance, curl, fixed-template null, spectra.
"""
from __future__ import annotations
import json,time,resource
from pathlib import Path
import numpy as np
import run_mock_validation as v
from mock import generate_mock,save_mock,sightlines_for_variant,mock_spectrum_check
from lowz import lowz_bundles,tracer_summary,save_tracer_maps,joint_amplitudes,optimal_combination,TRACERS
from validation_stats import absolute_statistics,slope_statistics

GRID_ROLES=('full',)
COMBINED_BOUND=0.5     # GATES v7: 95 % residual bound on the combined null and recovery, in A (chosen with N)


def slice_names(bundle): return sorted(k for k in bundle['templates'] if k.startswith('slice'))


def process_seed(seed,cfg,root,role,basis,fixed_maps,extras=False):
    t=time.perf_counter(); path=root/f'seed{seed:03d}.h5'; fits_path=root/f'fits{seed:03d}.h5'
    avals=[0,.5,1,2] if role in GRID_ROLES else [0,1]
    m=generate_mock(cfg,seed,A_true=1,response=True,magnification=True,completeness=True,real_mask=True,cmb_noise=True,
                    variant_A_values=avals,disjoint_selection=True,lowz=True)
    save_mock(m,path); b=lowz_bundles(m,cfg,fixed=fixed_maps); save_tracer_maps(fits_path,b['maps'])
    names=slice_names(b)+['combined','truth']+(['fixed_combined'] if fixed_maps is not None else [])
    pairs=v.find_pairs(m.sightlines,cfg.r_perp_max/m.sightlines.chi.min())
    diag={'seed':seed,'role':role,'fits':{},'joint':{},'xi_fit':{},'templates':tracer_summary(b['maps']),
          'config':vars(cfg),'mock_attrs':m.attrs,'files':{'mock':str(path),'fits':str(fits_path)}}
    todo=[(A,True) for A in avals]+([(0,False),(1,False)] if role in GRID_ROLES else [])
    for A,resp in todo:
        key=f'A{A:g}_R{int(resp)}'; sl=sightlines_for_variant(m,A,resp)
        ft=v.table_for(sl,cfg,basis); ft.save(fits_path,f'xi/{key}'); diag['xi_fit'][key]=ft.params
        cat=v.cat_for(sl,ft.table,cfg,pairs); cat.save(fits_path,f'catalogues/{key}')
        fits,jk,reg=joint_amplitudes(cat,b,cfg,sl,names)
        diag['fits'][key]=fits
        sn=slice_names(b); diag['joint'][key]=optimal_combination([fits[n]['A'] for n in sn],jk[[names.index(n) for n in sn]])
        import h5py
        with h5py.File(fits_path,'a') as f:
            g=f.require_group(f'jackknife/{key}'); g['samples']=jk; g.attrs['names']=json.dumps(names); g['regions']=reg
        print(key,{n:round(fits[n]['A'],3) for n in names},'joint',round(diag['joint'][key]['A'],3),flush=True)
    diag['spectra']=mock_spectrum_check(m)
    diag['wall_s']=time.perf_counter()-t; diag['peak_memory_gb']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2
    v.dump(root/f'diagnostics{seed:03d}.json',diag)
    print('SEED COMPLETE',seed,role,'seconds',diag['wall_s'],'A_combined',diag['fits']['A1_R1']['combined']['A'],flush=True)
    return diag


def build_rows(diags,extras,bound=COMBINED_BOUND):
    """GATES v7 rows. Every seed carries A_true in {0, 1} with the response on; the null and recovery rows share
    realisations (declared). Slope rows use the A-grid seeds."""
    rows=[]
    def row(name,tol,stats,ok,detail='',required=None): v.record(rows,len(rows)+1,name,tol,stats,ok,detail,required)
    every=list(diags); rec=[d for d in diags if d['role'] in GRID_ROLES]
    sn=sorted(k for k in every[0]['fits']['A0_R1'] if k.startswith('slice'))
    # Bias recovery per tracer (own auto-spectrum against the generator's bias).
    labels=sorted(every[0]['templates']['tracers'])
    for lab in labels:
        vals=[d['templates']['tracers'][lab]['bias_used']/d['templates']['tracers'][lab]['bias_true'] for d in every]
        st=absolute_statistics(vals,1); st['mean_sigma_b_over_b']=float(np.mean([d['templates']['tracers'][lab]['bias_fit']['sigma_b']/d['templates']['tracers'][lab]['bias_true'] for d in every]))
        row(f'bias from the angular auto-spectrum: {lab}','|mean-1| <= 2 SEM',st,abs(st['residual'])<=2*st['sem'],
            'b^2 from the masked auto-spectrum of the kernel-weighted map in 40 <= L <= 300, shot noise subtracted, mask coupling and pixel window on the theory.')
    # Per-slice null and recovery (own templates, response on).
    for n in sn:
        st=absolute_statistics([d['fits']['A0_R1'][n]['A'] for d in every]); row(f'null (A_true = 0), own template, {n}','|mean| <= 2 SEM',st,abs(st['mean'])<=2*st['sem'])
        st=absolute_statistics([d['fits']['A1_R1'][n]['A'] for d in every],1); row(f'recovery (A_true = 1), own template, {n}','|mean-1| <= 2 SEM',st,abs(st['residual'])<=2*st['sem'])
    # Combined template: the acceptance rows.
    st=absolute_statistics([d['fits']['A0_R1']['combined']['A'] for d in every])
    row('combined null (A_true = 0), own templates',f'|mean| <= 2 SEM and bound95 <= {bound:g} A',st,abs(st['mean'])<=2*st['sem'] and st['bound95']<=bound,
        'Sum of the Wiener-filtered slice maps; lognormal tracers of the same realisation; response on.')
    st=absolute_statistics([d['fits']['A1_R1']['combined']['A'] for d in every],1)
    row('combined recovery (A_true = 1), own templates',f'|mean-1| <= 2 SEM and residual bound95 <= {bound:g} A',st,abs(st['residual'])<=2*st['sem'] and st['bound95']<=bound,
        'Unbiased detection against a lognormal redshift tracer: the gate of iteration 7.')
    st=absolute_statistics([d['fits']['A1_R1']['combined']['A']-d['fits']['A0_R1']['combined']['A'] for d in every],1)
    row('paired combined response A(1) - A(0), same realisation','|mean-1| <= 2 SEM',st,abs(st['residual'])<=2*st['sem'])
    # Jackknife combination of the slice amplitudes against the combined template.
    diff=[d['joint']['A1_R1']['A']-d['fits']['A1_R1']['combined']['A'] for d in every if np.isfinite(d['joint']['A1_R1']['A'])]
    st=absolute_statistics(diff); st['joint_mean_error']=float(np.mean([d['joint']['A1_R1']['error'] for d in every if np.isfinite(d['joint']['A1_R1']['A'])]))
    st['combined_mean_jk_error']=float(np.mean([d['fits']['A1_R1']['combined']['jk_error'] for d in every]))
    row('jackknife-covariance combination of the slice amplitudes minus the combined-template amplitude (A_true = 1)','report only',st,np.isfinite(st['mean']),
        'Optimal combination with the Hartlap-corrected joint jackknife covariance; the two should agree within their errors.')
    # Truth template, fixed template, curl, covariance, slopes, spectra.
    st=absolute_statistics([d['fits']['A1_R1']['truth']['A']-d['fits']['A0_R1']['truth']['A'] for d in every],1)
    row('paired truth-template response A(1) - A(0)','report only',st,np.isfinite(st['mean']))
    if 'fixed_combined' in every[0]['fits']['A0_R1']:
        st=absolute_statistics([d['fits']['A1_R1']['fixed_combined']['A'] for d in every])
        row('fixed (other-realisation) combined template, A_true = 1','|mean| <= 2 SEM',st,abs(st['mean'])<=2*st['sem'])
    st=absolute_statistics([d['fits']['A1_R1']['combined']['curl'] for d in every])
    row('curl null, combined template (A_true = 1)','|mean| <= 2 SEM',st,abs(st['mean'])<=2*st['sem'])
    scatter=float(np.std([d['fits']['A1_R1']['combined']['A'] for d in every],ddof=1))
    rms=float(np.sqrt(np.mean([d['fits']['A1_R1']['combined']['jk_error']**2 for d in every])))
    jk={'nside':sorted({int(d['fits']['A1_R1']['combined']['nside_jk']) for d in every}),'regions_mean':float(np.mean([d['fits']['A1_R1']['combined']['nregion'] for d in every]))}
    row('absolute covariance, combined template','0.7 <= scatter/RMS jackknife <= 1.3',{'scatter':scatter,'rms_jk':rms,'ratio':scatter/rms,'jackknife':jk},.7<=scatter/rms<=1.3)
    if rec:
        for name in ('combined','truth'):
            st=slope_statistics([[d['fits'][f'A{A:g}_R1'][name]['A'] for A in (0,.5,1,2)] for d in rec])
            row(f'{name} normalisation slope (A grid, response on)','|slope-1| <= 2 SEM',st,abs(st['residual'])<=2*st['sem'])
        st=absolute_statistics([d['fits']['A0_R1']['combined']['A']-d['fits']['A0_R0']['combined']['A'] for d in rec])
        row('response on minus off, combined template, A_true = 0 (A-grid seeds)','|mean| <= 2 SEM',st,abs(st['mean'])<=2*st['sem'],
            'The forest response to its own long modes must not correlate with a low-redshift template.')
    sv=np.array([d['spectra']['ratio'] for d in every]); st={'ratios':sv.mean(axis=0).tolist(),'sem':(sv.std(axis=0,ddof=1)/np.sqrt(len(sv))).tolist()}
    row('foreground and total spectra klkl, klkc, kckc','each total ratio within 10%',st,np.all(abs(np.array(st['ratios'])-1)<=.1))
    for r in extras.get('rows',[]):
        r.setdefault('required',not str(r.get('tolerance','')).lower().startswith('report')); rows.append(r)
    return rows
