"""Stage-A mock acceptance runner (IMPLEMENTATION.md section 2.8)."""
from __future__ import annotations

import argparse,json,resource,time,sys
from functools import lru_cache
from pathlib import Path
import numpy as np
from scipy.integrate import trapezoid
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent; CODE=HERE.parent; ROOT=CODE.parent
if str(CODE) not in sys.path: sys.path.insert(0,str(CODE))
from cosmo import chi as chi_of_z,z_of_chi
from forest_power import ForestPower
from config import Config,kernel_product_g1,SightlineSet
from mock import generate_mock,save_mock,sightlines_for_variant,mock_spectrum_check
from xi_model import xi_from_model,xi_from_data
from pairs import find_pairs,accumulate,benchmark,pair_midpoint_regions
from templates import (Template,cosine_band,flat_sky_band_templates,
                       flat_sky_wiener_transfer,matched_template_flat)
from amplitude import amplitude,independent_response_prediction
from inject import injection_test
from three_tracer import spectra

SCIENCE=((40,100),(100,200),(200,300))


def record(rows,item,name,tolerance,metrics,passed,detail=""):
    rows.append({"item":item,"acceptance":name,"tolerance":tolerance,"measured":metrics,
                 "pass":bool(passed),"detail":detail})


def mean_sem(values):
    x=np.asarray(values,float)
    return float(np.mean(x)),float(np.std(x,ddof=1)/np.sqrt(len(x))) if len(x)>1 else np.nan


def cat_for(sl,xi,cfg,pairs=None):
    if pairs is None: pairs=find_pairs(sl,cfg.r_perp_max/max(float(sl.chi.min()),1))
    start=time.perf_counter(); cat=accumulate(sl,pairs,xi,cfg)
    cat.attrs["wall_s"]=time.perf_counter()-start
    return cat


def midpoint_regions(cat,sl):
    reg=pair_midpoint_regions(cat,sl,8); nside=8
    if len(np.unique(reg))<30:
        reg=pair_midpoint_regions(cat,sl,16); nside=16
    return reg,nside,len(np.unique(reg))


def solve(F,y):
    return np.linalg.solve(F,y)


def common_science(result,nscience=3):
    """Collapse independently fitted science bands to one common amplitude."""
    diag=np.diag(result.F)
    active=list(range(nscience)); nuisance=list(range(nscience,len(diag)))
    M=np.zeros((len(diag),1+len(nuisance))); M[active,0]=1
    for j,i in enumerate(nuisance,1): M[i,j]=1
    Fr=M.T@result.F@M; yr=M.T@(result.q-result.mf); pars=solve(Fr,yr)
    jk=[]
    for k in range(len(result.regions)):
        Fk=result.F-result.partial_F[k]
        yk=(result.q-result.partial_q[k])-(result.mf-result.partial_mf[k])
        jk.append(solve(M.T@Fk@M,M.T@yk)[0])
    jk=np.asarray(jk); nr=len(jk)
    jke=float(np.sqrt((nr-1)/nr*np.sum((jk-jk.mean())**2))) if nr>1 else np.nan
    return {"A":float(pars[0]),"jk_error":jke,
            "sigma_F":float(np.sqrt(max(np.linalg.inv(Fr)[0,0],0))),
            "jk":jk,"active":active}


def fit_bundle(cat,bundle,g1,sl):
    reg,nside,nreg=midpoint_regions(cat,sl); result=amplitude(cat,bundle,g1,reg)
    summary=common_science(result); summary.update(result=result,nside_jk=nside,nregion=nreg)
    return summary


def make_bundles(mock,template_margin=None,fixed_maps=None):
    """Common masked/pixelized maps; CXX measured from each actual noisy map.

    Signal spectra include the actual radial window and magnification coefficient.
    The same mask convolution is applied to S and to the observed map power.
    The ratio is frozen in annuli of width 40, with no fitted normalization.
    """
    from scipy.ndimage import zoom
    from scipy.signal import resample
    from templates import map_spectrum
    from cross_spectrum import kernel,Z_CMB
    side=np.deg2rad(20*float(mock.attrs['scale'])); n=128; shape=(n,n); pix=side/n
    margin=float(mock.attrs['template_margin']) if template_margin is None else template_margin
    c1=float(chi_of_z(mock.sightlines.attrs['zmin'])); c2=float(chi_of_z(mock.sightlines.attrs['zmax']))
    def selection(cat):
        keep=(cat['chi']>=c1-margin)&(cat['chi']<=c2+margin)
        return {k:np.asarray(v)[keep] for k,v in cat.items()}
    q=selection(mock.quasars); r=selection(mock.randoms)
    matched,mm,meta=matched_template_flat(q,r,lambda z:np.full_like(z,3.5),shape,pix)
    def resize(a): return resample(resample(a,n,axis=0),n,axis=1).real
    act=zoom(mock.maps['act_mask'],n/mock.maps['act_mask'].shape[0],order=1)
    mask=np.clip(act,0,1)*mm
    # CMB saved unmasked as well: common mask is applied exactly once.
    cmbmap=resize(mock.maps['kappa_CMB_unmasked'])*mask
    matched=matched*mask
    depmap=cmbmap-matched
    edges=np.arange(0,np.sqrt(2)*np.pi/pix+80,40.); L=.5*(edges[1:]+edges[:-1])
    S=spectra(z1=mock.sightlines.attrs['zmin'],z2=mock.sightlines.attrs['zmax'],zs=2.4,
              template_z1=float(z_of_chi(c1-margin)),template_z2=float(z_of_chi(c2+margin)),
              L_values=L,n_q_slab=len(q['ra'])/np.rad2deg(side)**2,mag=.5*bool(mock.attrs['magnification']))
    cc=np.linspace(c1-margin,c2+margin,400)
    magcoef=(.5/3.5)*trapezoid(kernel(cc,float(chi_of_z(Z_CMB))),cc) if mock.attrs['magnification'] else 0.
    Sm=S['skl']+magcoef*S['klkl']
    # NGP count pixel window and common-mask mode coupling on the full FFT grid.
    ellx=2*np.pi*np.fft.fftfreq(n,pix); ell=np.hypot(ellx[:,None],ellx[None,:])
    window=np.sinc(ellx[:,None]*pix/(2*np.pi))*np.sinc(ellx[None,:]*pix/(2*np.pi))
    maskpower=abs(np.fft.fft2(mask))**2
    def masked_cross(s,win=1.):
        grid=np.interp(ell,L,s,left=0,right=0)*win
        convolved=np.fft.ifft2(np.fft.fft2(maskpower)*np.fft.fft2(grid)).real/n**4/max(np.mean(mask**2),1e-30)
        count=np.histogram(ell,edges)[0]
        return np.divide(np.histogram(ell,edges,weights=convolved)[0],count,out=np.zeros_like(L),where=count>0)
    sc=masked_cross(S['klkc']); sm=masked_cross(Sm,window); sd=sc-sm
    used={'cmb':cmbmap,'matched':matched,'deprojected':depmap}
    if fixed_maps is not None:
        used=fixed_maps['used_maps']; mask=fixed_maps['common_mask']
        matched=used['matched']; mm=fixed_maps['matched_mask']; meta=fixed_maps['matched_meta']
    out={}; operators={}; spectra_used={}
    for name,signal in [('cmb',sc),('matched',sm),('deprojected',sd)]:
        cxx=map_spectrum(used[name],pix,edges,mask)
        transfer=flat_sky_wiener_transfer(shape,pix,L,signal,cxx)
        if fixed_maps is not None:
            transfer=fixed_maps['operators'][name]
            signal=fixed_maps['spectra_used'][name]['S']; cxx=fixed_maps['spectra_used'][name]['CXX']
        out[name],_=flat_sky_band_templates(used[name],mock.sightlines.ra,mock.sightlines.dec,pix,transfer=transfer,source=name)
        operators[name]=transfer; spectra_used[name]={'S':signal,'CXX':cxx}
    truthmap=mock.maps['kappa_lya'] if fixed_maps is None else fixed_maps['truth_map']
    truthpix=side/truthmap.shape[0]
    out['truth'],_=flat_sky_band_templates(truthmap,mock.sightlines.ra,mock.sightlines.dec,truthpix,source='unfiltered truth')
    out.update(matched_meta=meta,matched_map=matched,matched_mask=mm,common_mask=mask,
               used_maps=used,operators=operators,spectra_used=spectra_used,L=L,spectra=S,
               truth_map=truthmap,pixel_size=pix,magcoef=magcoef)
    return out


def shape_consistency(cat,bundle,g1,sl):
    vals=[]; errs=[]; reg,_,_=midpoint_regions(cat,sl)
    for ib in range(6):
        summary=common_science(amplitude(cat,bundle,g1,reg,bins=[ib]))
        vals.append(summary["A"]); errs.append(summary["jk_error"])
    return np.asarray(vals),np.asarray(errs)


def gaussian_flat_map(shape,side,L,Cl,rng):
    nx,ny=shape; pixarea=(side/nx)*(side/ny)
    lx=2*np.pi*np.fft.fftfreq(nx,side/nx); ly=2*np.pi*np.fft.rfftfreq(ny,side/ny)
    ell=np.hypot(lx[:,None],ly[None,:]); c=np.interp(ell,L,Cl,left=0,right=Cl[-1])
    z=np.fft.rfft2(rng.normal(size=shape)); f=z*np.sqrt(np.maximum(c,0)/pixarea); f[0,0]=0
    return np.fft.irfft2(f,s=shape).astype(np.float32)



def serializable(x):
    if isinstance(x,np.ndarray): return x.tolist()
    if isinstance(x,np.generic): return x.item()
    if isinstance(x,Path): return str(x)
    raise TypeError(type(x).__name__)

def dump(path,data):
    path=Path(path); tmp=path.with_name(path.name+'.tmp')
    tmp.write_text(json.dumps(data,indent=2,default=serializable)+'\n'); tmp.replace(path)

def save_bundles(path,bundles,group='templates'):
    import h5py
    with h5py.File(path,'a') as f:
        if group in f: del f[group]
        g=f.create_group(group)
        for name in ('truth','cmb','matched','deprojected'):
            tg=g.create_group(name)
            for i,t in enumerate(bundles[name]):
                tt=tg.create_group(str(i)); tt['alpha']=t.alpha; tt['phi_lm']=t.phi_lm
                for k in ('name','kind','Lmin','Lmax','filter','source'): tt.attrs[k]=getattr(t,k)
        for name in ('common_mask','matched_map','L'): g[name]=bundles[name]
        for category in ('operators','used_maps'):
            gg=g.create_group(category)
            for name,value in bundles[category].items(): gg[name]=value
        gg=g.create_group('spectra_used')
        for name,value in bundles['spectra_used'].items():
            ng=gg.create_group(name)
            for k,v in value.items(): ng[k]=v
        g.attrs['matched_meta']=json.dumps(bundles['matched_meta'],default=serializable)

def load_mock(path):
    import h5py
    from mock import MockResult,load_sightlines
    with h5py.File(path,'r') as f:
        ds={name:{k:v[()] for k,v in f[name].items()} for name in ('truth','quasars','randoms','maps')}
        attrs=dict(f.attrs)
    return MockResult(load_sightlines(path),ds['truth']['alpha_lya'],ds['truth'],ds['quasars'],ds['randoms'],ds['maps'],attrs)

def fit_save(cat,bundle,cfg,sl,path,group):
    result=fit_bundle(cat,bundle,cfg.g1,sl); result['result'].save(path,group)
    rr=result['result']; raw=common_science(__import__('dataclasses').replace(rr,mf=np.zeros_like(rr.mf)))['A']
    return {k:v for k,v in result.items() if k!='result'}|{'raw':raw,'bands':rr.A}

def development(root):
    """Numerical decisions use ONLY seeds 100--104, never acceptance seeds."""
    root=Path(root); root.mkdir(parents=True,exist_ok=True); t=time.perf_counter()
    cfg=Config(scale=1); results=[]
    prior_progress=root/"development_progress.json"
    known={d["seed"]:d for d in json.loads(prior_progress.read_text())} if prior_progress.exists() else {}
    for seed in range(100,105):
        if seed in known:
            results.append(known[seed]); continue
        path=root/f'dev{seed}.h5'
        if path.exists(): m=load_mock(path)
        else:
            m=generate_mock(cfg,seed,A_true=1,response=True,magnification=True,completeness=True,real_mask=True,cmb_noise=True,variant_A_values=[0,1,5,10])
            save_mock(m,path)
        sl=m.sightlines; pix=np.deg2rad(20)/m.maps['kappa_lya'].shape[0]
        tb,_=flat_sky_band_templates(m.maps['kappa_lya'],sl.ra,sl.dec,pix)
        pairs=find_pairs(sl,cfg.r_perp_max/sl.chi.min()); vals={}
        # Development criterion: grid-converged width closest to unit physical slope.
        for width in (.76,.9,1.1):
            aa=[]
            for A in (0,1,5,10):
                sample=sightlines_for_variant(m,A,False); c=cfg.copy(xi_smoothing=width)
                xi=xi_from_data(sample,c); aa.append(fit_bundle(cat_for(sample,xi,c,pairs),tb,c.g1,sample)['A'])
            vals[str(width)]=aa
        results.append({'seed':seed,'widths':vals}); dump(root/'development_progress.json',results)
        print('DEVELOPMENT',seed,vals,flush=True)
        del m
    from validation_stats import slope_statistics
    stats={str(w):slope_statistics([r['widths'][str(w)] for r in results]) for w in (.76,.9,1.1)}
    width=min(stats,key=lambda w:abs(stats[w]['mean']-1))
    result={'seeds':list(range(100,105)),'xi_smoothing':float(width),'xi_step':.25,'analytic_nk':6400,
            'results':results,'slopes':stats,'wall_s':time.perf_counter()-t,'peak_memory_gb':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2,
            'choice_rule':'among .76,.9,1.1 choose smallest absolute deviation of development mean physical slope from 1'}
    dump(root/'development.json',result); return result


def process_seed(seed,cfg,root,role,raw_xi,fixed_bundle):
    from response import prediction_catalogue
    t=time.perf_counter(); path=root/f'seed{seed:03d}.h5'; fits_path=root/f'fits{seed:03d}.h5'
    recovery=role in ('recovery','extension'); null=role in ('null','extension')
    avals=[0,1,5,10] if recovery else [0,1]
    m=generate_mock(cfg,seed,A_true=1,response=True,magnification=True,completeness=True,
                    real_mask=True,cmb_noise=True,variant_A_values=avals)
    save_mock(m,path); b=make_bundles(m); save_bundles(fits_path,b)
    pairs=find_pairs(m.sightlines,cfg.r_perp_max/m.sightlines.chi.min())
    diag={'seed':seed,'role':role,'fits':{},'slopes':{},'predictions':{},'baseline_fixed':[], 'shape':[],
          'config':vars(cfg),'mock_attrs':m.attrs,'files':{'mock':str(path),'fits':str(fits_path)}}
    xfixed=xi_from_data(sightlines_for_variant(m,0,False),cfg); xfixed.save(fits_path,'xi_fixed')
    todo=[(A,False) for A in avals]+[(0,True),(1,True)]
    for A,resp in todo:
        key=f'A{A}_R{int(resp)}'; sl=sightlines_for_variant(m,A,resp)
        xi=xi_from_data(sl,cfg); xi.save(fits_path,f'xi/{key}')
        cat=cat_for(sl,xi,cfg,pairs); cat.save(fits_path,f'catalogues/{key}')
        if A==0 and resp:
            pc=prediction_catalogue(sl,cat,m.truth['delta_L_pixel'],raw_xi,xi,cfg)
            pc.save(fits_path,'prediction_catalogue')
            for name in ('cmb','matched','deprojected'):
                diag['predictions'][name]=fit_save(pc,b[name],cfg,sl,fits_path,f'predictions/{name}')
        diag['fits'][key]={}
        for name in ('truth','cmb','matched','deprojected'):
            diag['fits'][key][name]=fit_save(cat,b[name],cfg,sl,fits_path,f'fits/{key}/{name}')
        if not resp and recovery:
            fixedcat=cat_for(sl,xfixed,cfg,pairs)
            diag['baseline_fixed'].append(fit_save(fixedcat,b['truth'],cfg,sl,fits_path,f'baseline_fixed/{A}')['A'])
        if A==1 and resp:
            for ib in range(6):
                reg,_,_=midpoint_regions(cat,sl); rr=amplitude(cat,b['deprojected'],cfg.g1,reg,bins=[ib])
                rr.save(fits_path,f'shape/{ib}'); diag['shape'].append(common_science(rr)['A'])
        if null and A==0 and not resp:
            fixed=make_bundles(m,fixed_maps=fixed_bundle); save_bundles(fits_path,fixed,'fixed_templates')
            diag['fixed_null']={name:fit_save(cat,fixed[name],cfg,sl,fits_path,f'fixed_null/{name}') for name in ('truth','cmb','matched')}
        if seed==0 and A==1 and resp:
            # Margin changes use the SAME density, sightlines, CMB map and forest.
            diag['margin']={}
            for margin in (0.,150.,300.):
                bm=make_bundles(m,template_margin=margin); save_bundles(fits_path,bm,f'margin_templates/{margin}')
                diag['margin'][str(margin)]=fit_save(cat,bm['deprojected'],cfg,sl,fits_path,f'margin/{margin}')['A']
            # Independent Gaussian maps receive the identical common mask and transfer.
            rng=np.random.default_rng(9912); diag['random']=[]
            L=b['L']; S=b['spectra']; shape=b['common_mask'].shape; side=np.deg2rad(20)
            for i in range(100):
                km=gaussian_flat_map(shape,side,L,S['kckc']+m.attrs['cmb_noise_level'],rng)*b['common_mask']
                tb,_=flat_sky_band_templates(km,sl.ra,sl.dec,b['pixel_size'],transfer=b['operators']['cmb'],source='independent masked Gaussian')
                diag['random'].append(fit_save(cat,tb,cfg,sl,fits_path,f'random/{i}'))
                import h5py
                with h5py.File(fits_path,'a') as f: f[f'random_maps/{i}']=km
    diag['spectra']=mock_spectrum_check(m)
    diag['wall_s']=time.perf_counter()-t; diag['peak_memory_gb']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2
    dump(root/f'diagnostics{seed:03d}.json',diag)
    print('SEED COMPLETE',seed,role,'seconds',diag['wall_s'],'Acombined',diag['fits']['A1_R1']['deprojected']['A'],flush=True)
    return diag


def build_rows(diags,extras):
    from validation_stats import absolute_statistics,slope_statistics,hotelling_shape
    rec=[d for d in diags if d['role'] in ('recovery','extension')]
    null=[d for d in diags if d['role'] in ('null','extension')]
    rows=[]
    def row(name,tol,stats,ok,detail=''): record(rows,len(rows)+1,name,tol,stats,ok,detail)
    for name in ('truth','cmb','matched'):
        vals=[[d['fits'][f'A{A}_R0'][name]['A'] for A in (0,1,5,10)] for d in rec]
        st=slope_statistics(vals)
        row(f'{name} normalization slope','|slope-1| <= 0.03',st,abs(st['residual'])<=.03,
            'Slope of absolute amplitudes at A_true={0,1,5,10}; per-seed intercept fitted freely.')
    fixed=slope_statistics([d['baseline_fixed'] for d in rec])
    refit=slope_statistics([[d['fits'][f'A{A}_R0']['truth']['A'] for A in (0,1,5,10)] for d in rec])
    diff=absolute_statistics(np.array(refit['slopes'])-np.array(fixed['slopes']))
    diff['absolute_A1_refit_minus_fixed']=absolute_statistics([d['fits']['A1_R0']['truth']['A']-d['baseline_fixed'][1] for d in rec])
    row('fixed versus refitted baseline','|slope difference| <= 0.03',{'fixed_absolute_slope':fixed,'refitted_absolute_slope':refit,'diagnostic_difference':diff},abs(diff['mean'])<=.03)
    for name in ('truth','cmb','matched'):
        for conditional in (False,True):
            vals=[(d['fixed_null'] if conditional else d['fits']['A0_R0'])[name]['A'] for d in null]
            raw=[(d['fixed_null'] if conditional else d['fits']['A0_R0'])[name]['raw'] for d in null]
            st=absolute_statistics(vals); st['raw']=absolute_statistics(raw)
            row(f'{"fixed" if conditional else "varying"} template mean field: {name}','|mean| <= 2 SEM',st,abs(st['mean'])<=2*st['sem'])
    st=absolute_statistics([d['fits']['A1_R0']['cmb']['A'] for d in rec],1)
    row('absolute stochastic recovery','|mean-1| <= 2 SEM (precision established separately by slope)',st,abs(st['residual'])<=2*st['sem'])
    for name in ('cmb','matched','deprojected'):
        obs=np.array([d['fits']['A0_R1'][name]['A'] for d in null]); pred=np.array([d['predictions'][name]['A'] for d in null])
        st=absolute_statistics(obs); st['prediction']=absolute_statistics(pred); st['prediction_difference']=absolute_statistics(obs-pred)
        ok=abs(st['prediction_difference']['mean'])<=2*st['prediction_difference']['sem']
        tol='|observed-predicted| <= 2 SEM'
        if name=='deprojected':
            tol+='; absolute null bound95 <= 0.3 A'; ok=ok and st['bound95']<=.3
        row(f'absolute response-only: {name}',tol,st,ok,
            'Prediction uses stored delta_L at every pixel and P_a (D_a C_ab D_b - C_ab) P_b^T. Discrete-grid, phase-averaged trilinear covariance truncated at saved table support.')
    st=absolute_statistics([d['fits']['A1_R1']['deprojected']['A'] for d in rec],1)
    row('absolute combined deprojected recovery','|mean-1| <= 2 SEM and residual bound95 <= 0.3 A',st,abs(st['residual'])<=2*st['sem'] and st['bound95']<=.3)
    st=hotelling_shape([d['shape'] for d in rec]); row('six-bin Hotelling shape','p > 0.01',st,st['p_value']>.01)
    scatter=np.std([d['fits']['A1_R1']['deprojected']['A'] for d in rec],ddof=1)
    rms=np.sqrt(np.mean([d['fits']['A1_R1']['deprojected']['jk_error']**2 for d in rec]))
    row('absolute covariance','0.7 <= scatter/RMS jackknife <= 1.3',{'scatter':scatter,'rms_jk':rms,'ratio':scatter/rms},.7<=scatter/rms<=1.3)
    sv=np.array([d['spectra']['ratio'] for d in rec]); st={'ratios':sv.mean(axis=0).tolist(),'sem':(sv.std(axis=0,ddof=1)/np.sqrt(len(sv))).tolist()}
    row('outside-box and total spectra','each total ratio within 10%',st,np.all(abs(np.array(st['ratios'])-1)<=.1))
    for r in extras.get('rows',[]): rows.append(r)
    return rows


def rebuild(root,out_root):
    root=Path(root); out_root=Path(out_root); out_root.mkdir(parents=True,exist_ok=True)
    diags=[json.loads(p.read_text()) for p in sorted(root.glob('diagnostics*.json'))]
    extras=json.loads((root/'extras.json').read_text()) if (root/'extras.json').exists() else {}
    rows=build_rows(diags,extras)
    result={'iteration':3,'acceptance':rows,'provenance':json.loads((root/'provenance.json').read_text()),
            'seeds':diags,'extras':extras,'timing':{'ensemble_s':sum(d['wall_s'] for d in diags),'extras_s':extras.get('wall_s',0)},
            'peak_memory_gb':max([d['peak_memory_gb'] for d in diags]+[extras.get('peak_memory_gb',0)])}
    if (root/'development_interruption.json').exists(): result['development_interruption']=json.loads((root/'development_interruption.json').read_text())
    if (root/'campaign_timing.json').exists(): result['timing'].update(json.loads((root/'campaign_timing.json').read_text()))
    dump(out_root/'mock_validation.json',result)
    lines=['# Stage A iteration 3 validation','', 'Every recovery/null estimate below is absolute: F^-1(q-mf). Paired differences appear only as labelled diagnostics.',
           '', '| Gate | Statistic | Tolerance | Result |','|---|---|---|---|']
    for r in rows:
        st=r['measured']; summary=(f"{st['mean']:.6g} ± {st['sem']:.3g} SEM; residual bound95={st['bound95']:.4g} A" if 'mean' in st else json.dumps(st,default=serializable))
        if 'prediction' in st:
            pre=st['prediction']; diff=st['prediction_difference']
            summary+=f"; prediction {pre['mean']:.5g} ± {pre['sem']:.3g} SEM; obs−prediction {diff['mean']:.5g} ± {diff['sem']:.3g} SEM (bound95 {diff['bound95']:.4g})"
        elif 'fixed_absolute_slope' in st:
            f=st['fixed_absolute_slope']; d=st['refitted_absolute_slope']
            summary=f"Fixed {f['mean']:.5g} ± {f['sem']:.3g} SEM (bound95 {f['bound95']:.4g}); refitted {d['mean']:.5g} ± {d['sem']:.3g} SEM (bound95 {d['bound95']:.4g})"
        elif 'p_value' in st: summary=f"Hotelling T²={st['T2']:.4g}, p={st['p_value']:.4g}"
        elif 'benchmark' in st: summary=f"{st['benchmark']['pixel_pairs_per_s']:.5g} pixel pairs/s; peak {st['peak_gib']:.3f} GiB"
        elif 'passed' in st: summary=f"{st['passed']} passed; {st['failed']} failed; {st['wall_s']:.2f} s"
        lines.append(f"| {r['acceptance']} | {summary} | {r['tolerance']} | {'PASS' if r['pass'] else 'FAIL'} |")
    lines+=['','## Protocol','',json.dumps(result['provenance'],indent=2),'','## Wall time and memory','',json.dumps(result['timing']),f"Peak RSS: {result['peak_memory_gb']:.3f} GiB.",'','## Open failures','']
    lines += [f"- {r['acceptance']}: {r['detail'] or 'The declared tolerance was not met; no seed or tolerance was changed in response.'}" for r in rows if not r['pass']]
    (out_root/'mock_validation.md').write_text('\n'.join(lines)+'\n')
    figs=out_root/'figures'; figs.mkdir(exist_ok=True)
    rec=[d for d in diags if d['role'] in ('recovery','extension')]
    fig,axes=plt.subplots(1,3,figsize=(11,3.6))
    for ax,name in zip(axes,('truth','cmb','matched')):
        ys=np.array([[d['fits'][f'A{A}_R0'][name]['A'] for A in (0,1,5,10)] for d in rec])
        ax.errorbar([0,1,5,10],ys.mean(axis=0),yerr=ys.std(axis=0,ddof=1)/np.sqrt(len(ys)),marker='o',capsize=3)
        ax.plot([0,10],[0,10],'k--'); ax.set(xlabel='A true',ylabel='Absolute fitted A ± SEM',title=name)
    fig.tight_layout(); fig.savefig(figs/'mock_physical_normalisation.pdf'); plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4)); xs=[d['seed'] for d in rec]; ys=[d['fits']['A1_R1']['deprojected']['A'] for d in rec]; errors=[d['fits']['A1_R1']['deprojected']['jk_error'] for d in rec]
    ax.errorbar(xs,ys,yerr=errors,fmt='o'); ax.axhline(1,color='k',ls='--'); ax.set(xlabel='Seed',ylabel='Absolute deprojected A (jackknife)'); fig.tight_layout(); fig.savefig(figs/'mock_combined_ensemble.pdf'); plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4)); sv=np.array([d['shape'] for d in rec]); ax.errorbar(np.arange(6),sv.mean(axis=0),yerr=sv.std(axis=0,ddof=1)/np.sqrt(len(sv)),fmt='o'); ax.set(xlabel='Shape bin',ylabel='Absolute A ± SEM'); fig.tight_layout(); fig.savefig(figs/'mock_shape.pdf'); plt.close(fig)
    fig,ax=plt.subplots(figsize=(6,4)); sv=np.array([d['spectra']['ratio'] for d in rec]); ax.errorbar(np.arange(3),sv.mean(axis=0),yerr=sv.std(axis=0,ddof=1)/np.sqrt(len(sv)),fmt='o'); ax.axhline(1,color='k',ls='--'); ax.set_xticks([0,1,2],['kl-kl','kl-kCMB','kCMB-kCMB']); ax.set_ylabel('Measured / theory ± SEM'); fig.tight_layout(); fig.savefig(figs/'mock_spectra.pdf'); plt.close(fig)
    return result


def extra_checks(cfg,root):
    """Fresh seed-0 dense physical controls and explicit diagnostic changes."""
    from scipy.stats import linregress
    from validation_stats import absolute_statistics
    t=time.perf_counter(); rows=[]; info={}
    model=xi_from_model(ForestPower(model='kaiser'),cfg,nk=cfg.analytic_nk); model.save(root/'analytic.h5')
    doubled=xi_from_model(ForestPower(model='kaiser'),cfg.copy(xi_max=30),nk=2*cfg.analytic_nk)
    selection=(model.r_perp>=10)&(model.r_perp<=30); radial=model.r_par<=30
    ref=doubled.xi_rp[doubled.r_perp>=10]
    change=np.linalg.norm(model.xi_rp[selection][:,radial]-ref)/np.linalg.norm(ref)
    doubled.save(root/'analytic_doubled.h5')
    record(rows,0,'analytic derivative convergence','doubling nk changes derivative norm by <1%',{'nk':cfg.analytic_nk,'relative_norm_change':change},change<.01)
    quad=json.loads((root/'grid_covariance_convergence.json').read_text())
    record(rows,0,'discrete-grid covariance quadrature','doubling angular samples changes derivative norm by <1%',quad,quad['derivative_norm_change']<.01)
    preflight=json.loads((root/'preflight.json').read_text())
    info['preflight']=preflight
    record(rows,0,'unit and regression tests','all tests pass',{k:preflight[k] for k in ('passed','failed','wall_s')},preflight['failed']==0)

    slopes={}; densevals={}
    for gon in (True,False):
        m=generate_mock(cfg,0,A_true=1,g_on=gon,response=False,n_los=100,pixel_noise_power=0,variant_A_values=[0,1,5,10])
        path=root/f'dense_g{int(gon)}.h5'; save_mock(m,path)
        pix=np.deg2rad(20)/m.maps['kappa_lya'].shape[0]; basis,_=flat_sky_band_templates(m.maps['kappa_lya'],m.sightlines.ra,m.sightlines.dec,pix)
        vals=[]; analytic=[]; omitted=[]
        pairs=find_pairs(m.sightlines,cfg.r_perp_max/m.sightlines.chi.min())
        for A in (0,1,5,10):
            sl=sightlines_for_variant(m,A,False); xi=xi_from_data(sl,cfg); xi.save(path,f'xi_A{A}')
            cat=cat_for(sl,xi,cfg,pairs); cat.save(path,f'catalogue_A{A}')
            c=cfg.copy(g1=cfg.g1 if gon else 0.)
            vals.append(fit_save(cat,basis,c,sl,path,f'fits_A{A}')['A'])
            if gon:
                analytic.append(fit_save(cat_for(sl,model,cfg,pairs),basis,cfg,sl,path,f'analytic_A{A}')['A'])
                omitted.append(fit_save(cat,basis,cfg.copy(g1=0),sl,path,f'omitted_A{A}')['A'])
                if A==1:
                    finer=xi_from_data(sl,cfg.copy(xi_step=.125)); finer.save(path,'xi_fine')
                    af=fit_save(cat_for(sl,finer,cfg,pairs),basis,cfg,sl,path,'fit_fine')['A']
                    change=abs(af-vals[-1])/max(abs(af),1e-30)
                    record(rows,1,'xi interpolation convergence','relative absolute A change < 0.005',{'A_025':vals[-1],'A_0125':af,'fractional_change':change},change<.005)
                    # Mixed response: estimator uses g=0, truth uses g=cfg.g1.
                    from amplitude import pair_scalars
                    d,s,_=pair_scalars(cat,basis); dt=d[:3].sum(axis=0)+d[-1]; st=s[:3].sum(axis=0)+s[-1]
                    acc=cat.accum.copy(); acc[:,0]=(cat.accum[:,3]+cfg.g1*cat.accum[:,4])*dt[:,None]+.5*cfg.g1*cat.accum[:,6]*st[:,None]
                    acc[:,1:3]=0; acc[:,8:]=0
                    from pairs import PairCatalogue
                    predcat=PairCatalogue(cat.a,cat.b,cat.thx,cat.thy,cat.theta,acc,cat.npair)
                    predicted=fit_save(predcat,basis,cfg.copy(g1=0),sl,path,'omitted_prediction')['A']
        lr=linregress([0,1,5,10],vals); slopes[str(gon)]=lr.slope
        densevals[str(gon)]={'A_true':[0,1,5,10],'A':vals,'slope':lr.slope,'regression_stderr':lr.stderr,
                             'ensemble_sem':None,'ensemble_bound95':None,'error_note':'single dense realization; regression stderr is descriptive, not an ensemble SEM'}
        record(rows,2,f'dense physical g={gon}','slope within 0.05 of 1',densevals[str(gon)],abs(lr.slope-1)<=.05,'One full scale-1 noiseless dense realization. Error is regression error, not an ensemble SEM.')
        if gon:
            sm=linregress([0,1,5,10],analytic).slope; sn=linregress([0,1,5,10],omitted).slope
            record(rows,3,'data versus analytic baseline','slope ratio within 3%',{'data_slope':lr.slope,'analytic_slope':sm,'ratio':lr.slope/sm,'analytic_A':analytic},abs(lr.slope/sm-1)<=.03,
                   'The empirical 0.903 factor is removed. A continuum- and grid-mismatched analytic baseline can fail this gate.')
            record(rows,4,'first moments','omitted/correct slope ratio within 5% of catalogue prediction',{'observed_ratio':sn/lr.slope,'predicted_ratio':predicted},abs(sn/lr.slope-predicted)<=.05)
        del m,cat
    # A fiducial sample, not the dense one, supplies injection and flag diagnostics.
    base=load_mock(root/'seed000.h5'); sl=sightlines_for_variant(base,1,True); xi=xi_from_data(sl,cfg); b=make_bundles(base)
    inj=injection_test(sl,xi,base.alpha_lya,[-2,-1,0,1,2],cfg,templates=b['truth'],output=root/'injection.h5')
    curlerr=float(np.sqrt(np.mean(inj['errors'][:,3:6]**2)))
    record(rows,5,'injection bookkeeping','odd slope within 0.05 of 1; curl within jackknife error',{'slope':inj['paired_slope'],'curl_slope':inj['curl_slope'],'curl_error':curlerr},abs(inj['paired_slope']-1)<=.05 and abs(inj['curl_slope'])<=curlerr)
    info['injection']=inj; info['flag_shifts']={}
    fc=fit_bundle(cat_for(sl,xi,cfg),b['deprojected'],cfg.g1,sl)['A']
    for flag in ('magnification','completeness','real_mask'):
        flags=dict(magnification=True,completeness=True,real_mask=True); flags[flag]=False
        m=generate_mock(cfg,0,A_true=1,response=True,cmb_noise=True,**flags); path=root/f'flag_{flag}_off.h5'; save_mock(m,path)
        bb=make_bundles(m); save_bundles(path,bb); xx=xi_from_data(m.sightlines,cfg); xx.save(path)
        fit=fit_save(cat_for(m.sightlines,xx,cfg),bb['deprojected'],cfg,m.sightlines,path,'fit')
        info['flag_shifts'][flag]={'absolute_A':fit['A'],'baseline_A':fc,'diagnostic_difference':fit['A']-fc}
    record(rows,6,'flag and common-realization margin diagnostics','all requested diagnostics persisted',{'flags':info['flag_shifts'],'margin_source':'diagnostics000.json: margin (same realization)'},True)
    d=json.loads((root/'diagnostics000.json').read_text()); rr=d['random']; st=absolute_statistics([r['A'] for r in rr]); sd=np.std([r['A'] for r in rr],ddof=1)
    st.update(scatter_over_rms_jk=sd/np.sqrt(np.mean([r['jk_error']**2 for r in rr])),scatter_over_rms_sigma_F=sd/np.sqrt(np.mean([r['sigma_F']**2 for r in rr])))
    record(rows,7,'100 random templates','finite absolute statistics and ratios (report only)',st,np.all(np.isfinite([st['mean'],st['sem'],st['scatter_over_rms_jk'],st['scatter_over_rms_sigma_F']])))
    bench=benchmark(sl,.1,xi,cfg,root/'benchmark.json'); peak=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2
    record(rows,8,'benchmark and memory','24 Numba threads; positive throughput; peak <40 GB',{'benchmark':bench,'peak_gib':peak},bench['threads']==24 and bench['pixel_pairs_per_s']>0 and peak*1024**3<40e9)
    info.update(rows=rows,wall_s=time.perf_counter()-t,peak_memory_gb=peak); dump(root/'extras.json',info)
    return info


def run_campaign(root,out_root):
    campaign_started=time.perf_counter()
    import hashlib,datetime,shutil,importlib.metadata
    from numba import set_num_threads
    from validation_stats import absolute_statistics
    root=Path(root); out_root=Path(out_root); set_num_threads(24)
    if not (root/'provenance.json').exists():
        import subprocess,re
        started=time.perf_counter()
        test=subprocess.run([sys.executable,'-m','pytest',str(HERE/'tests'),'-q'],text=True,capture_output=True)
        (root/'pytest.log').write_text(test.stdout+test.stderr)
        match=re.search(r'(\d+) passed',test.stdout)
        preflight={'passed':int(match.group(1)) if match else 0,'failed':int(test.returncode!=0),
                   'wall_s':time.perf_counter()-started,'output':test.stdout,'test_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (HERE/'tests').glob('*.py')}}
        dump(root/'preflight.json',preflight)
        if test.returncode: raise RuntimeError('preflight tests failed; fresh validation was not started')

    dev=json.loads((root/'development.json').read_text())
    cfg=Config(scale=1,xi_smoothing=dev['xi_smoothing'],analytic_nk=6400)
    frozen={'config':vars(cfg),'development':dev,'GATES_sha256':hashlib.sha256((ROOT/'GATES.md').read_bytes()).hexdigest()}
    config_hash=hashlib.sha256(json.dumps(vars(cfg),sort_keys=True,default=serializable).encode()).hexdigest()
    if not (root/'frozen.json').exists(): dump(root/'frozen.json',frozen)
    provenance={'UTC_start':datetime.datetime.now(datetime.timezone.utc).isoformat(),'recovery_seeds':list(range(20)),
                'null_seeds':list(range(20,40)),'extension_rule':'all seeds 40-59 once iff either required residual bound95 >0.3',
                'GATES_sha256':frozen['GATES_sha256'],'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(HERE.glob('*.py'))},
                'versions':{name:importlib.metadata.version(name) for name in ('numpy','scipy','numba','healpy','h5py','camb','fitsio')},
                'python':sys.version,'config_sha256':config_hash,'development_seeds':dev['seeds'],'xi_smoothing':cfg.xi_smoothing,'xi_step':cfg.xi_step,'analytic_nk':cfg.analytic_nk,'threads':24}
    if (root/'provenance.json').exists():
        prior=json.loads((root/'provenance.json').read_text())
        if prior['source_sha256']!=provenance['source_sha256'] or prior['GATES_sha256']!=provenance['GATES_sha256'] or prior['config_sha256']!=config_hash:
            raise RuntimeError('source or gates changed after launch; cannot mix validation products')
    else: dump(root/'provenance.json',provenance)
    shutil.copyfile(ROOT/'GATES.md',root/'GATES.md')
    from grid_covariance import xi_from_mock_grid,sample_covariance
    devmock=load_mock(root/'dev100.h5')
    gridshape=devmock.attrs['grid']
    if isinstance(gridshape,str): gridshape=json.loads(gridshape)
    raw,cube=xi_from_mock_grid(cfg,gridshape,return_cube=True); raw.save(root/'raw_prediction_xi.h5')
    _,gd=sample_covariance(cube,raw.r_perp,2.,.5,256)
    take=(raw.r_perp>=10)&(raw.r_perp<=30)
    convergence=np.linalg.norm(raw.xi_rp[take]-gd[take])/np.linalg.norm(gd[take])
    import h5py
    with h5py.File(root/'raw_prediction_xi.h5','a') as f: f['grid_covariance_cube']=cube
    dump(root/'grid_covariance_convergence.json',{'nangle':128,'doubled':256,'derivative_norm_change':convergence})
    if convergence>=.01: raise RuntimeError('development-grid angular covariance quadrature did not converge; validation has not begun')
    del devmock
    fixed=make_bundles(load_mock(root/'dev100.h5')); save_bundles(root/'fixed_templates.h5',fixed)
    diags=[]
    for seed in range(40):
        p=root/f'diagnostics{seed:03d}.json'
        diags.append(json.loads(p.read_text()) if p.exists() else process_seed(seed,cfg,root,'recovery' if seed<20 else 'null',raw,fixed))
    recovery=absolute_statistics([d['fits']['A1_R1']['deprojected']['A'] for d in diags[:20]],1)
    null=absolute_statistics([d['fits']['A0_R1']['deprojected']['A'] for d in diags[20:]])
    extension=recovery['bound95']>.3 or null['bound95']>.3
    dump(root/'extension_decision.json',{'recovery':recovery,'null':null,'extend':extension,'seeds':list(range(40,60)) if extension else []})
    if extension:
        for seed in range(40,60):
            p=root/f'diagnostics{seed:03d}.json'
            diags.append(json.loads(p.read_text()) if p.exists() else process_seed(seed,cfg,root,'extension',raw,fixed))
    if not (root/'extras.json').exists(): extra_checks(cfg,root)
    dump(root/'campaign_timing.json',{'this_invocation_wall_s':time.perf_counter()-campaign_started,
                                    'initial_development_wall_s':dev.get('initial_scan',dev)['wall_s'],
                                    'expanded_development_scan_wall_s':dev['wall_s']-dev.get('initial_scan',{}).get('wall_s',0),
                                    'development_note':'initial and expanded development jobs overlapped; their wall times are not additive'})
    return rebuild(root,out_root)




def finalize_development(root):
    """Complete the declared smoothing scan from saved development realizations."""
    from xi_model import xi_from_counts
    from validation_stats import slope_statistics
    root=Path(root); cfg=Config(scale=1); start=time.perf_counter(); results=[]
    widths=(.76,.9,1.1,1.3,1.5,1.7,2.0)
    for seed in range(100,105):
        while True:
            progress=root/'development_progress.json'
            if progress.exists() and any(d['seed']==seed for d in json.loads(progress.read_text())): break
            time.sleep(10)
        m=load_mock(root/f'dev{seed}.h5'); s=m.sightlines
        basis,_=flat_sky_band_templates(m.maps['kappa_lya'],s.ra,s.dec,np.deg2rad(20)/m.maps['kappa_lya'].shape[0])
        pairs=find_pairs(s,cfg.r_perp_max/s.chi.min()); vals={str(w):[] for w in widths}; changes={}
        for A in (0,1,5,10):
            sl=sightlines_for_variant(m,A,False); hist=xi_from_data(sl,cfg)
            hist.save(root/f'development_fits{seed}.h5',f'counts_A{A}')
            for width in widths:
                c=cfg.copy(xi_smoothing=width); xi=xi_from_counts(*hist.counts,c)
                cat=cat_for(sl,xi,c,pairs); fit=fit_save(cat,basis,c,sl,root/f'development_fits{seed}.h5',f'width{width}/A{A}')
                vals[str(width)].append(fit['A'])
                if A==1:
                    fine=xi_from_counts(*hist.counts,c.copy(xi_step=.125)); af=fit_bundle(cat_for(sl,fine,c,pairs),basis,c.g1,sl)['A']
                    changes[str(width)]={'A025':fit['A'],'A0125':af,'relative':abs(fit['A']-af)/max(abs(af),1e-30)}
        results.append({'seed':seed,'widths':vals,'convergence':changes}); print('FINAL DEVELOPMENT',seed,vals,flush=True)
        dump(root/'development_scan_progress.json',results)
    stats={str(w):slope_statistics([r['widths'][str(w)] for r in results]) for w in widths}
    chosen=min(stats,key=lambda w:abs(stats[w]['mean']-1))
    prior=json.loads((root/'development.json').read_text())
    final={'seeds':list(range(100,105)),'xi_smoothing':float(chosen),'xi_step':.25,'analytic_nk':6400,
           'results':results,'slopes':stats,'wall_s':prior['wall_s']+time.perf_counter()-start,
           'peak_memory_gb':max(prior['peak_memory_gb'],resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2),
           'choice_rule':'minimum abs(mean slope-1) among .76,.9,1.1,1.3,1.5,1.7,2.0 on seeds 100-104 only',
           'initial_scan':prior}
    dump(root/'development.json',final); return final


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--phase',choices=['development','freeze','validation','rebuild'],required=True)
    parser.add_argument('--mock-root',default=str(__import__('pathlib').Path(__import__('os').environ.get('LYALENSER_DATA','/data/LyaLenser'))/'mocks'/'iteration3'))
    parser.add_argument('--output',default=str(ROOT/'report'))
    args=parser.parse_args()
    if args.phase=='development': development(args.mock_root)
    elif args.phase=='freeze': finalize_development(args.mock_root)
    elif args.phase=='rebuild': rebuild(args.mock_root,args.output)
    else:
        result=run_campaign(args.mock_root,args.output)
        print('\nACCEPTANCE SUMMARY')
        for r in result['acceptance']: print('PASS' if r['pass'] else 'FAIL',r['acceptance'],r['measured'])
        print('Wall times:',result['timing'],'Peak GiB:',result['peak_memory_gb'])
