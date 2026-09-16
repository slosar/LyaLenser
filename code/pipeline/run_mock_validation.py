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
from mock import generate_mock,save_mock,sightlines_for_variant,mock_spectrum_check,patch_side_rad
from xi_model import xi_from_model,xi_from_data
from pairs import find_pairs,accumulate,benchmark,pair_midpoint_regions
from templates import (Template,cosine_band,flat_sky_band_templates,
                       flat_sky_wiener_transfer,matched_template_flat)
from amplitude import amplitude,independent_response_prediction
from inject import injection_test
from three_tracer import spectra

SCIENCE=((40,100),(100,200),(200,300))


def record(rows,item,name,tolerance,metrics,passed,detail="",required=None):
    """A row is required (blocks Stage B) unless its tolerance declares it report-only or ``required`` says so."""
    if required is None: required=not tolerance.lower().startswith("report")
    rows.append({"item":item,"acceptance":name,"tolerance":tolerance,"measured":metrics,
                 "pass":bool(passed),"required":bool(required),"detail":detail})


def mean_sem(values):
    x=np.asarray(values,float)
    return float(np.mean(x)),float(np.std(x,ddof=1)/np.sqrt(len(x))) if len(x)>1 else np.nan


def cat_for(sl,xi,cfg,pairs=None):
    if pairs is None: pairs=find_pairs(sl,cfg.r_perp_max/max(float(sl.chi.min()),1))
    start=time.perf_counter(); cat=accumulate(sl,pairs,xi,cfg)
    cat.attrs["wall_s"]=time.perf_counter()-start
    return cat


JK_MIN_REGIONS=30

def midpoint_regions(cat,sl):
    """Pair-midpoint HEALPix jackknife regions: nside 8, refined to 16 when fewer than JK_MIN_REGIONS regions are
    populated (the 20-degree scale-1 patch: ~27 regions at nside 16). The chosen nside and the region count are
    recorded with every fit; the campaign reports them (GATES v6)."""
    reg=pair_midpoint_regions(cat,sl,8); nside=8
    if len(np.unique(reg))<JK_MIN_REGIONS:
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
    # The generator's patch side (nx cells of dx at chi_ref), so that the 128-pixel operators, the quasar
    # pixelisation and the sightline coordinates share one angular scale (review 5, finding 2).
    side=patch_side_rad(mock); n=128; shape=(n,n); pix=side/n
    margin=float(mock.attrs['template_margin']) if template_margin is None else template_margin
    c1=float(chi_of_z(mock.sightlines.attrs['zmin'])); c2=float(chi_of_z(mock.sightlines.attrs['zmax']))
    def selection(cat):
        keep=(cat['chi']>=c1-margin)&(cat['chi']<=c2+margin)
        return {k:np.asarray(v)[keep] for k,v in cat.items()}
    q=selection(mock.quasars); r=selection(mock.randoms)
    matched,mm,meta=matched_template_flat(q,r,lambda z:np.full_like(z,3.5),shape,pix,radial_range=(c1-margin,c2+margin))
    def resize(a): return resample(resample(a,n,axis=0),n,axis=1).real
    act=zoom(mock.maps['act_mask'],n/mock.maps['act_mask'].shape[0],order=1)
    mask=np.clip(act,0,1)*mm
    # CMB saved unmasked as well: common mask is applied exactly once.
    cmbmap=resize(mock.maps['kappa_CMB_unmasked'])*mask
    matched_unmasked=matched; matched=matched*mask
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
    out.update(matched_meta=meta,matched_map=matched,matched_map_unmasked=matched_unmasked,matched_mask=mm,common_mask=mask,
               used_maps=used,operators=operators,spectra_used=spectra_used,L=L,spectra=S,
               truth_map=truthmap,pixel_size=pix,side_rad=side,magcoef=magcoef)
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
        for name in ('common_mask','matched_map','matched_map_unmasked','L'): g[name]=bundles[name]
        g.attrs['side_rad']=float(bundles['side_rad']); g.attrs['pixel_size']=float(bundles['pixel_size'])
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
        names=['truth','quasars','randoms','maps']+[n for n in ('lowz_catalogue','lowz_randoms') if n in f]
        ds={name:{k:v[()] for k,v in f[name].items()} for name in names}
        attrs=dict(f.attrs)
    if 'lowz' in attrs and isinstance(attrs['lowz'],str): attrs['lowz']=json.loads(attrs['lowz'])
    return MockResult(load_sightlines(path),ds['truth']['alpha_lya'],ds['truth'],ds['quasars'],ds['randoms'],ds['maps'],attrs,
                      ds.get('lowz_catalogue'),ds.get('lowz_randoms'))

def fit_save(cat,bundle,cfg,sl,path,group):
    result=fit_bundle(cat,bundle,cfg.g1,sl); result['result'].save(path,group)
    rr=result['result']; raw=common_science(__import__('dataclasses').replace(rr,mf=np.zeros_like(rr.mf)))['A']
    return {k:v for k,v in result.items() if k!='result'}|{'raw':raw,'bands':rr.A}

def table_for(sl,cfg,basis,counts=None):
    """Measured 1 Mpc/h counts of this sample -> fitted response table.

    ``cfg.xi_correction`` selects iteration 5 (two-parameter Kaiser fit of the projected model basis) or
    iteration 8 (the same fit plus a spline correction and a same-wavelength term, `xi_spline`)."""
    num,den=(xi_from_data(sl,cfg).counts if counts is None else counts)
    if getattr(cfg,'xi_correction','none')=='spline':
        from xi_spline import Envelope,XiCorrection,fit_corrected_table
        from xi_fit import FittedTable,BASIS
        from xi_model import XiTable
        b=basis['projected']; c=[1.,2*1.4,1.4**2]
        ref=XiTable(b[BASIS[0]].r_perp,b[BASIS[0]].r_par,sum(x*b[k].xi for x,k in zip(c,BASIS)),
                    sum(x*b[k].xi_rp for x,k in zip(c,BASIS)),{})
        KN={'default':((3.,6.,10.,16.,30.),(0.,4.,10.,30.)),
            'coarse':((3.,8.,16.,30.),(0.,6.,30.)),
            'small':((3.,30.),(0.,30.))}[getattr(cfg,'xi_knots','default')]
        corr=XiCorrection(Envelope(ref),rp_knots=KN[0],rz_knots=KN[1],
                          rz_sw=1. if getattr(cfg,'xi_same_wavelength',True) else 0.)
        tab,par,_=fit_corrected_table(num,den,b,basis['coarse'],cfg,corr,ridge=cfg.xi_correction_ridge)
        return FittedTable(tab,par)
    from xi_fit import fit_model_table
    return fit_model_table(num,den,basis['projected'],cfg,basis['coarse'])


GRID_ROLES=('recovery','extension','full')          # carry the A grid {0, .5, 1, 2} (slopes, fixed-baseline check)
NULL_ROLES=('null','extension','full','core')       # carry the fixed-template null fits
DEPROJECTED_BOUND=0.5                               # GATES v6: 95 % residual bound on the deprojected rows, in A


def process_seed(seed,cfg,root,role,raw_xi,fixed_bundle,basis,extras=None):
    """Roles: iteration 6 uses 'full' (A grid + every null extra) and 'core' (A in {0, 1} + null extras); every seed
    carries A0_R1/A1_R1 with the response on, which is what the deprojected rows use. 'recovery', 'null' and
    'extension' are the iteration 4-5 roles and stay valid."""
    from response import prediction_catalogue
    t=time.perf_counter(); path=root/f'seed{seed:03d}.h5'; fits_path=root/f'fits{seed:03d}.h5'
    recovery=role in GRID_ROLES; null=role in NULL_ROLES
    if not (recovery or null): raise ValueError(f'unknown seed role {role!r}')
    # Margin and 100-random-template diagnostics run on one designated seed (seed 0 in iterations 4-5).
    extras=(seed==0) if extras is None else bool(extras)
    avals=[0,.5,1,2] if recovery else [0,1]
    m=generate_mock(cfg,seed,A_true=1,response=True,magnification=True,completeness=True,
                    real_mask=True,cmb_noise=True,variant_A_values=avals,disjoint_selection=True)
    save_mock(m,path); b=make_bundles(m); save_bundles(fits_path,b)
    pairs=find_pairs(m.sightlines,cfg.r_perp_max/m.sightlines.chi.min())
    diag={'seed':seed,'role':role,'fits':{},'slopes':{},'predictions':{},'baseline_fixed':[], 'shape':[],'xi_fit':{},
          'config':vars(cfg),'mock_attrs':m.attrs,'files':{'mock':str(path),'fits':str(fits_path)}}
    ffixed=table_for(sightlines_for_variant(m,0,False),cfg,basis); ffixed.save(fits_path,'xi_fixed'); xfixed=ffixed.table
    todo=[(A,False) for A in avals]+[(0,True),(1,True)]
    for A,resp in todo:
        key=f'A{A:g}_R{int(resp)}'; sl=sightlines_for_variant(m,A,resp)
        ft=table_for(sl,cfg,basis); ft.save(fits_path,f'xi/{key}'); xi=ft.table; diag['xi_fit'][key]=ft.params
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
        if extras and A==1 and resp:
            # Margin changes use the SAME density, sightlines, CMB map and forest.
            diag['margin']={}
            for margin in (0.,150.,300.):
                bm=make_bundles(m,template_margin=margin); save_bundles(fits_path,bm,f'margin_templates/{margin}')
                diag['margin'][str(margin)]=fit_save(cat,bm['deprojected'],cfg,sl,fits_path,f'margin/{margin}')['A']
            # Independent Gaussian maps receive the identical common mask and transfer.
            rng=__import__('random_streams').seed_streams(seed)['random_templates']; diag['random']=[]
            L=b['L']; S=b['spectra']; shape=b['common_mask'].shape; side=b['side_rad']
            for i in range(100):
                km=gaussian_flat_map(shape,side,L,S['kckc']+m.attrs['cmb_noise_level'],rng)*b['common_mask']
                tb,_=flat_sky_band_templates(km,sl.ra,sl.dec,b['pixel_size'],transfer=b['operators']['cmb'],source='independent masked Gaussian')
                diag['random'].append(fit_save(cat,tb,cfg,sl,fits_path,f'random/{i}'))
                import h5py
                with h5py.File(fits_path,'a') as f: f[f'random_maps/{i}']=km
    from template_audit import audit_mock
    diag['template_coefficients']=audit_mock(m,fits_path,make_bundles)
    # Identical forest, maps and randoms; only catalogue selection changes.
    from dataclasses import replace
    q={k:m.truth['shared_'+k] for k in ('ra','dec','chi','qid')}
    q['z']=z_of_chi(q['chi'])
    shared=replace(m,quasars=q)
    sb=make_bundles(shared); save_bundles(fits_path,sb,'shared_templates')
    diag['shared']={}
    for A,resp in ((0,False),(0,True),(1,True)):
        sl=sightlines_for_variant(m,A,resp); xi=table_for(sl,cfg,basis).table; cat=cat_for(sl,xi,cfg,pairs)
        key=f'A{A}_R{int(resp)}'
        diag['shared'][key]={name:fit_save(cat,sb[name],cfg,sl,fits_path,f'shared/{key}/{name}')
                             for name in ('matched','deprojected')}
    diag['spectra']=mock_spectrum_check(m)
    diag['wall_s']=time.perf_counter()-t; diag['peak_memory_gb']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2
    dump(root/f'diagnostics{seed:03d}.json',diag)
    print('SEED COMPLETE',seed,role,'seconds',diag['wall_s'],'Acombined',diag['fits']['A1_R1']['deprojected']['A'],flush=True)
    return diag


def build_rows(diags,extras,bound=DEPROJECTED_BOUND):
    """Acceptance rows (GATES v6). Slope rows use the A-grid seeds; the mean-field rows the seeds with fixed-template
    null fits; the deprojected, response-only, covariance, shape and spectra rows use EVERY seed (all carry A0_R1
    and A1_R1), so null and recovery summaries share the realisations (declared in GATES v6)."""
    from validation_stats import absolute_statistics,slope_statistics,hotelling_shape
    rec=[d for d in diags if d['role'] in GRID_ROLES]
    null=[d for d in diags if d['role'] in NULL_ROLES]
    every=list(diags)
    rows=[]
    def row(name,tol,stats,ok,detail='',required=None): record(rows,len(rows)+1,name,tol,stats,ok,detail,required)
    for name in ('truth','cmb','matched'):
        vals=[[d['fits'][f'A{A:g}_R0'][name]['A'] for A in (0,.5,1,2)] for d in rec]
        st=slope_statistics(vals)
        row(f'{name} normalization slope','|slope-1| <= 0.03',st,abs(st['residual'])<=.03,
            'Slope of absolute amplitudes at A_true={0,0.5,1,2}; per-seed intercept fitted freely.')
    fixed=slope_statistics([d['baseline_fixed'] for d in rec])
    refit=slope_statistics([[d['fits'][f'A{A:g}_R0']['truth']['A'] for A in (0,.5,1,2)] for d in rec])
    diff=absolute_statistics(np.array(refit['slopes'])-np.array(fixed['slopes']))
    diff['absolute_A1_refit_minus_fixed']=absolute_statistics([d['fits']['A1_R0']['truth']['A']-d['baseline_fixed'][2] for d in rec])
    row('fixed versus refitted baseline','|slope difference| <= 0.03',{'fixed_absolute_slope':fixed,'refitted_absolute_slope':refit,'diagnostic_difference':diff},abs(diff['mean'])<=.03)
    for name in ('truth','cmb','matched'):
        for conditional in (False,True):
            vals=[(d['fixed_null'] if conditional else d['fits']['A0_R0'])[name]['A'] for d in null]
            raw=[(d['fixed_null'] if conditional else d['fits']['A0_R0'])[name]['raw'] for d in null]
            st=absolute_statistics(vals); st['raw']=absolute_statistics(raw)
            row(f'{"fixed" if conditional else "varying"} template mean field: {name}','|mean| <= 2 SEM',st,abs(st['mean'])<=2*st['sem'])
    st=absolute_statistics([d['fits']['A1_R0']['cmb']['A'] for d in every],1)
    row('absolute stochastic recovery','|mean-1| <= 2 SEM (precision established separately by slope)',st,abs(st['residual'])<=2*st['sem'])
    for name in ('cmb','matched','deprojected'):
        obs=np.array([d['fits']['A0_R1'][name]['A'] for d in every]); pred=np.array([d['predictions'][name]['A'] for d in every])
        st=absolute_statistics(obs); st['prediction']=absolute_statistics(pred); st['prediction_difference']=absolute_statistics(obs-pred)
        ok=abs(st['prediction_difference']['mean'])<=2*st['prediction_difference']['sem']
        tol='|observed-predicted| <= 2 SEM'
        if name=='deprojected':
            # GATES v6: the bound is chosen with the ensemble size (N = 400, per-seed scatter ~3.7 A -> t SEM ~0.37 A).
            tol+=f'; |mean| <= 2 SEM and absolute null bound95 <= {bound:g} A'; ok=ok and abs(st['mean'])<=2*st['sem'] and st['bound95']<=bound
        row(f'absolute response-only: {name}',tol,st,ok,
            'Prediction uses stored delta_L at every pixel and P_a (D_a C_ab D_b - C_ab) P_b^T with the production pair selection. Discrete-grid, phase-averaged trilinear covariance truncated at saved table support.')
    st=absolute_statistics([d['fits']['A1_R1']['deprojected']['A'] for d in every],1)
    row('absolute combined deprojected recovery',f'|mean-1| <= 2 SEM and residual bound95 <= {bound:g} A',st,abs(st['residual'])<=2*st['sem'] and st['bound95']<=bound,
        'Same seeds as the response-only null (every seed carries A_true = 0 and 1 with the response on).')
    st=absolute_statistics([d['fits']['A1_R1']['deprojected']['A']-d['fits']['A0_R1']['deprojected']['A'] for d in every],1)
    row('paired deprojected response A(1) - A(0), same realisation','|mean-1| <= 2 SEM',st,abs(st['residual'])<=2*st['sem'],
        'Sample variance of the deprojection cancels in the pair; this is the deprojected normalisation.')
    st=hotelling_shape([d['shape'] for d in every]); row('six-bin Hotelling shape','p > 0.01',st,st['p_value']>.01)
    scatter=np.std([d['fits']['A1_R1']['deprojected']['A'] for d in every],ddof=1)
    rms=np.sqrt(np.mean([d['fits']['A1_R1']['deprojected']['jk_error']**2 for d in every]))
    jk={'nside':sorted({int(d['fits']['A1_R1']['deprojected'].get('nside_jk',0)) for d in every}),
        'regions_mean':float(np.mean([d['fits']['A1_R1']['deprojected'].get('nregion',np.nan) for d in every]))}
    row('absolute covariance','0.7 <= scatter/RMS jackknife <= 1.3',{'scatter':scatter,'rms_jk':rms,'ratio':scatter/rms,'jackknife':jk},.7<=scatter/rms<=1.3,
        'Midpoint HEALPix jackknife; nside and populated region count as recorded per fit.')
    sv=np.array([d['spectra']['ratio'] for d in every]); st={'ratios':sv.mean(axis=0).tolist(),'sem':(sv.std(axis=0,ddof=1)/np.sqrt(len(sv))).tolist()}
    row('outside-box and total spectra','each total ratio within 10%',st,np.all(abs(np.array(st['ratios'])-1)<=.1))
    for r in extras.get('rows',[]):
        r.setdefault('required',not str(r.get('tolerance','')).lower().startswith('report')); rows.append(r)
    return rows


if __name__=='__main__':
    from campaign4 import main
    main()
