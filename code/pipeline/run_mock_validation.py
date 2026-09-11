"""Stage-A mock acceptance runner (IMPLEMENTATION.md section 2.8)."""
from __future__ import annotations

import argparse,json,resource,time,sys
from functools import lru_cache
from pathlib import Path
import numpy as np
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
    return accumulate(sl,pairs,xi,cfg)


def midpoint_regions(cat,sl):
    reg=pair_midpoint_regions(cat,sl,8); nside=8
    if len(np.unique(reg))<30:
        reg=pair_midpoint_regions(cat,sl,16); nside=16
    return reg,nside,len(np.unique(reg))


def solve(F,y):
    try: return np.linalg.solve(F,y)
    except np.linalg.LinAlgError: return np.linalg.pinv(F,rcond=1e-12)@y


def common_science(result,nscience=3):
    """Collapse independently fitted science bands to one common amplitude."""
    diag=np.diag(result.F); threshold=max(float(np.max(np.abs(diag))),1e-300)*1e-11
    active=[i for i in range(nscience) if abs(diag[i])>threshold]
    nuisance=[i for i in range(nscience,len(diag)) if abs(diag[i])>threshold]
    if not active: return {"A":np.nan,"jk_error":np.nan,"sigma_F":np.nan,"jk":np.array([]),"active":[]}
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
            "sigma_F":float(np.sqrt(max(np.linalg.pinv(Fr,rcond=1e-12)[0,0],0))),
            "jk":jk,"active":active}


def fit_bundle(cat,bundle,g1,sl):
    reg,nside,nreg=midpoint_regions(cat,sl); result=amplitude(cat,bundle,g1,reg)
    summary=common_science(result); summary.update(result=result,nside_jk=nside,nregion=nreg)
    return summary


def subtract_bundles(a,b,name="deprojected"):
    return [Template(x.alpha-y.alpha,x.name,x.kind,Lmin=x.Lmin,Lmax=x.Lmax,
                     filter=x.filter,source=name) for x,y in zip(a,b)]


def radial_transfer_average(L,h,ratio):
    vals=[]
    for lo,hi in SCIENCE:
        w=(2*L+1)*cosine_band(L,lo,hi,10)**2
        vals.append(float(np.sum(w*h*ratio)/np.sum(w)) if np.sum(w)>0 else 0.)
    w=(2*L+1)*(((L<40)|(L>300)).astype(float))
    vals.append(float(np.sum(w*h*ratio)/np.sum(w)) if np.sum(w)>0 else 0.)
    return np.asarray(vals[:3]+[0.,0.,0.]+vals[3:])


@lru_cache(maxsize=8)
def _theory_grid(z1,z2,margin,dx,cref):
    lmax=int(np.ceil(np.sqrt(2)*np.pi/(dx/cref)))
    L=np.unique(np.r_[2.,np.linspace(2,lmax,640)])
    c1=float(chi_of_z(z1)); c2=float(chi_of_z(z2))
    S=spectra(z1=z1,z2=z2,zs=2.4,template_z1=float(z_of_chi(c1-margin)),
              template_z2=float(z_of_chi(c2+margin)),L_values=L)
    return L,S


def theory_for(mock):
    side=np.deg2rad(20*float(mock.attrs["scale"])); dx=float(mock.attrs["dx"])
    z1,z2=mock.sightlines.attrs["zmin"],mock.sightlines.attrs["zmax"]
    margin=float(mock.attrs["template_margin"])
    L,S=_theory_grid(float(z1),float(z2),margin,dx,float(mock.sightlines.attrs["chi_ref"]))
    return L,S,side


def make_bundles(mock):
    """Truth, CMB, matched, deprojected, and density template bases."""
    L,S,side=theory_for(mock); shape=mock.maps["kappa_CMB"].shape; pix=side/shape[0]
    noise=float(mock.attrs["cmb_noise_level"])
    hc=S["klkc"]/(S["kckc"]+noise); hm=S["skl"]/(S["ss_sig"]+S["shot_s"])
    denom=S["kckc"]+noise-2*S["skc"]+S["ss_sig"]+S["shot_s"]
    hdep=(S["klkc"]-S["skl"])/denom
    tc=flat_sky_wiener_transfer(shape,pix,L,S["klkc"],S["kckc"]+noise)
    td=flat_sky_wiener_transfer(shape,pix,L,S["klkc"]-S["skl"],denom)
    truth,_=flat_sky_band_templates(mock.maps["kappa_lya"],mock.sightlines.ra,mock.sightlines.dec,pix,
                                    science_bands=SCIENCE,taper=10,source="mock truth")
    cmb,_=flat_sky_band_templates(mock.maps["kappa_CMB"],mock.sightlines.ra,mock.sightlines.dec,pix,
                                  science_bands=SCIENCE,taper=10,source="masked noisy CMB",transfer=tc)
    density,_=flat_sky_band_templates(mock.maps["delta_L"],mock.sightlines.ra,mock.sightlines.dec,pix,
                                      science_bands=SCIENCE,taper=10,source="mock delta_L")
    nmatch=max(32,int(np.ceil(np.rad2deg(side)/0.0573))); mpix=side/nmatch
    matched_map,matched_mask,meta=matched_template_flat(
        mock.quasars,mock.randoms,lambda z:np.full_like(z,3.5,dtype=float),(nmatch,nmatch),mpix)
    tm=flat_sky_wiener_transfer((nmatch,nmatch),mpix,L,S["skl"],S["ss_sig"]+S["shot_s"])
    tmd=flat_sky_wiener_transfer((nmatch,nmatch),mpix,L,S["klkc"]-S["skl"],denom)
    matched,_=flat_sky_band_templates(matched_map,mock.sightlines.ra,mock.sightlines.dec,mpix,
                                      science_bands=SCIENCE,taper=10,source="clustered-quasar matched template",transfer=tm)
    cmb_dep,_=flat_sky_band_templates(mock.maps["kappa_CMB"],mock.sightlines.ra,mock.sightlines.dec,pix,
                                      science_bands=SCIENCE,taper=10,source="CMB part of deprojected template",transfer=td)
    matched_dep,_=flat_sky_band_templates(matched_map,mock.sightlines.ra,mock.sightlines.dec,mpix,
                                          science_bands=SCIENCE,taper=10,source="matched part of deprojected template",transfer=tmd)
    dep=subtract_bundles(cmb_dep,matched_dep)
    base_ratio=S["dkc"]/S["dd"]
    ratios={"cmb":radial_transfer_average(L,hc,base_ratio),
            "matched":radial_transfer_average(L,hm,base_ratio),
            "deprojected":np.zeros(7)}
    return {"truth":truth,"cmb":cmb,"matched":matched,"deprojected":dep,"density":density,
            "ratios":ratios,"matched_meta":meta,"L":L,"spectra":S,
            "matched_map":matched_map,"matched_mask":matched_mask}


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


def run(scale,seeds,out_root,mock_root,run_flag_shifts=True):
    if len(seeds)<20 and scale>=1: raise ValueError("full configuration requires at least 20 seeds")
    cfg=Config(scale=scale,xi_step=.5,seeds=tuple(seeds))
    out_root=Path(out_root); figs=out_root/"figures"; figs.mkdir(parents=True,exist_ok=True)
    mock_root=Path(mock_root); mock_root.mkdir(parents=True,exist_ok=True)
    rows=[]; timing={}; ttotal=time.perf_counter(); g1,mean_lens=kernel_product_g1()

    norm_scale=min(scale,.15); ncfg=cfg.copy(scale=norm_scale)
    t=time.perf_counter()
    base=generate_mock(ncfg,seeds[0],scale=norm_scale,A_true=0,response=False,n_los=100,
                       pixel_noise_power=0,variant_A_values=[0,.5,1,2])
    save_mock(base,mock_root/"normalisation_fiducial.h5"); timing["normalisation_mock_s"]=time.perf_counter()-t
    t=time.perf_counter(); xm=xi_from_model(ForestPower(model="kaiser"),ncfg,nk=1600,kmax=25); timing["xi_model_s"]=time.perf_counter()-t
    t=time.perf_counter(); xd=xi_from_data(sightlines_for_variant(base,0,False),ncfg); timing["xi_data_s"]=time.perf_counter()-t
    bench=benchmark(sightlines_for_variant(base,0,False),.15,xd,ncfg,mock_root/"benchmark.json")
    table_ok=np.isfinite(xd.xi).all() and np.isfinite(xm.xi).all() and bench["pixel_pairs_per_s"]>0
    record(rows,1,"Tables and benchmark","finite tables; benchmark reported",
           {"pixel_pairs_per_s":bench["pixel_pairs_per_s"],"threads":bench["threads"]},table_ok)
    bundles=make_bundles(base); p=find_pairs(base.sightlines,ncfg.r_perp_max/base.sightlines.chi.min())
    Agrid=[0,.5,1,2]; dataA=[]; modelA=[]; no_moment=[]
    for A in Agrid:
        sl=sightlines_for_variant(base,A,False)
        dataA.append(fit_bundle(cat_for(sl,xd,ncfg,p),bundles["truth"],g1,sl)["A"])
        modelA.append(fit_bundle(cat_for(sl,xm,ncfg,p),bundles["truth"],g1,sl)["A"])
        no_moment.append(fit_bundle(cat_for(sl,xd,ncfg,p),bundles["truth"],0.,sl)["A"])
    sd=float(np.polyfit(Agrid,dataA,1)[0]); sm=float(np.polyfit(Agrid,modelA,1)[0]); ratio=sd/sm
    record(rows,2,"Baseline absorption","data/model slopes agree to 3%",{"ratio":ratio,"data_slope":sd,"model_slope":sm},abs(ratio-1)<=.03)
    off=generate_mock(ncfg,seeds[0],scale=norm_scale,A_true=0,g_on=False,response=False,n_los=100,
                      pixel_noise_power=0,variant_A_values=[0,.5,1,2])
    offb=make_bundles(off); po=find_pairs(off.sightlines,ncfg.r_perp_max/off.sightlines.chi.min())
    offA=[]
    for A in Agrid:
        sl=sightlines_for_variant(off,A,False); offA.append(fit_bundle(cat_for(sl,xd,ncfg,po),offb["truth"],0.,sl)["A"])
    soff=float(np.polyfit(Agrid,offA,1)[0]); snom=float(np.polyfit(Agrid,no_moment,1)[0])
    predicted_nom=1+g1*(float(np.mean(base.sightlines.chi))-cfg.chi_ref)
    physical_ok=all(abs(x-1)<=.05 for x in (sd,soff)) and abs((snom/sd)-predicted_nom)<=.05
    record(rows,3,"Physical normalization","slopes 1 +/- 0.05; missing-moment shift predicted",
           {"g_on_moments":sd,"g_off":soff,"g_on_no_moments":snom,"predicted_no_moment_ratio":predicted_nom},physical_ok)
    inj=injection_test(sightlines_for_variant(base,0,False),xd,base.alpha_lya,[-2,-1,0,1,2],ncfg,
                       junk_alpha=bundles["truth"][-1].alpha)
    cerr=float(np.nanmedian(inj["errors"][:,1])); inj_ok=abs(inj["paired_slope"]-1)<=.05 and abs(inj["curl_slope"])<=cerr
    record(rows,5,"Injection bookkeeping","slope 1 +/- 0.05; curl 0 within error",
           {"slope":inj["paired_slope"],"curl_slope":inj["curl_slope"],"curl_error":cerr},inj_ok)

    null={"truth":[],"cmb":[],"matched":[]}; null_raw={"truth":[],"cmb":[],"matched":[]}
    stochastic=[]; response_obs={"cmb":[],"matched":[],"deprojected":[]}; response_pred={k:[] for k in response_obs}
    combined=[]; combined_jk=[]; shape_values=[]; spectrum_ratios=[]; nregions=[]; first_random=None
    t=time.perf_counter()
    for seed in seeds:
        m=generate_mock(cfg,seed,scale=scale,A_true=1,response=True,magnification=True,
                        completeness=True,real_mask=True,cmb_noise=True,variant_A_values=[0,1])
        save_mock(m,mock_root/f"mock_seed{seed:03d}.h5")
        b=make_bundles(m); pairs=find_pairs(m.sightlines,cfg.r_perp_max/m.sightlines.chi.min())
        sl00=sightlines_for_variant(m,0,False); sl01=sightlines_for_variant(m,0,True)
        sl10=sightlines_for_variant(m,1,False); sl11=sightlines_for_variant(m,1,True)
        c00=cat_for(sl00,xd,cfg,pairs); c01=cat_for(sl01,xd,cfg,pairs)
        c10=cat_for(sl10,xd,cfg,pairs); c11=cat_for(sl11,xd,cfg,pairs); c11.save(mock_root/f"catalogue_seed{seed:03d}.h5")
        fits00={}
        for name in null:
            f=fit_bundle(c00,b[name],g1,sl00); null[name].append(f["A"])
            fits00[name]=f
            rr=f["result"]; old=rr.mf.copy(); rr.mf[:]=0; null_raw[name].append(common_science(rr)["A"]); rr.mf[:]=old
        fits00["deprojected"]=fit_bundle(c00,b["deprojected"],g1,sl00)
        stochastic.append(fit_bundle(c10,b["cmb"],g1,sl10)["A"]-fits00["cmb"]["A"])
        fits01={}
        for name in response_obs:
            fitr=fit_bundle(c01,b[name],g1,sl01); response_obs[name].append(fitr["A"])
            fits01[name]=fitr
            response_obs[name][-1]-=fits00[name]["A"]
            pred=independent_response_prediction(c01,b[name],b["density"],m.truth["delta_L"],b["ratios"][name],cfg.response_delta,g1)
            rr=fitr["result"]; oldq=rr.q.copy(); oldmf=rr.mf.copy(); rr.q[:]=pred["scores"]; rr.mf[:]=0
            response_pred[name].append(common_science(rr)["A"]); rr.q[:]=oldq; rr.mf[:]=oldmf
        fc=fit_bundle(c11,b["deprojected"],g1,sl11); combined.append(fc["A"]-fits01["deprojected"]["A"])
        jkd=fc["jk"]-fits01["deprojected"]["jk"]; nr=len(jkd)
        combined_jk.append(float(np.sqrt((nr-1)/nr*np.sum((jkd-jkd.mean())**2))))
        vals1,_=shape_consistency(c11,b["deprojected"],g1,sl11)
        vals0,_=shape_consistency(c01,b["deprojected"],g1,sl01); shape_values.append(vals1-vals0)
        spectrum_ratios.append(mock_spectrum_check(m)["ratio"]); nregions.append(fc["nregion"])
        if first_random is None: first_random=(m,c00,b,sl00)
    timing["full_ensemble_s"]=time.perf_counter()-t

    mean_field_extra=[]; tmf=time.perf_counter()
    if scale>=1 and len(seeds)==20:
        seed=max(seeds)+1; mean_field_extra.append(seed)
        m=generate_mock(cfg,seed,scale=scale,A_true=0,response=False,magnification=True,
                        completeness=True,real_mask=True,cmb_noise=True)
        save_mock(m,mock_root/f"mock_seed{seed:03d}_meanfield.h5")
        b=make_bundles(m); c=cat_for(m.sightlines,xd,cfg)
        c.save(mock_root/f"catalogue_seed{seed:03d}_meanfield.h5")
        for name in null:
            f=fit_bundle(c,b[name],g1,m.sightlines); null[name].append(f["A"])
            rr=f["result"]; old=rr.mf.copy(); rr.mf[:]=0
            null_raw[name].append(common_science(rr)["A"]); rr.mf[:]=old
    timing["mean_field_extra_s"]=time.perf_counter()-tmf

    mfmetrics={}; mfok=True
    for name in null:
        mu,se=mean_sem(null[name]); rmu,_=mean_sem(null_raw[name]); mfmetrics[name]={"mean":mu,"sem":se,"raw_mean":rmu}; mfok &= abs(mu)<=2*se
    record(rows,4,"Mean field","all three means 0 within 2 SEM",mfmetrics,mfok,"raw means show the mf/F subtraction")
    ms,es=mean_sem(stochastic); record(rows,6,"Stochastic normalization","mean 1 within 2 SEM",{"mean":ms,"sem":es},abs(ms-1)<=2*es)
    respmetrics={}; respok=True
    for name in response_obs:
        obs=np.asarray(response_obs[name]); pre=np.asarray(response_pred[name]); diff=obs-pre
        mo,eo=mean_sem(obs); mp,_=mean_sem(pre); md,ed=mean_sem(diff)
        respmetrics[name]={"observed":mo,"observed_sem":eo,"predicted":mp,"difference":md,"difference_sem":ed}
        respok &= (abs(mo)<=2*eo) if name=="deprojected" else (abs(md)<=2*ed)
    record(rows,7,"Response-only null","predictions agree within errors; deprojected is 0",respmetrics,respok)
    mc,ec=mean_sem(combined)

    shifts={}; flag_time=time.perf_counter()
    if run_flag_shifts:
        base_flags=dict(magnification=True,completeness=True,real_mask=True,template_margin=150.)
        for label,change in (("magnification_off",{"magnification":False}),("completeness_off",{"completeness":False}),
                             ("real_mask_off",{"real_mask":False}),("no_margin",{"template_margin":0.})):
            mm=generate_mock(cfg,seeds[0],scale=scale,A_true=1,response=True,cmb_noise=True,
                             variant_A_values=[0,1],**{**base_flags,**change})
            bb=make_bundles(mm); pp=find_pairs(mm.sightlines,cfg.r_perp_max/mm.sightlines.chi.min())
            ss1=sightlines_for_variant(mm,1,True); ss0=sightlines_for_variant(mm,0,True)
            inc=fit_bundle(cat_for(ss1,xd,cfg,pp),bb["deprojected"],g1,ss1)["A"]-fit_bundle(cat_for(ss0,xd,cfg,pp),bb["deprojected"],g1,ss0)["A"]
            shifts[label]=inc-combined[0]
    timing["flag_shifts_s"]=time.perf_counter()-flag_time
    record(rows,8,"Combined recovery","mean 1 within 2 SEM",{"mean":mc,"sem":ec,"switch_off_shifts":shifts},abs(mc-1)<=2*ec)

    sv=np.asarray(shape_values); shape_mean=sv.mean(axis=0); shape_sem=sv.std(axis=0,ddof=1)/np.sqrt(len(sv)); ok=np.isfinite(shape_sem)&(shape_sem>0)
    if ok.sum()>1:
        from scipy.stats import chi2 as chi2dist
        wm=np.sum(shape_mean[ok]/shape_sem[ok]**2)/np.sum(1/shape_sem[ok]**2); chi2=float(np.sum(((shape_mean[ok]-wm)/shape_sem[ok])**2)); pshape=float(chi2dist.sf(chi2,ok.sum()-1))
    else: chi2=np.nan; pshape=np.nan
    record(rows,9,"Shape test","six-bin p-value > 0.01",{"p_value":pshape,"chi2":chi2,"values":shape_mean.tolist()},pshape>.01)

    m,cbase,b0,slbase=first_random; L,S,side=theory_for(m); rng=np.random.default_rng(9912); aval=[]; ajk=[]; af=[]
    noise=float(m.attrs["cmb_noise_level"]); pix=side/m.maps["kappa_CMB"].shape[0]
    transfer=flat_sky_wiener_transfer(m.maps["kappa_CMB"].shape,pix,L,S["klkc"],S["kckc"]+noise)
    for _ in range(100):
        km=gaussian_flat_map(m.maps["kappa_CMB"].shape,side,L,S["kckc"]+noise,rng)
        tb,_=flat_sky_band_templates(km,slbase.ra,slbase.dec,pix,science_bands=SCIENCE,taper=10,source="Gaussian CMB",transfer=transfer)
        fr=fit_bundle(cbase,tb,g1,slbase); aval.append(fr["A"]); ajk.append(fr["jk_error"]); af.append(fr["sigma_F"])
    std=float(np.std(aval,ddof=1)); medjk=float(np.nanmedian(ajk)); medf=float(np.nanmedian(af))
    randmetrics={"std":std,"std_over_jackknife":std/medjk,"std_over_sigma_F":std/medf}
    record(rows,10,"Random-template ensemble","ratios reported (no tolerance)",randmetrics,np.all(np.isfinite(list(randmetrics.values()))))
    scatter=float(np.std(combined,ddof=1)); rmscombinedjk=float(np.sqrt(np.nanmean(np.asarray(combined_jk)**2))); covratio=scatter/rmscombinedjk
    record(rows,11,"Covariance","mock scatter / jackknife = 1 +/- 0.3",{"ratio":covratio,"scatter":scatter,"rms_jk":rmscombinedjk},abs(covratio-1)<=.3)

    specmean=np.mean(np.asarray(spectrum_ratios),axis=0); specok=bool(np.all(np.abs(specmean-1)<=.10))
    rows[0]["measured"]["spectrum_ratios_klkl_klkc_kckc"]=specmean.tolist(); rows[0]["pass"] &= specok
    rows[0]["tolerance"] += "; mock spectra within 10% over 40<L<300"
    rows.sort(key=lambda row: row["item"])

    fig,ax=plt.subplots(figsize=(5.5,4)); ax.plot(Agrid,dataA,"o-",label="xi from data"); ax.plot(Agrid,modelA,"s-",label="xi model"); ax.plot(Agrid,offA,"^-",label="g off"); ax.set(xlabel=r"$A_{true}$",ylabel=r"$\hat A$",title="Physical normalization"); ax.legend(); fig.tight_layout(); fig.savefig(figs/"mock_physical_normalisation.pdf"); plt.close(fig)
    fig,ax=plt.subplots(figsize=(5.5,4)); ax.hist(combined,bins=max(5,int(np.sqrt(len(combined))))); ax.axvline(1,color="k",ls="--"); ax.set(xlabel="deprojected recovered amplitude",ylabel="mocks",title="Full ensemble"); fig.tight_layout(); fig.savefig(figs/"mock_combined_ensemble.pdf"); plt.close(fig)
    fig,ax=plt.subplots(figsize=(5.5,4)); ax.errorbar(np.arange(6),shape_mean,yerr=shape_sem,fmt="o"); ax.set(xlabel="shape bin",ylabel=r"$\hat A$",title=f"Shape consistency, p={pshape:.3g}"); fig.tight_layout(); fig.savefig(figs/"mock_shape.pdf"); plt.close(fig)
    fig,ax=plt.subplots(figsize=(5.5,4)); ax.bar(["kl-kl","kl-kCMB","kCMB-kCMB"],specmean); ax.axhspan(.9,1.1,color="0.8"); ax.set(ylabel="mock / three_tracer",title="Patch spectra, 40 < L < 300"); fig.tight_layout(); fig.savefig(figs/"mock_spectra.pdf"); plt.close(fig)

    timing["total_s"]=time.perf_counter()-ttotal; peak_gb=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2
    result={"configuration":{"scale":scale,"normalisation_scale":norm_scale,"seeds":list(seeds),"threads":24,"mock_root":str(mock_root),"g1":g1,"g1_mean_lens_chi":mean_lens,"mean_field_extra_seeds":mean_field_extra,"jackknife_nsides":[8,16],"jackknife_region_range":[int(min(nregions)),int(max(nregions))]},"timing":timing,"peak_memory_gb":peak_gb,"acceptance":rows,"raw":{"Agrid":Agrid,"data_A":dataA,"model_A":modelA,"g_off_A":offA,"combined":combined,"combined_jk":combined_jk,"shape_values":sv.tolist(),"spectrum_ratios":np.asarray(spectrum_ratios).tolist()}}
    (out_root/"mock_validation.json").write_text(json.dumps(result,indent=2,default=float)+"\n")
    lines=["# Stage A mock validation","",f"Run: `{len(seeds)}` seeds, scale `{scale}`, 24 Numba threads."]
    if mean_field_extra: lines.append(f"Mean-field null: `{len(seeds)+len(mean_field_extra)}` seeds (extra {mean_field_extra}).")
    lines += [f"Kernel-product `g1 = {g1:.8g} (Mpc/h)^-1` (`<chi_l> = {mean_lens:.3f} Mpc/h`).","","| Sec. 2.8 | Acceptance item | Tolerance | Measured | Result |","|---:|---|---|---|:---:|"]
    for r in rows:
        measured=json.dumps(r["measured"],default=float,separators=(",",":")); lines.append(f'| {r["item"]} | {r["acceptance"]} | {r["tolerance"]} | `{measured}` | {"PASS" if r["pass"] else "FAIL"} |')
    lines += ["","## Wall times and memory",""]+[f"- `{k}`: {v:.3f} s" for k,v in timing.items()]
    lines += [f"- Peak resident memory: {peak_gb:.3f} GB.","","## Remaining failures",""]
    failures=[r for r in rows if not r["pass"]]; lines += ([f'- Section 2.8 item {r["item"]}, {r["acceptance"]}: {r["detail"] or "measured tolerance was not met"}.' for r in failures] or ["- None."])
    (out_root/"mock_validation.md").write_text("\n".join(lines)+"\n"); return result


if __name__=="__main__":
    ap=argparse.ArgumentParser(); ap.add_argument("--scale",type=float,default=.15); ap.add_argument("--seeds",type=int,nargs="*",default=list(range(20))); ap.add_argument("--output",default=str(ROOT/"report")); ap.add_argument("--mock-root",default="/data/LyaLenser/mocks"); ap.add_argument("--skip-flag-shifts",action="store_true")
    args=ap.parse_args(); result=run(args.scale,args.seeds,args.output,args.mock_root,not args.skip_flag_shifts)
    print(json.dumps({"passed":sum(x["pass"] for x in result["acceptance"]),"failed":sum(not x["pass"] for x in result["acceptance"]),"wall_s":result["timing"]["total_s"],"peak_memory_gb":result["peak_memory_gb"]},indent=2))
