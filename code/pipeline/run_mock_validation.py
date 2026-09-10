"""Run and document the Stage-A mock acceptance suite.

Standalone usage from this directory::
  /home/anze/anaconda3/bin/python3 run_mock_validation.py --scale .03 --seeds 0 1
"""
from __future__ import annotations

import argparse,json,time,sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent; CODE=HERE.parent; ROOT=CODE.parent
if str(CODE) not in sys.path: sys.path.insert(0,str(CODE))
from forest_power import ForestPower
from config import Config
from mock import generate_mock,save_mock,alpha_from_map
from xi_model import xi_from_model,xi_from_data
from pairs import find_pairs,accumulate,benchmark
from templates import Template,curl
from amplitude import amplitude,shape_test,random_ensemble
from inject import injection_test


def cat_for(m,x,cfg):
    p=find_pairs(m.sightlines,cfg.r_perp_max/max(float(m.sightlines.chi.min()),1))
    return accumulate(m.sightlines,p,x,cfg)


def regions(cat,n=12):
    # Reduced patches occupy too few nside=8 pixels for a jackknife.  The full
    # geometry path uses HEALPix; reduced validation uses deterministic pair
    # groups and records this in the output.
    return np.arange(len(cat.a),dtype=np.int32)%min(n,max(len(cat.a),1))


def fit(cat,alpha,g1):
    ts=[Template(alpha,"gradient","truth"),Template(curl(alpha),"curl","curl")]
    return amplitude(cat,ts,g1,regions(cat))


def slope(x,y): return float(np.polyfit(np.asarray(x,float),np.asarray(y,float),1)[0])


def record(rows,item,name,tolerance,value,passed,detail=""):
    rows.append({"item":item,"test":name,"tolerance":tolerance,"measured":value,
                 "pass":None if passed is None else bool(passed),"detail":detail})


def run(scale,seeds,out_root,mock_root):
    cfg=Config(scale=scale,xi_step=.5,nside_alpha=64,lmax_alpha=128,seeds=tuple(seeds))
    norm_scale=min(scale,.05)
    norm_cfg=cfg.copy(scale=norm_scale)
    out_root=Path(out_root); figs=out_root/"figures"; figs.mkdir(parents=True,exist_ok=True)
    mock_root=Path(mock_root); mock_root.mkdir(parents=True,exist_ok=True)
    rows=[]; timing={}; ttotal=time.perf_counter()
    t=time.perf_counter(); base=generate_mock(norm_cfg,seeds[0],scale=norm_scale,A_true=0,response=False,n_los=100,pixel_noise_power=0)
    save_mock(base,mock_root/"fiducial_reduced.h5"); timing["fiducial_mock_s"]=time.perf_counter()-t
    t=time.perf_counter(); xm=xi_from_model(ForestPower(model="kaiser"),norm_cfg,nk=320,kmax=20); timing["xi_model_s"]=time.perf_counter()-t
    t=time.perf_counter(); xd=xi_from_data(base.sightlines,norm_cfg); timing["xi_data_s"]=time.perf_counter()-t
    timing["benchmark"]=benchmark(base.sightlines,1.0,xd,norm_cfg,mock_root/"benchmark.json")
    record(rows,1,"tables and benchmark","report only",timing["benchmark"]["pixel_pairs_per_s"],None,"pixel pairs/s, 24 threads")

    # Fixed-geometry physical slopes, same seed so cosmic/noise terms cancel.
    Agrid=[0,.5,1,2]; ahd=[]; ahm=[]; ahnog=[]; ahnom=[]
    for A in Agrid:
        m=generate_mock(norm_cfg,seeds[0],scale=norm_scale,A_true=A,g_on=True,response=False,n_los=100,pixel_noise_power=0)
        ahd.append(fit(cat_for(m,xd,norm_cfg),m.alpha_lya,norm_cfg.g1).A[0])
        ahm.append(fit(cat_for(m,xm,norm_cfg),m.alpha_lya,norm_cfg.g1).A[0])
        mn=generate_mock(norm_cfg,seeds[0],scale=norm_scale,A_true=A,g_on=False,response=False,n_los=100,pixel_noise_power=0)
        ahnog.append(fit(cat_for(mn,xd,norm_cfg),mn.alpha_lya,0).A[0])
        ahnom.append(fit(cat_for(m,xd,norm_cfg),m.alpha_lya,0).A[0])
    sd,sm,snog,snom=(slope(Agrid,v) for v in (ahd,ahm,ahnog,ahnom))
    ratio=sd/sm if sm else np.nan
    record(rows,2,"baseline absorption: data/model normalisation","|ratio-1| <= 0.03",ratio,abs(ratio-1)<=.03)
    record(rows,3,"physical slope, g on + moments","1 +/- 0.05",sd,abs(sd-1)<=.05)
    record(rows,3,"physical slope, g off + moments disabled","1 +/- 0.05",snog,abs(snog-1)<=.05)
    predicted_no_moment=1+cfg.g1*(np.mean(base.sightlines.chi)-cfg.chi_ref)
    record(rows,3,"g on, moments disabled deviation","deviates by predicted amount",snom,
           abs((snom-sd)-(predicted_no_moment-1))<.05,
           f"predicted multiplicative shift {predicted_no_moment-1:+.4f}")

    # Injection bookkeeping.
    inj=injection_test(base.sightlines,xd,base.alpha_lya,[-2,-1,0,1,2],norm_cfg)
    record(rows,5,"paired coordinate injection slope","1 +/- 0.05",inj["paired_slope"],abs(inj["paired_slope"]-1)<=.05)
    ce=np.nanmedian(inj["errors"][:,1])
    record(rows,5,"coordinate injection curl","0 within errors",inj["curl_slope"],abs(inj["curl_slope"])<=ce if np.isfinite(ce) else False)

    # Seed ensemble. Paired fixed-geometry differences isolate the lensing
    # increment; raw A=0 values test the exact mean-field subtraction.
    null=[]; null_cmb=[]; null_slab=[]; lens=[]; stochastic=[]; resp=[]; dep=[]; combined=[]; combined_raw=[]; jk=[]; shapes=[]
    t=time.perf_counter()
    for seed in seeds:
        m0=generate_mock(cfg,seed,scale=scale,A_true=0,response=False)
        m1=generate_mock(cfg,seed,scale=scale,A_true=1,response=False)
        c0=cat_for(m0,xd,cfg); c1=cat_for(m1,xd,cfg)
        r0=fit(c0,m0.alpha_lya,cfg.g1); r1=fit(c1,m1.alpha_lya,cfg.g1)
        null.append(r0.A[0]); lens.append(r1.A[0]-r0.A[0])
        raw0=alpha_from_map(m0,"kappa_CMB"); raw1=alpha_from_map(m1,"kappa_CMB")
        null_cmb.append(fit(c0,raw0,cfg.g1).A[0]); null_slab.append(fit(c0,alpha_from_map(m0,"kappa_slab"),cfg.g1).A[0])
        # Map-level Wiener scalar, estimated from the two fields (no fit data).
        h=np.sum(m0.alpha_lya*raw0)/max(np.sum(raw0*raw0),1e-30)
        rs0=fit(c0,h*raw0,cfg.g1); rs1=fit(c1,h*raw1,cfg.g1); stochastic.append(rs1.A[0]-rs0.A[0])
        mr0=generate_mock(cfg,seed,scale=scale,A_true=0,response=True)
        mr1=generate_mock(cfg,seed,scale=scale,A_true=1,response=True,magnification=True,completeness=True,real_mask=True)
        cr0=cat_for(mr0,xd,cfg); cr1=cat_for(mr1,xd,cfg)
        slab=alpha_from_map(mr0,"kappa_slab"); cmb=alpha_from_map(mr0,"kappa_CMB"); deproj=cmb-slab
        rr=fit(cr0,slab,cfg.g1); resp.append(rr.A[0])
        rd=fit(cr0,deproj,cfg.g1); dep.append(rd.A[0])
        rc0=fit(cr0,mr0.alpha_lya,cfg.g1); rc1=fit(cr1,mr1.alpha_lya,cfg.g1)
        combined.append(rc1.A[0]-rc0.A[0]); combined_raw.append(rc1.A[0]); jk.append(rc1.jk_error[0])
        sh1=shape_test(cr1,Template(mr1.alpha_lya,"truth","truth"),cfg.g1,regions(cr1))
        sh0=shape_test(cr0,Template(mr0.alpha_lya,"truth","truth"),cfg.g1,regions(cr0))
        shapes.append(np.asarray(sh1["A"])-np.asarray(sh0["A"]))
    timing["ensemble_s"]=time.perf_counter()-t
    def meanerr(v):
        v=np.asarray(v,float); return float(v.mean()),float(v.std(ddof=1)/np.sqrt(len(v))) if len(v)>1 else np.nan
    mn,se=meanerr(null); record(rows,4,"unlensed truth-template mean","0 within 2 sigma",mn,abs(mn)<=2*se,f"SEM={se:.4g}; mf/F subtraction included")
    for label,vals in (("filtered kappa_CMB",null_cmb),("matched slab",null_slab)):
        mm,ee=meanerr(vals); record(rows,4,f"unlensed {label} mean","0 within 2 sigma",mm,abs(mm)<=2*ee,f"SEM={ee:.4g}")
    msig,esig=meanerr(stochastic); record(rows,6,"stochastic Wiener-template normalisation","1 within 2 sigma",msig,abs(msig-1)<=2*esig,f"SEM={esig:.4g}")
    mr,er=meanerr(resp); record(rows,7,"response-only predicted amplitude","prediction within errors",mr,False,"required independent pre-fit prediction is not implemented")
    md,ed=meanerr(dep); record(rows,7,"deprojected response null","0 within errors",md,abs(md)<=2*ed,f"SEM={ed:.4g}")
    record(rows,7,"matched-slab response","consistent with prediction",mr,False,"no independent theoretical prediction was computed")
    mc,ec=meanerr(combined); record(rows,8,"combined deprojected recovery","1 within errors",mc,False,f"SEM={ec:.4g}; truth-template proxy, not the required deprojected fit")
    # Switch-off shifts are directly measured relative to the paired no-response ensemble.
    ml,el=meanerr(lens); record(rows,8,"flag/slab switch-off shifts","report only",mc-ml,None,f"all flags on minus all off; off mean={ml:.4g}")
    sha=np.asarray(shapes); shmean=sha.mean(axis=0); sherr=sha.std(axis=0,ddof=1)/np.sqrt(len(sha)) if len(sha)>1 else np.full(6,np.nan)
    ok=np.isfinite(sherr)&(sherr>0)
    if ok.sum()>1:
        from scipy.stats import chi2 as chi2dist
        wm=np.sum(shmean[ok]/sherr[ok]**2)/np.sum(1/sherr[ok]**2); ch=np.sum(((shmean[ok]-wm)/sherr[ok])**2); pmed=float(chi2dist.sf(ch,ok.sum()-1))
    else: pmed=np.nan
    record(rows,9,"six-bin shape consistency","p > 0.01",pmed,pmed>.01,"paired increments across seeds")

    # 100 independent random templates on one catalogue.
    rng=np.random.default_rng(912); sig=np.std(base.alpha_lya,axis=0); randoms=[rng.normal(size=base.alpha_lya.shape)*sig for _ in range(100)]
    cbase=cat_for(base,xd,norm_cfg); re=random_ensemble(cbase,randoms,norm_cfg.g1,regions(cbase)); std=float(np.std(re["A"],ddof=1))
    rjk=float(np.nanmedian(re["jk_error"])); rf=float(np.nanmedian(re["sigma_F"]))
    record(rows,10,"random-template std / jackknife","report only",std/rjk if rjk else np.nan,None)
    record(rows,10,"random-template std / sigma_F","report only",std/rf if rf else np.nan,None)
    scatter=float(np.std(combined_raw,ddof=1)); medjk=float(np.nanmedian(jk)); covratio=scatter/medjk if medjk else np.nan
    record(rows,11,"mock scatter / jackknife","1 +/- 0.3",covratio,abs(covratio-1)<=.3)

    # Required figures.
    fig,ax=plt.subplots(); ax.plot(Agrid,ahd,"o-",label="measured xi"); ax.plot(Agrid,ahm,"s-",label="model xi")
    ax.set(xlabel=r"$A_{true}$",ylabel=r"$\hat A$",title="Fixed-geometry physical normalisation"); ax.legend(); fig.tight_layout(); fig.savefig(figs/"mock_physical_normalisation.pdf"); plt.close(fig)
    fig,ax=plt.subplots(); ax.hist(combined,bins=max(5,int(np.sqrt(len(combined))))); ax.axvline(1,color="k",ls="--"); ax.set(xlabel="paired recovered amplitude",ylabel="mocks"); fig.tight_layout(); fig.savefig(figs/"mock_combined_ensemble.pdf"); plt.close(fig)
    fig,ax=plt.subplots(); ax.errorbar(np.arange(6),shmean,yerr=sherr,fmt="o"); ax.set(xlabel="shape bin",ylabel=r"paired $\Delta\hat A$"); fig.tight_layout(); fig.savefig(figs/"mock_shape.pdf"); plt.close(fig)
    timing["total_s"]=time.perf_counter()-ttotal
    result={"configuration":{"scale":scale,"normalisation_scale":norm_scale,"seeds":list(seeds),"threads":24,"mock_root":str(mock_root),
                              "reduced_patch_pair_regions":True},"timing":timing,"acceptance":rows,
            "raw":{"Agrid":Agrid,"data_A":ahd,"model_A":ahm,"null":null,"lensing_increment":lens,
                   "combined_increment":combined,"shape_values":sha.tolist(),"shape_p":pmed}}
    (out_root/"mock_validation.json").write_text(json.dumps(result,indent=2,default=float)+"\n")
    lines=["# Stage A mock validation","",f"Run: `{len(seeds)}` seeds, mock scale `{scale}`, 24 Numba threads.","",
           "| Sec. 2.8 | Acceptance test | Tolerance | Measured | Result |","|---:|---|---|---:|:---:|"]
    for r in rows:
        status="—" if r["pass"] is None else ("PASS" if r["pass"] else "FAIL")
        val=r["measured"]; val=f"{val:.6g}" if isinstance(val,(float,int)) else str(val)
        lines.append(f'| {r["item"]} | {r["test"]} | {r["tolerance"]} | {val} | {status} |')
    lines += ["","## Wall times",""]+[f"- `{k}`: {v if isinstance(v,dict) else f'{v:.3f} s'}" for k,v in timing.items()]
    lines += ["","## Notes","", "This run uses the reduced patch controlled by `scale`; `/data` was read-only, so large products were written under `/tmp/LyaLenser/mocks`. See `code/pipeline/NOTES.md` for deviations and diagnoses."]
    (out_root/"mock_validation.md").write_text("\n".join(lines)+"\n")
    return result


if __name__=="__main__":
    ap=argparse.ArgumentParser(); ap.add_argument("--scale",type=float,default=.03); ap.add_argument("--seeds",type=int,nargs="*",default=list(range(20)))
    ap.add_argument("--output",default=str(ROOT/"report")); ap.add_argument("--mock-root",default="/tmp/LyaLenser/mocks")
    args=ap.parse_args(); r=run(args.scale,args.seeds,args.output,args.mock_root)
    print(json.dumps({"timing":r["timing"],"passed":sum(x["pass"] is True for x in r["acceptance"]),
                      "failed":sum(x["pass"] is False for x in r["acceptance"])},indent=2,default=float))
