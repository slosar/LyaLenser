"""Template scores, response matrices, mean fields and jackknife errors."""
from __future__ import annotations

from dataclasses import dataclass, field
import numpy as np
from lyalenser.pairs import PairCatalogue


@dataclass
class AmplitudeResult:
    names: list
    q: np.ndarray
    F: np.ndarray
    mf: np.ndarray
    A: np.ndarray
    sigma_F: np.ndarray
    jk_samples: np.ndarray
    jk_cov: np.ndarray
    regions: np.ndarray
    partial_q: np.ndarray
    partial_F: np.ndarray
    partial_mf: np.ndarray
    attrs: dict = field(default_factory=dict)

    def save(self,path,group="fit"):
        import h5py,json
        with h5py.File(path,"a") as f:
            if group in f: del f[group]
            g=f.create_group(group)
            for k in ("q","F","mf","A","sigma_F","jk_samples","jk_cov","regions","partial_q","partial_F","partial_mf"):
                g[k]=getattr(self,k)
            g.attrs["names"]=json.dumps(self.names)
            g.attrs["config"]=json.dumps(self.attrs,default=lambda x: x.item())

    @property
    def jk_error(self):
        return np.sqrt(np.maximum(np.diag(self.jk_cov),0))


def n_science(result,default=3):
    """Number of leading science bands in a fit result; the curl partners follow them."""
    try: return int(result.attrs["n_science"])
    except Exception: return default


def curl_amplitude(result):
    """Mean curl amplitude and its jackknife error, for any number of bands."""
    n=n_science(result)
    return float(np.mean(result.A[n:2*n])),float(np.sqrt(np.mean(result.jk_error[n:2*n]**2)))


def _alpha_name(t,i):
    return (np.asarray(t.alpha,float),getattr(t,"name",f"template{i}")) if hasattr(t,"alpha") else (np.asarray(t,float),f"template{i}")


def pair_scalars(cat,templates):
    ds=[]; ss=[]; names=[]
    for i,t in enumerate(templates):
        a,name=_alpha_name(t,i)
        da=a[cat.a]-a[cat.b]; sa=a[cat.a]+a[cat.b]
        ds.append(cat.thx*da[:,0]+cat.thy*da[:,1])
        ss.append(cat.thx*sa[:,0]+cat.thy*sa[:,1]); names.append(name)
    return np.asarray(ds),np.asarray(ss),names


def _partials(cat,templates,g1,regions,bins=None):
    d,s,names=pair_scalars(cat,templates); nt=len(names)
    bins=np.arange(6) if bins is None else np.atleast_1d(bins)
    x=cat.accum[:,:,bins].sum(axis=2)
    v=x[:,0]+g1*x[:,1]; vc=g1*x[:,2]
    mm=x[:,3]+2*g1*x[:,4]+g1*g1*x[:,5]
    mc=g1*x[:,6]; mcc=g1*g1*x[:,7]
    beta=x[:,8]+g1*x[:,9]; betac=g1*x[:,10]
    qpair=d*(v[None,:])+s*(vc[None,:])
    bpair=d*(beta[None,:])+s*(betac[None,:])
    fpair=np.empty((nt,nt,len(cat.a)))
    for i in range(nt):
        for j in range(nt):
            fpair[i,j]=mm*d[i]*d[j]+.5*mc*(d[i]*s[j]+d[j]*s[i])+mcc*s[i]*s[j]
    regvals=np.unique(regions); nr=len(regvals)
    pq=np.zeros((nr,nt)); pb=np.zeros((nr,nt)); pf=np.zeros((nr,nt,nt))
    inv=np.searchsorted(regvals,regions)
    for k in range(nr):
        m=inv==k; pq[k]=qpair[:,m].sum(axis=1); pb[k]=bpair[:,m].sum(axis=1); pf[k]=fpair[:,:,m].sum(axis=2)
    return names,regvals,pq,pf,pb


def _solve(F,y):
    if not np.all(np.isfinite(F)) or np.linalg.matrix_rank(F) != len(F):
        raise ValueError("singular response matrix")
    return np.linalg.solve(F,y)


def _fit(cat: PairCatalogue,templates:list,g1:float=0.0,regions=None,bins=None):
    if not any(getattr(t,"kind","")=="junk" for t in templates):
        raise ValueError("amplitude fit requires a real junk-band template")
    bands=[t for t in templates if getattr(t,"Lmax",0)>0 and t.kind!="junk"]
    if bands:
        # every band present must carry both a science-like and a curl component (iteration 9: the band list is
        # no longer hard-coded, so the invariant is derived from the templates themselves)
        science_kinds={"signal","truth","injection","response","random"}
        wanted={(t.Lmin,t.Lmax) for t in bands}
        for lo,hi in sorted(wanted):
            for kind in ("signal","curl"):
                if not any(t.Lmin==lo and t.Lmax==hi and (t.kind=="curl" if kind=="curl" else t.kind in science_kinds) for t in bands):
                    raise ValueError(f"missing required {kind} band {lo}-{hi}")
    elif not any(getattr(t,"kind","")=="curl" for t in templates):
        raise ValueError("missing required curl component")
    if regions is None: regions=np.zeros(len(cat.a),np.int32)
    regions=np.asarray(regions)
    names,regvals,pq,pF,pmf=_partials(cat,templates,g1,regions,bins)
    q=pq.sum(axis=0); F=pF.sum(axis=0); mf=pmf.sum(axis=0)
    A=_solve(F,q-mf); Finv=np.linalg.inv(F)
    sigma=np.sqrt(np.maximum(np.diag(Finv),0))
    nr=len(regvals); jk=np.zeros((nr,len(names)))
    if nr>1:
        for r in range(nr):
            try: jk[r]=_solve(F-pF[r],(q-pq[r])-(mf-pmf[r]))
            except ValueError: jk[r]=np.nan
        avg=jk.mean(axis=0); dif=jk-avg
        cov=(nr-1)/nr*dif.T@dif
    else:
        jk[0]=A; cov=np.full_like(F,np.nan)
    kinds=[getattr(t,"kind","") for t in templates]
    nsci=sum(1 for k in kinds if k in {"signal","truth","injection","response","random"})
    return AmplitudeResult(names,q,F,mf,A,sigma,jk,cov,regvals,pq,pF,pmf,
                           {"g1":g1,"bins":list(np.arange(6) if bins is None else np.atleast_1d(bins)),
                            "has_junk":("junk" in kinds),"has_curl":("curl" in kinds),"n_science":int(nsci)})


def random_ensemble(cat,list_of_alpha,g1=0.0,regions=None):
    vals=[]; errs=[]; scales=[]
    for template_set in list_of_alpha:
        if not isinstance(template_set,(list,tuple)):
            raise ValueError("each random realization must include its science/curl/junk template set")
        r=amplitude(cat,list(template_set),g1,regions)
        vals.append(r.A[0]); errs.append(r.jk_error[0]); scales.append(r.sigma_F[0])
    return {"A":np.asarray(vals),"jk_error":np.asarray(errs),"sigma_F":np.asarray(scales)}


def shape_test(cat,templates,g1=0.0,regions=None,target=0):
    if not isinstance(templates,(list,tuple)): templates=[templates]
    vals=[]; errs=[]
    for b in range(6):
        r=amplitude(cat,list(templates),g1,regions,bins=[b]); vals.append(r.A[target]); errs.append(r.jk_error[target])
    vals=np.asarray(vals); errs=np.asarray(errs)
    ok=np.isfinite(errs)&(errs>0)
    if ok.sum()>1:
        mean=np.sum(vals[ok]/errs[ok]**2)/np.sum(1/errs[ok]**2)
        chi2=np.sum(((vals[ok]-mean)/errs[ok])**2)
        from scipy.stats import chi2 as chi2dist
        p=float(chi2dist.sf(chi2,ok.sum()-1))
    else: mean=np.nan; chi2=np.nan; p=np.nan
    return {"A":vals,"error":errs,"weighted_mean":mean,"chi2":chi2,"p_value":p}


def catalogue_modulation_scores(cat,templates,modulation,g1=0.0,bins=None):
    """Expected estimator scores per unit delta_F amplitude modulation.

    ``modulation`` is the long-density value at each sightline.  This uses the
    stored xi*G accumulators, not the measured fitted scores, so it is an
    independent pre-fit catalogue response calculation.
    """
    modulation=np.asarray(modulation,float)
    if len(modulation)<=max(np.max(cat.a,initial=-1),np.max(cat.b,initial=-1)):
        raise ValueError("modulation must contain one value per sightline")
    d,s,_=pair_scalars(cat,templates)
    use=np.arange(6) if bins is None else np.atleast_1d(bins)
    x=cat.accum[:,:,use].sum(axis=2)
    beta=x[:,8]+g1*x[:,9]; betac=g1*x[:,10]
    pair_mod=modulation[cat.a]+modulation[cat.b]
    return np.sum((d*beta[None,:]+s*betac[None,:])*pair_mod[None,:],axis=1)


def independent_response_prediction(cat,templates,density_template,modulation,
                                    cross_ratios,response_delta,g1=0.0,bins=None):
    """Legacy spectral toy retained for algebra regression, NOT a physical prediction.

    Production uses response.prediction_catalogue with per-pixel modulation.

    The density-template modulation score calibrates the catalogue estimator;
    ``cross_ratios[b]`` is C_L^{delta,target_b}/C_L^{delta,delta} in the
    corresponding band (zero for curl).  Nothing from the observed fit score
    enters this prediction.
    """
    if not any(getattr(t,"kind","")=="junk" for t in templates):
        raise ValueError("response prediction requires the same junk-inclusive fit basis")
    density_templates=list(density_template) if isinstance(density_template,(list,tuple)) else [density_template]*len(templates)
    if len(density_templates)!=len(templates):
        raise ValueError("density_template list must match the fitted template basis")
    q_density=catalogue_modulation_scores(cat,density_templates,modulation,g1,bins)
    ratios=np.asarray(cross_ratios,float)
    if ratios.shape != (len(templates),):
        raise ValueError("cross_ratios must have one entry per fitted template")
    dummy_regions=np.zeros(len(cat.a),np.int32)
    _,_,_,pF,_=_partials(cat,templates,g1,dummy_regions,bins)
    F=pF.sum(axis=0)
    predicted_scores=.5*float(response_delta)*q_density*ratios
    predicted_A=_solve(F,predicted_scores)
    fdd=np.empty(len(density_templates))
    for i,td in enumerate(density_templates):
        _,_,_,Fd,_=_partials(cat,[td],g1,dummy_regions,bins); fdd[i]=float(Fd.sum())
    return {"A":predicted_A,"scores":predicted_scores,
            "density_modulation_score":q_density,
            "catalogue_response":np.divide(q_density,fdd,out=np.full_like(q_density,np.nan),where=fdd!=0),
            "cross_ratios":ratios,"F":F}


def compress_score_per_sightline(cat,template,g1=0.0):
    """Return U_a and its exact dot(alpha) score (including mean field)."""
    alpha,_=_alpha_name(template,0); n=len(alpha)
    x=cat.accum.sum(axis=2)
    direction=np.column_stack((cat.thx,cat.thy))
    h=(x[:,0]+g1*x[:,1]-x[:,8]-g1*x[:,9])[:,None]*direction
    k=(g1*(x[:,2]-x[:,10]))[:,None]*direction
    U=np.zeros((n,2)); np.add.at(U,cat.a,h+k); np.add.at(U,cat.b,-h+k)
    score=float(np.sum(U*alpha))
    return U,score


def amplitude(cat,templates,g1=0.,regions=None,bins=None):
    """Public production fit: complete science/curl/junk basis is mandatory."""
    if not any(getattr(t,"kind","")=="junk" for t in templates):
        raise ValueError("amplitude fit requires a real junk-band template")
    # the band list is whatever the templates carry (iteration 13: alternative band sets); _fit checks that every
    # band present has both a science-like and a curl component
    return _fit(cat,templates,g1,regions,bins)


def common_science(result, nscience=None):
    """Collapse the independently fitted science bands of an `AmplitudeResult` to one common amplitude, with the
    nuisance components (curl partners, junk) marginalised: A, its jackknife error, the Fisher error and the
    jackknife samples."""
    if nscience is None: nscience = n_science(result)
    diag = np.diag(result.F)
    active = list(range(nscience)); nuisance = list(range(nscience, len(diag)))
    M = np.zeros((len(diag), 1 + len(nuisance))); M[active, 0] = 1
    for j, i in enumerate(nuisance, 1): M[i, j] = 1
    Fr = M.T @ result.F @ M; yr = M.T @ (result.q - result.mf); pars = np.linalg.solve(Fr, yr)
    jk = []
    for k in range(len(result.regions)):
        Fk = result.F - result.partial_F[k]
        yk = (result.q - result.partial_q[k]) - (result.mf - result.partial_mf[k])
        jk.append(np.linalg.solve(M.T @ Fk @ M, M.T @ yk)[0])
    jk = np.asarray(jk); nr = len(jk)
    jke = float(np.sqrt((nr - 1) / nr * np.sum((jk - jk.mean()) ** 2))) if nr > 1 else np.nan
    return {"A": float(pars[0]), "jk_error": jke,
            "sigma_F": float(np.sqrt(max(np.linalg.inv(Fr)[0, 0], 0))),
            "jk": jk, "active": active}
