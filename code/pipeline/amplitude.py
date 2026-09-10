"""Template scores, response matrices, mean fields and jackknife errors."""
from __future__ import annotations

from dataclasses import dataclass, field
import numpy as np
try:
    from .pairs import PairCatalogue
except ImportError:
    from pairs import PairCatalogue


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

    @property
    def jk_error(self):
        return np.sqrt(np.maximum(np.diag(self.jk_cov),0))


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
    try: return np.linalg.solve(F,y)
    except np.linalg.LinAlgError: return np.linalg.pinv(F,rcond=1e-12)@y


def amplitude(cat: PairCatalogue,templates:list,g1:float=0.0,regions=None,bins=None):
    if regions is None: regions=np.zeros(len(cat.a),np.int32)
    regions=np.asarray(regions)
    names,regvals,pq,pF,pmf=_partials(cat,templates,g1,regions,bins)
    q=pq.sum(axis=0); F=pF.sum(axis=0); mf=pmf.sum(axis=0)
    A=_solve(F,q-mf); Finv=np.linalg.pinv(F,rcond=1e-12)
    sigma=np.sqrt(np.maximum(np.diag(Finv),0))
    nr=len(regvals); jk=np.zeros((nr,len(names)))
    if nr>1:
        for r in range(nr): jk[r]=_solve(F-pF[r],(q-pq[r])-(mf-pmf[r]))
        avg=jk.mean(axis=0); dif=jk-avg
        cov=(nr-1)/nr*dif.T@dif
    else:
        jk[0]=A; cov=np.full_like(F,np.nan)
    kinds=[getattr(t,"kind","") for t in templates]
    return AmplitudeResult(names,q,F,mf,A,sigma,jk,cov,regvals,pq,pF,pmf,
                           {"g1":g1,"bins":list(np.arange(6) if bins is None else np.atleast_1d(bins)),
                            "has_junk":("junk" in kinds),"has_curl":("curl" in kinds)})


def random_ensemble(cat,list_of_alpha,g1=0.0,regions=None):
    vals=[]; errs=[]; scales=[]
    for a in list_of_alpha:
        r=amplitude(cat,[a],g1,regions)
        vals.append(r.A[0]); errs.append(r.jk_error[0]); scales.append(r.sigma_F[0])
    return {"A":np.asarray(vals),"jk_error":np.asarray(errs),"sigma_F":np.asarray(scales)}


def shape_test(cat,template,g1=0.0,regions=None):
    vals=[]; errs=[]
    for b in range(6):
        r=amplitude(cat,[template],g1,regions,bins=[b]); vals.append(r.A[0]); errs.append(r.jk_error[0])
    vals=np.asarray(vals); errs=np.asarray(errs)
    ok=np.isfinite(errs)&(errs>0)
    if ok.sum()>1:
        mean=np.sum(vals[ok]/errs[ok]**2)/np.sum(1/errs[ok]**2)
        chi2=np.sum(((vals[ok]-mean)/errs[ok])**2)
        from scipy.stats import chi2 as chi2dist
        p=float(chi2dist.sf(chi2,ok.sum()-1))
    else: mean=np.nan; chi2=np.nan; p=np.nan
    return {"A":vals,"error":errs,"weighted_mean":mean,"chi2":chi2,"p_value":p}


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
