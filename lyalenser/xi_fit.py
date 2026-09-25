"""Model-shaped correlation tables for the pair response (iteration 5).

The pair kernel needs xi(r_perp, r_par) and dxi/dr_perp at every pixel pair. Iteration 4 differentiated a
smoothed measurement of xi, whose normalisation depended on the smoothing width (a 4% bias). Here the SHAPE comes
from a model: P_F(k) = b_F^2 (1 + beta_F mu^2)^2 P_lin(k) F_NL(k, mu) is expanded in the three basis spectra
P_lin F_NL mu^{2i}, i = 0, 1, 2, each transformed to (xi_i, dxi_i/dr_perp) exactly as the pixels see it:
  * `provider='grid'` (mocks): discrete-grid covariance with the trilinear pixel window and aliases
    (`grid_covariance.xi_from_mock_grid`), i.e. the generator's own operators;
  * `provider='hankel'` (data): Hankel/cosine transform (`xi_model.xi_from_model`) with a line-of-sight pixel
    window W(k_par) = sinc(k_par dz/2) exp(-k_par^2 sigma_R^2 / 2).
Both are then pushed through the per-forest continuum projection (mean + slope removal) on the actual radial pixel
grid, which commutes with d/dr_perp, giving fine projected tables. The two physical parameters (b_F^2, beta_F) are
fitted by weighted least squares of the coarse-binned model to the measured 1 Mpc/h table of the same sample
(weights = the accumulated pair weights), over a declared fit range; the fitted table (xi and its analytic
derivative) is what the estimator uses. No smoothing, no width.
"""
from __future__ import annotations
from dataclasses import dataclass
import json, sys
import numpy as np
from scipy.optimize import least_squares
from lyalenser.forest_power import ForestPower, FID
from lyalenser.xi_model import XiTable, xi_from_model
from lyalenser.config import Config

BASIS=('mu0','mu2','mu4')


class BasisPower:
    """P_lin(k) F_NL(k, mu) mu^{2i} [W_los(k_par)]^2, with the Kaiser prefactor stripped."""
    def __init__(self,i,model='kaiser',los_pixel=0.,los_resolution=0.,**params):
        self.i=i; self.model=model; self.pf=ForestPower(model=model,**params); self.pf.p['b_F']=1.; self.pf.p['beta_F']=0.
        self.los_pixel=float(los_pixel); self.los_resolution=float(los_resolution)
        # Kaiser-type F_NL depends on k only: tabulate P_lin exp(-(k/kp)^2) once on a log grid (fast on 3D grids).
        self._kg=np.geomspace(1e-4,49,3000); self._pg=self.pf.plin(self._kg)*np.exp(-(self._kg/self.pf.p['kp'])**2)
    def __call__(self,kpar,kperp):
        kpar=np.abs(np.asarray(kpar,float)); kperp=np.asarray(kperp,float)
        k=np.maximum(np.sqrt(kpar**2+kperp**2),1e-6); mu2=(kpar/k)**2
        w=np.ones_like(k)
        if self.los_pixel>0: w=w*np.sinc(kpar*self.los_pixel/(2*np.pi))**2
        if self.los_resolution>0: w=w*np.exp(-(kpar*self.los_resolution)**2)
        base=np.interp(np.clip(k,self._kg[0],self._kg[-1]),self._kg,self._pg) if self.model=='kaiser' else self.pf(kpar,kperp)
        return base*mu2**self.i*w


def basis_tables(cfg:Config,provider:str,grid_shape=None,nk=None,model='kaiser',params=None,los_pixel=0.,los_resolution=0.):
    """Unprojected fine tables (xi, dxi/dr_perp) of the three basis spectra."""
    params=dict(params or {}); out={}
    for i,name in enumerate(BASIS):
        bp=BasisPower(i,model=model,los_pixel=los_pixel,los_resolution=los_resolution,**params)
        if provider=='grid':
            from grid_covariance import xi_from_mock_grid
            out[name]=xi_from_mock_grid(cfg,grid_shape,power=bp)
        elif provider=='hankel':
            out[name]=xi_from_model(bp,cfg,nk=nk or cfg.analytic_nk)
        else: raise ValueError(provider)
    return out


def _continuum_projector(cpix):
    c=np.asarray(cpix,float); u=np.column_stack((np.ones(len(c)),c-c.mean())); return u,u@np.linalg.inv(u.T@u)


def project_fine(table:XiTable,cpix,cfg:Config,subsamples=4):
    """Per-forest continuum projection P C P^T of a stationary table, averaged over pixel pairs at fixed r_par.

    Uniform weights within a forest (exact for the mocks, where the pixel weight is constant along a skewer;
    an approximation for data). Applied identically to xi and to dxi/dr_perp, which commute with the radial
    projection. Returns fine tables on the (r_perp, r_par) grid of ``cfg``.
    """
    from scipy.interpolate import RegularGridInterpolator
    c=np.asarray(cpix,float); n=len(c); u,v=_continuum_projector(c)
    rz=abs(c[:,None]-c[None,:]); step=float(np.median(np.diff(c))); kidx=np.rint(rz/step).astype(int)
    nlag=int(cfg.xi_max/step)+2; counts=np.bincount(kidx.ravel(),minlength=nlag)[:nlag]
    fine=np.arange(0,cfg.xi_max+.5*cfg.xi_step,cfg.xi_step)
    fx=RegularGridInterpolator((table.r_perp,table.r_par),table.xi,bounds_error=False,fill_value=0.)
    fg=RegularGridInterpolator((table.r_perp,table.r_par),table.xi_rp,bounds_error=False,fill_value=0.)
    lag_xi=np.zeros((len(fine),nlag)); lag_g=np.zeros_like(lag_xi)
    def project(M):
        return M-u@(v.T@M)-(M@v)@u.T+u@(v.T@M@v)@u.T
    for ir,rp in enumerate(fine):
        for f,acc in ((fx,lag_xi),(fg,lag_g)):
            M=project(f((np.full_like(rz,rp),rz)))
            acc[ir]=np.bincount(kidx.ravel(),weights=M.ravel(),minlength=nlag)[:nlag]/np.maximum(counts,1)
    lags=np.arange(nlag)*step; keep=counts>0
    xi=np.array([np.interp(fine,lags[keep],row[keep]) for row in lag_xi]); g=np.array([np.interp(fine,lags[keep],row[keep]) for row in lag_g])
    meta=dict(table.meta or {}); meta.update(projection='per-forest mean+slope, uniform weights',pixel_step=step)
    return XiTable(fine,fine,xi,g,meta)


def _weighted_projector(c,w,region=None):
    """u, v for the weighted mean+slope removal on pixel positions c with weights w (picca fits with weights).

    One block per delta region: a quasar that contributes both a Lya-region (A) and a Lyb-region (B) segment had
    TWO independent continuum fits, so the projector removes a mean and a slope from each segment separately.
    """
    c=np.asarray(c,float); w=np.asarray(w,float)
    reg=np.zeros(len(c),np.int8) if region is None else np.asarray(region,np.int8)
    cols=[]
    for r in np.unique(reg):
        m=reg==r
        if m.sum()<2: continue
        o=np.zeros(len(c)); o[m]=1.
        t=np.zeros(len(c)); t[m]=c[m]-np.average(c[m],weights=w[m])
        cols+= [o,t]
    u=np.column_stack(cols)
    return u,(u*w[:,None])@np.linalg.pinv(u.T@(u*w[:,None]))


def project_fine_sample(tables:dict,forests,cfg:Config,n_pairs=400,seed=11,step=None,rp_grid=None,verbose=False):
    """Continuum projection averaged over real forest PAIRS (iteration 8), replacing `project_fine`.

    The measured cell is sum_{ab} sum_{p in a, q in b} w_p w_q d_p d_q / sum w_p w_q, whose expectation is the
    pair-weight average of (P_a C P_b^T)_{pq}, with a DIFFERENT projector on each side. `project_fine` used one
    representative forest on both sides with uniform weights; on DR1 that forest spanned the whole slab (1295
    pixels) while the median forest has 481, and the resulting model over-predicts xi by 5-12 per cent across the
    kernel range. Here the average is estimated by Monte Carlo over sampled forest pairs, accumulating numerator
    and denominator in radial-lag bins exactly as the measurement does, with the real per-pixel weights.

    ``forests`` is a sequence of (chi, w[, in_range[, region]]): the projector uses the WHOLE picca forest
    (continuum fitting saw every pixel), while the pair sums use only the pixels the measurement kept, so a
    redshift cut shortens the pair range but not the projection. ``region`` tags each pixel with its delta region
    (0 = Lya, 1 = Lyb); the projector then has one mean-and-slope block per region, and B x B pixel pairs are
    dropped from the average exactly as the estimator drops them. Returns {name: XiTable} on (rp_grid, fine r_par).
    """
    forests=[(f[0],f[1],(np.asarray(f[2],float) if len(f)>2 else np.ones(len(f[0]))),
              (np.asarray(f[3],np.int8) if len(f)>3 else np.zeros(len(f[0]),np.int8))) for f in forests]
    rng=np.random.default_rng(seed)
    if step is None:
        step=float(np.median(np.diff(np.unique(np.round(np.concatenate([f[0] for f in forests]),4)))))
    fine=np.arange(0,cfg.xi_max+.5*cfg.xi_step,cfg.xi_step)
    rp_grid=fine if rp_grid is None else np.asarray(rp_grid,float)
    nlag=int(cfg.xi_max/step)+2; lags=np.arange(nlag)*step
    names=list(tables)
    from scipy.interpolate import RegularGridInterpolator as _RGI
    interp=[(_RGI((tables[k].r_perp,tables[k].r_par),a,bounds_error=False,fill_value=0.))
            for k in names for a in (tables[k].xi,tables[k].xi_rp)]
    rows=np.array([[f((np.full(nlag,rp),lags)) for rp in rp_grid] for f in interp])
    nsurf=rows.shape[0]
    numr=np.zeros((nsurf,len(rp_grid),nlag)); dend=np.zeros(nlag)
    if len(forests)<2: raise ValueError('need at least two forests to sample pairs')
    ia_all=rng.integers(0,len(forests),size=n_pairs)
    ib_all=(ia_all+1+rng.integers(0,len(forests)-1,size=n_pairs))%len(forests)   # never a forest with itself
    for it,(ia,ib) in enumerate(zip(ia_all,ib_all)):
        ca,wa,ma,ra=forests[ia]; cb,wb,mb,rb=forests[ib]
        ua,va=_weighted_projector(ca,wa,ra); ub,vb=_weighted_projector(cb,wb,rb)
        k=np.rint(np.abs(ca[:,None]-cb[None,:])/step).astype(np.int64)
        kfull=np.minimum(k,nlag-1); keep=(k<nlag); kf=kfull.ravel()
        # pair weight: zero beyond the table, on pixels the measurement does not keep, and on B x B pairs
        W=((wa*ma)[:,None]*(wb*mb)[None,:])*keep*(1.-(ra>0)[:,None]*(rb>0)[None,:])
        dend+=np.bincount(kf,weights=W.ravel(),minlength=nlag)[:nlag]
        for si in range(nsurf):
            for ir in range(len(rp_grid)):
                M=rows[si,ir][kfull]*keep           # the table is zero beyond xi_max, as in project_fine
                M=M-ua@(va.T@M)-(M@vb)@ub.T+ua@((va.T@M@vb)@ub.T)
                numr[si,ir]+=np.bincount(kf,weights=(W*M).ravel(),minlength=nlag)[:nlag]
        if verbose and it%20==0: print(f'  projection pair {it}/{n_pairs}',flush=True)
    prof=np.divide(numr,np.maximum(dend,1e-300),out=np.zeros_like(numr),where=dend>0); ok=dend>0
    out={}
    for i,k in enumerate(names):
        xi=np.array([np.interp(fine,lags[ok],row[ok]) for row in prof[2*i]])
        g =np.array([np.interp(fine,lags[ok],row[ok]) for row in prof[2*i+1]])
        meta=dict(tables[k].meta or {}); meta.update(projection=f'pair-weighted, {n_pairs} real forest pairs, weighted mean+slope',pixel_step=step)
        out[k]=XiTable(rp_grid,fine,xi,g,meta)
    return out


def refine_rp_grid(tables:dict,cfg:Config):
    """Cubic interpolation of coarse-in-r_perp projected tables onto the standard fine grid."""
    from scipy.interpolate import CubicSpline
    fine=np.arange(0,cfg.xi_max+.5*cfg.xi_step,cfg.xi_step); out={}
    for k,t in tables.items():
        xi=CubicSpline(t.r_perp,t.xi,axis=0)(fine); g=CubicSpline(t.r_perp,t.xi_rp,axis=0)(fine)
        out[k]=XiTable(fine,t.r_par,xi,g,dict(t.meta or {},rp_grid='cubic refinement'))
    return out


def coarse_bin(table:XiTable,cfg:Config,subsamples=4):
    """Average a fine table over 1 Mpc/h cells the way the measured num/den table is binned (uniform in r_perp)."""
    n=int(cfg.xi_max); out=np.zeros((n,n))
    from scipy.interpolate import RegularGridInterpolator
    f=RegularGridInterpolator((table.r_perp,table.r_par),table.xi,bounds_error=False,fill_value=0.)
    offs=(np.arange(subsamples)+.5)/subsamples
    for i in range(n):
        for j in range(n):
            rp=i+offs; rz=j+offs; out[i,j]=f((np.repeat(rp,subsamples),np.tile(rz,subsamples))).mean()
    return out


@dataclass
class FittedTable:
    table: XiTable
    params: dict
    def save(self,path,group='xi'):
        self.table.meta=dict(self.table.meta or {}); self.table.meta['fit']=self.params
        self.table.save(path,group)


def fit_model_table(num,den,projected_basis:dict,cfg:Config,coarse_basis:dict=None):
    """Fit (b_F^2, beta_F) of the projected basis to the measured coarse table; return the fitted fine table.

    Weighted least squares with weights = den (the accumulated pair weights, i.e. the inverse-variance weight of
    each cell up to a common factor); the fit range is declared in ``cfg`` (fit_rperp_min <= r_perp < r_perp_max,
    r_par < r_par_max). A free 3-coefficient linear fit is reported as a diagnostic of the model's adequacy.
    """
    raw=np.divide(num,den,out=np.zeros_like(num),where=den>0)
    coarse={k:(coarse_bin(projected_basis[k],cfg) if coarse_basis is None else coarse_basis[k]) for k in BASIS}
    n=raw.shape[0]; rp=np.arange(n)+.5; rz=np.arange(n)+.5
    use=(den>0)&(rp[:,None]>=cfg.fit_rperp_min)&(rp[:,None]<cfg.r_perp_max)&(rz[None,:]<cfg.r_par_max)
    w=np.sqrt(den[use]); y=raw[use]; B=np.column_stack([coarse[k][use] for k in BASIS])
    lin=np.linalg.lstsq(B*w[:,None],y*w,rcond=None)[0]
    def model(p): b2,beta=p; return b2*(B[:,0]+2*beta*B[:,1]+beta*beta*B[:,2])
    p0=np.array([max(lin[0],1e-6),lin[1]/(2*max(lin[0],1e-6))])
    res=least_squares(lambda p:(model(p)-y)*w,p0,x_scale=[max(abs(p0[0]),1e-6),1.])
    b2,beta=res.x; chi2=float(np.sum(res.fun**2)); dof=int(use.sum()-2)
    tab=projected_basis['mu0']
    coef=np.array([b2,2*b2*beta,b2*beta*beta])
    xi=sum(c*projected_basis[k].xi for c,k in zip(coef,BASIS)); g=sum(c*projected_basis[k].xi_rp for c,k in zip(coef,BASIS))
    params={'b_F2':float(b2),'beta_F':float(beta),'chi2':chi2,'dof':dof,'cells':int(use.sum()),
            'linear_coefficients':lin.tolist(),'linear_beta':float(lin[1]/(2*lin[0])) if lin[0]!=0 else None,
            'linear_beta2_consistency':float(lin[2]/lin[0]) if lin[0]!=0 else None,
            'fit_range':{'rperp_min':float(cfg.fit_rperp_min),'rperp_max':float(cfg.r_perp_max),'rpar_max':float(cfg.r_par_max)}}
    table=XiTable(tab.r_perp,tab.r_par,xi,g,{'provider':'model_fit','basis':list(BASIS),'accepted_weight':float(den.sum())},(num,den))
    return FittedTable(table,params)


def save_basis(path,basis:dict):
    for k,t in basis.items(): t.save(path,f'basis/{k}')

def load_basis(path):
    from lyalenser.tables import read_xi
    return {k:read_xi(path,f'basis/{k}') for k in BASIS}
