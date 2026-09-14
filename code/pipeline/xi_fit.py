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
from pathlib import Path
import json, sys
import numpy as np
from scipy.optimize import least_squares
HERE=Path(__file__).resolve().parent
for p in (HERE,HERE.parent):
    if str(p) not in sys.path: sys.path.insert(0,str(p))
from forest_power import ForestPower, FID
from xi_model import XiTable, xi_from_model
from config import Config

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
    from campaign4 import read_xi
    return {k:read_xi(path,f'basis/{k}') for k in BASIS}
