"""Independent pixel-covariance modulation propagated through continuum projectors.

No observed forest products or fitted scores enter this prediction. The raw
covariance provider, its radial support, and the measured response table are
saved separately. Projection includes all pixels in each forest, including
those outside the estimator's pair cuts; the final scores apply the production
selection (cfg.r_perp_min <= r_perp <= cfg.r_perp_max, |r_par| <= cfg.r_par_max).
"""
import numpy as np
from numba import njit,prange
from pairs import _interp, _interp_layer, _shape_bin, PairCatalogue

@njit(parallel=True,cache=True)
def _predict(starts,chi,w,slab,mod,a,b,theta,rpgrid,rzgrid,raw,rawgrad,
             grid,xi,grad,cref,rd,rpmin,rpmax,rzmax,chi0,dchi,nchi):
    """The final score sweep applies the production pair selection (rpmin <= r_perp <= rpmax, |r_par| <= rzmax);
    the continuum-projection sums keep the full covariance support."""
    out=np.zeros((len(a),11,6))
    for ip in prange(len(a)):
        ia,ib=starts[a[ip]],starts[b[ip]]
        na,nb=starts[a[ip]+1]-ia,starts[b[ip]+1]-ib
        ua=np.ones((na,2)); ub=np.ones((nb,2)); va=np.zeros((na,2)); vb=np.zeros((nb,2))
        ca=0.; cb=0.; wa=0.; wb=0.
        for p in range(na): ca+=np.float64(chi[ia+p])*np.float64(w[ia+p]); wa+=w[ia+p]
        for q in range(nb): cb+=np.float64(chi[ib+q])*np.float64(w[ib+q]); wb+=w[ib+q]
        ca/=wa; cb/=wb; da=0.; db=0.
        for p in range(na): ua[p,1]=chi[ia+p]-ca; da+=w[ia+p]*ua[p,1]**2
        for q in range(nb): ub[q,1]=chi[ib+q]-cb; db+=w[ib+q]*ub[q,1]**2
        for p in range(na): va[p,0]=w[ia+p]/wa; va[p,1]=w[ia+p]*ua[p,1]/da
        for q in range(nb): vb[q,0]=w[ib+q]/wb; vb[q,1]=w[ib+q]*ub[q,1]/db
        bv=np.zeros((na,2)); vt=np.zeros((2,nb)); middle=np.zeros((2,2))
        # First sweep: low-rank continuum terms of D C D - C.
        q0=0
        for p in range(na):
            cp=float(chi[ia+p]); mp=.5*rd*mod[ia+p]
            while q0<nb and chi[ib+q0]<cp-rzgrid[-1]: q0+=1
            q=q0
            while q<nb and chi[ib+q]<=cp+rzgrid[-1]:
                cq=float(chi[ib+q]); mq=.5*rd*mod[ib+q]
                v,_=_interp(.5*(cp+cq)*theta[ip],abs(cp-cq),rpgrid[0],rpgrid[1]-rpgrid[0],len(rpgrid),rzgrid[0],rzgrid[1]-rzgrid[0],len(rzgrid),raw,rawgrad)
                v*=mp+mq+mp*mq
                for j in range(2):
                    bv[p,j]+=v*vb[q,j]; vt[j,q]+=va[p,j]*v
                q+=1
        for p in range(na):
            for j in range(2):
                for k in range(2): middle[j,k]+=va[p,j]*bv[p,k]
        q0=0
        for p in range(na):
            cp=float(chi[ia+p]); mp=.5*rd*mod[ia+p]
            while q0<nb and chi[ib+q0]<cp-rzmax: q0+=1
            q=q0
            while q<nb and chi[ib+q]<=cp+rzmax:
                cq=float(chi[ib+q]); cm=.5*(cp+cq); dc=cp-cq; rp=cm*theta[ip]; rz=abs(dc)
                if rp<=rpmax and rp>=rpmin and slab[ia+p]>=0 and slab[ib+q]>=0:
                    v,_=_interp(rp,rz,rpgrid[0],rpgrid[1]-rpgrid[0],len(rpgrid),rzgrid[0],rzgrid[1]-rzgrid[0],len(rzgrid),raw,rawgrad)
                    mq=.5*rd*mod[ib+q]; v*=mp+mq+mp*mq
                    for j in range(2):
                        v-=ua[p,j]*vt[j,q]+bv[p,j]*ub[q,j]
                        for k in range(2): v+=ua[p,j]*middle[j,k]*ub[q,k]
                    _,g=_interp_layer(rp,rz,cm,grid[0],grid[1]-grid[0],len(grid),grid[0],grid[1]-grid[0],len(grid),chi0,dchi,nchi,xi,grad)
                    val=w[ia+p]*w[ib+q]*v*cm*g; binid=_shape_bin(rp,rz,rpmax,rzmax)
                    out[ip,0,binid]+=val; out[ip,1,binid]+=val*(cm-cref); out[ip,2,binid]+=val*dc*.5
                q+=1
    return out


def prediction_catalogue(sl,cat,pixel_long,raw_xi,measured_xi,cfg):
    if len(pixel_long)!=len(sl.chi): raise ValueError('one modulation value per forest pixel required')
    vals=_predict(sl.pix_start,sl.chi.astype(np.float64),sl.w.astype(np.float64),sl.slab,np.asarray(pixel_long),cat.a,cat.b,cat.theta,
                  raw_xi.r_perp,raw_xi.r_par,raw_xi.xi.ravel(),raw_xi.xi_rp.ravel(),
                  measured_xi.r_perp,measured_xi.xi.ravel(),measured_xi.xi_rp.ravel(),cfg.chi_ref,cfg.response_delta,
                  float(getattr(cfg,'r_perp_min',0.0)),float(cfg.r_perp_max),float(cfg.r_par_max),
                  *(measured_xi.layers() if hasattr(measured_xi,'layers') else (0.,1.,1)))
    # Same F, zero baseline: score is the independently calculated covariance perturbation.
    vals[:,3:8]=cat.accum[:,3:8]
    return PairCatalogue(cat.a,cat.b,cat.thx,cat.thy,cat.theta,vals,cat.npair,
                         {'prediction':'P_a [D_a C_ab D_b - C_ab] P_b^T','raw_support':raw_xi.r_par[-1],
                          'r_perp_min':float(getattr(cfg,'r_perp_min',0.0)),'r_perp_max':float(cfg.r_perp_max)})
