"""Low-redshift tracers of the lensing mass: the tracer table, the Limber spectra of a redshift slice, the
multipole range of the bias fit and the optimal combination of slice amplitudes.

A density tracer at z < 1.6 shares no modes with the forest at z = 1.96-3, so its correlation with the forest
pair products is lensing only (plus magnification of the sightline quasars, handled by the mean field, and shared
sky systematics). The per-slice HEALPix templates themselves are built by scripts/lowz_catalogues.py:
kernel-weighted tracer maps (forest-source lensing kernel W(chi; chi_ref), per-object 1/(b nbar)), the bias b of
every tracer from its own NaMaster auto-spectrum, the Wiener combination of the tracers of a slice into one
estimate of the slice's kappa_lya contribution, and the sum over slices as the combined estimate (slices are
independent in the Limber approximation; the class-wise, masked harmonic filtering only approximates the
conditional expectation of kappa_lya given all tracers, with normalisation checked against CMB lensing).
"""
from __future__ import annotations
from dataclasses import dataclass, asdict
from functools import lru_cache
import json
import numpy as np
from lyalenser.cosmo import chi as chi_of_z, z_of_chi
from lyalenser.lensing import kernel, Z_CMB
from lyalenser.lensing import DEG2

ANNULUS=40.
KMAX_BIAS=0.2        # h/Mpc: the bias is fitted on large scales only, ell <= KMAX_BIAS * chi(z_mid) of the slice (user, 2026-09-15)
LMIN_BIAS=40.


@dataclass(frozen=True)
class Tracer:
    name: str
    zmin: float
    zmax: float
    bias: float          # nominal bias (the pipeline measures b from the auto-spectrum)
    nbar_deg2: float     # nominal surface density in the slice
    smoothing_mpc: float = 3.0   # transverse smoothing of the Stage-A mock tracers; the data path uses none

    @property
    def chi_mid(self): return float(chi_of_z(.5*(self.zmin+self.zmax)))
    def filter(self,ell):
        sig=self.smoothing_mpc/self.chi_mid
        return np.exp(-.5*(np.asarray(ell,float)*sig)**2)

    @property
    def label(self): return f'{self.name}_{self.zmin:g}_{self.zmax:g}'
    @property
    def chi_range(self): return (float(chi_of_z(self.zmin)),float(chi_of_z(self.zmax)))


# Nominal DESI tracer table (bias, surface density per deg^2); scripts/lowz_catalogues.py adds BGS and BOSS and
# measures every bias. b_QSO(z) from the DESI fit 0.278((1+z)^2-6.565)+2.393.
TRACERS=(Tracer('LRG',.4,.6,1.9,200.),Tracer('LRG',.6,.8,2.1,250.),Tracer('LRG',.8,1.1,2.3,150.),
         Tracer('ELG',.8,1.1,1.3,350.),Tracer('ELG',1.1,1.6,1.4,500.),
         Tracer('QSO',.8,1.1,1.7,40.),Tracer('QSO',1.1,1.6,2.1,60.),Tracer('QSO',1.6,1.75,2.5,25.))


def slices_of(tracers=TRACERS):
    """Sorted redshift slices and the tracers in each; every tracer must cover exactly one slice."""
    bounds=sorted({(t.zmin,t.zmax) for t in tracers})
    for (a,b),(c,d) in zip(bounds[:-1],bounds[1:]):
        if c<b: raise ValueError(f'overlapping slices {(a,b)} and {(c,d)}')
    return [{'zmin':a,'zmax':b,'tracers':[t for t in tracers if (t.zmin,t.zmax)==(a,b)]} for a,b in bounds]


def _ell_grid(lmax):
    return np.unique(np.r_[2.,np.linspace(2,max(2.,lmax),640)])


@lru_cache(maxsize=32)
def slice_spectra(zmin,zmax,chi_source,lmax):
    """Limber 3 x 3 covariance per L of (slice-mean density d, kappa_lya part l, kappa_CMB part c) over [zmin, zmax]."""
    from lyalenser.lensing import limber
    c1,c2=float(chi_of_z(zmin)),float(chi_of_z(zmax)); ccmb=float(chi_of_z(Z_CMB))
    wd=lambda c: np.where((c>=c1)&(c<=c2),1./(c2-c1),0.)
    wl=lambda c: kernel(c,chi_source); wc=lambda c: kernel(c,ccmb)
    L=_ell_grid(lmax); W=(wd,wl,wc); C=np.zeros((len(L),3,3))
    for i in range(3):
        for j in range(i,3):
            C[:,i,j]=C[:,j,i]=limber(L,W[i],W[j],c1,c2,nchi=400,to_recombination=False)
    return L,C


# ----------------------------------------------------------------------------------------------------------------


def bias_band(zmin,zmax,kmax=KMAX_BIAS,lmin=LMIN_BIAS):
    """Multipole range of the bias fit: lmin <= ell <= kmax chi(z_mid) (Limber k = (ell + 1/2)/chi)."""
    cmid=float(chi_of_z(.5*(zmin+zmax)))
    return (float(lmin),float(kmax*cmid-.5))


# ----------------------------------------------------------------------------------------------------------------
# Amplitudes: per slice, combined, and the jackknife combination of the slice amplitudes.


def optimal_combination(A,jk):
    """A: slice amplitudes [k]; jk: jackknife samples [k, nregions]. C = (n-1)/n sum (jk - mean)(jk - mean)^T,
    Hartlap-corrected inverse; A_opt = (1^T C^-1 A)/(1^T C^-1 1)."""
    A=np.asarray(A,float); jk=np.asarray(jk,float); k,nr=jk.shape
    d=jk-jk.mean(axis=1,keepdims=True); C=(nr-1)/nr*d@d.T
    if nr<=k+2: return {'A':float('nan'),'error':float('nan'),'covariance':C.tolist(),'note':'too few regions'}
    hart=(nr-k-2)/(nr-1); Ci=hart*np.linalg.pinv(C); one=np.ones(k)
    var=1/float(one@Ci@one); Aopt=float(one@Ci@A*var)
    naive_w=1/np.maximum(np.diag(C),1e-30); naive=float(np.sum(naive_w*A)/np.sum(naive_w))
    corr=C/np.sqrt(np.outer(np.diag(C),np.diag(C)))
    return {'A':Aopt,'error':float(np.sqrt(var)),'naive_A':naive,'naive_error':float(1/np.sqrt(np.sum(naive_w))),
            'covariance':C.tolist(),'correlation':corr.tolist(),'hartlap':hart,'nregions':nr}
