"""Observed forest correlation tables used by the pair response."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import sys
import numpy as np
from scipy.special import j0, j1
from scipy.ndimage import gaussian_filter
from numba import njit, prange

HERE = Path(__file__).resolve().parent
CODE = HERE.parent
if str(CODE) not in sys.path:
    sys.path.insert(0, str(CODE))
from forest_power import ForestPower
try:
    from .config import Config, SightlineSet
except ImportError:
    from config import Config, SightlineSet


@dataclass
class XiTable:
    r_perp: np.ndarray
    r_par: np.ndarray
    xi: np.ndarray
    xi_rp: np.ndarray
    meta: dict | None = None
    counts: tuple | None = None

    def __post_init__(self):
        self.r_perp = np.ascontiguousarray(self.r_perp, dtype=np.float64)
        self.r_par = np.ascontiguousarray(self.r_par, dtype=np.float64)
        self.xi = np.ascontiguousarray(self.xi, dtype=np.float64)
        self.xi_rp = np.ascontiguousarray(self.xi_rp, dtype=np.float64)
        expected = (len(self.r_perp), len(self.r_par))
        if self.xi.shape != expected or self.xi_rp.shape != expected:
            raise ValueError(f"xi arrays must have shape {expected}")

    def save(self,path,group="xi"):
        import h5py,json
        with h5py.File(path,"a") as f:
            if group in f: del f[group]
            g=f.create_group(group)
            for k in ("r_perp","r_par","xi","xi_rp"): g[k]=getattr(self,k)
            g.attrs["meta"]=json.dumps(self.meta)
            if self.counts is not None:
                g["coarse_num"]=self.counts[0]; g["coarse_den"]=self.counts[1]

    def interp(self, rp, rz):
        return bilinear(rp, rz, self.r_perp, self.r_par,
                        self.xi.ravel(), self.xi_rp.ravel())


@njit(cache=True)
def bilinear(rp, rz, rp_grid, rz_grid, xi_flat, xirp_flat):
    """Numba-callable interpolation returning (xi, dxi/dr_perp)."""
    nrz = rz_grid.size
    if rp < rp_grid[0] or rz < rz_grid[0] or rp > rp_grid[-1] or rz > rz_grid[-1]:
        return 0.0, 0.0
    drp = rp_grid[1] - rp_grid[0]
    drz = rz_grid[1] - rz_grid[0]
    i = min(int((rp-rp_grid[0])/drp), rp_grid.size-2)
    j = min(int((rz-rz_grid[0])/drz), rz_grid.size-2)
    tx = (rp-rp_grid[i])/drp
    ty = (rz-rz_grid[j])/drz
    k00 = i*nrz+j; k10 = (i+1)*nrz+j
    def one(a):
        return ((1-tx)*(1-ty)*a[k00] + tx*(1-ty)*a[k10] +
                (1-tx)*ty*a[k00+1] + tx*ty*a[k10+1])
    return one(xi_flat), one(xirp_flat)


def _trap_weights_log(k):
    x = np.log(k)
    w = np.empty_like(x)
    w[0] = .5*(x[1]-x[0]); w[-1] = .5*(x[-1]-x[-2])
    w[1:-1] = .5*(x[2:]-x[:-2])
    return w


def xi_from_model(pf: ForestPower, cfg: Config, nk: int = 6400,
                  kmin: float = 1e-4, kmax: float = 35.0) -> XiTable:
    """Hankel/cosine transform of ``ForestPower`` on logarithmic k grids."""
    rp = np.arange(0, cfg.xi_max + .5*cfg.xi_step, cfg.xi_step)
    rz = rp.copy()
    kp = np.geomspace(kmin, kmax, nk)
    kz = np.geomspace(kmin, kmax, nk)
    wp = _trap_weights_log(kp) * kp**2/(2*np.pi)
    wz = _trap_weights_log(kz) * kz/np.pi  # 2 * int_0^inf dk/(2pi)
    J0 = j0(np.outer(rp, kp)); J1 = j1(np.outer(rp, kp))
    xi=np.zeros((len(rz),len(rp))); xirp=np.zeros_like(xi); variance=0.
    for i in range(0,nk,128):
        kzpart=kz[i:i+128]; wzpart=wz[i:i+128]
        P=pf(kzpart[:,None],kp[None,:])
        C=np.cos(np.outer(rz,kzpart))*wzpart
        xi+=C@((P*wp)@J0.T)
        xirp-=C@((P*(wp*kp))@J1.T)
        variance+=float(np.sum(P*wzpart[:,None]*wp))
    xi=xi.T; xirp=xirp.T
    return XiTable(rp, rz, xi, xirp,
                   {"provider": "model", "nk": nk, "variance_integral": variance,
                    "kmin": kmin, "kmax": kmax})


@njit(parallel=True,cache=True)
def _data_hist_kernel(pix_start, chi, delta, weight, a, b, theta,
                      rpmax, rzmax, step, slab, slab_edges, slab_index):
    nr = int(np.ceil(rpmax/step))
    nums = np.zeros((24,nr,nr),np.float64)
    dens = np.zeros((24,nr,nr),np.float64)
    for chunk in prange(24):
        for ip in range(chunk*a.size//24,(chunk+1)*a.size//24):
            aa, bb = a[ip], b[ip]
            q0=pix_start[bb]; qend=pix_start[bb+1]
            for p in range(pix_start[aa], pix_start[aa+1]):
                cp=float(chi[p])
                while q0<qend and float(chi[q0])<cp-rzmax: q0+=1
                q=q0
                while q<qend and float(chi[q])<=cp+rzmax:
                    rz = abs(float(chi[p])-float(chi[q]))
                    rp = .5*(float(chi[p])+float(chi[q]))*theta[ip]
                    cm=.5*(float(chi[p])+float(chi[q]))
                    selected=slab[p]>=0 and slab[q]>=0
                    if slab_index>=0: selected=selected and slab_edges[slab_index,0]<=cm<slab_edges[slab_index,1]
                    if rz < rzmax and rp < rpmax and selected:
                        i, j = int(rp/step), int(rz/step)
                        ww = float(weight[p])*float(weight[q])
                        nums[chunk,i,j] += ww*float(delta[p])*float(delta[q])
                        dens[chunk,i,j] += ww
                    q+=1
    return nums.sum(axis=0),dens.sum(axis=0)


def xi_from_data(sl: SightlineSet, cfg: Config) -> XiTable:
    """Measure, smooth and differentiate the observed cross-sightline xi."""
    try:
        from .pairs import find_pairs
    except ImportError:
        from pairs import find_pairs
    chi_min = max(1.0, float(np.min(sl.chi)))
    p = find_pairs(sl, cfg.xi_max/chi_min)
    num, den = _data_hist_kernel(sl.pix_start, sl.chi, sl.delta, sl.w,
                                 p[0], p[1], p[4], cfg.xi_max,
                                 cfg.xi_max, 1.0,sl.slab,
                                 np.asarray([[float(__import__("cosmo").chi(z)) for z in b] for b in cfg.slabs]),cfg.slab_index)
    return xi_from_counts(num,den,cfg,len(p[0]))


def xi_from_counts(num,den,cfg,pair_count=0):
    """Smooth saved 1-Mpc/h counts, allowing reproducible development scans."""
    raw = np.divide(num, den, out=np.zeros_like(num), where=den > 0)
    # Normalised convolution prevents empty cells being interpreted as zeros.
    valid = (den > 0).astype(float)
    # Width is chosen on seeds 100-104, then frozen for every validation sample.
    smooth = gaussian_filter(raw*valid, cfg.xi_smoothing, mode="reflect")
    norm = gaussian_filter(valid, cfg.xi_smoothing, mode="reflect")
    coarse = np.divide(smooth, norm, out=np.zeros_like(smooth), where=norm > 1e-8)
    centers = np.arange(coarse.shape[0]) + .5
    fine = np.arange(0, cfg.xi_max + .5*cfg.xi_step, cfg.xi_step)
    from scipy.interpolate import RectBivariateSpline
    # Correlation is even in both separation coordinates. Extend bin centers
    # across zero before differentiating; clipping all r<0.5 to the first
    # center produced a spurious nonzero gradient at r_perp=0.
    signed=np.r_[-centers[::-1],centers]
    even=np.block([[coarse[::-1,::-1],coarse[::-1,:]],
                   [coarse[:,::-1],coarse]])
    spl = RectBivariateSpline(signed, signed, even, kx=3, ky=3, s=0)
    evaluation=np.minimum(fine,centers[-1])
    xi = spl(evaluation,evaluation)
    xirp = spl(evaluation,evaluation,dx=1)
    return XiTable(fine, fine, xi, xirp,
                   {"provider": "data", "coarse_step": 1.0,"smoothing_bins":cfg.xi_smoothing,"boundary":"even reflection about zero",
                    "accepted_weight": float(den.sum()), "n_sightline_pairs": pair_count},(num,den))


def direct_variance(pf, nk=700, kmin=1e-4, kmax=35.0):
    kp = np.geomspace(kmin, kmax, nk); kz = np.geomspace(kmin, kmax, nk)
    wp = _trap_weights_log(kp)*kp**2/(2*np.pi)
    wz = _trap_weights_log(kz)*kz/np.pi
    return float(np.sum(pf(kz[:,None], kp[None,:])*wz[:,None]*wp[None,:]))
