"""Observed forest correlation tables used by the pair response."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import sys
import numpy as np
from scipy.special import j0, j1
from scipy.ndimage import gaussian_filter
from numba import njit

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

    def __post_init__(self):
        self.r_perp = np.ascontiguousarray(self.r_perp, dtype=np.float64)
        self.r_par = np.ascontiguousarray(self.r_par, dtype=np.float64)
        self.xi = np.ascontiguousarray(self.xi, dtype=np.float64)
        self.xi_rp = np.ascontiguousarray(self.xi_rp, dtype=np.float64)
        expected = (len(self.r_perp), len(self.r_par))
        if self.xi.shape != expected or self.xi_rp.shape != expected:
            raise ValueError(f"xi arrays must have shape {expected}")

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


def xi_from_model(pf: ForestPower, cfg: Config, nk: int = 420,
                  kmin: float = 1e-4, kmax: float = 35.0) -> XiTable:
    """Hankel/cosine transform of ``ForestPower`` on logarithmic k grids."""
    rp = np.arange(0, cfg.xi_max + .5*cfg.xi_step, cfg.xi_step)
    rz = rp.copy()
    kp = np.geomspace(kmin, kmax, nk)
    kz = np.geomspace(kmin, kmax, nk)
    wp = _trap_weights_log(kp) * kp**2/(2*np.pi)
    wz = _trap_weights_log(kz) * kz/np.pi  # 2 * int_0^inf dk/(2pi)
    P = pf(kz[:, None], kp[None, :])
    C = np.cos(np.outer(rz, kz))
    J0 = j0(np.outer(rp, kp))
    J1 = j1(np.outer(rp, kp))
    inner0 = (P * wp[None, :]) @ J0.T
    inner1 = (P * (wp*kp)[None, :]) @ J1.T
    xi = (C * wz[None, :]) @ inner0
    xirp = -((C * wz[None, :]) @ inner1)
    xi = xi.T; xirp = xirp.T
    variance = float(np.sum(P * wz[:, None] * wp[None, :]))
    return XiTable(rp, rz, xi, xirp,
                   {"provider": "model", "nk": nk, "variance_integral": variance,
                    "kmin": kmin, "kmax": kmax})


@njit(cache=True)
def _data_hist_kernel(pix_start, chi, delta, weight, a, b, theta,
                      rpmax, rzmax, step):
    nr = int(np.ceil(rpmax/step))
    num = np.zeros((nr, nr), np.float64)
    den = np.zeros((nr, nr), np.float64)
    for ip in range(a.size):
        aa, bb = a[ip], b[ip]
        q0=pix_start[bb]; qend=pix_start[bb+1]
        for p in range(pix_start[aa], pix_start[aa+1]):
            cp=float(chi[p])
            while q0<qend and float(chi[q0])<cp-rzmax: q0+=1
            q=q0
            while q<qend and float(chi[q])<=cp+rzmax:
                rz = abs(float(chi[p])-float(chi[q]))
                rp = .5*(float(chi[p])+float(chi[q]))*theta[ip]
                if rz < rzmax and rp < rpmax:
                    i, j = int(rp/step), int(rz/step)
                    ww = float(weight[p])*float(weight[q])
                    num[i,j] += ww*float(delta[p])*float(delta[q])
                    den[i,j] += ww
                q+=1
    return num, den


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
                                 cfg.xi_max, 1.0)
    raw = np.divide(num, den, out=np.zeros_like(num), where=den > 0)
    # Normalised convolution prevents empty cells being interpreted as zeros.
    valid = (den > 0).astype(float)
    # A 0.76-bin normalized kernel is the smallest stable width on the
    # predeclared 1 Mpc/h counts: narrower kernels visibly follow empty-cell
    # noise in dxi/dr_perp, while wider kernels bias the remapping response.
    smooth = gaussian_filter(raw*valid, 0.76, mode="nearest")
    norm = gaussian_filter(valid, 0.76, mode="nearest")
    coarse = np.divide(smooth, norm, out=np.zeros_like(smooth), where=norm > 1e-8)
    centers = np.arange(coarse.shape[0]) + .5
    fine = np.arange(0, cfg.xi_max + .5*cfg.xi_step, cfg.xi_step)
    from scipy.interpolate import RectBivariateSpline
    spl = RectBivariateSpline(centers, centers, coarse, kx=3, ky=3, s=0)
    xi = spl(np.clip(fine, centers[0], centers[-1]),
             np.clip(fine, centers[0], centers[-1]))
    xirp = spl(np.clip(fine, centers[0], centers[-1]),
               np.clip(fine, centers[0], centers[-1]), dx=1)
    return XiTable(fine, fine, xi, xirp,
                   {"provider": "data", "coarse_step": 1.0,"smoothing_bins":0.76,
                    "accepted_weight": float(den.sum()), "n_sightline_pairs": len(p[0])})


def direct_variance(pf, nk=700, kmin=1e-4, kmax=35.0):
    kp = np.geomspace(kmin, kmax, nk); kz = np.geomspace(kmin, kmax, nk)
    wp = _trap_weights_log(kp)*kp**2/(2*np.pi)
    wz = _trap_weights_log(kz)*kz/np.pi
    return float(np.sum(pf(kz[:,None], kp[None,:])*wz[:,None]*wp[None,:]))
