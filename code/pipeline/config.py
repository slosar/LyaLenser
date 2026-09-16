"""Configuration and lightweight data containers for the Stage-A pipeline."""
from __future__ import annotations

from dataclasses import dataclass, field, replace
from functools import lru_cache
from pathlib import Path
import sys
import numpy as np
from scipy.integrate import trapezoid

CODE = Path(__file__).resolve().parents[1]
if str(CODE) not in sys.path:
    sys.path.insert(0, str(CODE))
from cosmo import chi as chi_of_z


@lru_cache(maxsize=1)
def kernel_product_g1():
    """Linear source-distance coefficient from the kl x kCMB lens kernel.

    The effective lens distance is averaged with the same kernel product and
    matter-power weight that enters the three-tracer Limber cross spectrum,
    averaged over the science range 40 <= L <= 300.  With that lens
    distribution fixed, W(chi_s)=<1-chi_l/chi_s> and g1=W'/W at chi_ref.
    """
    from cosmo import z_of_chi, linear_pk_interp
    from cross_spectrum import kernel, Z_CMB
    cref = float(chi_of_z(2.4))
    ccmb = float(chi_of_z(Z_CMB))
    chis = np.linspace(1.0, cref * (1.0 - 1e-5), 1200)
    zs = z_of_chi(chis)
    pk = linear_pk_interp(zmax=6.0, kmax=200.0, nonlinear=True)
    weight = np.zeros_like(chis)
    for ell in (40.0, 70.0, 100.0, 150.0, 200.0, 250.0, 300.0):
        kval = (ell + 0.5) / chis
        weight += (2.0 * ell + 1.0) * kernel(chis, cref) * kernel(chis, ccmb) \
                  * pk.P(zs, kval, grid=False) / chis**2
    norm = trapezoid(weight, chis)
    mean_lens_chi = float(trapezoid(weight * chis, chis) / norm)
    wref = 1.0 - mean_lens_chi / cref
    return mean_lens_chi / (cref**2 * wref), mean_lens_chi


SHAPE_BINS = (((0.0, 10.0), (0.0, 10.0)),
              ((0.0, 10.0), (10.0, 30.0)),
              ((10.0, 20.0), (0.0, 10.0)),
              ((10.0, 20.0), (10.0, 30.0)),
              ((20.0, 30.0), (0.0, 10.0)),
              ((20.0, 30.0), (10.0, 30.0)))


@dataclass
class Config:
    r_perp_max: float = 30.0
    r_par_max: float = 30.0
    xi_step: float = 0.25
    xi_max: float = 40.0
    xi_smoothing: float = 0.9
    analytic_nk: int = 6400
    slab_index: int = -1
    shape_bins: tuple = SHAPE_BINS
    nside_jk: int = 8
    nside_alpha: int = 1024
    lmax_alpha: int = 2000
    chi_ref: float = field(default_factory=lambda: float(chi_of_z(2.4)))
    slabs: tuple = ((2.1, 3.0),)
    data_root: Path = field(default_factory=lambda: __import__("paths").MOCKS)
    report_root: Path = field(default_factory=lambda: CODE.parent / "report")
    seeds: tuple = tuple(range(20))
    g1: float = field(default_factory=lambda: kernel_product_g1()[0])
    response_delta: float = 2.0
    n_los: float = 22.0
    pixel_noise_power: float = 0.33
    scale: float = 1.0
    # Iteration 5: model-shaped xi table (xi_fit). The pair kernel accepts r_perp_min <= r_perp <= r_perp_max.
    forest_model: str = "kaiser"
    xi_correction: str = "none"          # "none" = iteration-5 two-parameter fit; "spline" = iteration-8 correction
    xi_correction_ridge: float = 1e-2
    xi_same_wavelength: bool = True      # False drops the first-bin term from the fit (diagnostic)
    xi_knots: str = "bicubic"            # "bicubic" (production) | "medium" | "fine" (diagnostics)
    fit_rperp_min: float = 2.0
    r_perp_min: float = 0.0
    los_pixel: float = 0.0
    los_resolution: float = 0.0

    def copy(self, **changes):
        return replace(self, **changes)


@dataclass
class SightlineSet:
    qid: np.ndarray
    ra: np.ndarray
    dec: np.ndarray
    zq: np.ndarray
    pix_start: np.ndarray
    chi: np.ndarray
    delta: np.ndarray
    w: np.ndarray
    slab: np.ndarray
    attrs: dict = field(default_factory=dict)
    region: np.ndarray | None = None      # per pixel: 0 = Lya region (A), 1 = Lyb region (B)

    def __post_init__(self):
        self.qid = np.asarray(self.qid, dtype=np.int64)
        self.ra = np.asarray(self.ra, dtype=np.float64)
        self.dec = np.asarray(self.dec, dtype=np.float64)
        self.zq = np.asarray(self.zq, dtype=np.float32)
        self.pix_start = np.asarray(self.pix_start, dtype=np.int64)
        self.chi = np.asarray(self.chi, dtype=np.float32)
        self.delta = np.asarray(self.delta, dtype=np.float32)
        self.w = np.asarray(self.w, dtype=np.float32)
        self.slab = np.asarray(self.slab, dtype=np.int8)
        self.region = (np.zeros(len(self.chi), np.int8) if self.region is None
                       else np.asarray(self.region, dtype=np.int8))
        if len(self.region) != len(self.chi):
            raise ValueError("region must have one entry per pixel")
        if len(self.pix_start) != len(self.qid) + 1:
            raise ValueError("pix_start must have Nq+1 elements")
        if self.pix_start[-1] != len(self.chi):
            raise ValueError("pix_start[-1] must equal Npix")

    @property
    def nq(self):
        return len(self.qid)
