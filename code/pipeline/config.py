"""Configuration and lightweight data containers for the Stage-A pipeline."""
from __future__ import annotations

from dataclasses import dataclass, field, replace
from pathlib import Path
import sys
import numpy as np

CODE = Path(__file__).resolve().parents[1]
if str(CODE) not in sys.path:
    sys.path.insert(0, str(CODE))
from cosmo import chi as chi_of_z


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
    shape_bins: tuple = SHAPE_BINS
    nside_jk: int = 8
    nside_alpha: int = 1024
    lmax_alpha: int = 2000
    chi_ref: float = field(default_factory=lambda: float(chi_of_z(2.4)))
    slabs: tuple = ((2.1, 3.0),)
    data_root: Path = Path("/data/LyaLenser/mocks")
    report_root: Path = field(default_factory=lambda: CODE.parent / "report")
    seeds: tuple = tuple(range(20))
    g1: float = field(default_factory=lambda: (float(chi_of_z(1.0))/float(chi_of_z(2.4))**2) /
                      (1.0-float(chi_of_z(1.0))/float(chi_of_z(2.4))))
    response_delta: float = 2.0
    n_los: float = 22.0
    pixel_noise_power: float = 0.33
    scale: float = 1.0

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
        if len(self.pix_start) != len(self.qid) + 1:
            raise ValueError("pix_start must have Nq+1 elements")
        if self.pix_start[-1] != len(self.chi):
            raise ValueError("pix_start[-1] must equal Npix")

    @property
    def nq(self):
        return len(self.qid)
