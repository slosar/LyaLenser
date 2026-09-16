"""Spline-corrected forest correlation table for the pair response (iteration 8).

Iteration 5 fitted a two-parameter Kaiser model (`xi_fit.fit_model_table`) to the measured 1 Mpc/h table and used
the analytic derivative of that fit as the response kernel. On DR1 that model misses the data by 5-12 per cent
exactly where the lensing information lives (small r_par), and by a factor of up to 1.6 in the first radial bin
(the same-wavelength excess of sky and calibration residuals). Because the estimator normalisation is
A_hat = A * <W g_used g_true> / <W g_used^2>, a shape error in the kernel biases the amplitude at first order.

This module adds a smooth, analytically differentiable correction fitted to the same measured table:

    xi(rp, rz)      = xi_ref(rp, rz) + sigma(r) S(rp, rz) + N(rp) 1[rz < rz_sw]
    dxi/drp |lensed = dxi_ref/drp + d[sigma(r) S(rp, rz)]/drp

with

  * `xi_ref`  the model table of `xi_fit` (Kaiser basis, pixel windows, continuum projection): a smooth,
    physically motivated shape that controls the extrapolation where the data are noisy;
  * `S`       a tensor-product cubic B-spline in (r_perp, r_par), the empirical correction;
  * `sigma`   a fixed, strictly positive envelope (the monopole of `xi_ref`) that removes the dynamic range of
    xi from the spline coefficients. A *multiplicative* correction xi_ref (1 + S) blows up at the zero crossing
    of the anisotropic model; an additive correction on a positive envelope does not;
  * `N`       the same-wavelength excess: a 1-D cubic B-spline in r_perp confined to the first radial bin. It is
    an additive instrumental correlation between pixels at the same observed wavelength on different sightlines,
    so it belongs in the mean field (which subtracts <delta delta>) but carries no lensing response and is
    therefore excluded from the kernel.

Every basis element is averaged over the 1 Mpc/h cell exactly as the measurement bins the pair counts before the
weighted least-squares fit, and the derivative is the analytic derivative of the fitted surface (no smoothing
width, no finite differences). The base parameters (b_F^2, beta_F) and the correction are fitted alternately.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import numpy as np
from scipy.interpolate import BSpline, RegularGridInterpolator

try:
    from .xi_model import XiTable
except ImportError:
    from xi_model import XiTable


# --------------------------------------------------------------------------------------- 1-D spline bases
def clamped_knots(breaks, k=3):
    """Clamped knot vector from a sorted list of breakpoints (first and last are the interval ends)."""
    b = np.asarray(breaks, float)
    return np.r_[[b[0]] * k, b, [b[-1]] * k]


class Basis1D:
    """Cubic B-spline basis on [breaks[0], breaks[-1]]. Outside, the value is clamped and the slope is zero."""

    def __init__(self, breaks, k=3):
        b = np.asarray(breaks, float)
        self.lo = float(b[0]); self.hi = float(b[-1]); self.k = int(k)
        self.t = clamped_knots(b, k)
        self.n = len(self.t) - self.k - 1
        eye = np.eye(self.n)
        self._val = [BSpline(self.t, eye[i], self.k, extrapolate=False) for i in range(self.n)]
        self._der = [b.derivative() for b in self._val]

    def __call__(self, x, deriv=0):
        x = np.asarray(x, float)
        xc = np.clip(x, self.lo, self.hi)
        src = self._der if deriv else self._val
        out = np.stack([np.nan_to_num(f(xc)) for f in src], axis=-1)
        if deriv:
            out = out * ((x >= self.lo) & (x <= self.hi))[..., None]
        return out


# --------------------------------------------------------------------------------------- envelope
class Envelope:
    """Strictly positive scale sigma(r), r = hypot(r_perp, r_par), from the monopole of a reference table.

    Interpolated in log sigma against log r, so sigma > 0 everywhere and d sigma / d r is smooth.
    """

    def __init__(self, table: XiTable, n_mu=41, r_min=0.5, n_r=60, floor=1e-3):
        from scipy.interpolate import CubicSpline
        f = RegularGridInterpolator((table.r_perp, table.r_par), table.xi, bounds_error=False, fill_value=None)
        r_max = float(min(table.r_perp[-1], table.r_par[-1]))
        self.r = np.geomspace(r_min, r_max, n_r)
        mu = np.linspace(0.0, 1.0, n_mu)
        m = np.array([f(np.column_stack((r * np.sqrt(1 - mu ** 2), r * mu))).mean() for r in self.r])
        m = np.maximum(m, floor * np.max(np.abs(m)))
        # enforce monotone decrease so the envelope never turns up on noise of the reference monopole
        m = np.minimum.accumulate(m)
        self.lnm = np.log(m)
        self._spl = CubicSpline(np.log(self.r), self.lnm)      # exact derivative, no finite differences
        self._dspl = self._spl.derivative()

    def __call__(self, rp, rz, deriv=False):
        rp = np.asarray(rp, float); rz = np.asarray(rz, float)
        r = np.hypot(rp, rz)
        rc = np.clip(r, self.r[0], self.r[-1])
        s = np.exp(self._spl(np.log(rc)))
        if not deriv:
            return s
        inside = (r > self.r[0]) & (r < self.r[-1])
        dln_dlnr = self._dspl(np.log(rc))
        return s, s * dln_dlnr * rp / np.maximum(r, 1e-12) ** 2 * inside


# --------------------------------------------------------------------------------------- correction model
DEFAULT_RP_KNOTS = (3.0, 6.0, 10.0, 16.0, 30.0)
DEFAULT_RZ_KNOTS = (0.0, 4.0, 10.0, 30.0)
DEFAULT_SW_KNOTS = (3.0, 8.0, 16.0, 30.0)


@dataclass
class XiCorrection:
    """Design of the correction: tensor spline on a positive envelope plus the same-wavelength term.

    The r_par direction uses breakpoints in r_par^2, so every basis function is an even function of r_par with
    zero slope at r_par = 0, as the correlation must be. Without that the smooth surface can mimic the narrow
    same-wavelength excess in the first radial bin and the two become degenerate.
    """
    envelope: Envelope
    rp_knots: tuple = DEFAULT_RP_KNOTS
    rz_knots: tuple = DEFAULT_RZ_KNOTS
    sw_knots: tuple = DEFAULT_SW_KNOTS
    rz_sw: float = 1.0
    brp: Basis1D = field(init=False)
    brz: Basis1D = field(init=False)
    bsw: Basis1D = field(init=False)

    def __post_init__(self):
        self.brp = Basis1D(self.rp_knots)
        self.brz = Basis1D(np.asarray(self.rz_knots, float) ** 2)
        self.bsw = Basis1D(self.sw_knots)
        self.rp_lo = float(self.rp_knots[0]); self.rp_hi = float(self.rp_knots[-1])
        self.rz_hi = float(self.rz_knots[-1])

    @property
    def n_s(self):
        return self.brp.n * self.brz.n

    @property
    def n_par(self):
        return self.n_s + self.bsw.n

    def _tensor(self, rp, rz, deriv=0):
        A = self.brp(rp, deriv); B = self.brz(np.asarray(rz, float) ** 2, 0)
        return (A[..., :, None] * B[..., None, :]).reshape(A.shape[:-1] + (self.n_s,))

    def design(self, rp, rz):
        """Columns of the correction to xi."""
        s = self.envelope(rp, rz)
        sw = (np.abs(np.asarray(rz)) < self.rz_sw)[..., None] * self.bsw(rp)
        return np.concatenate([s[..., None] * self._tensor(rp, rz), sw], axis=-1)

    def design_deriv(self, rp, rz):
        """Columns of d/dr_perp of the LENSABLE part (the same-wavelength term is not lensed)."""
        s, ds = self.envelope(rp, rz, deriv=True)
        cols = ds[..., None] * self._tensor(rp, rz) + s[..., None] * self._tensor(rp, rz, deriv=1)
        return np.concatenate([cols, np.zeros(np.shape(rp) + (self.bsw.n,))], axis=-1)

    def penalty(self):
        """Second-difference roughness penalty on the coefficient grid."""
        def d2(n):
            if n < 3:
                return np.zeros((0, n))
            D = np.zeros((n - 2, n))
            for i in range(n - 2):
                D[i, i:i + 3] = (1.0, -2.0, 1.0)
            return D
        nrp, nrz = self.brp.n, self.brz.n
        P = np.zeros((self.n_par, self.n_par))
        Drp = np.kron(d2(nrp), np.eye(nrz)); Drz = np.kron(np.eye(nrp), d2(nrz))
        P[:self.n_s, :self.n_s] = Drp.T @ Drp + Drz.T @ Drz
        Dsw = d2(self.bsw.n)
        P[self.n_s:, self.n_s:] = Dsw.T @ Dsw
        return P


# --------------------------------------------------------------------------------------- fitting
def _cell_points(use, sub=4):
    o = (np.arange(sub) + 0.5) / sub
    idx = np.argwhere(use)
    RP = np.repeat((idx[:, 0][:, None] + o[None, :])[:, :, None], sub, axis=2)
    RZ = np.repeat((idx[:, 1][:, None] + o[None, :])[:, None, :], sub, axis=1)
    return RP, RZ


def fit_corrected_table(num, den, basis_proj, basis_coarse, cfg, corr: XiCorrection,
                        ridge=0.0, sub=4, iterations=3):
    """Fit (b_F^2, beta_F) of the projected model basis and the spline correction to the measured cells.

    Returns (XiTable, params). The table's ``xi`` is the full model (correction and same-wavelength term
    included) and ``xi_rp`` is the derivative of the lensable part only.
    """
    from xi_fit import BASIS, coarse_bin
    from scipy.optimize import least_squares

    raw = np.divide(num, den, out=np.zeros_like(num), where=den > 0)
    n = raw.shape[0]
    rp_c = np.arange(n) + 0.5
    rz_c = np.arange(n) + 0.5
    use = ((den > 0) & (rp_c[:, None] >= cfg.fit_rperp_min) & (rp_c[:, None] < cfg.r_perp_max)
           & (rz_c[None, :] < cfg.r_par_max))
    coarse = {k: (coarse_bin(basis_proj[k], cfg) if basis_coarse is None else basis_coarse[k]) for k in BASIS}

    RP, RZ = _cell_points(use, sub)
    A = corr.design(RP, RZ).mean(axis=(1, 2))
    B = np.column_stack([coarse[k][use] for k in BASIS])
    y = raw[use]
    w = np.sqrt(den[use])
    P = corr.penalty()

    def base_model(p):
        b2, be = p
        return b2 * (B[:, 0] + 2 * be * B[:, 1] + be * be * B[:, 2])

    lin = np.linalg.lstsq(B * w[:, None], y * w, rcond=None)[0]
    p = np.array([max(lin[0], 1e-6), lin[1] / (2 * max(lin[0], 1e-6))])
    c = np.zeros(corr.n_par)
    for _ in range(iterations):
        r = least_squares(lambda q: (base_model(q) + A @ c - y) * w, p,
                          x_scale=[max(abs(p[0]), 1e-6), 1.0])
        p = r.x
        M = A * w[:, None]
        scale = np.maximum(np.linalg.norm(M, axis=0), 1e-30)
        Ms = M / scale
        b = (y - base_model(p)) * w
        if ridge > 0:
            c = np.linalg.solve(Ms.T @ Ms + ridge * (P / np.outer(scale, scale)), Ms.T @ b) / scale
        else:
            c = np.linalg.lstsq(Ms, b, rcond=None)[0] / scale

    model_cells = base_model(p) + A @ c
    resid = (y - model_cells) * w
    chi2 = float(np.sum(resid ** 2))

    b2, be = p
    coef = np.array([b2, 2 * b2 * be, b2 * be * be])
    ref = basis_proj[BASIS[0]]
    xi_ref = sum(ci * basis_proj[k].xi for ci, k in zip(coef, BASIS))
    g_ref = sum(ci * basis_proj[k].xi_rp for ci, k in zip(coef, BASIS))
    GRP, GRZ = np.meshgrid(ref.r_perp, ref.r_par, indexing='ij')
    xi = xi_ref + corr.design(GRP, GRZ) @ c
    g = g_ref + corr.design_deriv(GRP, GRZ) @ c

    sw_only = np.concatenate([np.zeros(corr.n_s), c[corr.n_s:]])
    params = dict(b_F2=float(b2), beta_F=float(be), chi2=chi2, cells=int(use.sum()),
                  n_correction=int(corr.n_par), dof=int(use.sum() - 2 - corr.n_par), ridge=float(ridge),
                  knots=dict(r_perp=list(corr.rp_knots), r_par=list(corr.rz_knots),
                             same_wavelength=list(corr.sw_knots), rz_sw=corr.rz_sw,
                             note='r_par breakpoints are applied in r_par^2 (even correlation)'),
                  fit_range={'rperp_min': float(cfg.fit_rperp_min), 'rperp_max': float(cfg.r_perp_max),
                             'rpar_max': float(cfg.r_par_max)},
                  coefficients=c.tolist(),
                  same_wavelength_rperp=[float(x) for x in
                                         (corr.design(np.array([3.5, 5.5, 10.5, 20.5, 29.5]),
                                                      np.zeros(5)) @ sw_only)])
    meta = dict(provider='model_fit+spline', basis=list(BASIS),
                accepted_weight=float(den.sum()),
                note='xi includes the same-wavelength term; xi_rp is the lensable derivative only')
    table = XiTable(ref.r_perp, ref.r_par, xi, g, meta, (num, den))
    table.meta['fit'] = params
    return table, params, dict(use=use, model_cells=model_cells, raw=raw, A=A, base=base_model(p))
