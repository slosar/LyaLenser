"""Redshift-evolving forest correlation table for the pair response (Stage B iteration 10).

The iteration-8 table (`xi_spline`) is one surface xi(r_perp, r_par) fitted to the pair-weighted cells of the whole
slab 2.1 < z < 3.0. The A x B pairs (one pixel in the Lyb region) sit at lower redshift than the A x A pairs and
their correlation amplitude is 0.84 of the A x A one; the natural reading is redshift evolution of the forest
correlation, not a property of the region-B deltas. Here the measured cells are split by the pair mean redshift
(`Config.xi_z_edges`) and one model with a few power laws in x = (1 + z)/(1 + z_ref) is fitted to all of them:

    xi(rp, rz, z) = b_F^2 [D(z)/D(z_ref)]^2 x^{2 gamma_b} [T0 + 2 beta(z) T1 + beta(z)^2 T2](rp, rz)
                    + x^{gamma_S} sigma(r) S(rp, rz) + N(rp) 1[rz < rz_sw],        beta(z) = beta_F x^{gamma_beta},

with (T0, T1, T2) the projected basis of `xi_fit` (P_lin at z_ref, pixel windows, continuum projection), D the
linear growth factor, S the tensor-spline correction of `xi_spline` on its positive envelope sigma, and N the
same-wavelength term (instrumental, not evolved, not lensed). The three exponents (gamma_b, gamma_beta, gamma_S)
are the only new parameters. The fitted table is LAYERED: one 2-D surface per node of a uniform chi grid spanning
the slab, and the pair kernels interpolate linearly between the layers at the pair's mean distance.
"""
from __future__ import annotations

from pathlib import Path
import sys
import numpy as np
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
for p in (HERE, HERE.parent):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))
from xi_model import XiTable
from xi_fit import BASIS, coarse_bin, FittedTable
from xi_spline import XiCorrection, _cell_points
from cosmo import growth, chi as chi_of_z, z_of_chi

PARAMS = ('b_F2', 'beta_F', 'gamma_b', 'gamma_beta', 'gamma_S')
DEFAULT_START = {'gamma_b': 2.9, 'gamma_beta': 0.0, 'gamma_S': 3.8}
BOUNDS = {'b_F2': (1e-8, np.inf), 'beta_F': (-5., 10.), 'gamma_b': (-10., 10.), 'gamma_beta': (-10., 10.), 'gamma_S': (-15., 15.)}


def zfactor(z, z_ref):
    return (1.0 + np.asarray(z, float)) / (1.0 + z_ref)


def growth_ratio2(z, z_ref):
    return (growth(z) / growth(z_ref)) ** 2


def base_coefficients(b2, beta, gamma_b, gamma_beta, z, z_ref):
    """(c0, c1, c2) multiplying the projected basis (mu0, mu2, mu4) at redshift z."""
    x = zfactor(z, z_ref)
    amp = b2 * growth_ratio2(z, z_ref) * x ** (2 * gamma_b)
    bz = beta * x ** gamma_beta
    return amp, amp * 2 * bz, amp * bz * bz


def chi_nodes_for(cfg, z_edges):
    """Uniform chi grid of the layered table: from one node step below the first z edge to one above the last."""
    step = float(getattr(cfg, 'xi_z_node_step', 0.05))
    zlo = float(min(z_edges)) - step; zhi = float(max(z_edges)) + step
    z_ref = float(getattr(cfg, 'xi_z_ref', 2.4))
    dchi = float(chi_of_z(z_ref + step) - chi_of_z(z_ref))
    c_lo = float(chi_of_z(zlo)); c_hi = float(chi_of_z(zhi))
    n = int(np.ceil((c_hi - c_lo) / dchi)) + 1
    return c_lo + dchi * np.arange(n)


def fit_evolving_table(counts_z, basis_proj, basis_coarse, cfg, corr: XiCorrection | None = None, ridge=0.0,
                       sub=4, iterations=4, fixed=None, start=None, chi_nodes=None):
    """Fit the evolving model to the z-resolved cells and return (layered XiTable, params, debug).

    ``counts_z`` = (num, den, chisum), each [n_z, n, n] (from `xi_model.xi_from_data(...).counts_z`); the cell
    redshift is the pair-weighted mean redshift of the cell. ``corr`` None fits the base model only (the evolving
    analogue of the two-parameter Kaiser fit). ``fixed`` maps parameter names to values held fixed (e.g. the
    exponents of a joint fit when refitting one pair type). Weighted least squares with weights = den.
    """
    num, den, csum = (np.asarray(a, float) for a in counts_z)
    nz, n, _ = num.shape
    z_ref = float(getattr(cfg, 'xi_z_ref', 2.4))
    raw = np.divide(num, den, out=np.zeros_like(num), where=den > 0)
    zc = np.full(den.shape, z_ref)
    ok = den > 0
    zc[ok] = z_of_chi(csum[ok] / den[ok])
    rp_c = np.arange(n) + 0.5; rz_c = np.arange(n) + 0.5
    use2 = (rp_c[:, None] >= cfg.fit_rperp_min) & (rp_c[:, None] < cfg.r_perp_max) & (rz_c[None, :] < cfg.r_par_max)
    use = ok & use2[None, :, :]
    coarse = {k: (coarse_bin(basis_proj[k], cfg) if basis_coarse is None else basis_coarse[k]) for k in BASIS}
    # cell quantities, flattened over the used (z, rp, rz) cells
    iz, ii, jj = np.nonzero(use)
    y = raw[use]; w = np.sqrt(den[use]); z = zc[use]
    x = zfactor(z, z_ref); D2 = growth_ratio2(z, z_ref)
    B = np.column_stack([coarse[k][ii, jj] for k in BASIS])
    if corr is not None:
        RP, RZ = _cell_points(use2, sub)
        A2 = corr.design(RP, RZ).mean(axis=(1, 2))          # [n_cells_2d, n_par] on the 2-D cell list
        idx2 = -np.ones((n, n), int); idx2[use2] = np.arange(int(use2.sum()))
        A = A2[idx2[ii, jj]]
        A_S = A[:, :corr.n_s]; A_sw = A[:, corr.n_s:]
        P = corr.penalty()
    else:
        A_S = np.zeros((len(y), 0)); A_sw = np.zeros((len(y), 0)); P = np.zeros((0, 0))

    fixed = dict(fixed or {}); start = dict(DEFAULT_START, **(start or {}))
    lin = np.linalg.lstsq(B * w[:, None], y * w, rcond=None)[0]
    start.setdefault('b_F2', max(lin[0], 1e-6)); start.setdefault('beta_F', lin[1] / (2 * max(lin[0], 1e-6)))
    free = [k for k in PARAMS if k not in fixed]
    full = np.array([fixed.get(k, start[k]) for k in PARAMS], float)

    def unpack(q):
        v = full.copy()
        for k, val in zip(free, q):
            v[PARAMS.index(k)] = val
        return v

    def base_model(v):
        b2, be, gb, gbe, _ = v
        amp = b2 * D2 * x ** (2 * gb); bz = be * x ** gbe
        return amp * (B[:, 0] + 2 * bz * B[:, 1] + bz * bz * B[:, 2])

    def corr_design(v):
        gS = v[4]
        return np.column_stack([x[:, None] ** gS * A_S, A_sw]) if corr is not None else np.zeros((len(y), 0))

    def linear_coefficients(v):
        """Spline coefficients at fixed nonlinear parameters (variable projection)."""
        if corr is None:
            return np.zeros(0)
        M = corr_design(v) * w[:, None]
        scale = np.maximum(np.linalg.norm(M, axis=0), 1e-30); Ms = M / scale
        b = (y - base_model(v)) * w
        if ridge > 0:
            return np.linalg.solve(Ms.T @ Ms + ridge * (P / np.outer(scale, scale)), Ms.T @ b) / scale
        return np.linalg.lstsq(Ms, b, rcond=None)[0] / scale

    def residual(qq):
        v = unpack(qq); c = linear_coefficients(v)
        return (base_model(v) + corr_design(v) @ c - y) * w

    q = np.array([full[PARAMS.index(k)] for k in free])
    lo = np.array([BOUNDS[k][0] for k in free]); hi = np.array([BOUNDS[k][1] for k in free])
    q = np.clip(q, lo + 1e-9, hi - 1e-9)
    if free:
        for _ in range(int(iterations)):        # restarts from the previous solution tighten the tolerances
            r = least_squares(residual, q, bounds=(lo, hi), x_scale=np.maximum(np.abs(q), 1e-3),
                              ftol=1e-12, xtol=1e-12, gtol=1e-12)
            q = r.x
    v = unpack(q); c = linear_coefficients(v)
    model_cells = base_model(v) + corr_design(v) @ c
    resid = (y - model_cells) * w; chi2 = float(np.sum(resid ** 2))
    per_z = [{'z_bin': int(k), 'cells': int(np.sum(iz == k)), 'chi2': float(np.sum(resid[iz == k] ** 2)),
              'mean_z': float(np.average(z[iz == k], weights=w[iz == k] ** 2)) if np.any(iz == k) else None}
             for k in range(nz)]

    # ---- layered table
    z_edges = getattr(cfg, 'xi_z_edges', None)
    if chi_nodes is None:
        chi_nodes = chi_nodes_for(cfg, z_edges if z_edges else (z_ref - 0.3, z_ref + 0.6))
    chi_nodes = np.asarray(chi_nodes, float)
    ref = basis_proj[BASIS[0]]
    GRP, GRZ = np.meshgrid(ref.r_perp, ref.r_par, indexing='ij')
    if corr is not None:
        Dv = corr.design(GRP, GRZ); Dg = corr.design_deriv(GRP, GRZ)
        cS = np.r_[c[:corr.n_s], np.zeros(corr.n_sw)]; csw = np.r_[np.zeros(corr.n_s), c[corr.n_s:]]
        corr_S = Dv @ cS; corr_sw = Dv @ csw; corr_gS = Dg @ cS          # the derivative of the sw block is zero
    xi = np.zeros((len(chi_nodes),) + ref.xi.shape); g = np.zeros_like(xi)
    node_z = z_of_chi(chi_nodes)
    for k, zk in enumerate(node_z):
        c0, c1, c2 = base_coefficients(v[0], v[1], v[2], v[3], float(zk), z_ref)
        xi[k] = c0 * basis_proj['mu0'].xi + c1 * basis_proj['mu2'].xi + c2 * basis_proj['mu4'].xi
        g[k] = c0 * basis_proj['mu0'].xi_rp + c1 * basis_proj['mu2'].xi_rp + c2 * basis_proj['mu4'].xi_rp
        if corr is not None:
            xk = float(zfactor(zk, z_ref)) ** v[4]
            xi[k] += xk * corr_S + corr_sw; g[k] += xk * corr_gS
    # the amplitude evolution actually implied: d ln(b^2 D^2) / d ln(1+z) at z_ref
    dz = 1e-3
    dlnD2 = float(np.log(growth_ratio2(z_ref + dz, z_ref) / growth_ratio2(z_ref - dz, z_ref)) / np.log((1 + z_ref + dz) / (1 + z_ref - dz)))
    params = dict(zip(PARAMS, (float(t) for t in v)))
    params.update(z_ref=z_ref, chi2=chi2, cells=int(len(y)), n_correction=int(len(c)),
                  dof=int(len(y) - len(free) - len(c)), free=free, fixed={k: float(val) for k, val in fixed.items()},
                  ridge=float(ridge), per_z=per_z, z_edges=(list(map(float, z_edges)) if z_edges else None),
                  amplitude_exponent_total=float(2 * v[2] + dlnD2), growth_exponent=dlnD2,
                  fit_range={'rperp_min': float(cfg.fit_rperp_min), 'rperp_max': float(cfg.r_perp_max), 'rpar_max': float(cfg.r_par_max)},
                  coefficients=c.tolist(), chi_nodes=chi_nodes.tolist(), node_z=[float(t) for t in node_z],
                  model=('base + spline correction (x^gamma_S) + same-wavelength term' if corr is not None else 'base only'))
    if corr is not None:
        params['knots'] = dict(r_perp=list(corr.rp_knots), r_par=list(corr.rz_knots), same_wavelength=list(corr.sw_knots), rz_sw=corr.rz_sw)
    meta = dict(provider='model_fit+zevol' + ('+spline' if corr is not None else ''), basis=list(BASIS),
                accepted_weight=float(den.sum()), z_ref=z_ref,
                note='layered in chi; xi includes the same-wavelength term; xi_rp is the lensable derivative only')
    table = XiTable(ref.r_perp, ref.r_par, xi, g, meta, (num.sum(axis=0), den.sum(axis=0)), chi_nodes=chi_nodes)
    table.meta['fit'] = params
    return table, params, dict(use=use, model_cells=model_cells, raw=raw, z=z, iz=iz, ii=ii, jj=jj, resid=resid)


def evolving_table_for(counts_z, cfg, basis, fixed=None, correction=None):
    """`run_mock_validation.table_for` analogue for the evolving model: ``correction`` None follows
    ``cfg.xi_correction`` ('spline' -> the production spline correction, 'none' -> base only)."""
    from run_mock_validation import correction_for
    use_corr = (getattr(cfg, 'xi_correction', 'none') == 'spline') if correction is None else bool(correction)
    corr = correction_for(basis, cfg) if use_corr else None
    tab, par, _ = fit_evolving_table(counts_z, basis['projected'], basis['coarse'], cfg, corr,
                                     ridge=(cfg.xi_correction_ridge if use_corr else 0.0), fixed=fixed)
    return FittedTable(tab, par)
