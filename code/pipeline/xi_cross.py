"""Quasar x forest cross-correlation for the pair lensing estimator (Stage B iteration 12).

For a quasar at (theta_q, chi_q) and a forest pixel at (theta_p, chi_p) on another sightline, lensing remaps both
positions, so the observed cross-correlation is xi_qF at the true separation,
    <delta_F(p) | quasar at q> = xi_qF(r_perp, r_par) + xi'_qF chi theta_hat_qp . (alpha_q - alpha_p),
exactly the pair form of the forest auto-correlation with (a, b) -> (q, p). Magnification bias modulates the
observed quasar density but, once the mean field is subtracted per OBSERVED quasar, contributes nothing at linear
order in a foreground template (its residual is second order and dropped).

This module provides
  * the measured cells: quasar-pixel pairs on different sightlines, 1 Mpc/h cells in (r_perp, r_par) with r_par
    SIGNED (pixel minus quasar; positive behind the quasar), split by the pair mean redshift;
  * the one-sided continuum projection of the model basis (the projector acts on the forest pixels only);
  * the evolving model fit: P_qF = b_q(z) b_F(z) (1 + beta_q mu^2)(1 + beta_F mu^2) P(k) with the forest side
    (b_F, beta_F and their power laws) FIXED from the auto-correlation fit, the quasar bias b_q x^gamma_q free
    and beta_q(z) = f(z)/b_q(z) set by the cosmology, plus the quasar redshift offset Delta r_par and Gaussian
    redshift-error smoothing sigma_par along the line of sight, and a spline correction (odd terms allowed) on a
    positive envelope. The fitted table is layered in chi like the auto one;
  * the quasar-pixel pair accumulator with the eleven accumulators of `pairs.accumulate`, indexing the quasar as
    position N_sightlines + q so that `amplitude.amplitude` runs unchanged on templates evaluated at the
    concatenated (sightline, quasar) positions.
"""
from __future__ import annotations
from pathlib import Path
import sys
import numpy as np
from numba import njit, prange
from scipy.spatial import cKDTree
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import least_squares
HERE = Path(__file__).resolve().parent
for p in (HERE, HERE.parent):
    if str(p) not in sys.path: sys.path.insert(0, str(p))
from xi_model import XiTable
from xi_fit import BASIS, _weighted_projector
from xi_spline import Basis1D, Envelope, _cell_points
from xi_zevol import zfactor, growth_ratio2, chi_nodes_for
from pairs import _unit_vectors, pair_geometry, _interp_layer, _shape_bin, PairCatalogue
from cosmo import chi as chi_of_z, z_of_chi, growth_rate


# ------------------------------------------------------------------------------------------- pairs
def find_cross_pairs(sl, qso, theta_max):
    """(a, b, thx, thy, theta) for quasar-sightline pairs within theta_max, a = N_sl + quasar index, b = sightline
    index, with the pair direction from the sightline to the quasar (theta_a - theta_b, as in `find_pairs`).
    The quasar's own sightline (same TARGETID) is excluded."""
    tree = cKDTree(_unit_vectors(sl.ra, sl.dec)); vq = _unit_vectors(qso.ra, qso.dec)
    lists = tree.query_ball_point(vq, 2 * np.sin(theta_max / 2))
    nq = np.array([len(l) for l in lists]); q = np.repeat(np.arange(qso.n), nq); s = np.concatenate([np.asarray(l, np.int64) for l in lists]) if nq.sum() else np.zeros(0, np.int64)
    keep = sl.qid[s] != qso.qid[q]
    q = q[keep]; s = s[keep]
    ra = np.r_[sl.ra, qso.ra]; dec = np.r_[sl.dec, qso.dec]; a = (sl.nq + q).astype(np.int32); b = s.astype(np.int32)
    thx, thy, theta = pair_geometry(ra, dec, a, b)
    m = theta <= theta_max
    return a[m], b[m], thx[m], thy[m], theta[m]


# ------------------------------------------------------------------------------------------- measured cells
@njit(parallel=True, cache=True)
def _cross_hist_kernel(pix_start, chi, delta, weight, region, chiq, a, b, theta, nsl, rpmax, rzmax, step, chi_zedges):
    """num, den, chisum [2 (pixel region), n_z, n_rp, 2 n_rz] over quasar-pixel pairs; r_par = chi_p - chi_q."""
    nr = int(np.ceil(rpmax / step)); nz = chi_zedges.size - 1
    nums = np.zeros((24, 2, nz, nr, 2 * nr)); dens = np.zeros((24, 2, nz, nr, 2 * nr)); csum = np.zeros((24, 2, nz, nr, 2 * nr))
    for chunk in prange(24):
        for ip in range(chunk * a.size // 24, (chunk + 1) * a.size // 24):
            cq = float(chiq[a[ip] - nsl]); s = b[ip]
            for p in range(pix_start[s], pix_start[s + 1]):
                cp = float(chi[p]); rz = cp - cq
                if rz < -rzmax or rz >= rzmax: continue
                cm = .5 * (cp + cq); rp = cm * theta[ip]
                if rp >= rpmax: continue
                iz = -1
                for k in range(nz):
                    if chi_zedges[k] <= cm < chi_zedges[k + 1]:
                        iz = k; break
                if iz < 0: continue
                i = int(rp / step); j = int((rz + rzmax) / step); reg = 1 if region[p] > 0 else 0
                w = float(weight[p])
                nums[chunk, reg, iz, i, j] += w * float(delta[p]); dens[chunk, reg, iz, i, j] += w; csum[chunk, reg, iz, i, j] += w * cm
    return nums.sum(axis=0), dens.sum(axis=0), csum.sum(axis=0)


def xi_qf_from_data(sl, qso, cfg, pairs=None):
    """Measured quasar-forest cells over r_perp < xi_max, |r_par| < xi_max, in the redshift bins of the evolving fit.
    Returns dict(all=(num, den, chisum), A=..., B=..., chi_edges, z_edges, pairs)."""
    if pairs is None: pairs = find_cross_pairs(sl, qso, cfg.xi_max / max(1., float(min(sl.chi.min(), qso.chi.min()))))
    zed = np.asarray(getattr(cfg, 'xi_z_edges', (2.1, 3.0)), float)
    chi_zedges = np.array([float(chi_of_z(z)) for z in zed]); chi_zedges[-1] = np.nextafter(chi_zedges[-1], np.inf)
    num, den, cs = _cross_hist_kernel(sl.pix_start, sl.chi, sl.delta, sl.w, sl.region, np.asarray(qso.chi, np.float32), pairs[0], pairs[1], pairs[4],
                                      sl.nq, float(cfg.xi_max), float(cfg.xi_max), 1.0, chi_zedges)
    out = {'all': (num.sum(axis=0), den.sum(axis=0), cs.sum(axis=0)), 'A': (num[0], den[0], cs[0]), 'B': (num[1], den[1], cs[1]),
           'chi_edges': chi_zedges, 'z_edges': zed, 'pairs': pairs, 'n_pairs': int(len(pairs[0])),
           'pair_weight_by_z': den.sum(axis=(0, 2, 3)).tolist(), 'pair_weight_by_region': [float(den[0].sum()), float(den[1].sum())]}
    return out


# ------------------------------------------------------------------------------------------- projection
def project_cross_sample(tables: dict, forests, cfg, n_forests=400, seed=11, step=None, rp_grid=None, q_step=4.0, verbose=False):
    """One-sided continuum projection of stationary tables xi(r_perp, r_par): for sampled real forests b (whole
    picca span, weights, in-range mask, region blocks) and quasar distances chi_q on a grid across and beyond the
    forest, average the projected column (M - u v^T M)_p over pixels p at lag chi_p - chi_q, weighted by w_p.
    Returns {name: XiTable} on (rp_grid, lags in [-xi_max, xi_max]) with the SIGNED r_par axis."""
    forests = [(f[0], f[1], (np.asarray(f[2], float) if len(f) > 2 else np.ones(len(f[0]))),
                (np.asarray(f[3], np.int8) if len(f) > 3 else np.zeros(len(f[0]), np.int8))) for f in forests]
    rng = np.random.default_rng(seed)
    if step is None: step = float(np.median(np.diff(np.unique(np.round(np.concatenate([f[0] for f in forests]), 4)))))
    fine = np.arange(0, cfg.xi_max + .5 * cfg.xi_step, cfg.xi_step); rp_grid = fine if rp_grid is None else np.asarray(rp_grid, float)
    nlag = 2 * int(cfg.xi_max / step) + 3; lags = (np.arange(nlag) - (nlag - 1) // 2) * step
    from scipy.interpolate import RegularGridInterpolator as _RGI
    names = list(tables)
    interp = [_RGI((tables[k].r_perp, tables[k].r_par), a, bounds_error=False, fill_value=0.) for k in names for a in (tables[k].xi, tables[k].xi_rp)]
    # the raw tables are even in r_par: evaluate at |lag|
    rows = np.array([[f((np.full(nlag, rp), np.abs(lags))) for rp in rp_grid] for f in interp]); nsurf = rows.shape[0]
    numr = np.zeros((nsurf, len(rp_grid), nlag)); dend = np.zeros(nlag)
    idx = rng.choice(len(forests), min(n_forests, len(forests)), replace=False)
    for it, ib in enumerate(idx):
        cb, wb, mb, rb = forests[ib]; ub, vb = _weighted_projector(cb, wb, rb)
        for cq in np.arange(cb.min() - cfg.xi_max, cb.max() + cfg.xi_max + 1e-6, q_step):
            k = np.rint((cb - cq) / step).astype(np.int64) + (nlag - 1) // 2
            keep = (k >= 0) & (k < nlag); kk = np.clip(k, 0, nlag - 1)
            W = wb * mb * keep
            dend += np.bincount(kk, weights=W, minlength=nlag)[:nlag]
            for si in range(nsurf):
                M = rows[si][:, kk] * keep[None, :]                 # [n_rp, n_pix]
                M = M - (M @ vb) @ ub.T                              # projector on the pixel side
                for ir in range(len(rp_grid)):
                    numr[si, ir] += np.bincount(kk, weights=W * M[ir], minlength=nlag)[:nlag]
        if verbose and it % 50 == 0: print(f'  cross projection forest {it}/{len(idx)}', flush=True)
    prof = np.divide(numr, np.maximum(dend, 1e-300), out=np.zeros_like(numr), where=dend > 0); ok = dend > 0
    lag_fine = np.arange(-cfg.xi_max, cfg.xi_max + .5 * cfg.xi_step, cfg.xi_step); out = {}
    for i, k in enumerate(names):
        xi = np.array([np.interp(lag_fine, lags[ok], row[ok]) for row in prof[2 * i]])
        g = np.array([np.interp(lag_fine, lags[ok], row[ok]) for row in prof[2 * i + 1]])
        out[k] = XiTable(rp_grid, lag_fine, xi, g, dict(tables[k].meta or {}, projection=f'one-sided, {len(idx)} real forests, quasar step {q_step}', pixel_step=step))
    return out


def refine_cross_grid(tables, cfg):
    from scipy.interpolate import CubicSpline
    fine = np.arange(0, cfg.xi_max + .5 * cfg.xi_step, cfg.xi_step); out = {}
    for k, t in tables.items():
        out[k] = XiTable(fine, t.r_par, CubicSpline(t.r_perp, t.xi, axis=0)(fine), CubicSpline(t.r_perp, t.xi_rp, axis=0)(fine), dict(t.meta or {}, rp_grid='cubic refinement'))
    return out


# ------------------------------------------------------------------------------------------- fit
class CrossCorrection:
    """Tensor cubic B-spline in (r_perp, r_par) on a positive envelope; odd terms in r_par allowed (the quasar
    redshift offset makes the correlation asymmetric)."""
    def __init__(self, envelope, rp_knots=(3., 30.), rz_knots=(-30., 0., 30.)):
        self.envelope = envelope; self.brp = Basis1D(rp_knots); self.brz = Basis1D(rz_knots); self.rp_knots = rp_knots; self.rz_knots = rz_knots
        self.n_par = self.brp.n * self.brz.n; self.n_s = self.n_par; self.n_sw = 0
    def _tensor(self, rp, rz, deriv=0):
        A = self.brp(rp, deriv); B = self.brz(rz, 0); return (A[..., :, None] * B[..., None, :]).reshape(A.shape[:-1] + (self.n_par,))
    def design(self, rp, rz): return self.envelope(rp, rz)[..., None] * self._tensor(rp, rz)
    def design_deriv(self, rp, rz):
        s, ds = self.envelope(rp, rz, deriv=True); return ds[..., None] * self._tensor(rp, rz) + s[..., None] * self._tensor(rp, rz, deriv=1)
    def penalty(self):
        def d2(n):
            D = np.zeros((max(n - 2, 0), n))
            for i in range(n - 2): D[i, i:i + 3] = (1., -2., 1.)
            return D
        Drp = np.kron(d2(self.brp.n), np.eye(self.brz.n)); Drz = np.kron(np.eye(self.brp.n), d2(self.brz.n)); return Drp.T @ Drp + Drz.T @ Drz


CROSS_PARAMS = ('b_q', 'gamma_q', 'dr_par', 'sigma_par', 'gamma_S')
CROSS_BOUNDS = {'b_q': (0.2, 10.), 'gamma_q': (-10., 10.), 'dr_par': (-10., 10.), 'sigma_par': (0., 15.), 'gamma_S': (-15., 15.)}
CROSS_START = {'b_q': 3.0, 'gamma_q': 0.0, 'dr_par': 0.0, 'sigma_par': 2.0, 'gamma_S': 5.0}


def _coarse_matrix(fine_axis, sub):
    """Averaging matrix from a fine axis (step h) to unit cells [k, k+1): each cell averages the ``sub`` fine points
    it contains (fine points at multiples of h = 1/sub). Returns (cell_lo_values, matrix [n_cells, n_fine])."""
    h = fine_axis[1] - fine_axis[0]; lo = np.floor(fine_axis[0] + 1e-9); n = int(np.rint((fine_axis[-1] - lo) * (1. / h))) // sub
    M = np.zeros((n, len(fine_axis)))
    for c in range(n):
        for s in range(sub):
            x = lo + c + s / sub; j = int(np.rint((x - fine_axis[0]) / h))
            if 0 <= j < len(fine_axis): M[c, j] = 1. / sub
    return lo + np.arange(n), M


def transform_along_rpar(arr, r_par, dr, sigma):
    """Shift by dr (the quasar redshift offset: the observed correlation is the true one at r_par - dr) and smooth
    with a Gaussian of width sigma along the SIGNED r_par axis of a fine table [n_rp, n_rz]."""
    h = r_par[1] - r_par[0]; out = arr
    if sigma > 0: out = gaussian_filter1d(out, sigma / h, axis=1, mode='nearest')
    if dr != 0: out = np.array([np.interp(r_par - dr, r_par, row) for row in out])
    return out


def fit_evolving_cross(counts_z, basis_proj, cfg, forest, corr=None, ridge=0.0, sub=4, fixed=None, start=None, chi_nodes=None, iterations=2):
    """Fit the evolving quasar-forest model to the z-resolved signed cells.

    ``forest`` = dict(b_F2, beta_F, gamma_b, gamma_beta) from the auto-correlation fit (b_F = -sqrt(b_F2) < 0).
    Free: b_q (at z_ref), gamma_q, dr_par, sigma_par, gamma_S (if corr); the spline coefficients are solved
    linearly inside the residual (variable projection). Returns (layered XiTable, params, debug)."""
    num, den, csum = (np.asarray(a, float) for a in counts_z); nz, nr, nz2 = num.shape
    z_ref = float(getattr(cfg, 'xi_z_ref', 2.4)); raw = np.divide(num, den, out=np.zeros_like(num), where=den > 0)
    zc = np.full(den.shape, z_ref); ok = den > 0; zc[ok] = z_of_chi(csum[ok] / den[ok])
    rp_c = np.arange(nr) + .5; rz_c = np.arange(nz2) - nz2 / 2 + .5
    use2 = (rp_c[:, None] >= cfg.fit_rperp_min) & (rp_c[:, None] < cfg.r_perp_max) & (np.abs(rz_c)[None, :] < cfg.r_par_max)
    use = ok & use2[None, :, :]; iz, ii, jj = np.nonzero(use)
    y = raw[use]; w = np.sqrt(den[use]); z = zc[use]; x = zfactor(z, z_ref); D2 = growth_ratio2(z, z_ref); fz = growth_rate(z)
    bF = -np.sqrt(forest['b_F2']) * x ** forest['gamma_b']; betaF = forest['beta_F'] * x ** forest['gamma_beta']
    ref = basis_proj[BASIS[0]]; rp_lo, Mrp = _coarse_matrix(ref.r_perp, sub); rz_lo, Mrz = _coarse_matrix(ref.r_par, sub)
    assert len(rp_lo) == nr and len(rz_lo) == nz2, (len(rp_lo), nr, len(rz_lo), nz2)
    fine = {k: basis_proj[k].xi for k in BASIS}
    def coarse_basis(dr, sig):
        return {k: Mrp @ transform_along_rpar(fine[k], ref.r_par, dr, sig) @ Mrz.T for k in BASIS}
    if corr is not None:
        RP, RZ = _cell_points(use2, sub); RZ = RZ - nz2 / 2          # signed r_par cell points
        A2 = corr.design(RP, RZ).mean(axis=(1, 2)); idx2 = -np.ones((nr, nz2), int); idx2[use2] = np.arange(int(use2.sum())); A = A2[idx2[ii, jj]]; P = corr.penalty()
    fixed = dict(fixed or {}); start = dict(CROSS_START, **(start or {}))
    free = [k for k in CROSS_PARAMS if k not in fixed and not (k == 'gamma_S' and corr is None)]
    full = np.array([fixed.get(k, start[k]) for k in CROSS_PARAMS], float)
    def unpack(q):
        v = full.copy()
        for k, val in zip(free, q): v[CROSS_PARAMS.index(k)] = val
        return v
    cache = {}
    def base_model(v):
        bq, gq, dr, sig, _ = v
        key = (float(dr), float(sig))                 # exact: the finite-difference steps of the optimiser must register
        if key not in cache:
            if len(cache) > 8: cache.clear()
            cache[key] = coarse_basis(dr, sig)
        B = cache[key]; b_q = bq * x ** gq; beta_q = fz / b_q; c0 = b_q * bF * D2
        return c0 * (B['mu0'][ii, jj] + (beta_q + betaF) * B['mu2'][ii, jj] + beta_q * betaF * B['mu4'][ii, jj])
    def corr_design(v): return (x[:, None] ** v[4]) * A if corr is not None else np.zeros((len(y), 0))
    def linear_coefficients(v):
        if corr is None: return np.zeros(0)
        M = corr_design(v) * w[:, None]; scale = np.maximum(np.linalg.norm(M, axis=0), 1e-30); Ms = M / scale; b = (y - base_model(v)) * w
        if ridge > 0: return np.linalg.solve(Ms.T @ Ms + ridge * (P / np.outer(scale, scale)), Ms.T @ b) / scale
        return np.linalg.lstsq(Ms, b, rcond=None)[0] / scale
    def residual(q):
        v = unpack(q); c = linear_coefficients(v); return (base_model(v) + corr_design(v) @ c - y) * w
    q = np.array([full[CROSS_PARAMS.index(k)] for k in free]); lo = np.array([CROSS_BOUNDS[k][0] for k in free]); hi = np.array([CROSS_BOUNDS[k][1] for k in free])
    q = np.clip(q, lo + 1e-9, hi - 1e-9)
    if free:
        for _ in range(int(iterations)):
            r = least_squares(residual, q, bounds=(lo, hi), x_scale=np.array([1., 1., 1., 1., 1.][:len(q)]), diff_step=1e-4, ftol=1e-10, xtol=1e-10, gtol=1e-10); q = r.x
    v = unpack(q); c = linear_coefficients(v); model_cells = base_model(v) + corr_design(v) @ c; resid = (y - model_cells) * w; chi2 = float(np.sum(resid ** 2))
    per_z = [{'z_bin': int(k), 'cells': int(np.sum(iz == k)), 'chi2': float(np.sum(resid[iz == k] ** 2)), 'mean_z': float(np.average(z[iz == k], weights=w[iz == k] ** 2)) if np.any(iz == k) else None} for k in range(nz)]
    # ---- layered table on the fine signed grid (offset and smoothing applied, so the kernel is read at OBSERVED r_par)
    if chi_nodes is None: chi_nodes = chi_nodes_for(cfg, getattr(cfg, 'xi_z_edges', (z_ref - .3, z_ref + .6)))
    chi_nodes = np.asarray(chi_nodes, float); node_z = z_of_chi(chi_nodes)
    T = {k: transform_along_rpar(basis_proj[k].xi, ref.r_par, v[2], v[3]) for k in BASIS}; G = {k: transform_along_rpar(basis_proj[k].xi_rp, ref.r_par, v[2], v[3]) for k in BASIS}
    GRP, GRZ = np.meshgrid(ref.r_perp, ref.r_par, indexing='ij')
    if corr is not None: Dv = corr.design(GRP, GRZ) @ c; Dg = corr.design_deriv(GRP, GRZ) @ c
    xi = np.zeros((len(chi_nodes),) + ref.xi.shape); g = np.zeros_like(xi)
    for k, zk in enumerate(node_z):
        xk = float(zfactor(zk, z_ref)); b_q = v[0] * xk ** v[1]; beta_q = float(growth_rate(zk)) / b_q
        bFk = -np.sqrt(forest['b_F2']) * xk ** forest['gamma_b']; betaFk = forest['beta_F'] * xk ** forest['gamma_beta']; c0 = b_q * bFk * float(growth_ratio2(zk, z_ref))
        coef = (c0, c0 * (beta_q + betaFk), c0 * beta_q * betaFk)
        xi[k] = sum(cc * T[b] for cc, b in zip(coef, BASIS)); g[k] = sum(cc * G[b] for cc, b in zip(coef, BASIS))
        if corr is not None: xi[k] += xk ** v[4] * Dv; g[k] += xk ** v[4] * Dg
    params = dict(zip(CROSS_PARAMS, (float(t) for t in v)))
    params.update(z_ref=z_ref, chi2=chi2, cells=int(len(y)), dof=int(len(y) - len(free) - len(c)), free=free, fixed={k: float(val) for k, val in fixed.items()},
                  forest_fixed=dict(forest), beta_q_at_zref=float(growth_rate(z_ref)) / v[0], n_correction=int(len(c)), ridge=float(ridge), per_z=per_z,
                  coefficients=c.tolist(), chi_nodes=chi_nodes.tolist(), node_z=[float(t) for t in node_z], z_edges=list(map(float, getattr(cfg, 'xi_z_edges', ()))),
                  fit_range={'rperp_min': float(cfg.fit_rperp_min), 'rperp_max': float(cfg.r_perp_max), 'rpar_max': float(cfg.r_par_max)},
                  b_q_at_z=[(float(zk), float(v[0] * zfactor(zk, z_ref) ** v[1])) for zk in (2.0, 2.2, 2.4, 2.6, 2.8, 3.0)])
    if corr is not None: params['knots'] = dict(r_perp=list(corr.rp_knots), r_par=list(corr.rz_knots))
    meta = dict(provider='cross_fit+zevol' + ('+spline' if corr is not None else ''), basis=list(BASIS), accepted_weight=float(den.sum()), z_ref=z_ref, signed_rpar=True,
                note='quasar x forest; r_par = chi_pixel - chi_quasar; offset and smoothing applied; xi_rp is d/dr_perp')
    table = XiTable(ref.r_perp, ref.r_par, xi, g, meta, (num.sum(axis=0), den.sum(axis=0)), chi_nodes=chi_nodes); table.meta['fit'] = params
    return table, params, dict(use=use, model_cells=model_cells, raw=raw, z=z, iz=iz, ii=ii, jj=jj, resid=resid, rz_c=rz_c, rp_c=rp_c)


def cross_correction_for(basis_proj, forest, cfg):
    """Envelope from the reference cross shape (b_q = 3, beta_q = f/3, forest side from the auto fit)."""
    ref = basis_proj[BASIS[0]]; z_ref = float(getattr(cfg, 'xi_z_ref', 2.4)); bq = 3.; betaq = float(growth_rate(z_ref)) / bq; bF = -np.sqrt(forest['b_F2']); betaF = forest['beta_F']
    c = (bq * bF, bq * bF * (betaq + betaF), bq * bF * betaq * betaF)
    shape = XiTable(ref.r_perp, ref.r_par, -sum(cc * basis_proj[k].xi for cc, k in zip(c, BASIS)), -sum(cc * basis_proj[k].xi_rp for cc, k in zip(c, BASIS)))
    env = Envelope(shape)
    # Envelope uses hypot(rp, rz), symmetric in the sign of r_par, so it serves the signed grid
    return CrossCorrection(env, rp_knots=(cfg.fit_rperp_min, cfg.r_perp_max), rz_knots=(-cfg.r_par_max, 0., cfg.r_par_max))


# ------------------------------------------------------------------------------------------- accumulator
@njit(parallel=True, cache=True)
def _accumulate_cross_kernel(pix_start, chi, delta, weight, chiq, pa, pb, theta, nsl, rp_grid, rz_grid, xi, xirp, rpmax, rzmax, chi_ref,
                             slab, slab_edges, slab_index, rpmin, theta_true, expectation, chi0, dchi, nchi):
    n = pa.size; out = np.zeros((n, 11, 6)); counts = np.zeros(n, np.int32)
    rp0 = rp_grid[0]; rz0 = rz_grid[0]; drp = rp_grid[1] - rp_grid[0]; drz = rz_grid[1] - rz_grid[0]; nrp = rp_grid.size; nrz = rz_grid.size
    for ip in prange(n):
        cq = float(chiq[pa[ip] - nsl]); s = pb[ip]
        for p in range(pix_start[s], pix_start[s + 1]):
            cp = float(chi[p]); rz = cp - cq
            if rz < -rzmax or rz > rzmax: continue
            cm = .5 * (cp + cq); rp = cm * float(theta[ip]); dc = cq - cp
            selected = slab[p] >= 0
            if slab_index >= 0: selected = selected and slab_edges[slab_index, 0] <= cm < slab_edges[slab_index, 1]
            if rp <= rpmax and rp >= rpmin and selected:
                ib = _shape_bin(rp, abs(rz), rpmax, rzmax)
                if ib >= 0:
                    xv, xg = _interp_layer(rp, rz, cm, rp0, drp, nrp, rz0, drz, nrz, chi0, dchi, nchi, xi, xirp)
                    G = cm * xg; dm = cm - chi_ref; ww = float(weight[p])
                    if expectation:
                        dd, _ = _interp_layer(cm * float(theta_true[ip]), rz, cm, rp0, drp, nrp, rz0, drz, nrz, chi0, dchi, nchi, xi, xirp)
                    else: dd = float(delta[p])
                    out[ip, 0, ib] += ww * dd * G; out[ip, 1, ib] += ww * dd * G * dm; out[ip, 2, ib] += ww * dd * G * dc * .5
                    gg = ww * G * G; out[ip, 3, ib] += gg; out[ip, 4, ib] += gg * dm; out[ip, 5, ib] += gg * dm * dm
                    out[ip, 6, ib] += gg * dc; out[ip, 7, ib] += gg * dc * dc * .25
                    bx = ww * xv * G; out[ip, 8, ib] += bx; out[ip, 9, ib] += bx * dm; out[ip, 10, ib] += bx * dc * .5
                    counts[ip] += 1
    return out, counts


def accumulate_cross(sl, qso, pairs, table: XiTable, cfg, shifted_positions=None, true_positions=None):
    """Quasar-pixel pair sums with the eleven accumulators of `pairs.accumulate`; ``a`` indexes the quasar as
    N_sl + q on the concatenated (sightline, quasar) position list. ``shifted_positions``/``true_positions`` are
    (ra, dec) arrays over that concatenated list (injection tests)."""
    a, b, thx, thy, theta = pairs
    if shifted_positions is not None:
        pos = np.asarray(shifted_positions); thx, thy, theta = pair_geometry(pos[:, 0], pos[:, 1], a, b)
    expectation = true_positions is not None
    theta_true = pair_geometry(np.asarray(true_positions)[:, 0], np.asarray(true_positions)[:, 1], a, b)[2] if expectation else theta
    chi0, dchi, nchi = table.layers()
    edges = np.asarray([[float(chi_of_z(z)) for z in bounds] for bounds in cfg.slabs])
    out, n = _accumulate_cross_kernel(sl.pix_start, sl.chi, sl.delta, sl.w, np.asarray(qso.chi, np.float32), a, b, theta, sl.nq, table.r_perp, table.r_par,
                                      table.xi.ravel(), table.xi_rp.ravel(), cfg.r_perp_max, cfg.r_par_max, cfg.chi_ref, sl.slab, edges, cfg.slab_index,
                                      float(getattr(cfg, 'r_perp_min', 0.)), theta_true, expectation, chi0, dchi, nchi)
    keep = n > 0
    attrs = {'chi_ref': cfg.chi_ref, 'kind': 'quasar x forest', 'pair_direction': 'theta_quasar - theta_sightline', 'quasar_index_offset': int(sl.nq),
             'r_perp_min': float(getattr(cfg, 'r_perp_min', 0.)), 'expectation': expectation, 'xi_layers': int(nchi)}
    return PairCatalogue(a[keep], b[keep], thx[keep], thy[keep], theta[keep], out[keep], n[keep], attrs)


class Positions:
    """(ra, dec) container over the concatenated sightline + quasar list, for the jackknife regions."""
    def __init__(self, sl, qso): self.ra = np.r_[sl.ra, qso.ra]; self.dec = np.r_[sl.dec, qso.dec]
