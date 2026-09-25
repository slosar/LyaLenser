"""Iteration 8: the pair-weighted continuum projection and the spline-corrected response table."""
import numpy as np
import pytest

from lyalenser.config import Config
from lyalenser.xi_model import XiTable
from lyalenser.xi_fit import basis_tables, project_fine, project_fine_sample, refine_rp_grid, coarse_bin, BASIS
from lyalenser.xi_spline import Basis1D, Envelope, XiCorrection, fit_corrected_table


@pytest.fixture(scope='module')
def cfg():
    return Config(xi_max=20.0, xi_step=0.5, analytic_nk=600, fit_rperp_min=3.0, r_perp_max=15.0, r_par_max=15.0)


@pytest.fixture(scope='module')
def raw(cfg):
    return basis_tables(cfg, 'hankel', nk=cfg.analytic_nk, los_pixel=0.55, los_resolution=0.4)


# ------------------------------------------------------------------ projection
def test_pair_projection_reproduces_single_forest(cfg, raw):
    """Identical, fully-kept, uniformly-weighted forests must give exactly `project_fine`."""
    c = 3700 + np.arange(200) * 0.55
    one = project_fine(raw['mu0'], c, cfg)
    rp = np.arange(0, cfg.xi_max + 0.25, 1.0)
    many = project_fine_sample({'mu0': raw['mu0']}, [(c, np.ones(len(c)))] * 3, cfg,
                               n_pairs=4, seed=1, rp_grid=rp)['mu0']
    i = np.searchsorted(one.r_perp, rp)
    for a, b in ((one.xi[i], many.xi), (one.xi_rp[i], many.xi_rp)):
        m = np.abs(a) > 1e-6 * np.max(np.abs(a))
        assert np.max(np.abs(a[m] / b[m] - 1)) < 2e-3


def test_shorter_forests_suppress_more(cfg, raw):
    """A shorter continuum-fitting span removes more of the correlation."""
    long_c = 3700 + np.arange(400) * 0.55
    short_c = 3700 + np.arange(150) * 0.55
    rp = np.array([0.0, 4.0, 8.0, 12.0])
    kw = dict(cfg=cfg, n_pairs=4, seed=1, rp_grid=rp, step=0.55)
    a = project_fine_sample({'mu0': raw['mu0']}, [(long_c, np.ones(len(long_c)))] * 3, **kw)['mu0']
    b = project_fine_sample({'mu0': raw['mu0']}, [(short_c, np.ones(len(short_c)))] * 3, **kw)['mu0']
    j = np.searchsorted(a.r_par, 2.0)
    assert np.all(b.xi[1:, j] < a.xi[1:, j])


def test_in_range_mask_changes_pairs_not_projector(cfg, raw):
    """Masking pixels out of the pair sums must not change the projector (the continuum saw them)."""
    c = 3700 + np.arange(300) * 0.55
    keep = np.zeros(len(c), bool); keep[50:250] = True
    rp = np.array([0.0, 6.0, 12.0])
    kw = dict(cfg=cfg, n_pairs=4, seed=1, rp_grid=rp, step=0.55)
    full = project_fine_sample({'mu0': raw['mu0']}, [(c, np.ones(len(c)))] * 3, **kw)['mu0']
    cut = project_fine_sample({'mu0': raw['mu0']}, [(c, np.ones(len(c)), keep)] * 3, **kw)['mu0']
    shortened = project_fine_sample({'mu0': raw['mu0']}, [(c[keep], np.ones(keep.sum()))] * 3, **kw)['mu0']
    j = np.searchsorted(full.r_par, 2.0)
    # the masked version keeps the long projector: closer to `full` than the truly shortened forest is
    assert abs(cut.xi[1, j] - full.xi[1, j]) < abs(shortened.xi[1, j] - full.xi[1, j])


def test_refine_rp_grid_preserves_values(cfg, raw):
    rp = np.arange(0, cfg.xi_max + 0.25, 1.0)
    t = {'mu0': XiTable(rp, raw['mu0'].r_par, raw['mu0'].xi[::2][:len(rp)], raw['mu0'].xi_rp[::2][:len(rp)])}
    out = refine_rp_grid(t, cfg)['mu0']
    i = np.searchsorted(out.r_perp, rp)
    assert np.allclose(out.xi[i], t['mu0'].xi, rtol=1e-10, atol=1e-14)


# ------------------------------------------------------------------ spline basis
def test_basis_partition_of_unity():
    b = Basis1D((3.0, 6.0, 10.0, 30.0))
    x = np.linspace(3.0, 30.0, 51)
    assert np.allclose(b(x).sum(axis=-1), 1.0, atol=1e-10)


def test_basis_clamped_outside():
    b = Basis1D((3.0, 10.0, 30.0))
    assert np.allclose(b(np.array([1.0])), b(np.array([3.0])))
    assert np.allclose(b(np.array([40.0])), b(np.array([30.0])))
    assert np.allclose(b(np.array([1.0, 40.0]), deriv=1), 0.0)


def test_basis_derivative_matches_finite_difference():
    b = Basis1D((3.0, 6.0, 12.0, 30.0))
    x = np.array([4.0, 8.0, 20.0]); h = 1e-5
    fd = (b(x + h) - b(x - h)) / (2 * h)
    assert np.allclose(b(x, deriv=1), fd, atol=1e-6)


# ------------------------------------------------------------------ correction design
@pytest.fixture(scope='module')
def corr(cfg, raw):
    proj = {k: project_fine(raw[k], 3700 + np.arange(200) * 0.55, cfg) for k in BASIS}
    ref = XiTable(proj['mu0'].r_perp, proj['mu0'].r_par,
                  sum(c * proj[k].xi for c, k in zip((0.021, 0.06, 0.043), BASIS)),
                  sum(c * proj[k].xi_rp for c, k in zip((0.021, 0.06, 0.043), BASIS)))
    return proj, ref, XiCorrection(Envelope(ref), rp_knots=(3., 6., 10., 15.), rz_knots=(0., 4., 15.),
                                   sw_knots=(3., 8., 15.))


def test_envelope_positive_and_decreasing(corr):
    _, ref, c = corr
    r = np.linspace(1.0, 14.0, 40)
    s = c.envelope(r, np.zeros_like(r))
    assert np.all(s > 0)
    assert np.all(np.diff(s) <= 1e-18)


def test_correction_is_even_in_rpar(corr):
    _, _, c = corr
    rp = np.array([4.0, 9.0]); rz = np.array([3.0, 3.0])
    assert np.allclose(c.design(rp, rz), c.design(rp, -rz))
    h = 1e-4                                            # zero slope in r_par at r_par = 0
    d = (c.design(rp, np.full(2, h))[:, :c.n_s] - c.design(rp, np.zeros(2))[:, :c.n_s]) / h
    assert np.max(np.abs(d)) < 1e-3 * np.max(np.abs(c.design(rp, np.zeros(2))[:, :c.n_s]))


def test_design_derivative_matches_finite_difference(corr):
    _, _, c = corr
    rp = np.array([4.0, 7.0, 12.0]); rz = np.array([2.0, 5.0, 8.0]); h = 1e-4
    fd = (c.design(rp + h, rz) - c.design(rp - h, rz)) / (2 * h)
    an = c.design_deriv(rp, rz)
    assert np.allclose(an[:, :c.n_s], fd[:, :c.n_s], rtol=1e-3, atol=1e-8)
    assert np.allclose(an[:, c.n_s:], 0.0)              # the same-wavelength term is not lensed


# ------------------------------------------------------------------ fitting
def _synthetic(cfg, proj, corr, coef, b2=0.021, beta=1.4):
    """Noise-free measured cells from a known (base + correction) table."""
    from scipy.interpolate import RegularGridInterpolator
    n = int(cfg.xi_max)
    c = np.array([b2, 2 * b2 * beta, b2 * beta * beta])
    fine = XiTable(proj['mu0'].r_perp, proj['mu0'].r_par,
                   sum(x * proj[k].xi for x, k in zip(c, BASIS)),
                   sum(x * proj[k].xi_rp for x, k in zip(c, BASIS)))
    f = RegularGridInterpolator((fine.r_perp, fine.r_par), fine.xi, bounds_error=False, fill_value=0.)
    o = (np.arange(4) + 0.5) / 4
    RP = np.repeat((np.arange(n)[:, None, None] + o[None, :, None]), 4, axis=2)
    RZ = np.repeat((np.arange(n)[None, :, None] + o[None, None, :])[:, :, :], 1, axis=0)
    RP = np.broadcast_to(RP[:, None, :, :], (n, n, 4, 4))
    RZ = np.broadcast_to((np.arange(n)[None, :, None, None] + o[None, None, None, :]), (n, n, 4, 4))
    base = f((RP, RZ)).mean(axis=(2, 3))
    extra = (corr.design(RP, RZ) @ coef).mean(axis=(2, 3))
    return base + extra


def test_fit_recovers_an_injected_correction(cfg, corr):
    proj, ref, c = corr
    rng = np.random.default_rng(0)
    coef = np.r_[rng.normal(0, 0.05, c.n_s), rng.normal(0, 1e-3, c.bsw.n)]
    num_over_den = _synthetic(cfg, proj, c, coef)
    n = num_over_den.shape[0]
    den = np.full((n, n), 1e9)
    coarse = {k: coarse_bin(proj[k], cfg) for k in BASIS}
    tab, par, dbg = fit_corrected_table(num_over_den * den, den, proj, coarse, cfg, c, ridge=0.0)
    truth = _synthetic(cfg, proj, c, coef)
    m = dbg['use']
    assert np.max(np.abs(dbg['model_cells'] - truth[m])) < 1e-3 * np.max(np.abs(truth[m]))


def test_same_wavelength_term_is_in_xi_but_not_in_the_kernel(corr):
    """A coefficient vector carrying only the same-wavelength block must move xi inside the first radial bin,
    leave xi untouched outside it, and leave the kernel untouched everywhere."""
    _, _, c = corr
    coef = np.r_[np.zeros(c.n_s), np.full(c.bsw.n, 1e-3)]
    rp = np.array([4.0, 8.0, 12.0])
    inside = c.design(rp, np.full(3, 0.4)) @ coef
    outside = c.design(rp, np.full(3, 3.0)) @ coef
    kernel = c.design_deriv(rp, np.full(3, 0.4)) @ coef
    assert np.all(inside > 5e-4)
    assert np.allclose(outside, 0.0)
    assert np.allclose(kernel, 0.0)
