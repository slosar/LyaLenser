import numpy as np
import pytest
from lyalenser.config import Config
from lyalenser.xi_model import XiTable
from lyalenser.xi_fit import BasisPower, basis_tables, project_fine, coarse_bin, fit_model_table, BASIS


def _hankel_basis(cfg, nk=4000):
    # r k reaches ~1000, so the log-k quadrature needs a few thousand points for percent-level tables.
    return basis_tables(cfg, 'hankel', nk=nk)


def test_hankel_basis_analytic_derivative_matches_finite_difference():
    cfg = Config(xi_max=32, xi_step=.25)
    b = _hankel_basis(cfg)['mu0']
    fd = np.gradient(b.xi, cfg.xi_step, axis=0)
    use = (b.r_perp > 3) & (b.r_perp < 28)
    assert np.linalg.norm((fd - b.xi_rp)[use]) / np.linalg.norm(b.xi_rp[use]) < .01


def test_kaiser_basis_recombines_to_forest_power():
    from lyalenser.forest_power import ForestPower
    pf = ForestPower(model='kaiser'); b2, beta = pf.p['b_F']**2, pf.p['beta_F']
    kpar = np.array([.05, .3, 1.]); kperp = np.array([.1, .5, 2.])
    parts = [BasisPower(i, model='kaiser')(kpar[:, None], kperp[None, :]) for i in range(3)]
    combined = b2 * (parts[0] + 2*beta*parts[1] + beta*beta*parts[2])
    # The basis tabulates P_lin exp(-(k/kp)^2) on a log grid: agreement to the interpolation accuracy.
    assert np.allclose(combined, pf(kpar[:, None], kperp[None, :]), rtol=1e-3)


def test_projection_commutes_with_transverse_derivative_and_kills_constants():
    cfg = Config(xi_max=32, xi_step=.25)
    b = _hankel_basis(cfg)['mu0']
    cpix = np.arange(3700., 4400., .55)
    p = project_fine(b, cpix, cfg)
    fd = np.gradient(p.xi, cfg.xi_step, axis=0)
    use = (p.r_perp > 3) & (p.r_perp < 28)
    assert np.linalg.norm((fd - p.xi_rp)[use]) / np.linalg.norm(p.xi_rp[use]) < .01
    # A table constant along r_par is removed entirely by the mean+slope projection when every pixel pair of the
    # forest lies inside the table range (forest shorter than xi_max); beyond that range the table is zero.
    const = XiTable(b.r_perp, b.r_par, np.ones_like(b.xi), np.zeros_like(b.xi))
    q = project_fine(const, np.arange(3700., 3728., .55), cfg)
    assert np.max(abs(q.xi)) < 1e-8


def test_fit_recovers_amplitude_and_beta_from_synthetic_table():
    cfg = Config(xi_max=32, xi_step=.5, fit_rperp_min=3.)
    basis = _hankel_basis(cfg, nk=2000)
    cpix = np.arange(3700., 4400., .55)
    proj = {k: project_fine(basis[k], cpix, cfg) for k in BASIS}
    coarse = {k: coarse_bin(proj[k], cfg) for k in BASIS}
    b2, beta = .02, 1.4
    truth = b2 * (coarse['mu0'] + 2*beta*coarse['mu2'] + beta*beta*coarse['mu4'])
    rng = np.random.default_rng(3); den = np.full(truth.shape, 1e6)
    # Small noise: (b_F^2, beta_F) are degenerate at the few-percent level on a single noisy table.
    num = (truth + 3e-5 * rng.normal(size=truth.shape)) * den
    fit = fit_model_table(num, den, proj, cfg, coarse)
    assert fit.params['b_F2'] == pytest.approx(b2, rel=.01)
    assert fit.params['beta_F'] == pytest.approx(beta, abs=.03)
    assert fit.table.xi.shape == proj['mu0'].xi.shape and np.all(np.isfinite(fit.table.xi_rp))


def test_pair_kernel_respects_r_perp_min():
    from lyalenser.pairs import find_pairs, accumulate
    from test_pairs import sample, table
    sl = sample(); t = table(); p = find_pairs(sl, 30/sl.chi.min())
    full = accumulate(sl, p, t, Config(chi_ref=3900))
    cut = accumulate(sl, p, t, Config(chi_ref=3900, r_perp_min=10.))
    # shape bins 0 and 1 are r_perp < 10: nothing may remain there; bins >= 2 are unchanged pair by pair.
    assert np.all(cut.accum[:, :, :2] == 0)
    index = {(int(a), int(b)): i for i, (a, b) in enumerate(zip(full.a, full.b))}
    rows = [index[(int(a), int(b))] for a, b in zip(cut.a, cut.b)]
    assert np.allclose(full.accum[rows][:, :, 2:], cut.accum[:, :, 2:])
    assert len(cut.a) < len(full.a)


