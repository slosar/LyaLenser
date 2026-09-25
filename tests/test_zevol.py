"""Iteration 10: the redshift-evolving correlation table (layered in chi) and its fit."""
import numpy as np
import pytest

from lyalenser.config import Config
from lyalenser.xi_model import XiTable, xi_from_data, layered_bilinear
from lyalenser.xi_fit import basis_tables, project_fine, coarse_bin, BASIS
from lyalenser.xi_spline import Envelope, XiCorrection
from lyalenser.xi_zevol import fit_evolving_table, base_coefficients, zfactor, chi_nodes_for
from lyalenser.pairs import find_pairs, accumulate
from lyalenser.cosmo import chi as chi_of_z, z_of_chi
from test_pairs import sample, table, brute


def _layered(t, nodes, scale):
    """Layered copy of a 2-D table whose layer k is scale[k] times the table."""
    return XiTable(t.r_perp, t.r_par, np.stack([s * t.xi for s in scale]), np.stack([s * t.xi_rp for s in scale]),
                   {}, chi_nodes=nodes)


def test_identical_layers_reproduce_the_2d_table():
    sl = sample(); cfg = Config(chi_ref=3900); p = find_pairs(sl, 30 / sl.chi.min()); t = table()
    ref = accumulate(sl, p, t, cfg).accum
    lay = _layered(t, 3800 + 50. * np.arange(5), np.ones(5))
    assert np.allclose(accumulate(sl, p, lay, cfg).accum, ref, rtol=1e-12, atol=1e-12)


def test_layer_interpolation_matches_bruteforce():
    """Layers scaled linearly in chi: the kernel must see s(chi_mean) x table, with s interpolated between nodes
    and clamped outside them."""
    sl = sample(); cfg = Config(chi_ref=3900); p = find_pairs(sl, 30 / sl.chi.min()); t = table()
    nodes = 3905 + 10. * np.arange(4); scale = np.array([1.0, 1.3, 0.8, 1.1])
    lay = _layered(t, nodes, scale)
    got = accumulate(sl, p, lay, cfg)
    # brute force with an explicit scale at the pair mean distance
    a, b, tx, ty, th = p; out = np.zeros((len(a), 11, 6)); count = np.zeros(len(a), int)
    for ip, (aa, bb) in enumerate(zip(a, b)):
        for i in range(sl.pix_start[aa], sl.pix_start[aa + 1]):
            for j in range(sl.pix_start[bb], sl.pix_start[bb + 1]):
                dc = float(sl.chi[i]) - float(sl.chi[j]); rz = abs(dc); cm = .5 * (float(sl.chi[i]) + float(sl.chi[j])); rp = cm * th[ip]
                if rz > 30 or rp > 30: continue
                s = np.interp(cm, nodes, scale)
                ir = 0 if rp < 10 else (1 if rp < 20 else 2); iz = 0 if rz < 10 else 1; ib = 2 * ir + iz
                xv, xg = t.interp(rp, rz); xv *= s; xg *= s
                G = cm * xg; dm = cm - cfg.chi_ref; ww = float(sl.w[i]) * float(sl.w[j]); dd = float(sl.delta[i]) * float(sl.delta[j])
                vals = (ww * dd * G, ww * dd * G * dm, ww * dd * G * dc / 2, ww * G * G, ww * G * G * dm, ww * G * G * dm * dm,
                        ww * G * G * dc, ww * G * G * dc * dc / 4, ww * xv * G, ww * xv * G * dm, ww * xv * G * dc / 2)
                out[ip, :, ib] += vals; count[ip] += 1
    keep = count > 0
    assert np.allclose(got.accum, out[keep], rtol=1e-9, atol=1e-11)
    v, g = lay.interp(5., 2., 3960.)
    v0, g0 = t.interp(5., 2.)
    assert np.isclose(v, 1.1 * v0) and np.isclose(g, 1.1 * g0)      # clamped above the last node
    v, g = lay.interp(5., 2., 3910.)
    assert np.isclose(v, 1.15 * v0)


def test_at_chi_and_save_roundtrip(tmp_path):
    t = table(); lay = _layered(t, 3900 + 20. * np.arange(3), np.array([1., 2., 3.]))
    assert np.allclose(lay.at_chi(3910.).xi, 1.5 * t.xi)
    lay.save(tmp_path / 'x.h5', 'xi')
    from lyalenser.tables import read_xi
    back = read_xi(tmp_path / 'x.h5', 'xi')
    assert back.layered and np.array_equal(back.chi_nodes, lay.chi_nodes) and np.array_equal(back.xi, lay.xi)


def test_zbinned_cells_sum_to_the_unbinned_ones():
    sl = sample(nq=30, npix=80); z = z_of_chi(sl.chi)
    edges = (float(z.min()) - 1e-3, 2.62, 2.64, float(z.max()) + 1e-3)
    cfg0 = Config(chi_ref=3900, xi_max=20., r_perp_max=20., r_par_max=20.)
    cfg1 = cfg0.copy(xi_z_evolution=True, xi_z_edges=edges)
    t0 = xi_from_data(sl, cfg0); t1 = xi_from_data(sl, cfg1)
    n1, d1, c1 = t1.counts_z['all']
    assert n1.shape[0] == 3
    assert np.allclose(n1.sum(axis=0), t0.counts[0]) and np.allclose(d1.sum(axis=0), t0.counts[1])
    assert np.allclose(t1.counts[0], t0.counts[0])
    ok = d1 > 0
    zc = z_of_chi(c1[ok] / d1[ok]); iz = np.nonzero(ok)[0]
    lo = np.asarray(edges)[iz]; hi = np.asarray(edges)[iz + 1]
    assert np.all(zc >= lo - 1e-6) and np.all(zc <= hi + 1e-6)


@pytest.fixture(scope='module')
def machinery():
    cfg = Config(xi_max=20.0, xi_step=0.5, analytic_nk=600, fit_rperp_min=3.0, r_perp_max=15.0, r_par_max=15.0,
                 xi_z_evolution=True, xi_z_edges=(2.1, 2.3, 2.5, 2.7, 3.0), xi_z_ref=2.4)
    raw = basis_tables(cfg, 'hankel', nk=cfg.analytic_nk, los_pixel=0.55, los_resolution=0.4)
    proj = {k: project_fine(raw[k], 3700 + np.arange(200) * 0.55, cfg) for k in BASIS}
    coarse = {k: coarse_bin(proj[k], cfg) for k in BASIS}
    ref = XiTable(proj['mu0'].r_perp, proj['mu0'].r_par,
                  sum(c * proj[k].xi for c, k in zip((0.021, 0.06, 0.043), BASIS)),
                  sum(c * proj[k].xi_rp for c, k in zip((0.021, 0.06, 0.043), BASIS)))
    corr = XiCorrection(Envelope(ref), rp_knots=(3., 15.), rz_knots=(0., 15.), sw_knots=(3., 8., 15.))
    return cfg, proj, coarse, corr


def _synthetic_cells(cfg, coarse, corr, truth, coef, zs):
    """Cells of the evolving model at the bin redshifts ``zs`` (exact on the coarse basis and the cell-averaged
    correction), with a huge weight so the fit is noise-free."""
    from lyalenser.xi_spline import _cell_points
    n = coarse['mu0'].shape[0]; use2 = np.ones((n, n), bool)
    RP, RZ = _cell_points(use2, 4); A = corr.design(RP, RZ).mean(axis=(1, 2)).reshape(n, n, -1)
    num = []; den = []; cs = []
    for z in zs:
        c0, c1, c2 = base_coefficients(truth['b_F2'], truth['beta_F'], truth['gamma_b'], truth['gamma_beta'], z, cfg.xi_z_ref)
        x = zfactor(z, cfg.xi_z_ref)
        cS = np.r_[coef[:corr.n_s], np.zeros(corr.n_sw)]; csw = np.r_[np.zeros(corr.n_s), coef[corr.n_s:]]
        cells = c0 * coarse['mu0'] + c1 * coarse['mu2'] + c2 * coarse['mu4'] + x ** truth['gamma_S'] * (A @ cS) + A @ csw
        d = np.full((n, n), 1e9); num.append(cells * d); den.append(d); cs.append(d * float(chi_of_z(z)))
    return np.array(num), np.array(den), np.array(cs)


def test_fit_recovers_injected_exponents(machinery):
    cfg, proj, coarse, corr = machinery
    rng = np.random.default_rng(1)
    truth = dict(b_F2=0.025, beta_F=1.3, gamma_b=2.6, gamma_beta=-0.8, gamma_S=3.1)
    coef = np.r_[rng.normal(0, 0.05, corr.n_s), rng.normal(0, 1e-3, corr.bsw.n)]
    zs = [2.2, 2.4, 2.6, 2.85]
    counts = _synthetic_cells(cfg, coarse, corr, truth, coef, zs)
    tab, par, dbg = fit_evolving_table(counts, proj, coarse, cfg, corr, ridge=0.0, iterations=6)
    for k, v in truth.items():
        assert abs(par[k] - v) < (2e-3 * abs(v) + 2e-3), (k, par[k], v)
    assert tab.layered and tab.xi.shape[0] == len(par['chi_nodes'])
    # every layer is the injected model at the node redshift (fine grid, correction included)
    GRP, GRZ = np.meshgrid(tab.r_perp, tab.r_par, indexing='ij')
    Dv = corr.design(GRP, GRZ); cS = np.r_[coef[:corr.n_s], np.zeros(corr.n_sw)]; csw = np.r_[np.zeros(corr.n_s), coef[corr.n_s:]]
    for k in (0, len(par['node_z']) // 2, -1):
        zk = par['node_z'][k]
        c0, c1, c2 = base_coefficients(truth['b_F2'], truth['beta_F'], truth['gamma_b'], truth['gamma_beta'], zk, cfg.xi_z_ref)
        want = c0 * proj['mu0'].xi + c1 * proj['mu2'].xi + c2 * proj['mu4'].xi + zfactor(zk, cfg.xi_z_ref) ** truth['gamma_S'] * (Dv @ cS) + Dv @ csw
        assert np.max(np.abs(tab.xi[k] - want)) < 2e-3 * np.max(np.abs(want))


def test_fixed_exponents_and_base_only(machinery):
    cfg, proj, coarse, corr = machinery
    truth = dict(b_F2=0.025, beta_F=1.3, gamma_b=2.6, gamma_beta=-0.8, gamma_S=0.0)
    coef = np.zeros(corr.n_par)
    counts = _synthetic_cells(cfg, coarse, corr, truth, coef, [2.2, 2.4, 2.6, 2.85])
    tab, par, _ = fit_evolving_table(counts, proj, coarse, cfg, None, fixed={'gamma_b': 2.6, 'gamma_beta': -0.8})
    assert par['free'] == ['b_F2', 'beta_F', 'gamma_S'] or set(par['free']) >= {'b_F2', 'beta_F'}
    assert abs(par['b_F2'] - 0.025) < 1e-5 and abs(par['beta_F'] - 1.3) < 1e-4
    assert par['n_correction'] == 0 and par['chi2'] < 1e-6 * 1e9


def test_chi_nodes_span_the_edges():
    cfg = Config(xi_z_edges=(2.1, 2.5, 3.0), xi_z_ref=2.4, xi_z_node_step=0.05)
    nodes = chi_nodes_for(cfg, cfg.xi_z_edges)
    assert nodes[0] <= chi_of_z(2.1) and nodes[-1] >= chi_of_z(3.0)
    assert np.allclose(np.diff(nodes), nodes[1] - nodes[0])
