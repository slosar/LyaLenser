"""Iteration 9: the Lyb region (B) as an extension of the Lya region (A) on the same sightline.

A region-B pixel carries Lya absorption from our slab and Lyb absorption from a much more distant one. For pixel
pairs inside 30 Mpc/h the Lya-Lyb cross term joins two fields ~1000 Mpc/h apart and is negligible, so A x B pairs
measure Lya-Lya alone; B x B pairs also carry Lyb-Lyb at a common redshift and are dropped everywhere: in the
pair accumulator, in the measured correlation and in the continuum-projection average.
"""
import sys
from pathlib import Path
import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
for p in (HERE.parent, HERE.parent.parent):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))
from config import Config, SightlineSet
from xi_model import XiTable, xi_from_data
from pairs import find_pairs, accumulate
from xi_fit import _weighted_projector, basis_tables, project_fine_sample


def _sightlines(region_of_second_half, n_per=24, nq=6, seed=3):
    """Two-column geometry: ``nq`` sightlines a few Mpc/h apart, each with 2 x n_per pixels."""
    rng = np.random.default_rng(seed)
    chi0 = 3700.0
    chi = np.tile(chi0 + np.arange(2 * n_per) * 1.0, nq)
    reg = np.tile(np.r_[np.zeros(n_per, np.int8), np.full(n_per, region_of_second_half, np.int8)], nq)
    npix = 2 * n_per
    ra = 180.0 + np.arange(nq) * 0.06
    dec = np.full(nq, 30.0)
    return SightlineSet(np.arange(nq), ra, dec, np.full(nq, 3.0, np.float32),
                        np.arange(nq + 1) * npix, chi,
                        rng.normal(size=nq * npix).astype(np.float32),
                        np.ones(nq * npix, np.float32), np.zeros(nq * npix, np.int8),
                        {}, reg)


@pytest.fixture(scope='module')
def cfg():
    return Config(xi_max=20.0, xi_step=0.5, r_perp_max=20.0, r_par_max=20.0, r_perp_min=0.0,
                  fit_rperp_min=3.0, chi_ref=3700.0)


@pytest.fixture(scope='module')
def table(cfg):
    g = np.arange(0, cfg.xi_max + .25, cfg.xi_step)
    x = np.exp(-np.hypot(g[:, None], g[None, :]) / 8)
    return XiTable(g, g, x, np.gradient(x, cfg.xi_step, axis=0))


# ------------------------------------------------------------------ pair accumulator
def test_bb_pairs_are_dropped(cfg, table):
    """Tagging half of every sightline as region B removes exactly the B x B pixel pairs."""
    allA = _sightlines(0)
    mixed = _sightlines(1)
    assert np.array_equal(allA.chi, mixed.chi) and np.array_equal(allA.delta, mixed.delta)
    pairs = find_pairs(allA, cfg.r_perp_max / float(allA.chi.min()))
    n_all = accumulate(allA, pairs, table, cfg).npair.sum()
    n_mix = accumulate(mixed, pairs, table, cfg).npair.sum()
    assert n_mix < n_all
    # count, directly, the B x B pixel pairs that pass the separation cuts
    bb = 0
    for a, b, th in zip(pairs[0], pairs[1], pairs[4]):
        for p in range(allA.pix_start[a], allA.pix_start[a + 1]):
            for q in range(allA.pix_start[b], allA.pix_start[b + 1]):
                if mixed.region[p] == 0 or mixed.region[q] == 0:
                    continue
                cp, cq = float(allA.chi[p]), float(allA.chi[q])
                rz = abs(cp - cq); cm = .5 * (cp + cq); rp = cm * float(th)
                if rz <= cfg.r_par_max and cfg.r_perp_min <= rp <= cfg.r_perp_max and rp <= 30 and rz <= 30:
                    bb += 1
    assert bb > 0                       # the geometry really does contain B x B pairs
    assert n_all - n_mix == bb          # and exactly those are the ones that disappear


def test_region_tag_defaults_to_A():
    sl = _sightlines(0)
    bare = SightlineSet(sl.qid, sl.ra, sl.dec, sl.zq, sl.pix_start, sl.chi, sl.delta, sl.w, sl.slab)
    assert bare.region.shape == sl.chi.shape and not bare.region.any()


# ------------------------------------------------------------------ measured correlation
def test_measured_counts_split_by_pair_type(cfg):
    mixed = _sightlines(1)
    t = xi_from_data(mixed, cfg)
    aa, ab = t.counts_by_type['AA'], t.counts_by_type['AB']
    assert aa[1].sum() > 0 and ab[1].sum() > 0
    assert np.allclose(t.counts[1], aa[1] + ab[1])          # the fitted table uses A x A + A x B
    allA = _sightlines(0)
    t0 = xi_from_data(allA, cfg)
    assert t0.counts[1].sum() > t.counts[1].sum()           # B x B weight is missing from the mixed run
    assert np.allclose(t0.counts_by_type['AB'][1], 0.0)     # an all-A sample has no A x B pairs


# ------------------------------------------------------------------ continuum projection
def test_block_projector_removes_a_mean_and_slope_per_region():
    c = 3700.0 + np.arange(40) * 1.0
    reg = np.r_[np.zeros(20, np.int8), np.ones(20, np.int8)]
    w = np.ones(40)
    u, v = _weighted_projector(c, w, reg)
    assert u.shape == (40, 4)
    P = lambda x: x - u @ (v.T @ x)
    for a0, b0, a1, b1 in ((1., 0., 0., 0.), (0., 2e-3, 0., 0.), (3., -1e-3, -2., 4e-3)):
        x = np.where(reg == 0, a0 + b0 * (c - c.mean()), a1 + b1 * (c - c.mean()))
        assert np.max(np.abs(P(x))) < 1e-9                  # anything linear per block is removed
    x = np.sin(c / 7.0)
    assert np.max(np.abs(P(x))) > 1e-3                      # but not everything


def test_single_region_projector_unchanged():
    c = 3700.0 + np.arange(30) * 1.0
    w = np.linspace(1, 2, 30)
    u0, v0 = _weighted_projector(c, w)
    u1, v1 = _weighted_projector(c, w, np.zeros(30, np.int8))
    assert np.allclose(u0, u1) and np.allclose(v0, v1)


def test_projection_average_drops_bb_pairs(cfg):
    """A two-region sightline projects to less correlation than the same pixels treated as one region."""
    raw = basis_tables(cfg, 'hankel', nk=400, los_pixel=0.55, los_resolution=0.4)
    c = 3700.0 + np.arange(120) * 0.55
    w = np.ones(120)
    reg = np.r_[np.zeros(60, np.int8), np.ones(60, np.int8)]
    rp = np.array([0.0, 4.0, 8.0])
    kw = dict(cfg=cfg, n_pairs=4, seed=1, rp_grid=rp, step=0.55)
    one = project_fine_sample({'mu0': raw['mu0']}, [(c, w, np.ones(120), np.zeros(120, np.int8))] * 3, **kw)['mu0']
    two = project_fine_sample({'mu0': raw['mu0']}, [(c, w, np.ones(120), reg)] * 3, **kw)['mu0']
    j = np.searchsorted(one.r_par, 2.0)
    assert not np.allclose(one.xi[1:, j], two.xi[1:, j])     # two continuum blocks, and no B x B pairs
    assert np.all(np.isfinite(two.xi))
