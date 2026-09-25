"""The tracer-slice helpers of `lyalenser.lowz`: bias-fit band, slice bookkeeping, optimal combination."""
import numpy as np
import pytest
from lyalenser.lowz import bias_band, KMAX_BIAS, slices_of, optimal_combination, Tracer
from lyalenser.cosmo import chi


def test_bias_band_is_kmax_times_chi():
    lo, hi = bias_band(.4, .6)
    assert lo == 40. and abs(hi - (KMAX_BIAS * float(chi(.5)) - .5)) < 1e-6
    assert bias_band(1.1, 1.6)[1] > bias_band(.4, .6)[1]


def test_slices_are_sorted_and_disjoint():
    tr = (Tracer('A', .4, .6, 2., 100.), Tracer('B', .1, .4, 1.5, 300.), Tracer('C', .4, .6, 1.2, 200.))
    sl = slices_of(tr)
    assert [(s['zmin'], s['zmax']) for s in sl] == [(.1, .4), (.4, .6)]
    assert [t.name for t in sl[1]['tracers']] == ['A', 'C']
    with pytest.raises(ValueError):
        slices_of(tr + (Tracer('D', .5, .8, 2., 10.),))


def test_optimal_combination_of_correlated_slices():
    rng = np.random.default_rng(2); k = 4; nr = 40
    common = rng.normal(size=nr); jk = np.array([1 + .05 * common + .03 * rng.normal(size=nr) for _ in range(k)])
    A = jk.mean(axis=1)
    out = optimal_combination(A, jk)
    # the optimal error never exceeds the best single slice (up to the Hartlap correction of the inverse)
    assert np.isfinite(out['A']) and out['error'] <= np.sqrt(np.diag(np.array(out['covariance'])) / out['hartlap']).min() * (1 + 1e-9)
    assert abs(out['A'] - 1) < .5 and 0 < out['hartlap'] < 1 and out['naive_error'] > 0
