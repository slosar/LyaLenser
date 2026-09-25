"""NaMaster conventions used by Stage B: spin-1 deflection E-mode, decoupled theory, amplitude fit."""
import numpy as np
import healpy as hp
import pytest

pytest.importorskip('pymaster')
from lyalenser.nmt_spectra import Spectra, deflection_maps, fit_amplitude


def test_gradient_template_cross_matches_2_over_sqrt_l_l1():
    nside = 128; lmax = 2 * nside; ell = np.arange(lmax + 1)
    cl = np.r_[0, 0, 1. / (ell[2:] * (ell[2:] + 1.)) ** 1.5]
    np.random.seed(3); phi = hp.synalm(cl, lmax=lmax, new=True)
    kap_alm = hp.almxfl(phi, ell * (ell + 1.) / 2); kap = hp.alm2map(kap_alm, nside, lmax=lmax)
    dth, dph = deflection_maps(kap_alm, nside, lmax)
    mask = np.ones(hp.nside2npix(nside)); mask[hp.query_disc(nside, (1, 0, 0), 0.8)] = 0
    S = Spectra(lmax, width=20)
    f1 = S.field(mask, [dth, dph], spin=1, key='m'); f0 = S.field(mask, [kap], spin=0, key='m')
    cE, cB = S.cross(f1, f0)
    ckk = hp.alm2cl(kap_alm, lmax=lmax); pred = S.theory(f1, f0, [2. / np.sqrt(np.maximum(ell * (ell + 1.), 1)) * ckk, 0 * ckk])[0]
    use = (S.ell_eff > 20) & (S.ell_eff < 200)
    assert np.allclose(cE[use] / pred[use], 1.0, atol=0.03)
    assert np.all(np.abs(cB[use]) < 0.01 * np.abs(cE[use]))


def test_fit_amplitude_recovers_scaling():
    rng = np.random.default_rng(0); T = rng.uniform(1, 2, 10); cov = np.diag(np.full(10, 0.01))
    y = 1.7 * T + rng.normal(0, 0.1, 10); use = np.ones(10, bool)
    b, sb, chi2, dof = fit_amplitude(y, T, cov, use)
    assert abs(b - 1.7) < 4 * sb and dof == 9 and chi2 < 30
