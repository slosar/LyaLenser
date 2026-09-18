"""NaMaster pseudo-C_ell helpers for the Stage B tracer maps (iteration 10).

Replaces the f_sky-scaled healpy pseudo-spectra and the 16-realisation Monte-Carlo mask transfer of iteration 7:
the mode-coupling matrix of each (binary) mask is computed analytically, the bandpowers are decoupled, the theory
is pushed through the same bandpower windows, and the bandpower covariance is the NaMaster Gaussian estimate
(mode coupling included). Bins are the width-40 annuli of the bias fit, starting at ell = 0.

Spin-1 convention (verified numerically, `tests/test_nmt.py`): with the deflection components ordered as
(d phi / d theta, d phi / d varphi / sin theta), i.e. (-north, east) of `templates.alpha_at`, the E-mode alm of
the spin-1 field is +sqrt(l(l+1)) phi_lm, so a gradient template phi = 2 kappa / [l(l+1)] has
C_l^{E kappa'} = 2 / sqrt(l(l+1)) C_l^{kappa kappa'} and no B-mode.
"""
from __future__ import annotations
import numpy as np
import healpy as hp
import pymaster as nmt


class Spectra:
    """Binning, fields, workspaces (cached per mask pair) and Gaussian covariances on one HEALPix resolution."""

    def __init__(self, lmax, width=40, n_iter=0):   # n_iter 0: the template alm are map2alm(iter=0), so the spectra must see the same alm
        self.lmax = int(lmax); self.width = int(width); self.n_iter = int(n_iter)
        edges = np.r_[np.arange(0, self.lmax, self.width), self.lmax + 1]
        self.edges = edges
        self.bins = nmt.NmtBin.from_edges(edges[:-1], edges[1:])
        self.ell_eff = self.bins.get_effective_ells()
        self.centres = .5 * (edges[:-1] + edges[1:] - 1)
        self._ws = {}; self._cws = {}

    # ---- fields
    def field(self, mask, maps, spin=0, key=None):
        f = nmt.NmtField(np.asarray(mask, float), [np.asarray(m, float) for m in np.atleast_2d(maps)],
                         spin=spin, lmax=self.lmax, n_iter=self.n_iter, masked_on_input=False)
        f._key = key if key is not None else id(mask)
        return f

    # ---- mode coupling
    def workspace(self, f1, f2):
        k = (f1._key, f1.spin, f2._key, f2.spin)
        if k not in self._ws:
            self._ws[k] = nmt.NmtWorkspace.from_fields(f1, f2, self.bins)
        return self._ws[k]

    def cross(self, f1, f2):
        """Decoupled bandpowers [n_cls, n_bins] of two fields."""
        w = self.workspace(f1, f2)
        return w.decouple_cell(nmt.compute_coupled_cell(f1, f2))

    def theory(self, f1, f2, cls):
        """Theory spectra [n_cls, lmax+1] pushed through the bandpower windows of the (f1, f2) workspace."""
        w = self.workspace(f1, f2)
        cls = np.atleast_2d(np.asarray(cls, float))[:, :self.lmax + 1]
        return w.decouple_cell(w.couple_cell(cls))

    def fsky_binned(self, cl):
        """The plain (2l+1)-weighted binning of a full-sky spectrum, for the mask-transfer diagnostic."""
        cl = np.asarray(cl, float)[:self.lmax + 1]; ell = np.arange(self.lmax + 1); out = np.zeros(len(self.centres))
        for i, (a, b) in enumerate(zip(self.edges[:-1], self.edges[1:])):
            m = (ell >= a) & (ell < b); w = 2 * ell[m] + 1.
            out[i] = np.sum(w * cl[m]) / w.sum()
        return out

    # ---- covariance
    def covariance_workspace(self, fa1, fa2, fb1, fb2):
        k = tuple((f._key, f.spin) for f in (fa1, fa2, fb1, fb2))
        if k not in self._cws:
            self._cws[k] = nmt.NmtCovarianceWorkspace.from_fields(fa1, fa2, fb1, fb2)
        return self._cws[k]

    def gaussian_covariance(self, fa1, fa2, fb1, fb2, cl_a1b1, cl_a1b2, cl_a2b1, cl_a2b2):
        """Covariance of the (fa1 x fa2) and (fb1 x fb2) bandpowers for spin-0 fields, given the (unmasked)
        spectra of the field pairs at every ell (signal plus noise)."""
        cw = self.covariance_workspace(fa1, fa2, fb1, fb2)
        wa = self.workspace(fa1, fa2); wb = self.workspace(fb1, fb2)
        f = lambda c: [np.asarray(c, float)[:self.lmax + 1]]
        return nmt.gaussian_covariance(cw, 0, 0, 0, 0, f(cl_a1b1), f(cl_a1b2), f(cl_a2b1), f(cl_a2b2), wa, wb)

    def spectrum_model(self, binned, floor=0.0):
        """Full-ell model of a spectrum from its bandpowers (linear interpolation, flat extrapolation), for the
        covariance inputs."""
        ell = np.arange(self.lmax + 1)
        return np.maximum(np.interp(ell, self.ell_eff, np.asarray(binned, float)), floor)


def fit_amplitude(y, T, cov, use):
    """Generalised least-squares amplitude b of y = b T over the bandpowers ``use``: returns (b, sigma_b, chi2, dof)."""
    y = np.asarray(y, float)[use]; T = np.asarray(T, float)[use]; C = np.asarray(cov, float)[np.ix_(use, use)]
    Ci = np.linalg.pinv(C)
    var = 1.0 / float(T @ Ci @ T); b = float(T @ Ci @ y) * var
    r = y - b * T
    return b, float(np.sqrt(var)), float(r @ Ci @ r), int(use.sum() - 1)


def deflection_maps(kappa_alm, nside, lmax=None, filt=None):
    """(d phi/d theta, d phi/d varphi / sin theta) of phi = 2 kappa / [l(l+1)] on the sphere: the spin-1 field of
    `Spectra.field(..., spin=1)`. ``filt[l]`` multiplies the alm first (a band window)."""
    alm = np.asarray(kappa_alm, complex); lm = hp.Alm.getlmax(len(alm)) if lmax is None else int(lmax)
    ell = np.arange(lm + 1); f = np.r_[0., 0., 2. / (ell[2:] * (ell[2:] + 1.))]
    if filt is not None: f = f * np.asarray(filt, float)[:lm + 1]
    phi = hp.almxfl(alm, f)
    _, dth, dph = hp.alm2map_der1(phi, nside, lmax=lm)
    return dth, dph


def gaussian_covariance_any(S: Spectra, fa1, fa2, fb1, fb2, cl_a1b1, cl_a1b2, cl_a2b1, cl_a2b2):
    """`Spectra.gaussian_covariance` for fields of any spin: each cl_* is [n_spectra, lmax+1] in NaMaster's order
    (spin-0 x spin-0: 1; spin-0 x spin-1: [E, B]; spin-1 x spin-1: [EE, EB, BE, BB]). Returns the covariance
    reshaped to [n_bins, n_cls_a, n_bins, n_cls_b]."""
    cw = S.covariance_workspace(fa1, fa2, fb1, fb2); wa = S.workspace(fa1, fa2); wb = S.workspace(fb1, fb2)
    f = lambda c: [np.asarray(x, float)[:S.lmax + 1] for x in np.atleast_2d(c)]
    n_a = len(f(cl_a1b1)) and (2 if (fa1.spin + fa2.spin) == 1 else (4 if fa1.spin + fa2.spin == 2 else 1))
    n_b = 2 if (fb1.spin + fb2.spin) == 1 else (4 if fb1.spin + fb2.spin == 2 else 1)
    cov = nmt.gaussian_covariance(cw, fa1.spin, fa2.spin, fb1.spin, fb2.spin, f(cl_a1b1), f(cl_a1b2), f(cl_a2b1), f(cl_a2b2), wa, wb)
    nb = len(S.ell_eff)
    return cov.reshape(nb, n_a, nb, n_b)
