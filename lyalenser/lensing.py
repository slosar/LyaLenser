"""Lensing kernel and Limber integrals (flat-sky, single source plane).

W_s(chi) = (3/2) Omega_m (H0/c)^2 (1+z) chi (chi_s - chi)/chi_s   [chi in Mpc/h, H0/c = 1/2997.9 h/Mpc]
C_L^{XY} = int dchi W_X(chi) W_Y(chi) / chi^2  P(k=(L+1/2)/chi, z(chi)).
"""
import numpy as np
from scipy.integrate import trapezoid
from lyalenser.cosmo import chi as chi_of_z, z_of_chi, linear_pk_interp, h, OMBH2, OMCH2

OM = (OMBH2 + OMCH2 + 0.06 / 93.14) / h ** 2
H0C = 1.0 / 2997.92458      # h/Mpc
Z_CMB = 1090.0
DEG2 = (np.pi / 180) ** 2   # steradians per square degree


def kernel(chi, chi_s):
    """Convergence kernel of a source at comoving distance chi_s (zero beyond the source)."""
    z = z_of_chi(chi)
    return 1.5 * OM * H0C ** 2 * (1 + z) * chi * np.clip(chi_s - chi, 0, None) / chi_s


def kernel_dsource(chi, chi_s):
    """dW/dchi_s: the derivative of the convergence kernel with respect to the source distance,
    (3/2) Omega_m (H0/c)^2 (1+z) chi^2 / chi_s^2 for chi < chi_s (zero beyond the source). A source at chi_s + d
    sees the lens at chi with kernel W + d dW/dchi_s; the derivative maps of the templates (lowz_catalogues.py)
    are built with this weight in place of W."""
    z = z_of_chi(chi)
    return 1.5 * OM * H0C ** 2 * (1 + z) * chi * chi * (np.asarray(chi) < chi_s) / chi_s ** 2


def limber(Ls, W1, W2, chimin=1.0, chimax=None, nchi=800, to_recombination=True):
    """Limber integral of two kernels W1(chi), W2(chi) over [chimin, chimax] (Mpc/h). Non-linear (halofit)
    P(k, z) for z < 6; linear P(k, z) from z = 6 to recombination when the range reaches that far."""
    pk = linear_pk_interp(zmax=6.0, kmax=200.0, nonlinear=True)
    chi6 = float(chi_of_z(6.0))
    chimax = chimax or float(chi_of_z(Z_CMB))
    chis = np.linspace(chimin, min(chimax, chi6), nchi); zs = z_of_chi(chis)
    w = W1(chis) * W2(chis) / chis ** 2
    out = np.zeros(len(Ls))
    for i, L in enumerate(Ls):
        k = (L + 0.5) / chis
        out[i] = trapezoid(w * pk.P(zs, k, grid=False), chis)
    if to_recombination and chimax > chi6:
        pkl = linear_pk_interp(zmax=1100.0, kmax=200.0, nonlinear=False)
        chis2 = np.linspace(chi6, chimax * 0.999, nchi); zs2 = z_of_chi(chis2)
        w2 = W1(chis2) * W2(chis2) / chis2 ** 2
        for i, L in enumerate(Ls):
            k = (L + 0.5) / chis2
            out[i] += trapezoid(w2 * pkl.P(zs2, k, grid=False), chis2)
    return out
