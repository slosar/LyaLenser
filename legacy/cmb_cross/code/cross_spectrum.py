"""Limber angular power spectra for kappa_lya (single source plane) and kappa_cmb,
and the cross-correlation signal-to-noise (report Sec. "Toy numbers").

C_L^{XY} = int dchi W_X(chi) W_Y(chi) / chi^2  P_NL(k=(L+1/2)/chi, z(chi)),
W_s(chi) = (3/2) Om_m (H0/c)^2 (1+z) chi (chi_s - chi)/chi_s   [chi in Mpc/h, H0/c = 1/2997.9 h/Mpc].
"""
import numpy as np
from scipy.integrate import trapezoid
from cosmo import chi as chi_of_z, z_of_chi, linear_pk_interp, h, OMBH2, OMCH2

OM = (OMBH2 + OMCH2 + 0.06 / 93.14) / h ** 2
H0C = 1.0 / 2997.92458  # h/Mpc
Z_CMB = 1090.0


def kernel(chi, chi_s):
    z = z_of_chi(chi)
    w = 1.5 * OM * H0C ** 2 * (1 + z) * chi * np.clip(chi_s - chi, 0, None) / chi_s
    return w


def limber(Ls, chi_s1, chi_s2, nonlinear=True, nchi=400):
    """Nonlinear P(k,z) to z=6, linear beyond up to the lower source plane (review 1, item 11)."""
    pk = linear_pk_interp(zmax=6.0, kmax=200.0, nonlinear=nonlinear)
    chi6 = float(chi_of_z(6.0))
    chimax = min(chi_s1, chi_s2)
    chis = np.linspace(1.0, min(chimax, chi6), nchi)
    zs = z_of_chi(chis)
    w = kernel(chis, chi_s1) * kernel(chis, chi_s2) / chis ** 2
    out = np.zeros(len(Ls))
    for i, L in enumerate(Ls):
        k = (L + 0.5) / chis  # h/Mpc
        out[i] = trapezoid(w * pk.P(zs, k, grid=False), chis)
    if chimax > chi6:
        pkl = linear_pk_interp(zmax=1100.0, kmax=200.0, nonlinear=False)
        chis2 = np.linspace(chi6, chimax * 0.999, nchi); zs2 = z_of_chi(chis2)
        w2 = kernel(chis2, chi_s1) * kernel(chis2, chi_s2) / chis2 ** 2
        for i, L in enumerate(Ls):
            k = (L + 0.5) / chis2
            out[i] += trapezoid(w2 * pkl.P(zs2, k, grid=False), chis2)
    return out


def spectra(Ls, z_lya=2.4):
    chi_l = float(chi_of_z(z_lya)); chi_c = float(chi_of_z(Z_CMB))
    return dict(ll=limber(Ls, chi_l, chi_l), lc=limber(Ls, chi_l, chi_c), cc=limber(Ls, chi_c, chi_c))


def calibrate_white_noise(Ls, Ccc, snr, lmin, lmax, fsky):
    """Constant N_L^{kk} that reproduces a quoted auto-spectrum detection S/N over [lmin,lmax]."""
    from scipy.optimize import brentq
    m = (Ls >= lmin) & (Ls <= lmax)
    def f(logN):
        N = 10 ** logN
        return np.sqrt(np.sum(fsky * (2 * Ls[m] + 1) / 2 * (Ccc[m] / (Ccc[m] + N)) ** 2)) - snr
    return 10 ** brentq(f, -12, -3)


def cross_snr(Ls, Clc, Cll, Nll, Ccc, Ncc, fsky):
    """Per-L (S/N)^2 contributions for C_L^{lya x cmb}, and cumulative S/N."""
    var = ((Cll + Nll) * (Ccc + Ncc) + Clc ** 2) / ((2 * Ls + 1) * fsky)
    s2 = Clc ** 2 / var
    return s2, np.sqrt(np.cumsum(s2))


if __name__ == "__main__":
    Ls = np.arange(2, 2001)
    S = spectra(Ls)
    for L in [10, 40, 100, 300, 1000]:
        i = L - 2
        print(f"L={L:5d}  C_ll={S['ll'][i]:.3e}  C_lc={S['lc'][i]:.3e}  C_cc={S['cc'][i]:.3e}  r={S['lc'][i]/np.sqrt(S['ll'][i]*S['cc'][i]):.3f}")
    N_act = calibrate_white_noise(Ls, S['cc'], 43.0, 40, 763, 0.23)
    N_pl = calibrate_white_noise(Ls, S['cc'], 42.0, 8, 400, 0.67)
    print(f"white-noise-equivalent CMB kappa noise: ACT-DR6-like {N_act:.2e}, Planck-PR4-like {N_pl:.2e}")
