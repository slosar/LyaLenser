"""Shared cosmology helpers (CAMB, Planck-2018-like flat LCDM).

Units: comoving distances in Mpc/h, wavenumbers in h/Mpc.
"""
import numpy as np
import camb
from functools import lru_cache

H0 = 67.66
h = H0 / 100.0
OMBH2, OMCH2, NS, AS = 0.02242, 0.11933, 0.9665, 2.105e-9
C_KMS = 299792.458


@lru_cache(maxsize=None)
def _pars(zs=(0.0,), nonlinear=False, kmax=50.0):
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=H0, ombh2=OMBH2, omch2=OMCH2, mnu=0.06, tau=0.054)
    pars.InitPower.set_params(As=AS, ns=NS)
    pars.set_matter_power(redshifts=sorted(set(zs), reverse=True), kmax=kmax)
    pars.NonLinear = camb.model.NonLinear_both if nonlinear else camb.model.NonLinear_none
    return pars


@lru_cache(maxsize=None)
def background():
    return camb.get_background(_pars())


def chi(z):
    """Comoving radial distance in Mpc/h."""
    z = np.asarray(z, dtype=float)
    return background().comoving_radial_distance(np.atleast_1d(z)).reshape(z.shape) * h


def hubble(z):
    """H(z) in km/s/Mpc (not divided by h)."""
    z = np.asarray(z, dtype=float)
    return background().hubble_parameter(np.atleast_1d(z)).reshape(z.shape)


def z_of_chi(chi_mpch):
    c = np.asarray(chi_mpch, float)
    return background().redshift_at_comoving_radial_distance(np.atleast_1d(c) / h).reshape(c.shape)


@lru_cache(maxsize=None)
def linear_pk_interp(zmax=4.0, kmax=50.0, nonlinear=False):
    """CAMB P(k,z) interpolator, k in h/Mpc, P in (Mpc/h)^3."""
    pars = _pars(zs=(0.0, zmax), nonlinear=nonlinear, kmax=kmax)
    return camb.get_matter_power_interpolator(
        pars, nonlinear=nonlinear, hubble_units=True, k_hunit=True,
        kmax=kmax, zmax=zmax, var1='delta_tot', var2='delta_tot')


@lru_cache(maxsize=None)
def _growth_interp():
    return linear_pk_interp(zmax=6.0, kmax=50.0, nonlinear=False)


def growth(z, k_ref=0.05):
    """Linear growth factor D(z)/D(0) from the ratio of the linear power spectrum at a fixed wavenumber."""
    z = np.asarray(z, dtype=float)
    pk = _growth_interp()
    z1 = np.atleast_1d(z)
    return np.sqrt(pk.P(z1, np.full(z1.shape, k_ref), grid=False) / pk.P(0.0, k_ref)).reshape(z.shape)


def sigma8():
    res = camb.get_results(_pars(zs=(0.0,), kmax=10.0))
    return float(res.get_sigma8_0())


def growth_rate(z):
    """f = dlnD/dlna, from Omega_m(z)^0.55 (adequate for the toy)."""
    om = (OMBH2 + OMCH2 + 0.06 / 93.14) / h ** 2
    omz = om * (1 + z) ** 3 / (om * (1 + z) ** 3 + 1 - om)
    return omz ** 0.55


if __name__ == "__main__":
    for z in [2.1, 2.4, 2.8, 3.2]:
        print(f"z={z}: chi={chi(z):.0f} Mpc/h, H={hubble(z):.1f} km/s/Mpc, "
              f"1 arcmin = {chi(z)*np.pi/180/60:.2f} Mpc/h, "
              f"1 Angstrom at Lya = {C_KMS/(1215.67*(1+z))*(1+z)/hubble(z)*h:.3f} Mpc/h")
    print("sigma8 =", sigma8())
