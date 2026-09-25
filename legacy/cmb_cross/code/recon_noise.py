"""Sample-variance-limited reconstruction noise for kappa from the Lya forest.

Implements Eqs. (N_kappa) of report/main.tex, Sec. "Continuum limit":

  1/N_kappa(L) = D int_{-inf}^{inf} dk_par/(2pi) I_kk(L, k_par),
  I_XY(L, k_par) = int d^2l/(2pi)^2  f_X(l, L-l) f_Y(l, L-l) / [2 P(l) P(|L-l|)],
  f_kappa(l, l') = (2/L^2) [ L.l P(l) + L.l' P(l') ],      (lensing response)
  f_a(l, l')     = P(l) + P(l'),                           (amplitude-modulation response)
  P(l) == P_F(k_par, k_perp = l/chi)  (the 1/chi^2 prefactor cancels in the ratios).

D is the comoving line-of-sight length over which the field is fully sampled
(for a survey with n_q sightlines/deg^2 each of mean forest length D_f, use
D = D_f together with the l_max set by n_q).
The bias-hardened noise against a local amplitude modulation a(x) is
  N_BH = N_kappa / (1 - rho^2),  rho^2 = F_ka^2 / (F_kk F_aa),
either with the F's summed over k_par first ("global" hardening: a is
k_par-independent) or per k_par slice ("per-slice": a may depend on k_par).
"""
import numpy as np
from scipy.integrate import trapezoid
from forest_power import ForestPower
from cosmo import chi as chi_of_z

DEG2 = (np.pi / 180) ** 2


def lmax_from_density(n_per_deg2, z=2.4):
    """Transverse Nyquist multipole for mean sightline separation 1/sqrt(n)."""
    theta_sep = np.sqrt(DEG2 / n_per_deg2)  # radians
    return np.pi / theta_sep


class ReconNoise:
    def __init__(self, pf=None, z=2.4, D=350.0, kpar_min=0.03, kpar_max=2.0,
                 lmax=1300.0, lmin=10.0, n_kpar=40, n_l=160, n_phi=64, noise_fn=None):
        """noise_fn(kpar) -> 3D noise power in the same units as P_F ((Mpc/h)^3), added to the
        denominators only (the lensing response is that of the signal)."""
        self.noise_fn = noise_fn
        self.pf = pf or ForestPower(z=z)
        self.chi = float(chi_of_z(z))
        self.D, self.lmax, self.lmin = D, lmax, lmin
        self.kpar = np.logspace(np.log10(kpar_min), np.log10(kpar_max), n_kpar)
        self.l = np.logspace(np.log10(lmin), np.log10(lmax), n_l)
        self.phi = np.linspace(0, 2 * np.pi, n_phi, endpoint=False)

    def P(self, kpar, l):
        return self.pf(kpar, l / self.chi)

    def fisher_slices(self, L):
        """Return arrays (n_kpar,) of I_kk, I_aa, I_ka at multipole L."""
        l = self.l[:, None]; ph = self.phi[None, :]
        lx, ly = l * np.cos(ph), l * np.sin(ph)
        l2x, l2y = L - lx, -ly
        l2 = np.hypot(l2x, l2y)
        Ldotl = L * lx; Ldotl2 = L * l2x
        mask = (l2 >= self.lmin) & (l2 <= self.lmax)
        Ikk = np.zeros(len(self.kpar)); Iaa = np.zeros_like(Ikk); Ika = np.zeros_like(Ikk)
        for i, kp in enumerate(self.kpar):
            P1 = self.P(kp, l * np.ones_like(ph))
            P2 = self.P(kp, l2)
            Nk = self.noise_fn(kp) if self.noise_fn is not None else 0.0
            fk = 2.0 / L ** 2 * (Ldotl * P1 + Ldotl2 * P2)
            fa = P1 + P2
            w = mask * l / (2 * (P1 + Nk) * (P2 + Nk)) / (2 * np.pi) ** 2  # d^2l = l dl dphi
            Ikk[i] = self._int(fk * fk * w)
            Iaa[i] = self._int(fa * fa * w)
            Ika[i] = self._int(fk * fa * w)
        return Ikk, Iaa, Ika

    def _int(self, g):
        # g has shape (n_l, n_phi); integrate dphi (uniform) then dl (log grid)
        dphi = 2 * np.pi / len(self.phi)
        gl = g.sum(axis=1) * dphi
        return trapezoid(gl, self.l)

    def _kpar_sum(self, I):
        # D int_{-inf}^{inf} dk/(2pi) I = 2 D int_0^inf dk/(2pi) I
        return 2 * self.D * trapezoid(I, self.kpar) / (2 * np.pi)

    def noise(self, Ls):
        """Return dict with N_kappa, N_BH_global, N_BH_perslice at each L."""
        out = {k: np.zeros(len(Ls)) for k in ("N", "N_bh_global", "N_bh_slice", "rho2_global", "R_ka")}
        for j, L in enumerate(Ls):
            if L >= 2 * self.lmax:   # no pixel pair can carry this multipole
                for k in out:
                    out[k][j] = np.inf if k != "rho2_global" else np.nan
                continue
            Ikk, Iaa, Ika = self.fisher_slices(L)
            Fkk, Faa, Fka = map(self._kpar_sum, (Ikk, Iaa, Ika))
            out["N"][j] = 1.0 / Fkk
            rho2 = Fka ** 2 / (Fkk * Faa)
            out["rho2_global"][j] = rho2
            out["R_ka"][j] = Fka / Fkk          # response of kappa-hat to a unit amplitude modulation
            out["N_bh_global"][j] = 1.0 / (Fkk * (1 - rho2))
            rho2s = Ika ** 2 / (Ikk * Iaa)
            out["N_bh_slice"][j] = 1.0 / self._kpar_sum(Ikk * (1 - rho2s))
        return out


def white_noise_check():
    """Mode-counting check: for P = const, 1/N = D k_max l_max^2/(2 pi^2) (report Eq. mode-count)."""
    class Flat:
        def __call__(self, kpar, kperp):
            return np.ones(np.broadcast(kpar, kperp).shape)
    rn = ReconNoise(pf=Flat(), D=100.0, kpar_min=1e-3, kpar_max=1.0, lmax=1000.0, lmin=1e-3, n_l=400)
    L = 50.0
    Ikk, _, _ = rn.fisher_slices(L)
    got = rn._kpar_sum(Ikk)
    exp = rn.D * (rn.kpar[-1] - rn.kpar[0]) * (rn.lmax ** 2 - rn.lmin ** 2) / (2 * np.pi ** 2)
    return got, exp


if __name__ == "__main__":
    g, e = white_noise_check()
    print(f"white-noise mode-count check: got {g:.4e}, expected {e:.4e}, ratio {g/e:.4f} "
          f"(ratio<1 because pairs with |L-l|>lmax are excluded)")
