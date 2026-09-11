"""Fiducial 3D Lyman-alpha forest flux power spectrum P_F(k_par, k_perp).

Model: linear-theory Kaiser form times the nonlinear correction of
Arinyo-i-Prats et al. (2015),
    P_F(k, mu) = b_F^2 (1 + beta_F mu^2)^2 P_lin(k, z) D_NL(k, mu),
    D_NL = exp{ [q1 Delta^2(k) + q2 Delta^4(k)] [1 - (k/k_v)^{a_v} mu^{b_v}] - (k/k_p)^2 },
    Delta^2(k) = k^3 P_lin / (2 pi^2).
Parameter values are APPROXIMATE fiducial numbers at z ~ 2.4 (to be refined
against the published tables / DESI fits later); they only set the shape of the
spectrum, which is all the sample-variance-limited noise depends on.
Units: k in h/Mpc, P in (Mpc/h)^3.
"""
import numpy as np
from scipy.integrate import trapezoid
from cosmo import linear_pk_interp

FID = dict(z=2.4, b_F=-0.13, beta_F=1.6, q1=0.6, q2=0.0, kv=1.0, av=0.55, bv=1.6, kp=16.0)


class ForestPower:
    def __init__(self, model="arinyo", **kw):
        self.p = dict(FID)
        self.p.update(kw)
        self.model = model
        self._pk = linear_pk_interp()

    def plin(self, k):
        k = np.atleast_1d(np.asarray(k, float))
        return self._pk.P(self.p["z"], np.clip(k, 1e-5, 49.0))

    def __call__(self, kpar, kperp):
        """P_F(k_par, k_perp) in (Mpc/h)^3; broadcasting inputs."""
        kpar = np.abs(np.asarray(kpar, float))
        kperp = np.asarray(kperp, float)
        k = np.sqrt(kpar ** 2 + kperp ** 2)
        k = np.maximum(k, 1e-6)
        mu = kpar / k
        p = self.p
        pl = self.plin(k).reshape(k.shape)
        kaiser = p["b_F"] ** 2 * (1 + p["beta_F"] * mu ** 2) ** 2 * pl
        if self.model == "kaiser":
            return kaiser * np.exp(-(k / p["kp"]) ** 2)
        d2 = k ** 3 * pl / (2 * np.pi ** 2)
        lnD = (p["q1"] * d2 + p["q2"] * d2 ** 2) * (1 - (k / p["kv"]) ** p["av"] * mu ** p["bv"]) - (k / p["kp"]) ** 2
        return kaiser * np.exp(lnD)

    def p1d(self, kpar, kperp_max=50.0, n=4000):
        """P_1D(k_par) = int d^2k_perp/(2pi)^2 P_F, for sanity checks."""
        kperp = np.logspace(-4, np.log10(kperp_max), n)
        out = []
        for kp in np.atleast_1d(kpar):
            integrand = self(kp, kperp) * kperp / (2 * np.pi)
            out.append(trapezoid(integrand, kperp))
        return np.array(out)


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    pf = ForestPower()
    kpar = np.logspace(-2, 0.7, 40)
    p1 = pf.p1d(kpar)
    print("k_par [h/Mpc]   k P1D/pi")
    for k, p in zip(kpar[::8], p1[::8]):
        print(f"{k:8.3f}  {k*p/np.pi:9.4f}")
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
    ax[0].loglog(kpar, kpar * p1 / np.pi)
    ax[0].set_xlabel(r"$k_\parallel$ [h/Mpc]"); ax[0].set_ylabel(r"$k_\parallel P_{1D}/\pi$")
    ax[0].set_title("fiducial P1D (z=2.4)")
    kperp = np.logspace(-2, 1.5, 200)
    for kp in [0.05, 0.2, 0.5, 1.0, 2.0]:
        ax[1].loglog(kperp, pf(kp, kperp), label=rf"$k_\parallel={kp}$")
    ax[1].set_xlabel(r"$k_\perp$ [h/Mpc]"); ax[1].set_ylabel(r"$P_F$ [(Mpc/h)$^3$]"); ax[1].legend(fontsize=8)
    ax[1].set_title(r"transverse shape at fixed $k_\parallel$")
    fig.tight_layout(); fig.savefig("../report/figures/forest_power.pdf")
    print("saved ../report/figures/forest_power.pdf")
