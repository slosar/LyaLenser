"""Driver: sample-variance-limited N_kappa(L) for DESI-like sightline densities,
comparison with C_L^{kappa kappa}(z_s=2.4), and toy S/N for kappa_lya x kappa_cmb.
Writes figures to ../report/figures and numbers to ../report/numbers.tex.
"""
import numpy as np, json
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from recon_noise import ReconNoise, lmax_from_density
from cross_spectrum import spectra, calibrate_white_noise, cross_snr
from cosmo import chi as chi_of_z

Z = 2.4
D_FOREST = 350.0      # Mpc/h, mean comoving forest length (1040-1200 A rest for z_q~2.8)
KPAR_MIN = 0.03       # h/Mpc, modes below this are removed by continuum fitting
Ls = np.unique(np.concatenate([np.arange(10, 100, 10), np.arange(100, 1000, 50), np.arange(1000, 3001, 250)])).astype(float)
Lfine = np.arange(2, 3001)

densities = [20, 50, 100, 400]   # sightlines per deg^2
kpar_maxs = [1.0, 2.0, 5.0]

import os, pickle
CACHE = "../report/recon_results.pkl"   # delete to force recomputation (~12 min)
results = pickle.load(open(CACHE, "rb")) if os.path.exists(CACHE) else {}
for n in densities:
    lmax = lmax_from_density(n, Z)
    for km in kpar_maxs:
        if (n, km) not in results:
            rn = ReconNoise(z=Z, D=D_FOREST, kpar_min=KPAR_MIN, kpar_max=km, lmax=lmax)
            results[(n, km)] = rn.noise(Ls)
            pickle.dump(results, open(CACHE, "wb"))
        print(f"n={n:4d}/deg2 lmax={lmax:6.0f} kpar_max={km}: N(L=100)={results[(n,km)]['N'][Ls==100][0]:.3e} "
              f"N_BH_global={results[(n,km)]['N_bh_global'][Ls==100][0]:.3e} N_BH_slice={results[(n,km)]['N_bh_slice'][Ls==100][0]:.3e}", flush=True)

S = spectra(Lfine, z_lya=Z)
Nc_act = calibrate_white_noise(Lfine, S['cc'], 43.0, 40, 763, 0.23)
Nc_pl = calibrate_white_noise(Lfine, S['cc'], 42.0, 8, 400, 0.67)

# ---- Figure 1: N_kappa(L) vs C_L^{kk}
fig, ax = plt.subplots(figsize=(6.8, 4.6))
ax.loglog(Lfine, S['ll'], 'k-', lw=2, label=r"$C_L^{\kappa\kappa}(z_s=2.4)$")
ax.loglog(Lfine, S['lc'], 'k--', lw=1.5, label=r"$C_L^{\kappa_{\rm Ly\alpha}\kappa_{\rm CMB}}$")
cols = plt.cm.viridis(np.linspace(0, 0.9, len(densities)))
for c, n in zip(cols, densities):
    r = results[(n, 2.0)]
    ax.loglog(Ls, r['N'], color=c, label=rf"$n_q={n}$/deg$^2$")
    ax.loglog(Ls, r['N_bh_slice'], color=c, ls=':')
ax.loglog(Ls, results[(50, 5.0)]['N'], color=cols[1], ls='-.', label=r"$n_q=50$, $k_{\parallel,\max}=5$")
ax.axhline(Nc_act, color='r', lw=0.8, ls='--', label=r"$N^{\rm CMB}$ ACT-like")
ax.axhline(Nc_pl, color='r', lw=0.8, ls=':', label=r"$N^{\rm CMB}$ Planck-like")
ax.set_xlabel(r"$L$"); ax.set_ylabel(r"$C_L$, $N_L$")
ax.set_ylim(1e-10, 1e-4); ax.set_xlim(10, 3000)
ax.set_title(r"$N_\kappa(L)$: naive (solid), amplitude-hardened (dotted); $k_{\parallel,\max}=2$", fontsize=10)
ax.legend(fontsize=7, ncol=2, loc="upper left"); fig.tight_layout(); fig.savefig("../report/figures/recon_noise.pdf")

# ---- Figure 2: cumulative S/N for the cross-correlation
fig, ax = plt.subplots(figsize=(6.4, 4.2))
snr_table = {}
for c, n in zip(cols, densities):
    r = results[(n, 2.0)]
    for key, ls in (('N', '-'), ('N_bh_slice', ':')):
        Nll = np.interp(Lfine, Ls, np.where(np.isfinite(r[key]), r[key], 1e30))
        for cmb, Nc, fsky in (('ACT', Nc_act, 0.15), ('Planck', Nc_pl, 0.22)):
            s2, cum = cross_snr(Lfine, S['lc'], S['ll'], Nll, S['cc'], Nc, fsky)
            snr_table[(n, key, cmb)] = float(cum[-1])
            if cmb == 'ACT':
                ax.semilogx(Lfine, cum, color=c, ls=ls, label=(rf"$n_q={n}$/deg$^2$" if key == 'N' else None))
ax.set_xlabel(r"$L_{\max}$"); ax.set_ylabel(r"cumulative S/N of $C_L^{\kappa_{\rm Ly\alpha}\kappa_{\rm CMB}}$")
ax.set_title(r"ACT-DR6-like $\kappa_{\rm CMB}$, $f_{\rm sky}=0.15$; dotted: amplitude-hardened")
ax.legend(fontsize=8); fig.tight_layout(); fig.savefig("../report/figures/cross_snr.pdf")

# ---- Figure 3: rho^2 (degeneracy of kappa with amplitude modulation)
fig, ax = plt.subplots(figsize=(6.0, 3.8))
for c, n in zip(cols, densities):
    ax.semilogx(Ls, results[(n, 2.0)]['rho2_global'], color=c, label=rf"$n_q={n}$/deg$^2$")
ax.set_xlabel(r"$L$"); ax.set_ylabel(r"$\rho^2_{\kappa a}$"); ax.set_ylim(0, 1); ax.legend(fontsize=8)
ax.set_title(r"squared correlation of $\kappa$ and amplitude responses")
fig.tight_layout(); fig.savefig("../report/figures/rho2.pdf")

# ---- numbers for the report
chi = float(chi_of_z(Z))
out = {"chi": chi, "Nc_act": Nc_act, "Nc_pl": Nc_pl,
       "lmax": {n: float(lmax_from_density(n, Z)) for n in densities},
       "N100": {f"{n}_{km}": float(results[(n, km)]['N'][Ls == 100][0]) for n in densities for km in kpar_maxs},
       "Nbh100": {f"{n}_{km}": float(results[(n, km)]['N_bh_slice'][Ls == 100][0]) for n in densities for km in kpar_maxs},
       "snr": {f"{n}_{k}_{c}": v for (n, k, c), v in snr_table.items()},
       "Cll100": float(S['ll'][Lfine == 100][0]), "Clc100": float(S['lc'][Lfine == 100][0]), "Ccc100": float(S['cc'][Lfine == 100][0])}
def _clean(o):
    if isinstance(o, dict): return {k: _clean(v) for k, v in o.items()}
    if isinstance(o, float) and not np.isfinite(o): return None
    return o
json.dump(_clean(out), open("../report/numbers.json", "w"), indent=1)
print(json.dumps(out, indent=1))
