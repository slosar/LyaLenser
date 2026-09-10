"""Step 3: three tracers (quasars q, forest delta_F, CMB kappa) on one footprint.

Squeezed-limit model for the local forest power modulation p_alpha(L) (Sec. 4 of the report):
    <p_alpha kappa_cmb> = R_alpha C^{d kc} + T_alpha C^{kl kc}
    <p_alpha q>         = R_alpha C^{d q}  + T_alpha C^{kl q}
with d = slab-projected matter density at the forest redshift, kl = kappa_lya, kc = kappa_cmb,
R_alpha the unknown forest response, T_alpha the known lensing template.  The deprojected statistic
    D_alpha = <p_alpha kc> - beta <p_alpha q>,   beta = C^{d kc}/C^{d q},
is response-free.  Equivalently one cross-correlates kappa-hat with kc_perp = kc - beta q.
This script evaluates the noise of that measurement with a realistic forest noise model
    N_3D(k_par) = [P_N,1D + P_1D(k_par)] chi^2 / nbar_eff          (McQuinn & White 2011 form)
(pixel noise + Poisson-sampling/aliasing noise), quasar shot noise and CMB-kappa noise.
"""
import os, pickle, json
import numpy as np
from cosmo import chi as chi_of_z, z_of_chi, linear_pk_interp
from forest_power import ForestPower
from recon_noise import ReconNoise, lmax_from_density, DEG2
from cross_spectrum import kernel as wkappa, Z_CMB, calibrate_white_noise

# ---------------- fiducial survey / model parameters ----------------
Z_S = 2.4                     # effective forest source plane
Z1, Z2 = 2.1, 3.0             # forest / quasar-template slab
B_Q = 3.5                     # quasar bias at z~2.4
MAG = 0.5                     # quasar magnification coefficient 5s-2 (toy value)
R_FID = 2.0                   # fiducial forest response dlnP_F/ddelta_L (order of magnitude, Chiang+2017)
N_Q_SLAB = 25.0               # template quasars per deg^2 in the slab
KPAR_MIN, KPAR_MAX = 0.03, 2.0
NEFFS = [22.0, 50.0, 100.0]   # unweighted available sightlines per deg^2 at a given redshift (45 x 350/712 = 22 for DR1)
PN1DS = [0.0, 0.17, 0.33]     # 1D pixel-noise power sigma_delta^2 dchi_pix [Mpc/h]; DESI DR1 ~0.33, cleaned ~0.17
SIGLNS = [1.0, 2.0]           # ln-scatter of the per-quasar noise power for the C^-1-weighted cases (median 0.33)
WEIGHTED = [(22.0, 0.33, 1.0), (22.0, 0.33, 2.0), (50.0, 0.33, 1.0), (50.0, 0.33, 2.0)]
FSKY = {"ACT": 0.15, "Planck": 0.22}

chi_s, chi1, chi2, chi_cmb = (float(chi_of_z(z)) for z in (Z_S, Z1, Z2, Z_CMB))
D_SLAB = chi2 - chi1
Ls = np.unique(np.concatenate([np.arange(10, 100, 10), np.arange(100, 1000, 50), np.arange(1000, 3001, 250)])).astype(float)
Lfine = np.arange(2, 3001)


def limber(Ls, W1, W2, chimin=1.0, chimax=None, nchi=800, to_recombination=True):
    """Limber integral. Nonlinear (halofit) P(k,z) for z<6; linear P(k,z) from z=6 to recombination
    (review 1, item 11: the earlier z=6 cap dropped 7-32% of C^{kappa kappa}_CMB)."""
    pk = linear_pk_interp(zmax=6.0, kmax=200.0, nonlinear=True)
    chi6 = float(chi_of_z(6.0))
    chimax = chimax or chi_cmb
    chis = np.linspace(chimin, min(chimax, chi6), nchi); zs = z_of_chi(chis)
    w = W1(chis) * W2(chis) / chis ** 2
    out = np.zeros(len(Ls))
    for i, L in enumerate(Ls):
        k = (L + 0.5) / chis
        out[i] = np.trapz(w * pk.P(zs, k, grid=False), chis)
    if to_recombination and chimax > chi6:
        pkl = linear_pk_interp(zmax=1100.0, kmax=200.0, nonlinear=False)
        chis2 = np.linspace(chi6, chimax * 0.999, nchi); zs2 = z_of_chi(chis2)
        w2 = W1(chis2) * W2(chis2) / chis2 ** 2
        for i, L in enumerate(Ls):
            k = (L + 0.5) / chis2
            out[i] += np.trapz(w2 * pkl.P(zs2, k, grid=False), chis2)
    return out


W_kc = lambda c: wkappa(c, chi_cmb)
W_kl = lambda c: wkappa(c, chi_s)
W_d = lambda c: np.where((c >= chi1) & (c <= chi2), 1.0 / D_SLAB, 0.0)


W_m = lambda c: np.where((c >= chi1) & (c <= chi2), W_kc(c), 0.0)   # matched template kernel (b_q divided out)


def spectra():
    S = {}
    # kernel-matched template spectra (signal parts); shot noise of the weighted template:
    #   N_shot = D^2 <W_CMB^2>_slab / (b_q^2 n_2D)  (uniform n(z) in the slab)
    S["ss_sig"] = limber(Lfine, W_m, W_m, chi1, chi2); S["skc"] = S["ss_sig"].copy()
    S["skl"] = limber(Lfine, W_m, W_kl, chi1, chi2)
    cg = np.linspace(chi1, chi2, 400); S["shot_s"] = D_SLAB ** 2 * np.mean(W_kc(cg) ** 2) / (B_Q ** 2 * N_Q_SLAB / DEG2)
    S["W_range"] = (float(W_kc(np.array([chi1]))[0]), float(W_kc(np.array([chi2]))[0]))
    S["klkl"] = limber(Lfine, W_kl, W_kl); S["klkc"] = limber(Lfine, W_kl, W_kc); S["kckc"] = limber(Lfine, W_kc, W_kc)
    S["dd"] = limber(Lfine, W_d, W_d, chi1, chi2); S["dkc"] = limber(Lfine, W_d, W_kc, chi1, chi2)
    S["dkl"] = limber(Lfine, W_d, W_kl, chi1, chi2)
    # quasar field q = b delta + MAG kappa_q, with kappa_q ~ kappa_lya (same redshift)
    S["qq"] = B_Q ** 2 * S["dd"] + 2 * B_Q * MAG * S["dkl"] + MAG ** 2 * S["klkl"]
    S["qkc"] = B_Q * S["dkc"] + MAG * S["klkc"]
    S["qkl"] = B_Q * S["dkl"] + MAG * S["klkl"]
    S["shot_q"] = DEG2 / N_Q_SLAB
    return S


def forest_noise_fn(pf, neff_deg2, pn1d, chi, sigma_ln=0.0):
    """Return noise_fn(kpar) in P_F units.
    sigma_ln = 0: homogeneous population, N = [P_N + P_1D(kpar)] chi^2 / nbar_sr.
    sigma_ln > 0: per-quasar pixel-noise power log-normally distributed with median P_N and
    ln-scatter sigma_ln, combined with inverse-variance (C^-1) weights per sightline, so that
    1/N = nbar_sr/chi^2 * < 1/(P_N,a + P_1D) >_a   (harmonic mean over the population)."""
    nbar_sr = neff_deg2 / DEG2
    kgrid = np.logspace(np.log10(KPAR_MIN) - 0.1, np.log10(KPAR_MAX) + 0.1, 30)
    p1d = pf.p1d(kgrid)
    if sigma_ln <= 0 or pn1d <= 0:
        return lambda kp: (pn1d + np.interp(kp, kgrid, p1d)) * chi ** 2 / nbar_sr
    x, w = np.polynomial.hermite_e.hermegauss(60)          # Gauss-Hermite for the log-normal average
    w = w / w.sum()
    pna = pn1d * np.exp(sigma_ln * x)
    def fn(kp):
        p1 = np.interp(kp, kgrid, p1d)
        inv = np.sum(w / (pna + p1))
        return chi ** 2 / (nbar_sr * inv)
    return fn


def recon_all(cache="../report/recon3_results.pkl"):
    res = pickle.load(open(cache, "rb")) if os.path.exists(cache) else {}
    pf = ForestPower(z=Z_S)
    configs = [(n, pn, 0.0) for n in NEFFS for pn in PN1DS] + WEIGHTED
    for neff, pn, sig in configs:
        lmax = lmax_from_density(neff, Z_S)
        key = (neff, pn) if sig == 0 else (neff, pn, sig)
        if key in res:
            continue
        rn = ReconNoise(pf=pf, z=Z_S, D=D_SLAB, kpar_min=KPAR_MIN, kpar_max=KPAR_MAX, lmax=lmax,
                        noise_fn=forest_noise_fn(pf, neff, pn, chi_s, sig))
        res[key] = rn.noise(Ls)
        pickle.dump(res, open(cache, "wb"))
        print(f"neff={neff:5.0f}/deg2 lmax={lmax:5.0f} PN1D={pn:.2f} sigma_ln={sig}: N(100)={res[key]['N'][Ls==100][0]:.3e} "
              f"N_bh_slice={res[key]['N_bh_slice'][Ls==100][0]:.3e} R_ka(100)={res[key]['R_ka'][Ls==100][0]:.2f}", flush=True)
    return res


def snr_curves(S, r, Nc, fsky, Lrange=(2, 3000)):
    """Cumulative S/N for: naive (kc), hardened-slice/global (kc), deprojected with a single slab-averaged beta
    (kc - beta q), and deprojected with the kernel-matched template (kc - kappa_s-hat).
    Review-1 corrections: the estimator's response to the field modulation a is R_ka, and a = (R_delta/2) delta_s for
    dlnP_F = R_delta delta_s (item 3); beta_cancel includes the magnification term (item 5); the reconstruction
    auto-power includes the response sample variance 2a C^{kl d} + a^2 C^{dd} (item 8)."""
    interp = lambda a: np.interp(Lfine, Ls, np.where(np.isfinite(a), a, 1e30))
    N, Nbs, Nbg = interp(r["N"]), interp(r["N_bh_slice"]), interp(r["N_bh_global"])
    Rka = np.interp(Lfine, Ls, np.where(np.isfinite(r["R_ka"]), r["R_ka"], 0.0))
    a = 0.5 * R_FID * Rka                                     # kappa-hat = kappa_lya + a delta_s + n
    resp_var = 2 * a * S["dkl"] + a ** 2 * S["dd"]            # response sample variance in <kappa-hat kappa-hat>
    beta = S["dkc"] / (B_Q * S["dd"] + MAG * S["dkl"])        # exact cancellation coefficient in this model
    Sdep = S["klkc"] - beta * S["qkl"]
    Cperp = S["kckc"] + Nc - 2 * beta * S["qkc"] + beta ** 2 * (S["qq"] + S["shot_q"])
    # kernel-matched template: kappa_perp = kc - kappa_s-hat; beta = 1; template = slab part of kappa_CMB + shot
    Smat = S["klkc"] - S["skl"]
    Cmat = S["kckc"] + Nc - 2 * S["skc"] + S["ss_sig"] + S["shot_s"]
    m = (Lfine >= Lrange[0]) & (Lfine <= Lrange[1])
    def cum(sig, Nl, Cc, extra=0.0):
        var = ((S["klkl"] + extra + Nl) * Cc + sig ** 2) / ((2 * Lfine + 1) * fsky)
        return np.sqrt(np.cumsum(np.where(m, sig ** 2 / var, 0.0)))
    bias = a * S["dkc"]                                        # response contamination of the naive statistic
    return dict(naive=cum(S["klkc"], N, S["kckc"] + Nc, resp_var), bh_slice=cum(S["klkc"], Nbs, S["kckc"] + Nc),
                bh_global=cum(S["klkc"], Nbg, S["kckc"] + Nc), deproj=cum(Sdep, N, Cperp, resp_var),
                matched=cum(Smat, N, Cmat, resp_var), resp=cum(bias, N, S["kckc"] + Nc, resp_var),
                beta=beta, Sdep=Sdep, Cperp=Cperp, bias_over_signal=bias / S["klkc"],
                signal_subtracted_frac=S["skl"] / S["klkc"])


if __name__ == "__main__":
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    S = spectra()
    res = recon_all()
    Nc = {"ACT": calibrate_white_noise(Lfine, S["kckc"], 43.0, 40, 763, 0.23),
          "Planck": calibrate_white_noise(Lfine, S["kckc"], 42.0, 8, 400, 0.67)}
    out = {"params": dict(Z_S=Z_S, Z1=Z1, Z2=Z2, D_SLAB=D_SLAB, B_Q=B_Q, MAG=MAG, R_FID=R_FID, N_Q_SLAB=N_Q_SLAB,
                          lmax={str(n): float(lmax_from_density(n, Z_S)) for n in NEFFS}), "spec": {}, "table": []}
    i100, i40, i300 = (int(np.where(Lfine == L)[0][0]) for L in (100, 40, 300))
    beta = S["dkc"] / (B_Q * S["dd"] + MAG * S["dkl"])
    r2_qk = S["qkc"] ** 2 / ((S["qq"] + S["shot_q"]) * S["kckc"])
    out["matched"] = dict(shot_s=float(S["shot_s"]), W_range=S["W_range"],
                          ss_sig100=float(S["ss_sig"][i100]), skl_over_klkc100=float(S["skl"][i100] / S["klkc"][i100]))
    # magnification-induced error in beta if MAG ignored (Mathematica check 3b, first order)
    dbeta_over_beta = (-2 * S["dkc"] * S["dkl"] + S["dd"] * S["klkc"]) * MAG / (B_Q * S["dd"] * S["dkc"])
    for L, i in (("40", i40), ("100", i100), ("300", i300)):
        out["spec"][L] = dict(klkl=S["klkl"][i], klkc=S["klkc"][i], kckc=S["kckc"][i], dd=S["dd"][i], dkc=S["dkc"][i],
                              dkl=S["dkl"][i], beta=beta[i], r2_qk=r2_qk[i], r_lc=S["klkc"][i] / np.sqrt(S["klkl"][i] * S["kckc"][i]),
                              resp_over_lens=0.5 * R_FID * S["dkc"][i] / S["klkc"][i], dbeta_over_beta=dbeta_over_beta[i],
                              shot_q=S["shot_q"], beta2_shot=beta[i] ** 2 * S["shot_q"])
    for key in [(n, pn) for n in NEFFS for pn in PN1DS] + WEIGHTED:
            neff, pn = key[0], key[1]; sig = key[2] if len(key) == 3 else 0.0
            r = res[key]
            Rka100 = float(r["R_ka"][Ls == 100][0])
            row = dict(neff=neff, pn=pn, sigma_ln=sig, N100=float(r["N"][Ls == 100][0]), Nbh100=float(r["N_bh_slice"][Ls == 100][0]),
                       Nbg100=float(r["N_bh_global"][Ls == 100][0]), Rka100=Rka100,
                       bias_over_signal100=0.5 * Rka100 * R_FID * S["dkc"][i100] / S["klkc"][i100])
            for cmb in ("ACT", "Planck"):
                c = snr_curves(S, r, Nc[cmb], FSKY[cmb])
                for k in ("naive", "bh_slice", "bh_global", "deproj", "matched", "resp"):
                    row[f"snr_{k}_{cmb}"] = float(c[k][-1])
                c2 = snr_curves(S, r, Nc[cmb], FSKY[cmb], Lrange=(40, 300))
                for k in ("naive", "deproj", "matched", "resp"):
                    row[f"snr{k}_{cmb}_L40_300"] = float(c2[k][-1])
            out["table"].append(row)
    pf_ = ForestPower(z=Z_S)
    out["neff_factor"] = {}
    for sig in SIGLNS:
        f0 = forest_noise_fn(pf_, 25.0, 0.33, chi_s, 0.0); f1 = forest_noise_fn(pf_, 25.0, 0.33, chi_s, sig)
        out["neff_factor"][str(sig)] = {str(k): float(f0(k) / f1(k)) for k in (0.05, 0.1, 0.3, 1.0)}
    out["Nc"] = Nc
    json.dump(out, open("../report/numbers3.json", "w"), indent=1, default=float)

    # ---- Figure: N_kappa with noise, and cumulative deprojected S/N
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
    ax[0].loglog(Lfine, S["klkl"], "k-", lw=2, label=r"$C_L^{\kappa\kappa}(z_s=2.4)$")
    ax[0].loglog(Lfine, S["klkc"], "k--", lw=1.5, label=r"$C_L^{\kappa_{\rm Ly\alpha}\kappa_{\rm CMB}}$")
    cols = plt.cm.viridis(np.linspace(0, 0.9, len(NEFFS)))
    for c, neff in zip(cols, NEFFS):
        for pn, ls in zip(PN1DS, ("-", "--", ":")):
            r = res[(neff, pn)]
            ax[0].loglog(Ls, r["N"], color=c, ls=ls, label=(rf"$n_{{\rm eff}}={neff:.0f}$/deg$^2$, $P_N={pn}$" if pn in (0.0, 0.33) else None))
    ax[0].set_ylim(1e-9, 1e-3); ax[0].set_xlim(10, 3000); ax[0].set_xlabel("$L$"); ax[0].set_ylabel("$N_\\kappa(L)$")
    ax[0].set_title("reconstruction noise with pixel + sampling noise", fontsize=10); ax[0].legend(fontsize=7)
    for c, neff in zip(cols, NEFFS):
        for pn, ls in zip(PN1DS, ("-", "--", ":")):
            cur = snr_curves(S, res[(neff, pn)], Nc["ACT"], FSKY["ACT"])
            ax[1].semilogx(Lfine, cur["matched"], color=c, ls=ls, label=(rf"$n_{{\rm eff}}={neff:.0f}$, $P_N={pn}$" if pn in (0.0, 0.33) else None))
    ax[1].set_xlabel(r"$L_{\max}$"); ax[1].set_ylabel(r"cumulative S/N, $\kappa_{\rm CMB}-\hat\kappa_{\rm s}$ (matched template)")
    ax[1].set_title("ACT-like CMB noise, $f_{\\rm sky}=0.15$", fontsize=10); ax[1].legend(fontsize=7); ax[1].set_ylim(0, None)
    fig.tight_layout(); fig.savefig("../report/figures/three_tracer_snr.pdf")

    # ---- Figure: bias/signal of naive estimator and beta, r^2
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.0))
    for c, neff in zip(cols, NEFFS):
        r = res[(neff, 0.33)]
        Rka = np.interp(Lfine, Ls, np.where(np.isfinite(r["R_ka"]), r["R_ka"], 0.0))
        ax[0].semilogx(Lfine, 0.5 * Rka * R_FID * S["dkc"] / S["klkc"], color=c, label=rf"$n_{{\rm eff}}={neff:.0f}$/deg$^2$")
    ax[0].axhline(1, color="k", lw=0.5); ax[0].set_xlabel("$L$"); ax[0].set_ylabel(r"response bias / lensing signal, $R_\delta=2$ (corrected)")
    ax[0].set_title(r"naive $\langle\hat\kappa\,\kappa_{\rm CMB}\rangle$: contamination", fontsize=10); ax[0].legend(fontsize=8); ax[0].set_ylim(0, 6); ax[0].set_xlim(2, 1000)
    ax[1].semilogx(Lfine, r2_qk, label=r"$r^2_{q\kappa_{\rm CMB}}$ (incl. shot noise)")
    ax[1].semilogx(Lfine, S["dkc"] ** 2 / (S["dd"] * S["kckc"]), label=r"$r^2_{\delta\kappa_{\rm CMB}}$ (no shot noise)")
    ax[1].semilogx(Lfine, beta ** 2 * S["shot_q"] / S["kckc"], label=r"$\beta^2/\bar n_q\,/\,C_L^{\kappa\kappa}$")
    ax[1].set_xlabel("$L$"); ax[1].legend(fontsize=8); ax[1].set_ylim(0, 1.2); ax[1].set_xlim(2, 1000); ax[1].set_title("cost of deprojecting $q$ from $\\kappa_{\\rm CMB}$", fontsize=10)
    fig.tight_layout(); fig.savefig("../report/figures/three_tracer_bias.pdf")
    print(json.dumps(out["spec"], indent=1, default=float))
    for row in out["table"]:
        print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in row.items()})
