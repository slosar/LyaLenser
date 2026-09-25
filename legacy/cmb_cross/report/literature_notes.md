# Literature notes: weak lensing of the Lyman-alpha forest (compiled 2026-09-10)

All entries verified against ADS (bibcode + arXiv ID). Citation sweeps of Croft+2018, Metcalf+2018, Metcalf+2020, Foreman+2018, Schaan+2018, Pourtsidou & Metcalf 2014 on ADS plus an INSPIRE title search were used to find follow-ups. Helper: `code/ads_query.sh` (needs `ADS_API_TOKEN`).

## A. Lensing OF the forest (reconstruction from the 3D flux field)
- **Croft, Romeo & Metcalf 2018**, MNRAS 477, 1814 (arXiv:1706.07870, 2018MNRAS.477.1814C). First proposal. Lenses simulated spectra (sources z=2-3, lens plane z~1) on 15x15 deg, 1.8 arcmin sightline grid (~1/arcmin^2, "extreme high end"). Configuration-space QE adapted from 21cm lensing, stacked over redshift slices treated as independent. Recovers degree-scale potential morphology. No per-survey S/N.
- **Metcalf, Croft & Romeo 2018**, MNRAS 477, 2841 (arXiv:1706.08939, 2018MNRAS.477.2841M). Gaussian noise forecasts. Noise per potential mode ~ sigma_tau^2 /[j_max (l_max-l_min)^2 l_min^2]. Cumulative kappa-power S/N (dz=0.1 / 0.5 bins): BOSS 1.3/6.0, eBOSS 3.7/15.9, DESI 7.0/28.8 (770k spectra, 14000 deg2, 8.1 arcmin separation, per-pixel S/N 3.3), CLAMATO 2.7/7.7, MSE 53/162. DESI -> 1.5% amplitude at z~2.5.
- **Metcalf, Tessore & Croft 2020**, A&A 642, A122 (arXiv:2005.04109, 2020A&A...642A.122M). Real-space minimum-variance pixel-pair estimator phi_mu = 1/2 F^-1 (d C^-1 P^nu d - tr C^-1 P^nu), potential in Legendre polynomials over a patch; handles irregular sightlines, holes, inhomogeneous noise. Needs >~0.5 sources/arcmin^2 and signal-dominated pixels; 200/deg2 fails. Scaling sigma ~ N_spec^-0.5 L^2.6 N_source^0.8.
- **Shaw, Croft & Metcalf 2023**, MNRAS 519, 5236 (arXiv:2209.04564, 2023MNRAS.519.5236S). Hydro sims + ray tracing. Non-Gaussianity of the forest reduces S/N by ~2.7x (noise-free), ~1.5x (S/N=1 spectra); ray-traced potentials another ~1.3x. Reconstructed potentials biased 5-25% (nonlinearity) + 20-30% (ray tracing); Gaussianization and bias correction shown.
- **Shaw, Croft & Metcalf 2025**, OJAp 8, 10 (arXiv:2410.20014, 2025OJAp....8E..10S). DESI-like: 50 QSO/deg2, 700k spectra, 14000 deg2, sigma_delta=0.3 per 2.8 A pixel after discarding noisiest 25%; Legendre order 5 (22 modes); detection defined as correlation of reconstructed potential with foreground galaxies (NOT CMB lensing). Forecast full DESI S/N ~3 (2-pixel), ~4 (4-pixel), ~9 noiseless. 200/deg2 -> 16 (37 noiseless); PFS-like 1150/deg2 -> 16/37; MSE-like 2000/deg2 -> 27/90.
- Not found: any other group; any eBOSS/DESI data application; any detection; any CMB-lensing cross-correlation forecast for forest lensing.

## B. Forest x CMB lensing (response/bispectrum, not deflection of the forest)
- Vallinotto, Das, Spergel & Viel 2009, PRL 103, 091304 (arXiv:0903.4171): <delta_F kappa_CMB>; S/N 9 for BOSS x Planck.
- Vallinotto, Viel, Das & Spergel 2011, ApJ 735, 38 (arXiv:0910.4125): S/N 30 (BOSS x Planck), 130 (BigBOSS x ACTPol); variance-kappa S/N 9.6 / 50.
- LoVerde, Marnerides, Hui, Menard & Lidz 2010, PRD 82, 103507 (arXiv:1004.1165): QSO magnification bias effects on forest statistics 0.1-1%.
- Doux et al. 2016, PRD 94, 103506 (arXiv:1607.03625): 5 sigma detection of <P_1D(local) kappa_CMB> BOSS DR12 x Planck 2015.
- Chiang & Slosar 2018, JCAP 01, 012 (arXiv:1708.07512): theory of the above, incl. mean-flux misestimation term.
- Karacayli et al. (DESI) 2024, PRD 110, 063505 (arXiv:2405.14988): 4.8 sigma DESI Y1 x Planck PR3 bispectrum, z_eff=2.4.
- La Posta & Schaan 2025, PRD 112, 063528 (arXiv:2405.01628): forecast S/N DESI x ACT 10.3, SO 14.8, S4 20.2.
- Alonso et al. 2018, JCAP 04, 053 (arXiv:1712.02738): DLA bias from CMB lensing (tangential).

## C. Response of forest power to long modes
- Chiang, Cieplak, Schmidt & Slosar 2017, JCAP 06, 022 (arXiv:1701.03375): separate-universe responses of P_F to density, tidal field, primordial amplitude. These responses mimic convergence/shear in a QE.

## D. Formalism: lensing reconstruction from 2D/3D fields
- Hu 2001 ApJL 557, L79 (astro-ph/0105424); Hu & Okamoto 2002 ApJ 574, 566 (astro-ph/0111606); Okamoto & Hu 2003 PRD 67, 083002.
- Zahn & Zaldarriaga 2006 ApJ 653, 922 (astro-ph/0511547): 3D source -> independent k_par screens, combined.
- Lu & Pen 2008 MNRAS 388, 1819 (arXiv:0710.1108): non-Gaussianity noise floor.
- Pourtsidou & Metcalf 2014 MNRAS 439, L36 (arXiv:1311.4484); 2015 MNRAS 448, 2368 (arXiv:1410.2533); Romeo, Metcalf & Pourtsidou 2018 MNRAS 474, 1787.
- Foreman, Meerburg, van Engelen & Meyers 2018 JCAP 07, 046 (arXiv:1803.04975): matter-bispectrum bias in 3D QE, bias hardening.
- Schaan, Ferraro & Spergel 2018 PRD 97, 123539 (arXiv:1802.05706); Schaan & Ferraro 2019 PRL 122, 181301 (arXiv:1804.06403) shear-only estimator.
- Jalilvand et al. 2019 JCAP 01, 020; 2020 PRL 124, 031101. Chakraborty & Pullen 2019 MNRAS 488, 1828. Zhu & Pen 2020 arXiv:2011.08251. Maniyar, Schaan & Pullen 2022 PRD 105, 083509; Fronenberg et al. 2024 PRD 109, 123518. Nistane et al. 2022 JCAP 06, 024. Lozano Torres & Schaefer 2022 MNRAS 512, 5135. Buncher, Holder & Hotinli 2025 OJAp 8, 30 (arXiv:2402.07988) lensing of galaxy clustering. Shen, Kokron & Schaan 2026 PRD 113, 023521 (arXiv:2507.17752).
- Bias hardening: Namikawa, Hanson & Takahashi 2013 MNRAS 431, 609; Sailer, Schaan & Ferraro 2020 PRD 102, 063517 (arXiv:2007.04325); Maniyar et al. 2021 PRD 103, 083524 (arXiv:2101.12193); Darwish et al. 2021 PRD 104, 123520 (arXiv:2007.08472).

## E. Data
- DESI 2024 IV (arXiv:2404.03001, 2025JCAP...01..124A): DR1 428,403 valid Lya-region forests (1040-1205 A) and 137,432 Lyb-region; >9,500 deg2; z_eff=2.33; ~45 forests/deg2.
- ACT DR6 lensing: Qu et al. 2024 ApJ 962, 112 (arXiv:2304.05202); Madhavacheril et al. 2024 ApJ 962, 113 (arXiv:2304.05203): 9,400 deg2, 43 sigma, A_lens=1.013+-0.023. Map-noise level and L-range not verified from abstracts.
- Planck PR4 lensing: Carron, Mirmelstein & Lewis 2022 JCAP 09, 039 (arXiv:2206.07773): GMV QE, 8<=L<=400 conservative, amplitude 1.004+-0.024 (~42 sigma equivalent; inferred, not quoted).

## Synthesis (short)
Only the Croft/Metcalf/Shaw group has worked on lensing of the forest (5 papers, 2018-2025). Gaussian forecasts: DESI kappa-power S/N ~7-29 (Metcalf+2018) but realistic DESI-like mocks with noise give S/N ~3-4 for potential x foreground galaxies (Shaw+2025), limited by pixel noise, sparse sampling and forest non-Gaussianity. No detection claimed; no data application; no CMB-lensing cross-correlation forecast or measurement; no bias-hardened forest estimator; no wide-field harmonic-space implementation. The separate P_F x kappa_CMB bispectrum line (Doux+2016, Karacayli+2024 at 4.8 sigma) measures the response of forest power to long modes, which is precisely the leading contaminant (mimics convergence) of a forest-lensing QE.
