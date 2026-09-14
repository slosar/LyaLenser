**Stage A remains unaccepted.** The products establish an own-template contaminant with the response switched off and strongly support a contribution from quadratic forest–density correlations. They do **not** establish that this contribution alone explains the complete bias, that all relevant configurations are non-squeezed, or that conditional-mean subtraction necessarily destroys useful S/N.

I read all requested files, including PDF §§4–5, and independently analyzed both JSON products with Python. For the campaign I used `sparse[*].fits`, `predictions`, `fixed_null`, `shared`, and `shape`, pairing realizations before calculating uncertainties. SEMs below use sample standard deviations with `ddof=1`. I also checked all 1,600 per-seed table entries in `mock_validation.md` against the JSON. The result is **46/54 rows passing, including 28/36 required rows**, with `Stage_B_allowed=false`.

The checkout advanced from `7be0c9f` to `78939eb` during the review. Its differences from the requested `a60fc89` are confined to `HANDOFF.md`, `PROGRESS.md`, and `MEMORY.md`; the reviewed science code and products are unchanged. All 39 source hashes recorded in the campaign JSON, plus the GATES and `numbers3.json` hashes, match the files. The subsequently added `longmode_diagnosis.py` is outside that frozen source inventory. I made no repository or data-directory writes, ran no campaign, and submitted no jobs. Seven focused iteration-6 tests passed with caches and scratch files confined to `/tmp`.

**1. [blocker] Question 1 — The contaminant is established; the proposed causal explanation is strongly supported but not quantitatively closed.**

The independently calculated campaign statistics are:

| Statistic, 400 paired seeds | Mean ± SEM |
|---|---:|
| Deprojected `A0_R0` | **−0.40785 ± 0.17602** |
| Deprojected `A0_R1` | −0.38495 ± 0.19216 |
| Deprojected `A1_R1` | 0.49364 ± 0.18673 |
| Deprojected response: `A0_R1 − A0_R0` | **0.02290 ± 0.07459** |
| Deprojected independent prediction | 0.10335 ± 0.04668 |
| Paired response minus prediction | **−0.08045 ± 0.05255** |
| Deprojected lensing, response off | **0.89133 ± 0.05750** |
| Deprojected lensing, response on | **0.87859 ± 0.05907** |

Thus the response-off additive term and the lensing response account arithmetically for the combined recovery. This is not evidence that the linear-response prediction itself is exact.

| Template | Own-template `A0_R0` | Fixed-template `A0_R0` | Paired response | Prediction | Paired response minus prediction |
|---|---:|---:|---:|---:|---:|
| CMB | 0.32685 ± 0.13551 | −0.01212 ± 0.13486 | 1.88247 ± 0.06162 | 1.63087 ± 0.03842 | **0.25160 ± 0.04236** |
| Matched | 10.32192 ± 1.58037 | 0.56500 ± 1.67472 | 28.72513 ± 0.66861 | 23.71771 ± 0.40416 | **5.00742 ± 0.49724** |

The response excesses are therefore **15.43% and 21.11%**, respectively. Their small deprojected residual supports effective cancellation at the current precision; it does not prove exact cancellation.

Evidence: per-seed records begin at [mock_validation.json:2060](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/report/mock_validation.json:2060), with predictions at [3303](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/report/mock_validation.json:3303), fixed nulls at [3759](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/report/mock_validation.json:3759), and the reported acceptance rows at [mock_validation.md:13](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/report/mock_validation.md:13). The implementation of the paired row is [run_mock_validation.py:345–350](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/run_mock_validation.py:345).

The requested zero-intercept regression, using `A0_R0` on all 400 seeds, gives

\[
A_{\rm dep}\simeq1.040213\,A_{\rm CMB}-0.0662338\,A_{\rm matched}.
\]

Its RMS residual is **0.59711 A**. Applied to the two mean fields, it predicts **−0.34366 A**, versus −0.40785 measured. With an intercept, the coefficients become \((1.042917,-0.065574)\), with intercept −0.071882. This reproduces NOTES’ numerical argument as an approximation. It is not an exact estimator identity: the three fits use different Wiener filters, response matrices, and nuisance marginalizations. Those differences are visible in [run_mock_validation.py:118–150](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/run_mock_validation.py:118). The regression is explanatory evidence, not an independently predicted bias subtraction.

For the long-mode diagnostic, I verified every `standard` mean and SEM against the corresponding 96 original campaign seeds, and checked all reported standard-minus-subtracted mean identities. The supplied summary contains aggregate statistics, not the 96 individual subtracted or m-only fits; those latter SEMs cannot be independently reconstructed from this file alone.

| Long-mode statistic, 96 seeds | Standard | Subtracted | Standard − subtracted |
|---|---:|---:|---:|
| Deprojected response-off null | −0.96726 ± 0.37163 | **−0.53871 ± 0.21759** | −0.42855 ± 0.27621 |
| Matched response-off null | 6.40862 ± 3.08902 | 3.65771 ± 1.94226 | 2.75091 ± 2.39025 |
| Truth lensing response, response on | 0.94571 ± 0.03887 | 0.52332 ± 0.03118 | 0.42239 ± 0.02061 |
| Deprojected lensing response, response on | 1.04590 ± 0.10906 | 0.59830 ± 0.09988 | 0.44761 ± 0.06775 |

The m-only scores are **matched \(11.54066\pm2.29134\)** and **deprojected \(-1.12431\pm0.29145\)**. This establishes that the constructed m-only statistic has a substantial own-template correlation. However:

- The deprojected subtraction-induced null change is only **1.55 SEM**, and the residual remains **2.48 SEM below zero**.
- The m-only score is not equal to the amount removed. Cross terms, the remaining field, and baseline estimation matter.
- The response excess has **not** been shown to disappear. The diagnostic reuses the original prediction. After subtraction, response-minus-old-prediction is **−1.05581 ± 0.07041** for CMB and **−15.28300 ± 0.82482** for matched. NOTES correctly acknowledges that these are not prediction tests.
- The deprojected response remaining near zero is an empirical result, not “unchanged by construction.”

Evidence: [longmode_iteration6_summary.json:306–410](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/report/longmode_iteration6_summary.json:306), [131–147](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/report/longmode_iteration6_summary.json:131), [216–253](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/report/longmode_iteration6_summary.json:216), and [longmode_diagnosis.py:95–107, 133–136](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/longmode_diagnosis.py:95).

**Alternative explanations that fixed-template rows would not expose remain concrete.** I found no `fit_save` mean-field bookkeeping error affecting the reported amplitudes: it computes and saves the proper fit before creating a separate zero-mean-field copy solely for the `raw` diagnostic. The unchanged partial mean fields in that temporary copy would make its jackknife inconsistent, but only its amplitude is retained. See [run_mock_validation.py:65–80, 213–216](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/run_mock_validation.py:213) and [amplitude.py:98–112](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/amplitude.py:98).

By contrast, the actual template uses a quasar-derived radial normalization, finite random-catalogue completeness, RSD transport, magnification, and a realization-dependent Wiener denominator. Sightline locations are another draw from the same lognormal density. These are correlations absent from the ideal fixed-geometry Gaussian argument and disrupted by fixing an independent template. They are supported implementation channels, although their individual contributions have not been measured. See [templates.py:159–174](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/templates.py:159) and [mock.py:434–480](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/mock.py:434).

A deterministic ACT mask or completeness pattern alone cannot create an odd Gaussian moment. Their interaction with selection, estimated normalization, and model error can matter. Likewise, independent Poisson sampling adds no mean to a statistic linear in a correctly normalized tracer at fixed geometry; that argument does not cover this complete nonlinear construction.

The shared-minus-disjoint response-off deprojected shift is **0.10924 ± 0.06896 A**. It limits the effect of adding the same objects under this construction, not the shared-density selection channel present in both versions. The six-bin result, independently reproduced as **\(p=0.00223718\)**, establishes shape inconsistency; no quantitative prediction of the m-only six-bin vector is supplied to establish unique shape-level closure.

**2. [significant] Question 2 — The missing contractions are real, but “non-squeezed bispectrum” is an incomplete description of why deprojection fails.**

For fixed geometry, weights, covariance kernel, and linear template operator, write the score schematically as

\[
y_T=\sum_{p<q}K_{pq}\,\mathcal L_{pq}[T]\,
       (d_pd_q-C_{pq}).
\]

Its mean is a weighted, projected forest–forest–template three-point function. The campaign amplitude additionally involves a random \(F^{-1}\), fitted \(\xi\), and estimated map filters; consequently, an ensemble score formula alone does not predict its exact mean amplitude. This distinction is explicit in PDF §5.2, [main.tex:404–417](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/report/main.tex:404).

**(i) Gaussian forest, ideal lognormal tracer, response off.** Let \(f_p,f_q,s_x\) be jointly zero-mean Gaussian and

\[
T_x=e^{b_qs_x-b_q^2\langle s_x^2\rangle/2}-1.
\]

Exponential tilting gives

\[
\left\langle f_pf_qT_x\right\rangle
=b_q^2 C_{ps}C_{qs}.
\]

This identity is exact for that ideal tracer. Fixed continuum projectors simply act on the two forest legs. For the kernel-matched map, one must also apply the radial integral, angular operators, and its **\(1/b_q\)** normalization; the displayed \(b_q^2\) coefficient is for the quasar overdensity before that division.

It is not an exact formula for the complete production template, which also has realization-estimated normalization, exponent clipping, RSD displacements, magnification, Poisson sampling, and an estimated Wiener filter. Nor is the diagnostic’s fitted \(cP\delta_L\) the exact Gaussian conditional mean. The relevant generation is [mock.py:43–89, 101–139, 434–449](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/mock.py:43); the diagnostic fits \(c\) from the same forest at [longmode_diagnosis.py:75–85](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/longmode_diagnosis.py:75).

**(ii) Modulation on.** Put \(a=R/2\), \(u=\delta_L(p)\), \(v=\delta_L(q)\), and use raw forest variables before continuum projection. For the ordered variables \((f_p,f_q,u,v)\), define \(h_i=b_q\langle y_i s_x\rangle\). Then

\[
\begin{aligned}
\langle f_pf_q(1+au)(1+av)T_x\rangle
={}&h_ph_q\\
&+a\sum_{t=u,v}
 \left(C_{pq}h_t+C_{pt}h_q+C_{qt}h_p+h_ph_qh_t\right)\\
&+a^2\left[
 \sum_{\{i,j\}\subset\{p,q,u,v\}}
 C_{ij}\!\!\prod_{k\notin\{i,j\}}\!\!h_k
 +h_ph_qh_uh_v
 \right].
\end{aligned}
\]

The sum contains six covariance–two-mean terms. This follows directly by shifting the mean of a Gaussian under exponential tilting. I checked the polynomial independently with five-dimensional Gaussian quadrature; the numerical agreement was \(4\times10^{-15}\) relative for a positive-definite test covariance.

For a *Gaussian* template \(\tau\), the response-off term vanishes at fixed geometry, but

\[
\langle f_pf_q(1+au)(1+av)\tau\rangle
=a\sum_{t=u,v}
\left[
 C_{pq}C_{t\tau}
 +C_{pt}C_{q\tau}
 +C_{qt}C_{p\tau}
\right].
\]

The independent prediction contains the first contraction. The other two are precisely the missing forest–long-mode contractions identified in NOTES. The \(a^2\) term vanishes here by Gaussian odd-moment symmetry. Projection must be applied to the complete expression; some spatially constant pieces can be annihilated.

An exact conditional formulation makes the error especially clear. For a specified long-mode vector \(\ell\),

\[
\mu=C_{f\ell}C_{\ell\ell}^{-1}\ell,\qquad
\Sigma=C_{ff}-C_{f\ell}C_{\ell\ell}^{-1}C_{\ell f}.
\]

The conditional second moment is \(M(\ell)=\Sigma+\mu\mu^T\). Modulation changes it by

\[
P[D M(\ell)D-M(\ell)]P^T,
\]

whereas production predicts \(P[D C_{ff}D-C_{ff}]P^T\). This explains both the response-off conditional term and additional response terms. It does not by itself establish that the entire measured 15–21% excess comes from them; radial-support and geometry approximations also remain. See [response.py:32–78](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/response.py:32).

**(iii) Real-universe tree level.** In real space, with the convention

\[
\delta_q=b_1\delta+\frac{b_2}{2}[\delta^2]+b_{K^2}[K^2]+\cdots,
\]

and analogous forest coefficients, the density-only contribution is

\[
B_{FFq}
=b_{F1}^2b_1B_m
+b_{F1}^2b_2P_1P_2
+b_{F1}b_{F2}b_1(P_1P_3+P_2P_3)
+\text{tidal terms},
\]

where

\[
B_m=2F_2(\mathbf k_1,\mathbf k_2)P_1P_2+\mathrm{cyc.},
\quad
F_2=\frac57+\frac{\mu}{2}
\left(\frac{k_1}{k_2}+\frac{k_2}{k_1}\right)+\frac27\mu^2.
\]

For observed redshift-space forests, the appropriate expression is

\[
\begin{aligned}
B_{FFq}={}&
2Z_{1F}(1)Z_{1F}(2)Z_{2q}(1,2)P_1P_2\\
&+2Z_{1F}(2)Z_{1q}(3)Z_{2F}(2,3)P_2P_3\\
&+2Z_{1F}(1)Z_{1q}(3)Z_{2F}(1,3)P_1P_3.
\end{aligned}
\]

Here the forest kernels must include flux nonlinearity, velocity-gradient dependence, and the relevant observational operators. A nonlinear power-spectrum fit alone does not specify them. The quasar-leg \(b_2P_1P_2\) and tidal structure are also given in [Lazeyras et al., Eq. 4.3](https://arxiv.org/pdf/1511.01096); the forest squeezed-response treatment is developed in [Chiang et al.](https://arxiv.org/abs/1701.03375).

**What cancels.** The report cancels contributions satisfying the assumed proportionality between the CMB and matched-tracer intrinsic bispectra, including arbitrary short-mode response coefficients when the radial kernels match. A quasar \(b_2\) or tidal contribution need not satisfy that proportionality.

However, if the matched map were an exact copy of the CMB slab-density contribution, their intrinsic correlations with *any* forest statistic would cancel, without requiring a squeezed expansion. Therefore gravitational \(F_2\) and forest nonlinearity do not automatically survive subtraction: their common matter-tracer part can cancel. Nonlinear/stochastic tracer mismatch, unequal projection operators, and nonlocal radial correlations are the relevant failure conditions. The report itself states the linear-tracer assumptions at [main.tex:295–339](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/report/main.tex:295).

**Scale separation.** With the recorded \(\chi_{\rm ref}=3955.72\,h^{-1}\mathrm{Mpc}\),

\[
k_L=L/\chi_{\rm ref}=0.0101\text{–}0.0758\,h\,\mathrm{Mpc}^{-1},
\qquad
k_Lr_\perp=0.030\text{–}2.28
\]

over the stated ranges. Some configurations are well squeezed; high-\(L\), wider pairs are not. At \(L=300,r_\perp=30\), the pair-response sinc factor is approximately 0.80, as the PDF states. “Not uniformly squeezed” is correct; “nothing is squeezed” is not. Full transverse wavelengths \(2\pi/k_L\) are approximately **83–621 Mpc/h**; a quoted “mode scale” needs its convention specified. See [main.tex:393](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/report/main.tex:393).

The quasar \(b_2P_1P_2\) term can also survive in a squeezed triangle: nonlinear small-scale structure contributes to the tracer’s long Fourier mode. Its existence is not synonymous with poor geometric scale separation.

**Size for DR1.** The stated \(b_2\simeq4.6\) is reproduced by the separate-universe halo fit

\[
b_2=0.412-2.143b_1+0.929b_1^2+0.008b_1^3,
\]

which gives **4.63475** at \(b_1=3.5\). This is a halo fitting relation, not a unique prediction from the phrase “peak-background split”; a Press–Schechter PBS calculation with Eulerian conversion gives approximately **5.016**. See [Lazeyras et al., Eq. 5.2 and Appendix B](https://arxiv.org/pdf/1511.01096).

Thus, holding everything else fixed, the isolated local-quadratic coefficient would be roughly **0.38–0.41** of the mock’s \(b_q^2=12.25\). **The total DR1 contaminant cannot be inferred by rescaling −0.41 A by this factor.** Quasar occupation, tidal bias, nonlinear forest kernels, redshift-space effects, selection, and estimator weighting are not supplied. The products establish a mechanism requiring treatment, not a quantitative DR1 bias prediction.

**3. [blocker] Question 3 — Prefer an explicit conditional pair-moment or bispectrum model, with nuisance marginalization; no proposed shortcut has yet demonstrated acceptance.**

My preferred next development step is to model the **conditional pair second moment** and its uncertainty, retain the existing linear matched subtraction where valid, and fit the additional contaminant shapes jointly with lensing. A quadratic nuisance is a sensible first component, but its sufficiency must be demonstrated.

The actual PDF Table 2 gives:

| Forecast | Relevant S/N |
|---|---:|
| DR1 empirical, matched, \(40\le L\le300\) | **1.0 ACT; 0.9 Planck** |
| DR1 empirical, response | **2.0** |
| DR1 empirical, all \(L\), naive / BH / matched | **1.3 / 0.8 / 1.2** |
| Complete-DESI example, 50 deg\(^{-2}\), \(P_N=0.17\), science band | **4.2 ACT; 3.9 Planck** |

These are the current table values, rather than older 0.9/1.8 figures in historical text. Evidence: PDF Table 2 and [numbers3_table.tex:12, 21](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/report/numbers3_table.tex:21).

For a lensing vector \(s\), nuisance columns \(N\), and covariance \(C\), the information after fitting nuisances is

\[
I_A=s^TC^{-1}s-
s^TC^{-1}N(N^TC^{-1}N)^{-1}N^TC^{-1}s.
\]

For one extra nuisance, the S/N retention is \(\sqrt{1-\rho^2}\), with \(\rho\) measured in this covariance metric. The current products do not contain the necessary new nuisance vectors and covariance to assign a numerical penalty.

| Candidate | Verdict, requirements, and decisive Stage A test | S/N consequence |
|---|---|---|
| **(a) Additional quadratic nuisance** | Promising. Construct a pair-covariance or bispectrum shape over separation, \(L\), and redshift. A deflection template made from a squared 2D quasar map is not automatically that shape. It requires tracer bias/cross-covariance information, radial selection, and shot-noise treatment, preferably split-catalogue products. Test recovery while independently varying its amplitude, quasar nonlinearity, and lensing. | Determined by the measured nuisance–lensing overlap above. No defensible fixed percentage is available. More shapes may be needed. |
| **(b) Subtract predicted forest means** | Not ruled out in principle. It requires a sufficiently accurate three-dimensional conditional mean, residual covariance, and recomputed lensing/response operators. Fit predictive coefficients outside the acceptance ensemble. The present scalar \(cP\delta_L\) subtraction is not such a model and fails its null. | Losing 45% of raw response does **not** imply losing 45% of S/N; see the calculation below. A data-usable predictor will also have quasar shot noise and reconstruction uncertainty absent from this oracle diagnostic. |
| **(c) Restrict \(r_\parallel\) or \(k_\parallel\)** | Useful diagnostic, not an established fix. A short-\(r_\parallel\) cut is not a high-\(k_\parallel\) filter. Explicit Fourier filtering requires propagating forest windows, continuum projection, and noise. Test the retained bias and response with the complete operator. | The recorded profile places about 80% of independent-pair information at \(r_\parallel<5\). Discarding it would retain only about \(\sqrt{0.2}=0.45\) of that approximate S/N. Keeping it does not remove the observed contaminant. A high-\(k_\parallel\) restriction has a different, currently unmeasured cost. |
| **(d) Condition the response prediction on the long field** | Required for a credible prediction, using \(M(\ell)\), not unconditional \(C\). It is a diagnosis/calibration improvement until the corresponding additive term is removed or marginalized in the estimator. Test on/off differences and the response-off conditional score separately. | A known correction has no automatic projection penalty; uncertainty in the correction adds variance or nuisance degeneracy. Correcting the prediction alone leaves the estimator bias intact. |
| **(e) Conditional pair-moment subtraction or improved density template** | Preferred development path. Subtract \(M_0(Q)=E[dd^T\mid Q,A=0]\), including its covariance component, or construct a density template corrected for nonlinear tracer bias. Keep its lensing response explicit. First test an oracle slab-density template, then continuous intensity, sampled quasars, and finally the data-usable construction. | Potentially retains information removed by field subtraction. Its actual cost depends on model uncertainty and additional tracer noise; it must be measured. |

The recorded signal-profile statement is at [PROGRESS.md:75](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/PROGRESS.md:75); its interpretation as independent-pair information, rather than full marginalized information, was already emphasized in [review 5:139](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/report/reviews/codex_review_5.md:139).

**The long-mode S/N argument needs correction.** Using the same 96 realizations and the response-on lensing difference:

| Quantity | Standard | Subtracted |
|---|---:|---:|
| Lensing response | 1.04590 | 0.59830 |
| Scatter of combined amplitude | 3.88932 | 2.16903 |
| Scatter divided by lensing response | **3.71863** | **3.62535** |
| Response-off null divided by that response | −0.92481 | −0.90041 |

These point estimates show approximately unchanged calibrated precision **and approximately unchanged calibrated additive bias**. They neither demonstrate an S/N improvement nor justify abandoning every conditional-mean approach. They show that this particular subtraction does not solve the problem. Evidence: [longmode_iteration6_summary.json:306–410](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/report/longmode_iteration6_summary.json:306).

For a genuinely jointly Gaussian forest and observed tracer, optimal linear conditioning gives

\[
C_{f\mid Q}=C_{ff}-C_{fQ}C_{QQ}^{-1}C_{Qf},
\]

or \(P_{\rm residual}=P_F(1-r^2)\) mode by mode in the homogeneous scalar case. That formula describes residual covariance. It is **not** a universal multiplicative correction to a nonlinear-tracer bispectrum bias. A lognormal/Poisson tracer needs the conditional second moment, not just a scalar correlation coefficient.

**Mock realism:** retain the present mock as a valuable analytic regression test, but add a scientifically relevant validation tier with independently variable quasar \(b_2\)/tidal bias and a nonlinear, redshift-space forest. Matching its two-point power is insufficient; validate its forest–forest–quasar/CMB bispectra. Merely lowering \(b_2\) until the current null becomes insignificant would weaken the test rather than fix the estimator. The present forest is explicitly a linear filter at [mock.py:73–76](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/mock.py:73).

**4. [significant] Question 4 — The sparse normalization deficit is suggestive, but neither its persistence nor its cause is established.**

The independently calculated response-off lensing differences are:

| Template | Mean ± SEM | Distance below unity |
|---|---:|---:|
| Truth | **0.95870 ± 0.01889** | 2.19 SEM |
| CMB | **0.88674 ± 0.04785** | 2.37 SEM |
| Deprojected | **0.89133 ± 0.05750** | 1.89 SEM |

Crucially, the same-seed difference is

\[
s_{\rm CMB}-s_{\rm truth}
=\mathbf{-0.07196\pm0.04575},
\]

only **1.57 SEM**. The products therefore do not establish a distinct additional CMB normalization failure at high confidence. The truth deficit warrants investigation, but two correlated 2-SEM discrepancies do not identify its mechanism.

The spectra check is reassuring but insufficient: it measures broad-band, unmasked **signal-map** spectra. Its ratios are **0.997964, 0.988202, 1.018194**; it does not test the masked noisy-map conditional transfer or the pair-weighted response. See [mock.py:264–280](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/mock.py:264).

The code provides these concrete diagnostic branches:

| Suspect | Assessment and discriminating check |
|---|---|
| **Estimated Wiener transfer** | `make_bundles` estimates \(C^{XX}\) separately from each actual noisy map and uses \(S/\widehat C^{XX}\). This is not the exact Gaussian conditional coefficient \(S/C^{XX}_{\rm ensemble}\). Nonlinear ratios and weighting across annuli can introduce a mean bias; its sign and size are not established here. Compare fixed ensemble transfer, independently estimated transfer, and production transfer using a known covariance response. |
| **Forest noise × fitted kernel** | Independent noise on distinct sightlines has zero mean cross-product for fixed weights/kernel. It cannot simply explain a persistent bias by being “noise.” However, the same noisy sample supplies fitted \(\xi\), its derivative, and a nonlinear normalization. Hold geometry and weights fixed while comparing generator, independent-fit, and same-sample kernels and paired noise realizations. |
| **The \(r_\perp\ge3\) cut** | Correctly implemented in both score paths for this campaign. A cut can change sensitivity to covariance-model errors; it does not inherently produce a multiplicative bias. Test the exact projected covariance derivative on the accepted pixels, particularly 3–10 Mpc/h. |
| **Projection/binning/geometry** | Still plausible. The fitted basis uses geometry-averaged projection, finite support, and uniform within-bin averaging. Those approximations are shared by successful tests and need direct accuracy bounds. |

Evidence: [run_mock_validation.py:90–147, 218–222, 247–263](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/run_mock_validation.py:90), [xi_fit.py:68–104, 116–140](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/xi_fit.py:68), and [pairs.py:124–147](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/pairs.py:124).

The dense comparison is not a controlled noise-only experiment: it also changes sightline density, weighting, mask/completeness/magnification settings, and template choice. See [campaign4.py:218–236](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/campaign4.py:218).

Dense fitted-versus-generator slopes differ by **−0.01007 ± 0.00163**, despite passing the 5% ratio gate. This is a small, resolved table dependence. Also, NOTES’ attribution of the change from iteration 5’s dense 0.980 to iteration 6’s 1.001 specifically to patch-side repair is not demonstrated: seeds changed, and review 5 measured that angular-scale effect’s common response as 0.99937. See [NOTES.md:1101–1103](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/NOTES.md:1101) and [review 5:25–31](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/report/reviews/codex_review_5.md:25).

**5. [significant] Question 5 — The principal production fixes are implemented correctly; “every review-5 fix is complete” is too broad.**

| Fix | Verdict and evidence |
|---|---|
| Single-mask sampled audit | **Correct.** The unmasked matched map enters the band regression, which applies the mask once. Recomputed margin-0 coefficient: **1.008444 ± 0.005514**. [template_audit.py:107–111](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/template_audit.py:107). |
| Production selection in prediction | **Correct for the campaign’s full-slab configuration.** Final accepted pixels obey the lower/upper transverse and LOS cuts; projection sums retain their separate covariance support. [response.py:32–64](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/response.py:32). |
| Generator patch side | **Correct in the reviewed template/audit/spectrum consumers.** [mock.py:256–261](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/mock.py:256), [run_mock_validation.py:101–103](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/run_mock_validation.py:101), [campaign4.py:227–229](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/campaign4.py:227). |
| Actual `zq` | **Correct.** Uses each selected sightline quasar’s observed displaced redshift. [mock.py:528–531](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/mock.py:528). |
| Seed roles and paired variants | **Correct.** 40 `full`, 360 `core`; all 400 contain both physical amplitudes with response on/off. Old acceptance ranges are refused. [run_mock_validation.py:225–260](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/run_mock_validation.py:225), [campaign4.py:261–275](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/campaign4.py:261). |
| Required flag and paired row | **Correct in collection.** Required rows alone control acceptance. The sparse slope tolerance is prospectively converted to the specified 2-SEM rule by `collect`. [campaign4.py:367–370, 408–416](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/campaign4.py:367). |
| Jackknife bookkeeping | **Correctly recorded.** All 400 use nside 16; 380 have 27 regions and 20 have 26, mean **26.95**. This is bookkeeping repair, not validation of DR1 region independence. [run_mock_validation.py:49–58, 352–357](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/run_mock_validation.py:49). |
| Fingerprint and basis dependency | **Substantially corrected.** Recursive source inventory, actual campaign configuration, seed ranges, and downstream basis digest are present. [campaign4.py:40–58, 129–143, 195–196](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/campaign4.py:40). |
| Injection expectation | **Valid as a bookkeeping gate.** True/unshifted separations supply the covariance while shifted geometry supplies cuts, kernel, and mean field. [pairs.py:124–147](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/pairs.py:124), [inject.py:39–64](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/inject.py:39). |

There are remaining implementation qualifications:

**Provenance migration is not fully repaired.** The utility requires exact source-hash pairs, verifies listed completed artifacts before rewriting, and preserves original marker provenance. But its float exception still accepts any relative change below \(10^{-9}\), without naming it in the manifest—not exact 12-significant-digit canonicalization. A scratch execution of its pure checking functions accepted `r_perp_min: 3.0 → 3.000000001` with an empty manifest.

It also omits basis completion/artifact verification from its migration inventory and does not migrate basis provenance. A genuine fingerprint change can therefore leave the basis inconsistent with the rewritten downstream products. These defects are not evidence that this campaign was improperly migrated, but the utility should not be certified for reuse. Evidence: [condor/reprovenance.py:33–49, 72–78, 89–108](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/condor/reprovenance.py:33).

**Selection equivalence is not general to tomography.** Production honors `cfg.slab_index` through a midpoint cut; the predictor checks valid pixel slabs but has no equivalent selected-slab argument. This does not affect iteration 6’s `slab_index=-1`, but it must be resolved before using the predictor for tomographic catalogs. Compare [pairs.py:127–130, 168](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/pairs.py:127) with [response.py:55, 70–73](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/response.py:55).

**Some generator uses of the nominal side remain.** Outside-box fields, CMB-noise pixel area, edge-distance scaling, and random-catalogue support still use nominal `side` or `20*scale`; this leaves a small internal geometric inconsistency, not an explanation for an 11% deficit. See [mock.py:413–419, 470, 543](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/mock.py:413).

**The expectation injection is an appropriate replacement for the old single-realization measured gate, with a limited claim.** Its smallest-amplitude slope is **0.978840**, with 0.978336, 0.976156, and 0.966307 at larger amplitudes. It tests signs, cuts, pair rebuilding, and basis representation against the supplied covariance model. It does not independently test whether that model is the physical covariance, nor whether conditional forest–template correlations are correctly modeled. The extrapolation toward small amplitude is near 0.979, not demonstrated exact unity.

The full-slab prediction and injection regression tests passed here: **7/7 in 64.82 seconds**. The campaign’s saved test result is **64 passed**, rather than the older 57-test count in the handoff. Evidence: [mock_validation.md:30–32](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/report/mock_validation.md:30), [test_iteration6.py:21–93](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/tests/test_iteration6.py:21). Review 5’s broader concerns about geometry-specific projection, support convergence, and realistic covariance remain relevant.

**6. [significant] Question 6 — The fixed-ensemble reading is legitimate; prospective power and acceptance logic still need improvement.**

Using all 400 seeds for both null and recovery is explicitly predeclared. The collector has no data-dependent extension and requires the complete fixed ranges. The paired normalization row is also predeclared. I found no stopping-rule mismatch or validation-informed change to the frozen round-6 acceptance result. Evidence: [GATES.md:36–50, 79–80](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/GATES.md:36), [campaign4.py:345–359](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/campaign4.py:345).

The eight failed required rows are correlated manifestations, not eight independent detections. Conversely, the additive null, linear-response excess, and multiplicative response should not be collapsed into a single proven mechanism. Post-campaign regressions and the 96-seed subtraction are legitimate development diagnostics, but their significance is exploratory and the 96 seeds are a subset of the 400, not independent confirmation.

The present 95% residual bounds reproduce as:

\[
\begin{array}{ll}
A0\_R0:&0.75390\,A,\\
A0\_R1:&0.76272\,A,\\
A1\_R1-1:&0.87346\,A.
\end{array}
\]

All exceed 0.5 A.

**N=400 is feasible but weakly powered for acceptance.** With per-seed scatter 3.7 and zero true bias,

\[
P\!\left(|\bar A|+t\,{\rm SEM}\le0.5\right)
\approx
2\Phi\!\left(\frac{0.5\sqrt N}{3.7}-t\right)-1,
\]

which gives approximately **54% at N=400** for a single such gate. Rough normal-approximation requirements are **576 seeds for 80%** and **712 for 90%** power, before accounting for the full correlated acceptance procedure or additional nuisance costs. Merely observing that \(t\,{\rm SEM}<0.5\) does not establish adequate power. The implemented bound is [validation_stats.py:5–9](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/validation_stats.py:5).

GATES v7 should require:

- A frozen estimator, nuisance model, transfer estimator, cuts, bias conventions, and development-derived calibration, followed by genuinely fresh acceptance realizations.
- Separate absolute own-template null tests, paired response-minus-complete-prediction tests, paired lensing normalization for **truth, CMB, and deprojected templates**, and absolute combined recovery. Neither a null nor a paired slope substitutes for the other.
- Scientifically justified **equivalence bounds with uncertainty**, rather than only failure to reject unity within 2 SEM. The matched slope currently passes despite a negative central value and an uncertainty of 1.47; that is consistency, not normalization certification.
- Mechanism tests that vary tracer nonlinearity independently of lensing and response, plus oracle-density, continuous-intensity, sampled-template, and selection controls. Preserve the difficult current mock.
- Shape tests against the complete predicted lensing-plus-nuisance vector, using joint covariance; scalar closure alone is insufficient.
- Known-covariance tests of the actual production selection, heterogeneous projection, radial-support convergence, and noisy-map transfer. Keep expectation injection as a bookkeeping test.
- A power calculation for the **whole correlated gate set**, including the revised estimator’s scatter. Define the handling of diagnostic p-values prospectively. An intersection of equivalence tests does not automatically require Bonferroni correction, but its probability of accepting a good pipeline still needs calculation.
- Immutable production provenance and separately fingerprinted diagnostic provenance, with complete per-seed outputs sufficient to reproduce every quoted uncertainty.

These requirements make another iteration decisive by separating estimator adequacy from an underpowered acceptance lottery.

**7. [blocker] Question 7 — Useful DR1 statistics are possible, but a trustworthy lensing upper limit is also blocked by the current estimator bias.**

Under the current protocol, Stage B remains blocked: [IMPLEMENTATION.md:11](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/IMPLEMENTATION.md:11) and [GATES.md:79–87](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/GATES.md:79).

Scientifically, a separately scoped diagnostic analysis could measure the raw forest-pair × CMB statistic, the matched-quasar statistic, band/separation dependence, and the per-skewer \(P_{1D}\) × CMB response statistic. However, the matched-quasar amplitude is now demonstrably a mixture of response, tracer nonlinearities, and other intrinsic terms. It cannot be labeled a pure measurement of \(R_\delta\) using the present prediction.

The current Table 2 suggests only **approximately 1σ lensing sensitivity and 2σ response sensitivity** for empirical DR1 inputs. ACT and Planck comparisons are useful, but their measurements share the forest, quasar template, and much of the sky signal; agreement does not validate away the same intrinsic contaminant.

**A lensing upper limit must wait for a calibrated estimator or an explicit contaminant likelihood with credible uncertainty.** A negative additive bias can make an upper limit artificially tight. Calling the result an upper limit does not remove the normalization or bias requirement, and a scalar subtraction of the mock’s −0.41 A is not transferable to DR1.

Additional release blockers are the realistic survey response and covariance checks already required by the report: varying forest lengths and weights, continuum and resolution operations, curved-sky map/template operators, selection and magnification, tomographic disjoint selections, and actual ACT/Planck transfer functions. The present predictor’s selected-slab omission is relevant here. Random CMB templates test a conditional component of the error; they cannot validate the real quasar-correlated mean field or full deprojected covariance. These limitations are stated at [main.tex:332–339, 437–450](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/report/main.tex:332).

**Final verdict:** The response-off own-template contaminant is established, and quadratic forest–density correlations are a demonstrated contributor. The diagnosis is not yet a complete quantitative explanation of the campaign bias or a prediction for DR1. The preferred fix is a conditional pair-moment/bispectrum model with validated nuisance marginalization or a correspondingly improved density template, followed by full response and covariance recalibration. The present long-mode subtraction fails, but its loss of raw response does not prove an S/N prohibition on all conditional methods. **Stage A could pass after one more substantive development-and-validation iteration; the current evidence does not justify promising that it will. Stage A is not passed, and Stage B remains blocked.**

Codex session ID: 01a0a0c0-a6b4-7d30-ac4d-cfb34bcdcda3
Resume in Codex: codex resume 01a0a0c0-a6b4-7d30-ac4d-cfb34bcdcda3
