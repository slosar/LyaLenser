# Codex review 5 (gpt-6-astra, high, 2026-09-14) — verbatim final output of task-mu15fudq-spqemj; request: codex_review_5_request.md

**Stage A is not ready.** The four recorded failures remain binding, and this review found two additional validation defects: inconsistent masking in the template audit and a missing separation cut in the independent response prediction.

I read the requested context, including PDF §5, and reviewed the specified code and tests. Checks were read-only; all scratch writes stayed under `/tmp/claude-10505/codex-review/`. The non-generating test subset passed **18 tests**, with three generation tests excluded. Iteration-5 products were still partly syncing, so I used the complete report JSON, available iteration-5 products, and saved iteration-4 products. I verified the report’s included source hashes and GATES hash against the current files, but could not certify the complete iteration-5 artifact chain.

1. **[significant] Model-shaped kernel — algebra largely correct; the claimed operator exactness is too strong.**

   The Kaiser expansion is correct:
   \[
   P_F=b_F^2(P_0+2\beta_F P_1+\beta_F^2P_2).
   \]
   The LOS pixel and resolution factors are correctly squared power windows. The discrete-grid B3 covariance and its derivative correctly describe **phase-averaged** trilinear interpolation, including aliases; the Hankel derivative uses the correct negative \(J_1\) factor and Fourier normalization. Evidence: [xi_fit.py:33–60](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/xi_fit.py:33), [grid_covariance.py:16–73](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/grid_covariance.py:16), [xi_model.py:83–104](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/xi_model.py:83).

   The implemented low-rank projection is \(PCP^T\), and a fixed projector commutes with differentiating its input covariance. However, `project_fine()` fixes one transverse separation throughout each covariance matrix, averages over radial-pair locations, assumes identical uniformly weighted forests, and sets covariance outside the 40 Mpc/h table to zero. Actual angular differentiation involves \(\chi_{\rm mid}\xi'\), with varying separation and potentially varying \(g(\chi)\) inside the projection. Multiplying the averaged projected derivative by those factors afterward is an approximation. Neither commutation nor the passing derivative test establishes its accuracy for heterogeneous DR1 forests. The test explicitly avoids testing constant removal on forests longer than the table support. Evidence: [xi_fit.py:64–92](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/xi_fit.py:64), [test_iteration5.py:31–43](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/tests/test_iteration5.py:31). The PDF itself correctly calls the geometry-averaged replacement an approximation: [main.tex:395–401](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/report/main.tex:395).

   WLS does minimize \(\sum \mathrm{den}_{ij}(\xi_{\rm model}-\mathrm{num}/\mathrm{den})^2\) over the declared \(3\le r_\perp<30,\ r_\parallel<30\) range. But accumulated pair weights are **not the true inverse covariance of correlation bins** when pairs share pixels and sightlines. The reported `chi2/dof` therefore has no calibrated goodness-of-fit interpretation. The model’s uniform within-bin averaging also differs from actual pair-weighted averaging. Evidence: [xi_fit.py:95–139](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/xi_fit.py:95).

   The same fitted table supplies both \(\xi\) for the mean field and \(\xi'\) for the response. No hidden empirical amplitude correction was found. Same-sample fitting remains a statistical dependence, but the fixed/refitted slope difference, **−0.00379 ± 0.00359**, is reassuring for this campaign. It does not establish the same bound for DR1. Evidence: [run_mock_validation.py:208–242](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/run_mock_validation.py:208), [pairs.py:129–139](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/pairs.py:129), [mock_validation.md:12](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/report/mock_validation.md:12).

2. **[significant] Dense normalization — acceptable under the frozen point-estimate gate, but not causally explained.**

   The fitted slope **0.98044 ± 0.01352** is only about 1.45 SEM below unity. The generator-table slope **0.98835 ± 0.01462** is compatible with unity. Thus the products do not establish a persistent −2% physical bias requiring a correction.

   They are also not “identical” or fully table-independent. Using the paired per-seed slopes, I obtain
   \[
   s_{\rm fitted}-s_{\rm generator}=-0.00791\pm0.00238.
   \]
   This is a small but resolved table dependence. Both tables share the projection, support truncation, and geometry approximations identified above. Evidence: [campaign4.py:186–219](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/campaign4.py:186), [mock_validation.md:48–60](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/report/mock_validation.md:48).

   I found a concrete angular-scale inconsistency: the generator uses pixel size `dx/chi_ref`, whereas truth templates use `radians(20*scale)/nx`, despite `nx` being rounded upward. On dense seed 200, summing science and junk deflections differs from the saved generator deflection by **1.57% in norm**; using the generator’s pixel size reduces this to approximately **\(4\times10^{-8}\)**. However, projecting the exact generator deflection through the saved catalogue and nominal template basis gives a common science response **0.99937**. This particular discrepancy does **not** explain −2%. Evidence: [mock.py:370–378](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/mock.py:370), [mock.py:210–226](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/mock.py:210), [campaign4.py:211–212](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/campaign4.py:211).

   The junk filter correctly completes the science windows, and gradient/curl components are fitted jointly. Trilinear sampling is already included in the phase-averaged reference; invoking it again without isolating phase or geometry effects is not an explanation. The spherical rays through an anisotropic Cartesian field remain a plausible approximation error, not a demonstrated cause. Evidence: [templates.py:74–84](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/templates.py:74), [mock.py:316–329](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/mock.py:316).

   **For DR1:** do not transfer a scalar correction from this patch. The local small-pair metric does not become a 20-degree approximation merely because the survey is large, but the global flat-map construction and Cartesian mock depths cannot validate spherical survey operators. Their response needs a separate deterministic or mock check.

3. **[significant] Parameter offsets — lag pattern supported; explanation of the parameter bias incomplete.**

   There is a real binning mismatch. For the saved pixel grid, the mean absolute lags in the first bins are approximately **0.367, 1.375, 2.475, 3.575, 4.675 Mpc/h**, rather than the model’s uniform-bin means. That can produce the reported oscillating residual pattern. Evidence: [mock.py:378](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/mock.py:378), [xi_fit.py:95–104](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/xi_fit.py:95).

   My read-only refit of dense seed 200 changed:

   | Model binning | \(b_F^2\) | \(\beta_F\) | Weighted residual sum |
   |---|---:|---:|---:|
   | Implemented | 0.0186516 | 1.46302 | 45,096 |
   | Actual-lag averaging | 0.0185777 | 1.46126 | 33,226 |

   Lag-aware averaging substantially improves residuals while barely changing parameters. Consequently, it **does not explain the full +10%/−9% offset**. This agrees with the numerical statement in [NOTES.md:1013–1020](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/NOTES.md:1013), but not with treating “degeneracy” as a complete causal diagnosis.

   The fit is correlated, not singular: for this seed, the weighted derivative-column cosine is **0.97145**, and the column-normalized Jacobian condition number is **8.31**. A correlated fit can amplify model mismatch; degeneracy alone does not explain a reproducible ensemble offset. The free three-coefficient fit barely reduces the original objective, from 45,096 to 45,076.

   The **0.8%** normalization sensitivity is established for the particular fitted-versus-generator comparison, not for arbitrary model errors. A third coefficient is not automatically needed. For DR1, actual lag/weight averaging, forest projection, resolution, and nonlinear-shape adequacy must be tested first. The Arinyo basis is callable, but the fitter varies only \(b_F^2,\beta_F\); nonlinear and LOS parameters remain fixed inputs, with approximate fiducial Arinyo values. Evidence: [xi_fit.py:116–141](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/xi_fit.py:116), [forest_power.py:8–17](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/forest_power.py:8).

4. **[significant; minor metadata defect] Mock changes — defensible toy-model repairs, not demonstrated DR1 realism.**

   Drawing sightline quasars behind the forest removes the obviously artificial placement of a quasar inside its own forest. Independent Poisson draws and named RNG streams correctly remove shared-object shot noise and the previous stream-prefix reuse. The shared diagnostic retains the same forest and random catalogue. Node-centred sampling fixes the half-cell displacement. Evidence: [mock.py:135–137](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/mock.py:135), [mock.py:445–472](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/mock.py:445), [random_streams.py:4–10](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/random_streams.py:4), [run_mock_validation.py:267–278](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/run_mock_validation.py:267).

   The 1.6× oversampling is a practical way to obtain enough candidates; random subsampling preserves their conditional positional distribution. However, masked runs can reduce the target count when candidates run short, so realized density remains selection-dependent.

   Radial smoothing before exponentiation is a reasonable regularization of this lognormal toy model. It is **not equivalent to applying redshift errors or FoG displacements to an already formed quasar population**:
   \[
   \exp(b\,S\delta)\ne S[\exp(b\delta)].
   \]
   Nor does preserving \(k_\parallel=0\) modes prove that all transverse lognormal statistics are unchanged. The accompanying test demonstrates one particular projected-field case. Evidence: [mock.py:426–438](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/mock.py:426), [test_iteration4.py:105–115](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/tests/test_iteration4.py:105).

   Disjoint catalogues remain correlated through their common density field and magnification. Moving quasars behind the slab reduces local proximity coupling; it does not mathematically eliminate every selection–forest correlation. Identical full-slab forests also do not reproduce DR1’s varying forest lengths and weights. The shared-minus-disjoint combined amplitude, **−0.136 ± 0.203**, is a useful diagnostic but not a tight bound on all such effects.

   **Minor:** `SightlineSet.zq` is filled with the slab’s upper redshift, rather than each selected quasar’s actual redshift. This is inconsistent with the newly implemented behind-slab selection and would mislead later redshift-based checks. Evidence: [mock.py:520–523](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/mock.py:520).

5. **[blocker] Template gate — the sampled statistic is incorrectly double-masked.**

   A science-band cross/auto coefficient is a sensible necessary transfer check when both maps receive identical operators. Mask mixing means it is not a complete per-band or estimator-weighted deprojection test.

   Here the operators differ. `make_bundles()` already multiplies `matched_map` by the common mask. `audit_mock()` passes that map to `band_regression()`, which multiplies it by the mask again. Continuous maps and truth receive only one mask. A soft ACT mask makes this consequential. Evidence: [run_mock_validation.py:103–108](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/run_mock_validation.py:103), [template_audit.py:78–84](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/template_audit.py:78), [template_audit.py:100–109](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/template_audit.py:100).

   Recomputing all **60 iteration-4 saved audits**, applying the common mask exactly once, gives:

   | Margin, Mpc/h | Stored sampled coefficient | Consistently masked coefficient |
   |---:|---:|---:|
   | 0 | 0.97171 ± 0.01284 | **0.99397 ± 0.01304** |
   | 150 | 0.99386 ± 0.01290 | **1.01525 ± 0.01315** |
   | 300 | 0.99716 ± 0.01305 | **1.01831 ± 0.01326** |

   All three corrected values satisfy the v5 two-SEM criterion. These are diagnostic recalculations, not a retroactive campaign acceptance.

   Therefore, the claimed Poisson-sampled edge explanation is unsupported; the masking error explains most of the margin-0 deficit. Random-catalogue normalization and completeness estimation can still introduce finite-sample effects, and the continuous audit approximates their expectation, but they have not been isolated as the cause here. The repeated coefficient in iterations 4 and 5 is also the same-seed result, not independent evidence for a systematic. Evidence: [templates.py:159–175](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/templates.py:159), [NOTES.md:998–1002](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/NOTES.md:998).

6. **[blocker] Protocol and injection corrections — neither proposed precision rule resolves the present failures; validation reuse must stop.**

   The two rows fail their absolute precision clauses: their \(t\,\mathrm{SEM}\) contributions alone are **1.240 A** and **1.195 A**. This establishes inadequate precision for the 1 A requirement, not that unbiasedness has been demonstrated.

   A **1.5 A** bound can be adopted prospectively if scientifically justified, but the current bounds **1.681 A and 1.706 A still fail it**. The alternative
   \[
   |\bar A-A_{\rm target}|+t\,\mathrm{SEM}\le t\,\mathrm{SEM}+0.3A
   \]
   reduces to an observed-residual requirement \(|\bar A-A_{\rm target}|\le0.3A\). It supplies no fixed upper limit on uncertainty; the current residuals **0.441 A and 0.510 A also fail it**. Choose the acceptable systematic bound and campaign power together, before new validation. Evidence: [validation_stats.py:5–9](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/validation_stats.py:5), [mock_validation.md:22–23](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/report/mock_validation.md:22).

   There is an actual stopping-rule mismatch: GATES specifies extension above **1 A**, while `collect()` uses **0.3 A**. Both trigger for iteration 5’s initial results, so the realized seed selection is unchanged. A single predeclared extension is preferable to repeated peeking, but ordinary fixed-\(N\) t bounds do not automatically retain their nominal coverage after data-dependent stopping. Fixed final \(N\), or calibrated sequential coverage, would be cleaner. Extension seeds enter both ensembles; those final null/recovery summaries are consequently correlated. Evidence: [GATES.md:29–34](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/GATES.md:29), [campaign4.py:325–334](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/campaign4.py:325).

   The injection does rebuild pairs after shifting coordinates, so boundary crossings are real. Removing the lower cut can provide a useful auxiliary bookkeeping test. It must not replace testing the production cut: the upper cut remains, and odd differences do not cancel every boundary contribution. Retain a production-cut test using a deterministic covariance expectation and/or controlled small-amplitude convergence, with the expected selection response specified. Evidence: [inject.py:27–54](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/inject.py:27), [campaign4.py:290–293](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/campaign4.py:290).

   Fixing demonstrated bugs is justified. Changing tolerances or test definitions after inspecting results is still validation-informed development. Renaming a campaign while reusing seeds 0–59 and 200–209 does not create a fresh holdout. Those seeds should now be regression/development material; acceptance needs previously unexamined realizations.

   **Minor protocol bookkeeping:** there is no required/report-only field, and `Stage_B_allowed` requires every row to pass. This contradicts GATES’ report-only exemption. The stated 34/38 count includes three diagnostic control rows that GATES marks report-only; literal classification gives **31/35 required passes**, with the same four failures. Evidence: [GATES.md:58–59](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/GATES.md:58), [campaign4.py:386–387](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/campaign4.py:386).

7. **[significant] Provenance — useful safeguards, but several claims exceed the implementation.**

   Completion markers check provenance and artifact hashes; interrupted phases reject changed provenance. The current report’s included source hashes and GATES hash match. I also inspected commit `560ab95`: its `campaign4.py` change is confined to fingerprint serialization, supporting the documented historical migration. Evidence: [campaign4.py:33–75](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/campaign4.py:33).

   Remaining holes are concrete:

   - The source search is three shallow globs, **not every `code/**/*.py`**. Newly imported nested code could escape hashing.
   - The serialized configuration is `Config(scale=scale)`, not `campaign_config(scale)`. It records `r_perp_min=0` and `fit_rperp_min=2`, although the campaign uses 3 and 3. Source hashing currently captures those overrides, but the configuration record is misleading.
   - Basis artifacts are checked locally but their completion digest is not pinned into downstream provenance. This dependency chain is weaker than the explicitly pinned development/freeze chain.
   - Cross-site merging requires the same freeze marker, not merely identical numerical source fingerprints: that marker includes timestamps, threads, and development digests.

   Evidence: [campaign4.py:33–47](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/campaign4.py:33), [campaign4.py:84–128](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/campaign4.py:84).

   The re-provenance utility is too permissive for reuse: it allows **any change anywhere in `campaign4.py`**, rather than one specific old/new hash pair, and tolerates small changes to any floating configuration field. It also rebuilds freeze artifact hashes without first verifying all original artifacts. The audited historical change is defensible; the utility does not enforce that narrow justification. Restrict it to an explicit migration manifest and preserve the original provenance. Evidence: [reprovenance.py:18–42](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/condor/reprovenance.py:18), [reprovenance.py:58–73](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/condor/reprovenance.py:58).

8. **[significant] Errors — encouraging campaign evidence, insufficient DR1 validation.**

   The covariance ratio **0.9617** is a standard-deviation ratio, not a full covariance-matrix comparison. Keeping `sigma_F` as a scale is correct. Evidence: [run_mock_validation.py:323–325](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/run_mock_validation.py:323), [mock_validation.md:25](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/report/mock_validation.md:25).

   Contrary to the question’s premise, the reported sparse checks use **nside 16**, with **27 regions for 58 seeds and 26 for two seeds**. The helper silently falls back from nside 8 when fewer than 30 regions exist; it does not enforce 30 regions afterward. Evidence: [run_mock_validation.py:47–51](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/run_mock_validation.py:47).

   Midpoint jackknife is a reasonable spatial resampling statistic, but adjacent regions share correlated sightlines and modes. Its simple equal-region formula also needs checking when DR1 regions have unequal information. The current jackknife holds the globally fitted correlation and templates fixed; the mock scatter check assesses that approximation only for these mocks. DR1 needs region-size/information-balance checks and covariance validation with its actual selection and continuum operations. Evidence: [run_mock_validation.py:58–74](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/run_mock_validation.py:58).

   The 100 random templates are independent Gaussian CMB-like maps on **one fixed forest**, using a fixed transfer. They are useful conditional diagnostics, not joint deprojected mocks. The planned 400 ACT simulations improve that conditional test but cannot by themselves establish the forest–quasar contribution to the deprojected covariance or the real-template mean field. This limitation is explicitly acknowledged in PDF §5.4. Evidence: [run_mock_validation.py:256–264](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/run_mock_validation.py:256), [main.tex:437–444](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/report/main.tex:437).

9. **[blocker] Additional acceptance defect and shortest path — the independent response prediction violates the production cut.**

   `response._predict()` accepts final score contributions whenever `rp <= 30`; it never applies `cfg.r_perp_min`. `prediction_catalogue()` then copies \(F\) from the production catalogue, which **does** apply the lower cut. Retaining only accepted sightline pairs does not fix this: a sightline pair can have accepted pixels above 3 Mpc/h and rejected pixels below it. Evidence: [response.py:45–72](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/response.py:45), [pairs.py:121–139](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/pairs.py:121).

   An explicit two-forest covariance check, with diagonal transverse separations **2.846, 3.000, 3.154 Mpc/h**, returned a predicted score **−99.7624** instead of the correctly cut **−84.1142**; the implementation exactly reproduced the uncut score. This demonstrates the defect, not its campaign-wide magnitude. The final score sweep must match production selection; the internal continuum-projection sums must retain the covariance support needed for projection. Existing response-agreement passes cannot certify the intended cut statistic until recomputed.

   The shortest credible route is to fix the audit mask, prediction selection, angular pixel-scale mismatch, and protocol/provenance bookkeeping; add focused regression checks; and use existing products to diagnose their impact without claiming renewed independent acceptance. Bound the remaining projection/geometry/model approximations with deterministic response calculations, including the populated \(3\!-\!10\) Mpc/h range and radial-support convergence. Then freeze scientifically justified tolerances and a powered stopping design, and run a genuinely fresh acceptance ensemble.

   The signal profile supports prioritizing these checks, but “the cut is free” is too strong: its **0.73% below 3 Mpc/h** is unmarginalized independent-pair \(F\), not measured information after correlated-pair covariance and nuisance marginalization. Evidence: [signal_profile.py:22–38](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/signal_profile.py:22), [signal_profile.py:54–65](/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser/code/pipeline/signal_profile.py:54).

**Final verdict: Stage A remains unaccepted, and Stage B should remain blocked.** The model-shaped kernel is a substantial improvement, and the corrected template audit is encouraging. Acceptance nevertheless requires repairing the validation defects, a defensible fresh-ensemble protocol, and explicit limits on transferring this mock calibration to DR1.
