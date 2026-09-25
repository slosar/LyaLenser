# Codex adversarial review 2 (gpt-6-astra, high), estimator only — 2026-09-10

Scope: report Sec. 5, Sec. 4.2 kernel matching, PLAN.md Sec. 0/6, mathematica/pair_template_checks.wls. Verbatim output follows.

---

I read the requested sections and previous findings, ran checks 1–8 with Mathematica, and checked several counterexamples independently. No files were modified.

The central response formula is correct: **within the stated flat-sky, homogeneous, single-source-plane model**, Eq. `Rij` is exact at first order for arbitrary \(L\). The finite-separation sinc error is removed. There is also no missing factor of two in Eq. `Ahat` with unique, distinct-sightline pairs. The problems concern its interpretation, calibration, covariance, and deprojection.

1. **Blocker — The specified normalization estimates a template coefficient, not the claimed cross-spectrum amplitude.**  
   [main.tex:394](/home/anze/Dropbox/work/LyaLenser/report/main.tex:394), Eqs. `Rij`, `Ahat`; [main.tex:405](/home/anze/Dropbox/work/LyaLenser/report/main.tex:405).

   Define \((U,V)_w=\sum_{i<j}w_iw_jU_{ij}V_{ij}\). With fixed weights, correct baseline covariance, and an arbitrary true deflection,
   \[
   \mathbb E[\hat A\mid\alpha_{\rm true},\alpha_T]
   =\frac{(R_T,R_{\rm true})_w}{(R_T,R_T)_w}
   \]
   at first order. This is the weighted projection of the true covariance perturbation onto the template response. It equals \(A\) when \(R_{\rm true}=A R_T\); Gaussianity is unnecessary for that identity.

   For a stochastic tracer \(X=\kappa_\perp\), suppose the desired parameter scales \(C_L^{\kappa_{\rm Ly\alpha}X}=A S_L\). Under a Gaussian conditional model,
   \[
   \mathbb E[\kappa_{\rm Ly\alpha}\mid X]
   =A\,h_LX_L,\qquad h_L=S_L/C_L^{XX}.
   \]
   If the fitted template is \(f_LX_L\), its conditional expectation is instead
   \[
   A\,\frac{(R[fX],R[hX])_w}{(R[fX],R[fX])_w}.
   \]
   The report’s extra \(1/(C^{\kappa\kappa}+N_\kappa)\) in \(f_L\) therefore invalidates the stated unit normalization.

   The Fisher-optimality claim is also wrong as written. In the homogeneous, weak-cross-correlation limit, let the **unnormalized** pair score have response \(G_L\), so \(s_L=G_L\hat\kappa_L\). The usual optimal cross-spectrum weighting translates into a template filter proportional to
   \[
   \frac{S_L}{G_L\,C_L^{\hat\kappa\hat\kappa}C_L^{XX}}.
   \]
   The report omits \(G_L^{-1}\). In the ideal noise-dominated diagonal model, \(C_L^{\hat\kappa\hat\kappa}\simeq G_L^{-1}\), and these factors cancel: the appropriate template is \(h_LX_L\).

   **Fix:** choose explicitly between a conditional template-coefficient estimator and a cross-spectrum amplitude estimator. For the latter, normalize the numerator by its response to the specified fiducial cross-spectrum, including the survey window. Division by a realization-dependent \(F(X)\) additionally requires justification: cancellation of an ensemble numerator does not alone prove cancellation of the ratio. The normalized statistic is not literally a fixed-weight cubic estimator.

2. **Significant — Separate \(F_b\) values do not establish unbiased bandpowers.**  
   [main.tex:405](/home/anze/Dropbox/work/LyaLenser/report/main.tex:405), [main.tex:418](/home/anze/Dropbox/work/LyaLenser/report/main.tex:418).

   For band templates, the response generally contains
   \[
   F_{bc}=\sum_{a<b}
   \Delta\alpha_{ab,b}^{T}M_{ab}\Delta\alpha_{ab,c}.
   \]
   Disjoint harmonic bands do not make this matrix diagonal on an irregular footprint. Dividing each score by \(F_{bb}\) leaves responses to other bands. Cross-spectrum bandpowers similarly need a window/response matrix, although its definition differs from this deterministic-template matrix.

   **Fix:** specify the measured band windows and full response matrix, including out-of-band leakage. These contractions can use the existing pair catalogue; correcting this does not itself require another pixel pass.

3. **Significant — The mean-field oddness claim is false.**  
   [main.tex:405](/home/anze/Dropbox/work/LyaLenser/report/main.tex:405); [pair_template_checks.wls:35](/home/anze/Dropbox/work/LyaLenser/mathematica/pair_template_checks.wls:35), check 4.

   Although \(\xi\nabla\xi\) is odd, \(\Delta\alpha_{ab}\) also changes sign when the pair is reversed. Thus
   \[
   \xi\,\nabla\xi\cdot\Delta\alpha
   \]
   is **even**. For constant convergence, it is proportional to \(-\kappa r\,\xi\xi'\), which survives an isotropic orientation average. Check 4 verifies a product-rule identity, not cancellation of \(b\).

   Computing the stated sum exactly removes the modeled mean field, provided \(\xi_{ij}\) is the correct observed baseline covariance. It does not make covariance-model errors disappear.

   **Fix:** remove the oddness argument; retain template-specific conditional subtraction and validate it against unlensed mocks with the actual geometry, projection, and selection.

4. **Blocker — Coordinate relabelling does not establish the claimed model-independent injection calibration.**  
   [main.tex:429](/home/anze/Dropbox/work/LyaLenser/report/main.tex:429); [PLAN.md:20](/home/anze/Dropbox/work/LyaLenser/PLAN.md:20).

   The **sign is correct**: carrying samples from source coordinates to \(\theta'=\theta-A\alpha_{\rm inj}\) gives source positions \(\theta'+A\alpha_{\rm inj}\) at first order. This is a valid restricted mock construction. It changes the observed sampling geometry, however, and is not automatically equivalent to lensing a field sampled at fixed observed positions.

   More seriously, recomputing a model-dependent statistic can measure the derivative of the model rather than the unknown physical response. Write the unnormalized statistic as
   \[
   Q(g,d)=\sum_e K_e(g)[d_id_j-C^{\rm mod}_e(g)],
   \]
   where \(g\) denotes geometry. Coordinate shifts give
   \[
   \frac{d\,\mathbb E Q}{dA}
   =\sum_e\dot K_e(C^{\rm true}_e-C^{\rm mod}_e)
     -K_e\dot C^{\rm mod}_e.
   \]
   The desired fixed-geometry lensing response is instead \(\sum_eK_eR^{\rm true}_e\). Equality needs additional assumptions.

   A concrete counterexample: take one retained pair, constant relative injected displacement, and a locally linear model \(\xi_{\rm mod}(r)=c+gr\). The shifted score has slope \(g^2\), independently of the data. Real lensing has score response \(g\,\xi'_{\rm true}\). For \(g=2\), \(\xi'_{\rm true}=5\), the injection returns 4 while the physical response is 10.

   **Fix:** withdraw “any approximation … affects only optimality.” Establish calibration with fixed-observed-geometry simulations or a validated projected response operator. Use paired \(+A/-A\) injections and subtract the zero-injection result. Specify whether pair membership, cuts, masks, template evaluation positions, and \(b\) are recomputed. Re-associating pixel pairs is necessary for the relabelled catalogue, but also introduces selection/boundary responses that must be separated from physical lensing.

5. **Significant — Check 5 does not validate the continuum response used by the estimator.**  
   [main.tex:430](/home/anze/Dropbox/work/LyaLenser/report/main.tex:430); [pair_template_checks.wls:38](/home/anze/Dropbox/work/LyaLenser/mathematica/pair_template_checks.wls:38).

   For fixed, sightline-specific projectors and a common source plane, the correct identity is
   \[
   R^{\rm obs}_{ab}
   =\left[P_a(\nabla_\theta C_{ab})P_b^T\right]\cdot\Delta\alpha_{ab}
   =\nabla_\theta(P_aC_{ab}P_b^T)\cdot\Delta\alpha_{ab}.
   \]
   This remains exact when \(P_a\) and \(P_b\) have different lengths and weights. Their variation does not itself make the pair-specific identity approximate.

   What is unverified is replacing that pair-specific derivative by the derivative of one geometry-averaged, distorted \(\xi(r_\perp,r_\parallel)\). Check 5 never performs that averaging. It also assumes fixed projectors and weights; dependencies on fitted data require separate treatment.

   With redshift-dependent deflections, the displacement cannot generally be pulled through the radial projection or compressed into one \(\Delta\alpha_{ab}\). The underlying discrete-estimator paper explicitly identifies the common-deflection approximation. [Metcalf, Tessore & Croft](https://arxiv.org/pdf/2005.04109)

   **Fix:** specify the projected derivative actually accumulated in the pair pass, and quantify the source-plane approximation per sub-slab. Review-1 item 14 remains unresolved.

6. **Significant — The pair compression is exact for the specified statistic, but “no per-sightline compression exists” is demonstrably wrong.**  
   [main.tex:408](/home/anze/Dropbox/work/LyaLenser/report/main.tex:408), Eq. `Vab`; [main.tex:421](/home/anze/Dropbox/work/LyaLenser/report/main.tex:421); check 3.

   Let \(H_{ab}=V_{ab}-B_{ab}\), using stored orientation \(a<b\), and define
   \[
   U_a=\sum_{b>a}H_{ab}-\sum_{b<a}H_{ba}.
   \]
   Then exactly,
   \[
   Q_T=\sum_{a<b}H_{ab}\cdot(\alpha_a-\alpha_b)
       =\sum_a U_a\cdot\alpha_a.
   \]
   Thus **one vector per sightline suffices for the numerator**. Check 3’s local isotropic cancellation does not invalidate this algebra.

   The normalization still contains pair couplings: \(F=\alpha^TL_M\alpha\), with \(L_M\) the block graph Laplacian built from \(M_{ab}\). A vector per sightline cannot encode this arbitrary quadratic form.

   Also, \(V,M,B\) preserve only the selected statistic with fixed kernels, weights, cuts, and source-plane treatment. They do not preserve arbitrary short-mode bins or the information needed for the promised \(\mathcal D/T\) versus \(\mathbf k\) consistency test.

   **Fix:** distinguish numerator compression from normalization storage. Retain additional bins/basis responses for the advertised diagnostics.

7. **Significant — Diagonal weights do not imply independent pair products or \(\mathrm{Var}(\hat A)=1/F\).**  
   [main.tex:405](/home/anze/Dropbox/work/LyaLenser/report/main.tex:405), [main.tex:434](/home/anze/Dropbox/work/LyaLenser/report/main.tex:434).

   For Gaussian pixels,
   \[
   \mathrm{Cov}(d_id_j,d_kd_l)=C_{ik}C_{jl}+C_{il}C_{jk}.
   \]
   Pairs \((a,b)\) and \((a,c)\) therefore correlate through terms such as \(C_{aa}C_{bc}\). Different pixel pairs within the same sightline pair also contain \(C_{aa}(p,p')C_{bb}(q,q')\), including strong radial correlations.

   Excluding \(a=b\) correctly removes direct same-sightline products from this estimator; their lensing response vanishes for a common transverse remapping. It **does not remove same-sightline covariance contractions** from its variance.

   **Fix:** specify the sandwich covariance of the quadratic scores, or a validated joint-mock estimate retaining shared sightlines and connected forest terms. The claim that the remaining suboptimality is only between neighboring sightlines overlooks correlations along each sightline. The asserted “under 20 per cent” gain is unverified.

8. **Significant — The curl statistic is a useful diagnostic, not an automatically clean null on the implemented geometry.**  
   [main.tex:427](/home/anze/Dropbox/work/LyaLenser/report/main.tex:427); [pair_template_checks.wls:45](/home/anze/Dropbox/work/LyaLenser/mathematica/pair_template_checks.wls:45), check 6.

   Check 6 returns
   \[
   \pi(g_1g_{2t}-g_{1t}g_2),
   \]
   after an isotropic orientation integral. It does not prove zero response for a fixed template and actual pair catalogue. For one horizontal pair and the symmetric lens Jacobian
   \[
   \Psi=\begin{pmatrix}1&1\\1&0\end{pmatrix},
   \]
   the same-field gradient/curl response product is \(-1\), not zero.

   A zero ensemble mean requires the appropriate symmetry of the fields **and the effective measurement operators**. Independent geometry with statistically isotropic fields can preserve the null; masking/filtering and conditioning on the actual template require explicit checks.

   **Fix:** calculate or inject the gradient-to-curl response through the footprint, preferably using a joint response matrix. A passing curl null cannot exclude parity-even density/tidal contamination, magnification, shared sampling, or higher-order lensing. Review-1 item 15 is only partly resolved.

9. **Significant — The CMB simulation ensemble remains ambiguous and is not a physical deprojection null.**  
   [main.tex:428](/home/anze/Dropbox/work/LyaLenser/report/main.tex:428); [PLAN.md:25](/home/anze/Dropbox/work/LyaLenser/PLAN.md:25).

   If simulations use \(X_r=\kappa_{{\rm CMB},r}-\hat\kappa_s^{\rm data}\), the mean **numerator** contains the fixed subtraction \(-Q[\hat\kappa_s^{\rm data}]\); it is not centered at zero. If they use only independent CMB simulations, they test a different statistic with different template power and no quasar subtraction.

   Every realization needs its own \(b[X_r]\) and normalization/response. Since an independent zero-mean template also averages an omitted mean field to zero, this ensemble can pass despite a wrong mean field for the real template.

   Fixed-forest simulations retain the realized forest structure but do not independently sample its full non-Gaussian covariance or its joint correlations with quasars and CMB lensing.

   **Fix:** define separately the conditional random-template diagnostic and the joint physical covariance ensemble. Moving joint mocks into phase 1 is a real improvement, but 10–20 mocks give approximately 32–47% fractional uncertainty on a Gaussian scalar variance and cannot establish a stable multidimensional covariance without further structure.

10. **Significant — Kernel matching resolves the arbitrary-radial-weight problem only within its stated transfer model.**  
    [main.tex:319](/home/anze/Dropbox/work/LyaLenser/report/main.tex:319), Eqs. `kernelresid`, `matched`; [main.tex:334](/home/anze/Dropbox/work/LyaLenser/report/main.tex:334); check 8.

    The pointwise cancellation is correct in Eq. `kernelresid`, and check 8 correctly verifies it for arbitrary \(A_\alpha(\chi)\). Review-1 item 4 is therefore resolved **for the local Limber density model**.

    Beyond Limber, even a perfect subtraction of the slab density leaves
    \[
    \int_{\rm slab}d\chi\,A_\alpha(\chi)
    \int_{\rm outside}d\chi'\,W_{\rm CMB}(\chi')
    C_L^{\delta(\chi)\delta(\chi')},
    \]
    because densities inside and outside the slab remain correlated. Exact matching removes the identical slab contribution; it does not imply that the remaining map is orthogonal to every slab response.

    Quasar RSD adds \(f\mu^2\delta/b_q\) to \(\delta_q/b_q\), and tidal responses have their own angular/radial transfer. The statement \(k_\parallel\lesssim0.01\) does not establish small \(\mu^2\): at \(L=40,\chi=3956\,h^{-1}\mathrm{Mpc}\), \(k_\perp\simeq0.010\,h\,\mathrm{Mpc}^{-1}\). Narrower slabs also broaden the radial Fourier window. Narrow-bin and small-overlap projections particularly require checking Limber accuracy. [LoVerde & Afshordi](https://arxiv.org/html/0809.5112v1)

    **Fix:** specify non-Limber density/tidal/RSD response tests and a residual-bias tolerance over the actual bands and sub-slabs.

11. **Significant — The quasar weight is correct, but the map normalization and general shot noise are underspecified.**  
    [main.tex:327](/home/anze/Dropbox/work/LyaLenser/report/main.tex:327), Eq. `matched`; line 330.

    The weight
    \[
    u(\chi)=\frac{W_{\rm CMB}(\chi)}
    {b_q(\chi)\bar n_{3D}(\chi)\chi^2}
    \]
    correctly converts a point catalogue into the desired angular field. The sum must mean \(\sum_i u_i\delta_D^{(2)}(\theta-\theta_i)\), or pixel sums divided by pixel solid angle. “Minus mean” must include the actual angular/radial selection through randoms.

    The general Poisson angular noise is
    \[
    N_L^{ss}=\int_{\rm slab}d\chi\,
    \frac{W_{\rm CMB}^2(\chi)}
    {b_q^2(\chi)\bar n_{3D}(\chi)\chi^2}.
    \]
    The quoted \(D^2\langle W^2\rangle/(b_q^2\bar n_q)\) follows for constant bias and a flat **angular count distribution per unit \(\chi\)**, \(d\bar N/(d\chi\,d\Omega)=\bar n_q/D\). It is not the general result for evolving quasar counts.

    **Fix:** specify pixel normalization, completeness/random weights, random-catalogue noise, and use the measured radial selection in the noise integral. Geometry and \(b_q\) are not the only implementation inputs.

12. **Significant — Magnification leaves intrinsic-response leakage, not just a modified lensing signal.**  
    [main.tex:336](/home/anze/Dropbox/work/LyaLenser/report/main.tex:336); [main.tex:426](/home/anze/Dropbox/work/LyaLenser/report/main.tex:426).

    With redshift-dependent coefficients, define
    \[
    M(\theta)=\int_{\rm slab}d\chi\,
    W_{\rm CMB}(\chi)\frac{m(\chi)}{b_q(\chi)}
    \kappa_q(\theta,\chi).
    \]
    Then \(\hat\kappa_s=\kappa_{\rm slab}+M+\epsilon_s+\cdots\). After matching the density term, the intrinsic forest response still contributes \(-C^{p_{\rm int}M}\). Including \(-C^{\kappa_{\rm Ly\alpha}M}\) in the signal template does not remove this separate contaminant.

    Likewise, the \(\hat\kappa_s\)-only control measures intrinsic response **plus real forest-lensing correlation** with the slab and magnification components. It is not a pure response measurement.

    **Fix:** model or remove magnification using the actual quasar redshift distribution and \(m(z)\), propagating its uncertainty and intrinsic-response leakage. The scalar-\(\beta\) discussion does not specify the correction for the new matched estimator.

13. **Significant — The \(1-r_{q\delta}^2\) stochasticity claim was carried over from a different estimator.**  
    [main.tex:335](/home/anze/Dropbox/work/LyaLenser/report/main.tex:335); its cited [three_tracer_checks.wls:12](/home/anze/Dropbox/work/LyaLenser/mathematica/three_tracer_checks.wls:12).

    That check derives the residual for \(\beta=C^{q\kappa}/C^{qq}\). The redesign uses known bias weighting and \(\beta=1\). If \(q=b\delta+\epsilon\), with \(\epsilon\) uncorrelated with the intrinsic forest response, additional stochasticity adds variance without producing a universal \(1-r^2\) residual. If it correlates with that response, the bias depends on the relevant cross-spectrum.

    **Fix:** replace the inherited formula with explicit assumptions about \(C^{p_{\rm int}\epsilon}\), scale-dependent bias, and how \(b_q\) is inferred.

14. **Significant — Shared-quasar sampling remains an estimator bias requirement, and the implementation phases contradict the adopted remedy.**  
    [main.tex:337](/home/anze/Dropbox/work/LyaLenser/report/main.tex:337); [PLAN.md:23](/home/anze/Dropbox/work/LyaLenser/PLAN.md:23), [PLAN.md:70](/home/anze/Dropbox/work/LyaLenser/PLAN.md:70).

    Exact subtraction of \(b\) handles pure geometry multiplying the correct baseline covariance. It does not prove
    \[
    \mathbb E[d_id_j\mid\text{selected quasars, weights, template}]=C^0_{ij}.
    \]
    Shared objects can correlate template shot noise, forest coverage, source properties, continuum errors, and flux-dependent weights. Disjoint object selections remove direct shared-object terms but do not remove shared-LSS or magnification correlations.

    A sightline-density–\(\kappa\) spectrum diagnoses one part of this; neither a passing curl test nor independent CMB templates bounds the resulting parity-even bias.

    **Fix:** make joint selection/response mocks and the intended source/template selections part of the estimator specification. Section 0 says the broad slab is forecasts-only, while phase 1 explicitly measures it and postpones tomography and magnification to phase 2. Resolve this before assigning a clean deprojected amplitude.

15. **Significant — First-order exactness and the proposed nulls do not bound higher-order cross-response bias.**  
    [pair_template_checks.wls:9](/home/anze/Dropbox/work/LyaLenser/mathematica/pair_template_checks.wls:9), check 1; [PLAN.md:79](/home/anze/Dropbox/work/LyaLenser/PLAN.md:79).

    The second-order covariance term is
    \[
    \tfrac12\Delta\alpha_\mu\Delta\alpha_\nu
    \partial_\mu\partial_\nu\xi.
    \]
    Its cross with the external tracer requires \(\langle\alpha\alpha X\rangle\): zero for jointly zero-mean Gaussian lensing/tracer fields under the independent-forest assumptions, but potentially nonzero with a lensing bispectrum or additional correlations. Gaussian cross-response corrections arise at third remapping order.

    The script still calls the second-order term “needed for the N1-type bias,” reviving the interpretation rejected in review 1. Independent-template and curl nulls can miss these physical, parity-even corrections.

    **Fix:** correct the check’s description and define a quantitative higher-order response test. Small independent injections alone do not calibrate the correlated physical correction.

16. **Minor — \(O(10^7)\) is a plausible contraction count, not the complete cost per requested operation.**  
    [main.tex:421](/home/anze/Dropbox/work/LyaLenser/report/main.tex:421), [main.tex:434](/home/anze/Dropbox/work/LyaLenser/report/main.tex:434); [PLAN.md:49](/home/anze/Dropbox/work/LyaLenser/PLAN.md:49).

    A fixed template with prepared sightline deflections costs \(O(N_{\rm pair})\) for \(F\); the numerator can use the sightline compression above. Template construction, harmonic filtering, interpolation, and band generation are additional costs. \(B\) bands require \(O(BN_{\rm pair})\) diagonal contractions; a direct full band-response matrix requires \(O(B^2N_{\rm pair})\).

    Injections explicitly require a new pixel pass, contradicting their inclusion in the “\(10^7\) operations” sentence. Other changes can also force new passes: changing \(\nabla\xi\), continuum treatment, pixel masks/weights, unstored radial or separation bins, or source-plane response bases. Template-only filtering changes do not.

    Seven float32 values for \(10^7\) pairs occupy 280 MB before pair indices and metadata; float64 doubles this. The claimed storage is therefore precision-dependent.

    **Fix:** state contraction costs separately from map construction and catalogue rebuilding; document retained bins, precision, and metadata. Review-1 item 18’s arithmetic is corrected, but runtime remains a benchmark requirement.

Relative to review 1: **13 is resolved within the explicit source-plane model; 4 is resolved within the Limber density model; 14 is not resolved; 15 is partly resolved; 17 improves the validation schedule but remains incomplete; 18 corrects the arithmetic but awaits measurement.**

**Verdict: the estimator is not ready to be implemented as specified.** The pair-response and compression machinery are sound foundations. First fix the physical amplitude definition and normalization, band response, projected covariance derivative, and injection-calibration argument. Then specify valid deprojection/null ensembles and how magnification, radial-transfer mismatch, and shared sampling are bounded. Otherwise the implementation could recover proportional-template injections and pass both proposed nulls while returning a biased physical cross-correlation amplitude.

