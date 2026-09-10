(* Checks for the three-tracer (quasar q, forest power modulation p, CMB kappa) separation.
   Squeezed-limit model for the cross-spectra of the local-forest-power field p_alpha with long-mode fields:
     C^{p k} = R Cdk + T Clc        (R = unknown forest response, T = known lensing template,
     C^{p q} = R Cdq + T Clq         Cdk = <delta_2.4 kappa_cmb>, Clc = <kappa_lya kappa_cmb>, ...)
   Deprojected statistic: Dst = C^{pk} - beta C^{pq} with beta = Cdk/Cdq.                                  *)

(* Check 1: exact cancellation of the unknown response R *)
Cpk = R Cdk + T Clc; Cpq = R Cdq + T Clq; beta = Cdk/Cdq;
Dst = Simplify[Cpk - beta Cpq];
Print["Check 1: deprojected statistic = ", Dst, " ; d/dR = ", Simplify[D[Dst, R]]];

(* Check 2: stochastic tracer q = b delta + eps, <eps delta>=<eps kappa>=0, Cee = <eps eps> (non-Poisson part);
   beta estimated as C^{qk}/C^{qq,no shot} *)
Cqk = b Cdk; Cqq = b^2 Cdd + Cee; Cdq = b Cdd;
betaHat = Cqk/Cqq;
resid = Simplify[R (Cdk - betaHat Cdq)];
r2 = b^2 Cdd/(b^2 Cdd + Cee);
Print["Check 2: residual response leak / (R Cdk) = ", Simplify[resid/(R Cdk)], "  (expect 1 - r^2 = ", Simplify[1 - r2], ")"];

(* Check 3: quasar magnification, q = b delta + m kappa_q (m = 5s-2), Cee = 0.
   Cdkq = <delta kappa_q>, Ckqk = <kappa_q kappa_cmb>, Ckqkq = <kappa_q kappa_q>. *)
Clear[Cqk, Cqq, Cdq];
Cqk = b Cdk + m Ckqk; Cqq = b^2 Cdd + 2 b m Cdkq + m^2 Ckqkq; Cdq = b Cdd + m Cdkq;
betaHat = Cqk/Cqq;
leak = Series[R (Cdk - betaHat Cdq), {m, 0, 1}] // Normal // Simplify;
Print["Check 3: response leak from ignoring magnification, first order in m: ", leak];
Print["Check 3b: fractional error in beta at first order in m: ", Simplify[Normal[Series[betaHat/(Cdk/(b Cdd)) - 1, {m, 0, 1}]]]];

(* Check 4: dilation response integrates to a boundary term, so the per-skewer (transverse-integrated) power is
   lensing invariant: Integral_0^lmax (2P + l P') l dl = lmax^2 P(lmax). *)
Print["Check 4: l (2P + l P') - d(l^2 P)/dl = ", Simplify[l (2 P[l] + l P'[l]) - D[l^2 P[l], l]],
  "  => int_0^lmax (2P+lP') l dl = lmax^2 P(lmax), a pure boundary term"];

(* Check 5: for k_par = 0 long modes the tidal tensor is fixed by delta_L: K_ij = (Lhat_i Lhat_j - delta_ij/3) delta_L,
   so k_hat.K.k_hat for a short mode with angle phi to L in the transverse plane and mu along the LOS
   is a linear combination of {1, cos 2 phi} with mu-dependent coefficients only. *)
Lhat = {1, 0, 0}; khat = {Sqrt[1 - mu^2] Cos[phi], Sqrt[1 - mu^2] Sin[phi], mu};
K = Outer[Times, Lhat, Lhat] - IdentityMatrix[3]/3;
proj = Simplify[khat.K.khat];
Print["Check 5: khat.K.khat = ", TrigReduce[proj], "  (only a constant and a cos 2phi harmonic)"];
