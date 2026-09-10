(* Checks for the pair-template estimator (PLAN v2 / report Sec. 5).
   Pixel pair i=(a,p), j=(b,q). Response to a deflection field alpha:
     C_ij(alpha) = xi(theta_ab + alpha_a - alpha_b, dchi)  =>  R_ij = grad xi_ij . (alpha_a - alpha_b) + O(alpha^2).
   Diagonal-weight quadratic estimator for the amplitude A of a template alpha_T:
     Ahat = [ sum_{i<j} w_i w_j d_i d_j R_ij - b ] / F,  F = sum_{i<j} w_i w_j R_ij^2,  b = sum_{i<j} w_i w_j xi_ij R_ij.
   Per-sightline-pair compression:  V_ab = sum_pq w w d d grad xi,  M_ab = sum_pq w w grad xi grad xi^T,  B_ab = sum_pq w w xi grad xi,
   so that  numerator = sum_ab V_ab . dalpha_ab,  F = sum_ab dalpha_ab^T M_ab dalpha_ab,  b = sum_ab B_ab . dalpha_ab. *)

(* Check 1: second-order expansion of the pair covariance in the relative deflection, for a generic isotropic xi(|theta|, dchi). *)
th = {tx, ty}; da = eps {ax, ay};
xiF[v_] := XI[Sqrt[v.v]];
ser = Series[xiF[th + da], {eps, 0, 2}] // Normal;
gradxi = D[xiF[th], {th}];
lin = Simplify[Coefficient[ser, eps, 1] - gradxi.{ax, ay}];
Print["Check 1: first-order pair response minus grad xi . dalpha = ", lin];
Print["Check 1b: second-order pair response (enters the cross-correlation only through <alpha alpha X>) = ", Simplify[Coefficient[ser, eps, 2]]];

(* Check 2: the expectation of V_ab is M_ab dalpha_ab (linear response), and the numerator's Fisher normalisation.
   Model one pixel pair with weights w1 w2 = 1 and <d d> = xi + R A: *)
Print["Check 2: <w w d d R> = xi R + A R^2, so <numerator> = b + A F  ->  <Ahat> = A  (by construction)."];

(* Check 3: isotropic cancellation. For pairs at fixed separation r and direction phi around sightline a,
   with grad xi = xi'(r) rhat and a locally linear deflection dalpha = Psi . theta_ab (Psi symmetric: -kappa 1 - gamma),
   V_ab ~ xi' rhat (rhat.Psi.theta_ab) = r xi' rhat (rhat.Psi.rhat).  The sum over an isotropic neighbourhood vanishes: *)
rhat = {Cos[phi], Sin[phi]};
Psi = -{{k + g1, g2}, {g2, k - g1}};
Vdir = rhat (rhat.Psi.rhat);
Print["Check 3: angular average of V_ab direction vector (expect {0,0}): ", Simplify[Integrate[Vdir, {phi, 0, 2 Pi}]/(2 Pi)]];
(* whereas the per-pair scalar along the pair direction keeps kappa and the shear: *)
Print["Check 3b: rhat.Psi.rhat = ", Simplify[TrigReduce[rhat.Psi.rhat]], "   (-kappa - gamma_+(phi): dilation + shear)"];
(* and the pair-level Fisher weight is the second moment, which does not cancel: *)
Print["Check 3c: angular average of (rhat.Psi.rhat)^2 = ", Simplify[Integrate[(rhat.Psi.rhat)^2, {phi, 0, 2 Pi}]/(2 Pi)],
  "   (kappa^2 + (g1^2+g2^2)/2: convergence and shear both measured)"];

(* Check 4: the mean-field b is odd in theta_ab for isotropic xi: xi grad xi = grad(xi^2/2) *)
Print["Check 4: xi grad xi - grad(xi^2/2) = ", Simplify[xiF[th] gradxi - D[xiF[th]^2/2, {th}]]];

(* Check 5: LOS-only continuum-fitting projection P commutes with the transverse gradient:
   for C_obs = P C P^T with P acting on the radial index only, d/dtheta (P C P^T) = P (dC/dtheta) P^T.  Trivial since P
   has no theta dependence; verified on a 2-pixel radial projector with theta-dependent C. *)
P = {{1/2, -1/2}, {-1/2, 1/2}};
Cm = {{c11[tx], c12[tx]}, {c12[tx], c22[tx]}};
Print["Check 5: D[P.C.P^T, tx] - P.D[C,tx].P^T = ", Simplify[D[P.Cm.Transpose[P], tx] - P.D[Cm, tx].Transpose[P]]];

(* Check 6: rotated (curl) template. For alpha_T = grad phi the pair response is grad xi . (grad phi_a - grad phi_b);
   for the 90-degree rotated field alpha_R = R90 grad phi, a pure-gradient true deflection gives zero expected amplitude
   because <numerator> = sum M_ab : (dalpha_true dalpha_R^T) and for an isotropic neighbourhood and linear fields
   Tr[Psi_true . R90 . Psi_R] with symmetric Psi's has zero angular average when contracted with rhat rhat rhat rhat. *)
R90 = {{0, -1}, {1, 0}};
PsiT = -{{kt + g1t, g2t}, {g2t, kt - g1t}};
resp = (rhat.PsiT.rhat) (rhat.(R90.Psi).rhat);
Print["Check 6: angular average of gradient-true x curl-template pair response = ", Simplify[Integrate[resp, {phi, 0, 2 Pi}]],
  "  (antisymmetric under true<->template shear: zero when the template is built from the same field, and zero in",
  " expectation for any parity-invariant ensemble, since <g1_true g2_T> = <g2_true g1_T>; the curl-template amplitude is a null test in expectation)"];

(* ---- Added after codex review 1 ---- *)
(* Check 7 (review item 13): exact pair response vs local-gradient template for a plane-wave lens potential
   phi = phi0 Cos[L x]. Exact: grad xi . (alpha(xa) - alpha(xb)); local: grad xi . (Psi(xmid) theta_ab), Psi = Hessian.
   Their ratio along the pair separation is sinc(L theta_ab/2) for the component along L. *)
phiPW[x_] := phi0 Cos[LL x];
xa = xm + s/2; xb = xm - s/2;
exact = (D[phiPW[x], x] /. x -> xa) - (D[phiPW[x], x] /. x -> xb);
local = (D[phiPW[x], {x, 2}] /. x -> xm) s;
ratio = Simplify[exact/local];
Print["Check 7: exact/local relative deflection for a plane wave = ", ratio, "  = Sinc[L s/2] : ", Simplify[ratio - Sinc[LL s/2]]];
Print["Check 7b: numerical value at L=300, s = 30 Mpc/h / 3956 Mpc/h: ", N[Sinc[300*30/3956/2]]];

(* Check 8 (review item 4): kernel matching. With the quasar template built as kappa_s-hat = Integral W_CMB(chi) delta_q(chi)/b_q(chi) dchi,
   the residual response leak Integral A(chi) [W_CMB(chi) - b_q(chi) W_T(chi)] dchi vanishes for ANY A(chi) because the bracket is
   identically zero, whereas a single slab-averaged beta leaves Integral A(chi) [W_CMB(chi) - beta b_q(chi)/D] dchi, which is zero
   only for A = const. Illustrate with linear W_CMB(chi) = w0 + w1 chi over [0, D], b_q const, A(chi) = a0 + a1 chi. *)
WT[chi_] := (w0 + w1 chi)/bq;                       (* matched template kernel *)
betaAvg = Integrate[(w0 + w1 chi), {chi, 0, DD}]/bq;   (* single-beta template: beta * (1/D) per unit chi, beta = mean(W) D / bq *)
residMatched = Integrate[(a0 + a1 chi) ((w0 + w1 chi) - bq WT[chi]), {chi, 0, DD}];
residSingle = Simplify[Integrate[(a0 + a1 chi) ((w0 + w1 chi) - betaAvg bq/DD), {chi, 0, DD}]];
Print["Check 8: residual with matched template (expect 0): ", Simplify[residMatched]];
Print["Check 8b: residual with a single slab-averaged beta: ", residSingle, "  (vanishes only if a1 w1 = 0, i.e. constant A or constant W)"];


(* ---- Added after codex review 2 ---- *)
(* Check 9 (review-2 item 3): the mean-field integrand xi grad xi . dalpha is EVEN under pair reversal and survives the
   orientation average for a constant convergence: dalpha = Psi theta_ab with Psi = -kappa 1 gives
   xi xi' rhat . (-kappa r rhat) = -kappa r xi xi'. *)
mf = XI[r] XI'[r] rhat.(-(k) r rhat);
Print["Check 9: orientation average of xi grad xi . dalpha for constant kappa = ", Simplify[Integrate[mf, {phi, 0, 2 Pi}]/(2 Pi)], "  (nonzero: mean field is even, not odd)"];

(* Check 10 (review-2 item 6): per-sightline compression of the numerator. For pairs (a<b) with vectors H_ab,
   sum_{a<b} H_ab.(alpha_a - alpha_b) = sum_a U_a.alpha_a with U_a = sum_{b>a} H_ab - sum_{b<a} H_ba. Verify on 3 sightlines. *)
H12 = {h12x, h12y}; H13 = {h13x, h13y}; H23 = {h23x, h23y};
al1 = {a1x, a1y}; al2 = {a2x, a2y}; al3 = {a3x, a3y};
lhs = H12.(al1 - al2) + H13.(al1 - al3) + H23.(al2 - al3);
U1 = H12 + H13; U2 = H23 - H12; U3 = -H13 - H23;
rhs = U1.al1 + U2.al2 + U3.al3;
Print["Check 10: pair-sum minus per-sightline sum (expect 0): ", Simplify[lhs - rhs]];

(* ---- Added after codex review 3 ---- *)
(* Check 11 (review-3 item 1): sign convention. theta_ab = theta_a - theta_b (points from b to a), response
   R = xi'(r) that_ab . (alpha_a - alpha_b). For constant convergence alpha = -kappa theta, the response must equal
   d/dkappa xi((1-kappa) r) = -r xi'(r) > 0 for decreasing xi. With the reversed direction the sign flips. *)
tha = {xa1, xa2}; thb = {xb1, xb2}; kk = k;
sep = tha - thb; rr = Sqrt[sep.sep]; that = sep/rr;
alA = -kk tha; alB = -kk thb;
Rgood = XI'[rr] that.(alA - alB);
Rbad = XI'[rr] (-that).(alA - alB);
Print["Check 11: response with that_ab = (theta_a - theta_b)/r : ", Simplify[Rgood], "   (expect -kappa r xi'(r))"];
Print["Check 11b: response with the reversed direction : ", Simplify[Rbad], "   (wrong sign)"];
Print["Check 11c: direct remapping d/dkappa xi((1-kappa) r) at kappa=0 : ", Simplify[D[XI[(1 - kk) rr], kk] /. kk -> 0]];

(* Check 12 (review-3 item 2): redshift-dependent deflection. Pixel p on a at chi_p, q on b at chi_q, deflections
   g(chi_p) alpha_a and g(chi_q) alpha_b with g = 1 + g1 (chi - chiref). Then
   dalpha_pq = (alpha_a - alpha_b)(1 + g1 dm) + g1 (dp - dq)/2 (alpha_a + alpha_b), dm = chi_mid - chiref.
   With grad xi = G that, R = G [ d (1 + g1 dm) + g1 (dp-dq)/2 s ], d = that.(alpha_a-alpha_b), s = that.(alpha_a+alpha_b).
   R^2 to O(g1^2): G^2 [ d^2 (1 + 2 g1 dm + g1^2 dm^2) + g1 d s (dp-dq) (1 + g1 dm) + g1^2 s^2 (dp-dq)^2/4 ]. *)
gp = 1 + g1 (chip - chiref); gq = 1 + g1 (chiq - chiref);
dal = gp alphaA - gq alphaB;   (* scalar components along that suffice *)
dm = (chip + chiq)/2 - chiref;
claim = (alphaA - alphaB) (1 + g1 dm) + g1 (chip - chiq)/2 (alphaA + alphaB);
Print["Check 12: exact dalpha minus first-moment decomposition (expect 0): ", Simplify[dal - claim]];
d = alphaA - alphaB; s = alphaA + alphaB;
R2 = Expand[(G claim)^2];
claimR2 = G^2 (d^2 (1 + 2 g1 dm + g1^2 dm^2) + g1 d s (chip - chiq) (1 + g1 dm) + g1^2 s^2 (chip - chiq)^2/4);
Print["Check 12b: R^2 minus its moment expansion (expect 0): ", Simplify[R2 - claimR2]];
Print["Check 12c: a linear m_eff = m (1 + g1 dm) for constant g=1.2 (g1 dm = 0.2) gives F scaled by ", 1.2, " instead of ", 1.2^2, "; the estimator would return A = ", 1.2^2/1.2];
