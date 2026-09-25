(* Checks for the first-order lensing response of a 2D power spectrum.
   Conventions: dtilde(x) = d(x + grad phi(x)),  kappa = -Laplacian[phi]/2,
   <dtilde(l) dtilde(l')> = (2pi)^2 delta(l+l') P(l) + f(l,l') phi(l+l'),
   f(l,l') = L.l P(l) + L.l' P(l'),  L = l + l'.
   Check A: squeezed expansion of f(l, L-l) in powers of L.
   Check B: exact power spectrum under constant kappa, gamma1, gamma2
            equals the squeezed-limit prediction 2 kappa P + (kappa + g1 cos2phi + g2 sin2phi) l P'. *)

(* ---------- Check A ---------- *)
l = {lx, ly}; Lv = eps {Lx, Ly};
lp = Lv - l;
norm[v_] := Sqrt[v.v];
f = Lv.l P[norm[l]] + Lv.lp P[norm[lp]];
ser = Series[f, {eps, 0, 2}] // Normal // Simplify;
(* expected: eps^2 [ L^2 P(l) + (L.l)^2 P'(l)/l ] *)
expected = eps^2 (Lv.Lv/eps^2 P[norm[l]] + (Lv.l)^2/eps^2 P'[norm[l]]/norm[l]);
Print["Check A (squeezed expansion residual, should be 0): ", Simplify[ser - expected]];

(* angle average over direction of l relative to L, for L along x *)
lA = lr {Cos[t], Sin[t]}; LA = {LL, 0};
fk = 2/LL^2 (LL^2 P[lr] + (LA.lA)^2 P'[lr]/lr);   (* f_kappa = 2 f / L^2 *)
avg = Integrate[fk, {t, 0, 2 Pi}]/(2 Pi) // Simplify;
Print["Check A2 angle-averaged kappa response f_kappa (expect 2P + l P'): ", avg];
Print["Check A2 cos(2t) part (expect l P'): ", Simplify[Integrate[fk Cos[2 t], {t, 0, 2 Pi}]/Pi]];
Print["Check A3: for P = l^n, angle-averaged response/P = ", Simplify[(avg/P[lr]) /. {P -> (#^n &)}]];

(* ---------- Check B ---------- *)
(* phi = -(k/2)(x^2+y^2) - (g1/2)(x^2-y^2) - g2 x y  ;  grad phi = -Psi.x with Psi = [[k+g1, g2],[g2, k-g1]] *)
phi = -(k/2) (x^2 + y^2) - (g1/2) (x^2 - y^2) - g2 x y;
Print["Check B0: kappa from -Laplacian/2 = ", Simplify[-(D[phi, x, x] + D[phi, y, y])/2]];
M = IdentityMatrix[2] + {{D[phi, x, x], D[phi, x, y]}, {D[phi, x, y], D[phi, y, y]}};  (* dtilde(x) = d(M.x) *)
(* If dtilde(x) = d(M x), then dtilde(l) = d(M^-T l)/|det M| and Ptilde(l) = P(|M^-T l|)/|det M| *)
lvec = lr {Cos[t], Sin[t]};
Mit = Transpose[Inverse[M]];
Pt = P[norm[Mit.lvec]]/Det[M];
dP = Series[Pt /. {k -> e k, g1 -> e g1, g2 -> e g2}, {e, 0, 1}] // Normal;
dP1 = Simplify[Coefficient[dP, e, 1], Assumptions -> lr > 0];
pred = 2 k P[lr] + (k + g1 Cos[2 t] + g2 Sin[2 t]) lr P'[lr];
Print["Check B (dilation/shear power response minus squeezed prediction, should be 0): ", Simplify[dP1 - pred, Assumptions -> lr > 0]];
