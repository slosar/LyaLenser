(* Check C: the "pair displacement" (real-space) picture reproduces the Fourier mode-coupling f.
   1D periodic toy on [0,2pi): delta(x) = Sum_q a_q e^{iqx}, <a_q a_q'^*> = P_q delta_qq'.
   Lensed field: dtilde(x) = delta(x + phi'(x)), phi(x) = phi0 Cos[L0 x]  (so phi_{L0} = phi0/2).
   Pair covariance: Ctilde(xa,xb) = xi(xa - xb + phi'(xa) - phi'(xb)) ~ C + xi'(xa-xb) (phi'(xa)-phi'(xb)).
   Claim: <atilde_l atilde_l'^*> = phi_{L0} L0 (l P_l - l' P_l')  for l - l' = L0,
   i.e. with l'' = -l':  phi_{L0} [ L0 l P_l + L0 l'' P_l'' ]  = f(l, l'') phi(l+l'').  *)
qmax = 4;
xi[s_] := Sum[P[q] Exp[I q s], {q, -qmax, qmax}];
phi[x_, L0_] := phi0 Cos[L0 x];
dC[xa_, xb_, L0_] := xi'[xa - xb] (D[phi[xa, L0], xa] - D[phi[xb, L0], xb]);
coupling[l_Integer, lp_Integer, L0_Integer] := Simplify[(1/(2 Pi))^2 Integrate[Expand[dC[xa, xb, L0] Exp[-I l xa + I lp xb]], {xa, 0, 2 Pi}, {xb, 0, 2 Pi}]];
res = Table[{L0v, lv, lv - L0v,
    Simplify[coupling[lv, lv - L0v, L0v] - (phi0/2) L0v (lv P[lv] - (lv - L0v) P[lv - L0v])]},
   {L0v, 1, 2}, {lv, -1, 2}];
Print["Check C residuals {L0, l, l', residual} (all residuals should be 0): ", Flatten[res, 1]];
(* also check that a pair with l - l' != L0 has no coupling *)
Print["Check C2 (no coupling when l-l' != L0, should be 0): ", coupling[2, 2, 1], ", ", coupling[3, 1, 1]];
