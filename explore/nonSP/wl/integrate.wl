(* usage: wolframscript -file integrate.wl base
   F(g, r) = P/Q^10 on Delta^2 x [0,1]^3, Q = gA rhoA + gB rhoB + gC rhoC (rho = 1 - r), gC = 1 - gA - gB.
   Homogenize P to degree Dg = max(deg_g P, 7) with G = gA + gB + gC (= 1 on the simplex);
   b = Dg - 7.  b = 0: Feynman formula  int_{Delta^2} g^beta/(c.g)^10 = beta!/9! prod c^(-beta-1), c = rho.
   b > 0:  1/(Q^10 G^b) = Gamma(10+b)/(Gamma(10)Gamma(b)) int_0^1 xi^9 (1-xi)^(b-1) (xi Q + (1-xi) G)^(-10-b),
           c_v = 1 - xi r_v;  then r-integrals factorize:  phi[gam, m] = int_0^1 r^gam (1 - xi r)^(-m) dr. *)
base = $ScriptCommandLine[[2]];
apow = If[Length[$ScriptCommandLine] >= 3, ToExpression[$ScriptCommandLine[[3]]], 10];
P = Get[base <> "_P.m"]; Qp = Get[base <> "_Q.m"];
t0 = AbsoluteTime[];
PP = Expand[P /. {gA -> lam gA, gB -> lam gB}];
degG = Exponent[PP, lam];
a = apow; Dg = Max[degG, a - 3]; b = Dg - (a - 3);
Print["joint g-degree ", degG, "  Dg ", Dg, "  b ", b];
parts = CoefficientList[PP, lam];
Ph = Expand[Sum[parts[[d + 1]] (gA + gB + gC)^(Dg - d), {d, 0, Length[parts] - 1}]];
check = Expand[(Ph /. gC -> 1 - gA - gB) - P];
Print["homogenization check (0 expected): ", check];
rules = CoefficientRules[Ph, {gA, gB, gC}];
Print["g-monomials ", Length[rules]];
If[b == 0,
  G1 = Sum[With[{be = rl[[1]], c = rl[[2]]},
        Factorial[be[[1]]] Factorial[be[[2]]] Factorial[be[[3]]]/Factorial[a - 1] *
        (c /. {rA -> 1 - pA, rB -> 1 - pB, rC -> 1 - pC}) pA^(-be[[1]] - 1) pB^(-be[[2]] - 1) pC^(-be[[3]] - 1)],
      {rl, rules}];
  G1 = Expand[G1];
  minexp = Min[Flatten[Exponent[#, {pA, pB, pC}, Min] & /@ (List @@ G1)]];
  Print["min exponent of rho after g-integration: ", minexp];
  res = Integrate[G1, {pA, 0, 1}, {pB, 0, 1}, {pC, 0, 1}];
  Print["RESULT c = ", InputForm[res], "  N = ", N[res, 30]];
  ,
  (* b > 0 *)
  phi[gam_, m_] := phi[gam, m] = Simplify[Integrate[r^gam (1 - xi r)^(-m), {r, 0, 1}, Assumptions -> 0 < xi < 1]];
  terms = Flatten[Table[With[{be = rl[[1]], cr = CoefficientRules[rl[[2]], {rA, rB, rC}]},
        Table[{be, cc[[1]], cc[[2]]}, {cc, cr}]], {rl, rules}], 1];
  Print["(beta, gamma) terms ", Length[terms]];
  pref = 1/(Gamma[a] Gamma[b]);
  H = Sum[With[{be = tm[[1]], ga = tm[[2]], c = tm[[3]]},
        c Factorial[be[[1]]] Factorial[be[[2]]] Factorial[be[[3]]] *
        phi[ga[[1]], be[[1]] + 1] phi[ga[[2]], be[[2]] + 1] phi[ga[[3]], be[[3]] + 1]], {tm, terms}];
  Print["H built, time ", AbsoluteTime[] - t0];
  Hx = pref xi^(a - 1) (1 - xi)^(b - 1) H;
  Put[Hx, base <> "_Hxi.m"];
  hn = NIntegrate[Hx, {xi, 0, 1/2, 1}, WorkingPrecision -> 60, PrecisionGoal -> 40, MaxRecursion -> 30];
  Print["NUMERIC c = ", N[hn, 45]];
  hL = Collect[Expand[Hx /. Log[1 - xi] -> LL], LL];
  Put[hL, base <> "_HL.m"];
  res = Integrate[hL /. LL -> Log[1 - xi], {xi, 0, 1}];
  Print["RESULT c = ", InputForm[Simplify[res]], "  N = ", N[res, 30]];
];
Print["total time ", AbsoluteTime[] - t0];
