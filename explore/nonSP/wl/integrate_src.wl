(* usage: integrate_src.wl base : symbolic sources SA, SB present.  Loads base_HL.m (= Collect[Hx, LL], LL = Log[1-xi]),
   splits into monomials SA^j SB^k, integrates each over xi in (0,1) exactly; also exports R_k(xi) per monomial. *)
base = $ScriptCommandLine[[2]];
hL = Get[base <> "_HL.m"];
t0 = AbsoluteTime[];
kmax = Exponent[hL, LL];
mons = Flatten[Table[SA^j SB^k, {j, 0, 4}, {k, 0, 4 - j}]];
cf[expr_, j_, k_] := Coefficient[Coefficient[expr, SA, j], SB, k];
total = 0; s = OpenWrite[base <> "_HLexport_src.txt"];
Do[
  {j, k} = Exponent[m, {SA, SB}];
  pieces = Table[Together[cf[Coefficient[hL, LL, q], j, k]], {q, 0, kmax}];
  If[Union[pieces] === {0}, Continue[]];
  WriteString[s, "mon ", j, " ", k, "\n"];
  Do[Rq = pieces[[q + 1]]; num = Expand[Numerator[Rq]]; den = Expand[Denominator[Rq]];
     p = Exponent[den, xi]; c = Coefficient[den, xi, p];
     If[Expand[den - c xi^p] =!= 0, Print["DENOMINATOR NOT MONOMIAL ", Factor[den]]];
     WriteString[s, "k ", q, " p ", p, "\n", StringRiffle[ToString[#, InputForm] & /@ CoefficientList[num/c, xi], " "], "\n"],
     {q, 0, kmax}];
  val = Integrate[Sum[pieces[[q + 1]] Log[1 - xi]^q, {q, 0, kmax}], {xi, 0, 1}];
  val = Simplify[val];
  Print["COEFF SA^", j, " SB^", k, " : ", InputForm[val], "   N = ", N[val, 25]];
  total += val m,
  {m, mons}];
Close[s];
Print["RESULT K(SA,SB) = ", InputForm[Expand[total]]];
Print["time ", AbsoluteTime[] - t0];
