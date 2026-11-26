(* usage: export_HL.wl base : writes base_HLexport.txt with, for k = 0..3, the rational function
   R_k(xi) = N_k(xi)/xi^p_k  (after Together) of  Hx = sum_k R_k(xi) Log[1-xi]^k  *)
base = $ScriptCommandLine[[2]];
hL = Get[base <> "_HL.m"];
out = {};
Do[Rk = Together[Coefficient[hL, LL, k]];
   num = Expand[Numerator[Rk]]; den = Expand[Denominator[Rk]];
   (* denominator must be c xi^p *)
   p = Exponent[den, xi]; c = Coefficient[den, xi, p];
   If[Expand[den - c xi^p] =!= 0, Print["DENOMINATOR NOT A MONOMIAL: ", Factor[den]]];
   AppendTo[out, {k, p, CoefficientList[num/c, xi]}], {k, 0, Exponent[hL, LL]}];
s = OpenWrite[base <> "_HLexport.txt"];
Do[WriteString[s, "k ", o[[1]], " p ", o[[2]], "\n"]; WriteString[s, StringRiffle[ToString[#, InputForm] & /@ o[[3]], " "], "\n"], {o, out}];
Close[s];
Print["exported ", Length[out], " log powers; degrees ", Length[#[[3]]] & /@ out, " p ", #[[2]] & /@ out];
