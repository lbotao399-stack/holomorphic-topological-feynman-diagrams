# hCS on C^{3|2}: Q_4 and (Q²)_{ħ⁴} = 0

All arithmetic is modulo p₁ = 2³¹−1 or p₂ = 2147483629; the denominators (factorials ≤ 599!, spred rationals) are invertible.

## 1. Operator

Q = Σ_n ħⁿ Q_n, Q_0 = i Q̂_0, Q_n = i/(8π³)ⁿ Σ_{Γ∈R_n} κ_Γ Q̂_Γ, κ_Γ = σ_Γ/|Aut Γ| [src six-d-Q3/README §Ordered Wick weights; `wick_weights.py`]

Q̂_Γ O = ∫d²θ I_Γ(∂_out; ∂_in) Π_v f_{legs(v)} Π_{v∈out} C^{D_v}(θ) δ_{A_1}(θ)…δ_{A_{n+1}}(θ) O
- C = c + θ₁g₁ + θ₂g₂ + θ₁θ₂b
- δ/δC = θ₁θ₂∂_c + θ₂∂_{g₁} − θ₁∂_{g₂} + ∂_b; δ_{n+1} acts first
- legs(v) = (in-edges, output, out-edges), each in edge order; outputs in vertex order [src verify/vops.py]

I_Γ = (−1)^{n(n+1)/2} ∫_{ΠΔ} g_Γ(ξ;λ) Π_j Π_{a∈α_j} (p_j)_a dξ   [(−1)^{10} = +1 at n = 4]
- p_j = y_{e_j} = Σ_v c_{jv}(ξ) λ_v, e_j = (v_j, root), Bᵀy = λ
- Σ_j c_{jv} = 1 [drv htk CONVENTIONS §3, §5; chk `kern.moment_polys`]

All terms of (Q²)_{ħ⁴} carry i²/(8π³)⁴:
(Q²)_{ħ⁴} = 0 ⟺ Σ_{a+b=4} Σ κ_{Γa} κ_{Γb} Q̂_{Γa} Q̂_{Γb} = 0, κ_{Q0} = 1 [drv]; 15 + 15 + 5 + 5 + 4 = 44 compositions.

Pairing: ⟨Φ,X⟩ = [ε₁…ε_N] X|_{letter(t,D,α) ↦ Σ_{i:t_i=t} ε_i A_i^D u_i^α} [src vops.py].

Composition kernels at slot momenta u_{S(s)} [drv; `comp.key_kernel`]:
- Q̂_Γ Q̂_0: I_Γ(λ_v = u_{S(v)})
- Q̂_0 Q̂_Γ: I_Γ(λ_{v*} = u_{S(q₁)} + u_{S(q₂)})
- Q̂_{Γa} Q̂_{Γb}: ∫dξ_a g_a(ξ_a;u) · I_{Γb}(λ^b_v = p^a_{j(v)}(ξ_a) or u_{S(v)}); symbolic inner kernel, or joint ∫dξ_a dξ_b when the inner graph carries jets.

## 2. The 15 four-loop graphs

Root = vertex 9; graphs from `rooted_trees.graphs(4)`; probe = output momenta cycling e₁, e₂, e₃ [chk test_kern_n4.log].

| Γ | tree edges | marks | \|Aut\| | σ | κ | plan | I_Γ(probe) |
|---|---|---|---|---|---|---|---|
| 1 | 01 12 23 34 45 36 27 18 | 05678 | 2 | + | 1/2 | tight | −683/447897600 |
| 2 | 01 12 23 34 45 46 27 18 | 05678 | 4 | − | −1/4 | tight | 7/3359232 |
| 3 | 01 12 23 34 25 56 57 18 | 04678 | 8 | − | −1/8 | tight | 427/44789760 |
| 4 | 01 12 23 34 25 56 67 18 | 04578 | 2 | + | 1/2 | tight | −1441/134369280 |
| 5 | 01 12 23 34 25 16 07 78 | 04568 | 1 | − | −1 | tight | −3139/2239488000 |
| 6 | 01 12 23 34 35 16 07 78 | 04568 | 2 | + | 1/2 | tight | 7/3359232 |
| 7 | 01 12 23 34 35 26 07 78 | 04568 | 2 | − | −1/2 | tight | −7/3359232 |
| 8 | 01 12 23 14 45 56 07 78 | 03468 | 2 | − | −1/2 | tight | −2743/223948800 |
| 9 | 01 12 23 24 45 56 07 78 | 03468 | 1 | + | 1 | tight | 5843/3359232000 |
| 10 | 01 12 03 34 45 56 67 78 | 02468 | 2 | − | −1/2 | tight | 12541/6718464000 |
| 11 | 01 12 23 14 05 56 67 58 | 03478 | 2 | + | 1/2 | non-tight | −37/19906560 |
| 12 | 01 12 23 14 05 56 67 68 | 03478 | 2 | − | −1/2 | tight | 787/895795200 |
| 13 | 01 12 23 14 05 56 67 78 | 03468 | 1 | + | 1 | tight | 1613/4478976000 |
| 14 | 01 12 23 24 05 56 67 68 | 03478 | 8 | + | 1/8 | tight | −77/44789760 |
| 15 | 01 12 23 24 05 56 67 78 | 03468 | 2 | − | −1/2 | tight | 77/44789760 |

Σ 1/|Aut| = 8; multiplicities of |Aut| = 1, 2, 4, 8: 3, 9, 1, 2 [chk; = src README].

## 3. Kernels

Tight plans (14 graphs): g = ±Π_m det L_m [src htk spred]; vectorized over plane-wave assignments mod p; ∫_{Δ_m} Π b_j^{a_j} = Π a_j!/(m+Σa_j)!; last two factors integrated bilinearly [`kern.integrate_factors`].

Non-tight plan (Γ₁₁), block lower-triangular density [drv `blocktri.py`]:
- ℓ_i = w_i y_red^{(m)} − m̂_i + Q_{i,j_i}
- E_{>m} is a cotree of G_m; T_m = E(G_m) ∖ E_{>m} spanning tree [chk assert]
- y_red^{(m)} = Σ_{E>m} C ℓ^{(>m)} + T^{(m)}, T^{(m)} = Σ_{v∈p-side of e_m in T_m} λ_v^{G_m}
- ⋀ dℓ = ⋀_m R^{(m)}, R^{(m)}_{i,a} = y_red,a dw_i + w_i dT^{(m)}_a − dm̂_{i,a} + dQ_{i,j_i,a} ∈ span{dξ_{≤m}}
- g = σ Σ_{(U_m)} ε(U) Π_m det R_m[U_m]
- Γ₁₁: (|R_m|) = (3,3,6), (|C_m|) = (4,2,6), g = det R₃[C₃] Σ_{a∈C₁} ε_a det R₁[C₁∖a] det R₂[{a}∪C₂], σ = +1

| check | result | log |
|---|---|---|
| n = 1–3 all graphs: XP vs spred at 3 probes; moment (p₀)₁(p_last)₃ vs spred at 2 probes | equal | test_kern_n3.log |
| n = 4: 15/15 at probe; 14 tight at 2 more probes; moments 14 tight | equal | test_kern_n4.log |
| Γ₁₁ moments (p₀)₁(p₄)₃, (p₁)₂(p₂)₂ at 2 probes | equal | test_moment411.log |
| block-triangular det vs spred local factors (3.1, 3.4, 4.1, 4.10), vs Kirchhoff density at 6 points per graph | equal; ratio ±1 | test_blocktri.log |
| Γ₁₁ block-triangular vs spred generic: −37/19906560, −4898417/5184000 | equal; 0.1 s vs 26 s | test_blocktri.log |
| n = 5: 49/49 kernels = five_loop_periods.jsonl; κ = wick_weight; 6 non-tight by block triangulation | equal | test_kern5.log |

## 4. Identity checks

E(Γa, Γb) = ⟨Φ, Q̂_{Γa} Q̂_{Γb} O⟩ [`run_nilp.py`]; residual = Σ κκ′E, κ = σ/|Aut|.

| n | g | O | functionals | F | p | residual | rank | κ_n solved | log |
|---|---|---|---|---|---|---|---|---|---|
| 2 | so5 | bbb | c³g₁g₂ | 3 | p₁ | 0 | 2/2 | (−1/2, 1/2) = σ/\|Aut\| | nilp2_so5_bbb(_joint).log |
| 2 | so5 | ∂₁b·b·∂₂∂₃g₁ | c⁴g₁ | 3 | p₁ | 0 | 2/2 | same | nilp2_so5_bbg1_jets.log |
| 2 | su2 | ∂₁b·b·g₁ | c⁴g₁ | 3 | p₁ | 0 | 2/2 | same | nilp2_su2_bbg1_jet.log |
| 3 | so5 | bbbb | c⁴g₁g₂ / c⁵b | 7 / 6 | p₁ | 0 | 5/5 | (1/2, −1/8, 1, −1/2, −1/2) = σ/\|Aut\| | nilp3_so5_bbbb*.log |
| 3 | so5 | g₁g₂bb | c⁶ | 6 | p₁ | 0 | 5/5 | same | nilp3_so5_g1g2bb_c6.log |
| 3 | su2 | g₁g₂·∂₁b·b (joint) | c⁶ | 6 | p₁ | 0 | 5/5 | same | nilp3_su2_g1g2bb_jet.log |
| 4 | so5 | bbbbb | c⁵g₁g₂ | 2 | p₁ | 0, 0 | — | — | nilp4_so5_bbbbb_F2.log |
| 4 | so5 | bbbbb | c⁵g₁g₂ | 16 | p₂ | 0 (16/16) | 15/15 | κ_4 = σ/\|Aut\| for all 15 Γ | nilp4_so5_bbbbb_F16_p2.log |
| 4 | so5 | g₁g₂bbb | c⁷ | 1 | p₁ | 0 | — | — | nilp4_so5_g1g2bbb_c7_F1.log |
| 4 | so5 | bbbbc | c⁷ | 1 | p₁ | 0 | — | — | nilp4_so5_bbbbc_c7_F1.log |

- Sensitivity [chk analyze_*.log]: 704/704 composition values nonzero (44 × 16); flipping or deleting any one κ_{Γ₄} gives nonzero residual (15/15); κ = +1/|Aut| fails.
- Pure-ghost sign [drv; chk compare_bbbbc_g1g2bbb.log]: Π_{j=1}^{4}(θ_{j1}−θ_{51})(θ_{j2}−θ_{52}) has coefficient +1 for b⁴c, −1 for g₁g₂b³; E(bbbbc) = −E(g₁g₂bbb) for 44/44 compositions; n = 3: E(bbbc) = −E(g₁g₂bb) for 14/14 (= README C10, independent engine).
- su(2) [chk]: every composition vanishes on jet-free O at n = 2 (g₁g₂b, bbb), n = 3 (g₁g₂bb), n = 4 (g₁g₂bbb, F = 2), and at n = 4 with one jet (g₁g₂·∂₁b·bb, F = 1); no nontrivial order-4 su(2) test [logs nilp2_su2_g1g2b_jetfree, nilp3_su2_g1g2bb_jetfree, nilp4_su2_g1g2bbb_F2_jetfree, nilp4_su2_g1g2bbb_jet_F1].
- Timing: order 3 6–8 s (verifier engine 966 s [src VERIFY C7]); order 4 608 s (F = 2), 4192 s (F = 16).

## 5. Physical content

- [chk single_op_n4_so5.log] so(5), cochain level: every Q̂_Γ (n = 4) ≠ 0 on bbbbb (paired with c⁵b, c⁴g₁g₂), on g₁g₂bbb and bbbbc (with c⁶); all 15 = 0 on g₁bbbb (with c⁵g₁).
- [chk single_op_n4_su2.log, single_op_jets_su2.log, single_op_scan_su2.log] su(2): Q̂_Γ (n = 4) = 0 on the jet-free monomials tested and on g₁g₂·∂₁b·bb; n = 3: 0 on all jet-free monomials, ≠ 0 with one jet.
- [src lifting/README] all orders: E₁ = E_∞ for gl_N, sl_N at field length p ≤ N; zero-jet sl₂ sector survives; field-length-2, -3 sectors survive for simple g; strict scalar cocycles killed by every Q_n; Q_n ≠ 0 on sl₂ cochains at every n.
- [open] lifted classes at ħ⁴ outside these ranges (E₄ page at low weight, so(5) or su(2) with jets).

## 6. Q_5 and order 5

- [chk test_kern5.log] 49/49 kernels assembled; κ = σ/|Aut| = `wick_weight`.
- [open] order-5 identity: O = b⁶, 8 final letters (8! assignments per functional set); Q̂_0Q̂_5 ≈ 2.5×10⁵ histories, Q̂_5Q̂_0 ≈ 1.5×10⁶; colour tensor 8⁸ entries per history; 15 coordinates exceed the 12-coordinate int64 key packing at 5 bits.

## 7. 10d (hCS on C⁵)

[open] no 10d operator layer (field content, leg rules, odd-free pairing); kernel data: 28 three-loop graphs [src circuit-general/assembly/results/cone5_o2.jsonl].

## 8. Failures and limits

1. spred with symbolic λ (18 parameters, Γ₄.₁): 8.1 GB after 6 min, killed; replaced by numeric evaluation vectorized over assignments.
2. Symbolic three-loop kernels with two moment insertions > 4.5 GB; replaced by joint ξ_a ∪ ξ_b integration (validated n = 2, 3).
3. Nontrivial su(2) test at n = 4 needs ≥ 2 jets; Γ₄.₁₀ greedy integration has a 4.3×10⁵-monomial intermediate (5.2 GB at 32 assignments); not run [open].
4. Order-3 so(5) with two jets (bbbg₁, F = 3) too slow under the current key structure; killed.
5. Q_4 has no closed-form colour-reduced expression; it exists as `comp.apply_graph` + `kern.GraphKernel`.

## 9. Files

- Engine: `kern.py`, `blocktri.py`, `comp.py`, `run_nilp.py`.
- Analysis: `analyze.py`, `compare_terms.py`, `single_op_scan.py`, `single_op_n4.py`, `single_op_jets.py`.
- Kernel tests: `test_kern.py`, `test_blocktri.py`, `test_moment411.py`, `test_kern5.py`.
- Copied unchanged from six-d-Q3: `rooted_trees.py`, `wick_weights.py`.
- Logs `*.log`; `nilp4_so5_bbbbb_F2_oldcode_partial.log` = partial run with the first integration code, agrees where it overlaps.

## 10. Spot-check (main session)

- `run_nilp.py 3 so5 b,b,b,b c,c,c,c,g1,g2 3 13`: residual [0, 0, 0], rank 3 [chk spot_nilp3_so5_bbbb_F3_s13.log].
- `run_nilp.py 4 so5 b,b,b,b,b c,c,c,c,c,g1,g2 2 23 --prime 2147483629`: 44 compositions, residual [0, 0], rank 2, 474 s [chk spot_nilp4_so5_bbbbb_F2_s23_p2.log].
