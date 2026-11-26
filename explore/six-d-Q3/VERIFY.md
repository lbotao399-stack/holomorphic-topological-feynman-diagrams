# VERIFY: six-d-Q3

Scripts and logs: `P/explore/six-d-Q3/verify/`. No explorer script is run; nothing is imported from `q3lib.py`, `opcheck.py`, `s2code/`, `htk`. Explorer files were read for conventions only (graph labels, move order, sources of `numerics_simplex.log`, colour/functional set-up of `q3check.py`).

`P/explore/six-d-Q3/REPORT.md` does not exist (directory listing at verification time). The logical gaps below are read off the structured report and the scripts.

Notation: $\langle ijk\rangle=\det(\lambda_i,\lambda_j,\lambda_k)$; $I_\Gamma$ = S2 eq-kernel-def,
$$I_\Gamma=\varepsilon_T\det B_T\int_{\omega_m}e^{-\sum_ey_e\cdot z_e}\bigwedge_{e\notin T}d^3y_e ,\qquad g(r)=\varepsilon_T\det B_T\,\det\Big[\frac{\partial y^a_c}{\partial r_k}\Big]_{(c,a),k}\operatorname{sgn}\det\Big[t;\frac{\partial t}{\partial r}\Big],$$
cotree $c$ increasing, $a=1,2,3$; $\hat Q_0=-iQ_0$, $Q_1=\frac{i}{16\pi^3}\hat Q^I_{\rm box}$, $Q_2=\frac{i}{128\pi^6}(\hat Q^I_A+\hat Q^I_B)$ (S2 (1.4), (1.5), (7.15)), $Q_3=\frac{i}{512\pi^9}\sum_\Gamma\kappa_\Gamma\hat Q^I_\Gamma$.

---

## C1 (graphs of $Q_3$) — confirmed

`v_enum.py` (all edge subsets of $K_7$ that are forests, all spoke multisets with $\deg_F(v)+s_v\le3$, $\sum s_v=n+p$; bound $2|I_S|\le3|V_S|-4$ on all $2^8-1$ vertex subsets, induced with multiplicity; connectivity; classes by VF2 with $\mathcal O$ coloured; $|\mathrm{Aut}|$ by brute force over $S_7$).
```
n=3 p=1: forests 16807, candidates (deg<=3) 679770, passing 13230
n=3 p=2..7: passing 0
n=3: 5 iso classes (O fixed)
  ~G1: degO 4, G-O tree, degrees [3,3,2,2,2,2,2], outputs 5, labelled 2520, |Aut| 2, chains [1,1,2,3,3]
  ~G2: ... labelled 2520, |Aut| 2, chains [1,2,2,2,3]
  ~G3: ... labelled 5040, |Aut| 1, chains [1,1,2,3,3]
  ~G4: ... labelled  630, |Aut| 8, chains [2,2,2,2,2]
  ~G5: ... labelled 2520, |Aut| 2, chains [1,2,2,2,3]
sum labelled = 13230 = sum 7!/|Aut|
explorer labelled graphs: |Aut| 2,2,1,8,2; stated generators are automorphisms and generate groups of orders 2,2,1,8,2;
  degree-2 vertices [1,2,4,6,7],[1,2,4,6,7],[1,2,5,6,7],[1,2,4,6,7],[1,2,5,6,7]
```
Same code: $n=1$: box only; $n=2$: two classes, chains $[1,3,3]$ (A), $[2,2,3]$ (B), $|\mathrm{Aut}|=2$ (S2 eq-AB).

## C2 ($p=1$ forced) — confirmed

Component $i$ of $\Gamma\setminus\mathcal O$: $k_i$ vertices, $k_i-1$ edges, $s_i$ spokes; $S$ = component $+\,\mathcal O$, $I_S\neq\emptyset$ (connectivity):
$$2(k_i-1+s_i)\le3(k_i+1)-4\iff2s_i\le k_i+1,\qquad\sum_is_i=\deg\mathcal O=n+p,\ \sum_ik_i=2n+1\ \Rightarrow\ 2(n+p)\le2n+1+p\iff p\le1 .$$
Brute force (C1 code, no use of the derivation): $p\ge2$ gives 0 graphs for $n=1,2,3$.

## C3 (kernels) — confirmed; one intermediate number convention-dependent

`vk.py`, `v_check_kernels.py`: nodal solver $L=B^{\mathsf T}t^{-1}B$, $\varphi=L^{-1}\lambda$, $y=t^{-1}B\varphi$; $\partial y/\partial t_j$ in closed form; chart $t(r)$ built backwards from the final box ($T_{\rm box}=1$) through S2 eq-chain, eq-bubble with forward-mode exact derivatives; 40 random rational points per graph; sources at the outputs (compact $J$) and at all vertices (product of local factors through own reduced networks).
```
G1..G5 (and box, A, B): 40 points; g tree-independent (5 trees) True; g == I-sign * J True;
  general sources: ratios recomputed from t(r) == r True, reduced-network currents == full currents True,
  g == s * prod(local factors) True    (s = +1,-1,+1,+1,-1;  box +1, A -1, B +1)
S2: eps_T det B_T = -1 (box, tree e2 e3 e4), +1 (A, cotree [0,4]), -1 (B, cotree [0,4])
```
$I_\Gamma=+\!\int J_1,+\!\int J_2,+\!\int J_3,-\!\int J_4,-\!\int J_5$ reproduced pointwise including sign.

Move-cotree values $\varepsilon_T\det B_T$:
```
cotree        S2 order (increasing)   order as listed   explorer
G1 [0,4,7]          -1                     -1             -1
G2 [0,6,4]          +1                     -1             -1
G3 [0,4,7]          +1                     +1             +1
G4 [0,6,4]          -1                     +1             +1
G5 [0,4,7]          -1                     -1             -1
```
The stated values for G2, G4 hold for the cotree ordered $(0,6,4)$; with S2's $\varepsilon_T$ (tree edges moved in front, cotree increasing) they are $+1,-1$. Densities are unaffected.

## C4 (displacement maps) — confirmed

$-\sum_ey_e\cdot z_e=\sum_v\lambda_v\cdot X_v$ exact at 40 random rational $(r,z,\lambda)$ per graph with own currents (`v_check_kernels.log`: `X map True` for G1..G5, box, A, B with S2 eq-res-X).
Convexity (`v_convex.py`): $\sum_jc_{v,j}=1$ identically (sympy), $c_{v,j}\ge0$ at 2000 rational points of $D_\Gamma$ (0 negative; minima $1.5\cdot10^{-5},5.1\cdot10^{-6},1.8\cdot10^{-6},9.5\cdot10^{-6},1.5\cdot10^{-5}$); by inspection every map is a nested combination with weights $r,1-r\in(0,1)$.

## C5 (zero-displacement values) — confirmed

`v_integrals.py`: multilinear expansion of $J_\Gamma$ in formal brackets, each coefficient integrated by sympy iterated integrals over $D_\Gamma$ (no monomial formula).
```
G1..G5: 36 positional triples, denominators 622080, 622080, 622080, 622080, 207360;
        equal to json intJ as polynomials in independent bracket symbols: True (5/5)
        G1 list of the statement == json: True
I(lambda;0) at the explorer sources: 52924313/311040, -85183/25920, 58733/25920, 7757/1440, -18409/51840  (5/5 equal)
factor_list in independent brackets: one irreducible factor of degree 3 (5/5)
```
Numerical (`v_numeric.py`, own float Kirchhoff, brute-force $g$, Sobol $8\times2^{20}$ in ratio coordinates):
```
G1 170.1504455 +- 0.0022   exact 170.1527553   4.9 digits, 1.06 sigma
G2 -3.284684213 +- 0.0013  exact -3.286381173  3.3 digits, 1.29 sigma
G3 2.266798598 +- 0.00047  exact 2.265933642   3.4 digits, 1.83 sigma
G4 5.386803288 +- 4.2e-05  exact 5.386805556   6.4 digits, 0.05 sigma
G5 -0.354322914 +- 0.00066 exact -0.3551118827 2.7 digits, 1.19 sigma
```
The values quoted as "htk.kernel_exact ... 7757/1440, −18409/51840" are $I_\Gamma=-\!\int J$ for G4, G5 ($\int J_4=-7757/1440$, $\int J_5=18409/51840$).

## C6 (Wick tensors, $Q_3$) — confirmed (signs through C7)

`v_tensors.py`: $T_\Gamma$ regenerated from the edge lists with legs (in-edges in edge order, output, out-edges), spoke $j\to A_j$ in edge order: equal to the explorer strings 5/5; the same rule gives S2's $T_{\rm box}$, $T_A$, $T_B$ ((1.4), (1.6), (1.7)).
$$\frac{i}{512\pi^9}\kappa_\Gamma I_\Gamma=\frac{i}{1024\pi^9}\,2\kappa_\Gamma\,\mathrm{sgn}_I\!\int J_\Gamma:\quad2\kappa\,\mathrm{sgn}_I=(-1,1,-2,-\tfrac14,1),\qquad|\kappa_\Gamma|=\tfrac1{|\mathrm{Aut}\,\Gamma|}=(\tfrac12,\tfrac12,1,\tfrac18,\tfrac12).$$
4 inputs, 5 outputs, three brackets per term: $(\partial_1,\partial_2,\partial_3)^3$. $\kappa$ from C7 (independent engine).

## C7 (nilpotency at $\hbar^3$) — confirmed

$$\{Q_0,Q_3\}+\{Q_1,Q_2\}=0\iff\sum_\Gamma\kappa_\Gamma\{\hat Q_0,\hat Q^I_\Gamma\}=-\tfrac14\{\hat Q^I_{\rm box},\hat Q^I_A+\hat Q^I_B\}$$
($i\cdot\frac{i}{512\pi^9}=-\frac1{512\pi^9}$, $\frac{i}{16\pi^3}\frac{i}{128\pi^6}=-\frac1{2048\pi^9}$).

Engine `vops.py` (own code): tagged letters in a Grassmann algebra with $\theta_{1,2}$; $C=c+\theta_1g_1+\theta_2g_2+\theta_1\theta_2b$; $\delta/\delta C=\theta_1\theta_2\partial_c+\theta_2\partial_{g_1}-\theta_1\partial_{g_2}+\partial_b$ (S2 eq-fd); $\hat Q_0c=\frac12fcc$, $\hat Q_0g_i=fcg_i$, $\hat Q_0b=fcb+fg_1g_2$ (from S2 eq-Q0, $Q_0\theta=-\theta Q_0$); derivation in insertion form; colours by einsum over $\mathfrak{so}(5)$ built from $[L_{ij},L_{kl}]$ (Jacobi checked); jets by momenta; kernels $J_{\rm box}=\langle123\rangle$, $J_A$, $-J_B$ (S2 (1.6),(1.7)), $\pm J_\Gamma$ (C3), $X$ maps S2 (1.8) and C4; exact polynomial integration; test functionals: homomorphism letter$\mapsto\sum_i\epsilon_ia_i^Du_i^m$, top coefficient; test operators with random integer colour vectors.

Validation against S2 (`v_validate.log`, `v_validate_seed2.log`):
```
Q0^2 = 0 on b, g1 b, c g2 (4 functional sets)
{Q0,Q_box} = 0 on g1g2, bb (2 sets), g1b, d_1 g1 g2, b d_3 g2   (each term nonzero, sum 0)
n=2: O = g1g2b, bbb, bbg1, d_1g1 g2 b, bb d_2g2:  rank 2/2, (kappa_A, kappa_B) = (1/2, 1/2)  [= S2 eq-kAkB k = (2,-2)]
```
Order $\hbar^3$, operators not used by the explorer (`v_nilp3_main.log`):
```
so(5) O = bbbb          (6 eq.)  rank 5, augmented 5, kappa = (-1/2, 1/2, -1, 1/8, -1/2)
so(5) O = bbbg1         (4 eq.)  rank 4, augmented 4, residual of (-1/2,1/2,-1,1/8,-1/2): 0,0,0,0
so(5) O = bbbg2         (3 eq.)  rank 3, augmented 3, residual 0,0,0
so(5) O = (d_2 b)bbb    (4 eq.)  rank 4, augmented 4, residual 0,0,0,0
so(5) O = bbb(d_3 g1)   (3 eq.)  rank 3, augmented 3, residual 0,0,0
all 20 equations: rank 5, augmented 5, kappa = (-1/2, 1/2, -1, 1/8, -1/2), residual 0 (20/20)   [966 s]
```
Explorer's operator with own engine (`v_nilp3_c10.log`): $O=g_1g_2bb$, 6 equations, rank 5/5, $\kappa=(-\frac12,\frac12,-1,\frac18,-\frac12)$.
Second algebra (`v_nilp3_so6.log`): $\mathfrak{so}(6)$, $O=bbbb$, 6 equations, rank 5/5, $\kappa=(-\frac12,\frac12,-1,\frac18,-\frac12)$.

## C8 (three-letter operators) — confirmed

`v_nilp3_c8.log`, `v_c8_detail.log`, $O=bbb$, five-$c$ functionals (3 sets) and type-mismatch controls: all 14 terms ($\hat Q_\Gamma\hat Q_0$, $\hat Q_0\hat Q_\Gamma$, box$\cdot X$, $X\cdot$box) are 0.
```
G1 Q0(bbb): part f c b -> 13654980383/3240,   part f g1 g2 -> -13654980383/3240,   both -> 0
G2:          53472198043/3888,                 -53472198043/3888
G3:         -244160257597/12960,               244160257597/12960
G4:          70041013/10,                      -70041013/10
G5:         -37115538109/2160,                 37115538109/2160
box.A, A.box, box.B, B.box on bbb: every binding cancels at the Grassmann level (0 surviving histories)
```

## C9 (sign rule, fitted) — confirmed as a labelling-consistent rule

`v_relabel.py`: (a) $\sigma$ on the reference labellings: box $+1$, A $+1$, B $+1$, G1..G5 $(-1,+1,-1,+1,-1)$; with C7 and the S2 values: $\kappa_\Gamma=\sigma_\Gamma/|\mathrm{Aut}\,\Gamma|$ for box ($\frac12$), A, B ($\frac12,\frac12$), G1..G5.
(b) Operator-level invariance without nilpotency: 6 random relabellings of each of box, A, B, G1..G5 (vertex permutation, edge reversals, edge-order permutation), $\varepsilon=I_{\Gamma'}/I_\Gamma$ by brute force at the same Schwinger point, $\hat Q^I_{\Gamma'}$ built from the relabelled leg pattern:
```
48 relabellings, 96 comparisons (O = g1g2 / g1g2b / g1g2bb, with and without one input derivative):
sigma' <Q_G' O> == sigma <Q_G O> in all: True; comparisons with value != 0: 92
```
Not derived from the Feynman rules (as stated).

## C10 ($bbbc$ vs $g_1g_2bb$) — confirmed

`v_nilp3_c10.log`: same colour vectors $w_{1..4}$ in $g_1^{w_1}g_2^{w_2}b^{w_3}b^{w_4}$ and $b^{w_1}b^{w_2}b^{w_3}c^{w_4}$, 6 own functional sets (six $c$):
```
term by term (all 14 compositions) T(bbbc) == -T(g1g2bb): True (6/6);  == +T: False
O = bbbc: rank 5, augmented 5, kappa = (-1/2, 1/2, -1, 1/8, -1/2)
```
The $\theta$-translation explanation is not derived (as stated).

---

## Logical gaps

1. C1, C2: S2 eq-count (5.11) (including the subgraph bound as the criterion for nonvanishing) and eq-theta-rule (5.13) are inputs, not re-derived. The five surviving graphs have nonzero kernels (C5).
2. C3, C4: pointwise identities at finitely many random rational points (explorer 12, here 40 per graph and 5 spanning trees) are identity tests of rational functions of $r$; the brute-force side is not reduced symbolically.
3. C3: the move-cotree values $\varepsilon_T\det B_T$ for G2, G4 refer to the cotree order $(0,6,4)$, not to S2's increasing order.
4. C5: the explorer's direct Schwinger-simplex integral (0.2–1.5 digits) is not a confirmation; the 3.3–5.2 digit check (explorer) and the 2.7–6.4 digit check (here) both use the ratio chart.
5. C6, C7: $Q_3$ is the unique solution inside the ansatz $\frac{i}{512\pi^9}\sum_\Gamma\kappa_\Gamma\hat Q^I_\Gamma$ (S2 leg pattern, $|\kappa|$ predicted by S2 (5.15)); the Koszul signs of $W_\Gamma$ are not derived from the Feynman rules (S2 derives only the box, (6.7)). The order-$\hbar^3$ relation is linear in $Q_1$; the overall sign of $\kappa$ is tied to the sign of $Q_1$ in S2 (1.4), while the order-$\hbar^2$ check fixes $Q_2$ relative to $Q_1^2$ only (S2 Section 8.2 item 4: $\theta_1\to-\theta_1$ gives $Q_n\to(-1)^nQ_n$).
6. C7: random test functionals and finitely many test operators; an operator identity is supported, not proved. Explorer: $\mathfrak{so}(5)$ only; the explorer's operator code is validated against S2's own code (salg.py copied), which does not test S2's code itself.
7. C7, C9, C10: the explorer's runs share one operator code (opcheck.py); C9 and C10 are not independent of it.
8. C9, C10: fitted rule and unexplained sign relation, as stated in the report.
