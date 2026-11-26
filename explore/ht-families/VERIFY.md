# VERIFY: ht-families

Scripts and outputs: `explore/ht-families/verify/` (run from that folder with `python3 -u <script>`; outputs `<script>.out`).
No explorer script (`hf_*.py`, `a_*.py`, `b*.py`, `c_*.py`, `d_*.py`, `t_*.py`) was imported or run.
`REPORT.md` does not exist in `explore/ht-families/` (nor anywhere under P); the gaps below are read off the structured report.

## 0. Independent layer

`vcore.py` (written from CONVENTIONS.md 1-2 only):

$$B_{ev}=\delta_{v,e(0)}-\delta_{v,e(1)},\quad \Lambda=B^{\mathsf T}t^{-1}B,\quad \varphi=\Lambda^{-1}\lambda,\quad y=t^{-1}B\varphi,$$

DFS tree from the root, fundamental cycles, first row flipped so that $\det\binom{B^{\mathsf T}}{C}>0$;
exterior algebra on bitmasks; $\mathrm{Pf}(A)=\sum_{j\ge2}(-1)^jA_{1j}\mathrm{Pf}(A_{\hat1\hat\jmath})$;

$$D(t)=\frac{(-1)^{e_0}}{t_{e_0}}\Big[U^{-H}\,U^{-1/2}\,\mathrm{Pf}\big(-\tfrac12dM\,M^{-1}dM\big)\wedge\textstyle\prod_{i}\prod_r\psi^i_r\Big]_{dt_{\hat e_0}} .$$

Float version: $A_{rs}=-\frac12\sum_k\beta_{rk}\wedge\beta_{sk}$, $\beta_{rk}=\sum_q dM_{rq}L_{qk}$, $M^{-1}=LL^{\mathsf T}$ (eigh), top coefficient as $\sum$ of determinants.

`vsp.py`: own SP decomposition (parallel merge of equal end points, series merge at a degree-2 vertex; deterministic or random order, flattened or not); ratio map

$$R_{c_j}=x_{s,j}R_s\ (\text{S-node}),\qquad R_{c_i}=R_p/w_{p,i}\ (\text{P-node}),\qquad R_{\rm root}=1,\qquad g(\xi)=D(t(\xi))\,\big|\det[t;\partial_\xi t]\big| ;$$

QMC: scrambled Sobol, Dirichlet($\alpha$) per block, $\alpha_b=\max(\frac12,e_b+1)$ with $e_b$ the leading exponent of $b$ in $g$ probed by the own float density; points with $\max t/\min t>10^7$ ($|I|=9$) or $>10^6$ ($|I|=13$) re-evaluated with the exact density.
`vdir.py`: $\int_{\Delta^m}\prod_jb_j^{a_j}=\prod_j\Gamma(a_j+1)/\Gamma(m+1+\sum a_j)$, exact half-integer Gamma.
`v_rec.py`: own implementation of the rule of a3 from its statement.

Engine (htk) use: (i) `spred` ratio density $g$ at numeric sources in V1/W1, accepted only after the pointwise identity with `vcore` (sign and square) at rational points, then integrated by `vdir`; (ii) `spred` plans as input structures for `v_rec`, `v_a6comp`; (iii) one cross-check S6.

Self-tests (`v0_selftest.out`):

```
S1 G2: D(e0=N-1) == D(e0=0) == D(e0=3): True ; S1 G4: True
S2 G2, G4: pullback of Omega to own ratio coordinates == D det[t; dt/dxi] at 3 points: True
S3 float/exact: rel 1.3e-13 (G2), 7.5e-16 (G4)
S6 own density == htk.brute.density (sign and square): True (4 points)
S7 G2, G4: own ratio map bijective (5 random decompositions, exact): True
```

(S4/S4', $\int|\det[t;\partial t]|(\sum t)^{-N}=1/(N-1)!$ by QMC: $10^{-3}$ to $3\cdot10^{-3}$, heavy-tailed integrand; not used.)
Float error of the batched density versus $\max t/\min t$ (measured, 60 points per graph): $\le7\cdot10^{-8}$ below $10^{9}$ and $10^{-2}$ at $10^{12}$ for G4; $1.4\cdot10^{-4}$ at $10^{7}$ and $4.5\cdot10^{-2}$ at $10^{8}$ for nested theta(1,3,3); hence the exact re-evaluation thresholds above.

## 1. Verdicts

| id | verdict | method |
|---|---|---|
| a1 | confirmed | re-derivation; exact Wick-moment evaluation of the $a$-integral vs own Pf |
| a2 | confirmed | exact: $C_{re}C_{se}=C_{rf}C_{sf}$ on chains; Pf invariant under $(t_e,t_f)\to(t_e+\delta,t_f-\delta)$ and $(dt_e,dt_f)\to(0,dt_e+dt_f)$; no $dx$ of plain chains in own ratio coordinates |
| a3 | confirmed | own implementation of the rule vs own Pf($C_{\rm rec}$), 20 plans, 2 points each |
| a4 | confirmed | own recursion insertion vs closed form, $k=2..6$, exact |
| a5 | refuted | $\rho_n$ product formula and $\int=\pi^{n/2}$ confirmed; $K_{\rm banana}$ closed form has the right modulus but the wrong sign for $H=1$ (theta(1,2,2), theta(2,2,1), banana(2,1,2,2,2), banana(2,2,2,2,2,2,1)) and $H=3$ (theta(1,4,4)) |
| a6 | refuted | count is 19 of 32 merges, not 20 of 30; the rule fails on 3 merges of the nested (2,1) graphs; component counts 4 and 10 confirmed |
| a7 | confirmed | exact 4-form identity in the 9 $dt_e$; second form to $2\cdot10^{-16}$ |
| a8 | confirmed | proof re-checked; exact $\rho\wedge\rho=0$ on Sym$_4$, a 12-dim family in Sym$_6$, two $d=2$, $n=6$ graphs and a nested (2,1) graph (the last three not in the explorer's set) |
| b1 | confirmed | exact (engine density verified pointwise, own Dirichlet), 20+3 source vectors per graph; own QMC |
| b2 | confirmed | exact at the 8 printed source sets; own QMC; banana closed form at 4 fresh sets |
| c1 | confirmed | exact $U(t(\xi))$ = monomial, all SP graphs of the six $(D,n)$, 5 decompositions each |
| c2 | confirmed (conditional) | Dirichlet $\pi$-count re-derived; exponent parity $\equiv-e_b(U)/2$ and $\min e_b=-\frac12$ on 9 graphs; hypothesis not proved |
| c3 | confirmed | re-derivation of $\phi$; identity on 25 decompositions of every SP graph; E distribution |
| c4 | confirmed | from b1, b2, d1, c3 |
| c5 | unverified | status open; not in scope |
| d1 | confirmed | exact Jacobian identity, $\Pi$ on the central square, 2D quadrature of the central integral at fixed ears, QMC |

## 2. Claims

### a1

$$t_e\hat A_e^2=t_eA_e^2+A_e\,dt_e\Theta_e\quad\big((dt_e\Theta_e)^2=0\big),\qquad \sum_eA_e\,dt_e\Theta_e=a^{\mathsf T}dM\,\theta=\sum_ra_rZ_r,\quad Z_r=\sum_sdM_{rs}\theta_s,$$

$$\pi^{-n/2}\!\int d^na\,e^{-a^{\mathsf T}Ma-a\cdot Z}=\det M^{-1/2}\sum_{k}\frac{(-1)^k}{k!}\big\langle(a\cdot Z)^k\big\rangle,\qquad\langle a_ra_s\rangle=\tfrac12(M^{-1})_{rs},$$

$\theta$-top part only from $k=n$. Check (`v_a.py` A1): $\big[\frac{(-1)^n}{n!}\sum_{r_1..r_n}\langle a_{r_1}\cdots a_{r_n}\rangle Z_{r_1}\cdots Z_{r_n}\big]_{\theta_1\cdots\theta_n}=\mathrm{Pf}(-\frac12dM M^{-1}dM)$, exact:

```
A1 theta(1,2,2) n=2: [True, True] ; G2 n=4: [True, True] ; G4 n=4: [True, True] ; banana(2,2,3,3,3) n=4: [True, True]
```

### a2

On a chain $\{e,f\}$: $C_{\cdot e}=\pm C_{\cdot f}\Rightarrow M=\dots+(t_e+t_f)C_{\cdot f}C_{\cdot f}^{\mathsf T}$, $dM=\dots+(dt_e+dt_f)C_{\cdot f}C_{\cdot f}^{\mathsf T}$.

```
A2 G0..G4: Pf(t) depends on chain only via (R, dR): True (all degree-2 chains) ;
           own ratio coords: #components with dx(plain) = 0 ; Pf, U unchanged under x(plain) change: True
```

Second sentence of a2 (x enters only through insertions of reduced edges) checked only in the form above plus a3.

### a3

`v_rec.py`, rule exactly as stated (variables $A_i=w_iJ+j_i$, $\Theta_i=w_i\Theta+\vartheta_i$; $\langle j_aj_b\rangle=(w_a\delta_{ab}-w_aw_b)/2t_m$; $(1+dt_mJ\Theta)$; Berezin at the right end; consumed insertions in creation order with the sign of moving past unconsumed odd ones), $C_{\rm rec}$ built from "loop $(m,a)$ = chain $a$ − chain $k$, new edge routed through chain $k$, rows by decreasing $m$":

```
A3 G0..G4 (large, small), theta(1,2,2) (large, small), banana(1,3,3,3,3), banana(2,2,3,3,3),
   nested theta(1,3,3), nested theta(2,2,3) (large, small): own recursion == Pf(C_rec), U == det M(C_rec):
   [(True, True), (True, True)] in all 20 plans
```

The four nested-graph plans are not in the explorer's test set.
Derivation re-done: $C_{\rm rec}$ flows $A_a=a_{(m,a)}$, $A_k=J-\sum_aa_{(m,a)}$; shift $j_a=a_{(m,a)}-w_aJ$ (unimodular), $\sum_iR_iA_i^2=t_mJ^2+\sum_iR_ij_i^2$ since $R_iw_i=t_m$, $\sum w_i=1$; $dt_m=t_m^2\sum dR_i/R_i^2=\sum_iw_i^2dR_i$; constrained covariance $\frac12\big(\mathrm{diag}(R^{-1})-R^{-1}R^{-\mathsf T}/\sum R^{-1}\big)=\frac{w_a\delta_{ab}-w_aw_b}{2t_m}$; $\det M_{\rm loc}=\prod_iR_i\sum_iR_i^{-1}=t_m^{k-1}/\prod w_i$.

### a4

```
A4 k=2..6 plain merge: own recursion insertion == s_k t^{k-1}/prod w dw_1..dw_{k-1} :Jhat^{k-1}:_t (2 points each): True
```

with $:X^p:_t=\sum_m\frac{(-1)^mp!}{m!(p-2m)!}(4t)^{-m}X^{p-2m}$, $\hat J^q=J^q+qJ^{q-1}\frac{dt}{2t}\Theta$.

### a5

$\rho_n$ part (`v_a.py` A5), $C_{\rm rec}$ rows $e_a-e_{n+1}$, $t_i=1/w_i$:

$$\mathrm{Pf}=\frac{(n-1)!!}{2^{n/2}}\frac1{\prod_iw_i},\quad U=\frac1{\prod_iw_i}\ \Longrightarrow\ \rho_n=\frac{(n-1)!!}{2^{n/2}}\prod_{i=1}^{n+1}w_i^{-1/2}\,dw_1\cdots dw_n ,$$

```
A5 plain banana n=2, 4, 6: True ; ((n-1)!!/2^{n/2}) int prod w^{-1/2} = pi^1, pi^2, pi^3
A5 banana(2,2,3,3,3) in (w, x) coordinates: Pf == (3/4)/prod w dw_1..dw_4 (no dx), U == 1/prod w: True
```

$K$ part. Closed form $F=\frac{(n-1)!!}{2^{n/2}}D_\Gamma\prod_P\frac{\Gamma(H+\frac32-a_P)}{(a_P-1)!}\big/\Gamma(H+\frac{n+1}2)$, $D_\Gamma$ as in b21.py (rows $(i,a)$, $i$ outer; $+\Lambda_P$ for $P=a$, $-\Lambda_P$ for $P=n+1$), graph `banana(lens)` of b21.py ($p=0$, $q=$ root, chains in order, edges $p\to q$). Exact $K$ (engine density, pointwise check, own Dirichlet), 2 source sets each (`v_a5_theta.out`, `v_a5_order.out`, `v_b2.out`):

```
H=2, n=4: banana(1,3,3,3,3), banana(2,2,3,3,3): K/F = 1 at 2 printed + 2 fresh sets each
H=2, n=4: banana(3,3,1,3,3), (3,3,3,3,1): K/F = 1, 1
H=2, n=2: theta(1,3,3), (2,2,3): K/F = 1
H=3, n=2: theta(2,3,4), (3,3,3): K/F = 1 ;  theta(1,4,4): K/F = -1, -1
H=1, n=4: banana(2,2,2,2,1), (1,2,2,2,2): K/F = 1 ;  banana(2,1,2,2,2): K/F = -1, -1
H=1, n=2: theta(2,1,2): K/F = 1 ;  theta(1,2,2): K/F = -1, -1 ;  theta(2,2,1): K/F = -1
```

Discrepancy, $H=1$, theta(1,2,2) (= `htk.theta(1,2,2)`): $D_\Gamma=\lambda_1\lambda_2$, $F=+\frac\pi3\lambda_1\lambda_2$, exact $K=-\frac\pi3\lambda_1\lambda_2$. Own-coordinate QMC without engine input (`v_a5_theta_qmc.out`): $\lambda=(2,3,-5)$: $K=15.70796327\pm2\cdot10^{-13}$ $(=-\frac\pi3\lambda_1\lambda_2)$, $F=-15.70796327$; $\lambda=(-1,4,7)$: $K=-29.32153143\pm4\cdot10^{-13}$, $F=+29.32153143$. The exact $K$ equals the CONVENTIONS theta formula $s\frac{D_\theta}2\prod_P\Gamma(H+\frac32-n_P)/(\Gamma(H+\frac32)\prod(n_P-1)!)$, $s=(-1)^{\lceil H/2\rceil+(H+1)(q+1)}$, in all six thetas of $H=1,2,3$ tested; $F$ differs from it by $-1$ for theta(1,2,2) ($H=1$) and theta(1,4,4) ($H=3$). Hence "at $n=2$ this is the theta formula of CONVENTIONS.md" is false, and $F$ without a sign factor is not the banana kernel for $H$ odd ($D=H+1$ even, where the kernel depends on the edge order through $\mathrm{sgn}(\pi)^{D+1}$). All cases with $H=2$ agree.

$n=6$, $H=1$, own QMC without engine input (`v_a5_n6_qmc.out`, $2^{10}\times8$ points, Dirichlet($e_b+1$)), $F=D_\Gamma\pi^3/7$:

```
banana(2,2,2,2,2,2,1) lambda [2,-1,3,1,-2,1,2]: K/pi^3 = -1.71428570 +- 8.3e-09 ; F/pi^3 = +1.71428571 (= 12/7)
banana(2,2,2,2,2,2,1) lambda [1,3,-2,2,1,-1,3]: K/pi^3 = -5.14285711 +- 2.5e-08 ; F/pi^3 = +5.14285714 (= 36/7)
explorer n=6 graph 1 (poles 5, 6; root 7 interior), lambda [3,1,-3,1,2,-3,1]: K/pi^3 = -5.14285711 +- 2.5e-08 ; |D_Gamma|/7 = 36/7
```

$|K|=|F|$ to 8 digits; $K=-F$ in the b21.py construction. For graph 1 the explorer's $-5.118\pm0.028$ is consistent with $-36/7$; $|D_\Gamma|=36$ only with the root source $-\sum_v\lambda_v=-2$ among the interior sources.

### a6

Polynomial-in-$\hat J$ test: $I=\sum_pc_pJ^p+\sum_pp\,c_pJ^{p-1}\frac{dt_m}{2t_m}\Theta$ with $c_p$ the $\Theta$-free coefficients (`v_rec.py`, `v_rec_a6.out`):

```
A6 restricted to G0..G4 (both plans): 32 merges, 19 polynomial in Jhat, rule holds: True
A6 all 57 merges (with thetas, (2,1) bananas, nested): 34 polynomial ; exceptions to the rule:
   nested(1,3,3) small, move 2: chains [[7,8,9],[14]] (plain, single reduced edge): not polynomial
   nested(1,3,3) small, move 3: chains [[10,11,12],[15]]: not polynomial
   nested(2,2,3) small, move 3: chains [[10,11,12],[15]]: not polynomial
   (the single reduced edges 14, 15 carry insertions that are themselves not polynomial in their Jhat)
```

The explorer's own `t_hyp.out` lists 32 merges with 19 `True`, not "20 of 30". "In every merge with a chain holding a reduced edge together with another edge, it is not" holds in all 57.
Component counts (`v_a6comp.out`, engine plan coordinates, own Pf and Jacobian):

```
G3 (fan) small plan: #nonzero components of rho_4 = 4 ; G4 small plan: 10   (large plans: 6, 6)
```

### a7

`v_a.py` A7, $C_{\rm rec}$ of the small plan (rows `[[0,1,0,-1,-1,0,1,0,1],[0,0,0,0,0,1,-1,1,0],[0,0,1,-1,0,0,0,0,1],[1,-1,0,0,1,0,0,0,0]]`), all 9 $dt_e$, no gauge:

$$\mathrm{Pf}(C_{\rm rec})=-\prod_{m=1}^3\frac{t_m\,dw_{m,1}}{w_{m,1}w_{m,2}}\wedge\frac1{4T}\,d\log\frac{t_9t_{10}t_{11}}{T^3},\qquad U=\prod_{m=1}^4\frac{t_m}{w_{m,1}w_{m,2}} :$$

```
A7: True at 3 rational points ; U: True
A7 second form -(1/2) prod (w w')^{-1/2} dw_m1 ^ d sqrt(Pi) vs rho at t_4 = 1: max rel. deviation 1.9e-16
```

Algebra between the two forms: $\frac1{4T}d\log\Pi=\frac{d\sqrt\Pi}{2T\sqrt\Pi}$, $\prod_mt_m^{1/2}/\sqrt\Pi=T^{3/2}$, $T\,w_{4,1}w_{4,2}=t_4$.

### a8

Proof steps re-checked: entries of $A$ are even, so $\mathrm{Pf}(A)^2=\det A$; $A=-\frac12B^2M$, $B=dM\,M^{-1}$; $\mathrm{tr}B^{2j}=-\mathrm{tr}B^{2j}$ (moving one odd entry past $2j-1$ odd entries); Newton identities over the commutative $\mathbb Q$-algebra of even forms give $\det B^2=0$; $\rho\wedge\rho=\det M^{-1}\det A=0$.

```
A8 Sym_4 (10 generators): #components of Pf = 141, of Pf^Pf = 0
A8 Sym_6, random 12-dim linear family (top form): 924, 0
A8 d=2 n=6 SP graphs 0 and 6, slice (12 generators, top form): 176, 0 ; 188, 0
A8 nested theta(1,3,3), (2,1) n=4 (12 generators): 145, 0
A8 G2: #components of tr B^2, tr B^4 = [0, 0]
```

### b1

V1 (`v_b1.out`, `v_b1_poly.out`): exact $K$ = claimed formula

```
G0..G4: 3 source vectors in [-7,7]^5 each (pointwise check True at 4 points each) and
        20 source vectors in [-100,100]^5 each: 20/20 for every graph
```

(a nonzero quartic vanishes at a uniform point of $\{-100..100\}^5$ with probability $\le4/201$.)
V2, own coordinates, own float density, no engine input:

```
G0 [-6,-2,2,2,4]: 94.748202 +- 5.4e-12, claim 94.74820225 (rel 1.6e-13) ; [2,2,-3,-3,1]: 71.061152, claim 71.06115169 (rel 1.7e-13)
G1 [1,2,6,7,1]: -145.41168 +- 0.00082, claim -145.4121715 (rel 3.4e-6) ; [6,2,-3,4,6]: -1677.8291 +- 0.0077, claim -1677.832748 (rel 2.2e-6)
G2 [-5,-5,4,7,-5]: 0.062 +- 0.057 ; [-4,-6,3,5,-3]: 0.024 +- 0.035 ; claim 0
G3 [-6,-6,-7,-5,2]: -9493.791 +- 0.024 (2^15 x 16); rerun 2^16 x 16, 3 seeds: -9493.7481 +- 0.020, -9493.7583 +- 0.024, -9493.7324 +- 0.014 ; claim -9493.736967 ; [6,5,-2,-7,4]: -237.9708 +- 0.0065, claim -237.9671283 (rel 1.5e-5)
G4 [2,-5,5,2,-7]: 32.0971 +- 0.0078 (2^15 x 16); rerun 2^16 x 16, 3 seeds: 32.07180 +- 0.0037, 32.07455 +- 0.0033, 32.08001 +- 0.0034 ;
   claim 32.0762143 ; [-6,6,1,6,-3]: 2158.165 +- 0.017, claim 2158.153496 (rel 5.3e-6)
```

Agreeing digits: 12 (G0), 5-6 (G1, G3), 4-5 (G4).

### b2

W4: $|I|=3n+1$, $|V|=2n+2$ and $2|I_S|\le3|V_S|-4$ for every vertex subset: True for the four graphs (own enumeration of subsets).
W1 (`v_b2.out`), printed source sets, exact $K$ = claim, pointwise check True:

```
banana(1,3,3,3,3): 1044576/5 pi^2, 77064/7 pi^2
banana(2,2,3,3,3): 19327113/280 pi^2, -3136 pi^2
nested theta(1,3,3): 139072 pi^2, -58327/200 pi^2
nested theta(2,2,3): -353210/3 pi^2, -587279/300 pi^2
```

W3 own QMC (first printed set): banana(1,3,3,3,3) 2061910.377 (rel 1.6e-13); banana(2,2,3,3,3) 681253.4269 (rel 2.7e-14); nested theta(1,3,3) 1372639 +- 31, claim 1372585.623 (rel 3.9e-5); nested theta(2,2,3) -1162030 +- 32, claim -1162014.324 (rel 1.3e-5).
Banana closed forms $K=\frac3{560}\pi^2D_\Gamma$, $K=\frac{\pi^2}{280}D_\Gamma$: $K/F=1$ at the 4 printed and 4 fresh source sets.

### c1

`v_c.py` C1: $U(t(\xi))=\prod_bb^{e_b}$, $e(x\text{ of child }c\text{ of an S-node})=L(c)$, $e(w\text{ of child }c\text{ of a P-node})=-1-L(c)$, own cycle basis, exact, 2 points, 5 decompositions (1 deterministic, 2 random flattened, 2 random non-flattened):

```
(D,n) = (2,2): 1 graph ; (2,4): 5 ; (2,6): 39 ; (3,2): 2 ; (3,3): 10 ; (4,2): 3 : True in every case
```

### c2

$\Gamma(k+\frac12)\in\mathbb Q\sqrt\pi$; a block with $h$ half-integer exponents gives $\pi^{h/2}$ from the numerator and $\sqrt\pi$ in the denominator iff $h$ odd: $\pi^{\lfloor h/2\rfloor}$, the same for every term of the Laurent polynomial. $g=U^{-H-1/2}\times$(rational in $t$ with denominators powers of $U$) $\times|\det[t;\partial t]|$; in ratio coordinates $t_e$ and $U$ are monomials (c1), so $g\,U^{1/2}$ is a Laurent polynomial in the ratios.

```
C2 G0..G4, banana(1,3,3,3,3), banana(2,2,3,3,3), nested(1,3,3), nested(2,2,3), own decompositions:
   probed exponent == -Uexp/2 mod 1 for all ratios: True ; min exponent -1/2 ; p = 2 (E=0), G2: p = 4 (E=1)
```

Every exact $K$ above (b1, b2, thetas of a5) is $\pi^{n/2}\times\mathbb Q[\lambda]$.

### c3

$\sum_j\lceil L_{ij}/2\rceil=\lceil L_i/2\rceil+\lfloor h_i/2\rfloor$ (since $h_i\equiv L_i$ mod 2); with $X=\sum_iL_i+k$, $X\equiv h_w$ (mod 2), $\sum_i\lceil L_i/2\rceil=(X-h_w)/2$, $\lceil L_{\rm new}/2\rceil=\lceil(X-1)/2\rceil$:

$$\sum_i\lceil L_i/2\rceil+\lfloor h_w/2\rfloor-\lceil L_{\rm new}/2\rceil=0\quad(h_w\text{ even or odd}),\qquad \phi_{\rm new}=\sum\phi_{ij}+2\sum_i\lfloor h_i/2\rfloor .$$

Valid for any decomposition tree (also nested P in P, S in S). `v_c.py` C3, 25 decompositions per graph:

```
p == ceil(n/2) + 2E: True for all graphs of the six (D,n)
E distribution: (2,4) {0: 4, 1: 1} ; (2,6) {0: 33, 1: 5, 2: 1} ; (2,2), (3,2), (4,2): all 0 ; (3,3): min E = 0 for all 10
decomposition-dependent E: none for even n ; (3,3): graphs 1,2,3,4,6,7,8,9 take E in {0,1}
```

(The explorer reports plan dependence at (3,3) for 2 graphs with its 2 plans; with 25 decompositions it is 8 of 10. $n$ odd, $K=0$.)

### c4

$K\in\pi^{n/2+2E}\mathbb Q[\lambda]$ and $K\in\pi^{n/2}\mathbb Q[\lambda]$ with $E\ge1$ give $\pi^{2E}P=Q$ coefficientwise with $P,Q$ rational, hence $K=0$. Computed cases: G0, G1, G3, G4 and the four (2,1) graphs $E=0$, $K\in\pi^2\mathbb Q[\lambda]$ (b1, b2); G2 $E=1$, $K=0$ (b1 V1, d1).

### d1

Hand-built decomposition: ear $A=P((0,1),S((0,2),(1,2)))$, $B=P((0,4),S((0,5),(4,5)))$, $C=P((1,4),S((1,3),(3,4)))$, centre $P(A,S(B,C))$; $X$ = ratio of $C$ in $S(B,C)$, $W$ = ratio of $S(B,C)$ in the centre; $l=y_0+y_1$ (current from vertex 0 into ear $A$).

Z1 (exact, forward-mode derivatives over Fractions, unit sources at each vertex, so the statement holds as a quadratic form in the sources):

$$\partial_X y_e\,\partial_Wl-\partial_Wy_e\,\partial_Xl=0\ \ (e=0..8),\qquad \Pi=\frac{t_At_Bt_C}{(t_A+t_B+t_C)^3}=X(1-X)\,W(1-W)^2 ,$$

```
Z1 3 ear points x 3 central points: True ; dl != 0: True ; edges with dy_e/dW != 0: [0..8]
```

Z2, central integral at fixed ear point, $X=\sin^2\phi$, $W=\sin^2\chi$, $80\times80$ Gauss-Legendre, $\lambda=(3,-2,5,1,-4)$, ratio $\int g F/\int|gF|$:

```
ear point 0: z=0: -1.1e-11/7.5 ; z random: -2.6e-11/13.1 ; z random, F=y0^2 y3: -6.5e-11/17.2 ; control F=X: 1.52/2.62
ear point 1: -3.6e-13/2.0 ; -1.5e-11/4.0 ; -1.5e-10/22.2 ; control 0.51/0.76
ear point 2: -9.6e-9/211 ; -3.0e-8/907 ; -1.1e-7/2.7e4 ; control 48.8/76.7
ear point 3: -1.2e-11/6.5 ; -1.8e-11/8.8 ; -2.2e-11/6.4 ; control 0.64/2.71
control convergence (ear point 3): nq = 40, 80, 120: 0.6369498329843, ...794, ...731
```

Z3, full 8-dim own QMC: $K(\lambda;0)=-0.006\pm0.067$, $K(\lambda;z)=0.025\pm0.065$ (random $z$), $\int|g|=38.89\pm0.07$.
Stokes argument re-checked: $|y_e|\le\sum_v|\lambda_v|$ (acyclic flow decomposition), so Fubini applies to $g\,e^{-y\cdot z}F(y)$; $y_e$ polynomial in the ratios, $\sqrt\Pi\,F'(l)\,dl$ continuous on the closed square and $\sqrt\Pi=0$ on its boundary.

## 3. Logical gaps in the explorer's arguments

- a1: none. The sign convention $\int d^n\theta\,\theta_1\cdots\theta_n=1$ is unambiguous because the top terms carry $n$ (even) $dt$'s.
- a2: the second sentence ("$x$ enters only through insertions of reduced edges") is a description of the recursion, not a separately stated or proved property.
- a3: the report gives the rule and exact checks, no derivation. The derivation (Fubini over the loop variables of each move, the unimodular shift $j_a=a_{(m,a)}-w_aJ$, the constrained Gaussian) is not written; the sign bookkeeping for consumed odd insertions is fixed only by the checks. Independence of $C_{\rm rec}$ from the routing of reduced edges (triangular unimodular change, det 1 with rows by decreasing $m$) is not stated.
- a4: full $J,\Theta$ dependence not proved (stated by the explorer); still only checked here.
- a5: no derivation of the $K$ formula is given (the fibre-determinant step is only named); the sign of $D_\Gamma$ is fixed only by code; the formula fails by a sign for $H$ odd (section 2). The $n=6$ evidence is "up to sign" and $36/7$ requires the root (an interior vertex of the banana) to carry $-\sum\lambda$; this is not stated.
- a6: the counts do not match the explorer's own output; the rule is an observation that fails outside the tested set.
- a7: intermediate formula for the central Berezin integral not checked separately (the final identity is).
- c1: the factorization $M=S^{\mathsf T}\mathrm{diag}(M_{\rm red},M_{\rm loc})S$, $\det S=1$, is asserted, not shown.
- c2: (i) the hypothesis that every term has exponents $>-1$ is unproved; (ii) "the recursion a3 and the polynomial currents give no other denominators" is an argument via a3, whose sign rules are only checked; the direct argument (denominators of $y$, $M^{-1}$ are powers of $U$, $U$ and $t$ monomial in the ratios) does not need a3.
- c3: complete. Plan independence of $E$ for even $n$ is observed (here on 25 decompositions per graph), not proved; the headline "Theorem: $K\in\pi^{n/2+2E}\mathbb Q[\lambda]$" inherits the hypothesis of c2.
- c4: none beyond c2.
- d1: "every current is affine in $l$": only functional dependence on $l$ at fixed ears is needed and is what is verified here; affinity was not checked. Absolute integrability (Fubini) and continuity of $\sqrt\Pi\,F'(l)\,dl$ up to the boundary are used without being stated.

## 4. Not done

- a5 at $n=6$: no exact value; the engine density of the $n=6$ banana at numeric sources (`v_a5_n6.py`) ran about 9 CPU-minutes without output and was stopped; only the own QMC above exists.
- c5: not examined (status open).
- a4, a6: no proofs attempted; only independent checks.
- No symbolic-source exact kernels were computed independently; b1 polynomials were confirmed at 23 integer source vectors per graph.
