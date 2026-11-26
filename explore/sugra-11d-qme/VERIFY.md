# VERIFY — sugra-11d-qme

Scripts and logs: `verify/` (`vcommon.py`, `v_structure.py`, `v_flows.py`, `v_density.py`, `v_integrand.py`, `v_degree.py`, `v_timing.py`, `v_jets.py`; each `*.py` writes the `*.log` of the same name). Exact `Fraction` arithmetic throughout. None of the verify scripts imports `htk` or the explorer's code, except `v_integrand.py`. That script calls `code/flows.py` (`patterns`, `flow_integrand`) and `htk.graph.theta` to get the value being tested. The evidence is that this value equals an evaluation by a different algorithm.

`P/explore/sugra-11d-qme/REPORT.md` does not exist (the folder contains only `code/`, `structure.log`, `checks.log` and `count_jets.log`). The gaps listed below come from the structured report, the code and the logs.

Graph: $\theta_{pqr}$ with $N=0$, interior vertices numbered along $A,B,C$, $N'=11$ (root), edges $N\to N'$ along the paths, edge order $(A,B,C)$ (S3 §sec-Q2). `vcommon.theta_graph` gives the same edge list as `htk.graph.theta` (assert in `v_integrand.py`).

---

## a1 — no graph at $g^0$

$$\#\{\text{incoming letters }\gamma,\beta\}=n_a+n_b+3m=12+2m\ \ge\ |I|=13\ \Rightarrow\ m\ge1 .$$

Method: transfer matrices along the three paths. The edge labels are (field at $e(0)$, field at $e(1)$) $\in\{(\mu,\gamma),(\gamma,\mu),(\nu,\beta),(\beta,\nu)\}$ (S3 eq-Delta). The vertex weight is the number of injective letter maps of the types $a,b_0,b_1,(b_2),c$ (S3 R-action, $f=\sum_k\nu^k$). This replaces the explorer's DFS over vertices.

```
(1, 6, 6) no Dtheta rule, kmax=1: total maps 198893259456 ; with m=0: 0 ; with n_a>0: 76826430144 ; min m = 1
(2, 5, 6) ... with m=0: 0 ...   (3, 4, 6) ... with m=0: 0 ...   (3, 5, 5) ... with m=0: 0 ...   (4, 4, 5) ... with m=0: 0
```

Even without the $D_\theta$ rule, no letter map exists at $m=0$. **confirmed**

## a2 — no vertex (a)

$$D_\theta=\det\begin{pmatrix}\Lambda_A&\Lambda_B&0\\0&\Lambda_B&\Lambda_C\end{pmatrix}:\ \text{column }v=(\lambda_v;0),(\lambda_v;\lambda_v),(0;\lambda_v)\ \Rightarrow\ \deg_{\lambda_v}D_\theta=1\ \ \forall v\ \text{interior};$$

$K_\theta[F]=\tfrac{D_\theta}2\int\cdots$ (S3 eq-insertion). An interior vertex with no output letter gives $\lambda_v\mapsto\partial_{z_v}$ acting on a constant, so the term is $0$. Type (a) has 2 letters, so it is neither interior with an external letter nor trivalent.

```
(1, 6, 6) Dtheta rule, kmax=1: maps 122066829312 ; m=0: 0 ; n_a>0: 0 ; ext beta>0: 0 ; m values [1, 2, 3, 4, 5]
(1, 6, 6) Dtheta rule, kmax=2: maps 4458050224128 ; m=0: 0 ; n_a>0: 0 ; ext beta>0: 0
```

The same holds for all five graphs, with kmax $=1$ and $2$. Without the rule, $n_a>0$ maps exist (above), so the rule alone removes (a). **confirmed**

## a3 — sectors

$$\#\mu_{\rm ext}=2(12-m)-13=11-2m,\qquad\#\gamma_{\rm ext}=(12-m)+3m-13=2m-1,\qquad\deg_\lambda=\underbrace{(12-m)+2m}_{\rm vertices}+\underbrace{10}_{D_\theta}=22+m .$$

```
all five graphs: keys violating (mu,ga)=(11-2m,2m-1): 0 ; max K by m: {1: 11, 2: 10, 3: 9, 4: 8, 5: 7}   (kmax=1)
   K=0 maps by m: {1: 491520, 2: 17694720, 3: 111476736, 4: 159252480, 5: 39813120}
   comparison with structure.log lines for this graph: 50 lines, all equal          (166 and 445)
v_degree.log: (pqr) integrand x D_theta homogeneous of degree 23 in lambda (2 flows): True     (all five, m=1)
```

The degree $23$ comes from the integrand, evaluated by the grouped formula of b1, at $\lambda$ and $3\lambda$. At $K=0$ every interior vertex has exactly one external letter, and $N,N'$ have none (a4 roles). **confirmed**. The explorer's evidence says "K up to 12" for $b_0,b_1$; the enumeration (its own log and mine) gives $K_{\max}=12-m$, i.e. $11$ at $m=1$.

## a4 — flows

Trivalent $(\mathrm{in},\mathrm{out})\in\{(1,2),(3,0)\}$ has probability $\tfrac12$ per trivalent vertex, also for $p=1$. Hence $2^{13}/4=2^{11}$.

```
(pqr) admissible flows 2048 ; by (m, ext mu, ext ga): {(1, 9, 1): 40, (2, 7, 3): 480, (3, 5, 5): 1008, (4, 3, 7): 480, (5, 1, 9): 40} ; equals 4*binom(10,2m-1): True
    m=1: roles {"c@N'": 10, 'c@interior': 20, 'c@N': 10} ; letter maps per flow {12288: 40} ; a4 structure at m=1: True
```

The same holds for all five graphs. "Structure" checks two cases. With c interior on $P$, the edges of $P$ point into the sink and the other two paths are directed one $N\to N'$ and one $N'\to N$. With c at $N$ or $N'$, there is exactly one interior source and no sink. **confirmed**

## b1 — explicit form of the $g\,\mu^9\gamma$ sector

Independent evaluation of $W_\Gamma$. $D_e=\sum_i\partial_L/\partial\mu^i_{t(e)}\,\partial/\partial\gamma_{i,h(e)}$ is odd and acts on the single even factor $W_{t(e)}$. By the graded Leibniz rule,

$$D_{e_1}\cdots D_{e_{13}}\bigl(W_0\cdots W_{11}\bigr)=(-1)^{\#\{k<l:\ t(e_k)>t(e_l)\}}\prod_{v=0}^{11}\Bigl[\Bigl(\prod_{e:\,t(e)=v}^{\rm incr.}\frac{\partial_L}{\partial\mu_e}\Bigr)\Bigl(\prod_{e:\,h(e)=v}\frac{\partial}{\partial\gamma_e}\Bigr)W_v\Bigr],$$

with the external odd letters in increasing vertex order. Each bracket is computed by literal left derivatives on the ordered words $(\mu^a,\mu^b,\partial_a\gamma_b)$ and $\tfrac16\varepsilon^{ijklm}(\gamma_i,\partial_j\gamma_k,\partial_l\gamma_m)$. A derivative on the contracted letter inserts $\mp y_e$ at $e(0)$ and $e(1)$. The result is contracted along the three paths. There is no sum over letter maps and no per-map sign.

```
vertex factor pass-through = u^j J_i - delta^j_i (J.u) : True
vertex factor trivalent b0: dL/dmu_{i'} dL/dmu_{i''} = -(delta^j_{i''} J_{i'} - delta^j_{i'} J_{i''}) : True
vertex factor source: = -(lam_{i'} ga_{i''} - lam_{i''} ga_{i'}) : True
vertex factor c sink: (1/6) eps, J on contracted, lam on external (sum over 6 maps) : True
(2, 5, 6) c@interior flow 22 | grouped (mine) == explorer flow_integrand: True | nonzero: True | letter maps 12288 | sec mine 0.0, explorer 35.9
(2, 5, 6) c@N flow 32 | ... True | nonzero: True | ... explorer 39.6
(4, 4, 5) c@N' flow 3 | ... True | nonzero: True | ... explorer 38.7
(4, 4, 5) c@interior flow 29 | ... True | nonzero: True | ... explorer 40.7
(3, 4, 6) c@interior flow 5 | ... True | nonzero: True | ... explorer 42.0
v_timing.log: flows evaluated 200 ; nonzero 200 ; total sec 7.68 ; per flow 0.038 s
```

These are exact equalities at random rational $(x,w,\lambda,u,\gamma)$ with Kirchhoff currents from my own solver. The graphs $256$, $445$ and $346$ were not in the explorer's pointwise checks. Other factors of the formula:

- $|\mathrm{Aut}|=4,2,2,4,4$ (networkx; `v_flows.log`);
- $40$ flows per graph and $12288=2^{11}\cdot6$ maps per flow (a4);
- $\tfrac{D_\theta}2\prod_Pw_P^{11/2-n_P}$ with its sign (b2);
- the currents $y^P_j=w_P\sum_{P'}\hat m_{P'}-\hat m_P+Q^P_j$ at $\lambda_N=0$, a special case of the check in b2.

**confirmed** as an evaluation of S3's $\cK_\Gamma(\partial;z)\circ W_\Gamma$ (§sec-assembly, eq-insertion), with the reading of the bullet tensors in gap G4 and the conventions in G5 and G6 assumed as stated.

## b2 — reduced density

$$\Omega=U^{-5}\rho_2\wedge\prod_{i=1}^5\psi^i_1\psi^i_2=D(t)\,\iota_E\,dt_0\cdots dt_{12},\qquad D=\frac{(-1)^{e_0}}{t_{e_0}}[\Omega]_{\widehat{e_0}},\qquad \mathrm{Pf}(-\tfrac12dM\,M^{-1}dM)=A_{12},$$

$$D(t(\xi))\,\bigl|\det[t;\partial_\xi t]\bigr|=-\frac{D_\theta}2\prod_Pw_P^{\frac{11}2-n_P},\qquad t^P_j=\frac{x^P_j}{w_P},\quad \xi=(x^A_{\ge2},x^B_{\ge2},x^C_{\ge2},w_B,w_C).$$

The cycle basis is $A-B$, $B-C$; for $D=6$ the result does not depend on the basis. The exterior algebra, Jacobian, Kirchhoff solver and $D_\theta$ are my own code. $\lambda_N\neq0$.

```
(pqr) poly currents == own Kirchhoff (+KCL,KVL), 3 pts: True | Pf == eq-rho2: True | D(t) same for e0=0,12: True | D|J| = -(Dtheta/2) prod w^(11/2-n) (square & sign): True
```

All five graphs give this. The sign of S3 eq-Ktheta follows from $\varsigma_\theta$. I recomputed $\varsigma^{\rm ext}$ of S3 eq-extraction as the $78\times78$ determinant of $\bigwedge_e(\theta^\tau_e\theta^1_e\cdots\theta^5_e)$, $\theta^\tau_e=d\tau_e-w_edt_e$, $\theta^i_e=d\bar z^i_e-4y^i_edt_e$, divided by $4^{5n}\bigl(\prod_r\psi^\tau_r\prod_{i,r}\psi^i_r\bigr)_S$:

```
(pqr) ratios ['1', '1', '1', '1'] -> varsigma^ext = 1 , varsigma_theta = -1      (all five)
control C_7: varsigma^ext ratios ['-1', '-1', '-1', '-1']                        (S3: -1)
```

$$\varsigma_\theta=(-1)^n(-1)^{\binom62\binom{11}2}\varsigma^{\rm ext}=(+1)(-1)^{825}(+1)=-1,\qquad(-1)^{|I|(|I|-1)/2}=(-1)^{78}=+1=(-1)^n,$$

$$\cK^{\rm S3}_\theta=\varsigma_\theta\int\Omega=+\frac{D_\theta}2\int\prod_Pw_P^{\frac{11}2-n_P}\cdots\quad(\text{S3 eq-Ktheta}).$$

**confirmed**

## b3 — jet count

$$N=\sum_{b}\sum_{n_1+\dots+n_5=9}\bigl[t^{(3,..,3)+n-e_b}\bigr]\prod_aE_{n_a}(t)\prod_i\frac1{1-t_i},\qquad E_n=\sum_{\lambda\vdash n}\frac{\varepsilon_\lambda}{z_\lambda}\prod_{k\in\lambda}\prod_i\frac1{1-t_i^k}.$$

This organization by $n_a=\#\{\mu^a\text{-jets}\}$ differs from the explorer's (Newton for $e_9$ of the whole alphabet, with parts distributed over directions). It is checked against explicit enumeration:

```
n_mu=1 n_der=5 target=(1, 1, 1, 1, 1): brute 640, generating function 640, equal True
n_mu=2 n_der=6 target=(1, 1, 1, 1, 1): brute 23330, generating function 23330, equal True
n_mu=4 n_der=3 target=(0, 0, 0, 0, 0): brute 5870, generating function 5870, equal True
n_mu=3 n_der=2 target=(0, 0, 0, 0, 0): brute 660, generating function 660, equal True
n_mu=3 n_der=5 target=(1, 1, 1, 0, 0): brute 41370, ... True
n_mu=5 n_der=6 target=(1, 1, 0, 0, 0): brute 1086934, ... True
n_mu=5 n_der=4 target=(0, 0, 0, 0, 0): brute 45420, ... True
n_mu=9 n_der=6 target=(0, 0, 0, -1, -1): brute 356196, ... True
n_mu=9 n_der=7 target=(1, 0, 0, -1, -1): brute 3001889, ... True
n_mu=6 n_der=6 target=(2, 0, -1, 0, 0): brute 834970, ... True
mu^9 gamma, weight (3,3,3,3,3) [23 derivatives]: 63414968175297950 ; explorer 63414968175297950: True
```

The weight $(3,3,3,3,3)$ was checked on the integrand itself. Under $\lambda_i,y_i,\gamma_i\mapsto s_i(\cdot)$ and $\mu^i\mapsto\mu^i/s_i$, the integrand times $D_\theta$ scales by $(s_1\cdots s_5)^3$:

```
(pqr) torus weight (3,3,3,3,3) of integrand x D_theta (2 flows): True     (all five)
```

**confirmed**

## c1, c2

c1 has status open and c2 has status derived-unchecked; neither was verified. Measured fact: the flow integrand evaluates in $0.038$ s per flow per point (`v_timing.log`), against the explorer's figure of about 30 s.

---

## Logical gaps

- **G1** (a3 evidence). "K up to 12 for 166 and 445 with b_1 allowed" does not match `structure.log`, which gives $K_{\max}=12-m$ ($11$ at $m=1$).
- **G2** (b3 evidence). The small-case brute force (640, 23330, 5870, 660) is not in `code/count_jets.py`. `count_jets.log` prints "check small case ... skipped". The four numbers are reproduced here by enumeration.
- **G3** (b1 evidence). In `checks.log` the comparison "integrand from polynomial currents = from Kirchhoff currents" uses the same tensors and signs on both sides, so it tests only the currents. The edge-swap and vertex-order tests check only that $s_\pi$ is internally consistent. Neither tests the vertex tensors or $s_\pi$ against an independent evaluation of $W_\Gamma$. The grouped formula above provides that test on 5 flows.
- **G4** (b1 statement). The bullet tensors (pass-through $\mu^jJ_i-\delta^j_i J\cdot\mu$, trivalent $b_0$, source) are each already the sum over that vertex's letter maps, with the relative Grassmann sign of the two $\mu$ letters. Read literally, "$\sum_\pi s_\pi N_\pi$ with $N_\pi$ the contraction of these tensors" makes $N_\pi$ independent of $\pi$. The two maps at a pass-through vertex have opposite $s_\pi$, so the literal reading gives $0$. `flows.vertex_tensor` uses single-map tensors, and only that reading is the one confirmed. The displayed trivalent-$b_0$ and source tensors also carry the sign of the derivative order $\partial_{i''}\partial_{i'}$; with $\partial_{i'}\partial_{i''}$, $e'<e''$, they change sign (`v_integrand.log`).
- **G5** (b1). The sign between the forms $[P_e]$ and the odd contractions $\Delta_e$ is not fixed by S3. It is taken as uniform over the five graphs, since they have equal $|I|$ and $|V|$ (explorer's open issue). S3 §sec-assembly states $\cK_\Gamma\circ W_\Gamma$ directly, and b1 follows that statement.
- **G6** (b1). The normalization $12!\sum|\mathrm{Aut}|^{-1}$ reads S3 def-bracket literally, for 12 identical vertices. It also uses the invariance of $\cK_\Gamma W_\Gamma$ under the choice of which vertex sits at the origin, together with $\lambda_{N'}=0$ (the quotient by total derivatives). S3 does not normalize $\{\cI^{\otimes12}\}_0$ against $\Delta I+\frac12\{I,I\}$ (explorer's open issue).
- **G7** (b2). The explorer's sign statement relies on $\varsigma^{\rm ext}_\theta=+1$ taken from S3 and `htk`, without recomputing it. It is recomputed here.
- **G8** (c2). The claim that the order-$g$ obstruction is $Q_0$-closed is asserted from the vanishing at $g^0$. No BV-Jacobi argument within S3 is given.
