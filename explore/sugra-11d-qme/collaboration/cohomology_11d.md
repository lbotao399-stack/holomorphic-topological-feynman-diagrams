# Local cohomology for the order-$g$ two-loop obstruction (Claude, 2026-11-21)

Tags: [drv] derived here; [chk] exact computation (scripts in `../cohomology/`, pending where marked).

## 1. Strict model of the tree differential [drv; chk pending: T1 of `../cohomology/REPORT.txt`]

R-Q0 in the variables $X^i=f(\nu)\mu^i=\mu^i/(1-\nu)$, i.e. $\mu^i=(1-\nu)X^i$:
$$
QX^i=-X^j\partial_jX^i,\qquad Q\nu=(1-\nu)\,\partial_iX^i-X^i\partial_i\nu,\qquad
Q\gamma_i=\partial_i\beta+X^k(\partial_i\gamma_k-\partial_k\gamma_i),\qquad Q\beta=X^iX^j\partial_i\gamma_j .
$$
Derivation: $Q\mu^i=-\partial_j\big((1-\nu)X^jX^i\big)$ and $Q\mu^i=-(Q\nu)X^i+(1-\nu)QX^i$ give $(1-\nu)QX^i=-(1-\nu)X^j\partial_jX^i$ (sign fixed by $Q^2=0$; T1 of `../cohomology/REPORT.txt`). Also $f\mu=X$ and $f^2\mu\mu=XX$. The factor $f(\nu)$ disappears at $g=0$. The $g$-term becomes $QX^i\supset\frac g2 f\,\varepsilon^{ijklm}\partial_j\gamma_k\partial_l\gamma_m$.

$\nu$ enters only $Q\nu$. $\rho=1-\nu$ is a density: $Q\rho=-\partial_i(\rho X^i)$.

## 2. Gradings [drv]

- $w=\#\gamma+\#\beta$: preserved by $Q|_{g=0}$; the $g$-term raises it by $2$. Sector $m$ has $w=2m-1$.
- $G=\sum|\alpha|-\#X-\#\beta$: preserved, including the $g$-term.
- Torus weight $t$ ($X^i:+e_i$, $\gamma_i:-e_i$, $\partial_k:-e_k$): preserved. At $g=0$, $Q$ is $GL_5$-equivariant.

Order $g^1$: $O=gX_1+O(g^2)$ and $I_2=I_2^{(0)}+gI_2^{(1)}+\dots$, with $I_2^{(0)}$ in $w=-1$, hence $=0$. So
$$[X_1]=0\iff X_1\in Q|_{g=0}\big(C_{w=1}\big).$$
$X_1$ slot: $w=1$, $G=23-9=14$, $t=(-3,-3,-3,-3,-3)$ (upper $9$, lower $24$), odd.

## 3. Constant modes and divergence [drv; chk pending: T2, T3]

Linear part per letter: $\nu_\alpha\mapsto\partial^\alpha(\partial\cdot X)$ and $\gamma_{i,\alpha}\mapsto\partial^\alpha\partial_i\beta$. Single-letter cohomology:
$$X\ \mathrm{mod}\ p\cdot X\ \ \oplus\ \ \{\text{functions of }p\wedge\gamma\}\ \ \oplus\ \ \mathbb C\,\beta_\varnothing .$$

- $\nu$ is eliminated together with $\partial^\alpha(\partial\cdot X)$. The ideal $(\nu_\alpha,\partial^\alpha\partial\cdot X)$ is a dg ideal, because $Q(\partial\cdot X)=X\cdot\partial(\partial\cdot X)$ (using $\mathrm{tr}((\partial X)^2)=0$ for odd $X$).
- Quotient: $C_{\rm quot}=C^\bullet\big(S_5^{\rm formal};M\big)$ with $M$ the $(\beta,\gamma)$ complex $\Omega^0\xrightarrow{\partial}\Omega^1$, so $H(M)=\Omega^2_{cl}\oplus\mathbb C\,1_\beta$.
- Constant sector retained: $\beta_\varnothing$ ($w=1$, $G=-1$, $t=0$) is closed for the linear part, but $Q\beta_\varnothing=X^i_\varnothing X^j_\varnothing\gamma_{j,e_i}\neq0$. Its class is to be read off from the T5 tables.

Split of $w=1$ along $\beta_\varnothing$. Write a $w=1$ cochain as $\beta_\varnothing\,c(X)+b(X;\gamma)$. Then $Q(\beta_\varnothing c)=\kappa\,c\pm\beta_\varnothing\,Qc$, with
$$\kappa:=Q\beta_\varnothing=X^i_\varnothing X^j_\varnothing\gamma_{j,e_i}\in C^2\big(S_5;(\Omega^2_{cl})^\vee\big),\qquad\kappa(x,y;F)=F(x,y)|_0,\qquad G(\kappa)=-1 .$$
So $C_{\rm quot,w=1}=\mathrm{Cone}\big(\kappa\cdot:C^\bullet(S_5)\to C^{\bullet+2}(S_5;(\Omega^2_{cl})^\vee)\big)$, and
$$\cdots\to H^q(S_5)_{G+1}\xrightarrow{\ \kappa\ }H^{q+2}\big(S_5;(\Omega^2_{cl})^\vee\big)_G\to H_{w=1,G}\to H^{q+1}(S_5)_{G+1}\xrightarrow{\ \kappa\ }\cdots$$
At the $X_1$ slot the inputs are $H(S_5)$ at $G=15$ and $H\big(S_5;(\Omega^2_{cl})^\vee\big)$ at $G=14$, both at $t=(-3)^5$, shifted by the translation factor of §4.

Sources [src, web extraction]:
- RSW 2111.03049, Thm 2.1(2) and eq. (33): the flat-space symmetry algebra is the non-trivial $L_\infty$ central extension of $E(5,10)$ by the even cocycle $\varphi(\mu,\mu',\alpha)=\langle\mu\wedge\mu',\alpha\rangle|_{z=0}$. This is $\kappa$, so the constant-$\beta$ sector does not split at the $L_\infty$ level.
- Hahner–Paquette–Raghavendran 2506.02246, §4.3.2: the all-order obstruction is "controlled by the twelfth Gelfand–Fuchs cohomology of $E(5|10)$ – a group which remains unknown". The $X_1$ slot is the order-$g$, $w=1$, $G=14$ piece of that group.

## 4. Translations [drv; chk pending: T4]

With $x=e_k$ a constant vector field, $\iota_k:X^i_\alpha\mapsto\delta^i_k\delta_{\alpha0}$, and $h_k:\beta_\alpha\mapsto-\gamma_{k,\alpha}$ (T4: $D_k=-[Q,\iota_k+h_k]$):
$$[Q,\iota_k]=D_k-\Delta_k,\qquad \Delta_k\gamma_i=\pm\partial_i\gamma_k=\pm[Q,h_k]\gamma_i,\qquad\Rightarrow\qquad D_k=[Q,\iota_k+h_k],$$
$$[\iota_k+h_k,\iota_l+h_l]=0,\qquad[\iota_k+h_k,D_l]=0 .$$
Conjugating the Koszul complex of $V=C_{\rm quot}/\sum_kD_k$ by $\exp\sum_k(\iota_k+h_k)\otimes\iota_{e_k^*}$, using freeness over $\mathbb C[D]$ in positive arity ($w=1$ has arity $\ge1$):
$$H(V)_{w=1}\ \cong\ \bigoplus_{K\subseteq\{1..5\}}H(C_{\rm quot})_{w=1,\ G-|K|,\ t+\sum_{k\in K}e_k}\ [\text{parity}+|K|] .$$
$SL_5$-invariant form: $\big[\bigoplus_jH(C_{\rm quot})_{G-j}\otimes\Lambda^j\mathbb C^5\big]^{SL_5}$, at $\det^{-3}$.

Non-holomorphic directions: as in flagship Thm 3.1, translation invariance alone gives $\bigoplus_{S\subseteq\{\tau,\bar z_1,\dots,\bar z_5\}}$, which is $64$ summands. The $\eta$-invariant part is $S=U$. Equivariance of the 11d scheme: $[P]_\Omega\propto(\tfrac\tau2d^5\bar z-d\tau\wedge\vartheta)/r^{11}$ is a form in $d(x-y)$, and the argument of flagship Prop. 3.6 carries over [drv]. All 64 summands are kept until that is checked.

## 5. Abelian probes cannot detect $X_1$ [drv]

Configurations $\mu_a=v_ae^{p_a\cdot z}$, $\gamma=c\,e^{p\cdot z}$ with $v_a\cdot p_b=0$, $v_a\cdot c=0$, $v_a\cdot p=0$ satisfy $Q\phi=0$ identically in the formal parameters. Every coboundary therefore vanishes on them.

By the first fundamental theorem for $SL_5$, each monomial of $X_1$ is a product of pairings $v\cdot p$, $v\cdot c$ and brackets $[\,5\text{ covectors}\,]$, $[\,5\ v\text{'s}\,]$. Since $9\not\equiv0\pmod5$, every monomial contains a pairing. So $X_1$ vanishes on all such probes, and the probe route is closed.

## 6. Size of the $X_1$ slot [drv; exact counts pending: T7]

$\dim S_5^{(k)}=5\binom{5+k}{k+1}-\binom{4+k}{k}$:

| $k$ | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| $\dim S_5^{(k)}$ | 70 | 160 | 315 | 560 | 924 | 1440 |

Nine vector arguments with $\sum k=13-k_F$ ($k_F\ge0$ the Taylor degree of $\partial\gamma$) span spaces of order $10^{10}$ before $SL_5$ reduction. Direct rank computation at $G=14$ is not feasible.

Order of work:
1. Exact tables for $n=2,3$ (T5, T6), to see where $H_{w=1}$ lives.
2. A structural computation of $H\big(S_5;M^\vee\big)$ at $G=14-j$, via Hochschild–Serre relative to $\mathfrak{sl}_5$ and the Koszul resolution of $\partial\cdot X=0$ over $W_5$, where the Euler field is inner.
3. Codex's antisymmetrized $X_1$ on generic fixtures: if $X_1\equiv0$, the cohomology question is void.

## 7. Exact results so far ([chk], `../cohomology/REPORT.txt`; $n$ = holomorphic dimension)

- **T1.** $Q^2=0$ with $QX=-X\partial X$ ($n=1,2,3$, $|\alpha|\le3$; $n=5$, $|\alpha|\le1$). Substitution $\mu=(1-\nu)X$ reproduces R-Q0. The ideal $(\partial^\alpha\partial\cdot X)$ is a dg ideal.
- **T3.** $H(C_{\rm full}^{\le A})=H(C_{\rm quot}^{\le A})$ for every truncation $A$ ($n=2$, $w=0,1$). This is the $\nu$-elimination.
- **T4.** $D_k=-[Q,\iota_k+h_k]$ with $h_k\beta=-\gamma_k$. The $H(V)$ formula of §4 holds on $494$ slots ($n=2$, $w=1$), with explicit representatives $\eta_Kc$.
- **T4, $\mathfrak{sl}_n$ Cartan.** $[Q,\iota_\xi+h_\xi]=E_{ba}$ for $\xi=z^b\partial_a$, and $[Q,\iota_c+h_c]=(c\cdot t)\,\mathrm{Id}$ for $\sum c=0$. Hence $H(C_{\rm quot})$ is supported on $t=(k,\dots,k)$ and is $SL_n$-trivial. For $n=5$:
$$H(V)_{(-3)^5,N}=H(C_{\rm quot})_{(-3)^5,N}\oplus H(C_{\rm quot})_{(-2)^5,N+5}\qquad(G=14\ \text{and}\ G=9).$$
- **T5, T6, $n=2$, $w=1$.** $H(C_{\rm quot})\neq0$ only at $G=-3,-1,1,7$ (for $G\le7$), each in two parities. Every class has a $\beta_\varnothing$-component $\beta_\varnothing c$ with $c$ a nonzero class of $H(C^{w=0})$, or is $c\,F_0$ with $F_0=\gamma_{2,e_1}-\gamma_{1,e_2}$.
  - $G=7$: $\beta_\varnothing c$ with $c\in H^7(C^{w=0})_{G=8,(-4,-4)}$. For $n=2$, $S_2$ is the Hamiltonian algebra and $c$ sits in degree 7, weight 8: the place of the Gelfand–Kalinin–Fuks class [identification drv].
  - The constant-$\beta$ sector carries classes, so it cannot be dropped.
- **T5, $n=3$, $w=1$**, $H(C_{\rm quot})$ for $G=-4..-1$; slots with $t\ne(k,k,k)$ vanish:
  - $G=-4$, $t=(1,1,1)$: $N=4,7,9,12$, one class each, $\beta_\varnothing X^1_0X^2_0X^3_0\cdot(1,x_3,x_5,x_3x_5)$ with $x_3,x_5$ the $\mathfrak{sl}_3$ cocycles; $c(z)$ is not exact in each case.
  - $G=-1$, $t=0$: $N=5,7,8,10,12,13,15$ with $\dim=1,1,1,2,1,1,1$. $c(z)=0$ at $N=5$ and for one class at $N=10$; $c(z)$ is exact at $N=7$; $c(z)$ is not exact at $N=8,12,13,15$ and for the other class at $N=10$.
  - $G=-3,-2$: $0$. $G=0,1$: $0$ by the Cartan support. $G=2$: not computed.
  - $H(C^{w=0})$: at $t=(1,1,1)$, $G=-3$, classes at $N=3,6,8,11$; at $t=0$, $G=0$, classes at $N=0,7,9,11,12,14$.
  - So the $\beta$-type classes at $G=-1$ have $N-1\in\{7,9,11,12,14\}=\{N>0:H(C^{w=0})_{t=0,G=0,N}\ne0\}$, and $\beta_\varnothing\cdot1$ is removed by $\kappa$.
  - $n=3$, $H(V)$: not computed.
- **T7, $n=5$ target**, $SL_5$-invariant cochain dimensions ($\det^{-3}$, $G=14$):

| $N$ | $C_{\rm quot}$ | $V$ |
|---|---|---|
| 8 | 20,799,946 | 1,313,542 |
| 9 | 182,332,617 | 14,656,912 |
| 10 | 1,290,504,415 | 125,272,048 |

So:
$$[X_1]\ \text{is decided by}\ H(C_{\rm quot})_{w=1}\ \text{at}\ (G,t)=(14,(-3)^5),\,(9,(-2)^5),$$
which on the $n=2,3$ evidence is built from $H^\bullet(S_5)$ (trivial coefficients) at $\det^{-3},\det^{-2}$ through $\beta_\varnothing$ and $\kappa$, and from $F$-type products. Direct ranks at $10^7$–$10^8$ are beyond the present code.
