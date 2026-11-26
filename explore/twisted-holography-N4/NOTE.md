# Twisted holography: holomorphic twist of $\mathcal N=4$ $U(N)$ vs twisted IIB on $\mathbb C^5\setminus\mathbb C^2$

2026-11-21. Tags: [src] read in the source; [drv] derived here; [chk] exact computation (`code/single_particle.py`); [open].

Sources:
- [src] [LLC manuscript in this project](../../sources/LLC_main.tex):
  - Thm b-00088: all-loop $Q_n$, $Q_n:(\lambda)_Y\mapsto(\lambda-n\mathsf d)_{Y+1-n}$ with $\mathsf d=(1,1,1;-1,-1)$;
  - eq. (2-60) and Thm b-00630: $Q_nG=0$ for single gravitons $G$, the components of $\mathrm{Tr}(dC)^p$;
  - Thm b-00096: $Q_n\mathcal O=0$ for $n>\#(\mathsf b,\beta)$.
- [src] Costello–Li arXiv:1606.00365 §9: twisted IIB on $AdS_5\times S^5$ is BCOV on $X=\mathbb C^5\setminus\mathbb C^2$ with the $PV^{2,2}$ field equal to $NF$, $\bar\partial F=\delta_{\mathbb C^2}$; differential $\bar\partial+t\partial+N\{F,-\}$; the $Q$-cohomology of $\mathfrak{psl}(4|4)$ is $\mathfrak{psl}(3|3)$.
- [src] standard: Loday–Quillen–Tsygan (LQT), Hochschild–Kostant–Rosenberg (HKR).

## 1. Boundary: planar single-trace classes, all loops [drv]

$$\text{local operators}=C^\bullet(\mathfrak{gl}_N\otimes A),\qquad A=\mathbb C[[z_1,z_2]]\otimes\Lambda[\theta_1,\theta_2,\theta_3],\qquad Q_0=d_{CE}$$
$$\lim_{N\to\infty}H\big(C^\bullet(\mathfrak{gl}_N\otimes A)\big)=\mathrm{Sym}\big(HC^{\bullet-1}(A)\big)\qquad\text{(LQT)}$$
$$HC_n(A)=\Omega^n_A/d\Omega^{n-1}_A\oplus\mathbb C[n\ \text{even},\,n\ge2]\cong B^{n+1}\oplus\mathbb C[n\ \text{even}]\qquad\text{(HKR, Poincaré on }\mathbb C^{2|3})$$
[chk] $HC_n=B^{n+1}\oplus\mathbb C[n\ \text{even}][\beta=0]$ on 290 slots ($|\beta|\le5$, $n\le4$); reduced $\overline{HC}_n=B^{n+1}$; cyclic-word CE model: relative $=B^Y$, absolute $=B^Y\oplus\mathbb C[\beta=0,\,Y\ \text{odd}]$ (`dictionary/REPORT.txt` A3).
$B^p$ = exact $p$-forms on $\mathbb C^{2|3}$, with $z$, $d\theta$ even and $\theta$, $dz$ odd.
$$E_1^{\rm single}=\bigoplus_{p\ge1}B^p\ \ (\text{components of }\mathrm{Tr}(dC)^p,\ Y=p)\ \oplus\ \bigoplus_{k\ge0}\mathbb C\,\mathrm{tr}\,\mathsf c^{2k+1}.$$
All loops:
- $Q_nG=0$;
- $Q_n\,\mathrm{tr}\,\mathsf c^{2k+1}=0$ for $n\ge1$, since there are no $\mathsf b,\beta$ letters;
- at planar order the single-trace sector is closed.

This displayed $E_1^{\rm single}$ is the absolute CE statement. The LLC physical complex consists of constant-gauge-invariant cochains independent of the undifferentiated ghost at the insertion. Algebraically this is the relative complex for constant $\mathfrak{gl}_N$. The augmentation $A\to\mathbb C$ splits. In the stable cyclic-word calculation its constant cyclic homology consists of the odd ghost traces, and the reduced cyclic homology leaves $\bigoplus_{p\ge1}B^p$. Equivalently, at each finite flavor degree split the positive-jet cochain complex equivariantly for the reductive constant gauge algebra. The CE complex of its split semidirect extension by the constant algebra then has cohomology equal to the relative invariants tensored with the cohomology of the constant algebra. This removes the constant primitive ghost traces on taking relative primitive classes. The $p=1$ representatives that survive for $U(N)$ are differentiated or matter letters in the decoupled center; $SU(N)$ has $p\ge2$. [Relative physical operators and constant-ghost condition](../../sources/LLC_main.tex#L151)

## 2. Bulk: single-particle states [drv; chk at $t^0$]

$\bar\partial$-cohomology of $PV(X)$:
$$A_0=\mathbb C[z,w]\otimes\Lambda[\xi_z,\xi_w],\qquad A_2=\mathbb C[z]\otimes L\otimes\Lambda[\xi_z,\xi_w],\qquad L=H^2(\mathbb C^3\setminus0,\mathcal O)=\mathrm{span}\{w^{-a-1}\}.$$
The normal Dolbeault cohomology has rows zero and two. Products of two degree-two normal representatives vanish. Their primitives under multiplication by $w_i$ enter higher operations: three commuting tangent shears give a nonzero secondary action on a flux mode. [Explicit descent](README.md#simultaneous-tangent-shears-of-a-laurent-scalar)
$$F=f\,\xi_{z_1}\xi_{z_2},\quad f=(w_1w_2w_3)^{-1},\qquad D=t\partial+N\{F,\cdot\},\qquad \{F,\cdot\}:A_0\to A_2,\qquad \{F,\cdot\}|_{A_2}=0,\qquad\{F,F\}=0 .$$
The last two hold because $H^2\cdot H^2=0$, i.e. $F\wedge F=0$ in degree $(0,4)$ on $\mathbb C^3$. Therefore
$$0\to\operatorname{coker}\varphi_*\to H(D)\to\ker\varphi_*\to0,\qquad \varphi_*=[N\{F,\cdot\}]:H(A_0[[t]],t\partial)\to H(A_2[[t]],t\partial).$$

Dictionary $A_2\cong\Omega(\mathbb C^{2|3})$ as complexes:
- $w^{-a-1}\leftrightarrow(-1)^{|a|}(d\theta)^a/a!$;
- $\xi_{w_I}\leftrightarrow\theta_I$;
- $\xi_z\leftrightarrow$ the complementary $dz$;
- $\partial\leftrightarrow d$;
- weights $\beta(W)=(W_z+1,\,-W_w-1)$.

$$\{F,\ker\partial_{A_0}\}=\partial(F\cdot\ker\partial_{A_0}),\qquad F\cdot\ker\partial_{A_0}\leftrightarrow\Omega^0\ \Rightarrow\ \{F,\ker\partial_{A_0}\}\leftrightarrow B^1 .$$
$$\operatorname{coker}\varphi_*=B^{\ge2}\ \oplus\ \bigoplus_{k\ge0}\mathbb C\,F\,t^k .$$
[chk] At $t^0$, $\dim\operatorname{coker}(W)=\dim B^{\ge2}(\beta)+[\beta=0]$ on all $60$ slots $\beta\le(2,2;2,2,2)$ (up to the $S_2\times S_3$ symmetry), with $0$ mismatches. The bulk side is computed from the divergence on the Laurent module and multiplication by $F$; the boundary side from the de Rham complex of $\mathbb C^{2|3}$ ($d^2=0$ checked).

## 3. Comparison

In the $t^0$ comparison, after fixing the flux and removing its one-dimensional variation, the proposed physical dictionary reads
$$
H^{\rm planar\ single,relative}_{\rm bdy}(U(N))
=\big(\operatorname{coker}\varphi_*|_{t^0}/\mathbb C F\big)\oplus B^1,\qquad
B^1=\{d\,\mathrm{Tr}\,C\}\quad(U(1)\text{ single letters}).
$$

- For $SU(N)$ the decoupled $B^1$ sector is omitted.
- The assignment $\mathrm{tr}\,\mathsf c^{2k+1}\leftrightarrow F\,t^k$ compares the absolute CE tower at charge zero. The physical relative complex omits these undifferentiated ghosts.
- $\ker\varphi_*$ (regular modes) $\leftrightarrow$ sources dual to single-trace operators; this is the dictionary assumption of Costello–Gaiotto type [open].

## 4. Loop orders [drv]

For one input fundamental trace the $n$-loop LLC fan has $n+2$ ribbon vertices, including the input, and $2n+1$ internal edges. For handle number $h$, $b$ nonempty output traces and $f$ closed color faces, Euler's relation gives
$$
(n+2)-(2n+1)+f=2-2h-b,\qquad
\hbar^nN^f=(\hbar N)^n N^{1-b-2h}.
$$
Use the ordinary fundamental traces and keep $\hbar N$ fixed. There are two differentiated output letters at the endpoints, so every nonzero term has $b\ge1$. The leading power at each order retains precisely the one-boundary planar sector. If every trace is rescaled by $1/N$, a coefficient with $b$ outputs acquires the extra factor $N^{b-1}$; in that convention the same complex requires projection to one boundary. The fan outputs contain differentiated ghosts, and the classical undifferentiated ghost acts by the total constant gauge transformation, zero on invariant words. The relative complex is consequently preserved. [Fan rule in the user's paper](../../sources/LLC_main.tex#L201)

**H3 as a filtered linear statement.** Take the stable single-trace limit at fixed charges. Choose a projection $\pi$ and inclusion $\iota$ of the classical cohomology representatives above, and a contracting homotopy $s$ with $Q_0s+sQ_0=1-\iota\pi$. They can be chosen separately in each finite charge space: split off representatives, then pair each boundary with a preimage. Thus $s$ preserves the flavor and derivative counts and lowers the letter count $Y$ by one. Let $\Delta=\mathbf Q-Q_0$ for the leading one-boundary differential. The single-graviton identity gives $\Delta\iota=0$; the linear $U(1)$ representatives also obey it because a positive-loop fan contracts at least two letters. From $\mathbf Q^2=0$ one obtains
$$
\mathbf Q(1+\Delta s)-(1+\Delta s)Q_0
=\Delta(1-Q_0s-sQ_0)=\Delta\iota\pi=0 .
$$
The operator $1+\Delta s$ has invertible constant part. Its order-$n$ term shifts the charges by $-n\mathbf d$ and $Y$ by $-n$, so the inverse also exists in the finite charge filtration. At the level of individual planar operators the resulting two-loop identity is
$$
Q_2=[Q_2s,Q_0]-Q_1^2s .
$$
Thus P1 holds with $Q'_1=0$ for filtered linear changes of basis. This argument establishes a vector-space conjugacy. A local field redefinition compatible with products would require additional identities. [Single-graviton cancellation](../../sources/LLC_main.tex#L4125) · [LLC charge shifts](../../sources/LLC_main.tex#L443)

**A planar and one-handle $Q_2$ on the same cochain.** Take the physical gauge-invariant word $\operatorname{tr}(\beta_1\beta_2\beta_3\gamma^I)$ with all fields undifferentiated. Only the three $\beta$ slots can contract: replacing a $\beta$ by the $\gamma$ in the three spokes leaves relative Grassmann degree four under a three-variable integral. Write $\operatorname{tr}(t_A t_B)=\tau\delta_{AB}$ to retain the given trace normalization. In the two-loop formula the color tensor is
$$
\sum_{\sigma\in S_3}\sum_{A,B,C,R,S}
\operatorname{tr}(t_A t_B t_C t_F)
C_{\sigma(A)DR}C_{\sigma(B)RS}C_{\sigma(C)SE}
=-\tau^3(N^2+6)C_{DEF}.
$$
Here $\sigma$ permutes the three input slots. To obtain the identity replace the chain by $i/\tau$ times the trace of $[t_{\sigma(C)},[t_{\sigma(B)},[t_{\sigma(A)},t_D]]]t_E$. Contract the three internal generator pairs in the two fundamental traces. Each contraction contributes $\tau$; the center added when using $U(N)$ completeness drops out of the commutators for $SU(N)$. The $6\cdot2^3$ terms sum, after extracting $\tau^3$, to $(N^2+6)\operatorname{tr}([t_D,t_E]t_F)$. For $t_i=\sigma_i/2$ in $SU(2)$, the $(D,E,F)=(1,2,3)$ component computed directly with three Levi-Civita tensors is $-5/4$. It follows with the manuscript's dotted-index contractions that
$$
Q_2\operatorname{tr}(\beta_1\beta_2\beta_3\gamma^I)=
\frac{ig^2\tau^3(N^2+6)}{48(32\pi^2)^2}
C_{DEF}(\gamma^I)^F
(\partial_{+\dot a}\partial_{+\dot b}c^D)
(\partial_+{}^{\dot a}\partial_+{}^{\dot b}c^E).
$$
The one-handle to planar ratio is $6/N^2$. Removing the spectator $\gamma^I$ instead gives zero at every rank: an invariant two-index color tensor on the simple adjoint is symmetric, whereas the two odd ghosts with two dotted contractions are antisymmetric in their colors. The four-letter cochain need not represent a classical class. For $I=1$, its classical variation has a term proportional to $\operatorname{tr}([\gamma^2,\gamma^3]\beta_2\beta_3\gamma^1)$ that is nonzero on Pauli generators even with $\beta_1=0$. [Two-loop three-fermion fan](../../sources/LLC_main.tex#L4079) · [Exact ribbons and independent Pauli computation](code/two_loop_trace.py)

**Bulk comparison.** In the state dictionary of §2, with the three odd charges written first, $N\{F,-\}$ shifts $(\theta_1,\theta_2,\theta_3;n_+,n_-)$ by $+(1,1,1;-1,-1)$, since the weight of $F$ is $(-1,-1;-1,-1,-1)$ and the normal weights in $\beta$ are negated. LLC $Q_n$ shifts them by $-n(1,1,1;-1,-1)$. A direct identification of these arrows requires a duality reversing the weight shift. The vanishing of $F^2$ does not imply the vanishing of a boundary fan on cochains.

One explicit $t^0$ state is $2w_1^{-3}w_2^{-1}w_3^{-1}\xi_{z_1}\xi_{z_2}$, the second normal derivative of the unit flux representative. Its weight is $W=(-1,-1;-3,-1,-1)$. Here $A_2$ has two monomials, the other being $w_1^{-2}w_2^{-1}w_3^{-1}\xi_{z_1}\xi_{z_2}\xi_{w_1}$ with nonzero divergence. The regular sector at $W-W(F)$ is empty. Thus the cokernel in §2 has dimension one and its representative maps to $(d\theta_1)^2=d(\theta_1d\theta_1)$ by the stated factorial dictionary. The corresponding boundary graviton is the explicit $\delta^{XY}(\gamma^1)_X(\gamma^1)_Y$ in the user's paper. Its charges are $(2,0,0;0,0)_2$. An outgoing $Q_n$ would require negative second and third flavor counts and an incoming $Q_n$ would require negative derivative counts. Hence every positive-loop boundary differential into or out of this class vanishes, including every handle contribution. This identifies the linear bulk state and supplies a protected boundary control for the quantum dictionary. [Exact mode count](code/single_particle.py) · [Observable and its superfield component](../../sources/LLC_main.tex#L4145)

For a concrete H4 selection test take the two outgoing ghost momenta to be dual coordinate covectors in the boundary $\mathbb C^2$ and keep the scalar spectator constant. The contraction of their two dotted derivatives gives their squared two-row determinant, equal to one. On embedding just these momenta into the tangent $\mathbb C^2\subset\mathbb C^5$, the flat fivefold positive-loop heat map has rank at most $2n$ in its $5n$ cycle-current coordinates; its local top form vanishes before integration. This excludes the flat elementary or product kernels with only these tangent injections as the source of the nonzero boundary handle term. The punctured bulk also has transverse Laurent modes and the $NF$ background, which were absent from that restricted calculation. The $K3$ and $T^4$ computations concern a threefold effective bulk with its own cohomology-ring labels. The four-dimensional $\mathcal N=4$ observable above specifies neither those labels nor a boundary map for either compactification. [Flat fivefold and compactification kernels](../iib-bcov/README.md) · [Punctured fivefold background, §9](https://arxiv.org/abs/1606.00365)

The two protected-product computations also constrain an attempted identification. The fivefold level-$r$ point class $[E_r]_5$ has tree coordinate degree $2r-1$, and each loop correction preserves its parity. Hence $[E_0]_5\star[E_1]_5$ is even. In the manuscript's explicit marked-product example $A=\delta^{XY}(\gamma^1)_X(\gamma^1)_Y$ is even and $B=(\gamma^2)_X(\beta_1)^X$ is odd, making $A\star B$ odd. The two Ward cancellations therefore cannot be matrix elements of the same pair under a parity-preserving multiplicative identification. [Fivefold seven- and eight-field tests](../iib-bcov/operator/README.md#quantum-products-and-the-constant-field-test) · [Boundary pair](../../sources/LLC_main.tex#L4145)

## Signed component of the boundary fan

For the literal input $\operatorname{tr}(\beta_1\beta_2\beta_3 b^{n-2}\gamma^I)$ fix the even spectator $\gamma^I$ last. The three operators $-\theta_I\partial_L/\partial\beta_I$ have even parity; every assignment yields $-\theta_1\theta_2\theta_3$. Each contracted $b$ gives $-ig$ and the two endpoint ghosts contribute $(ig)^2$. Converting $n+1$ structure constants into a nested matrix commutator and using fundamental completeness gives $i^{-(n+1)}\tau^n$. Thus the phase and couplings multiplying the color traces are
$$
(-1)^ni^{n+1}\big[-(-ig)^{n-2}(ig)^2\big]i^{-(n+1)}\tau^n
=-(ig\tau)^n .
$$
Multiply this by the color polynomial in the all-order calculation, divide by $2(32\pi^2)^n(n!)^3(n+1)$, and contract
$$
\left[\operatorname{tr}(t_Dt_Et_F)+(-1)^{n+1}\operatorname{tr}(t_Dt_Ft_E)\right]
(\gamma^I)^F
(\partial_{+\dot a_1}\cdots\partial_{+\dot a_n}c^D)
(\partial_+{}^{\dot a_1}\cdots\partial_+{}^{\dot a_n}c^E).
$$
All labeled assignments of the identical $b$ fields are included in that polynomial. For the physical-letter word $\operatorname{tr}(\psi_+^1\psi_+^2\psi_+^3 f_{++}^{n-2}\bar\phi_I)$ replace the two differentiated ghosts by the two symmetrized $(n-1)$-derivative antichiral gauginos and multiply the coefficient by $-8\sqrt2\,i\,4^{n-2}$, using the manuscript's physical-letter map. [Zero-displacement fan](../../sources/LLC_main.tex#L4019) · [Component and physical-letter conventions](../../sources/LLC_main.tex#L4076) · [Color polynomial](README.md)

At $n=3$, take $D=E=t_1$ in the Gell-Mann convention $t_i=\lambda_i/2$, and take the scalar spectator color matrix $\operatorname{diag}(1,1,-2)$. A direct contraction with the six off-diagonal matrix units, two traceless diagonal generators, and their trace-dual basis gives the four-structure-constant color tensor $63/8$. The twisted coefficient is $+63ig^3/[8\cdot1728(32\pi^2)^3]$, and the physical-letter coefficient is $+7\sqrt2\,g^3/[48(32\pi^2)^3]$. These refer to the summed dotted contraction as displayed. [Direct traceless three-by-three computation](code/fan_component_check.py)

## 5. Further tests [open]

- **T-b.** The filtered linear statement and a nonzero planar $Q_2$ cochain are computed in §4. A bulk-induced local change of variables remains open.
- **T-c.** $1/N^2$: one-loop BCOV on $X$, using the $T=0$ kernels at $(H,T)=(5,0)$, against non-planar $Q_n$.
- **T-d.** Finite $N$: the LLC census ($E_1$, fortuitous classes, $d_r$) against the giant-graviton quantization conjecture for the finite-$N$ BPS Hilbert space, Budzik–Gaiotto arXiv:2611.18539 (body not read).

## 6. Results of the parallel checks (Claude, 2026-11-21)

Sources: `dictionary/REPORT.txt` (tools in `../../code/holography/`, tests pass) and `verify/VERIFY.txt`. All [chk] are exact over $\mathbb Q$.

**Bulk complex.**
- [drv] The Čech transfer leaves $D=t\partial+N\{F,-\}$ exact on $T=A_0[[t]]\oplus A_2[[t]]$. Finite grading: $\hat W=W-(1^5)[A_0]$, $h=k-[A_0]$, with $T(\hat W,h)=t^hA_{2,\hat W}\oplus t^{h+1}A_{0,\hat W+1}$ of dimension $\le64$.
- [chk] $D^2=0$ on $206143$ basis elements, and $\Phi\circ\partial=d\circ\Phi$ on $16807$ monomials.
- [chk] On $6300$ pieces through $t^3$:
$$\operatorname{coker}\varphi_*=B^{\ge2}(\beta)\oplus\textstyle\bigoplus_{k\ge0}\mathbb C\,Ft^k,\qquad \ker\varphi_*=\{\alpha\in A_0:\partial\alpha=0,\ \alpha|_{w=0,\xi_z=0}\in\mathbb C\}\oplus\textstyle\bigoplus_{k\ge1}\mathbb C\,\xi_{\rm top}t^k .$$

**Index.**
- [chk] Through $t^{30}$, with the sector monomial $(-1)^{Y+\Theta}t^{2\Theta+3n}y^{n_+-n_-}v^{\theta_1-\theta_2}w^{\theta_2-\theta_3}$:
$$\mathrm{sch}\,B^1=f\ \ (\text{letter index}),\qquad \textstyle\sum_{p\ge1}\mathrm{sch}\,B^p=I_{\rm sp}\ \text{(KMMR (4.14))},\qquad \sum_{p\ge2}\mathrm{sch}\,B^p=I_{\rm sp}-f .$$
- The bulk normalizable spectrum therefore reproduces the supergravity single-graviton index with the $n=1$ multiplet removed.

**Sources vs operators.**
- [chk] Weight rule: $W_r=\beta(W_n)-(1^5)$, with $t^k\leftrightarrow t^k$.
- [chk] On 525 weights:
$$\dim\ker\varphi_*(\beta-1^5)-\dim\operatorname{coker}\varphi_*(\beta)=[\beta\in M_1]+[\beta=1^5],\qquad M_1=\{\beta\ne0:\beta_\theta\in\{0,1\}^3,\ \min\beta_z=0\}.$$
- At $\beta=1^5$ the extra regular mode lies in the span of the constant $1$ and the Cartan vector fields.
- [drv] A $D$-invariant pairing of this weight type has $B^1$ in its radical, so no perfect source/operator pairing exists.

**Product-compatible P1** (version (M): $e^{-\xi}\mathbf Qe^{\xi}=Q_0+\hbar Q_1'$ with $\xi$ an algebra derivation).
- [chk] (M) fails at $\hbar^2$ at finite $N$. Algebra automorphisms do not raise the Koszul–Akman order, so (M) requires $\Phi^3_{Q_2}=0$, and
$$\Phi^3_{Q_2}(\mathrm{tr}\,\beta_1\gamma,\mathrm{tr}\,\beta_2\gamma,\mathrm{tr}\,\beta_3\gamma)\neq0,$$
with $SU(2)$ value $ig^2(x\cdot x)(x_1B^{23}-x_2B^{13}+x_3B^{12})/(32768\pi^4)$.
- At fixed $\hbar N$ the obstruction is of order $N^{-2}$: three traces joined into one. At $N=\infty$, (M) reduces to the linear statement, which holds.
- [drv] At order $1/N$ the joining obstruction vanishes because multi-gravitons are $\star$-closed.
- Bulk counterpart [open]: tree-level four-point exchange with two flux insertions on separate vertices ($F^2=0$ only per vertex).

**Codex's values, recomputed independently.**
- [chk] Colour polynomials $N^2+6$, $N^3+12N$, $N^4+36N^2+120$, $N^5+63N^3+360N$.
- [chk] Further polynomials:
  - $n=6$: $N^6+127N^4+1996N^2+5040$;
  - $n=7$: $N^7+199N^5+4960N^3+20160N$;
  - $n=8$: $N^8+331N^6+16732N^4+169632N^2+362880$.
- [chk] One-handle formula for $n\le300$; values $-5/4$ and $63/8$.
- [chk] $\sum_hB_{m,h}N^{m-2h}=[\prod_{j=0}^m(N+j)-\prod_{j=0}^m(N-j)]/(m(m+1))$ for $m\le9$.
- [chk] Leading coefficients $1/(2^{2h}(2h+1)!)$ and $4h/(2^{2h}(2h+1)!)$ for $h\le3$; the $\sinh(x/2)/(x/2)$ limit numerically.
- [chk] Flux hexagon: $N/12$, $-N/6$, determinant $8$.
- [chk] Genus-two output $-N/(1296\,w_1w_2w_3^3)\,\partial_{z_2}$, with the source's BV sign conventions taken as given.

**Finite $N$: Budzik–Gaiotto against the census** (`finiteN/REPORT.txt`; tools `bg_words`, `bg_complex`, `bg_complex_k`, `uN_q0`, `rank_mod`, `giant`, `census_io` in `../../code/holography/`).
- [src] Budzik–Gaiotto (1.7) = (4.51), in PDF numbering:
$$\lim_{k\to\infty}H_\bullet\big(\mathfrak{gl}'_{k|k}[A^!]\,\big|\,\mathfrak b_{k|k},\mathrm{Ber}^N\big)\cong H^\bullet\big(\mathfrak{gl}_N[A]\,\big|\,\mathfrak{gl}_N\big)=E_1^{U(N)},$$
  with $A=\mathbb C[\theta_{1,2,3},z_{1,2}]$, $A^!\simeq\mathbb C[x_{1,2,3},\zeta_{1,2}]$ (Koszul-dual D5-branes), tree level only. The paper's own tables cover only $A=\mathbb C[\theta]$ and $\mathbb C[z]$.
- [chk] The meson complex reproduces BG Tables 1–5 ($k=1$) and Tables 1–2 ($k=2$).
- [chk] For the $1/16$-BPS case at $k=1$, $H^{BG}=E_1^{U(N)}$ sector by sector:

| $N$ | levels | sectors agreeing |
|---|---|---|
| 1 | $\le9$ | 23/23 |
| 2 | $\le12$ | 54/54 |
| 3 | $\le10$ | 32/32 |
| 4 | $\le10$ | 32/32 |

  At $k=2$, $U(2)$ agrees through level 6. The agreement includes the finite-$N$ trace relations: $U(2)$ differs from $U(\infty)$ in 27 of 38 representatives. Largest sector: $(2,2,2;0,0)$, 115,596 chains, $H=1,9,29,35,14$ at $Y=2..6$.
- [chk] Index: $Z_2=Z_\infty(1+G^{(1)}+G^{(2)}+G^{(3)})$ through level 33, with onsets at levels $6,16,30=2k(N+k)$. The census index equals the Molien–Weyl $SU(2)$ index for levels $2..33$.
- [drv, chk to level 11] $E_r^{U(N)}=E_r^{SU(N)}\otimes\mathrm{Sym}(V_1)$ and $d_r^{U(N)}=d_r^{SU(N)}\otimes1$.
- Limits:
  - $k=1$ fails from level 16: a spurious ghost-number $-3$ class in $(8,0,0;0,0)$.
  - The $SU(2)$ fortuitous $\mathcal O_0$ (level 24) needs $k\ge2$; P39, X40 and $\mathcal O_1$ need $k\ge3$. None is tested.
- [drv] Criterion in BG language: a class is fortuitous iff every representative needs an antifield ($u,\bar u$) tail, i.e. it closes only through $\Delta_{BV}$.

**Finite $N$ at $k=2$** (`finiteN_k2/REPORT.txt`; tools `bg_fast`, `bg_sizes`, `native/bgk.cpp`; ranks mod $p=2^{31}-1$ and $p'=2^{31}-19$).
- [chk] `bg_fast` $=$ `bg_complex(_k)` on the 151 sectors of `finiteN/REPORT.txt` ($p,p'$), entry by entry for $k=2,3$. Reduced $b$-set $\{P0_{12},MM0_{12},R0_{12}\}$ $=$ all eight $b$-operators.
- [chk] $\mathrm{res}(\lambda)=[Z_{U(2)}]_\lambda-[Z_\infty(1+G_2^{(1)})]_\lambda$: at level 16 nonzero exactly for $\Theta=8$, $n=0$, with orbit sum $126=[G_2^{(2)}]_{x^{16}}$. In every computed sector, $H^{BG}_{k=1}=E_1^{U(2)}\oplus(\text{classes with }\chi=-\mathrm{res}(\lambda))$.
- [chk, $p,p'$] $U(2)$, $k=2$:

| sector | level | $g$ | $\dim H^{BG}_{k=1}-\dim E_1^{U(2)}$ | $\dim H^{BG}_{k=2}-\dim E_1^{U(2)}$ |
|---|---|---|---|---|
| $(8,0,0;0,0)$ | 16 | all | 1 | 0 |
| $(7,1,0;0,0)$ | 16 | $-3$ | 1 | 0 |
| $(7,0,0;1,0)$ | 17 | $-2$ | 1 | 0 |
| $(9,0,0;0,0)$ | 18 | $-4..-2$ | 2 | 0 |
| $(10,0,0;0,0)$ | 20 | $-3$ | 5 | 0 |

  $(8,0,0;0,0)$: $H^{BG}_{k=2}=\{Y=8:5\}=\mathbb C[\mathrm{Tr}X,\mathrm{Tr}X^2]_8$ over $g=-16..0$ ($g\le-9$: `finiteN_k2/e1_lowg_check.log`, $p$).
- [chk, $p$] One-matrix sector $(n,0,0;0,0)$, $k=2$:
  - $N=0$: $H=0$ for $n=4..8$; $H_{-8}=1$ at $n=9$ (also $p'$);
  - $N=1$: $H=\mathbb C[\mathrm{Tr}X]_n$ for $n=6,7,8$;
  - $N=2$: $H=\mathbb C[\mathrm{Tr}X,\mathrm{Tr}X^2]_n$ for $n=5,6,7$, over $g=-2n..0$.
- Onset rule. The first class of $H^{BG}_k$ outside $E_1^{U(N)}$ sits at the onset of $G_N^{(k+1)}$, at level $2(k+1)(N+k+1)$.
  - [chk] for $(k,N)=(1,0),(1,1),(1,2),(2,0)$; no excess below the onset at $(2,1)$, $(2,2)$ where computed.
  - $(k,N)=(2,2)$: $30>24$ $\Rightarrow$ [conj] $\mathcal O_0\otimes1\in H^{BG}_{k=2}$.
- [chk] $\mathcal O_0$ sector $(4,4,4;0,0)$, $g=-5$:
  - $\dim V_T$: $9.39\cdot10^8$ ($k=1$); $1.055\cdot10^{12}$ ($k=2$, $\det^N$-isotypic part $3.42\cdot10^{10}$); $>2.3\cdot10^{14}$ ($k=3$, [conj]);
  - gauge cochains: $\dim C^{U(2)}_{Y=7}=1821$, $\dim C^{SU(2)}_{Y=7}=153$;
  - not computed (storage $\ge6\cdot10^{13}$ B).
- [chk] Strict $HA^!$ model ($\mathbb C[x]$ commutative, $HC_{\ge1}$ as spectators): $E_2\ne H^{BG}$ from level 8. At $(2,1,1;0,0)$: $E_2=\{2:5,3:12,4:6\}$, $H=\{2:3,3:10,4:6\}$ $\Rightarrow$ $d_{r\ge2}\ne0$.
- Not done: $\mathcal O_0$; the antifield-tail criterion; 100 of the 105 $k=1$-failing representatives at levels 16–20.

**Bulk symmetry algebra $\mathfrak g=\ker\varphi_*|_{t^0}$** (`symmetry/REPORT.txt`; tool `symmetry.py`, tests pass).
- [drv, chk] Dictionary to polyvectors:
$$\Psi:A_0\xrightarrow{\sim}PV(\mathbb C^{2|3}),\quad z\mapsto z,\ w_I\mapsto-\partial_{\theta_I},\ \xi_{z_a}\mapsto\partial_{z_a},\ \xi_{w_I}\mapsto\theta_I;\qquad \Psi\partial=\Delta\Psi,\quad \Delta|_{PV^1}=\mathrm{div}_{\rm Ber},$$
$$\Phi\{\alpha,\zeta\}=(L_{\Psi\alpha}-\iota_{\Psi\partial\alpha})\Phi\zeta.$$
- [drv, chk] The flux map is $\Phi\varphi(\alpha)=d(\Psi\alpha)_0$, with $(\Psi\alpha)_0$ the function part. Hence
$$\mathfrak g=\mathbb C\cdot1\times\Psi^{-1}PV^{\ge1}_\Delta(\mathbb C^{2|3}),\quad[\mathfrak g_k,\mathfrak g_l]\subset\mathfrak g_{k+l-1},\quad\mathfrak g_1=S(2|3),$$
$$\dim\mathfrak g(W)=\dim Z(A_0)(W)-\dim B^1(W_z,-W_w).$$
- [chk] $\mathfrak{psl}(3|3)\subset\mathfrak g_1$:
  - $\psi(E_{PQ})=-V_{PQ}$, with $V_{PQ}=y_P\partial_{y_Q}$ and $V_{P0}=-y_PE$ on $\mathbb P^{2|3}$;
  - $\mathrm{div}_{\rm Ber}(y_PE)=-y_P+y_P=0$, since $\mathrm{div}E=2-3$;
  - homomorphism on 1156 pairs; same list as Costello–Li Lemma 9.2.1 [src];
  - finite part of $\mathfrak g$: $\mathfrak{psl}(3|3)\times\mathbb C\cdot1$.
- [drv] Excess:
$$\dim\mathfrak g(\beta-1^5)-\dim\operatorname{coker}(\beta)=B^1(\beta)-B^1(\gamma(\beta)),\quad\gamma(\beta)=(\beta_z-(1,1),(1,1,1)-\beta_\theta),$$
$$\Rightarrow\quad\dim\mathfrak g(\beta-1^5)-\dim\operatorname{coker}(\beta)=[\beta\in M_1]+[\beta=1^5].$$
  [chk, independent code `verify/c5_symmetry_spotcheck.py`, 980 weights] $\mathrm{rank}\,\varphi|_{Z(A_0)(W)}=B^1(W_z,-W_w)$, and the excess formula.
- [drv] Higgs pairs: $\varphi:O(\mathbb C^{2|3})(\gamma)/\mathbb C\xrightarrow{\sim}B^1(\gamma)$, $f\mapsto df$. The $U(1)$ shifts $C\mapsto C+f\,1_N$ are broken by $NF$; their partners are the $U(1)$ letters $d\,\mathrm{tr}\,C$.
- [chk] Sources of the $U(1)$ letters:
  - $\mathrm{Ann}(B^{\ge2}(\beta))\cap Z(A_0)(\beta-1^5)=\mathbb Cu_\beta$, with $u_\beta=\theta^{(1,1,1)-\beta_\theta}H_{z^{\beta_z}}$, $H_h=\partial_2h\,\partial_{z_1}-\partial_1h\,\partial_{z_2}$, $H_1=\partial_{z_1}\wedge\partial_{z_2}$;
  - $u_\beta\in\mathfrak g$ (199 weights);
  - $U=\Lambda[\theta]\otimes\mathrm{Ham}(\mathbb C^2)\oplus(\Lambda[\theta]/\mathbb C\theta^{123})\,\partial_{z_1}\wedge\partial_{z_2}$ is a subalgebra; its stabilizer in $\mathfrak{psl}(3|3)$ is a 26-dimensional parabolic.
- [chk] Excess modes:
  - $\beta\in M_1$: $u_\beta$ with $z_1z_2\nmid z^{\beta_z}$;
  - $\beta=1^5$: $1$;
  - $\mathfrak{psl}(3|3)\cap\text{excess}=\{\partial_{z_a},\theta_I\partial_{z_a},z_b\partial_{z_a}\ (a\ne b)\}$.
- [chk] Action on $\operatorname{coker}$: $\Phi\,\rho(\psi A)\,\Phi^{-1}=L_{V_A}$ on $B^{\ge2}$ (34 generators, 150 forms). The only part of $\mathfrak g$ acting trivially is $\mathbb C\cdot1$ (60 weights).
- Corrected dictionary:
$$\dim\mathfrak g(\beta-1^5)=\dim B^{\ge1}(\beta)+[\beta=0]-[\beta\in R'],\qquad R'=\{\beta_\theta\in\{0,1\}^3,\ \beta_z\ge(1,1),\ \beta\ne1^5\}.$$
- [drv] Bracket on $H(D)$: $\{\alpha,\beta\}+\tfrac{N^2}{2}\ell_4(F,F,\alpha,\beta)$. On $\mathfrak{psl}(3|3)$, 7 weight-allowed pairs lie among $\{z_aE,z_b\partial_{\theta_I}\}$. Costello–Li Prop. 9.2.2 requires $[\ell_4(F,F,\cdot,\cdot)|_{\mathfrak{psl}(3|3)}]=0$ in $H^2(\mathfrak{psl}(3|3);\operatorname{coker}\varphi_*)$.
- [open] Values of $\ell_4(F,F,\cdot,\cdot)$; a natural pairing for $(u_\beta,O_\beta)$ at $\beta\in R'$.
