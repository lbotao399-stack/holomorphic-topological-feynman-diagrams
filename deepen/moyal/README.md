# Laurent counterterms for the Moyal vertex

Use the flagship's Moyal action $s_e$ on $\mathfrak{gl}_{N\ge3}$ with fundamental trace, its loop parameter $\delta$, and the choice of zero first-order quantum deformation in $\delta$. The order-$\delta^2$ counterterms are determined coefficient by coefficient in $e$ at the graded bidegrees of Stage M2. Vacuum terms are omitted. The cochains and the signs $\varepsilon,U$ are those in the [flagship](../../notes/flagship/flagship.tex).

Write $f=[p_1p_2]^2\operatorname{tr}(X_1[X_2,X_3])$ for the flagship cochain denoted there by $\mathsf f$ with its momentum factor (the color tensor alone has no momentum factor). The first terms of the classical vertex are $s_0+e\,j_M+e^2 i_M$, with $j_M=D/2$, $i_M=f/8$, $d=-\{s_0,-\}$ and $df=\beta$. The order-$e^{-2}$ term in the undressed equivariant sector has bidegree $(0,0)$. In degree three the only invariant alternating tensor on $\mathfrak{gl}_N$ is the commutator trace. Thus the leading contribution is a multiple of $s_0$, which extends to the closed cochain $e^{-2}s_e$ since $\{s_e,s_e\}=0$. Differentiation gives the second closed direction $e^{-1}\partial_e s_e$. These remove the free multiples of $s_0$ at $e^{-2}$ and of $j_M$ at $e^{-1}$. In particular the allowed $e^{-2}\delta^2$ term cannot change the coefficient below.

After these two closed directions are subtracted, the $e^{-1}$ coefficient must be closed for $d$. Its mixed trace components follow by polarizing the flagship quadratic map at $j_M$, whose coefficients in the basis $(\pi^*D_{\mathfrak{sl}},M_1,S_0)$ are $(1/2,1/N,1/N^2)$. Pairing it with coefficients $(\alpha,\chi,\gamma)$ gives
$$
[\{j_M,\alpha\pi^*D_{\mathfrak{sl}}+\chi M_1+\gamma S_0\}]
=\left(\frac{\alpha}{12N^2}-\frac{\chi}{24N}\right)[F_{223}]
 -\frac{N\gamma-\chi}{48}[F_1].
$$
The bare two-loop class has no $[F_1]$ term, so $\chi=N\gamma$. Removing the free $j_M$ direction leaves the coefficient of $\pi^*D_{\mathfrak{sl}}$ equal to
$$
t=-12N^2\varepsilon c_O=\frac{\varepsilon U N^2}{256\pi^4},
\qquad c_O=\frac{\pi\nu}{48}=-\frac U{3072\pi^4}.
$$

For the finite term use $\pi^*D_{\mathfrak{sl}}=D-2S_1/N+4S_0/N^2$. The formal trace bracket gives
$$
\{j_M,\pi^*D_{\mathfrak{sl}}\}
 =\frac{F_{223}}{12N^2}+
   \left(\frac12+\frac2{N^2}\right)\beta ,
\qquad
o_2=\pi\nu\left(\frac{F_{223}}{48}-\frac{N^2}{40}\beta\right).
$$
In the cochain equation, $\varepsilon o_2+t\{j_M,\pi^*D_{\mathfrak{sl}}\}=k\beta$ with $k=-\varepsilon\pi\nu(3N^2+10)/20$. Choose the finite counterterm $kf$.

At order $e$ the full two-loop graph correction comes from one marked $j_M$ vertex, so the flagship one-mark calculation gives $\pi\nu\,d(-N^2D_3/120+NS_{1,3}/60)$. The remaining tree terms are $\{j_M,kf\}+\{i_M,t\pi^*D_{\mathfrak{sl}}\}$. The identities $\{D,f\}=dD_3/3$, $\{S_1,f\}=dS_{1,3}/3$, $\{S_0,f\}=0$ give the following counterterm coefficients, with $\delta^2$ removed:
$$
\begin{array}{c|l}
e^{-2} & 0\\
e^{-1} & t\pi^*D_{\mathfrak{sl}}\\
e^0 & -\varepsilon\pi\nu(3N^2+10)f/20\\
e^1 & \varepsilon\pi\nu\big(-(21N^2+40)D_3/480+3N S_{1,3}/80\big).
\end{array}
$$
Adding any constant multiples of $e^{-2}s_e$ and $e^{-1}\partial_e s_e$ gives their complete classical Moyal extensions. The choices in the table fix these directions to zero; higher homogeneous counterterms may have further closed components. [Earlier matching equation](../../notes/interpretation/J1_and_5d_CS_level.md) · [Exact Laurent coefficient check](check_laurent.py)

### Formal-rank counterterms through $e^6$

Divide the quantum cubic by $\delta^2\varepsilon\pi\nu$. At power $e^n$ its momentum factor is $[p_1p_2]^{n+2}$. The columns below multiply the flagship color tensors $\mathsf f=\operatorname{tr}(1[23])$, $\mathsf d=\operatorname{tr}(1\{23\})$, $s_1=\operatorname{tr}(1)\operatorname{tr}(23)+\operatorname{tr}(2)\operatorname{tr}(31)+\operatorname{tr}(3)\operatorname{tr}(12)$ and $s_0=\operatorname{tr}(1)\operatorname{tr}(2)\operatorname{tr}(3)$. The triple-central coefficient at $e^1$ left free in the preceding table is fixed by the equation at $e^2$.

| $n$ | $\mathsf f$ | $\mathsf d$ | $s_1$ | $s_0$ |
|---|---|---|---|---|
| $1$ | $0$ | $-(21N^2+40)/480$ | $3N/80$ | $-1/40$ |
| $2$ | $-N^2/120-517/30240$ | $0$ | $0$ | $0$ |
| $3$ | $0$ | $-3N^2/2560-1/576$ | $151N/241920$ | $-1/4480$ |
| $4$ | $-N^2/7680-4813/25945920$ | $0$ | $0$ | $0$ |
| $5$ | $0$ | $-11N^2/921600-1751/129729600$ | $1241N/276756480$ | $-1/967680$ |
| $6$ | $-N^2/1075200-4877479/4940103168000$ | $0$ | $0$ | $0$ |

At each degree, [the graph calculation](theta_moyal.py) expands every cyclic vertex and integrates its Schwinger currents by Dirichlet moments. [The replay](solve_laurent.py) solves and checks the full BV equation as a polynomial in all external momenta, all formal trace words and $N$. Odd momentum powers at cubic arity multiply symmetric tensors and even powers multiply alternating tensors. This follows from conservation of the three momenta and the fact that $[p_1p_2]$ generates their $SL(2)$-invariant polynomial ring. Each odd degree admits the undetermined triple-central term $[p_1p_2]^{n+2}s_0$ because its classical CE differential vanishes. Its next bracket with $j_M$ enters at the following even degree.

### Planar single trace for every power of $e$

The disk ribbon contractions have one external boundary and two closed index loops. Expanding the two cyclic orders at each of six vertices gives $2^6$ terms per graph. Precisely two terms of $(1,3,3)$ and zero terms of $(2,2,3)$ have that boundary and those index loops. In both terms the sum of the six signed vertex brackets is the sum of $[p_ip_j]$ over the ordered pairs on the external trace. This is a polynomial identity before Schwinger integration, checked with all path parameters indeterminate in [the replay](solve_laurent.py). Thus the disk part of the entire graph is the Moyal dressing of its $e=0$ value $-N^2\beta/40$. The dressing multiplies each cyclic word by its external Moyal exponential. On a tree the two vertex phases concatenate to the external cyclic phase by momentum conservation, so it takes the primitive $[p_1p_2]^2\mathsf f$ of $\beta$ to $[p_1p_2]^2s_e=4\partial_e^2s_e$.

In the pole convention above the disk counterterm, divided by $\delta^2\varepsilon\pi\nu$, is
$$
-\frac{N^2}{10}\left(\partial_e^2s_e+\frac5e\partial_e s_e\right).
$$
The derivative term is closed and fixes the disk pole $-N^2 D/(4e)$. Consequently the disk coefficient at every $n\ge-1$, with the cyclic color tensor and $[p_1p_2]^{n+2}$ removed, is $-N^2(n+2)(n+6)/(10\cdot2^{n+2}(n+2)!)$.

### Abelian and decoupled central Poisson cubic

At rank one on the external matrix $1$, the Moyal vertex is $2\sinh(e[p_1p_2]/2)$. All six vertices start at $e$. The first bare two-loop coefficient occurs at $\delta^2 e^6$. On momenta $(1,0),(0,1),(-3,2),(2,-3)$ the two graph values, in units of $\pi\nu$, give
$$
\frac{4823}{56100}-\frac{1483}{336600}=\frac{323}{3960}.
$$
With the cubic cochain identified with its polynomial on three scalar waves, [the direct integration](abelian_poisson.py) gives the full four-input polynomial
$$
\frac{323}{259459200}\,\big\{[p_1p_2],[p_1p_2]^7\big\}.
$$
The bracket is the flagship $(2,2)$ shuffle with the unit scalar pairing. Its value on these waves is $65520$. The rank-one part of the table at $e^1,e^3$ equals the restriction of the closed term $-2e^{-2}(e\partial_e-1)s_e$. Subtracting this term, the $e^5\delta^2$ cubic left in the same units $\varepsilon\pi\nu$ is $-323[p_1p_2]^7/259459200$.

For a decoupled center of $\mathfrak{gl}_N$ with $\kappa=\operatorname{tr}$ and first cubic $\hbar\gamma S_0$, each internal contraction of two one-letter traces gives $N$. Seven internal edges give $N^7$ and the four external legs give the product of their traces. At rank one $N=1$ with the external matrix $1$ all four traces equal one. The two-loop contact at $\hbar^8$ is the polynomial above times $\hbar^8\pi\nu N^7\gamma^6$ and the four external traces. With only the Poisson cubic at first order, choose the cancelling cubic at order seven to be
$$
-\hbar^7\varepsilon\pi\nu\,
\frac{323N^6\gamma^5}{259459200}\,
[p_1p_2]^7\operatorname{tr}(X_1)\operatorname{tr}(X_2)\operatorname{tr}(X_3).
$$
Its bracket with $\hbar\gamma S_0$ supplies the seventh power of $N$. The scalar Poisson cubic squares to zero by the Jacobi identity. This calculation uses no higher classical Moyal vertex. For the fully Moyal theory the classical $e^7$ cubic belongs to the separate tree action.

Recompute the results with `python3 abelian_poisson.py` and `python3 solve_laurent.py 6` in this directory.

### Finite formula for every formal coefficient at two loops

Use the classical vertex $s_e$ and write its two-loop graph cochain in units of $\pi\nu$. Its BV bracket with $s_e$ vanishes by curvature consistency at order $\delta^2$: the Moyal product is associative, all one-loop contacts vanish, and the free differential vanishes on full-wedge holomorphic local vertices. Subtracting the equation solved through $e^{n-1}$ leaves a closed coefficient at $e^n$ of momentum bidegree $(n+2,n+2)$. The invariant classical cohomology at this weight, and its one-order-higher consistency test, are computed in [the matrix calculation](../o3/higher_orders/center_classes.md). Apply that argument with first classical derivative $j_M$, whose mixed coefficient is $1/N$. At even $n\ge2$, the previous free triple-central cubic removes the two-central class. At odd $n$ this class vanishes. After taking a primitive of the other noncentral components, the bracket of $j_M$ with the residual central four-cochain is a classical boundary by the next coefficient of consistency. Two equal commuting traceless matrices and three scalar matrices detect this bracket and force the central polynomial to vanish. This proves existence of a primitive at every $n$. The only free cubic at a fixed $n\ge2$ has three central legs and odd $n$; the next even equation fixes it. All statements concern the formal Laurent expansion with the two pole redefinitions fixed above.

The [two-probe calculation](probe_laurent.py) gives the coefficients in finite rational arithmetic at every degree. At order $n\ge2$, form the bare graph plus the brackets of the positive powers of $s_e$ with the counterterms already obtained, taking the undetermined preceding central coefficient as zero. Put $k=n+2$. The following extractions on this four-cochain give all cubic coefficients at this step.

| trace word | successive momenta | coefficient extracted |
|---|---|---|
| $(0123)$ | $(1,0),(0,1),(1,1),(-2,-2)$ | divide by $2^k+(-1)^k-1$ for $\mathsf f$ at even $k$ or $\mathsf d$ at odd $k$ |
| $(0)(123)$, odd $k$ | same | divide by $2^k-2$ for $s_1$ |
| $(01)(2)(3)$, even $k$ | $(1,0),(0,1),(-1,0),(0,-1)$ | negate for the $s_0$ coefficient at $n-1$ |

The single-trace divisor follows from merging neighbors in the ordered four-cycle. The double-trace divisor is $-([p_0p_1]^k+[p_0p_2]^k+[p_0p_3]^k)$. The triple-trace divisor before the coefficient is moved across the equation is $[p_0p_1][p_2p_3]^{k-1}=1$. The disjoint trace words isolate these three columns. For an odd requested truncation, evaluate one extra even order to fix its final triple-central cubic.

For each probe, expand the $2^6$ cyclic ribbon assignments per theta with the external labels alternated. Each assignment contributes the $n$-th power of its summed signed vertex brackets divided by $2^n n!$, multiplied by the bare Schwinger density, trace contraction and graph symmetry factor. For path lengths $n_P$ the normalized moments of the three conductance weights have Dirichlet parameters $7/2-n_P$. The path ratios have parameters $(1,\ldots,1)$ independently on each path. For example a monomial $w_A^i w_B^j$ has mean $(7/2-n_A)_i(7/2-n_B)_j/(7/2)_{i+j}$, where $(x)_i=x(x+1)\cdots(x+i-1)$. On a path of length $n_P$ the corresponding moment of $\prod_j(x^P_j)^{i_j}$ is $\prod_j i_j!/(n_P)_{\sum_j i_j}$. Thus every coefficient reduces to explicit factorials, trace shuffles and the nonzero integer divisors in the table. The recurrence leaves no period to be evaluated at a specified finite order. Its output through $n=6$ agrees with the full four-input polynomial calculation above.

### Triple trace for every power of $e$

For the $(01)(2)(3)$ opposite-pair probe above, multiply the bare graph coefficient by $e^4$. With shares $w_A+w_B+w_C=1$, its path-ratio integral, expressed in the symmetric Dirichlet measure with parameters $(1/2,1/2,1/2)$, is
$$
\frac{\prod_{P\in\{A,B,C\}}(e^{e w_P}-1)+\prod_{P\in\{A,B,C\}}(e^{-e w_P}-1)}
{w_Aw_Bw_C}.
$$
Indeed the density changes from the symmetric measure by $5w_A^2$ for $(1,3,3)$ and $15w_Aw_B$ for $(2,2,3)$. Each signed vertex phase in the three-trace channel is affine in the four path ratios at fixed shares. Integration of $\exp(e\sum_j y_j x_j)$ over the normalized uniform simplex of length $l$ gives $(l-1)!e^{1-l}$ times the divided difference of $y\mapsto\exp(ey)$ at its $l$ vertex values. Repeated values use the corresponding derivatives. [The symbolic replay](central_phase.py) expands both graph shapes and symmetrizes the three shares. For each exponential frequency and power of $e$ it verifies the displayed identity over $\mathbb Q(w_A,w_B)$.

For even $m\ge4$ the coefficient of $e^m$ after the remaining Dirichlet expectation is
$$
\frac{2}{(3/2)_{m-3}}
[t^m]\left(\sum_{\substack{i\ge1\\\text{Taylor degrees}}}\frac{(1/2)_{i-1}}{i!}\,t^i\right)^3
=\frac{16}{(3/2)_{m-3}}[t^m](1-\sqrt{1-t})^3
=\frac{12(m-2)}{m!}.
$$
Here $[t^m]$ extracts the Taylor coefficient. Odd $m$ cancel between the two exponentials. The bare probe is consequently $12\sinh(e)/e^3-24(\cosh(e)-1)/e^4$. In this trace word only the two-letter matrix trace at the Moyal vertex and two remaining one-letter traces contribute. The Moyal vertex supplies $2\sinh(e/2)$ since $[p_0p_1]=1$. Dividing the negative bare probe by this factor gives the complete coefficient of $s_0$ on the three scalar trace covectors at $[p_1p_2]=1$:
$$
\frac{24\sinh(e/2)}{e^4}-\frac{12\cosh(e/2)}{e^3}.
$$
Equivalently, for odd $k\ge1$, the term of degree $e^{k-2}$ is $-6[p_1p_2]^k s_0/((k+2)2^k k!)$ in units of $\delta^2\varepsilon\pi\nu$. Even $k$ give zero.

### Double trace and genus-one single trace without recursion

Take $p_0=(1,0)$, $p_1=p_2=(0,1)$, $p_3=(-1,-2)$ and evaluate the word $(0)(123)$. Write $\Theta(e)$ for its two-loop graph coefficient divided by $N\pi\nu$. The same [path-ratio identity check](central_phase.py) gives the conditional integrand of $e^3\Theta(e)$ with the symmetric Dirichlet measure. The contribution indexed by the share $w_A$, before adding the other two cyclic choices, is
$$
\frac{w_B^2+w_C^2}{3w_A^2(1-w_A)^2}
\big(\cosh(e w_A)-1\big)\big(\cosh(e(1-w_A))-1\big).
$$
Conditioned on $w_A$, the mean of $w_B^2+w_C^2$ is $3(1-w_A)^2/4$. The marginal density of $w_A$ is $1/(2\sqrt{w_A})$. Thus
$$
\Theta(e)=\frac{3}{8e^3}\int_0^1
\frac{(\cosh(et)-1)(\cosh(e(1-t))-1)}{t^{5/2}}\,dt,
\qquad
[e^{m-3}]\Theta(e)=
\frac{3}{8\Gamma(m-1/2)}
\sum_{\substack{2\le i\le m-2\\ i\text{ even Taylor degree}}}
\frac{\Gamma(i-3/2)}{i!}
$$
for even $m\ge4$. All other coefficients vanish. In the tree shuffle the pair $(12)$ has zero determinant and joins to a cubic at determinant $-2$. The other two nonzero terms join to the cubic at determinant $1$ with opposite Moyal phases. Hence, denoting the coefficient of $s_1$ at determinant $1$ divided by $N\delta^2\varepsilon\pi\nu$ by $B(e)$, the full equation on this word is
$$
4B(2e)=2\cosh(e/2)B(e)+\Theta(e),\qquad [e^{-1}]B(e)=1/2.
$$
Division by $\sinh e$ and the double-angle identity diagonalize dilation on the resulting Taylor coefficients:
$$
B(e)=\sinh(e/2)\left(e^{-2}
+\sum_{\substack{j\ge0\\ j\text{ even Taylor degree}}}
\frac{e^j}{2^{j+2}-1}\,[x^j]\frac{\Theta(x)}{\sinh x}\right).
$$

For the genus-one single trace let $\Theta(e)$ now denote the rank-independent part of the bare $(0123)$ coefficient at the same momenta, divided by $\pi\nu$. Its coefficient at $e^r$ is given by the finite signed-ribbon sum in the preceding calculation. Each shape has 40 distinct genus-one phases with nonzero weight after identical phases are combined. The ordered three-letter trace has two equal positive cyclic contributions at determinant $1$ and one merged term at determinant $-2$. Their tree equation has the same dilation as the double trace with $\cosh(e/2)$ replaced by $\exp(e/2)$, and the genus-one pole vanishes. The rank-independent coefficient of the ordered cubic at determinant $1$, in units $\delta^2\varepsilon\pi\nu$, is exactly
$$
\frac{\exp(e/2)}2
\sum_{\substack{j\ge0\\\text{Taylor degrees}}}
\frac{e^j}{2^{j+1}-1}\,[x^j]\big(\exp(-x)\Theta(x)\big).
$$
Each $[x^j]\Theta(x)$ is a finite polynomial moment with the Dirichlet factorial weights specified above. Together with the disk and triple-trace formulas this gives [direct coefficients in all four sectors](boundary_coefficients.py). At general cubic momenta restore the factor $[p_1p_2]^{n+2}$ in degree $e^n$.
