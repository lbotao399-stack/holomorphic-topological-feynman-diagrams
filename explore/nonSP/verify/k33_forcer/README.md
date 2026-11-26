# Four-source coefficient of $K_{3,3}$ on $\mathbb C^2$

For the edge order and source bracket $\langle i,j\rangle=\det(\lambda_i,\lambda_j)$ of the [four-source calculation](../k33_positive/README.md), with $\lambda_4=\lambda_5=0$ and total momentum zero,

$$
K_{K_{3,3}}=\frac{2160\zeta(5)-1414\zeta(3)-539}{14400}\,
\langle0,2\rangle\langle1,3\rangle
\big(\langle1,3\rangle^2-\langle0,2\rangle^2\big).
$$

At $\lambda_0=(1,0),\lambda_1=(0,1),\lambda_2=(1,1)$ the bracket polynomial is $6$. The [Kirchhoff replay](../k33_positive/check.py) computes the polynomial $U^6D$ from the density $D$ and the first Symanzik polynomial $U$. Summing the numerator over the 72 graph automorphisms gives 175 monomial orbits. Divide each orbit coefficient by $6\cdot72$ to obtain the scalar coefficient above.

For a monomial with exponents $\alpha_e-1$, all $\alpha_e$ are positive integers and their sum is $24$. Delete edge $0$ and place external Euclidean momentum $Q$ with $Q^2=1$ at its ends. Write $U=t_0U_0+F_0$. In any projective chart fixing another time, the integral over $t_0$ is

$$
\int_0^\infty \frac{t_0^{\alpha_0-1}\,dt_0}{(t_0 U_0+F_0)^6}
=\frac{\Gamma(\alpha_0)\Gamma(6-\alpha_0)}{\Gamma(6)}
  U_0^{-\alpha_0}F_0^{\alpha_0-6}.
$$

The remaining three-loop propagator integral $I_{12}$, with powers $\alpha_1,\ldots,\alpha_8$ and measure $d^{12}p/\pi^6$ per loop, has the same parameter integrand with prefactor $\Gamma(6-\alpha_0)/\prod_{e=1}^8\Gamma(\alpha_e)$. Hence the projective monomial contributes $I_{12}\prod_{e=0}^8\Gamma(\alpha_e)/120$. The strict convergence inequalities for each monomial are checked in the [positive-period certificate](../k33_positive/check.py).

Edges $1,\ldots,5$ form a spanning tree of the cut graph. The remaining three momenta $p_6,p_7,p_8$ are independent. For the same propagator powers let $I_d[\det\operatorname{Gram}(Q,p_6,p_7,p_8)]$ denote the integral with a Gram determinant in the numerator. At fixed Schwinger parameters, the independent projections perpendicular to $Q$ are Gaussian with loop covariance $M^{-1}/2$, where $M$ is the cut graph's loop matrix. Cauchy–Binet gives $\binom{d-1}{3}$ minors; the mean square of each $3\times3$ minor is $3!/(8\det M)$. Since increasing $d$ by two multiplies the normalized Gaussian integral by $1/\det M$,

$$
I_{d+2}=\frac{8\,I_d[\det\operatorname{Gram}(Q,p_6,p_7,p_8)]}
{(d-1)(d-2)(d-3)}.
$$

The Gaussian identity continues meromorphically in $d$ by the Schwinger formula. Apply it four times from $12-2\epsilon$ to $4-2\epsilon$. At each step FORCER reduces ordinary integrals and integrals with the Gram numerator to the same five masters. The ordinary graph basis omits internal edge sets $\varnothing,\{5,8\},\{5,6\},\{1,3\},\{1,5\}$ from the graph with unit powers. The rational basis matrix is triangular with nonzero diagonal. Its inverse expresses the physical dimension shift entirely in ordinary graph integrals, cancelling the common master normalization at each step. At $4-2\epsilon$ the finite limit is unchanged by the common MS-bar factor, whose value at zero is one. [FORCER, Ruijl–Ueda–Vermaseren](https://arxiv.org/abs/1704.06650)

The [executable replay](replay.py) uses the lexicographically least and greatest exponent vectors in each orbit. The two dimensional continuations have different linear terms in $\epsilon$. Both give zero for every negative Laurent coefficient and the same constant displayed above.

```sh
python3 replay.py --form /path/to/form --forcer /path/to/forcer
```

The FORCER directory contains `forcer.h` and its `forcer/` subdirectory. The graph calculation uses C++17, NumPy, NetworkX and SymPy. Installing `python-flint` enables SymPy's fast univariate exact arithmetic.

## Arbitrary conserved sources

At zero displacement the kernel for six arbitrary conserved sources is the linear combination in the table. Each row denotes the **unnormalized** sum of its bracket monomial over all 72 graph automorphisms, weighted by the sign of the induced permutation of the nine edges. Write $[ij]=\langle i,j\rangle$ in this table.

| Bracket monomial | Coefficient of its signed orbit sum |
|---|---|
| $[02]^3[13]$ | $193/138240+\zeta(3)/315-5\zeta(5)/1008$ |
| $[01]^2[02][03]$ | $11189/806400+3583\zeta(3)/100800-23\zeta(5)/420$ |
| $[01]^2[02][23]$ | $-1739/96768-1007\zeta(3)/20160+19\zeta(5)/252$ |
| $[01]^2[02][34]$ | $-637/13824-281\zeta(3)/2240+4\zeta(5)/21$ |
| $[01][02][12][34]$ | $1/3456$ |

**Completeness.** Source conservation leaves the five-dimensional quotient of the six-vertex permutation representation by its constant line. The determinant of the eight current rows transforms by $\det(g)^4$ under $g\in GL_2$. The Cauchy decomposition therefore selects the Schur representation of partition $(4,4)$ on the five-dimensional vertex-label space. On a vertex permutation of cycle lengths $\ell$, let $h_k$ be the coefficient of $t^k$ in $(1-t)\prod_{\ell}(1-t^\ell)^{-1}$, with one factor per cycle. The Jacobi–Trudi formula gives character $h_4^2-h_3h_5$. Its sum weighted by the edge-permutation sign is $360$, so the dimension is $360/72=5$. The five displayed signed orbit sums evaluated on the first five momentum assignments in the replay give

$$
\det\begin{pmatrix}
-24&-6&6&0&0\\
48&44&-4&0&-32\\
-24&-16&2&2&-12\\
96&-24&-24&0&0\\
-48&60&20&8&0
\end{pmatrix}=2654208.
$$

The sixth assignment in vertex order is $((0,-2),(-1,0),(1,-1),(1,1),(1,1),(-2,1))$. Its evaluation row is the linear combination of those five rows with coefficients $(6,-3/2,3,1/4,-3)$. Direct integration gives $-123421/33600-6417\zeta(3)/700+1492\zeta(5)/105$, equal to that combination of the five integrated probes. Each interpolation tests the uniform-time value against a separate grounded Kirchhoff solve and verifies strict convergence for every monomial orbit. [General-source replay](all_sources.py) · [Seven momenta checked in both cuts](all_sources.txt)

```sh
python3 all_sources.py --form /path/to/form --forcer /path/to/forcer
```

**Irrational period and higher dimensions.** For the conserved integer momenta in vertex order

$$
\big((-3,-1),(-2,2),(0,-1),(1,-1),(3,1),(1,0)\big),
$$

the five signed bracket orbit sums are $(432,-1800,540,-720,-1152)$. Substitution in the table gives

$$
K_{K_{3,3}}=-\frac{789}{640}+\frac{59}{80}\zeta(3)\ \notin\mathbb Q
$$

by Apéry's theorem. The coefficient of $\zeta(5)$ cancels. This resolves the rationality question in E6.1 for $K_{3,3}$. The four-loop cotree subdivision identity expresses this period as a finite signed sum of child periods with integral sources, increasing the holomorphic dimension by one. An irrational summand exists at each iteration. Thus every $H\ge2$ admits an admissible four-loop subdivision of $K_{3,3}$ on $\mathbb C^H$ with an irrational zero-displacement kernel at integral momenta. [Signed subdivision formula and convergence](../../../ht-families/dimension_raise.md)

## A displacement derivative at a zero of the bare graph

Take the six conserved momenta $((1,0),(0,1),(1,0),(-1,0),(-1,1),(0,-2))$. The odd edge permutation induced by $(02)$ fixes the momenta, so $K(\lambda;0)=0$. In the oriented edge order of the project, the two components of $-\partial_{z_e}K(\lambda;0)$ on edges $01$ and $03$ have the following coefficients.

| Edge and component | $1$ | $\zeta(3)$ | $\zeta(5)$ |
|---|---:|---:|---:|
| $01$, first | $35039171/65318400$ | $316373/226800$ | $-32/15$ |
| $01$, second | $1820521/10886400$ | $60881/113400$ | $-493/630$ |
| $03$, first | $16083553/65318400$ | $141/224$ | $-407/420$ |
| $03$, second | $-4695731/13063680$ | $-153887/151200$ | $107/70$ |

The linear combination of these rows with coefficients $(0,-3,4,1)$ gives $2695393/21772800-16711\zeta(3)/151200$, irrational by Apéry's theorem. The current on $05$ is minus the sum on $01,03$ by Kirchhoff conservation. Reflection gives the currents on $12,23,25$ as those on $01,-03,-05$ respectively, where a minus sign reverses orientation. The currents on $14,34,45$ vanish.

The exact replay interpolates the degree-fifteen density numerator on a downward-closed integer grid, multiplies each term by the tree-current numerator, and cuts the resulting $U^{-7}$ projective integral in fourteen dimensions. The least and greatest representatives of each time orbit give the same finite part, with every Laurent pole canceled. [Current-insertion evaluator and replay](../prism_forcer/insertions.py)
