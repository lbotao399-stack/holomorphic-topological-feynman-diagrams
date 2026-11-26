# Four-loop prism on $\mathbb C^2$

For the ordered edges $01,02,03,14,15,23,25,34,45$, the full zero-displacement kernel is in [kernel.json](kernel.json). Each row gives four vertex pairs and the coefficients of $1,\zeta(3),\zeta(5)$. Multiply the four determinants $[ij]=\det(\lambda_i,\lambda_j)$ and sum over all twelve prism automorphisms with the sign of the induced edge permutation. Source conservation is $\lambda_0+\cdots+\lambda_5=0$.

For the six momenta $((-1,0),(1,1),(-1,-1),(0,-1),(1,0),(0,1))$, direct integration gives

$$
K_{\rm prism}(\lambda;0)=-\frac{577}{40320}-\frac{2}{35}\zeta(3)+\frac5{63}\zeta(5).
$$

The edge-signed character sum on the conserved-source Schur module of shape $(4,4)$ is $490-3(18)+3(18)-18+2(4)=480$. Division by twelve gives dimension $40$. The first forty integer probes have rank $40$ modulo $65521$. The displayed forty-first probe uses all forty interpolation values with nonzero rational coefficients. All fifteen four-mark coefficients at $0,1,2,3$ and all four signed pair projections, including the rung projection, agree with the independent star-chart and Ward calculations. [Interpolation](sources.py) · [Four-mark calculation](../four_mark_relations.md)

The Kirchhoff polynomial $U$ has 75 spanning-tree monomials. For the Schwinger density $D$, interpolate the degree-fifteen polynomial $U^6D$ on the time grid $\{1,2,3,4\}^9$ and check it against grounded Kirchhoff solves at equal and unequal times. For each monomial $t^k$ the period bank integrates $t^k/U^6$ with projective measure. Its 1,481 rows cover every prism orbit with $0\le k_e\le3$, total degree fifteen and $\sum_{e\in S}(1+k_e)>6h(S)$ for every nonempty proper edge subset $S$, where $h(S)$ is its cycle rank. Every monomial in the tested source densities lies in this list. [Exact period bank](periods.json)

Cut either edge $01$ or $02$ at $Q^2=1$. The projective integral equals the normalized three-loop propagator integral in dimension twelve times $\prod_e k_e!/5!$. For the Gram dimension shift use independent momenta $p_5,p_7,p_8$ in the $01$ cut. The ordinary master basis deletes edge sets $(6,7),(1,3),(1,4),(1,8),(5,8)$ from the eight propagators. Its triangular reductions transport the five masters through four dimension shifts. Every Laurent pole cancels and the two inequivalent cuts give identical finite parts separately for all 1,481 monomials. [Dimension shift](shift.py) · [Full replay](reconstruct.py) · [Projective cut normalization](../k33_forcer/README.md)

~~~sh
python3 evaluate.py '[[-1,0],[1,1],[-1,-1],[0,-1],[1,0],[0,1]]'
python3 reconstruct.py --form /path/to/form --forcer /path/to/forcer
~~~

## First displacement derivative

At $\lambda=((-2,0),(0,-2),(1,0),(1,0),(0,1),(0,1))$ the reflection $(23)(45)$ fixes the sources and has odd edge sign. Hence $K_{\rm prism}(\lambda;0)=0$. For edge $ij$, the two components of $-\partial_{z_{ij}}K$ at zero displacement are its integrated Kirchhoff current. The two components on edge $23$ and the first on $25$ are respectively

$$
\begin{aligned}
&-\frac{4219729}{2721600}-\frac{68831}{16200}\zeta(3)+\frac{809}{126}\zeta(5),\\
&-\frac{520883}{340200}-\frac{18917}{4536}\zeta(3)+\frac{284}{45}\zeta(5),\\
&\phantom{-}\frac{1935809}{907200}+\frac{173}{30}\zeta(3)-\frac{306}{35}\zeta(5).
\end{aligned}
$$

Contracting these three entries with $(28,20,35)$ cancels $\zeta(5)$ and gives $588751/907200-5123\zeta(3)/9450$, irrational by Apéry's theorem. Thus a derivative vertex can contribute on a six-source configuration where the graph with no insertion vanishes. Each inserted current raises the denominator to $U^7$ and the cut dimension to fourteen. Both cuts agree on all three entries. Reflection, interchange of the triangles and current conservation give the other fifteen components. The replay checks all eighteen identities on the orbit sums of the Schwinger numerators before integration. The independent [sparse interpolation evaluator](insertions.py) recovers all three period values with its required Laurent depth calculated from the actual reductions. [Current calculation](moment.py)

## Periods of higher current insertions

For $m$ current insertions the prism numerator is the original degree-fifteen polynomial times $m$ tree-current numerators, with denominator $U^{6+m}$. Each tree factor has degree at least the cycle rank on every edge subset, preserving strict convergence term by term. Cutting a monomial gives an integer-power three-loop propagator integral in dimension $12+2m$. Equivalently, before the cut integrate its first time using the beta formula. The only non-coordinate divisor before this integration is the prism Kirchhoff polynomial; changing the positive integer powers introduces no further divisor. Suppressing any subdivided paths contributes only factorial ratios. [Polynomial and convergence proof](../k4_periods/README.md#schwinger-numerator-and-convergence)

Panzer identifies the two prism cuts as $L,M$ and the $K_{3,3}$ cut as $N$. His polynomial reductions use ordinary multiple zeta values for $L,M$ and allow alternating Euler sums for $N$. They apply to the convergent rational projective integrands above because their divisors are unchanged. Thus every finite displacement coefficient of the prism at rational sources belongs to the multiple zeta value algebra, and every such coefficient of $K_{3,3}$ belongs to the alternating Euler sum algebra. At fixed real sources the displacement series converges everywhere by boundedness of the tree currents. [Panzer, Fig. 2, Thm. 1.1, §§2–3 and 6.1](https://arxiv.org/abs/1305.2161)

## Four-loop catalog on $\mathbb C^2$

In catalog order the five series–parallel graphs reduce to a single edge. Four graphs reduce to $K_4$, two to the four-spoke wheel, one to the prism and one to $K_{3,3}$. For the wheel, contract the prism rung $25$ and omit its propagator. There are eight core edges and one bivalent subdivision in either catalog graph. For a bare numerator each original exponent is at most three. For denominator $U^6$ the eight positive core propagator powers sum to 24, each is at most five, and at most one can be five. Checking the strict cycle inequalities and dividing by the eight wheel automorphisms gives 1,421 scalar integrals. Both cuts give the same finite part individually; each result lies in the rational span of $1,\zeta(3),\zeta(5)$. [Wheel bank](wheel.json) · [Reconstruction of both cuts](wheel.py)

Thus for **all thirteen graphs**, the entire bare kernel for arbitrary conserved rational sources belongs to $\mathbb Q+\mathbb Q\zeta(3)+\mathbb Q\zeta(5)$. The [catalog evaluator](catalogue.py) returns the three rational coefficients for any of the thirteen graphs without running FORM. It uses the signed bracket tables for the prism and $K_{3,3}$, the wheel bank, two-loop triangle integration by parts for $K_4$ and Dirichlet integrals for series–parallel graphs. The [thirteen graph values](catalogue.json) use the same six sources $((-3,-1),(-2,2),(0,-1),(1,-1),(3,1),(1,0))$. A separate direct reduction evaluates the entire weighted propagator sum in two cuts before taking its limit. The five series–parallel results also agree with the original ratio-coordinate engine. [Direct graph reduction](four_loop.py)

~~~sh
python3 catalogue.py
python3 catalogue.py 10 '[[-3,-1],[-2,2],[0,-1],[1,-1],[3,1],[1,0]]'
~~~

The same core classification holds in every holomorphic dimension at four loops. Suppress bivalent chains, and merge parallel edges until neither remains. A simple core of loop number four has $|E|=|V|+3$ and minimum valence three, hence at most six vertices. At five vertices the two missing edges of $K_5$ must be disjoint, giving the wheel. At six vertices the graph is cubic. If triangle-free, a vertex's three neighbors exhaust the common neighborhood of the other two vertices, giving $K_{3,3}$. If it has a triangle, its three remaining vertices must form another triangle and the intertriangle edges form a matching, giving the prism. Four or fewer vertices cannot carry seven simple edges. A parallel reduction lowers loop number; the only nonsplittable lower-loop core is $K_4$. The series and parallel beta factors retain positive integer arguments for every convergent current monomial. Therefore at any finite current order the wheel and prism give ordinary multiple zeta values and the bipartite core admits the alternating sums of Panzer's Theorem 1.1. [Three-loop reduction and convergence](../k4_periods/README.md) · [Five primitive cut types](https://arxiv.org/abs/1305.2161)

### Field generated by each graph's bare periods

The field generated over $\mathbb Q$ by the six-source values is $\mathbb Q$ for the five series–parallel graphs, $\mathbb Q(\zeta(3))$ for the four with $K_4$ core, and $\mathbb Q(\zeta(3),\zeta(5))$ for the two wheels, prism and $K_{3,3}$. The four $K_4$ coefficients of $\zeta(3)$ in the catalog fixture are all nonzero. The four remaining graphs each have a rational-source value with a nonzero $\zeta(3)$ coefficient and zero $\zeta(5)$ coefficient:

| Graph | Six momenta in vertex order | Value |
|---|---|---:|
| Wheel, catalog 8 | $(1,0),(0,1),(0,0),(1,1),(0,0),(-2,-2)$ | $-239/75600+3\zeta(3)/1400$ |
| Wheel, catalog 10 | $(1,0),(0,1),(0,0),(1,1),(-2,-2),(0,0)$ | $-3209/604800+13\zeta(3)/4200$ |
| Prism | $(1,0),(0,1),(-1,0),(0,-1),(0,0),(0,0)$ | $169/151200-3\zeta(3)/2800$ |
| $K_{3,3}$ | $(-3,-1),(-2,2),(0,-1),(1,-1),(3,1),(1,0)$ | $-789/640+59\zeta(3)/80$ |

These values recover $\zeta(3)$ by rational arithmetic. The catalog fixture also has nonzero $\zeta(5)$ coefficients $1433/45,-3323/105,-467/630$ on the two wheels and prism respectively. On $K_{3,3}$ the four-mark fixture $((1,0),(0,1),(1,1),(-2,-2),(0,0),(0,0))$ has coefficient $9/10$. Subtracting the known rational and $\zeta(3)$ terms recovers $\zeta(5)$. Thus all eight non-series–parallel topologies have explicit irrational periods, using only the irrationality of $\zeta(3)$.
