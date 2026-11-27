# Exact reduction for every product-forest topology

The starting point is the user's existing [manuscript](../../../sources/LLC_main.tex):
the product-forest theorem `b-00249`, the cubic integration theorem `b-00405`,
and the sourced current and full current-area rules
`local-B-sourced-current-map` and `local-B-sourced-operator`.
The contribution here makes the sourced rule terminate on every classified
forest, composes the existing integration methods, and supplies an independent
convergent integration route directly from the complete currents.

Throughout, the two insertions are the ordered vertices A,B. Original vertex
labels, edge order, edge orientations, source spinors, and free lower spinor
indices are retained. The geometric normalization is
[htk/CONVENTIONS.md](../../../code/htk/CONVENTIONS.md), sections 2 and 2.3.

## 1. Why the sourced rule terminates on every forest

Remove A,B. Each remaining connected component is a tree. Choose any root
in each tree and eliminate its descendants first. Once a vertex's descendants
have disappeared, its only possible retained neighbours are its parent,A,B.
Elimination creates edges only between those neighbours. It therefore cannot
create an edge between two distinct unprocessed interaction branches.

This proves the invariant inductively: every step has at most three retained
neighbours, and every interaction vertex is removed exactly once. The final
network is the marked edge AB. A high-valence superpotential vertex causes no
exception: when it is chosen as root, its descendants disappear before it,
leaving only A,B at the last step of that component.

Intermediate graphs need not obey the original Laman edge count. The circuit
identity remains valid on those graphs. The implementation retains the
original current-area form throughout; it does not replace such an intermediate
network by a lower-loop Laman kernel with a scalar prefactor.

## 2. Complete source transfer and reverse reconstruction

The following sums run over all currently retained neighbours of the removed
vertex v. An absent edge has zero conductance. The exact local rule is

\[
\frac1{t_{uw}^{\mathrm{red}}}
=\frac1{t_{uw}}
+\frac{1}{t_{vu}t_{vw}
\displaystyle\sum_{x\text{ adjacent to }v}1/t_{vx}},
\]

\[
\lambda^{\mathrm{red}}_{u,+\dot a}
=\lambda_{u,+\dot a}
+\frac{1/t_{vu}}{\displaystyle\sum_{x\text{ adjacent to }v}1/t_{vx}}
\lambda_{v,+\dot a}.
\]

The transfer coefficients sum to one, so the full source conservation law is
preserved at every step. All original sources are carried as exact coefficient
rows; a source may be nonzero at any vertex.

Use the repository's electric potential, defined by
\(t_{uv}y_{uv,+\dot a}=\varphi_{u,+\dot a}-\varphi_{v,+\dot a}\).
After the retained network has been solved, reverse substitution gives

\[
\varphi_{v,+\dot a}
=\frac{\lambda_{v,+\dot a}
+\displaystyle\sum_{u\text{ adjacent to }v}\varphi_{u,+\dot a}/t_{vu}}
{\displaystyle\sum_{u\text{ adjacent to }v}1/t_{vu}},
\qquad
y_{vu,+\dot a}=\frac{\varphi_{v,+\dot a}-\varphi_{u,+\dot a}}{t_{vu}}.
\]

The source in this reverse formula is the complete source accumulated at that
stage, with its dependence on the original sources and resistances retained.
The full current on every original edge follows, including edges removed many
stages earlier. All pivots are positive for positive original resistances.

The weighted Laplacian is the repository's
\(\Lambda=B^{\mathsf T}t^{-1}B\), grounded at the original last vertex.
The exact Schur-complement determinant identity yields

\[
U_\Gamma(t)
=\left(\prod_{e\in I}t_e\right)
\frac1{t_{AB}^{\mathrm{final}}}
\prod_{v\text{ eliminated}}
\left(\sum_{u\text{ retained at the removal of }v}\frac1{t_{vu}}\right).
\]

Thus the algorithm produces every source-linear current coefficient and U
without expanding a Kirchhoff polynomial or enumerating spanning trees.
There are linearly many local steps, at most three retained neighbours per
step, and quadratically many exact coefficient operations for the complete
source-linear response. These are arithmetic-operation counts; the bit sizes
of exact rational numbers are not bounded by this statement.

## 3. The complete current differential

The resistance-dependent transferred sources must be differentiated.
Freezing them would change the current-area form in the manuscript's sourced
operator identity.

The implementation avoids that mistake by differentiating the reconstructed
**original** network at fixed original sources:

\[
\Lambda\varphi_{+\dot a}=\lambda_{+\dot a},
\qquad
\frac{\partial\varphi_{+\dot a}}{\partial t_f}
=\Lambda^{-1}B_f^{\mathsf T}\frac{y_{f,+\dot a}}{t_f}.
\]

Here \(B_f\) is the original incidence row on edge f. Since
\(y_e=B_e\varphi/t_e\),

\[
\boxed{
\frac{\partial y_{e,+\dot a}}{\partial t_f}
=\left(\frac{B_e\Lambda^{-1}B_f^{\mathsf T}}{t_e}
-\delta_{ef}\right)\frac{y_{f,+\dot a}}{t_f}.
}
\]

Every shifted-source differential is included in this original-network
identity. The source-linear response already contains the required factors,
so a separate cubic-size tensor is unnecessary when only the current
differential or area density is requested.

The full original chord-current form
\(\bigwedge_{e\in I\setminus T}d^2y_e\)
then gives the H=2,T=0 density with the repository's exact orientation factor.
Changing the gauge edge leaves the projective density unchanged. The tests
also verify changes of component root, vertex labels, original edge order,
and edge orientations.

## 4. Integrate complete current polynomials before individual periods

`ForestReduction` accepts an arbitrary finite polynomial in the original
currents, including terms of different degrees. Its first operation uses the
exact tree-flow identity from the original graph:

\[
y_e=y_e^{\mathrm{tree}}+
\sum_{c\in I\setminus T}(C_T)_{ce}\,y_c.
\]

This is an equality for every original edge. Every coefficient is retained.
It reduces source-conservation insertions, including mixed-degree polynomial
identities, before any integration. For example,

\[
\int\left(
\sum_{e:e(0)=v}y_{e,+\dot a}
-\sum_{e:e(1)=v}y_{e,+\dot a}
\right)\Omega_\Gamma
=\lambda_{v,+\dot a}\int\Omega_\Gamma.
\]

The remaining complete polynomial follows the applicable exact routes:

1. Cubic forests with their physical root sources use the manuscript's
   local current-message recursion and finite Dirichlet factorials.
2. Generic series-parallel graphs retain all sources in the full density.
3. Source-free single and paired stars are integrated locally; sourced ears
   are then integrated on the retained graph, including source transfers.
4. Full source-dependent residual numerators are combined before scalar
   period lookup. A missing finite expression keeps its complete convergent
   rational series and its rational remainder bound.

The certified five-loop full-source polynomial is now usable under marked
graph isomorphisms, arbitrary original edge orientations, and the correct
edge-permutation sign. Its source-support hypotheses and all recorded
dependency fingerprints are checked before use. It applies to an insertion
which becomes constant after the **complete** Kirchhoff simplification.
Unmatched source data continue through the general exact route.

For the supplied five-loop example, both the bare kernel and the complete
current-conservation insertion at vertex 1 evaluate to

\[
\boxed{\frac{1187-816\zeta(3)}{435456000}.}
\]

This result is returned directly by the executable interface; it is no longer
an isolated research formula unavailable to the general forest solver.

## 5. A compact cube cover for direct current integration

The following construction applies to the bridgeless connected graphs under
consideration. Choose a spanning tree T and order its n complementary edges
by increasing resistance. Gauge the last complementary resistance to one.
For every tree edge, its upper resistance is the first complementary edge
in that ordering whose fundamental cycle contains it. Its resistance divided
by that upper resistance is one independent unit-interval coordinate.
The n-1 successive ratios of complementary resistances supply the remaining
coordinates. There are

\[
(|V|-1)+(n-1)=2n
\]

coordinates in the product-forest case.

The minimum-spanning-tree cycle criterion says that these inequalities are
exactly the cone in which T is the minimum spanning tree. Away from equality
walls, that tree and the complementary-edge order are unique. The cubes
therefore cover the positive projective resistance space once, apart from
measure-zero common walls. The number of cubes is

\[
\#\{\text{spanning trees}\}\,n!,
\]

instead of \((2n+1)!\) fully ordered Hepp chambers. This compares explicit
chart counts; the older Hepp-series implementation already avoids enumerating
all fully ordered chambers by a subset recursion.

Divide U by the monomial belonging to the complementary edges of T. Every
other spanning-tree monomial becomes a monomial with nonnegative exponents
in the cube coordinates, and the chosen tree contributes the constant 1.
To see nonnegativity for a successive complementary ratio, take the subgraph
formed by the corresponding complementary prefix and all tree edges whose
upper resistance belongs to that prefix. Its cycle rank is the prefix length.
Every alternative tree complement meets this subgraph in at least that many
edges. These are precisely the required nonnegative exponents. Tree-edge
fraction exponents are already zero or one.

Consequently the normalized denominator is at least 1 on the **closed** cube.
The original current on any edge is the weighted sum of its exact tree flows,
divided by this positive denominator. It is a rational function without poles
on the closed cube. Its full differential, area form, and every finite current
polynomial are therefore continuous there, with all sources unchanged.

## 6. An exact rational series without a global numerator

At dyadic level k, divide every coordinate of every cube into \(2^k\) equal
intervals. Evaluate the complete current-area integrand at every cell centre;
its coordinates are \((2j+1)/2^{k+1}\). Multiply by the exact cell volume
\(2^{-2nk}\) and sum over every cell and cube, with the original orientation.
Every operation and every resulting sum is rational for rational sources.

The exact integral is the level-zero sum plus the infinite sum of the
successive level differences. This is a specified convergent rational series,
not an unevaluated resistance integral. The implementation also encloses the
complete integrand on each closed cell using rational interval arithmetic.
The sum of all cell lower bounds and the sum of all cell upper bounds enclose
the full integral. The larger distance from the midpoint sum to these two
bounds is its certified absolute tail bound.

Uniform refinement makes the interval widths tend to zero: the normalized
denominator stays at least 1, and the finite rational expression for every
current and derivative is continuous on a compact cube. Thus every positive
rational requested width is attained after finitely many levels, provided
the explicit chart, evaluation and storage budgets permit those levels.
Budget exhaustion returns no partial result as an exact value.

This route never needs the expanded global current numerator or a finite
period table. It complements the finite rational/zeta reductions; it does
not assert that every higher-superpotential integral has a finite expression
in Riemann zeta values.

## 7. Complete checks and computed values

The independent graph census examines every simple original graph through
four loops, before using the production forest recognizer. It agrees with
the constructive unrooted-tree/spoke enumeration.

| Loops | All marked forests | Cubic forests | Higher-valence forests |
|---:|---:|---:|---:|
| 1 | 1 | 1 | 0 |
| 2 | 3 | 3 | 0 |
| 3 | 11 | 8 | 3 |
| 4 | 46 | 23 | 23 |
| 5 | 200 | 65 | 135 |
| 6 | 874 | 194 | 680 |

Every one of the 1,135 topologies through six loops passes the full
source-response comparison at two resistance assignments. All 190,326
source-response coefficients agree with the original weighted-Laplacian
calculation. The independent original cycle-form density agrees as well,
including two gauge choices. The plans contain 1,721 intermediate stages
which fail the Laman edge count; all are retained by the general circuit rule.

At the explicitly recorded admissible source assignment on **each** topology
through four loops, both implementations compute all 61 bare kernels and all
1,018 original-edge first-current entries. The results include every edge
and both entries of each free spinor. Their complete spinor expressions in
original-source bases are printed in [FULL_RESULTS.md](FULL_RESULTS.md).
All 692 integrated current-conservation identities hold exactly.

For example, the four-loop graph with ordered edges

\[
(01,03,04,12,14,25,34,35,45)
\]

and marked vertices 4,5 has the source rows explicitly recorded in
`integral_certificate.json`, forest 14. Its complete bare value is

\[
\boxed{-\frac{127957}{896}+\frac{25553}{280}\zeta(3).}
\]

Individual forest periods can therefore remain present. Their cancellation
in a fully assembled product is a separate Wick/color identity; the existing
all-jet fixed-input cancellation theorem remains applicable in its stated
field family.
