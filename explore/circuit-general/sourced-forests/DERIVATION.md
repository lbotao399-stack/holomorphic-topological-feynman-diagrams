# Source-complete product forests: local contractions and exact periods

The source-free single and paired star identities are theorems of
[the manuscript, Section 7](../../../sources/LLC_main.tex). Their hypotheses
require zero source at the removed star centres. They do **not** require
the retained parent to be cubic or to have zero source. This distinction
extends the existing cubic-tree implementation to higher-superpotential
forests without imposing a cubic source pattern on the graph.

The new computation has three parts: integrate any permitted local stars,
evaluate the remaining complete current polynomial using an exact sourced
core, and retain an explicit convergent series if the finite period library
does not close that core. The old identities and period algorithms are
credited below; the general source-preserving planner and denominator
cancellation connect them into a complete computational route.

The implementation accepts conserved rational source vectors and arbitrary
finite current words, subject to explicit resource budgets. All finite
fractions and zeta coefficients in this note use **htk's kernel
`K`**, including the original edge order and orientations. The manuscript's
current-area convention has a separate loop-dependent sign.

## 1. Local operators with arbitrary retained sources

Write the spinor contraction as in the manuscript,

\[
\langle y,z\rangle=y_{+\dot a}z_+{}^{\dot a}.
\]

For a single zero-source star, use the local edges

\[
(vu,vw,vx,ux,wx)\longrightarrow(uw,ux,wx).
\]

Let `y,l,h` be the retained currents on `(uw,ux,wx)`. The inverse currents
on the five original edges, in that order, are

\[
-y-(1-r_1)l,\quad y-(1-r_2)h,\quad
(1-r_1)l+(1-r_2)h,\quad r_1l,\quad r_2h.
\]

Substitute these in the requested current polynomial, multiply by
`r1*r2*<l,h>`, and integrate both intervals. In the ordinary monomial
basis, a local factor `r1^a*r2^b` integrates with weight

\[
\frac1{(a+2)(b+2)}.
\]

For paired zero-source stars use

\[
(pu,pA,pB,qu,qA,qB,AB)\longrightarrow(uA,uB,AB).
\]

For retained currents `l=y_uA`, `h=y_uB`, `y=y_AB`, their inverse map is

\[
\begin{array}{lll}
y_{pu}=-r_1l-r_2h,&y_{pA}=r_1l-r_3y,&y_{pB}=r_2h+r_3y,\\
y_{qu}=-(1-r_1)l-(1-r_2)h,&
y_{qA}=(1-r_1)l-r_4y,&
y_{qB}=(1-r_2)h+r_4y,\\
y_{AB}^{\rm old}=(1-r_3-r_4)y.
\end{array}
\]

The ratio density is

\[
(r_1-r_2)(1-r_3-r_4)\langle l,y\rangle\langle h,y\rangle,
\]

on two intervals and `r3,r4>0`, `r3+r4<1`. Thus the monomial
`r1^a*r2^b*r3^c*r4^d` has scalar weight

\[
\frac{a-b}{(a+1)(a+2)(b+1)(b+2)}
\frac{c!d!}{(c+d+3)!}.
\]

Both operators are exact on every finite current word. All sources on
retained vertices remain unchanged and independent of these ratios.
Consequently the local ratio integration can be performed before the
remaining graph integration. This remains true if the parent has degree
four or larger, and if any retained source is nonzero. The retained graph
may require a different exact method at the next step.

These statements fail for a sourced star centre: the source-transfer
coefficients in the manuscript's sourced-star appendix depend on the
remaining resistances. The implementation never substitutes the passive
formula in that case.

### Signs

Removing `d` loops changes the current-area normalization by

\[
(-1)^{n(n+1)/2-(n-d)(n-d+1)/2}.
\]

For a single move this is `(-1)^n`; for a paired move it is `-1`.
Multiply this by the permutation that brings the displayed original
local edges to the front, with all remaining edges in their original
order. The reduced graph uses the displayed three new edges first,
followed by those remaining edges. Reversing an original current edge
negates every occurrence of that current. Bare `H=2` density is unchanged
by an edge reversal. These signs are retained separately in the code.

The paired move on the five-loop example below has positive full
resistance-coordinate Jacobian and kernel factor

\[
-(r_1-r_2)(1-r_3-r_4)\langle l,y\rangle\langle h,y\rangle.
\]

An independent differentiated resistance map in
[tests_exact.py](tests_exact.py) checks all original currents and the
complete ten-dimensional density identity at four exact rational points.

## 2. One Kirchhoff denominator per antisymmetric current bracket

Use the existing circuit notation: `C` is the oriented cycle matrix,
`M=C diag(t) C^T`, and `U=det(M)`. Choose any constant particular current
`y0` satisfying the source equations and put `b=C diag(t)y0`. The exact
current is

\[
y=y_0-C^{\mathsf T}M^{-1}b.
\]

For two original edges `e,f`, the Schur-complement determinant gives

\[
\boxed{
U\,\langle y_e,y_f\rangle
=\det\begin{pmatrix}
M&b^1&b^2\\
C_e^{\mathsf T}&y_{0,e}^1&y_{0,e}^2\\
C_f^{\mathsf T}&y_{0,f}^1&y_{0,f}^2
\end{pmatrix}.}
\]

The superscripts here select the two coordinates of the existing spinor;
the determinant retains their complete antisymmetric contraction. Every
entry in its first `n` rows is linear in the resistances, and its last two
rows are constant. Its polynomial is homogeneous of degree `n`.
Therefore

\[
\langle Uy_e,Uy_f\rangle
\quad\text{is exactly divisible by }U.
\]

This is an immediate inference from the circuit equations, not a new
assumption about period values. It applies to arbitrary conserved sources.
For a subset of resistances, the coefficient of their common scaling
parameter factors through the corresponding cycle columns. The bordered
determinant consequently has degree on that subset at most the rank of
those columns, exactly as a single numerator current `Uy_e` does.

Every passive single contraction removes one loop and adds one bracket;
every paired contraction removes two loops and adds two brackets. Hence
if `d` loops are removed from an original `n`-loop graph carrying a word
of degree `k`, the remaining insertion has ordinary current degree
`k+2*d`, but its numerator has at least `d` common factors of `U`.
After exact cancellation its denominator power is bounded by

\[
(n-d)+2+(k+2d)-d=n+2+k.
\]

Thus passive elimination preserves the original denominator level.
Counting each bracket as two unrelated current denominators would
unnecessarily raise it. The backend combines the complete current
polynomial before cancelling exact factors of `U`; further cancellations
from a particular source assignment or a conservation identity are also
valid. A failed exact division leaves the polynomial unchanged.

The persistent tests construct the bordered determinants independently,
without division or the production determinant algorithm, at generic
sources on three- and five-loop graphs. They verify the complete
polynomial identity on several edge pairs.

## 3. A complete sourceful five-loop local example

Use the ordered edges

\[
(01,02,03,04,15,26,35,36,45,46,56)
\]

with marked insertions `(5,6)` and conserved sources

\[
((1,2),(3,-1),(-2,4),(0,0),(0,0),(5,3),(-7,-8)).
\]

The retained parent `0` has full degree four and a nonzero source.
Vertices `3,4` satisfy the paired-star hypotheses. Removing them leaves
the source-complete three-loop series-parallel graph with ordered
vertices `(0,1,2,5,6)`. Its generic density retains every remaining source.

The bare kernel is zero, but its full first-current spinors are nonzero:

\[
\begin{aligned}
K[y_{35}]&=\left(\frac{82364363}{2721600},
                  \frac{294587393}{5443200}\right),\\
K[y_{36}]&=\left(\frac{153917681}{5443200},
                  \frac{35278553}{5443200}\right),\\
K[y_{03}]&=\left(\frac{106215469}{1814400},
                  \frac{164932973}{2721600}\right).
\end{aligned}
\]

The `4`-branch moments are the negatives of the corresponding `3`-branch
moments, and every first moment on the untouched edges vanishes. In
particular `K[y03]=K[y35]+K[y36]`, with all other original vertex
conservation equations also satisfied exactly.

The local first-current operators in the htk convention are

\[
-\frac{l}{72}\langle l,y\rangle\langle h,y\rangle,
\qquad
+\frac{h}{72}\langle l,y\rangle\langle h,y\rangle.
\]

An additional quadratic marked-current word changes `1/72` to `1/240`
in the first expression, with the requested marked-current factor
retained. Exact integration gives, with zero-based current coordinates,

\[
\begin{aligned}
K[y_{35}^{0}(y_{56}^{1})^2]&=-\frac{2542321651}{47628000},\\
K[(y_{35}^{0})^2y_{01}^{1}]&=\frac{79343557121}{6858432000}.
\end{aligned}
\]

These values are raw current moments, before the displacement
exponential's signs and factorials.

### An all-loop higher-valence comb

Keep any retained sourceful core with the three currents `l,h,y` of the
paired identity. Nest `k` passive paired branches on that boundary and
insert the current on the deepest `A` spoke. With the fixture's original
edge order, the complete integrated boundary polynomial is

\[
\frac{2(-1)^k}
 {k!((k+1)!)^2(k+2)!(2k+1)!}
\,l^{\dot a}\langle l,y\rangle^k\langle h,y\rangle^k.
\]

This follows directly from the local ratio operator. At the next paired
step the current contributes `r1*l-r3*y`, while the existing brackets
contribute `r1^k*r2^k*(1-r3-r4)^(2*k)`. The `y` contribution vanishes
because its interval integral is antisymmetric under `r1↔r2`. The `l`
contribution has interval integral

\[
\int_0^1\!\int_0^1
 r_1^{k+1}r_2^k(r_1-r_2)\,dr_1dr_2
=\frac1{(k+1)(k+2)^2(k+3)}
\]

and simplex integral

\[
\int_{r_3,r_4>0,\ r_3+r_4<1}
 (1-r_3-r_4)^{2k+1}\,dr_3dr_4
=\frac1{(2k+2)(2k+3)}.
\]

The paired htk sign is negative. Multiplying these two factors at each
step, starting with coefficient one at `k=0`, gives the displayed
factorials. No condition is imposed on the retained core sources or
parent valence. At `k=8`, a three-loop sourceful core produces a genuine
nineteen-loop higher-valence forest; the existing generic series-parallel
core integrator can complete its remaining integral. This formula
compresses that family to a factored current polynomial, rather than
asserting that an arbitrary forest has the same shape.

## 4. Source-free moves do not close every physical forest

A counterexample already occurs at five loops. Its **ordered** edges are

\[
(01,13,02,24,05,06,35,36,45,46,56),
\]

its marked pair is `(5,6)`, and its sources are

\[
((0,0),(1,2),(3,-1),(0,0),(0,0),(5,3),(-9,-4)).
\]

The interaction tree is the path `3–1–0–2–4`. Vertex `0` has degree four
and zero source, as a fully contracted quartic interaction may require.
Vertices `1,2` have degree two and carry the two surviving sources.
Vertices `3,4` are source-free cubic leaves. Neither leaf has two edges
among its neighbours, so no single move applies; their neighbour sets
are different, so no paired move applies. Series removal of `1,2` is
valid but transfers sources to the remaining stars.

Suppressing those two series vertices gives `K5` with one edge deleted.
Its planar dual is the triangular prism. This duality and the prism
period table are existing repository results, used here with the full
physical source polynomial. The exact value is

\[
\boxed{K=
\frac{3327759589}{46656000}
+\frac{50942017}{324000}\zeta(3)
-\frac{45619}{180}\zeta(5).}
\]

[independent_duality.py](independent_duality.py) constructs the face dual
from a planar embedding, independently of the period router. It verifies
the complete 75-monomial Kirchhoff inversion identity, performs the two
series beta integrals, and contracts all 2,474 nonzero prism orbits by
direct table lookup. Every weight vector satisfies every strict scalar
convergence inequality before orbit cancellation; the smallest margin is
one. The result does not use a floating-point fit or a stabilization rule.

Physical forests are not all planar either. The six-loop graph

\[
(01,13,02,24,05,06,36,37,46,47,56,57,67)
\]

with marked pair `(6,7)` and sources

\[
((0,0),(1,2),(3,-1),(0,0),(0,0),(0,0),(5,3),(-9,-4))
\]

is a valid higher-valence product forest and is nonplanar. Its zero-source
vertex `5` admits one single move, leaving the preceding hard five-loop
graph with a quadratic current bracket. The denominator cancellation in
Section 2 places that residual integral at level eight, rather than nine.
In the preceding fixed five-loop edge order its insertion is
`-<y05,y56>/4`. The complete six-loop value is

\[
\boxed{K=
-\frac{2553462100217}{369515520000}
-\frac{107288658589}{1710720000}\zeta(3)
+\frac{7666393}{95040}\zeta(5).}
\]

[independent_six.py](independent_six.py) constructs the inserted numerator
directly from the bordered determinant, without dividing by `U` and
without calling either production contraction or core integration. Its
287,184 numerator terms give 7,784 nonzero prism orbits; the independent
face-dual table replay gives the displayed three coefficients and verifies
every strict convergence inequality, again with smallest margin one.

This example proves that a source-free planner can reduce the graph while
retaining a genuine sourced core, and that the remaining complete period
can still be evaluated. No claim of universal closure under local passive
moves is made.

## 5. Polynomiality and separate convergence of every scalar term

The established numerator and support theorem is proved in
[the series-fibre derivation](../series-fibres/DERIVATION.md). For a
connected admissible `H=2` graph with `n` loops and a word of `k` current
factors, let the polynomial numerator have denominator `U^(n+2+k)`.
Its homogeneous degree is `n*(n+k)-1`, and on every nonempty edge
subset its degree is at most `(n+k)*rank(C_subset)-1`.

For a nonempty proper edge subset `T`, let `h(T)` denote the number of
independent cycles in its edge-induced subgraph, counting isolated
vertices consistently. The complementary cycle-column rank is
`n-h(T)`. Complementing the support inequality therefore shows that
every numerator monomial has degree on `T` at least `(n+k)*h(T)`.
If its positive Schwinger powers are `alpha_e=exponent_e+1`, then

\[
\sum_{e\in T}\alpha_e-(n+2+k)h(T)
\ge |T|-2h(T)>0.
\]

The last strict inequality is the admissibility condition, applied to
each connected component of the subgraph. Thus **each monomial period
converges separately**. No cancellation of divergent integrals is needed.

Each bordered bracket numerator has the same subset-degree bound as one
numerator current. Replacing two current denominators by its single
denominator consequently preserves this argument. Exact division by
further common factors of `U` lowers numerator subset degree by the
corresponding rank and lowers denominator power by one, preserving the
strict inequalities. Grouping a whole series chain adds its Schwinger
powers and contributes its exact beta factor; the lifted edge subset has
the same cycle number, so convergence also descends to the core.

## 6. An exact series when finite period reduction is unavailable

This part reuses and exposes
[the existing Hepp-series algorithm](../../nonSP/verify/hepp_series/period.py).
It gives a convergent rational series with a rational remainder bound.
It is a complete exact representation of the period, **not a claim that
the period has been reduced to a finite list of zeta values**.

The persistent tests also exercise an actually unavailable finite core,
without replacing the table lookup. The product-forest edges
`(01,12,03,04,25,26,35,36,45,46,56)` with sources
`((0,0),(1,0),(0,0),(0,1),(0,0),(0,0),(-1,-1))` and marked pair `(5,6)`
have no permitted passive move. Their series core is `K3,3` with one added
edge at denominator power seven. The public adapter retains its non-root
source and returns `ExactForestSeries`, with the missing finite table
explicitly identified. Its 10,577 original numerator terms reduce to
7,368 core terms. A partial sum is never substituted for this exact object.

For a core graph write its Kirchhoff polynomial as the sum of its `tau`
spanning-tree-complement monomials. Let `Umax` be the largest of those
monomials at the given positive resistances. Then

\[
1\le U/U_{\max}\le\tau.
\]

Every chamber with an ordering of edge resistances has a greedy dominant
cotree. For positive powers `alpha` and a denominator exponent `q`, the
integral of one monomial divided by `Umax^q` over that chamber is

\[
\prod_{j=1}^{|E|-1}
\frac1{\sum_{e\in S_j}\alpha_e-qh(S_j)},
\]

where `S_j` is the set of its `j` smallest edges. This follows by setting
each smaller resistance equal to a ratio times the next one and fixing
the largest resistance to one. Every denominator is strictly positive
by Section 5. Summing all chambers is a finite subset recursion: at a
proper subset sum over all choices of its largest edge and divide by the
displayed denominator; at the full edge set perform the sum without
division. This replaces explicit enumeration of all edge orders.

For a numerator `P`, let `mu_j` be the resulting exact rational integral
of `P*U^j/Umax^(q+j)`. Multiplication by a cotree monomial contributes at
least `h(T)` powers to every subset `T`; hence all these moments continue
to converge separately. Also let `Habs` be the rational zeroth moment
with each coefficient of `P` replaced by its absolute value.

Set `c=(tau+1)/2` and `rho=(tau-1)/(tau+1)`. Uniformly on the whole
positive projective domain, `|1-U/(c*Umax)| <= rho < 1`. The original
period equals the explicit series

\[
\sum_{m=0}^{\infty}
c^{-q}\binom{q+m-1}{m}
\sum_{j=0}^{m}(-1)^j\binom mj c^{-j}\mu_j.
\]

After order `N`, an absolute remainder bound is the rational number

\[
H_{\rm abs}
\left[1-c^{-q}\sum_{m=0}^{N}\binom{q+m-1}{m}\rho^m\right].
\]

This follows from the negative-binomial series and `1-rho=1/c`. The
bound decreases to zero. A finite sum of numerator components is treated
by adding their terms and their absolute bounds.

### The available Gegenbauer acceleration

The existing implementation also permits a rational number `s` with
`s^2>=tau`, sets `rho=(s-1)/(s+1)` and `c=((s+1)/2)^2`, and writes

\[
U/U_{\max}=c(1-2\rho z+\rho^2),\qquad -1\le z\le1.
\]

Expand `(1-2*rho*z+rho^2)^(-q)` by its polynomial generating function.
Its degree-`m` coefficient has absolute value at most
`binom(2*q+m-1,m)` on this interval. To see the bound, put `z=cos(theta)`
and factor the generating function into
`(1-rho*exp(i*theta))^(-q)*(1-rho*exp(-i*theta))^(-q)`; the triangle
inequality on its positive binomial coefficient convolution gives the
stated bound. All actual calculations use rational polynomial
coefficients, with no trigonometric or floating-point evaluation.

The resulting absolute tail is

\[
H_{\rm abs}
\left[1-(1-\rho)^{2q}
\sum_{m=0}^{N}\binom{2q+m-1}{m}\rho^m\right].
\]

The exact series is valid at every finite requested current order.
Computing its terms can be expensive: the subset recursion has `2^|E|`
states, coefficient polynomials grow, and rational bit lengths grow.
Explicit budgets fail with an error before a partial answer is returned
as an exact value. Increasing those budgets changes computational
resources, not the mathematical statement. Finite zeta completion of
every sourced forest remains a separate problem.
