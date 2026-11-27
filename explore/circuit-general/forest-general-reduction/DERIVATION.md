# Source-complete ear bundles on arbitrary cores

The series and parallel circuit identities and the source transfer are
already present in [the manuscript](../../../sources/LLC_main.tex) and
[`htk` conventions, Sections 2.3 and 5](../../../code/htk/CONVENTIONS.md).
The contribution here is an exact local implementation that retains a
non-series-parallel remainder, removes **sourced** degree-two vertices,
supports all original current words, and batches an arbitrary number of
ears through a fixed two-dimensional source moment calculation.

## 1. One sourced ear

Let a degree-two vertex `v` join `u` to `w`, with an existing direct edge
`u--w`. Orient the local currents as `u→v`, `v→w`, `u→w`. Write the
retained resistance as `T` and its current as `J`. The exact positive
coordinates are

\[
t_{uv}=\frac{rT}{z},\qquad
t_{vw}=\frac{(1-r)T}{z},\qquad
t_{uw}=\frac{T}{1-z},\qquad 0<r,z<1.
\]

The complete retained source rows are

\[
\lambda_u^a\longmapsto\lambda_u^a+(1-r)\lambda_v^a,
\qquad
\lambda_w^a\longmapsto\lambda_w^a+r\lambda_v^a.
\]

All other rows are unchanged. Every original current is retained:

\[
y_{uv}^a=zJ^a-(1-r)\lambda_v^a,\qquad
y_{vw}^a=zJ^a+r\lambda_v^a,\qquad
y_{uw}^a=(1-z)J^a.
\]

These identities follow by solving the two series resistors and then
merging their branch with the direct resistor. In particular
`y_vw-y_uv=lambda_v`; no zero-source hypothesis is used.

Choose a reduced spanning tree containing `u--w`. At fixed reduced
cotree currents, the source transfer changes `J` by `-lambda_v dr`.
The new original cotree current is `(1-z)J`. After wedging its complete
area form with every reduced cotree area form, terms containing a repeated
reduced-current differential vanish. The surviving local factor is exactly

\[
(1-z)\langle\lambda_v,J\rangle\,dr\wedge dz.
\]

The time Jacobian in the local edge order `(uv,vw,uw)` is

\[
\det\frac{\partial(t_{uv},t_{vw},t_{uw})}{\partial(r,z,T)}
=-\frac{T^2}{z^3(1-z)^2}.
\]

The ratio of the `htk` current-volume signs at `n` and `n-1` loops is
`(-1)^(n-1)` by conventions Section 2.3. Combining it with this negative
Jacobian gives `(-1)^n` for the positive ratio measure. For arbitrary
original edge order, also multiply by the parity of the permutation that
places `(uv,vw,uw)` before the other edges. Reversing an edge leaves the
bare H=2 kernel unchanged and reverses that current in every inserted word.

This proves the density identity used by `EarBundle.check`, including its
sign, for arbitrary retained sources and arbitrary retained topology.

## 2. Any number of sourced ears on one edge

Let `v_1,...,v_k` be degree-two vertices with the same neighbours `u,w`.
Keep their arbitrary source spinors `lambda_v_i`. Use one interval
coordinate `r_i` for each series branch, and splitting weights
`w_0,w_1,...,w_k>0`, with their sum equal to one. The direct edge has
weight `w_0`. Then

\[
\begin{aligned}
t_{uw}&=T/w_0, &
t_{uv_i}&=r_iT/w_i, &
t_{v_iw}&=(1-r_i)T/w_i,\\
y_{uw}^a&=w_0J^a, &
y_{uv_i}^a&=w_iJ^a-(1-r_i)\lambda_{v_i}^a, &
y_{v_iw}^a&=w_iJ^a+r_i\lambda_{v_i}^a.
\end{aligned}
\]

The full transferred sources are

\[
\lambda_u^a+\sum_i(1-r_i)\lambda_{v_i}^a,
\qquad
\lambda_w^a+\sum_i r_i\lambda_{v_i}^a.
\]

Apply the one-ear identity successively. Its intermediate current on the
merged direct edge scales at each later parallel merge. If `z_i` denotes
the successive splitting coordinates, the resulting density has factor
`product_i (1-z_i)^i`. The simplex coordinates are

\[
w_i=z_i\prod_{j>i}(1-z_j),\qquad
w_0=\prod_i(1-z_i),\qquad
\left|\frac{\partial(w_1,\ldots,w_k)}{\partial(z_1,\ldots,z_k)}\right|
=\prod_{j=2}^k(1-z_j)^{j-1}.
\]

Consequently the complete local density, against the reduced kernel at
the transferred sources, is

\[
(-1)^{kn-k(k-1)/2}\,
\operatorname{sgn}(\text{edge permutation})\,
w_0\prod_{i=1}^{k}\langle\lambda_{v_i},J\rangle.
\]

Here the permutation places `(uw,uv_1,v_1w,...,uv_k,v_kw)` before the
remaining edges. The reduced graph places `u→w` first and keeps all other
edge directions and their inherited order. The positive affine Jacobian is

\[
\frac{T^{2k}}{w_0^2\prod_{i=1}^k w_i^3}.
\]

Both current and density identities are checked directly against the
original Kirchhoff solve and determinant. If any removed source spinor is
zero, the displayed density is identically zero, including after arbitrary
finite current insertions.

## 3. Every original current word

Expand only the requested current powers using the affine current maps.
For each selected term, write the powers of splitting weights as
`j_0,...,j_k`. Their complete integral, including the density factor `w_0`,
is

\[
\int_{\Delta^k}w_0^{j_0+1}\prod_{i=1}^kw_i^{j_i}\,d^kw
=\frac{(j_0+1)!\prod_{i=1}^k j_i!}
       {(k+1+\sum_{i=0}^k j_i)!}.
\]

The source offsets on removed currents contribute only finite factors
`r_i^alpha_i(1-r_i)^beta_i`. Thus the remainder is a polynomial in the
retained currents, multiplied by the full source-transfer functional.
The two ends of a removed edge are never identified in the current word:
their different source offsets are both retained.

Let the retained core have `m=n-k` loops and let the original current word
have degree `q`. At fixed positive retained times, its kernel is polynomial
of degree `2m` in all source spinors; the Kirchhoff denominator depends only
on those retained times. The particular transferred-source dependence is
smaller. Changing the aggregate source displacement between `u,w` changes
every core current by the same spinor displacement times an edge-dependent
scalar. In the projective current determinant, its varying block consists
of two copies of the same `m`-row cycle matrix, multiplied by the two
coordinates of that spinor displacement. Its rank over the rational
function field is at most `m`. Every determinant term using more than `m`
varying columns vanishes, so the complete kernel has displacement degree
at most `m`. This is an exact determinant-rank bound, not a truncation.

The added bracket product has retained-current degree `k`, and the retained
part of the word has degree at most `q`. Therefore the remaining source
polynomial has total degree at most

\[
m+k+q.
\]

This is a bound on the degree in the **entire source displacement**
`sum_i r_i lambda_v_i`, which has two spinor coordinates. It does not grow
the interpolation dimension with the number of ears. Every removed source
also remains explicitly in the bracket product and current-offset factors.

## 4. Exact source moments and finite interpolation

There is also a covariant finite differential-operator form. For a retained
kernel/current polynomial of source-displacement degree at most `m+k+q`,
start at the full source rows
`lambda_u+sum_i lambda_v_i, lambda_w`. The weighted source integral is the
product, over every removed vertex, of

\[
\sum_{d=0}^{m+k+q}
\frac{(d+\alpha_i)!\,\beta_i!}
     {d!\,(d+\alpha_i+\beta_i+1)!}
\left[
\lambda_{v_i}^{a}
\left(\frac{\partial}{\partial\lambda_w^{a}}
     -\frac{\partial}{\partial\lambda_u^{a}}\right)
\right]^d.
\]

The repeated spinor index is contracted. These derivatives are tangent
to the full source-conservation constraint. Terms whose combined order
exceeds the proved displacement degree annihilate the polynomial exactly.
The implementation evaluates this complete finite operator by rational
source samples, so it does not need a globally expanded source polynomial.

Use a basis of the span of the actual removed source spinors. If two are
independent, choose those two; every other source is reconstructed by its
two exact bracket ratios. If all are parallel, one source is a basis and
all sources are reconstructed by their exact scalar coefficients. This
is an invertible coordinate description of the source displacement on its
actual span; it removes no source information.

For one interval factor, the exponential generating coefficient of a
source-coordinate monomial of total degree `d` is its directional monomial
divided by the component factorials, multiplied by

\[
\int_0^1 r^{d+\alpha}(1-r)^\beta\,dr
=\frac{(d+\alpha)!\,\beta!}{(d+\alpha+\beta+1)!}.
\]

Multiplying these finite generating series gives every source moment
through degree `m+k+q`. Terms of greater degree cannot affect a retained
coefficient and need not be constructed in subsequent products.

The implementation transposes the exact Newton interpolation formula on
the integer total-degree grid. In two source coordinates, that formula is

\[
F(x,y)=\sum_{i+j\le m+k+q}
 (\Delta_x^i\Delta_y^jF)(0,0)
 \binom{x}{i}\binom{y}{j}.
\]

Integrate its binomial polynomials using the complete source moments just
computed. This produces rational weights on at most

\[
\frac{(m+k+q+1)(m+k+q+2)}2
\]

fixed-core source evaluations. In one source coordinate the corresponding
bound is `m+k+q+1`. Some weights cancel exactly. Grid points need not lie
inside the positive source-ratio box: interpolation is an exact identity
of source polynomials, while every retained Schwinger integration still
uses positive times.

No global polynomial in the `2k` original ratios is expanded. For every
source sample, all original source rows are redistributed exactly and the
residual graph is fixed before interpolation. A special source value never
triggers an extra node-specific graph simplification.

## 5. Aggregate before period integration

All samples have identical retained time variables. Sum their exact
Schwinger numerator polynomials first, keeping each complete denominator
power. Apply series-chain beta factors and graph automorphisms only once
to that sum. The isolated read-only period worker is then called once per
remaining denominator power.

If a finite period is unavailable, `moment` returns the complete
`ExactForestSeries` of the aggregated polynomial; `series` requests that
object explicitly. No partial finite-period answer is substituted for the
missing expression. The time, source-sample, polynomial, matrix, and
weighted-period budgets remain explicit.

Different ear bundles may be composed along a geometric plan. Each local
interface still has at most two source coordinates. The total number of
samples for many distinct interfaces can grow, and the explicit global
sample budget is enforced. Forests with no applicable sourced ear retain
their complete core-period problem. There is no claim of polynomial-time
finite zeta evaluation of every product forest.

## 6. Physical higher-superpotential examples

The fixtures start with the three-loop fork and sources

\[
\lambda_0=(2,-4),\quad\lambda_1=(1,2),\quad\lambda_2=0,
\quad\lambda_A=(5,3),\quad\lambda_B=(-8,-1).
\]

Attach `k` sourced cubic degree-two ears along the edge `1--A` and retain
their sources `(i+2,3i-3)` for `i=0,...,k-1`. The balancing source at `B`
is decreased by their full sum. Vertex `1` has contracted degree `k+3`
and one surviving field, so a superpotential of degree `k+4` supports this
assignment. Vertex `0` and each ear have one surviving cubic field;
vertex `2` is a fully contracted cubic vertex and has zero source.

The quintic and sextic bare kernels are respectively

\[
-\frac{966893}{8400}+\frac{5148}{175}\zeta(3),
\qquad
-\frac{833635741}{2016000}+\frac{262369}{1750}\zeta(3).
\]

Their nonzero zeta coefficients are retained. These are individual graph
kernels, not complete Wick/color sums. The certificate also includes the
bare kernels and complete first-current spinors through eleven loops,
where the higher vertex is supported by a degree-twelve superpotential.

## 7. A genuinely sourced higher-loop Ward lift

The earlier complete three-loop all-current-word Ward identity can be
acted on by this source-complete operator. Place the ears on its common
`1--A` edge. Source merging at each of its three surviving path-edge
contractions commutes with the affine ear transfer. Its rational boundary
is therefore the base rational spider boundary with insertion
`product_i <lambda_v_i,y_1A>`, integrated over every source ratio.

For four and five loops the new script enumerates all tree and proper
faces and integrates all three actual higher-loop quotient kernels. The
boundary is checked independently with tensor Lagrange interval rules,
without using the new two-variable source weights. At the physical base
assignment, the four-loop relation is

\[
\begin{aligned}
&\left(-\frac{966893}{8400}+\frac{5148}{175}\zeta(3)\right)
-\left(\frac{45421121}{483840}-\frac{394}{7}\zeta(3)\right)
+\left(\frac{6615901}{268800}-\frac{14998}{175}\zeta(3)\right)\\
&\hspace{4em}=-\frac{398239}{2160}.
\end{aligned}
\]

The five-loop relation is

\[
\begin{aligned}
&\left(-\frac{833635741}{2016000}+\frac{262369}{1750}\zeta(3)\right)
-\left(-\frac{90661753}{216000}+\frac{295809}{1000}\zeta(3)\right)
+\left(-\frac{5063094713}{18144000}+\frac{1021187}{7000}\zeta(3)\right)\\
&\hspace{4em}=-\frac{70717559}{259200}.
\end{aligned}
\]

Every removed ear has nonzero source. The three zeta coefficients cancel
by rational arithmetic in each identity; no independence assumption about
zeta values is needed. These are complete graph-level Ward relations.
Full operator assembly remains a separate, explicitly identified task.
