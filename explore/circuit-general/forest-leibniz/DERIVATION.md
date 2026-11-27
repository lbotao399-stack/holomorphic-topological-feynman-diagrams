# Forest periods and the quantum Leibniz identity

This note explains the compatibility of the repository's forest integrals
with the manuscript's quantum product. The chain identity and its
higher-bracket proof are existing results of
[the manuscript](../../../sources/LLC_main.tex), Section
“Quantum multiplication and the Leibniz rule”, especially labels
`eq-1-1`, `b-00504`, `eq-chain`, and `eq-2-55`. No new proof of that
identity is claimed here. The exact graph example below replays existing
four-loop Ward data while retaining every source.

## 1. Which operation obeys the identity

The manuscript defines

\[
\mathbf Q=\sum_{n\ge0}\hbar^nQ_n,
\qquad
Q_n\mathcal O=
\frac1{(n+1)!}
\bigl\{\underbrace{\mathcal I,\ldots,\mathcal I}_{n+1},
\mathcal O\bigr\}.
\]

For inputs independent of `hbar`, write `[hbar^n]` for the coefficient
of its indicated power. The marked product is

\[
[\hbar^n](A\star B)=
\frac1{n!}
\bigl\{\underbrace{\mathcal I,\ldots,\mathcal I}_{n},A,B
\bigr\}^{\gamma,0},
\qquad [\hbar^0](A\star B)=AB.
\]

The marked auxiliary edge joins the ordered insertions `A,B`; it does
not supply a physical Wick contraction between their original fields.
The `n`-loop product graph has `n` interaction vertices and `2n`
physical edges, while a graph for `Q_n` has `n+1` interaction vertices
and `2n+1` physical edges. Both become `n`-loop current graphs with
`2n+1` edges after the marked edge is included in the product.

When the interaction obeys the quantum master equation, and with `|A|`
the fermion parity of a homogeneous input, the exact cochain identity is

\[
\boxed{
\mathbf Q(A\star B)
=(\mathbf QA)\star B+(-1)^{|A|}A\star(\mathbf QB).}
\]

The ordinary product already has a generally nonzero one-loop defect

\[
Q_1(AB)-(Q_1A)B-(-1)^{|A|}A(Q_1B),
\]

given explicitly in manuscript label `eq-2-29`. The marked product
cancels this defect. On cochains it need not be associative or graded
commutative; those properties hold on quantum cohomology through the
manuscript's specified homotopies. None of these statements requires
rational graph coefficients.

## 2. The exact equation at every loop order

Expanding the boxed identity at order `hbar^n` gives

\[
\begin{aligned}
0=\sum_{j=0}^{n}\Big(&
Q_j\bigl([\hbar^{n-j}](A\star B)\bigr)
-[\hbar^{n-j}]\bigl((Q_jA)\star B\bigr)\\
&-(-1)^{|A|}[\hbar^{n-j}]\bigl(A\star(Q_jB)\bigr)
\Big).
\end{aligned}
\]

Thus the equation determining the classical variation of the complete
order-`n` product is

\[
\begin{aligned}
&Q_0\bigl([\hbar^n](A\star B)\bigr)
-[\hbar^n]\bigl((Q_0A)\star B\bigr)
-(-1)^{|A|}[\hbar^n]\bigl(A\star(Q_0B)\bigr)\\
&\quad=-\sum_{j=1}^{n}\Big(
Q_j\bigl([\hbar^{n-j}](A\star B)\bigr)
-[\hbar^{n-j}]\bigl((Q_jA)\star B\bigr)\\
&\hspace{43mm}
-(-1)^{|A|}[\hbar^{n-j}]\bigl(A\star(Q_jB)\bigr)
\Big).
\end{aligned}
\]

All factorials are already included in the displayed coefficient
definitions. There is no extra `n!` or loop shift in this equation.
Each composed term has `n+1` interaction entries before their collision
is resolved. In particular, order five couples `Q_0` to the five-loop
product, `Q_1` to the four-loop product, through `Q_5` to the ordinary
product. Rationality of the five-loop supercharge alone places no
rationality requirement on the five-loop product.

For classical cocycles the two terms containing `Q_0A,Q_0B` vanish.
At one loop the equation then reads

\[
Q_0\bigl([\hbar](A\star B)\bigr)
=-Q_1(AB)+(Q_1A)B+(-1)^{|A|}A(Q_1B),
\]

which is manuscript label `eq-2-55`. At higher orders the right side
also contains variations of every lower product coefficient.

For full quantum cocycles write `A=sum_i hbar^i A_i` and
`B=sum_j hbar^j B_j`. Manuscript label `eq-2-32` gives the complete
coefficient of their product:

\[
[\hbar^r](A\star B)
=\sum_{i+j+n=r}\frac1{n!}
\bigl\{\underbrace{\mathcal I,\ldots,\mathcal I}_{n},A_i,B_j
\bigr\}^{\gamma,0}.
\]

These coefficients solve
`sum_{j=0}^r Q_j([hbar^(r-j)](A star B))=0` at every order.
Replacing the quantum inputs by only their classical leading terms
would omit required terms in this equation.

## 3. Rational supercharges leave homogeneous freedom

The fan and spider reductions give rational Taylor coefficients of
their normalized geometric kernels. This statement keeps the common
loop normalization, interaction couplings, and gauge or flavor tensors
explicit. It does not assert that arbitrary numerical choices of those
couplings are rational. The manuscript classifies the nonzero graphs
of `Q_n` as fans or spiders (`b-00247`).

Suppose at a given order all lower assembled product coefficients and
the normalized supercharge coefficients are defined over the rational
coefficient field of the chosen theory. Then the preceding equation is
a rational linear system for the next bilinear operation, with rational
forcing. At any finite charge and derivative truncation, rational row
reduction proves the following precise statement:

* If the system has a solution over a larger field of period values,
  it has a rational particular solution: set the free coordinates of
  its reduced rational system to zero.
* Every other solution is that particular solution plus a linear
  combination of rational null vectors. The coefficients of those
  null vectors may be the actual period values of the geometric
  construction.

The proof uses only rational row operations. It requires neither a
choice of independent zeta values nor any claim about their arithmetic.
In the difference of two such solutions, the complete classical
bilinear chain equation has zero right side. Hence a newly introduced
period can multiply a `Q_0`-closed bilinear map. On classical closed
inputs, its output is `Q_0`-closed. This argument says nothing about
whether that output is `Q_0`-exact.

At later orders, the lower product coefficients can themselves contain
periods. The right side then transports those same coefficients through
the rational supercharge operators. It is therefore incorrect to treat
every higher-order period term as an independent homogeneous term with
zero forcing. The first appearance and the later propagation have
different equations.

The marked-graph prescription selects a particular product among the
solutions of the chain equation. Further homotopy identities, locality,
field statistics, and physical symmetries constrain that choice. The
chain equation alone does not fix it. Manuscript label
`eq-product-contact-equivalence` makes the relevant cohomological
distinction explicit: two products induce the same cohomology product
when their difference on every closed pair is locally `Q`-exact.
Closedness alone does not establish this equivalence.

## 4. A scalar forest is one part of the assembled product

Manuscript label `eq-2-31` multiplies each complete current integral by
its ordered physical Wick contractions, sums all labelled graphs and
interaction choices, includes `1/n!`, and then applies the source
derivatives to the surviving output fields. A current insertion is
obtained from the displacement exponential with its prescribed minus
sign (`eq-product-current-insertion`). Raw moments in the repository
precede these displacement signs and factorials.

A source at a vertex represents derivatives on its surviving fields.
For a cubic product component there is exactly one such field, at its
unique full-degree-two vertex. Every other cubic vertex has zero source.
The manuscript's all-order polynomial reduction then proves rationality
for every finite current word of these cubic forests.

A component with one superpotential monomial has instead

\[
\begin{aligned}
\#\{\text{degree-two cubic vertices}\}
 &=\deg_\Gamma\bar W-2,\\
\#\{\text{fields surviving at }\bar W\}
 &=\deg\bar W-\deg_\Gamma\bar W,\\
\#\{\text{fields surviving in the component}\}
 &=\deg\bar W-2.
\end{aligned}
\]

These are manuscript label `eq-1-23`. A fully contracted quartic vertex
can therefore have zero source while two different degree-two cubic
vertices carry nonzero sources. The exact five- and six-loop examples
in [sourced-forests](../sourced-forests/DERIVATION.md) obey precisely that
pattern. Setting either source to zero would change the operator
coefficient being calculated. A higher-degree superpotential can also
leave a field at a vertex whose graph degree is three or larger.

The same counting gives, for the complete order-`n` product graph,

\[
\#(b,\beta)_{\rm out}
=\#(b,\beta)_A+\#(b,\beta)_B
-n-\#(\bar W\text{ vertices}).
\]

This controls finite loop termination on polynomial inputs. Increasing
the permitted superpotential degree does not remove the source or
antifield constraints.

An exact scalar expression with nonzero coefficients of `zeta(3)` and
`zeta(5)` establishes an evaluated graph kernel. It does not by itself
prove that its numerical value is irrational. It also does not prove
that those terms survive the complete theory-specific Wick and coupling
sum. Gauge invariance, superpotential tensor symmetries, and signs of
surviving fermionic letters can relate different graph contributions.
Conversely, the Leibniz equation does not imply that all such terms must
cancel. Determining their survival requires a complete assembled product
for a specified theory and specified input operators. That calculation
is distinct from the source-complete scalar kernel certificates.

## 5. An exact existing graph Ward cancellation

The full-source four-loop calculation in
[full_ward](../../nonSP/verify/full_ward/README.md) certifies sixty
polynomial graph Ward identities. Its rational edge-contraction matrix
has rank `2344` on `2359` independent source coefficients, leaving
fifteen homogeneous directions. Integer null vectors and a modular
rank bound certify the rational rank in `homogeneous.py`.
`verify_identities.py` compares all source coefficients with the proper
lower-loop faces. These are existing repository results, not new
operator-level calculations in this directory.

A short complete example uses the ordered sliding-graph edges

\[
(01,05,06,12,15,16,23,34,45,46)
\]

and all seven conserved sources

\[
((1,2),(-2,1),(3,-1),(-1,-2),(2,3),(-4,1),(1,-4)).
\]

Only contraction of edges `12`, `23`, and `34` gives nonzero admissible
quotients. Restoring each quotient to its recorded edge order gives
signs `-1,+1,-1`. The full three quotient values, in that order, are

\[
\begin{aligned}
&-\frac{105245}{1152}
-\frac{286719}{1120}\zeta(3)
+\frac{5373}{14}\zeta(5),\\
&\frac{61350797}{80640}
+\frac{2377057}{1120}\zeta(3)
-\frac{22570}{7}\zeta(5),\\
&\frac{11841847}{13440}
+\frac{83243}{35}\zeta(3)
-\frac{50513}{14}\zeta(5).
\end{aligned}
\]

Their signed sum is exactly the independently stored lower face
`{0,1,4,5,6}`:

\[
\begin{aligned}
\text{rational term}
&=\frac{7367150+61350797-71051082}{80640}
=-\frac{66661}{2304},\\
\text{coefficient multiplying }\zeta(3)
&=\frac{286719+2377057-2663776}{1120}=0,\\
\text{coefficient multiplying }\zeta(5)
&=\frac{-5373-45140+50513}{14}=0.
\end{aligned}
\]

[replay_existing_ward.py](replay_existing_ward.py) reconstructs these
values from the full graph source-polynomial table, computes each
contraction and source merge, and compares the result to the independent
lower-face certificate. No source is removed. The cancellation is proved
by the displayed rational arithmetic before the named constants are
evaluated. It is consequently valid regardless of any relation between
odd zeta values.

Across all sixty identities, the coefficient multiplying `zeta(5)`
vanishes in the contraction sum, while the lower faces can carry
`zeta(3)` at the same order. This is an explicit instance of homogeneous
period information alongside transported lower-period information.
The graph identities enter the higher-bracket proof of the quantum
Leibniz rule; this single geometric example is not being identified with
a complete higher-superpotential operator product.

## 6. The three-loop relation for every original-edge current word

[ward_three.py](ward_three.py) evaluates a different, three-loop sliding
graph with ordered edges

\[
(01,05,13,14,23,24,25,34).
\]

This extension uses the existing graph Stokes identity with an arbitrary
finite current polynomial. The following argument explains why the
complete original currents, including the contracted-edge currents,
obey the same identity.

### Closedness with every finite current insertion

At fixed original vertex sources, Kirchhoff conservation expresses every
edge-current differential as a constant linear combination of the
independent chord-current differentials. There are two spinor
differentials for each independent cycle. Consequently, for any finite
polynomial `F` in all original currents,

\[
\mathrm d\left(
F(y)\bigwedge_{e\in I\setminus T}(-\mathrm d^2y_e)
\right)=0.
\]

Indeed, every term of `dF` is a linear combination of differentials
already occurring in the complete chord volume, and its wedge product
with that volume vanishes. This proves closedness before any current is
integrated. It is the current-form version of differentiating the full
displacement-dependent Ward identity. A word of degree `k` obtained by
displacement differentiation has the common sign `(-1)^k` on every
face, so the same identity holds for raw current moments.

The positive-resistance spanning-tree formula expresses each current
as a convex combination of finitely many constant tree currents.
Thus every fixed finite current polynomial is bounded on Schwinger
space. Multiplying by it introduces no new singularity at the faces of
the existing graph Stokes formula. The source-complete convergence
argument in [sourced-forests](../sourced-forests/DERIVATION.md) also proves
separate convergence of the scalar monomials used to compute the
resulting face integrals. Faces whose current volume vanishes continue
to vanish after every finite insertion.

### Contracted edges retain their affine source term

For a contracted original edge `e=(a,b)` oriented from `a` to `b`, let
`inc(a,j)` be `+1` on an edge leaving `a`, `-1` on an edge entering `a`,
and zero otherwise. The limiting current is fixed by the original
Kirchhoff equation:

\[
\boxed{
y_e=\lambda_a-\sum_{j\ne e}\operatorname{inc}(a,j)y_j.}
\]

The quotient source is `lambda_a+lambda_b`; all other sources remain
unchanged. The current on the removed edge is generally nonzero. Its
source term is required when the current word differentiates that edge.

There is no additional displacement contact term. The complete phase
on this face is identically

\[
\left.e^{-\sum_j y_{j,+\dot a}z_j^{+\dot a}}\right|_{t_e=0}
=e^{-\lambda_{a,+\dot a}z_e^{+\dot a}}
e^{-\sum_{j\ne e}y_{j,+\dot a}
 (z_j^{+\dot a}-\operatorname{inc}(a,j)z_e^{+\dot a})}.
\]

Differentiating the affine translated phase produces exactly the boxed
current, including its source. Keeping only the shifted quotient
displacements would incorrectly omit the first exponential.

For the specified graph, only edge positions `0,1,6` survive. Their
signs in the inherited quotient orders are `+1,-1,+1`, and their
contracted currents are respectively

\[
y_{01}=\lambda_0-y_{05},\qquad
y_{05}=\lambda_0-y_{01},\qquad
y_{25}=\lambda_2-y_{23}-y_{24}.
\]

Each equality is used on its own quotient. The other five contractions
produce parallel edges; their current volume is identically zero,
independently of the inserted polynomial. Each nonzero quotient is a
once-subdivided complete graph on four vertices.

### The proper face retains the full shifted inner sources

The only nonzero proper tight face has vertices `{1,2,3,4}`. Its inner
edge positions, in inherited order, are `(2,3,4,5,7)`; its outer positions
are `(0,1,6)`. Call the currents on these three original outer edges
`J0,J1,J2`, respectively. The outer graph is the triangle with edges
`(0,S),(0,5),(S,5)`, where `S` is the collapsed face. Its sources are
the original sources at `0,5` and the sum of all four face sources.

The complete sources of the inner diamond, in original vertex order
`(1,2,3,4)`, are

\[
\boxed{
(\lambda_1+J_0,\ \lambda_2-J_2,\ \lambda_3,\ \lambda_4).}
\]

They sum to zero by the outer Kirchhoff equation. Original edges
`0,1,6` retain their outer currents; every other original edge retains
its complete inner current evaluated at these shifted sources. A word
containing both inner and outer currents is substituted before the
remaining outer integration.

Although the shifted inner sources depend on the outer currents, no
extra Jacobian or derivative term has been omitted. Each inner current
differential has a part along the inner ratios and a part along the
outer currents. The latter wedges to zero against the complete outer
current volume. The surviving form is exactly the inner current volume
at fixed shifted sources, followed by the outer current volume.

The edge split `(2,3,4,5,7,0,1,6)` has odd permutation parity. The inner
and outer loop numbers are two and one, so moving this face to the
right side gives

\[
-(-1)^{2\cdot1}(-1)=+1.
\]

The two-loop diamond and the outer triangle both have exact polynomial
series-parallel charts with all their sources retained. Inner current
integration therefore gives a polynomial in the shifted sources;
substituting the outer polynomial currents and integrating the triangle
gives finite factorial sums for every finite original-edge word.
At conserved rational sources this proper face is rational at every
current order. The complete sum of the three contracted core moments
equals that rational expression, with their affine current lifts and
all inherited edge signs included. This supplies a concrete family of
Ward relations in which individual core periods coexist with rational
forcing, without deleting any period term or any source.
