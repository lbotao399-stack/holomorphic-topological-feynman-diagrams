# Fixed-core period stabilization at every integer pole order

## 1. Statement and normalization

Fix one finite connected graph with at least two edges. A one-edge
connected graph has either no positive-weight scalar family, for a
bridge, or the constant scalar period one, for a self-loop; its
finite-span conclusion is immediate. The original ordered edge set is
denoted by `E`; `|E|` is its number of edges. Its loop number is
`|E| - |V| + 1`. The notation `U` is the Kirchhoff polynomial in the
repository's resistance convention:

\[
U(t)=\sum_{\text{spanning trees }T}\prod_{e\notin T}t_e.
\tag{1}
\]

Every coefficient in (1) is one, and its degree is the loop number. A
cotree means the complement of a spanning tree in the original edge set.
Repeated edges, when present, remain distinct labeled edges.

Use the existing scalar-period convention

\[
I_q(G;\alpha)=
\int_{\mathbb{RP}_{>0}^{|E|-1}}
\frac{\prod_{e\in E}t_e^{\alpha_e-1}}{U(t)^q}\,\Omega,
\qquad
\sum_{e\in E}\alpha_e=q\bigl(|E|-|V|+1\bigr).
\tag{2}
\]

Here `q` and all `alpha_e` are positive integers. The positive projective
domain is the set of positive Schwinger parameters modulo common positive
rescaling. The projective integration form is

\[
\Omega=\sum_{j=1}^{|E|}(-1)^{j-1}t_j\,
dt_1\wedge\cdots\wedge\widehat{dt_j}\wedge\cdots\wedge dt_{|E|}.
\tag{3}
\]

The hat means that the indicated differential is omitted. Orient the
positive domain so that its integral is the ordinary positive measure in
the gauge `sum_e t_e = 1`. This fixes all sector orientations consistently.
All convergence assertions below mean **absolute convergence**.

The more general input is a homogeneous polynomial with rational
coefficients, including a signed complete numerator:

\[
\int \frac{P(t)}{U(t)^q}\,\Omega,
\qquad
\deg P=q\bigl(|E|-|V|+1\bigr)-|E|.
\tag{4}
\]

The homogeneity condition makes (4) a projective form. A zero numerator
has zero integral. If the displayed degree is negative, the polynomial
numerator space is zero. Finite sums with several pole orders are also
allowed. No additional variable divisor, regularization, or changing
kinematical denominator is included in this statement.

**Theorem.** For this fixed graph, the rational vector space spanned by
all absolutely convergent integrals (4), at every positive integer `q`,
is finite dimensional. If the full scalar span at a given `q` contains
every convergent weight vector in (2), these spans are nested as `q`
increases and are eventually equal. The theorem supplies an existential
finite stabilization threshold. It does not supply its value, the
dimension of the numerical period span, or a chosen finite list of named
constants.

The proof below derives a fixed finite family of relative pairs. It is
compatible with Brown's fixed-graph treatment of polynomial numerators
in [Brown, Section 6.7, Lemma 6.12](#references). The sector construction
uses the greedy-tree structure of Hepp sectors described in
[Panzer, Lemma 2.8 and Proposition 2.9](#references). Neither source is
invoked to identify a numerical basis of the periods studied here.

## 2. Every sector has one fixed positive denominator

Order the edges in one permutation `sigma` and take the largest Schwinger
parameter to be one. Introduce the `|E|-1` ratios of successive parameters:

\[
0\le t_{\sigma(1)}\le\cdots\le t_{\sigma(|E|)}=1,
\qquad
t_{\sigma(i)}=\prod_{j=i}^{|E|-1}x_j,
\qquad 0\le x_j\le1.
\tag{5}
\]

For each prefix `S_j = {sigma(1), ..., sigma(j)}`, let `b_1(S_j)` be
the number of independent cycles of the subgraph on the original vertex
set with precisely those edges. Its maximum forest has
`|S_j| - b_1(S_j)` edges.

The intersection of any spanning tree with `S_j` is a forest, so

\[
|T\cap S_j|\le |S_j|-b_1(S_j),
\qquad
|(E\setminus T)\cap S_j|\ge b_1(S_j).
\tag{6}
\]

Run the greedy algorithm in the order `sigma`: accept an edge precisely
when it joins two currently distinct connected components. After each
prefix, the accepted edges form a maximal forest of that prefix. They
therefore have exactly `|S_j| - b_1(S_j)` edges. The final accepted edges
form a spanning tree, and its cotree attains equality in every inequality
(6) simultaneously.

The exponent of `x_j` in a cotree monomial after (5) is its number of edges
in `S_j`. Therefore every term of `U(t(x))` is divisible by the same
monomial, and the division is an exact polynomial identity:

\[
U(t(x))=
\left(\prod_{j=1}^{|E|-1}x_j^{b_1(S_j)}\right)
U_\sigma(x),
\tag{7}
\]

where the sector polynomial is explicitly

\[
U_\sigma(x)=
\sum_{\text{spanning trees }T}
\prod_{j=1}^{|E|-1}
x_j^{|(E\setminus T)\cap S_j|-b_1(S_j)}.
\tag{8}
\]

All exponents are nonnegative. The greedy cotree contributes a constant
term. It is the only one: the complete sequence of prefix counts
determines membership of each of the first `|E|-1` edges by differences,
and membership of the last edge follows from the fixed cotree size.
Thus the constant coefficient of `U_sigma` is exactly one. On the entire
closed cube,

\[
1\le U_\sigma(x)\le
\#\{\text{spanning trees of }G\}.
\tag{9}
\]

Equations (7)--(9) hold for every sector, independently of any numerator
or pole order. No monomial of `U` has been discarded or substituted for
the full denominator. The number of sectors is `|E|!`; equality faces
have zero interior measure, and their boundary roles will be retained.

## 3. Absolute convergence removes every coordinate pole

The positive Jacobian of (5) is

\[
\prod_{j=1}^{|E|-1}x_j^{j-1}.
\tag{10}
\]

Substituting one scalar monomial from (2), using (7) and (10), gives

\[
\frac{
\prod_{j=1}^{|E|-1}
x_j^{\sum_{e\in S_j}\alpha_e-qb_1(S_j)-1}
}{U_\sigma(x)^q}
\,dx_1\cdots dx_{|E|-1}.
\tag{11}
\]

By (9), this is absolutely integrable precisely when every coordinate
exponent is greater than `-1`. Since the exponents are integers, this
means that every exponent is nonnegative. Every nonempty proper edge
subset occurs as a prefix in some ordering. Thus the exact convergence
criterion is

\[
\boxed{
\sum_{e\in S}\alpha_e>q\,b_1(S)
\quad\text{for every }\varnothing\ne S\subsetneq E.
}
\tag{12}
\]

### Complete signed numerator

After (5), (7), and the Jacobian, a general numerator (4) produces a
finite Laurent polynomial divided by `U_sigma^q`. A Laurent polynomial
means a finite sum of integer powers, allowing negative powers at this
intermediate step. Collect identical monomials exactly before considering
their exponents.

Suppose a nonzero collected numerator has least exponent `r < 0` in one
coordinate `x_i`. It has the form

\[
x_i^r\bigl(L(x_1,\ldots,\widehat{x_i},\ldots)+x_i R(x)\bigr),
\tag{13}
\]

with `L` a nonzero Laurent polynomial in the remaining coordinates. Choose
a compact positive patch of the face `x_i = 0`, away from its other
coordinate faces, on which the absolute value of `L` is bounded below by
a positive number. Such a patch exists because a nonzero polynomial
cannot vanish on an open real set. For sufficiently small positive
`x_i`, the absolute value of (13) is bounded below by a positive constant
times `x_i^r`. Equation (9) bounds its denominator above. The integral of
`x_i^r` then diverges, since `r` is an integer at most `-1`.

It follows that every absolutely convergent complete numerator has a
polynomial, rather than Laurent-polynomial, pullback on each closed cube:

\[
\frac{P_\sigma(x)}{U_\sigma(x)^q}\,
dx_1\cdots dx_{|E|-1},
\qquad P_\sigma\in\mathbb Q[x_1,\ldots,x_{|E|-1}].
\tag{14}
\]

There is no cancellation of divergent absolute integrals in this argument.
For homogeneous numerators, there is also no concealed identification of
distinct original monomials: successive differences of their prefix
exponents recover the first `|E|-1` original exponents, and homogeneity
recovers the last. Hence every surviving original numerator monomial
satisfies (12) when the complete numerator is absolutely convergent.
Conversely, finitely many absolutely convergent monomials give an
absolutely convergent complete polynomial. This proves the equivalence
between the polynomial and scalar formulations of the theorem.

## 4. One fixed relative pair for every numerator and pole order

For the remainder of this section, `m = |E|-1` is the number of coordinates
in one fixed sector. Define its rational-function ring by

\[
A_\sigma=\mathbb Q[x_1,\ldots,x_m,U_\sigma^{-1}].
\tag{15}
\]

The notation means that only powers of this one fixed polynomial are
allowed as denominators. Let `X_sigma` be affine `m`-space with the zero
set of `U_sigma` removed. Let `B_sigma` be the union, inside `X_sigma`, of
all the faces `x_i = 0` and `x_i = 1`. These are algebraic varieties over
the rational numbers. Each nonempty intersection of distinct compatible
faces is smooth and affine. Locally, the boundary components meet as
distinct coordinate hyperplanes; this is the simple normal crossings
property. Removing the zero set of `U_sigma` preserves that property.

The real closed cube lies in `X_sigma` by (9), and its boundary lies in
`B_sigma`. Equation (14) is a regular algebraic top-degree differential
form on `X_sigma`. The relevant invariant is therefore the relative de
Rham cohomology of this fixed pair, which remembers boundary values of
primitives.

Here is a direct description of the required relative complex. In degree
`p`, take the regular differential forms on `X_sigma` whose restrictions
to every boundary component vanish. For a coordinate wedge with index
set `J`, its coefficient must be divisible by

\[
\prod_{i\notin J}x_i(1-x_i).
\tag{16}
\]

No factor for `i` is needed if `dx_i` is already in the wedge, because
that wedge restricts to zero on either face normal to `x_i`.

To verify that this complex computes the relative cohomology, augment
the restrictions to single faces by the alternating restrictions to
their double intersections, triple intersections, and so forth. The
result is a finite face complex. In local normal crossing coordinates,
it is exact in each form degree after adjoining the kernel (16). On
the topological side the same alternating complex resolves the constant
sheaf on the boundary. Apply the algebraic de Rham comparison theorem
to `X_sigma` and to every smooth affine face intersection. This
identifies the kernel complex with the relative cohomology of the pair.
The comparison theorem used here is [Grothendieck, Theorem 1](#references).

In conventional notation, the kernel complex is
`Omega^bullet(log B_sigma)(-B_sigma)`: the logarithmic differential complex
twisted to vanish on the boundary. The preceding argument specifies it
without requiring that notation as an extra hypothesis. All its terms
are locally free sheaves on an affine variety, so their higher coherent
cohomology vanishes and global sections compute its cohomology.

In top degree, (16) imposes no factors. In degree `m-1`, a term omitting
`dx_i` has a coefficient divisible by `x_i(1-x_i)`. Consequently,

\[
\boxed{
H^m_{\mathrm{dR}}(X_\sigma,B_\sigma;\mathbb Q)
=
\frac{A_\sigma}
{\displaystyle\sum_{i=1}^{m}
\partial_{x_i}\bigl(x_i(1-x_i)A_\sigma\bigr)}.
}
\tag{17}
\]

The denominator in (17) is a **rational vector subspace**, generally not
an ideal. Multiplication by an arbitrary polynomial does not preserve an
exact integration-by-parts relation without the product-rule correction.
Ordinary commutative ideal reduction is therefore not justified by (17).

Every variety in the finite face complex has finite-dimensional
cohomology. Equivalently, these algebraic pairs have finite topological
type, and comparison gives finite-dimensional de Rham cohomology.
The complex is defined over the rational numbers; extending its scalars
to the complex numbers commutes with its cohomology. Thus (17) is finite
dimensional **over the rational numbers**.

## 5. Boundary-safe relations and the finite period span

For every polynomial `P`, every coordinate, and every nonnegative integer
`q`, a directly checkable relation in (17) is

\[
\begin{aligned}
\partial_{x_i}\left(
\frac{x_i(1-x_i)P}{U_\sigma^q}\right)
={}&
\frac{(1-2x_i)P+x_i(1-x_i)\partial_{x_i}P}{U_\sigma^q}\\
&-q\frac{x_i(1-x_i)P\,\partial_{x_i}U_\sigma}{U_\sigma^{q+1}}.
\end{aligned}
\tag{18}
\]

Its integral over the full cube is exactly zero. The primitive is regular
on the closed cube, vanishes at the two faces normal to `x_i`, and its
differential wedge restricts to zero on the other faces. All terms in
(18) are individually integrable. This avoids an operation on divergent
monomials or a discarded boundary term.

Allowing every polynomial and every pole order in (18) gives exactly the
vector subspace in the denominator of (17). These are therefore complete
relations for **relative exactness** in a sector. This statement is about
the full family of relations, not a particular finite list of generated
rows. It also does not assert that every numerical identity between
period values comes from relative exactness.

Choose a rational basis of each finite-dimensional space (17), and choose
rational top-form representatives. Their integrals over the cube form a
finite spanning list of numbers for all integrals on that sector.
Integration need not be injective on cohomology, so this list need not
be numerically independent. Sum over the finitely many sectors. Equations
(5)--(14) place every original convergent numerator in this fixed direct
sum of relative spaces. This proves the fixed-core finite-span theorem.

Different sectors introduce artificial internal faces. The finite direct
sum retains those faces, and may be larger than the image of the original
graph family. No equality between these two spaces is assumed.

## 6. Nesting and the existential threshold

Multiplication by the complete polynomial `U` gives the exact identity
already compatible with the repository's scalar normalization:

\[
\boxed{
I_q(G;\alpha)=
\sum_{\text{cotrees }C}
I_{q+1}(G;\alpha+\mathbf 1_C).
}
\tag{19}
\]

The vector `1_C` adds one precisely on the original edges of that cotree.
The new weights sum to `(q+1)` times the loop number. Each new summand is
absolutely convergent separately: its positive integrand is the old
positive integrand multiplied by

\[
\frac{\prod_{e\in C}t_e}{U(t)},\qquad
0\le\frac{\prod_{e\in C}t_e}{U(t)}\le1.
\tag{20}
\]

The same conclusion follows directly from (6) and (12). Thus the full
scalar rational span at order `q` is contained in that at order `q+1`.
All these spaces lie in the one finite-dimensional period span proved
above. They therefore eventually stabilize. Explicitly, choose a finite
basis of their union from its actual scalar periods and take the maximum
of the finitely many pole orders at which those basis elements occur.
Every larger level contains that basis by (19).

This is an existence argument. It does not compute a threshold, prove a
period-independence assertion, or show that a chosen implementation's
finite rows generate the full relative relations.

### Consecutive equal levels do not certify stabilization

For the cycle with `E` edges, the loop number is one and
`U = sum_e t_e`. Positive integer weights with sum `q` do not exist when
`q < |E|`. At `q = |E|`, however,

\[
I_{|E|}(G;1,\ldots,1)=\frac1{(|E|-1)!}.
\tag{21}
\]

More generally, the gauge `sum_e t_e = 1` gives the exact simplex moment

\[
I_q(G;\alpha)
=\frac{\prod_e(\alpha_e-1)!}{(q-1)!},
\qquad \alpha_e\ge1,\quad \sum_e\alpha_e=q.
\tag{22}
\]

The scalar span is zero for every `q < |E|` and is the rational numbers
for every `q >= |E|`. This family has arbitrarily long initial plateaux.
Finite dimension bounds the number of strict increases, not the position
of the final increase. A row count or two equal consecutive ranks is not
a general stopping certificate.

## 7. An exactly solved relative example and a strict upper bound

Consider one coordinate with `U_sigma = 1+x`. The ring in (15) consists
of rational functions whose only finite possible pole is `x=-1`.
Polynomial division and partial fractions write any such function as a
polynomial plus a finite sum of powers of `(1+x)^(-1)`. Every summand
except the simple pole has a rational primitive in the same ring.

Let `a` be the simple-pole coefficient and let `R` be a rational primitive
of the other terms. Subtract its endpoint difference times the constant
function:

\[
f(x)=
\bigl(R(1)-R(0)\bigr)+\frac{a}{1+x}
+\partial_x\left[
R(x)-R(0)-x\bigl(R(1)-R(0)\bigr)\right].
\tag{23}
\]

The expression in square brackets vanishes at both endpoints. Its
denominator is a power of `1+x`, so its numerator is divisible by
`x(1-x)`. The last term is a relation of (17).

The two classes represented by `1` and `1/(1+x)` are independent in
relative cohomology. A rational derivative has zero residue at `x=-1`,
so an exact linear combination must have zero simple-pole coefficient.
The integral of the remaining constant derivative is zero because its
primitive vanishes at both endpoints, forcing that constant to be zero.
Thus this relative quotient has dimension exactly two, with integrals

\[
\int_0^1 1\,dx=1,
\qquad
\int_0^1\frac{dx}{1+x}=\log 2.
\tag{24}
\]

For example, the complete identity

\[
\frac1{(1+x)^2}
=\frac12+
\partial_x\left(\frac{x(1-x)}{2(1+x)}\right)
\tag{25}
\]

gives integral `1/2` without a boundary correction.

This sector polynomial occurs for the graph with two parallel edges,
which is the two-edge cycle. Yet (22) shows that its entire original
convergent scalar family has rational span of dimension one. The extra
relative class in (24) belongs to the larger sector space, not to an
asserted numerical master basis of that graph. This explicitly exhibits
the strict-upper-bound possibility in Section 5. No general statement
of numerical period independence is needed.

## 8. Source spinors, currents, dimensions, and computation

The theorem applies coefficientwise to a finite numerator polynomial
carrying all original free spinor indices. Every original current and
source transfer may be retained. Multiplication by a finite current
polynomial or a finite number of derivatives is covered whenever the
existing exact reduction yields the normalized rational forms (4) with
the same fixed `U` and absolute convergence. No external spinor index
has to be contracted to establish the scalar period statement.

Varying an integer holomorphic dimension is covered under that same
condition: it changes integer pole orders and numerator polynomials,
not the fixed divisor. The coefficients lie in one finite rational
period span, while the source polynomial itself may have arbitrarily
large degree. This does not turn the set of all arbitrary complex source
evaluations into a finite-dimensional rational vector space. Known
normalization factors such as dimension-dependent powers of `pi` must
be separated before making the normalized fixed-core assertion.

The proof does not cover noninteger pole orders, logarithmic regulator
coefficients, additional changing divisors, or an unbounded set of core
graphs. Nor does it prove nesting for a restricted physical family at
fixed source degree; (19) establishes nesting for the full scalar family.

There is a constructive route in principle. General algorithms compute
de Rham cohomology of affine complements and finite complexes of
algebraic differential modules; see [Oaku--Takayama and Walther](#references).
Apply such an algorithm to the finite relative face complex of Section 4.
Once a rational basis with representatives `f_j` is available, a requested
sector integrand `f` can be reduced by finding the exact rational identity

\[
f-\sum_j c_j f_j
=\sum_i\partial_{x_i}\bigl[x_i(1-x_i)g_i\bigr],
\qquad c_j\in\mathbb Q,\quad g_i\in A_\sigma.
\tag{26}
\]

Search increasing polynomial-degree and pole-order bounds for the `g_i`.
Each search is finite rational linear algebra, and a spanning cohomology
basis guarantees eventual success. The rational-function identity (26)
and the explicit boundary factors provide a directly checkable
certificate. Basis construction is the substantial preliminary step;
it is not implemented by the verifier accompanying this note.

The existing [geometric cube series](../forest-cube-expansion/README.md)
provides an independent exact representation of rational cube integrals.
For the sector polynomials here, (9) also gives the same centered
geometric-series construction mathematically. No finite evaluation of
its constants in zeta values follows from finite cohomology alone.
In particular, this theorem does not close the all-order specified
three-zeta-span question in [the core report](../cores/REPORT.md) or the
all-order selected-relation question in [the K5 report](../k5/REPORT.md).

## 9. What the executable checks establish

The local verifier uses Python's standard library only and imports no
production period or graph modules. It reconstructs Kirchhoff polynomials
independently as determinants of reduced weighted graph Laplacians and
compares them with spanning-tree complements. It checks every edge
ordering of its recorded small graphs, the exact factorization (7),
the unique constant term, convergence inequalities, and every recorded
convergence-preserving cotree shift.

It separately checks signed complete-numerator pullbacks, a divergent
numerator witness, exact cycle simplex moments and their nesting, the
complete one-variable reductions (23), and polynomial cube relations
with explicit vanishing primitives. These finite controls audit the
formulas and boundary conditions; the unbounded theorem follows from
Sections 2--6 and the cited comparison theorem, not from extrapolating
the tested graph sizes or pole orders. The JSON certificate records the
exact coverage and hashes of every replayed input and source file.

## References

1. Francis Brown, *Feynman amplitudes, coaction principle, and cosmic
   Galois group*, [arXiv:1512.06409](https://arxiv.org/abs/1512.06409),
   Section 6.7, especially Lemma 6.12. This is the fixed-graph
   polynomial-numerator and resolved-boundary reference.
2. Alexander Grothendieck, *On the de Rham cohomology of algebraic
   varieties*, Publications Mathématiques de l'IHÉS **29** (1966),
   95--103, [DOI:10.1007/BF02684807](https://numdam.org/articles/10.1007/BF02684807/),
   Theorem 1. This supplies the smooth affine algebraic de Rham comparison
   used on the finite face complex.
3. Erik Panzer, *Hepp's bound for Feynman graphs and matroids*, Annales
   de l'Institut Henri Poincaré D **10** (2023), 31--119,
   [DOI:10.4171/AIHPD/126](https://ems.press/journals/aihpd/articles/8830429),
   Lemma 2.8 and Proposition 2.9. These identify the greedy cotree on a
   Hepp sector and the positive power-counting domain.
4. Toshinori Oaku and Nobuki Takayama, *An algorithm for de Rham cohomology
   groups of the complement of an affine variety via D-module computation*,
   [arXiv:math/9801114](https://arxiv.org/abs/math/9801114), Theorem 1.1
   and Algorithm 1.2.
5. Uli Walther, *Algorithmic Computation of de Rham Cohomology of
   Complements of Complex Affine Varieties*,
   [arXiv:math/9807176](https://arxiv.org/abs/math/9807176), Theorem 4.3
   and the algorithm for complexes of differential modules.
