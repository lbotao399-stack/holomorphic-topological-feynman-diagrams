# Explicit rational coefficients and a finite tail for every forest

This continuation uses the manuscript's full sourced current rule and the
[proved compact cube atlas](../all-product-forests/DERIVATION.md). It evaluates
every coefficient of a general geometric series by factorial-free rational
polynomial moments and gives a sufficient series order before those
coefficients are computed. The original ordered graph, all balanced sources,
and the full original current polynomial remain the input throughout.

The normalization and edge-order signs are the repository's
[H=2,T=0 conventions](../../../code/htk/CONVENTIONS.md). No Wick, color,
or interaction coefficient is inferred from an individual geometric kernel.

## 1. Complete original currents on each compact cube

There are `number of spanning trees` times `n!` cubes at loop order `n`.
The normalized Kirchhoff polynomial on each cube has a constant term one
and positive integer coefficients. Its value at the all-one corner is the
number of spanning trees. Hence, including every boundary face,

\[
1\le U\le \text{number of spanning trees}.
\]

Every original current is its complete weighted tree-flow numerator divided
by this same U. The tree-flow weights contain every original source and are
computed with the original edge orientations. The original chord currents
give `2n` scalar rows of the current-area form; this is a storage arrangement
of the full spinor, not a restriction to selected tensor components.

## 2. A smaller exact denominator

For the original chord-current rows, form the augmented polynomial matrix
whose top row is U and all its cube derivatives, and whose remaining rows
are each `U times the original current` and all its derivatives. The matrix
has size `2n+1`. Subtracting each current times the top row from its
corresponding lower row gives a first column with only the top entry U,
and multiplies every entry of the remaining current Jacobian by U.
Therefore

\[
\boxed{
\det\!\left(\frac{\partial y}{\partial\text{cube coordinates}}\right)
=
\frac{
\det\!\begin{pmatrix}
U & \partial U\\
Uy & \partial(Uy)
\end{pmatrix}}
{U^{2n+1}}.
}
\]

Here `y` contains every original chord-current spinor entry, and the
derivatives hold every **original** source fixed. The derivatives of all
resistance-dependent source transfers are consequently retained. The
existing cube orientation supplies the remaining overall sign.

This identity avoids the denominator power `4n` obtained by multiplying
`2n` separately differentiated current denominators. A complete current
word of degree q adds exactly q powers of U before cancellations. Words of
different degrees are first grouped by their degree. Every exact factor U
in a complete numerator is divided out, with the zero remainder checked.
No numerator term is discarded.

The implementation uses fraction-free polynomial elimination with exact
division checks. It subsequently adds complete numerators with identical
denominator polynomials and powers. A denominator power zero is integrated
as a polynomial immediately.

## 3. The complete rational series

In this section `centre` and `ratio` are the explicit rational numbers

\[
\text{centre}=\frac{\text{number of spanning trees}+1}{2},
\qquad
\text{ratio}=\frac{\text{number of spanning trees}-1}
{\text{number of spanning trees}+1}.
\]

The two endpoint deviations from this centre have the same absolute value.
Thus, on the entire closed cube,

\[
\left|\frac{\text{centre}-U}{\text{centre}}\right|
\le\text{ratio}<1.
\]

For a retained denominator power `power`, the binomial theorem gives

\[
\frac{\text{numerator}}{U^{\text{power}}}
=
\sum_{k=0}^{\infty}
\frac{\binom{\text{power}+k-1}{k}}
{\text{centre}^{\text{power}+k}}
\,\text{numerator}\,(\text{centre}-U)^k.
\]

Convergence is absolute and uniform on the closed cube. Every term on the
right is a **finite complete polynomial**, including its signs. Its integral
is obtained by the exact rule

\[
\int_{[0,1]^{2n}}
\prod_{j=1}^{2n}x_j^{a_j}\,\prod_{j=1}^{2n}dx_j
=\prod_{j=1}^{2n}\frac1{a_j+1}.
\]

The indices and variables in this display are the ordinary monomial
exponents and coordinates of the already constructed cube. They do not
alter any field or source index. The term at order k requires only sparse
polynomial multiplication, rational arithmetic, and these finite moments.
Summing all the retained components and the already integrated polynomial
constant gives the original forest integral exactly.

The complete answer can be stored finitely as the list of component
numerators, denominator polynomials, their positive integer powers, the
centre and ratio, and this coefficient rule. That representation contains
no unevaluated resistance integral. It does not require a finite zeta basis.

## 4. A closed rational bound for the entire remaining series

For each complete numerator, take its coefficientwise absolute value and
integrate that polynomial by the same monomial rule. This positive rational
number, called `absolute_moment` in the executable output, bounds the
integral of the absolute value of the original numerator.

After order K, the full absolute remainder is bounded by

\[
\boxed{
\text{ratio}^{K+1}
\sum_{\text{retained components}}
\text{absolute\_moment}
\sum_{j=0}^{\text{power}-1}
\frac{\binom{K+j}{j}}{\text{centre}^{j}}.
}
\]

For completeness, the finite negative-binomial tail identity used here is

\[
\frac1{\text{centre}^{\text{power}}}
\sum_{k=K+1}^{\infty}
\binom{\text{power}+k-1}{k}\,\text{ratio}^{k}
=
\text{ratio}^{K+1}
\sum_{j=0}^{\text{power}-1}
\frac{\binom{K+j}{j}}{\text{centre}^{j}},
\]

because `1-ratio=1/centre`. One proof subtracts the finite initial sum
from the full binomial identity and uses Pascal's identity successively
in the positive integer denominator power. Equivalently, both sides have
the same K=0 value and their differences under K increasing by one are
the same omitted binomial term. Their limit is zero.

The bound is rational, nonnegative, and decreasing. It is a geometric
factor times a polynomial in K of finite degree. Therefore it tends to
zero, and every positive rational requested interval width has a finite
sufficient order. `required_order(width)` finds the least order for which
twice this bound is at most the requested width, without computing the
series terms. `enclose(width)` then computes every term through that order
and returns the exact rational enclosing interval.

Intersecting the intervals from successive orders gives nested certified
intervals. The displayed `partial_sum` alone is never called the exact
period. Arbitrarily small width can require large polynomial support;
the resource limits remain explicit and a budget exception returns no
incomplete answer as a value.

## 5. Independent all-order triangle calculation

Take the ordered edges `(01,02,12)` and the balanced source spinors
`((1,0),(0,1),(-1,-1))`. The three complete cube integrands sum to

\[
-\frac3{(1+x+y)^3},\qquad 0\le x,y\le1.
\]

Here the number of spanning trees is three and the exact centre is two.
Two integrations give

\[
\int_0^1\!\int_0^1(1-x-y)^k\,dx\,dy
=\frac{1+(-1)^k}{(k+1)(k+2)}.
\]

Thus every odd series term vanishes, and every even term is

\[
\boxed{\text{term}_{2j}=-\frac3{2^{2j+3}},
\qquad\text{term}_{2j+1}=0.}
\]

The complete sum is

\[
-\frac38\sum_{j=0}^{\infty}\frac1{4^j}=-\frac12.
\]

Its actual absolute tail after order K is

\[
\frac1{2^{2\lfloor K/2\rfloor+3}}.
\]

The general algorithm independently reproduces every tested term. For
requested width `1/1000`, its sufficient order is 18 and it returns

\[
\text{partial sum}=-\frac{1048575}{2097152},
\qquad
\text{certified absolute error}=\frac{87}{262144}.
\]

The sharper actual error in this solvable control is `1/2097152`.
The bound is deliberately valid without using that separately known sum.

The same source assignment gives every complete first-current spinor:

\[
\begin{aligned}
\int y_{01,+\dot a}\,\Omega
&=-\frac16\lambda_{0,+\dot a}+\frac16\lambda_{1,+\dot a},\\
\int y_{02,+\dot a}\,\Omega
&=-\frac13\lambda_{0,+\dot a}-\frac16\lambda_{1,+\dot a},\\
\int y_{12,+\dot a}\,\Omega
&=-\frac16\lambda_{0,+\dot a}-\frac13\lambda_{1,+\dot a}.
\end{aligned}
\]

Both entries of each of these full spinors are independently enclosed by
the executable. The mixed polynomial control integrates to
`-1/2-1/6+2/24=-7/12` by the original simplex rule.

## 6. Arbitrary source spinors and finite derivative order

The current numerator and every current derivative are linear in the
original balanced sources. The full area determinant contains `2n` such
rows. Hence, at inserted-current degree q, the source dependence of every
coefficient is a homogeneous polynomial of degree `2n+q`. Source
conservation leaves one source spinor fixed by all the others.

The preceding construction therefore applies coefficient by coefficient
to this complete finite source polynomial. Uniform convergence applies to
each coefficient, or uniformly on every fixed bounded source domain.
Finite exact interpolation on rational source assignments also recovers
the same polynomial using the degree bound. No source is set to zero by
this argument. The delivered numeric interface accepts complete rational
source assignments; it does not pretend to return an already reconstructed
symbolic source polynomial when only one assignment was requested.

Every finite original-current polynomial, including terms of mixed degree
and arbitrary choices of original edges, is admitted. `spinor(edge)` returns
both entries with the original free lower spinor index represented in the
same basis as the input. Reversing an inserted edge supplies the expected
sign; permuting original edges preserves the complete orientation rule.

## 7. What has and has not been completed

This is a terminating exact coefficient formula, a finite representation
of the entire convergent answer, and a proved finite order-to-error rule
for every valid product forest at any fixed finite current order. The
proof is general; the certificates test concrete complete instances.

It is not a theorem that all forests have finite expressions in Riemann
zeta values, that all individual forests are rational, or that periods
cancel in every fully assembled interacting quantum product. Polynomial
support and the complete atlas may grow rapidly. The existing finite
rational and period reductions remain preferable whenever they close.
