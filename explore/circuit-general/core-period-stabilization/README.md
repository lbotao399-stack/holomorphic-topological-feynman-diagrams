# Every integer pole order belongs to one finite fixed-core period span

For one fixed graph, all absolutely convergent normalized projective
integrals with rational polynomial numerator and a positive integer power
of its Kirchhoff polynomial lie in a finite-dimensional rational period
span. This covers the scalar coefficients of arbitrary finite original
current and source polynomials whenever their reduction has this fixed
denominator. All free spinor indices can be retained.

The [complete proof](DERIVATION.md) gives the missing all-order statement:
the full scalar period spans at successive integer pole orders are nested
and eventually equal. The stabilization threshold is existential. No
numerical threshold, general core master basis, or numerical period
independence is asserted.

## The boundary is part of the reduction

Every edge ordering gives a Hepp cube. A greedy cotree supplies a monomial
dividing every term of the Kirchhoff polynomial. The exact quotient has
constant term one and nonnegative integer coefficients. It stays between
one and the number of spanning trees on the entire closed cube.

Absolute convergence of the complete signed rational numerator forces
all surviving coordinate exponents to be nonnegative. The resulting
forms are regular on one fixed affine complement for each sector,
independently of the pole order or numerator degree.

For a sector, write `A = Q[x_1,...,x_m,1/U_sector]`. Here `m` is the
number of cube coordinates, and `U_sector` is its fixed denominator.
The relative cohomology in the proof is the finite-dimensional rational
vector-space quotient

\[
\frac{A}{\sum_i\partial_{x_i}\bigl(x_i(1-x_i)A\bigr)}.
\]

Every displayed primitive vanishes on the appropriate cube faces. The
denominator is a vector subspace, not a polynomial ideal. The proof uses
the full family of these relative relations; it does not establish that
the repository's currently selected graph identities generate all of
them.

Multiplying the complete numerator by the full Kirchhoff polynomial
embeds every scalar pole level into the next. Each cotree summand stays
convergent because its positive integrand is bounded by the original
integrand. Finiteness then proves eventual stabilization.

## Exact controls and two distinctions

The [independent verifier](verify.py) uses the Python standard library
and imports no production modules. It reconstructs Kirchhoff polynomials
from reduced Laplacian determinants and from explicit spanning trees,
checks every sector of its recorded graphs, and checks convergence and
pole shifts using exact integer and rational arithmetic.

The recorded controls include all 908 sectors of eight graphs, 12,700
cotree factorizations, 667 scalar weight vectors, 501 explicit divergence
witnesses, 430,388 nonnegative coordinate-exponent checks, and 2,082
convergence-preserving cotree shifts. Further controls retain a complete
signed numerator in all 720 sectors of the four-vertex complete graph,
integrate 490 cycle weights independently, and verify 169 complete
one-variable relative identities and 150 polynomial cube identities.

Two complete analytic examples prevent stronger interpretations of the
theorem. The cycle with `E` edges has no positive integer scalar weights
at pole orders below `E`, followed by the nonzero period `1/(E-1)!`.
Thus arbitrarily many consecutive equal levels do not certify the final
plateau. The one-variable sector with denominator `1+x` has two relative
classes, represented by `1` and `1/(1+x)`, with integrals `1` and `log(2)`.
The original two-edge cycle family nevertheless has rational period
span of dimension one. The larger sector cohomology can be a strict
upper bound on the original graph's period space.

```bash
python3 explore/circuit-general/core-period-stabilization/run_checks.py
```

[checks.json](checks.json) records every tested graph and pole range,
all convergent tested weight vectors, exact complete example identities,
hashes of the exhaustive sector traces, and source-file hashes. The replay
recomputes every control. `verify.py --write` explicitly regenerates the
certificate after a reviewed change. The finite checks audit the formulas
and boundary factors; they do not replace the all-order proof.

## What this adds to the existing solution

The [exact geometric cube series](../forest-cube-expansion/README.md)
already supplies convergent rational coefficients and certified tails
for every finite sourced current polynomial. The present result proves
that, for a fixed core, the normalized scalar constants needed across all
integer pole orders fit in one finite rational span. General algebraic
de Rham algorithms provide a route to computing a finite relative basis;
that substantial basis computation is not performed here.

The result does not identify the span with any finite set of zeta values,
prove closure of a selected finite relation bank, or decide cancellation
in a full Wick/color sum. The all-order specified three-zeta-span question
in [the core report](../cores/REPORT.md) and the selected-relation
completeness question in [the K5 report](../k5/REPORT.md) remain separate.
The theorem concerns one fixed normalized core and integer rational
forms; changing divisors, logarithmic regulators, noninteger powers,
dimension-dependent normalization constants, and arbitrarily varying
numerical source values require separate treatment.

The proof cites primary sources directly: Brown's fixed-graph numerator
construction, Grothendieck's algebraic de Rham comparison, and Panzer's
Hepp-sector and greedy-tree results. Their precise roles and links are
listed in [the derivation](DERIVATION.md#references).
