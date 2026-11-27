# Complete bulk kernels with one topological direction

This directory contains new exact scalar tables, complete source and
current polynomials, and a general bulk evaluator. No article credit is
claimed. The setting is zero topological momentum; defect insertions and
inserted topological currents are outside the calculation.

## Four-loop cores: every source and every first current

For both the triangular prism and `K33`, the holomorphic dimension is one
and the topological dimension is one. Five independent original vertex
sources are retained; the sixth is minus their sum. The bare polynomial
has degree four and each of the nine first-current polynomials has degree
five. Thus each complete answer contains

\[
\binom{8}{4}+9\binom{9}{4}=1204
\]

original monomial slots, for each scalar constant. No source is set to
zero in the resulting polynomials. The files
`prism_complete_polynomials.json` and `K33_complete_polynomials.json`
contain all ten outputs, with every nonzero coefficient written exactly.
An absent monomial is zero. Their scalar constants are \(\pi^2,\pi^4\).

`source_covariance.py` imposes all graph automorphisms, with their edge
order and current orientation signs, and all five independent original
Kirchhoff identities. Exact elimination leaves 49 and 5 coefficients,
respectively, from the complete 1,204-slot arrays. Fourteen prism source
assignments and three `K33` source assignments determine those coefficients
uniquely. The full rational interpolation matrices, ranks, and all unused
equation residuals are retained. This is finite polynomial reconstruction
with an exact degree bound, not numerical fitting of a value space.

`evaluate_sources.py` constructs the complete original-edge Pfaffian
density at every assignment, integrates all scalar orbit coefficients,
and reconstructs every original source/current coefficient. Independent
grounded-Laplacian controls check the density before integration. An
additional source assignment is calculated using the distinct Gaussian
linear-minor algorithm, including all nine original currents.

The bare `K33` polynomial is identically zero, but its current polynomials
are nonzero. For the ordered edges

```text
(0,3),(0,4),(0,5),(1,3),(1,4),(1,5),(2,3),(2,4),(2,5)
```

and sources \((1,-1,0,1,0,-1)\), one exact result is

\[
K_{K_{3,3}}=0,\qquad
\int D\,y_{03}
=\frac{2243\pi^2}{4725}-\frac{121\pi^4}{2520}.
\]

Here \(D\) is the complete Schwinger density, and \(y_{03}\) is the
original oriented edge current. With the repository's
\(\exp(-\sum_e y_ez_e)\) convention, the first derivative with respect to
\(z_{03}\) is minus this current integral.

## New scalar tables

| Core | Pole order | All convergent orbits | Original exact relation rows |
|---|---:|---:|---:|
| `K33` | \(7/2\) | 16 | 125 |
| `K33` | \(9/2\) | 154 | 1,415 |
| `K33` | \(11/2\) | 936 | 9,758 |
| prism | \(11/2\) | 3,948 | 43,776 |
| Total | | **5,054** | **55,074** |

The existing exact integration-by-parts system generates these equations.
Every original row is retained in deterministic gzip files; the JSON
tables bind both compressed and uncompressed SHA-256 hashes. Verification
substitutes the complete rational answers into all 55,074 rows. A separate
sparse rank calculation modulo 2147483647 records a full-rank minor for
each system, proving uniqueness over the rationals. The values themselves
are rational dictionaries, never reconstructed from modular values.

There are also **785 independent dimension-shift checks**: multiplication
by the Kirchhoff polynomial compares successive tables, with 81 terms per
`K33` identity and 75 terms per prism identity. Full source covariance,
all source conservation identities, every Gaussian-minor holdout output,
and three older prism values with different vertex/edge orders are checked
by `verify_complete.py`.

`verify_source_convergence.py` closes the analytic step before using
source interpolation. It rederives every original bare numerator and
tests each nonzero monomial against all 510 proper nonempty edge subsets,
before any orbit cancellation. The prism source points and their graph
relabelings give 162 points and exact degree-four rank 70. For `K33`, two
additional source assignments give 252 points and the same full rank 70.
Thus every formal bare source coefficient is a rational linear combination
of absolutely convergent integrands. No evaluation of the extra integrals
is needed for this convergence argument.

Every original current numerator is supported on the positive Kirchhoff
cotree monomials. The source assignments also have full linear rank five,
so this support statement holds coefficientwise in all five sources.
Dividing by the Kirchhoff polynomial expresses each current coefficient as
a bounded linear combination of positive weights whose sum is one.
Multiplication by any finite word of these bounded coefficients preserves
absolute convergence of the bare coefficients. This proves convergence
for every finite holomorphic-current word at arbitrary sources. The
`K33` degree-five evaluation rank is recorded as 124 of 126; the
bounded-current argument is what closes that case.

`cached_half_tables.py` adds bounded memoization at fixed pole order and
Laurent precision to the unchanged scalar recurrence. Each completed
high-pole table records cold-cache comparisons. A change in Laurent
precision clears the gamma cache. No equation or pole prescription is
changed.

## One evaluator in arbitrary holomorphic dimension

The entry point is `htk.bulk_kernels.evaluate_bulk`. It accepts a connected
tight admissible bulk graph, exact rational source vectors, and any finite
word of original holomorphic currents. It verifies all local Laman
inequalities after exact bivalent-chain contraction. For at most five
loops there are at most twelve nontrivial chains, independent of the
holomorphic dimension; all their subsets are tested.

Its routes are exact bulk vanishing, stored complete source polynomials,
series-parallel integration, a complete numerator with scalar-period
reduction, and a complete convergent rational series if finite period data
are unavailable. The finite scalar worker reads tables without generating
or changing them. Parallel scalar edges retain their original labels:
the inherited automorphism helper only supports simple graphs.

```python
from htk.graph import Graph
from htk.bulk_kernels import evaluate_bulk

graph = Graph(6, [(0,3),(0,4),(0,5),(1,3),(1,4),(1,5),
                  (2,3),(2,4),(2,5)])
answer = evaluate_bulk(graph, [[1],[-1],[0],[1],[0]],
                       H=1, T=1, word=((0,0),))
assert str(answer.as_sympy()) == '-121*pi**4/2520 + 2243*pi**2/4725'
```

Run with `PYTHONPATH=code` from the repository root. Current slots are
`(original_edge_index, original_spinor_index)`, both zero based. Repeated
slots are retained with their multiplicities. A finite value has exact
`coefficients` and `as_sympy()`. Otherwise the returned `ExactBulkSeries`
has `terms(order)`, `interval(order)`, and `required_order(width)` with
exact rational remainder bounds. Explicit resource limits raise errors;
no incomplete polynomial, chart atlas, or partial sum is called a value.

`topological_polynomial.py` constructs every original current as its
Kirchhoff numerator divided by the full graph polynomial. The complete
adjugate, Kirchhoff, and voltage equations are checked symbolically.
For one topological direction it computes the full exterior Gaussian
moment. At two and four loops there is a second algorithm: add one full
holomorphic source direction and perform its quadratic or quartic
Gaussian average by exact polarization. This avoids large intermediate
Pfaffian minors in higher holomorphic dimensions. All source directions
remain present. At positive loop number, two or more topological
directions vanish in this bulk sector; one topological direction vanishes
at odd loop number. These statements retain arbitrary holomorphic-current
words and do not concern topological momenta or defect kernels.

The scalar pushforward integrates every original bivalent-chain source
monomial with its exact Dirichlet factor. It acts after the complete
source/current numerator has been formed, so it loses neither source
dependence nor correlations between inserted currents.

## Rational series for half-integer poles

`bulk_series.py` retains every dominant-tree chart. Squaring each positive
cube coordinate makes all integrable half-integer boundary powers
nonnegative integers. The complete transformed density is polynomial
divided by the normalized Kirchhoff polynomial to the original pole
order. Its denominator lies between one and the number of spanning trees
throughout the closed cube.

Choose a positive integer square whose double is at least one plus that
tree count. The binomial expansion about this square has rational
coefficients even at a half-integer pole. Its absolute geometric ratio is
strictly less than one. Replacing the pole order by its integer ceiling
bounds the omitted binomial tail by a finite rational expression. The
code integrates every polynomial coefficient by exact cube moments and
multiplies that tail by the complete absolute polynomial moment.
Consequently `interval(order)` encloses the entire original integral.
This construction does not claim finite closure in any proposed set of
zeta values.

The focused tests include a rational enclosure of \(\pi\) of width less
than \(1/10000\), all sixteen original-current chart identities for a
half-integer two-loop kernel, a nonzero three-holomorphic-dimensional
density checked against the original grounded Laplacian, and a mixed
spinor-current integral checked against independent series-parallel
fibres. They also check edge reversal, edge permutation, exact input
validation, all-source bulk vanishing, and resource-limit failures.

## Reproduction

```bash
python3 code/tests/test_bulk_kernels.py
python3 explore/circuit-general/topological-completion/verify_complete.py
python3 explore/circuit-general/topological-completion/verify_source_convergence.py
```

The second command checks the full scalar tables, all complete source
polynomials, and independent integral controls. The third checks every
original monomial's absolute convergence and the complete source-space
ranks. `independent_api_review.json` also records an actual missing-table
fallback for a generic prism with two original current insertions, and
strict rejection of invalid source vertices or extra root components.
To regenerate the source
certificates and scalar tables, see the `--help` options of
`source_covariance.py`, `evaluate_sources.py`, `build_half_tables.py`, and
`cached_half_tables.py`. `compact_tables.py` compresses all original
relation rows without discarding any equation.

Remaining scope: defect kernels, nonzero topological momentum, and
efficient finite scalar closure for every dimension and insertion order.
The general evaluator supplies a convergent exact series when finite
closure has not been established.
