# Sourced ears on general product-forest cores

`htk.forest_ears` adds an exact reduction of **nonzero-source degree-two
vertices** whose neighbours already share an edge. The remainder may be
non-series-parallel. Any number of ears on one edge is integrated through
at most two source coordinates, and every finite word in the original
currents is retained.

The old zero-source single/paired-star planner cannot remove these
vertices. This module is additive: import it explicitly; the existing
dispatcher and package exports are unchanged.

## API

```python
from htk.forest_ears import EarBundle, ForestEarMoments, find_ear_bundles

with ForestEarMoments(graph, full_sources, marked=(A, B)) as engine:
    bare = engine.moment()
    first_component = engine.moment({(original_edge, 0): 1})
    second_component = engine.moment({(original_edge, 1): 1})
    complete_series = engine.series({(original_edge, 0): 1})
    statistics = engine.last_stats
```

`moment` returns a complete `Fraction`, `ExactPeriod`, or `ExactForestSeries`.
Set `allow_series=False` to require an available finite period expression.
`integrate(polynomial)` accepts a mixed-degree sparse current polynomial
with one exponent per original edge/component.

For a single local operation,

```python
bundle = EarBundle(graph, edge=(u, w), marked=(A, B))
bundle.check(full_sources, count=2)
result = bundle.pushforward(full_sources, powers)
```

`result.samples` contains the exact lower-loop full-source assignments and
current polynomials. Their integrals sum to the original answer. The
geometric reduction plan is fixed before source interpolation.

## What is reduced

For `k` ears and a degree-`q` original current word, a retained core with
`m` loops requires source moments only through degree `m+k+q`. A bundle
uses at most

\[
\frac{(m+k+q+1)(m+k+q+2)}2
\]

source evaluations; parallel source spinors require only `m+k+q+1`.
The complete source-integrated Schwinger polynomials are summed before
series beta factors, graph automorphisms, and period lookup. The algorithm
never expands one global polynomial in the original `2k` ratios.

The [derivation](DERIVATION.md) gives the full current/source maps, the
edge-order sign, Dirichlet factors, source-degree bound, and exact Newton
moment functional. Other ear bundles can be composed, with an explicit
global source-sample budget. Irreducible residual periods remain complete
finite expressions or explicit convergent series.

## Complete physical examples

The source-complete fixtures have a higher-superpotential vertex with one
surviving field, a fully contracted cubic vertex with zero source, and
source-bearing cubic vertices of contracted degree two. Every graph,
source, word, and field count is recorded in
[`ear_certificate.json`](ear_certificate.json).
[FULL_RESULTS.md](FULL_RESULTS.md) prints all seven example kernels and
their complete first-current spinors in the original source basis.

| Loops | Superpotential degree | Sourced ears | Result recorded |
|---:|---:|---:|---|
| 4 | 5 | 1 | Bare kernel and complete first-current spinor |
| 5 | 6 | 2 | Bare kernel and complete first-current spinor |
| 6 | 7 | 3 | Bare kernel and complete first-current spinor |
| 8 | 9 | 5 | Bare kernel and complete first-current spinor |
| 11 | 12 | 8 | Bare kernel and complete first-current spinor |

At four and five loops, the physical bare values are

\[
-\frac{966893}{8400}+\frac{5148}{175}\zeta(3),\qquad
-\frac{833635741}{2016000}+\frac{262369}{1750}\zeta(3).
\]

The eleven-loop physical bare value is

\[
-\frac{478851806802137300846107313237}{48331644641280000}
+\frac{2422391296201519477646891}{332972640000}\zeta(3).
\]

Its full first-current spinor uses 91 fixed-core source evaluations and
at most 163,571 terms in a residual numerator. The certificate stores both
spinor components and their equivalent coefficients in the original source
basis. No coordinate is omitted from the reported tensor.

These are individual geometric kernels. Their values do not establish
irrationality of a fully assembled quantum product.

## Sourced higher-loop Ward cancellation

[`sourced_ward.py`](sourced_ward.py) evaluates all three actual four/five-loop
quotients obtained by attaching nonzero-source ears to the earlier sliding
Ward graph. An independent tensor interval integration gives their rational
spider boundary. All tree and proper faces are enumerated.

For the physical base sources, the complete signed sums are

\[
-\frac{398239}{2160},\qquad -\frac{70717559}{259200},
\]

although all three individual quotient kernels in each relation contain
nonzero `zeta(3)` coefficients. Their complete expressions and source/edge
maps are in [`sourced_ward_certificate.json`](sourced_ward_certificate.json).
This is a graph-level identity; Wick/color/operator assembly is not asserted
by this script.

## Reproduce

From the repository root:

```bash
python3 -m pip install -r explore/circuit-general/forest-general-reduction/requirements.txt
python3 explore/circuit-general/forest-general-reduction/run_checks.py --examples
```

The complete replay log is [`CHECKS.log`](CHECKS.log). To run the stages
separately:

```bash
python3 code/tests/test_forest_ears.py
python3 explore/circuit-general/forest-general-reduction/certify.py
python3 explore/circuit-general/forest-general-reduction/sourced_ward.py
```

The first command covers independent current/density checks, reversed and
permuted graphs, exact source moments, complete original-graph period
controls, multistage source transfer, explicit budgets, and exact-series
fallback. The certificate replay also independently integrates the complete
original four/five-loop control graphs. Its eleven-loop first-current
calculations are the longest stage. `--write` regenerates either certificate.
