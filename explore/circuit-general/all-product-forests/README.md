# Solve product forests with the existing exact sourced rule

The manuscript already contains the exact source-transfer and current-area
rules. This contribution makes them executable on **every classified
product-forest topology**, including higher superpotentials and arbitrary
conserved rational sources. It also joins the previously separate integration
methods and adds a direct exact series which avoids the expanded global
current numerator.

The [complete derivation](DERIVATION.md) proves termination of the network
reduction, retains the differentials of all transferred sources, constructs
the compact cube cover, and proves convergence of its rational series.
Every requested finite current polynomial and every free spinor index is
retained. Explicit budgets bound computation and storage.

**Further exact integration:** the [geometric cube expansion](../forest-cube-expansion/README.md)
integrates every series coefficient by a finite polynomial-moment formula.
Its centered binomial expansion has a closed rational remainder bound and
an explicit sufficient order. The [direct JSON entry point](../../../code/solve_product_forest_cubes.py)
retains complete current words, mixed-degree inputs, and free-index tensors.

**Complete finite five-loop continuation:** [All 200 additional topologies](../five-loop-forest-census/README.md)
have finite bare kernels and all 4,400 original first-current entries at the
recorded admissible full source assignments. Together with this directory's
61 cases, the finite census contains 261 topologies and 5,418 current entries
through five loops. [Every five-loop spinor](../five-loop-forest-census/FULL_RESULTS.md)
is written in original-source directions with its free lower index retained.

## What is now executable

| Task | Exact implementation | Output |
|---|---|---|
| Any forest topology, arbitrary sources | [ForestKronPlan](../../../code/htk/forest_kron.py) | Every source-linear current, source transfer, current differential, and full area density |
| Any finite original-current polynomial | [ForestReduction](../../../code/htk/forest_reduction.py) | Complete rational value, finite period expression, or complete convergent rational series |
| Full finite spinor tensor | `ForestReduction.tensor` | Every entry and every named free lower index |
| General direct integration without the global numerator | [ForestCubes](../../../code/htk/forest_cubes.py) | Complete rational telescoping series with rigorous rational intervals |
| Exact JSON input and result | [solve_product_forest.py](../../../code/solve_product_forest.py) | Complete machine-readable values or series specifications |

The cubic theorem still gives finite Dirichlet-factorial answers at every
loop and finite current order. The general sourced algorithm supplies an
exact answer through finite periods or a specified convergent series.
A universal finite Riemann-zeta expression, or a theorem that every complete
nonrenormalizable Wick sum is rational, is not established by this work.

## Immediate calculations

From the repository root:

```bash
python3 -m pip install -r explore/circuit-general/all-product-forests/requirements.txt
python3 code/solve_product_forest.py explore/circuit-general/all-product-forests/five_loop.json --finite-only
python3 code/solve_product_forest.py explore/circuit-general/all-product-forests/sourced_fork.json --finite-only
```

The first input evaluates the complete five-loop kernel

\[
\frac{1187-816\zeta(3)}{435456000}.
\]

Its second request inserts the complete conserved current at vertex 1.
Kirchhoff simplification gives exactly the same value. The old public
dispatcher could return only the general series for this kernel; the new
solver consumes the existing all-source coefficient certificate, verifies
its dependencies, and transports it through marked graph isomorphisms with
the original edge-order sign.

The second input computes the complete three-loop first-current spinor.
Both free-index entries are returned, with their rational and zeta(3)
coefficients. The results of both commands are archived as
[five_loop_solution.json](five_loop_solution.json) and
[sourced_fork_solution.json](sourced_fork_solution.json).

JSON source values and coefficients must be exact integers or fraction
strings. Each request supplies exactly one of `powers`, `polynomial`, or
`tensor_edges`. A current-power row specifies the original edge, a stored
spinor entry, and its nonnegative exponent. `tensor_edges` requests the
entire free-index tensor. An unavailable finite period is returned as the
whole series specification, including every numerator coefficient.
`--finite-only` requires a finite answer. `--series-order` additionally
computes certified finite series terms; those terms do not replace the
full series as the result.

## Python interface

With the repository's `code` directory on `PYTHONPATH`:

```python
from htk import Graph
from htk.forest_kron import ForestKronPlan
from htk.forest_reduction import ForestReduction

graph = Graph(5, [(0,1), (0,2), (1,3), (1,4), (2,3), (2,4), (3,4)])
sources = [(1,2), (2,-3), (0,0), (9,6), (-12,-5)]

plan = ForestKronPlan(graph, marked=(3,4))
response = plan.response([1,2,3,4,5,6,7])
all_original_currents = response.currents(sources)
all_original_current_derivatives = response.current_jacobian(sources)
full_current_area_density = response.density(sources)

with ForestReduction(graph, sources, marked=(3,4)) as forest:
    bare = forest.moment()
    spinor = forest.tensor((0,), indices=("+dot_a",))
    full_tensor = forest.tensor((0,4), indices=("+dot_a", "+dot_b"))
    # A separate direct-current route; no global current numerator is built.
    complete_series = forest.cube_series()
```

`integrate(polynomial)` accepts a sparse polynomial with one exponent for
each original edge and each stored spinor entry. Terms of different degrees
are kept in the same request. Exact source-conservation identities are
simplified before passive stars, sourced ears, or scalar periods are used.
The result remains in the original source and orientation conventions.

The direct cube interface is also available independently:

```python
from htk.forest_cubes import ForestCubes

atlas = ForestCubes(graph, sources, marked=(3,4))
series = atlas.series()
first_complete_level = next(series.terms(0))
assert first_complete_level.lower <= first_complete_level.upper
```

`series.terms(k)` includes **every** chart and cell through level k.
`series.enclose(width)` refines until the complete rational interval has at
most the requested positive rational width, subject to explicit evaluation
budgets. This method gives rigorous enclosures; generic initial bounds can
be broad because the current determinant has cancellations. The cubic and
finite-period methods remain the efficient routes where they apply.

## Every topology through six loops

The census keeps A,B ordered and identifies only interaction-vertex
isomorphisms. It enumerates every unrooted interaction tree, every admissible
spoke assignment, and every unordered multiset of components.

| Loops | All forests | Cubic forests | Higher-valence forests |
|---:|---:|---:|---:|
| 1 | 1 | 1 | 0 |
| 2 | 3 | 3 | 0 |
| 3 | 11 | 8 | 3 |
| 4 | 46 | 23 | 23 |
| 5 | 200 | 65 | 135 |
| 6 | 874 | 194 | 680 |

All 1,135 topologies pass the universal sourceful reduction. At two exact
positive resistance assignments per graph, all **190,326** source-response
coefficients agree with the original grounded-Laplacian calculation.
The full current-area density also agrees with an independent original
cycle-form determinant at two gauge choices. The complete plans include
**1,721** intermediate stages outside the Laman edge count.

The [census certificate](kron_certificate.json) retains every graph, original
source assignment, resistance assignment, density, and full-response digest.
Independent exhaustive enumeration of all simple original graphs through
four loops gives the same topology counts.

## Complete integrated census through four loops

At the explicitly recorded admissible full source assignment for each graph,
all **61 bare kernels** and **1,018 first-current entries** have been evaluated
twice: by direct integration of the complete original graph and by the new
compositional solver. All **692** integrated source-conservation identities
hold exactly. The census includes every original edge and both entries of
every current spinor. No zero bare kernel is used to infer vanishing of
higher current moments.

The [complete formulas](FULL_RESULTS.md) retain the free lower spinor index
and express each first-current spinor in two original-source directions.
The [integral certificate](integral_certificate.json) contains every exact
coefficient. These are fixed-source geometric kernels; arbitrary-source
interpolation and a full Wick/color assembly are separate calculations.

The recorded complete runs took 28.975 seconds for the original-graph
reference and 7.613 seconds for the compositional solver. These timings
describe the supplied run, with three workers; exact equality of all
coefficients is the validation criterion.

## A complete physical five-loop composition

The [five-loop composition certificate](composed_five_loop_certificate.json)
has a degree-five superpotential at vertex 1, source-bearing degree-two
vertices 0,5, and source-free cubic vertices 2,6. A passive star followed by
a sourced ear gives a complete three-loop integral with the source transfer
integrated exactly. The original ordered edges and every source are recorded
in the certificate and replayed by [composed_five_loop.py](composed_five_loop.py).

The complete bare value is

\[
-\frac{877164533}{3456000}+\frac{1409387}{7200}\zeta(3).
\]

The full original current spinor on edge \((6,1)\) is

\[
\begin{aligned}
\int y_{(6,1),+\dot a}\,\Omega_\Gamma
={}&\left(\frac{4175805823}{12441600}
-\frac{6328007}{24000}\zeta(3)\right)\lambda_{0,+\dot a}\\
&+\left(-\frac{3104462473}{20736000}
+\frac{8478617}{72000}\zeta(3)\right)\lambda_{3,+\dot a}.
\end{aligned}
\]

The complete mixed-degree spinor insertion also closes:

\[
\begin{aligned}
\int\left(\lambda_{0,+\dot a}-\frac27y_{(6,1),+\dot a}\right)\Omega_\Gamma
={}&\left(-\frac{5438599621}{15552000}
+\frac{4879469}{18000}\zeta(3)\right)\lambda_{0,+\dot a}\\
&+\left(\frac{443494639}{10368000}
-\frac{1211231}{36000}\zeta(3)\right)\lambda_{3,+\dot a}.
\end{aligned}
\]

All five scalar storage entries of these three complete objects agree with
an independent composition of the existing passive-star and sourced-ear
integrators, with no residual period. The two original-source directions
express each entire spinor; no free index is contracted away.

## Compact direct-current integration

The new cubes are indexed by a spanning tree and an ordering of its
complementary edges. The normalized tree denominator has constant term 1
and nonnegative coefficients on every closed cube.

| Complete graph | Spanning trees | Compact cubes | Fully ordered Hepp sectors |
|---|---:|---:|---:|
| Three-loop sourced K4 subdivision | 24 | 144 | 5,040 |
| Five-loop irreducible sourced forest | 195 | 23,400 | 39,916,800 |

These are exact chart counts, not runtime estimates. The new series works
directly with the complete source-dependent currents and their derivatives.
It avoids the global current numerator and period-table construction.

## Reproduce

```bash
python3 explore/circuit-general/all-product-forests/run_checks.py
```

This replays the focused module checks, the independent graph census, the
1,135-topology source-response/density certificate, the complete integrated
four-loop census, the direct-current cube certificate, and the command-line
examples. The original-graph integral control is independently replayed with:

```bash
python3 explore/circuit-general/all-product-forests/certify_integrals.py --backend original
```

The committed logs and certificates record the actual completed checks.
