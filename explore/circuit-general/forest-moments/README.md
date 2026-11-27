# Exact cubic forests at arbitrary loop order

**27 November 2026 — computational contribution, H=2, T=0.**

The new evaluator computes every finite original-current word on every
cubic marked product forest at supplied exact rational sources. It generates
subtree moments recursively, without precomputed tree charts or a polynomial
in all graph ratios. Sources at every tree root and both insertion vertices
are retained, with their complete conservation law. Results use the
[htk kernel convention](../../../code/htk/CONVENTIONS.md).

The all-loop forest reduction theorem is already in
[Li–Liu–Chang, Theorems b-00249 and b-00405](../../../sources/LLC_main.tex).
The contribution here is its computational realization through local current
messages, exact degree grading, and source translation. The
[derivation](DERIVATION.md) states the local identities, proves the recurrence,
and separates the broader sourced-star problem.

## What changes computationally

| Method | Representation | Limitation addressed here |
|---|---|---|
| Existing [forest ratio evaluator](../../nonSP/verify/full_ward/forest_jets.py) | Original currents in all `2*n` ratios, with factor contraction | Ratio dependence grows with the graph |
| Companion [product-kernel evaluator](https://github.com/lbotao399-stack/All-loop-holomorphic-kernels/blob/main/product_kernels/code/product_kernels.py) | Cached boundary moments and forest assembly | Reads tree charts tabulated through seven interactions |
| New [local forest moments](../../../code/htk/forest_moments.py) | Three boundary currents; only two or four temporary local ratios | Generates arbitrary tree shapes on demand; no chart table |
| New [batched fan Taylor transfer](../../../code/htk/fan_moments.py) | At most `binom(order+3,3)` exact coefficient states | Computes a whole directional Taylor polynomial without enumerating its current words |

Each cubic tree is rooted at its unique degree-two interaction. A leaf has
two boundary spokes, a single-child vertex uses one star reduction, and a
two-child vertex uses paired stars. Local ratios are integrated immediately.
The retained message is a polynomial in three boundary spinors, represented
exactly by six scalar coordinates.

Two details are essential. The marked-current degree of the first child
shifts the second child's degree argument, preserving their correlations.
When a component root is removed, its source translates the insertion source
inside the remaining forest polynomial. This translation is evaluated
completely. Repeated subtrees with identical inserted words share messages.

For one current word of degree `k` on an `n`-loop graph, a conservative bound
on stored boundary coefficients is
`n*(2*n+k+1)*binom(2*n+k+6,6)`.
This is a bound on rational coefficient counts. Integer bit lengths, source
polynomial output, and the number of distinct requested words have separate
costs. There is no table-imposed loop cutoff. The implementation has explicit
term and state budgets and raises an error if they are exceeded.

## Public use

Add the repository's `code` directory to `PYTHONPATH`. Cubic forest messages
require `python-flint`; the fan transfer uses Python's exact `Fraction` type.

```python
from fractions import Fraction
from htk import ForestMoments, forest_graph, fan_taylor

# Three-loop fork, canonical original graph and orientation.
graph, roots = forest_graph(((("A", "B"), ("A", "B")),))
moments = ForestMoments(graph, {roots[0]: (0, 1), 3: (1, 0)}, marked=(3, 4))
assert moments.moment() == 0
assert moments.moment({(2, 0): 1}) == Fraction(1, 8640)
assert moments.moment({(2, 1): 1}) == -Fraction(1, 4320)

# Complete finite displacement jet on all 257 edges of a 128-loop fan.
displacements = [(Fraction(e % 5 + 1, 7), Fraction(e % 7 + 1, 11))
                 for e in range(257)]
jet = fan_taylor(128, (2, -1), (3, 4), displacements, order=8)
assert len(jet.coefficients) == 9
assert jet.peak_states == jet.state_bound == 165
```

`ForestMoments.moment(powers)` accepts a dictionary from an original edge
index and current coordinate to its nonnegative exponent. It returns a raw
current moment, without Taylor signs or factorials. `moments(words)` shares
subtree work across a finite list. `source_polynomial(powers)` gives the full
polynomial in the insertion source with all root sources fixed; the other
insertion source is determined by conservation. Linear combinations of
these values evaluate any finite current polynomial.

`fan_taylor` returns coefficients of the common displacement parameter,
including every sign and factorial from the exponential. It supports
displacements on the original path edges and spokes, including their exact
constant source phase. `fan_current_moment` provides arbitrary original-edge
words, and `fan_endpoint_moment` exposes the existing ordered-simplex moment
identity with unbounded integer exponents. The batched Taylor algorithm,
rather than that already known identity, supplies the new state bound.

`CubicForest` and `recognize_product_forest` accept arbitrary vertex labels,
original edge ordering, and edge orientations. They retain the corresponding
kernel and insertion signs. The existing topology counts and automorphism
factors are exposed through `htk.forest_model`; enumeration is inherently
larger than evaluating one chosen graph.

## Exact results and replay

The [forest certificate](forest_certificate.json) records four complete
integrals and their full insertion-source polynomials:

| Fixture | Loops | Current degree | Largest local polynomial | Source-polynomial coefficients |
|---|---:|---:|---:|---:|
| Alternating single tree | 20 | 3 | 3,509 terms | 86 |
| Asymmetric branched tree | 9 | 4 | 540 terms | 54 |
| Three nonparallel root sources | 9 | 4 | 538 terms | 99 |
| Repeated fork, known-formula regression | 36 | 12 | 4 terms | 1 |

The twenty-loop integral is exactly

\[
\frac{51063405172108037818499}
{12808478037547423258484132867205600038721945600000000000}.
\]

The two nine-loop values are respectively

\[
\frac{330274357}{292188041773056000},\qquad
-\frac{228258973157261949061}{2192323400928460800000}.
\]

Both nine-loop values independently agree with the old global chart.
The thirty-six-loop control reproduces the existing value

\[
\frac1{72^{12}49!}
=\frac1{11805783789483870322440822060126121595397476067369105346800411100447600803840000000000}.
\]

Each certificate specifies the complete original graph, balanced source
assignment, inserted word, all polynomial coefficients, convention, and
implementation hashes. [certify_forests.py](certify_forests.py) regenerates
these fixtures from their tree shapes and checks the independent controls.

The complete [128-loop certificate](fan_certificate.json) contains all nine
coefficients through order eight, all sources, and a deterministic assignment
of nonzero displacements to every original edge. Its recurrence uses 165
states and 384,549 exact transfer additions. The normalized first coefficient
is `19896208/215215`; all higher coefficients are stored as complete rational
numbers in the certificate.

The independent [forest checks](tests_exact.py) compare with the pre-existing
global ratio charts on **all 100 cubic forest topologies through five loops**:
**2,292 exact scalar comparisons**, covering the bare kernel, every original
first current moment, and mixed words of degrees three and five. They also
check every original vertex conservation equation, nonparallel root sources,
source-polynomial translations, orientation changes, asymmetric branched
trees, and known all-order fan and fork controls.

The [fan checks](../../../code/tests/test_fan_moments.py) independently compare
all-edge Taylor coefficients through order four with the existing
series-parallel integrator through four loops. They check mixed original
path/spoke words, the exact source-phase identity, and exponents 255, 256,
and 257. The [graph-model checks](../../../code/tests/test_forest_model.py)
compare signs and currents against the original circuit density, validate
sources, and compare product-forest sparsity recognition with an exhaustive
small-graph test.

```bash
python3 -m pip install -r explore/circuit-general/forest-moments/requirements.txt
python3 explore/circuit-general/forest-moments/run_checks.py --examples
```

Omit `--examples` to run the independent controls alone. Certificates store
every returned rational number; their replay recomputes the full answer and
checks exact equality. The low-loop oracle comparisons and analytic controls
provide independent validation beyond certificate replay.

The recorded [complete check run](CHECKS.log) passes 31 validation groups
and all five higher-loop examples. An additional independent read-through
checked cache keys, degree shifts, source handling, orientation factors,
and exact source-preserving dispatch.

## Scope of “all product forests”

The local polynomial engine covers **all cubic product-forest shapes** with
the source support required by the cubic theorem, at every finite current
order. A zero bare kernel does not justify discarding the tree: the fork
above already has nonzero first moments.

Higher-degree superpotential vertices can leave multiple sources in a
component. The sourced-star transfer then depends rationally on the remaining
resistances. The cubic three-current polynomial recursion cannot be applied
to that transfer. The structural recognizer distinguishes these graphs and
the source validator identifies additional surviving sources.

`ProductForestMoments`, also exported from `htk`, recognizes the broad family
and selects the appropriate implemented method. Cubic root-source forests
use local messages. Higher-valence forests and cubic topologies with extra
surviving sources use the existing series-parallel integrator when its
generic reduction succeeds. This route computes the complete generic
current-Jacobian density and checks it against the original circuit at two
exact rational points before integration. It retains every source and
reuses the density for subsequent current words. Its `method` attribute
records which route was chosen. State and term budgets apply to the cubic
method; the existing series-parallel integrator has its own computational
costs.

The higher-valence three-loop control in
[test_product_forests.py](../../../code/tests/test_product_forests.py) has
complete bare value `55625/432` and a mixed degree-three moment
`4369801/8400`, with nonzero sources at the higher-valence vertex and both
degree-two vertices. Both are exact generic-source integrals; this routing
is existing-method coverage, not a new sourced-star formula.

The [source-complete continuation](../sourced-forests/README.md) now handles
these remaining sourced non-series-parallel forests. It contracts only
permitted passive stars, retains the full residual source polynomial,
cancels exact Kirchhoff factors in current brackets, and integrates the
remaining sourced core. Available periods give finite rational/zeta
expressions; otherwise the adapter returns an explicit convergent rational
series with rigorous remainder bounds. `UnsupportedProductForestReduction`
now occurs only when a finite-only result was requested and its period is
unavailable. The new certificates include physical five/six-loop quartic
forests and higher-valence moments through nineteen loops. Finite zeta
closure for every forest and efficient evaluation at arbitrary size remain
open.
