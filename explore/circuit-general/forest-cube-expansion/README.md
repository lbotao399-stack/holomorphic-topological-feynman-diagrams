# Every product forest as an explicit geometric rational series

The existing sourced reduction already reaches every product-forest topology.
This addition supplies a complete algebraic answer at every finite original
current order: a finitely stored convergent rational series, explicit exact
coefficients, and a closed rational bound for all remaining coefficients.
It applies to arbitrary balanced rational source assignments, including
higher-superpotential sources. The proof also applies coefficientwise to
the complete finite polynomial in arbitrary source spinors.

The new result uses the existing compact spanning-tree cubes. It does not
require a finite scalar-period table or the expanded original Schwinger
numerator. Each cube's full current determinant is evaluated polynomially,
with every original source and every derivative of every source transfer
retained. [Full proof](DERIVATION.md).

## Exact reduction and termination

On each closed cube the normalized Kirchhoff polynomial lies between one
and the number of spanning trees. Expanding about the exact midpoint of
those endpoints gives absolute geometric ratio

\[
\frac{\text{number of spanning trees}-1}
{\text{number of spanning trees}+1}<1.
\]

Every series coefficient is a finite polynomial cube integral. Each monomial
is integrated by the exact product of reciprocal positive integers.
The augmented full-current determinant gives initial denominator power
`2n+1`, and exact polynomial divisions reduce it further whenever possible.

The full tail is a geometric power times an explicitly given finite
polynomial in the truncation order. Thus every positive rational requested
interval width has a finite sufficient order. `required_order` computes
that order before expanding the corresponding coefficients. Explicit
storage and order budgets can stop a request; a partial answer is never
reported as the exact integral.

This is a universal exact series evaluation, with a finite representation
of the entire answer. A universal finite expression in Riemann zeta values
and universal period cancellation in a complete Wick/color sum are separate
questions. Polynomial support and atlas size can still grow rapidly.

## Executable interface

```python
from fractions import Fraction
from htk.graph import Graph
from htk.forest_cube_expansion import ForestCubeExpansion

graph = Graph(3, ((0, 1), (0, 2), (1, 2)))
engine = ForestCubeExpansion(
    graph, ((1, 0), (0, 1), (-1, -1)), marked=(1, 2))

bare = engine.series()
complete_formula = bare.specification()
sufficient_order = bare.required_order(Fraction(1, 1000))
certified_interval = bare.enclose(Fraction(1, 1000))
complete_original_current = engine.spinor(0)
```

`specification()` stores every complete numerator, denominator, power,
coefficient rule, and exact tail rule, together with the original graph,
balanced sources and requested current polynomial. `series(polynomial=...)`
accepts mixed degrees with one exponent for every original edge/current
entry. `spinor(edge)` returns both entries of the complete original current.

The [JSON command](../../../code/solve_product_forest_cubes.py) supports
bare kernels, complete free-spinor tensors, repeated-edge tensors, and
mixed current polynomials. The default backend is this geometric expansion;
`--backend midpoint` retains the earlier independent current-cube route.
Every requested tensor entry must succeed before an output file is replaced.

```bash
python3 code/solve_product_forest_cubes.py \
  explore/circuit-general/forest-cube-expansion/triangle_input.json \
  --series-order 6 --width 1/1000 --output /tmp/triangle_result.json
```

The complete [input](triangle_input.json) and [result](triangle_result.json)
contain the bare kernel, all three full first-current spinors, and one
mixed-degree polynomial. The three original source rows are retained.

## Fully summed independent control

For the triangle and sources in the input, every odd coefficient is zero
and every even coefficient is `-3/2^(2j+3)`. Therefore the complete series
has exact sum `-1/2`. At requested width `1/1000`, the sufficient order
is 18, the exact finite sum is `-1048575/2097152`, and the certified
absolute tail is `87/262144`. The separately derived actual tail is
`1/2097152`. The finite sum is not substituted for the exact result.

The full first-current spinors and the independent mixed simplex value
`-7/12` are written out in the derivation and the certificate. Every entry
is enclosed at the requested width.

## Complete sourced three-loop control

The ordered graph `(01,03,04,12,14,23,34)`, marked at `3,4`, has 144
compact cubes. At the recorded balanced higher-source assignment,
the new augmented density and every original current agree with direct
original-network calculations in **every cube**.

The certificate stores the complete infinite rational specifications of
the bare kernel and both entries of the first current on edge `03`.
Their first three coefficients are computed exactly. Separate existing
finite-period values are

\[
\frac{2023}{96},\qquad
\frac{5477}{864}-\frac{304}{15}\zeta(3),\qquad
-\frac{8205703}{302400}-\frac{6289}{2100}\zeta(3).
\]

The two latter entries together are the complete free spinor in the
recorded source basis; the source rows are all included in the output.
Rigorous rational zeta bounds verify compatibility with every displayed
finite series enclosure. The output also records a sufficient order for
width `1/1000`; coefficients through that much higher order have **not**
been computed in this control. No precise three-loop enclosure is claimed.

## Verification

```bash
python3 explore/circuit-general/forest-cube-expansion/run_checks.py
```

The [readable certificate summary](geometric_certificate.json) links to
[the complete compressed certificate](geometric_certificate.json.gz),
which includes every three-loop chart control, all complete series
coefficient data, all triangle outputs, and source-file digests. The summary
records the SHA256 of its uncompressed canonical JSON. `certify.py --write`
regenerates both; the default command replays and compares the entire result.

The focused tests verify the finite tail identity by independent finite
subtraction, the all-order triangle coefficient formula, full original
source and density maps, exact edge signs, complete mixed-degree inputs,
zero-loop values, sufficient-order bounds, and atomic recovery after a
budget failure. The separate command tests verify both backends and every
entry of repeated-edge tensors. These are internal computational checks,
not peer review.
