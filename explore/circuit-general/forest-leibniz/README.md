# Forest periods and the quantum Leibniz rule

Rational supercharge kernels do not require every product-forest kernel
to be rational. The Leibniz identity constrains the variation of the
complete quantum product and couples all lower loop orders. The
manuscript already proves that identity. This continuation turns it into
an exact current-word reduction with a separately integrated rational
boundary, and explains precisely what the period cancellations imply.

## The equation being solved

Use the manuscript's full differential and marked product. For
homogeneous inputs independent of the loop parameter, the coefficient
of order `n` obeys

\[
\begin{aligned}
0=\sum_{j=0}^{n}\Bigl(
&Q_j\bigl([\hbar^{n-j}](A\star B)\bigr)
-[\hbar^{n-j}]\bigl((Q_jA)\star B\bigr)\\
&-(-1)^{|A|}[\hbar^{n-j}]\bigl(A\star(Q_jB)\bigr)
\Bigr).
\end{aligned}
\]

Here `[hbar^n]` means the indicated coefficient and `|A|` is the fermion
parity already used in the manuscript. Thus the new product coefficient
is acted on by `Q0`, while every higher `Qj` acts on an earlier product
coefficient. It is not an equality between a single forest integral and
a single spider integral.

At a first order with rational forcing, the equation is a rational
linear system. Rational row reduction gives a rational particular
solution and rational null vectors. Period values can multiply the
null vectors. The resulting bilinear correction is closed for the
classical differential; exactness is a separate condition. At subsequent
orders, periods in earlier product coefficients also enter the forcing.
No arithmetic independence of odd zeta values is assumed.

Moreover, **physical cubic product forests already have rational
coefficients at every finite current order**: their only interaction
sources are at their degree-two roots. The period-valued sourced cores
are needed when additional fields survive, as in higher-degree
superpotentials. A scalar graph value still requires Wick contractions,
couplings, statistics, and the complete labelled graph sum before it
becomes an operator coefficient.

See [DERIVATION.md](DERIVATION.md) for the exact coefficient equation,
source counting, homogeneous freedom, and the distinction between
closed and exact corrections. The original marked-product theorem and
its geometric higher-operation inputs are credited there.

## New complete current-word Ward reduction

The ordered sliding graph has edges

\[
(01,05,13,14,23,24,25,34),
\]

with marked insertions `A=3,B=4`. Its surviving edge-contraction terms
are the three three-loop forests obtained by contracting `01`, `05`,
and `25`, with signs `+,-,+`. All five other contractions have parallel
edges and vanish. There is one nonvanishing proper tight face: the
two-loop diamond on vertices `1,2,3,4`, composed with an outer triangle.
Both pieces have exact rational spider charts.

[ward_three.py](ward_three.py) evaluates **every finite word in all
eight original edge currents**. A word may include the edge being
contracted. Its limiting current is reconstructed from the complete
Kirchhoff equation, including the source at its tail. Removing that
affine term would give an incorrect Ward identity.

The proof works for every such word: multiplying the top current-volume
form by any current polynomial preserves closedness, because its
differential is a linear combination of the already present chord-current
differentials. Restriction to a boundary therefore substitutes the full
limiting currents. The diamond sources are shifted by the outer currents
before its integration, preserving every original source.

### An explicit period cancellation

Write `lambda_v` for the retained source at vertex `v`. The recorded
source assignment can be stated without expanding the spinors:

\[
\lambda_1=\lambda_2=0,\qquad
\lambda_3=\frac{13}{4}\lambda_0+\frac78\lambda_5,\qquad
\lambda_4=-\frac{17}{4}\lambda_0-\frac{15}{8}\lambda_5,\qquad
\langle\lambda_0,\lambda_5\rangle=-8.
\]

The three quotient source patterns comprise two higher-superpotential
patterns and one physical cubic-root pattern. Keep the free spinor
index of the first current insertion on `13`. The three complete values,
in the contraction order `01,05,25`, are

\[
\begin{aligned}
&\left(\frac{3376337}{120960}-\frac{10285}{336}\zeta(3)\right)
  \lambda_0^{\dot a}
 +\left(\frac{171081}{44800}-\frac{12881}{5600}\zeta(3)\right)
  \lambda_5^{\dot a},\\[2mm]
&-\frac{6859}{6912}\lambda_0^{\dot a}
 +\frac{6859}{7680}\lambda_5^{\dot a},\\[2mm]
&\left(-\frac{7617139}{241920}+\frac{10285}{336}\zeta(3)\right)
  \lambda_0^{\dot a}
 +\left(\frac{5034451}{2419200}+\frac{12881}{5600}\zeta(3)\right)
  \lambda_5^{\dot a}.
\end{aligned}
\]

With signs `+,-,+` their sum is exactly

\[
-\frac{1115}{432}\lambda_0^{\dot a}
+\frac{21629}{4320}\lambda_5^{\dot a}.
\]

The right side is independently obtained from the diamond–triangle
boundary. Both spinor components and all zero terms are recorded in
[ward_three_certificate.json](ward_three_certificate.json).
The implementation also gives a universal first-spoke boundary formula
using spinor brackets and exact simplex moments, with no coordinate
expansion in that formula.

The generic source example gives `-219793/1080` for the same insertion.
Its three individual zeta coefficients cancel by the integer identity

\[
-324999+2(157526)+7(1421)=0.
\]

This is a graph Ward calculation with original currents and sources.
The complete Wick and coupling assembly for a specified interacting
theory is a separate calculation. Neither this identity nor a nonzero
scalar period establishes that a period survives that assembly.

## Higher-loop lifts and reconstruction

[lifted_ward.py](lifted_ward.py) replaces the two spokes at vertex `1`
by `k` nested source-free paired stars. The original full graph then
has `3+2k` loops. Insert the deepest first spoke current and eliminate
only the newly added passive stars. In the inherited base edge order,
the exact transferred current polynomial is

\[
\frac{2(-1)^k\,y_{13}^{\dot a}
 \langle y_{13},y_{34}\rangle^k
 \langle y_{14},y_{34}\rangle^k}
 {k!((k+1)!)^2(k+2)!(2k+1)!}.
\]

Here `y` denotes the original edge current already used in the circuit
formalism; the brackets retain its spinor contractions. This polynomial
has current degree `4k+1`. The factorial coefficient is the existing
[paired-star comb identity](../sourced-forests/DERIVATION.md#an-all-loop-higher-valence-comb),
proved there by exact interval and simplex integration at every step.
This continuation applies that transfer to all three Ward terms and
their independently integrated boundary. The implementation compares the complete
factorial polynomial with the actual passive contractions, including
edge permutation signs and every source row. Exact bracket factors
cancel `2k` Kirchhoff denominators before the three-loop core is
integrated.

The separate rational boundary integrates the lifted diamond by explicit
beta moments, then integrates the outer triangle. It does not call the
core period engine. Its full spinor agrees with the independent cubic
forest message engine through eighteen inner loops.

| Full graph loops | Paired insertions | Transferred current terms | Outer-core denominator power |
|---:|---:|---:|---:|
| 5 | 1 | 4 | 8 |
| 7 | 2 | 9 | 10 |
| 9 | 3 | 16 | 12 |
| 13 | 5 | 36 | 16 |
| 19 | 8 | 81 | 22 |

[lifted_ward_certificate.json](lifted_ward_certificate.json) records
every rational and period coefficient for both spinor components:
generic sources at five and seven loops, and the higher-superpotential
source pattern at five, seven, nine, thirteen and nineteen loops.
These are fourteen exact scalar checks of seven complete spinor
identities. Five- and seven-loop graphs have exhaustive face censuses,
including dense-subgraph witnesses for every vanishing quotient.
The larger examples are proved as exact inductive lifts of the base
all-current-word identity; no exhaustive larger face census is asserted.

At nineteen loops, the full signed boundary at the displayed source
assignment is

\[
\begin{aligned}
&\frac{17819061863763165565499114948621}
 {107165802777757324493110193307833004167402618880000000000}
 \lambda_0^{\dot a}\\
&\quad+
\frac{18863367752895416954141643575833}
 {82435232905967172687007841006025387821078937600000000000}
 \lambda_5^{\dot a}.
\end{aligned}
\]

The middle quotient is rational. The two outer quotients have opposite
full `zeta(3)` spinor coefficients, so the relation reconstructs either
outer forest from the other, the rational middle forest and the
rational boundary. `reconstruct_last` performs this exact coefficient
reconstruction and preserves every supplied period label. The
certificates evaluate all three independently to verify it.

## Complete verification through quadratic current order

[certify_words.py](certify_words.py) computes all **153** monomials
through degree two in the sixteen original current coordinates, retaining
every edge and every source. Of these, **152** have nonzero period
coefficients in individual quotient terms; every signed sum agrees
exactly with its independent rational boundary.

It also verifies four mixed cubic words, all bare and first-current
words at the higher-superpotential source pattern, and four nonintegral
rational source controls. In total, **178 complete Ward identities**
are stored in
[current_words_certificate.json](current_words_certificate.json).
Every quotient coefficient and every exact rational boundary is
recomputed during replay. The period cancellations are literal rational
coefficient equalities, not inferences from conjectural transcendence.

[test_ward_local.py](test_ward_local.py) independently checks the
diamond chart against the original generic series-parallel integrator,
all current limits at exactly zero resistance on each contraction
face, and spinor covariance. No numerical limiting extrapolation is used.

The earlier four-loop Ward computation is also replayed explicitly in
[replay_existing_ward.py](replay_existing_ward.py): three wheel kernels
with both `zeta(3)` and `zeta(5)` combine to `-66661/2304`. These older
coefficients and their independent lower-face data are attributed to
the existing [full_ward](../../nonSP/verify/full_ward/README.md)
calculation; they are not presented as newly evaluated periods.

## Run

```python
import sys
sys.path.insert(0, "explore/circuit-general/forest-leibniz")

from ward_three import PHYSICAL_SOURCES, WardThree

with WardThree(PHYSICAL_SOURCES) as calculation:
    result = calculation.verify({(2, 0): 1})
    assert result["rational_boundary"].numerator == 8027
    assert result["rational_boundary"].denominator == 1080
```

Run the scripts from the repository root:

```bash
python3 -m pip install -r explore/circuit-general/forest-leibniz/requirements.txt
python3 explore/circuit-general/forest-leibniz/run_checks.py --examples
```

The example option replays all 178 stored current-word identities,
all fourteen higher-loop spinor-coordinate checks, and the complete
weighted-period combination below.
[CHECKS.log](CHECKS.log) records the verified commands and outcomes.

## Reducing complete period combinations

The new [weighted combination reducer](../symmetry/combination.py)
adds coefficients of equivalent completed periods before requesting
their values. It uses the proved uniqueness, inverse uniqueness,
Fourier and twist transformations with their exact rational factors.
Thus some unknown periods cancel without ever being evaluated.
The terminal engine also reads the existing `K5` and `K5-minus-edge`
certificates.

The [production certificate](../symmetry/forest_combination_certificate.json)
starts with every one of the 7,368 scalar terms of the previously
unresolved five-loop sourced forest. At the declared search bounds it
finitely resolves 1,007 input completed orbits and cancels five unresolved
components. The returned expression retains all **5,694 remaining
weighted completed periods**, together with its finite contribution.
Their complete weighted-graph and coefficient digest is replayed from
the original graph and sources. The finite contribution is not exposed
as a complete value. See the [exact counts and coefficients](../symmetry/REPORT.md#31b-reproducible-full-sourced-forest-control).

## Remaining work

Ward relations need not determine every period. The existing four-loop
source-polynomial system has fifteen homogeneous directions, so
independent period data remain necessary. The general forest evaluator
continues to return its complete convergent series when finite reduction
is unavailable. Universal finite zeta closure and the complete
theory-specific higher-superpotential product are not claimed here.
