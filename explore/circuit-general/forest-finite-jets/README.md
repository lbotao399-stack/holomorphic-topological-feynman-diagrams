# Complete five-loop kernel with arbitrary supported source spinors

This calculation completes the **entire bare source polynomial** on the ordered graph

```text
(01,12,03,04,25,26,35,36,45,46,56),  marked vertices (5,6).
```

The sources at vertices 1, 3 and 5 are arbitrary spinors, the source at vertex 6 is
`lambda_6 = -lambda_1-lambda_3-lambda_5`, and the sources at vertices 0, 2 and 4
vanish structurally. This is the complete kernel for that source support, with
`H=2`, `T=0` and no current insertion. It does not assert arbitrary sources on the
three structurally empty vertices, nor arbitrary current moments.

## Full covariant answer

The spinor bracket is
\(\langle\lambda_i,\lambda_j\rangle=\epsilon^{\dot a\dot b}\lambda_{i,+\dot a}\lambda_{j,+\dot b}\), with
\(\epsilon^{\dot+\dot-}=1\). The kernel is the sum of the ten rows below: multiply the
monomial by the rational entry, by the second entry times \(\zeta(3)\), and by
the third entry times \(\zeta(5)\), and add every contribution. These are reduced
kernels in the repository's oriented-edge normalization; no Wick or coupling
factor has been attached.

| Full spinor monomial | Rational coefficient | Coefficient of \(\zeta(3)\) | Coefficient of \(\zeta(5)\) |
|---|---:|---:|---:|
| \(\langle\lambda_1,\lambda_3\rangle^4\langle\lambda_1,\lambda_5-\lambda_6\rangle\) | \(-456409/10450944000\) | \(-20101/217728000\) | \(1/6720\) |
| \(\langle\lambda_1,\lambda_3\rangle^4\langle\lambda_3,\lambda_5-\lambda_6\rangle\) | \(-99833/2090188800\) | \(-1549/17418240\) | \(1/6720\) |
| \(\langle\lambda_1,\lambda_3\rangle^2\langle\lambda_1,\lambda_5-\lambda_6\rangle^3\) | \(18001/387072000\) | \(3571/27216000\) | \(-1/5040\) |
| \(\langle\lambda_1,\lambda_3\rangle^2\langle\lambda_3,\lambda_5-\lambda_6\rangle\langle\lambda_1,\lambda_5-\lambda_6\rangle^2\) | \(79091/696729600\) | \(127/362880\) | \(-1/1920\) |
| \(\langle\lambda_1,\lambda_3\rangle^2\langle\lambda_3,\lambda_5-\lambda_6\rangle^2\langle\lambda_1,\lambda_5-\lambda_6\rangle\) | \(18733/139345920\) | \(2887/7257600\) | \(-1/1680\) |
| \(\langle\lambda_1,\lambda_3\rangle^2\langle\lambda_3,\lambda_5-\lambda_6\rangle^3\) | \(1433/19906560\) | \(7613/43545600\) | \(-11/40320\) |
| \(\langle\lambda_3,\lambda_5-\lambda_6\rangle\langle\lambda_1,\lambda_5-\lambda_6\rangle^4\) | \(-29/9331200\) | \(-79/17418240\) | \(1/120960\) |
| \(\langle\lambda_3,\lambda_5-\lambda_6\rangle^2\langle\lambda_1,\lambda_5-\lambda_6\rangle^3\) | \(-18811/2090188800\) | \(-17/806400\) | \(1/30240\) |
| \(\langle\lambda_3,\lambda_5-\lambda_6\rangle^3\langle\lambda_1,\lambda_5-\lambda_6\rangle^2\) | \(-27217/2090188800\) | \(-361/14515200\) | \(1/24192\) |
| \(\langle\lambda_3,\lambda_5-\lambda_6\rangle^4\langle\lambda_1,\lambda_5-\lambda_6\rangle\) | \(-587/149299200\) | \(-479/43545600\) | \(1/60480\) |

The two pure fifth powers have zero coefficient. Every other possible degree-five
bracket monomial vanishes by the symmetries proved below. **No scalar period is
left unresolved in any source coefficient.**
The focused check also proves that each of the ten listed real coefficients
is nonzero: exact sums of the first 64 terms of each zeta series, with
integral-test rational bounds on both tails, give an interval excluding zero
for each coefficient. No independence assumption about zeta values is used.

The coefficient multiplying \(\zeta(5)\) also has the useful exact form

\[
\begin{aligned}
\frac{\langle\lambda_5,\lambda_6\rangle}{60480}\bigl(&
18\langle\lambda_1,\lambda_3\rangle^4
-33\langle\lambda_1,\lambda_3\rangle^2
   \langle\lambda_3,\lambda_5-\lambda_6\rangle^2\\
&-39\langle\lambda_1,\lambda_3\rangle^2
   \langle\lambda_3,\lambda_5-\lambda_6\rangle
   \langle\lambda_1,\lambda_5-\lambda_6\rangle
-24\langle\lambda_1,\lambda_3\rangle^2
   \langle\lambda_1,\lambda_5-\lambda_6\rangle^2\\
&+2\langle\lambda_3,\lambda_5-\lambda_6\rangle^3
   \langle\lambda_1,\lambda_5-\lambda_6\rangle
+3\langle\lambda_3,\lambda_5-\lambda_6\rangle^2
   \langle\lambda_1,\lambda_5-\lambda_6\rangle^2\\
&+\langle\lambda_3,\lambda_5-\lambda_6\rangle
   \langle\lambda_1,\lambda_5-\lambda_6\rangle^3\bigr).
\end{aligned}
\]

Here source conservation gives
\(\langle\lambda_3,\lambda_5-\lambda_6\rangle+
\langle\lambda_1,\lambda_5-\lambda_6\rangle
=2\langle\lambda_5,\lambda_6\rangle\).
The factorization is checked by exact expansion. It explains why the earlier
specialization with zero source at vertex 5 had no \(\zeta(5)\) term. Its complete
value for arbitrary remaining supported sources is
\[
\frac{1187-816\zeta(3)}{435456000}\langle\lambda_1,\lambda_3\rangle^5.
\]
The earlier numerical fixture has \(\langle\lambda_1,\lambda_3\rangle=1\).
At the source sample named `(3,1)` in the certificate the \(\zeta(5)\) coefficient
is exactly \(1/280\).

## Why these coefficients determine the complete source dependence

The bare five-loop numerator has source degree ten. A simultaneous linear
transformation of the two spinor coordinates multiplies the kernel by the fifth
power of its determinant. After conservation there are three unrestricted
spinors: the sources at vertices 1 and 3, and the difference of the marked
sources. Their three brackets have no polynomial relation, so the kernel is a
homogeneous polynomial of bracket degree five. This gives 21 initial monomials.

Exchanging the marked vertices 5 and 6 swaps the edge pairs 25/26, 35/36 and
45/46. The remaining edge 56 reverses orientation, which contributes +1 for
`H=2`; the three edge transpositions contribute -1. The kernel is therefore odd
in the difference of the marked sources. Only total degrees one, three and five
in the two brackets involving that difference survive: twelve monomials.

When the source at vertex 1 vanishes, its two series currents coincide, and the
bare projective determinant has dependent source columns. The kernel is zero;
this excludes the pure fifth power of the bracket of vertex 3 with the marked
source difference. When the source at vertex 3 vanishes, exchanging vertices 3
and 4 preserves all sources and swaps the three edge pairs 03/04, 35/45 and
36/46. Its odd sign excludes the other pure fifth power. Ten monomials remain.

`source_family.py` uses an exact eleven-node interpolation containing one of
these excluded powers as a zero control. The resulting ten nonzero scalar
numerator dictionaries have exact rational rank ten. A fresh source frame,
with a nonunit determinant transformation and rational marked sources, agrees
in **all 14,916 scalar coefficients**. The numerator identities extend from the
open set where the first two source spinors are independent to every source
specialization by polynomiality; degenerate source spinors are included.

## Computational improvement

`../symmetry/tree_minor.py` jointly chooses a tree-supported particular source
flow and the positively oriented fundamental cycle basis of that tree before
forming the existing homogeneous linear minor. Both choices preserve every
source and every original edge convention. Changing only the source-flow tree
is less effective on the measured five-loop example.

The planner ranks linear-matrix support within an explicit finite tree scan,
then tries a bounded number of exact determinants. Each attempt retains the
caller's polynomial term and matrix size limits. The support score is a
heuristic; a low score can still hit the exact arithmetic budget. Only that
specific budget failure triggers a retry. Exhaustion raises explicitly, leaving
the numerator cache empty. No source coefficient is set to zero. Successful
numerators pass independent exact Kirchhoff density checks before installation.

For the source sample `(3,1)`, the observed default numerator time was
**109.748 seconds**. Joint tree/basis selection produced the 37,514-term
numerator in **3.796 seconds**, or **4.732 seconds including comparison of all
12,014 scalar-core coefficients**. These are elapsed times from one controlled
exploratory run, not machine-independent performance guarantees.
`benchmark_trees.py --baseline` reruns both exact computations; without that
flag it checks the optimized output against the full certificate digests.

`../symmetry/batch_witness.py` shares one sparse elimination history across
all source targets, then reverses that history separately for each target to
produce rational weights on the **original complete rows**. Direct summation of
those rows checks every coefficient and every boundary constant. The bounded
exploratory solve used **65,226 scalar variables and 258,089 inherited rows**, and
finished in **43.922 seconds including loading**. The clean certificate records
its own deterministic row and operation counts. Both runs close all eleven
candidate source coefficients, including the structural zero.

## Direct use

With this folder and `code` on the Python import path:

```python
from source_family import five_loop_source_kernel
value = five_loop_source_kernel(
    ((0,0),(1,0),(0,0),(0,1),(0,0),(0,0),(-1,-1))
)
# ExactPeriod({"1": Fraction(1187,435456000),
#              "z3": Fraction(-17,9072000)})
```

The evaluator accepts every exact balanced source specialization on the stated
support, including degenerate spinors. It rejects nonzero sources at 0, 2 or 4.

## Exact replay and scope of verification

```bash
python3 explore/circuit-general/forest-finite-jets/run_checks.py
```

This runs the focused tree/basis, orientation, parallel-edge, exact-zero,
allocation-budget, stopped-target and factorization checks, then:

1. Reconstructs all eleven complete original source numerators and every scalar
   coefficient after exact series-chain integration and graph automorphisms.
2. Checks the full fresh rational source frame.
3. Regenerates the bounded local, dimension-shift and polynomial-annihilator
   identities from the unchanged earlier implementations and existing period
   tables. Every accepted row retains all coordinate-face terms; invalid
   convergence or unavailable-boundary rows are rejected as whole identities.
4. Sums the stored original-row weights and verifies zero residual separately
   for every source coefficient, without rerunning elimination.
5. Checks immutable input/table/implementation hashes and the full previously
   completed source specialization.

The independent polynomial-annihilator and lower-level boundary controls from
`../forest-finite-target/check_annihilators.py` remain applicable: the same
verified generators and row formulas are used unchanged. This contribution
constructs no complete scalar period table and modifies no earlier hashed
engine or certificate.

To regenerate the shared elimination and witnesses explicitly:

```bash
python3 explore/circuit-general/forest-finite-jets/certify_sources.py --write
```

`source_certificate.json` holds source choices, complete numerator hashes,
finite coefficients, construction bounds and provenance. The deterministic
`source_witnesses.json.gz` contains only exact original-row indices and rational
weights; all source numerators and rows are reconstructed from the stated inputs.
It contains no truncated numerator or omitted residual. `CHECKS.log` distinguishes
certificate generation from direct stored-witness replay.
