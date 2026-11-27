# Five-loop bulk3: complete supported-source polynomials and full currents

This calculation evaluates **92 of the 410 bare bulk3(3,5) source
polynomials**, including **89 previously unevaluated non-SP graphs**.
It also evaluates **every original first-current vector on 37 graphs**:
592 complete vectors, with all 1,776 coefficients of the three independent
source vectors retained. Thirty-four of these full current tensors are
non-SP. These are evaluations of the original five-loop graphs in
`fiveloop/enum_d3.json`, beyond the separate product-forest census.

The source support is precise: the four bivalent vertices, in the order
`S0` of each original graph, have arbitrary three-component sources with
zero sum. Sources at all other vertices vanish. These results include
singular source frames. They do not replace the arbitrary-source D=2
calculations, the twistor census, higher current orders, or Wick and color
assembly.

| Original core | Complete bare polynomials | Complete first-current tensors |
|---|---:|---:|
| SP | 3 | 3 |
| K4 | 11 | 11 |
| prism | 32 | 0 |
| K33 | 12 | 12 |
| V8_0, a two-sum of K4 cores | 11 | 11 |
| V8_1 | 4 exact signed zeros | 0 |
| V8_2 | 16 exact signed zeros | 0 |
| V8_4 | 3 exact signed zeros | 0 |
| **Total** | **92** | **37** |

All complete answers are in [`certificate.json`](certificate.json), with
the original ordered oriented edges and source vertices. The individual
[`cases`](cases) contain the complete finite-target digests and evaluation
records. The aggregate explicitly lists the 318 unresolved bare indices
and 373 unresolved first-current indices. The 32 prism current tensors
require denominator level nine; the available prism table stops at eight.
Historical logs record a complete 223,549-orbit prism level-nine build and
77,322 successful dimension checks, but its JSON table is absent from this
checkout and its current Git tree. Restoring that declared input would
remove this particular data limitation.
The remaining irreducible V8 families require additional scalar period
evaluation. A signed zero of a bare kernel never implies a zero current
tensor here.

## Reconstructing the entire source polynomial

The source vertices are `s0,s1,s2,s3`, with
$\lambda_{s3}=-\lambda_{s0}-\lambda_{s1}-\lambda_{s2}$.
For each graph the complete bare result is

$$
K_\Gamma(\lambda;0)=
\det\!\begin{pmatrix}\lambda_{s0}\\\lambda_{s1}\\\lambda_{s2}\end{pmatrix}^{5}
\;\texttt{bare\_coefficient}.
$$

For every original oriented edge $e$, the complete first-current moment,
with its free source index $a$ retained, is

$$
\int D_\Gamma(t;\lambda)y_e^a(t;\lambda)\,\Omega
=\det\!\begin{pmatrix}\lambda_{s0}\\\lambda_{s1}\\\lambda_{s2}\end{pmatrix}^{5}
\sum_{j=0}^{2}\lambda_{sj}^{a}
\;\texttt{all\_original\_first\_current\_coefficients}[e][j].
$$

The stored convention is the positive current moment. With the original
kernel exponential $\exp(-\sum_e y_e\cdot z_e)$, the derivative with
respect to $z_e$ at zero is the negative of this displayed vector.

The coefficient dictionaries use `1`, `z3`, `z5`, and `z3*z3` for
$1,\zeta(3),\zeta(5),\zeta(3)^2$. For example, in the canonical source frame,

$$
\begin{aligned}
K_{3}&=\frac{311}{246903552000}-\frac{\zeta(3)}{952560000},\\
K_{14}&=\frac{47501}{9258883200000}-\frac{61\zeta(3)}{14288400000},\\
K_{46}&=-\frac{61437767}{18517766400000}
-\frac{189281\zeta(3)}{21432600000}
+\frac{4\zeta(5)}{297675},\\
K_{62}&=\frac{2804909}{1866590853120}
-\frac{9019\zeta(3)}{3600676800}
+\frac{29\zeta(3)^2}{27783000}.
\end{aligned}
$$

These displayed numbers multiply the full determinant power above; they
are not isolated source samples.

**Proof of source completeness.** Start with the three unit source
vectors and their negative sum. Any invertible matrix of three source
vectors takes that frame to the requested balanced sources. Each of the
five loop-current rows for a source direction transforms linearly by this
matrix. The complete density therefore acquires its determinant to the
fifth power. Every original current transforms by the same matrix once.
This proves both displayed formulas on the open set of invertible source
frames. The original numerators are polynomials in all source entries, so
the identities extend to determinant-zero frames. No interpolation in the
source parameters or omission of a current direction is involved.

This source argument also holds at every loop number in the connected
cubic bulk3 family, provided every original vertex is bivalent or
trivalent and only the bivalent vertices carry sources. There are
$3n+1$ edges and $2n+2$ vertices, so the number of bivalent vertices is
exactly $3(2n+2)-2(3n+1)=4$. The density transforms by the determinant
to the power $n$, and a full original current transforms once more by
the source matrix. The calculation below evaluates the resulting scalar
and vector coefficients specifically at $n=5$.

## Exact numerator and scalar evaluation

For these graphs $n=5$, $H=3$, and there are sixteen original edges. The
unchanged linear-minor construction in
`explore/nonSP/verify/k4_periods/current_minor.py` produces the entire
homogeneous polynomial

$$P(t)=U(t)^8 D_\Gamma(t),\qquad \deg P=24,$$

as a $24\times24$ homogeneous linear determinant. This replaces the old
whole-graph interpolation requirement of 4,855,374,080 points in the
specific four-source family. Exact alternative source/cycle spanning
trees control intermediate polynomial growth; their complete resulting
polynomials agree. The positively oriented cycle matrix is used for the
density. The unmodified fundamental-cycle matrix is used for affine
reconstruction of original currents.

For a maximal series chain and an original numerator monomial
$\prod_e t_e^{p_e}$, set its total scalar power to
$\alpha_{\rm chain}=\sum_{e\in\rm chain}(p_e+1)$.
The complete positive simplex integration contributes

$$
\prod_{\rm chains}
\frac{\prod_{e\in\rm chain}p_e!}{(\alpha_{\rm chain}-1)!}.
$$

Every surviving scalar term is absolutely convergent: the original
physical density is integrable, and the integer-power monomial
convergence criterion applies to its collected numerator. The same
statement holds with the full bounded original-current insertion. No
divergent term is discarded or assigned zero. The period worker performs
read-only exact evaluation at the required level, using the inherited
wheel formulas, prism/K33 tables, and normalized K4 two-sums. It refuses
unavailable periods and never builds a new table implicitly.

### Integrating the full current before polynomial expansion

Let the original particular tree flow be $y_0$ and use the existing
cycle matrix $C$. The exact current identity is

$$
y=y_0-C^{\mathsf T}
\bigl(C\,\operatorname{diag}(t)C^{\mathsf T}\bigr)^{-1}
C\,\operatorname{diag}(t)y_0.
$$

The matrix between parentheses depends only on the total chain
resistances. Its determinant is $U$, and its full adjugate has degree four.
Thus one can first integrate both $P$ and every $P t_e$ over the chain
fractions, then multiply by the adjugate in total-chain variables.
For the monomial above, insertion of $t_e$ raises the appropriate chain
power by one and multiplies its exact simplex coefficient by

$$\frac{p_e+1}{\alpha_{\rm chain}}.$$

[`chain_currents.py`](chain_currents.py) retains the complete affine
$y_0$ term and all loop/source entries. It verifies all twenty-five
polynomial matrix-adjugate identities. This rearranges finite rational
sums and avoids a multi-million-term intermediate polynomial on the
sixteen original edges. The resulting first-current scalar denominator
is nine, with only exact common-$U$ cancellations allowed.

### Complete original-current relations

Every original vertex equation

$$\sum_{e\text{ leaving }v}\int D_\Gamma y_e^a\Omega
-\sum_{e\text{ entering }v}\int D_\Gamma y_e^a\Omega
=\lambda_v^a K_\Gamma$$

is retained. In addition, a graph automorphism transports each original
edge with its orientation sign and permutes all four source vertices.
The density multiplier is the sign of the full vertex permutation,
times $(-1)$ to the number of reversed original edges, times the fifth
power of the source-permutation sign. The source permutation acts on the
three independent source vectors by its complete three-by-three matrix.
`current_relations.py` records the resulting exact equations on all
forty-eight current coefficients and the bare coefficient. Rational row
reduction selects independent integrals; every original current entry is
then reconstructed and checked against every original equation. A
nonzero bare term is never dropped from these affine equations.

The same orientation law proves the twenty-seven exact bare zeros in
[`symmetry_certificate.json`](symmetry_certificate.json). Twenty-three
lie beyond the available denominator-eight families. Its independent
checker classifies all 410 graphs and checks 54 full original density
transformations and 2,592 original current-entry transformations.

## Independent exact controls

[`independent_controls.py`](independent_controls.py) performs the
following complete checks and compares its result to
[`independent_certificate.json`](independent_certificate.json):

* All three SP graphs are integrated using the original polynomial ratio
  Jacobian and simplex integrals, without scalar period tables: three bare
  values and all 144 original current coefficients.
* Graphs 3 and 46 use different source/cycle trees and reintegrate every
  original edge and all three source directions: two bare values and
  all 96 original current coefficients.
* Before integration, 63 entire scalar dictionaries are compared between
  original-edge multiplication and exact integration over chain
  fractions: **101,517 exact rational coefficients** agree.
* Source matrices with determinants $4,-2,0$ give six complete
  coefficientwise numerator comparisons, including the degenerate
  polynomial identity, and 288 original-current source-transform checks.
* Every stored full current tensor also satisfies the complete source and
  automorphism equations in the quick census replay.

The original-Laplacian point checks establish normalization independently
of the linear minor; they are not used to infer a polynomial identity.
The non-SP period values still depend on the declared inherited exact
period tables and wheel formulas. The complete dictionary controls test
the new reduction before those period lookups.

## Reproduction

Run from the repository root with Python 3, python-flint, SymPy, and
NetworkX installed. No network access is used.

**Quick stored-result replay** checks every complete result, every
reconstruction equation, 136 original-density control points, all declared
dependencies, and the precise unresolved sets. It does not repeat all
scalar integrations:

```bash
python3 explore/circuit-general/five-loop-completion/verify.py
```

**Independent reintegration** recomputes the five complete control
tensors and all coefficient/source checks. The zero proof has its own
independent original-network replay:

```bash
python3 explore/circuit-general/five-loop-completion/independent_controls.py
python3 explore/circuit-general/five-loop-completion/symmetry_controls.py
```

For a full fresh computation, the bounded original route writes complete
targets and explicit failures. Alternative trees and exact current
fraction integration then finish the required graphs. These commands
regenerate the result files; the first command can return a nonzero exit
status for intermediate-budget failures which the following commands
resolve. The final verifier requires every declared answer to be complete.

```bash
python3 explore/circuit-general/five-loop-completion/certify.py --indices 0-68 --currents lower --jobs 2 --seconds 600
python3 explore/circuit-general/five-loop-completion/planned_cases.py --indices 19,25,39 --seconds 180
python3 explore/circuit-general/five-loop-completion/finish_currents.py --indices 47,48,51,53-55 --seconds 600
python3 explore/circuit-general/five-loop-completion/resume_currents.py --indices 56-68 --seconds 900
python3 explore/circuit-general/five-loop-completion/independent_controls.py --write
python3 explore/circuit-general/five-loop-completion/verify.py --write
python3 explore/circuit-general/five-loop-completion/verify.py
```

`resume_currents.py` saves every complete independent scalar target before
continuing, and validates its dependencies on resumption. A timeout leaves
an explicitly partial record; it does not label the whole current tensor
complete. `--write` is required to refresh comparison certificates.

## Remaining five-loop work

The twelve historical D=2 gaps remain:

| Core | Original D=2 indices |
|---|---|
| Z5 | 49, 50, 51, 53, 54, 58 |
| K33+e | 57, 61, 65 |
| V7_0 | 66 |
| V7_2 | 68 |
| V7_3 | 67 |

Their generic source numerators and denominator-seven periods are a
different problem. The historical report describes large numerator
caches for three cases, but those particular caches are absent from this
checkout. The bounded reconstruction experiment did not produce a new
D=2 value. This directory claims no completion of those twelve kernels,
the remaining 318 bare bulk3 polynomials, or the non-SP twistor families.
