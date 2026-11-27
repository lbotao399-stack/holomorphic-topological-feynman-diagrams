# Reduce the whole sourced numerator before solving scalar periods

The complete missing five-loop product-forest fixture reduces from **7,368
scalar terms to 19 retained periods by local identities, then to 17 with
687 exact lower-level constraints**. Every original source and every signed
coefficient is retained. The best finite contribution is accompanied by all
17 periods below; this is an incomplete finite evaluation.
The existing public `ExactForestSeries` remains an exact representation of
the complete original kernel.

## Exact fixture and identity

The original ordered oriented edges, balanced sources, and marked vertices are

```python
edges = ((0,1),(1,2),(0,3),(0,4),(2,5),(2,6),
         (3,5),(3,6),(4,5),(4,6),(5,6))
sources = ((0,0),(1,0),(0,0),(0,1),(0,0),(0,0),(-1,-1))
marked = (5,6)
```

For the scalar core, `I_7(powers)` denotes the existing unnormalized
projective period `integral_Omega product(t_e**(powers_e-1))/U**7`.
Its edge order is

```python
((0,1),(0,2),(0,3),(1,4),(1,5),(2,4),(2,5),(3,4),(3,5),(4,5))
```

There is no `product Gamma(powers_e)` division in this convention. The full
original kernel equals

\[
-\frac{3706491467747}{27935373312000}
-\frac{1906264441}{5388768000}\zeta(3)
+\frac{57527}{106920}\zeta(5)
+\sum_{\text{all rows below}}\text{coefficient}\;I_7(\text{powers}).
\]

| Powers, in the displayed core edge order | Coefficient |
|---|---:|
| `(1,3,4,3,4,4,4,3,4,5)` | `-1918/45` |
| `(1,3,4,3,4,4,5,3,4,4)` | `-959/45` |
| `(1,3,4,4,4,3,4,3,4,5)` | `-4438/45` |
| `(1,3,4,4,5,3,4,3,4,4)` | `-2219/45` |
| `(1,3,4,4,5,4,3,3,4,4)` | `133/15` |
| `(1,3,5,4,5,4,3,3,3,4)` | `21/2` |
| `(1,4,4,4,4,2,4,4,4,4)` | `-7/3` |
| `(1,4,4,4,5,3,3,3,4,4)` | `-77/15` |
| `(2,2,4,3,5,4,4,4,3,4)` | `154/5` |
| `(2,3,3,3,4,4,4,5,3,4)` | `1673/180` |
| `(2,3,3,4,5,3,4,4,3,4)` | `-7994/45` |
| `(2,3,4,2,4,4,3,4,4,5)` | `-14` |
| `(2,3,4,3,3,4,4,3,4,5)` | `56/3` |
| `(2,3,4,3,4,4,4,3,3,5)` | `28/15` |
| `(2,3,4,4,4,3,3,3,4,5)` | `308/15` |
| `(2,3,4,4,4,3,4,3,3,5)` | `182/15` |
| `(3,3,3,2,4,4,4,4,4,4)` | `77/9` |

These coefficients do not prove that the full kernel is zero or nonzero.
They provide a compact complete identity for its remaining finite reduction.

## Reusable method

[`family_target.py`](../symmetry/family_target.py) builds only scalar orbits
reached from the supplied complete numerator. It reuses the exact diamond
and star--triangle rules in [`msys5.py`](../fiveloop/msys5.py), the inverse
star--triangle and position-space vertex rules in
[`rel5.py`](../fiveloop/rel5.py), and read-only finite terminals in
[`engine.py`](../symmetry/engine.py). It does not build a full period table.

The diamond moves preserve the extra-edge power and the power sums in the
three cubic columns `(0,3,4)`, `(1,5,6)`, `(2,7,8)`. The original source
polynomial occupies sixteen blocks of these invariants. Restricting the
vertex identities to extra-edge powers at most four keeps their shifted
terms in those blocks. Every accepted relation is complete; identities that
would need an untreated zero-power boundary contribution are rejected by
the inherited boundary guard.

The resulting local system has **47,092 variables and 163,907 rows**:

| Relation | Rows |
|---|---:|
| Known finite input periods | 2,277 |
| Diamond | 115,850 |
| Star--triangle | 10,008 |
| Inverse star--triangle | 11,398 |
| Position-space vertex identities | 24,374 |

Sparse exact elimination chooses a short row and a low-incidence pivot,
then updates both the remaining equations and the entire target. It takes
46,879 pivots and exhausts every active equation. The target still contains
19 terms, all stored in [target_certificate.json](target_certificate.json).
The identity above additionally uses the dimension constraints below.
Neither residual list is asserted to be a minimal independent basis under
every possible identity.

Input cores must be simple. Repeated undirected edges are rejected before
the inherited simple-graph orbit canonicalization; transformed multigraphs
inside accepted identities still receive the existing exact reductions.
Construction limits raise explicitly. Elimination limits return the whole
finite-plus-residual expression and its stop reason. Accessing `.value` on
an incomplete result raises `IncompleteFamilyTarget`.

## Checked dimension shift

For the same graph and projective measure, multiplying numerator and
denominator by `U**2` gives the exact identity

\[
I_{q-2}(\beta)
=\sum_m [t^m]U^2\;I_q(\beta+m).
\]

No gamma prefactor occurs. Each shifted term is integrable: it is a positive
summand of the original convergent scalar integrand.

The tracked `K33+e` table contains 687 suitable `q=5` entries: the extra-edge
power is three and every cubic-column sum is at most eight. Their `U**2`
shifts stay in the same sixteen source blocks. The optional dimension
certificate, [dimension_certificate.json](dimension_certificate.json),
records 61,130 variables, 215,623 rows, and 60,844 pivots. Every active
equation is exhausted within the declared bounds. This reduces the target
from 19 to **17 retained scalar periods**, with every coefficient displayed
above. The lift does not establish a complete finite value. Both compact
certificates preserve the full expression and make that remaining
calculation explicit.

## Reproduction and independent controls

```bash
python3 explore/circuit-general/forest-target-reduction/run_checks.py
python3 explore/circuit-general/forest-target-reduction/run_checks.py --examples
python3 explore/circuit-general/forest-target-reduction/run_checks.py --dimension-shift
```

The seven focused tests check all **23 `q=3` and 871 `q=4` local relation
rows** against the tracked complete tables, all **18 `q=3` to `q=5`
double-dimension identities**, exact bounded-target preservation, duplicate
input cancellation, canonical constant labels, explicit storage limits,
contradictory equations, and rejection of parallel input cores.

`certify_target.py` reconstructs all 10,577 original numerator coefficients
and all 7,368 scalar coefficients through the public source-complete kernel
pipeline. Their SHA256 references are fixed to the preceding contribution.
The compact certificates contain the original graph, every balanced source,
the full residual identity, exact counts and budgets, and implementation and
period-table hashes. Timings are confined to `CHECKS.log`; no intermediate
matrix or multi-megabyte numerator dump is needed. Default execution compares
every deterministic certificate field. `--write` explicitly refreshes the
chosen certificate.

The previous conformal-combination certificate, period engine, and public
forest dispatcher remain unchanged. This additional target solve is an
explicit computational route, with no automatic expensive table build.
