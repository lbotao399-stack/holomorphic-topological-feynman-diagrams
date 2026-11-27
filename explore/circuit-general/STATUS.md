# Product forests and core periods: current evidence and remaining gates

Updated 27 November 2026. Scope: the product-forest and scalar-core work in
this repository, not its separate physical-theory and holography projects.

The repository's [100-point rubric](../../SCORE.md) now records **75 / 100**.
The accepted changes from 69.5 are C3, 3.5 to 4, for uniform integer-pole
period statements, and D2, 0 to 5, for the complete methods draft. The
[acceptance audit](../../notes/methods-paper/ACCEPTANCE.md) and
[final internal review](../../notes/methods-paper/verification/INDEPENDENT_REVIEW.md)
record the requirements and evidence. This ledger preserves the exact scope
of each result and its remaining gates. The admission rule and P1--P9
pipeline in [PLAN.md](../../PLAN.md) continue to govern the calculations.

## Completed, with exact scope

| Result | Accepted scope | Evidence and replay |
|---|---|---|
| Universal sourced network elimination | All classified product forests; original source and current transfer retained; 1,135 topologies through six loops checked | [All-product-forest reduction](all-product-forests/README.md) |
| Complete finite census through five loops | 261 topologies and 5,418 original first-current entries, at the explicitly recorded source assignments | [Five-loop census](five-loop-forest-census/README.md) |
| Universal exact scalar representation | Complete convergent rational series with rational remainder bounds, at rational sources and every finite current order; explicit resource failures | [Geometric cube expansion](forest-cube-expansion/README.md) |
| Arbitrary supported sources on an irreducible five-loop bare kernel | Ten full bare source coefficients; source vertices 0,2,4 structurally zero; no first-current claim for this topology | [Five-loop source polynomial](forest-finite-jets/README.md) |
| Arbitrary supported-source bare and first-current polynomials on the physical quartic five-loop forest | All eleven original spinors, 468 nonzero current/source monomials and 1,334 nonzero rational period coefficients; sources 0,3,4 structurally zero; denominator-free degenerate extension and complete fresh-frame edge integrations | [Full covariant source formulas](forest-source-currents/README.md) |
| Fixed-core stabilization | All normalized, absolutely convergent rational numerators and positive integer pole orders on a fixed graph; existential threshold and finite ambient period span | [Theorem and independent exhaustive controls](core-period-stabilization/README.md) |
| Explicit all-loop wheel scalar span | Every $n\ge3$, every integer $q\ge\max(3,n-1)$, full convergent scalar span exactly $\mathbb Q+\sum_{j=2}^{n-1}\mathbb Q\zeta(2j-1)$ | [Constructive wheel proof and exact certificate](core-period-closure/README.md) |
| Complete relative target certificates | General relation generator and bounded exact search; every successful specified target carries all boundary primitives; no numerical master independence claim | [Production engine and independent certificate replay](core-sector-reduction/README.md) |
| Complete methods draft | 35-page manuscript and 13-page complete coefficient supplement; normalization, derivations, exact controls, provenance, reproducible build and final internal review | [Methods paper and replay](../../notes/methods-paper/README.md) |
| Six-loop $K_5$ bare examples | All three recorded admissible subdivision shapes at specified balanced source assignments; complete core integration | [Complete reciprocal reconstructions](series-fibres/CORE_RESULTS.md) |

The stabilization theorem uses fixed relative cohomology and the complete
Kirchhoff identity. Its finite checks support the formulas; they are not an
extrapolation from consecutive equal ranks. The wheel theorem additionally
computes a specific stable span. A restricted physical source family can
span a proper subspace of the full scalar space.

## Next acceptance gates

| Open task | Required completion evidence |
|---|---|
| Extend the new full supported-source first-current result to other irreducible forests | Complete covariant tensor, all original edges and free spinor indices, every source coefficient, explicit degenerate-source extension, independent fresh-source checks for each additional family |
| Current insertion on the existing `forest-finite-jets` topology | Its suppressed core is $K_{3,3}$ plus an edge, not $K_5$ minus an edge; the first-current level-eight scalar family needs a complete new reduction or witness |
| Numerical stable basis of a general core | An explicit spanning proof or completed relative-cohomology calculation; finite tables or a rank plateau do not certify saturation |
| Prism, $K_{3,3}$ and $K_5$ three-zeta closure at all integer powers | Uniform reduction to the stated constants, including boundary values and exceptional powers; the general finite-span theorem alone does not close this gate |
| Scalable execution on large cores | Complete target identities with declared memory/time bounds and reproducible nonzero higher-loop targets; no residual dropped at a resource limit |
| Additional interacting products | Every Wick/color/coupling contribution in the stated theory and input family, followed by a complete independent contraction check |
| General $T=1$ current/contact transfer | Full topological contact terms and convergence control; no inference from the purely holomorphic integer-power theorem |

## Historical reports and superseded computational limits

The reports record their original calculation dates and are not silently
reinterpreted as current universal claims.

- The old level-eight insertion entry in [cores/REPORT.md](cores/REPORT.md)
  is superseded for the actual built families by
  [assembly/REPORT.md](assembly/REPORT.md) and the five-loop census. It does
  not mean every core has a level-eight table.
- The uncomputed six-loop numerator described in
  [k5/REPORT.md](k5/REPORT.md) is superseded for the three recorded shapes by
  [series-fibres/CORE_RESULTS.md](series-fibres/CORE_RESULTS.md).
- The earlier sufficient wheel threshold `q >= n+2` remains a valid theorem;
  the new wheel result proves the stronger sufficient range
  `q >= max(3,n-1)`.
- A complete geometric forest kernel is distinct from the full interacting
  Wick/color sum. None of the new scalar statements supplies an omitted
  interaction contribution or proves universal period cancellation.
