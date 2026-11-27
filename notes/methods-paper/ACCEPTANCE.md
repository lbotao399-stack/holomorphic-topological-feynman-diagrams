# Methods-paper acceptance and uniform-period audit

Review baseline: repository revision
`fac864229e7990310e9069244b5d2484bfc0a5e0`, 27 November 2026.
This note records an independent audit and prospective acceptance checks.
It does not change [SCORE.md](../../SCORE.md) or award points.

## The existing rubric

[SCORE.md](../../SCORE.md) records a total of **69.5**, gives credit only
for results with recorded exact checks, assigns half credit to `[drv]`
and none to `[open]`, and contains these two relevant rows:

| Item | Recorded score | Remaining requirement, verbatim |
|---|---:|---|
| C3, period theorems: value fields, reach criteria, vanishing | 3.5 / 4 | `uniform-in-q statements` |
| D2, methods paper: circuits and beyond | 0 / 5 | `draft` |

The [project admission rule](../../PLAN.md#1-admission-rule) additionally
requires an exact check independent of the derivation. These requirements
support a complete mathematical manuscript with reproducible exact examples.
They do not prescribe a page count, journal submission, journal acceptance,
new physical application, or the solution of every remaining graph family.
Those would be additional requirements, not the recorded D2 gate.

## C3: what the existing result establishes

Commit `75a9218` contains two actual uniform statements, with different
strengths and explicitly different scopes.

1. The [fixed-core theorem](../../explore/circuit-general/core-period-stabilization/DERIVATION.md)
   proves finite-dimensional rational span for every absolutely convergent
   normalized rational polynomial numerator on a fixed graph, across all
   positive integer powers of its Kirchhoff polynomial. The full scalar
   spans are nested and eventually equal. The threshold is existential.
2. The [wheel theorem](../../explore/circuit-general/core-period-closure/README.md)
   identifies the full scalar span on every wheel with `n >= 3` spokes as
   `Q + sum_{j=2}^{n-1} Q zeta(2j-1)` for every integer
   `q >= max(3,n-1)`. It gives a sufficient threshold and constructive
   spanning periods, without asserting arithmetic independence.

The independent audit found no blocking defect in the following proof steps:

- Greedy cotree factorization gives a fixed polynomial denominator between
  one and the spanning-tree count on each closed Hepp cube.
- Absolute convergence of a collected signed homogeneous numerator removes
  every negative integer coordinate exponent. Homogeneity makes the
  monomial substitution injective, so this also controls every surviving
  original scalar monomial.
- The finite complex of the cube's coordinate faces computes relative
  cohomology. Its top-degree quotient is by the vector subspace of all
  derivatives with primitives divisible by `x_i(1-x_i)`, retaining the
  complete boundary condition. Finite-dimensionality concerns this full
  relative quotient, not an implementation's selected finite relation bank.
- Multiplication by the complete Kirchhoff polynomial preserves convergence
  term by term and gives nesting. Finite-dimensionality then proves an
  eventual threshold; it does not identify that threshold from a rank plateau.
- The wheel lower inclusion retains the odd-zeta radial periods and supplies
  a nonzero rational direction. The finite four-vertex base range is joined
  to the exact uniform inequality for every `q >= 21`, rather than being
  extrapolated from sampled powers.

Both independent audit replays completed with exit status zero:

```bash
python3 explore/circuit-general/core-period-stabilization/run_checks.py
python3 explore/circuit-general/core-period-closure/verify.py
```

The first replay checked all 908 recorded sectors, 12,700 cotree
factorizations, 430,388 coordinate-exponent checks, 2,082 convergence-preserving
shifts, signed-numerator controls and complete relative identities. The second
checked 18 four-vertex base values, the uniform inequality, 10 rational wheel
controls, nine independent star--triangle embeddings and 10 spanning
determinants. These finite checks audit the exact formulas; the unbounded
claims follow from the written proofs.

**C3 recommendation:** the literal remaining requirement is satisfied by
these uniform-in-integer-pole statements. Recording the final 0.5 under C3
is defensible with precisely this evidence and scope. It would not be
defensible to describe that increment as all-power three-zeta closure for
the triangular prism, the complete bipartite graph on three plus three
vertices, or the complete graph on five vertices. Those questions remain
open. The complete-graph all-power problem is also explicitly a separate
remaining requirement under C2.

Neither theorem includes arbitrary noninteger powers, regulator expansions,
changing denominator divisors, unboundedly varying cores, or an asserted
numerical basis for every physical source/insertion subfamily.

## D2: acceptance of a complete draft

The following checks explain what constitutes the requested **draft** in
the existing exact-check standard. They do not create a new point schedule.

| Check | Evidence needed before accepting the draft |
|---|---|
| Complete manuscript | A readable paper, with abstract, introduction, consistent conventions, connected methods/results sections, proofs or precise primary references, worked examples, limitations and bibliography. No missing section or unresolved placeholder needed for a principal claim. |
| Circuits and beyond | Explain the sourced current equations and exact local/network reduction, then explain the treatment of remaining scalar periods. State how these methods compose and what each computational route returns. A list of links to separate repository reports is insufficient as the paper itself. |
| Exact mathematical scope | State graph and convergence assumptions, scalar normalization, source support, original edge orientations and free-index conventions before using them. Separate complete network reduction, a convergent exact series, a finite evaluated period, and a complete interacting Wick/color sum. |
| Proof integrity | Include the principal arguments claimed by the paper, including retained affine source transfer and current insertions, relevant boundary terms, and the hypotheses of any uniform period theorem. Identify every inherited theorem and its source. A methods draft may synthesize established results; the rubric does not require an additional new theorem. |
| Exact worked results | Present fully evaluated examples with a recorded exact check independent of their production derivation. Preserve the entire tensor and every declared coefficient. If an example has restricted sources or is only a scalar control, say so explicitly. Report unresolved residuals and computational limits honestly. |
| Reproducibility | Supply manuscript source, an executable build command and a successfully compiled PDF. Supply a documented check command, recorded successful output and dependency/source provenance sufficient to bind the checked data to the manuscript's examples. Distinguish cheap stored-certificate replay from expensive regeneration. |
| Independent review | Record a separate review of the assembled manuscript's mathematical statements, normalization, citations, example-to-certificate agreement and computational scope. Resolve every finding that invalidates a principal claim. Explicitly stated open problems are not defects requiring closure for draft acceptance. |

A successful build alone does not establish these checks. Conversely,
new numerical experiments, a required page count or a journal decision
are not necessary once a complete rigorous draft and its exact evidence
meet them. All existing finite computations need not be repeated from
scratch merely to prepare the manuscript; their preserved certificates,
independent replays and declared provenance may be used.

## Scoring interpretation and current decision

Under the recorded wording, a complete independently reviewed draft meeting
the above checks can justify **D2 = 5 / 5**. A lower score would be appropriate
if the artifact remains an outline, lacks an essential proof, has unreviewed
principal claims, or relies only on derived assertions without the required
exact checks. This note invents no fractional sub-rubric for those incomplete
states. The `[drv]` rule remains applicable when that is the actual evidence
status.

At the time this acceptance note was written, the assembled methods manuscript
had **not yet passed final independent review**. D2 is therefore **pending**;
no draft credit is awarded by this note. The proposed arithmetic
`69.5 + 0.5 + 5 = 75` is a conditional consequence of accepting both results,
not a reason to accept either one. Any later score change should cite the
completed manuscript, build/check outputs and final review explicitly.

Finally, at the audit baseline the introductory paragraph in
[the circuit status ledger](../../explore/circuit-general/STATUS.md)
stated that no numerical rubric existed. That description is superseded by
the present `SCORE.md`; it must not be reused as the current scoring policy.

## Final disposition, 27 November 2026

**The methods-paper draft is accepted.** This final disposition supersedes
the pending status recorded above, while retaining the original criteria
and scoring interpretation unchanged. The separate
[independent review](verification/INDEPENDENT_REVIEW.md) found no unresolved
blocking mathematical or provenance error and closed the build, layout
and exact-artifact gates.

The completed paper and supplement have 35 and 13 pages. The final full
replay passed all 12 commands and all 23 tests without skips. The reviewer
independently checked all 93 computational input hashes, all 24 command
log hashes, both PDFs and their build records, all 48 final render hashes,
and the complete layout record, with a separately disclosed visual sample.
The exact worked examples retain every declared coefficient and all eleven
original current spinors. The accepted scope and inherited-computation
provenance are explicit; this acceptance does not claim fresh reintegration
of every stored census result or external peer review.

The final review records the complete artifact hashes. In particular,
the accepted full replay manifest has SHA-256
`b3cea0fde11a88d2943d3e5747a0736f9498ab92911e6578cbac598df0cb3f15`,
the build record
`06f7e876291a6d13846d10742d474bb8e0fc8b51b0ba96d306446373313d53f0`,
and the layout record
`0b29bc23ca135fa7736089607d071a1f2adf106131db23c6cb5a387c665baf2d`.

Under the existing literal gates, the completed evidence supports
**D2 = 5 / 5** for the draft and the remaining **0.5 under C3** for the
scoped uniform integer-pole statements. The corresponding arithmetic is

```text
69.5 + 0.5 + 5 = 75.
```

This is an acceptance recommendation grounded in the completed artifacts;
the score ledger is a separate edit. No requirement was removed or weakened
to obtain the target, and the open general-core, unsupported-source and
physical-assembly questions remain open as stated in the manuscript.
