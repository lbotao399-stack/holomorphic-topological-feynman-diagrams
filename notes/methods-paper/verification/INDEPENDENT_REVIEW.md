# Independent mathematical and provenance review

Date: 27 November 2026. Repository baseline:
`fac864229e7990310e9069244b5d2484bfc0a5e0`.

Reviewer: the separate `audit_reduction` artificial-intelligence agent.
The reviewer did not draft the assembled manuscript or its network and
period sections. The reviewer previously implemented the underlying
supported-source current study and its verifier, so this review is not
an independent reintegration of that entire earlier study. Independent
formula checks, original-network controls and the manuscript's separate
printed-table parser are identified below. This is an internal technical
review, not external peer review.

## Final decision

**Accepted as a complete methods-paper draft under the existing D2 gate.
No unresolved blocking mathematical, provenance, build or layout finding
remains in the exact version identified below.**

The manuscript meets the acceptance criteria in
[ACCEPTANCE.md](../ACCEPTANCE.md). The initially pending build, layout and
final-manifest checks were completed on 27 November 2026; the final gate
record below supersedes that initial checkpoint. Full D2 draft credit is
supported by this review, which does not itself edit or award a score.
The remaining C3 requirement is satisfied only in the stated uniform
integer-pole scope; this does not close the separate all-power three-zeta
questions for general cores.

## Actual findings and their disposition

| Finding | Consequence | Disposition |
|---|---|---|
| Several TeX commands in the first reviewed period section lacked their initial backslash, including the de Rham subscript and spacing commands in the radial/embedding formulas. | They would render as unintended mathematical variables even if compilation succeeded. | Reported directly to the section author. The missing escapes were corrected; the wheel certificate was refreshed and its independent replay passed afterward. |
| The first introduction described the starting normalization as that of Li--Liu--Chang without immediately distinguishing the repository's loop-dependent orientation sign. | A reader could incorrectly apply the same sign to the two kernel conventions. | Root revised the introduction to name the repository convention and point to the explicit conversion. The conversion was independently checked against the original source below. |
| At the initial checkpoint, final compiled-document/layout records and the final manuscript manifest had not been reviewed. | Buildability, final cross-references, visual legibility and the exact accepted source/PDF binding were not yet established by this review. | Resolved by the final gate below: both PDFs, all 93 input hashes, all 24 command logs, the complete layout record and all 48 final page-render hashes were checked. |

The review did not find an incorrect period value, omitted source direction,
dropped original current, missing convergence assumption, or unjustified
finite-table extrapolation in the reviewed version.

## Mathematical audit

### Original currents and normalization

The incidence convention is positive at the initial vertex and negative
at the terminal vertex. With resistances rather than conductances as
the edge variables, the stated grounded solution and its derivative give

\[
\frac{\partial y_e}{\partial t_f}
=\left(\frac{B_e\Lambda^{-1}B_f^{\mathsf T}}{t_e}
             -\delta_{ef}\right)\frac{y_f}{t_f}.
\]

The source is held fixed in this differentiation. The Schur transfer
retains its resistance dependence and the affine source term. The chain
currents and both endpoint source transfers satisfy their original
Kirchhoff and voltage equations; no derivative of an accumulated source
is silently frozen by the general rule.

Regrouping the two source blocks of the cycle form into complete chord
spinors contributes `(-1)^(n(n-1)/2)`. In
[the local original manuscript](../../../sources/LLC_main.tex), equation
`eq-operatope` has one additional minus sign per chord, while
`eq-UV-boundary` fixes the same positive radial orientation. Therefore

\[
\mathcal I^{\rm LLC}_\Gamma
=(-1)^{n+n(n-1)/2}K_\Gamma
=(-1)^{n(n+1)/2}K_\Gamma.
\]

The paper states this conversion explicitly. The triangle bare value,
all three triangle spinors and the negative displacement derivative
convention agree with it. Reversal of an original edge changes that
current's sign and leaves the bare two-dimensional kernel unchanged.
The full fork calculation independently checks the determinant sign
against two positive original-resistance controls, including a nonunit
source frame.

### Compact cubes and the complete rational series

The minimum-spanning-tree inequalities give the stated atlas up to
measure-zero equality walls. For each fixed chord order the tree-edge
upper bounds are exactly those in the cube map. The prefix cycle-rank
argument makes every cotree ratio a polynomial with nonnegative powers.
Its selected cotree is the unique constant term; the upper bound is the
spanning-tree count. Every original current consequently extends to the
closed cube, which establishes absolute convergence of every finite
current word in the stated two-dimensional current-area setting.

The augmented determinant has denominator exponent `2n+1`, before
insertions or verified exact cancellations. Its row operation supplies
the complete current Jacobian. The centered negative-binomial expansion
has ratio bounded by `(tau-1)/(tau+1)`. The finite tail formula follows
from the stated negative-binomial tail identity and the absolute
polynomial moments. The resulting finite truncation is correctly
described as a partial sum with an explicit rational bound; the exact
answer remains the full series specification.

The polynomial-numerator denominator `U^(n+2)` follows from the rank-one
adjugate along each simple irreducible divisor and the resulting `n`
block-row cancellations. The apparent omitted-edge factor is removed by
the alternative cofactor identities. The stated time and source degrees
then follow by homogeneity. These arguments concern the complete density;
they do not substitute an intermediate nontight graph's kernel.

### Fixed-core finiteness and stabilization

The proof keeps a fixed denominator divisor and all coordinate faces.
For a collected signed numerator, the least negative integer coordinate
power causes absolute divergence on a positive boundary patch. The
homogeneous monomial substitution is injective, so the argument also
justifies termwise scalar convergence.

The finite face complex gives the top relative quotient by the **vector
subspace** of derivatives of `x_i(1-x_i)` times regular functions. The
smooth affine face intersections permit the cited de Rham comparison and
finiteness argument. No ordinary ideal reduction is substituted for this
vector-space quotient. Multiplication by the full Kirchhoff polynomial
embeds a scalar level into the next with every summand separately
convergent. The finite-dimensional containing space then proves eventual
stabilization.

The examples correctly distinguish an existential threshold from a
finite observed plateau and relative cohomology from the numerical image
of the original graph family. The proof does not identify a general
numerical basis, assert independence of zeta values, or cover variable
divisors, regulator logarithms or noninteger powers.

### Uniform wheel theorem

The expanded upper-inclusion proof was reviewed, including the
position-space gamma normalization and squared-radius gauge, angular
multiplicities, complete link identity, parity of the raw residue nodes,
the common strip, and oddness of the convolution function. In particular,
the squared-radius delta supplies the stated factor `1/2`; the link
rescaling cancels the remaining powers of two. The pole locations lie
strictly below the harmonic summation range. Multiplicity zeros reduce
the maximal pole order, oddness cancels the even zeta coefficients, and
convergence excludes both the polynomial part and a harmonic divergence.

The lower-inclusion proof retains the triangular radial odd-zeta
directions and constructs a nonzero rational direction. Its finite
four-vertex base range is complete through `q=20`; the explicit lower
bound and positive growth polynomial prove every `q>=21`. Hence the
threshold `q>=max(3,n-1)` is genuinely uniform. The proof does not require
numerical independence of the displayed constants.

### Complete source-current family

The graph has the eleven stated ordered edges. The support restriction
is explicit: sources at `0,3,4` vanish, the three independent spinors
are at `1,2,5`, and vertex `6` carries their negative sum. The scalar
denominator levels are seven and eight after complete chain integration
to the five-vertex complete graph with one edge absent. Planar duality
with complemented edge powers gives the existing prism tables.

Source degrees ten and eleven, the full `GL(2)` covariance, and polynomial
regularity prove the frame degree bound before interpolation. The
degree-six identity `b f_6 = a g_6` gives the denominator-free lift in
all three original source directions. Its equality extends to degenerate
frames as a polynomial identity. The 28 triangular nodes, sixteen direct
nodes, twelve symmetry images and 644 frame coefficient dictionaries
are the correct complete counts.

All eleven currents are kept. Tree currents include the bare kernel
times the original source flow. The fresh and degenerate controls
integrate both entries of every edge separately. The text correctly
discloses their shared numerator formula and period tables. It also
distinguishes this family from the other unresolved five-loop core and
from a fully assembled physical product.

## Exact checks performed during this review

The following commands completed successfully in the reviewer's session:

```bash
python3 notes/methods-paper/checks/printed_tables.py
python3 notes/methods-paper/checks/fork_control.py
python3 notes/methods-paper/checks/wheel_residue.py
```

- The printed-table parser checked all **1,452 rational entries** in the
  complete bare/current supplement, including every zero column. It
  verified 16 nonzero bare monomials, 468 nonzero current monomials and
  1,334 nonzero rational current coefficients; all eleven specialized
  spinors; 42 specialized coefficient-conservation checks; all 18 wheel
  base rows; and all 16 printed boundary primitives with 32 vanishing
  end faces. It imports neither the table generator nor the production
  current/period implementation.
- The fork checker recomputed the full six-variable chord-current
  determinant, all seven spinors, ten polynomial and ten integrated
  source identities, and two original-network orientation controls. Its
  symbolic moment calculation is separate from the production forest
  reduction; the original-network point controls use the existing
  brute-force Kirchhoff engine, which is disclosed in its dependencies.
- The wheel checker imports no production residue evaluator. It checked
  18 complete bases, 1,116 proper-subset convergence inequalities, 396
  principal coefficients, 1,458 spectral coefficients and all 18 complete
  polynomial-part cancellations, together with the uniform bound.

The preceding rubric audit also replayed both existing fixed-core and
wheel certificates successfully. Finite checks support the exact formulas;
the uniform theorems rely on their written proofs, not extrapolation.

The finite-census counts were independently read from their complete
certificates: 1,135 network topologies; 61 integrated topologies and 1,018
current entries through four loops; 200 five-loop topologies and 4,400
current entries. Thus the stated totals are 261 and 5,418, with 2,200
complete five-loop spinors. Their recorded-source scope is preserved.

## Provenance and attribution

The cited source revision resolves exactly to
`75a9218178da47e7926df75ef29ee86f4616e421`. The complete source-current
certificate used by the printed parser has SHA-256
`d16bc6c84383ba67d4dc80ecf4eee0b4c9efa339f235c2439643cde98dc74b98`.
The paper presents these as inherited repository calculations, not new
integrations completed while drafting.

The attribution to the all-order Type I-prime recursion was checked
directly in [Gui--Wang--Williams, version 2, Theorem 4.13 and Remark 4.14](https://arxiv.org/html/2510.26608v2).
The cited theorem is not restricted to the three explicitly calculated
loop orders. The present introduction acknowledges that correctly.
The manuscript also credits the inherited current-area construction,
holomorphic finiteness, relative graph-period framework and algebraic
de Rham comparison. It distinguishes its complete source retention and
specified computational examples from physical Wick/colour assembly.

The local Li--Liu--Chang source used for the sign check has SHA-256
`97c470530a1965a7bb1fcf7dfa28e5c02dccf19d1879e459be54775af05032bb`.
The accepted final manifests bind the exact manuscript and generated
table inputs, as checked below. The primary-source record distinguishes
locally supplied source material from sources accessed through their
public original versions.

## Completed final gate

The reviewer read the final `README.md`, `run_checks.py`, `build.py`,
`PRIMARY_SOURCES.md`, build record, replay manifest and layout record.
The documented full command is

```bash
python3 notes/methods-paper/run_checks.py
```

The final full replay was run by the coordinating agent after the final
page-break changes. This reviewer independently audited its actual
recorded outputs and current files without rerunning the command and
thereby replacing the frozen manifest:

- All **93 computational input hashes** match the current files; all
  **24 standard-output and standard-error log hashes** match their logs.
  All **12 commands** completed successfully. The three test-suite logs
  contain exactly **10 + 7 + 6 = 23 passing tests**, with no skips.
- The runner checks the expected test counts, removes an inherited
  `PYTHONOPTIMIZE`, enforces bounded subprocess runs, and compares input
  hashes before and after replay. Its mathematical-only mode cannot
  claim a full build pass. The documented full mode includes PDF build
  verification. Missing inputs, failed commands and incomplete test
  counts fail the replay.
- Every final build-source hash and each recorded build artifact hash
  matches. The main PDF has **35 pages** and the supplement **13 pages**,
  independently confirmed from the actual files. Every page has
  nonempty extractable text. The final compiler logs contain no
  overflowing boxes, unresolved references or citations, or
  multiply-defined labels. The build runs the compiler and bibliography
  passes and rejects those failures; verification binds the existing
  PDFs to their sources and page counts.
- The layout record binds the same final manifest, build record and
  PDFs. All **48 final rendered-page hashes** were independently checked.
  The coordinating agent records inspection of all 48 rendered pages,
  including the final changed pages. This reviewer separately viewed
  main-paper pages **2 and 31--35** and supplement pages **3 and 13**.
  No clipped formula, overlapping fraction row or stranded table label
  was observed in that separate sample. This report does not represent
  that sample as an independent full-page visual review of all 48 pages.
- The final page-break and display-layout corrections preserve the
  mathematical claims, exact table contents and proof scope reviewed
  above. The bounded replay explicitly does **not** claim to repeat the
  expensive original five-loop integrations or regenerate the inherited
  period tables. Their certificates, separate controls and computational
  provenance remain identified.

The accepted files have the following SHA-256 hashes; paths are relative
to `notes/methods-paper/`.

| Artifact | SHA-256 |
|---|---|
| `verification/checks.json` | `b3cea0fde11a88d2943d3e5747a0736f9498ab92911e6578cbac598df0cb3f15` |
| `verification/build.json` | `06f7e876291a6d13846d10742d474bb8e0fc8b51b0ba96d306446373313d53f0` |
| `verification/layout.json` | `0b29bc23ca135fa7736089607d071a1f2adf106131db23c6cb5a387c665baf2d` |
| `verification/PRIMARY_SOURCES.md` | `6c49d9f11719bdea5aebf23c88c23316a51785284a87005f47a862646cb0fe43` |
| `output/pdf/paper.pdf` | `3c72970c8ceb94a321aaa8c0d7edcbfe8b28e468cef1e60849e66a9e740e30e1` |
| `output/pdf/supplement.pdf` | `e39beb85f0b65e4203cfeb06a986f95c5b3f1813266482931bf9553f4a35716d` |

The mutable review, acceptance, layout and score reports are excluded
from the computational-input manifest intentionally, avoiding circular
hashes. This review binds the final layout report explicitly above.

All seven recorded D2 acceptance checks are now met. A complete draft
therefore supports **D2 = 5 / 5** under the unchanged rubric. The scoped
recommendation for the remaining **C3 = 0.5** is also unchanged. Neither
conclusion asserts external peer review, a universal finite-zeta basis,
unrestricted source support for the five-loop polynomial family, or a
completed Wick/colour assembly. No score file was changed by this review.
