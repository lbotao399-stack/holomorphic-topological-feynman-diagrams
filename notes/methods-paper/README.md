# Exact sourced circuits and fixed-core periods

This directory contains a complete methods draft, its complete coefficient
supplement, manuscript-specific independent checks, and reproducible build
records. The source files are [paper.tex](paper.tex) and
[supplement.tex](supplement.tex); the built documents are
[paper.pdf](output/pdf/paper.pdf) and
[supplement.pdf](output/pdf/supplement.pdf).

## Mathematical and computational scope

The paper derives source-preserving elimination on every product forest,
reconstructs every original current with its affine source terms, fixes the
current-area orientation, and proves a complete convergent rational cube
expansion at every finite current order with an explicit rational tail.
It then proves eventual stabilization of the full scalar period span at a
fixed graph, with an existential threshold. For wheels with $n\ge3$ spokes,
the full convergent scalar span is

$$
\mathbb Q+\sum_{j=2}^{n-1}\mathbb Q\,\zeta(2j-1)
\qquad\text{for every integer }q\ge\max(3,n-1).
$$

No arithmetic independence of these constants is assumed. The fixed-core
theorem does not determine a general numerical basis or stabilization
threshold, and no all-power three-zeta closure is claimed for the prism,
the complete bipartite graph on three plus three vertices, or the complete
graph on five vertices.

The supported five-loop example retains arbitrary spinors at vertices
1, 2 and 5, with $\lambda_6=-\lambda_1-\lambda_2-\lambda_5$ and zero
sources at vertices 0, 3 and 4. The supplement prints **484 complete source
monomial rows: 16 bare and 468 current rows**. The current rows contain
**1,334 nonzero rational coefficients** in the channels
$1,\zeta(3),\zeta(5)$ and give **all eleven original current spinors**,
including degenerate source frames. The displayed free spinor index is
retained. These are normalized geometric kernels before Wick, color,
representation, coupling and propagator-normalization factors.

The finite census, the complete supported-source polynomial, the all-order
rational series, and the uniform scalar-period theorems have distinct scopes.
The paper states those scopes separately and supplies full triangle, fork
and relative-boundary controls.

## Dependencies

Use Python 3.11 or later with **SymPy, python-flint and pypdf** installed.
The bounded runner uses the standard-library unittest interface for the
existing focused tests; **pytest is optional**. An equivalent optional
invocation for the relative-engine tests is:

~~~bash
python3 -m pytest -q code/tests/test_core_sector_periods.py
~~~

Building also requires **pdfLaTeX, BibTeX and kpsewhich**, with the installed
TeX packages used by [preamble.tex](preamble.tex), including Latin Modern.
The build script uses installed resources only. If the TeX tree is present
but its generated format or font map is missing, it creates local files
under the manuscript's temporary directory. It does not download packages
or change a system-wide TeX installation.

## Build and full bounded replay

Run these commands from the repository root:

~~~bash
python3 notes/methods-paper/build.py
python3 notes/methods-paper/build.py --verify
python3 notes/methods-paper/run_checks.py
~~~

The build compiles both documents and rejects unresolved references,
overflowing boxes and empty pages. Its
[build record](verification/build.json) records the exact manuscript
source hashes, compiler version, PDF hashes and page counts. The verify
command checks an existing build without compiling it again.

The bounded runner executes eleven mathematical and data checks, followed
by verification of the two existing PDFs. It records each command,
return code, elapsed time, structured coverage counts, stdout and stderr
logs, installed versions, and SHA256 hashes of every selected source and
certificate dependency. Its result is
[verification/checks.json](verification/checks.json), with separate logs
under [verification/check_logs](verification/check_logs).
Inputs are hashed before and after replay; a changed input is an error.
The output manifest and logs are never included in their own source hashes.

Each command has a 300-second timeout. Any failed command, missing required
input, incomplete focused-test count, skipped focused test, changed input,
missing PDF or failed build verification makes the runner exit nonzero.
An explicitly requested timeout can be supplied as, for example,
<code>--timeout 600</code>.

To check the mathematics while a document is being prepared:

~~~bash
python3 notes/methods-paper/run_checks.py --math-only
~~~

This writes separate
[verification/math_checks.json](verification/math_checks.json) and
<code>verification/math_check_logs/</code>. Its success status is
<code>math-only-passed</code>, with <code>full_replay_passed=false</code>;
it cannot satisfy the PDF/build gate. A successful full replay establishes
the listed computational and document bindings. Mathematical review,
primary-source review and visual page inspection are recorded separately;
the runner neither awards a score nor substitutes for those reviews.

## Individual exact checks

| Command, from the repository root | What is checked |
|---|---|
| <code>python3 notes/methods-paper/generate_tables.py</code> | All four generated TeX tables agree with their complete exact data. |
| <code>python3 notes/methods-paper/checks/printed_tables.py</code> | Independent parsing of all 484 printed monomials, all 1,452 displayed rational channel entries, eleven specialized spinors, eighteen wheel bases and sixteen complete boundary primitives. |
| <code>python3 notes/methods-paper/checks/fork_control.py</code> | Full six-variable current Jacobian, all seven fork spinors, ten polynomial and ten integrated Kirchhoff identities, and two original-network orientation controls. |
| <code>python3 notes/methods-paper/checks/wheel_residue.py</code> | Independent rational residue calculation of every base from $q=3$ through $20$, including the complete polynomial part, plus the exact inequality supporting every $q\ge21$. |
| <code>python3 explore/circuit-general/core-period-stabilization/run_checks.py</code> | Independent exact finite controls of the fixed-core proof, including 908 sectors, 12,700 cotree factorizations, 430,388 coordinate-exponent checks and 2,082 convergence-preserving shifts. |
| <code>python3 explore/circuit-general/core-period-closure/verify.py</code> | Eighteen wheel base values, the uniform inequality, ten rational-wheel controls, nine star--triangle embeddings and ten spanning determinants. |
| <code>python3 explore/circuit-general/core-sector-reduction/run_checks.py</code> | Ten focused production tests and eight stored complete target identities, independently differentiated with all 82 primitive end faces retained. The manuscript runner executes these two constituent scripts separately. |
| <code>python3 explore/circuit-general/forest-source-currents/run_checks.py</code> | Complete stored supported-source reconstruction, all 644 frame coefficients, fourteen polynomial source identities, all eleven spinors, fresh and degenerate controls, dependency hashes and the full printed result. |
| <code>python3 code/tests/test_forest_cube_expansion.py</code> | Seven focused tests of complete rational coefficients, certified tails, original sources, orientations and explicit budgets. |
| <code>python3 code/tests/test_solve_product_forest_cubes.py</code> | Six focused tests of the complete tensor interface, original edge reversal, exact enclosures, failed requests and atomic output. |

The uniform statements follow from the written proofs. These finite controls
check the formulas and certificates used in those proofs; no observed rank
plateau or finite sample range is treated as an all-order argument.

## Regeneration and expensive original integration

The normal runner replays stored source and scalar certificates. It freshly
recalculates the small manuscript controls, but **does not reintegrate the
five-loop interpolation samples or rebuild the original scalar period
tables**. The current manuscript replay therefore makes no claim that
those expensive calculations were repeated during draft preparation.

To regenerate the complete printed tables after an intentional data change,
and then refresh their independent parsed record:

~~~bash
python3 notes/methods-paper/generate_tables.py --write
python3 notes/methods-paper/checks/printed_tables.py --write
~~~

The fork and wheel-residue checkers also accept <code>--write</code> to
regenerate their exact control records. Inspect a mismatch before replacing
a certificate. Rebuild the PDFs and rerun the full checks after any
manuscript or generated-table change.

The explicitly expensive supported-source command is:

~~~bash
python3 explore/circuit-general/forest-source-currents/run_checks.py --deep --jobs 4
~~~

It recomputes all sixteen directly integrated interpolation nodes, constructs
the complete 28-node data, and reintegrates both independent original-edge
controls. It still uses the existing exact scalar period tables. This
command is intentionally absent from the bounded manuscript runner.

The complete earlier cube-certificate regeneration/replay is also separate:

~~~bash
python3 explore/circuit-general/forest-cube-expansion/run_checks.py
~~~

That command recalculates every recorded three-loop cube and its recorded
series coefficients. The bounded manuscript runner instead executes the
two focused cube test files listed above.

See [ACCEPTANCE.md](ACCEPTANCE.md) for the review criteria and
[verification/PRIMARY_SOURCES.md](verification/PRIMARY_SOURCES.md) for
primary-reference provenance. Build and check records describe what was
actually replayed; they do not assert completion of the open physical or
all-core classification problems.
