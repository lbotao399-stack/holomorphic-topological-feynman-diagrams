# Complete finite five-loop product-forest census

This continuation evaluates every one of the **200** classified five-loop
marked product-forest topologies at the admissible full source assignments
recorded in [certificate.json](certificate.json). It retains the complete bare
kernel and the first-current spinor on every original edge, including the
marked auxiliary edge: **4,400 exact stored entries**, or **2,200 full
spinors**. Every answer is a finite rational combination of
$1,\zeta(3),\zeta(5)$, evaluated from the existing exact tables.
There are no residual period symbols in this census.

| Complete bare result | Number of topologies |
|---|---:|
| Zero | 10 |
| Nonzero rational | 107 |
| Rational and $\zeta(3)$ coefficients | 77 |
| Rational, $\zeta(3)$ and $\zeta(5)$ coefficients | 6 |

All ten zero bare kernels have nonzero first-current spinors. The complete
current calculation therefore retains information absent from the bare
kernel. These displayed coefficient channels require no arithmetic
independence assumption about the zeta values.

[FULL_RESULTS.md](FULL_RESULTS.md) prints every bare kernel and every current
spinor in two original-source directions, with the lower spinor index retained.
These are geometric kernels in the original ordered-edge convention. The
source assignments are fixed inputs; arbitrary-source interpolation and full
Wick/color/coupling assembly are separate computations.

## Exact completion method

The existing sourced forest reduction, passive star moves, sourced ear moves,
and exact scalar-period tables are used unchanged. Six topologies, numbered
42, 47, 50, 117, 120 and 157 in the complete census, have no applicable local
forest move at these sources. Exact integration of their degree-two chains
leaves the graph obtained by deleting one edge from the complete graph on
five vertices. Its scalar integrals are the planar-dual triangular-prism
integrals already stored in `assembly/tables/prism_q8.json` for first-current
insertions, and the existing level-seven table for bare kernels.

Only the five original cotree currents need to be integrated independently.
The original grounded Kirchhoff equations express every other current as
the complete tree-supported source flow plus the signed fundamental-cycle
sum of those five currents. Multiplication by the full current-area form and
integration preserves this exact linear identity. Thus every original
first-current spinor is reconstructed from the ten independently integrated
cotree entries and the bare kernel. In particular, the affine source term
is never dropped, and a zero integrated bare value is not used to infer that
a current insertion vanishes.

The complete resulting data satisfy all **2,800** integrated source-conservation
identities, separately in each retained period channel. Coefficientwise
verification does not assume algebraic independence of zeta values.

The additive [ChainFirstPeriods](../../../code/htk/forest_chain_periods.py)
backend first sums the exact Dirichlet-factorial contributions with equal
chain powers, then applies the same graph automorphisms as the original
engine. This reassociates finite rational sums and retains every coefficient.
The independent [chain-first audit](chain_first_audit.json) compares the
complete old and new scalar dictionaries for a hard five-loop bare kernel
and both entries of its first-current spinor. Original hashed engines and
period tables remain unchanged.

## Independent controls and replay

The same current reconstruction agrees with all **1,018** entries in the
earlier complete original-graph reference through four loops, including an
additional check after reversing two currents and an odd edge permutation.
For **60 five-loop graphs**, [direct_reference.json](direct_reference.json)
also records all **1,320** entries integrated edge by edge without using the
cotree reconstruction. Those controls reuse the sourced reduction engine;
they are not described as independent original-graph integrations.

Each computation runs in a separate bounded subprocess. A timeout or a
resource failure leaves that graph unfinished and raises explicitly. A
complete census certificate is assembled only after all 200 full records
are present and every source-conservation check passes.

From the repository root, with the dependencies in
`../all-product-forests/requirements.txt` installed:

```bash
python3 explore/circuit-general/five-loop-forest-census/run_checks.py
python3 explore/circuit-general/five-loop-forest-census/certify.py \
  --checkpoint /tmp/five-loop-forest-replay --jobs 3 --timeout 900
```

The first command validates the complete stored census, independently
integrated hard-core values, every recorded audit dependency hash, and all
printed spinors without recomputing all 200 graphs. A missing full census
certificate is an error. Add `--deep` to regenerate both independent hard-core
audits. The second command recomputes every finite coefficient. The work directory
holds complete per-graph checkpoints; it can be reused after an interruption.
Pass `--direct --indices 0,6,7` with a fresh work directory to integrate every
original edge directly on the indicated graphs. The executable renderer
rechecks both entries of each printed original-source spinor.

The two immediate stored-data checks can also be run separately:

```bash
python3 explore/circuit-general/five-loop-forest-census/check_reconstruction.py
python3 explore/circuit-general/five-loop-forest-census/independent_core.py --compare-only
```

This finite census complements the universal convergent-series solution.
It does not assert a universal finite zeta basis for arbitrary loop number,
current degree, retained sources or interacting theory.
