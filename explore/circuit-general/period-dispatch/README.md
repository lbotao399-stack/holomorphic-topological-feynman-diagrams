# Existing K5 periods in the public reduction pipeline

`k5_bridge.py` exposes the tracked `../k5/tables/K5_q*.json` certificates in
the same `I_q` normalization as `per5.period`. It does not construct tables
or solve new relations. `per5.lower` and `per5.value` now use this bridge for
complete five-vertex cores at `q = 3, ..., 8`.

| q | Available weighted K5 orbits |
|---|---:|
| 3 | 2 |
| 4 | 30 |
| 5 | 312 |
| 6 | 2,107 |
| 7 | 10,512 |
| 8 | 41,428 |
| Total | 54,391 |

These are existing certificates, including the generic orbits. Their
derivation and exact relation checks are recorded in
[`../k5/REPORT.md`](../k5/REPORT.md), sections 1--6. The bridge checks the
file format, exact rational coefficients, strict convergence, and distinct
isomorphism orbits. It records the source file and SHA-256 fingerprint; it
does not repeat the original relation solve.

## Exact lookup and bulk contraction

```python
from fractions import Fraction
from k5_bridge import E5, K5Lookup, k5_period_with_provenance

lookup = K5Lookup(E5, 5)
assert lookup.period([3] * 10) == {
    "1": Fraction(-13, 15120),
    "z3": Fraction(19, 720),
    "z5": Fraction(-5, 168),
}

# Accumulate numerator coefficients before any period arithmetic.
grouped = {}
for powers, coefficient in weighted_numerator.items():
    if coefficient:
        key = lookup.key(powers)
        grouped[key] = grouped.get(key, 0) + coefficient

answer = {}
for key, coefficient in grouped.items():
    if coefficient:
        for label, value in lookup.values[key].items():
            answer[label] = answer.get(label, Fraction(0)) + coefficient * value

certificate = k5_period_with_provenance(E5, [3] * 10, 5)
```

`E5` is the order `01,02,03,04,12,13,14,23,24,34`. Any edge ordering,
orientation, and consistently ordered vertex labels are accepted. Input
weights are positive integers or exact integral `Fraction` values. Labels
`1`, `z3`, and `z5` mean the rational part, zeta(3), and zeta(5).

`lookup.key(powers)` validates convergence and returns a collision-free
orbit key; `lookup.values` contains immutable exact coefficient mappings.
`lookup.metadata` identifies the table. `k5_period_with_provenance` also
returns its original orbit representative and certificate kind
(`T1`, `T2`, or `G`, defined in the K5 report).

The three failure cases are explicit:

- `K5NotConvergent`: the positive projective integral fails a strict
  convergence or homogeneity condition.
- `K5Unavailable`: the selected table or convergent orbit is absent.
- `K5TableError`: a present table is malformed or inconsistent.

None returns an empty coefficient dictionary. The public `per5.period`
adapter preserves its existing `ValueError` for divergence and raises
`UnknownCore` for missing convergent periods; `per5.value` returns `None`
when the integral diverges. Outside the available K5
levels it retains the existing reduction attempts. The bridge never
starts an automatic table build; this does not change the older builders
for other core families.

## Canonicalization and convergence

A vertex signature is the sorted four-tuple of its incident weights.
Sort vertices by these signatures, then minimize the edge-weight tuple
only over permutations inside equal-signature blocks. An isomorphism
preserves every signature and therefore maps these blocks bijectively.
The minimization is consequently a complete canonical form. It uses one
permutation when all five signatures differ; the full 120 permutations
remain available only for tied signatures. Tables are indexed once with
this convention. Their original full-permutation representatives are
retained as provenance, and `canonical_powers` remains an independent
full-permutation reference implementation.

For K5, strict projective convergence is equivalent to

```text
sum of all ten weights = 6q;
0 < each edge weight < q;
each triangle has total weight > q;
each four-vertex clique has total weight > 3q.
```

To see sufficiency, complete each connected proper subgraph on at most
four vertices to its clique. Adding an edge increases `q * loops - sum`
by `q - edge_weight > 0`, so the listed clique inequalities imply the
inequality for the original subgraph. For a connected proper subgraph on
all five vertices, the complementary missing edges and `sum = 6q` imply
the inequality directly. Disconnected subgraphs add these bounds.
Necessity follows by testing triangles, four-vertex cliques, and K5 with
one edge removed. Thus fifteen clique checks replace enumeration of all
1,022 proper nonempty edge subsets in each bulk lookup.

## Verification

```bash
python3 explore/circuit-general/period-dispatch/test_k5_bridge.py
python3 explore/circuit-general/period-dispatch/benchmark_lookup.py --count 100000
```

The tests compare every homogeneous positive weight vector at `q=3,4`
with an independent all-subgraph convergence check, compare the new
canonical form with all 120 vertex permutations, test random relabelings
at `q=8`, and check seven exact reference periods across all available
levels. They also exercise public dispatch with table construction
disabled, missing-data failures, and exact contraction after grouping
integer coefficients by orbit. Period availability is distinct from
completion of a source-dependent graph numerator.
