# K5: an exact dimension-lowering relation and the missing input identity

This directory adds three mathematical results to the existing finite K5
tables.

1. A downward dimension identity, valid at every integer `q >= 3` as one
   convergent integral with its complete Gram numerator retained. Its
   expansion into separate completed periods is permitted precisely under
   the explicit strict-cut guard proved below.
2. A derivation of the previously unexplained relation between the four
   tower inputs. It comes from a uniform radial-wheel identity, one exact
   finite boundary term, and a certificate consisting of eleven vertex
   relations and one upward dimension relation.
3. Complete exact values at the new denominator level `q = 9` for
   `3 <= beta_e <= 7`, `sum beta_e = 36`, where `alpha_e = 9-beta_e` on
   the ten original K5 edges. This is **4,905 labeled inputs and 73
   isomorphism orbits: 13 generic, 7 T1, and 53 T2**. All 9,490 terms in
   their downward relations are evaluated against the existing complete
   `q = 8` tables.

**The all-integer-q numerical closure of arbitrary K5 periods remains
open.** The guard fails on explicit convergent inputs. Neither the new
values nor the general identity is presented as a proof that the existing
four-input reduction closes at every level.

## Exact result beyond the previous tables

Use the existing edge order

```text
01,02,03,04,12,13,14,23,24,34.
```

One of the thirteen newly evaluated generic orbits is

\[
\begin{aligned}
I_9(K_5;3,5,6,6,5,6,6,6,6,5)
={}&\frac{1637270809}{3295339407360000}\\
&+\frac{50430067}{22884301440000}\zeta(3)
-\frac{551}{181621440}\zeta(5).
\end{aligned}
\]

Every one of the 73 complete answers is in
[`q9_interior.json`](q9_interior.json). The file includes the kind of the
original completed orbit and the digest of every complete 130-term
lowering row. It does not replace the existing public dispatch tables.

## The previously missing relation

Retain the numbering of the fifteen level-three completed classes in
[`../k5/REPORT.md`](../k5/REPORT.md). The exact identity is

\[
\boxed{
J_2=4J_3^{(14)}+16J_3^{(6)}+\frac{16}{3}J_3^{(5)}.
}
\]

Its derivation does not fit already-known zeta coefficients. The explicit
vertex and dimension rows first give

\[
J_2=24J_3^{(9)}+20J_3^{(14)}.
\]

The class `9` is the radial four-spoke wheel in six dimensions. The
uniform radial identity and the retained rational boundary then give

\[
J_2+6J_3^{(9)}=20J_3^{(6)}+\frac{20}{3}J_3^{(5)}.
\]

Eliminating `J_3^(9)` gives the boxed relation. The graph maps, all
normalization factors, rational telescoping function, and complete row
witness appear in [`DERIVATION.md`](DERIVATION.md) and
[`certificate.json`](certificate.json).

## What the guard checks

For the completed K6 integer weights, require vertex sums `2q`. The
complete determinant expansion is termwise convergent at `q-1` exactly
when

\[
\beta_e\le q-2\quad\text{for every edge},\qquad
\sum_{e\subset S}\beta_e\le2q-3
\quad\text{for every three-vertex set }S.
\]

Nonpositive position-space weights are allowed in the completed relation;
they represent polynomial numerators. The positive-weight K5 scalar
interface additionally checks every shifted original edge weight and
every lower scalar cut before consulting any table.

The obstruction is concrete. At `q = 5`, the convergent vector

```text
alpha = (2,2,3,4,2,3,3,4,3,4)
```

has a nonzero formal lowering term at `q = 4` with

```text
alpha' = (1,1,3,3,1,3,3,3,3,3).
```

All these lower edge weights are positive, but the triangle `012` has
weight sum `3`, below the required strict bound `4`. Its determinant
coefficient is `1/8`, and its falling-factor coefficient is `1`.
The software refuses this termwise reduction. A second example at
`q = 3` shows that discarding zero falling factors loses an entire
nonzero boundary contribution.

## Reproduction

From the repository root:

```bash
python3 explore/circuit-general/k5-uniform/verify.py
```

The replay independently constructs the determinant in both its six by
six distance form and its four by four Gram form, checks five values
against the previously independent upward relation tables, rebuilds all
73 new answers, scans all 54,391 existing K5 orbits for the proposed
lowering obstruction, and verifies the complete extra-input witness.
The existing tables are read as declared inputs; no production graph
reducer or period evaluator is imported. Rational arithmetic uses the
standard library; NumPy batches the finite boundary scan.

`--write` explicitly refreshes both complete certificates. The recorded
digests bind the implementation, proof, and every input table. New input
normalizations or stale source hashes cause the default replay to fail.

## Files

| File | Content |
|---|---|
| [`DERIVATION.md`](DERIVATION.md) | Full proofs, boundary hypotheses, extra-input origin, and failure examples |
| [`k5_dimension.py`](k5_dimension.py) | Complete determinant relation with explicit refusal of missing boundaries |
| [`verify.py`](verify.py) | Independent determinant, graph, relation, table, and obstruction replay |
| [`q9_interior.json`](q9_interior.json) | All 73 exact values in the stated new region |
| [`certificate.json`](certificate.json) | Complete row witness, five controls, boundary counts, and source digests |
| [`RUN.log`](RUN.log) | Successful default replay |

There is no claimed score change, arithmetic independence theorem,
complete `q=9` table, or unconditional iteration to the old base levels.
