# Finite completion of the sourced five-loop product forest

The previously unresolved source-complete kernel is now evaluated exactly:

\[
\boxed{\frac{1187-816\zeta(3)}{435456000}}.
\]

There are **zero remaining scalar periods**. The convention is the existing
`htk` kernel, with the original ordered oriented graph and full balanced
sources below; physical coupling and propagator prefactors are not included.

```python
edges = ((0,1),(1,2),(0,3),(0,4),(2,5),(2,6),
         (3,5),(3,6),(4,5),(4,6),(5,6))
sources = ((0,0),(1,0),(0,0),(0,1),(0,0),(0,0),(-1,-1))
marked = (5,6)
```

Both certificate generation and the fast default verification reconstruct
all **10,577 original numerator coefficients** and all **7,368 scalar
coefficients**. Their hashes must equal the immutable references in the
[preceding complete 17-term identity](../forest-target-reduction/dimension_certificate.json).
No source is removed or changed during this calculation.

## Exact closure

The preceding identity contained the finite contribution

\[
-\frac{3706491467747}{27935373312000}
-\frac{1906264441}{5388768000}\zeta(3)
+\frac{57527}{106920}\zeta(5)
\]

plus seventeen explicit periods. The new calculation evaluates the **entire
seventeen-term remainder** as

\[
\frac{529509659453}{3990767616000}
+\frac{272322049}{769824000}\zeta(3)
-\frac{57527}{106920}\zeta(5).
\]

Adding these expressions gives

\[
\begin{aligned}
-\frac{3706491467747}{27935373312000}
+\frac{529509659453}{3990767616000}
&=\frac{1187}{435456000},\\
-\frac{1906264441}{5388768000}
+\frac{272322049}{769824000}
&=-\frac{17}{9072000},\\
\frac{57527}{106920}-\frac{57527}{106920}&=0.
\end{aligned}
\]

Positivity follows without any assumption about independence of zeta values.
The integral test after the first two terms gives

\[
\frac{85}{72}\leq\zeta(3)\leq\frac54,
\qquad
0<\frac{167}{435456000}
\leq\frac{1187-816\zeta(3)}{435456000}
\leq\frac{671}{1306368000}.
\]

## What changes the reduction

The earlier local equations preserved sixteen blocks of power sums. The
degree-two polynomial annihilators couple those blocks: the new equations
reach 4,096 scalar variables outside the previous closed system. Increasing
the old graph-search budget was unnecessary.

For the graph polynomial `U`, the archive supplies homogeneous vector fields
`a_i` and a polynomial `b` satisfying

\[
\sum_i a_i\,\partial_i U=bU.
\]

Every coefficient of this polynomial identity is verified exactly. The
scalar convention is `integral_Omega product(t_e**(powers_e-1))/U**q`, with
no product-Gamma normalization; `q` is the denominator power. Positive base
powers `beta` satisfy `sum(beta) = (loop number)*q - degree(a_i) + 1`. Then

\[
\sum_i\partial_i\!\left(
 a_i\,\frac{\prod_j t_j^{\beta_j-1}}{U^q}
\right)
=\frac{\sum_i\partial_i(a_i\prod_j t_j^{\beta_j-1})
-q b\prod_jt_j^{\beta_j-1}}{U^q}.
\]

The inherited [`syz5.relation`](../fiveloop/syz5.py) includes every
coordinate-face contribution. In particular, terms with `beta_i=1` retain
the restriction of `a_i` to `t_i=0` and the corresponding contracted graph.
An identity is rejected if its complete derivative terms are not absolutely
integrable or a required boundary period cannot be classified. No shifted
term or boundary contribution is dropped.

The method is the parametric-annihilator approach of
[Bitoun, Bogner, Klausen and Panzer, arXiv:1712.09215](https://arxiv.org/abs/1712.09215).
The concrete three degree-one and fourteen degree-two generators already
existed in the repository's `fiveloop/syz_cache`. The new rational JSON
archive records that file's hash and its producer. Independent controls
verify exact equality with the historical generators after transport to
the actual ordered core. This contribution adds source-target generation,
bounded exact elimination, and a directly summable proof of the full kernel.

## Small proof, reproducible search

| Quantity | Exact count |
|---|---:|
| Polynomial generators | 17 |
| Verified Kirchhoff monomials | 135 |
| Candidate base/generator pairs | 44,382 |
| Complete nontrivial rows | 42,472 |
| Convergence-guard rejections | 1,898 |
| Identically zero relations | 12 |
| Scalar variables reached | 14,553 |
| Elimination pivots before the entire target is finite | 2,580 |
| Nonconstant rational coefficient products | 1,686,085 |
| Rows in the final exact witness | 1,039 |
| Residual periods in the kernel | **0** |

The solver stops as soon as the complete target is finite. It need not solve
every other scalar variable in the generated system.

[`relation_witness.py`](../symmetry/relation_witness.py) records row
dependencies and propagates proof weights backward through the elimination.
The resulting 1,039 rational weights are attached to original complete
identities. Direct exact summation checks that they reproduce every one of
the seventeen target coefficients and the displayed finite remainder.
The [65,399-byte certificate](finite_certificate.json) stores all witness
bases and weights, not a digest in place of the proof.

[`annihilator_target.py`](../symmetry/annihilator_target.py) accepts another
complete scalar target and a verified polynomial archive. Candidate,
variable, row, arithmetic, fill, and proof-history bounds are explicit.
A stopped computation retains its complete finite-plus-residual expression;
`.value` is available only after all residual coefficients vanish. Archives
are rational JSON, and no annihilator basis or period table is built
automatically. Input cores must be simple; the existing classifier handles
transformed multigraphs within complete identities.

## Run

Dependencies are listed in [requirements.txt](requirements.txt), shared
with the existing exact forest integrator.

```bash
# Independent controls and fast exact proof verification:
python3 explore/circuit-general/forest-finite-target/run_checks.py

# Also repeat the bounded search and compare every recomputed field:
python3 explore/circuit-general/forest-finite-target/run_checks.py --recompute

# Explicitly regenerate the certificate:
python3 explore/circuit-general/forest-finite-target/certify_finite.py --recompute --write
```

Default proof verification reconstructs the full original source numerator,
checks implementation/data hashes, verifies every polynomial annihilator,
rebuilds the witnessed identities including their faces, sums the complete
target exactly, and checks the positive rational interval. It does not
repeat the larger search. `--recompute` additionally reproduces the search
statistics and witness selection.

Independent controls cover **917 exact lower-level period identities**,
**310 convergence/face guard agreements**, **34 fresh rational checks of
the `q=7` differential identity**, and contracted
coordinate-face polynomials. Five additional regressions cover whole-target
preservation under bounds, proof tampering, contradictory equations,
malformed generators, and unavailable graph archives. Timings and observed
commands are in [CHECKS.log](CHECKS.log).

This finishes the listed source-complete five-loop kernel. It does not
claim a complete `K33+e, q=7` scalar table or a universal finite zeta basis
for every product forest. Earlier hashed implementations and certificates
remain unchanged.
