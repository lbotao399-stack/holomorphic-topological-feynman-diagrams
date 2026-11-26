# Exact chain integration and graph-rank numerator bounds

**26 November 2026.** Computational contribution for the full $T=0$ current
determinant, with arbitrary conserved sources and finite original-current
monomials. The [derivation](DERIVATION.md) proves the identities. The
[certificates](examples.json) specify every edge, source, resistance and value.

## What changes

1. **Contract constant source columns once.** On a chain, successive current
   columns differ by $C_c\otimes\lambda_v$. With $N$ original edges and $m$
   maximal chains, exact elimination reduces the determinant from $N-1$ to
   $m-1$. Every source and original edge current is retained. The core inverse
   is reused at all chain-fraction points.
2. **Integrate the chain fractions first.** At fixed core resistances, the
   determinant with $k$ current factors is polynomial of total degree at most
   $m-1+k$. Rational Grundmann–Möller rules combined by sparse tensor
   differences integrate it exactly, including the full series Jacobian.
3. **Prune the numerator before sampling.** For the existing numerator
   $P=U^{n+H}D$, every nonempty edge subset satisfies

   $$\deg_S P\le n\operatorname{rank}C_S-1,$$

   where $C_S$ is its cycle-column matrix. In particular, the total degree
   on an entire series chain is at most $n-1$. Graph-complement bridges
   enumerate the relevant rank constraints. A bounded enumeration may return
   fewer constraints without compromising correctness.
4. **Remove setup enumeration.** A single tree traversal constructs the
   particular current. Enumerating all spanning trees is no longer needed by
   `FastKernel`; zero-source initialization also works.

The existing [assembly interpolator](../assembly/numerator.py) and
[five-loop evaluator](../fiveloop/fk5.py) now use the rank constraints by
default. `BareNumerator(..., rank_bounds=False)` retains the old support as a
reference. Sparse grids have an explicit allocation budget, and singular
modular samples raise an error. Period-weight grouping now compares complete
integer tuples, avoiding overflow in the former packed keys.

## Exact computational evidence

The [consolidated replay](verify.log) passes all 33 focused test groups and
six existing engine controls. The longer full-period and higher-loop fibre
replays are recorded separately below.

The restricted and unrestricted interpolators give identical **complete
coefficient dictionaries** at primes $2147483629$ and $2147483587$:

| Fully sourced graph | Old samples per prime | New samples per prime | Nonzero coefficients |
|---|---:|---:|---:|
| Three-loop $K_4$ subdivision, $H=3$ | 8,272 | 4,141 | 1,138 |
| Four-loop prism, $H=2$ | 56,598 | 48,475 | 11,039 |

The updated five-loop fan replay uses **3,906,400** samples per prime instead
of **8,209,410** and agrees with the independent series–parallel integral:

$$K=-\frac{5713177}{576000}.$$

The full four-loop controls also agree, in the repository's conventions:

$$
K_{\mathrm{prism}}=-\frac{577}{40320}-\frac{2}{35}\zeta(3)
                    +\frac{5}{63}\zeta(5),\qquad
K_{K_{3,3}}=-\frac{789}{640}+\frac{59}{80}\zeta(3).
$$

Sources for these complete-period controls are in
[`check_fk5.py`](../fiveloop/check_fk5.py); the replay is
[`check_fk5.log`](check_fk5.log). The fan uses deterministic nonzero sources
at every vertex. These are regression controls, not new period evaluations.

### Higher-loop chain integrals

Both graphs below pass the graph-admissibility condition. At the retained
core resistances specified in [examples.json](examples.json), the complete
chain-fraction integrals are:

| Graph | Determinant size | Exact integration nodes | Remaining core dimension |
|---|---:|---:|---:|
| Five-loop cube subdivision, $H=3$ | $15\to11$ | 945 | 11 |
| Six-loop $K_5$ subdivision, $H=2$ | $12\to9$ | 165 | 9 |

$$
\text{Cube fibre density}=
-\frac{1190741613927780610249}{109165759454560087711309761228960},
$$

$$
\text{$K_5$ fibre density}=
\frac{1759441639913}{1678343272235712000000}.
$$

All **1,110** integration nodes are independently replayed with the original
grounded vertex Laplacian and the original full determinant. These values
are exact functions of the retained core resistances evaluated at the
specified points; the remaining projective core integrals are **unevaluated**.

The integration tests also compare against old whole-graph Fraction
numerators integrated monomial by monomial using factorial Dirichlet moments,
including a generic $H=3$ tetrahedron and current insertions. Cubature tests
check all 3,003 independent monomials through degree five on ten intervals
using 221 rational nodes.

Separate nonzero pointwise controls compare all currents and the complete
scalar density against the original vertex solver at eight loops ($H=3$,
determinant $24\to15$) and twelve loops ($H=2$, $24\to23$). Both graphs
pass admissibility; these are pointwise controls, not completed periods.

## Use and replay

With the repository's `code` directory on `PYTHONPATH`:

```python
from htk.series import SeriesReduction
from htk.fibres import integrate_fibre

# graph and sources describe the original graph, before any contraction.
reduction = SeriesReduction(graph, sources)
result = integrate_fibre(reduction, core_resistances,
                         powers={(original_edge_id, component): 2})
print(result.value)  # Fraction, including the chain-coordinate Jacobian
```

Run from the repository root:

```sh
python3 explore/circuit-general/series-fibres/run_checks.py
python3 explore/circuit-general/series-fibres/run_examples.py
python3 explore/circuit-general/fiveloop/check_fk5.py
```

The first command runs the five focused test suites and six existing engine
controls. The second regenerates the complete five/six-loop certificates.
The third is the longer complete-period regression. The combined form is
`run_checks.py --examples --periods`.

## Boundary of the result

The implementation covers $T=0$ and arbitrary finite current monomials.
The $T=1$ Gaussian contact terms remain a separate extension. The full
displacement exponential is handled coefficient by coefficient, not by a
finite-degree rule for the entire exponential.

Core period tables and subsequent core integration still limit complete
higher-loop answers. Fibre integration is especially useful for a small core
with many subdivisions; it is not uniformly faster for graphs with few
subdivisions. The support bound and the determinant identity hold at every
loop order; finite replays test their implementation rather than replace
their proofs.
