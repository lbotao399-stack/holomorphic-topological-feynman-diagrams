# Complete six-loop K5 kernels by reciprocal core reconstruction

**27 November 2026.** The three admissible $K_5$ subdivision shapes at six
loops and $H=2$ now have complete bare $T=0$ kernel evaluations at explicit,
conserved source assignments. Every original vertex source is nonzero and
retained. There are **no remaining projective integrals** in these results.
The evaluator accepts other exact integer or rational source assignments;
the certificates below are fixed-source evaluations, not a symbolic
polynomial in all source variables.

The [family certificate](k5_family.json) records every original oriented
edge, source vector, exact result, reconstruction bound and table hash.
The [derivation](RECIPROCAL.md) proves the reduction and stopping criterion.

## 1. The three completed results

Number the five core vertices $0,1,2,3,4$ and order the ten core edges
lexicographically. Subdivide the displayed edges once, assigning the new
vertices $5,6,7$ in that edge order. Orient each resulting path from the
smaller original endpoint to the larger endpoint, and list its edges in
path order. The conserved source vectors, in vertex order, are

```text
0: (-1,  3)    1: (-1, -2)    2: ( 1,  3)    3: ( 0, -2)
4: ( 3,  0)    5: ( 1, -1)    6: ( 0, -3)    7: (-3,  2)
```

Their sum is $(0,0)$. In the normalization
$K_\Gamma(\lambda;z=0)$ of
[CONVENTIONS.md](../../../code/htk/CONVENTIONS.md), with $\zeta$ the Riemann
zeta function, the results are:

| Subdivided core edges | Shape | Exact complete kernel |
|---|---|---|
| $01,23,34$ | Two-edge path and a disjoint edge | $\dfrac{823413972469}{182891520000}-\dfrac{10923017647}{635040000}\zeta(3)+\dfrac{37445089}{2352000}\zeta(5)$ |
| $01,02,12$ | Triangle | $\dfrac{362193894187}{18289152000}+\dfrac{1985715463}{42336000}\zeta(3)-\dfrac{12884579}{176400}\zeta(5)$ |
| $01,12,23$ | Three-edge path | $-\dfrac{482542848959}{15240960000}-\dfrac{210499252657}{3810240000}\zeta(3)+\dfrac{20153461}{201600}\zeta(5)$ |

Individual certificates:
[disjoint-edge case](six_loop_k5.json),
[triangle](six_loop_k5_triangle.json), and
[three-edge path](six_loop_k5_three_edge_path.json).

All three values are strictly positive. This is certified by exact
rational partial sums and integral-test bounds for the zeta tails; no
linear-independence assumption about $1,\zeta(3),\zeta(5)$ is used:

| Shape | Exact enclosing interval |
|---|---|
| Two-edge path and disjoint edge | $[83149129/250000000,\ 13457827/40000000]$ |
| Triangle | $[440635169/1000000000,\ 451155691/1000000000]$ |
| Three-edge path | $[1395803389/250000000,\ 5595616629/1000000000]$ |

### Completeness of the topology list

There are $\binom{12}{3}=220$ ways to distribute three inserted vertices
among ten core edges. The original five core vertices must induce at most
$2\cdot5-3=7$ edges, so three **distinct** edges must be subdivided. This
excludes all 100 allocations with repeated subdivisions.

Every four original vertices must induce at most $2\cdot4-3=5$ edges.
Consequently the subdivided edges must meet every four-vertex subset.
The three-edge star fails this condition; its 20 labeled allocations are
excluded. The remaining $10+60+30=100$ allocations form the three shapes
above. The script checks every vertex subset of all 220 graphs, agrees
with `catalog.tight_ok`, and checks complete-graph isomorphism rather than
only the degree sequences. See [run_k5_family.py](run_k5_family.py).

## 2. What makes the computation feasible

The original numerator has degree 35 in thirteen edge resistances. Chain
integration leaves a polynomial of degree 35 in ten core resistances,
with degree at most five in each one. Inverting each resistance and
multiplying by the fifth power of every conductance changes its total
degree to $10\cdot5-35=15$.

The graphic-rank inequalities of [RECIPROCAL.md](RECIPROCAL.md) then give:

| Computation | Exact size |
|---|---:|
| Unconstrained original-edge interpolation grid | 1,788,632,625 samples per prime |
| Reciprocal grid before graphic-flat constraints | 877,844 samples per prime |
| Reciprocal grid with graphic-flat constraints | 683,834 samples per prime |
| Allowed full homogeneous support | 518,124 exponent vectors |
| Original current determinant | $12\times12$ |
| Contracted current determinant | $9\times9$ |
| Rational sparse cubature used in the previous fibre certificate | 165 nodes per core point |
| Exact finite-field Gauss integration | $3^3=27$ nodes per core point |
| Contracted determinants per reconstruction prime | $683834\cdot27=18463518$ |

The compiler evaluates one core solve per conductance point, reuses the
affine determinant coefficients at all 27 nodes, and batches the modular
pivot inversions. Every operation contributing to the answer is integer
or rational arithmetic. The square roots in the three-node rule are exact
elements of the selected prime fields.

The resulting integer-scaled coefficient counts are:

| Shape | Nonzero coefficients | Proved coefficient bound, in bits | Reconstruction modulus, in bits |
|---|---:|---:|---:|
| Two-edge path and disjoint edge | 458,342 | 195 | 217 |
| Triangle | 446,594 | 196 | 217 |
| Three-edge path | 459,743 | 197 | 217 |

## 3. Deterministic reconstruction certificate

Multiplication by $60^3=216000$ clears every chain-integration denominator.
The coefficients of this scaled polynomial are integers. A permanent of
the exact tree-current coefficient norms, combined with Mahler's
coefficient inequality, bounds their absolute values. The bound and its
inputs are stored in each certificate; the proof is in
[RECIPROCAL.md, section 5](RECIPROCAL.md).

For the first fixture the bound is exactly

$$25384925169104535878250163645356935162028057805002375168000.$$

The seven distinct primes are

```text
2147483647  2147483629  2147483587  2147483579
2147483563  2147483477  2147483399
```

Their product exceeds twice each fixture's bound. Centered Chinese
remainder reconstruction therefore determines **every coefficient
uniquely**. Agreement of successive reconstructions is not used to stop.
All excluded homogeneous-support coefficients are checked to be zero.

For every interpolation point, the grounded core Laplacian determinant is
an integer between $1$ and $125\cdot6^4=162000$. Thus none of the seven
large primes creates a singular core denominator anywhere on this grid.
An additional nonunit-gauge sample is checked at every prime.

Each reconstructed polynomial is also evaluated exactly at two positive
conductance points and compared with independent rational 165-node fibre
cubature. For the first fixture the two complete fibre densities are

$$\frac{271759}{195312500000},\qquad
-\frac{1194839063}{4246732800000},$$

at conductances $(1,1,1,1,1,1,1,1,1,1)$ and
$(1,2,3,1,2,3,1,2,3,1)$, respectively. These replays test the intermediate
polynomial. The period sum in section 1 completes the remaining integration.

## 4. Reusing the existing period calculations

The new [K5 bridge](../period-dispatch/README.md) makes all **54,391**
tracked convergent period orbits at levels 3 through 8 available through
`per5.period` and a dedicated bulk lookup. Level 8 contains 41,428 orbits.
It checks homogeneity and strict convergence and retains exact rational
coefficients of $1,\zeta(3),\zeta(5)$.

Weighted vertex signatures identify an exact graph orbit. Integer
numerator coefficients are added within that orbit before multiplication
by its rational period. The first fixture uses 26,119 nonzero accumulated
orbits. Missing data and divergence raise distinct errors in the bridge.
Its calls never generate missing period tables.

The existing level-8 table is unchanged. Its SHA-256 is

```text
7c3028f06348032b4d814f5dc05c9d377ed98254cf0d6bc6aa3a022055620bec
```

This contribution supplies the previously missing practical numerator
calculation and connection to those period certificates.

The [verification log](reciprocal_verify.log) records 58 focused test
groups and six existing engine controls. They include complete coefficient
comparisons with the original exact whole-graph numerator, all projective
gauges, edge order and orientation changes, rational sources, the 21 Gauss
moment identities, singular-input rejection, corrupted-checkpoint
rejection, and exhaustive small-level K5 period convergence checks. The
[family log](k5_family.log) records the three complete computations and
their independent rational replays.

## 5. Run and extend

Dependencies: Python with NumPy, SymPy and NetworkX, plus `g++` with OpenMP
support. The sampler builds its shared library in the system temporary
directory. Modular arrays are optional local checkpoints, protected by
implementation fingerprints and checksums; they are not committed.

```sh
# Original checks plus all reciprocal, sampler, reconstruction and bridge tests.
python3 explore/circuit-general/series-fibres/run_checks.py --reciprocal

# Recompute the first complete six-loop fixture.
python3 explore/circuit-general/series-fibres/run_core.py --threads 4

# Recompute all three sourced fixtures and the exhaustive topology census.
python3 explore/circuit-general/series-fibres/run_k5_family.py --threads 4 --recompute
```

With `code` and this directory on `PYTHONPATH`, arbitrary exact source
assignments use the same public computation:

```python
from htk.series import SeriesReduction
from core_numerator import certified_core_numerator, integrate_core

reduction = SeriesReduction(graph, full_conserved_sources)
polynomial = certified_core_numerator(reduction, threads=4)
kernel, period_certificate = integrate_core(polynomial)
```

The compiled sampler currently handles bare $T=0$ numerators through six
loops, with chains of length one or two. These include every admissible
six-loop $H=2$ graph with $K_5$ core. Exact source scaling, edge orientations,
edge ordering and alternate projective gauges are preserved. The general
reciprocal bound and existing rational-simplex fibre integration have wider
scope than this compiled backend. Explicit point, subset, cubature and
height-certificate budgets stop unsupported allocations.

Current insertions in the compiled reconstruction, higher-order interval
rules, longer-chain simplex rules, the $T=1$ contact terms, and period
tables for additional core families remain extensions. These three
completed bare kernels do not constitute all six-loop amplitudes or the
full displacement-dependent generating function.
