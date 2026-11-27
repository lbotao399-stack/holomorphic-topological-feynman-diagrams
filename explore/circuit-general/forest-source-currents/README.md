# Eleven complete five-loop current spinors at arbitrary supported sources

This calculation gives the entire bare source polynomial and **all eleven
original first-current spinors** for the physical quartic forest with ordered
edges

\[
(01,13,02,24,05,06,35,36,45,46,56).
\]

The marked vertices are 5 and 6. The three independent source spinors are
arbitrary rational values at vertices 1, 2 and 5:

\[
\lambda_6=-\lambda_1-\lambda_2-\lambda_5,
\qquad \lambda_0=\lambda_3=\lambda_4=0.
\]

Every coefficient is evaluated in the rational span of
\(1,\zeta(3),\zeta(5)\), with zero residual scalar periods. The result is a
polynomial in the three original spinors, with no inverse source brackets.
It includes every degenerate source frame. [FULL_RESULTS.md](FULL_RESULTS.md)
prints every coefficient, retaining the free spinor index on every original
edge. [certificate.json](certificate.json) records the complete exact data.

These are the normalized `H=2,T=0` geometric kernels, before Wick, color and
coupling factors. The original current is inserted with the sign of the
negative displacement derivative. Sources outside the displayed support are
not part of this family.

## A vanishing bare kernel retains nonzero currents

When \(\lambda_5=0\), the bare kernel vanishes for every pair of root sources.
Nevertheless the complete current on the original edge 01 is

\[
\begin{aligned}
\langle\lambda_1,\lambda_2\rangle^5\bigg[
&\left(\frac{923}{1828915200}-\frac{\zeta(3)}{8467200}\right)
\lambda_{1,+\dot a}\\
&+\left(\frac{17}{114307200}-\frac{\zeta(3)}{8467200}\right)
\lambda_{2,+\dot a}\bigg].
\end{aligned}
\]

Both scalar coefficients are strictly positive: the elementary exact bound
\(\zeta(3)\le5/4\) gives respective lower bounds
\(653/1828915200\) and \(2/1828915200\). Thus the bare zero cannot be used
to omit the first-current spinor.

## Complete reconstruction and checks

The source degree is ten for the bare kernel and eleven for a first current.
Polynomial covariance proves that, after choosing the two root spinors as
the coordinate frame, each current entry has degree at most six in the two
marked-source coordinates. The 28 triangular integer nodes therefore
determine every coefficient exactly. An original-graph automorphism supplies
12 of those nodes from 16 complete directly integrated nodes. Every direct
node retains both entries of all five chord currents and reconstructs all
tree currents with their affine source flow intact.

The degree-six frame terms obey an exact compatibility identity. It gives
a denominator-free expression in the three original source directions,
proving validity even when the root-source bracket vanishes. The
[derivation](DERIVATION.md) proves this before integration, so no arithmetic
independence of the zeta values is assumed.

The independent verifier uses rational finite differences rather than the
production matrix inverse. It checks all reconstructed frame coefficients,
every coefficient of all fourteen original source-conservation equations,
nonunit source transformations, and the earlier independent physical-quartic
bare value.

Two further controls integrate every original edge separately, without
source interpolation or chord reconstruction:

- [Fresh determinant-three frame](fresh_nonunit_frame.json): original root
  sources `(2,1)`, `(1,2)` and marked source `(-1,4)`, outside the sampling grid.
- [Degenerate root frame](degenerate_frame.json): root sources `(1,0)`, `(2,0)`
  and marked source `(0,1)`.

Both use a different source-flow tree and cycle basis. The complete original
numerator agrees coefficient by coefficient with a fresh construction in the
sampling basis. The controls share the established numerator formula and
scalar period tables; they are independent source, basis and full-edge
integration checks.

## Replay and executable result

From the repository root, with the dependencies in
`../all-product-forests/requirements.txt` installed:

```bash
python3 explore/circuit-general/forest-source-currents/run_checks.py
python3 explore/circuit-general/forest-source-currents/run_checks.py --deep --jobs 5
```

The first command verifies the complete stored data, exact dependency hashes,
independent controls, and printed result. Missing data is an error. The second
recomputes all direct interpolation nodes and both complete original-edge
controls. Exact node checkpoints may be used during construction; the published
certificate contains every required sample and does not depend on those files.

`source_currents.five_loop_source_currents(sources)` evaluates the stored
polynomial at any balanced rational assignment on the displayed support. It
returns the bare kernel and every full original current spinor as exact
rational coefficient dictionaries. Unsupported nonzero sources are rejected.

## Relation to the earlier five-loop source polynomial

This extends the [physical quartic forest](../sourced-forests/README.md) from
one evaluated source assignment to its complete supported-source bare and
first-current polynomials. Its two sourced chains leave the graph obtained
by deleting one edge from the complete graph on five vertices. Planar duality
then uses the complete existing prism tables at levels seven and eight.

The [other five-loop bare source family](../forest-finite-jets/README.md) has a
different scalar core, the complete bipartite graph on three plus three
vertices with an extra edge. Its level-eight current periods require further
relations. This calculation makes no completion claim for that other family,
for arbitrary sources on vertices 0, 3 or 4, or for a complete interacting
Wick/color product.
