# Complete first-current spinors with arbitrary supported sources

## 1. Fixed graph, full source support, and normalization

The original ordered edges are

\[
(01,13,02,24,05,06,35,36,45,46,56).
\]

The marked vertices are 5 and 6. The original sources at vertices 1, 2,
and 5 are three arbitrary rational spinors; conservation fixes

\[
\lambda_6=-\lambda_1-\lambda_2-\lambda_5,
\qquad \lambda_0=\lambda_3=\lambda_4=0.
\]

This is the physical quartic forest already evaluated at one source assignment
in [sourced-forests](../sourced-forests/README.md). The result here includes its
entire supported-source bare polynomial and the entire first-current spinor
on every one of its eleven original edges. The quartic vertex retains its
specified zero source. No conclusion is made for sources outside this support.

The normalization is the existing `H=2,T=0` current-area form with the original
ordered edges. A first-current moment is the integral of that full form times
the complete original spinor current, equivalently the negative first
displacement derivative. No free spinor index is removed. Coupling, color,
Wick and displacement factorial factors are not attached to these geometric
kernels.

## 2. A priori finite polynomial reconstruction

The original five-loop density is polynomial of source degree ten, divided
by the source-independent Kirchhoff polynomial to power seven. An original
current has a source-linear polynomial numerator and one Kirchhoff
denominator. Thus a first-current moment has source degree eleven. Under a
simultaneous invertible two-by-two source transformation, the bare density
acquires the fifth power of the determinant, and each current additionally
acquires that same source transformation on its free spinor index.

The source degree can also be seen directly in the complete current-area
form: each original Kirchhoff current is linear in all original sources,
each of the ten current differentials is therefore source-linear, and their
wedge has degree ten. Multiplying by one full current raises that degree
to eleven. None of the denominators depends on a source.

Integration of these source polynomials is coefficientwise legitimate.
At every rational balanced source assignment, the established forest
current-area form is absolutely integrable; multiplying it by an original
current preserves absolute convergence because that current is a bounded
weighted average of its finitely many tree flows. Choose a finite rational
unisolvent grid for polynomials of degree at most eleven in the six
independent source coordinates. Each coefficient integrand is a finite
rational linear combination of these absolutely convergent evaluations.
The triangle inequality therefore proves absolute convergence of every
source coefficient separately. Finite linearity then permits its
coefficientwise integration. This uses the complete original forms and
does not infer convergence from a vanishing integrated bare value.

These are polynomial identities before integration. Integration, exact chain
factorials and the finite period lookup are rational linear maps, so they
preserve the identities coefficientwise. No independence of the zeta values
is needed.

Set the independent root spinors to `(1,0)` and `(0,1)`, and write the marked
spinor as `(a,b)`. Each full current then has two ordinary frame polynomials,
denoted here by `f(a,b)` and `g(a,b)`. These are a storage arrangement of the
complete spinor. Their total degree is at most six, as follows directly from
polynomial covariance.

To prove that bound, temporarily take

\[
\lambda_1=(1,0),\qquad \lambda_2=(s,t),\qquad\lambda_5=(0,1).
\]

For nonzero `t`, covariance makes the two entries of the same current

\[
t^5\bigl(f(-s/t,1/t)+s\,g(-s/t,1/t)\bigr),
\qquad t^6g(-s/t,1/t).
\]

The original current numerator is polynomial in `s,t`, including at `t=0`.
If the highest homogeneous frame degree exceeded six, the second expression
would force its highest `g` part to vanish; the first would then force the
corresponding `f` part to vanish. This contradiction proves the degree bound.
At degree six the first expression instead requires

\[
b f_6(a,b)=a g_6(a,b).
\]

Since `a,b` are relatively prime polynomials, there is a homogeneous
degree-five polynomial with `f_6` equal to `a` times that polynomial and
`g_6` equal to `b` times it. The same argument for the scalar density gives
frame degree at most five for the bare kernel.

There are exactly 28 bivariate monomials of total degree at most six.
Their values at the 28 integer nodes `a>=0, b>=0, a+b<=6` determine every
coefficient. In the falling-factorial basis the evaluation matrix is
triangular with nonzero diagonal `i!j!`; hence this is exact interpolation,
not a numerical fit. The production reconstruction inverts the rational
monomial matrix. The verifier independently uses bivariate finite
differences and expands the falling-factorial basis.

## 3. Fewer direct evaluations by an exact source exchange

Exchanging vertices 1 and 2 and simultaneously vertices 3 and 4 induces
the edge permutation

\[
(0\ 2)(1\ 3)(6\ 8)(7\ 9).
\]

Every orientation is preserved and the permutation is even. Following this
graph isomorphism by exchange of the two source coordinates restores the
root frame and takes `(a,b)` to `(b,a)`. The latter matrix has determinant
minus one, so the bare kernel changes sign, and the complete current on
the corresponding original edge changes sign and exchanges its two entries.
Thus the 16 nodes with `a<=b` supply all 28 values exactly. Every directly
evaluated node includes the bare kernel and both entries of each of the five
original chord currents. No zero bare value is used to discard a current.

## 4. Exact finite integration and every original edge

The sourced chains are `0-1-3` and `0-2-4`. Suppressing them leaves five
vertices and nine edges: the complete graph on five vertices with the edge
between vertices 3 and 4 absent. Exact simplex integration along each chain
contributes only factorial ratios, retaining the original sources and the
current inserted on its original edge.

Reciprocal Schwinger variables identify this scalar core with the planar-dual
triangular prism, replacing each scalar edge power by `q` minus that power.
The complete existing prism banks at levels seven and eight evaluate the
bare kernel and first-current insertions. Every period monomial is checked
for strict convergence. No table is generated and no absent value is replaced
by zero.

The complete original source minor is computed in a declared source-flow
tree and cycle basis. A grounded-network evaluation independently checks the
density, and complete tree-current polynomials are checked against every
original current and all fourteen coefficientwise source identities.

At each node the tree currents follow from the complete tree-supported
source flow plus the signed fundamental-cycle sum of the independently
integrated chord currents. The affine source term is multiplied by the bare
kernel and retained. All eleven resulting spinors satisfy the fourteen
integrated source equations. Their coefficientwise polynomial versions are
also checked after reconstruction.

## 5. Denominator-free result at every degenerate frame

For arbitrary supported source spinors, the frame coordinates are

\[
a=\frac{\langle\lambda_5,\lambda_2\rangle}
        {\langle\lambda_1,\lambda_2\rangle},
\qquad
b=\frac{\langle\lambda_1,\lambda_5\rangle}
        {\langle\lambda_1,\lambda_2\rangle}
\]

where the bracket is the determinant of the two displayed spinors. For a
frame monomial `a^i b^j` of degree at most five, its first or second entry
lifts by multiplication with the corresponding original spinor and

\[
\langle\lambda_1,\lambda_2\rangle^{5-i-j}
\langle\lambda_5,\lambda_2\rangle^i
\langle\lambda_1,\lambda_5\rangle^j.
\]

For the degree-six part, divide `f_6` by `a` (equivalently divide `g_6` by
`b`) as an exact polynomial operation. Lift the resulting degree-five
monomials with the original spinor `lambda_5` and the same bracket formula.
There are therefore no negative bracket powers in the final result.

The bare result uses the same bracket monomials without an extra spinor.
The resulting polynomial agrees with the integrated original graph on the
dense open set where the two root spinors are independent. Both expressions
are polynomials, so they agree for every source assignment, including
degenerate frames. The certificate retains all three original-source
directions and every free spinor index.

## 6. Independent controls and exact scope

The new-source control uses roots `(2,1)` and `(1,2)` and marked source
`(-1,4)`. Its determinant is three and its frame point is `(-2,3)`, outside
the interpolation grid. Every original edge is integrated separately with
a different source/cycle tree; interpolation and chord reconstruction are
not used by that control. A second complete original-graph control has
parallel root spinors, checking the denominator-free specialization.
The entire original numerator is also compared across the two source/cycle
bases. These controls share the established numerator formula and scalar
period tables; their independent parts are the fresh sources, source/cycle
construction and all-edge integrations.

The earlier bare physical quartic fixture supplies a further exact value
at a nonunit source determinant. A separate verifier reconstructs the
interpolation by rational finite differences, checks every coefficient of
every integrated source equation, and compares every stored full spinor.

This is a new complete supported-source current family on the physical
quartic core. It is not the different `K3,3 plus one edge` topology of
[forest-finite-jets](../forest-finite-jets/README.md). A direct probe on that
older topology produces 66,430 scalar powers for one first-current entry
and encounters a missing level-eight table. Its current completion remains
a separate relation problem.
