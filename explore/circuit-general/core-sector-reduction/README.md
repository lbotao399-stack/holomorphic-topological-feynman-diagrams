# Exact relative-period reductions with complete boundary witnesses

The fixed-core [stabilization theorem](../core-period-stabilization/DERIVATION.md)
has a computational relation system: on every Hepp cube, use all derivatives

\[
\partial_{x_i}\left[x_i(1-x_i)\frac{x^a}{U_\sigma^q}\right].
\]

The primitive vanishes at both faces in its differentiation coordinate.
The complete denominator has a positive constant coefficient and nonnegative
coefficients, so the primitive has no pole on the closed cube. Its integral is
exactly zero. No face contribution is discarded or estimated.

The new [production module](../../../code/htk/core_sector_periods.py) implements
this full relation schema, graph-to-sector pullback, and a bounded exact search
for a reduction to supplied master forms. Every successful reduction contains
the entire list of rational master coefficients and boundary primitives.
Clearing the common denominator verifies the resulting rational-function
identity exactly. A second verifier here differentiates the primitives directly
with SymPy, independently of the production relation generator and elimination.

This is a **certificate-producing finite search**, not a general master-basis
algorithm. A finite search failure retains the whole residual and returns no
period value. An observed rank plateau is not used to certify completeness.

## Complete relation

Here `a` is a tuple of nonnegative integer numerator powers, `q` is a
nonnegative integer, and `i` labels a cube coordinate. The expanded identity is

\[
\begin{aligned}
0=\int_{[0,1]^m}\Bigg[&
\frac{\bigl((a_i+1)-(a_i+2)x_i\bigr)x^a}{U_\sigma^q}
\\[-1mm]
&-q\frac{x_i(1-x_i)x^a\,\partial_{x_i}U_\sigma}
{U_\sigma^{q+1}}\Bigg]\,dx_1\cdots dx_m.
\end{aligned}
\]

The complete infinite family presents the top relative de Rham quotient in the
theorem. The implementation enumerates only the explicitly requested finite
primitive powers and numerator degrees. The certificate asserts the complete
target identity, never that these finite rows generate the universal quotient.

## Complete evaluated controls

### A graph period, including every sector

For the three-edge banana, with positive scalar edge powers `(2,2,2)` and
denominator power `3`,

\[
U=t_0t_1+t_0t_2+t_1t_2,
\qquad
I_3(2,2,2)=\frac12.
\]

All six ordered sectors have the same transformed integral

\[
\int_0^1\!\int_0^1
\frac{x}{(1+x+xy)^3}\,dy\,dx=\frac1{12}.
\]

The certificate reduces each complete sector to `1/12` times the constant
cube form, and sums all six sectors. Independently, integrating in `y` gives

\[
\frac12\int_0^1\left(\frac1{(1+x)^2}-\frac1{(1+2x)^2}\right)dx
=\frac12\left(\frac12-\frac13\right)=\frac1{12}.
\]

Projective inversion also turns the complete original period into the unit
triangle Dirichlet integral, whose value is `1/2`.

### Nonrational cube controls

These are controls of the relative algorithm. They are not asserted to be
new irreducible Feynman core periods or bases of the original graph family.

\[
\int_0^1\frac{x^7}{(1+x)^5}\,dx
=\frac{1553}{64}-35\log2.
\]

The supplied master forms are `1` and `1/(1+x)`. The certificate retains all
six nonzero boundary primitives, including the terms at the higher pole orders.

The two-variable control is

\[
\int_0^1\!\int_0^1\frac{x^2y^2}{(1+xy)^3}\,dy\,dx
=\frac{\pi^2}{12}-\frac32\log2+\frac14.
\]

Its three supplied master forms are
`1`, `1/(1+xy)`, and `x/(1+xy)`, with integrals
`1`, `pi^2/12`, and `2 log2-1`. Their exact reduction coefficients are

\[
\left(-\frac12,\ 1,\ -\frac34\right).
\]

Direct integration independently gives

\[
\begin{aligned}
\int_{[0,1]^2}\frac{dx\,dy}{(1+xy)^2}&=\log2,\\
\int_{[0,1]^2}\frac{dx\,dy}{(1+xy)^3}&=\frac14+\frac12\log2,
\end{aligned}
\]

and the identity `(xy)^2/(1+xy)^3 = 1/(1+xy)-2/(1+xy)^2+1/(1+xy)^3`
reconstructs the same complete value.

## API and replay

```python
from htk.core_sector_periods import CubePeriodRelations

cube = CubePeriodRelations({(0,): 1, (1,): 1})
target = cube.monomial((7,), power=5)
masters = [cube.monomial((0,)), cube.monomial((0,), power=1)]
witness = cube.search(target, masters,
                      max_primitive_power=4, max_primitive_degree=8)
assert cube.verify(witness)
assert tuple(map(str, witness.coefficients)) == ('1553/64', '-35')
```

The graph-facing `CoreSector` preserves every original edge label, including
parallel edges. Its default pullback checks every proper edge-subset convergence
inequality. The optional sector-only mode is explicitly named and does not claim
convergence of the other sectors. Scalar orientation is the positive projective
measure; the original oriented source coefficient remains an external exact
scalar. Polynomial, tree-enumeration, relation-count and sparse-witness budgets
fail explicitly before a partial answer can be returned as a complete value.

```bash
python3 explore/circuit-general/core-sector-reduction/run_checks.py
```

[reductions.json](reductions.json) contains every master form, primitive and
rational coefficient for the eight complete controls. The replay checks the
stored certificates without rerunning elimination, and independently binds
every denominator, target density and supplied master to the stated control.
It differentiates all eight complete identities and checks all 82 primitive
end faces. `certify.py --write` rebuilds
them using the declared finite bounds. The focused production tests additionally
compare transformed densities against original Kirchhoff determinants and full
Jacobian matrices, reject divergent weights, retain unsuccessful-search
residuals, and reject deliberately corrupted witnesses. Ten tests include a
valid algebraic witness for a changed target: the target-specific replay rejects
it even though its integrated value agrees.

This module is additive. It changes no existing period table, hashed numerator
engine, current convention, or physical Wick/color assembly.
