# Complete products with a higher superpotential

This continuation treats nonrenormalizable superpotentials. It computes
complete Wick sums, including the representation coefficients, fermionic
signs, differentiated input legs, and all surviving source spinors.
Individual geometric forest periods and coefficients of the assembled
operator product are distinct quantities.

**Universal continuation:** [all color indices and arbitrary-degree pure
matter superpotentials](../universal-quartic-product/README.md) derives the
zeta cancellation from the representation identity before any color
contraction. It prints the complete quartic tensor and proves rationality
at every order for the same fixed input family and a neutral free spectator.

## What the Leibniz identity forces

For the full differential and marked product of
[the manuscript](../../../sources/LLC_main.tex), the quantum master
equation gives

\[
\mathbf Q(A\star B)
-(\mathbf QA)\star B-(-1)^{|A|}A\star(\mathbf QB)=0.
\]

This holds with a higher superpotential. Its proof uses the complete
marked higher operations, not the rationality of a particular graph.
The underlying geometric identities are those of
[Gaiotto, Kulp and Wu, arXiv:2403.13049](https://arxiv.org/abs/2403.13049).

Thus the entire period-valued contribution to the Leibniz defect must
vanish. This does not, by itself, set the period coefficients of
`A star B` to zero. At an order with rational lower-order forcing,
the difference from a rational particular solution obeys the
homogeneous classical bilinear chain equation. Closedness alone is
insufficient to remove that difference by an exact change of product.
The precise coefficient equation is in
[forest-leibniz/DERIVATION.md](../forest-leibniz/DERIVATION.md).

Conversely, a nonzero period in a single source-compatible scalar
forest does not prove survival in the complete product. The calculations
below explicitly test this additional cancellation.

## A complete abelian quartic product

Take a `U(1)` gauge theory with four matter fields whose charges are
\[
(1,-1,0,0),\qquad
\bar W=\gamma_1\gamma_2\gamma_3\gamma_4.
\]
Both the sum of charges and the sum of their cubes vanish. The
superpotential is gauge invariant and has degree four. Use the
manuscript's fields and spinor derivatives, with inputs
\[
A=\beta^1\partial_{+\dot a}\beta^2,\qquad B=b\beta^3.
\]
Both inputs have zero gauge charge. They are not assumed to be
`Q`-closed.

At order `g^2 hbar^3`, the complete coefficient is rational. With
`partial_1` acting on `gamma_4(x_1)` and `partial_2` on `c(x_2)`,
its physical normalization is

\[
\begin{aligned}
&[\hbar^3]\bigl(A\star B\bigr)\\
&\quad=
\frac{g^2}{552960(8\pi^2)^3}
\left[
\left(\partial_{1,+\dot b}\partial_{2,+}{}^{\dot b}\right)^3
\left(12\partial_{1,+\dot a}-19\partial_{2,+\dot a}\right)
\gamma_4(x_1)c(x_2)
\right]_{x_1=x_2=x}.
\end{aligned}
\]

Every term differentiates the surviving ghost. There is no
undifferentiated `c` output in this expression.

The stronger calculation in [quartic_wick.py](quartic_wick.py)
retains a free neutral spectator `gamma_5` in the first input,
`A=beta^1 (partial beta^2) gamma_5`, and computes the complete
source polynomial on `gamma_4,c,gamma_5`. Setting `gamma_5=1`
gives the displayed result. The four interacting charges are kept
symbolic until the end: all seven occurring quadratic charge monomials
have identically zero complete `zeta(3)` coefficient.

[FULL_QUARTIC_RESULT.md](FULL_QUARTIC_RESULT.md) prints every term
of the covariant source polynomial for the gauge-invariant charge
choice. [quartic_wick_certificate.json](quartic_wick_certificate.json)
contains every graph value, source assignment, matching and sign.

### Complete Wick enumeration

Fix the superpotential at interaction slot zero and retain the two
ordered matter-interaction slots. The physical edges are all allowed
pairings of the five `beta` fields with their corresponding `gamma`
fields, together with the contraction of the external `b` with one
of the two interaction ghosts. The differentiated `beta^2` is tracked
through its actual edge.

The enumeration gives twenty nonvanishing admissible matchings.
Twenty-eight further matchings have explicit parallel-edge or Laman
subset witnesses for vanishing. All three placements of the
superpotential agree after vertex, edge and current transport.
The input matter vertices are not replaced by a preselected forest.

The Wick signs are obtained by literal left differentiation using
manuscript equation `eq-A-20`. A separate test applies those
differential operators to every candidate physical edge set.

### Why finitely many source evaluations prove the full polynomial

At three loops, the current-area determinant has source degree six.
The inserted current raises the total degree to seven.
Varying only the spectator source changes the current determinant
through at most three independent cycle forms, so its degree in
that source is at most three. The current insertion raises this
bound to four.

For independent `lambda_gamma4` and `lambda_c`, an invertible
spinor change of basis fixes them to the two coordinate basis
vectors. The entire remaining dependence is a polynomial of total
degree at most four in the two spectator-source entries. The full
triangular interpolation set has
\[
\binom{4+2}{2}=15
\]
points. Exact Newton differences recover every coefficient.
Restoration uses the full spinor transformation law
`det(M)^3 M` and produces a polynomial in the original spinors
and their brackets. Polynomial continuation includes dependent
source spinors. Independent rational and dependent-source checks
reproduce the restored answer.

This is a complete polynomial calculation. It does not infer
an all-source statement from a collection of numerical tests.

### Normalization

For this three-loop coefficient:
\[
\begin{array}{c|c}
\text{factor}&\text{value}\\ \hline
\text{superpotential placements divided by }3!&1/2\\
\text{manuscript kernel divided by htk kernel}&1\\
\text{marked-product graph sign}&-1\\
\text{two matter-interaction factors of }i&(-i)^2=-1\\
\text{physical derivative on the oriented input head}&-2
\end{array}
\]

Their product is `-1`. The raw current polynomial therefore carries
`-g^2/(8*pi^2)^3`. Each retained source becomes
`-partial_+/2` on its own surviving field. Its total source degree
is seven, so the physical differential polynomial carries
`g^2/[128(8*pi^2)^3]`. Its compact raw coefficient `1/4320`
then gives `1/552960` in the displayed operator.

## A complete nonabelian quartic product

Take `SU(2)` with four adjoint matter fields, the invariant adjoint
metric used to contract repeated color indices, and

\[
\bar W=(\gamma_1^D\gamma_{3,D})(\gamma_2^E\gamma_{4,E}),\qquad
A=\beta^{1D}\partial_{+\dot a}\beta^2_D,\qquad
B=b_E\beta^{3E}.
\]

These are color-singlet cochains. No `Q`-closure is assumed.
The adjoint representation has vanishing cubic gauge-anomaly trace.

Retaining a neutral free `gamma_5` multiplying the first input,
the twenty matter–matter Wick matchings alone have a nonzero
`zeta(3)` polynomial. The eight matter–gauge matchings give its
negative. The full cancellation holds for all three retained
source spinors. In the compact specialization `gamma_5=1`, each
of these two period polynomials already vanishes separately.
All nine entries of the output color matrix are computed and
give the required invariant metric.

The complete three-loop result is

\[
\begin{aligned}
&[\hbar^3]\bigl(A\star B\bigr)\\
&\quad=
-\frac{g^2}{92160(8\pi^2)^3}
\left[
\left(\partial_{1,+\dot b}\partial_{2,+}{}^{\dot b}\right)^3
\left(8\partial_{1,+\dot a}+5\partial_{2,+\dot a}\right)
\gamma_4^D(x_1)c_D(x_2)
\right]_{x_1=x_2=x}.
\end{aligned}
\]

As in the abelian calculation, the stronger
[nonabelian_wick.py](nonabelian_wick.py) computation retains a
neutral free `gamma_5` multiplying the first input and proves the
full source polynomial before setting that spectator to one.
Every matching, every complete color contraction, every graph
period and the reconstructed covariant polynomial are in
[nonabelian_wick_certificate.json](nonabelian_wick_certificate.json).
[FULL_NONABELIAN_RESULT.md](FULL_NONABELIAN_RESULT.md) prints the
entire covariant source polynomial and the two opposite period
contributions, keeping the free spinor index.

The three-loop raw current polynomial after removing the spectator is
\[
\frac{\langle\lambda_{\gamma_4},\lambda_c\rangle^3}{720}
\left(8\lambda_{\gamma_4,+\dot a}+5\lambda_{c,+\dot a}\right).
\]
The implementation absorbs the adjoint matrices' factors of `i`
into real epsilon tensors in the matter interaction. With the
gauge-vertex coefficient `-g/2` retained, the common multiplier
of this raw expression is `+g^2/(8*pi^2)^3`.
The degree-seven source substitution supplies `-1/128`,
giving the displayed physical sign and denominator `92160`.

## Completeness at three loops and finite termination

For both displayed pairs of inputs there are four antifields:
three `beta` letters and one `b`. A product graph at loop order
`n` with `r` superpotential vertices loses `n+r` antifields.
Consequently `n+r<=4` is necessary.

If there is no superpotential vertex, the `n` cubic interaction
vertices and the input `b` supply exactly `n+1` contractible
Grassmann-even fields. Each of the `2n` physical contractions
requires one even field, so `2n<=n+1` is also necessary.
Thus cubic-only graphs vanish for every `n>=2`.

At `n=3` the only surviving possibility is exactly one
superpotential vertex and two cubic vertices. The displayed
calculations therefore give the **entire three-loop coefficient**,
not a selected graph contribution. At `n>=4` the antifield
bound excludes every graph with a superpotential vertex, while
the even-field count excludes the cubic-only graphs. Hence

\[
[\hbar^n](A\star B)=0\qquad(n\ge4)
\]

for these specific inputs in the stated theories. The neutral
spectator does not change this argument because it is never
contracted.

## Reproduce

From the repository root:

```bash
python3 explore/circuit-general/nonrenormalizable-product/run_checks.py --examples
```

The runner checks both independent Wick/color suites. `--examples`
also replays the two complete source-polynomial certificates.
Each default certificate command regenerates and compares the
complete deterministic calculation. `--write` refreshes its
certificate explicitly.

The complementary [sourced-ear reduction](../forest-general-reduction/README.md)
retains nonzero sources at removed cubic vertices and evaluates
individual higher-superpotential forests through eleven loops.
Its four- and five-loop identities exhibit exact period cancellation
in complete Leibniz combinations independently of this Wick assembly.
For a larger irreducible core, the [whole-numerator period solver](../forest-target-reduction/README.md)
reduces the earlier five-loop control to a complete seventeen-term residual
identity. It preserves the unresolved terms explicitly; a finite closure
of that control is still open.

## Scope

This proves cancellation in the specified complete abelian and
nonabelian quartic products, including every retained source.
The inputs are invariant under constant gauge transformations;
ordinary derivative jets are not asserted invariant under every
local gauge transformation. This is the constant-gauge relative
cochain calculation, with all positive ghost jets retained.
It does not establish
rationality of every product in every higher-superpotential theory.
No claim about nontriviality of a quantum cohomology class is made.
The inputs need not be cocycles for their cochain product to be
calculated.
