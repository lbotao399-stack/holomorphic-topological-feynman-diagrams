# Derivation of the complete fourth-order coefficient

The conventions are those of [the source manuscript](../../../sources/LLC_main.tex),
particularly `eq-A-20`, `eq-W-jet-edge`, `eq-D-lambda`, and `eq-2-31`,
and [the kernel convention](../../../code/htk/CONVENTIONS.md).

## 1. Exact field count

The interaction is

\[
-igc\sum_{i=1}^{4}q_i\gamma_i\beta^i
+\gamma_1\gamma_2\gamma_3\gamma_4.
\]

Here \(q_i\) denotes the charge of the indicated interacting field.
The original inputs are

\[
\beta^1(\partial_{+\dot a}\beta^2)\gamma_5,
\qquad b^2\beta^3,
\]

where \(\gamma_5\) is free, neutral, and absent from the interaction.
There are five original input antifields: three matter antifields and
two antighosts. If the interaction order is \(n\), and \(r\) interaction
vertices are pure superpotential vertices, the total number of
antifields is \(n-r+5\). The product has \(2n\) physical contractions,
so

\[
2n\le n-r+5,\qquad n+r\le5.
\]

With no superpotential vertex, the contractible even fields are the
\(n\) interacting matter fields and the two original antighosts; the
free spectator cannot contract. Thus

\[
2n\le n+2,\qquad n\le2.
\]

At fourth order there is therefore exactly one superpotential vertex
and three matter vertices. All eight antifields are consumed. The
remaining fields are exactly \(c\gamma_4\gamma_5\). Fifth and higher
orders vanish: a positive superpotential count violates \(n+r\le5\),
and a zero superpotential count violates \(n\le2\).

## 2. Complete original Wick sum

For a fixed superpotential position, the original vertex order is the
four interaction positions, the first input, and the second input.
Within a matter interaction the order is \(c\gamma_i\beta^i\).
Within the second input the two original \(b\)'s are temporarily
distinguished to record the ordinary derivative multiplicity of
\(b^2\). They are not divided by a factorial.

Every antifield is matched to its compatible original field at a
different vertex. Parallel contractions vanish; all vertex subsets
are checked against the four-dimensional Laman inequality. This leaves
324 complete matchings and 732 exclusions with explicit witnesses.
Each of the four possible superpotential positions gives the same
complete sum after edge orientation and edge-order signs are retained.

There are 14 distinct sourced geometries after relabelling the three
matter vertices. The differentiated \(\beta^2\) remains attached to
its original edge and contributes \(-2y_{+\dot a}\), with the
orientation sign separately retained under graph relabelling.

The independent calculation in [check_wick.py](check_wick.py) enumerates
the admissible graphs first and then acts with the manuscript's ordered
\(\Delta\) differential operators. Its complete polynomial agrees
with the matching enumeration. Every nonzero coefficient has the
expected factor two from differentiating \(b^2\).

## 3. All-source reconstruction is finite

The graph has four loops. Its bare holomorphic current determinant has
source degree eight; the original differentiated input adds one current,
so the full spinor has degree nine. The source of each output field
stays at the vertex where that field survives. Only the second marked
input is used as the balancing root.

The current form has four pairs of spinor-valued chord rows. The part
due to the spectator source has rank at most four: in each chord pair,
the two spinor components multiply the same scalar current differential.
Multilinearity therefore bounds the spectator degree in the bare
determinant by four. Its degree in the complete first-current spinor is
at most five.

Choose independent \(\lambda_{\gamma_4}\) and \(\lambda_c\) as a
temporary interpolation frame. The complete triangular grid of 21
nodes in [quartic_b_squared.py](quartic_b_squared.py) determines every
polynomial of spectator degree at most five by exact rational Newton
interpolation. The six highest-degree vector coefficients must combine
into the original spectator source spinor; the implementation checks
this condition before restoring every coefficient as products of

\[
\lambda_{\gamma_4,+\dot a},\quad
\lambda_{c,+\dot a},\quad
\lambda_{\gamma_5,+\dot a},\qquad
\langle\lambda_{\gamma_4},\lambda_c\rangle,\quad
\langle\lambda_{\gamma_5},\lambda_c\rangle,\quad
\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle.
\]

There are exactly four brackets in every term. A simultaneous invertible
linear change of spinor frame acts with the fourth power of its
determinant and with its usual action on the remaining spinor index.
The restored formulas are polynomial, so their identities extend to
dependent and zero source assignments by polynomial continuation.

The exact rational and every period channel are interpolated
separately. Every period channel vanishes at every unisolvent node,
for each symbolic cubic charge monomial. This proves its vanishing as
a full source polynomial. Four further direct source evaluations check
the restored formulas, including dependent and zero source cases.

## 4. Exact integrations and independent geometry

The original certificate integrates every inserted original current
with [CoreCurrentPeriods](../../../code/htk/forest_periods.py), using the
full current determinant and read-only exact period reduction. No
missing period is replaced by zero; convergent-series fallback is
disabled for this finite calculation.

At a generic full-source control, seven individual geometries contain
\(\zeta(3)\). The 16 nonzero symbolic cubic charge coefficients are
nevertheless rational after the complete sum. All 20 possible cubic
charge monomials, including the four zero coefficients, are written in
[the complete tensor table](FULL_SOURCE_TENSOR.md).

The independent [geometry check](check_geometries.py) evaluates all 54
required original current components using the product-forest
dispatcher: ten graphs contract source-free vertices first, one has a
direct series-parallel reduction, and three go through sourced core
periods. These complete results agree with the direct numerator route.
The last three share the same core-period library and are not claimed
as independent period evaluations.

## 5. Physical normalization

The fourth-order graph sign is positive, and the conversion from the
kernel convention to the manuscript is also positive:

\[
(-1)^4=1,\qquad (-1)^{4(4+1)/2}=1.
\]

The four superpotential placements, the three matter-vertex phases, the
head derivative, and the degree-nine source substitution give

\[
\frac{4}{4!}\,(-i)^3\,(-2)\,
\left(-\frac12\right)^9
=\frac{i}{1536}.
\]

Consequently the table of full source spinors multiplies

\[
\frac{ig^3}{1536(8\pi^2)^4}
\]

in the coefficient of \(\hbar^4\), with the source of each surviving
field replaced by the derivative acting on that field.

For \((q_1,q_2,q_3,q_4)=(1,-1,0,0)\),

\[
\sum_{i=1}^{4}q_i=0,\qquad
\sum_{i=1}^{4}q_i^3=0.
\]

The two input monomials and the quartic superpotential are invariant
under constant gauge transformations. No assertion that either input
is \(Q\)-closed is used in this calculation.

The exact source restriction
\(\lambda_{\gamma_5}=\lambda_{\gamma_4}\) simplifies the entire
physical spinor to

\[
\frac{\langle\lambda_{\gamma_4},\lambda_c\rangle^4}{9600}
\bigl(10\lambda_{\gamma_4,+\dot a}-3\lambda_{c,+\dot a}\bigr).
\]

Setting \(\lambda_{\gamma_5}=0\) gives zero. Both statements are
checked as exact identities in the complete covariant table.
