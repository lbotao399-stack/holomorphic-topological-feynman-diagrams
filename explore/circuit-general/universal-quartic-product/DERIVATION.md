# A complete representation identity in the three-loop product

The new calculation keeps every input color index open. It does not choose
a gauge group, a representation, an invariant tensor for the
superpotential, or an input color contraction. Its geometric input is the
complete, fingerprinted first-current certificate in
[`nonrenormalizable-product`](../nonrenormalizable-product/).

The result is a factorization of the complete nonrational contribution by
the ordinary representation relation. The remaining ten color tensors have
explicit rational spinor coefficients for arbitrary source spinors.
[`FULL_SOURCE_TENSOR.md`](FULL_SOURCE_TENSOR.md) prints these coefficients;
[`quartic_tensor_certificate.json`](quartic_tensor_certificate.json) retains
the complete unreduced and reduced tensors, all Wick matchings, all source
routings, and the independent full-source replays.

## 1. Inputs and the meaning of the tensor indices

Initially polarize the pure-matter quartic superpotential. Its coefficient
tensor is defined by

\[
\bar W(\gamma)=\bar W^{i_1i_2i_3i_4}
\gamma_{1i_1}\gamma_{2i_2}\gamma_{3i_3}\gamma_{4i_4}.
\]

The four slots may carry different representations. They may also be
identified after the computation, with the exact polarization
multiplicities. Distinct polarization slots are not an assumption of
distinct physical matter species. The separate field \(\gamma_5\) below
is free, gauge neutral, and absent from \(\bar W\).

The inputs are

\[
A=\beta^{1i_1}(\partial_{+\dot a}\beta^{2i_2})\gamma_5,
\qquad B=b_A\beta^{3i_3}.
\]

Every one of \(i_1,i_2,A,i_3\) remains open. The output is
\(c^B\gamma_{4i_4}\gamma_5\), and \(B,i_4\) also remain open.
The free spinor index \(\dot a\) is never contracted. No quantum-closure or
observable-class statement is made for these inputs.

The interaction, contraction signs, derivative routing, and product
normalization are precisely those in
[`sources/LLC_main.tex`](../../../sources/LLC_main.tex), equations
`eq-B-10`, `eq-A-20`, `eq-W-fibre-explicit`, `eq-W-jet-edge`, `eq-D-lambda`,
and `eq-2-31`.

## 2. Complete Wick census and generator order

At third order there is one superpotential vertex and two cubic vertices.
Fixing the superpotential at vertex 0, the cubic vertices are 1 and 2, and
the marked inputs are 3 and 4. The complete census is

| Ordered cubic vertices | Valid matchings |
| --- | ---: |
| matter, matter | 20 |
| matter, gauge | 4 |
| gauge, matter | 4 |
| gauge, gauge | 0 |

All 44 excluded matchings retain their vanishing witnesses. In particular,
with two gauge vertices both input beta fields at vertex 3 must contract
directly to vertex 0, giving two parallel edges and a zero Laman kernel.

This census also accounts for matter species absent from the displayed
superpotential and input beta lines. At third order every antifield is
consumed. A single matter vertex of such a species has no partner for its
beta field. Two such vertices can only pair their beta and gamma fields
mutually, producing two parallel edges; self contractions are excluded.
Thus they add no admissible matching. Intermediate color components within
the interacting representations are all retained by the open matrix
products and are not excluded species.

For each valid matching the implementation applies the ordered contraction
operators from right to left. A beta--gamma contraction contributes its
explicit minus sign, and every left odd derivative contributes the parity
of the preceding odd fields. It then retains the derivative on the
original edge incident to \(\partial_{+\dot a}\beta^{2i_2}\).
The canonical graph relabeling only exchanges vertices 1 and 2. Its edge
permutation and the current orientation are included exactly.

Let the ghost paired to the original \(b_A\) carry the generator label
\(A\), and let the surviving ghost carry \(B\). A color path starting at
the superpotential and meeting the \(A\) interaction followed by the
\(B\) interaction gives

\[
(t_B)^{i}{}_{j}(t_A)^{j}{}_{k}\bar W^{\cdots k\cdots}.
\]

Thus the matrix order is the reverse of the path order. This follows
directly from the manuscript interaction
\(-igc^A(t_A)^j{}_i\gamma_j\beta^i\), without a choice of representation.
The matter--gauge matchings similarly give
\(iC^C{}_{AB}(t_C)^i{}_j\bar W^{\cdots j\cdots}\).
Exchanging the two ordered gauge ghosts changes the sign by their
antisymmetry; the factor \(-1/2\) of the gauge vertex is retained.

The result consists of ten two-generator words and two
\(iC^C{}_{AB}t_C\) words. The extraction visits each matter vertex exactly
once, retains every color path, and checks its endpoints against the
original Wick matching.

## 3. The three geometries in the commutator identity

Each displayed pair is an oriented edge, in the displayed edge order.
The final edge is the marked edge. The first and third source columns
give the locations of the surviving superpotential field and the free
spectator; the middle column gives the surviving ghost.

| Ordered graph | Source locations \((\gamma_4,c,\gamma_5)\) | Complete routed color weights |
| --- | --- | --- |
| \((01,03,04,12,14,23,34)\) | \((0,2,3)\) | \(+2y_{03,+\dot a}(t_Bt_A)\) on slot 1; \(-2y_{23,+\dot a}(t_Bt_A)\) on slot 2 |
| \((01,03,04,12,23,24,34)\) | \((0,1,3)\) | \(+2y_{03,+\dot a}(t_At_B)\) on slot 1; \(-2y_{23,+\dot a}(t_At_B)\) on slot 2 |
| \((01,03,04,12,13,24,34)\) | \((0,2,3)\) | \(-2y_{03,+\dot a}iC^C{}_{AB}t_C\) on slot 1; \(+2y_{13,+\dot a}iC^C{}_{AB}t_C\) on slot 2 |

For the first and third graphs the full conserved source rows are

\[
(\lambda_{\gamma_4},0,\lambda_c,\lambda_{\gamma_5},
-\lambda_{\gamma_4}-\lambda_c-\lambda_{\gamma_5}).
\]

For the second graph they are

\[
(\lambda_{\gamma_4},\lambda_c,0,\lambda_{\gamma_5},
-\lambda_{\gamma_4}-\lambda_c-\lambda_{\gamma_5}).
\]

The free source and current indices are the lower \(+\dot a\) indices
of the manuscript's `eq-D-lambda` and `eq-C-current-delta-integral`.
The current on \(03\), \(23\), or \(13\) is the actual original
propagator current. No source is removed when the ghost changes location.
The other canonical graphs give rational coefficients: one is generically
series-parallel, while the other two have a source-free cubic single-star
contraction to a series-parallel graph.

## 4. Why fifteen source assignments prove the full spinor identity

The bare three-loop density is homogeneous of source degree six. The
current insertion adds one, so the complete spinor has degree seven.
For a variation of only the spectator source, with the balancing source
at vertex 4 varied at the same time, every current changes by a scalar
edge response times the same spinor. The source-dependent part of the
six-row current determinant consequently has rank at most three. Its
degree in that spectator spinor is at most three, and the inserted
current raises this bound to four.

For two independent source spinors \(\lambda_{\gamma_4},\lambda_c\),
an invertible linear coordinate change sends them to \((1,0),(0,1)\).
The complete triangular grid

\[
\lambda_{\gamma_5}=(a,b),\qquad
a,b\in\{0,1,2,3,4\},\qquad a+b\leq4
\]

therefore determines both entries of the complete spinor polynomial
exactly. The result transforms by the cube of the determinant and one
copy of the spinor transformation. Polynomial continuation includes the
cases in which \(\lambda_{\gamma_4}\) and \(\lambda_c\) are dependent.

The certificate performs exact Newton interpolation, restores the
covariant expression in the original source spinors, and replays every
formal color word at five further complete source assignments, including
dependent and zero source spinors. It retains the full rational-plus-zeta
expression at every step. The geometric period values and their linear
dependence on the source polynomial are already exact in the reused
certificate; no assumption about independence of period values is needed.

## 5. The complete cancellation and the rational answer

For each of slots 1 and 2, the complete coefficient of the ordered word
\(t_At_B\) has a zeta term opposite to the coefficient of \(t_Bt_A\)
and opposite to the coefficient of \(iC^C{}_{AB}t_C\). No other color word
has a zeta term. Thus the entire nonrational expression has tensor factors

\[
\bigl([t_A,t_B]-iC^C{}_{AB}t_C\bigr)^{i_1}{}_{j_1}
\bar W^{j_1i_2i_3i_4},
\qquad
\bigl([t_A,t_B]-iC^C{}_{AB}t_C\bigr)^{i_2}{}_{j_2}
\bar W^{i_1j_2i_3i_4}.
\]

Both factors vanish by the representation relation. The two full spinor
polynomials multiplying them are printed in
[`FULL_SOURCE_TENSOR.md`](FULL_SOURCE_TENSOR.md). They are generally
nonzero when the spectator source is retained.

Substituting

\[
iC^C{}_{AB}t_C=t_At_B-t_Bt_A
\]

in the **complete** expression leaves ten explicitly computed rational
spinor coefficients. Gauge invariance of \(\bar W\) was not used. In
particular, the same identity holds before imposing the physical
gauge-invariance restriction on the superpotential tensor or contracting
the input color indices.

This computation uses the existing exact scalar periods, then establishes
a new universal tensor identity in their complete Wick sum. It does not
derive rationality of an arbitrary product from the Leibniz rule alone.

## 6. Physical normalization

The fixed-superpotential sum has the factor \(3/3!=1/2\). At three loops
the manuscript kernel equals the repository kernel. The product graph sign
is \(-1\), and the derivative at the head of the original propagator is
\(-2y_{+\dot a}\). For the two-matter words,

\[
(-1)\frac{3}{3!}(-2)(-i)^2=-1,
\qquad
\left(-\frac12\right)^7=-\frac1{128}.
\]

Consequently the rational tensor in the displayed manuscript
\(t_A,t_B\) convention, with each source spinor replaced by the derivative
on its own surviving field, has the common physical multiplier

\[
\boxed{\frac{g^2}{128(8\pi^2)^3}}.
\]

The same multiplier applies to the matter--gauge words when they are
written as \(iC^C{}_{AB}t_C\). Their original gauge factor and ghost-slot
signs are already in the coefficient table. The coefficient is that of
\(\hbar^3\); any superpotential coupling is part of the original
\(\bar W\) tensor. This convention reproduces both complete previously
computed abelian and SU(2) products after their input color contractions.

## 7. Arbitrary pure-matter superpotential degree

This extension has a separate exact certificate,
[`higher_W_certificate.json`](higher_W_certificate.json), generated by
[`higher_W.py`](higher_W.py). All additional surviving fields and their
original vertices are included.

Polarize a homogeneous term of \(\bar W(\gamma)\) of degree at least three, and fix the three
distinct tensor slots at which the three original input beta color lines
end. The slot labels are bookkeeping for the multilinear tensor. The
physical free spectator \(\gamma_5\) is a separate field absent from
\(\bar W\); the implementation reserves its label and uses 4,6,7,... for
the other polarization slots.

Terms of degree below three cannot supply the three distinct superpotential
slots at which the three input beta color lines must end in this
third-order channel.

There are only two finite cases to classify. If neither matter vertex
acts on a superpotential slot outside those first three, all contractions
are in the complete degree-three ledger. If one does, relabel that slot
as 4, giving the complete degree-four ledger. Two such matter vertices
are impossible in an admissible graph: both beta inputs at vertex 3
would then contract directly to vertex 0, producing a parallel pair.

The complete classification per fixed ordered choice of the three input
slots is

| Class | Number of matchings |
| --- | ---: |
| Same representation-commutator identity | 16 |
| Source-free cubic single-star contraction to an SP graph | 8 |
| Generic source-complete SP graphs | \(4(\deg\bar W-3)\) |

In the sixteen commutator matchings, every remaining superpotential field
is still at vertex 0. The already proved spinor identity therefore applies
with \(\lambda_{\gamma_4}\) replaced by the **sum of the independent
sources of all those surviving fields**. The source remains attached to
the full product of those fields. In the other eight non-SP matchings,
the source-free centre is cubic vertex 1 or 2; adding surviving
superpotential fields at vertex 0 cannot change that zero-source fact.
The generic SP graphs permit arbitrary retained sources. Their complete
finite current words integrate rationally by the existing exact local
ratio integrations.

The code independently enumerates degrees three through nine, preserves
every original sign and current in all sixteen commutator matchings, and
checks the complete counts. The all-degree conclusion follows from the
two finite cases above, rather than from extrapolating those replays.

To undo polarization, for a homogeneous term of degree \(\deg\bar W\),
divide the fully multilinear coefficient by \((\deg\bar W)!\). Its
diagonal is the original homogeneous polynomial. Differentiation of the
diagonal sums every differentiated-slot assignment, with its integer
multiplicity. The identities hold for the entire multilinear tensor
before that sum; their linear combination therefore holds for arbitrary
physical species and repeated fields. A formal power series is treated
coefficient by coefficient. Gauge invariance of the original
superpotential may then be imposed, but is not an ingredient of the
commutator cancellation.

This proves rationality of the complete third-order channel for arbitrary
pure \(\bar W(\gamma)\) in this input family. It does not give the same
claim for \(\bar W(\gamma,w)\), derivative superpotentials, other input
field contents, or higher input jet orders.

## 8. All-order rationality and termination for these inputs

Here the interaction is exactly the manuscript interaction in
`eq-B-10`, with pure \(\bar W(\gamma)\) and the free neutral spectator.
At order \(\hbar^n\) there are \(n\) interaction vertices and \(2n\)
Wick contractions. If \(r\) interaction vertices are superpotential
vertices, the available antifields are the \(n-r\) cubic-vertex
antifields plus the four input antifields. Since each contraction
consumes one antifield,

\[
2n\leq n-r+4,
\qquad n+r\leq4.
\]

If \(r=0\), each cubic vertex supplies exactly one contractible
Grassmann-even field, and the input \(b_A\) supplies one. The free
\(\gamma_5\) has no propagator to the interacting theory. Each
contraction consumes one Grassmann-even field, so

\[
2n\leq n+1,
\qquad n\leq1.
\]

Therefore the third-order coefficient has exactly one superpotential and
two cubic vertices, and every coefficient of order at least four is
zero. At orders zero, one, and two, every nonzero geometric graph is
series-parallel and has rational source-current coefficients. Together
with the complete third-order commutator calculation, this proves
rationality at every order for these fixed first-derivative inputs.

The coefficient ring is the rational ring generated by the original
interaction coefficients and superpotential tensor, with the universal
\((8\pi^2)^{-n}\) factors restored. Arbitrary numerical coupling values
are not asserted to be rational. This is an all-order rationality and
termination statement; it is distinct from an explicit evaluation of all
lower-order product formulas.

## 9. Independent controls and replay

The word extraction is checked against direct contraction of every
original color index in every matching, using arbitrary noncommuting
matrices of dimensions 2,3,2,3 and an arbitrary non-invariant
superpotential tensor. This check does not impose a representation
identity. A separate SU(2) calculation uses all three invariant quartic
pairings and all 729 open input/output color entries at every source
interpolation node. Another check uses a non-invariant W tensor with
honest representation matrices. Finally, contracting the two input color
pairs reproduces every full-source SU(2) result in the earlier certificate.

```bash
python explore/circuit-general/universal-quartic-product/quartic_tensor.py
python explore/circuit-general/universal-quartic-product/higher_W.py
python explore/circuit-general/universal-quartic-product/test_quartic_tensor.py
```

The first two commands regenerate and compare the entire deterministic
certificates. They reuse the published exact geometric data after checking
its implementation hashes; they do not build a period table. The complete
geometric certificate remains independently replayable with its original
command in the earlier directory.
