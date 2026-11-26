# An exact flux-derived bulk kernel reaching the glueball channel

26 November 2026. Independent calculation in the fixed-flux linear bulk model.

## Result and scope

There is an independently computed, nonzero bulk operation with one flux insertion,
two tangent scalar test insertions, one tangent area contraction, and an output in
the stable glueball channel. Its coefficient in the repository's mathematical
flux convention is exactly fixed below. A local primitive also identifies its
global relative Ward obstruction.

This is **not yet an independent evaluation of the N=4 Q1 amplitude**. The test
insertions are not full bulk cocycles, their identification with physical composite
sources has not been constructed, and the action/source normalization producing
the boundary factor g²/(64π²) is unresolved. No direct Q2 coefficient is obtained.

This distinction is particularly important here: the test insertions lie in the
regular scalar sector whose single flux images correspond to the removed B¹
central-shift sector. They cannot simply be called physical SU(N) external states.
The calculation supplies an explicit relative/source operation, not an operation
on the cohomology of those input states.

## 1. Existing conventions

Use the repository's tangent coordinates z1,z2 and normal coordinates w1,w2,w3.
The tangent bivector and flux class are

$$
P=\partial_{z_1}\wedge\partial_{z_2},\qquad
F=\frac{P}{w_1w_2w_3}.
$$

The bracket is the Schouten bracket, with [vector,function] equal to differentiation
of the function. Equivalently, for a polyvector of exterior degree p, it is

$$
\{\alpha,\beta\}=(-1)^{p-1}
\left(\partial_\Omega(\alpha\wedge\beta)
-\partial_\Omega\alpha\wedge\beta
-(-1)^p\alpha\wedge\partial_\Omega\beta\right).
$$

The flux has even normal Dolbeault degree two. Thus no additional odd coefficient
sign occurs in the calculations with ordinary scalar test functions below.
The full linear bulk differential is

$$
D=\bar\partial+t\partial_\Omega+N\{F,-\}.
$$

For computations on normal cohomology, F is the indicated Čech class. In smooth
formulas below it is P multiplied by the unit angular form U.

## 2. A precise obstruction to using only the linear flux map

For every regular polyvector α,

$$
\boxed{\{F,\alpha\}\text{ contains at least one tangent vector in every term}.}
$$

Proof: if the bracket differentiates the coefficient of α using P, one of the
two tangent vectors of P remains. If it differentiates the normal coefficient
1/(w1w2w3), the two tangent vectors of P remain. Differentiating constant vector
bases gives no extra terms. This proves the statement for arbitrary coefficients
and exterior polarizations, not just for monomials.

The scalar Laurent mode

$$
\frac{1}{w_1w_2w_3}
$$

therefore cannot be the image of a single linear flux bracket. At descendant
level zero it is not a boundary from t∂Ω either, because negative descendant
levels are absent. Its normal Dolbeault class is nonzero. The existing factorial
dictionary maps it to the exact two-form dz1∧dz2, the charge sector of the
single-trace gaugino bilinear. The fact that this two-form is exact does not make
it a boundary in the bulk descendant complex.

Thus identifying Q1's local glueball output with the linear operation N{F,-}
alone is a type mismatch. The two-source operation below removes both tangent
vectors and reaches the correct exterior type.

## 3. The nonzero two-source kernel

Let φ and ψ be ordinary holomorphic scalar test functions of z1,z2, independent
of the normal coordinates. Direct differentiation gives

$$
\{NF,\phi\}
=\frac{N}{w_1w_2w_3}
\left((\partial_{z_2}\phi)\partial_{z_1}
-(\partial_{z_1}\phi)\partial_{z_2}\right),
$$

and hence

$$
\boxed{
\{\{NF,\phi\},\psi\}
=-\frac{N}{w_1w_2w_3}
\left(\partial_{z_1}\phi\,\partial_{z_2}\psi
-\partial_{z_2}\phi\,\partial_{z_1}\psi\right).
}
$$

The unit example is

$$
\boxed{\{\{NF,z_1\},z_2\}=-\frac{N}{w_1w_2w_3}.}
$$

For φ=z1^a z2^b and ψ=z1^c z2^d, the complete answer is

$$
-N(ad-bc)\frac{z_1^{a+c-1}z_2^{b+d-1}}{w_1w_2w_3}.
$$

If a displayed exponent is negative, its coefficient is zero and the term is
absent. For formal plane waves exp(p1z1+p2z2), exp(q1z1+q2z2), the result is

$$
-N(p_1q_2-p_2q_1)
\frac{\exp((p_1+q_1)z_1+(p_2+q_2)z_2)}{w_1w_2w_3}.
$$

This has one tangent area and one factor of N. It was derived from the bulk
Schouten bracket and quantized flux input; the N=4 fan formula was not used.
The local output is D-closed: its divergence is zero, its normal Dolbeault
class is closed, and its bracket with F multiplies two degree-two normal
classes, which vanishes in normal degree four. When its coefficient is nonzero,
it is a nonzero primary scalar class by the type argument in Section 2.

### Coordinate comparison with the manuscript's spinors

For this comparison alone, identify z1=x^{+dot+}, z2=x^{+dot-} and retain the
repository flux P/(w1w2w3). The manuscript uses
∂_{+dot a}=−2∂/∂x^{+dot a} and epsilon^{dot+ dot-}=−1. Therefore, for test
waves exp(−p_{+dot a}x^{+dot a}/2), exp(−q_{+dot a}x^{+dot a}/2),

$$
\boxed{
\{\{NF,e^{-p_{+\dot a}x^{+\dot a}/2}\},
e^{-q_{+\dot a}x^{+\dot a}/2}\}
=\frac{N}{4}
(p_{+\dot a}q_+{}^{\dot a})
\frac{e^{-(p+q)_{+\dot a}x^{+\dot a}/2}}{w_1w_2w_3}.
}
$$

This factor 1/4 is the exact coordinate conversion for the stated bulk flux
convention. It is **not** an identification of action normalizations. The
boundary paper's holomorphic volume, physical field normalizations, ordered
odd ghost probes, and the bulk composite source map still have to be transported.

## 4. A local primitive and a global Ward obstruction

Use the already established unit angular representative on the punctured normal
space. Set r_i=bar(w_i)/|w|², and for cyclic (i,j,k) set

$$
b_i=r_j\bar\partial r_k-r_k\bar\partial r_j,\qquad
U=2\sum_{\mathrm{cyclic}(i,j,k)}r_i\bar\partial r_j\wedge\bar\partial r_k.
$$

These normalized forms obey

$$
\bar\partial b_i=w_iU,\qquad U\wedge b_i=0.
$$

On the chart wi≠0, the flux potential P b_i/wi satisfies
bar∂(P b_i/wi)=P U. Its twice-derived source operation is

$$
\{\{NPb_i/w_i,\phi\},\psi\}
=-N\frac{b_i}{w_i}
\left(\partial_{z_1}\phi\,\partial_{z_2}\psi
-\partial_{z_2}\phi\,\partial_{z_1}\psi\right).
$$

Applying D gives the smooth representative of the boxed kernel in Section 3.
Indeed bar∂ gives U, the divergence of the scalar is zero, and the extra flux
bracket has the factor U∧b_i=0. Thus this is a full linearized primitive on
each chart, including the descendant and background-flux terms.

The primitives do not combine into a global primitive when the tangent area is
nonzero. Their obstruction is the class −N(ad−bc) z1^{a+c−1}z2^{b+d−1}/(w1w2w3).
Equivalently, its normalized normal residue is precisely the tangent coefficient.

This identifies a concrete relative/source Ward building block. It does not yet
identify a heat-kernel-regulated bulk Feynman amplitude: the local primitive is
a patchwise Dolbeault potential, not a globally specified propagator with the
external boundary conditions and contact prescription required for Q1.

## 5. The normal Euler relation and its limitation

For the normal scaling vector E=Σ_i wi∂wi,

$$
\partial_\Omega E=3,\qquad \{F,E\}=3F,
\qquad \boxed{DE=3t+3NF.}
$$

More generally, for any diagonal vector with weights on all five coordinates,
its divergence is the sum of those weights, and its flux bracket is that same
sum times F. Traceless scalings are closed; the trace part relates the scalar
descendant and the flux variation.

Although the N-dependent term resembles an anomalous scaling Ward identity,
it is not the local Konishi identity. F contains the tangent bivector P; the
glueball mode in Section 2 is a scalar. The Euler relation belongs to the
flux/descendant sector, including the absolute flux variation removed when the
physical relative background is fixed. Calling DE=3(t+NF) the desired Q1
answer would skip precisely the source and descent identification at issue.

## 6. Why this does not establish higher Q_n

The scalar tests φ,ψ are not D cocycles unless their tangent derivatives vanish.
The existing regular-mode dictionary identifies their flux images with exact
one-forms dφ,dψ, the broken central-shift sector. Therefore the operation in
Section 3 is not a cohomology product of physical SU(N) single-graviton states.

Moreover, operations that multiply two flux representatives at the same point
vanish by their normal Dolbeault degree. Higher Q_n cannot be recovered merely
by taking powers of this local flux bivector. Separated flux vertices with
propagators, quantum vertices, and the actual composite source transport have
to be calculated.

No value has been assigned to the remaining normalization. In particular,
the boundary g²N/(64π²) coefficient is not used to define any term above.

## 7. Exact verification

The standalone standard-library script flux_derived_kernel.py imports no
boundary graph evaluator. It computes exterior products, holomorphic divergence,
and Schouten brackets directly with exact rational Laurent coefficients.

Its saved certificate flux_derived_kernel.json records:

- 1296 polynomial source pairs, including 1126 nonzero pairs;
- 3125 diagonal scaling vectors;
- 7776 regular polyvectors checking the linear flux type obstruction.

The unrestricted formulas are proved in the text. The finite tests independently
check their signs, exterior types, and coefficients.

## Primary sources and prior exact inputs

- Costello–Li, *Twisted supergravity and its quantization*,
  https://arxiv.org/html/1606.00365v1, Section 9.1: proposed flux background and
  linear differential.
- Costello–Li, *Quantization of open-closed BCOV theory I*,
  https://arxiv.org/pdf/1505.06703, Sections 3 and 7: open/closed disk maps and
  annulus cancellation in the space-filling model. This theorem does not by
  itself evaluate the punctured D3 composite amplitude.
- Li–Wang, *On the exact quantization of BCOV Theory on Calabi-Yau Manifolds*,
  https://arxiv.org/html/2610.21927v1: the stated higher-genus bulk interaction
  construction does not supply the D3 external/source map used in this task.
- Repository at d209471: explore/twisted-holography-N4/NOTE.md Sections 2 and 6
  for the Laurent dictionary and the B¹ image; README.md, “Simultaneous tangent
  shears of a Laurent scalar,” for the normalized U,b_i descent identities.

