# Further twisting: the exact normal-momentum complex and its locality test

Date: 2026-11-26.

This calculation concerns the strict field-strength model of twisted
eleven-dimensional supergravity. It determines the complete linear
further-twist complex at every holomorphic momentum, its polynomial
cohomology near the residual plane, and the dependence on the cubic
coupling. The distinction between fixed momentum and formal position
dependence is essential for using the further twist in an anomaly argument.

## 1. Starting point and the precise question

The established two-loop calculation has the nonzero integrated pairing

$$
-\frac{63g}{256\pi^{10}}.
$$

The existing further-twist calculation shows that this vertex vanishes when
evaluated on the residual planar Hamiltonian representatives. It does not
construct a primitive for the entire shifted local density.
[Source calculations](../collaboration/current_complex/README.md#constant-two-form-background-and-the-further-twist).

The published further-twist background is

$$
\gamma_{\mathrm{nm}}
=\frac12(z_1\,dz_2-z_2\,dz_1),
\qquad
\partial\gamma_{\mathrm{nm}}=dz_1\wedge dz_2.
$$

Raghavendran–Saberi–Williams identify the resulting full classical theory
with Poisson Chern–Simons theory on
$\mathbb C^2\times\mathbb R^7$; their formal symmetry calculation uses
primitives with dependence on the three transverse coordinates.
[Sections 4.2–4.3, equation (45), Theorem 4.1](https://arxiv.org/html/2111.03049v1#S4).
Hahner–Paquette–Raghavendran also describe the three additional holomorphic
and three additional antiholomorphic directions made topological by this
further twist. [Section 4.3](https://arxiv.org/html/2506.02246v3#S4.SS3).

Here we retain every vector and closed-two-form polarization at a fixed
momentum. This gives a direct test of whether the formal simplification
can be implemented by constant-coefficient differential operators with
finite derivative order. The answer requires the full normal Koszul
complex below.

## 2. The differential on the complete wave spaces

Let $p$ be the holomorphic momentum covector. An even element of the
exceptional super Lie algebra is a vector wave with polarization $X$;
an odd element is a two-form wave with polarization $F$. Their constraints
are

$$
p(X)=0,\qquad p\wedge F=0.
$$

Even and odd in this note refer to the super Lie algebra grading,
before the shift used for physical superfields. The common wave factor
is $\exp(p\cdot z)$.

Use the five-volume
$\Omega=dz_1\wedge dz_2\wedge dz_3\wedge dz_4\wedge dz_5$.
For a four-form $\alpha$, the vector $\alpha^\sharp$ is defined by
$\iota_{\alpha^\sharp}\Omega=\alpha$.
The further-twist differential induced by the background
$dz_1\wedge dz_2$ is

$$
\boxed{
\delta X=-p\wedge\iota_X(dz_1\wedge dz_2),
\qquad
\delta F=g\bigl(dz_1\wedge dz_2\wedge F\bigr)^\sharp .
}
$$

The first sign is fixed by
$[dz_1\wedge dz_2,X]=-\mathcal L_X(dz_1\wedge dz_2)$.
The coefficient $g$ is retained from the original cubic interaction.
Both maps preserve the stated constraints. For example,

$$
p\wedge dz_1\wedge dz_2\wedge F=0
\quad\Longrightarrow\quad
p\!\left((dz_1\wedge dz_2\wedge F)^\sharp\right)=0.
$$

The two square-zero identities are immediate: the one-form
$\iota_X(dz_1\wedge dz_2)$ is planar, while
$(dz_1\wedge dz_2\wedge F)^\sharp$ is transverse.

### Theorem 1. Complete fixed-momentum cohomology

Assume $g\ne0$. Write $p\vert_{\mathbb C^3}$ for the restriction of $p$
to the three transverse coordinate directions.

| Momentum | Even wave dimension | Odd wave dimension | Even-to-odd rank | Odd-to-even rank | Even cohomology | Odd cohomology |
|---|---:|---:|---:|---:|---:|---:|
| $p\vert_{\mathbb C^3}\ne0$ | 4 | 4 | 2 | 2 | 0 | 0 |
| $p\vert_{\mathbb C^3}=0$, $p\ne0$ | 4 | 4 | 0 | 0 | 4 | 4 |
| $p=0$ | 5 | 10 | 0 | 3 | 2 | 7 |

**Proof.** For $p\ne0$, the first constraint has a four-dimensional
kernel. Every closed two-form is $p\wedge a$, with the one-form $a$
defined modulo the line spanned by $p$, so the odd wave space is also
four-dimensional.

If $p\vert_{\mathbb C^3}\ne0$, the condition
$p\wedge\iota_X(dz_1\wedge dz_2)=0$ forces the planar part of $X$ to
vanish. Its remaining transverse polarization satisfies one nonzero
linear condition, leaving dimension two. Conversely the transverse
two-form part of $p\wedge a$ has the form
$p\vert_{\mathbb C^3}\wedge a|_{\mathbb C^3}$ and supplies precisely these
two vector polarizations. Both kernel/image equalities follow.

If $p$ is nonzero and planar,

$$
p\wedge\iota_X(dz_1\wedge dz_2)
=p(X)\,dz_1\wedge dz_2=0.
$$

Every closed $F=p\wedge a$ has no purely transverse two-form part.
Consequently $\delta F=0$ as well.

At $p=0$, the constraints vanish. The map on vectors is zero.
The map on two-forms identifies the three-dimensional transverse
two-form space with the three transverse vector directions.
Its kernel has dimension $10-3=7$ and its cokernel has dimension
$5-3=2$. This proves every entry. $\square$

At $g=0$ the odd-to-even map is zero. The corresponding three cohomology
pairs are $(2,2)$, $(4,4)$ and $(5,10)$.

## 3. Exact primitives away from zero transverse momentum

These formulas establish acyclicity without discarding any polarization.
Let
$\Omega_\perp=dz_3\wedge dz_4\wedge dz_5$.
Choose a transverse vector $n$ with $p(n)=1$.

For a closed even polarization $X$, Theorem 1 makes $X$ transverse.
The closed two-form

$$
F=\frac1g\,p\wedge\iota_n\iota_X\Omega_\perp
$$

is an exact primitive:

$$
p\wedge F=0,\qquad
p\vert_{\mathbb C^3}\wedge\iota_n\iota_X\Omega_\perp
=\iota_X\Omega_\perp,\qquad
\delta F=X.
$$

For a closed odd polarization $F$, its transverse two-form part
vanishes. Thus $a=\iota_nF$ is a planar one-form, and the closure
condition gives $F=p\wedge a$. Let $Y$ be the unique planar vector
satisfying $\iota_Y(dz_1\wedge dz_2)=-a$. Then

$$
X=Y-p(Y)n,\qquad
p(X)=0,\qquad
\delta X=F.
$$

For instance $n=\partial_{z_j}/p_j$ when a transverse component
$p_j$ is nonzero. These exact formulas expose inverse transverse
momenta, and the even-cycle primitive also uses $g^{-1}$.
They are valid localized identities; locality must be checked before
using them to construct counterterms.

## 4. The whole polynomial normal complex

We now compute the dependence on the transverse momenta, instead of
only evaluating individual fibers. On a nonzero planar momentum
chart, an area-preserving planar change of frame puts

$$
p=dz_1+q_3dz_3+q_4dz_4+q_5dz_5.
$$

Thus $q=p-dz_1$ is the normal momentum covector.
Write $R=\mathbb C[q_3,q_4,q_5]$ for the ring of polynomial transverse
momenta. For the moment set $g=1$.

The vector wave space is free over $R$, with the basis

$$
\partial_{z_2},\qquad
\partial_{z_j}-q_j\partial_{z_1}\quad(j=3,4,5).
$$

The two-form wave space is free over $R$, with the basis

$$
p\wedge dz_2,\qquad p\wedge dz_j\quad(j=3,4,5).
$$

Every vector and form displayed satisfies its full constraint.
For the three transverse coefficients use the volume
$dz_3\wedge dz_4\wedge dz_5$ to identify two-forms with vectors.
The differential, after ordering these four summands by exterior
degree, is exactly

$$
\boxed{
0\longrightarrow R
\xrightarrow{\ q\ } R^3
\xrightarrow{\ q\wedge\ } R^3
\xrightarrow{\ q\cdot\ } R
\longrightarrow0 .
}
$$

The first scalar basis element is $-\partial_{z_2}$ and the last is
$p\wedge dz_2$; this accounts for the first differential's sign.
The successive parities are even, odd, even, odd. Thus this is the
three-variable Koszul complex with its integer degrees folded by
parity. It describes all eight wave polarizations.

### Theorem 2. Normal polynomial cohomology

$$
\boxed{
H^{\mathrm{even}}=0,\qquad
H^{\mathrm{odd}}
=R/(q_3,q_4,q_5).
}
$$

Here $H$ denotes the kernel of the outgoing map divided by the image
of the incoming map. The surviving class is the last scalar
$p\wedge dz_2$ at polynomial degree zero.

**Constructive proof.** Identify the displayed sequence with
$R\otimes\Lambda^\bullet\mathbb C^3$, with differential given by
wedge multiplication by $q$.
On a homogeneous polynomial of degree $m$ and an exterior tensor of
degree $k$,

$$
(q\wedge)\left(\sum_{j=3}^5
  \iota_{\partial_{z_j}}\partial_{q_j}\right)
+\left(\sum_{j=3}^5
  \iota_{\partial_{z_j}}\partial_{q_j}\right)(q\wedge)
=(m+3-k)\,\mathrm{id}.
$$

For a closed element, division by the positive integer $m+3-k$
therefore supplies a polynomial primitive. The only exception is
$m=0,k=3$. In top degree every positive-degree polynomial is in
the ideal $(q_3,q_4,q_5)$, and the constant survives.
This proves the claimed cohomology in every polynomial degree.

The displayed homotopy differentiates with respect to momentum:
it is linear over $\mathbb C$, but is not linear over $R$.
In position space it introduces coordinate dependence. Consequently
this constructive polynomial proof is compatible with the
nonexistence of the uniform constant-coefficient contraction discussed
next.

## 5. Why the zero-momentum limit retains eight modes

Setting all three normal momenta to zero in the free complex
annihilates every differential. The resulting exterior dimensions
are

$$
1,\quad3,\quad3,\quad1.
$$

Even and odd dimensions are therefore $1+3=4$ and $3+1=4$.
By contrast, setting normal momenta to zero only after taking the
polynomial cohomology would retain a single odd line.
These operations differ.

The precise algebraic correction is the derived tensor product.
In this context it means retaining a free resolution when imposing
$q_3=q_4=q_5=0$; its extra homology records the failure of ordinary
tensor product to preserve exactness. The Koszul resolution gives

$$
\operatorname{Tor}^{R}_j
\bigl(R/(q_3,q_4,q_5),\,R/(q_3,q_4,q_5)\bigr)
=\Lambda^j\mathbb C^3,\qquad 0\le j\le3,
$$

where $\operatorname{Tor}_j$ is that additional homology in resolution
degree $j$. The four dimensions are exactly $1,3,3,1$.
At the level of the linear symbol, these multiplicities supply
the three holomorphic directions that become topological.
Combining them with the three transverse
antiholomorphic directions and the original real direction gives
the seven topological directions in the published maximal twist.

There is also an exact one-parameter test. Put
$q_3=t,\ q_4=q_5=0$. Over $\mathbb C[t]$, the four nonzero arrows
are

$$
\begin{aligned}
\delta(\partial_{z_2})&=-t(p\wedge dz_3),&
\delta(\partial_{z_3}-t\partial_{z_1})&=t(p\wedge dz_2),\\
\delta(p\wedge dz_4)&=t\partial_{z_5},&
\delta(p\wedge dz_5)&=-t\partial_{z_4}.
\end{aligned}
$$

Consequently

$$
H^{\mathrm{even}}=(\mathbb C[t]/t)^2,\qquad
H^{\mathrm{odd}}=(\mathbb C[t]/t)^2.
$$

All four contractible pairs at $t\ne0$ require division by $t$.
No matrix of polynomials, nor a matrix regular at $t=0$, can give
an identity $\delta h+h\delta=\mathrm{id}$ on these pairs:
the left side vanishes at $t=0$ and the right side does not.
Equivalently a constant-coefficient finite-order differential
operator cannot perform that uniform elimination.

This is a precise limit on the use of the formal residual
Hamiltonian algebra. It does not exclude more general local
constructions with position-dependent coefficients or the full
classical equivalence proved in the source.

More generally, restrict normal momentum to an $r$-dimensional
linear subspace, with $r=0,1,2,3$. A change of transverse basis
leaves $r$ independent momentum variables and $3-r$ zero entries
in the Koszul sequence. Its cohomology is the top cohomology of
the $r$-variable Koszul complex tensored with the exterior algebra
on the remaining $3-r$ directions. Over the quotient by the
$r$ remaining momentum variables, the free generator counts are

| Independent normal momenta $r$ | Even generators | Odd generators |
|---:|---:|---:|
| 3 | 0 | 1 |
| 2 | 1 | 1 |
| 1 | 2 | 2 |
| 0 | 4 | 4 |

This gives every linear normal-momentum specialization, including
the one-parameter family and the zero-normal-momentum fiber.

## 6. Retaining the cubic coupling

The anomaly computation treats $g$ explicitly. To compare the
categories, let now
$R=\mathbb C[g,q_3,q_4,q_5]$.
The complex of free modules has the middle map multiplied by
$g$:

$$
0\longrightarrow R
\xrightarrow{\ q\ }R^3
\xrightarrow{\ g\,q\wedge\ }R^3
\xrightarrow{\ q\cdot\ }R
\longrightarrow0 .
$$

Let $\ker(q\cdot)$ mean the polynomial transverse vectors whose
contraction with $q$ is zero. The unscaled Koszul calculation
identifies this kernel with the image of $q\wedge$.
Because multiplication by $g$ is injective on the free modules,

$$
\boxed{
H^{\mathrm{even}}
=\frac{\ker(q\cdot)}{g\,\ker(q\cdot)},
\qquad
H^{\mathrm{odd}}=R/(q_3,q_4,q_5).
}
$$

Thus even cohomology has a sector annihilated by $g$. It disappears
after inverting $g$, but remains in the polynomial or formal
power-series category before that inversion.
Its normal coefficient degree $m$ has dimension

$$
3\binom{m+2}{2}-\binom{m+3}{2}=m(m+2),\qquad m\ge0.
$$

Indeed contraction with $q$ maps the three degree-$m$ vector
coefficients onto all degree-$(m+1)$ scalar polynomials. This
provides the exact multiplicity in every normal coefficient
degree; $m=0$ gives zero.
On the one-parameter normal line,

$$
\boxed{
H^{\mathrm{even}}=(\mathbb C[g,t]/(gt))^2,
\qquad
H^{\mathrm{odd}}=(\mathbb C[g,t]/t)^2.
}
$$

These are complete cohomology modules, not numerical rank
extrapolations. At fixed $g\ne0$ they reduce to the preceding
answer. The primitive formulas in Section 3 display exactly where
the division by $g$ enters.

The existing anomaly is itself proportional to $g$, so this
torsion statement alone does not prove that a further-twisted
counterterm is impossible. It specifies a divisibility check that
an actual counterterm construction must satisfy.

## 7. Consequence for the anomaly problem

The further-twist consistency test has two separate outcomes.

* The full normal wave complex reproduces the exterior multiplicities
  required by three additional topological holomorphic directions.
  This provides an exact linear check of the maximal-twist
  interpretation.
* The vanishing of the original two-loop vertex on planar
  Hamiltonian inputs cannot by itself supply a
  constant-coefficient finite-jet primitive. The complete normal
  Koszul resolution, its zero-momentum classes, and its coupling
  torsion must be retained in a local calculation.

An exactness theorem for the shifted two-loop density still
requires the nonlinear local differential and total-derivative
relations. The calculation here neither proves nor disproves
that theorem. It does establish the full linear symbol input
and the exact failure of the naive uniform contraction.

This note does not calculate the spectrum of the complete kinetic
operator $d_\tau+\bar\partial+\delta$, an on-shell scattering
amplitude, or a new loop integral. Its statements concern
the further-twist differential of the strict field-strength
model and polynomial dependence on holomorphic momentum.

## 8. Independent certificate

Run from the repository root:

```bash
python explore/sugra-11d-qme/developments_2026_11/further_twist_symbols.py
```

The script uses an independent exterior-algebra implementation and
SymPy exact arithmetic, with no imports from the existing graph or
cohomology code. It derives the maps from the five-volume and the
background two-form, checks both constraints and both square-zero
identities, and obtains the full polynomial matrices before taking
any momenta.

The saved [exact output](further_twist_symbols.json) records:

* all 243 integer momenta in $\{-1,0,1\}^5$, plus two rational-primitives
  controls with larger integer momenta;
* all three strata and their complete dimensions;
* 472 even-cycle and 472 odd-cycle primitive identities;
* the full symbolic normal complex and the localized contraction,
  including the coupling factor;
* 672 polynomial/exterior basis identities through polynomial
  degree six for the explicit homotopy.

The all-degree module statements are proved analytically above.
The bounded certificate independently checks the formulas and signs.
No PDF is generated.

## References and provenance

1. Surya Raghavendran, Ingmar Saberi, Brian R. Williams,
   *Twisted eleven-dimensional supergravity*,
   [arXiv:2111.03049v1](https://arxiv.org/html/2111.03049v1),
   Sections 2.2 and 4.2–4.3.
2. Fabian Hahner, Natalie M. Paquette, Surya Raghavendran,
   *Higher Symmetries in Twisted Eleven-Dimensional Supergravity*,
   [arXiv:2506.02246v3](https://arxiv.org/html/2506.02246v3),
   Sections 2.2.4 and 4.3.
3. Repository input: [current-complex calculations](https://github.com/lbotao399-stack/holomorphic-topological-feynman-diagrams/blob/aa0d2f1725dff135494eef132f5e92a8c0401a42/explore/sugra-11d-qme/collaboration/current_complex/README.md).
   The previously integrated two-loop coefficient is used only
   to state the anomaly question; it is not an input to the new
   exterior-algebra certificate.
