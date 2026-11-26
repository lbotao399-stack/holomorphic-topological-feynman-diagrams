# Linear bulk states for the all-order joining test, and the exact primary-tree obstruction

## Result

The all-order boundary test in
[the boundary report](../boundary/REPORT.md) uses $n+1$ individually
protected quadratic traces. Their complete linear Laurent state spaces
are determined explicitly below. Every field in these four charge
sectors contains the same tangent bivector.

**Every rooted tree built from the local primary bulk vertices vanishes
on these directly attached Laurent representatives.** This remains true
with any number of separate flux insertions and arbitrary linear
propagators on the internal edges. A local positive-genus vertex does
not change this conclusion when the sewing graph remains a tree.

This is a concrete obstruction to identifying the nonzero boundary
joining coefficient with that bare bulk tree. It does not establish
vanishing of the physical string amplitude. External bulk-to-boundary
kernels, the map of composite observables, and relative open/closed
contact terms are additional data. Graphs with geometric cycles are
also outside the tree argument.

We also determine, for every boundary order $n\ge2$, the unique
**stable protected one-particle channel** at the output charge of this
test and its exact Laurent representative. Thus the input fields and
the possible protected output state are explicit at every order.
The raw boundary joining output is a cochain; it is not automatically
that protected state.

The new calculation uses the dictionary and local bulk vertices at
repository revision
[29ef4be3c92b8055ac51067357524c1048e6cda0](https://github.com/lbotao399-stack/holomorphic-topological-feynman-diagrams/tree/29ef4be3c92b8055ac51067357524c1048e6cda0).
The previous notes already proved the corresponding cubic-tree
vanishing for the three-input, two-loop example. The result here
supplies the whole input basis and proves the statement at every order.

## 1. Protected boundary constituents and their charges

Retain the paper's scalar fields $\gamma^I$, fermionic fields $\beta_I$,
field $b$, and coupling $g$. The four boundary operators are

$$
\begin{aligned}
G_1&=\operatorname{Tr}(\beta_1\gamma^2),&
G_2&=\operatorname{Tr}(\beta_2\gamma^3),\\
G_3&=\operatorname{Tr}(\beta_3\gamma^1),&
H&=\operatorname{Tr}(b\gamma^1-ig\beta_2\beta_3).
\end{aligned}
$$

The first three are odd and the fourth is even. They are components
of the quadratic graviton $\operatorname{Tr}(dC)^2$ and are individually
annihilated by every positive $Q_n$, as well as by $Q_0$.
Here $d$ is the exterior derivative on the existing boundary
superspace; it acts on the generating superfield $C$.
The component identities, with left odd derivatives, are

$$
\begin{aligned}
G_1&=-\frac12[\theta_3]\,
 \operatorname{Tr}\bigl((\partial_{\theta_2}C)^2\bigr),\\
G_2&=-\frac12[\theta_1]\,
 \operatorname{Tr}\bigl((\partial_{\theta_3}C)^2\bigr),\\
G_3&=-\frac12[\theta_2]\,
 \operatorname{Tr}\bigl((\partial_{\theta_1}C)^2\bigr),\\
H&=\frac{ig}{2}[\theta_2\theta_3]\,
 \operatorname{Tr}\bigl((\partial_{\theta_1}C)^2\bigr).
\end{aligned}
$$

The square brackets extract the indicated coefficient. The component
conventions are fixed in
[the manuscript, equation (A.13)](https://github.com/lbotao399-stack/holomorphic-topological-feynman-diagrams/blob/29ef4be3c92b8055ac51067357524c1048e6cda0/sources/LLC_main.tex).

Their three odd-coordinate charges are

$$
\begin{array}{c|cccc}
 &G_1&G_2&G_3&H\\ \hline
(\theta_1,\theta_2,\theta_3)&
(0,2,1)&(1,0,2)&(2,1,0)&(2,1,1).
\end{array}
$$

Every constituent has letter number two and zero tangent jet counts.
The order-$n$ ordinary product is $G_1G_2G_3H^{n-2}$, in this order.
Its maximal connected joining coefficient differentiates one letter
in each of its $n+1$ factors; the boundary report computes that
coefficient without replacing the factors by nonclosed traces.

## 2. Exact two-form bases

In the linear dictionary, a length-two single-particle class is an
exact two-form on $\mathbb C^{2|3}$. The odd coordinates are
$\theta_1,\theta_2,\theta_3$, while their differentials are even.
All four charges have zero tangent counts, so their form bases have
no tangent differential.

Choose the following exact basis forms, each with the first displayed
coefficient equal to one:

$$
\begin{aligned}
d(\theta_2\theta_3\,d\theta_2)
 &=\theta_3(d\theta_2)^2-\theta_2d\theta_2d\theta_3,\\
d(\theta_3\theta_1\,d\theta_3)
 &=\theta_1(d\theta_3)^2-\theta_3d\theta_1d\theta_3,\\
d(\theta_1\theta_2\,d\theta_1)
 &=\theta_2(d\theta_1)^2-\theta_1d\theta_1d\theta_2,\\
d(\theta_1\theta_2\theta_3\,d\theta_1)
 &=\theta_2\theta_3(d\theta_1)^2
 -\theta_1\theta_3d\theta_1d\theta_2
 +\theta_1\theta_2d\theta_1d\theta_3.
\end{aligned}
$$

At each of the first three charges the one-form space has dimension
one and the two-form space dimension two; the derivative has rank one.
At the fourth charge these dimensions are one and three, again with
rank one. Therefore each exact two-form space is one-dimensional.
The basis choice above fixes the fields below exactly. It does not
silently fix the physical source/operator pairing or the nonlinear map
of products. The vanishing theorem applies to every normalization
within these complete one-dimensional spaces.

## 3. The four Laurent states and the full linear differential

Let $z_1,z_2$ be the tangent coordinates and $w_1,w_2,w_3$ the normal
coordinates. As in the existing flux calculation, write

$$
P=\partial_{z_1}\wedge\partial_{z_2}.
$$

The Laurent/form dictionary sends a monomial with the two tangent
vectors present to

$$
\Phi\left(
\frac{P\wedge\partial_{w_{j_1}}\wedge\cdots\wedge\partial_{w_{j_k}}}
{w_1^{a_1+1}w_2^{a_2+1}w_3^{a_3+1}}
\right)
=\frac{(-1)^{a_1+a_2+a_3}}{a_1!a_2!a_3!}
\theta_{j_1}\cdots\theta_{j_k}
(d\theta_1)^{a_1}(d\theta_2)^{a_2}(d\theta_3)^{a_3}.
$$

Here $\Phi$ is the existing dictionary, and the indices $j_1<\cdots<j_k$
give the exterior ordering. Its inverse on the four chosen form
bases is

$$
\begin{aligned}
G_1\text{ sector}:\quad&
P\wedge\left(
\frac{2\partial_{w_3}}{w_1w_2^3w_3}
-\frac{\partial_{w_2}}{w_1w_2^2w_3^2}\right),\\
G_2\text{ sector}:\quad&
P\wedge\left(
\frac{2\partial_{w_1}}{w_1w_2w_3^3}
-\frac{\partial_{w_3}}{w_1^2w_2w_3^2}\right),\\
G_3\text{ sector}:\quad&
P\wedge\left(
\frac{2\partial_{w_2}}{w_1^3w_2w_3}
-\frac{\partial_{w_1}}{w_1^2w_2^2w_3}\right),\\
H\text{ sector}:\quad&
P\wedge\left(
\frac{2\partial_{w_2}\wedge\partial_{w_3}}{w_1^3w_2w_3}
-\frac{\partial_{w_1}\wedge\partial_{w_3}}{w_1^2w_2^2w_3}
+\frac{\partial_{w_1}\wedge\partial_{w_2}}{w_1^2w_2w_3^2}
\right).
\end{aligned}
$$

The Laurent coefficients represent degree-two normal Dolbeault
cohomology, or equivalently the three-chart Čech cocycles.
Let $U$ be the unit normal angular representative of
$1/(w_1w_2w_3)$. Use the actual smooth representative

$$
\frac{1}{w_1^{a_1+1}w_2^{a_2+1}w_3^{a_3+1}}
\quad\longleftrightarrow\quad
\frac{(-1)^{a_1+a_2+a_3}}{a_1!a_2!a_3!}
\partial_{w_1}^{a_1}\partial_{w_2}^{a_2}\partial_{w_3}^{a_3}U.
$$

For instance the first field becomes
$P\wedge((\partial_{w_2}^2U)\partial_{w_3}
-(\partial_{w_2}\partial_{w_3}U)\partial_{w_2})$.
Its divergence vanishes by commutation of derivatives, as an
identity of smooth forms on the puncture. The other three
cancellations have the same form. Thus the following closure
statements hold on these explicit representatives, as well as
on their Čech classes.

Let $\partial_\Omega$ be divergence for the holomorphic five-volume.
The linear descendant differential at flux strength $N$ is

$$
D=\bar\partial+t\partial_\Omega+N\{F,-\},
\qquad
F=\frac{P}{w_1w_2w_3}
$$

in the cohomological Laurent notation. The parameter $t$ records
descendant level, and $\{F,-\}$ is the Schouten bracket.
Every displayed field is divergence free. For example, on the
first field its two normal divergences are

$$
-\frac{2P}{w_1w_2^3w_3^2}
+\frac{2P}{w_1w_2^3w_3^2}=0.
$$

The certificate checks the whole exterior expression in each case.
The bracket with the flux is zero: the product of two normal
degree-two representatives has degree four, and the same result
also follows from their common tangent bivector and tangent
independence. Therefore all four are closed under the full displayed
linear differential for every $N$ and $t$.

They are nonzero classes at descendant level zero. At that level
the incoming flux image is the space of exact one-forms under
$\Phi$, whereas these are nonzero exact two-forms. A free descendant
divergence has a factor of $t$ and cannot kill a level-zero class.
Calling the boundary forms “exact” thus does not make these bulk
states boundaries.

## 4. All-order primary-tree vanishing

Consider a rooted sewing graph with no geometric cycle. Attach
any number of the four displayed Laurent fields directly to its
external inputs. Also allow any number of separate copies of $NF$
as external background leaves. The internal edges may be arbitrary
linear propagators.

At a local primary bulk vertex, the scalar momentum polynomial
acts by holomorphic derivatives on the coefficients of the
incoming polyvectors. Their exterior product is then formed,
followed by the root divergence when the vertex is written as a
Taylor operation. This is the structure of the fixed local
worldsheet vertices used in the repository.

Every holomorphic derivative of every input still has $P$ as an
exterior factor:

$$
\partial_z^{\,a}\partial_w^{\,b}(P\wedge Y)
=P\wedge\partial_z^{\,a}\partial_w^{\,b}Y.
$$

The integers in the multi-indices $a,b$ are nonnegative derivative
orders. For two such differentiated inputs,

$$
(P\wedge Y_1)\wedge(P\wedge Y_2)=0,
\qquad
\partial_\Omega\bigl((P\wedge Y_1)\wedge(P\wedge Y_2)\bigr)=0.
$$

This is the full product, before selecting any output component.

Every finite rooted tree whose nonlinear vertices have at least
two incoming flags has a nonlinear vertex furthest from the root.
All of that vertex's incoming flags are external leaves. At least
two of them carry $P$, so the vertex is identically zero. Its image
under any subsequent linear propagator is zero, and hence the
whole graph is zero.

Explicit flux insertions satisfy the same hypothesis because $F$
contains $P$. Treating the flux bracket as a unary background
operator does not evade the argument: it vanishes on the displayed
external fields, and acts linearly on zero after the first nonlinear
vertex.

This proves the statement at every external valence and every
choice of local vertex genus. Higher-genus local contact vertices
can occur in a sewing graph with no cycles; their holomorphic
coefficient derivatives leave $P$ unchanged, so they obey the
same zero product. At genus zero with primary external fields the
vertices are cubic, giving the usual binary rooted-tree case.

### What the statement excludes, and what it does not compute

The nonzero connected boundary coefficient at order $n$ joins
$n+1$ protected trace factors into one trace and has the topology
of an $(n+2)$-point bulk tree. The linear state dictionary by
itself gives precisely the fields above, and direct attachment to
the local bulk vertices gives zero.

Therefore the equality cannot be obtained by that direct attachment.
A cochain map for the ordinary boundary product must include the
appropriate composite-operator dressing, or the relevant physical
amplitude must use its actual external kernels and relative
open/closed contact terms. The present proof places no restriction
on such kernels. It does not evaluate sewing graphs with cycles,
for which the furthest-vertex argument is unavailable.

## 5. The possible protected output channel at every order

The total odd-coordinate charge of $G_1G_2G_3H^{n-2}$ is

$$
(3,3,3)+(n-2)(2,1,1)=(2n-1,n+1,n+1).
$$

The boundary $Q_n$ removes $(n,n,n)$ from these charges, adds
$(n,n)$ tangent derivative counts, and changes the letter number
by $1-n$. Hence its output sector is exactly

$$
(\theta_1,\theta_2,\theta_3;n_+,n_-)_Y
=(n-1,1,1;n,n)_{n+3}.
$$

In the stable single-particle dictionary, this is the highest
possible form degree at its charge. Both tangent differentials
must occur, and no undifferentiated odd coordinate may occur.
The only form monomial is

$$
z_1^{n-1}z_2^{n-1}\,
dz_1\wedge dz_2\,
(d\theta_1)^{n-1}d\theta_2d\theta_3.
$$

It is a nonzero exact form, with explicit primitive

$$
\frac{z_1^n z_2^{n-1}}{n}\,
dz_2\,(d\theta_1)^{n-1}d\theta_2d\theta_3.
$$

For a scalar Laurent field the dictionary has the additional
complementary-tangent sign minus. Its inverse is consequently

$$
\boxed{
(-1)^n(n-1)!\,
\frac{z_1^{n-1}z_2^{n-1}}{w_1^n w_2^2w_3^2},
\qquad n\ge2.
}
$$

Its explicit smooth angular representative is

$$
-z_1^{n-1}z_2^{n-1}
\partial_{w_1}^{n-1}\partial_{w_2}\partial_{w_3}U.
$$

This determines an exact all-order list of protected output
**channels**, including their factorial and sign. It is not the
value of the boundary joining cochain. At finite rank, trace
relations can further identify or remove the single-trace channel;
the stated one-dimensional result is the stable linear dictionary.

The distinction is essential because the ordinary product
$G_1G_2G_3H^{n-2}$ is not its complete quantum product.
For a fixed parenthesization let $\mathcal O_j$ be the coefficient
of $\hbar^j$ in the marked quantum product of these same
constituents. Their individual protection and the derivation
identity for the quantum product give the complete order-$n$
equation

$$
Q_n\mathcal O_0
+\sum_{j=1}^{n-1}Q_j\mathcal O_{n-j}
+Q_0\mathcal O_n=0.
$$

Thus a nonzero raw $Q_n\mathcal O_0$ neither produces a new
single-particle differential nor identifies a nonzero physical
bulk output by itself. Its effect is tied to the lower-order
operations and the composite representative. The explicit bulk
tree obstruction shows why that representative data cannot be
omitted from the holographic test.

## 6. Independent exact certificate

Run

```bash
python explore/twisted-holography-N4/boundary_bulk_all_orders/bulk/primary_tree_test.py
```

The [script](primary_tree_test.py) imports no repository evaluator
and uses only the Python standard library. Its
[saved output](primary_tree_test.json) records:

* the four complete Laurent fields and all input-sector dimensions;
* 32 independent checks of $\Phi\partial_\Omega=d\Phi$ on general
  polyvector polarizations;
* 420 normal-derivative checks through total derivative order six;
* 2,100 complete zero exterior products and their root divergences;
* the exact output charge, form dimension, Laurent sign, and
  factorial at each boundary order two through twenty-four;
* the exact counts of ordered cubic rooted trees through twelve
  input leaves.

The proof in Section 4 covers every finite tree and every
holomorphic derivative order; the bounded certificate verifies
its algebraic inputs independently.

## Sources and provenance

1. [Dictionary report, Sections A1 and A3](https://github.com/lbotao399-stack/holomorphic-topological-feynman-diagrams/blob/29ef4be3c92b8055ac51067357524c1048e6cda0/explore/twisted-holography-N4/dictionary/REPORT.txt):
   Laurent/form normalization, full linear differential, and the
   stable single-trace exact-form space.
2. [Earlier holography note, Sections 1–4](https://github.com/lbotao399-stack/holomorphic-topological-feynman-diagrams/blob/29ef4be3c92b8055ac51067357524c1048e6cda0/explore/twisted-holography-N4/NOTE.md):
   protected quadratic constituents and the distinction between
   cochains, states, and source/operator pairing.
3. [Earlier three-trace test](https://github.com/lbotao399-stack/holomorphic-topological-feynman-diagrams/blob/29ef4be3c92b8055ac51067357524c1048e6cda0/explore/twisted-holography-N4/README.md#the-three-trace-joining-symbol-and-the-four-point-bulk-tree):
   the previously proved two-loop cubic-tree vanishing and the
   explicit lower-order Ward cancellation.
4. K. Costello and S. Li, *Twisted supergravity and its quantization*,
   [arXiv:1606.00365](https://arxiv.org/abs/1606.00365), Section 9:
   punctured fivefold with the D3 flux background.
5. S. Li and M. Wang, *On the exact quantization of BCOV Theory on
   Calabi–Yau Manifolds*,
   [arXiv:2610.21927v1](https://arxiv.org/abs/2610.21927v1), Section 1:
   the announced coefficient-derivative/wedge form of the local
   interaction. Applying the tree statement to those positive-genus
   local vertices retains this prescription as an input.
