# Exact Bulk Tests in Twisted Supergravity

## Locality of the further twist, delayed quantum responses, and Ward identities

Research development, 26 November 2026.

### Abstract

We develop three exact tests of the twisted supergravity calculations in
this repository. In eleven dimensions, the complete linear further-twist
complex becomes a three-variable Koszul complex in normal momentum.
Its cohomology, singular fibers, and dependence on the cubic coupling
identify the obstruction to a uniform translation-invariant local
elimination of the transverse modes. In the flat type-IIB model, an
admissible polarization of three scalar descendants cancels every donor
fusion and has its first nonzero elementary scalar response at genus
four. Finally, two complete three-wave compositions eliminate unknown
four-, five-, and six-loop binary coefficients from Ward tests. These
calculations give exact inputs for the nonlinear anomaly problem and
for bulk-to-boundary transfer. They do not perform that transfer.

The type-IIB calculations use the flat vertex prescription announced by
Li–Wang and the complete period tables already in this repository.
The new work consists of algebraic proofs, exact recombinations,
and independent certificates; it does not claim a new evaluation of
those inherited worldsheet integrals.

## 1. Questions and results

The existing calculations leave two physically consequential questions.
Does the classical further twist supply a local primitive for the
eleven-dimensional two-loop defect? Which flat bulk quantum terms can
enter a consistent holographic symmetry operation?

The results below settle concrete preliminary questions without assuming
the desired anomaly cancellation or holographic identification.

| Question | Exact result | Detailed proof |
|---|---|---|
| Which transverse modes survive the further-twist differential? | A full Koszul complex, singular momentum fibers, and coupling torsion | [Further-twist report](../../explore/sugra-11d-qme/developments_2026_11/REPORT.md) |
| Can allowed scalar descendants delay the first bulk response? | A finite polarization with zero response at genera zero through three and a nonzero genus-four response | [Delayed-response report](../../explore/twisted-holography-N4/developments_2026_11/REPORT.md) |
| Can unknown higher-genus periods repair the isolated binary Jacobi identity? | Explicit combinations eliminate those periods and retain nonzero lower-genus products | [Ward-identity report](../../explore/iib-bcov/developments_2026_11/REPORT.md) |
| Is the compact-momentum support bound sharp? | An exact routing criterion and momentum-conserving saturation in every compact dimension | [Compact-momentum calculation](compact_momenta.md) |

The starting repository revision is
[`aa0d2f1725dff135494eef132f5e92a8c0401a42`](https://github.com/lbotao399-stack/holomorphic-topological-feynman-diagrams/tree/aa0d2f1725dff135494eef132f5e92a8c0401a42).
Each detailed report identifies its inherited inputs separately from
its new derivations.

## 2. Eleven dimensions: the full further-twist symbol

The published further twist uses the constant two-form
$dz_1\wedge dz_2$. In the strict field-strength model, retain a vector
polarization $X$ and a two-form polarization $F$ at holomorphic
momentum $p$. Their full constraints are

$$
p(X)=0,\qquad p\wedge F=0.
$$

Let $\Omega=dz_1\wedge dz_2\wedge dz_3\wedge dz_4\wedge dz_5$ be the
holomorphic volume. For a four-form $\alpha$, the vector
$\alpha^\sharp$ is defined by $\iota_{\alpha^\sharp}\Omega=\alpha$.
The further-twist differential, denoted by $\delta$, is

$$
\delta X=-p\wedge\iota_X(dz_1\wedge dz_2),\qquad
\delta F=g(dz_1\wedge dz_2\wedge F)^\sharp,
$$

where $g$ is the cubic coupling. Even and odd below refer to the
exceptional super Lie algebra before the physical field shift.

For nonzero normal momentum and $g\ne0$, both differentials have
rank two on four-dimensional spaces; their kernel/image equalities
give zero cohomology. At nonzero planar momentum, both maps vanish.
At zero momentum, the two-form map has rank three. Thus the complete
cohomology dimensions are

$$
\begin{array}{c|cc}
\text{momentum stratum}&\text{even}&\text{odd}\\ \hline
p\vert_{\mathbb C^3}\ne0&0&0\\
p\vert_{\mathbb C^3}=0,\ p\ne0&4&4\\
p=0&2&7
\end{array}
$$

Here $p\vert_{\mathbb C^3}$ is the covector evaluated on the three
transverse directions.

### 2.1 Polynomial cohomology and the formal coupling

On a chart with nonzero planar momentum, let $q_3,q_4,q_5$ be the
three normal momentum components and use the frame
$p=dz_1+q_3dz_3+q_4dz_4+q_5dz_5$. The vector $q$ in the following
maps consists of these three components; the wedge and contraction
maps use the transverse volume. For the polynomial ring
$R=\mathbb C[g,q_3,q_4,q_5]$, the complete free complex is

$$
0\longrightarrow R
\xrightarrow{\ q\ }R^3
\xrightarrow{\ g\,q\wedge\ }R^3
\xrightarrow{\ q\cdot\ }R
\longrightarrow0.
$$

The successive parities are even, odd, even, odd.
With $H$ denoting kernel modulo image, exactness of the unscaled
Koszul sequence gives

$$
\boxed{
H^{\mathrm{even}}
=\frac{\ker(q\cdot)}{g\,\ker(q\cdot)},\qquad
H^{\mathrm{odd}}=R/(q_3,q_4,q_5).
}
$$

The even module is annihilated by $g$. Inverting $g$ removes it.
Setting $g=1$ leaves a single odd polynomial cohomology module
supported at zero normal momentum.

This one module does not mean that the zero-momentum fiber has
only one mode. Taking that fiber in the free complex annihilates
all arrows and retains exterior dimensions $1,3,3,1$.
These are the additional homology groups produced by imposing the
three normal-momentum equations on a free resolution. They give
four even and four odd modes, in agreement with the direct calculation.

### 2.2 The locality test

Restrict to $q_3=t,\ q_4=q_5=0$, where $t$ is a single normal momentum.
The exact cohomology modules become

$$
H^{\mathrm{even}}=(\mathbb C[g,t]/(gt))^2,\qquad
H^{\mathrm{odd}}=(\mathbb C[g,t]/t)^2.
$$

At fixed nonzero $g$ and $t$, every pair is contractible.
Its inverse uses $t^{-1}$, and two of the pairs also require
$g^{-1}$. If a matrix $h$ of constant-coefficient finite-order
differential operators contracted the full normal complex uniformly,
its momentum entries would be polynomial and it would satisfy

$$
\delta h+h\delta=\mathrm{id}.
$$

At $t=0$, the left side is zero and the right side is the identity.
Therefore such a uniform contraction does not exist.

The established anomaly pairing is $-63g/(256\pi^{10})$, and its
restriction to the residual Hamiltonian inputs vanishes. The new
calculation explains why that restriction alone supplies no
translation-invariant finite-jet primitive. A primitive for the
particular nonlinear density could still exist. It must be tested in
the full local complex, including total derivatives and the formal
coupling. The published classical equivalence to Poisson
Chern–Simons theory on $\mathbb C^2\times\mathbb R^7$ remains compatible
with this distinction.

## 3. Type IIB: an admissible response beginning at genus four

Use the minimal fivefold scalar descendants $A_k$ with
$0\le k\le4$. Three designated scalar insertions carry planar momenta
$e_1,e_2,-e_1-e_2$, where $e_1,e_2$ are the unit planar covectors.
Their total descendant level is four.

The following coefficients apply to each distinct permutation of the
listed levels:

$$
\begin{array}{c|rrrr}
\text{levels}&(4,0,0)&(3,1,0)&(2,2,0)&(2,1,1)\\ \hline
\text{coefficient}&2&-1&1&0\\
\text{number of permutations}&3&6&3&3
\end{array}
$$

More explicitly, the full three hard momenta are

$$
e_1-4e_3-4e_4-4e_5,\qquad e_2,\qquad -e_1-e_2.
$$

Add the primary scalar $e^{z_3}+e^{z_4}+e^{z_5}$ and the primary
vector $e^{z_4}\partial_{z_3}$. Take the part linear in each of the
three separately labeled hard waves and in that vector, and sum
its descendant assignments with the table's weights.
The report supplies the equivalent explicit differential operator
on independent source amplitudes. The construction uses only the
allowed descendant fields. The statement below concerns this
selected multilinear coefficient; other coefficients of the source
are not asserted to vanish.

### 3.1 Cancellation of every donor fusion

Write $C_{k_1,k_2,k_3}$ for precisely the coefficients in the table,
extended symmetrically. For every triple of nonnegative levels
$\ell_1+\ell_2+\ell_3=3$,

$$
C_{\ell_1+1,\ell_2,\ell_3}
+C_{\ell_1,\ell_2+1,\ell_3}
+C_{\ell_1,\ell_2,\ell_3+1}=0.
$$

A donor is a primary insertion absorbed into one of the three
descendant insertions in the string-equation expansion. Repeatedly
applying the displayed identity cancels every nonempty set of such
donors, including the labeled multinomial multiplicities.
The classical multinomial coefficient also cancels:

$$
3\cdot2\cdot\frac{4!}{4!}
-6\cdot\frac{4!}{3!\,1!}
+3\cdot\frac{4!}{2!\,2!}
=6-24+18=0.
$$

After cancellation, genus $g$ has $2g+4$ vertical marked points.
The three-dimensional transverse Abel map has complex rank at most
$3g$. For $g=1,2,3$,

$$
2g+4-3g=4-g>0.
$$

The required vertical top-degree form therefore vanishes. The
classical cancellation above supplies genus zero separately.

### 3.2 The first nonzero coefficient

The complete inherited planar genus-four periods, with the three
unit-area momenta and the displayed descendant levels, are

$$
\begin{aligned}
V_{4;(4,0,0)}&=\frac{731}{1024192512000},&
V_{4;(3,1,0)}&=\frac{1493}{438939648000},\\
V_{4;(2,2,0)}&=\frac{157}{29262643200},&
V_{4;(2,1,1)}&=\frac{433}{37623398400}.
\end{aligned}
$$

Here $V$ denotes the scalar worldsheet vertex in the existing
normalization. The new polarization gives

$$
6V_{4;(4,0,0)}
-6V_{4;(3,1,0)}
+3V_{4;(2,2,0)}
=-\frac{31}{1024192512000}.
$$

The transverse fiber contributes $(4!)^3=13824$. The contributing
fixed-source scalar multiplicities are $(3,3,4)$, giving the divisor
$3!3!4!=864$. The root-vector pairing is one.
The elementary-field convention contributes the sign
$-(-1)^4=-1$. Consequently, for $Q_g\phi(0)$ denoting the selected
genus-$g$ scalar response at the origin,

$$
\boxed{
Q_0\phi(0)=Q_1\phi(0)=Q_2\phi(0)=Q_3\phi(0)=0,\qquad
(2\pi i)^{20}Q_4\phi(0)
=\frac{31}{64012032000}.
}
$$

The response is a flat local bulk coefficient in the fixed field
coordinates. Identifying a corresponding boundary observable still
requires the punctured D3-brane background, the complete differential,
and the bulk-to-boundary map.

The full finite linear-algebra classification is also computed:
symmetric three-insertion polarizations obeying the same donor
identity and $0\le k_j\le4$ exist only at total levels
$0,2,3,4,6$, with one dimension at each level.
The level-six polarization kills lower genera, but its first possible
genus-six period has not been evaluated here.

## 4. Ward identities without unknown higher-loop periods

In a fixed normal cubic background $B$, the inherited positive-genus
primary tangent operations are binary. Let $q_g$ denote their
coefficients after multiplication by $(-1)^g(2\pi i)^{5g}$, and use
$\hbar$ as the formal genus variable. The multiplier of the planar
determinant $s$ is

$$
F(s)=\sum_{g\ge0}q_g\hbar^g s^{2g+1},
\qquad q_0=1,\quad q_1=0,\quad q_2=\frac{s(B)}{51840}.
$$

The symbol $s(B)$ here is the existing degree-four invariant of the
normal cubic; it is distinct from the determinant argument of $F$.
Its normalization and the known coefficient $q_3$ are given in the
linked report.

Use the complete cyclic Jacobiator of these binaries. Denote its
genus-$g$ coefficients by $J_g^{(1)}$ and $J_g^{(2)}$ on the two
three-wave inputs

$$
\begin{aligned}
(p,q,r)&=((1,0),(0,1),(1,2)),\\
(p,q,r)&=((1,0),(0,1),(2,3)).
\end{aligned}
$$

At genus four, every ordered genus split is included:

$$
J_4^{(1)}=17640q_4-780q_2^2,\qquad
J_4^{(2)}=408240q_4-26280q_2^2.
$$

The following combination cancels the unknown $q_4$ exactly:

$$
\boxed{
7J_4^{(2)}-162J_4^{(1)}
=-57600q_2^2
=-\frac{s(B)^2}{46656}.
}
$$

Similarly,

$$
\boxed{
32J_5^{(2)}-1389J_5^{(1)}=-5040000q_2q_3,
}
$$

and the full genus-six identity is

$$
\boxed{
13J_6^{(2)}-1038J_6^{(1)}
=4838400q_2q_4-20180160q_3^2.
}
$$

On $s(B)=0$, the genus-six right side is
$-20180160q_3^2$, independent of the unknown $q_6$ and $q_4$.
Other contributions in the complete Ward identity must cancel the
displayed targets if its total vanishes. Choosing the higher-loop
coefficient of the isolated binary cannot do so.

### 4.1 The all-order algebraic explanation

For three planar momenta, write $a,b,c$ for their cyclic pairwise
determinants. Strict Jacobi for an alternating determinant-only
binary is equivalent to

$$
F(a)F(b-c)+F(b)F(c-a)+F(c)F(a-b)=0.
$$

The normalized odd series satisfies $F'(0)=1$. Differentiating in
$c$ at zero gives

$$
F(a-b)=F(a)F'(b)-F(b)F'(a).
$$

Two further derivatives in $b$ at zero yield

$$
F''(a)=F'''(0)F(a).
$$

The exact formal solution is therefore

$$
F(a)=\sum_{g\ge0}
\frac{(6q_1\hbar)^g}{(2g+1)!}\,a^{2g+1},
\qquad
q_g=\frac{(6q_1)^g}{(2g+1)!}.
$$

These are the standard Moyal-type coefficients. Since the inherited
flat bulk binary has $q_1=0$, strict Jacobi in this isolated class
would force every $q_g$ with $g>0$ to vanish.
The new higher-loop combinations quantify the required completion.
They do not establish an anomaly of the full type-IIB quantum theory
or identify the isolated binary with a boundary charge bracket.

## 5. Compact directions: exact support and routing

For the eleven-dimensional theta graph, the three paths carry
momentum columns $(p,0)$, $(p,p)$, and $(0,p)$.
When momenta are evaluated on $r$ compact directions, let
$\operatorname{span}A,\operatorname{span}B,\operatorname{span}C$
be the spans of the compact covectors on the corresponding paths.
The compact rank equals

$$
\dim\operatorname{span}A+\dim\operatorname{span}B
+\dim\operatorname{span}C
-\dim(\operatorname{span}A\cap\operatorname{span}B
\cap\operatorname{span}C).
$$

This follows from the kernel of
$(a,b,c)\mapsto(a+b,b+c)$, which consists of $(-b,b,-b)$ in the
triple intersection.

The noncompact rank is at most $2(5-r)$, so a nonzero full kernel
requires at least $2r$ slots with compact momentum.
This bound is attained simultaneously for every $1\le r\le5$:
on shape $(1,6,6)$, place the five unit covectors on the second path
and their negatives on the third. Momentum conservation holds, the
full determinant is $-1$, and exactly $2r$ slots carry compact momentum.
In the inherited orientation the geometric kernel is

$$
K_{166}=-\frac{\pi}{158400}.
$$

This is sharp support for the geometric theta kernel. The complete
supergravity vertex also includes tensor contractions and all graph
shapes; their sum is not calculated for this new fixture.

## 6. Reproduction and remaining calculations

From the repository root, with Python 3.10 or later:

```bash
python -m pip install -r notes/twisted-sugra-bulk-tests/requirements.txt
python notes/twisted-sugra-bulk-tests/run_checks.py
```

The five certificates use exact integers, rational numbers, symbolic
polynomials, and exact zeta constants. The replay compares saved
outputs byte for byte and records source/output hashes in
[`verification/replay.json`](verification/replay.json).

The checks include the complete fixed-momentum maps, explicit cycle
primitives, all allowed three-scalar descendant levels, labeled donor
assignments, Fourier multiplicities, full nested binary compositions,
an independent Leibniz-rule evaluator, and compact routing ranks.
The all-degree conclusions are proved in the reports; finite
certificates check their formulas and normalization independently.

Three specific calculations remain necessary for a holographic
interpretation:

1. Solve the further-twisted nonlinear local primitive equation for
   the actual eleven-dimensional two-loop density, retaining the
   normal complex and the formal cubic coupling.
2. Transport the explicit genus-four scalar coefficient through the
   punctured D3-brane background and identify its boundary image or
   demonstrate its cancellation there.
3. Evaluate the complete Ward contributions complementary to the
   isolated binary on the two displayed momentum inputs. Their
   combinations must cancel the exact targets above.

## References

1. S. Raghavendran, I. Saberi, B. R. Williams,
   *Twisted eleven-dimensional supergravity*,
   [arXiv:2111.03049](https://arxiv.org/abs/2111.03049),
   especially Sections 4.2–4.3.
2. F. Hahner, N. M. Paquette, S. Raghavendran,
   *Higher Symmetries in Twisted Eleven-Dimensional Supergravity*,
   [arXiv:2506.02246v3](https://arxiv.org/abs/2506.02246v3),
   especially Section 4.3.2.
3. S. Li, M. Wang,
   *On the exact quantization of BCOV Theory on Calabi–Yau Manifolds*,
   [arXiv:2610.21927v1](https://arxiv.org/abs/2610.21927v1),
   especially the flat vertex prescription in Section 1.
   This source announces the construction; the results using its
   prescription retain that assumption.
4. Exact repository inputs and their pinned revisions are listed
   in the three sector reports linked in Section 1.
