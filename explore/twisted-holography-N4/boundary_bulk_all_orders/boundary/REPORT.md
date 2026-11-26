# Every boundary loop order on protected multiparticle constituents

## Result

Within the all-loop fan construction in the repository's manuscript, every
$Q_n$ has a nonzero, independently isolated matrix element on an ordinary
product of single-graviton representatives. Each constituent is annihilated
by every $Q_k$, while its ordinary product requires quantum dressing.

The family for $n\ge2$ is

$$
\boxed{
\operatorname{Tr}(\beta_1\gamma^2)\,
\operatorname{Tr}(\beta_2\gamma^3)\,
\operatorname{Tr}(\beta_3\gamma^1)\,
\left[\operatorname{Tr}(b\gamma^1-ig\beta_2\beta_3)\right]^{n-2}.
}
\tag{1}
$$

Set $\beta=b=0$ and $\gamma^1=\gamma^2=\gamma^3=X$, with $X$ a fixed
traceless matrix. Use two independent odd amplitudes $\eta_1,\eta_2$ and
the ghost

$$
c(x)=
\eta_1D\!\left(e^{-p_{+\dot a}x^{+\dot a}/2}-1\right)
+\eta_2E\!\left(e^{-q_{+\dot a}x^{+\dot a}/2}-1\right).
\tag{2}
$$

The subtraction makes the undifferentiated ghost vanish at the insertion.
The matrices $D,E$ are traceless. The momenta $p,q$ are the
eigenvalues of the manuscript's spinor derivatives
$\partial_{+\dot a}=-2\partial/\partial x^{+\dot a}$ on their
respective exponentials. This fixes the factor of four in each
paired derivative. Retain the manuscript's fundamental trace normalization
$\operatorname{Tr}(t_A t_B)=\tau\delta_{AB}$ and
$[t_A,t_B]=iC_{AB}{}^C t_C$.

For the full formal boundary differential
$\mathbf Q=\sum_{k\ge0}\hbar^kQ_k$, the coefficient of
$\eta_1\eta_2$ in its action on (1), evaluated at the origin, is exactly

$$
\boxed{
-\frac{\hbar^n(ig\tau)^n}{(32\pi^2)^n(n!)^2}
\left(p_{+\dot a}q_+{}^{\dot a}\right)^n
\operatorname{Tr}\!\left((\operatorname{ad}_X)^{n+1}(D)\,E\right).
}
\tag{3}
$$

Here $\operatorname{ad}_X(Y)=[X,Y]$. **Every loop order other than $n$
has zero coefficient in this same measurement.** Thus (3) distinguishes
all the boundary $Q_n$ without identifying boundary loop order with
bulk genus.

For every matrix rank $N\ge2$, choose the upper-left block

$$
X=\operatorname{diag}\!\left(\frac12,-\frac12,0,\ldots,0\right),
\qquad D=E_{12},\qquad E=E_{21},
\qquad p_{+\dot a}q_+{}^{\dot a}=1.
\tag{4}
$$

The notation $E_{ij}$ denotes the matrix with its sole nonzero entry
equal to one in position $(i,j)$. Then

$$
(\operatorname{ad}_X)^{n+1}(D)=D,\qquad
\operatorname{Tr}(DE)=1.
$$

At $\tau=1/2$, the measured answer is

$$
\boxed{
-\frac{(ig\hbar)^n}{(64\pi^2)^n(n!)^2}\ne0
\qquad(n\ge2).
}
\tag{5}
$$

The separate protected two-factor test in Section 5 supplies $Q_1$.

## 1. Each factor is individually protected

The only inherited dynamical input is the manuscript's all-loop fan
formula and its single-graviton identity. In its conventions,

$$
C=igc+\theta_I\gamma^I
-\frac12\varepsilon^{IJK}\theta_I\theta_J\beta_K
+\frac{\theta_1\theta_2\theta_3}{ig}b.
\tag{6}
$$

The components of $\operatorname{Tr}(dC)^2$ obey

$$
Q_0G=0,\qquad Q_kG=0\quad(k\ge1).
\tag{7}
$$

The derivative $d$ in this inherited expression is the super de Rham
derivative in the two tangent and three odd coordinates, as defined in
the manuscript. The new family uses only its odd-coordinate part.
For example,

$$
\partial_{\theta_1}C
=\gamma^1-\theta_2\beta_3+\theta_3\beta_2
+\frac{\theta_2\theta_3}{ig}b.
$$

Let $[\theta_{i_1}\cdots\theta_{i_r}]$ mean extraction of the indicated
ordered odd monomial. The four required factors are

$$
\begin{aligned}
\operatorname{Tr}(\beta_1\gamma^2)
&=-\frac12[\theta_3]\operatorname{Tr}(\partial_{\theta_2}C)^2,\\
\operatorname{Tr}(\beta_2\gamma^3)
&=-\frac12[\theta_1]\operatorname{Tr}(\partial_{\theta_3}C)^2,\\
\operatorname{Tr}(\beta_3\gamma^1)
&=-\frac12[\theta_2]\operatorname{Tr}(\partial_{\theta_1}C)^2,\\
\operatorname{Tr}(b\gamma^1-ig\beta_2\beta_3)
&=\frac{ig}{2}[\theta_2\theta_3]
 \operatorname{Tr}(\partial_{\theta_1}C)^2.
\end{aligned}
\tag{8}
$$

For the last equality, graded multiplication gives

$$
[\theta_2\theta_3]\operatorname{Tr}(\partial_{\theta_1}C)^2
=\frac{2}{ig}\operatorname{Tr}(b\gamma^1)
-2\operatorname{Tr}(\beta_2\beta_3).
$$

All four factors are therefore protected individually. Since $Q_0$
is a derivation, their ordinary product (1) is also classically closed.
There is no appeal to a nonclosed constituent in the argument.

## 2. Exact isolation of one loop order

In a chosen monomial of (1), suppose that $j$ of its $n-2$ even
factors contribute the term $-ig\operatorname{Tr}(\beta_2\beta_3)$.
The monomial then contains

$$
3+2j\quad\text{beta letters},\qquad
n-2-j\quad\text{$b$ letters},\qquad
n+1-j\quad\text{scalar letters}.
\tag{9}
$$

At the evaluation $\beta=b=0$, every beta and $b$ letter must be
contracted. The functional derivative in the manuscript is

$$
\frac{\delta}{\delta C^A(\theta)}
=\frac{\theta_1\theta_2\theta_3}{ig}
 \frac{\partial_L}{\partial c^A}
-\frac12\varepsilon^{IJK}\theta_J\theta_K
 \frac{\partial}{\partial(\gamma^I)_A}
-\theta_I\frac{\partial_L}{\partial(\beta_I)^A}
-ig\frac{\partial}{\partial b_A}.
\tag{10}
$$

Each beta contraction supplies one odd coordinate, each scalar
contraction supplies two, and each $b$ contraction supplies none.
To extract the two-ghost output, both endpoint fields in the fan
must supply their $\theta^0$ components $igc$.

If $a$ additional scalar letters contract, the number of spokes and
the odd degree are

$$
k+1=(3+2j)+(n-2-j)+a=n+1+j+a,
$$

$$
\text{odd degree}=3+2j+2a=3+2(k-n).
\tag{11}
$$

The integral over three odd coordinates requires degree three.
Therefore $j=a=0$ and $k=n$. Fewer than $n+1$ spokes leave a beta
or $b$ letter uncontracted; additional spokes have excessive odd
degree. This proves the exact all-order isolation statement.

In particular, every surviving spoke differentiates exactly one
factor of (1), leaving its scalar partner. All $n+1$ scalar partners
equal $X$ at the measurement. Every assignment has the same
Grassmann sign. There are $(n+1)!$ labeled assignments; when the
last $n-2$ factors are written as a power, their ordinary derivative
factorials reproduce the same count.

The argument uses the complete fan action on (1). It does not retain
only the highest derivative term of an otherwise nonzero answer.

## 3. Complete coefficient, phase, and endpoint exchange

For an undifferentiated input, the inherited fan prefactor is

$$
\frac{(-1)^ni^{n+1}}
 {2(32\pi^2)^n(n!)^3(n+1)}.
\tag{12}
$$

The three beta derivatives give $-\theta_1\theta_2\theta_3$ in every
assignment. The $b$ contractions and two endpoint ghosts supply,
respectively,

$$
(-ig)^{n-2},\qquad (ig)^2.
\tag{13}
$$

For each differentiated quadratic trace the remaining scalar gives
$\operatorname{Tr}(t_A X)=\tau X^A$. The full structure-constant
chain is consequently

$$
\begin{aligned}
&\tau^{n+1}X^{A_1}\cdots X^{A_{n+1}}\,
C_{A_1DR_1}C_{A_2R_1R_2}\cdots C_{A_{n+1}R_nE}\\
&\qquad=\tau^n i^{-(n+1)}
\operatorname{Tr}\!\left((\operatorname{ad}_X)^{n+1}(t_D)t_E\right).
\end{aligned}
\tag{14}
$$

Combining (12)–(14) and the $(n+1)!$ assignments gives the coefficient
before substituting the two odd ghost waves:

$$
\begin{aligned}
&\frac{(-1)^ni^{n+1}}
 {2(32\pi^2)^n(n!)^3(n+1)}
\left[-(-ig)^{n-2}(ig)^2\right]
(n+1)!\tau^n i^{-(n+1)}\\
&\qquad=
-\frac{(ig\tau)^n}{2(32\pi^2)^n(n!)^2}.
\end{aligned}
\tag{15}
$$

Swapping the endpoints changes the adjoint-chain trace by
$(-1)^{n+1}$, the contracted spinor factor by $(-1)^n$, and the
ordered odd amplitudes by $-1$. Their product is one. The two
endpoint assignments therefore add, canceling the denominator two
in (15) and proving (3).

For a diagonal $X$ whose first two eigenvalues have difference
$\Delta$, the complete result specializes to

$$
-\frac{\hbar^n(ig\tau)^n}{(32\pi^2)^n(n!)^2}
\Delta^{n+1}\left(p_{+\dot a}q_+{}^{\dot a}\right)^n.
\tag{16}
$$

For reference, with $\Delta=1$, unit spinor contraction and
$\tau=1/2$, the first five nonzero coefficients of the corresponding
families are

$$
\begin{array}{c|c}
n&\text{coefficient of }\hbar^n\eta_1\eta_2\\ \hline
2&g^2/(16384\pi^4)\\
3&ig^3/(9437184\pi^6)\\
4&-g^4/(9663676416\pi^8)\\
5&-ig^5/(15461882265600\pi^{10})\\
6&g^6/(35624176739942400\pi^{12})
\end{array}
\tag{17}
$$

The genus-two value agrees with the independent Pauli-matrix
three-trace calculation already in the repository. The formulas
for arbitrary $n$, their protected factors, and the exact isolation
of $Q_n$ are the new extension.

## 4. Differential order and product-compatible changes

Multiplication by a homogeneous invariant operator $G$ is denoted
by $M_G$. The graded commutator with multiplication is
$[Q,M_G]=QM_G-(-1)^{|Q||G|}M_GQ$.
A differential operator has order at most $r$ if every
$(r+1)$-fold such commutator vanishes.

The fan formula contains $n+1$ functional derivatives, so $Q_n$
has order at most $n+1$ on the invariant local-operator algebra.
Apply $n+1$ multiplication commutators using the factors of (1).
At the measurement locus, every factor itself vanishes. All
subtraction terms that retain an undifferentiated factor vanish,
and the remaining term is exactly (3). Hence

$$
\boxed{\operatorname{ord}(Q_n)=n+1\qquad(n\ge2).}
\tag{18}
$$

Section 5 gives order two for $Q_1$.

This order is preserved by a formal change of variables that is an
automorphism of ordinary multiplication and is the identity at
$\hbar=0$. Its formal logarithm is a series of derivations.
Commuting an operator with a derivation never raises differential
order. At coefficient $\hbar^n$, every conjugation term involving
$Q_k$ with $k<n$ therefore has order at most $k+1\le n$.
It cannot cancel the nonzero order-$n+1$ symbol of $Q_n$.
An invertible constant change of variables conjugates that symbol
and also cannot annihilate it.

Thus no product-compatible formal coordinate change removes any
member of this tower, or turns the full boundary differential into
an operator of bounded differential order. This does not contradict
the previously established filtered *linear* conjugacy on stable
single traces: a general linear change of basis need not preserve
products. Nor does it prohibit a description with a deformed
quantum product or higher multiparticle operations.

## 5. The one-loop member and the completed products

The manuscript already supplies the protected pair

$$
\operatorname{Tr}((\gamma^1)^2),\qquad
\operatorname{Tr}(\beta_1\gamma^2).
\tag{19}
$$

At $\beta=b=0$, $\gamma^1=\gamma^2=X$, with the ghost (2), its
ordinary product has the exact coefficient

$$
\boxed{
\frac{2\hbar\tau g^2}{32\pi^2}
\left(p_{+\dot a}q_+{}^{\dot a}\right)
\operatorname{Tr}\!\left((\operatorname{ad}_X)^2(D)E\right).
}
\tag{20}
$$

The complete classical differential vanishes. Every $Q_{k\ge2}$
vanishes because the input has odd-coordinate degree five and
$Q_k$ lowers that degree by $3k$. At the unit test (4) and
$\tau=1/2$, (20) is $\hbar g^2/(32\pi^2)$.
The coupling is $g^2$ in this component sector, while (1) gives
$g^n$; these factors come from the physical component normalization
in (6), not from identifying $g$ with the formal loop parameter.

For every $n$, replacing the ordinary multiplication in (1) by
the manuscript's marked quantum products gives a perturbatively
closed operator. The equation at order $\hbar^n$ is the complete
sum

$$
Q_0\mathcal O_n+Q_1\mathcal O_{n-1}
+Q_2\mathcal O_{n-2}+\cdots+Q_n\mathcal O_0=0,
\tag{21}
$$

where $\mathcal O_j$ denotes the coefficient of the chosen iterated
marked product. Formula (3) fixes one nonzero exact coefficient
that this sum must cancel. It does not compute the individual
completion coefficients $\mathcal O_j$, nor imply that the
classical multiparticle class is lifted. The isolated matrix
element is a differential on the ordinary representative, not a
nonzero differential between protected cohomology classes.

### Nonexactness of the completed family

For this family one can also prove nonexactness; quantum closure
alone would not establish it. Use a second evaluation, distinct
from the off-diagonal ghost experiment:
all holomorphic jets and $c$ vanish, and

$$
\gamma^1=\gamma^2=\gamma^3=X,\qquad
\beta_I=\chi_I X,\qquad b=vX,
$$

where the $\chi_I$ are independent odd scalars and $v$ is an even
scalar. The matrix $X$ is the diagonal block in (4).
Every classical color commutator vanishes, and every positive-loop
fan has differentiated output fields. Thus this evaluation
annihilates $\mathbf Q$ on every input cochain.

The ordinary product (1) has the exact value

$$
\begin{aligned}
&\left(\operatorname{Tr}X^2\right)^{n+1}
\chi_1\chi_2\chi_3\left(v-ig\chi_2\chi_3\right)^{n-2}\\
&\qquad=\frac{\chi_1\chi_2\chi_3v^{n-2}}{2^{n+1}}\ne0.
\end{aligned}
\tag{22}
$$

A marked completion has this nonzero constant coefficient in
$\hbar$. Its full evaluation is therefore nonzero, whereas every
quantum boundary evaluates to zero. The completion is nonexact
both over formal power series in $\hbar$ and after allowing
finitely many negative powers of $\hbar$.
For the pair (19), the corresponding value is $\chi_1/4\ne0$.
These arguments work at every finite rank $N\ge2$ by using the
same embedded two-by-two block.

## 6. The bulk target and its limits

Every nonzero term in (3) touches all $n+1$ input traces and leaves
one output trace. Its ribbon graph has

$$
\text{vertices}=2n+2,\qquad
\text{edges}=2n+1,
$$

so it is a tree. It has no empty color face and no handle. In the
general trace-count formula derived in the
[matching report](../matching/REPORT.md), this is
$r=n+1$, $b=1$, $h=0$, and $F=0$.
With unnormalized fundamental traces and the formal parameter
$\hbar N$ fixed, its coefficient scales exactly as $N^{-n}$.
The relevant closed-string target is therefore a genus-zero
$(n+2)$-leg operation with $n+1$ incoming particles and one
outgoing particle.

This is a topological and algebraic specification of the bulk test.
It is not a calculation of a punctured BCOV amplitude. The
dictionary must carry the marked product and the full relative
bulk/brane operations, not only the linear single-particle state
map. In particular one cannot infer that the genus-$n$ elementary
flat BCOV Hamiltonian is the image of this boundary $Q_n$.
The [bulk report](../bulk/REPORT.md) tests the ordinary primary
bulk-tree interpretation of these same constituent charge sectors.

## 7. Independent exact certificate

Run from the repository root:

```bash
python explore/twisted-holography-N4/boundary_bulk_all_orders/boundary/all_order_join.py
```

The script imports no inherited evaluator. Its saved output is
[all_order_join.json](all_order_join.json). It verifies:

* all four single-graviton component formulas by graded exterior
  multiplication and graded cyclic traces;
* 10,920 ordered Berezin-derivative arrangements at $2\le n\le14$,
  retaining the factorials of the identical even factors;
* 24 direct structure-constant chains against nested Pauli
  commutators;
* 540 further exact root-pair checks at matrix ranks two through
  seven, with positive and negative root gaps;
* loop isolation for every possible monomial/derivative count in
  the 39 families $2\le n\le40$;
* the full fan phase, trace normalization, endpoint exchange and
  the separate one-loop Berezin calculation and independent two-loop
  Pauli normalization;
* the exact nonzero constant-Cartan values of the completed-family
  classical terms.

The proofs in Sections 2–4 hold at arbitrary loop order. These
finite checks independently verify the signs and integer formulas.

### Provenance

Inherited inputs are pinned to repository revision
[`29ef4be3c92b8055ac51067357524c1048e6cda0`](https://github.com/lbotao399-stack/holomorphic-topological-feynman-diagrams/tree/29ef4be3c92b8055ac51067357524c1048e6cda0):

1. [`sources/LLC_main.tex`](https://github.com/lbotao399-stack/holomorphic-topological-feynman-diagrams/blob/29ef4be3c92b8055ac51067357524c1048e6cda0/sources/LLC_main.tex):
   equations `eq-Qnzero`, `eq-A-13`, `eq-A-16` and the
   single-graviton theorem `b-00630`.
   Git blob `db6b3994c61a469268e57b71bf426cc1281b4a10`.
2. [Existing holography note](https://github.com/lbotao399-stack/holomorphic-topological-feynman-diagrams/blob/29ef4be3c92b8055ac51067357524c1048e6cda0/explore/twisted-holography-N4/NOTE.md):
   stable single-trace linear conjugacy and the earlier $Q_2$
   three-trace test.
   Git blob `439f5c157fb587742c954131f4a6beaa0e784de6`.

The new conclusions concern the boundary operator construction and
the exact targets it imposes. No complete boundary/bulk quantum
equivalence is assumed.
