# Integrate the series fibres before reconstructing the core

Status: **[drv]** algebraic derivations below; **[src]** the existing current
and density conventions and the Grundmann–Möller simplex formula. Exact
implementation checks live in `code/tests/test_dirichlet.py` and the tests
of `code/htk/series.py`. A check at finitely many points is not used as a
proof of an identity in this note.

The contribution is a reordering of the exact computation. Bivalent
subdivision vertices become constant columns of the current determinant.
After eliminating those columns, integrate a polynomial over the series
fractions. Only the remaining core resistances need reconstruction or
period reduction. All original vertex sources and edge currents remain in
the calculation.

## 1. Conventions and scope

Use [`code/htk/CONVENTIONS.md`](../../../code/htk/CONVENTIONS.md), Section 2,
and [`code/htk/brute.py`](../../../code/htk/brute.py). The graph is connected,
loopless and oriented; parallel edges are allowed. Let its edge count be
$N$, loop number $n\ge1$, and holomorphic dimension $H\ge1$, with

$$N-1=nH,\qquad T=0.$$

The conserved sources are the original $H$-component vertex vectors
$\lambda_v$. The positive cycle matrix is $C$, the edge currents are $y_e$,
and

$$
M=C\operatorname{diag}(t)C^{\mathsf T},\qquad U=\det M,\qquad
\Psi_{(i,r),e}=C_{re}y_e^i.
$$

Kirchhoff's voltage law and the density formula are

$$
\Psi t=0,\qquad
D(t)=\frac{(-1)^{e_0}\det\Psi_{\widehat{e_0}}}{t_{e_0}U^H}.
$$

Edge and row indices in this formula start at zero. The hat deletes a
column; the row order is $(i,r)$, with $i$ outermost. The choice of $e_0$
does not change the result. No restriction is placed on the sources at
bivalent vertices.

## 2. An exact column operation for every chain

Partition the original edge IDs into maximal chains. Write $m$ for their
number and $L_c$ for the length of chain $c$. Fix a traversal

$$p_c\longrightarrow u_{c,1}\longrightarrow\cdots
\longrightarrow u_{c,L_c-1}\longrightarrow q_c.$$

For the $j$th edge, $s_{c,j}$ is $+1$ when its stored orientation agrees
with traversal and $-1$ otherwise. The signed cycle columns are constant
along the chain:

$$s_{c,j}C_{\cdot,e_{c,j}}=C_c.$$

Introduce the chain resistance and its barycentric fractions immediately
before using them:

$$
R_c=\sum_jt_{e_{c,j}},\qquad
x_{c,j}=\frac{t_{e_{c,j}}}{R_c},\qquad
\sum_jx_{c,j}=1.
$$

The traversal current on edge $j$ is $s_{c,j}y_{e_{c,j}}$. Current
conservation gives

$$
s_{c,j+1}y_{e_{c,j+1}}-s_{c,j}y_{e_{c,j}}
=\lambda_{u_{c,j}}.
$$

Consequently the corresponding determinant columns, denoted here by
$\Psi_{c,j}$, obey

$$
\Psi_{c,j+1}-\Psi_{c,j}=C_c\otimes\lambda_{u_{c,j}}.
$$

Replace the $L_c$ columns of this block by

$$
\left(
\sum_jx_{c,j}\Psi_{c,j},\quad
\Psi_{c,2}-\Psi_{c,1},\quad\ldots,\quad
\Psi_{c,L_c}-\Psi_{c,L_c-1}
\right).
$$

The column-change matrix has determinant $\sum_jx_{c,j}=1$. Define
$Y_c=\sum_jx_{c,j}s_{c,j}y_{e_{c,j}}$. The transformed block is exactly

$$
\boxed{\left(C_c\otimes Y_c,\quad
C_c\otimes\lambda_{u_{c,1}},\quad\ldots,\quad
C_c\otimes\lambda_{u_{c,L_c-1}}\right).}
$$

Thus $N-m$ columns contain no Schwinger variables. Dependence among these
fixed columns proves $D\equiv0$ for the specified sources. Otherwise,
ordinary exact row elimination makes the fixed block triangular and
records every row exchange and pivot. The remaining determinant has size
$m-1$; its multiplier is the recorded signed product of pivots. This is a
determinant identity, with no discarded source or current component.

### Orientation

Let $\pi$ be the permutation that groups the original columns in chain
traversal order. Let $p_g$ be the position of the averaged column of a
chosen chain $g$ in that grouped order. The transformed null vector is

$$
(R_1,0,\ldots,0;\ R_2,0,\ldots,0;\ \ldots;\ R_m,0,\ldots,0),
$$

because the old chain time vector is $R_c(x_{c,1},\ldots,x_{c,L_c})$.
Deleting the averaged column of chain $g$ therefore divides its cofactor
by $R_g$. If $\tau$ moves the fixed columns of this minor to the front,
then, writing $F$ for that fixed block and $A_{\ne g}$ for the retained
averaged columns,

$$
\boxed{
D(t)=\frac{\operatorname{sgn}(\pi)(-1)^{p_g}
\operatorname{sgn}(\tau)}{R_gU^H}
\det[\,F\mid A_{\ne g}\,].
}
$$

This fixes the sign before row elimination. Different $g$ give the same
scalar density. For a pure polygon, cut its single closed chain at one
vertex: $m=1$, and the residual determinant is the empty determinant $1$.
The implementation also retains closed chains attached to a larger
graph; their averaged current is zero by the voltage law.

## 3. Affine current transport and the fibre degree

The source-transport formula already appears in
[`assembly/REPORT.md`, Section 1](../assembly/REPORT.md). With

$$\sigma_{c,k}=\sum_{j>k}x_{c,j},$$

it reads

$$
\begin{aligned}
s_{c,j}y_{e_{c,j}}
&=Y_c-\sum_k\sigma_{c,k}\lambda_{u_{c,k}}
  +\sum_{k<j}\lambda_{u_{c,k}},\\
\lambda_w^{s}
&=\lambda_w+\sum_{c,k}\lambda_{u_{c,k}}
\left([p_c=w]\sigma_{c,k}+[q_c=w](1-\sigma_{c,k})\right).
\end{aligned}
$$

All original sources occur. At fixed $R$, solving the core circuit is
linear in $\lambda^s$; hence every $Y_c$ and every original $y_e$ is affine
in all chain fractions together. The same fact follows without creating
a separate core graph. Choose any tree current $y_0$ with the original
divergence. Then

$$
y=y_0-C^{\mathsf T}M(R)^{-1}C\operatorname{diag}(Rx)y_0,
\qquad M(R)=\sum_c R_c C_cC_c^{\mathsf T}.
$$

The inverse depends only on $R$ and can be reused at every fibre point.
There are $m-1$ variable columns after the fixed block is eliminated.
For any monomial consisting of $k$ original edge-current factors, counting
repeated factors with multiplicity,

$$
\boxed{\deg_x\left(D(Rx)\prod_{a=1}^{k}y_{e_a}^{i_a}(Rx)\right)
\le m-1+k.}
$$

This is a bound on total degree, not a separate bound in every fraction.
It is valid with arbitrary internal sources. Source specializations may
lower it; no such lowering is required by the algorithm.

The absolute Jacobian of the series coordinates is

$$
\prod_e dt_e
=\left(\prod_cR_c^{L_c-1}\right)
\prod_c dR_c\,d^{L_c-1}x_c.
$$

Therefore the fibre coefficient to integrate is

$$
\boxed{
\int_{\prod_c\Delta^{L_c-1}}
\left(\prod_cR_c^{L_c-1}\right)
D(Rx)\prod_{a=1}^{k}y_{e_a}^{i_a}(Rx)\,dx.
}
$$

The factors $R_c^{L_c-1}$ do not change its fraction degree. Fixing one
positive $R_g=1$ gives the remaining projective core chart. Its dimension
is $m-1$, whereas the fibre dimension is $N-m$.

## 4. Exact rational integration on the product of simplices

For a simplex of dimension $b$, let $\beta$ run over nonnegative integer
tuples of length $b+1$. The Grundmann–Möller rule of order $s$ is

$$
Q_s^{(b)}f=
\sum_{i=0}^{s}
\frac{(-1)^i(b+2s-2i+1)^{2s+1}}
{2^{2s}i!\,(b+2s-i+1)!}
\sum_{|\beta|=s-i}
f\left(\frac{2\beta_0+1}{b+2s-2i+1},\ldots,
\frac{2\beta_b+1}{b+2s-2i+1}\right).
$$

It integrates every polynomial of total degree at most $2s+1$. Every
node is rational and strictly inside the simplex. Its weights may have
either sign and sum to $1/b!$, because the measure is ordinary simplex
Lebesgue measure. Attribution: A. Grundmann and H. M. Möller,
*Invariant Integration Formulas for the n-Simplex by Combinatorial
Methods*, SIAM Journal on Numerical Analysis **15** (1978), 282–290,
[doi:10.1137/0715019](https://doi.org/10.1137/0715019). **[src]**

For the product of $a$ positive-dimensional chain simplices, define
$\Delta_s^{(b)}=Q_s^{(b)}-Q_{s-1}^{(b)}$, with $Q_{-1}=0$. The sparse rule
is

$$
\sum_{s_1+\cdots+s_a\le s}
\Delta_{s_1}^{(b_1)}\otimes\cdots\otimes\Delta_{s_a}^{(b_a)}.
$$

**Proof of the total-degree guarantee.** For $s_j\ge1$, the $j$th
difference annihilates every block polynomial of degree below $2s_j$.
A monomial of total degree at most $2s+1$ cannot have block degree at
least $2s_j$ for every factor of a term with $\sum_js_j>s$. Thus all
omitted tensor differences vanish. The full tensor product telescopes to
a product rule exact on that monomial. Linearity proves the claim.

Equivalently, the implemented combination formula is

$$
\sum_{\max(0,s-a+1)\le|\alpha|\le s}
(-1)^{s-|\alpha|}\binom{a-1}{s-|\alpha|}
\bigotimes_{j=1}^{a}Q_{\alpha_j}^{(b_j)}.
$$

For the chain determinant with $k$ current factors, choose

$$s=\left\lfloor\frac{m-1+k}{2}\right\rfloor.$$

This rule integrates the complete fibre exactly. Duplicate rational
nodes are combined and zero weights removed. Zero-dimensional chain
factors contribute their single barycentric coordinate $1$ and volume
$1$. The implementation has explicit node and enumeration limits; a
request above those limits fails before returning an incomplete rule.

The independent reference moments in the tests are

$$
\int_{\Delta^b}\prod_{j=0}^b x_j^{\beta_j}\,dx
=\frac{\prod_j\beta_j!}{(b+|\beta|)!}.
$$

## 5. Graph-rank bounds for the Schwinger numerator

The existing numerator convention is

$$P=U^{n+H}D,\qquad\deg P=n^2-1,$$

as in [`cores/numer.py`](../cores/numer.py) and
[`assembly/numerator.py`](../assembly/numerator.py). The support bound
below strengthens the existing individual-edge bound and can constrain
interpolation before the fibre step or after it.

### Statement

For every nonempty edge set $S$, let $C_S$ be the submatrix of its cycle
columns. For a polynomial, $\deg_S$ means its largest total exponent in
the variables indexed by $S$. Then

$$
\boxed{\deg_S P\le n\,\operatorname{rank}C_S-1.}
$$

If $\operatorname{rank}C_S=0$, the assertion means $P=0$. The empty set
is excluded. In particular, the columns of one non-bridge series chain
have rank one, so

$$\sum_{e\in c}\text{exponent}_e\le n-1.$$

This is a bound on the sum over the chain, with no source-free assumption.
The cycle-column rank can also be computed without a symbolic matrix:

$$
\operatorname{rank}C_S=|S|+1-\#\text{connected components of }\Gamma\setminus S,
$$

where isolated vertices count as components.

### Proof of the rank bound

Use the tree current $y_0$ from Section 3 and set

$$b=C\operatorname{diag}(t)y_0,\qquad
Uy=Uy_0-C^{\mathsf T}\operatorname{adj}(M)b.$$

For each component and edge there is an augmented determinant identity:

$$
(Uy)_e^i=
\det\begin{pmatrix}
M&b^i\\ C_{\cdot e}^{\mathsf T}&(y_0)_e^i
\end{pmatrix}.
$$

Scale $t_e\mapsto zt_e$ for $e\in S$ and let
$r=\operatorname{rank}C_S$. The coefficient of $z$ in this augmented
matrix factors through $C_S$:

$$
\begin{pmatrix}C_S\\0\end{pmatrix}
\operatorname{diag}(t_S)
\begin{pmatrix}C_S^{\mathsf T}&(y_0)^i_S\end{pmatrix}.
$$

Its rank is at most $r$, so determinant multilinearity gives
$\deg_z(Uy)_e^i\le r$. Moreover $\deg_zU=r$: the positive
Cauchy–Binet expansion of $U$ contains a term using a maximal independent
subset of $S$, extended to a basis of the full cycle-column matrix.

For proper $S$, choose the omitted edge $e_0$ outside $S$. Every column
of the numerator matrix $C_{re}(Uy)_e^i$ has degree at most $r$. The
coefficient of $z^{r+1}$ in its voltage identity is

$$
\sum_{e\in S}t_e C_{re}(Uy)_e^{i,[r]}=0.
$$

Thus the highest-degree coefficient matrix of the chosen square minor
has a nonzero column dependence, supported on $S$. Its determinant has
degree at most $rnH-1$. Since

$$
P=\frac{(-1)^{e_0}
\det[C_{re}(Uy)_e^i]_{e\ne e_0}}
{t_{e_0}U^{n(H-1)}},
$$

division gives $rnH-1-rn(H-1)=nr-1$. For $S$ equal to the entire edge
set, the result is the homogeneous degree $n^2-1$. If $r=0$, each
selected cycle column is zero, and a minor omitting an edge outside $S$
vanishes identically.

For $k$ current factors, the polynomial numerator is
$P\prod_a(Uy)_{e_a}^{i_a}$ and obeys

$$
\deg_S\left(P\prod_a(Uy)_{e_a}^{i_a}\right)
\le(n+k)\operatorname{rank}C_S-1.
$$

After multiplying the series Jacobian, the pushed-forward core numerator
therefore has degree at most $n+k+L_c-2$ in a non-bridge chain resistance
$R_c$. Its total homogeneous degree is
$n^2-1+nk+(N-m)$.

### Algebraic justification of polynomiality

The formulas above use polynomiality already assumed by the existing
interpolators. For completeness it follows under the stated loopless
graph assumptions as follows. This argument is algebraic; it does not
infer divisibility from sampled values.

The Kirchhoff polynomial $U$ is multi-affine and hence squarefree.
For every irreducible factor $f$ of $U$, the matrix
$\operatorname{adj}(M)$ has rank at most one over the residue field at
$f$. Modulo $f$, the edge-by-component matrix $Uy$ consequently has
rank at most one. The matrix with entries $C_{re}(Uy)_e^i$ then has rank
at most $n$. A square minor has size $nH$, so exact elimination over the
discrete valuation ring at $f$ makes its determinant divisible by
$f^{n(H-1)}$. Squarefreeness gives divisibility by $U^{n(H-1)}$.

The voltage identity also shows that the minor omitting $e_0$ vanishes
when $t_{e_0}=0$: its remaining columns annihilate the nonzero remaining
time vector. Hence the minor is divisible by $t_{e_0}$. A non-loop edge
belongs to some spanning tree, so $U$ has a term independent of
$t_{e_0}$ and $\gcd(U,t_{e_0})=1$. Both denominators divide
simultaneously. Each $(Uy)_e^i$ has total degree $n$; subtracting the
denominator degree gives

$$n^2H-n^2(H-1)-1=n^2-1.$$

## 6. What is and is not computed

- `code/htk/series.py` performs the exact chain contraction, retains the
  complete original currents, and returns the density with its original
  sign. At fixed core resistances its remaining fibre is polynomial.
- `code/htk/dirichlet.py` integrates that polynomial with rational nodes
  and weights. Signed-weight cancellation is carried out over the
  rationals.
- The rank constraints remove impossible numerator exponents before
  polynomial reconstruction. They do not assert that every permitted
  coefficient is nonzero.
- This implementation treats $T=0$ and finite current monomials. A full
  displacement exponential is not a polynomial; its individual Taylor
  coefficients fit the method. The $T=1$ Gaussian contact terms require
  a separate extension.
- The remaining projective integral over the core resistances is still
  a period problem. Exact fibre integration does not supply missing
  high-level period tables or prove that the remaining periods belong
  to a prescribed transcendental number field.
- The best reconstruction route depends on the graph. A small fibre
  degree does not by itself establish an overall speedup: count both
  core samples and exact fibre evaluations, and compare with the
  constrained whole-graph grid for the actual source assignment.
