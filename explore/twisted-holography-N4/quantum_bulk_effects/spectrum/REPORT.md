# Direct higher-loop sensitivity and the physical bulk gap

## Result and scope

The boundary fan formula implies the exact selection rule

$$
\boxed{Q_n\mathcal O=0\quad\text{when}\quad 2J_L(\mathcal O)<2n-1.}
$$

This holds on the physical relative operator complex: gauge-invariant polynomials
in the retained letters and their holomorphic derivatives, independent of the
undifferentiated ghost. It applies before taking cohomology. It is a statement in
the manuscript's fan presentation of the quantum differential; an arbitrary
change of quantum operator scheme can redistribute individual higher-order
terms.

Consequently, the level-34 SO(13) BMN example, its level-40 mixed descendant,
and the listed pure glueball three-loop examples all have **only direct
$Q_0+\hbar Q_1$** on every cochain in their perturbative source chains.
Their nonzero third spectral differential is an effect of three successive
one-loop operations. This extends the paper's explicit argument for the pure
glueball examples to the BMN and mixed examples without needing their enormous
representatives.

A bulk reproduction of these examples would establish state recombination from
iterated one-loop dynamics. A direct boundary $Q_2$ test requires input
$J_L\ge3/2$; a direct $Q_3$ test requires $J_L\ge5/2$.
The existing all-order family saturates these bounds at every order and remains
the appropriate direct-operation test.

The separate Hilbert-space calculation below identifies the exact bulk observable
for a genuine lift: the positive excess of global anti-de Sitter energy over the
supersymmetry bound. It also shows precisely why a nonzero cochain coefficient
does not supply a normalized energy shift. No bulk action, bulk propagator, or
bulk quantum amplitude is evaluated here.

## 1. Proof of the spin selection rule

The physical letters in Table A.1 of the manuscript have the following
nonnegative undotted spin:

| Letter | Odd degree | Holomorphic derivative degree | $2J_L$ |
|---|---:|---:|---:|
| $\gamma^I$ | 1 | 0 | 0 |
| $\beta_I$ | 2 | 0 | 1 |
| $b$ | 3 | 0 | 2 |
| $\partial_{+\dot a}c$ | 0 | 1 | 0 |

Every additional $\partial_{+\dot a}$ raises $2J_L$ by one. For a
sector $(\theta_1,\theta_2,\theta_3;n_+,n_-)_Y$, equation (A.10) gives

$$
2J_L=\theta_1+\theta_2+\theta_3+n_++n_--Y.
\tag{1}
$$

An $n$-loop fan differentiates $n+1$ input letters. All functional
derivatives have the same three odd coordinates. A functional derivative of
an input letter of odd degree $0,1,2,3$ has odd degree $3,2,1,0$,
respectively. The total degree of these derivatives is therefore

$$
3(n+1)-\deg_\theta(\text{contracted input}).
$$

The two endpoint superfields add a nonnegative odd degree. Integration over
three odd coordinates forces

$$
\deg_\theta(\text{contracted input})\ge3n.
\tag{2}
$$

Applying (1) to the $n+1$ contracted letters gives

$$
\begin{aligned}
2J_L(\text{contracted input})
&=\deg_\theta(\text{contracted input})-(n+1)
  +\deg_\partial(\text{contracted input})\\
&\ge 3n-(n+1)=2n-1.
\end{aligned}
\tag{3}
$$

All remaining physical letters have nonnegative $J_L$, so (3) proves the
claim. Holomorphic jets, colour contractions, finite-rank trace identities,
and multi-trace products cannot weaken this necessary degree condition.
No statement is made about the unrestricted complex containing the
undifferentiated ghost, whose spin is negative.

For fixed input spin, the maximal possible direct loop order is

$$
n\le\left\lfloor J_L+\frac12\right\rfloor.
\tag{4}
$$

At every $n\ge2$, the previously calculated nonzero detector

$$
\operatorname{Tr}(\beta_1\gamma^2)
\operatorname{Tr}(\beta_2\gamma^3)
\operatorname{Tr}(\beta_3\gamma^1)
\left[\operatorname{Tr}(b\gamma^1-ig\beta_2\beta_3)\right]^{n-2}
\tag{5}
$$

has exactly

$$
J_L=\frac32+(n-2)=n-\frac12.
$$

Thus the bound is attained by a family with a proved nonzero direct
$Q_n$ matrix element. Each constituent of (5) is protected, and the
ordinary product is replaced by its quantum product. This proves sensitivity
to the quantum operation, not loss of the completed protected state.

## 2. The SO(13) BMN example contains no direct higher-loop term

The source in the current manuscript has sector

$$
(7,5,5;0,0)_{16},\qquad 2J_L=7+5+5-16=1.
$$

The cochain multiplying $\hbar^j$ lies in

$$
(7-j,5-j,5-j;j,j)_{16-j}.
$$

For every such correction,

$$
2J_L=(7-j)+(5-j)+(5-j)+j+j-(16-j)=1.
\tag{6}
$$

Hence

$$
Q_n\mathcal O_j=0\qquad(n\ge2)
\tag{7}
$$

on the entire correction chain. In particular, the complete order-three
equations reduce exactly to

$$
\begin{aligned}
Q_0\mathcal O_0&=0,\\
Q_0\mathcal O_1&=-Q_1\mathcal O_0,\\
Q_0\mathcal O_2&=-Q_1\mathcal O_1,\\
Q_0\mathcal O_3+Q_1\mathcal O_2&=\text{nonzero target class}.
\end{aligned}
\tag{8}
$$

The same identity holds for the mixed source:

$$
(5-j,3-j,3-j;3+j,3+j)_{16-j},\qquad 2J_L=1,
$$

and for every pure glueball source in this family:

$$
(3-j,3-j,3-j;p-3+j,p-3+j)_{2p+2-j},\qquad 2J_L=1.
$$

Multiplication by $J_L=0$ graviton factors does not change this conclusion.
This is a source-selection result, not a proof that higher spectral
differentials vanish: (8) is precisely how a nonzero third spectral
differential survives despite (7).

## 3. The physical quantity a bulk calculation must reproduce

In the convention of Chang and Yin, equation (2.1), the physical radial
quantization adjoint satisfies

$$
\boxed{2\{Q,Q^\dagger\}=E-2J_L-H_1-H_2-H_3.}
\tag{9}
$$

The manuscript writes the three R charges as $R_I$; they are the
$H_I$ in (9). For the SO(13) source and its order-three target, equation
(A.10) gives exactly

$$
\begin{array}{c|c|c|c}
&E&(H_1,H_2,H_3)&(J_L,J_R)\\\hline
\text{source}&33/2&(13/2,9/2,9/2)&(1/2,0)\\
\text{target}&17&(7,5,5)&(0,0).
\end{array}
\tag{10}
$$

Both have zero classical right-hand side in (9), and both have
$2(E+J_L)=34$. The target differs from the source by the supercharge
weights

$$
\left(\frac12;\frac12,\frac12,\frac12;-\frac12,0\right).
$$

Consequently their pairing has equal weights under every charge commuting with
$Q$, and opposite fermion parity. It cancels from the supersymmetric index.
An independent bulk computation can instead measure the positive eigenvalue
of the left-hand side of (9), or the corresponding shift of global
anti-de Sitter energies.

For an isolated physical pair, suppose

$$
Q\mathcal O_{\mathrm{source}}
=\kappa\hbar^r\mathcal O_{\mathrm{target}},\qquad
Q\mathcal O_{\mathrm{target}}=0,
$$

where $\kappa$ is the matrix element in the displayed, potentially
unnormalized state basis. If the physical adjoint preserves this pair and
annihilates the source, its action on the target is fixed by the inner product:

$$
Q^\dagger\mathcal O_{\mathrm{target}}
=\overline\kappa\hbar^r
\frac{\|\mathcal O_{\mathrm{target}}\|^2}
     {\|\mathcal O_{\mathrm{source}}\|^2}
\mathcal O_{\mathrm{source}}.
$$

Both states then have the exact common BPS excess

$$
\boxed{
E-2J_L-H_1-H_2-H_3
=2|\kappa|^2\hbar^{2r}
\frac{\|\mathcal O_{\mathrm{target}}\|^2}
     {\|\mathcal O_{\mathrm{source}}\|^2}.}
\tag{11}
$$

Rescaling either state changes $\kappa$ and the norm ratio by reciprocal
factors, leaving (11) invariant. A cohomological elimination is generally not
unitary; its effective two-state matrix cannot be inserted into (11) using
untransformed norms. The next calculation makes this issue explicit.

## 4. Exact finite Hilbert complex: high-page lifting from only one-loop terms

Take $r$ even and $r$ odd orthonormal states. The supercharge vanishes
on odd states. Its even-to-odd matrix is

$$
Q\big|_{\mathrm{even}\to\mathrm{odd}}=
\begin{pmatrix}
\hbar&1&0&\cdots&0\\
0&\hbar&1&\cdots&0\\
0&0&\hbar&\ddots&0\\
\vdots&\vdots&\ddots&\ddots&1\\
0&0&\cdots&0&\hbar
\end{pmatrix}.
\tag{12}
$$

This has only $Q_0+\hbar Q_1$, and its oddness makes $Q^2=0$.
At $\hbar=0$, its cohomology consists of the first even state and last
odd state. The exact completion identity is

$$
Q
\begin{pmatrix}1\\-\hbar\\\hbar^2\\\vdots\$-\hbar)^{r-1}\end{pmatrix}
=(-1)^{r-1}\hbar^r
\begin{pmatrix}0\\0\\\vdots\\0\\1\end{pmatrix}.
\tag{13}
$$

For every real $\hbar>0$, the smallest eigenvalue of
$2\{Q,Q^\dagger\}$, denoted below by the physical phrase “BPS gap”, obeys
the exact bounds

$$
\boxed{
\frac{2\hbar^{2r}}{\displaystyle\sum_{j=0}^{r-1}(j+1)\hbar^{2j}}
\le \text{BPS gap}\le
\frac{2\hbar^{2r}}{\displaystyle\sum_{j=0}^{r-1}\hbar^{2j}}.}
\tag{14}
$$

The upper bound is the Rayleigh quotient of (13). The inverse of (12) has
entry $(-1)^{j-i}\hbar^{-(j-i+1)}$ for $j\ge i$, and zero otherwise.
Its squared Frobenius norm is

$$
\hbar^{-2r}\sum_{j=0}^{r-1}(j+1)\hbar^{2j}.
$$

The smallest squared singular value is the reciprocal of the largest squared
singular value of that inverse; the Frobenius norm bounds the latter from
above. This gives the lower bound in (14).

For the three-step case, $\lambda$ denotes the smallest eigenvalue of
$Q^\dagger Q$ on the even states. Its complete exact equation is

$$
\boxed{
\lambda^3-(2+3\hbar^2)\lambda^2
+(1+2\hbar^2+3\hbar^4)\lambda-\hbar^6=0,
\qquad\text{BPS gap}=2\lambda.}
\tag{15}
$$

The root is unambiguously the smallest positive real root for
$\hbar>0$. The completion-vector quotient
$2\hbar^6/(1+\hbar^2+\hbar^4)$ is an upper bound, not the exact
answer. Its failure to give the exact gap is the effect of mixing with the
classically nonzero pairs.

For a fully exact coefficient check, the formal root at $\hbar=0$ obeys

$$
\begin{array}{c|rrrrrrrrrr}
\text{power of }\hbar&6&8&10&12&14&16&18&20&22&24\\\hline
\text{coefficient in the BPS gap}
&2&-4&2&12&-40&44&98&-520&882&640.
\end{array}
\tag{16}
$$

Equation (15), rather than a truncated series, specifies the result at every
value of the parameter.

## 5. Exact dependence on the Gram matrix

Replace the diagonal entries of (12) by
$a_1\hbar,\ldots,a_r\hbar$, and its superdiagonal entries by
$b_1,\ldots,b_{r-1}$, all nonzero. Let
$u_1,\ldots,u_r$ and $v_1,\ldots,v_r$ be the positive squared norms
of the orthogonal even and odd basis states, respectively. These quantities
are defined here only to exhibit the dependence on physical normalization.

The generalized characteristic polynomial for the squared singular values is

$$
\det\left(
Q^\dagger_{\mathrm{matrix}}\operatorname{diag}(v_1,\ldots,v_r)Q
-\lambda\operatorname{diag}(u_1,\ldots,u_r)
\right),
\tag{17}
$$

where $Q^\dagger_{\mathrm{matrix}}$ is the conjugate transpose of the
displayed matrix, not the adjoint in the non-unit Gram metric. Its constant
term is exactly

$$
\hbar^{2r}\prod_{i=1}^{r}|a_i|^2\prod_{i=1}^{r}v_i,
$$

and at $\hbar=0$ its coefficient of $\lambda$ is exactly

$$
-u_1\prod_{i=1}^{r-1}|b_i|^2\prod_{i=1}^{r-1}v_i.
$$

Consequently the first nonzero formal coefficient of the physical gap is

$$
\boxed{
[\hbar^{2r}]\,\text{BPS gap}
=2\left|\frac{\prod_{i=1}^{r}a_i}{\prod_{i=1}^{r-1}b_i}\right|^2
\frac{v_r}{u_1}.}
\tag{18}
$$

The bracket means the coefficient of the displayed power. Equation (18) is
valid for the specified block with parameter-independent positive metric.
If the physical norms carry powers of $\hbar$, their powers must be included
before assigning any physical perturbative order. In particular, a third
spectral differential in the twisted formal parameter does not, by itself,
fix a six-loop physical anomalous dimension or a power of the string coupling.

## 6. What this adds to the bulk program

1. **A complete exclusion for the existing lifting examples.** The direct
   $Q_{n\ge2}$ vanish on the SO(13) BMN and mixed chains as well as the
   pure glueball chains. Matching those lifts is a meaningful quantum bulk
   recombination test with a precise limitation.
2. **A sharp choice of direct-operation sources.** The first direct
   $Q_2$ probe has $J_L=3/2$; the existing three-factor detector has
   exactly that spin. At every order, (5) achieves the minimal input spin.
   A bulk construction that fails already on this source does not reproduce
   the complete boundary differential.
3. **A physical observable beyond the index.** The excess energy in (9)
   and the Gram-normalized matrix element in (11) distinguish actual pairing.
   Equations (14)–(18) are exact calibration checks for any proposed finite
   bulk state block.
4. **The remaining dynamical input.** Neither the source/target bulk states
   nor their physical inner products and interaction matrix elements have
   been constructed from the bulk action. Equations (11) and (18) therefore
   do not furnish a numerical anti-de Sitter energy shift for SO(13).

## Sources and reproducibility

- [Li–Liu–Chang, *Holomorphic Feynman diagrams in four dimensions
  and the quantum BPS ring*, arXiv:2611.24160v2](https://arxiv.org/pdf/2611.24160v2):
  the fan operation (8.3), the relative-letter conventions and charges in
  Appendix A, equation (A.9), equation (A.10), the explicit three-loop cases
  in Section 2.3 and Table 2.6, and the quantum-comparison qualification in
  Appendix B. The spin theorem above is derived from these formulas; it is
  not presented as a separately quoted theorem of that paper.
- [C.-M. Chang and X. Yin, *1/16 BPS states in N=4 super-Yang–Mills theory*,
  arXiv:1305.6314](https://arxiv.org/pdf/1305.6314), equation (2.1): the
  normalization of the physical anticommutator.
- Existing project calculation:
  `explore/twisted-holography-N4/boundary_bulk_all_orders/boundary/REPORT.md`,
  equations (1), (6), and (10): the saturating all-order detector and its
  functional derivatives.

Run `python verify.py` in this directory. Exact checks cover 66,040 contracted
letter inventories through loop order 32, 18 correction sectors, the complete
characteristic determinant and inverse identities for matrix sizes 1 through
8, the three-step gap coefficients through degree 24, and a symbolic
three-step nontrivial Gram metric. The unbounded selection theorem and matrix
bounds are proved above; these finite checks verify their implementation.
