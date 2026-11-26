# All boundary supercharges as bulk multiparticle tests

Research development, 26 November 2026.

**Scope clarification:** the target is an independent bulk derivation of the
$Q_n$ in the four-dimensional $\mathcal N=4$ paper. This note supplies
boundary targets and restricted bulk calculations; it does not establish
that derivation. See [Reconstructing the quantum supercharge from the
bulk](../independent_bulk_reconstruction/README.md) for the precise matching
equations, current-paper lifting targets, and a nonzero classical source
residue. The local primary-tree zero below is not a physical amplitude zero.

## Result

Every positive boundary correction $Q_n$ now has an explicit nonzero
test on products of individually protected single-graviton operators.
The test isolates exactly that boundary order, has a closed formula
for every $n$, and fixes its bulk particle number and ribbon genus.
The corresponding linear bulk input states and a possible protected
output channel are also determined at every order.

There is a definite remaining step: the direct local bulk trees on
these linear states vanish, while the ordinary boundary product has
the nonzero coefficient below. The quantum dictionary for composite
operators, including the appropriate external interaction kernels,
must reproduce the marked-product completion. The calculations here
give exact all-order targets for that construction; they do not claim
that the full punctured bulk amplitude has already been evaluated.

| Calculation | Result and proof |
|---|---|
| Every boundary $Q_n$ | [Protected constituents, exact isolation, full coefficient](boundary/REPORT.md) |
| Bulk particle number and genus | [Complete ribbon theorem and sewing](matching/REPORT.md) |
| Actual linear bulk states | [Laurent fields, full linear closure, tree vanishing, output channels](bulk/REPORT.md) |
| Entire tower and dictionary corrections | [Generating function and all-order Ward sum](resummation/REPORT.md) |

All new results start from repository revision
[`29ef4be3c92b8055ac51067357524c1048e6cda0`](https://github.com/lbotao399-stack/holomorphic-topological-feynman-diagrams/tree/29ef4be3c92b8055ac51067357524c1048e6cda0).
The boundary dynamical input is the all-loop fan and marked-product
construction of Li–Liu–Chang. The punctured fivefold and its flux are
the proposed twisted D3-brane bulk. Their normalization and scope
are retained explicitly.

## 1. Boundary order and bulk genus

Let $n$ be boundary perturbative order. A connected $Q_n$ fan
touches $r$ input traces and produces $b$ nonempty output traces.
Let $h$ be its ribbon genus and $f$ its number of empty color faces.
With ordinary fundamental traces, Euler's formula gives

$$
\boxed{
f=n+2-r-b-2h,\qquad
\hbar^n N^f=(\hbar N)^nN^{2-r-b-2h}.
}
$$

Here $N$ is matrix rank and $\hbar$ records boundary order.
The component coupling $g$ remains explicit in the formulas below;
no additional identification of it with $\hbar N$ is assumed.
There are $r+n+1$ ribbon vertices and $2n+1$ internal edges, so

$$
(r+n+1)-(2n+1)+(f+b)=2-2h.
$$

Thus, for a direct connected fan,

$$
n\ge r+b+2h-2,\qquad
r\le n+1,\qquad
h\le\left\lfloor\frac{n+1-r}{2}\right\rfloor.
$$

| Bulk channel under the single-trace dictionary | Genus | Lowest allowed boundary order |
|---|---:|---:|
| $1\to1$ | $0,1,2,3,4$ | $0,2,4,6,8$ |
| $2\to1$ or $1\to2$ | $0,1,2,3$ | $1,3,5,7$ |
| $3\to1$ or $2\to2$ | $0,1,2$ | $2,4,6$ |
| $(n+1)\to1$ | $0$ | $n$ |

The bounds are necessary; tensor or fermionic identities can still
annihilate an allowed coefficient. In particular the previous flat
genus-four scalar contact is not a boundary $Q_4$ matrix element
merely because both carry the integer four.

Sewing must preserve the same grading. For a connected graph of
vertex surfaces, its total genus equals the sum of their genera plus
the number of independent graph cycles. Two genus-zero pieces sewn
along two circles therefore contribute at genus one.

## 2. A protected family detecting every $Q_n$

Use the manuscript's fields, coupling $g$, fundamental trace
$\operatorname{Tr}(t_A t_B)=\tau\delta_{AB}$, and spinor derivatives
$\partial_{+\dot a}$. Define the four specified graviton components
by their full expressions:

$$
\begin{aligned}
G_1&=\operatorname{Tr}(\beta_1\gamma^2),&
G_2&=\operatorname{Tr}(\beta_2\gamma^3),\\
G_3&=\operatorname{Tr}(\beta_3\gamma^1),&
G_4&=\operatorname{Tr}(b\gamma^1-ig\beta_2\beta_3).
\end{aligned}
$$

They are components of $\operatorname{Tr}(dC)^2$, where $C$ is the
manuscript's superfield and $d$ is its superspace exterior derivative.
Consequently each satisfies $Q_kG_i=0$ for every $k\ge0$.
The first three are odd, and $G_4$ is even.

For $n\ge2$ take the ordinary product

$$
\mathcal O^{(n)}_0=G_1G_2G_3G_4^{n-2}.
$$

Set $\beta_I=b=0$ and
$\gamma^1=\gamma^2=\gamma^3=X$ with constant $X$.
Probe the output with two independent odd amplitudes
$\varepsilon_+,\varepsilon_-$ and traceless colors $D,E$.
The ghost waves are

$$
c(x)=
\varepsilon_+D\left(e^{-p_{+\dot a}x^{+\dot a}/2}-1\right)
+\varepsilon_-E\left(e^{-q_{+\dot a}x^{+\dot a}/2}-1\right).
$$

The exponent matches
$\partial_{+\dot a}=-2\,\partial/\partial x^{+\dot a}$.
The subtraction makes the undifferentiated ghost zero at the
insertion. The momenta $p,q$ are thus precisely the eigenvalues of
the spinor derivatives; no extra factor of four is implicit.

Write $\operatorname{ad}_X(Y)=[X,Y]$ for matrix commutation.
In the ordered $\varepsilon_+\varepsilon_-$ coefficient at the
origin, the complete supercharge gives

$$
\boxed{
\mathbf Q\mathcal O^{(n)}_0
\ \longmapsto\
-\frac{\hbar^n(ig\tau)^n}{(32\pi^2)^n(n!)^2}
\left(p_{+\dot a}q_+{}^{\dot a}\right)^n
\operatorname{Tr}\left((\operatorname{ad}_X)^{n+1}(D)E\right).
}
$$

The arrow here means exactly the stated evaluation and
coefficient extraction, not a bulk map. Every order other than $n$
has zero value in this same measurement.

### Exact order isolation

Choose $j$ factors from the $-ig\beta_2\beta_3$ term of $G_4$.
To survive at $\beta=b=0$, a fan must differentiate all
$3+2j$ beta letters and all $n-2-j$ copies of $b$.
If it also differentiates $a$ scalars, then

$$
k+1=n+1+j+a,\qquad
\text{Berezin degree}=3+2j+2a.
$$

The two output ghosts carry no odd coordinates.
The three-dimensional Berezin integral therefore requires
$j=a=0$, hence $k=n$. This proves the result for the complete
fan, rather than for a selected derivative term.

The $(n+1)!$ surviving assignments combine with the inherited
$1/((n!)^3(n+1))$ to give $1/(n!)^2$.
The two endpoint assignments give equal contributions because
the ghost exchange, spinor exchange, and adjoint-chain reversal
signs multiply to one. The boundary report derives the remaining
phase and every normalization factor.

### A nonzero value at every matrix rank

Let $E_{ij}$ be the matrix unit. For any $N\ge2$, use the upper-left
two-by-two block

$$
X=\operatorname{diag}(1/2,-1/2,0,\ldots,0),\qquad
D=E_{12},\quad E=E_{21},
$$

Then $[X,D]=D$ and $\operatorname{Tr}(DE)=1$.
At $\tau=1/2$ and $p_{+\dot a}q_+{}^{\dot a}=1$,

$$
\boxed{
\mathbf Q\mathcal O^{(n)}_0
\ \longmapsto\
-\frac{(ig\hbar)^n}{(64\pi^2)^n(n!)^2}\ne0.
}
$$

For the omitted order one, use the protected ordinary product
$\operatorname{Tr}((\gamma^1)^2)\operatorname{Tr}(\beta_1\gamma^2)$.
Its same normalized measurement is exactly
$\hbar g^2/(32\pi^2)$, and every higher order vanishes.

### Intrinsic order and survival of the completed states

Each $Q_n$ has at most $n+1$ functional derivatives. The nonzero
$(n+1)$-factor measurement above shows that its differential order
with respect to ordinary multiplication is exactly $n+1$.
No multiplicative formal change of variables can remove that highest
symbol: at order $n$, transformations of lower $Q_k$ have
differential order at most $k+1<n+1$.

This does not lift the protected multiparticle states. Replace
ordinary multiplication by a fixed bracketing of the marked quantum
product. The resulting finite completion is closed to every order.
Its nonexactness follows from a separate constant-Cartan evaluation:
take $\gamma^I=X$, $\beta_I=\chi_I X$, $b=vX$, and $c=0$, with
independent odd $\chi_I$. Every differential image evaluates to
zero, whereas the leading product evaluates to

$$
(\operatorname{Tr}X^2)^{n+1}
\chi_1\chi_2\chi_3v^{n-2}\ne0.
$$

Thus these tests detect required quantum dressing of surviving
multiparticle states. The nonzero fan measurement uses off-diagonal
ghosts; the nonexactness test puts every color on the Cartan.

## 3. Their explicit bulk states

Use tangent coordinates $z_1,z_2$, normal coordinates
$w_1,w_2,w_3$, and the tangent bivector
$P=\partial_{z_1}\wedge\partial_{z_2}$.
The linear dictionary identifies each of the four input sectors with
a one-dimensional space of exact two-forms on
$\mathbb C^{2|3}$. For example, the $G_1$ sector has the exact basis

$$
d(\theta_2\theta_3\,d\theta_2)
=\theta_3(d\theta_2)^2-\theta_2d\theta_2d\theta_3
$$

and the corresponding Laurent field

$$
P\wedge\left(
\frac{2\partial_{w_3}}{w_1w_2^3w_3}
-\frac{\partial_{w_2}}{w_1w_2^2w_3^2}
\right).
$$

The bulk report lists all four complete fields. They are closed
under the linear differential

$$
D=\bar\partial+t\partial_\Omega+N\{F,-\},\qquad
F=\frac{P}{w_1w_2w_3},
$$

where $t$ is descendant level, $\partial_\Omega$ is holomorphic
divergence, and $\{F,-\}$ is the Schouten action of the flux.
They are nonzero level-zero classes. Their exact-form description
does not make them bulk boundaries: the incoming flux image
corresponds to one-forms, and descendant divergence has a factor $t$.

Every input field and every separate flux insertion contains $P$.
Holomorphic derivatives preserve this exterior factor. At the
nonlinear vertex furthest from the root of any finite sewing tree,
at least two incoming leaves therefore multiply to zero:

$$
(P\wedge Y_1)\wedge(P\wedge Y_2)=0.
$$

Root divergence and subsequent internal linear propagators preserve
zero. This proves the all-order vanishing of local primary trees
with these linear fields attached directly. Positive-genus local
vertices obey the same argument when their sewing graph is a tree.
The statement does not cover external bulk-to-boundary kernels or
sewing graphs with cycles.

The stable protected output channel is also explicit. Its charges
are $(n-1,1,1;n,n)$ and its letter number is $n+3$.
Its unique exact form is

$$
z_1^{n-1}z_2^{n-1}\,
dz_1\wedge dz_2\,
(d\theta_1)^{n-1}d\theta_2d\theta_3,
$$

whose Laurent representative is

$$
\boxed{
(-1)^n(n-1)!
\frac{z_1^{n-1}z_2^{n-1}}{w_1^n w_2^2w_3^2}.
}
$$

This is a possible protected channel at the required charges.
The raw joining cochain is not identified with this state.
Finite-rank trace relations must also be imposed when using this
stable single-particle description.

## 4. Exact requirements on the composite map

For fixed $n$, write $\mathcal O^{(n)}_j$ for the coefficient of
$\hbar^j$ in the marked product of the same protected constituents.
Its classical term is $\mathcal O^{(n)}_0$.
Quantum closure gives

$$
Q_n\mathcal O^{(n)}_0+
\sum_{j=1}^n Q_{n-j}\mathcal O^{(n)}_j=0.
$$

Consequently the normalized two-ghost measurement fixes the exact
sum of the required corrections at every order:

$$
\boxed{
\sum_{j=1}^n Q_{n-j}\mathcal O^{(n)}_j
\ \longmapsto\
\frac{(ig)^n}{(64\pi^2)^n(n!)^2},
\qquad n\ge2.
}
$$

It does not fix each summand separately.
The [resummation report](resummation/REPORT.md) sums the whole family
using $G_1G_2G_3/(1-uG_4)$, where $u$ counts the even constituent.
Its exact coefficient series is

$$
-\sum_{n\ge2}
\frac{(ig\hbar)^n u^{n-2}}{(64\pi^2)^n(n!)^2},
$$

with a closed expression in the modified Bessel function and an
opposite-sign series for the required corrective Ward terms.

These are concrete tests for the nonlinear holographic map.
The topology is fixed: every member is maximal joining,
$(n+1)\to1$, with $h=f=0$, and hence has bulk tree weight $N^{-n}$
after factoring out $(\hbar N)^n$.
The input states are fixed, the candidate protected output sector
is fixed, and the boundary coefficient including phase and factorial
is fixed. The remaining calculation is the actual relative
open/closed composite-observable map or external kernel that
realizes this Ward-compatible result in the punctured bulk.

## 5. Reproduction

From the repository root:

```bash
python -m pip install -r explore/twisted-holography-N4/boundary_bulk_all_orders/requirements.txt
python explore/twisted-holography-N4/boundary_bulk_all_orders/run_checks.py
```

Four independent certificates verify the boundary Berezin and color
calculation, complete ribbon topology, the bulk exterior complex,
and the all-order coefficient series. The proofs establish the
unbounded statements. Finite exact checks independently validate
their coefficients, signs, and constructions.

The replay records source/output hashes and requires saved outputs
to reproduce byte for byte. No inherited graph engine is imported
by these four certificates.

## References

1. B. Li, K. Liu, C.-M. Chang,
   *Holomorphic Feynman diagrams in four dimensions and the quantum
   BPS ring*, [arXiv:2611.24160v2](https://arxiv.org/abs/2611.24160v2).
   The pinned repository manuscript supplies the all-loop fan,
   single-graviton closure, and finite marked-product completion.
   Its interpretation in the untwisted physical theory retains the
   stated assumption that twisting commutes with quantization of
   composite operators.
2. K. Costello, S. Li,
   *Twisted supergravity and its quantization*,
   [arXiv:1606.00365](https://arxiv.org/abs/1606.00365), Section 9.
   This supplies the proposed punctured fivefold bulk and D3 flux.
3. The four detailed reports give the pinned repository input
   files, exact conventions, standard antecedents, and scope of
   every new calculation.
