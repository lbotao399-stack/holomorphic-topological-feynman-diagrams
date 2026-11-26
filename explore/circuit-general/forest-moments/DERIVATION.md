# Exact moment recursion for cubic product forests

This is a coefficient algorithm for the complete polynomial current integrals in
[the manuscript, Section 7](../../../sources/LLC_main.tex). It implements the same
series, single-star, and paired-star changes of variables as
[the existing ratio-chart evaluator](../../nonSP/verify/full_ward/forest_jets.py),
but integrates the newly introduced ratios before moving to the next vertex.
It never constructs all original currents as polynomials in all `2*n` ratios.

The implemented domain is a connected simple cubic product forest on
`C^2`: deleting the endpoints `A,B` of the marked edge leaves trees, each with
one vertex of full degree two and all other vertices of full degree three.
Only that degree-two root can carry a source inside its tree. Sources at all
roots and both insertions are retained. They may be arbitrary rational
two-component spinors satisfying conservation. The construction accepts a
finite monomial in currents on **any original edges**, including internal
edges and the marked edge. Linearity gives arbitrary finite polynomial
insertions.

Throughout, an oriented current from `v` to `w` is denoted by `y_vw`, and

\[
\langle y,z\rangle=y_{+\dot a}z_+{}^{\dot a}
\]

is the manuscript's antisymmetric spinor contraction. Scalar polynomial
coordinates in the implementation are an exact representation of these
spinors; no source direction or current direction is discarded.

## 1. Local moment operators

Let `F` be a polynomial in the original currents. Currents on unlisted edges
are left unchanged in each of the following formulas.

### A single source-free star

Use the manuscript's local edge orders

\[
(vu,vw,vx,ux,wx)\longrightarrow(uw,ux,wx).
\]

Write the three retained currents as
`y = y_uw`, `l = y_ux`, and `h = y_wx`. The exact operator is

\[
\begin{aligned}
&\int_0^1\!\int_0^1
r_1r_2\langle l,h\rangle\,
F\bigl(-y-(1-r_1)l,\ y-(1-r_2)h,\
       (1-r_1)l+(1-r_2)h,\ r_1l,\ r_2h\bigr)
\,dr_1\,dr_2.
\end{aligned}
\]

In particular, the ratio-free factor for the constant insertion is
`<l,h>/4`. A coefficient carrying
`r1^a*(1-r1)^b*r2^c*(1-r2)^d` is multiplied by

\[
\frac{(a+1)!b!}{(a+b+2)!}
\frac{(c+1)!d!}{(c+d+2)!}.
\]

This formula includes the two ratio factors from the density. All four
exponents are nonnegative integers.

### Paired source-free stars

Use local edge orders

\[
(pu,pA,pB,qu,qA,qB,AB)\longrightarrow(uA,uB,AB),
\]

and retained currents `l = y_uA`, `h = y_uB`, `y = y_AB`.
Set `s = 1-r3-r4` only in this local formula. The inverse map is

\[
\begin{array}{lll}
y_{pu}=-r_1l-r_2h,&y_{pA}=r_1l-r_3y,&y_{pB}=r_2h+r_3y,\\
y_{qu}=-(1-r_1)l-(1-r_2)h,
 &y_{qA}=(1-r_1)l-r_4y,
 &y_{qB}=(1-r_2)h+r_4y,\\
y_{AB}^{\rm old}=s\,y.
\end{array}
\]

Multiply `F` evaluated on these currents by

\[
(r_1-r_2)s\,\langle l,y\rangle\langle h,y\rangle
\]

and integrate over the two intervals and the triangle
`r3,r4 > 0`, `r3+r4 < 1`. For the Bernstein monomial

\[
r_1^a(1-r_1)^b r_2^c(1-r_2)^d r_3^e r_4^f s^g
\]

the complete scalar moment is

\[
\left[
\frac{(a+1)!b!}{(a+b+2)!}\frac{c!d!}{(c+d+1)!}
-\frac{a!b!}{(a+b+1)!}\frac{(c+1)!d!}{(c+d+2)!}
\right]
\frac{e!f!(g+1)!}{(e+f+g+3)!}.
\]

The constant insertion therefore vanishes. The first current moments are
already nonzero before the root is removed:

\[
\begin{aligned}
F=y_{pA}:&\qquad
 \frac{l}{72}\langle l,y\rangle\langle h,y\rangle,\\
F=y_{pB}:&\qquad
-\frac{h}{72}\langle l,y\rangle\langle h,y\rangle,\\
F=y_{pu}:&\qquad
 \frac{h-l}{72}\langle l,y\rangle\langle h,y\rangle.
\end{aligned}
\]

The moments on the corresponding `q` edges are their negatives. These
identities explain why a zero bare kernel cannot justify dropping a tree
from the displacement expansion.

### Series removal of a root

Let the root source be `lambda`, the retained marked current be `J`, and
let `s,r` be the series and parallel fractions. The inverse map is

\[
l=s\lambda-(1-r)J,\qquad
h=(1-s)\lambda+(1-r)J,\qquad
y_{AB}^{\rm old}=rJ.
\]

The density contributes `r*<lambda,J>`. The reduced insertion source is

\[
\lambda_A^{\rm red}=\lambda_A+s\lambda,
\qquad
\lambda_B^{\rm red}=\lambda_B+(1-s)\lambda.
\]

The source shift must be applied to the remaining forest kernel before
integrating `s`. Treating the reduced source as independent of `s` would
give an incorrect answer even for two disconnected one-vertex trees.

## 2. Why only three boundary currents are needed

Root each interaction tree at its unique source vertex. A non-root vertex
has two downward ports. Each port is an insertion spoke or a child subtree.
A leaf has one `A` and one `B` spoke. A vertex with one child is reduced by
the single-star rule, and one with two children by the paired rule.

After its descendants are removed, a subtree is seen by its parent only
through the three currents on its two root spokes and the marked edge.
Its parent edge is unchanged during these descendant reductions. At a
source-free child its current in the parent-to-child orientation is
`l+h`, by current conservation. This also retains insertions on internal
edges: they are included as powers of that exact current when the child
is attached to its parent.

Every reduction rescales the old marked current by a scalar fraction;
its spinor direction is unchanged. Consequently an arbitrary homogeneous
polynomial of marked-current degree `q` factors out after a subtree is
integrated. Define its subtree message to be the polynomial multiplier
in the three retained currents. The multiplier depends on `q`, but not
on which degree-`q` polynomial was factored out.

There is an essential degree shift when two disjoint descendants are
composed. First decompose the first message into homogeneous degrees `h`
in its marked current. The second descendant then sees marked-current
degree `q+h`, since it rescales both the outside degree-`q` polynomial
and the marked current inside the first message. Thus the composition is

\[
\sum_{h\ge0}
\bigl(\text{marked-degree-}h\text{ part of the first message at }q\bigr)
\bigl(\text{second message at }q+h\bigr).
\]

The displayed factors still have their own boundary spokes; the final
single or paired current map places both in the same three retained
currents. Multiplying messages at the same `q` without this shift is
incorrect. Exact degree grading is what replaces a dense polynomial in
all forest ratios.

## 3. Exact assembly of all component trees

For an exponent pair `alpha`, consider the remaining forest moment with
the marked-current monomial `J^alpha` inserted. Its dependence on the two
coordinates of `lambda_A` is polynomial; all root sources are held fixed
and `lambda_B` is determined by conservation.

For the next root, take its descendant message at `q = |alpha|` and form

\[
r^{|\alpha|+1}\langle\lambda,J\rangle J^\alpha
\,\text{message}_{|\alpha|}
\bigl(s\lambda-(1-r)J,(1-s)\lambda+(1-r)J,rJ\bigr).
\]

Integrate `r` exactly and expand the result as a finite sum of monomials
`c_(p,beta)*s^p*J^beta`. The required forest moment is then

\[
\sum_{p,\beta}c_{p,\beta}
\int_0^1s^p\,
\text{remaining-forest moment}_{\beta}
       (\lambda_A+s\lambda)\,ds.
\]

For an empty forest the moment is simply `lambda_A^alpha`.
This recursion accounts for every source transfer. It involves only the
two symbolic insertion-source coordinates and one temporary shift
variable. Equivalently, if the remaining polynomial is `f`, its integral
is the finite identity

\[
\int_0^1s^p f(\lambda_A+s\lambda)\,ds
=\sum_{j=0}^{\deg f}
\frac{(\lambda_{+\dot a}\,\partial/\partial\lambda_{A,+\dot a})^j
      f(\lambda_A)}{j!(p+j+1)}.
\]

The derivatives here describe exact translation of a polynomial. There
is no series cutoff or omitted derivative order.

## 4. Correctness, signs, and finite size

Each local identity is the pullback of the full current-area form in the
manuscript. Its ratio domain is a product of compact intervals and
simplices. After fixing the sources and a finite current word, all
integrands are polynomials. Fubini's theorem therefore permits each
fresh ratio block to be integrated immediately. Induction over the
rooted tree gives the boundary message, and induction over the component
trees gives the source-translation recursion. This proves equality to
the original complete ratio-chart integral for every finite current
word, at every loop order.

The single-star resistance change reverses orientation and cancels the
minus sign in its local current-area wedge. The paired-star resistance
change preserves orientation. The formulas in Section 1 use exactly the
local edge orders displayed there. The implementation also transports
the original edge ordering and current orientations. Reversing an edge
negates a current inserted on that edge; permuting edges changes the
kernel by the appropriate permutation sign. In the historical sorted
chart convention the accumulated factor is

\[
(-1)^{n(n+1)/2}
\prod_{\text{reduction steps}}(\text{local edge-order sign}).
\]

The graph model converts between its canonical rooted representation and
the original graph before returning a moment. In particular the bare
kernel sign and the signs of inserted oriented currents are separate.

For an `n`-loop forest and a requested current word of total degree `k`,
the full current degree is `2*n+k`. A descendant message removing `m`
loops and carrying `h` original insertions is homogeneous of degree
`2*m+h` in its three boundary spinors. The outside marked-current degree
and this message degree together never exceed `2*n+k` in the recursion.
There are six scalar boundary coordinates, so a conservative bound on
stored rational boundary coefficients is

\[
n(2n+k+1)\binom{2n+k+6}{6}.
\]

This bound includes all marked-degree states and allows inhomogeneous
messages; it is deliberately larger than the actual homogeneous count.
Local substitution introduces at most four new ratio variables, and
integrates them immediately. Thus the number of variables per local
operation is bounded independently of the loop order. The root-assembly
states have two source variables and two marked-current exponents.
These are bounds on coefficient counts, not on integer bit lengths or
on the cost of requesting every possible word. All words through degree
`k` on the `2*n` physical edges already have
`binom(4*n+k,k)` scalar component monomials.

The two source coordinates at `A` have degree at most `n+k`: each of the
`n` spinor brackets in the density has degree at most one in this source,
since the quadratic term contracts the source with itself; each current
factor adds at most one. This also bounds the exact translation above.
Keeping other root sources symbolic is mathematically allowed, but
their coefficient ring then has its own unavoidable output size.

The implementation uses explicit state and term budgets. A budget
failure is an error, not a partial polynomial or an approximate result.
There is no built-in loop cutoff or current-order cutoff.

## 5. Batched Taylor transfer for a fan

The ordered-uniform moment identity used here is already implemented and
independently checked in
[scalar_fan_joint.py](../../iib-bcov/scalar_fan_joint.py), function
`ordered_uniform_moment`. The new
[fan_moments.py](../../../code/htk/fan_moments.py) provides a batched
coefficient transfer and accepts displacements on every original graph
edge. It does not claim a new identity for the underlying scalar fan.

Use the graph order of `htk.graph.fan(n)`: path edge `(i,i+1)` followed by
spoke `(i,hub)`, for `i=0,...,n-1`, and then the final spoke. Only the two
path endpoints carry interaction sources. Write them as `lambda_first`
and `lambda_last`; the hub source is their negative sum. The bare kernel
in the `htk` convention is

\[
\frac{(-1)^{n(n+1)/2}
\langle\lambda_{\rm first},\lambda_{\rm last}\rangle^n}
{(n!)^3(n+1)}.
\]

On spoke `i`, let its first- and last-endpoint weights be `a_i,b_i`, so
its current is `a_i*lambda_first+b_i*lambda_last`. For nonnegative
exponents `p_i,q_i`, the existing normalized moment identity is

\[
\begin{aligned}
\left\langle\prod_{i=0}^{n}a_i^{p_i}b_i^{q_i}\right\rangle
={}&\frac{(n!)^2(n+1)!\prod_{i=0}^{n}(1+p_i+q_i)!}
{(n+\sum_i p_i)!(n+\sum_i q_i)!}\\
&\times\int_{0<x_n<\cdots<x_0<1}
\prod_{i=0}^{n}x_i^{p_i}(1-x_i)^{q_i}\,dx_0\cdots dx_n.
\end{aligned}
\]

The angle brackets in this one formula denote normalized scalar moments;
all spinor contractions retain their two arguments.

First remove path displacements by the exact conservation identity

\[
y_{(i,i+1)}=\lambda_{\rm first}
             -\sum_{j=0}^{i}y_{(j,\mathrm{hub})}.
\]

Thus the displacement exponential acquires the constant linear phase

\[
-\lambda_{\rm first,+\dot a}
  \sum_{i=0}^{n-1}z_{(i,i+1)}^{+\dot a},
\]

while the effective displacement on spoke `j` is its original
displacement minus the sum of path displacements from `j` to the final
endpoint. The last spoke has no such sum. This identity retains all
original path-current moments and their source phase.

For each effective spoke displacement, denote its contractions with
`-lambda_first` and `-lambda_last` by `A_i,B_i`. The factorial in the
ordered moment combines with the exponential's factorials into

\[
\sum_{p,q\ge0}\frac{(p+q+1)!}{p!q!}
 (sA_i x)^p(tB_i(1-x))^q
=\frac{1}{(1-sA_i x-tB_i(1-x))^2}.
\]

Here `s,t` are formal coefficient tags for total first- and last-endpoint
degree. They introduce no new integration variables. After `m` ordered
integrations, a state `(A,B,j)` carries the monomial `x^(m+A+j)`, with
`0 <= j <= B`. Appending local exponents `(p,q)` and choosing the term
`ell` in `(1-x)^q` multiplies its coefficient by

\[
\frac{(p+q+1)!}{p!q!}
A_i^pB_i^q(-1)^\ell\binom q\ell,
\]

and then divides by

\[
(m+1)+(A+p)+(j+\ell).
\]

The new state is `(A+p,B+q,j+ell)`. Start with `(0,0,0)` of coefficient
one and process the `n+1` spokes in reverse order. At the final upper
endpoint `x=1`, sum over `j` and multiply each `(A,B)` coefficient by

\[
\frac{(n!)^2(n+1)!}{(n+A)!(n+B)!}.
\]

Collect terms with `A+B=k`, multiply by the bare kernel, and convolve
with the exact Taylor coefficients of the constant source phase.
This yields the complete coefficient of the common displacement scale
through any requested finite order.

For requested order `k`, discard a state only if `A+B>k`: subsequent
moves can only increase this degree, so no requested coefficient is
affected. The exact state count is

\[
\sum_{A=0}^{k}\sum_{B=0}^{k-A}(B+1)=\binom{k+3}{3},
\]

independent of the loop number. The direct transfer uses at most a
quadratic number of these states per spoke. This count concerns exact
rational coefficient operations; integer bit lengths still grow with
the calculation. The finite coefficient request does not assign an
approximate value to the kernel at a nonzero displacement.

[test_fan_moments.py](../../../code/tests/test_fan_moments.py) compares
all-original-edge directional Taylor coefficients and mixed raw moments
to the pre-existing series-parallel integrator through four loops. It
also checks the original ordered-uniform identity, exponents beyond
255, and the complete gradient-displacement phase identity. For example,

\[
\langle a_0^{256}\rangle_{n=1}=\frac1{129},\qquad
\langle b_0^{256}\rangle_{n=1}=\frac1{33153}.
\]

## 6. Independent controls and the noncubic boundary

The one-loop tree `(AB)` gives

\[
\mathcal I_{(AB)}=-\frac{\langle\lambda_0,\lambda_A\rangle}{2}.
\]

At two loops the two chains `(A(AB))`, `(B(AB))` in lexicographic edge
order give respectively

\[
\frac{\langle\lambda_0,\lambda_A\rangle^2}{24},\qquad
-\frac{\langle\lambda_0,\lambda_A\rangle^2}{24}.
\]

For the two separate roots `(AB),(AB)`, the complete source polynomial is

\[
\begin{aligned}
\mathcal I={}&
-\frac{\langle\lambda_0,\lambda_A\rangle
        \langle\lambda_1,\lambda_A\rangle}{6}\\
&+\frac{(\langle\lambda_0,\lambda_A\rangle
         -\langle\lambda_1,\lambda_A\rangle)
         \langle\lambda_0,\lambda_1\rangle}{12}
+\frac{\langle\lambda_0,\lambda_1\rangle^2}{24}.
\end{aligned}
\]

The last two terms are an explicit control on the source translation:
they survive for nonparallel root sources.

For the three-loop fork with root `u=0`, leaves `p=1,q=2`, insertions
`A=3,B=4`, and lexicographic edge order, direct evaluation of the local
moments above gives

\[
\begin{aligned}
\mathcal I[y_{pA,+\dot a}]&=
\frac{\langle\lambda_u,\lambda_A\rangle^3}{8640}
(2\lambda_u-\lambda_A)_{+\dot a},\\
\mathcal I[y_{pB,+\dot a}]&=
-\frac{\langle\lambda_u,\lambda_A\rangle^3}{8640}
(\lambda_A+3\lambda_u)_{+\dot a},\\
\mathcal I[y_{up,+\dot a}]&=
-\frac{\langle\lambda_u,\lambda_A\rangle^3}{8640}
(2\lambda_A+\lambda_u)_{+\dot a}.
\end{aligned}
\]

Interchanging `p,q` negates these expressions, and the bare fork is zero.
For example `lambda_A=(1,0)`, `lambda_u=(0,1)` gives the complete spinors

\[
\mathcal I[y_{pA}]=(1/8640,-1/4320),\qquad
\mathcal I[y_{pB}]=(1/8640,1/2880),\qquad
\mathcal I[y_{up}]=(1/4320,1/8640).
\]

[tests_exact.py](tests_exact.py) compares the recursion to the pre-existing
full ratio charts on all cubic forest topologies through five loops,
with independent nonparallel sources at every root and insertions on
the original edges: 100 topologies and 2,292 exact scalar comparisons.
For each graph this includes the bare kernel, every scalar component of
every original edge current, a repeated-current word of degree three,
and a mixed physical/marked-current word of degree five. All original
vertex conservation equations are also checked on the first moments.
Additional tests evaluate source polynomials at rational source values,
transport arbitrary words under edge reversals and permutations, and
compare asymmetric branched trees at six and nine loops with the old
global chart. The all-order fan and repeated-fork formulas already
recorded in the repository are controls, not newly inferred formulas.

For a higher-degree superpotential a component can carry several source
vertices. The manuscript's sourced-star appendix then transfers sources
with coefficients such as

\[
\frac{t_{wx}^{\rm red}}{(1-r_2)
 [t_{uw}^{\rm red}+t_{ux}^{\rm red}/(1-r_1)
                       +t_{wx}^{\rm red}/(1-r_2)]}.
\]

Those coefficients couple the local ratios to the remaining resistances.
The source-free polynomial operators above do not apply to that case.
General product-forest recognition and any further reduction must retain
this distinction; an arbitrary nonzero source cannot be silently erased.

The public [product-forest adapter](../../../code/htk/product_forests.py)
uses the cubic recursion only after validating this source condition.
For a broader forest that is series-parallel it instead uses the existing
generic current-Jacobian integrator with every original source retained.
This adds coverage through an already established method. A sourced
non-series-parallel remainder raises `UnsupportedProductForestReduction`;
the rational sourced-star coefficients above have not been integrated by
the cubic message engine.
