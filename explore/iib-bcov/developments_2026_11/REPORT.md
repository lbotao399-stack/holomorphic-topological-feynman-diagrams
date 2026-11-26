# Exact Ward targets for planar bulk currents in twisted type IIB theory

## Results and scope

We calculate exact differences of complete cyclic compositions of the planar local binary operations in the normal cubic background. At total genus four, seven times the coefficient at momenta

$$
p=(1,0),\qquad q=(0,1),\qquad r=(2,3)
$$

minus 162 times the coefficient at the same $p,q$ and $r=(1,2)$ is

$$
\boxed{-\frac{s(B)^2}{46656}}.
$$

Here $s(B)$ is the ratio of the measured genus-two binary coefficient at $B$ to its value at $xyz$. The entire genus-four coefficient is multiplied by $(2\pi i)^{20}$, and the common output wave is removed separately from each matrix element. The normal cubic invariant $s(B)$ and the heat convention are specified below. This equality holds for every value of the previously uncomputed genus-four binary coefficient: it cancels exactly. We obtain corresponding exact targets at genus five, and at genus six on $s(B)=0$. The genus-four and genus-six targets together are nonzero on every smooth normal projective cubic.

These are identities for the isolated local binary operations on specified cochains. They constrain the additional terms required in a complete quantum Ward identity. They do not establish that this binary induces an operation on the cohomology of the full theory or on D3-brane charges. The genus-four through genus-six results below involve no new integration of a higher-genus period.

The restricted Moyal rigidity statement used to organize the test is standard deformation theory. Its application, the new integer-momentum witness, the complete cyclic compositions, and the differences that eliminate unknown coefficients are the development in this note. The earlier repository already established nonzero first-order Jacobi matrix elements at genus two and three.

## 1. Prior input and the exact field sector

We use the local fivefold Bershadsky–Cecotti–Ooguri–Vafa theory as developed in the repository at commit [`aa0d2f1725dff135494eef132f5e92a8c0401a42`](https://github.com/lbotao399-stack/holomorphic-topological-feynman-diagrams/tree/aa0d2f1725dff135494eef132f5e92a8c0401a42). The source is the section [“Planar primary brackets and the normal elliptic curve”](https://github.com/lbotao399-stack/holomorphic-topological-feynman-diagrams/blob/aa0d2f1725dff135494eef132f5e92a8c0401a42/explore/iib-bcov/README.md#planar-primary-brackets-and-the-normal-elliptic-curve), with its [existing graph evaluator](https://github.com/lbotao399-stack/holomorphic-topological-feynman-diagrams/blob/aa0d2f1725dff135494eef132f5e92a8c0401a42/explore/iib-bcov/planar_cubic_primary.py).

The construction of the all-genus local interaction is given by Si Li and Minghao Wang, [*On the exact quantization of BCOV Theory on Calabi-Yau Manifolds*](https://arxiv.org/html/2610.21927v1), equations (1.9)–(1.11), Lemma 1.2, and Theorem 1.11. That paper is an announcement, explicitly reserving a detailed account of its proofs for a subsequent paper. The exact algebraic results here require the stated local binary coefficients and their determinant form; a complete physical Ward interpretation additionally uses the announced quantum master equation.

Let $x,y,z$ be the three normal coordinates and $z^1,z^2$ the two tangent coordinates. The background $B(x,y,z)$ is a homogeneous cubic scalar primary field. A positive-genus primary local operation with $k$ tangent inputs uses $2g+2-k$ background slots. Its $2g$ five-dimensional determinants require $6g$ normal derivatives, whereas the cubics supply at most

$$
3(2g+2-k).
$$

Consequently $k\leq 2$. Rank in the tangent plane requires $k\geq2$. Thus $k=2$, and every one of the $2g$ normal slots is differentiated exactly three times. The resulting coefficient is constant in the normals. The isolated binary consequently closes on tangent Hamiltonian waves, and on their action on tangent scalar waves, at the level of these representatives.

Writing the primary field as $A_0$ and its holomorphic volume trace as $\operatorname{Tr}$, the classical primary interaction is cubic:

$$
\frac{1}{3!}\operatorname{Tr}(A_0^3).
$$

In particular its binary contains no background insertion and retains the ordinary linear Poisson determinant. For tangent momenta, the repository notation is

$$
[p,q]=p_1q_2-p_2q_1,
\qquad
X_p=e^{p\cdot z}\left(-p_2\partial_{z^1}+p_1\partial_{z^2}\right).
$$

The classical vector and scalar outputs are therefore

$$
[X_p,X_q]=[p,q]X_{p+q},
\qquad
[X_p,e^{q\cdot z}]=[p,q]e^{(p+q)\cdot z}.
$$

At genus $g$, the original heat-convention scalar matrix element $Q_g\vert_B$ multiplies $[p,q]^{2g+1}$, with the same output type as the corresponding classical operation. Its reduced coefficient is

$$
(-1)^g(2\pi i)^{5g}Q_g\vert_B.
$$

The genus-zero coefficient is $Q_0\vert_B=1$, and $Q_1\vert_B=0$. Heat reduction is multiplicative under composition: two vertices of genera $i,j$ have precisely the common factor $(-1)^{i+j}(2\pi i)^{5(i+j)}$. All differences below use this factor for the **total** genus.

For the Hesse cubic, the previously measured coefficients are

$$
\begin{aligned}
B&=\frac{x^3+y^3+z^3}{6}+bxyz,\\
s(B)&=b^4-b,\qquad t(B)=1-20b^3-8b^6,\\
(2\pi i)^{10}Q_2\vert_B&=\frac{s(B)}{51840},\\
-(2\pi i)^{15}Q_3\vert_B
&=t(B)\left(
\frac{1163706689}{25609494822912000}
-\frac{74897}{1976041267200}\zeta(3)\right).
\end{aligned}
$$

Here $\zeta(3)=\sum_{k=1}^{\infty}k^{-3}$. On general ternary cubics, $s(B)$ is the ratio of the measured genus-two coefficient to its value at $xyz$; $t(B)$ is the ratio of the measured genus-three coefficient to its value at the Fermat cubic. The repository proves

$$
t(B)^2-64s(B)^3
=4096\operatorname{Res}(\partial_xB,\partial_yB,\partial_zB),
$$

where the resultant vanishes exactly when the three normal gradients have a common nonzero zero. Smoothness of the projective cubic therefore implies that $s(B),t(B)$ cannot vanish simultaneously.

The local linearized differential receives no positive-genus correction on a single tangent variation, since even with descendants all momenta span at most four directions. At the normal origin, the cubic and all its first derivatives vanish. A regular primary boundary from descendant level $r$ is

$$
\partial_\Omega\left(\frac{B^{r+1}Y}{(r+1)!}\right),
$$

where $Y$ is the polyvector input and $\partial_\Omega$ is divergence with respect to the holomorphic volume form. It vanishes at the normal origin. The tangent scalar waves and divergence-free Hamiltonian waves used here are consequently closed and nonexact in this **regular local linear complex**. This observation does not establish representative independence of the isolated binary, and it does not extend the assertion to the punctured bulk with flux.

## 2. Standard rigidity in the determinant ansatz

For the proof, $\hbar$ records genus. Suppose that an ordinary Lie bracket on the Hamiltonian waves has multiplier

$$
\sum_{g\geq0}\hbar^g Q_g\vert_B\,[p,q]^{2g+1},
\qquad Q_0\vert_B=1.
$$

This is the determinant ansatz fixed by the local binary. The normalized sine or hyperbolic-sine solution and its rigidity belong to the classical Moyal deformation problem; see Guy Patissier, [*Déformations différentielles à coefficients constants et produits de Moyal généralisés*](https://www.numdam.org/item/AIHPA_1979__30_4_275_0/), Section 7, and Paul Fletcher, [*The uniqueness of the Moyal algebra*](https://doi.org/10.1016/0370-2693(90)90300-U). We give the restricted argument needed here to make its assumptions explicit.

As in the existing planar test, set

$$
a=[p,q],\qquad b'=[q,r],\qquad c=[r,p].
$$

The three cyclic factors after fusion are $b'-c,c-a,a-b'$. The Jacobi identity is therefore the formal polynomial identity

$$
\sum_{i,j\geq0}\hbar^{i+j}Q_i\vert_B Q_j\vert_B
\left[
a^{2i+1}(b'-c)^{2j+1}
+(b')^{2i+1}(c-a)^{2j+1}
+c^{2i+1}(a-b')^{2j+1}
\right]=0.
$$

These are independent determinant variables: for $a\ne0$, choose

$$
p=(1,0),\qquad q=(0,a),\qquad r=(-b'/a,-c).
$$

Polynomiality extends the resulting identity to all $a$. Differentiate once in $c$, set $c=0$, then differentiate twice in $b'$ and set $b'=0$. Since the multiplier is odd with derivative one at zero, the result is

$$
\sum_{g\geq1}(2g)(2g+1)\hbar^g Q_g\vert_B\,a^{2g-1}
=6\hbar Q_1\vert_B
\sum_{g\geq0}\hbar^g Q_g\vert_B\,a^{2g+1}.
$$

Coefficient comparison gives the entire recursion and its solution:

$$
\boxed{
(2g)(2g+1)Q_g\vert_B=6Q_1\vert_B Q_{g-1}\vert_B,
\qquad
Q_g\vert_B=\frac{6^g(Q_1\vert_B)^g}{(2g+1)!}.}
$$

Conversely these coefficients satisfy the addition identity obtained after the first differentiation, hence the full Jacobi identity. This can also be verified by multiplying their formal exponential series; no analytic convergence assumption is needed.

Since the actual local binary has $Q_1\vert_B=0$, strict Jacobi within this ansatz would force $Q_g\vert_B=0$ for every positive genus. The measured genus-two and genus-three coefficients violate this requirement on every smooth normal cubic. This conclusion concerns the isolated binary and its scalar module action. A linear invertible change of its basis conjugates the Jacobiator and cannot turn a nonzero Jacobiator into zero. Adding further operations, changing the field content, or forming the full relative bulk and brane theory is a different question.

## 3. Universal first nonzero genus and a new witness

At the first nonzero positive genus $g$, divide the Jacobiator by its reduced binary coefficient and by its output wave. Its universal value is

$$
\begin{aligned}
&a^{2g+1}(b'-c)+a(b'-c)^{2g+1}\\
&\quad +(b')^{2g+1}(c-a)+b'(c-a)^{2g+1}\\
&\quad +c^{2g+1}(a-b')+c(a-b')^{2g+1}.
\end{aligned}
$$

For genus two this polynomial is exactly

$$
10ab'c(a-b')(a-c)(b'-c).
$$

For genus three it is exactly

$$
7ab'c(a-b')(a-c)(b'-c)
\left(3a^2+3(b')^2+3c^2-2ab'-2ac-2b'c\right).
$$

The new fixture $p=(1,0),q=(0,1),r=(2,3)$ has $(a,b',c)=(1,-2,-3)$. At every genus its first-order coefficient is

$$
\boxed{2\left(1+2^{2g+2}+4^{2g+1}-3^{2g+2}\right).}
$$

It is strictly positive for $g\geq2$: the ratio $4^{2g+1}/3^{2g+2}$ exceeds one at $g=2$, where it equals $1024/729$, and increases by $16/9$ when $g$ increases by one. The old fixture $r=(1,2)$ has coefficient

$$
3^{2g+1}-2^{2g+3}+5.
$$

The corresponding exact integer values are:

| Genus | $r=(1,2)$ | $r=(2,3)$ |
|---:|---:|---:|
| 2 | 120 | 720 |
| 3 | 1680 | 20160 |
| 4 | 17640 | 408240 |
| 5 | 168960 | 7333920 |
| 6 | 1561560 | 124684560 |

Consequently the new fixture gives the first nonzero reduced Jacobi matrix element

$$
\frac{s(B)}{72}
$$

when $s(B)\ne0$, and

$$
t(B)\left(
\frac{1163706689}{1270312243200}
-\frac{74897}{98017920}\zeta(3)\right)
$$

when $s(B)=0$ and $t(B)\ne0$.

The nonvanishing of the Fermat coefficient has a short rational proof. Since $x^{-3}$ is strictly decreasing,

$$
\zeta(3)>
\sum_{k=1}^{8}\frac1{k^3}+\int_9^\infty\frac{dx}{x^3}
=\frac{2136104771}{1778112000}.
$$

The difference from the exact zero threshold is

$$
\frac{2136104771}{1778112000}
-\frac{1163706689}{970665120}
=\frac{327281302787}{133175254464000}>0.
$$

It follows that

$$
\frac{1163706689}{25609494822912000}
-\frac{74897}{1976041267200}\zeta(3)
< -\frac{327281302787}{3513622689703526400000}<0.
$$

This proof needs no transcendence or irrationality theorem.

## 4. Complete compositions and elimination of unknown periods

At total genus $g$, the complete heat-convention coefficient is the sum over **all** ordered genus splits:

$$
\sum_{i=0}^{g}Q_i\vert_B Q_{g-i}\vert_B
\left[
a^{2i+1}(b'-c)^{2g-2i+1}
+(b')^{2i+1}(c-a)^{2g-2i+1}
+c^{2i+1}(a-b')^{2g-2i+1}
\right].
$$

The numerical coefficient in each summand simplifies to

$$
1+3^{2g-2i+1}-2^{2g+2}
\quad\text{at }r=(1,2),
$$

and to

$$
1+2^{4g-2i+3}-3^{2g+2}
\quad\text{at }r=(2,3).
$$

These formulas hold for every genus and include the classical vertex on either side. The following are the complete unreduced coefficients through genus six, using $Q_1\vert_B=0$:

| Genus | $r=(1,2)$ | $r=(2,3)$ |
|---:|---|---|
| 4 | $17640Q_4\vert_B-780(Q_2\vert_B)^2$ | $408240Q_4\vert_B-26280(Q_2\vert_B)^2$ |
| 5 | $168960Q_5\vert_B-5760Q_2\vert_B Q_3\vert_B$ | $7333920Q_5\vert_B-407520Q_2\vert_B Q_3\vert_B$ |
| 6 | $1561560Q_6\vert_B-12840Q_2\vert_B Q_4\vert_B-14196(Q_3\vert_B)^2$ | $124684560Q_6\vert_B-653040Q_2\vert_B Q_4\vert_B-2685816(Q_3\vert_B)^2$ |

Each row is multiplied by $(-1)^g(2\pi i)^{5g}$ when reporting a reduced Ward coefficient. The table is derived directly from the full cyclic sum; no part of a given matrix element has been discarded.

### 4.1 Genus four

Seven times the second matrix element minus 162 times the first gives

$$
\begin{aligned}
&(2\pi i)^{20}\left[
7\left(408240Q_4\vert_B-26280(Q_2\vert_B)^2\right)
-162\left(17640Q_4\vert_B-780(Q_2\vert_B)^2\right)
\right]\\
&\qquad=-57600(2\pi i)^{20}(Q_2\vert_B)^2
=\boxed{-\frac{s(B)^2}{46656}}.
\end{aligned}
$$

This fixes an exact target for the same difference of the remaining Ward contributions, with the opposite sign. It does not determine the individual connected genus-four binary period.

### 4.2 Genus five

Thirty-two times the second matrix element minus 1389 times the first gives

$$
\begin{aligned}
&-(2\pi i)^{25}\left[
32\left(7333920Q_5\vert_B-407520Q_2\vert_B Q_3\vert_B\right)
-1389\left(168960Q_5\vert_B-5760Q_2\vert_B Q_3\vert_B\right)
\right]\\
&\qquad=5040000(2\pi i)^{25}Q_2\vert_B Q_3\vert_B\\
&\qquad=\boxed{s(B)t(B)\left(
-\frac{1163706689}{263411946749952}
+\frac{374485}{101624979456}\zeta(3)\right)}.
\end{aligned}
$$

The result is nonzero whenever $s(B)t(B)\ne0$. The unknown genus-five binary period cancels before any lower-genus value is substituted.

### 4.3 Genus six on the equianharmonic locus

Thirteen times the second matrix element minus 1038 times the first gives, for every cubic,

$$
(2\pi i)^{30}
\left(4838400Q_2\vert_B Q_4\vert_B-20180160(Q_3\vert_B)^2\right).
$$

On $s(B)=0$, the first term vanishes without any assumption about $Q_4\vert_B$. The full result is then

$$
\boxed{
-20180160\,t(B)^2
\left(
\frac{1163706689}{25609494822912000}
-\frac{74897}{1976041267200}\zeta(3)\right)^2.}
$$

With every coefficient expanded, it is

$$
t(B)^2\left(
-\frac{1354213258023342721}{32499555260451677798400000}
+\frac{87158139886033}{1253840866529771520000}\zeta(3)
-\frac{5609560609}{193493960884224000}\zeta(3)^2
\right).
$$

It is nonzero for $t(B)\ne0$. At the Fermat cubic $t(B)=1$, it is strictly negative by the rational bound above. Thus the genus-four target covers $s(B)\ne0$, and the genus-six target covers the remaining smooth locus. Neither target requires a period beyond the already measured genera two and three.

## 5. What the full physical test requires

The regular local primary binary is a genuine subcalculation of the twisted bulk theory. Its nonzero Jacobiator is not by itself an anomaly of that theory: the quantum master equation includes loop contractions and the rest of the effective vertices. In the relative theory it additionally includes brane and flux terms. The original genus-two and genus-three graph calculations already exhibit compensating contributions for the old fixture. The new results specify exact higher-genus comparisons that any full calculation must meet in the same convention.

For example, the sum of the non-binary Ward terms at total genus four must have the difference

$$
7\left[\text{coefficient at }r=(2,3)\right]
-162\left[\text{coefficient at }r=(1,2)\right]
=\frac{s(B)^2}{46656},
$$

provided the complete quantum Ward identity is evaluated on these same regular inputs with the same heat reduction. This is a necessary condition, not a computation of the individual non-binary terms.

The cohomological distinction matters for holography. In the existing punctured, flux-dependent construction, smooth cubic backgrounds admit homotopies that make the relevant scalar loop outputs exact; singular critical directions behave differently. Closure and nonexactness in the regular local complex consequently cannot be transported to the D3 relative charge complex without an additional argument. The present algebraic targets remain well-defined on the specified cochains and can be used when that transfer is calculated.

## 6. Independent exact certificate

Run from the repository root:

```bash
python explore/iib-bcov/developments_2026_11/planar_ward_rigidity.py
python explore/iib-bcov/developments_2026_11/planar_ward_rigidity.py --check
```

The [certificate program](planar_ward_rigidity.py) uses Python and SymPy, imports no existing repository evaluator, and performs no floating-point arithmetic. Its [complete output](planar_ward_rigidity.json) records every exact coefficient and the pinned source blob hashes.

| Verification | Exact scope |
|---|---|
| Universal factorizations | Full polynomials in three independent determinant variables at genera two and three |
| Independent nested differential evaluation | Five integer-momentum fixtures; every ordered genus split at total genera zero through six; 140 complete cyclic comparisons |
| Closed fixture formulas | Both fixtures at genera zero through twelve |
| Complete compositions and unknown-coefficient cancellation | Full genus-four, genus-five and genus-six expressions before taking their stated differences |
| Moyal positive control | Complete Jacobi polynomial vanishes at every genus zero through eight |
| Deliberately perturbed control | Adding one to the Moyal genus-two coefficient gives exactly 720 at the new fixture |
| Fermat nonvanishing | Exact rational integral-test bound |
| Earlier graph normalization | Arithmetic checks of the existing genus-two and genus-three sums; no reintegration |

The independent differential evaluator distributes every outer derivative over the two factors of the inner vertex using the Leibniz rule. It never substitutes a fused momentum or a determinant in that calculation. Its results agree with the determinant formula for each ordered genus split. The finite checks verify the displayed coefficients; the proof in Section 2 establishes the all-order rigidity statement.

The inherited genus-two normalization control is

$$
\frac{41}{90720}+\frac{11}{11340}+\frac1{2016}+\frac1{2520}
=\frac1{432}=\frac{120}{51840}.
$$

The inherited genus-three control is

$$
\begin{aligned}
&-\frac{3361}{27433728000}
+\frac{13904672281}{177843714048000}
-\frac{74897}{1176215040}\zeta(3)
-\frac1{580608}\\
&\qquad=1680\left(
\frac{1163706689}{25609494822912000}
-\frac{74897}{1976041267200}\zeta(3)\right).
\end{aligned}
$$

The three terms correspond respectively to the double-square graph, the tetrahedral graph, and the hexagon with a genus-two insertion. These are prior results retained solely as convention checks.

## References

1. Si Li and Minghao Wang, *On the exact quantization of BCOV Theory on Calabi-Yau Manifolds*, [arXiv:2610.21927v1](https://arxiv.org/abs/2610.21927), 15 October 2026. [Full HTML](https://arxiv.org/html/2610.21927v1). The announcement supplies the all-genus local interaction and states that a detailed account of the proofs will follow.
2. Guy Patissier, *Déformations différentielles à coefficients constants et produits de Moyal généralisés*, Annales de l'Institut Henri Poincaré, Section A **30** (1979), 275–293, especially Section 7. [Primary publication](https://www.numdam.org/item/AIHPA_1979__30_4_275_0/).
3. Paul Fletcher, *The uniqueness of the Moyal algebra*, Physics Letters B **248** (1990), 323–328. [DOI:10.1016/0370-2693(90)90300-U](https://doi.org/10.1016/0370-2693(90)90300-U).
4. Repository source at commit [`aa0d2f1725dff135494eef132f5e92a8c0401a42`](https://github.com/lbotao399-stack/holomorphic-topological-feynman-diagrams/tree/aa0d2f1725dff135494eef132f5e92a8c0401a42): `explore/iib-bcov/README.md`, `explore/iib-bcov/planar_cubic_primary.py`, and `explore/iib-bcov/operator/README.md`. Exact Git blob hashes are recorded in [the certificate](planar_ward_rigidity.json).
