# A finite descendant polarization with first nonzero scalar response at genus four

## Result and scope

The finite field profile and coefficient extraction below satisfy

$$
Q_0\phi(0)=Q_1\phi(0)=Q_2\phi(0)=Q_3\phi(0)=0,
\qquad
(2\pi i)^{20}Q_4\phi(0)=\frac{31}{64012032000}
\tag{1}
$$

**Every equality in (1) refers to the specified polarized coefficient.** The full, unpolarized response of the field profile is not asserted to vanish at lower genus. All nonzero scalar field slots have descendant levels between zero and four, so the experiment lies in the minimal fivefold field space used by this repository.

This calculation is a local elementary-field test in the fixed flat BCOV coordinates. It gives a new genus-four benchmark for twisted supergravity and for a future bulk-to-boundary comparison. It does not establish a transferred observable, invariance under general field redefinitions, or a boundary correlation function. No D3 flux is present; no power of its flux parameter is inferred.

The worldsheet input is the flat vertex prescription in Si Li and Minghao Wang, [*On the exact quantization of BCOV Theory on Calabi–Yau Manifolds*, arXiv:2610.21927v1](https://arxiv.org/html/2610.21927v1), equations (1.9)–(1.11), Lemmas 1.2 and 1.6, and Theorem 1.11. This source announces the quantization theorem and states that detailed proofs will appear subsequently. The result here is conditional on that prescription and on the repository's complete genus-four planar periods, listed explicitly below. The new calculation does not recompute those master-period integrals.

## 1. Exact field profile and the measured coefficient

Use the existing fivefold coordinates $z^1,\ldots,z^5$, with dual covectors $e^1,\ldots,e^5$. The three labeled hard scalar waves have momenta

$$
p_1=e^1-4e^3-4e^4-4e^5,
\qquad p_2=e^2,
\qquad p_3=-e^1-e^2.
\tag{2}
$$

Their independent formal amplitudes are $\alpha_1,\alpha_2,\alpha_3$. The variables $t_1,t_2,t_3$ only record the descendant level on these three waves. For each scalar descendant slot $A_k$, set

$$
A_k^{\mathrm{scalar}}(z)
=\sum_{j=1}^3\alpha_j t_j^k e^{p_j\cdot z},
\qquad 0\leq k\leq4.
\tag{3}
$$

Add $e^{z^3}+e^{z^4}+e^{z^5}$ to the primary scalar $A_0^{\mathrm{scalar}}$, and set the primary vector to

$$
A_0^{\mathrm{vector}}(z)=\lambda e^{z^4}\partial_{z^3}.
\tag{4}
$$

Here $\lambda$ is its independent formal amplitude. All remaining slots are zero. The vector is divergence free since $e^4(\partial_{z^3})=0$. These are finite descendant profiles. The scalar cutoff agrees with `elementary_one_loop` in [the existing elliptic evaluator](../../iib-bcov/elliptic_descendants.py).

In $Q_g\phi(0)$, take the coefficient of $\alpha_1\alpha_2\alpha_3\lambda$, apply

$$
\frac1{12}
\left(
\partial_{t_1}^2+\partial_{t_2}^2+\partial_{t_3}^2
-\partial_{t_1}\partial_{t_2}
-\partial_{t_1}\partial_{t_3}
-\partial_{t_2}\partial_{t_3}
\right)^2,
\tag{5}
$$

and then set $t_1=t_2=t_3=0$. This completely defines the measurement in (1). In particular, only total descendant degree four enters. Formula (5) weights each **distinct permutation** of the three levels as follows.

| Descendant levels | Weight per distinct permutation | Number of permutations |
|---|---:|---:|
| $(4,0,0)$ | $2$ | $3$ |
| $(3,1,0)$ | $-1$ | $6$ |
| $(2,2,0)$ | $1$ | $3$ |
| $(2,1,1)$ | $0$ | $3$ |

No scalar level-five field occurs. If descendant fields are written using the alternate convention $t=-\zeta$, the total degree four contributes $(-1)^4=1$.

## 2. All-degree cancellation of fused primary marks

The following algebraic identity explains why (5) isolates a late response without introducing genus-dependent source choices.

For a homogeneous polynomial $P(t_1,t_2,t_3)$ of total degree $K$, write its coefficients in divided-power form:

$$
P(t_1,t_2,t_3)
=\sum_{k_1+k_2+k_3=K}
C_{k_1,k_2,k_3}
\frac{t_1^{k_1}t_2^{k_2}t_3^{k_3}}{k_1!k_2!k_3!}.
\tag{6}
$$

Thus the numbers $C_{k_1,k_2,k_3}$ are precisely the weights with which $P(\partial_{t_1},\partial_{t_2},\partial_{t_3})$ evaluates the ordinary monomials in (3). Coefficients outside the stated degree or allowed field range are zero. Suppose

$$
(\partial_{t_1}+\partial_{t_2}+\partial_{t_3})P=0.
\tag{7}
$$

For any residual levels of total degree $K-1$, this is the finite identity

$$
C_{\ell_1+1,\ell_2,\ell_3}
+C_{\ell_1,\ell_2+1,\ell_3}
+C_{\ell_1,\ell_2,\ell_3+1}=0.
\tag{8}
$$

If a fixed set of $d$ distinct primary marks is fused into the three hard marks, its complete labeled assignment coefficient, at any residual levels of total degree $K-d$, is

$$
\sum_{a_1+a_2+a_3=d}
\frac{d!}{a_1!a_2!a_3!}
C_{\ell_1+a_1,\ell_2+a_2,\ell_3+a_3}=0,
\qquad 1\leq d\leq K.
\tag{9}
$$

Equation (9) is the corresponding coefficient of
$(\partial_{t_1}+\partial_{t_2}+\partial_{t_3})^dP$, so it follows from (7) at **every degree**. It includes the exact multinomial for distinct donor labels; no unlabeled replacement is made.

To apply it to the scalar vertex, let all forgettable primary momenta lie in the span of $e^3,e^4,e^5$, and let the three retained planar momenta be $e^1,e^2,-e^1-e^2$. For any fixed set of fused donors, two assignments differ by transverse shifts on the three hard momenta whose sum is zero. A determinant-one shear fixing $e^3,e^4,e^5$ realizes any such triple of shifts: choose its images on $e^1,e^2$, and its image on $-e^1-e^2$ is then forced. The scalar vertex is unchanged by this shear. Therefore, at fixed residual descendant levels, all assignments of this donor set have the same scalar vertex.

The cotangent pullback formula for forgetting primary marks expands into ordinary descendant vertices, with a minus sign and one lost cotangent power for each fused donor. Equations (8)–(9) cancel every nonempty fusion layer. Consequently the ordinary weighted scalar vertex is exactly the vertex with all three cotangent classes pulled back from the retained three-mark moduli space. This extends the two-pair cancellation previously used in [the genus-three experiment](../README.md#a-fixed-elementary-scalar-polarization-whose-first-response-is-at-three-loops).

The absorbed root and the primary vector can be included among the forgettable marks of the **symmetric scalar trace**. Their roles in the elementary-field operation do not change: the root–vector pairing is held outside this identity. The fused intermediate traces are only terms of the cotangent identity, not new physical source insertions.

For (5), the divided-power polynomial is

$$
P(t_1,t_2,t_3)
=\frac1{12}
\left(t_1^2+t_2^2+t_3^2-t_1t_2-t_1t_3-t_2t_3\right)^2.
\tag{10}
$$

It obeys (7), so the cancellation is exact for all genera. The certificate independently sums every labeled assignment for all $1\leq d\leq4$: 20 residual-level tests and 246 individual assignments.

## 3. Why genera zero through three vanish

For target dimension five the worldsheet degree condition is

$$
\sum_i k_i=-2g+n-3.
\tag{11}
$$

At total descendant degree four, $n=2g+7$. Retaining the three hard marks therefore forgets $2g+4$ primary marks, including the root and vector. At fixed curve and hard marks, the forgotten momenta enter only the three transverse Abel maps. Each Jacobian has complex dimension $g$, so their combined vertical rank is at most $3g$.

After the cancellation proved above, all cotangent classes are pulled back from the retained moduli space. If

$$
2g+4>3g,
\tag{12}
$$

there are more forgotten complex marking directions than the Betti form can saturate. Its required vertical top-degree coefficient is zero. For $g=1,2,3$, the deficits are respectively $3,2,1$. This proves the three positive-genus vanishings in (1) for every ordinary-primary Fourier multiplicity, before their sum. The argument uses the Betti-form rank property and the cotangent pullback prescription, with the convergence and extension statements of Li–Wang justifying integration.

At genus zero the complete classical check is

$$
3(2)\frac{4!}{4!}
+6(-1)\frac{4!}{3!1!}
+3(1)\frac{4!}{2!2!}
=6-24+18=0.
\tag{13}
$$

The free differential also vanishes on the chosen holomorphic scalar waves and divergence-free primary vector.

As an independent positive-genus control, the certificate evaluates ordinary genus-one stable partitions at primary scalar multiplicities $(1,1,2)$. Their balancing root is $3e^3+2e^4+2e^5$, and the complete fivefold input has full momentum rank. The three ordinary positive scalar vertices are

$$
V_{1;(4,0,0)}=\frac73,
\qquad V_{1;(3,1,0)}=\frac{28}3,
\qquad V_{1;(2,2,0)}=14.
\tag{14}
$$

Every distinct level permutation is evaluated separately at the fixed, transversely shifted hard momenta: all twelve weighted tuples reproduce the values in (14), with $12,48,72$ nonzero stable partitions for each tuple in its respective orbit. This check does not assume the shear or forgetful identity. Each vertex is nonzero, whereas

$$
6\frac73-6\frac{28}3+3(14)=14-56+42=0.
\tag{15}
$$

This independently checks genuine descendant cancellation on a full-rank input, not a vanishing caused by deficient total momentum rank.

## 4. Complete genus-four coefficient and normalization

At genus four the twelve forgotten directions exactly saturate the three transverse Jacobians. The planar retained vertex uses the unit triangle $(e^1,e^2,-e^1-e^2)$. The inherited complete periods are

| Levels | Genus-four planar vertex |
|---|---:|
| $(4,0,0)$ | $731/1024192512000$ |
| $(3,1,0)$ | $1493/438939648000$ |
| $(2,2,0)$ | $157/29262643200$ |
| $(2,1,1)$ | $433/37623398400$ |

These are the constants `LEVELS` and `EXPECTED[4]` in [the existing complete planar evaluator](../code/planar_four_loop.py); its derivation and original graph checks are described in [the genus-four planar section](../README.md#three-point-genus-four-vertex-in-two-target-dimensions). Mark permutations preserve the eighth power of the planar area determinant, so the weights in Section 1 give

$$
\begin{aligned}
6\frac{731}{1024192512000}
-6\frac{1493}{438939648000}
+3\frac{157}{29262643200}
&=\frac{4386-20902+16485}{1024192512000}\\
&=-\frac{31}{1024192512000}.
\end{aligned}
\tag{16}
$$

For four unit transverse momenta on each of the three axes, the ordered degree-four Abel maps have product degree

$$
(4!)^3=13824.
\tag{17}
$$

Hence the complete **labeled scalar vertex** is

$$
13824\left(-\frac{31}{1024192512000}\right)
=-\frac{31}{74088000}.
\tag{18}
$$

The fixed source in (3)–(4) requires the ordinary-primary Fourier multiplicities and their Taylor factorials. Let $a,b,c$ count primary scalar copies of $e^{z^3},e^{z^4},e^{z^5}$. At genus four, $a+b+c=10$. The vector contributes another momentum $e^4$, and momentum conservation gives the absorbed root

$$
p_{\mathrm{root}}=(4-a)e^3+(3-b)e^4+(4-c)e^5.
\tag{19}
$$

For twelve arbitrary transverse momenta, the fiber degree is the sum over their ordered partitions into four triples of the product of four squared three-by-three determinants. To derive it, expand the vertical Abel determinant and integrate the paired evaluations of normalized holomorphic one-forms at every labeled point. The bilinear relations identify the two form labels. Each of the four form labels occurs in a triple, giving exactly the stated ordered partitions. Equivalently it is the coefficient of the product of twelve independent weights in the fourth power of the determinant of the weighted sum of their outer products. This also proves that transverse recoil at the hard marks does not change the saturated coefficient.

Enumeration of all $\binom{12}{2}=66$ choices in (19) gives exactly the following nonzero fiber degrees.

| $(a,b,c)$ | Root | Fiber degree | $p_{\mathrm{root}}(\partial_{z^3})$ |
|---|---|---:|---:|
| $(3,3,4)$ | $e^3$ | $13824$ | $1$ |
| $(4,2,4)$ | $e^4$ | $13824$ | $0$ |
| $(4,3,3)$ | $e^5$ | $13824$ | $0$ |

Indeed, a nonzero triple must receive one direction from each transverse axis. Apart from the one root, all vectors are axis vectors. A surviving partition therefore has four copies of two axes and three of the third before the root is included, which gives exactly these rows. Thus the full root-weighted, factorial-divided sum is

$$
\frac{13824}{3!3!4!}=16.
\tag{20}
$$

There is no extra $3!$ from the hard waves: their amplitudes are separately labeled. There is no factorial for the single vector or the absorbed root. The repeated ordinary-primary waves carry exactly $a!b!c!$.

Finally the repository's elementary scalar convention is

$$
Q_g\phi(0)
=-\frac{(-1)^g}{(2\pi i)^{5g}}
\,p_{\mathrm{root}}(\partial_{z^3})\,
V_{g;\mathbf k}
\tag{21}
$$

on a labeled scalar trace with the single primary vector. The sign comes from the root divergence and the heat normalization; it is the same convention used by the existing elementary-field experiments in [the BCOV notes](../../iib-bcov/README.md). At $g=4$, its sign is minus, while (16) is negative. Combining all factors gives

$$
\begin{aligned}
(2\pi i)^{20}Q_4\phi(0)
&=-(-1)^4\,
\frac{(4!)^3}{3!3!4!}\,
\left(-\frac{31}{1024192512000}\right)\\
&=\frac{31}{64012032000},
\end{aligned}
\tag{22}
$$

for the measurement (5). Before the repeated-source factorials, the labeled elementary response is $+31/74088000$.

## 5. Classification within the minimal scalar field space

Require that the polarization be symmetric in the three hard marks, obey (7), and use only individual levels $0\leq k_i\leq4$. The space is one dimensional at each total degree

$$
K=0,2,3,4,6,
\tag{23}
$$

and zero at every other degree. Polynomial representatives are

| Total degree | Representative, with irrelevant overall scale suppressed |
|---|---|
| $0$ | $1$ |
| $2$ | $t_1^2+t_2^2+t_3^2-t_1t_2-t_1t_3-t_2t_3$ |
| $3$ | $(2t_1-t_2-t_3)(2t_2-t_3-t_1)(2t_3-t_1-t_2)$ |
| $4$ | $(t_1^2+t_2^2+t_3^2-t_1t_2-t_1t_3-t_2t_3)^2$ |
| $6$ | $(t_1-t_2)^2(t_2-t_3)^2(t_3-t_1)^2$ |

Here is a finite algebraic proof, independent of genus. Declare distinct monomials orthogonal and give the monomial of levels $(k_1,k_2,k_3)$ the positive squared norm

$$
\prod_{i=1}^3 k_i!(4-k_i)!.
\tag{24}
$$

On the space with individual degrees at most four, the adjoint of $\sum_i\partial_{t_i}$ is $\sum_i(4t_i-t_i^2\partial_{t_i})$, and their commutator is $12-2\sum_i t_i\partial_{t_i}$. A homogeneous polynomial satisfying (7) therefore obeys

$$
\left\|\sum_i(4t_i-t_i^2\partial_{t_i})P\right\|^2
=(12-2K)\|P\|^2.
\tag{25}
$$

Positivity gives $K\leq6$. Symmetric polynomials invariant under common translation are generated by the degree-two and degree-three representatives in the table: center the three variables by their average, then use the two nontrivial elementary symmetric polynomials. Up to degree six this gives one candidate in degrees zero, two, three, four and five, and two candidates in degree six. The degree-five product has individual degree five, so it is inadmissible. In degree six the unique admissible linear combination is fixed by

$$
\begin{aligned}
&\left[(2t_1-t_2-t_3)(2t_2-t_3-t_1)(2t_3-t_1-t_2)\right]^2\\
&\quad-4\left(t_1^2+t_2^2+t_3^2-t_1t_2-t_1t_3-t_2t_3\right)^3\\
&=-27(t_1-t_2)^2(t_2-t_3)^2(t_3-t_1)^2.
\end{aligned}
\tag{26}
$$

The last polynomial has individual degree four. This proves (23). Exact rational row reduction independently verifies every total degree from zero through twelve, the full range permitted by three level-four slots.

For the degree-six representative, divided-power coefficients can be normalized to the following integers per distinct permutation:

| Levels | Weight |
|---|---:|
| $(4,2,0)$ | $2$ |
| $(4,1,1)$ | $-2$ |
| $(3,3,0)$ | $-3$ |
| $(3,2,1)$ | $1$ |
| $(2,2,2)$ | $-2$ |

The same theorem cancels all nonempty donor sets. At total degree $K$, the number of forgotten marks is $2g+K$, so the selected positive-genus response vanishes whenever $g<K$. At genus zero, applying the polarization to the classical multinomial gives a derivative of the common-translation invariant polynomial and is zero for $K>0$. In particular the degree-six polarization vanishes through genus five. **A nonzero genus-six coefficient is not established here.** The theorem identifies its first possible genus, and the evaluation at that genus remains an explicit next calculation.

## 6. Reproduction and audit trail

Run from the repository root:

```sh
python explore/twisted-holography-N4/developments_2026_11/delayed_bulk_response.py
```

The certificate uses only the Python standard library and exact integers and rational numbers. Its recorded output is [delayed_bulk_response.json](delayed_bulk_response.json). It performs these separate checks:

1. Exact kernel dimensions for every total level zero through twelve, plus independent construction of each polynomial representative and its translation derivative.
2. Explicit labeled donor sums at total degrees two, three, four and six: respectively 18, 72, 246 and 2412 individual assignments, covering all nonempty donor counts and all residual levels.
3. Classical multinomial cancellation and an independently implemented genus-one stable-partition calculation of all twelve distinct weighted level tuples, with nonzero ordinary controls (14).
4. Exact recombination (16) of the inherited complete genus-four periods, read through a restricted syntax tree without executing their graph engine.
5. All 66 primary Fourier multiplicities, with determinant partition evaluation of every transverse fiber degree.
6. Independent comparison of the determinant-partition algorithm with finite differences of the Gram determinant on twelve generic integer vectors. Both give the exact coefficient $2671280640$.
7. Repeated-source factorials, root–vector pairing, and the elementary heat sign in (22).

The inherited period file is recorded at repository revision `aa0d2f1725dff135494eef132f5e92a8c0401a42`; its SHA-256 is `a17bd30096fa14d7c61e5a0e9a57c0675e42ae37f45f360f5f15ae784719481e`. The graph master-period integrations themselves are an explicit input, not an independent check performed by this certificate.

The new results are the admissible fixed profile, the all-degree donor-cancellation criterion, its complete symmetric minimal-field classification, and the nonzero genus-four elementary response with every multiplicity and sign fixed. A holographic comparison still requires the corresponding state continuation, the full punctured-background transfer and its Ward terms, and a boundary computation in matching conventions.
