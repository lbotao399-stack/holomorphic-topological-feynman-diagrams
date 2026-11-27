# A uniform wheel-period stabilization theorem

For every wheel with $n\ge3$ spokes, the rational span of **all absolutely
convergent positive-integer Schwinger monomials** is

$$
\boxed{
\operatorname{span}_{\mathbb Q}
\{I_q(W_n;\alpha):\alpha\text{ convergent}\}
=\mathbb Q+\sum_{j=2}^{n-1}\mathbb Q\zeta(2j-1),
\qquad q\ge\max(3,n-1).
}
$$

This improves the preceding sufficient threshold $q\ge n+2$ in the
[all-power wheel proof](../../nonSP/verify/full_ward/README.md#arbitrary-integer-powers-on-an-n-spoke-wheel).
The proof below supplies a rational direction at $q=n-1$ and retains every
odd-zeta direction. It does not assume arithmetic independence of zeta values,
and does not infer uniform closure from a finite table or a rank plateau.
The threshold is sufficient; optimality is not asserted.

```sh
python3 explore/circuit-general/core-period-closure/verify.py
```

The [certificate](certificate.json) records all finite base values, exact graph
reductions, coefficient determinants, the uniform inequality, and dependency
hashes. `--write` regenerates it. The executable replays every recorded check
with rational arithmetic and the existing exact residue algorithm.

## 1. Normalization, powers, and convergence

The hub is vertex $0$ and the cyclic rim vertices are $1,\ldots,n$. The
ordered edges are the spokes $0\to1,\ldots,0\to n$, followed by
$1\to2,\ldots,(n-1)\to n,n\to1$. Scalar periods do not depend on these
orientations. The projective measure is positive, and the convention is

$$
I_q(G;\alpha)=\int_{t_e>0,\ \sum_e t_e=1}
 \frac{\prod_e t_e^{\alpha_e-1}}{U_G(t)^q}\,
 dt_1\cdots dt_{|E|-1},\qquad \sum_e\alpha_e=nq.
$$

There is no division by $\prod_e\Gamma(\alpha_e)$. Kirchhoff's polynomial
$U_G$ is the sum of the spanning-tree-complement monomials. Every power is a
positive integer. Writing $h(S)$ for the cycle rank of an edge subset,
absolute convergence is exactly

$$
\sum_{e\in S}\alpha_e>q\,h(S)
\quad(\varnothing\ne S\subsetneq E).
$$

The [existing wheel residue theorem](../../nonSP/verify/full_ward/README.md#arbitrary-integer-powers-on-an-n-spoke-wheel)
already proves the upper inclusion in the boxed statement, at every integer
$q$ for which the monomial converges. Only the reverse inclusion needs to be
established here.

Multiplication by $U_G$ gives

$$
I_q(G;\alpha)=\sum_{T\text{ spanning tree}}
 I_{q+1}(G;\alpha+\mathbf1_{E\setminus T}).
$$

Every term converges: $|(E\setminus T)\cap S|\ge h(S)$ preserves the strict
inequality. Thus the full scalar spans are nested in $q$.

## 2. Exact graph operations in this normalization

At a trivalent vertex, let the incident Schwinger powers be $a_1,a_2,a_3$.
If $a_1+a_2+a_3=q$, remove the vertex and insert the triangle with the edge
opposite $a_i$ carrying power $q-a_i$. The value of $I_q$ is unchanged.

To see the normalization directly, write the three new resistances as
$x_1,x_2,x_3$ and set

$$
r_i=\frac{x_jx_k}{x_1+x_2+x_3},\qquad
U_{\triangle}=(x_1+x_2+x_3)U_{\star},\qquad
\left|\frac{\partial(r_1,r_2,r_3)}
{\partial(x_1,x_2,x_3)}\right|
=\frac{x_1x_2x_3}{(x_1+x_2+x_3)^3}.
$$

The transformed density is exactly

$$
\frac{\prod_i x_i^{q-a_i-1}}{U_{\triangle}^{q}},
$$

because the remaining power of $x_1+x_2+x_3$ is $q-\sum_i a_i=0$.
This positive, homogeneous change of variables respects the positive
projective measure. Its inverse is
$x_i=(r_1r_2+r_2r_3+r_3r_1)/r_i$.

If two parallel edges have powers $a,b$, combine them into power $a+b-q$.
The exact prefactor is

$$
\frac{\Gamma(q-a)\Gamma(q-b)}{\Gamma(2q-a-b)}.
$$

All uses below satisfy $0<a,b<q$ and $a+b>q$, so these are ordinary
convergent beta integrals. The checker implements these two operations
directly, independently of the wheel residue evaluation.

## 3. An explicit rational wheel at $q=n$

First consider $q=n$. Give the spokes at rim vertices $1,n$ power $q-1$,
every other spoke power $q-2$, the rim path $1\!-!2\!-!\cdots\!-!n$
powers $1$, and the closing rim edge $n\!-!1$ power $q-1$.

Eliminate vertices $2,3,\ldots,n-1$ in order. Every eliminated star has
total power $q$. At step $k$, the two nontrivial parallel factors are

$$
\frac1k,\qquad
\begin{cases}
1/[k(k+1)],&1\le k\le q-3,\\
1/(q-2),&k=q-2.
\end{cases}
$$

The last rim merge has factor one. Their product is $1/((q-2)!)^3$.
The remaining one-loop triangle has powers $(1,1,q-2)$ and value
$\Gamma(q-2)/\Gamma(q)$. Consequently this complete wheel period is

$$
\boxed{\frac{1}{(n-1)(n-2)((n-2)!)^3}}.
$$

All transformations are positive and the final triangle converges, proving
convergence of the starting monomial as well. At $n=3$ this gives $1/2$ at
$q=3$; together with the $q=2$ radial value $6\zeta(3)$ lifted once by $U$,
it proves the theorem for $W_3$.

## 4. A four-vertex period with a nonzero rational coefficient

For each integer $q\ge3$, consider $K_4=W_3$ with spoke powers
$(q-1,1,1)$ and cyclic rim powers $(1,q-1,q-1)$. Its three heavy edges form
a path, and its three light edges form the complementary path. All triangle
sums exceed $q$, all four-cycle sums exceed $q$, and deleting one edge
leaves power sum greater than $2q$. Thus it is absolutely convergent.

Denote its rational coefficient by $R_q$. The exact residue formula gives

$$
I_q(K_4;q-1,1,1;1,q-1,q-1)
=R_q+\frac{12((q-2)!)^2}{(2q-2)!}\zeta(3).
\tag{1}
$$

### The coefficient of $\zeta(3)$ for every $q$

For this calculation put $r=q-2$. The occupied residue multiplicities are
$(2,2,\ldots,2,3)$ at $0,1,\ldots,r$, and the factorial prefactor is
$(r!)^4$. The univariate function before convolution is

$$
\frac{1}{(w-r)^3\prod_{j=0}^{r-1}(w-j)^2}.
$$

Its highest principal coefficients are $1/(r!)^2$ at $w=r$ and
$-1/((j!)^2((r-j)!)^2(r-j))$ at $w=j<r$.
The spherical multiplicity is

$$
\frac{2s^2\prod_{j=1}^{r}(s^2-j^2)}{(2r+2)!}.
$$

Only the order-five pole at shift zero and the order-four poles at shifts
$\pm1,\ldots,\pm r$ can contribute to $\zeta(3)$. The multiplicity has a
double zero at zero and simple zeros at those other shifts. Substitution
of the displayed principal coefficients gives, including $(r!)^4$,

$$
\frac{12(r!)^2}{(2r+2)!}
\sum_{k=0}^{r}(-1)^{r-k}\binom rk\binom{r+k}{k}
=\frac{12(r!)^2}{(2r+2)!}.
$$

The binomial sum is one: it is the coefficient of $x^r$ in
$(1+x)^r((1+x)-1)^r=x^r(1+x)^r$.
The existing wheel theorem bounds the weight by three, proving (1) with
$R_q\in\mathbb Q$.

### Exact finite base range

The residue algorithm calculates the following rational coefficients in
full. The checker verifies all entries and the coefficient in (1).

| $q$ | $R_q$ |
|---:|---:|
| 3 | $-1/6$ |
| 4 | $7/360$ |
| 5 | $1459/60480$ |
| 6 | $41597/2721600$ |
| 7 | $16605397/1796256000$ |
| 8 | $52657567/9081072000$ |
| 9 | $58321119121/15256200960000$ |
| 10 | $1540257744143/583549686720000$ |
| 11 | $140089601530061/73916293651200000$ |
| 12 | $480053046971681/341493276668544000$ |
| 13 | $2221011870323291011/2073547175931399168000$ |
| 14 | $7813118641826693161/9359761558023676800000$ |
| 15 | $60998079785613798244717/91987736592256695590400000$ |
| 16 | $31234328161374002908531/58330416060687554784000000$ |
| 17 | $16113559532439076284373/36738903322033050378240000$ |
| 18 | $3279951793970779750019213/9017104584101486802209280000$ |
| 19 | $58910656224526265294836325869/193146380191453847303322777600000$ |
| 20 | $216462973067870397915716294869/838159909102420090458246374400000$ |

### A uniform proof for every $q\ge21$

On the six-variable simplex, $U_{K_4}$ contains sixteen of the twenty
square-free triple products. The sum of all twenty is at most $5/54$:
replacing any two variables by their average preserves their sum and
increases the third elementary symmetric polynomial. Its maximum is
therefore attained at six equal variables, giving $20/6^3=5/54$.
Consequently the positive integral in (1) satisfies

$$
I_q(K_4;q-1,1,1;1,q-1,q-1)
\ge \left(\frac{54}{5}\right)^q
\frac{((q-2)!)^3}{(3q-1)!}.
$$

Dividing this lower bound by the positive coefficient of $\zeta(3)$ in
(1) gives

$$
\left(\frac{54}{5}\right)^q
\frac{(q-2)!(2q-2)!}{12(3q-1)!}.
\tag{2}
$$

At $q=21$ its exact value is

$$
\frac{2778146948724569113512799974260736}
{2099879042434213173389434814453125}>\frac54,
$$

because four times the numerator minus five times the denominator is
$613192582727210587104025824777319>0$.
The ratio of (2) at $q+1$ to its value at $q$ is

$$
\frac{36(q-1)(2q-1)}{5(3q+1)(3q+2)}>1\qquad(q\ge6),
$$

since its numerator minus denominator equals

$$
27(q-6)^2+171(q-6)+80>0.
$$

Finally,
$\zeta(3)<1+1/8+\int_2^\infty x^{-3}dx=5/4$.
Thus (1)--(2) imply $R_q>0$ for every $q\ge21$. Combined with the finite
table, this proves

$$
R_3=-\frac16,\qquad R_q>0\ (q\ge4),\qquad
\boxed{R_q\ne0\text{ for every integer }q\ge3.}
$$

## 5. The rational direction at $q=n-1$

Let $n\ge4$ and put $q=n-1$. On $W_n$ assign spokes

$$
(q-1,q-1,\underbrace{q-2,\ldots,q-2}_{n-3},q-1)
$$

and cyclic rim powers $(1,\ldots,1,q-1)$. Eliminate vertices
$3,4,\ldots,n-1$. The same beta factors as in Section 3 multiply to
$1/((q-2)!)^3$. The remaining vertices are $0,1,2,n$, and their six powers
are

$$
\alpha_{01}=q-1,\quad \alpha_{02}=\alpha_{0n}=\alpha_{12}=1,
\quad \alpha_{2n}=\alpha_{1n}=q-1.
$$

This is exactly the period (1). The original wheel is therefore

$$
\frac{R_q}{((q-2)!)^3}
+\frac{12}{(q-2)!(2q-2)!}\zeta(3).
\tag{3}
$$

It converges because the final $K_4$ and every beta integral converge.
The coefficient of $1$ is nonzero by Section 4.

## 6. Every odd-zeta direction and completion of the proof

For $2\le q'<n$, the radial wheel with every spoke power $q'-1$ and every
rim power $1$ has exact value

$$
\binom{2n-2}{n-1}
\sum_{s=1}^{\infty}
\frac{2s^2\prod_{j=1}^{q'-2}(s^2-j^2)}
{(2q'-2)!\,s^{2n-1}}.
$$

Its lowest zeta weight is $2n-2q'+1$, with nonzero coefficient
$2\binom{2n-2}{n-1}/(2q'-2)!$.
As $q'$ runs from $2$ to $n-1$, the successive lowest weights are
$2n-3,2n-5,\ldots,3$. These values triangularly span all the odd zeta
values in the theorem. Multiplying numerator and denominator by
$U^{n-1-q'}$ places every radial value in the scalar span at $q=n-1$;
all resulting monomials converge individually.

Subtract the known $\zeta(3)$ contribution from (3) and divide by its
nonzero rational coefficient. This supplies $1$ at the same level.
The upper inclusion and nesting from Section 1 finish the proof.

The verifier additionally compares complete rational coefficient vectors,
graph reductions and radial masters through twelve loops. Its nonzero
determinants are reproducible controls of the constructive proof, without
any assertion that the displayed real zeta values are arithmetically
independent.

## 7. Scope

This is a theorem about the full scalar period space for each fixed wheel.
For a holomorphic kernel the denominator exponent is $q=n+H+m$, with $m$
current insertions. Exact series-chain integration has rational factorial
coefficients at rational sources, so every such wheel-cored coefficient
remains in the displayed space. A specified physical source or insertion
family may span a proper subspace. No equality for that restricted family
is asserted.

For triangular prism, $K_{3,3}$ and $K_5$, existing finite tables are not
promoted here to an all-$q$ closure theorem. The stronger conjectural
threshold $q=3$ for every wheel also remains open in this contribution.
