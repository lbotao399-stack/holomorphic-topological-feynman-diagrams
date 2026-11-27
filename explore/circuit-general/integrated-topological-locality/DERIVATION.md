# Complete integration of a sourced three-terminal topological response

## 1. Original normalization and retained variables

Use the original graph and conventions of
[`htk/CONVENTIONS.md`](../../../code/htk/CONVENTIONS.md), with

\[
V=\{0,1,2,3,4,5\},\qquad
I=(01,02,12,34,35,45,03,14,25).
\]

There are four independent cycles, nine original edge times, and five
independent vertex sources. Vertex 5 is grounded and
\(\lambda_5=-\sum_{v=0}^4\lambda_v\). All sources below are retained.
The holomorphic and topological dimensions are \(H=T=1\), with zero
topological source. The original-edge current insertion has the original
edge orientation; it is not a current on a contracted graph.

For the original unimodular cycle matrix \(C\), the cycle forms and
Kirchhoff polynomial are

\[
\psi_r[\lambda]=\sum_e C_{re}y_e[\lambda]\,dt_e,
\qquad U=\sum_{\text{spanning trees }T}\prod_{e\notin T}t_e.
\]

Let the grounded Gaussian vertex field have covariance
\(\tfrac12\Lambda^{-1}\), where \(\Lambda\) is the conductance
Laplacian. Write \(w_e=(X_{e(0)}-X_{e(1)})/t_e\) for its original
edge current. The density convention is

\[
U^{-3/2}\,
\mathbb E\!\left[\bigwedge_{r=1}^4\sum_e C_{re}w_e\,dt_e\right]
\wedge\bigwedge_{r=1}^4\psi_r[\lambda]
=D(t;\lambda)\,\iota_E\prod_e dt_e.
\tag{1}
\]

Here \(E=\sum_e t_e\partial_{t_e}\). In the code, the coefficient
of the wedge omitting the last original edge is divided by its time.
The sign is positive because that edge has zero-based index 8. This
fixes the normalization of every integral below.

## 2. The three-terminal source response

Delete the star at vertex 0. Its terminals are \(1,2,3\); its original
edge indices are \(0,1,6\). Put

\[
g_1=t_{01}^{-1},\qquad g_2=t_{02}^{-1},\qquad
g_3=t_{03}^{-1},\qquad G=g_1+g_2+g_3,
\qquad h_a=\frac{g_a}{G}.
\]

The three by three response on the terminals is

\[
L_r=\operatorname{diag}(g_1,g_2,g_3)
-\frac1G(g_1,g_2,g_3)^{\mathsf T}(g_1,g_2,g_3).
\tag{2}
\]

Add it to the six untouched outer edges. Grounding vertex 5 gives the
four by four matrix \(\Lambda'\) on vertices \(1,2,3,4\). With
\(g_{ab}=t_{ab}^{-1}\) on the untouched edges, it is exactly

\[
\Lambda'=
\begin{pmatrix}
g_{12}+g_{14}+g_1-g_1^2/G&-g_{12}-g_1g_2/G&-g_1g_3/G&-g_{14}\\
-g_{12}-g_1g_2/G&g_{12}+g_{25}+g_2-g_2^2/G&-g_2g_3/G&0\\
-g_1g_3/G&-g_2g_3/G&g_{34}+g_{35}+g_3-g_3^2/G&-g_{34}\\
-g_{14}&0&-g_{34}&g_{14}+g_{34}+g_{45}
\end{pmatrix}.
\tag{3}
\]

The reduced sources and recovered potential are

\[
\begin{pmatrix}\Phi_1\\\Phi_2\\\Phi_3\\\Phi_4\end{pmatrix}
=(\Lambda')^{-1}
\begin{pmatrix}
\lambda_1+h_1\lambda_0\\
\lambda_2+h_2\lambda_0\\
\lambda_3+h_3\lambda_0\\
\lambda_4
\end{pmatrix},
\qquad
\Phi_0=\sum_{a=1}^3h_a\Phi_a+\frac{\lambda_0}{G},
\qquad \Phi_5=0.
\tag{4}
\]

Every original current is then

\[
y_{ab}=g_{ab}(\Phi_a-\Phi_b),\qquad ab\in I.
\tag{5}
\]

Equations (2)--(5) follow by solving the vertex-0 Kirchhoff equation
first and substituting it into every remaining vertex equation. Thus
they retain the source at the eliminated vertex and its transfer to all
three terminals. They also determine the currents on the three deleted
star edges; no original current slot is lost.

The implementation checks, as complete polynomials, all six original
Kirchhoff equations and every original cycle equation separately for
each of the five independent sources.

## 3. Markov decomposition and all contact terms

Completing the square in the Gaussian quadratic form gives

\[
X^{\mathsf T}\Lambda X
=G\left(X_0-\sum_{a=1}^3h_aX_a\right)^2
+(X_1,X_2,X_3,X_4)\Lambda'
(X_1,X_2,X_3,X_4)^{\mathsf T}.
\]

Consequently

\[
X_0=\sum_{a=1}^3h_aX_a+\xi,
\qquad \mathbb E\xi^2=\frac1{2G},
\qquad
\xi\ \text{is independent of the reduced Gaussian field}.
\tag{6}
\]

Let \(w_e^R\) denote the current recovered from the harmonic extension
\(X_0=\sum_a h_aX_a\). On the star the edges are oriented out of
vertex 0, so

\[
w_{0a}=g_a\left(\sum_bh_bX_b-X_a\right)+g_a\xi,
\qquad w_e=w_e^R\quad(e\notin\{01,02,03\}).
\tag{7}
\]

If an original star orientation is reversed, both its current and its
coefficient of \(\xi\) are reversed. The independent calculation
deleting vertex 1 tests precisely this case.

For each original cycle, define its contact one-form by the actual
coefficient of \(\xi\):

\[
\alpha_r=\sum_{a=1}^3 C_{r,0a}\,g_a\,dt_{0a},
\qquad
\psi_r^\tau=\psi_{r,R}^\tau+\xi\alpha_r.
\tag{8}
\]

The original cycle conservation law at vertex 0 says
\(\sum_a C_{r,0a}=0\). Therefore all \(\alpha_r\) lie in the span
of the two one-forms

\[
g_1dt_{01}-g_3dt_{03},\qquad
g_2dt_{02}-g_3dt_{03}.
\]

Every wedge containing three or more contact one-forms vanishes.
Odd powers of \(\xi\) have zero expectation. At any even loop number
the complete contact formula is therefore

\[
\begin{aligned}
\mathbb E\bigwedge_{r=1}^{n}\psi_r^\tau
={}&\mathbb E_R\bigwedge_{r=1}^{n}\psi_{r,R}^\tau\\
&+\frac1{2G}\sum_{1\le r<s\le n}(-1)^{s-r-1}
\alpha_r\wedge\alpha_s\wedge
\mathbb E_R\bigwedge_{k\notin\{r,s\}}\psi_{k,R}^\tau.
\end{aligned}
\tag{9}
\]

The remaining wedge on the second line is in its original cycle order.
The sign counts the transpositions needed to move the two selected
contact forms to the front. Equation (9) is not a truncation: every
higher contact contribution vanishes by exterior rank. It extends the
adapted two-cycle contact expression in
[`ht-d2/REPORT.md`, §1.5](../ht-d2/REPORT.md) to an arbitrary
three-terminal star and an arbitrary original cycle basis.

For computation, (9) is equivalently the complete covariance

\[
2\mathbb E(w_ew_f)
=2\mathbb E_R(w_e^Rw_f^R)
+\frac{g_eg_f}{G}
\quad(e,f\text{ on the outward-oriented star}),
\tag{10}
\]

with the contact term zero if either edge is outside the star. The
code assembles both terms **before** their common denominator is
cancelled. It checks the full result against the independent recovered
source response. For four distinct exterior edge slots, a product of
two contact covariances vanishes because the star has only three edges.

The certificate records three exact controls of the coefficient of
\(U^{-1/2}\). In each case the single-contact contribution is nonzero,
the double-contact contribution is zero, and their sum equals the
original grounded-Laplacian density. One control sets \(\lambda_0=0\),
so the Gaussian contact cannot be confused with the sourced transfer
term in (4).

## 4. Exact polynomial construction from the reduced response

For this paragraph the variables are conductances. Let
\(K_\Gamma=\det\Lambda\) be the original conductance tree polynomial.
It is computed through the reduced matrix, using

\[
\det\Lambda=G\det\Lambda',\qquad
\det(G\Lambda')=G^3K_\Gamma,
\qquad
\frac{\operatorname{adj}(G\Lambda')}{G^2}
=K_\Gamma(\Lambda')^{-1}.
\tag{11}
\]

Only four by four determinants and three by three cofactors are needed.
All displayed divisions are exact in the conductance polynomial ring.
Equations (4) give each potential with denominator \(K_\Gamma\).
Equations (7),(10) then give every off-diagonal covariance with the same
denominator, after the reduced and Dirichlet contributions are added.

For any collected multiaffine conductance polynomial, substitute
\(g_e=t_e^{-1}\) and multiply by \(\prod_e t_e\). Equation (11)
then gives exactly the original \(U\), with 75 spanning-tree terms.
The code independently enumerates all original spanning trees to check
this normalization. The recovered original currents and covariances
take the form

\[
y_e=\frac{\text{degree-four polynomial, linear in every original source}}{U},
\qquad
2\mathbb E(w_ew_f)=\frac{\text{degree-three polynomial}}{U}
\quad(e\ne f).
\tag{12}
\]

The remaining step is the finite exterior expansion of (1). Each
selection of four of the eight non-gauge edges supplies an integer
cycle minor, a complete three-pair Wick moment, and the complementary
holomorphic minor. Twenty-eight of these exterior splits are nonzero.
The complete common denominator is initially \(U^{15/2}\). The Euler
factor and exactly three factors of \(U\) divide the full numerator,
giving

\[
D(t;\lambda)=\frac{\text{degree-nine time polynomial, degree four in the sources}}{U^{9/2}}.
\tag{13}
\]

Every first-current moment uses the full product of the numerator in
(13) with the original current numerator in (12); its scalar denominator
is \(U^{11/2}\), with time degree 13 and source degree 5. The verifier
never extracts only a chosen current component or a chosen source
coefficient from these products.

At \((\lambda_0,\ldots,\lambda_5)=(1,2,3,4,5,-15)\), the complete
bare time numerator has 2,220 terms. It and all nine current numerators
agree as entire polynomials with the separate direct Pfaffian engine.
Deleting vertex 1, including its incoming original edges, produces the
same entire polynomials. Neither comparison supplies any coefficient
to the response construction.

## 5. The actual positive response chart and its measure

The triangle response resistances are

\[
R_{12}=\frac{t_{01}t_{02}+t_{01}t_{03}+t_{02}t_{03}}{t_{03}},
\quad
R_{13}=\frac{t_{01}t_{02}+t_{01}t_{03}+t_{02}t_{03}}{t_{02}},
\quad
R_{23}=\frac{t_{01}t_{02}+t_{01}t_{03}+t_{02}t_{03}}{t_{01}}.
\]

Writing \(P_\Delta=R_{12}+R_{13}+R_{23}\), the inverse is

\[
t_{01}=\frac{R_{12}R_{13}}{P_\Delta},\qquad
t_{02}=\frac{R_{12}R_{23}}{P_\Delta},\qquad
t_{03}=\frac{R_{13}R_{23}}{P_\Delta}.
\tag{14}
\]

Both maps preserve strict positivity and are inverse on the full
positive cone. Thus the response coordinates have no omitted star
fibre. The untouched edge times remain unchanged. Direct differentiation
gives

\[
\det\frac{\partial(t_{01},t_{02},t_{03})}
{\partial(R_{12},R_{13},R_{23})}
=-\frac{R_{12}R_{13}R_{23}}{P_\Delta^3}.
\tag{15}
\]

The code verifies the cleared determinant identity as a polynomial in
all three resistances. Because (14) is homogeneous of degree one, the
same factor converts the original projective measure. Integration in
the positively oriented response coordinates uses its absolute value;
the original density orientation remains the one fixed by (1).

The complete transformed denominator is

\[
P_\Delta^2U(t(R))=P_\Delta U_{O+\triangle}(R,t_O),
\tag{16}
\]

where the response triangle has its three distinct edges, including
the triangle edge parallel to the original edge 12. Equation (16)
does not insert an additional triangle circulation into the kernel.
The circuit response uses only its two terminal potential differences;
the determinant identity accounts for the redundant triangle cycle.

For each monomial of the complete bare numerator, let its exponents
on the three star times sum to \(s\). Equations (14),(15) give its
exact transformed numerator over the denominator
\((P_\Delta^2U(t(R)))^{9/2}\) by multiplying its substituted
numerator by

\[
R_{12}R_{13}R_{23}\,P_\Delta^{6-s}.
\tag{17}
\]

All exponents in this computation are nonnegative. The full response
numerator and denominator are constructed and hashed in the certificate.
Their degrees are 18 and 6, so projective homogeneity is preserved:
\(18-6(9/2)=-9\). For a first-current moment the corresponding factor
is \(R_{12}R_{13}R_{23}P_\Delta^{8-s}\). The same construction applies
to the entire finite current product.

Integration is carried out by the exact inverse chart (14): its full
pulled-back response numerator is reduced to the original scalar
periods in §6. This is a complete integration of the reduced response,
with the full positive domain and its Jacobian. It does not assert a
lower-dimensional fibre formula or a product-of-simplices evaluation of
the unmerged response coordinates.

## 6. Exact scalar integration and every boundary

The scalar convention is

\[
I_q(\Gamma;\alpha)
=\int_{\mathbb P^{8}_{>0}}
\frac{\prod_{e=0}^8t_e^{\alpha_e-1}}{U^q}\,
\iota_E\prod_{e=0}^8dt_e,
\qquad \sum_e\alpha_e=4q.
\tag{18}
\]

If a complete numerator monomial has exponent vector \(a\), its
scalar weight is exactly \(\alpha_e=a_e+1\); there is no factorial
renormalization. Equations (13),(12) require only \(q=9/2\) and
\(q=11/2\). The declared exact prism tables contain 615 and 3,948
isomorphism orbits at these levels.

For every nonzero monomial, before any automorphism-orbit cancellation,
the code checks positivity, projective degree, and

\[
\sum_{e\in F}\alpha_e>q\,b_1(F)
\quad\text{for every nonempty proper original edge subset }F.
\tag{19}
\]

Here \(b_1(F)=|F|-\operatorname{rank}(F)\), computed by exact graph
connectivity. Forest subsets satisfy (19) automatically because all
weights are positive. The strict inequalities are precisely the
positive projective convergence conditions: on every ordered Hepp
sector, successive coordinate faces contribute the corresponding
power \(\sum_{e\in F}\alpha_e-qb_1(F)-1\); the remaining Kirchhoff
factor is positive and bounded on the resolved closed cube. Thus every
individual monomial used in (18) is absolutely integrable. The proof
and complete relative-boundary discussion for this criterion are also
recorded in
[`core-period-stabilization/DERIVATION.md`](../core-period-stabilization/DERIVATION.md).

Each finite polynomial sum can therefore be integrated term by term.
Automorphisms identify equal absolutely convergent scalar periods; no
integration-by-parts boundary is discarded. All resulting values are
exact rational combinations of \(\pi^2\) and \(\pi^4\), using the
declared scalar input tables. Their complete source dependencies and
the level-\(11/2\) relation certificates are retained by the input
bindings.

### Arbitrary-source convergence before reconstruction

Convergence at fourteen source points alone would not justify an
arbitrary-source integral. The finite-dimensional covariance equations
in the next section concern integrated polynomials; the pointwise
density is not asserted to belong to that smaller space.

The separate check [`source_orbits.py`](source_orbits.py) closes this
issue before any source reconstruction. It independently enumerates all
twelve prism vertex automorphisms and applies them to the fourteen
balanced source points. There are exactly 162 distinct resulting points
in the five independent original sources. Their complete homogeneous
evaluation matrices have exact ranks

\[
\operatorname{rank}(\text{degree-four evaluation matrix})=70,
\qquad
\operatorname{rank}(\text{degree-five evaluation matrix})=126.
\tag{20}
\]

[`source_orbit_certificate.json`](source_orbit_certificate.json) retains
all 162 rows, every source monomial, a selected square minor of each
matrix, and its complete nonzero integer determinant. No integrated
symmetry equation is used in this rank calculation.

Choose the seventy source rows forming the recorded degree-four minor
and invert their evaluation matrix over the rationals. Every coefficient
function of the degree-four source polynomial in (13) is then a fixed
rational linear combination of the seventy specialized density
functions. Apply the same argument with the 126 degree-five rows to
each original-current integrand.

The same evaluation inverses can be applied to the common numerator
before integration. Since \(U\) is invariant under graph automorphisms,
the three exact \(U\) divisions in every checked specialization imply
divisibility by \(U^3\) coefficient by coefficient in the complete
formal-source numerator. Thus (13) holds uniformly in all source
coefficients as well.

Each orbit specialization is obtained from one of the fourteen checked
specializations by a graph automorphism. That operation only permutes
the positive edge times, permutes or reverses the original current
labels, and contributes the established orientation sign. Its absolute
Jacobian is one. Since all nine current slots were checked at every base
source, every orbit specialization is absolutely integrable by (19).
The finite rational combinations just described are consequently
absolutely integrable. This proves coefficientwise absolute convergence
for every bare and first-current source coefficient. It also proves
absolute convergence for arbitrary complex sources and justifies the
polynomial reconstruction after integration in the next section.

## 7. Complete source and current reconstruction

The full bare polynomial has degree 4 in five independent sources, so
it has \(\binom84=70\) possible coefficients. Each of the nine
original first-current polynomials has degree 5 and
\(\binom95=126\) possible coefficients. Thus the original joint
coefficient vector has

\[
70+9\cdot126=1204
\]

entries. All graph automorphism equations, with their exact original
edge-orientation signs, and all five independent integrated Kirchhoff
identities are imposed on this whole vector:

\[
\sum_{e(0)=v}K[y_e]-\sum_{e(1)=v}K[y_e]
=\lambda_vK,\qquad v=0,1,2,3,4.
\tag{21}
\]

The sixth equation follows from balance. These equations are valid
before integration by (4),(5), so they also hold for the integrated
response.

The declared complete covariance calculation has 13,874 equation rows,
exact rank 1,155, and a 49-dimensional solution space. The verifier
replays that complete rank and every basis vector. Fourteen fixed
rational source points have exact evaluation rank 49 in this entire
space. Each point supplies the bare kernel and all nine original
currents, making 140 separately integrated outputs. Exact elimination
uses all 140 rows, checking every surplus equation. No coefficient is
removed by a source specialization: unisolvence determines every entry
of the original 1,204-coefficient vector.

The reconstructed ten polynomials contain respectively

```text
114, 230, 237, 237, 238, 237, 237, 236, 236, 228
```

nonzero rational coefficients of the named constants \(\pi^2,\pi^4\).
The independent generic source \((1,2,3,4,5,-15)\) supplies ten more
fully integrated controls. Only after this calculation does the verifier
read the independent direct source-polynomial file; every coefficient
agrees. The displayed special value in the README follows by evaluating
the complete bare polynomial at \((1,1,0,-1,0,-1)\).

## 8. Scope of the closure

This is a finite exact integrated \(T=1\), three-terminal response for
a non-series-parallel four-loop graph, with arbitrary original rational
sources and all original first-current insertions. The equality of the
source polynomials extends formally to arbitrary complex sources.

The calculation preserves the full response domain, the source
transfer, the Gaussian contact pair, the determinant normalization, and
all projective boundaries. It establishes the stated prism milestone.
It does not claim that arbitrary composite response elements have been
integrated, that arbitrary higher current words close in the same two
scalar levels, or that the general star-merge fibre is elementary.

