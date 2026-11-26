# Marked operators in minimal BCOV theory on $\mathbb C^5$

The minimal fields are $t^k\mathrm{PV}^{i,j}$ with $i+k\le4$, differential $\bar\partial+t\partial$, and degree $i+j+2k-2$. Here $i$ is polyvector rank, $j$ is Dolbeault form degree, and $k$ is descendant level. A closed-string propagator has two primary ends of total polyvector rank four. The free BV kernel $\pi_L=(\partial\otimes1)K_L^C$ also has primary ends of total rank four, with total Dolbeault degree five. On a local monomial its order-two contact contracts a pair of letters. In the flat heat kernel the scalar coefficient of the term missing polyvector direction $i$ is proportional to
$$
(\bar z^i-\bar w^i)e^{-|z-w|^2/4L}.
$$
A holomorphic derivative of any order multiplies it by further antiholomorphic differences. Evaluation at the coincident point gives zero at every $L>0$. Thus this free contact vanishes on point monomials of holomorphic jets. The contraction of two smeared observables still pairs their test sections by $\int\alpha\partial\beta$. Pure closed loops are counted by $\lambda^2$ in the open-closed theory's string coupling. [Costello–Li, minimal field space and cubic interaction, pp. 13–15](https://arxiv.org/pdf/1606.00365#page=13) · [Closed heat kernels and BV differential, §§4.1, 7](https://arxiv.org/pdf/1505.06703#page=65)

## Graph degree

A connected geometric marked graph at loop number $n$ has $5n+1$ internal edges by Dolbeault form degree. Euler's identity gives $4n+1$ interaction vertices and one operator vertex. Suppose $k$ primary letters contract at an operator of length $r$. At a classical interaction of valence $l$ the sum of descendant levels is $l-3$. Since descendants cannot contract, let $D$ be their total level on the surviving interaction legs. Counting those legs, then adding the $r-k$ operator spectators, gives
$$
3(4n+1)+D-2(5n+1)+k+(r-k)=r+2n+1+D .
$$
A point holomorphic primary field coordinate of rank $i$ has dual degree $2-i$. A descendant coordinate has degree $2-i-2k$ at Dolbeault degree zero. The arity and rank counts give operator-degree change $2(2n+1+D)-1-2D=4n+1$. The source coupling has $|\lambda|=3-5=-2$, so multiplication by $\lambda^{2n}$ gives total degree one. [Costello–Li, §4.2](https://arxiv.org/pdf/1505.06703#page=40)

The polyvector ranks on the surviving interaction legs sum to $5(4n+1)-4(5n+1)$ plus the sum of the ranks of the $k$ contracted operator letters, namely one plus that sum. The same counting applies to any selection of primary indices and any spectator jets. [Geometric density degree](../../../code/htk/CONVENTIONS.md) · [Classical descendant weights](https://arxiv.org/pdf/1606.00365#page=14)

At one loop a nonzero graph must be a six-cycle. A bridge would leave a cycle on fewer than six vertices; the five current rows on such a cycle have rank below five after projectivization. The operator consumes two letters. With two scalar primary letters contracted and no descendants there are five outgoing interaction legs of total polyvector rank one, necessarily four scalars and one vector.

## The local operator $\phi^2/2$

Write $\phi$ for the primary scalar coefficient and take four scalar waves and one vector wave with polarization $v$. The ordered five bulk momenta $p_1,\ldots,p_5$ are holomorphic covectors; the marked point is fixed. Put the operator at vertex $0$, use edges $0\to1\to\cdots\to5\to0$, and set $\sum_{e=0}^5t_e=1$. The currents satisfy $y_i-y_{i-1}=p_i$ and $\sum_e t_e y_e=0$. The determinant of five ordered momenta uses the source holomorphic volume form.

To contract the polyvector indices use the five Grassmann coordinates $\theta_i$ at each bulk vertex. Set $\theta_0=0$ at the two scalar evaluation slots. For edge $e$ put $z_e=\theta_e-\theta_{e+1}$ cyclically. Its polyvector tensor is $\iota_{y_e}(z_e^1\cdots z_e^5)$ with the exterior contraction sign $(-1)^{a-1}$ for omission of direction $a$. At the vertex $j$ carrying the vector, $\theta_j(v)=\sum_{e=j}^5z_e(v)$. In each term $z_e(v)\iota_{y_e}(z_e^1\cdots z_e^5)=y_e(v)z_e^1\cdots z_e^5$. The remaining vertex-top coefficient is $(-1)^{e+1}\det(y_0,\ldots,\widehat y_e,\ldots,y_5)$. The signs in edge order are $-,+,-,+,-,+$, by the five rooted incidence minors. Current conservation gives
$$
(-1)^e\det(y_0,\ldots,\widehat y_e,\ldots,y_5)
=t_e\det(p_1,\ldots,p_5).
$$
The ordered polyvector tensor is therefore $-\det(p_1,\ldots,p_5)\sum_{e=j}^5t_e y_e(v)$. All six edge polyvector factors have even degree. For the antiholomorphic factor the six heat markings give $-t_e\det(p_1,\ldots,p_5)$ each before the projective facet integration. Their sum is the reduced density $-\det(p_1,\ldots,p_5)$. This calculation orders the vertex holomorphic volume forms on the left, the edge forms in cyclic order, and the external vector after the six polyvector factors. [Scalar heat facets](../../type-i/README.md#the-hexagon)

Let $T_i=\sum_{e=i}^5t_e$. In the remaining Schwinger integral the coefficient of $p_i(v)$ is $T_{\max(i,j)}-T_iT_j$. The normalized simplex means are $\mathrm E[t_e]=1/6$ and $\mathrm E[t_et_f]=(1+\delta_{ef})/42$. Including the simplex volume yields the reduced ordered result
$$
\frac{\det(p_1,\ldots,p_5)^2}{5!\,42}
\sum_{i=1}^5\min(i,j)(6-\max(i,j))p_i(v).
$$
For the Hessian of $\phi^2/2$ the root factor is one. Each polarized cubic $\mathrm{Tr}(\alpha^3)/3!$ gives its trilinear trace; the root-fixing six-cycle has two automorphisms. Sum over the $5!$ bulk label assignments and divide by two. If $p_v$ denotes the momentum of the vector wave, the reduced result and the result with the normalized heat kernels, respectively, are
$$
\frac{\det(p_1,\ldots,p_5)^2}{24}(p_1+\cdots+p_5+p_v)(v),
\qquad
-\frac{\det(p_1,\ldots,p_5)^2}{24(2\pi i)^5}(p_1+\cdots+p_5+p_v)(v).
$$
Indeed the pure closed and open kernels have the same antiholomorphic heat forms; each closed edge adds exactly one holomorphic divergence factor. Its Gaussian mean is the Kirchhoff covector. Holomorphic Gaussian insertions have no pair contractions with each other. All six divergence signs multiply to $+1$. With the normalized Dolbeault delta form, Gaussian integration over five relative points gives $(-32i)^{-1}(4\pi)^{-5}$, and the five gauge-fixing contractions give $4^5$. Their product is $-1/(2\pi i)^5$. The odd vector is polarized with its parameter on the left. This is the bare marked coefficient before multiplication by $\lambda^2$ and before any finite redefinition of the operator. For $p_i=e^i$, $p_v=e^5$, $v=e_1$, the reduced value is $1/24$. Here $e_i,e^i$ are dual coordinate bases and the vector wave is divergence-free. [Closed propagator, §7](https://arxiv.org/pdf/1505.06703#page=65) · [Heat factor](../../type-i/README.md#the-hexagon) · [Exterior coefficients, exact degree-seven cubature and polarizations](marked_hexagon.py)

For an integrated bilinear with all outgoing primary fields and cubic vertices, momentum conservation gives $p_1+\cdots+p_5=0$. The geometric determinant vanishes before integration for every choice of bilinear jets and polyvector indices. At a marked point, the incoming operator may absorb nonzero momentum.

The classical scalar variation contains $\partial(\phi\mu)$ from $\partial(\alpha^2)/2$, with $\mu$ the primary vector. Already $\mu=\partial_{z^1}$, $\phi=1+z^1$ gives a nonzero classical variation of $\phi^2/2$ at the origin. A cohomology class of this operator therefore does not exist in the unrestricted classical local complex. At the completely differentiated five-primary projection, a genus-one bulk counterterm entering through the one-edge bracket always leaves an undifferentiated scalar spectator from the bilinear. It cannot alter that projection. A finite genus-one redefinition of the marked operator changes its differential by the classical differential of that redefinition. The elementary evaluation above computes the bare operation; the vanishing of a six-field bulk anomaly does not determine it.

### Equivariant redefinitions of this marked operator

Project onto holomorphic waves with four scalar outputs and a transverse vector output, $p_v(v)=0$. Keep total derivative order eleven and at least two holomorphic derivatives on each of the five outputs. The displayed determinant square has exactly two powers of each covector. Its remaining factor increases the degree at one output by one. It therefore belongs to this projection and evaluates to $1/24$ on the indicated transverse input in reduced units.

A translation- and $SL_5$-equivariant finite change $\phi^2/2\mapsto\phi^2/2+\lambda^2 O_1$ adds $Q_0 O_1$. For the primary output restriction the two possible contributions are as follows.

- In the nonlinear contribution the contracted end at $O_1$ is primary. The classical genus-zero vertex with only primary external legs is cubic by its descendant-level rule. Replacing one letter by the classical quadratic variation increases both arity and total polyvector rank by one, and applies one holomorphic derivative. Hence the contributing part of $O_1$ would be a point scalar with four primary scalar letters of total derivative degree ten. Its polarized symbol would be an invariant homogeneous polynomial on four covectors. An element of $SL_5$ scales their span by any common factor and compensates the determinant on a complementary line on the dense set of independent quadruples. A positive homogeneous invariant must therefore vanish there, and polynomiality extends the vanishing everywhere.
- In the linear contribution the Dolbeault differential vanishes on holomorphic waves. The $t\partial$ term could only replace a single level-one scalar by the divergence of the primary vector, leaving four primary scalar spectators. Every such term contains $p_v(v)=0$. Higher descendant levels, antiholomorphic jets and positive Dolbeault components do not survive this projection.

A genus-one bulk counterterm acting by the tree bracket has one contracted operator letter and one undifferentiated scalar spectator. It vanishes after the fully differentiated projection. Consequently the bare coefficient in this projection is invariant under the specified marked and bulk changes. This uses the given genus-zero interaction and preserves its heat normalization and the normalization of the marked operator. A change of basis on the entire complex of composite operators also acts on $Q_0(\phi^2/2)$, which is nonzero; invariance under such a transformation is a separate question. An elementary-field change acting as a derivation has commutator with $Q_0$ again a derivation, so its value on $\phi^2/2$ retains the scalar spectator. On transverse waves the displayed symbol has the factor $(p_1+\cdots+p_5)(v)$, which accounts for its vanishing after spacetime integration. [Type and invariance checks](marked_redefinitions.py)

### Signs of the closed heat edges at two loops

In the normalized heat convention above, write the uncontracted closed kernel in the published order as the scalar Gaussian times its five antiholomorphic differences followed by its five polyvector differences. Acting by $\partial_\Omega=\sum_i\partial_{z^i}\partial_{\theta_z^i}$, the left odd derivative crosses five antiholomorphic differences with sign $-1$. The holomorphic derivative of the Gaussian supplies $-(\bar z^i-\bar w^i)/(4t)$. The product gives $+y_i$ times the usual signed interior contraction of the polyvector five-form, for heat and propagator edges alike. In the normalized heat coordinate the gauge fixing is $-4\sum_j\partial_{z^j}\iota_{\bar j}$. Its holomorphic derivative annihilates the antiholomorphic linear prefactor already present and differentiates the Gaussian with a second minus sign. Thus each of the ten propagators contributes $+4y_j$ in its signed antiholomorphic contraction. The eleventh closed divergence is on the BV heat edge. All eleven contracted polyvector factors have degree four. Gathering nine integration volumes has parity $25\binom92=900$. The raw marked heat minors can consequently be multiplied by the raw closed-index tensor with the signs used in the square calculation. The conversion is
$$
(-32i)^{-2}(4\pi)^{-10}4^{10}=(2\pi i)^{-10}.
$$
In this convention the reduced two-loop evaluation $1/405$ gives $\lambda^4/[405(2\pi i)^{10}]$. This specifies the closed heat normalization on the marked probe; no independent conversion to a ten-dimensional supergravity coupling has been used. [Five BV and 25 propagator signs](closed_heat_sign.py) · [Direct raw graph minors](../point_graph.py) · [Reduced theta sum](../marked_theta.py) · [Eleven heat factors and nine Gaussian integrations](../../type-i/README.md#two-loops-eight-external-gauge-fields)

## An integrated class with cubic scalar leading term

Let $\eta$ be the constant holomorphic top polyvector with $\Omega(\eta)=1$. Take an odd formal parameter of degree $-(2r+3)$ multiplying the full BCOV field $t^r\eta$, for $r\ge0$. This is a closed constant background. Translation in this direction commutes in the graded sense with $\bar\partial+t\partial$ and preserves the constant Poisson tensor. The coefficient linear in the parameter in the translated full classical interaction is therefore a classical cocycle. Restriction to minimal fields is a cochain map for this operation: the minimal field space is a subcomplex containing both legs of the full Poisson tensor. [Full potential and minimal field space](https://arxiv.org/pdf/1606.00365#page=14)

The top polyvector forces all other inputs to be scalar. Write $\phi_k$ for the level-$k$ scalar superfield, with $\phi_0=\phi$ and $0\le k\le4$. Products under the integral select Dolbeault degree five and use the fixed holomorphic volume form. The completed functional, of degree $2r-1$, is
$$
F_r=\sum_{m\ge r+2}\frac{(m-2)!}{m!r!}
\sum_{\substack{k_1+\cdots+k_m=m-r-2\\0\le k_i\le4\;\mathrm{descendant\ levels}}}
\int\Omega\,\prod_{i=1}^m\frac{\phi_{k_i}}{k_i!}.
$$
At $r=0$ this gives the integrated scalar bilinear, whose genus-two lift is derived in the [fivefold calculation](../README.md#the-integrated-bilinear-and-simultaneous-scalar-deformations). Its leading holomorphic-derivative-free term is $\int\Omega\,\phi^{r+2}/(r+2)!$. It is nonexact in the translation-invariant classical local functional complex: every term in the classical differential and every total derivative has positive jet order. The displayed leading polynomial has nonzero zero-momentum symbol after the top Dolbeault polarization.

The scalar/Beltrami check can be performed without any higher-rank fields. Solving the formal equation $u=\sum_{k=0}^4\phi_k u^k/k!$ in field degree, Lagrange inversion gives
$$
\frac{\partial F_r}{\partial\phi_j}
=\frac{u^{r+j+1}}{(r+j+1)r!j!}
$$
for the derivative of the integrand with the common density removed. In the primary scalar Hamiltonian variation the coefficient of a level-$j$ vector $\mu_j$ inside the divergence is $u^{j+1}/(j+1)!$, obtained by choosing that vector in one slot of the differentiated interaction. After integration by parts the two contributions cancel by
$$
d\left(\frac{\partial F_r}{\partial\phi_{j+1}}\right)
=\frac{u^{j+1}}{(j+1)!}\,
d\left(\frac{\partial F_r}{\partial\phi_0}\right),\qquad 0\le j\le3 .
$$
Here $d$ differentiates the scalar coefficients. In particular $F_1$ at $\phi_2=\phi_3=\phi_4=0$ is $\int\Omega\,\phi^3/[6(1-\phi_1)^2]$. The background argument proves closure on the full minimal complex. [Descendant coefficients and twenty closure identities](integrated_descendants.py)

### Closed point classes from scalar source insertions

The flat Li–Wang interaction supplies the full descendant master equation used in the source differentiation below. Consequently all five $[E_r]_5$, $0\le r\le4$, lift to every loop order over the formal string coupling. Here $[E_r]_5$ is the polyvector-rank-five component of the local Euler derivative with respect to the level-$r$ field. A quantum primitive would give a classical primitive for its leading term, excluded by the constant-field test below. On the projected saturated planar probe, the $[E_0]_5$ coefficient at each genus $n\ge2$ is $-n!(n-1)!(n-1)^2(2n-1)J_n/(3n-1)>0$ in reduced units. The exact-improvement calculation below applies at every such genus with the field coordinates and bulk action fixed. Thus this closed local operator actually requires a loop-dependent representative. [All-genus input and trace comparison](../README.md#all-genus-primary-vertex-from-the-worldsheet)

The restriction of every closed holomorphic point observable to primary scalar jets is constant. Indeed keep arbitrary jets of $\phi$, turn on a constant primary vector with a single odd parameter, and set all other fields to zero. To first order in that parameter only $\phi$ varies, by the directional derivative of $\phi$ along the constant vector. The free scalar-descendant variation is the divergence of a constant, and the nonlinear vector variation requires two vector factors or a higher polyvector. The same restrictions apply after differentiating the field jets. Thus the scalar restriction has zero total derivative in the first coordinate direction. If its highest jet order is $m$, the coefficients of the order-$(m+1)$ jets obtained by appending this first derivative are precisely its partial derivatives with respect to the distinct order-$m$ jets. They must all vanish. Descending in $m$ gives the claim, field degree by field degree in the completed jet algebra. In particular extra tensor and descendant terms cannot complete a nonconstant pure-primary-scalar restriction to a classical point cocycle.

In the preceding integrand $F_r$, replace each scalar $\phi_j$ by the full even holomorphic polyvector superfield $\alpha_j$ and take its polyvector-rank-five coefficient at a fixed point. The resulting point observable is classically closed for every $r\ge0$. Indeed, with odd parameters ordered on the left as in the bulk hexagon, the component rules are $Q_0\alpha_0=-\partial_\Omega F_0$ and $Q_0\alpha_{i+1}=\partial_\Omega\alpha_i$. They give $Q_0\phi=\operatorname{div}(u\mu)$ and $Q_0\phi_1=-\operatorname{div}\mu$ on the primary scalar/vector restriction. Subscripts after a comma below mean differentiation in $\alpha_i$. The same Lagrange formula, now for $U=\sum_{j=0}^4\alpha_j U^j/j!$, gives

$$
F_{r,i}=\frac{U^{r+i+1}}{(r+i+1)r!i!},\qquad
F_{r,0}F_{0,i}-F_{r,i+1}=F_{r+1,i},\qquad
F_{r,0}F_{0,ij}=F_{r+1,ij}.
$$

The middle equation is needed for $i<4$ since $\partial_\Omega\alpha_4=0$. In the chain rule for the second-order operator $\partial_\Omega$, the linear-divergence and bilinear-bracket coefficients are exactly these two identities. Hence $Q_0F_r=-\partial_\Omega F_{r+1}$. The right-hand side has no polyvector-rank-five coefficient. Exact holomorphic point observables vanish on constant fields because every classical variation contains a spacetime derivative. On constant primary inputs consisting of $r$ scalars, a bivector and a trivector whose wedge has unit volume, the polarized rank-five coefficient is one. This proves nonexactness in the local point complex.

For $0\le r\le4$ a quantum representative follows directly from the solved bulk master equation. Translate the level-$r$ field by an arbitrary compactly supported scalar of Dolbeault degree five with an odd test parameter. Both $\bar\partial$ and $t\partial_\Omega$ annihilate this source. Translation preserves the constant BV kernel and commutes with heat contraction. Differentiate the master equation in the parameter: the marked equation is the linearized quantum master equation for the derivative of the interaction. At finite scale the same differentiation of each connected graph is the sum over placements of the source. In the local interaction integrate all derivatives off the test form; its arbitrary coefficient then selects the Euler derivative point by point. The tree term is the rank-five coefficient above. Thus each of these five classes lifts through genus three in this bulk solution. The same argument applies at any genus where the bulk master equation is solved. For $1\le r\le4$, apply the differentiated equation instead to a level-$(r-1)$ vector-valued top Dolbeault test form. Its free differential is its scalar divergence at level $r$. Integrating that derivative by parts expresses each holomorphic derivative of the level-$r$ class as a quantum differential. Rank-$k$ tests at level $r-k$ give the further descent equations for $k\le r$.
 
The descent continues to all five holomorphic degrees. Let $E_j$ denote the local Euler derivative of the chosen bulk interaction with respect to $\alpha_j$, including the solved genus corrections, with tree part $F_j$. Its allowed polyvector ranks are $j+1,\ldots,5$ since the trace pairs complementary ranks and $\alpha_j$ has ranks $0,\ldots,4-j$. In the holomorphic bottom-component complex the graded source identity gives $Q E_j=-\partial_\Omega E_{j+1}$ on these ranks for $j<4$. The rank-five coefficient of $E_4$ is closed. The tree identity fixes the overall sign, and the linear free differential on the test field is the same at every genus. The next identities are
$$
Q\alpha_0=-\partial_\Omega E_0,\qquad
Q\alpha_{j+1}=+\partial_\Omega\alpha_j.
$$
For the first identity the local mark has one primary end. A graph of positive loop number with this mark has a bridge. Its cycle column vanishes and the Kirchhoff relation on all other edges kills its current cofactor. The surviving one-edge tree is the Hamiltonian variation of the local action. For the second identity a descendant mark has no propagating end at all. Both statements use holomorphic bottom components, so the Dolbeault term vanishes.

In the tower starting with the rank-five component of $E_r$, contract $\Omega$ with polyvector rank $5-k$ to obtain a holomorphic $k$-form. For $0\le k\le r$ use the rank-$(5-k)$ component of $E_{r-k}$. The differentiated source has rank $k$ and level $r-k$, totaling $r\le4$. For $r+1\le k\le5$ use the rank-$(5-k)$ component of $\alpha_{k-r-1}$. Its rank plus level equals $4-r$. At $k=r+1$ the rank-$(5-r)$ component of $E_0$ is followed by the rank-$(4-r)$ component of $\alpha_0$. At $k=r+2$ the free positive-sign identity acts on $\alpha_1$ whenever $r\le3$. These component identities give successive holomorphic exterior derivatives with the exterior contraction signs of $\Omega$. In particular every holomorphic derivative of each of the five starting classes is quantum exact through the solved bulk order. The terminal linear representatives are $\phi_{4-r}\Omega$. Their Dolbeault top components integrate to closed functionals, since the terminal variation is a holomorphic total derivative modulo the Dolbeault differential.

For $r=4$, the terminal pair is the rank-one coefficient of $E_0$ and the primary scalar $\phi$. In the three-loop term, on seven scalar waves and a transverse vector $v$ with odd coefficient $\beta$ on the left, $E_0$ contributes $-C\,\beta\theta(v)$ times the exponential at total external momentum $-p_{\rm root}$. Thus $-\partial_\Omega E_0$ has scalar coefficient $(-C)(-p_{\rm root})(v)=p_{\rm root}(v)C$ after extracting $\beta$. This fixes the elementary-field orientation. The three-loop elementary variation is the terminal holomorphic descent of this protected level-four point class.

On the level-one point class, differentiating the chosen one-descendant bulk counterterm gives the rank-five component of $\det(\partial_i\partial_j\alpha_0)/24$ at one loop, in reduced units. On the level-zero point class, include the scalar source with momentum minus the sum of the external momenta and differentiate the primary bulk action. The six- and eight-field polarized corrections are the trace product times $-P$ and $-C$ at two and three loops. Take external momenta $(e^1,e^2,e^3,e^4,e^5)$ at one loop, $(e^1,e^1,e^2,e^3,e^4,e^5)$ at two loops and $(e^1,e^1,e^1,e^2,e^2,e^3,e^4,e^5)$ at three loops. On each list give the first slot the transverse polarization $e_2\wedge e_3\wedge e_4\wedge e_5$ and the last slot $e_1$, with scalar polarizations elsewhere. The trace product is one. In the chosen bulk action the three reduced quantum corrections to these closed point observables are

$$
\frac1{24},\qquad \frac1{31104},\qquad
\frac{49188779}{2715939072000}-\frac{7783}{523908000}\zeta(3).
$$

These are coefficients of quantum representatives of nonzero classical point classes. Their obstruction to lifting through genus three vanishes. [Classical coefficient identities and three primary insertions](point_descendants.py)

The four-loop seven-planar-wave test fixes a necessary coefficient of the level-zero point representative in every bulk completion of the solved trace restriction. Keep its absorbed momentum $e^4$. At momentum $e^5$ choose the rank-four polarization $e_1\wedge e_2\wedge e_3\wedge e_4$, and at $e^3$ choose the vector $e_5$. Both are transverse and their wedge is the unit top polyvector. Differentiating the negative eleven-field action trace gives $-9072J_4/11$ in the rank-five component of $E_0$. On the elementary scalar, the same trace with vector polarization $e_4$ gives $+9072J_4/11$ because $e^4(e_4)=1$. This identifies the two opposite coefficients on that probe. Extension of the full closed point representative to genus four depends on the unrestricted bulk equation. [Moment sums and wave configuration](../README.md#a-four-loop-measurement-on-the-elementary-scalar)

On the whole real saturated seven-planar-wave locus the same polarization gives $(2\pi i)^{20}[E_0]_5=-(2\pi i)^{20}Q_4\phi(0)$ in the genus-four coefficient of the source-derived representative. The polarization change keeps all ten momenta fixed. The four-vector and vector wedge to the oriented volume polyvector, whereas applying the elementary divergence contracts $e^4$ with $e_4$. The complete period signs thus give a strictly positive coefficient whenever $A>0$. This holds for the specified representative in every completion by a doubly soft primary trace.

The two integer measurements distinguish this representative from its cohomology class. Consider a local primitive with a trivector $\xi$ of momentum $p_\xi$, a vector $v$ of momentum $p_v$, and seven scalar momenta $p_i$. Take the polarized symbol
$$
\begin{aligned}
p_\xi(v)\!
\sum_{\substack{i,j,k,l\ {\rm distinct}\\{\rm scalar\ slots}}}
&\xi^{abc}(p_v)_a(p_i)_b(p_j)_c\,
\det(p_\xi,p_v,p_i,p_j,p_k)^{h_1}\\
&\qquad\cdot\det(p_\xi,p_v,p_i,p_j,p_l)^{h_2}
\det(p_\xi,p_v,p_i,p_k,p_l)^{h_3}.
\end{aligned}
$$
The determinant uses the source volume. The two exponent triples $(h_1,h_2,h_3)=(1,2,4),(1,3,3)$ both give $39$ derivatives. These translation- and $SL_5$-equivariant local polynomials have the field degree of a primitive for the ten-field probe. Add $\lambda^8 Q$ of either primitive to the representative. At genus four only the classical variation enters. On primary waves only a fusion of the four-vector $W=e_1\wedge e_2\wedge e_3\wedge e_4$ with a scalar can leave a trivector and a vector. If that scalar is the balancing wave, all determinant rows span at most four directions. For fusion with a planar scalar $p$, the new momentum is $e^5+p$ and the induced trivector is $\iota_{e^5+p}W$, up to the common cochain sign. All three determinants then require the balancing scalar at slot $i$. With $[p,s]$ the ordered planar bracket, the two exact changes reduce to
$$
\sum_{\substack{p,s,t,u\\{\rm distinct\ planar\ slots}}}
[p,s]\,[s,t]^{h_1}[s,u]^{h_2}[t,u]^{h_3}.
$$
On the two integer configurations above, taking the same exterior order in both columns gives the matrix and determinant
$$
\begin{pmatrix}191032&-236568\\-1207384&-30648\end{pmatrix},
\qquad -291483166848 .
$$
Thus exact local improvements independently shift both genus-four measurements of $[E_0]_5$. The canonical source representative has the specified positive values, the underlying point class is unchanged by these improvements, and the two elementary differential measurements keep their fixed-mark invariance. [Five-dimensional contractions and independent bracket sum](point_exact_improvements.py)

There is an invariant projection of these representative measurements. Multiply each planar wave momentum by an independent scalar and rebalance the last wave. Retain only monomials of degree at least two in every multiplier. To evaluate the retained component at unit multipliers, apply $f(1)-f(0)-f'(0)$ in each multiplier separately, with all other multipliers fixed. A doubly soft source trace already obeys this constraint, since each planar momentum enters at least twice before substituting the balancing wave. This projection kills exact improvements in the translation- and $SL_5$-equivariant holomorphic jet complex with the bulk action and output coordinates fixed.

Indeed, the classical variation adds one holomorphic derivative and one field. A contributing primary primitive has $39$ momentum slots and four polyvector indices. In an $SL_5$-invariant tensor the four polyvector indices pair with four covectors; the other $35$ covectors fill seven volume determinants. A determinant of five polyvector indices is impossible with only four indices. The free descendant terms give the divergences of the two transverse polarizations and vanish. Fusion of two nonplanar waves leaves at most four independent momenta. A planar scalar has zero contraction with $e_5$, and two scalars have zero primary fusion. Thus the only surviving case fuses a planar scalar of momentum $p$ with $W$. Each volume determinant must use the three nonplanar covectors $e^5+p,e^3$ and the balancing wave, together with two other planar covectors. Its value is independent of the multiplier of $p$, since its two planar slots already span the plane. Of the covectors at the primitive, $e^5+p$ annihilates $\iota_pW$. A nonzero contraction of the trivector requires $e^3$, the balancing covector supplying its $e^4$ component, and a remaining planar covector $s$. This yields $\pm[p,s]$. The remaining vector pairs only with $e^5+p$ or the balancing covector, giving $\pm1$. Every surviving term is exactly linear in its fused planar momentum. Hence the projected genus-four value, positive on the nontrivial saturated real locus, is insensitive to these exact improvements. At genus $n\ge2$ the identical argument uses $10n-1$ covectors, four polyvector indices and $2n-1$ determinants for $2n-1$ planar waves, wherever the source trace and bulk completion exist. For an earlier primitive of genus $s$, the marked arity shift $2(n-s)+1$ gives $2s+1$ primary arguments and the rank shift gives total rank four. A genus-$h$ primary trace carries $10h$ derivatives, and a geometric graph of loop number $\ell$ inserts $5\ell$ heat momenta and $5\ell+1$ holomorphic propagator momenta. Thus the required primitive has $10s-1$ derivatives. At $s=0$ this is impossible. For $s>0$ there are $2s-1$ volume determinants, so at most $2s-2$ external planar scalars can enter the mark. Each geometric cubic taking a planar wave contributes an edge-column difference in the $2\ell$-dimensional tensor product of cycles with the plane; at most $2\ell$ such columns survive. A positive-genus bulk trace of genus $h$ takes at most $2h-1$ planar waves by momentum conservation in the transverse quotient. Their sum at total genus $n=s+\ell+\sum_{h\in{\rm bulk\ vertices}}h$ is at most $2n-2$, below the $2n-1$ planar inputs. At tree level a positive-genus trace obeys the same bound. The remaining classical cubic was treated above. Descendant legs cannot propagate, free divergences vanish on the transverse waves, and the coincident free BV contact vanishes on holomorphic jets. This includes lower-genus exact improvements of the representative.

On the elementary-scalar configuration with $t^2=-2(n-1)^2(3n-2)/(n^2(3n-1))$, applying this projection to the genus-$n$ rank-five coefficient of $E_0$ gives
$$
-\frac{n!(n-1)!(n-1)^2(2n-1)}{3n-1}J_n>0
$$
in the reduced units $(-1)^n(2\pi i)^{5n}$. The factors are $6/5$, $30$, $9072/11$ at genera two, three and four. The first two values hold on the designated quantum source class in the constructed flat quantization, with fixed field coordinates and its representative taken modulo exact terms. The genus-four value holds in every bulk completion satisfying the stated trace hypotheses. [All-genus elementary configuration](../README.md#an-elementary-scalar-period-at-each-genus)

### Quantum source classes at every descendant index

Use the full flat interaction (1.9)–(1.10) of Li–Wang on $\operatorname{PV}(\mathbb C^5)[[\zeta]][2]$. Their Theorem 1.11 is announced with proof sketches. For any integer $r\ge0$, take its local Euler derivative $E_r$ in a level-$r$ source direction **before** setting the external fields to their minimal components of rank plus level at most four. Write $[E_r]_5$ for its rank-five coefficient in the affine volume frame. For $r>4$ this constructs a composite of the minimal fields from an auxiliary source in the full theory. [Full field complex, equations (1.9)–(1.10) and Theorem 1.11](https://arxiv.org/html/2610.21927v1)

Indeed $\bar\partial-\zeta\partial_\Omega$ preserves the rank-plus-level cutoff. Every endpoint of the heat BV kernel or propagator is primary; its two ranks sum to four. Thus restriction commutes with the free differential, BV contraction and effective graph sum, including graphs with one distinguished external source. Translate the full fields by an arbitrary compactly supported scalar-valued top Dolbeault test form at level $r$. Its free differential is zero. Differentiating the full effective master equation gives the marked equation, and restriction leaves its kernel unchanged. Integrating derivatives off the arbitrary test form yields $Q[E_r]_5=0$ in the minimal point complex. The same argument with a rank-$k$ test at level $r-k$, retaining its free divergence, gives the successive descent equations. For $r\ge5$ all six components are Euler derivatives: contract the rank-$(5-k)$ part of $E_{r-k}$ into $\Omega$ for $k=0,\ldots,5$. In particular the top integrated descendant of $[E_{s+5}]_5$ is the quantum integrated scalar deformation $F_s$ for every $s\ge0$.

Equivalently, translate the full local action by the constant top polyvector of unit volume at each chosen level $s$, with the existing odd coefficients of degrees $-(2s+3)$. Every local vertex has total polyvector rank five, so terms with two such backgrounds vanish. Translation commutes with the scale-dependent contractions, which generate the mixed effective graphs at positive scale. This gives the simultaneous all-genus lifts of the integrated $F_s$ in the flat worldsheet quantization. For the point classes take a constant primary scalar $b$ and a constant bivector and trivector whose exterior product has unit volume, extracting once in each polarization. The tree coefficient is $b^r/r!$; all positive-genus terms vanish on constant jets. Evaluation here annihilates the full local differential by the zero-current argument below. The different powers of $b$ therefore detect a linearly independent family over the formal loop parameter. [Classical source normalization](#an-integrated-class-with-cubic-scalar-leading-term) · [Constant-field chain map](#quantum-products-and-the-constant-field-test)

#### A finite loop expansion for each protected source

Take the primary scalar $b+B$ with $B=((z^2)^2+(z^3)^2+(z^4)^2+(z^5)^2)/2$ and constant $b$. Extract the coefficient linear in each of the primary bivector $f(z^1)e_2\wedge e_3$ and primary trivector $e_1\wedge e_4\wedge e_5$, where $e_i=\partial/\partial z^i$. Their wedge has unit volume and each is divergence-free. In worldsheet action units, the genus-$g$ coefficient of $[E_r]_5(0)$ on these fields is
$$
\begin{cases}
\displaystyle \frac{b^{r-2g}}{(r-2g)!}\,
\frac{\partial_1^{2g} f(0)}{4^g(2g+1)!},& 2g\le r,\\[4pt]
0,&2g>r .
\end{cases}
$$
These are coefficients of the representative from the source derivative. The constant-field test above establishes the cohomology statement independently.

To compute the coefficient of $b^m$, there are $2g+r-m$ copies of $B$, by the fivefold valence $2g+3+r$ including the absorbed source and the two polyvectors. Each quadratic uses two normal derivatives at zero, and the four transverse target directions supply $8g$ normal derivatives. Eliminating the absorbed momentum by conservation preserves that homogeneity. Thus $2(2g+r-m)=8g$, or $m=r-2g$. The constant scalar factorial is $m!$. To pull back the source cotangent class after forgetting the transverse fields, allow the test source at level $r-j$ to be $(-B)^j/j!$ times the original test form, for $0\le j\le r$. After variation and integration by parts the point insertion is the sum of $(-B)^j E_{r-j}/j!$. At the origin $B=0$, so it is exactly $E_r(0)$. The correlated source then gives the transverse Hessian factor one by the cotangent pullback. For $f=e^{p z^1}$, the only nonzero axial charge besides that of the source is $p$. The one-mark double-ramification formula gives $[x^{2g}]\sinh(px/2)/(px/2)=p^{2g}/[4^g(2g+1)!]$. Zero-charge primary marks contribute factors one, or equivalently lower the source exponent by the string equation. Polarizing the finite jet gives the displayed formula. At $g=0$ it follows directly from the classical multinomial. [Cotangent pullback and Hessian](../README.md#the-transverse-hessian-on-a-general-holomorphic-background) · [One-mark formula, Theorem 1](https://arxiv.org/abs/1211.5273)

For instance at $r=6$ the four coefficients through genus three are $b^6 f(0)/720,\ b^4\partial_1^2 f(0)/576,\ b^2\partial_1^4 f(0)/3840,\ \partial_1^6 f(0)/322560$. At $b=0$ each even source index measures a single loop order. For a fivefold graph check use $\cosh z^i-1$ on each normal axis. At the selected valence each occurrence must supply exactly its quadratic jet. At genus two the four axes must each occur twice; the signed and permuted Fourier assignments reduce to $14$ nonzero orbits. After string at the constant primary mark, the ordinary fivefold Ward graphs sum to $1/1920$ at $p=1$. The one-loop Ward graphs give $1/6$ at $p=-2$. The scalar stable-graph expansion checks all five positive genera through $g=5$, and the full elliptic vertex checks source levels $5,6$ before applying string. [Exact graphs and marked-source replay](../protected_source_tower.py)

#### All tangent primary fields on a transverse quadratic

Keep the transverse quadratic $B$ above. Allow $A_0=B+\alpha(z^1)$, where $\alpha$ is any even primary superfield of minimal polyvector rank zero through four, and put every positive-level field equal to zero. Its coefficients can vary on the line and its polarizations can point in any of the five directions. Use $x$ to count the source index and weight genus $g$ by $\lambda^{2g}$ in worldsheet action units. At the transverse origin the full Euler derivatives satisfy
$$
\sum_{\substack{r\ge0\\\mathrm{source\ levels}}} x^{r+2}E_r
=\exp\!\left(x\int_{-1/2}^{1/2}
       \alpha(z^1+x\lambda s)\,ds\right)-1-x\alpha(z^1).
$$
The integral means termwise Taylor expansion. Integration over the symmetric interval keeps only even powers of $x\lambda$. Exterior products on the right use the superfield order of the action. Projection to polyvector rank five removes the last two terms and gives the entire family of quantum point cocycles on these fields. General primary inputs here are off shell.

At genus $g$ polarize $k$ tangent waves with axial momenta $p_i$. The number of background insertions is $2g+2+r-k$. After eliminating the test momentum by conservation, every monomial has normal derivative degree $8g$. Each background insertion must receive exactly two derivatives at the transverse origin, hence there are $4g$ such insertions and $k=r+2-2g$. The source cotangent pullback used above leaves the unit Hessian factor. On the retained curve the one-mark formula in its top charge degree gives
$$
[y^{2g}]\prod_{\substack{i=1,\ldots,k\\\mathrm{tangent\ marks}}}
       \frac{\sinh(p_i y/2)}{p_i y/2}.
$$
The quotient is one at zero momentum. At positive genus a term with no tangent marks has zero axial charge degree. At genus zero the multinomial gives $E_r=\alpha^{r+2}/(r+2)!$. Each factor for a tangent mark is the Fourier symbol of averaging its field on $[-1/2,1/2]$ as in the display. Taylor division by $k!$ and the equality $k+2g=r+2$ give the exponential. This also proves the formula for polynomial jets by polarization. With a general homogeneous transverse quadratic the genus-$g$ coefficient acquires $(\det\partial_\perp^2B)^g$ by the same fiber integral. In particular a singular quadratic kills every positive genus in this sector. [Vertical Hessian calculation](../README.md#the-transverse-hessian-on-a-general-holomorphic-background)

Rank five requires at least two tangent primary fields, so $[E_r]_5$ terminates at $g=\lfloor r/2\rfloor$ here. For a test on free-closed waves take the primary scalar $e^{-p z^1}$, bivector $e^{p z^1} e_2\wedge e_3$, and constant trivector $e_1\wedge e_4\wedge e_5$. Extract once in each of these three fields. Their axial momenta sum to zero. Only $r=2g+1$ contributes at genus $g$, with coefficient
$$
[y^{2g}]\left(\frac{\sinh(py/2)}{py/2}\right)^2
=\frac{2p^{2g}}{(2g+2)!}.
$$
The identity uses $2\sinh^2(py/2)=\cosh(py)-1$. This is a point insertion with conserved axial momentum. The conserved five-dimensional integrated test uses the constant top-polyvector source separately.

For the independent fivefold Ward test, after string reduction the labeled tangent momenta $(1,2,3,4,5)$ give $55/24$ at one loop. At two loops $(1,2)$ and $(1,1,1)$ give $91/5760$ and $13/1920$. The last test uses an unrestricted level-five source and twelve waves in the actual fivefold graph recursion. [Exact graph and source replay](../protected_source_tower.py)

#### One protected operator with a nonzero coefficient at every genus

Turn on the level-three scalar amplitude $a$ in the same field coordinates as the elementary differential. The scalar background at levels zero through four is $(B-aB^3/6,\ aB^2/2,\ -aB,\ a,\ 0)$. Add the three labeled primary waves in the preceding test, with axial momenta $-p,p,0$ and scalar, bivector and trivector polarizations respectively. All inputs are annihilated by the free differential. For the fixed closed operator $[E_1]_5$, the coefficient linear in each of the three waves at genus $g$ is
$$
a^g p^{2g}\frac{2(3g+1)!}{6^g g!(2g+2)!}.
$$
Indeed the level-three field and its correlated lower components pull back its cotangent class. For each occurrence of this field, a lowering by $j$ trades $j$ independent quadratic insertions for $j$ attached ones. The source class pulls back after Euler integration by parts as above. Normal saturation leaves $4g$ quadratic insertions altogether and requires exactly $g$ copies of the level-three field. In the remaining one-target integral the absorbed level-one mark has zero axial momentum. Forgetting it by the dilaton equation gives the factor $2g-2+(3+g)=3g+1$. The other $g$ marks have level three and zero momentum. The retained primary charges are $-p,p,0$. Their one-target value, divided by $g!$, is $2p^{2g}(3g)!/[6^g g!(2g+2)!]$ by the marked-point formula. [Primary specializations in the one-target formula](../README.md#an-exact-all-loop-scalar-response-on-one-quadratic-background)

With $x=a p^2\lambda^2/6$ and the branch $C(0)=1$ of $C=1+xC^3$, the full measured coefficient is $C^2$. Formal inversion gives $[x^g]C^2=2(3g+1)!/[g!(2g+2)!]$. For $a=p=1$ the tree through genus three give $1,1/3,7/36,5/36$ in worldsheet action units. At genus one the fivefold graph terms with zero through three transverse donors attached to the level-three field are $4/3,-1,0,0$. Independent elliptic partitions give the same four entries. Scalar stable graphs verify the protected insertion through genus three. These are matrix elements of the source-derived representative of a nonexact class; its closure is an identity before specializing the input fields. [Fivefold graphs and scalar stable graphs](../protected_source_tower.py)

More generally, allow constant scalar amplitudes $a_j$ at levels $1\le j\le4$, supplying each amplitude also as $a_j(-B)^l/l!$ at level $j-l$ for $1\le l\le j$. Retain the primary quadratic $B$ and $k\ge2$ labeled primary tangent waves of axial momenta $p_i$ with total zero and ordered wedge of unit rank-five trace. Write $U(s)$ for the classical scalar inverse determined by $U=s+\sum_{j=1,\ldots,4\,;\,\mathrm{levels}}a_jU^j/j!$. The complete genus-$g$ coefficient of $[E_r]_5$ on these waves is
$$
\frac{(2g+k-1)!}{(r+1)!}\,[s^{2g+k-1}]U(s)^{r+1}
\ [y^{2g}]\prod_{\substack{i=1,\ldots,k\\\mathrm{tangent\ marks}}}
          \frac{\sinh(p_i y/2)}{p_i y/2}.
$$
To see this, take $m_j$ copies of each amplitude and let $m$ be their total number. Normal saturation and valence give $r+\sum_j j m_j=2g+k+m-2$. After cotangent pullback the retained source and all $m$ descendants have zero charge. In the permutation formula of Buryak–Shadrin–Spitz–Zvonkine, setting the variable of each primary tangent mark to zero removes it with a factor $X\sinh(p_i X/2)/(p_i X/2)$, where $X$ is the sum of the remaining cotangent variables. The remaining $m+1$ zero-charge marks have numerator $X^{m-1}$, including $1/X$ for $m=0$. Only the leading reciprocal $1/X$ of the universal denominator contributes in top charge degree. Thus the cotangent coefficient is the multinomial $(2g+k+m-2)!/[r!\prod_j(j!)^{m_j}]$ times the axial charge coefficient above. Divide by $\prod_j m_j!$ for the repeated amplitudes. Lagrange inversion in $U(s)$ sums exactly these multinomials and gives the displayed formula. The same multinomial at genus zero follows from the classical action. With $a_1$ retained, the expansion is formal at $a_1=0$. [Marked-point formula, Theorem 2](https://arxiv.org/abs/1211.5273)

#### A normal D3 residue of the third holomorphic descendant

In the inherited flat local vertex, take the rank-two Euler derivative $[E_{2g-1}]_2$ for $g\ge1$. Its contraction with the fivefold volume is the third holomorphic descendant of the protected point class $[E_{2g+2}]_5$. Use $4g$ copies of the same primary quadratic scalar $B$ and one primary polyvector with scalar coefficient $f$. The Taylor coefficient here uses the ordinary background $B$, with the identical-field factorial included. If $H_{ij}=\partial_i\partial_j B$ is its constant Hessian, the scalar factor of the answer in worldsheet action units is
$$
\frac{\big((\operatorname{adj}H)^{ij}\partial_i\partial_j\big)^g f}
     {4^g(2g+1)!}.
$$
Here $\operatorname{adj}H$ is the matrix of cofactors in the given affine volume frame, with contravariant slots. To prove the formula take $f=e^{p_i z^i}$. The genus-$g$ vertex has holomorphic derivative degree $10g$ and determinant character $2g$. After elimination of the source momentum its polynomial is a sum of products of $2g$ five-row minors of input momenta. Each minor takes the distinguished wave at most once. Thus the $4g$ quadratics supply at least $8g$ derivatives, exhausting their available degree. Each contributes its Hessian and $f$ receives precisely $2g$ derivatives. On a nonsingular Hessian the special orthogonal group fixes the answer and its invariants of a single covector are polynomials in $(H^{-1})^{ij}p_i p_j$. Degree $2g$ in $p$ and $4g$ in $H$ force a constant multiple of $\big((\operatorname{adj}H)^{ij}p_i p_j\big)^g$. This is a polynomial identity and extends to singular $H$. The coefficient is fixed by the rank-four unit transverse Hessian and one tangent wave in the source calculation above. For a polyvector the same symbol multiplies its unchanged polarization.

In the D3 model use brane coordinates $z_1,z_2$, normal coordinates $w_1,w_2,w_3$, and the oriented volume $dz_1\wedge dz_2\wedge dw_1\wedge dw_2\wedge dw_3$. The flux is $N(w_1w_2w_3)^{-1}\partial_{z_1}\wedge\partial_{z_2}$ in the Laurent model. Take $B=z_1w_1+z_2w_2$. Its Hessian on $(z_1,w_1,z_2,w_2)$ has two exchange blocks of determinant $-1$, so the adjugate differentiates only along $w_3$, with coefficient $+1$. At first order in the flux, a rank-one correction to a regular scalar lift cannot supply rank two with only other scalar inputs. Contracting the local answer into the volume gives the normal three-form
$$
\frac{N}{4^g(2g+1)}
\frac{dw_1\wedge dw_2\wedge dw_3}{w_1w_2w_3^{2g+1}}.
$$
Its primitive is $-N\,dw_1\wedge dw_2/[2g\,4^g(2g+1)w_1w_2w_3^{2g}]$. The unweighted linking residue is zero. Pairing with a holomorphic normal test function gives $N\partial_{w_3}^{2g}f_{\rm test}(0)/[4^g(2g+1)!]$ with the unit Laurent residue convention. At test $w_3^{2g}$ the first three values are $N/12,N/80,N/448$. Replacing $w_1,w_2$ in $B$ by two independent linear forms and changing normal coordinates by $SL_3$ supplies every direction of this derivative. Pure powers span all symmetric normal derivatives of order $2g$. An additional pure-normal quadratic in $B$ leaves this answer unchanged: all normal-normal cofactors delete a normal row and column from the Hessian, and the two brane rows must use the two remaining normal columns. Hence these cofactors depend only on the two mixed linear forms.

The genus-one descent uses a minimal source. From genus two onward the full source differentiation defines its restriction. Direct fivefold Ward graphs on the two exchange blocks verify the genus-one and genus-two contact, including a flux wave with momentum in three normal directions. Only the flat local vertex has been evaluated here. Puncture-dependent propagators can still change the transferred operation and its weighted residue. [Fivefold replay](../protected_source_tower.py) · [Laurent flux and boundary conventions](../../twisted-holography-N4/NOTE.md)

#### All quadratic backgrounds and the unweighted D3 linking residue

Keep the same source representative and worldsheet action units. On a primary scalar quadratic $B$ of constant Hessian $H$, retain the term linear in a primary field with coefficient $f$. Write $E_r(B;f)$ for this linear term, keeping its original polyvector polarization. At any point on the fivefold,
$$
\sum_{\substack{r\ge0\\\mathrm{source\ levels}}}x^{r+2}E_r(B;f)
=x\exp\!\left(xB+\frac{x^3\lambda^2}{24}\det H\right)
 \sum_{\substack{k\ge0\\\mathrm{even\ derivatives}}}
 \frac{(x\lambda)^{2k}}{4^k(2k+1)!}
 \big((\operatorname{adj}H)^{ij}\partial_i\partial_j\big)^k f-xf .
$$
In particular the value, Hessian and the derivatives of $f$ suffice even at a point where $dB$ is nonzero. The genus-$g$, level-$r$ coefficient has one term for each $0\le j\le g$ with $m=r+1-2g-j\ge0$:
$$
\frac{B^m(\det H)^j}{m!\,24^j j!}
\frac{\big((\operatorname{adj}H)^{ij}\partial_i\partial_j\big)^{g-j}f}
     {4^{g-j}(2g-2j+1)!}.
$$
The tree coefficient is $B^{r+1}f/(r+1)!$.

To justify dropping the linear Taylor term of $B$, vary a single primary momentum by a scalar $t$ and balance it entirely at the distinguished source. Forget that primary mark. In the pulled-back cotangent term, a holomorphic and an antiholomorphic fiber direction must both be supplied by the Betti form. Each vertical differential of the Abel map contains the charge of the forgotten mark, giving a factor $t^2$. In the ordinary cotangent term the only additional contribution is fusion onto the source with level lowered by one. Its fused momentum is the original sum of the two, independent of $t$. At source level zero this additional term is absent. This proves that every primary slot in the local Euler derivative has zero coefficient of its first holomorphic jet on arbitrary primary fields, at every positive genus. The constant slot lowers the source level by the string equation. Multiple constant slots therefore give the factor $\exp(xB)$. The four fivefold graph tests in the replay check the constant term, vanishing first derivative and nonzero second derivative at one and two loops. [Cotangent pullback and single-mark fusion](../README.md#an-all-genus-scalar-target-measurement-of-the-fivefold-differential)

For the remaining quadratic at the point take a plane wave $f=e^{p_i z^i}$. At generic $H,p$ the kernel of $p$ has nonsingular restricted Hessian. Choose a volume frame with one axial coordinate and the other four in that kernel, orthogonal to the axial direction for $H$. Write $h$ for the axial Hessian entry and $D$ for the determinant of the transverse block. The four transverse directions multiply the genus-$g$ axial response by $D^g$ as in the tangent formula above. The average of the axial quadratic in its exponential contributes $x^3\lambda^2 hD/24$, since the centered unit interval has second moment $1/12$. Here $hD=\det H$ and the square of the axial momentum times $D$ equals $(\operatorname{adj}H)^{ij}p_i p_j$. This gives the first display. Coefficientwise polynomiality extends it to singular $H$ and all finite jets of $f$. At genus one the independent elliptic action on primary fields gives $E_r(u)=u^{r-1}\det(\partial_i\partial_j u)/[24(r-1)!]$ for $r\ge1$; linearization gives exactly both terms.

For the D3 test use every homogeneous quadratic vanishing on the brane,
$$
B=z_1L_1(w)+z_2L_2(w)+\tfrac12 H^\perp(w,w).
$$
Here $w$ is the normal position vector, $L_1,L_2$ are linear normal forms and $H^\perp$ is the normal Hessian. Use the oriented normal vector $v$ defined by $\iota_v(dw_1\wedge dw_2\wedge dw_3)=dL_1\wedge dL_2$, and let $\delta=\det H$. The normal block of $\operatorname{adj}H$ is $v\otimes v$ by its two mixed rows, and expansion of $\det H$ gives $\delta=H^\perp(v,v)$. In particular $B(z,tv)=\delta t^2/2$ for all brane coordinates. Use a single D3 flux of strength $N$ and the ordinary repeated quadratic Taylor coefficient. Pair the resulting rank-two Euler derivative, contracted into the volume, with a holomorphic normal test $\varphi$ by the normalized Laurent residue.

For $0\le d\le g$ the level $r=3g-1-d$ gives
$$
\boxed{\displaystyle
\frac{N\,\delta^{\,g-d}}{(2d)!\,4^d24^{g-d}(g-d)!}
\left(\int_0^1 t^{2d}(1+3t^2)^{g-d}\,dt\right)
\partial_v^{2d}\varphi(0).}
$$
All other nonnegative source levels give zero at this genus on these tests. Indeed Laurent integration by parts moves the $2g-2j$ derivatives of the flux onto $B^m\varphi$. Evaluating along $tv$ leaves $(\delta/2)^m(2g-2j)!\partial_v^{2d}\varphi(0)/(2d)!$ if $m\ge0$ and $d=g-j-m\ge0$. The remaining binomial with $1/[2(g-j)+1]=\int_0^1t^{2(g-j)}dt$ gives the displayed coefficient. This includes the singular mixed Hessian of the preceding calculation.

For $B=z_1w_1+z_2w_2+w_3^2/2$ the constant test has a **nonzero** local linking residue at source level $3g-1$. It is $N/12,N/240,29N/181440$ in the first three genera. At genus one the two contributions to the three-form at $z_1=z_2=0$ are $Nf/24$ and $N(w_3^2/2)\partial_{w_3}^2f/24$ times the normal volume, with $f=1/(w_1w_2w_3)$. Each has residue $N/24$. At genus two the tests $1,w_3^2,w_3^4$ on levels $5,4,3$ give $N/240,7N/720,N/80$. The constant residue detects a nonexact normal three-form of this local source insertion. Identifying it with a quantum D3 charge requires adding the punctured transfer graphs. [Ward coefficients and independent Laurent expansions](../quadratic_source_residues.py)

#### An all-genus linking test with cancellation of every multivertex graph

Use the primary scalar $B=z_1w_1+z_2w_2+c w_3^2/2$ and the D3 flux of strength $N$. For each $g\ge1$ take the third descendants of the protected point sources at the levels corresponding to $E_{3g-1-j}$. At tangent origin, let $\operatorname{Res}$ contract the rank-two coefficient into the oriented volume and take the normalized normal linking residue. Let $E_r^{(g)}$ denote the genus-$g$ source measurement, including connected transfer diagrams from the flat bulk vertices. Keep exactly $5g$ regular scalar inputs counting the $j$ moved into the external test, so the $j$th measurement uses $5g-j$ scalars. In the term linear in $N$ the following combination is independent of the internal propagator kernels under the hypotheses below:
$$
\boxed{\displaystyle
[c^g]\sum_{\substack{j=0,\ldots,g\\\mathrm{source\ shifts}}}
\frac{(-c/2)^j}{j!}\,
\operatorname{Res}\!\left(w_3^{2j}[E_{3g-1-j}^{(g)}]_2\right)
=\frac{N}{24^g g!}.}
$$
The weights are prescribed external tests. The first two checks on the flat vertices are $N(1/12-1/24)=N/24$ and $N(1/240-7/1440+1/640)=N/1152$. The genus-one puncture pentagons cancel with the same relative multiplier at every regulated proper time. [Labeled pentagon calculation](../../twisted-holography-N4/README.md#genus-one-pentagons-of-the-third-source-descendant) · [Exact signed moments and root counts](../quadratic_source_residues.py)

Here the vertices are those of the flat descendant action, with two primary internal ends per propagator, each of polyvector rank at most four. Allow arbitrary tensor kernels on those ends. Use a common regulator for which the pushforward to the marked point exists at tangent origin, respects tangential translations, and introduces no negative powers of $c$. Internal kernels are fixed while supplying the external scalar slots; extra background insertions at bulk vertices are counted as vertices. These conditions admit tangent delta derivatives acting in relative coordinates. The calculation concerns the inherited connected diagrams; an independent local term attached to the brane or a change of source representative can modify the measured combination.

Put $r=3g-1$. Before taking $[c^g]$, one may use source weights $(-B)^j/j!$ up to $j=r$, since $B(0,w)=cw_3^2/2$. Consider one connected graph with geometric loop number $\ell$ and marked-vertex genus $h$. At a fixed distribution away from the marked vertex let $k$ be its internal degree, counting both ends of a self-loop, and let $s=1$ if it also carries the exceptional flux leg or the rank-one correction to a mixed scalar, otherwise $s=0$. The axial scalar has zero flux bracket. The source level $r$ and genus $h$ require $b=2h+2+r-k-s$ regular scalar slots at the marked vertex. In the level-$r-j$ term, $j$ of these have been moved into the external test. At fixed incident legs the string equation identifies every constant scalar there with a lowering of the source. Thus multiplication by $\exp(-xB)$ in the source variable cancels all slots contributing their undifferentiated value. Each identical-field split has precisely the factor $1/(j!(b-j)!)$. The single-primary soft calculation above forbids first derivatives independently at every remaining scalar slot. For $h>0$, the marked vertex can survive only if $2b\le10h$. At $h=0$ and $k+s\ge2$ the classical Taylor sum gives the same condition $b=0$ directly. Hard internal fields can carry arbitrary jets in this argument.

Every nonroot vertex needs polyvector rank five. With only scalar external legs it requires at least two internal ends. At order $N$ at most one nonroot vertex can carry an exceptional external leg, in which case it can have internal degree one. For a connected graph with at least two vertices, $\sum(\deg-2)=2\ell-2$ over vertices therefore gives $k+s\le2\ell+1$. The same bound holds for a single vertex with self-loops, whose degree is $2\ell$. For $\ell>0$ it follows that
$$
b-5h=3(g-h)+1-k-s\ \ge\ 3(g-h)-2\ell\ \ge\ \ell.
$$
The remaining unstable classical possibility is $h=0,k=1,s=0$. Its truncated Taylor sum leaves $(-1)^r B^{r+1}/(r+1)!$ at the root. This has degree $r+1=3g$ in $c$ at tangent origin and disappears in $[c^g]$. Thus positive geometric loops give zero already at the common regulator. A contribution with $N$ in a change of the propagator and no exceptional external leg satisfies the same degree bound.

For completeness, with $V$ vertices and $E=V+\ell-1$ internal edges, summing the flat vertex valences $2h_v+3+\text{total descendant level}$ gives $V=4\ell+1$ with an explicit flux, and $V=4\ell$ when the exceptional contribution is a scalar lift in place of a regular scalar or is in the propagator. A tree thus leaves just one genus-$g$ vertex with the explicit flux. In the local generating function the factor $\exp(-xB)$ removes $\exp(xB)$. At $r=3g-1$ only $(\det H)^g f/(24^g g!)$ survives, where $f=(w_1w_2w_3)^{-1}$ and $\det H=c$. Its constant linking residue gives the boxed result. The normal-trivector tests multiplied by $B^j$ can have free divergence, so an identification with a closed charge of the boundary theory requires a compatible boundary test and its contact prescription.

**A compact-source test for the boundary interpretation.** Represent the normal residue by an annular Dolbeault test, with a compact tangent test of unit integral. Suppose a free-closed test in the full descendant theory had the same pairings as the prescribed tests on tangent-bivector Laurent modes at each level. In the free equation at the next level the divergence of its preceding coefficient must be Dolbeault exact. Pair that equation with the Dolbeault-closed normal Laurent class represented on the triple overlap by
$$
Y_j=\frac{\partial_{z_1}\wedge\partial_{z_2}\wedge\partial_{w_3}}
             {w_1w_2w_3^{2j}},\qquad 1\le j\le g.
$$
The Dolbeault-exact term pairs to zero. Moving the holomorphic divergence onto $Y_j$ gives the tangent bivector with coefficient $-2j/(w_1w_2w_3^{2j+1})$ in the stated volume convention. At level $r-j$ its pairing with the prescribed weight is
$$
-\frac{2j}{j!}\left(-\frac{c}{2}\right)^j .
$$
This is nonzero already at $j=1$. The calculation uses only the pairing with this tangent bivector, so adding compact components invisible on that mode cannot restore free closure. For the test concentrated at tangent origin the terms with $z_1w_1+z_2w_2$ give zero; equivalently one extracts the coefficient of $c^j$ before testing. Under the established linear dictionary the primary Laurent polarization at $j=1$ is a nonzero multiple of $(d\theta_3)^2$, the protected single-graviton component $\operatorname{Tr}((\gamma^3)^2)$. Multiplication by the descendant variable turns the displayed divergence into a free boundary. Consequently the nonzero signed measurement and its compact-source obstruction refer to different levels of the complex. [Linear state and boundary charge test](../../twisted-holography-N4/NOTE.md) · [Single-graviton supercharge in the user's paper](../../../sources/LLC_main.tex)

The unweighted test itself admits an explicit compact free completion. Use the project's variable $t=-\zeta$ and $D=\bar\partial+t\partial_\Omega$. Let $\tau_z$ be the compact tangent Dolbeault top form of unit integral, $\eta_w=\partial_{w_1}\wedge\partial_{w_2}\wedge\partial_{w_3}$ and $\chi$ a normal cutoff constant one near the deleted plane and zero at large normal distance. Choose the annular orientation to give the unit residue. At level $r$ take
$$
D\!\left(t^r\chi\,\tau_z\,\eta_w\right).
$$
Both terms have compact support in the punctured space since only derivatives of $\chi$ occur. The source is closed by $D^2=0$. Its leading term pairs to one with the divergence-free unit flux, excluding a compact primitive. Inserting $w_3^{2j}$ into the potential for $j>0$ additionally gives $\chi\,\partial_\Omega(w_3^{2j}\eta_w)$, which extends to the deleted plane. The preceding Laurent pairing proves that there is no compact substitute preserving the same test values. Thus among the independently prescribed level weights in the boxed sum only the constant weight admits a compact free completion, up to an overall scalar.

The extra term in the compact constant-weight source has rank two and level $r+1$. At $r=3g-1$, on $5g$ regular scalars and one explicit flux, the valence count for a graph using this term gives $V=4\ell$ vertices and $E=5\ell-1$ edges. It has no local single vertex. On componentwise polyvector propagators each edge has degree four in odd differences, but $4E=5(V-1)+1$, so their product vanishes. With one rank-one scalar lift instead, $V=4\ell-1$ and $4E=5(V-1)+2$ gives the same zero. This establishes the local value and removes these additional graphs in the componentwise transfer. The original unweighted graphs, and possible additional brane interactions, still enter the completed boundary comparison.

#### Normal derivatives required of a relative D3 source

Let $P(w)$ be the holomorphic weight of one tangent-bivector residue test at source level $k$. At the D3 plane, where the cutoff is one, the extension of its compact potential to the filled normal ball has differential
$$
D(t^k\chi\tau_zP\eta_w)\big|_{\chi=1}
=t^{k+1}\tau_z\,\partial_\Omega(P\eta_w),\qquad
\iota_{\partial_\Omega(P\eta_w)}(dw_1\wedge dw_2\wedge dw_3)=d_wP .
$$
Here $d_w$ differentiates only the normal variables. To verify what part of this germ is forced by the prescribed measurements, take any nonnegative normal monomial exponent $\beta$ and any coordinate $w_i$. Pair the divergence of $w^{-\beta-1}\partial_{z_1}\wedge\partial_{z_2}\wedge\partial_{w_i}$ with $P$. In the orientation of the preceding test its residue is
$$
-(\beta_i+1)[w^{\beta+e_i}]P=-[w^\beta]\partial_{w_i}P .
$$
Here $[w^\beta]$ is ordinary monomial coefficient extraction and $e_i$ increments exponent $i$. Thus all three derivatives are measured coefficientwise. A constant is the full kernel on holomorphic germs and has the compact completion above. At nonzero $N$ the same obstruction holds for the flux-linearized differential on the punctured space. The flux and each displayed Laurent test have the two tangent vector factors in common, so their Schouten bracket vanishes (one bracket term removes at most one repeated factor). Equivalently one can pair with smooth Dolbeault representatives and move the flux bracket onto the Laurent test.

Under a volume-preserving change of normal frame the axial degree-$2j$ measurements use powers of an arbitrary linear normal form. These powers span all homogeneous polynomials of degree $2j$ by polarization. Their germs are precisely the gradients of that degree, with dimension $\binom{2j+2}{2}$. For the pure weight $w_3^{2j}$ the germ is $2j w_3^{2j-1}\partial_{w_1}\wedge\partial_{w_2}$. Its first nonzero normal derivative at the plane has order $2j-1$ and value $(2j)!$ in that component. Multiplying by the prescribed $(-c/2)^j/j!$ restores the weight of the linking combination. At genus $g$, the levels $3g-1-j$ distinguish $1\le j\le g$. Extracting each power of $c$ independently gives $g(4g^2+15g+17)/6$ independent polynomial gradients in the polarized family. To reproduce the highest-degree germ from its Taylor coefficients on the reduced plane requires its normal jet through order $2g-1$. No fixed finite truncation of these coefficient jets reproduces this family for every genus. This statement specifies the data of this cutoff extension and its Laurent pairings. A dynamical disk vertex requires a separate bulk-to-brane map. [Independent regular divergences, Laurent residues and exact ranks](../../twisted-holography-N4/code/normal_source_jets.py)

For $p=2j$, the project's Laurent dictionary sends the axial primitive to $-\theta_3(d\theta_3)^{p-1}/(p-1)!$ and its divergence to $-(d\theta_3)^p/(p-1)!$. Multiplication by $t$ makes the second mode a free boundary of the first. The primary mode instead corresponds in the linear dictionary to the polarization of $\operatorname{Tr}((\gamma^3)^p)$. At boundary charges $(0,0,p;0,0)$ the physical relative cochains contain only $\gamma^3$; the classical variation vanishes on all of them. A positive-loop $Q_n$ leaving this sector would require negative first and second flavor counts; a $Q_n$ entering it would require negative tangent derivative counts. The traceless diagonal $(1,-1,0,\ldots,0)$ gives trace $2$ for every even $p\ge2$ at rank $N\ge2$. Hence this boundary polarization survives at every loop order. Its Laurent descendant exactness is reproduced by $t d$ in the linear form model. An identification of the nonzero radial obstruction with an ordinary boundary supercharge variation on this polarization is excluded. [Single-particle dictionary](../../twisted-holography-N4/NOTE.md#2-bulk-single-particle-states-drv-chk-at-t0) · [Charge rule in the paper](../../../sources/LLC_main.tex#L443)

#### All-genus measurement on a cubic normalizable D3 state

Take the Laurent **scalar** $-\partial_{w_I} f$ for $f=(w_1w_2w_3)^{-1}$. Its linear state is $+dz_1dz_2d\theta_I$, a polarization of the cubic boundary graviton: the sign from the first normal pole derivative and the sign from the complementary tangent form in the dictionary cancel for $-\partial_{w_I} f$. It is closed on the punctured fivefold, including the flux bracket: its tangent derivatives vanish and the product of two normal Dolbeault two-forms is zero. A normal test $w_J$ has normalized residue $\delta_{IJ}$ on these three scalar modes. [Signed Laurent map](../../../code/holography/derham.py)

Use $B=z_1w_1+z_2w_2+c w_3^2/2$ at tangent origin, in the worldsheet action units above. In the $j$th term take exactly $5g-j$ regular scalar copies of $B$ and one copy of $-\partial_{w_I}f$. Use $[E_{3g-1-j}^{(g)}]_0$, the rank-zero Euler coefficient, obtained with a rank-five source of the full theory and minimal dynamical fields. Its volume form is the fifth holomorphic descendant of $[E_{3g+4-j}]_5$. Strip the common tangent volume and use the oriented normal Laurent residue. For $g\ge1$ the inherited connected diagrams give
$$
 [c^g]\sum_{\substack{0\le j\le g\\\mathrm{source\ shifts}}}
 \frac{(-c/2)^j}{j!}\,
 \operatorname{Res}\!\left(w_J w_3^{2j}[E_{3g-1-j}^{(g)}]_0\right)
 = \frac{\delta_{IJ}}{24^g g!}.
$$
The signed Laurent coefficient before pairing is the original scalar $-\partial_{w_I}f$ divided by $24^g g!$. The assumptions on kernels and common regulation are the ones in the all-genus linking test above. Extra brane-local interactions require their own calculation.

At zero explicit flux, half-edge counting with the marked level gives $V=4\ell+1$ and $E=5\ell$ for geometric loop number $\ell$. A tree is the one genus-$g$ vertex. At every nonroot vertex the external fields now have polyvector rank zero; at least two internal ends are needed for the trace. The marked internal degree $k$ is therefore at most $2\ell$ in a multivertex graph. Write $s=1$ if the distinguished Laurent wave is at the mark, and zero otherwise. The signed root string equation removes the undifferentiated value of every regular scalar. At marked genus $h$, the required number $b$ of remaining regular slots in the unshifted term obeys
$$
 b-5h=3(g-h)+1-k-s\ge\ell .
$$
Every such slot now requires at least two holomorphic derivatives and the local budget is $10h$. For $\ell>0$ the graph vanishes. The unstable classical leaf has the previously computed factor $B^{3g}$ and cannot enter the coefficient of $c^g$. For explicit additional fluxes and lifted external leaves, each additional flux raises the output Dolbeault degree by one. Indeed $m$ explicit flux leaves change $V$ to $4\ell+1+m$, while adding $2m$ external degrees; $4E-5(V-1)$ then equals $-m$. A flux bracket and its homotopy in a leaf add one degree directly. The Laurent output has degree two already at zero flux. This excludes further powers when the flux expansion uses the stated homogeneous kernels.

For the remaining local vertex the signed source shifts remove $\exp(xB)$ in the generating function above. At level $3g-1$ only the undifferentiated input multiplied by $(\det\partial^2 B)^g/(24^g g!)$ remains; $\det\partial^2B=c$. The axial choice $I=J=3$ tests the double pole independently of the simple-pole calculation. Before the shift weights, setting $c=1$ gives $(1/6,1/4)$ at genus one and $(1/72,1/24,1/16)$ at genus two. Their signed sums are $1/24$ and $1/1152$. Fivefold Ward coefficients give these same two rows. [Laurent and Ward replay with the ward option](../quadratic_source_residues.py)

The first compact-source obstruction is present already at $j=0$. The Laurent vector $f\partial_{w_I}$ has divergence $\partial_{w_I}f$, whose pairing with $w_J$ is $-\delta_{IJ}$. A compact free-closed source in the full descendant tower cannot preserve this pairing. The regular normal divergence of $w_1\partial_{z_1}\wedge\partial_{z_2}\wedge\partial_{w_1}\wedge\partial_{w_2}\wedge\partial_{w_3}$ is the regular rank-four tensor in the cubic disk calculation. In this obstruction it occurs at level $3g$. More generally the normal germ of the $j$th prescribed test differentiates $w_I w_3^{2j}$ and has normal degree $2j$, at level $3g-j$. In the normal Koszul symbol the two surviving normal vector factors become odd coefficients, the $2j$ normal powers become derivatives with respect to odd coordinates, and both tangent derivatives remain. Its open polyvector rank is $2j+2$. The cyclic disk arity rule adds two slots per descendant level and one trace slot, giving $6g+3$ open slots for every $j$. In the full closed theory the corresponding positive descendant is the free differential of the top polynomial primitive one level lower; both polyvectors wedge to zero with the flux even after one Schouten contraction. These measurements distinguish the punctured source response from the primary regular disk coupling to three open fields. [Regular disk and primary cohomology](../../twisted-holography-N4/README.md#which-two-ghost-measurements-survive-the-boundary-differential)

For these axial germs the ordinary open propagator cannot contract the single descendant disk down to a cubic. At least one of the three odd directions is absent from the odd derivatives of every monomial of $\partial_w(w_Iw_3^{2j})$. Along that direction each open edge between vertices $u,v$ supplies the undifferentiated difference $\theta_u-\theta_v$. On a connected open graph with a cycle there are more edge factors than independent differences; their exterior product vanishes, also for a self-loop. Attaching ordinary cubic open vertices in a tree only increases the external arity. This rules out the open sewings of the isolated inherited descendant vertex. [Exterior graph check](../../twisted-holography-N4/code/regular_cubic_disk.py)

Putting the regular backgrounds back on the open–closed surface changes the degree test. The source $t^{3g-j}\partial_\Omega(w_Iw_3^{2j}\eta)$, where $\eta$ is the ordered top polyvector, has degree $2+6g-2j$. The $5g-j$ regular primary scalars contribute $-2(5g-j)$, leaving $2-4g$ including the Dolbeault-smeared two ghosts and scalar. On a single nonempty trace the fivefold axiom therefore requires
$$
 2g_{\rm s}+f=2g ,
$$
where $g_{\rm s}$ is the surface genus and $f$ the number of empty boundaries. Consider sewing only the known disks with one interior closed insertion and the ordinary cubic open vertices. The mixed background $z_iw_i$ gives a two-flag open disk, and $w_3^2/2$ gives a three-flag disk. After the external weight $c^j$, extracting $c^g$ uses $4g$ mixed disks and $g-j$ axial disks. The open ribbon graph has cycle number $2g_{\rm s}+f=2g$. The external valence equals two minus twice the cycle number plus the sum of vertex valences minus two. Its exact value is $3+3g-j$ plus the number of ordinary cubic vertices, always greater than three. Hence this entire sector of regular open sewings vanishes even with all the prescribed backgrounds. Multi-interior disks, closed propagators, independent higher-surface vertices and the singular puncture map are outside this count.

There is a stronger vanishing for separated graphs on the flat regular test. Let $S$ vertices be in the fivefold and $D$ on its D3 plane. The two external antiholomorphic degrees require total propagator degree $5S+2D-2$. Translations in the two common tangent directions annihilate every ordinary flat propagator, leaving exactly $5S+2D-2$ relative antiholomorphic directions. Each open or closed edge is also annihilated by simultaneous antiholomorphic dilation: its angular form is contraction of its edge volume by the antiholomorphic edge displacement. This remains true after restriction to the plane at either end and after holomorphic jet differentiation. The dilation vector is nonzero modulo the translations on a dense set for $S+D\ge2$. Its horizontal forms cannot fill the relative top degree. Thus every such separated graph vanishes at finite positive proper times. A mixed graph with two bulk and three brane vertices requires degree $14$ and has rank $13$ in a direct incidence calculation. This leaves the intrinsic local marked vertex in that regular projection. A heat-edge Ward term or a singular puncture propagator has different form data. [Independent open and mixed graph tests](../../twisted-holography-N4/code/regular_cubic_disk.py)

#### A separating closed edge in the Ward equation for the regular cubic

Use the positive cotangent coefficients of the worldsheet vertices here. At genus $g$ and shift $j$ the regular primitive of the free-exact source has level $3g-j-1$, normal polynomial $w_I w_3^{2j}$, and top polyvector $\eta$. Its free divergence has the level and polarization of the regular source above. A homogeneous free Ward solution at fixed arity has zero evaluation on that exact input and free-closed companions. The primitive belongs to the full external field complex. In a local cyclic vertex one can smear the open scalar by a compact Dolbeault top form and use two holomorphic ghost germs on its support. The inhomogeneous Ward equation also has separating brackets and nonseparating BV contractions.

Consider the particular separating edge joining a closed genus-$h$ vertex, $0\le h\le g$, to the open surface. Write $n=g-h$, and $e=1$ for $I=1,2$, $e=0$ for $I=3$. Put $B_0=z_1w_1+z_2w_2$. After the original source weights and through axial order $c^g$, the Euler derivative on the primary end of the closed vertex is
$$
 \frac{c^h w_I B_0^{3n}\eta}{(3n)!}
 \sum_{\substack{0\le k\le h\\\mathrm{derivative\ pairs}}}
 \frac{(3n+k+1+e)!}{(3n+1+e)!\,24^{h-k}(h-k)!\,2^k(2k+1)!}.
$$
The marked primary Hamiltonian variation in the coordinate Ward convention applies $-\partial_\Omega$ to this expression. This sign is before moving a term to the opposite side of the Ward equation. It follows by transposing the even derivatives in the quadratic Euler operator from the varying primary to the source. Its $k$ cofactor derivatives act on the binomial sum of $w_I(-cw_3^2/2)^j B^{3n+k-j}/[j!(3n+k-j)!]$, equal to $w_I B_0^{3n+k}/(3n+k)!$ through the required order. On this polynomial the axial second derivative vanishes. The remaining derivative is $2c(\partial_{z_1}\partial_{w_1}+\partial_{z_2}\partial_{w_2})$. On $w_I B_0^m$ the operator inside parentheses gives $m(m+1+e)w_I B_0^{m-1}$, proving the sum. Of the $2n$ residual regular scalar slots exactly $n$ contain $B_0$ and $n$ contain $w_3^2/2$ at the extracted order, with the repeated-field factor $1/(n!)^2$.

At $h=g$ there are no residual scalar slots. The other vertex is the known disk on $\partial_\Omega(w_I\eta)$ on the genus-$g$ one-boundary surface with no empty boundary. The Euler coefficients after extracting $c^g$, for a mixed and an axial direction respectively, are $(7/24,5/24)$ at one loop and $(209/5760,13/640)$ at two. The corresponding sums above are strictly positive in every genus. The other separating and nonseparating terms and the intrinsic D3 contact are needed for the full boundary Ward equation. [Direct differentiation before summing the source shifts](../quadratic_source_residues.py) · [Regular disk in its matrix trace units](../../twisted-holography-N4/README.md#regular-d3-product-and-integrated-disk-coefficient)

### Classical multiplication of the closed point observables

Write $[\,]_5$ for the rank-five coefficient in the chosen volume frame. For $r+s\ge1$ the ordinary coincident product satisfies the cochain identity
$$
[F_r]_5[F_s]_5=
\frac{(s-r)(r+s-1)!}{r!s!}[F_0]_5[F_{r+s}]_5 .
$$
Every triple product vanishes. Each bilinear $[F_0]_5[F_m]_5$ for $m\ge1$ is nonexact.

To calculate the products, use the formal polynomial integral
$$
F_r=\frac1{r!}\int_0^U x^r\left(\sum_{j=0}^4\alpha_j x^j/j!-x\right)\,dx .
$$
Differentiating in $\alpha_i$ gives $F_{r,i}$ above since the integrand vanishes at $U$. The integral through the scalar component $u$ has polyvector rank at most four. In the remaining integral shift $x=u+y$ and let $B=\sum_{j=0}^4\alpha_j u^j/j!-u$. This has positive polyvector rank. Let $H=(1-\sum_{j=1}^4\alpha_j u^{j-1}/(j-1)!)^{-1}$. The shifted expression in parentheses has constant coefficient $B$ and linear coefficient $-H^{-1}$. Implicit inversion gives $U-u=HB$ modulo $B^2$. If $R_j$ is the coefficient of $z^j$ in $e^{-uz}\sum_{r\ge0}[F_r]_5z^r$, the shifted integral therefore gives
$$
R_j=\left[\frac{H^{j+1}B^{j+2}}{(j+2)!}+\text{a multiple of }B^{j+3}\right]_5 .
$$
Thus $R_j=0$ for $j>3$. A product of two rank-five coefficients has ten polyvector slots. If its monomials contain at least eight factors of $B$, at least six of those factors must supply rank one. There are only five odd components of its rank-one part. Such a monomial vanishes, so $R_i R_j=0$ for $i+j\ge4$.

For $i+j=3$ only the displayed leading terms contribute. Put $D=HB$ and use $H^{-1}$ as the common weight in their rank-five coefficients. In the product there are seven factors of $D$. At least four must have rank one. With four rank-one factors the other three have rank two and fill all ten slots. With five, the two others have rank at least two and exactly one further slot is available. Hence it suffices to include the rank-one and rank-two parts of $D$, and terms linear in either its rank-three part or the rank-one part of $H^{-1}$. The rank-two part has normal form $e_1\wedge e_2+e_3\wedge e_4$ on a dense set. Take all five independent odd coefficients in the rank-one part. Exterior multiplication on this frame, with the ten trivector choices and five additional weight vector choices, gives
$$
\left[H^{-1}D^3/3!\right]_5\left[H^{-1}D^4/4!\right]_5
=\left[H^{-1}D^2/2!\right]_5\left[H^{-1}D^5/5!\right]_5 .
$$
The weight scalar is invertible and can be divided out of both factors. Multilinearity in the remaining choices and polynomiality in the bivector give $R_1R_2=R_0R_3$ on all fields. [Exact exterior products](point_products.py)

The product of generating series at $x,y$ now equals
$$
e^{u(x+y)}(y-x)\left(R_0R_1+(x+y)R_0R_2+(x+y)^2R_0R_3\right).
$$
Comparing with the series at $x=0$ gives the claimed factor. Since $[F_0]_5$ is odd its square vanishes, which kills every triple product. For nonexactness use constant primary fields with scalar value $\phi$, bivector $e_1\wedge e_2+e_3\wedge e_4$, vector $e_5$ and trivector $e_3\wedge e_4\wedge e_5$, the last two with separate odd coefficients. Here $R_2=R_3=0$ and $R_0R_1$ pairs to unit magnitude. Consequently the bilinear at index $m$ pairs to $\phi^{m-1}/(m-1)!$ up to their common orientation sign. Every classical boundary vanishes on this constant configuration. The calculation concerns the ordinary classical multiplication of the point classes.

### Quantum products and the constant-field test

Evaluation on constant fields kills the full local differential through every order of the chosen flat bulk solution. The free variation and the coincident free BV contact vanish. In a connected graph on a mark and interaction vertices, setting all external momenta to zero sets every Kirchhoff current to zero. At positive geometric loop number the top current Jacobian then vanishes identically. Further holomorphic derivatives multiply by currents or external momenta, since two holomorphic Gaussian differentiations have no mutual contraction. At geometric loop number zero the one physical heat edge has a holomorphic divergence, also zero at these momenta. This argument includes insertions of quantum interaction vertices. It gives a chain map from translation-invariant holomorphic point jets to constant-field polynomials with zero differential. Every lift with nonzero leading restriction survives incoming differentials, including each of its coupling multiples through the solved order.

For two holomorphic point marks use the paper's auxiliary Dolbeault heat contraction between them, normalized by its unit tree contact. It has the same antiholomorphic form as a physical edge and carries no polyvector divergence. Equivalently, attach one field of a decoupled Dolbeault conjugate pair to each mark and project the dressed two-input bracket onto its unique auxiliary contraction. In five dimensions its two bottom components can both have even parity. The second mark stays at the output point. At geometric loop number $n$ the holomorphic heat degree and Euler relation give $5n+1$ edges, $4n$ interaction vertices, two marks, and $5n$ physical contractions. With classical vertices of total descendant level $D$, the total output length grows by $2n+D$, the sum of polyvector ranks is unchanged, and the coordinate degree grows by $4n$. Thus the factor $\lambda^{2n}$ gives a product of degree zero. Inserted quantum vertices carry their own coupling powers.

These collision graphs give local jet coefficients. For fixed real external momenta the spanning-tree formula bounds every edge current by the largest tree flow. The heat form is the pullback of the ordinary volume form on the $5n$ cycle-current coordinates. On its full-rank locus the signed maximal minors of the cycle voltage equations recover the positive times with sum one. Therefore the absolute heat integral is the volume of bounded current regions. The physical contractions and input jet derivatives multiply it by polynomials in bounded currents and external momenta. At fixed numbers of fields and derivatives the integral is a polynomial in the momenta, by homogeneity and coefficientwise interpolation in the multilinear Gaussian formula. Complex coefficients follow by polarization. This is the local collision operation in the completed jet algebra. [Current inverse and Jacobian](../../nonSP/README.md#signed-current-coordinates) · [Absolute convergence with numerator insertions](../../nonSP/verify/k4_periods/README.md#schwinger-numerator-and-convergence)

Apply the local heat-bracket quadratic identity in the enlarged theory. A collapsed subgraph contains both marks, the first alone, the second alone, or neither. In the first case the auxiliary edge contracts inside the subgraph. In the next two it survives to the quotient. In the last case the inner operation is the bulk master-equation defect. Independent assignments of labeled interaction vertices give the two vertex factorials, and the first Betti numbers of the subgraph and its quotient add. The free $t\partial_\Omega$ commutes with the contractions by the free-kernel identity; the auxiliary pair is decoupled. Projecting to its one contraction therefore gives the graded tensor differential on the inputs of this degree-zero operation. With $|A|$ the parity of the first mark, this reads through genus three
$$
Q(A\star B)=(QA)\star B+(-1)^{|A|}A\star(QB)
\pmod{\lambda^8},\qquad (A\star B)|_{\lambda=0}=AB .
$$
The unit heat edge alone gives the stated tree value. At positive geometric loop number the current Jacobian vanishes on constant external fields, just as for the one-mark differential. Thus constant-field evaluation takes the local product to ordinary multiplication coefficient by coefficient. [Two-mark boundary partition in the user's paper](../../../sources/LLC_main.tex#L2849)

Take the five rank-five components of $E_0,\ldots,E_4$ above as the two inputs. Their ten products for $0\le r<s\le4$ lift all seven independent classical symbols $[F_0]_5[F_m]_5$, $1\le m\le7$. For $1\le m\le4$ use the pair $(0,m)$. For $m=5,6,7$ use $(1,4),(2,4),(3,4)$ and divide respectively by $3,5,5$. On the constant witness in the preceding section the seven tree values have common orientation and distinct powers $\phi^{m-1}/(m-1)!$. Hence each has zero outgoing obstruction and survives every incoming differential through genus three. Whether products with proportional tree terms acquire independent derivative classes remains open. [Exact ten-pair test](point_products.py)

One genuine product graph is measured on the scalar marks $A=\phi^2/2$ and $B=\phi$. At genus one, fix $B$ at vertex zero, put four cubic vertices in order, and put $A$ at vertex five. The auxiliary edge goes from five to zero. Each physical end at a mark evaluates at zero polyvector coordinate. With edge times normalized to sum one, the four bulk traces of the five physical contractions give $\det(y_0,\ldots,y_4)=-t_5\det(p_1,\ldots,p_5)$, where the $p_i$ are the four external cubic momenta followed by the momentum of the spectator of $A$. The reduced heat form is $-\det(p_1,\ldots,p_5)$ in this edge order. The oriented current and polyvector factors multiply to $t_5\det(p_1,\ldots,p_5)^2$. Its integral is $\det(p_1,\ldots,p_5)^2/6!$. The moving spectator and four bulk slots admit $5!$ assignments and the two labeled marks fix the cycle orientation. Thus the squared reduced polarized coefficient is $\det(p_1,\ldots,p_5)^4/36$. With the two marks interchanged the moving linear field has no spectator and the fifth injected momentum vanishes. This measures a local product on two scalar cochains. [Exact exterior contraction and order sum](point_product_hexagon.py)

There is also a one-loop seven-field coefficient of the product of the protected classes $[E_0]_5\star[E_1]_5$. Give two constant waves polarizations $e_1\wedge e_2$ and $e_3\wedge e_4\wedge e_5$. On the other five waves set $p_i=e^i$. The first carries $e_2\wedge e_3\wedge e_4\wedge e_5$, the fifth carries $e_2$, and the middle three are scalar. Every polarization is transverse. On primaries the moving tree mark has two letters, and the root has three. At all four cubic vertices and the moving mark there must be a nonzero injection for the five-momentum heat determinant. Thus both constants occupy the root and their trace leaves a scalar end.

In an ordering of the five active waves along the path let $j$ be a site's position from $1$ to $5$, and put $T_j=t_j+\cdots+t_5$. For the coordinate $i$ whose momentum enters at position $j$, the edge current is $y_{e,i}=\mathbf 1_{e\ge j}-T_j$. A polarization at position $j$ has Grassmann variable $\theta_j=-\sum_{e<j}(\theta_e-\theta_{e+1})$. Each physical edge leaves one direction to fill. Consequently the odd integral is the determinant of the matrix whose column for a polarization vector $v$ at $j$ has entries $\mathbf 1_{e<j}y_e(v)$, $0\le e\le4$, with the vertex/exterior ordering of the preceding scalar test. Multiply by the negative of the determinant of the five ordered momenta. All $5!$ orders with a single auxiliary contraction enter with unit weight. For descending $j_1,\ldots,j_m$ the normalized simplex mean of the tails in the column expansion is
$$
\mathbb E[T_{j_1}\cdots T_{j_m}]
=\prod_{h=1}^m\frac{6-j_h+h-1}{6+h-1}.
$$
The sum of the oriented current and exterior factors integrates to $1/45$. The level-zero representative has no primary one-loop term. The level-one representative's Hessian correction differentiates all five active waves, whose polyvector product vanishes because both contain $e_2$. Hence its ordinary product with the tree mark contributes zero, and the whole reduced one-loop coefficient on this fixture has square $1/2025$. The protected class has zero obstruction and acquires a nonzero jet coefficient in this representative. Replacing the fifth polarization by $e_1$ gives a genuine hexagon factor $1/60$, while the ordinary Hessian contribution then has magnitude $1/24$. [Independent degree-five cubature, simplex tails and exterior coefficients](point_product_hexagon.py)

The one-loop correction to multiplication is required by closure. Take the ordinary product of the same two individually corrected observables. On eight primary waves its genus-one variation consists of the terms in which $Q_1$ contracts one letter from each mark. The two terms in which both letters belong to the same mark cancel by the individual cocycle equations; the Hamiltonian bulk tree satisfies the ordinary Leibniz rule. For an exact cross-contraction test take the five active momenta $p_i=e^i$, with scalar polarizations for $i=1,\ldots,4$ and bivector $e_1\wedge e_2$ for $i=5$. The other three momenta are zero, with polarizations $e_2\wedge e_3\wedge e_4\wedge e_5$, $e_1\wedge e_2$, $e_3\wedge e_4\wedge e_5$. A zero-momentum bulk vertex kills the heat determinant. Thus the three constants are mark spectators. Only the four-vector can occupy the single spectator position of the level-zero mark, since the remaining two must have total rank at most five at the level-one mark. Their rank is exactly five, leaving a scalar end there. The two assignments of the contracted marks to the two cycle ends cancel the root-fixing reflection divisor.

Split the two Grassmann coordinates at the common spatial root. Start with the scalar end, pass the five cubic vertices, and finish at the end with the four-vector spectator. All six edges are physical. In the determinant rule above the two columns for $e_1,e_2$ at the bivector vertex $j$ have entries $\mathbf 1_{e<j}y_e(e_1),\mathbf 1_{e<j}y_e(e_2)$, $0\le e\le5$. The four final columns are $y_e(e_2),y_e(e_3),y_e(e_4),y_e(e_5)$. The oriented integral of this determinant times the heat determinant, summed over the $5!$ orders, is $-1/24$. Thus the ordinary product of the two quantum cocycles has nonzero reduced one-loop variation of squared magnitude $1/576$ on this transverse test. Altering the two first-order representatives by individual classical cocycles preserves that variation because $Q_0$ obeys the ordinary Leibniz rule. The marked product supplies the compensating classical variation. [Six-edge exterior tensor, resistance-tail moments and independent degree-seven cubature](point_product_hexagon.py)

The classical variation of the genuine marked correction can also be evaluated on these same eight waves. Transversality kills the free descendant variation. Fusing two active waves leaves only four independent injected momenta, and two constants have zero fusion momentum. When a constant fuses with an active wave, the remaining constants must fit together at the level-one root. This forces the four-vector to participate in the fusion. Its contraction with $e^1$ vanishes; its wedge with the fifth wave's $e_1\wedge e_2$ also vanishes. Fusions with the second, third and fourth scalar waves leave the oriented trivectors $e_3\wedge e_4\wedge e_5$, $-e_2\wedge e_4\wedge e_5$, $e_2\wedge e_3\wedge e_5$, respectively. For each, put its three columns at the fused wave's position and the two bivector columns at the fifth wave's position in the five-edge determinant. In the common orientation fixed by the cross term the three signed moments give the order-$\lambda^2$ Ward identity
$$
-\frac1{24}+\frac1{72}+\frac1{72}+\frac1{72}=0.
$$
[Three separately integrated fusion moments](point_product_hexagon.py)

### All-primary tests

Any descendant letter at the marked vertex is uncontracted because all propagators and the BV kernel have primary ends. Therefore the genuine positive-loop marked graphs on exclusively primary outgoing fields see just the leading $\phi^{r+2}/(r+2)!$. The free $t\partial$ variation can replace a single scalar descendant at level one by a primary vector divergence. Higher levels or two remaining descendants cannot contribute to this primary projection. On transverse holomorphic vector waves that contact is zero.

For $F_1$ with five scalar momenta $p_i$ and one vector wave $(p_v,v)$, impose $\sum_i p_i+p_v=0$. At one loop choose the spectator among the five scalar labels. The other four scalar momenta and $p_v$ have the same squared determinant as all five scalar momenta, by conservation. The ordered reduced result above gives
$$
\frac{\det(p_1,\ldots,p_5)^2}{24}
\sum_{\substack{i=1,\ldots,5\\\mathrm{scalar\ spectators}}}(p_v-p_i)(v)
=\frac14\det(p_1,\ldots,p_5)^2 p_v(v).
$$
In coordinate components, let $G=\int\Omega\,\phi_1\det(\partial_i\partial_j\phi)$, with the determinant taken in the volume-normalized basis. The polarized Hessian has coefficient $\det(p_1,\ldots,p_5)^2$ by determinant multiplicativity. The free variation of $G$ on these inputs is its negative times $p_v(v)$. Indeed $|G|=5$ gives the graded-dual sign $-(-1)^5=+1$. With the vector input last and of Dolbeault degree $j$, the preceding inputs give $(-1)^{5-j}$ and the holomorphic divergence on the form coefficient gives $(-1)^j$. Thus $+G/4$ cancels the displayed reduced primary projection. With the normalized heat coefficient it is $-G/[4(2\pi i)^5]$. This verifies the primary projection of the genus-one equation. Its descendant-output equations may require additional terms.

For the bare genus-two equation at the integrated cubic, one classical nonlinear variation of a possible marked correction producing seven scalar waves and one vector requires seven primary scalar inputs and twenty holomorphic derivatives. The only possible free descendant term is longitudinal in the vector wave. A one-loop hexagon acting on a genus-one marked correction would require five primary scalar inputs and ten derivatives. With only five conserved scalar covectors every positive-degree $SL_5$-invariant polynomial vanishes. In particular $G$ has a scalar descendant spectator in any such hexagon. These statements concern holomorphic jet functionals in the fixed classical action. General smooth translation-invariant terms with antiholomorphic derivatives on contracted letters require a reduction to that model; holomorphic external waves alone cannot eliminate internal antiholomorphic derivatives. A genus-one bulk interaction inserted in a wheel, the genus-one tree bracket on a genus-one marked correction, and a genus-two bulk tree bracket on the classical mark also enter the expanded differential when the bulk action itself is corrected. [Two-source comparison for seven primary scalar inputs](../marked_boundaries.py)

### Quantum background dependence

For a full background-dependent effective action satisfying the quantum master equation, translation in $t^r\eta$ commutes gradedly with the free differential, regularized BV Laplacian and propagator contraction. Differentiating its quantum master equation and its renormalization group equation gives a closed, compatible family of quantum integrated observables. This assertion applies to formal directional derivatives of local functionals about compactly supported perturbations on flat space. The constant background may be treated as a coefficient in that local expansion. [Constant-translation proof, Lemma 3.8.1](https://arxiv.org/pdf/1201.4501) · [Effective closed-string BV equation, §4](https://arxiv.org/pdf/1505.06703)

Identifying its classical term with this specific $F_r$ requires identification of the full genus-zero action and of the background direction through the equivalence. Costello–Li's open-closed quantization produces a classical closed action; equivalence to the explicit descendant potential in general dimension is a conjecture in §8 of their paper. The pure-primary cubic comparison alone does not fix the $t^r\eta$ derivative. A bare positive two-loop primary value and its exclusion from the marked primitive image do not decide the corresponding full quantum lifting question. [Scope of the classical-action comparison](https://arxiv.org/pdf/1505.06703#page=84)

### One nonconstant scalar descendant at genus one

The double primary derivative of the genus-zero interaction at a hexagon bulk vertex is multiplication by the even polyvector superfield solving $U=\sum_{k=0}^4\alpha_k U^k/k!$. This follows by applying Lagrange inversion to the coefficient $(l-3)!/(l-2)!$ left when two primary slots are removed from the valence-$l$ vertex. Its scalar and vector components are
$$
u=\sum_{k=0}^4\phi_k\frac{u^k}{k!},
\qquad
w=\frac{\sum_{j=0}^3\mu_j u^j/j!}{1-\sum_{k=1}^4\phi_k u^{k-1}/(k-1)!}.
$$
For a scalar mark the rank count allows only this scalar or vector component at each bulk vertex. Thus the previous five-slot hexagon applies to four copies of $u$ and one of $w$. The marked spectator coefficient is the ordinary second primary derivative of the scalar integrand.

Restrict external fields to $\phi$, the level-one scalar $b=\phi_1$, and $\mu=\mu_0$. Put $h=(1-b)^{-1}$. Then $u=h\phi$, $w=h\mu$ and the marked coefficient for $F_1$ is $hu$. On scalar coordinate functions the nonlinear differential gives $Q_0\phi=\partial_i(u\mu^i)$, while the free dual differential gives $Q_0 b=-\partial_i\mu^i$. Consequently
$$
Q_0 u=w^i\partial_i u,\qquad Q_0 h=w^i\partial_i h-h\partial_iw^i.
$$
All evaluations below use rational wave momenta and the reduced heat convention.

For spatially constant $b$, the four scalar slots contribute $h^4$, the vector slot $h$, and the primary spectator $h^2$. Thus on arbitrary $\phi,\mu$ the result is $h^7\int\Omega\,\det(\partial_i\partial_j\phi)\partial_l\mu^l/4$. The local correction that solves this entire constant-$b$ sector is
$$
\int\Omega\,\frac{(1-b)^{-6}-1}{24}\det(\partial_i\partial_j\phi).
$$
Its coefficient has derivative $h^7/4$ in $b$. For a constant coefficient the first variation of the integral of a Hessian determinant is zero by integration by parts and the divergence-free cofactor. This leaves only the displayed free descendant variation.

For one nonconstant $b$ write $x=h-1$ and first use the independent variables $x,u,w$. Take five scalar $u$ waves of momenta $p_i=e^i$, one vector $w$ wave of momentum $k=e^1$ and polarization $v=e_2$, and one $x$ wave of momentum $q=-k-\sum_i p_i$. The vector is transverse. The coefficient linear in $x$ in the marked factor $hu$ is $xu$. Sum over the choice $i$ of the scalar at this marked point. The reduced amplitude is
$$
\frac1{24}\sum_i\det(p_1,\ldots,\widehat p_i,\ldots,p_5,k)^2
                 (k-q-p_i)(v)=\frac1{24};
$$
only $i=1$ has a nonzero determinant. The sums here are over the five scalar labels.

For an equivariant local marked change, terms quadratic in $x$ could contribute linearly through the free part $-\partial_iw^i$ of $Q_0x$ only, which vanishes. At first order $x$ the candidate primitive has five $u$ inputs and one $x$ input of total derivative degree ten. Six conserved covectors leave five independent covectors; thus every holomorphic $SL_5$ invariant at derivative degree ten is a multiple of their determinant squared. A representative is $\int\Omega\,x\det(\partial_i\partial_j u)$. On these transverse inputs its variation is
$$
q(v)\det(p_1,\ldots,p_5)^2
+\sum_i p_i(v)\det(p_1,\ldots,p_i+k,\ldots,p_5)^2=-1+1=0 .
$$
An $x$-independent term cannot generate $x$ since $Q_0u=w^i\partial_i u$. Higher scalar descendants have free partners of positive vector descendant level that vanish on this restriction. Positive polyvector rank in the primitive cannot yield total output rank one through the nonlinear variation. Antiholomorphic jets are excluded in the holomorphic jet model.

Independently, in the original variables $\phi,b,\mu$ with these same individual momenta, the linear-$b$ terms in the marked coefficient, bulk scalar and bulk vector are respectively $2b\phi,b\phi,b\mu$. Their reduced values are $1/12,-1/24,-1/4$, hence the genuine marked wheel gives $-5/24$. The necessary primary-projection correction $\int\Omega\,b\det(\partial_i\partial_j\phi)/4$ contributes $+1/4$ in these variables; the result is again $1/24$. The comparison rules out a marked local solution on this one-descendant sector for the fixed bare differential. A genus-one bulk action term acting by its classical bracket on $F_1$ has not been included. [Direct sums in both sets of fields](one_descendant_hexagon.py)

### The genus-one bulk term on the one-descendant sector

In reduced heat units the leading bulk correction can be chosen as
$$
\frac1{24}\int\mathrm{Tr}\,\alpha_1\det(\partial_i\partial_j\alpha_0).
$$
Here the fields are the full even polyvector superfields at the indicated descendant levels. The determinant uses their commutative superfield product, and the trace selects the top polyvector and top Dolbeault form. Its six-linear polarized symbol with one distinguished descendant is the trace product times the squared determinant of the five primary momenta. On conserved momenta any five of the six give the same square. The leading bulk QME fixes the displayed choice up to free cocycles.

Indeed the ungrounded bulk hexagon pairs one external polyvector contraction with its vertex potential, with a common factor equal to the squared momentum determinant. Averaging the grounded Green matrix over the six placements leaves only the diagonal minus the off-diagonal entry, since the covectors sum to zero. The common vertex potential also drops out since a simultaneous odd translation of the six factors leaves their Berezin integral invariant. The sum over labeled orders, the dihedral symmetry factor and the simplex volume give $6!/(12\cdot 12\cdot 5!)=1/24$. The contraction at slot $i$ is precisely its holomorphic divergence. The free differential of the displayed functional supplies all six distinguished choices.

The relative sign can be fixed in the full polyvector space before restriction. Put each odd polarization parameter on the left and extract the coefficient with all parameters to the left of the odd polyvector coordinates. With $\partial_\Omega=\partial_{z^i}\partial_{\theta_i}$, the free action on the even generating field is $Q_0\alpha_1=\partial_\Omega\alpha_0$; on a wave $\beta\theta(v)$ with odd $\beta$ this gives $-\beta k(v)$. Translate the interaction by the odd parameter $\epsilon$ in direction $t\eta$, in the order $\epsilon\beta$. At the translated hexagon vertex the primary Hessian contributes $\epsilon\eta\phi$ with scalar momentum $p$. Taking the divergence of this slot contributes $-p(v)$ to the trace polynomial. Taking the divergence of the vector contributes $+k(v)$. The resulting trace polynomial has coefficient $(k-p)(v)$, in agreement with the marked graph. The bulk master-equation defect with one odd insertion has coefficient *minus* the marked differential: the BV operation is odd and $\{I,\epsilon F\}=-\epsilon\{I,F\}$ for the even bulk action $I$. Thus the free variation of the positive bulk term above cancels the six-primary bulk defect. The same sign fixes its induced scalar variation to be the positive divergence of its vector variational derivative, as for the classical cubic. [Grounded Green function and exterior factor](../hexagon.py)

Keep five scalar $u$ momenta $p_i$, one $w$ of momentum $k$ and polarization $v$, and one $x$ of momentum $q=-k-\sum_i p_i$. In this sector $\alpha_1$ supplies $x$ in the bulk correction and its remaining inputs can be taken at $x=0$. The primary derivative of $F_1$ is $u^2/2$. Integration by parts in its bracket with the bulk correction puts a derivative on this scalar pair. The reduced bulk contribution is
$$
-\frac1{24}\sum_{i<j}(p_i+p_j)(v)
 \det(p_i+p_j,k,p_1,\ldots,\widehat p_i,\ldots,\widehat p_j,\ldots,p_5)^2 .
$$
On the preceding transverse probe only the pair $(1,2)$ survives. This gives $-1/24$ and cancels the entire residual.

For a general probe with $\Delta=\det(p_1,\ldots,p_5)\ne0$, write $k=\sum_i a_i p_i$, $c_i=p_i(v)$ and $d=k(v)$. In units $\Delta^2/24$, the genuine hexagon, the variation of $\int\Omega\, x\det(\partial_i\partial_j u)/4$, and the bulk bracket respectively give
$$
\begin{aligned}
 &(\textstyle\sum_i c_i+2d)\sum_i a_i^2-\sum_i c_i a_i^2,\\
 &6\big[q(v)-d+\sum_i c_i(1+a_i)^2\big]=6\sum_i c_i a_i^2,\\
 &-\sum_{i<j}(c_i+c_j)(a_i-a_j)^2
 =-\big(\textstyle\sum_i c_i\big)\sum_i a_i^2-5\sum_i c_i a_i^2+2d\sum_i a_i .
\end{aligned}
$$
What remains is $2d\sum_i(a_i^2+a_i)$. Let $M_i(k)$ be the determinant of the five scalar rows with row $i$ replaced by $k$. A local marked correction with two $x$ slots of momenta $k,q$ has the polarized symbol
$$
-\frac1{12}\sum_i M_i(k)M_i(q)
=\frac{\Delta^2}{12}\sum_i(a_i^2+a_i).
$$
This polynomial is symmetric in the two descendant waves and in the five scalar waves, so it defines a translation-invariant local functional by polarization. Its free variation replaces one $x$ by $-\partial_i w^i$. This cancels the last term for all polarizations. Polynomiality extends the identity to $\Delta=0$. Thus the primary and one-scalar-descendant equations for the integrated cubic admit simultaneous local solutions with this bulk choice. The higher arities follow from the free jet complex below. [Exact Green factor, full polarization and singular momenta](bulk_descendant_contact.py)

### Direct all-arity genus-one master equation on the fivefold

Let $F_0$ be the primary Euler derivative of the classical action. Use the even superfields $A_j$, primary Hamiltonian variation $Q_0A_0=-\partial_\Omega F_0$, and free variations $Q_0A_{j+1}=\partial_\Omega A_j$ in the positive classical descendant convention. Write $T(x)=\sum_{j=0}^4 A_jx^j/j!$, $U=T(U)$, and $H=(1-T'(U))^{-1}$. The two field derivatives required here are
$$
\frac{\partial F_0}{\partial A_i}=\frac{U^{i+1}}{(i+1)!},\qquad
\frac{\partial U}{\partial A_j}=H\frac{U^j}{j!}.
$$
For even superfields $X,Y$, their Schouten bracket in this calculation is $\{X,Y\}=\partial_\Omega(XY)-(\partial_\Omega X)Y-X\partial_\Omega Y$, with left odd derivatives in $\partial_\Omega=\partial_i\partial_{\theta_i}$. Expanding $Q_0U$ by the ordinary chain rule and $\partial_\Omega F_0$ by the second-order chain rule cancels all terms linear in $\partial_\Omega A_j$. The missing terminal free term is $\partial_\Omega A_4=0$ since $A_4$ is scalar. The remaining terms give
$$
Q_0U=-\tfrac12\{U,U\}.
$$
Differentiate in a spatially constant scalar shift of $A_0$. This shift commutes with spatial and odd derivatives. Its commutator with $Q_0$ sends $A_0$ to $-\partial_\Omega U$ and annihilates $A_1,\ldots,A_4$. Acting on $U$ gives $-H\partial_\Omega U$. Since the constant derivative of $U$ is $H$,
$$
Q_0H=H\partial_\Omega U-\{U,H\},
\qquad Q_0\log H=\partial_\Omega U-\{U,\log H\}.
$$
Put $M_{ij}=\partial_i\partial_j U$ and let $C^{ij}$ denote its cofactor. The explicit coordinate bracket on two even fields is $\{X,Y\}=(\partial_i X)\partial_{\theta_i}Y+(\partial_{\theta_i}X)\partial_i Y$. The identity $C^{ij}M_{ik}=\det(M)\delta^j_k$ gives the pointwise polynomial equality
$$
C^{ij}\{\partial_i U,\partial_j U\}=2\det(M)\partial_\Omega U .
$$
Differentiating $Q_0U$ twice in space therefore gives, for the even density $W=(\log H)\det M$,
$$
Q_0W=(\partial_\Omega U)\det M-\{U,W\}-2(\partial_\Omega U)W .
$$
Integration by parts in one spatial and one odd derivative gives $\int\mathrm{Tr}\,U\partial_\Omega W=\int\mathrm{Tr}\,(\partial_\Omega U)W$. Since $\int\mathrm{Tr}\,\partial_\Omega(UW)=0$, the integrated bracket is $-2\int\mathrm{Tr}\,(\partial_\Omega U)W$. It follows that
$$
Q_0\,\frac1{24}\int\mathrm{Tr}\,(\log H)\det M
=\frac1{24}\int\mathrm{Tr}\,(\partial_\Omega U)\det M.
$$

Each of the six vertices of the bare loop contributes the two-primary derivative of the classical action, exactly $U$. Sum the six-point voltage with these six even superfields. The integrated diagonal minus off-diagonal Green response supplies $1/24$, and its bulk master-equation orientation is the negative of the right-hand side, fixed above on a longitudinal six-wave test. A bridge kills the heat cofactor; the other genus-one diagrams are exactly the classical one-edge term already in $Q_0$. The coincident free heat trace is zero on the holomorphic jets. This proves the compact action solves the genus-one bulk master equation at every arity directly from graphs. The minimal restriction is preserved: $\partial_\Omega F_0$ has rank at most four in five odd variables, and $\partial_\Omega A_j$ has rank at most $3-j$ as required at level $j+1$. The constant shift of $A_0$ also preserves this subspace. All products and integration by parts were in the full exterior algebra before restricting the external superfields. [Six-point graph orientation](../hexagon.py) · [Left Grassmann derivatives, cofactor and signed trace tests](../elliptic_superfield.py)

### Seven-wave bulk graph equation with moving descendants

Use the coefficient $A_k$ of $\zeta^k$ in Li–Wang's field, with positive classical descendant weights. On seven conserved waves let $d$ be the label of the single level-one field in the master-equation defect. For distinct external labels $i,j$, put $B_{ij}=\det(p_{\widehat{i,j}})^2$, the square of the determinant of the remaining five momenta in label order, and put $B_{ii}=0$. Write $B_i=\sum_{\substack{j\\\mathrm{wave\ labels}}}B_{ij}$. The loop is a hexagon with an external pair $d,j$ at the quartic vertex. The tree has a classical cubic with a pair of primary waves and the six-wave genus-one vertex at the other end. Measured from that other end, its one-edge voltage is the momentum sum of its classical pair. In units of $1/24$ the two potentials at slot $r$ are
$$
\begin{array}{c|cc}
 &\text{six quartic hexagons}&\text{classical trees}\\ \hline
 r=d & B_d p_d+\displaystyle\sum_{\substack{j\\\mathrm{wave\ labels}}}B_{dj}p_j &0\\
 r\ne d &B_d p_r+B_{dr}p_d&
 \displaystyle\sum_{\substack{j\ne d,r\\\mathrm{primary\ labels}}}B_{rj}(p_r+p_j).
\end{array}
$$
The averaged hexagon voltage at a fused event uses the sum of both wave momenta for either exterior polarization there. The squared momentum cofactor is $B_{dj}$. The diagonal minus off-diagonal Green integral, including the $720$ labeled orders and dihedral factor $12$, is $1/24$.

The kernel of a rank-five seven-wave momentum matrix contains the all-ones vector by conservation. Complete it to a basis of its two-dimensional kernel by a vector with components $h_j$ in wave-label order. The signed complementary five-minors are Plücker coordinates of this kernel; hence $B_{ij}$ is a common constant times $(h_i-h_j)^2$. Since both $\sum_{\substack{j\\\mathrm{wave\ labels}}}p_j$ and $\sum_{\substack{j\\\mathrm{wave\ labels}}}h_jp_j$ vanish, $\sum_{\substack{j\\\mathrm{wave\ labels}}}B_{rj}p_j$ is independent of $r$. Subtracting this common exterior translation gives the free divergences of
$$
V_{1;2_d}=\frac{B_d}{24},\qquad
V_{1;1_d,1_r}=\frac{B_d+B_r-B_{dr}}{24}\quad(r\ne d).
$$
The formulas extend polynomially to deficient rank. Within the trace ansatz these coefficients are unique: the difference between two solutions supplies a multiple of each wave's own momentum and the same covector at every slot. Two nonparallel waves force that covector to vanish.

At higher arity, place the external waves in six nonempty groups around the loop. In each group $B$ require $\sum_{\substack{i\in B\\\mathrm{wave\ labels}}} k_i=|B|-1$, and use the classical multinomial in those levels. The potential at a wave is its group momentum, times the square of any five group momenta and the common $1/24$. For a one-edge tree, a classical group $S$ instead has total level $|S|-2$. Sew it to the lower-arity genus-one vertex evaluated on the complementary waves and the level-zero wave of momentum $\sum_{\substack{i\in S\\\mathrm{wave\ labels}}}p_i$. Its voltage has this value on $S$ and zero on the complement. All cycle graphs with a bridge have zero heat cofactor. Recovering the individual free divergences from differences of these potentials reproduces the independent moduli-space partition sum in $6,49,288,1485$ labeled tests at respectively six, seven, eight and nine waves. [Graph-only recursion and exact comparisons](../genus_one_graph_ward.py)

### The complete one-loop differential of a scalar point polynomial

Use the fixed normal-ordered coordinate $\phi=(A_0)_{\rm scalar}$ and a polynomial $f(\phi)$ at the marked point. In the existing equation $U=\sum_{j=0}^4A_j U^j/j!$ denote its scalar and vector components by $u,w^a$. The superscript $a$ is the tangent-vector slot. Let $C^{ij}$ be the cofactor of the ordinary five-coordinate Hessian $\partial_i\partial_j u$. With the one-letter term computed from the Euler derivative of the genus-one bulk action, the full genus-one coefficient is
$$
Q_1 f(\phi)=
f'(\phi)\,Q_1\phi-
\frac{f''(\phi)}{24(2\pi i)^5}
\left[
 \partial_a\!\left(C^{ij}\partial_i\partial_j w^a\right)
 + C^{ij}\partial_i\partial_j\partial_a w^a
\right].
$$
The string expansion multiplies $Q_1$ by $\lambda^2$. A one-edge graph consumes one letter and gives the first term. A positive geometric loop consumes two letters on the marked hexagon. Its five other vertices are the two-primary derivatives of the classical interaction, with external fields exactly $U$. Four scalar components give the cofactor; a single vector component gives $\partial_i\partial_j w^a$. In the rooted hexagon the Green average pairs the vector with the total of all five group momenta plus the momentum of its own group, accounting separately for the two derivatives in brackets. The marked Hessian contributes $f''(\phi)$, and the coincident free heat contact is zero on these holomorphic jets.

In the one-edge term the superfield operation in reduced units is $-\partial_\Omega$ on the Euler derivative with respect to $A_0$. On a vector with odd coefficient on the left this is the positive ordinary divergence. For root $e^1$, vector momentum $e^2$ with polarization $e_1$, and scalar momenta $e^3,e^4,e^5,e^1+e^2,-2e^1-2e^2-e^3-e^4-e^5$, the momentum at that divergence evaluates to $-1$. A level-two wave at $e^1+e^2$ gives reduced tree $-1/2$. Assigning level one to each of the final two waves gives $-7/12$. The genus-one heat factor $-1/(2\pi i)^5$ fixes both signs. Independent finite-jet products of the compact action, its variational derivative, the rooted hexagon and the bulk graph recursion reproduce these matrix elements and the moving-descendant second-order term. [Taylor and divergence replay](../elliptic_hessian.py)

For the protected point classes obtained by differentiating the action, write $E_r^{(1)}$ for the reduced genus-one Euler derivative with respect to $A_r$. With $H=\partial U/\partial A_0$, varying the implicit equation gives $\partial U/\partial A_r=HU^r/r!$. The explicit dependence of $\log H$ on $A_r$ at fixed $U$ supplies the second term in
$$
E_r^{(1)}=\frac{U^r}{r!}E_0^{(1)}
 +\frac{HU^{r-1}}{24(r-1)!}\det(\partial_i\partial_jU),
 \qquad 1\le r\le4.
$$
On primary superfields alone, $E_0^{(1)}=0$ and $E_r^{(1)}=A_0^{r-1}\det(\partial_i\partial_jA_0)/(24(r-1)!)$. In the rank-five component take primary waves of momenta $e^1,\ldots,e^5$. Give the first wave polarization $e_2\wedge e_3\wedge e_4\wedge e_5$, the fifth polarization $e_1$, and include $r-1$ additional labeled constant scalars. Their exterior trace is $+1$ and the Hessian determinant gives one. Thus all four differentiated classes have the same polarized coefficient $1/24$ on these probes in the $A_r$ source convention. The compact formula also includes their moving-descendant components. [Euler derivatives at each source level](../elliptic_hessian.py)

### Eliminating higher descendant equations in the holomorphic jet complex

Fix $N$ labeled external legs and fix $m_i$, the sum of polyvector rank and descendant level on leg $i$. The minimal theory has $0\le m_i\le4$. In the polynomial ring of its five momentum components the free dual differential is exterior multiplication by the leg covector $p_i$. Its complex consists of the exterior powers of the coordinate cotangent space from $0$ through $m_i$. Exterior power $m_i$ is at descendant level zero. The full five-variable Koszul complex is the tensor product of the five injective one-variable multiplication complexes; it is exact below exterior power five. The truncated complex therefore resolves its primary cokernel in exactly $m_i$ steps. That cokernel embeds into the next exterior power by exterior multiplication by $p_i$. It is locally free where a component of $p_i$ is invertible, by contraction against a vector pairing to one with $p_i$.

Tensor these complexes on the $N$ separate legs and impose the five linear momentum-conservation equations from total derivatives. The resulting free modules are over a polynomial ring in $5(N-1)$ variables. The length of the complex is $s=\sum_i m_i\le4N$. Its positive homology is supported at the common origin. Indeed, away from it, select a leg with an invertible momentum component and eliminate that leg momentum by conservation. That leg complex splits into a locally free primary module and contractible pairs. The other legs have independent momenta, so their tensor product is exact in positive descendant level. Their primary cokernels embed into free modules, and their tensor product over the ground field preserves these injections.

For $N\ge5$ we have $s\le5(N-1)$. Descend through the complex from its top free module, assuming exactness through positions strictly above $j>0$. The quotient of the free module in position $j$ by the incoming image then has a free resolution of length at most $s-j<5(N-1)$. At the origin no such module contains a nonzero finite-length submodule. To see this, the Koszul resolution of the residue field gives zero Ext into a free module in degrees below $5(N-1)$. Apply its long exact Ext sequence along the indicated resolution. The group of maps from the residue field to the quotient is zero, whereas any nonzero finite-length submodule contains a copy of the residue field. Since the homology at $j$ is such a submodule, it vanishes. This proves exactness at every positive total descendant level. For $N\ge6$ we also have $s<5(N-1)$ at position zero. The same argument excludes primary torsion at the origin, and the preceding injections exclude it away from the origin. The primary quotient is therefore torsion-free. Polynomial primary symbols can be detected by evaluation on generic conserved waves with each primary polyvector annihilated by its own momentum. Slot permutation symmetries and $SL_5$ invariance can be restored on primitives by averaging over the finite group and a compact real form at fixed polynomial degree.

In the one-loop bulk problem the defect at arity $N$ has external descendant level $N-6$. At $N=6$ the explicit free primitive above cancels the defect in all ranks. Having solved the equation at lower arity, the arity-$N$ residual is free-closed by the BV consistency identity: for the odd square-zero effective BV operator applied to the exponential of the action, applying it a second time gives zero, and the lowest uncanceled coefficient is killed by the free differential. At $N>6$ the exactness just proved supplies a local primitive with descendant level $N-5$. Its nonlinear classical variation only enters higher arities. This gives a formal genus-one bulk solution at every arity in the translation-invariant holomorphic jet model. For the integrated cubic the marked defect likewise has level $N-6$, including the bulk bracket. Its arity-six primary output consists of five scalars and one vector and is canceled above. The marked consistency identity follows by differentiating the bulk identity in a square-zero parameter for the classical mark. The same induction therefore gives a formal genus-one lift of the completed integrated cubic through every descendant level in this model. No equality with an open-closed effective action or with the untwisted IIB action is needed for this induction. At genus two the eight-primary bulk and cubic identities are solved by the seven-slot symbols $-P\,\mathrm{Tr}$ and $-9P$. The same induction applies above arity eight. [Genus-two trace and cubic calculation](../README.md#the-genus-two-bulk-action-and-integrated-cubic) [Flat effective BV framework, §4](https://arxiv.org/pdf/1505.06703)

The same genus-one lifting argument applies to every $F_r$ with $r\ge1$. There are $r+4$ primary scalar momenta $p_i$ on the first marked equation, together with one vector $(k,v)$. Put $c_i=p_i(v)$ and $d=k(v)$, so $\sum_i c_i=-d$. For a five-element subset $B$ of the scalar slots write $D_B=\det(p_B)$, using increasing order. For a four-element subset $A$ put $H_A=\det(p_A,k)^2$, and put $F=\sum_B D_B^2$. The scalar second primary derivative at this root is $hu^r/r!$. At zero descendant order the marked hexagon times 24 is
$$
\sum_A H_A\big(2d+\sum_{i\in A}c_i\big).
$$
Take the primary marked correction $\int\Omega\,u^{r-1}\det(\partial_i\partial_j u)/[24(r-1)!]$. Its symbol times 24 is $F$. At $Q_0u=w^i\partial_i u$, fuse $k$ with each scalar $i$ in turn, weighted by $c_i$. If $i\in B$ the determinant changes from $D_B$ to $D_B+M_i(k)$, with $M_i(k)$ the row-replacement determinant on $B$. Cramer's identity gives $\sum_{i\in B}c_i M_i(k)=dD_B$. Therefore the variation times 24 is
$$
-dF+2dF+\sum_A H_A\sum_{i\notin A}c_i
=dF+\sum_A H_A\big(-d-\sum_{i\in A}c_i\big).
$$
The residual $d(F+\sum_A H_A)/24$ is the negative free variation of a single-$x$ primitive with polarized symbol $(F+\sum_A H_A)/24$ and descendant momentum $k$. The numerator sums the squared determinants of all five-element subsets of the scalar and descendant momenta, so it defines a local polynomial. The bulk correction retains its descendant in a bracket on exclusively primary fields. Higher arities of the marked equation have positive descendant level, and the preceding argument gives a genus-one lift of every $F_r$.

### Finite transverse tests of the genus-two bulk and marked tensors

Take eight conserved primary momenta and a slot covector $A_i$ at each external leg. The genus-two heat and odd minors each have degree ten, and the potential adds one power of momentum. Hence the integrated tensor is an $S_8$-equivariant polynomial array of covectors of degree 21. It transforms by the fourth determinant power times the defining covector action. The signed slot contraction of these covectors with the polyvectors is insensitive to a common shift of all $A_i$ because a product of total polyvector rank six vanishes in five dimensions. A shift of slot $i$ by a multiple of $p_i$ also vanishes on transverse primary fields.

A converse follows on generic momenta. Put $\iota_{p_l}\eta$ at one fixed leg $l$, two transverse vectors $v_i,v_j$ at distinct legs $i,j$, and scalars elsewhere. The odd determinant evaluates on this input to
$$
 (A_i-A_l)(v_i)p_l(v_j)-(A_j-A_l)(v_j)p_l(v_i).
$$
If all such tests vanish, choose $v_j$ pairing to one with $p_l$. Then $(A_i-A_l)$ restricts to a fixed multiple of $p_l$ on $\ker p_i$, with the same multiple at all legs by the pair equations. Thus $A_i$ is a common covector plus a multiple of $p_i$ at each leg. This also proves sufficiency of the chosen test field ranks.

For the dimension count, the seven independent momentum slots transform as the sum-zero representation $W$ of $S_8$. The degree-21 Cauchy decomposition on the matrix of five rows and seven columns selects the Schur functor of partition $(5,4,4,4,4)$ on $W$: subtracting four determinant powers leaves one defining covector. Scalar coefficients of a longitudinal term have partition $(4,4,4,4,4)$. The eight labeled output slots transform as the permutation representation. The multiplicities in the two Schur functors are
$$
\begin{array}{c|rr}
 &\text{trivial}&\text{permutation}\cr
 (5,4,4,4,4)&6&32\cr
 (4,4,4,4,4)&4&12 .
\end{array}
$$
For cycle lengths $L$ the complete symmetric coefficients $h_j$ are extracted from $(1-t)/\prod_L(1-t^L)$. The Jacobi–Trudi determinant $\det(h_{\lambda_i-i+j})$ is multiplied by 1 or by the number of fixed points and averaged over $S_8$. The replay performs the rational character sums. The common shifts and longitudinal arrays intersect trivially since two generic covectors $p_i,p_j$ are independent. The physical quotient has dimension $32-6-12=14$. For the stronger condition that $10A_i-\sum_s A_s$ is longitudinal, the map on slot space has eigenvalues 2 and 10. Its quotient has dimension $32-12=20$.

The first certificate uses the five coordinate momenta, $(1,-1,1,0,-1)$, $(0,1,1,-1,-1)$, and their balancing eighth momentum. Fourteen transverse rows at these momenta suffice for the physical quotient. The twenty-row certificate uses eighteen rows here and two at a second fixture, with sixth and seventh momenta $(1,2,-1,0,-1)$ and $(-2,0,1,-1,2)$. For independence, symmetrize the covector of leg 1 placed at leg 0 times four five-row determinants over all $8!$ permutations. The replay gives the four five-subset labels for every such column and every chosen evaluation. The square matrices have determinants $96$ and $927$ modulo the prime 1009, respectively. Each is therefore injective on the entire indicated rational quotient. [Fixtures, exact orbit sums and both certificates](genus_two_probes.py)

### The integrated quartic at genus two

For $F_2$, take the eight labeled primary scalar momenta and one vector of momentum $k$ and polarization $v$. For a subset $B$ of five scalar labels let $D_B$ be their ordered determinant. On eight conserved scalar slots the signed genus-two marked primitive has polarized symbol
$$
\sum_{d=2}^5 c_d
\sum_{\substack{\{B,C\}\\ |B|=|C|=5,\ |B\cap C|=d}}
 D_B^2D_C^2,
\qquad
(c_2,c_3,c_4,c_5)=
\left(\frac{43}{207360},\frac1{69120},-\frac1{311040},\frac1{124416}\right).
$$
The inner sum runs over unordered pairs with repetition. The physical multiplier is $\lambda^4/(2\pi i)^{10}$. Specializing any scalar momentum to zero gives the signed cubic primitive $-9P$, where $P$ is the seven-leg bulk quartic with the coefficients in the [eight-primary calculation](../bulk_eight.py).

The bivalent theta contribution is the sum over pairs of scalar spectators. Fuse each pair at the root and evaluate $\Phi_v-\Phi_{\rm root}$ from the eight balanced bivalent sources. At a trivalent mark there is one scalar spectator. Ground that tip, omit the spectator from the eight bivalent inputs, and evaluate $2\Phi_v$ on the resulting generally unbalanced eight covectors. The tip carries the balancing charge. The factor two fixes the tip of the theta. In the genus-two bulk tree bracket fuse each triple of scalar spectators. The contribution is their total momentum paired with $v$, times $P$ evaluated on the fused triple, the other five scalars and the vector. These signs use the signed bulk action $-P$.

The genus-one primary marked correction is $\int\Omega\,u\det(\partial_i\partial_j u)/24$. In a hexagon rooted at this vertex the four scalar spectator momenta and the two incident currents sum to zero with their incoming signs. Every one of their six five-minors has the same square. For ordered bulk momenta $z_1,\ldots,z_5$, with the vector at $j$, put $T_i=\sum_{e=i}^5 t_e$ in the unit resistance simplex. The initial current is $-\sum_iT_i z_i$. For each set of four scalar spectators $S$, put $L_i=\det(p_S,z_i)$. The reduced genus-one insertion into the genus-two equation is
$$
\frac{\det(z_1,\ldots,z_5)^2}{8}\,
\operatorname{Av}_{S_5}\operatorname{E}_{t}
 \left[(\sum_i T_i L_i)^2
 \sum_l(T_{\max(j,l)}-T_jT_l)z_l(v)\right],
$$
summed over the spectator choices. The average permutes the five active labels. The factor is $6/(24\cdot2)$, including the root-preserving graph symmetry; the simplex volume cancels the number of permutations. For descending positions $i_1\ge\cdots\ge i_m$ the mean of the tails is $\prod_{h=1}^m(6-i_h+h-1)/(6+h-1)$ by elementary simplex integration. Both the genuine two-loop graph and this genus-one insertion carry the same normalized multiplier.

To evaluate the primary classical variation of each candidate marked primitive, fuse the vector momentum with each scalar and multiply by that scalar's pairing with $v$, in the variable $u$ with $Q_0u=w^i\partial_i u$. Adding the two theta contributions, the bulk tree contact, the genus-one marked insertion, and the signed primitive above gives zero on all 108 transverse rows in the replay. For nine slots the relevant covector and longitudinal scalar multiplicities are 75 and 24, so the whole transverse space has dimension 51. The multiplicity of scalar invariants at eight conserved slots is four. The four displayed columns of classical variations have rank four. The test rows on three integer fixtures contain a square 51-row matrix of determinant $714\pmod{1009}$, computed on $S_9$-averaged products of four minors times a slot momentum. Consequently the cancellation holds as a polynomial identity for arbitrary conserved momenta. The remaining longitudinal covector equals $k$ times a polynomial: the five components of a generic $k$ have no common divisor. This gives a local free primitive with one scalar descendant; the positive-descendant exactness above supplies the remaining arities. Thus the integrated quartic admits a genus-two lift in this holomorphic jet model. [Three exact graph fixtures](higher_genus_two.py) · [Character computation and 51-row certificate](quartic_unisolvence.py)

### Finite reduction for every scalar background at genus two

At the primary equation for $F_r$ put $n=r+6$. There are $n$ scalar waves and one vector. The bivalent theta fuses $r$ spectators at its root; the trivalent theta fuses $r-1$ there. The bulk tree fuses $r+1$ spectator momenta into $P$. For the one-loop marked insertion, take four active scalars and one vector, and $r+2$ spectators at the mark. The polarized genus-one primary symbol times 24 is the sum of the squared five-minors of these $r+2$ spectators and the two signed incident currents. The currents differ up to incoming signs by the fixed sum of the five active momenta, so even a minor involving both currents is at most linear in the resistance tails. Thus the integration uses the same fourth-order moments as for $F_2$.

Inserting a zero-momentum scalar specializes every term in this equation to the corresponding one at $r-1$. A zero at an unmarked bivalent vertex kills its current determinant. A zero spectator removes a single undifferentiated scalar with unit polarized coefficient. A minor using a zero row vanishes. In the genus-two bulk term a zero outside the spectator group also vanishes: on six balanced momenta and one zero the centered seven-entry relation is a multiple of $(-1,-1,-1,-1,-1,-1,6)$. The two quartics in $P$ have values $1764$ and $1302$ on this vector, and their signed coefficients obey $-31\cdot1764+42\cdot1302=0$. Advection of a scalar by the vector commutes with zero-slot restriction since the momentum pairing of that slot vanishes.

Eliminate the vector momentum as $k=-\sum_i p_i$. A symmetric covariant of degree 21 in the independent scalar momenta, transforming with four determinant powers, is a sum of expressions with four alternating five-row determinants and one uncontracted covector. In any term symmetric in the scalar labels, each single-occurrence scalar label must occur in a different determinant or in the free slot: interchanging two single-occurrence labels in the same determinant kills its symmetrization. At most five labels therefore occur once; all other active labels occur at least twice. There are at most $(21+5)/2=13$ active labels. For scalar coefficients of degree 20 the corresponding bound is $(20+4)/2=12$.

Symmetrized monomials that sum each choice of distinct labels once extend to any larger number of scalar waves. Restriction to a zero scalar deletes every summand using that label. On independent scalar momenta, orbit sums of ordinary component monomials with fixed nonzero exponent multisets are linearly independent. Hence restriction from more than 13 scalar waves is an isomorphism on this covector space, and similarly above 12 on the scalar longitudinal coefficients. Since $k$ also restricts compatibly, the transverse quotient stabilizes. A primary scalar counterterm on 13 *conserved* slots likewise extends: lift its symmetric polynomial representative to larger arities by the same distinct-label sums, then impose conservation. Different such lifts give the same restricted primitive and agree on its class wherever the corresponding restriction is injective.

| Total primary legs | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Transverse marked tensor dimension | 20 | 51 | 94 | 126 | 142 | 146 | 147 | 147 |
| Scalar primitive dimension on one fewer leg | 2 | 4 | 7 | 10 | 12 | 13 | 13 | 13 |

The dimensions follow by the same cycle character formula as above. Consequently it suffices to solve the single finite primary equation for $F_7$, with a 147-dimensional transverse target and 13 scalar primitive symbols. Restricting that solution gives every lower $r$. Extending its primitive gives every higher $r$ by the stable injectivity just proved. Each primary solution then lifts through all positive descendant arities. After the cubic and quartic, the next equation is $F_3$: after a lift of the quartic primitive the transverse residual vanishes when one scalar is zero. The new target dimension is $94-51=43$ and the new scalar primitive dimension is $7-4=3$. The restriction maps here are surjective by extension of symmetric orbit sums, including after imposing conservation and taking the longitudinal quotient. The quintic calculation below uses this difference of dimensions.

### The integrated quintic at genus two

For $F_3$, extend the four quartic overlap sums to its nine scalar marked slots. Three additional signed orbit sums suffice. Number the five-element subsets of $\{0,\ldots,8\}$ lexicographically starting at zero. A four-tuple in the table defines the product of the four corresponding ordered determinants. Sum its signed orbit under all permutations of the nine slots, counting each distinct monomial once with its orientation sign. The three products use all nine slots, so their restriction with one zero scalar is zero.

| Four subset numbers | Signed coefficient |
|---|---:|
| $(17,26,85,104)$ | $17/3110400$ |
| $(16,32,103,119)$ | $97/933120$ |
| $(50,102,106,119)$ | $-1/622080$ |

The coefficients multiply the same positive normalized factor $\lambda^4/(2\pi i)^{10}$. In the remaining genus-one marked hexagon the root has five scalar spectators. For four fixed spectator rows define $L_i=\det(p_S,z_i)$. The two choices of one internal row contribute $(\sum_i T_iL_i)^2+(\sum_i(1-T_i)L_i)^2$. For three fixed spectator rows set $L_i=\det(p_S,z_i,\sum_j z_j)$; choosing both internal rows contributes $(\sum_iT_iL_i)^2$. Choosing neither gives the squared five-spectator determinant. These three contributions with the same Green insertion give the complete one-loop term.

After including the two theta markings, the signed seven-leg bulk contact and the primary variations of the marked primitive, all 120 transverse rows on three ten-wave fixtures vanish. The soft limit was proved by the quartic. To test its 43-dimensional kernel, symmetrize products of four determinants and one covector on the nine independent scalar labels, requiring each label to occur. Each gives a covariant vanishing on a zero scalar. Canonical patterns record in which of the four determinants a scalar occurs, and which slot supplies the free covector. Two identical odd multiplicity patterns away from the free slot annihilate a symmetrization. Exact orbit sums on the first two fixtures give a 43-row minor equal to $175\pmod{1009}$. The quantum primary equation therefore holds for all momenta. The longitudinal and positive-descendant arguments above complete the formal genus-two lift of $F_3$ in the holomorphic jet model. [Signed primitive and graph replay](quintic_marks.py) · [Soft unisolvence certificate](quintic_unisolvence.py)

The next new equation, $F_4$, has $126-94=32$ soft transverse coefficients and $10-7=3$ new scalar primitive symbols.

### Genus-two lift for every constant scalar background

The primary correction is a single stable formula on the $n=r+6$ scalar slots of $F_r$. Start with the four overlap terms for $F_2$ and the three signed orbit terms for $F_3$ above. For a term displayed on $m$ differentiated labels, sum it over every $m$-element subset of the $n$ scalar slots. The other labels are undifferentiated. Add the following six terms; a term with $m>n$ contributes zero.

For each listed row, the integer at position $i$ specifies which of the four ordered determinants contains scalar $i$. Bits $1,2,4,8$ select the first, second, third, fourth determinant. In each determinant use the scalar rows in their displayed position order. Sum the product over injective assignments to the $n$ labels, requiring the assigned labels to increase among positions with identical integers. This matches the coefficient normalization of the replay.

| Differentiated slots $m$ | Incidence integers | Signed coefficient |
|---:|---|---:|
| 10 | $(3,3,3,3,3,12,12,12,12,12)$ | $5/1536$ |
| 10 | $(1,2,3,3,3,12,12,12,12,15)$ | $-5/15552$ |
| 10 | $(1,2,3,3,4,8,12,12,15,15)$ | $-149/37324800$ |
| 11 | $(1,2,3,3,3,3,12,12,12,12,12)$ | $1/1280$ |
| 11 | $(1,2,3,3,3,4,8,12,12,12,15)$ | $-133/2985984$ |
| 12 | $(1,2,3,3,3,3,4,8,12,12,12,12)$ | $1/23040$ |

All thirteen scalar terms carry the common multiplier $\lambda^4/(2\pi i)^{10}$. The zero-slot restriction argument reduces the successive required tests to the subspace of transverse covectors vanishing with one zero scalar. To construct columns in this subspace, use the four determinant bits and a fifth bit for the one free covector, and require every independent scalar label to occur. Polarization gives a symmetric polynomial. Evaluating it on the transverse basis at each vector slot gives the following exact certificates. The prime is 1009 throughout.

| $r$ | New target dimension | New scalar coefficients | Determinant modulo 1009 |
|---:|---:|---:|---:|
| 3 | 43 | 3 | 175 |
| 4 | 32 | 3 | 41 |
| 5 | 16 | 2 | 203 |
| 6 | 4 | 1 | 122 |
| 7 | 1 | 0 | 865 |

At $r=4$ the signed primary equations vanish on two 44-row fixtures. At $r=5$ they vanish on two 48-row fixtures. At $r=6$ they vanish on the 52-row fixture giving the four independent evaluations. At $r=7$ take two copies of each coordinate covector, an extra first and second coordinate covector, and $(2,-1,-2,1,0)$, followed by the negative sum. Use the first leg as the vector with polarization along the second coordinate. Before the twelve-slot correction its reduced value is $419/40$; the multiplier of that correction is $-241344$. The only new equation is
$$
\frac{419}{40}-\frac{241344}{23040}=0.
$$
The nonzero determinants certify all transverse primary equations up to the stable threshold. Beyond thirteen independent scalar waves, the zero-scalar restriction is injective on the entire transverse quotient. Thus all the genus-two transverse primary equations hold for every $r\ge1$. For each $r$ the remaining polynomial primary defect is a multiple of the vector momentum, which is canceled by a single scalar-descendant primitive. Positive-descendant exactness completes every higher arity. This proves a formal genus-two lift of each $F_r$ for the specified bulk counterterms in the flat translation-invariant holomorphic jet complex. [All primary graph equations and signed coefficients](stable_lifts.py) · [Four stable transverse certificates](stable_certificates.py) · [Polarized incidence sums](stable_orbits.cpp)

### Simultaneous sources and the weak-coupling effective action

The unweighted classical fivefold interaction has degree $-4$: a cubic trace uses ten total polyvector and Dolbeault degrees and three field shifts by $-2$. Thus for $|F_r|=2r-1$ the constant background parameter multiplying $F_r$ has degree $-(2r+3)$, in agreement with translation by $t^r\eta$. The common source weight $\lambda^{-2}$ makes the action degree zero. Use a finite set of independent odd background parameters, or work coefficientwise in a completed parameter algebra.

For a connected local heat graph of loop number $n$ with $k$ scalar insertions and $4n+2-k$ bulk vertices, the external polyvector rank is $5(4n+2-k)-4(5n+1)=6-5k$. Bulk counterterms can be chosen of rank five and marked counterterms of rank zero throughout the free resolution. Holomorphic jet decorations preserve the rank of a letter. Thus every term of parameter degree at least two vanishes in the local master equation. Together with the individual genus-two lifts this yields a simultaneous formal deformation by all the indicated scalar directions.

At positive scale a connected ordinary propagator graph can have a different number of edges. For $E$ edges and $V=E-n+1$ vertices its external polyvector and Dolbeault degrees are $E-5n+5-5k$ and $E-5n+5$, respectively. Already the tree chain of two scalar marked endpoints and four intermediate cubics has $E=5$, $V=6$, $k=2$, external rank zero and external Dolbeault degree ten. The incoming and outgoing ranks at its four cubics can be $(4,1),(3,2),(2,3),(1,4)$. For independent edge displacements, the five propagator four-forms in the twenty intermediate odd coordinates have a nonzero top coefficient. Taking their contraction covectors to be the five coordinate covectors omits edge $i$ in coordinate $i$; the four remaining differences in that coordinate form an invertible incidence matrix on the four intermediate vertices anchored at the endpoints. The five antiholomorphic four-forms are nonzero on five independent edge-coordinate differences and can be completed by the external forms. Such separated graph weights are included by BV renormalization-group flow. Local mixed source counterterms through genus two vanish by the preceding rank calculation. [Local integrated square and bulk rank count](../README.md#the-integrated-bilinear-and-simultaneous-scalar-deformations)

### Soft restriction on the genus-three trace and marked tower

At genus $g$, a primary marked covector of holomorphic derivative degree $10g+1$ uses $2g$ alternating volume tensors and one free slot. After eliminating the vector momentum, at most $2g+1$ singleton scalar gradients can survive symmetrization. There are at most $6g+1$ active scalar labels. Since $F_r$ has $r+2g+2$ primary scalar outputs, the soft-stability threshold is $r=4g-1$ for a strictly zero-slot-compatible source.

At genus three the nine-primary pure-trace candidate has 50 scalar coefficients. Put one scalar external field at constant momentum in the ten-leg bulk equation. It kills the bare geometric cofactor at the bivalent vertex. The one-loop correction involving the seven-primary genus-two trace vanishes when the constant scalar enters an ordinary cubic, or when it enters that trace, by the seven-slot identity used above. The surviving tree equation on the other nine waves is the transverse boundary of the eight-slot restriction of the candidate. When the cubic fuses the zero scalar with a nonzero leg, the latter receives only its own momentum and the scalar slot has rank zero. Such terms vanish on transverse fields. The restricted eight-slot trace has dimension 16. On actual four-vector and two-vector polarizations the matrix of its boundary on sixteen symmetric sextic minor products has determinant $243$ modulo 1009. Therefore the boundary is injective. Any nine-primary pure-trace solution must vanish on a zero scalar, leaving at most $50-16=34$ independent coefficients. The restriction map on invariant scalar symbols is surjective by extension of symmetric monomial sums. [Exact restricted-boundary certificate](three_loop_soft.py)

For the genus-three bulk solution constructed below, the known genus-one and genus-two local terms and the nine-slot pure trace respect the zero-scalar recursion of the marked problem. In the positive graph-loop contributions a zero scalar at a bivalent cubic kills the heat cofactor; at the genus-two primary bulk trace it kills its symbol. A zero scalar at the marked vertex restricts its polarized symbol to that of the previous $r$. The advective primary boundary commutes with this restriction. Thus a single formal lifting equation at $F_{11}$ suffices for this entire genus-three family. At its 19 scalar momenta the transverse covector space has dimension $31268$, and the scalar primitive space has dimension $637$, by the cycle-character calculation with partitions $(7,6,6,6,6)$ and $(6,6,6,6,6)$. The local jet calculation below gives a primary primitive at every arity, including $F_{11}$.

### First soft derivative of the genus-three bulk equation

The seven-primary genus-two trace has the shorter representation
$$
P=-\frac1{622080}\sum_{\substack{\{B,C\}\\|B|=|C|=5,\ |B\cap C|=3}}D_B^2D_C^2 .
$$
The two five-subsets cover all seven inputs. At seven slots the invariant trace space has dimension two and its one-zero-slot restriction has dimension one. Both the expression above and $P$ lie in the kernel. Evaluating on the five coordinate rows, $(1,2,3,4,5)$ and their balancing row fixes the coefficient. Hence a soft scalar at this vertex gives at least two powers of its momentum. The genuine three-loop heat cofactor and odd cofactor also each contain the bivalent scalar momentum, so the full primary bulk source vanishes to first order in a soft scalar.

For the candidate nine-slot trace, differentiate the tree equation at a zero scalar momentum $\epsilon q$. On the subspace already vanishing at zero, the terms fusing this scalar with a live slot give the same $\epsilon q$ times the unperturbed trace at every live slot. Their common translation drops out of the top exterior coefficient. Terms fusing two live legs have the usual eight-live-slot trace boundary with its ninth scalar slot differentiated in direction $q$. To impose conservation, perturb one live covector by $-\epsilon q$. Changing the balancing live slot changes a derivative tangent to the zero-slot restriction, which vanishes on the stated subspace.

The joint constant and first-derivative equations on the 50 symmetric trace coefficients have rank at least 32. An exact minor is $868\pmod{1009}$. The replay uses scalar evaluations for the constant restriction, whose physical boundary was proved injective above, and actual four-vector/two-vector contractions for the differentiated restriction. Conversely 18 explicitly specified $S_9$ orbit sums of six five-row determinants contain every scalar slot at least twice. They obey all the soft conditions and have independent evaluation determinant $840\pmod{1009}$. The permitted trace space is exactly their span. In particular every pure-trace solution admits a representative with two holomorphic derivatives on every primary letter.

At the first ten-wave quintic fixture, the transverse tree boundaries of these 18 traces have rank 18, with a minor $904\pmod{1009}$. Their coefficients are uniquely fixed by the full three-loop source at that fixture. The consistency calculation below proves that this determines the global solution. The list of six-minor words and an exact Chinese-remainder evaluator for the 18 boundary tensors are in the replay. [Soft-coordinate certificates and exact boundary arrays](three_loop_soft_jet.py)

On the planar family with triples $a_i e^1,b_j e^2$ and four retained momenta, every one of these 18 trace boundaries restricts to a multiple of the same transverse tensor. To see this, subtract the mean of the four retained covectors from each. Their span carries the three-dimensional standard representation of $S_4$, while the first two coordinate directions are fixed. Permuting the four retained labels is realized by these linear transformations; their determinant sign drops out of the sixth determinant power. Equivariance gives a common plane covector plus a common coefficient times the centered covector at every retained row. Subtracting a common shift and the longitudinal directions removes all four retained potentials. The collinear rows then lie in the plane, leaving only the $e^2$ component on an $a$-row and the $e^1$ component on a $b$-row.

Each remaining scalar component has two factors of every collinear strength. For a soft source distinct from the vector leg this is the common-translation calculation above. For its own soft row, expansion of the fused formula before projecting gives $-qP+\sum_j p_j(d_jP)(q)$ to first order, where $d_j$ differentiates the $j$th slot. The defining determinant weight gives $\sum_j p_j(d_jP)(q)=6qP$ by differentiating a linear change of coordinates. This first-order covector is longitudinal. The subtraction for the retained rows can be performed polynomially: the common standard-representation coefficient is read from a difference of two retained rows along their fixed difference, and the common plane part from their average. After division by $\prod_i a_i^2b_i^2$, the $e^2$ component on an $a$-row has separate degrees zero and one in $a,b$, respectively. Symmetry among the three $b$'s makes it a multiple of $\sum_j b_j$. Exchanging the first two axes fixes its partner. The explicit nonzero planar trace in the main calculation shows this restriction has rank one. Thus 17 of the 18 trace coefficients remain after imposing that family. [Polynomial example attaining the restriction](three_loop_trace.py)

For six arbitrary planar waves $u_1,\ldots,u_6$, the same $S_4$ argument removes the retained slots and the remaining potential at $i$ is detected by the binary determinant $[u_i,A_i]$. It has determinant weight seven on the two-plane and total planar degree fourteen. The quadratic soft property of the transverse boundary gives at least degree three at the distinguished wave $i$ and degree two at each other wave. The possible multidegrees are therefore $(4,2,2,2,2,2)$ with the four at $i$, or $(3,3,2,2,2,2)$ with one three at $i$. Let $V$ be the binary defining representation of $SL_2$. In the first case the other five undistinguished inputs carry
$$
\mathrm{Sym}^5(\mathrm{Sym}^2V)=\mathrm{Sym}^{10}V\oplus\mathrm{Sym}^6V\oplus\mathrm{Sym}^2V,
$$
which has no copy of the distinguished $\mathrm{Sym}^4V$. In the second case $\mathrm{Sym}^4(\mathrm{Sym}^2V)$ has summands of orders $8,4,0$, while $\mathrm{Sym}^3V\otimes\mathrm{Sym}^3V$ has orders $6,4,2,0$. The two common orders $4,0$ give precisely two possible invariant tensors. These binary decompositions follow by enumerating weights in a symmetric product of weights $2,0,-2$ and subtracting adjacent weight multiplicities. The traces $P_{223},P_{233}$ attain both independent components in the project calculation. Thus this family fixes two of the eighteen trace coefficients and leaves sixteen. [Full planar period and two polarized trace identities](../README.md#six-arbitrary-planar-waves)

In the zero-based order of the 18 orbit sums in the replay, write $c_j$ for their coefficients in $\Phi=C$. Using the first planar slot as a vector paired with the second coordinate, the two independent restrictions read
$$
\begin{aligned}
576(-4c_2-2c_7+4c_8+c_{17})&=K,\\
192(-12c_2-16c_7+12c_8+5c_{17})&=16K/9-4J .
\end{aligned}
$$
The first fixture is the two unit triples. The second uses $(1,0),(0,1),(1,1),(1,-1),(1,0),(0,1)$. Both retain the same three coordinate covectors and balancing final wave. Thus
$$
c_7=-\tfrac43(c_8-c_2)-K/10368+J/288,\qquad
c_{17}=-\tfrac{20}3(c_8-c_2)+K/648+J/144 .
$$

### The full quadratic-soft three-loop bulk target

Write $R$ for the polynomial ring of ten conserved five-component covectors and $I_i$ for the ideal generated by the components of slot $i$. Let $M=\prod_i I_i^2$. Before imposing conservation, the ten factors use independent variables. Each five-variable ideal of first or second powers of the maximal ideal has a free resolution of length at most four: its finite-length quotient has a filtration by residue fields, each resolved by the five-element Koszul complex. Tensoring the ten resolutions gives length at most 40. After conservation the ring has dimension 45. Away from the common origin one can eliminate a live leg whose ideal is locally the unit ideal, retaining independent ideals on the other nine slots. Exactness away from that point, followed by the descending free-resolution argument above, proves exactness also at the origin. The resulting zeroth module embeds in $R$ since its only possible kernel has finite length and its resolution length is below 45. Moreover $M=\bigcap_i I_i^2$: off the origin this is the product-intersection identity for independent slots; at the origin $R/M$ has resolution length at most 41 and has no finite-length submodule. These arguments also allow degree one in a specified factor.

The physical potential is a ten-tuple $A_i$ modulo $g+p_i h_i$, with $g$ a common covector and $h_i$ scalar. If all entries of this gauge array are in $M$, two generic nonparallel live rows force $g$ to vanish to second order on each soft locus. Thus $g\in M$ and $h_i\in I_i\prod_{j\ne i}I_j^2$. Common and longitudinal gauges intersect trivially over $R$. Quotienting double-soft representatives by exactly these polynomial gauges gives the full physical quadratic-soft target, as follows.

Let $C$ denote all covector arrays modulo gauge. At a prime containing exactly $t<10$ soft slot ideals, the local module $C/MC$ has projective dimension at most $4t+2$. Indeed it is the cokernel of the injection of the five common gauge copies of $R/M$ and the ten scalar gauge quotients just computed into fifty copies of $R/M$. Each appearing quotient has resolution length at most $4t+1$ after eliminating a live leg. With $t\ge3$ the prime has height at least $5t>4t+2$; the module has positive depth there. The same argument works at the common origin using $45>42$. For $t=1,2$ one need only check the remaining primes of height at most $4t+2$. There two live momenta are nonparallel: imposing collinearity on all live covectors has strictly larger codimension. The common gauge splits off using a unit minor of its map into these two live rows, after separately removing their longitudinal directions. The remaining module splits as free summands and copies of a one-slot quotient $R^5/Rp_i$ at the soft slots. In these local coordinates the soft slot variable blocks are independent, so imposing the powers on different blocks gives their product. With a single soft block these quotients are constant in all other variables and have no associated prime properly above the soft ideal. With two blocks the product quotient embeds in the sum of the two individual quotients, so has no associated prime on their intersection. Thus $C/MC$ has no submodule invisible at every generic single-soft locus. At such a locus the own row of a soft leg can have a first-order term in its momentum: a linear map taking every covector into its own line is scalar, exactly the longitudinal gauge. This proves the claimed saturation for transverse amplitudes.

For a permutation cycle of length $L$, the character of one independently assigned momentum slot of minimum degree $m$ is the complete-symmetric series starting in degree $m$, with every coordinate eigenvalue raised to the $L$th power. Multiply these factors around the cycles and multiply by $\prod_{a=1}^5(1-x_a)$ for the five regular total-momentum relations. The exterior sign alternation of the Weyl denominator extracts the multiplicity of $(7,6,6,6,6)$ for covectors and $(6,6,6,6,6)$ for the scalar longitudinal multipliers. The character average is over $S_{10}$ for common shifts and over the distinguished-slot stabilizer $S_9$ for a longitudinal coefficient. The three exact counts are
$$
\begin{array}{c|rrrr}
\text{minimum slot degree}&\text{array}&\text{common}&\text{longitudinal}&\text{physical}\\
0&4013&493&795&2725\\
1&2774&335&532&1907\\
2&770&96&156&518 .
\end{array}
$$

An explicit basis consists of 518 signed permutation sums. Each places slot 1's momentum at slot 0 and multiplies it by the six indicated ordered five-row determinants before summing over all $10!$ permutations. Both index conventions start at zero; five-subsets are in lexicographic order. For any one generic rational fixture the fiber has dimension $50-5-10=35$. Take the 15 integer fixtures in the replay. Subtract the last potential, project each of the nine remaining rows onto its four transverse directions, then eliminate the last longitudinal covector shared by these differences. The first 518 of the resulting 525 scalar tests on the listed basis have determinant $356\pmod{1009}$. Thus fifteen fixtures suffice and fourteen cannot test the entire soft target. The native replay groups permutations by their ordered pair of distinguished slots. The remaining $8!$ triple-minor products form a Gram matrix; reducing its entries modulo 1009 before the product bounds every unreduced Gram sum by $8!\,1008^2<2^{53}$, giving exact binary64 integer arithmetic. [Character calculation](soft_target.py) · [Basis words](soft_target_basis.json) · [Fifteen-fixture certificate](soft_target_basis.py)

The same conserved-momentum character calculation applies to the spaces with two derivatives on every slot at higher genus. At genus $g$ the trace has $2g+3$ slots and determinant weight $2g$; the potential has $2g+4$ slots and one additional covector. The dimensions are

| Genus | Trace | Potential array | Common shifts | Longitudinal shifts | Quotient |
|---:|---:|---:|---:|---:|---:|
| 2 | 1 | 6 | 1 | 2 | 3 |
| 3 | 18 | 770 | 96 | 156 | 518 |
| 4 | 2398 | 280142 | 24666 | 41599 | 213877 |
| 5 | 1005099 | 150267930 | 11041838 | 20458535 | 118767557 |

The integer cycle-character calculation has a signed 128-bit bound through genus five. An independent replay at 1009 uses Newton's recursion for symmetric powers of the single-slot series in place of cycle enumeration. These higher-genus dimensions count cochains; the ranks of the four- and five-loop consistency differentials are open. [Dimensions and independent recursion](soft_trace_dimensions.py)

### Injectivity of the transverse trace variation at every positive genus

At genus $g\ge1$ let $m=2g$. A primary quadratic-soft trace symbol $C$ is symmetric on $n=m+3$ conserved five-dimensional momenta and is homogeneous of degree $5m$, with determinant character $m$. It is a linear combination of products of $m$ five-row minors. Its transverse cubic variation is injective. In particular the image dimensions in genera four and five are exactly $2398$ and $1005099$. Since the next variation kills these images by the Jacobi identity, its respective ranks are at most $213877-2398=211479$ and $118767557-1005099=117762458$.

On $n+1$ boundary waves write $C_{ij}$ for the trace evaluated after fusing $i,j$. Its potential at slot $i$ is $A_i=\sum_{j\ne i}p_j C_{ij}$, with the sum over fusion partners. A zero physical boundary implies $A_i-A_j\in\operatorname{span}(p_i,p_j)$ for generic nonparallel $p_i,p_j$. Indeed place transverse trivectors at $i,j$, both containing an arbitrary common vector $w\in\ker p_i\cap\ker p_j$, and put scalars elsewhere. Choose the other four vectors to span the quotient by $w$. Only contraction of the repeated $w$ yields a top component, with coefficient $(A_i-A_j)(w)$ up to the nonzero oriented volume.

Insert $k$ soft waves $t_s q_s$ in $C$ and shift one hard wave to balance their sum. Define its leading joint symbol by extracting $\prod_s t_s^2$. Softness in each independent parameter makes this coefficient independent of the balancing hard slot. It is quadratic in each $q_s$ and symmetric in the soft slots. Polarizing in $q_sq_s$ therefore gives a symmetric $k$-linear function of symmetric two-tensors. Write $F_k(H;p)$ for its diagonal value at a symmetric matrix $H$, with $N=n-k$ conserved hard waves $p$. The zeroth symbol is $F_0=C$.

Suppose $F_{k+1}=0$. Apply the same extraction to a hard row difference of the boundary with $k+1$ soft waves, balancing at a third hard slot. All hard–hard fusions contain the vanishing joint symbol. In a fusion with soft slot $s$, its explicit prefactor $t_s q_s$ requires just the first derivative in $t_s$ of the remaining $k$-soft symbol. Subtracting hard rows $i,j$ gives the directional derivative transferring momentum $q_s$ from $j$ to $i$. The resulting sum over $s$ is the polarization of
$$
H(\partial_i-\partial_j)F_k(H;p)\ \in\ \operatorname{span}(p_i,p_j).
$$
Here $\partial_i-\partial_j$ denotes the derivative on the conserved hard configuration with increments of opposite sign at $i,j$. The inclusion holds for every $H$ by multilinearity, since rank-one symmetric tensors span. No null-boundary condition is imposed on the individual fused graphs.

Set $H$ to the identity in a complex coordinate frame. Represent the hard differential of $F_k$ by gradients at each slot, defined up to a common vector. Differences of gradients lie in the corresponding pair spans. On any generic independent triple, the sum of the three differences is zero. Independence identifies the coefficient along each $p_i$ in these differences with a single slot coefficient. Consequently the differential on the conserved configuration is a linear combination of the differentials of the $N$ squared lengths $\ell_i=p_i^Tp_i$.

To factor the polynomial through these lengths, use a complex Witt frame with components $(a_i,b_i,c_i,d_i,z_i)$, so $\ell_i=2a_i c_i+2b_i d_i+z_i^2$. Fix generic conserved $a,b,z$. The length equations have rank $N$ on the $2N-2$ independent $c,d$ entries. A linear relation on their differentials would require the row multiples of $(a_i,b_i)$ to equal one common vector; the nonparallel nonzero rows force every multiple to vanish. Hence one may use the independent $\ell_i$, the free $a,b,z$, and $N-2$ free dual entries as rational coordinates. The differential criterion in characteristic zero implies that $F_k(I;p)$ is rational in the lengths alone. Fixing generic constant $a,b$, zero $z$, and free dual entries zero gives a polynomial section of the length map. Substitution proves $F_k(I;p)=P(\ell_1,\ldots,\ell_N)$ with $P$ a polynomial homogeneous of degree $r=5m/2-k$.

Every invertible complex symmetric $H$ is congruent to the identity. Determinant covariance, with $m$ even, now gives exactly
$$
F_k(H;p)=(\det H)^{m/2}P\big(p_1^TH^{-1}p_1,\ldots,p_N^TH^{-1}p_N\big).
$$
At $H=\operatorname{diag}(t,1,1,1,1)$ the right hand side is $t^{-h}P(x_i^2+t b_i)$, where $h=r-m/2=2m-k$, the first coordinates $x_i$ sum to zero, and $b_i$ are squared lengths in the other four coordinates. The two Witt pairs there realize arbitrary $b_i$ by the same linear section, independently of the $x_i$. Polynomiality of the left hand side in $t$ forces all derivatives of $P$ of orders below $h$ to vanish at $(x_i^2)$ along $\sum_i x_i=0$. Thus $P(x_1^2,\ldots,x_N^2)$ is divisible by the $h$th power of this linear form. Independent sign flips of the $x_i$ give $2^{N-1}$ distinct linear factors, each with that multiplicity. For $k\le m-2$, one has $N\ge5$, $h\ge m+2$, and $2r=2h+m<2^{N-1}h$. The degree forces $P=0$, hence $F_k=0$.

For the initial step take $k=m-1$. The four conserved hard momenta span at most three dimensions. Under simultaneous soft scaling every five-row minor then vanishes through first order, including minors containing the balancing hard wave. A product of $m$ minors starts in degree at least $2m$, above the requested joint degree $2m-2$. Thus $F_{m-1}=0$. Descent to $k=0$ proves the claim. This gives uniqueness when a higher-genus transverse bulk solution exists. The equality of the consistency kernel with the trace image above genus three is a further condition.

The injection also holds over prime fields of characteristic greater than $5m(m+4)$, including $1009$ in genera four and five. Work first over the algebraic closure. In the Witt coordinates eliminate $c_N,d_N$ by conservation and solve an $N$-by-$N$ minor of the length equations. Cramer's rule expresses each original hard coordinate with common denominator of degree at most $N$ and numerator of degree at most $N+1$ in the rational coordinates. Therefore numerator and denominator of the transformed $F_k$ have degrees at most $5m(N+1)$. In any one free coordinate a rational function with vanishing derivative is a rational function of its characteristic-th power. This degree bound forces independence of that coordinate. All polarization orders are smaller than the characteristic. The rank-four degeneration can equivalently compare Taylor coefficients rather than normalized derivatives, and the signed hyperplanes are distinct in odd characteristic. Finally the components of $(A_i-A_j)\wedge p_i\wedge p_j$ have degree $5m+3$ in the independent hard momenta after eliminating the balancing wave. Vanishing of all these polynomials on the base prime field implies their identical vanishing, since this degree is below the field size. Thus finite-field momentum fixtures suffice to separate all trace symbols over that field.

The relation between arity and determinant weight is essential. With six conserved trace slots in dimension five, the square of their common signed maximal minor is a nonzero doubly soft trace of weight two with zero boundary. For seven boundary waves the kernel of their row matrix is spanned by the constant vector and a second vector with entries $z_i$. The oriented complementary five-row minors are a common scalar times $z_i-z_j$. The fused traces are their squares. Each potential row is the same common scalar squared times $\sum_j p_j z_j^2$, since both $\sum_j p_j$ and $\sum_j z_jp_j$ vanish. [Soft fusion, length coordinates and exceptional arity replay](trace_injectivity.py)

### A soft bivector in the bulk consistency equation

Let $A_i$ be the covector at hard slot $i$ of a symmetric doubly soft primary potential. On conserved hard momenta choose a balancing slot $k$. Write $D_i(q)$ for the derivative transferring the covector $q$ from $k$ to $i$, with $D_k(q)=0$. The constant and first-order soft bivector equations imply
$$
\sum_{\substack{i\\\mathrm{hard\ slots}}}p_i\wedge A_i=0,\qquad
q\wedge\sum_{\substack{i\\\mathrm{hard\ slots}}}p_i\wedge D_i(q)A_i=0 .
$$

Insert a constant bivector $H$ at zero momentum. Each hard–hard fusion puts the zero wave in the potential, so it vanishes. In the fusion with slot $i$, the new factor is the vector $\iota_{p_i}H$ wedged with the old field. Contracting an old field cancels after summing $i$ by conservation. Contracting the new vector gives $A_i(\iota_{p_i}H)$ times the common wedge of the old fields. Transverse vectors on five hard slots can have nonzero wedge, with scalars on the remaining slots. Varying $H$ proves the first equation by polynomiality.

At a fixed full-rank configuration let $P$ and $A$ be the matrices with rows $p_i$ and $A_i$. The first equation is precisely the symmetry of $P^{\mathsf T}A$. Choose any left inverse $L$ of $P$. The symmetric matrix
$$
W=AL+L^{\mathsf T}A^{\mathsf T}-L^{\mathsf T}(P^{\mathsf T}A)L
$$
satisfies $WP=A$. Diagonal entries of $W$ contribute longitudinal terms. A trace gives such a matrix with each off-diagonal entry evaluated on the fused pair. This last dependence imposes an additional condition.

Give the extra bivector momentum $tq$, require $\iota_qH=0$, and balance at $k$. Hard–hard fusions begin at $t^2$. In a soft–hard fusion the two new terms are $(\iota_{p_i}H)\wedge F_i$ and $tH\wedge\iota_q F_i$ for transverse hard fields $F_i$. The balancing shift leaves $\iota_{p_k}H$ unchanged. At first order, contracting the new vector gives the scalar $\sum_i (D_i(q)A_i)(\iota_{p_i}H)$ multiplying the wedge of the hard fields. Terms acting on an old field in slot $j$ cancel by determinant covariance. Indeed the infinitesimal frame change $p_i\mapsto p_i+t p_i(w)q$ for a vector $w$ gives
$$
\sum_{\substack{i\\\mathrm{hard\ slots}}}p_i(w) D_i(q)A_j
=A_j(w)q+ m q(w) A_j ,
$$
where $m$ is the determinant weight of $A_j$, equal to eight at genus four. After contracting $w$ against $H$, the second term vanishes and the first cancels the contraction of $A_j$ into $H$ in $tH\wedge\iota_q F_j$. The remaining common term is $H\wedge\iota_q$ of a hard four-vector. This is zero because $H$ wedges that four-vector to degree six and $\iota_q H=0$. Differentiating the old transverse polarizations gives zero by the constant identity. Testing all $H$ annihilated by $q$ proves the second equation. Both equations descend under common and longitudinal changes of the homogeneous equivariant representative. Changing $k$ in the second equation gives the derivative of the first equation wedged with $q$.

For a trace boundary $A_i=\sum_{j\ne i}p_j C_{ij}$, its fused scalar $C_{ij}$ has $D_i(q)C_{ij}=D_j(q)C_{ij}$. The pair $i,j$ therefore cancels in the second equation after wedging with $q$. This condition is independently measurable at genus four. For each six-slot subset $S$ write $d_i$ for the squared determinant of the five rows $S\setminus\{i\}$. On the complementary six slots take the sum of $d_a d_b$ over distinct unordered omissions, using their own six rows. For every unordered $i,j\in S$ add $d_i d_j$ times that complementary sum to $W_{ij}=W_{ji}$, with the diagonal zero. Sum all six-slot subsets. Every summand covers all twelve waves at least twice. The resulting $A=WP$ has degree 41, determinant weight eight and satisfies the constant identity everywhere. At the conserved fixture in the replay with $q=e^1$, the $e^1\wedge e^2\wedge e^3$ component of the first-derivative identity is $-8729370490059451776$, equal to $782$ modulo $1009$. Direct interpolation of the thirteen-field cubic variation on five transverse vectors of oriented volume $1332$ gives $-11627521492759189765632$ at first order. Thus a pointwise symmetric solution for the fusion coefficients does not ensure a local trace primitive. The example is a test potential constructed from minors. The physical graph source satisfies the consistency equation when the lower-genus master equations hold. [Exact fixed-fixture and cubic-variation replay](soft_bivector_descent.py)

### Exact global rank of the constant-bivector test

At determinant weight $m=2g$ let $F$ be any symmetric doubly soft two-form polynomial of momentum degree $5m+2$ on $m+4$ conserved waves. For $g\ge2$ it is the reciprocity tensor $\sum_{i\in\mathrm{slots}}p_i\wedge A_i$ of a doubly soft covariant potential. The common and longitudinal gauges both drop out of this tensor. Consequently the rank of the constant-bivector equation on the full potential space equals the entire two-form character multiplicity.

To construct a preimage, lift $F$ before momentum conservation using the product of the squared slot ideals above. One can take an equivariant homogeneous lift since every finite-degree space is a semisimple $GL_5$-module and permutation averaging is exact. Fix its multidegree component. Call a slot heavy if its degree is at least three; every other slot has degree exactly two. Apply the operator
$$
L=\sum_{i\in\mathrm{heavy\ slots}} p_i\wedge\operatorname{div}_{p_i},
\qquad
(\operatorname{div}_{p_i}F)_k=\frac{\partial F_{jk}}{\partial(p_i)_j}.
$$
Here divergence contracts the momentum derivative with the first index of the two-form. For a Schur component of highest weight $\lambda$ of the heavy coefficient variables, tensoring with the output $\Lambda^2 V^*$ subtracts one box in two distinct rows $i<j$; rows of terminal weight $-1$ are allowed. The operator is minus the mixed quadratic Casimir on these two factors: the matrix generator on the polynomial coefficients is multiplication by a momentum component followed by differentiation, and on the exterior covectors the dual action has opposite sign. With $c(\lambda)=\sum_{k\in\mathrm{weight\ rows}}\lambda_k(\lambda_k+6-2k)$ and $c(\Lambda^2V^*)=8$, its eigenvalue at weight $\lambda-e_i-e_j$ is
$$
\frac{c(\lambda)+8-c(\lambda-e_i-e_j)}2=\lambda_i+\lambda_j+9-i-j .
$$
Only $i=4,j=5$ and $\lambda_4=\lambda_5=0$ could give zero. For that summand to pair to determinant weight $m$, the light slots must have highest weight $(m+1,m+1,m-\lambda_3,m-\lambda_2,m-\lambda_1)$. The light slots are symmetric among themselves and quadratic individually. The symmetric powers of $\operatorname{Sym}^2 V$ have only even rows. Indeed, highest-weight polynomials on symmetric matrices restrict, on the dense triangular diagonalization chart, to products of diagonal pivots; a diagonal basis change multiplies each pivot by a square. Thus $m+1$ odd rules out the zero eigenvalue.

Apply $L^{-1}$ on the required multidegree component and take $A_i=\operatorname{div}_{p_i}L^{-1}F$ in heavy slots and zero in light slots. Differentiation leaves every slot of degree at least two. Permuting the multidegree components reassembles the covariant array, which then restricts to conserved waves. This proves surjectivity. Repeating the conserved character calculation with highest weight $(m+1,m+1,m,m,m)$ gives

| Genus | All potentials modulo gauge | Constant-bivector rank | Potentials satisfying reciprocity |
|---:|---:|---:|---:|
| 2 | 3 | 0 | 3 |
| 3 | 518 | 193 | 325 |
| 4 | 213877 | 77199 | 136678 |
| 5 | 118767557 | 37945455 | 80822102 |

An independent symmetric-power recursion reproduces every rank modulo $1009$. The explicit first-order control above violates an additional equation within the reciprocal space at genus four. [Exact character and modular replay](reciprocity_rank.py)

### Exactness of the four-loop fusion equations on a generic thirteen-wave configuration

Fix thirteen conserved covectors. Each first fusion joins two wave labels and gives twelve momentum groups. A second fusion gives either a triple or two disjoint pairs. Allow a scalar value for each of these $\binom{13}{3}+3\binom{13}{4}=2431$ partitions. Its trace variation at any first fusion is the twelve-wave array with entry at group $i$ equal to the sum of $p_j$ times the scalar assigned to the further fusion of $i,j$, summed over the eleven possible partners $j$. Each of the $78$ first fusions has $60-12-5=43$ potential coordinates after the common and longitudinal gauges. On this generic configuration the scalars span **every** collection of such potentials satisfying reciprocity and thirteen-wave consistency.

For the complete consistency test, distribute seven vector factors among the original slots and wedge the factors assigned to each slot. Each factor at slot $i$ belongs to $\ker p_i$. Their span of polarizations is the exterior algebra of the direct sum of these kernels. In the variation one divergence removes a vector from one of the fused slots and the potential removes a second vector. The remaining five vectors are paired with the volume. Thus the full tensor is the pullback of the volume under vector addition wedged with a two-form. Vector addition is surjective on generic polarizations. Vanishing of this tensor is equivalent to vanishing of the two-form on tuples $(v_i)$ with $p_i(v_i)=0$ and $\sum_{i\in\mathrm{slots}}v_i=0$. Each summand of the two-form is the wedge of the following linear forms:
$$
p_i(v_j)+p_j(v_i),\qquad
A_{ij}(v_i+v_j)+
\sum_{k\notin\{i,j\},\ k\in\mathrm{slots}} A_k(v_k).
$$
Here the second form uses the array on the corresponding twelve-wave fusion. Changing its common shift gives zero on the sum-zero space. Changing the longitudinal representative at the merged slot adds a multiple of the first form. This proves the gauge independence of the wedge. The seven-vector expression gives all allowed decomposable polyvectors: a nonzero-momentum slot has a four-dimensional transverse space. The replay independently checks the exterior sign against the signed Grassmann fusion on an eleven-wave configuration. Jacobi gives zero on trace variations, as can also be seen by canceling both orders of disjoint fusions and the three orders inside each triple.

The rank of this two-form test on reciprocal arrays has a closed expression. For $N$ generic conserved waves let $L_{ij}(v)=p_i(v_j)+p_j(v_i)$ on the sum-zero transverse vector space. Each $L_{ij}$ kills rotations $v_i=Jp_i$ for skew $J$. Conversely a motion annihilated by every $L_{ij}$, including the zero diagonal from transversality, is such a rotation. Indeed every linear dependence among the $p_i$ holds also among the $v_i$ after pairing with every $p_j$. The $p_j$ span, so $p_i\mapsto v_i$ is linear, and its symmetric bilinear form is zero. The space spanned by the $L_{ij}$ therefore has dimension $d=4N-15$.

At any first fusion a reciprocal array is a symmetric matrix times the momenta, modulo gauge. To see this, column-reduce the full-rank momentum matrix to an identity block with zero remaining rows, applying the dual change of basis to the array. Reciprocity makes the square leading block of the transformed array symmetric, and its other entries extend freely to a symmetric matrix. Transform back. Thus the entire image of the two-form map is spanned by $L_{ij}\wedge L_{kl}$ for disjoint pairs and $L_{ij}\wedge(L_{ik}+L_{jk})$ for triples of distinct labels.

An alternating bilinear form on the span of the $L_{ij}$ annihilating this image vanishes on disjoint pairs. On overlapping pairs its value on $(L_{ij},L_{ik})$ gives a totally alternating tensor in $i,j,k$ by the triple relation. For every relation $\sum_{k\in\mathrm{slots}} z_k p_k=0$, the edge weights $z_i z_j$ give a relation among the $L_{ij}$: their weighted sum reduces to a contraction with $\sum_k z_kp_k$ after omitting the zero diagonal. Testing this relation against $L_{ij}$ gives
$$
(z_i-z_j)\sum_{k\in\mathrm{slots}}z_k\,\omega_{ijk}=0 .
$$
At a generic $N$-tuple with $N\ge7$ the space of momentum relations contains a vector with $z_i\ne z_j$ for every fixed pair. Polynomiality in $z$ gives contraction zero for every relation. Thus $\omega_{ijk}$ is the value of one alternating three-tensor on $p_i,p_j,p_k$. Conversely every such tensor defines an annihilator. Indeed the relation with edge weights $w_{ij}=w_{ji}$ means $\sum_{k\ne i,\ k\in\mathrm{slots}}w_{ik}p_k=g+h_i p_i$, since the common and longitudinal potentials annihilate the sum-zero transverse space. The difference for slots $i,j$ belongs to their span. Evaluating the three-tensor on $p_i,p_j$ and this difference gives zero, so the bilinear form on edge labels descends to their span. The annihilator has precisely $\binom53=10$ dimensions. On a generic configuration where every first fusion has rank five and nonparallel momenta, the thirteen-wave rank $656$ is the instance $N=13$ of
$$
\operatorname{rank}=\binom{4N-15}{2}-10 .
$$

At one integral configuration the exact matrix ranks modulo $1009$ are

| Initial waves | Scalar values | Potential coordinates | Trace image | Reciprocity | Remaining two-form equations |
|---:|---:|---:|---:|---:|---:|
| 9 | 462 | 972 | 412 | 360 | 200 |
| 11 | 1155 | 1925 | 979 | 550 | 396 |
| 13 | 2431 | 3354 | 1918 | 780 | 656 |
| 15 | 4550 | 5355 | 3325 | 1050 | 980 |

The last column is the rank after restricting every first-fusion array to the kernel of reciprocity. Each modular rank gives a characteristic-zero lower bound at a generic configuration: the local gauge and reciprocity pivots remain invertible on an open set. The composition of the two maps vanishes identically. At thirteen waves $1918+780+656=3354$, so every lower bound is exact. The space of scalar assignments giving any one consistent collection has dimension $2431-1918=513$.

These assignments can be made compatible on overlapping configurations. On a specified first-fusion face there are $66$ scalar values, with $33$ independently fixed by the potential. After deleting those columns from the full thirteen-wave trace map, its rank is $1885=1918-33$. Hence the $513$-dimensional kernel maps onto **all** $33$ homogeneous choices on the common face. The same test works if prescribing the scalar values on two or three faces at once. For two faces the deleted-column rank is $1852$ and their union has $131$ scalar values. For three faces it is $1819$ and their union has $195$ scalar values except for the triangle of first fusions, which has $196$. Every simple graph type of up to three first fusions is included. In each case the restriction map onto the kernel of the prescribed face equations is surjective. At a generic configuration these equalities again follow from the modular lower bounds and the rank of the full map.

The same restriction-surjectivity holds for $k$ faces sharing one label, for every $1\le k\le9$. Their local rank is $33k$ and deletion leaves rank $1918-33k$. At $k=9$ their union prescribes $66k-\binom{k}{2}=558$ scalar entries with $261$ homogeneous choices. This also gives compatibility around cycles: start from a generic fourteen-wave configuration and select any ten distinct first fusions involving one fixed label. At the $j$th step, the new thirteen-wave tuple meets the preceding $j-1$ tuples on a star of $j-1$ twelve-wave faces. Indeed equality of sums at a generic tuple identifies the summands, so any shared eleven-wave partition must coarsen the join of the two first fusions. The existing scalar values satisfy the face equations and agree on intersections by induction. Surjectivity extends them to the new tuple. This constructs one assignment on the union of all ten thirteen-wave tests, including their triangular overlaps.

In particular one may attach a new thirteen-wave configuration by splitting a group of an existing twelve-wave face with a freely chosen new five-component momentum. Any old eleven-wave evaluation shared with the new tuple has to put the two new pieces back into the same group: otherwise an eleven-wave momentum depends on the free split. All shared values therefore occur on the parent face. The restriction-surjectivity just proved glues the scalar choices recursively along any finite tree of such generic splits. This gives a rational local solve in all the input momenta of each tree, by linear elimination on a pivot chart. A scalar depending **only** on its eleven final momenta throughout the full configuration space still requires compatibility beyond that tree. [Complete matrices, composition and gluing replay](fusion_fiber.py)

### A third-order test for compatibility across momentum splits

Choose a vector $v$, two slots $i,j$ with $p_i(v)=p_j(v)=0$, and fix all transverse components $p_s(v)$. Allow transfers within the four-dimensional covector hyperplane $H=\ker v$. If there is a scalar trace primitive, the gauge-independent measurement is
$$
(A_i-A_j)(v)=\sum_{\substack{s\\p_s(v)\ne0}}p_s(v)(C_{is}-C_{js}).
$$
Each fused scalar $C_{rs}$ is constant under every transfer $D_{rs}(q)$ of $q\in H$ from slot $s$ to $r$. Thus on the active slots this measurement must be a sum of polynomials, each independent of the split for its own edge. Here the momenta of all other slots and the total of the active slots are held fixed.

Before using frame covariance there is a complete ideal test for this decomposition. For $h$ outside slots, write $U=D_{ij}$ and $X_s=D_{is}$ for four-component blocks of commuting constant derivatives, evaluated along a basis of $H$. The other edge derivative is $D_{js}=X_s-U$. Write $(X_s)$ for the ideal generated by all four components. On each homogeneous polynomial degree the pairing by differentiation at zero is nondegenerate. The annihilator of the sum of the split-independent spaces is therefore the intersection of the $2h$ edge ideals. This intersection is the ideal of all polynomials $b$ such that
$$
 U_1^{\,h} b \ \in\ \prod_{\substack{s\\\mathrm{outside\ slots}}}(X_s)(X_s-U).
$$
Indeed expand $U_1^h$ as the product of $(X_s)_1-(X_s-U)_1$. Each term selects one generator for each $s$. The opposite ideals, one per $s$, are independent variable blocks and their intersection equals their product. Conversely $U_1$ belongs to none of the prime edge ideals. Testing all the resulting differential operators is sufficient for this one measurement to decompose into separate pair functions. Agreement of those functions with a single symmetric eleven-wave trace is a further constraint.

Frame covariance reduces the physical test. With only two outside waves their components are opposite. The reflection fixing $H$ pointwise and exchanging the two outside covectors fixes all inside rows at even determinant weight, so the measurement is identically zero. Take three outside slots $k,l,m$ with components $1,2,-3$. The shear $p_s\mapsto p_s+p_s(v)q$ with $q\in H$ has determinant one, fixes all spectators, and leaves the measurement and every trace term invariant. Hence $D_{ik}+2D_{il}-3D_{im}=0$ on these polynomials. Modulo this shear the six transfer blocks are
$$
X+tZ,\quad Y+tZ,\qquad t\in\{0,3,2\},\qquad
X=D_{ik},\quad Y=D_{jk},\quad Z=(D_{il}-D_{ik})/3 .
$$
The three blocks $X,Y,Z$ are independent on the reduced space. The saturated product test here has exponent four: its intersection ideal consists of $b$ such that
$$
Z_1^4 b\ \in\ \prod_{t\in\{0,3,2\}}(X+tZ)(Y+tZ).
$$
Interpolation expresses $Z_1^2$ as a constant linear combination of the three products of two distinct $(X+tZ)_1$, and likewise for $Y$. In their product expansion the two missing edge ideals, one $X+tZ$ and one $Y+t'Z$, are independent variable blocks. This proves the forward inclusion; $Z_1$ belongs to none of the six prime ideals.

The first independent equations on this physical slice occur at order three. The blocks $D_{ik},D_{il},D_{jk}$ span the three independent flavor directions. A degree-three operator in their intersection is trilinear. The remaining $D_{im}$ face enforces skewness between the first two arguments, and $D_{jl}=D_{jk}+D_{il}-D_{ik}$ enforces skewness under interchange with the third. The four independent alternating three-row minors kill all six faces, so they are exactly the degree-three equations. For example,
$$
\operatorname{Alt}_{a,b,c}\,
D_{ik}(a)D_{il}(b)D_{jk}(c)\,(A_i-A_j)(v)=0\qquad(a,b,c\in H).
$$
Here the alternation is the signed sum over the six orders of its three directions. At finite increments, applying one difference on each of the six edges also gives zero exactly. This tests 64 split configurations directly.

These equations impose additional conditions beyond the constant and first soft-bivector identities. For a six-slot set $S$ let $H(S)$ be the sum of products of two distinct squared five-row minors obtained by omitting one row of $S$. For each unordered partition into two six-slot sets $S,T$, add $H(S)H(T)$ to the symmetric weight of every pair separated by the partition. The resulting potential has determinant weight eight and two soft powers on every slot. Its reciprocity tensor vanishes by symmetry of the weights. In the first-bivector test, differentiating $H(S)$ in slot $r\in S$ in direction $q$, multiplying by $p_r$, and summing over $S$ gives $4qH(S)$ by determinant covariance, and likewise for $T$. Therefore the wedged first derivative vanishes. On the twelve-wave fixture in the replay use active inside slots $0,1$, outside slots $9,10,11$, and $v=e_5$. The alternating derivative in $e^1,e^2,e^3$ gives $-4726500469225893262080$, equal to $142$ modulo $1009$. The sixth finite difference with increment $e^1$ on every edge gives $-140342433697054124544$, equal to $776$ modulo $1009$. A single symmetric eleven-wave trace passes the finite test for both equal and unequal increments. These are obstruction tests on a model potential. [Exact differential and finite-split tests](fusion_transfer_box.py)

For the physical tree $\tfrac12\{I_2,I_2\}$ on that fixture, let $P$ be the seven-field genus-two symbol. Each partition $S,T$ with six labels on each side and $0\in S$, $1\in T$ contributes $-p_S(e_5)P(S,-p_S)P(T,p_S)$, where $p_S$ is the total momentum on $S$. The undifferentiated sum is $7698906551931026804317/403107840$. The sixth finite difference and alternating third derivative in the order used above give
$$
-\frac{219328932226387589}{2519424},
\qquad
\frac{1460732123756158891}{466560}.
$$
The remaining physical graphs must give the opposite values for a local eleven-wave primitive to exist. [Direct seven-field and cofactor evaluations](quartic_two_vertex_transfer.py)

A second test distinguishes an irrational period. Take the seven planar waves $e^1,e^1,e^1,e^2,e^2,e^2,-3e^1-3e^2$. In coordinates $e^3,e^4,e^5$ take five further waves $(1,2,0),(3,-1,0),(-2,1,1),(0,2,2),(-2,-4,-3)$. Measure the difference of the potentials at the first two of these five waves on $e_5$. Transfer $e^3$ along each of the six edges from the last three waves to the first two and take the product of the six forward differences. Every fused eleven-wave trace is constant under at least one of the six transfers. The tree with two seven-field vertices vanishes at every transfer configuration because each vertex admits at most three planar external waves.

In the one-loop hexagon with a nine-field vertex, exactly two planar waves $u_i,u_j$ and three normal waves $q_1,q_2,q_3$ must occur at the five ordinary vertices. More planar singletons kill their determinant, fewer leave at least six planar rows at the nine-field vertex. Eliminate one incident current by conservation. In each five-row minor of the vertex the three remaining normal rows then give the common factor $\det(r,s,\sum_k t_kq_k)$, where $r,s$ are the two normal external waves at the vertex and $t_k$ are the normal singleton positions in $[0,1]$. The six minors in the vertex give the sixth power of this factor. Their planar factor is the five-wave restriction of the integrated genus-three symbol. Its degree patterns $(4,2^4)$ and $(3,3,2^3)$ give at most three coefficients. Indeed $\mathrm{Sym}^4(\mathrm{Sym}^2)$ contains $\mathrm{Sym}^4$ once, and $\mathrm{Sym}^3(\mathrm{Sym}^2)$ contains each of $\mathrm{Sym}^6,\mathrm{Sym}^2$ once, the two summands of the symmetric square of binary cubics. Three evaluations on $A,B,D$ have nonzero determinant. The replay obtains the three coefficients from all eighteen nine-wave trace coefficients and verifies them on a fourth configuration with general mixed retained rows.

The heat and odd cofactors supply $[u_i,u_j]^2\det(q_1,q_2,q_3)^2$ times the grounded voltage. On the three-variable cube the Green kernel is $\min(t_i,t_j)-t_it_j$. Its monomial moments through homogeneous degree six are checked independently on the six ordering simplices. In the reduced conventions the full hexagon and its sixth finite difference are respectively
$$
-\frac{2994484730491}{6158592000}+
 \frac{4742213891}{11880000}\zeta(3),
\qquad
 \frac{52408071712547}{3302208000}-
 \frac{82996010347}{6370000}\zeta(3).
$$
[Integrated nine-field vertex insertion](quartic_seven_plane_hexagon.py)

The two-loop diagrams with one $I_2$ insertion have rational reduced periods on rational external data. On the simplex of edge times of total length one, the currents are bounded by the tree formula, the voltage by a finite sum of edge drops, and the second projective cofactor by deleting an edge of maximal time. The rational polynomial $P$ in the external and incident currents preserves absolute convergence against the first heat cofactor. A bridge has a zero cycle column, and the other columns have a positive time relation. Two blocks meeting at one vertex have independent time rescalings leaving all currents fixed, so their top heat Jacobian also vanishes. Suppressing bivalent vertices leaves three parallel paths. The degree sum above two is two on a two-loop core, allowing precisely two trivalent branch vertices here. For respective path lengths $x,y,z$, the only Kirchhoff denominator is a power of $xy+xz+yz$. Integrating segment monomials along each path gives rational simplex moments. Collect the resulting numerator. Scaling any pair of path lengths to zero on an open interval of positive ratios shows every monomial of the absolutely integrable polynomial density satisfies the strict convergence inequality there. In particular, for its integer exponents $a,b,c$ with $a+b+c=2d$, one has $0<a,b,c<d$. The substitution $x,y,z\mapsto1/x,1/y,1/z$ on the projective positive cone gives
$$
\int_{x+y+z=1}\frac{x^{a-1}y^{b-1}z^{c-1}}{(xy+xz+yz)^d}\,dx\,dy
=\frac{(d-a-1)!(d-b-1)!(d-c-1)!}{(d-1)!}.
$$
There are four ordinary planar jumps in the two-loop graphs and three planar waves at $I_2$. The four jumps fill the tensor product of the two cycle and two planar coordinates. Summing the placements of the four labeled jumps gives the squared Kirchhoff polynomial times twice the sum over unordered pairings of the squared planar brackets. The direct sum over all $3^4$ path assignments gives the same factor at unequal path lengths. At the quantum vertex the remaining four arguments give the fourth power of a three-row normal determinant. The vertex occupies either a trivalent endpoint with one external normal, or a bivalent site with two external normals. The reduced theta then has four bivalent positions, with path occupancies $(2,2,0)$ or $(2,1,1)$ up to permutation. Integrate their ordered position monomials and the conductance simplex using their factorial moments. The raw endpoint and bivalent normal integrals are $4213502549/170100000$ and $1864769/81000$. Their respective sixth differences are $825565581019/2546775000$ and $-175476439/220500$. The common signed planar factor is $17/160$. Thus the full two-loop insertion and its sixth difference are
$$
\frac{138201796633}{27216000000},
\qquad
-\frac{20420183920327}{407484000000}.
$$
[Both placements and unreduced five-dimensional checks](quartic_seven_plane_theta.py)

The self-tree vanishes. The genus-three hexagon and the genus-two theta exhaust the lower-vertex sources. Therefore any local solution of the genus-four master equation requires the **bare four-loop sixth difference** on the stated twelve waves to equal
$$
-\frac{502834915887479369}{31783752000000}
+\frac{82996010347}{6370000}\zeta(3).
$$
The integral tail bound $\zeta(3)<1+1/8+1/27+1/64+1/32$ gives the strict upper bound $-4408037849009413/63567504000000$. The required bare value is negative and irrational. Direct evaluation of the bare sum on this configuration remains open. [Projective cofactor convergence](../../nonSP/verify/k4_periods/README.md#schwinger-numerator-and-convergence)

A paired probe eliminates the entire theta sector before integration. Replace the seventh planar momentum by $t(e^1+e^2)$ and subtract $(3+t)(e^1+e^2)$ from the last normal wave. The changes to the tangential parts of normal waves cancel from every reduced normal cofactor. With the $-P$ vertex sign the planar theta factor is $t^2(18+4t+5t^2)/4320$. For the hexagon its two planar coefficients, in the order $1,\zeta(3)$, are
$$
-\frac{t^2(9133295639+2195940786t+2961866876t^2)}{23538138624000},
\qquad
\frac{t^2(2890121+695214t+938306t^2)}{9081072000}.
$$
The quadratic soft behavior gives the factor $t^2$. At most four powers can occur in the slot since all the other planar slots occur at least twice. The four nonzero integer tests thus fix each identity. Taking $27$ times the bare sixth difference at $t=-1$ minus $19$ times that at $t=1$, any local completion requires the exact value $9538548955331/37149840000-(3033930713/14332500)\zeta(3)$. [Independent planar placements and transfer evaluation](quartic_seven_plane_hexagon.py)

For an individual twelve-wave test keep the preceding five transverse projections and choose either of the [two integer seven-wave configurations](../README.md#a-four-loop-measurement-on-the-elementary-scalar). Give the five remaining waves arbitrary planar components satisfying total momentum conservation. Denote their three-dimensional projections by $q_0,\ldots,q_4$ and test the difference between slots $0,1$ on $e_5$ in this five-slot numbering. Fusing either with $s=2,3,4$ leaves seven planar rows and four retained rows. Each determinant of the fused trace acquires the determinant of the three untouched projected rows, by the change of frame fixing the planar rows. Polynomiality includes the singular cases. Thus the required total potential is the measured seven-wave trace times
$$
\sum_{\substack{s=2,3,4\\ \mathrm{outside\ slots}}}q_s(e_5)
\left(\det(q_{\widehat{0,s}})^8-\det(q_{\widehat{1,s}})^8\right)
=-4548608064 ,
$$
where a hat omits the two indicated rows out of five. Both probes annihilate the trace ambiguity. The two exact total potentials required by a completion are respectively
$$
-4075552825344K_4-419199719178240J_4,\qquad
-10261659792384K_4-330702000685056J_4.
$$
The two-vertex tree vanishes. Subtracting the integrated hexagon and theta insertions fixes the full required **bare** four-loop potentials on these individual configurations. In reduced units they lie strictly between $18550$ and $18551$, and between $36651$ and $36652$, respectively. These are necessary values; an independent unrestricted bare integral is still outstanding. The replay gives all six exact rational coefficients in the $1,\zeta(3),\zeta(5)$ expressions and checks the fused boundaries with the original eleven-slot determinant words for two different planar lifts. [Individual predictions and replay](quartic_seven_plane_target.py)

For a direct bare calculation on the paired family, remove all seven planar bivalent vertices while holding the other five positions fixed. There are fourteen remaining edge segments. Stack the constant projective row above the twelve three-component cycle-current rows. The fourteen alternating maximal minors of this $13\times14$ matrix give a kernel vector $w$. At a vertex $v$ and segment $f$, let $J_f(v)$ be the flow of a unit source at $v$ with a fixed sink, and $C_f$ the four-component cycle column. Set $K_v=\sum_f w_f C_f J_f(v)$. A planar insertion at fractional position $s$ on segment $e=(v,v')$ contributes the interpolated vector $(1-s)K_v+sK_{v'}$. Indeed its first segment current has the additional jump term $-(1-s)$ on that segment, whose cycle column is already one of the seven determinant columns. This term drops out by alternation.

For the three copies of $e^1$ assigned to segments $f,g,h$, write $\alpha(x)=\det(C_f,C_g,C_h,x)$. Let $a=\alpha(C_e)$ at the segment of the seventh insertion. Let $b$ be $\alpha$ evaluated on the sum of the three interpolated $K$ vectors and $t$ times the interpolated seventh vector, minus $(3+t)K$ at the tangential recoil vertex. The three copies of $e^2$ give an independent pair $a',b'$ with the same distribution. Laplace expansion gives the signed full determinant, up to orientation, as $t(ab'-ba')$. With the seventh placement fixed, the integral of its square therefore equals
$$
2t^2\left(\int a^2\,\int b^2-\left(\int ab\right)^2\right).
$$
Each triple integral sums $\binom{14}{3}=364$ distinct triples with six label orders and the product of segment lengths. The interpolated means and variances use only $\int_0^1 s\,ds=1/2$ and $\int_0^1(s-1/2)^2ds=1/12$. The final integral sums fourteen possible segments for the seventh wave, and the cofactor square includes the divisor $(\sum_f T_f)^2$. This replaces $14^7=105413504$ assignments and integrates their continuous positions exactly. The replay restores the original 21 edges at unequal times, compares independent seven-position moments on five allowed segments, and checks projective homogeneity after multiplying by the voltage and inverse fifth power of the Kirchhoff polynomial. On a $K_{3,3}$ subdivision the bare sixth-difference integrand is strictly positive at the recorded times. The integration over the remaining normal positions and core times and the topology sum are open. [Bare seven-wave integration and direct Kirchhoff checks](quartic_seven_plane_bare.py)

For arbitrary seven planar waves one can integrate the same determinant using exterior powers of the eight cycle-planar coordinates. After $k$ insertions let $W$ be the exterior product of their jumps, $V$ their total interpolated current including the fixed tangential lifts, and store $(W,W\wedge V)$. At the next position on edge $e$, the planar momentum $u$ gives jump $x=u\otimes C_e$ and residual current $z=u\otimes((1-s)K_v+sK_{v'})$. The exact update is
$$
(W,W\wedge V)\longmapsto
\big(x\wedge W,\,x\wedge(W\wedge V)+(-1)^{k+1}z\wedge x\wedge W\big).
$$
Its matrix is affine in $s$. If its midpoint and slope are $M,D$ and the previous second moment is $G$, integrating this edge gives $T_e(MGM^{\mathsf T}+DGD^{\mathsf T}/12)$. Sum over edges at each of the seven steps. The square of the final eight-form is the required numerator. The largest stored second-moment matrix has order $\binom84+\binom85=126$. On seven distinct planar directions the replay agrees with all $4^7$ individual assignments on four allowed segments. On the fourteen segments it agrees with the independent triple calculation. Arbitrary conserved lifts, normal momenta, normal test slots and polarizations are supported. [Exact exterior-covariance integration](quartic_seven_plane_bare.py)

The quadratic response to a tangential change of the five retained waves survives ordinary graph integration on individual labeled diagrams. In the seven planar slots take three momenta along $e^1$ and four along $e^2$, each set summing to zero. Lift the first two retained waves by opposite multiples of $e^1$. Four $e^2$ jumps give $4!\det M$ and the three $e^1$ jumps give $3!\operatorname{adj}M$ by Cauchy–Binet, with $M$ the cycle Gram matrix. The coefficient of the squared lift is nonnegative before the voltage. In the repository's lexicographic $K_{3,3}$ core put retained labels $2,0,1$ in that order on path $0$, label $3$ on path $1$ and label $4$ on path $3$. Use projections $(1,1,0),(2,-1,0),(-1,2,1),(1,0,-1),(-3,-2,0)$ and measure the voltage between $0,1$ in the last component. The unique positive source is $2$ and the unique sink is $3$. The branch reached from $2$ through $0,1$ lies strictly below the source potential by the resistor-network maximum principle. Linear interpolation gives a strictly positive voltage drop between $0$ and $1$ for all positive times. The replay checks a nonzero squared cofactor at unequal resistances, giving a strictly positive integral for this labeling. On a separate prism, permuting all five retained labels still gives a nonzero sum at fixed times.

The signed quadratic-lift determinant on the normal graph has five jump columns and nine initial-current columns. After using the constant projective row, at most eight columns depend linearly on the five normal position fractions. Moving the positions on one path changes the columns in a space of dimension at most four. The jump columns on that path already contain one of these directions, giving joint position degree at most three on the path. Four-point interpolation in each of the five fractions recovers the four components of the unsquared determinant. The squared polynomial then integrates by ordinary cube moments, with the same-edge Green minimum integrated separately. On the prism with core path resistances $(2,3,5,7,11,13,17,19,23)$, put the five retained waves on paths $(0,2,3,6,8)$. Use normal momenta $(1,2,0),(3,-1,0),(-2,1,1),(0,2,2),(-2,-4,-3)$ and add opposite tangential components at the first two slots. The sum over all $5!$ labelings after all five normal position integrations gives
$$
-\frac{948679641062255329222385597741966332326829569020906127211}
{541430684480898430785257421232153201165649075923292160000}
$$
for the quadratic cofactor contracted with $\operatorname{adj}M$ and the voltage, with the common positive planar and heat factors removed. A separate shared-path test integrates a nonzero Green minimum. Integration over the nine core resistances and the remaining marked graphs is still required. [Five-position integral and full permutation replay](quartic_seven_plane_lift_positions.py)

A doubled core edge with $m$ retained normals has a local cycle supported on $2+m$ segments. Its three normal current rows are orthogonal to the positive vector of segment times by Kirchhoff's voltage law. A nonzero minor requires $3\le 2+m-1$. Five normals therefore kill the core with three doubled edges. On the core with two doubled edges they permit either three normals on one bubble and two on the other, or two on each and one elsewhere.

For the three-plus-two case number the doubled paths $(0,1)$ and $(7,8)$. The four ordered path multisets $(0,0,0,7,7)$, $(0,0,0,7,8)$, $(0,0,1,7,7)$ and $(0,0,1,7,8)$ have stabilizers of orders $2,4,2,4$ in the core group of order $16$. At normal momenta $(1,2,0),(3,-1,0),(-2,1,1),(0,2,2),(-2,-4,-3)$, with opposite $e^1$ lifts at slots $0,1$ and the $e_5$ voltage difference, their integrated contributions including these symmetry factors are
$$
-\frac{17678261}{4400827200000},\qquad
-\frac{3279091}{1320248160000},\qquad
\frac{55782439}{264049632000000},\qquad
-\frac{430110173}{1188223344000000}.
$$
Their sum is $-3161360627/475289337600000$ before the common planar tensor. The numerator separates by the middle-path powers and the two arm degrees. For each arm the two bubble children have integer powers $a,b$. The projective beta factors require $a,b<13$ and $a+b>13$. A separated summand may violate these inequalities. Collect its coefficient as a polynomial in every segment time of the opposite arm. Every forbidden coefficient vanishes identically in the sum, after which each surviving monomial integrates by factorial ratios. Two different integral cycle bases give the same four answers, enumerating all $480$ label orders or their stabilizer orbits. For the first two shapes at the original label order direct integration of the $9611$ and $62601$ collected monomials gives $-3247/27725211360000$ and $-3247/9241737120000$. [Exact three-plus-two projection](quartic_seven_plane_double_three_two.py)

Adding the seven shapes with two normals on each bubble gives $367378814983/7129340064000000-3161360627/475289337600000=22854171827/509238576000000$ on this frame, before the planar factor. With three copies of $e^1$, three of $e^2$ and one $(2,3)$ the planar factor for opposite $e^1$ lifts is $1296$. The complete two-double-core quadratic response on the original fixture is therefore $22854171827/392931000000$. [Seven two-plus-two shapes](quartic_seven_plane_double.py) · [Four three-plus-two shapes](quartic_seven_plane_double_three_two.py)

For the complete quadratic tensor take two further normal frames
$$
(1,0,0),\ (0,1,0),\ (0,0,-1),\ (-1,0,1-c),\ (0,-1,c),\qquad c=2,3 .
$$
All eleven shapes give the stripped coefficients $2693233/354882705408000000$ and $51018137/638788869734400000$, respectively. The signed complementary triple minors give the five Gale coordinates $(-1,0,c-1,-1,0)$ up to a common shift, and $(4,8,14,-7,0)$ on the original frame. The cross product of the first two normal vectors equals the measured third normal basis vector on the two sparse frames and $-7$ times this vector on the original one. On each frame the seven-wave matching evaluated at the induced lift $Y=(z_0-z_1)e^1$ gives $5184(z_0-z_1)^2$. In the three normal generators of degree seven, ordered as in the [normal tensor calculation](planar_normal_tensors.py), the three evaluation rows have determinant $-64012032/5$. Dividing the cross-product measurements by the matching gives the complete coefficients
$$
\left(
-\frac{129013}{5110310957875200000},
-\frac{777251}{1022062191575040000},
-\frac{19109}{69686058516480000}
\right).
$$
Direct integration on the fourth frame $c=4$ gives $1860953/4655895552000000$, agreeing with the tensor determined by the other three frames. This replay reconstructs the projective cofactor from twenty integral source configurations. Label assignments with the same Gale coordinates have the same squared cofactor, so their voltage polynomials can be added before angular integration. Every seven-planar quadratic lift of the two-double topology follows from this tensor. The bare sum of the singly doubled, prism and bipartite cores must supply its negative if a genus-four trace solution exists. [Integrated tensor and angular replay](quartic_seven_plane_double_tensor.py)

All $120$ retained-label assignments of the quadratic lift can be evaluated from one current nullspace. Take four unit source vectors and their negative sum at the five positions. Delete the first edge column of the canonical four-cycle, four-source current array and transpose, giving a $13\times16$ matrix. At full row rank collect a basis of its kernel into entries $N_{r i \mu}$, where $r$ labels one of the four cycle blocks, $i$ one of the four unit sources and $\mu$ one of the three kernel vectors. For the first four physical normal momenta let $m_i=(-1)^i\det(q_{\widehat i})$, omitting row $i$ in zero-based order, and put $m_4=0$. The physical determinant selects all four source rows in its distinguished cycle block and three normal rows in each of the other blocks. The full block supplies the determinant of the four lift-plus-normal rows. Each other block supplies the four cofactors $m_i$. Complementary maximal minors of the canonical array are the three-row minors of $N$ up to one common factor. Thus the four physical determinants are the signed three-row minors of the $4\times3$ matrix with entries $\sum_{i=0}^3 m_i N_{r i \mu}$, multiplied by the common lifted determinant and the canonical factor. When the canonical array has smaller rank all four determinants vanish. Relabeling by $g$ replaces $m_i$ by $m_{g(i)}-m_{g(4)}$; the permutation signs cancel between the lifted determinant and the cubic.

The cubic has twenty monomials. Their span is attained separately by the prism and the bipartite core on paths $(0,2,3,6,8)$, with five sets of resistances giving square minors $990338,896731$ modulo $1000003$. In the first double-pair core put three retained events on the left bubble, two on the right and none on the middle branch. Set $x=m_0-m_2$, $y=m_1-m_2$, $c=m_2$, $d=m_3$. With cycle rows in repository order, component zero uses only the four cubics in $x,y$. Component one uses the six monomials obtained by multiplying one of $x^2,xy,y^2$ by $c$ or $d$; it vanishes if all three left events occupy the same path. The other two components vanish. Interpolation at twenty integral external-source frames checks these support statements as polynomial identities in all fourteen segment resistances. Squaring therefore leaves at most $34$ normal monomials, or seven if the left triple shares one path. In the voltage multiplying each monomial the ten pair-resistance weights sum to zero at every slot, by conservation and transversality of the two tested waves. The exact ranks of these weight arrays on the five normal momenta above are five and three respectively. [Nullspace, symbolic support and rank certificates](quartic_seven_plane_gale.py)

Three independent quadratic responses are invisible in the eight-planar sector even with full external symmetry, reciprocal potentials, determinant weight eight and double softness. Assign seven labels to three pairs and a singleton. Use each group in one squared five-row minor, selecting the other rows from five retained labels $0,\ldots,4$. Weight a symmetric pair of retained slots by the product of these four squared minors, then sum over all external permutations. In each model below every retained label occurs in at least two minors before squaring. Eight planar rows would therefore occupy at least nine entries, exceeding the four minors' total capacity of eight. For seven planar rows, a quadratic lift along $e^1$ requires all seven singly occurring labels in the plane. The singleton minor supplies a common squared four-row lift determinant and the planar labels supply a common matching factor. Divide by these nonzero factors. Take projected retained waves $(1,0,0),(0,1,0),(0,0,-3),(-4,-11/3,8),(3,8/3,-5)$. The columns measure the differences at pairs $01,02,03$ on their respective cross products.

| Omitted pairs of the three pair minors | Omitted label of the singleton minor | Weighted pair | $01$ | $02$ | $03$ |
|---|---:|---|---:|---:|---:|
| $01,01,23$ | 2 | $01$ | 6083520 | −3024384 | −133056 |
| $01,01,23$ | 2 | $02$ | −1750800 | 8443776 | −4291056 |
| $01,02,12$ | 3 | $01$ | −1704960 | −3091968 | 6386688 |

The determinant is $188775118544044032000$. Thus none of these three coefficient directions is determined by polarization of the eight-planar data under the listed kinematic conditions. [Exact minors, unsuppressed graph, and rank replay](quartic_seven_plane_lift.py)

### Twenty-five measurements for seven planar waves and five mixed waves

Keep seven momenta $u_1,\ldots,u_7$ in a plane. Write the remaining five as $(\ell_i,q_i)$ in the chosen planar and three-dimensional normal frames, with $i=0,\ldots,4$. Momentum conservation gives $\sum_i q_i=0$ and $\sum_i\ell_i=-\sum_j u_j$. In this section the pair measurement is $(A_i-A_j)(q_i\times q_j)$, where contraction with the normal vector $q_i\times q_j$ means $\det(q_i,q_j,\,\cdot\,)$. All six lower-vertex coefficient functions on this sector have been reconstructed from the integrated theta and hexagon. The remaining bare measurement belongs to a space of dimension at most 25 on the whole family, including arbitrary planar components on the five mixed waves.

In the oriented normal frame choose five numbers $z_i$ up to common shift by
$$
z_i-z_j=(-1)^{i+j-1}\det(q_{\widehat{i,j}})\qquad(i<j).
$$
The signed complementary minors form the exterior square of the kernel of the three-by-five normal momentum matrix. One kernel vector is the constant vector, so these equations hold together and $\sum_i z_i q_i=0$. Subtract the mean $\bar z$ when forming the planar covector $Y=\sum_i(z_i-\bar z)\ell_i$. At rank three this covector parametrizes the five planar components at fixed sum modulo all shears $\ell_i\mapsto\ell_i+Lq_i$ with $L$ a linear map from the normal three-space of covectors to the planar two-space. The measurement is unchanged by the shear, which fixes the testing normal vector in the dual action. Each bare heat determinant has eight planar rows, and seven are filled by the planar bivalent jumps. This leaves at most one factor of the planar lifts per determinant. The tested voltage depends only on the normal charges. The whole measurement is consequently of degree at most two in $Y$.

The coefficients are polynomial even where the normal rank drops. Before factoring $Y$, invariant polynomial dependence on the normal covectors is generated by their three-row determinants. They are differences of the $z_i$. At generic $z$, shear invariance factors the polynomial in lift variations through $Y$. In each lift degree $k=1,2$, specialize the zero-sum lift variations to a planar vector and its opposite on slots $r,s$. The factor is then $(z_r-z_s)^k$ times the corresponding homogeneous form on that vector. These factors have greatest common divisor one as the pair varies, excluding denominators in the coefficients. Polarization recovers all vector arguments. At $Y=0$ take all $\ell_i$ equal to minus one fifth of the total planar momentum.

In lift degrees $0,1,2$ there are respectively $3,1,1$ independent planar contractions. Every planar slot contributes degree at least two, and the total planar degree including the lifts is sixteen. At degree zero use the preceding $A,B,D$. At degree two all seven slots are quadratic. In $\operatorname{Sym}^7(\operatorname{Sym}^2\mathbb C^2)$ the binary quadratic representation occurs once. Its contraction with $Y^2$ is the matching polynomial of $(u_1,\ldots,u_7,Y)$, with the same normalization as in the planar graph code. At degree one exactly one of the seven slots is cubic. The symmetric sixth power of binary quadratics has highest weights $12,8,4,0$ under $SL_2$. Only weight four couples with a cubic to give weight one, with multiplicity one. The degree-one contraction is the bilinear polarization of the same matching polynomial in its last argument on $Y$ and $\sum_j u_j$. In particular **the linear lift term vanishes whenever the seven planar waves sum to zero**, also before the core integration after summing the planar labels.

After contracting with $q_i\times q_j$ the normal momentum degree is 27, or nine signed minors. In lift degree $k$ the remaining normal polynomial therefore has degree $9-k$ in $z$. Permuting the five labels permutes the centered $z_i$ with the additional permutation sign. The measurement is symmetric in the test pair since exchanging it reverses both the voltage and the testing vector. The residual polynomial thus alternates separately on the pair and on its complementary triple. Factor out their two Vandermonde products, with the sign of ordering the pair before the triple. The remaining translation-invariant polynomial has degree $5-k$, separately symmetric on the pair and the triple. Its generators have degrees $1,2,2,3$: the difference of the two means, the pair variance, and the second and third centered symmetric polynomials of the triple.

For each such polynomial extend the pair value to a symmetric matrix $M_{ij}$ by permutations, with zero diagonal. The identity $\sum_j M_{ij}=0$ follows from momentum conservation and the symmetric graph Green matrix, since the normal projection satisfies $\sum_j q_j\wedge A_j=0$. In this identity the seven planar momenta have zero normal projection. The relation $\sum_{i<j}z_i z_j M_{ij}=0$ follows by contracting twice with the other normal momentum relation. Expanding these identities coefficient by coefficient in four independent $z$ differences gives the following exact counts.

| Lift degree | Normal monomials before the identities | Rank of identities | Normal coefficients | Planar contractions | Number of readings |
|---:|---:|---:|---:|---:|---:|
| 0 | 9 | 3 | 6 | 3 | 18 |
| 1 | 7 | 3 | 4 | 1 | 4 |
| 2 | 4 | 1 | 3 | 1 | 3 |

The replay constructs an invertible rational 25-by-25 evaluation matrix using the two frames $z=(0,1,4,9,12),(3,-5,0,1,2)$. An independent fused scalar trace reconstructed from these evaluations passes ten further geometries and shear changes. For the physical equation the entire lower-vertex source is independent of $Y$. Its fused trace, hexagon and theta have linearly independent six-component normal coefficient vectors. Six original integrated readings determine these three vectors; two further normal geometries test all of them directly. Three planar configurations determine the rational and $\zeta(3)$ hexagon coefficients and the rational theta coefficient in $A,B,D$, with a fourth independent check. The function `required_bare` returns all three period channels and the coefficient of the undetermined eleven-field trace for arbitrary rational inputs with this normal polarization.

The allowed bare coefficient vectors thus lie on one specified affine line **if the local master equation holds**. All four linear and all three quadratic lift coefficients must cancel in the full integrated bare sum. The remaining eighteen constant-in-lift coefficients have one common unresolved trace coefficient, multiplying $22A-5B+18D$. This gives 24 independent compatibility equations in the 25-dimensional interpolation space. Direct topology integration is still needed to decide these equations. [Exact interpolation, integrated coefficient tables and checks](quartic_seven_plane_space.py)

### Seven-planar tests with coincident normal minors

Li–Wang (1.18) now supplies the analytic cancellation for all seven lift-dependent coefficients. At genus four and twelve labeled inputs its left side is zero on the transverse pair test. A fused eleven-field vertex has seven planar and four mixed arguments; the projections of the four sum to zero. At normal rank three any redistribution of their planar components with fixed sum is a shear. In the hexagon precisely two planar jumps fill the two planar cycle rows; the nine-field vertex retains five planar and four mixed arguments. In the theta precisely four planar jumps fill its four planar cycle rows; the seven-field vertex retains three planar and four mixed arguments. More planar jumps kill the heat determinant and fewer leave at most three arguments outside the plane at the quantum vertex. The two-genus-two tree cannot distribute seven planar waves between vertices which each admit at most three. The paired voltage in every case uses only normal charges. The fused and lower-vertex terms are therefore independent of the planar lifts. The genus-one primary vertex vanishes at rank at most four. The coefficient of each positive lift degree in the remaining bare sum is zero. On the original normalized probe this fixes the sum of the single-double, prism and complete bipartite cores at $-22854171827/392931000000$. The single-double term has now been integrated on all 36 shape orbits. This fixes the joint prism and bipartite reading. Their independent term-by-term integration remains open. [Single-double coefficients and the primitive target](../README.md#four-loop-bare-quadratic-variation-on-the-two-double-core) The constant-in-lift equation also has an analytic solution given by the eleven-point moduli integral; its unrestricted explicit period coefficients remain to be evaluated. [Trace, normalization and flux comparison](../README.md#all-genus-primary-vertex-from-the-worldsheet) · [Li–Wang, Propositions 1.8–1.10](https://arxiv.org/abs/2610.21927)

All 25 four-loop pair coefficients can be measured with at least two equal $z_i$. In the notation above use
$$
(0,1,2,0,1),\quad (0,1,3,0,1),\quad (0,1,4,0,1),\quad
(0,1,2,-1,0),\quad (0,1,3,-1,0).
$$
Each frame has integral conserved normal momenta of rank three. For the first three frames the $120$ normal labelings give $5!/(2!2!)=30$ distinct lists of $z_i$. The last two give $5!/2!=60$. For the quadratic term the signed current minors depend on those lists; the voltage weights of all represented labelings must still be summed.

At $Y=0$ use six normal readings with the three planar configurations $3e^1,3e^2,t$, where $t=e^1,e^1+e^2,2e^1+2e^2$. At $t=e^1+e^2$ add four readings at $Y=e^1$ and three at $Y=-e^1$. The three paired differences recover three linear and three quadratic readings. The last one-sided difference then supplies the fourth linear reading. In the four-generator basis below the normal determinants in degrees $7,8,9$ are
$$
-\frac{1728}{125},\qquad \frac{1410048}{3125},\qquad
\frac{216710774784}{1953125}.
$$
The resulting $25$-row matrix has full rank and yields the same change of basis as the original unrestricted measurements.

For a single coincident pair of $z_i$, alternation forces the pair measurement to vanish when the equal labels both belong to the test pair or both to its complement. In the remaining case permute the test pair to $(0,1)$ and translate and scale to $(0,1,x,y,0)$. Homogeneity and polynomiality recover every coincident configuration from this two-variable restriction. In normal degrees $4$ through $14$ its exact polynomial coefficient matrix has full column rank. In degree $15$ its rank is $23$ of $24$; the $(0,1)$ reading on $(0,1,2,3,5)$ restores rank $24$. Consequently coincident frames detect every normal coefficient on this family through genus six. At genus seven one constant-in-lift tensor escapes them. [Exact polynomials, interpolation nodes and reconstruction](quartic_coincident_frames.py)

On the prism and $K_{3,3}$ every subgraph of cycle rank $h$ needs at least $3h+1$ segments: the $3h$ normal current rows are supported there and annihilate its positive resistance vector. Enumerating the five marked positions modulo core automorphisms leaves respectively $80$ and $26$ path multisets up to symmetry, representing $786$ and $1134$ multisets before quotienting. Every one of these $106$ configurations has a nonzero full-label quadratic density at positive integral resistances on the first sparse frame. Grouping the $120$ orders by their $30$ distinct Gale lists retains each order's voltage, with independent literal determinant checks on both cores. These are fixed-resistance witnesses. [Orbit enumeration and exact modular witnesses](quartic_seven_plane_primitive.py)

### Four normal tensors at every genus

With $2n-1$ planar waves and five remaining waves, the same measurement uses normal polynomials of degree $m=2n+1-k$ in lift degree $k=0,1,2$. Center the five $z_i$ to sum to zero. The generating function for reciprocal symmetric pair measurements obeying the two normal momentum relations is
$$
\frac{x^4+x^5+x^6+x^7}{(1-x^2)(1-x^3)(1-x^4)(1-x^5)} .
$$
The denominator counts symmetric polynomials in the four nonconstant centered elementary symmetric functions. A proof using the row equations starts before imposing reciprocity. The alternating factors on the pair and complementary triple have total degree four. Their residual separately symmetric polynomials have generators of degrees $1,2,2,3$. The row-sum target on five slots alternates on the complementary four labels at each slot. Its numerator begins at the degree-six Vandermonde, with generator degrees $1,2,3,4$. The row-sum map is onto: row data $r_i$ lift off diagonal to $(r_i+r_j)/3-(\sum_s r_s)/12$.

In the row-sum kernel the remaining map is $M\mapsto z^{\mathsf T}Mz$, increasing degree by two, with alternating scalar target. Every alternating scalar is divisible by the five-slot Vandermonde $\Delta=\prod_{i<j}(z_i-z_j)$. Write $H_{ij}=\partial_i\partial_j\Delta$ using five independent arguments before centering. For $i\ne j$ set $M_{ij}=3H_{ij}+H_{ii}+H_{jj}$ and $M_{ii}=0$. The row sums vanish by translation invariance and $\sum_iH_{ii}=0$, the latter alternating of degree eight. The Vandermonde monomials have exponents $0,1,2,3,4$ in some order, giving $\sum_i z_i^2 H_{ii}=20\Delta$. Euler homogeneity and centering then give $z^{\mathsf T}Mz=(3\cdot10\cdot9-5\cdot20)\Delta=170\Delta$. Multiplication by symmetric polynomials proves surjectivity. Subtracting the two target series gives the displayed generating function.

For an explicit basis take the four symmetric weights $W_{is}=z_i^r z_s^h+z_i^h z_s^r$ with $(r,h)=(1,2),(1,3),(1,4),(2,4)$. Form their pair readings by summing $(W_{is}-W_{js})\det(q_i,q_j,q_s)$ over $s$. They have degrees $4,5,6,7$ and satisfy both normal momentum relations. At $z=(-2,-1,0,1,2)$ the determinant on the pairs $(0,1),(0,2),(0,3),(1,2)$ is $82944$. The four columns are therefore independent over the field of symmetric rational functions. The symmetric-polynomial module they generate has exactly the displayed series, proving it exhausts the kernel in every degree.

For $n\ge3$ the planar multiplicities are again $3,1,1$. The quadratic contraction uses the unique quadratic in the symmetric $(2n-1)$st power of binary quadratics. The degree-one contraction uses the single weight four in their symmetric $(2n-2)$nd power coupled with one cubic. At genus two the constant-in-lift multiplicity is two. Thus at genera $2,3,4,5,6,7,8$ the numbers of independent scalar readings sufficient on this family are $3,12,25,44,71,108,155$. These are counts for the stated polynomial tests on the bare source. The four-generator calculation uses no assumption that the bulk quantum master equation has a solution. [Polynomial basis, Hessian identities and finite-genus checks](planar_normal_tensors.py)

The simultaneous test can retain all thirteen choices of a first fusion with the anchor. The $513$ free scalar choices at each parent restrict to $33$ free scalar coordinates on each shared twelve-wave face. Take the difference of the two restrictions. At a triple overlap there is one shared scalar component corresponding to joining the anchor and all three labels into one group. Its three differences must sum to zero. These conditions are independent: on a shared face the eleven components joining a fourth label to its triple can be prescribed independently in the face kernel, by a rank-eleven test. Hence with $k$ parent tuples the compatible differences form a space of dimension $33\binom{k}{2}-\binom{k}{3}$.

| Parent tuples | Compatible differences | Certified rank of the difference map modulo $1009$ | Generic number of additional conditions, upper bound |
|---:|---:|---:|---:|
| 10 | 1365 | 1365 | 0 |
| 11 | 1650 | 1650 | 0 |
| 12 | 1958 | 1957 | 1 |
| 13 | 2288 | 2284 | 4 |

Indeed first choose scalar assignments separately on every parent using thirteen-wave exactness. Their differences belong to the indicated face kernels because the potentials agree on shared momenta. Adding a homogeneous solution at each parent changes these differences by the tested map. Every compatible difference can arise from a collection of local scalar assignments: for each shared scalar independently choose its value at a reference parent and add the prescribed differences. The triple equations give path independence among its parents. Thus the cokernel measures the remaining gluing conditions on this family. Eleven parents always glue generically. On all thirteen, four further tests suffice generically. The positive upper bounds are finite-field rank certificates; they give no nonzero evaluation on the physical graph sum. [Coupled overlap matrices and triple checks](fusion_overlap.py)

### Polynomiality of rational trace primitives at every genus

Suppose a symmetric rational trace of determinant weight on conserved five-dimensional momenta has a transverse boundary represented by a polynomial array. The trace is polynomial at every arity. The implication applies before imposing soft conditions. If the rational trace already has generic double zeros at the soft slots, its polynomial extension satisfies the same soft identities everywhere.

Work over the complex numbers and suppose an irreducible denominator has maximal pole order among all prime divisors of the trace. The connected group $SL_5$ fixes each of the finitely many pole divisors. Its irreducible equation is invariant since $SL_5$ has no characters. Group two slots together when their opposite momentum shifts by every covector preserve that equation. This is an equivalence relation since transfers compose, and the polynomial depends only on the sums in its equivalence classes. Translation invariance of the divisor is enough here, since a translate of its irreducible polynomial has the same highest homogeneous part. There must be at least six classes. With five or fewer, conservation leaves at most four independent covectors, and $SL_5$ acts transitively on their generic frames.

Split any wave into two and pull the divisor back to this boundary configuration. The linear fusion is a projection with affine fibers, so the lifted divisor remains prime with the same pole order. Its transfer classes are precisely the old classes with the two split slots in the class of their sum. For each pair $i,j$ take the leading residue of its fused trace at the chosen maximal order, denoted $R_{ij}=R_{ji}$. Here the residue means the coefficient after multiplying by the maximal power of the defining equation and restricting to its zero set. A nonzero $R_{ij}$ requires the pair to lie in the same class, because a function of their sum is invariant under their opposite transfers.

The leading potential row is $\sum_{j\ne i}p_j R_{ij}$. Polynomiality in the transverse quotient makes this a common covector plus a multiple of $p_i$. The common covector must be parallel to every class sum. Indeed, wedging the row with $p_i$ and summing over its class gives zero by the symmetry of $R_{ij}$. The class sums are pairwise nonparallel at the generic point of the divisor. Collinearity of any two of the independent class coordinates has codimension four, while the divisor has codimension one. Hence the common covector is zero. These residue equations pass to the quotient at the generic point, since the longitudinal and common gauge map has a unit maximal minor there, as seen from any two nonparallel live rows.

Keep the old momenta at a generic point of the pole divisor. Take the new slot $i$ at $tq$ and subtract $tq$ from the old slot that was split. For every old $j$ in its class, fusing $i,j$ leaves exactly the same old class sums. The defining denominator polynomial is therefore identical under these fusions, with the fused wave assigned to old slot $j$. Each $R_{ij}$ is the pullback of the same old trace residue. At $t=0$ the fused configurations all equal the original one, where that residue is finite and nonzero. Outside the class $R_{ij}$ vanishes identically. Wedge the longitudinal row of $i$ with $q$ and set $t=0$. This gives
$$
q\wedge\left(\sum_{\substack{j\\\mathrm{old\ slots\ in\ the\ class}}} p_j\right)
\left.\operatorname{Res}(C)\right|_{\mathrm{old\ configuration}}=0 .
$$
The class sum is nonzero at a generic point of the divisor. The soft direction $q$ is arbitrary, which contradicts a nonzero residue. Removing the balancing momentum identifies the coordinate ring with a polynomial ring, where a rational function without a divisorial pole is polynomial. Existence of a single rational primitive compatible with all momentum splits is a separate condition.

It also suffices to construct a symmetric determinant-covariant algebraic trace with compatible fused branches. Average its conjugates over the rational function field of the trace momenta. On the boundary take the normalized field trace from a common finite extension containing all the chosen fused values. Each pair fusion adds five independent transcendental coordinates to the function field of the trace inputs, so the minimal polynomial of the algebraic trace stays irreducible. Transitivity of the field trace therefore gives precisely the pullback of the same average at each pair. The transverse equation is linear and its target is rational, so averaging produces a rational solution. The pole argument makes this average polynomial.

A single-valued holomorphic primitive on a Zariski-open dense domain also extends. Among the finitely many exceptional divisors consider those across which the primitive fails to extend generically. These are preserved by $SL_5$ as above. Near a generic point of one lifted divisor, every cross-class fusion is holomorphic. Inside each class the identity obtained by wedging its rows and summing shows that the common gauge wedged with the class sum is holomorphic. Two nonparallel class sums recover the common gauge holomorphically by a unit minor.

At a generic old point of this divisor use $z$ for a transverse local coordinate given by its equation. Split a slot by $tq$ as above, and consider $q\wedge\sum_{j\in B}p_j C_{ij}$ with old partners in its class $B$. For $z\ne0$, the internal fusions preserve the class sums and remain near the old point, where there is only this exceptional divisor. For $z=0$ and small nonzero $t$, the partial sum extends holomorphically because the new row is longitudinal modulo holomorphic terms. Indeed a generic $q$ avoids every exceptional cross-class coincidence identically along the split line. Each such coincidence is algebraic and therefore has only finitely many nonzero values of $t$ on that line. Shrinking the bidisk leaves at most the point $(0,0)$ excluded. Hartogs extension across this point shows at $t=0$ that the original trace times $q\wedge\sum_{j\in B}p_j$ extends along its normal disk. This prefactor has a nonzero component. At generic smooth points of the divisor all negative Laurent coefficients in a normal coordinate consequently vanish. Analytic continuation in the tangential coordinates and Hartogs extension across the remaining codimension-two subsets give an entire function on the conserved momentum space. Its Taylor expansion and determinant homogeneity make it polynomial.

The fixture in the replay gives six nonparallel class sums spanning four dimensions. For six classes with two rows each, the off-diagonal residue and gauge equations have zero kernel. The same holds with five rows per class. At six rows per class, each block has one linear relation on its rows. The square of its vector of coefficients gives a symmetric matrix in the block kernel after a diagonal gauge, and the replay finds the six independent kernels. On a five-wave class split to six, the soft test detects a nonzero common fused residue directly. [Fixed-divisor and soft-split residue calculations](trace_pole.py)

### Genus-three bulk consistency and the trace counterterm

A lower bound on the transverse trace variation at arbitrary arity follows by imposing three powers of momentum at each soft slot. Suppose the symmetric trace symbol $C$ vanishes to this order on the conserved locus, and its transverse boundary vanishes. For an $(n+1)$-tuple of wave momenta the potential at slot $i$ is $\Phi_i=\sum_{j\ne i,\,j\ {\rm slot}}p_j C_{ij}$, with $C_{ij}$ the symbol after fusing that pair. At independent slots $i,j$, vanishing implies $\Phi_i-\Phi_j$ lies in their two-dimensional momentum span. Indeed use an adapted frame with $p_i=e^1,p_j=e^2$ and polarize only these two slots by the four-vector $e_2\wedge e_3\wedge e_4\wedge e_5$ and the bivector $e_1\wedge e_k$ for $k=3,4,5$. The three pairings measure the difference on the common kernel of the two covectors.

Starting from $n$ conserved momenta, append $tq$ and subtract $tq$ from the last old momentum to restore conservation. At every old slot $i$, fusion with another old slot contributes only at order $t^3$ or higher by softness. Fusion with the new slot contributes $tq C_{i,n+1}$. Thus $\Phi_i-\Phi_j$ starts at $t^2$ with coefficient $q$ times the derivative of $C$ transferring momentum $q$ from slot $j$ to slot $i$. Wedge the difference with the two actual old momenta, including the balancing shift if needed. Its quadratic coefficient is the transfer derivative times $q\wedge p_i\wedge p_j$. This forces the derivative to vanish for generic $q,p_i,p_j$, hence identically. All transfers generate the tangent space to conservation. The constant that remains is zero by softness. Consequently the variation is injective on the cubic-soft subspace; any doubly-soft trace in its kernel is determined by the quadratic Taylor coefficient in a single soft slot.

For the cubic-soft character remove also the fifteen quadratic monomials from the five-variable single-slot series. The conserved product/intersection and free-resolution dimension arguments above hold for a third power of the maximal ideal by the same filtration of its finite-length quotient. At genera $2,3,4,5$ the cubic-soft trace dimensions are respectively $0,2,122,31671$. Newton recursion followed by Weyl coefficient extraction works modulo $10000019$, greater than the already established exact dimensions of the containing doubly-soft spaces. In particular the transverse trace maps have ranks at least $122$ and $31671$ at genera four and five. [Cubic-soft character and rank bound](cubic_soft_trace.py)

More generally allow at most $h$ of the $n=2g+3$ slots to have degree exactly two before imposing conservation, with degree at least three on the others. These spaces form an increasing filtration of ideals. At the step with precisely $h$ quadratic slots the associated graded has one summand for each choice of those slots. Each chosen factor is the quadratic monomial space annihilated by its slot variables, hence has a length-five free resolution. Each unchosen cubic ideal has a resolution of length four. Extensions give the entire filtered ideal resolution length at most $4n+h$. Away from the common momentum origin eliminate a live slot by conservation. Localizing its variables turns both the quadratic and cubic ideals of that slot into the unit ideal, leaving the corresponding filtration in the independent remaining slots. As in the preceding resolution argument, $4n+h<5(n-1)$ gives exactness after imposing conservation and injection into the conserved polynomial ring. Symmetrization is exact in characteristic zero. Thus the contribution from precisely $j$ quadratic slots, $0\le j\le h$, is computed before the conservation factor by taking the symmetric power of order $n-j$ of the cubic-soft single-slot series and multiplying it by the symmetric power of order $j$ of its quadratic monomial space. The cumulative dimensions are

| Genus | $h=0$ | $h=1$ | $h=2$ | $h=3$ | Entire doubly-soft space |
|---:|---:|---:|---:|---:|---:|
| 3 | 2 | 6 | 12 | 17 | 18 |
| 4 | 122 | 649 | 1527 | 2149 | 2398 |
| 5 | 31671 | 198576 | 525868 | 827474 | 1005099 |

For transverse primary polyvectors $F_i$ of momenta $p_i$, with $\iota_{p_i}F_i=0$, the ten-wave cochain of a potential $A$ is
$$
T_A(F_1,\ldots,F_{10})=
\sum_{\substack{1\le i\le10\\\rm slots}}
(-1)^{|F_1|+\cdots+|F_{i-1}|}
\big(F_1\wedge\cdots\wedge\iota_{A_i}F_i\wedge\cdots\wedge F_{10}\big)_{\mathrm{top}} .
$$
Here $|F_i|$ is its polyvector rank and the subscript pairs the resulting five-vector with the volume form. A longitudinal shift kills its own transverse input. A common shift contracts the wedge of all ten inputs, a six-vector in dimension five.

The primary cubic variation $d_{\mathrm{CE}}$ sums over all pairs of slots. Move each pair $i<j$ to the end with its Koszul sign and replace the pair by $\iota_{p_i+p_j}(F_i\wedge F_j)$ of momentum $p_i+p_j$. The odd contraction gives one additional sign $(-1)^{\sum_{h\ne i,j}|F_h|}$ from the preceding slots. The fused field is transverse since contraction squares to zero. This is the Chevalley–Eilenberg differential in the symmetric, suspended convention. For two successive fusions on disjoint pairs the two odd contractions anticommute. For three input waves the remaining three signed terms give the Jacobi identity of the derived Schouten bracket. Explicitly, let the divergence on a wave of momentum $p$ act as $\iota_p$. It squares to zero, and its deviation from the graded derivation rule on a wedge has exactly the two terms contracting either field by the other's momentum. There is no third-order deviation on a triple wedge. Expanding the square of the divergence on that wedge therefore gives the three Jacobi terms on transverse fields. Hence $d_{\mathrm{CE}}^2=0$. The 18 injective trace boundaries belong to its kernel.

For an integer test take the eleven momenta in the replay and sum all 55 fusions. Put $\iota_{p_i}$ of the coordinate volume five-vector in one slot, three independent transverse vectors in three other slots, and unit scalars elsewhere. The eleven-slot variation takes ranks $4+1+1+1=7$ to ten-slot rank six. The three consecutive polarization seeds 5432, 5433 and 5434 give 598 nonempty rows on the 518 explicit soft potential classes above. The first 500 rows on the first 500 columns give
$$
\det =27\pmod {1009}.
$$
This proves rank at least 500 in characteristic zero. Jacobi and the boundary minor $904\pmod{1009}$ give rank at most $518-18=500$. A second integral eleven-momentum fixture, polarization seed 19875, also has rank 500 modulo 1009. Thus after any characteristic-zero scalar extension,
$$
\ker d_{\mathrm{CE}}\big|_{\text{ten-primary soft potentials}}
=\operatorname{span}\{d_{\mathrm{CE}}P:P\text{ in the listed eighteen nine-primary traces}\}.
$$

Insert the chosen bulk action through genus two in the effective BV master equation. Its first possible genus-three defect has arity ten, descendant level zero and rank six. The three-loop bare term and the hexagon with the genus-two trace each have the potential form. At a bivalent leg the heat and odd cofactors each contain its external momentum. At a leg incident to the genus-two trace its determinant symbol is itself quadratic-soft. Thus the source defines a class in the 518-dimensional space above. Applying the square-zero effective BV operator again gives, at genus three and arity eleven, its ten-slot classical cubic variation plus the free variation of the arity-eleven defect. All lower-genus defects vanish. On transverse holomorphic waves the free variation is zero: the Dolbeault derivative vanishes and $t\partial$ contracts an input with its own momentum. Thus the ten-slot class belongs to the displayed kernel. In the convention $\Phi=C$, the unique doubly soft pure trace with this boundary enters the action with a minus sign. Once it is included, the arity-ten remainder has zero transverse symbol and has a local free primitive by primary torsion-freeness above. At arities $N>10$ the residual has positive descendant level $N-10$, so the same free resolution and BV arity induction complete the formal genus-three bulk action in the translation-invariant fivefold holomorphic jet complex. The periods of its eighteen coefficients are determined by the full genus-three diagram sum at the ten-wave fixture of minor 904. [Eleven-wave consistency replay](soft_consistency.py) · [Nine-primary trace boundary](three_loop_soft_jet.py)

### Global symmetric polynomial weights for reciprocal potentials

Write $E$ for the vector space of momentum labels with the single conservation relation. In a basis of independent momenta take potential row differences relative to the eliminated leg. At determinant weight $m=2g$ let $S$ send symmetric scalar weights to the potential $(SW)_i=\sum_jW_{ij}p_j$. Reciprocity is $\sum_i p_i\wedge A_i=0$. Over characteristic zero every reciprocal polynomial potential admits polynomial symmetric weights on the full conserved locus.

Indeed the three flavor representations for sewing followed by reciprocity are
$$
\operatorname{Sym}^2E\otimes S_{(m^5)}E
\ \longrightarrow\ E\otimes S_{(m+1,m^4)}E
\ \longrightarrow\ S_{(m+1,m+1,m,m,m)}E ,
$$
where $S_\lambda E$ is the Schur representation of highest weight $\lambda$, and an exponent repeats a row length. The middle Pieri summands besides the rightmost shape are $(m+2,m^4)$ and $(m+1,m^4,1)$. Put $D=\det(p_1,\ldots,p_5)$. For the first, the weight $e_1^2D^m$ sews to $e_1\otimes p_1D^m$. For the second let $\Delta_{\widehat a}$ be the five-row determinant omitting $a$ from the first six labels. Take $e_1\odot e_a=e_1\otimes e_a+e_a\otimes e_1$. The alternating weight
$$
\sum_{a=1}^6(-1)^{a-1}(e_1\odot e_a)\Delta_{\widehat a}D^{m-1}
$$
is the highest vector obtained by alternating six labels and multiplying by the other highest factors $e_1$ and $D^{m-1}$. Upon sewing, the term with output $e_1$ summed against $\sum_a(-1)^{a-1}p_a\Delta_{\widehat a}$ vanishes. The other term has sixth row $-p_1D^m$. Reciprocity hits the rightmost component on $e_1\otimes p_2D^m$. This proves exactness in the middle. For five independent labels the six-row summand is absent; for $m=0$ this is ordinary symmetrization and antisymmetrization of $E\otimes E$.

The horizontal two-box rule gives just one further summand on the left, of shape $(m^5,2)$ for $m\ge2$ and $\dim E\ge6$. It is the whole kernel of sewing. A highest vector has weights $c_ic_j D^{m-2}$ with $c_i=(-1)^{i-1}\Delta_{\widehat i}$ on the six labels and zero elsewhere. Each row sews to zero by $\sum_jc_jp_j=0$; the $(6,6)$ entry is $D^m$. Thus all symmetric weight changes fixing an abstract potential are generated by flavor polarization of this squared six-row determinant relation.

For an explicit inverse in any independent label basis use $(B A)_{ij}=(\operatorname{div}_{p_j}A_i+\operatorname{div}_{p_i}A_j)/2$, where divergence contracts a derivative in that momentum with the spatial covector index. Differentiating a polynomial splits off one label tensor and the displayed average symmetrizes the two labels, so this operation is $GL(E)$-equivariant. On the two sewn summands $SB$ has eigenvalues $\alpha=m+5$ and $\beta=(m+4)/2$, respectively. For the first, the only nonzero divergence is $(m+5)D^m$ in entry $(1,1)$. For the second evaluate its sixth row at the coordinate frame for the first five momenta, with $p_6=0$. Differentiating the sixth row in slot one contributes $-(m+5)p_1$; differentiating the first row in slot six contributes $p_1$ by the alternating determinant formula. Therefore one abstract symmetric weight tensor is
$$
\widehat W=\frac{B\big((\alpha+\beta)A-SBA\big)}{\alpha\beta},
\qquad S\widehat W=A .
$$

If a representative of $A$ vanishes twice in every physical momentum, its abstract tensor has total degree at least two in the direction of each label $e_k$: choose any independent flavor basis containing this label and scale it alone. Total degree counts both polynomial momentum occurrences and explicit flavor indices. Equivariance gives the same weight bound for $\widehat W$, simultaneously for all $k$. Return to the redundant physical labels. Any two lifts of $\widehat W$ differ by $z_i+z_j$, whose sewing is a common potential. The unique lift of zero diagonal is
$$
W^\circ_{ij}=-\tfrac12\,\widehat W(e^i-e^j,e^i-e^j)
$$
with $e^i-e^j$ the dual difference, which descends to $E$. For $k\ne i,j$ it annihilates the flavor vector $e_k$, so both required occurrences are in the momentum polynomial. This gives simultaneous double zeros on every spectator slot by a polynomial formula.

If there are at least eight physical slots the endpoints can also be required to vanish once in every off-diagonal entry. At $p_i=0$ the weight bound forces $\widehat W$ to be a multiple of $e_i^2$. Hence $H_i=W^\circ_{ij}|_{p_i=0}$ is independent of $j\ne i$. For any remaining slot $k$ choose $j\ne i,k$; spectator softness gives $H_i\in I_k^2$ on the conserved remaining slots. For at least seven remaining slots the intersection of these ideals is their product by the free-resolution argument above. Thus $H_i$ lifts to a polynomial $Z_i\in\prod_{k\ne i} I_k^2$ in the original conserved ring, by lifting the coefficients of the product ideal. Then $W_{ij}=W^\circ_{ij}-Z_i-Z_j$ has the claimed endpoint and spectator orders. The diagonal $-2Z_i$ is twice soft in every spectator. Averaging the lifts over permutations and projecting to determinant weight in their finite polynomial degree preserves the conditions. These symmetric weights impose reciprocity; requiring them to arise from one trace symbol evaluated at every fused pair is the remaining bulk equation. [Highest vectors, differential inverse and zero-diagonal lift](reciprocal_sewing.py)

### Rank-four sources and the genus-four bulk consistency test

The twelve-wave bulk potential is completely detected by one rank-four wave, two vectors and scalars. With the four-vector $\iota_{p_k}\Omega^{-1}$ in slot $k$ and transverse vectors $v,w$ in slots $i,j$, respectively, take the slots in the order $k,i,j$. The contraction formula gives
$$
T_A=(A_i-A_k)(v)\,p_k(w)-(A_j-A_k)(w)\,p_k(v).
$$
For generic momenta take $v\in\ker p_i\cap\ker p_k$ and choose $w\in\ker p_j$ with $p_k(w)=1$. This measures the full pair difference on the common kernel. Its vanishing for all pairs means $A_i-A_k\in\operatorname{span}(p_i,p_k)$. On triples of independent momenta comparison of coefficients gives $A_i=a+c_ip_i$ with common covector $a$. Such an array represents zero. A rank-four source with independent scalar potential $h$ has polarization $\iota_{dh}\Omega^{-1}$ and can perform these tests.

The primary Hamiltonian insertion of the integrated $F_r$ instead has $h=u^{r+1}/(r+1)!$ at scalar field $u$. Its momentum is the sum of a group of $r+1$ scalar momenta and polarization sums the formula above over all those groups. In particular $r=0$ acts on each transverse input by wedge multiplication by $\iota_{p_i}\Omega^{-1}$, which kills positive-rank inputs since $\ker p_i$ has dimension four. The terms where the potential contraction misses this new factor sum to zero by conservation. The remaining term contracts $\sum_i p_i\wedge A_i$ into the volume five-vector. A symmetric Green matrix $A_i=\sum_jG_{ij}p_j$ makes this zero, including after adding common or longitudinal gauge. For $r>0$ the nonlinear substitution still defines a classical differential on scalar and divergence-free vector fields together with the rank-four insertion. Indeed $\iota_{du}\iota_{d h}\Omega^{-1}=0$ and a volume-preserving vector acts on $h$ by transport. A rank-six bulk cocycle pulled back to first order in this insertion therefore yields a rank-two cocycle in the integrated scalar and low-polyvector complex. Equivalently the classical differential commutes, with its graded sign, with the Hamiltonian variation of the closed functional $F_r$. [Twelve-wave direct exterior test and nonlinear-source control](rank_two_variation.py)

Every such rank-two cocycle of positive scalar scaling degree and positive derivative degree has a local primary primitive on transverse waves. Here is a direct calculation for the two types of inputs. For a single divergence-free bivector $B$, integrate its derivatives by parts and denote its dual coefficient by the horizontal two-form $L[u]$, paired with $B$. On a single divergence-free trivector the consistency equation sets the pairing with $du\wedge L$ to zero. Taking that trivector to be a divergence and integrating by parts gives $du\wedge dL=0$. The five independent scalar first jets give the Koszul homotopy below, now on the three-form $dL$. Iterating $dL_j=du\wedge L_{j-1}$ lowers scalar degree. At degree zero the positive derivative order kills the field-independent term. At every positive scalar degree horizontal closed forms of degrees one and two are exact by the independent total-momentum Koszul complex. Ascending the chain gives $L=dV+du\wedge J$ with local one-forms $V,J$. The first term vanishes against $B$. The second is the bivector variation of a one-vector functional pairing the vector with $J$. Subtract its classical differential.

The residual $R(v,w)$ is alternating in its two divergence-free vectors. Testing consistency on a bivector and a vector gives $R(\iota_{du}B,w)=0$. Integrating derivatives off $v$ and repeating the one-form descent below gives $R(v,w)=\int v(u)F[w]$, where $F$ is a scalar differential operator linear in $w$ with coefficients in scalar jets. The descent treats $w$ as a spectator. At positive scalar degree the total momentum still contains an independent scalar momentum. At zero scalar degree its horizontal one-form is closed because wedging its derivative with an independent $du$ gives zero. A closed one-form linear in transverse $w$ is horizontally exact: the symbol module in the vector momentum $\xi$ is the quotient of five free components by the single relation $\xi$. Its length-one free resolution and the five-coordinate Koszul complex give zero horizontal cohomology in degrees below four.

Skew symmetry also gives $R(v,\iota_{du}B)=0$, hence $dF[\iota_{du}B]\wedge du=0$. Write $B$ as the divergence of an arbitrary trivector test. On its highest derivatives this equality multiplies the symbol of $F[\iota_{du}B]$ by $\xi\wedge du$. A generic $\xi$ is independent of $du$. Descent in derivative order proves $F[\iota_{du}B]=0$. At the highest order of $F$ its covector symbol must then annihilate $\ker \xi\cap\ker du$, since contraction by $du$ of $\wedge^2\ker\xi$ surjects onto this intersection. Thus it equals $a(\xi)\xi+b(\xi)du$. These coefficients are polynomial: Cramer's rule puts each denominator among the two-row minors of the independent five-component variables $\xi,du$, whose common divisor is one. Subtract the corresponding differential operators on $\operatorname{div}w$ and $w(u)$. Both kill $w=\iota_{du}B$, so descending in order gives
$$
R(v,w)=\int v(u)\,K[w(u)]
$$
for a local scalar differential operator $K$.

This bilinear functional defines a variational two-form on unrestricted scalar variations. Indeed a multilinear scalar differential functional vanishing whenever its arguments are $v(u)$ vanishes identically, for two or more arguments. Integrate by parts off the first argument. Its coefficient $f$ must satisfy $df\wedge du=0$. Express the other vectors as divergences of arbitrary bivectors. The highest-derivative argument just given proves $f=0$, by linearity in those fields. In any remaining argument the highest symbol of $v\mapsto v(u)$ is evaluation on $du$, which surjects onto scalars for $v\in\ker\xi$ at generic $\xi$. Induction removes the restriction in each argument. Apply this to the skew-symmetry identity and then to the three-vector consistency identity. The latter is the variational exterior derivative of $\int h_1 K[h_2]$ evaluated at $h_i=v_i(u)$: variations of the vectors' action on $u$ cancel the vector-bracket terms. It follows that this two-form is closed. If its scalar scaling degree, counting the two variations as scalar fields, is $n>0$, contraction with the Euler variation $h=u$ divided by $n$ gives the local primitive $\int uK[h]/n$. Evaluate this at $h=v(u)$ and restore the one-vector term previously subtracted.

Consequently for any closed genus-four bulk potential every projection with $h=u^{r+1}/(r+1)!$ is a classical coboundary in this scalar-source sector. This result also holds for a potential representing a hypothetical nonzero bulk obstruction. Existence of all of these scalar-source primitives therefore cannot establish exactness of the twelve-wave bulk class. The independent $h$ detector above retains the information lost when $dh$ is constrained to the line of $du$.

### Genus-three lift of every integrated scalar observable

Fix the bulk solution through genus three and a genus-two lift of $F_r$. Write $u$ for the primary scalar, $v$ for its primary vector, and $D_i$ for total holomorphic differentiation in coordinate $i$. At the first possible genus-three marked arity the defect has rank one. Integration by parts on its unique vector puts it in the form $\mathcal A(v)=\int v^i A_i[u]$, with the same integration and volume as the marked functional. Each $A_i$ is a homogeneous polynomial in scalar jets with $n=r+8$ scalar factors and 31 derivatives. The scalar momenta are independent after eliminating the vector momentum by conservation. At the next arity the square-zero BV identity restricts on transverse primary waves to $d_{\mathrm{CE}}\mathcal A=0$. The free variation on the one-descendant defect vanishes there. Variations of lower-genus defects and of the bulk master defect vanish by the chosen lifts. Rank two allows one bivector or two vectors.

**One bivector.** With a divergence-free bivector $H$, fusing a scalar into it gives $\iota_{du}H$, where $du$ is the horizontal one-form with components $D_i u$. This includes each scalar choice once by polarization. Taking $H$ to be the divergence of any three-vector and integrating by parts shows $d(du\wedge A)=0$, with $d$ the horizontal exterior derivative. Thus $du\wedge dA=0$ as an identity in the unrestricted scalar jets. The five first-jet variables $D_i u$ are independent. Their Koszul complex for exterior multiplication by $du$ is exact in form degrees below five. Explicitly, on a $k$-form homogeneous of degree $m$ in those five variables, the anticommutator of exterior multiplication by $du$ with $\sum_i \iota_{e_i}\partial/\partial(D_i u)$ is multiplication by $m+5-k$. Therefore one can solve successively
$$
A_n=A,\qquad dA_j=du\wedge A_{j-1}\quad(1\le j\le n)
$$
with $A_j$ a one-form of scalar degree $j$ and the same derivative degree as $A$. At every step $d^2=0$ gives the next Koszul kernel condition. The zero-field term $A_0$ vanishes by positive derivative degree.

For $j\ge1$ a horizontal closed one-form of scalar degree $j$ is exact. Polarization identifies its differential with exterior multiplication by the sum of the $j$ independent scalar momenta. The five components of this sum are independent linear coordinates. The same Koszul argument applies and symmetrization restores the scalar-slot symmetry. Put $V_0=0$. Since $d(A_j+du\, V_{j-1})=0$ whenever $A_{j-1}=dV_{j-1}-du\,V_{j-2}$, ascending the chain gives
$$
A=dV_n+f\,du,\qquad f=-V_{n-1}.
$$
The gradient drops out on a divergence-free vector. The scalar polynomial $f$ has $n-1$ factors and 30 derivatives.

**Two vectors.** Take the ordered inputs $v,w$, followed by scalars. The fusion of $v,w$ has sign $+$ and gives $p_w(v)w-p_v(w)v=[v,w]$. The fusion of $v$ with a scalar has sign $+$ since moving $v$ past $w$ and moving the odd divergence past $w$ each gives a minus sign. The fusion of $w$ with a scalar has sign $-$ from the preceding $v$. Their Ward equation is
$$
\delta_v\mathcal A(w)-\delta_w\mathcal A(v)+\mathcal A([v,w])=0,
\qquad \delta_v(D_Ju)=D_J(v(u)).
$$
Here $v(u)=v^iD_i u$ and $J$ is the differentiated coordinate multi-index of the scalar letter. The bracket of transverse divergence-free waves is divergence-free, so the gradient part drops out in all three terms. Write $\ell_f$ for the Fréchet derivative of the jet polynomial $f$, applying $\ell_f$ to a variation $h$ by replacing each scalar jet $D_Ju$ in turn with $D_Jh$. Varying the distinguished scalar in $w(u)$ and $v(u)$ adds $f\big(w(v(u))-v(w(u))+[v,w](u)\big)=0$. Varying the other scalar factors and integrating by parts gives
$$
\int w(u)(\ell_f-\ell_f^\dagger)v(u)=0 ,
$$
where the dagger is the formal differential adjoint under the same integral.

Use divergence-free vectors with components $(-D_2 h,D_1 h,0,0,0)$ and $(-D_2 k,D_1 k,0,0,0)$ for arbitrary test jets $h,k$. Then $v(u)=Dh$, $w(u)=Dk$ for the scalar operator $D=(D_2 u)D_1-(D_1 u)D_2$. Integrating by parts against arbitrary $k$ forces the operator identity $D^\dagger(\ell_f-\ell_f^\dagger)D=0$. Equivalently the independent jets of $h$ test every differential coefficient after the derivatives on $k$ have been removed by momentum conservation. The polynomial jet ring is an integral domain. Products of nonzero scalar differential operators over it have nonzero principal symbol, and the symbol of $D$ is nonzero. Consequently $\ell_f=\ell_f^\dagger$.

**Primitive.** Denote by $\delta/\delta u$ the Euler derivative of a local scalar density, obtained by integrating every derivative on its variation by parts. Homogeneity of $f$ gives $\ell_f(u)=(n-1)f$. Hence
$$
\frac{\delta}{\delta u}\left(\frac{u f}{n}\right)
=\frac{f+\ell_f^\dagger(u)}{n}=f,\qquad
\delta_v\int\frac{u f}{n}=\int f\,v(u)=\mathcal A(v).
$$
Subtract this primary scalar functional from the action linear in the marked source. One may average its finite-degree coefficient over a compact real form of $SL_5$ to preserve invariance. The remaining gradient has a one-descendant free primitive. Subsequent arities have positive descendant level and the free resolution above completes the lift. This constructs a genus-three lift of each $F_r$, including $F_0$. Rank $6-5k$ removes mixed local defects of source degree $k\ge2$, giving the simultaneous genus-three deformation. The proof also applies at every positive genus for which the chosen bulk solution and the lower-genus marked lifts exist, with respective derivative degrees $10g+1$ and $10g$. The construction gives the marked primitive directly from its defect. [Local jet inversion and an exact replay](marked_consistency.py)

At positive derivative degree the scalar primitive is unique modulo total derivatives. Indeed, if a scalar density $B$ has zero boundary, its Euler derivative $f$ obeys $f\,du=dV$ by transversality, so $du\wedge df=0$. If $f$ has highest jet order $k\ge1$, the coefficients of the independent order-$(k+1)$ jets in $(D_1u)D_2 f-(D_2u)D_1 f$ have polynomial symbol
$$
\big((D_1u)\xi_2-(D_2u)\xi_1\big)
\sum_{\substack{|J|=k\\\rm jet\ multi\text{-}indices}}
\frac{\partial f}{\partial(D_Ju)}\xi^J=0 .
$$
The first factor is nonzero in the polynomial ring of the indeterminates $\xi_i$ and all jets through order $k$. Descent on $k$ therefore gives $f=f(u)$, which is zero in positive derivative degree. Integration by parts of the homogeneous field count gives $n\int B=\int u f=0$. At $F_{11}$ the dimensions of the full one-vector target, longitudinal quotient already taken, and of the scalar density space are $36722-5454=31268$ and $637$. The consistency map has kernel dimension $637$ and rank $30631$.

For a normalization check the known genus-one quartic correction $u\det(\partial_i\partial_j u)/24$ has Euler derivative $\det(\partial_i\partial_j u)/4$, by the divergence-free Hessian cofactor. On the five coordinate scalar momenta followed by $(1,-2,2,-1,1)$ choose the balancing vector momentum $k=(-2,1,-3,0,-2)$. The bare hexagon and the vector contraction of this Euler derivative have respective covectors
$$
\frac{(-180,72,-288,0,-180)}{24},\qquad
\frac{(12,12,36,0,12)}{24}.
$$
Their sum is $84k/24$ and vanishes on every transverse vector for that fixture. [Marked hexagon and analytic rank replay](marked_consistency.py)

For $F_0$ the genus-three primary counterterm can in fact be chosen zero. The genus-one and genus-two bilinear corrections each retain a descendant. On a positive-loop primary graph the classical bilinear has either two internal ends and zero injected momentum, or a single bridge. The first possibility gives equal heat columns and the second gives the vanishing bridge cofactor, also in the presence of another quantum vertex. Only its tree with the nine-primary bulk trace remains. With the vector in slot $i$ the eight choices of the other bilinear end give the same nine-slot trace symbol $C$, multiplied by their respective scalar momenta. Their sum is $-p_i C$ by conservation. In the sign convention of the genus-two bilinear, a scalar descendant of symbol $-C$ supplies $+p_i C$. Hence the integrated bilinear's leading genus-three correction is fixed by the same eighteen bulk periods. [Bilinear cofactor and two-loop sign](integrated_square.py)

### Exact partial cotree sums for sparse three-loop bulk probes

There are two integral ten-wave fixtures on which the eighteen trace boundary coefficients are detected by single coordinate differences between transverse potentials. The first has three parallel pairs, the second two. Seventeen rows from the first and one from the second have determinant $963\pmod{1009}$. [Momenta, projections and replay](sparse_trace_rows.py)

Fix the six cubic core path resistances and the positions of the retained waves, and insert two waves of strengths $a,b$ in one coordinate on reduced edges $e,f$. Let $C_e,C_f$ be their three-component cycle columns. Subtract consecutive columns of the sixteen-edge current determinant at each insertion. The two jump columns in this coordinate are $a C_e,b C_f$. Their exterior product removes two of the three coordinate rows. The residual row at a reduced edge $l$ is
$$
ab\,\det(C_e,C_f,C_l)\,I_l .
$$
Here $I_l$ is the current in that coordinate. The residual row vanishes at $e,f$. On every other edge use the original reduced network, replacing a source at fraction $s$ by endpoint charges with weights $1-s,s$. Thus the row is affine in the two insertion positions and independent of all other deleted coordinates. A single wave of strength $a$ leaves two residual rows instead. Choose two vectors perpendicular to $C_e$ with exterior product dual to $C_e$ up to sign, and multiply their two-row minors by $a$. The cycle columns have entries $0,1,-1$, so such rows are obtained by elementary pivoting at any nonzero component.

On deleting a final reduced edge column $l$, the projective cofactor equals the resulting square determinant divided by its resistance $T_l$, up to sign. This follows from the weighted cycle-current relation, which puts the vector of edge times in the current matrix kernel. Each independently removed color supplies the outer product of its row or of its two-row minors. Average that matrix over its own edge assignments with measure $T_e\,ds$ for each labeled insertion. The squared determinant is obtained by wedging these matrices on both index sets and extracting the top coefficient. The scalar Green kernel between fractions $s,t$ is its bilinear interpolation from the endpoint Green matrix plus $T_e(\min(s,t)-st)$ if both points are on edge $e$. For a transverse coordinate probe the target and driven colors are distinct. The coupled bridge contribution therefore uses just the two distinguished single-position conditional outer products. Their entries have degree at most four in the distinguished fraction, so integrating their product against this kernel uses a fixed five-by-five matrix of rational moments.

Thirty placements on a four-pair source give the original sixteen-edge cofactor times all four independently computed voltage differences, including fourteen nonzero placements. On a three-pair source twenty further placements pass the same original cofactor check. At the fixed reduced edge times $2,3,\ldots,11$ the three-pair position integral in the replay is
$$
-\frac{6604126224963280696502022854122480}{22071023253620227236325964673}.
$$
Its residue is $180\pmod{1009}$. [Exact rational placement and voltage tests](partial_cotrees.py)

The native evaluator sums every remaining path assignment and integrates the retained positions. At fixed core times the squared cofactor has total position degree at most four on any core path by the path determinant bound. The voltage adds at most two after integrating a removed coordinate against the same-edge bridge. For two removed targets the difference of potentials is therefore a polynomial of degree at most six in each retained position on the whole cube: the only minimum involving a retained coordinate pairs it with the integrated target. One may ground the voltage at a fixed branch to see this directly. A four-point finite-field rule is obtained from four distinct roots of $P_4+aP_3$, with $P_j$ the shifted Legendre polynomials. Division by this polynomial and its orthogonality to quadratics give exactness through degree six. Distinct choices of $a$ give disjoint grids for different labels. The evaluator checks every prescribed moment. A retained target with a retained source on the same path is instead integrated on the ordered chambers.

The executable takes a prime, a core, a fixture name and the zero-based triples $(i,j,k)$ denoting $(\Phi_i-\Phi_j)(e_k)$. Each input line has six path resistances. Each output line has their unaveraged core numerators, one per probe, with the common denominator $U^{12}$. Equivalently it outputs $U^7$ times the position integral of the squared cofactor and voltage. The fivefold heat measure supplies $U^{-5}$. The core symmetry divisors are $24$ for the tetrahedron and $16$ for the double square. The genus-two trace insertion in the hexagon is evaluated separately. With the four-pair fixture on the tetrahedron and times $(1,2,3,4,5,6)$, the projections $(2,j,0)$ for $j=3,4,5,6$ have residues $(840,211,390,284)\pmod{1009}$. [Complete finite-field core sum](partial_cotrees_mod.cpp)

For the last projection one may instead use a leading soft coefficient. Take the five coordinate covectors in their usual order and use the ten inputs
$$
e^1,2e^1,e^2,3e^2,e^3,e^4,e^5,
\epsilon(e^3+e^4),e^3+2e^4+e^5,p_9 ,
$$
where $p_9$ balances the first nine. Measure the coefficient of $\epsilon^2$ in $(\Phi_2-\Phi_3)(e_1)$. Each current determinant has a jump column $\epsilon C_e\otimes(e^3+e^4)$. At degree two in its square all currents and the voltage are evaluated at $\epsilon=0$. The soft insertion's position disappears. Its remaining edge sum is the contraction with $\sum_e T_e C_e C_e^{\mathsf T}$, evaluated by a rank-three symmetric factorization. A colored current perturbation lies in the three cycle rows of its direction, of which the jump already uses one, so the original cofactor has degree at most three in the strength. Squaring and multiplying by the voltage gives degree at most seven. At core times $(1,2,3,4,5,6)$ the six finite strengths $1,\ldots,6$ give the respective residues $102,290,534,921,183,750$ modulo 1009. Dividing out the quadratic zero and interpolating gives the constant $97$, also obtained by the soft edge sum.

In the trace variation the soft wave pairs to zero with the observed $e_1$. It therefore remains among the nine trace arguments in every contributing fusion. Keep precisely the terms in which this slot occurs twice, substitute $e^3+e^4$ there, and set $\epsilon=0$ in the balancing covector. The resulting integer trace row, together with the first eleven exact rows and the chosen four and two ordinary projections, has determinant $515\pmod{1009}$. The seven-field genus-two insertion contributes $-254/6075$ to this row.

An independent normalization control takes the earlier $3,2,2,1,1$ multiplicities. For a color with three waves the jump block fills its three cycle rows, giving its weighted cotree determinant directly. Both complete projections then reproduce the previously integrated rational and $\zeta(3)$ coefficients modulo 1009. In table order the tetrahedral pairs are $(471,778),(576,965)$ and the double-square pairs are $(271,0),(1005,0)$. [Exact soft row and hexagon](soft_trace_measurement.json) · [Independent complete-core controls and checkpointed integration](soft_trace_measurement.py)

For $|\epsilon|=1$, a tree flow across a cut is its value at zero plus the soft contribution. In each of coordinates three and four its absolute value is at most the hard positive-entry sum plus one. The other coordinates have their hard bounds. These are precisely the five positive-entry sums at $\epsilon=1$. The observed voltage is in coordinate one and hence independent of $\epsilon$. The existing determinant and two-forest bound therefore holds uniformly on the unit circle. Cauchy extraction of the quadratic coefficient preserves this bound, by the triangle inequality under the contour integral. Path integration acts coefficientwise on the integral strength polynomial and introduces the same factorial denominators as at fixed strength. In rational and $\zeta(3)$ channel order, the resulting numerator bounds have bit lengths $(392,354)$ on the tetrahedron and $(365,0)$ on the double square. [Denominators and height calculation](../period_bounds.py)

For a single removed wave, its two-row minor is quadratic in its position $s$. The three polynomials $1$, $s-1/2$, $(s-1/2)^2-1/12$ have diagonal Gram matrix $\operatorname{diag}(1,1/12,1/180)$ under uniform integration on $[0,1]$. The averaged outer product is therefore the sum of three coefficient outer products with these weights. For a parallel pair the residual row is affine in its two positions. Holding either position fixed gives the outer product of the conditional mean plus $1/12$ times the slope outer product for the other position. The conditional quadratic coefficients can be summed over partner edges before the voltage integration. [Implementation and graph controls](partial_cotrees_mod.cpp)

For a sharper rational reconstruction bound use the original integral $24$-column determinant. At quadratic soft order its two soft jump columns supply both powers of the strength. All particular flows, the balancing jump and the voltage therefore use the hard source at strength zero. Each component of total positive charge $H_i$ admits a particular flow of absolute size at most $H_i/2$. Match the positive charges to the negative charges and split every transfer equally between two edge-disjoint paths, which exist on either bridgeless subdivided core. Changing the particular flow adds a cycle flow and preserves the determinant.

For fixed core-path occupancies $m_e$ with total ten, subtract consecutive constant columns along each path. The remaining six base columns have the form $C_e\otimes(C_e,y)$ with $y$ the particular flow. The ten jumps have the form $C_e\otimes p$ in the source rows. Divide the three rows of source component $i$ by $\sqrt L\,H_i/2$ for any positive $L$. In the squared determinant the row scalings and ten jump norms give the hard-source factor
$$
L^5\prod_{\substack{i=1,\ldots,5\\\rm components}}(H_i/2)^6
\prod_{\substack{p\\\rm jumps}}\left(\sum_{\substack{i=1,\ldots,5\\\rm Euclidean\ components}}\frac{|p_i|^2}{(H_i/2)^2}\right)
\prod_{\substack{e\\\rm paths}}\|C_e\|^{2m_e}.
$$
Each base column has squared norm at most $\|C_e\|^2(\|C_e\|^2+5/L)$. For a fixed cycle row $r$, the variable columns indexed by $s$ are sums of $t_l C_{l,s}$ times the vector with $(C_l,y_l)$ in block $r$. Expand their exterior product. Its coefficient on a selection of segments is the corresponding cycle minor. For the full three columns or either remaining pair, its coefficient one-norm is bounded by the sum over distinct paths of the absolute minor times the product of $(m_e+1)\sqrt{\|C_e\|^2+5/L}$ on these paths. Cauchy–Schwarz bounds the square by the product of two integer-weighted sums, with zero and two powers of each square root. Two blocks use three columns and the third uses two. Permuting the cycle basis allows any one diagonal column to be omitted. Apply Hadamard to each selected tuple of segment vectors and then sum their absolute contributions. This bounds the full coefficient one-norm by the product of the three block bounds and the constant-column norms. The implementation minimizes over this choice and $L=1,2,3,4,5,8$ for each occupancy.

A single segment parameter in the current determinant contributes variable columns in a three-dimensional space containing its constant column. The determinant thus has segment degree at most two. The square times the separating two-forest has segment degree at most five. If the resulting core exponent is $k_e$, the segment exponents on path $e$ sum to $k_e-m_e$. The product of their factorials is at most $(5!)^{\lfloor(k_e-m_e)/5\rfloor}((k_e-m_e)\bmod5)!$, by concentrating exponents up to five. Multiply the bound by the absolute core moment divided by $\prod k_e!$, maximizing over the feasible exponent tuples for each occupancy. There are $10!$ orders at each occupancy and at most $\binom{16}{4}$ two-forests per voltage pair. Finally divide by the core symmetry order. Each squared jump also contributes the square of the common integer divisor of its five components. Its product over jumps divides the scaled integer numerator and reduces the proven common denominator accordingly. [Exact occupancy, minors and denominator calculation](../ordered_period_bounds.py)

The resulting rational and $\zeta(3)$ numerator bounds have bit lengths $(235,199)$ on the tetrahedron and $(194,0)$ on the double square. The independent prime $1009$ gives $(911,830)$ and $(981,0)$ respectively. On the double square ten further primes give a 200-bit modulus and the exact pair $(-159426073/3919104000,0)$.

The full core numerator admits a smaller reconstruction bound. On a core path with several insertions, replace its segment columns by their resistance-weighted mean and the consecutive jumps, with determinant of the change of columns equal to one up to sign. The current in the mean column solves the core flow problem with each source of fractional position $s$ split between its branch endpoints with weights $1-s,s$. This follows by summing the voltage drops and the endpoint current conservation equations along that path. The mean columns obey the core-time kernel relation. Dropping one mean column leaves five mean columns and all ten jumps. The mean currents are affine in all insertion positions together. For a single varying position its perturbation lies in the three cycle rows of its source direction, one of which is already used by that source's jump. Consequently the determinant has total position degree at most five and single-position degree at most two. Changing the ordering of the jumps changes only its sign. Multiplying the square by $U^6$ gives a polynomial with integral coefficients on each ordered chamber by the original current determinant, and the mean-column formula gives the same polynomial on the whole cube.

For a transverse voltage, the two observed slots carry zero charge in its component. Grounding at a branch expresses $U$ times the voltage by the polynomial core Green numerator, bilinearly interpolated in two distinct slots. On a shared path the additional term is the core-time factor times $\min(s,t)-st$, with the same integral $U$. Thus the squared cofactor times the polynomial voltage term has total position degree at most twelve and at most two slots of degree five. For the remaining minimum term the polynomial has total degree at most ten and slot degree at most four. Its only coupled position integral is
$$
\int_0^1\int_0^1 s^a t^b\min(s,t)\,ds\,dt
=\frac{1/(a+2)+1/(b+2)}{a+b+3}.
$$
Taking the least common multiple of the denominators of these monomial integrals gives $14370048000$. Each core coefficient divided by the squared integer jump content therefore has denominator dividing both this integer and twice the six core-exponent factorials multiplied together. The occupancy bound above applies coefficientwise by maximizing the segment-factorial product over the given exponent orbit. After using the intersection of the two denominator bounds, the largest scaled numerator on either core has at most 155 bits for the quadratic-soft probe and 148 bits for the ordinary three-pair probes. Eight primes suffice in both cases. The replay reconstructs every orbit coefficient, compares them all at the independent prime 1009 and integrates both period channels exactly. [Polynomial reconstruction and divisibility replay](../core_numerator_reconstruction.py)

### Exact nine-field bulk trace at genus three

In the basis order of the eighteen listed sextic words, write each coefficient of $C$ as its rational entry plus its second entry times $\zeta(3)$. The action contains $-C$ in the convention $\Phi=C$. The complete $18\times18$ rational system multiplies back to both period columns. Its determinant reduces to $515$ modulo $1009$, where each right-hand side and each coefficient agrees with the independent cubature. The first pair is
$$
\left(-\frac{185964611864833}{23664402472122187776000},\;
       \frac{28870230551}{4564892452184064000}\right).
$$
For $C$ the signs in basis order are $(-,+,+,+,-,-,-,+,-,+,+,+,+,+,+,+,-,+)$. Ten terms of Euler's transform for the alternating zeta series certify them with the rational endpoints in the coefficient file. The residual of the $N$-term lower bound is less than $4/(3\,2^N)$ because the $k$th finite difference of $1/(j+1)^3$ equals the integral of $(1-x)^k(\log x)^2/2$ on $[0,1]$. [All basis words, diagram measurements, exact matrix, coefficients and signs](three_loop_coefficients.json) · [Exact solve and readback](../trace_reconstruction.py)

Across the seven new transverse projections, the greatest common divisor of the unaveraged degree-thirty core numerators is $U^2$ for the tetrahedron and
$$
U^2(T_0T_1T_4T_5)^2(T_0+T_1)(T_4+T_5)
$$
for the double square, with the zero-based path order of the native evaluator. Dividing each numerator by its displayed factor gives a polynomial over the rationals. Two quotient polynomials on each core have greatest common divisor one modulo $1009$. A common nonconstant homogeneous primitive integer factor would retain its degree under that reduction. Integrating all fourteen numerators after canceling $U^2$ gives the same fourteen ordered period pairs with Kirchhoff denominator $U^{10}$. [Core factors and lower-denominator period replay](../core_numerator_reconstruction.py)

### Conserved symbols of the integrated scalar deformations

Take the integrated $F_s$ defined by a constant level-$s$ top-polyvector source in the full flat theory, in worldsheet action units. The regular scalar waves are free-closed. Use their translation-invariant polarized symbols, with a top Dolbeault polarization for the scalar integration. A density that is a total derivative has zero symbol at conserved momenta.

These scalar symbols descend to the quantum cohomology of the translation-invariant holomorphic local functional complex in the fixed flat theory. To see this, polarize a candidate primitive of total input polyvector rank $R\ge0$. In a local graph for its differential, write $E$ for the number of internal edges, $V$ for the number of bulk vertices in addition to the marked primitive, and $h=E-V$ for the geometric loop number. With one heat edge of Dolbeault degree five and $E-1$ propagators of degree four, relative position integration requires $4E+1=5V$. Thus $E=5h+1$, $V=4h+1$. Each bulk vertex has total rank five and every edge consumes rank four. The external fields must carry rank
$$
5V+R-4E=R+1 .
$$
Consequently this term vanishes on scalar fields. This includes vertices of positive worldsheet genus. A coincident free BV self-contraction vanishes on holomorphic jets by the diagonal heat-kernel calculation. The free divergence raises the required input rank by one, and the free Dolbeault differential vanishes on the scalar tests. Total derivatives vanish by momentum conservation. This argument uses the ordinary exterior product of Dolbeault components, with the top form on the external scalar fields. In particular a nonzero conserved scalar coefficient of $F_s$ cannot be removed by a quantum-exact local change of representative, even using a primitive with lower-genus terms. The statement holds at fixed bulk action and field coordinates. [Local graph rank count](#simultaneous-sources-and-the-weak-coupling-effective-action)

At one loop use the primary density $K(u)\det H$, where $u=A_0$, $H_{ij}=\partial_i\partial_j u$ and $K(u)=u^{s-1}/[24(s-1)!]$ for $s\ge1$. Each divergence of the Hessian cofactor vanishes by antisymmetry of the determinant and commutation of the three derivatives. Integrating twice by parts and using $\operatorname{adj}(H)^{ij}H_{ij}=5\det H$ therefore gives
$$
\frac{\delta F_s^{(1)}}{\delta u}
=6K'(u)\det H+K''(u)(\partial_i u)\operatorname{adj}(H)^{ij}(\partial_j u).
$$
For $s=1$ the density is a total derivative, since $\partial_i((\partial_j u)\operatorname{adj}(H)^{ij})=5\det H$. For $s=2$, on $u=B+f(z^1)$ with $B=((z^2)^2+(z^3)^2+(z^4)^2+(z^5)^2)/2$, the Euler derivative is $\partial_1^2f/4$. Pairing with a second wave of opposite axial momentum gives $p^2/4$ for $f=e^{p z^1}$, as a coefficient per unit volume. More generally the coefficient of $b^{s-2}$ linear in $f$ at $u=b+B+f$ is $\partial_1^2 f/[4(s-2)!]$.

For an independent integrated test take the six primary momenta $pe^1,t_2e^2,t_3e^3,t_4e^4,t_5e^5$ and their negative sum. The level-two source has momentum zero. A nonzero elliptic partition must pair this source with one of the six primaries. The resulting genus-one cotangent integral is $1/24$, and each complementary squared five-row minor equals $p^2(t_2t_3t_4t_5)^2$. The sum is $p^2(t_2t_3t_4t_5)^2/4$. Its quadratic coefficient in all four transverse momenta survives momentum conservation. At unit strengths the independent fivefold Ward rows, after deforming the source by $j(e^2-e^3)$ and compensating at the first scalar, give $1/3,7/12,1$ for $j=1,2,3$. Quadratic interpolation gives $3/3-3\cdot7/12+1=1/4$.

Keeping all six hard momenta, the first hard scalar can also carry a descendant. At unit strengths its level-one six-field elliptic vertex is $1/24$: the six singleton blocks place the genus-one cotangent class on the first hard mark. For $F_0$ put level two on that hard scalar. The string equation at the zero-momentum primary source gives $1/24$. For $F_1$ put level one on the hard scalar. The zero-momentum level-one source gives the dilaton factor $2-2+6=6$ and the value $1/4$. Both are nonzero on conserved momenta. Scaling the five simplex directions multiplies all three integrated results by $p^2(t_2t_3t_4t_5)^2$.

At every $g\ge2$, choose the $2g+3$ conserved momenta of a nonzero all-primary fivefold vertex $V_{g;0}$, and replace all its regular polarizations by scalars; the constant source already supplies the full polyvector trace. For $s\ge1$ add $s-1$ labeled zero-momentum primary scalars alongside the zero-momentum level-$s$ source. Each string equation lowers only the source. The last level-one zero-momentum insertion gives the dilaton factor $2g-2+(2g+3)=4g+1$. On a repeated constant scalar background $b$ the integrated coefficient is
$$
\frac{b^{s-1}}{(s-1)!}(4g+1)V_{g;0}.
$$
For $F_0$ take the same hard momenta, give exactly one hard scalar level one and use the zero-momentum primary source. The string equation gives $V_{g;0}$. Thus both measurements survive the total-derivative quotient at every such genus. At genus two the unit fivefold basis, $e^1+e^2$ and their negative sum give $V_{2;0}=17/155520$. The level-one source gives $17/17280$, and putting the level instead on the first hard scalar gives $17/155520$. These are coefficients of the chosen integrated deformation. [Nonzero planar family in every genus](../README.md#an-elementary-scalar-period-at-each-genus) · [Elliptic and genus-two replay with the integrated option](../protected_source_tower.py)
