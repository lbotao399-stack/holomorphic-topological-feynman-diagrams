# Local loop corrections in minimal BCOV theory on the fivefold

Use the flat minimal BCOV model proposed for twisted type IIB supergravity. The primary scalar field is $\phi$. The propagating fields have descendant level zero. The primary interaction is $\mathrm{Tr}(\alpha^3)/3!$ and the closed propagator applies holomorphic divergence to the Dolbeault kernel. At a vertex of valence $l$ the descendant levels sum to $l-3$. Thus only the cubic contributes on primary external fields. All currents below use the repository's orientations and incident vertex momenta. [Costello–Li, §§4–5](https://arxiv.org/pdf/1606.00365#page=13)

## All-genus primary vertex from the worldsheet

For the resulting quantum cohomology, the protected point operators $[E_r]_5$ exist at every descendant index. A fixed level-one point operator has a nonzero measured coefficient at every genus on the stated quadratic and level-three background. The integrated scalar deformations $F_s$ have nonzero conserved responses at every positive genus; the scalar rank count makes these responses invariant under exact local changes of representative at fixed bulk action and coordinates. [Point-source calculation](operator/README.md#one-protected-operator-with-a-nonzero-coefficient-at-every-genus) · [Integrated classes and exact changes](operator/README.md#conserved-symbols-of-the-integrated-scalar-deformations)

Li and Wang announce the flat descendant BCOV quantization in Theorem 1.11. The genus-zero interaction in their Lemma 1.3 is the descendant action used here. Their $Q_0=\bar\partial-\zeta\partial$ preserves the minimal subcomplex $i+k\le4$ after $t=-\zeta$, and their BV kernel contracts two primary ends of total polyvector rank four. Hence restriction commutes with every contraction and with the effective master equation. For the primary genus-$g$ vertex on $\mathbb C^5$, equations (1.9) and (1.10) give $n=2g+3$ conserved input momenta $P$ and the polarized coefficient
$$
V_{g;0}(P)=\int_{\mathcal M_{g,2g+3}}\frac{\Theta_g(P)^{5g}}{(5g)!}.
$$
Here $\mathcal M_{g,n}$ is the moduli space of smooth genus-$g$ curves with $n$ labeled marked points. Their Lemma 1.2 writes $\Theta_g$ as the sum of the symplectic forms on the five sets of Abel–Jacobi period coordinates. The action divides the polarized vertex by $n!$. The preprint gives proof sketches and announces a fuller treatment. [Li–Wang, §§1.1–1.6](https://arxiv.org/abs/2610.21927) · [Costello–Li, minimal fivefold theory, §4](https://arxiv.org/abs/1606.00365)

The ratio of our output trace to theirs is $-(2\pi i)^5$ since $\Omega\wedge\bar\Omega=-32i\,d^{10}x$ and their trace factor in dimension five is $\pi^{-5}$. Their heat time is four times ours. Its Jacobian cancels the inverse rescaling of the gauge-fixing operator. Each edge uses the inverse trace; each vertex supplies one trace. Removing the output trace gives $[-(2\pi i)^5]^{-g}$ on a bare connected $g$-loop graph. The same conversion holds at total genus after assigning the factor at each positive-genus vertex.

In (1.17) write $Y_{e,(r,i)}=C_{re}j_e^i$ for the cycle rows $C$ and minimizing current $j$. At projective lengths $\lambda_e$ summing to one let $M=C\operatorname{diag}(\lambda) C^{\mathsf T}$. The gradient of the edge energies with respect to the positive lengths is $-Y(1_5\otimes M^{-1})Y^{\mathsf T}$ in component-first column order. Restrict the derivatives to the columns $e_j-e_0$, $j>0$. Kirchhoff's law gives $Y^{\mathsf T}\lambda=0$, so the signed maximal minors of $Y$ are $c\lambda_e$ with common scalar $c$ (all zero at deficient rank). Expanding the flux determinant twice gives $c^2(\det M)^{-5}$ times the voltage paired with the odd sources. The $E-1$ minus signs from the derivatives cancel the orientation sign in (1.17). This includes the one-edge tree. Literal differentiation of the grounded Laplacian and an independent cycle determinant give ratio $+1$ on nonzero graphs at genera one through four. [Exact flux check](worldsheet_primary.py)

For real $P$, normalize the curve's period matrix $\tau$ by the $A$-cycles. Its imaginary part is positive definite. For each target coordinate the Abel–Jacobi lift $x+\tau y$ is locally holomorphic in the curve and its markings. For nonintegral real sources this follows by real linearity from the divisors $z_i-z_n$. Thus $dx+\tau\,dy=d(x+\tau y)-(d\tau)y$ is of type $(1,0)$. The canonical form in Lemma 1.2 is
$$
\sum_{\substack{\mu=1,\ldots,g\\\mathrm{period\ coordinates}}} dx_\mu\wedge dy_\mu
=\frac{i}{2}(dx+\tau dy)^{\mathsf T}(\operatorname{Im}\tau)^{-1}
                 \wedge\overline{(dx+\tau dy)}.
$$
Its top power on the moduli space is a squared absolute determinant times the positive volume density. Hence $V_{g;0}(P)\ge0$ for real momenta. If one momentum is scaled by $\varepsilon$ and another balances it, both real tangent directions moving only the soft marking carry $\varepsilon$ in the top determinant. This proves quadratic softness pointwise; polynomiality and convergence in Lemma 1.6 give the integrated identity. At $g=1$ five conserved waves have rank at most four, so Lemma 1.2 gives $V_{1;0}=0$.

On a single scalar letter only the one-edge tree with this vertex survives. In the repository's conventions $V_{2;0}=-P$ and $V_{3;0}=-C$: the lower-genus primary vertices agree inductively, and the transverse tree variation is injective on doubly soft primary traces. In general, with $p_{\rm root}$ the absorbed momentum and $v$ the vector polarization, the primary matrix element is
$$
Q_g\phi(0)=-\frac{(-1)^g}{(2\pi i)^{5g}}\,p_{\rm root}(v)\,V_{g;0}(P).
$$
The string expansion multiplies this coefficient by $\lambda^{2g}$. Together with the nonzero saturated planar period at each $g\ge2$, the announced quantization gives a nonzero elementary-field correction at every such genus. Moreover the reduced trace has one sign on all real conserved momenta and is strictly positive in the $V_{g;0}$ convention on a dense open set. [Elementary-field probes and strict period sign](#an-elementary-scalar-period-at-each-genus) · [Uniqueness of the primary trace](operator/README.md#injectivity-of-the-transverse-trace-variation-at-every-positive-genus)

### Complete one-loop bulk interaction with moving descendants

Use the descendant coefficients $A_k$ of (1.10), with primary field $A_0$, and the existing formal classical field $U$ with zero constant term,
$$
U=\sum_{k\ge0\,;\,\mathrm{descendant\ levels}}A_k\frac{U^k}{k!},
\qquad
\frac{\partial U}{\partial A_0}
=\left(1-\sum_{k\ge0\,;\,\mathrm{descendant\ levels}}A_{k+1}\frac{U^k}{k!}\right)^{-1}.
$$
Here the derivative holds every positive-level field fixed. In flat affine holomorphic coordinates normalized by the volume form, $i,j=1,\ldots,5$ label the coordinate derivatives. In worldsheet units the full flat genus-one bulk interaction of the even generating superfield is
$$
I_1=\frac1{24}\operatorname{Tr}\!\left[
\log\!\left(\frac{\partial U}{\partial A_0}\right)
\det\!\big(\partial_i\partial_j U\big)\right].
$$
The trace integrates the top polyvector and Dolbeault component. Restrict to the minimal theory by retaining components of polyvector rank plus descendant level at most four. The coefficient linear in $A_1$ is the previously verified $\operatorname{Tr}[A_1\det(\partial_i\partial_j A_0)]/24$. [Six-field diagram](operator/README.md#the-genus-one-bulk-term-on-the-one-descendant-sector)

For target dimension $N$, take $n$ labeled conserved momenta $p_i$, with descendant levels $k_i$ satisfying $\sum_{i\,;\,\mathrm{marks}}k_i=n-N$. Partition the marks into one distinguished nonempty block $B_0$ and $N$ unordered nonempty blocks $B_1,\ldots,B_N$. In $B_0$ the levels sum to $|B_0|$; in every other block they sum to $|B_j|-1$. Let $p_{B_j}$ be the sum of its marked momenta. The labeled vertex is
$$
V_{1;\boldsymbol k}(P)=
\sum_{\mathrm{admissible\ partitions}}
\left(\int_{\overline{\mathcal M}_{1,|B_0|}}\prod_{i\in B_0}\psi_i^{k_i}\right)
\prod_{j=1}^N\frac{(|B_j|-1)!}{\prod_{i\in B_j}k_i!}\,
\det(p_{B_1},\ldots,p_{B_N})^2 .
$$
All block weights are positive. For real momenta this vertex vanishes exactly when every admissible partition has dependent ordinary block momenta, including the case of no admissible partition. Every genus-one primary bulk vertex vanishes by this criterion. For a single level-one descendant the result is the squared determinant of the other $N$ momenta divided by $24$.

To derive the partition formula use $\psi_i=\lambda+\sum_{i\in S,\ |S|\ge2}\delta_{0:S}$ on $\overline{\mathcal M}_{1,n}$, with $\lambda$ the first Chern class of the Hodge line and $\delta_{0:S}$ the rational tail carrying precisely $S$. This follows by forgetting to the elliptic curve with mark $i$ and adding the tails contracted in stabilization. The Hodge class is pulled back from $\overline{\mathcal M}_{1,1}$, where its degree is $1/24$. A tail contributes the genus-zero descendant integral with zero exponent at the new node and momentum $p_S$ on the elliptic component. To justify the divisor substitution with Betti factors, use polarization of Chern–Weil for the Abel section of the Jacobian. The separating-edge current is fixed by $p_S$. Its linear energy yields the ordinary separating-divisor trace; all nonlinear height jumps on an elliptic graph depend only on the lengths of its unique cycle. Subdivisions in those lengths leave the separating smoothing parameter free. Restriction of the divisor extensions to the tail therefore gives the fused Abel section on the elliptic factor and the constant Abel map on the rational factor. The node cotangent terms follow from the normal self-intersection $-\psi_\star-\psi_{\star'}$. The single-node formula is also (1.15) with the conserved separating current. [Biextension Chern–Weil, Theorems 5.20 and 6.3](https://arxiv.org/abs/2112.09007) · [BCOV vertex and boundary restriction, (1.9) and (1.15)](https://arxiv.org/html/2610.21927v1)

Choosing $\lambda$ fixes a smooth elliptic curve. With $N+1$ surviving marks and zero remaining levels, fix one as the origin and integrate the Betti top form over its $N$-fold product. On integral momenta the Abel homomorphism has degree $\det(p_{B_1},\ldots,p_{B_N})^2$ if it is an isogeny, and zero top form otherwise. Polynomiality extends this to arbitrary momenta. On the moving base with no descendants the $N$ conserved rows have rank at most $N-1$. Thus every nonzero reduction ends on one fixed elliptic fiber. Each ordinary cluster of $b$ labels has the fixed-fiber weight $24\int\lambda\prod\psi_i^{k_i}=(b-1)!/\prod k_i!$: use $12\lambda=\delta_{\mathrm{irr}}$ and the two interchangeable branches on $\overline{\mathcal M}_{0,b+2}$. In the exceptional block the single base degree instead gives the ordinary genus-one cotangent integral. More explicitly, reduce any positive $k_i$ first. A tail $S$ through $i$ contributes only when $\sum_{j\in S} k_j=|S|-1$, with its genus-zero multinomial after lowering $k_i$. For a fixed final partition these cuts take place inside the block containing $i$. In an ordinary block the base term has the wrong fiber degree; the remaining recursion is that of the fixed-fiber weight. In the exceptional block both terms give precisely the moving genus-one cotangent recursion. This proves the displayed product weights, including all multiplicities.

Lagrange inversion gives the ordinary block coefficients of $U$. The string operator $\partial_{A_0}-\sum_{k\ge0\,;\,\mathrm{descendant\ levels}}A_{k+1}\partial_{A_k}$ sends $U$ to one and annihilates $\log(\partial U/\partial A_0)$. At $A_0=0$ this logarithm is $-\log(1-A_1)$, whose $b$-linear coefficient is $(b-1)!$. By the string and dilaton equations its coefficients divided by $24$ are exactly the distinguished block integrals. Finally the Hessian of each ordinary block plane wave is a rank-one matrix with momentum $p_B$ in both slots. Cauchy–Binet gives the squared determinant once per unordered collection of ordinary blocks. [Exact evaluator and two orders of rational-tail reduction](elliptic_descendants.py)

On the elementary primary scalar $\phi(0)$, a one-letter mark forces a bridge in every positive-geometric-loop graph. The surviving one-edge bracket of $I_1$ contributes $p_{\rm root}(v)V_{1;\boldsymbol k}(P)/(2\pi i)^5$ in the $A_k$ convention, where $p_{\rm root}$ is the absorbed momentum and $v$ the transverse primary vector polarization. Use $p_{\rm root}=e^1$, vector momentum $e^2$ with $v=e_1$, and the five scalar momenta
$$
e^3,\quad e^4,\quad e^5,\quad e^1+e^2,\quad -2e^1-2e^2-e^3-e^4-e^5 .
$$
For scalar descendant levels in this order,
$$
(2\pi i)^5 Q_1\phi(0)=
\begin{cases}
1/2,&(0,0,0,2,0),\\
7/12,&(0,0,0,1,1).
\end{cases}
$$
Both tuples have even total level and retain their signs under $t=-\zeta$. The independent genus-one graph Ward equation recovers both values from hexagons and lower-valence trees. [Exact point evaluation](elliptic_descendants.py) · [Graph calculation](genus_one_graph_ward.py) · [Variation of the compact action](elliptic_hessian.py)

#### Scalar descendant backgrounds at one loop

Fix holomorphic scalar backgrounds for the positive-level fields and vary the primary superfield. The substitution from $A_0$ to $U$ is invertible as a formal series. With $U$ as variable, write the logarithm in the action as
$$
F(z,U)=-\log\!\left(1-\sum_{k\ge0\,;\,\mathrm{descendant\ levels}}A_{k+1}(z)\frac{U^k}{k!}\right).
$$
The bulk interaction vanishes as a functional of all primary fields, modulo total derivatives, exactly on the backgrounds
$$
A_k=0\quad(k\ge2),\qquad
\partial_i\partial_j\log(1-A_1)=0\quad(1\le i,j\le5).
$$
Indeed the Euler derivative of $24I_1$ with respect to $U$ is
$$
F_U\det(\partial^2 U)+\operatorname{cof}(\partial^2 U)^{ij}\partial_i\partial_j F(z,U(z)).
$$
The Hessian cofactor has divergence zero by symmetry of the third derivatives. At fixed value and gradient of $U$, scale its Hessian by a free parameter. The term of degree five is $6F_U\det(\partial^2 U)$, so vanishing on every primary jet forces $F_U=0$, hence $A_k=0$ for $k\ge2$. At degree four the cofactors of invertible symmetric Hessians span the symmetric matrices, forcing the displayed condition on $A_1$. Conversely an affine function multiplying a Hessian determinant integrates to a total derivative, by two integrations by parts. The necessity test also applies in the minimal theory: every monomial here has at least five primary factors, so two distinct labeled primary slots can carry complementary polyvector ranks one and four in the trace. The positive-level backgrounds are held fixed.

For constant $A_1,A_2$ and six hard primary waves, the dimension rule allows exactly one $A_2$ and any number of $A_1$ insertions. In the coefficient linear in $A_2$, terms multiplying $-\log(1-A_1)$ integrate to zero. The remaining six-field term is
$$
\frac{A_2}{24(1-A_1)^7}\operatorname{Tr}\big[A_0\det(\partial_i\partial_j A_0)\big].
$$
On the unit fivefold simplex put the root at $e^1$, the transverse vector at $e^2$ with polarization $e_1$, and the four scalar waves at $e^3,e^4,e^5,-e^1-e^2-e^3-e^4-e^5$. Each choice of the undifferentiated slot has squared determinant one. Thus the elementary one-loop response on these hard waves is
$$
(2\pi i)^5Q_1\phi(0)=\frac{A_2}{4(1-A_1)^7}.
$$
The replay checks the labeled zero-momentum coefficients with zero through eight additional level-one insertions. [Point evaluator](elliptic_descendants.py)

### A single moving-descendant probe of the elementary differential at every loop order

Use the coordinate covectors $e^1,\ldots,e^5$ and their dual vectors on the flat fivefold. In $Q_g\phi(0)$ extract the coefficient linear in each of the primary scalar wave $e^{z^1}$, the transverse primary vector wave $e^{z^2}e_1$ and the first-descendant scalar wave $e^{z^3}$ in the $A_1$ convention. Supply the additional primary scalar field $a e^{z^4}+b e^{z^5}$. Let $D_g$ be the unit threefold four-point vertex with one cotangent insertion and let $W_g$ be the fourfold primary on $g-1$ equal axial waves with four waves projecting to a unit transverse simplex; $W_1=0$ by rank. For this extracted coefficient,
$$
(-1)^g(2\pi i)^{5g}Q_g\phi(0)
=(ab)^gD_g+\frac{(ab)^{g-1}(a+b)^2}{(g-1)!}W_g.
$$
At $a=1,b=-1$ this gives $(2\pi i)^{5g}Q_g\phi(0)=D_g>0$ for every $g\ge1$ in the Li–Wang flat quantization. All five supplied Fourier modes and their relative amplitudes are fixed as $g$ varies. In particular the first three heat coefficients are
$$
\frac1{24},\qquad \frac1{9720},\qquad
\frac{411065203}{94152554496000}-\frac{3617}{1009008000}\zeta(3).
$$
Here $\zeta(s)$ is the Riemann zeta function. The fifty-graph four-loop coefficient lies strictly between $110899269775/10^{22}$ and $110899269776/10^{22}$. [Threefold normalization and positivity](../kodaira-spencer/README.md#the-one-descendant-four-field-interaction-on-the-flat-threefold) · [Exact four-loop period](../kodaira-spencer/README.md#complete-four-loop-descendant)

The dimension rule fixes $2g$ additional primary waves. If $m$ of these carry $e^4$, the others carry $e^5$. The absorbed root momentum evaluates to $-1$ on $e_1$ at each $m$. Pull the descendant's cotangent class back from the moduli space of the four retained marks. The two forgotten sets supply $g$ vertical period directions each, so only $m=g$ contributes to that term, with labeled coefficient $(g!)^2D_g$. The ordinary cotangent line also has its rational-tail divisors. With a single cotangent insertion, only fusions with one forgotten wave survive. At $m=g$, each of the $2g$ fusions leaves one saturated direction, of ordered Abel degree $g!$, and the primary vertex $W_g$ in the quotient. At $m=g\pm1$, exactly the $g+1$ fusions with the majority direction contribute, each giving $g!W_g$. Longitudinal components of the four retained momenta are removed by a shear on the three independent transverse covectors; the balancing momentum carries the residual recoil. All other multiplicities have more than $g$ unfused primary marks in one direction, forcing a zero determinant. Divide each labeled result by $m!(2g-m)!$ and insert $a^m b^{2g-m}$ to obtain the formula. The projective current integrand for $D_g$ is nonnegative and the explicit ladder gives a positive term at every $g$. The graph replay checks every multiplicity against the full one- and two-loop vertices and every genus-three fused primary in both period columns. [Graph and factorial replay](fixed_moving_scalar.py) · [Positive all-genus ladder](../kodaira-spencer/README.md#one-integrated-ladder-at-every-genus)

One can also fix the entire real primary scalar background to $e^{z^1}+e^{z^4}+e^{z^5}$ and extract only the term linear in the vector and the moving descendant. With $j$ copies of $e^{z^1}$ the root pairing is $-j$. Every real one-descendant vertex is nonnegative, by pairing the closed positive Betti current with the nef cotangent line, and each Taylor factorial is positive. The $j=1$ terms above therefore give the lower bound
$$
(-1)^g(2\pi i)^{5g}Q_g\phi(0)\ \ge\ D_g+\frac{4W_g}{(g-1)!}>0
$$
on this linear response at every $g\ge1$.

The scalar amplitude at $e^{z^1}$ may also be fixed at one. Use the primary scalar field $e^{z^1}+e^{z^4}+e^{z^5}$ and take the part linear in the vector and in the first descendant above. In each multiplicity the negative of the root momentum on $e_1$ counts the copies of $e^{z^1}$. Every real one-descendant vertex is nonnegative by pairing the positive Betti current with the nef cotangent class. Retaining exactly one copy of $e^{z^1}$ therefore gives the all-genus lower bound
$$
(-1)^g(2\pi i)^{5g}Q_g\phi(0)\ \ge\ D_g+\frac{4W_g}{(g-1)!}>0 .
$$
At genus one each of the three scalar momenta occurs once. At genus two all three axes must occur for full rank. The multiplicities are permutations of $(3,1,1)$ and $(2,2,1)$ with Taylor-weighted vertices $W_2$ and $D_2+2W_2$, respectively, by permutation of the three primary axes in the preceding calculation. Within either orbit the root multiplicities sum to five. The one- and two-loop coefficients in these units are respectively $1/24$ and $5D_2+15W_2=47/62208$. [Direct labeled graph sum](fixed_moving_scalar.py)

### Genus-two fivefold vertices with arbitrary moving descendants

Fix the worldsheet action and the original point coordinates. For $n$ conserved fivefold momenta $p_i$ with minimal levels $0\le k_i\le4$, the genus-two labeled vertex $V_{2;\boldsymbol k}$ has $\sum_i k_i=n-7$. It is a degree-twenty polynomial in the momentum components with rational coefficients. Its elementary-scalar matrix element on one transverse primary vector of polarization $v$ is $-p_{\rm root}(v)V_{2;\boldsymbol k}/(2\pi i)^{10}$ when the absorbed root is primary. The following finite recursion uses the flat graph Ward equation of Li–Wang and the seven-field primary vertex. [Ward equation, Proposition 1.8](https://arxiv.org/abs/2610.21927) · [Complete rational evaluator](genus_two_graph_ward.py)

Lower any positive level by one, obtaining a Ward source of total level $n-8$. Its geometric loop numbers are two, one and zero.

1. At geometric loop number two the graph has eleven edges and ten vertices. A bridge gives a zero current minor; rescaling the two cycles independently gives an extra null direction on a figure-eight. The surviving theta has two trivalent endpoints and eight bivalent sites. At a site carrying the labeled set $S$, its levels sum to $|S|-1$ at a bivalent site and to $|S|$ at a trivalent site. Its charge is $\sum_{i\in S}p_i$ and its classical weight is $(\sum_{i\in S}k_i)!/\prod_{i\in S}k_i!$. An empty endpoint has weight one. The four nonzero shapes have path lengths $(1,5,5),(2,4,5),(3,3,5),(3,4,4)$. With both endpoints decorated the ungrounded source adds a constant to the sum of marked currents; the two special columns of the theta determinant are still linear in that sum. Its square times the grounded potential integrates by rational cube moments and the three-conductance simplex integral. Sum endpoint assignments with the original reversal symmetry. [Two-endpoint integrator and direct electrical checks](theta_endpoint_sources.py)
2. At geometric loop number one, one genus-one vertex lies on a hexagon. The other five sites carry the classical weights just given. In every term of the elliptic partition formula the two internal momenta have opposite displacements from their common mean. Each determinant is affine in that displacement, because the perturbation is rank one. Its square integrates against the circle Green function by moments of degree at most two beyond the Green factor. At geometric loop number zero there is one edge, joining either the two known genus-one vertices or a classical vertex to a lower-arity genus-two vertex. At the classical end the sum of levels is $|S|-2$; the fused leg at the other end is primary. [Elliptic sewing](two_loop_elliptic_insertions.py) · [Seven-field primary seed](#the-genus-two-bulk-action-and-integrated-cubic)
3. The summed potential at slot $i$ is its momentum times the vertex with that level restored, plus a covector common to all slots. Subtract the potential at a nonparallel wave and project onto a vector annihilated by that wave to recover the coefficient. The only recursive calls have fewer waves. String and dilaton give the zero-momentum primary and level-one inputs. For a higher-level zero wave, shift it and a balancing primary in one direction. The genus-two Betti volume has two holomorphic and two antiholomorphic rows in that target direction, so this polynomial has degree at most four in the shift. Five evaluations give its value at zero. A deficient fivefold momentum rank gives zero. On a generic set the recursion uses only rational operations; polynomiality of the source vertex then gives rational coefficients everywhere. [Ward recursion and independent lowerings](genus_two_graph_ward.py)

There is a finite arity bound on independent input data. Eliminate one balancing momentum and take inclusion-exclusion over which of the remaining labeled momenta are set to zero. Every monomial of a genus-two symbol has total degree twenty, so all differences involving more than twenty independent labels vanish. Restoring the balancing slot leaves at most 21 nonzero waves. Apply string and dilaton to the zero waves of levels zero and one. If $m$ zero momenta of higher level remain among $n$ marks, then $2m\le n-7$, whereas $n-m\le21$. Thus $m\le14$ and $n\le35$. This bounds the data needed by the polynomial reconstruction; actual evaluation at a given arity can use the direct recursion.

All 81 raised entries on the nine-wave test satisfy the 36 exchanges of distinct lowering slots. Take root $e^1$, transverse vector $e_1$ at $e^2$, scalar waves $e^3,-e^1-e^2-e^3-e^4-e^5,e^4+e^5$, and four further primaries $e^4,-e^4,e^5,-e^5$. For level two at $e^4+e^5$ the graph vertex is $649/77760$; for level one at that wave and at $e^3$ it is $421/38880$. Forgetting the four additional primaries computes the same two values from threefold periods and rational-tail fusions. On two ten-wave fixtures, independent lowerings of levels two and one give $8449/38880$ and $2015618689/19440$. The rank-five source derivatives inherit these vertices by varying the action. [Independent two-axis calculation](level_two_scalar_reduction.py) · [Eight-wave first-descendant formula](level_one_scalar_reduction.py)

For a normal-ordered point polynomial $f(\phi)$ in the same coordinates, the two-loop differential on moving scalar descendants and one transverse primary vector has three coefficients multiplying $f',f'',f'''$. The first is the elementary vertex above. For $f''$ constrain the point to a bivalent theta site, or to an ordinary bivalent site on the hexagon with an elliptic insertion. For $f'''$ it is a trivalent theta endpoint; fixing it removes the exchange of the two endpoints. In each case contract the vector polarization with the potential difference from the point. On the moving level-one momentum $-e^2+e^3+e^5$, the specified second-derivative probe gives $457/2592+7/108=625/2592$, and the third-derivative probe gives $139/270$, in heat units. The full polarized coefficients on $\phi^2/2$ and $\phi^3/6$, including spectator fields, are $1745/5184$ and $20251/9720$ on their respective fixtures. [Rooted calculation, inputs and string/dilaton checks](two_loop_scalar_point.py)

### A two-loop fivefold measurement on a moving level-four scalar

Retain three labeled momenta $e^1-e^2,e^1,-2e^1+e^2$ in the flat cotangent basis. They carry, respectively, the level-four scalar, the absorbed root, and a primary vector polarized by $e_1+2e_2$. The eight remaining primaries carry the labeled scalar momenta $F=(e^2,-e^2,e^3,-e^3,e^4,-e^4,e^5,-e^5)$. Thus the vector is transverse and the absorbed momentum evaluates to one. For a subset $S$ of these labels, fuse its momenta into the distinguished scalar and lower its level to $4-|S|$. In the following sum $V$ denotes the usual labeled genus-two fivefold vertex on these fused inputs, with every other retained and unfused input left in place. Then
$$
\sum_{S\subseteq F,\ |S|\le4\,;\,\mathrm{labeled\ fusions}}(-1)^{|S|}V
=2^4\left(\frac{1+16}{1920}+\frac{4}{576}\right)
=\frac{91}{360}.
$$
In the $A_k$ convention $(2\pi i)^{10}Q_2\phi(0)$ on each such set of waves is $-V$. This is a signed combination of 163 local matrix elements. The original matrix element alone satisfies
$$
-(2\pi i)^{10}Q_2\phi(0)\ \ge\ \frac{91}{360}
+12\frac1{31104}+8\frac{29}{155520}+8\frac{53}{155520}
=\frac{10007}{38880}.
$$

For the first identity replace the distinguished cotangent class $\psi_i$ by its pullback $\bar\psi_i$ after forgetting $F$. On a fixed smooth genus-two curve each ordered pair in $F$ maps to its Abel difference in one Jacobian. This map has degree two: replacing the negative point by its hyperelliptic conjugate gives the ordinary ordered Abel sum up to translation. The hard $e^2$ components translate one target fiber. The four Jacobian volumes therefore multiply to $2^4$; the first-coordinate Betti form and $\bar\psi_i^4$ are pulled back from the retained three-mark base. The remaining one-dimensional integral has charges $(1,1,-2)$ with level four at the first mark. [One-dimensional charge formula](../twisted-holography-N4/README.md#genus-two-in-one-target-dimension-and-the-double-ramification-cycle)

To undo the pullback, $\psi_i-\bar\psi_i$ is the sum of stable rational-tail divisors containing $i$ and a nonempty subset $S$ of forgotten labels alone. In $\psi_i^k-\bar\psi_i^k$, the tail integral keeps just $\psi_i^{|S|-1}$, of integral one. The separating edge has fixed conserved current and its rational Jacobian is trivial; the Betti forms on each divisor are the fused forms on the main component. The result expresses the ordinary descendant as the sum of the pullback insertions for all $|S|\le k$, with level $k-|S|$. Subset inversion gives the alternating signs. Each forward term is nonnegative for real momenta: pair the closed positive Betti current with the nef pullback cotangent class, by an ample perturbation and continuity. At $|S|=4$ all remaining vertices are primary. The independent seven-field graph polynomial gives 42 zero terms and the 28 nonzero terms displayed in the bound. Eight genus-one tests of the same inversion, at distinguished levels one through four, agree with the complete elliptic partition calculation. [Fusion replay](ancestor_projection.py) · [Seven-field graph polynomial](longitudinal_theta.py)

#### An all-genus scalar-target measurement of the fivefold differential

Retain a list of hard momenta with axial components $p_i$ along $e^1$, total axial charge zero. A primary mark is the absorbed root and another is the transverse primary vector. At each remaining hard scalar mark choose a descendant level $k_i\le4$; the primary levels are zero. At genus $g$ require $\sum_i k_i=2g-3+n$ on the $n$ retained marks. Supply $g$ labeled primary scalar waves in each of the four covector directions $e^2,\ldots,e^5$, all with unit momentum on their axis. Choose transverse components on the hard marks to conserve momentum, keeping the vector transverse and the root evaluation equal to one.

Polarize the ordinary fivefold vertices by these labeled waves. In their signed sum each of the $4g$ additional marks either stays at its own site or fuses into one retained scalar, adding its momentum and lowering the scalar's descendant level by one. Discard negative levels and take sign $(-1)$ for each fusion. In the $A_k$ convention this signed sum of ordinary $Q_g\phi(0)$ matrix elements multiplied by $(-1)^{g+1}(2\pi i)^{5g}$ is
$$
(g!)^4[t^{2g}]\int_{\overline{\mathcal M}_{g,n}}
\mathrm{DR}_g(tp)\prod_{i\,;\,\mathrm{retained\ marks}}\psi_i^{k_i}.
$$
Here $[t^{2g}]$ extracts the indicated charge coefficient, $\mathrm{DR}_g$ is the double-ramification cycle, and $\psi_i$ is the cotangent-line class of mark $i$. Indeed the alternating sum replaces each marked cotangent class by its pullback from the retained-mark space. The relevant rational tail has just one retained mark. Distinct retained marks can therefore receive only disjoint sets of forgotten labels; on each tail the top one-mark cotangent integral equals one. Inverting the subset sum supplies the displayed sign. With pulled-back cotangent classes all fiber degree is supplied by the four maps to the Jacobian. Each ordered degree-$g$ Abel sum has degree $g!$, since a generic degree-$g$ line bundle has a unique effective divisor, whose points can be ordered in $g!$ ways. Hard transverse components translate these four maps fiberwise. The remaining scalar-target integral is the top homogeneous double-ramification term by graph deletion and contraction. [Scalar graph identity at every genus](../twisted-holography-N4/README.md#deletion-and-contraction-for-every-scalar-double-ramification-graph)

For an explicit minimal-field family let $m=\lceil(2g-1)/3\rceil$ be the number of descendant marks. Take their axial charges $1,\ldots,m$, the root charge $1$ and the vector charge $-1-m(m+1)/2$. Give the first $m-1$ scalars level four and the last level $2g+3-3m$. Choose the vector polarization $e_2$, the root transverse component $e^2$, and the first scalar transverse component $-e^2-g(e^2+e^3+e^4+e^5)$. The root evaluation is one. In order of increasing genus starting at one, the signed sums in the displayed units are
$$
\frac5{24},\qquad \frac{91}{360},\qquad
\frac{985827}{1120},\qquad \frac{19064476708}{175}.
$$
They also bound the **individual unfused elementary matrix element** from below in the same signed heat units. Undoing the cotangent pullbacks by the rational-tail formula expresses the ordinary vertex as the sum, with positive unit multiplicities, of the pulled-back vertices on all allowed fusions. Every summand is nonnegative, by the nef Betti product paired with nef pulled-back cotangent classes at real momenta.

The bound is strictly positive at every genus for the displayed family. For nonzero conserved scalar charges $p$, the top charge part with every cotangent power at mark $i$ is $[z^{2g}]\prod_{j\ne i}\sinh(p_j z/2)/(p_j z/2)>0$ by the one-mark double-ramification formula; a zero charge supplies the factor one. This implies that the scalar Betti top power pairs positively with a complementary power of an ample class: a sufficiently large ample multiple dominates $\psi_i$ by a nef class. Each other $\psi_j$ is nef and big, since the string equation and the one-point formula give $\int\psi_j^{3g-3+n}=1/(24^g g!)$. Hence $\psi_j-\epsilon H$ is pseudo-effective for a sufficiently small rational $\epsilon>0$ and a fixed ample class $H$. Its pullback has nonnegative top intersection with the remaining nef factors, also on every nef approximating model for the Betti divisor. Replacing each cotangent factor by $\epsilon H$ preserves a strictly positive lower bound. Thus the rational number on the right is positive for every allowed distribution of levels at the retained marks. In particular the genus-three and genus-four unfused elementary responses on their specified inputs are bounded below by $985827/1120$ and $19064476708/175$ in signed heat units. Their full values include the fusion terms. [One-mark coefficient, Theorem 1](https://arxiv.org/abs/1211.5273)

These are polarized, labeled combinations. The [fixed field profile](#a-single-moving-descendant-probe-of-the-elementary-differential-at-every-loop-order) gives individual responses in the same quantization.

As a direct two-loop graph test with two retained axes, take the two-dimensional hard charges $((1,0),(0,1),(-1,-1))$ with levels $(2,0,0)$ and six additional primaries at $\pm e^3,\pm e^4,\pm e^5$. The root is the second hard mark and the vector polarization is $e_2-e_1$. The graph totals with zero, one and two fusions give
$$
\frac{67}{12960}-\frac7{2160}+\frac1{2592}
=\frac1{432}=2^3\frac1{3456}.
$$
For hard charges $((1,0),(0,1),(1,-1),(-2,0))$ with levels $(1,2,0,0)$, root at the third mark and vector polarization $-e_2$, the two-moving-mark test gives
$$
\frac{11759}{12960}-\frac{1307}{1296}+\frac{4079}{12960}
-\frac2{81}=\frac{17}{90}.
$$
Both right sides come independently from the rank-two Betti intersection. On a genus-two curve the difference of an ordered pair of points has degree two by the hyperelliptic involution. [Fivefold Ward replay and scalar-target coefficients](scalar_ancestor_response.py)

#### Exact graph value and the protected source representative

In the positive $A_k$ convention the actual level-four matrix element on the preceding eleven waves is
$$
(2\pi i)^{10}Q_2\phi(0)=-\frac{19579}{12960}.
$$
In subset sizes $0,1,2,3,4$, the ordinary fused-vertex totals are respectively
$$
\left(\frac{19579}{12960},\frac{9697}{4860},\frac{311}{360},
      \frac{85}{648},\frac{179}{38880}\right).
$$
Their alternating sum is $91/360$. After subset inversion the five nonnegative contributions to the unfused vertex are $91/360,347/540,3227/6480,137/1215,179/38880$.

The Ward calculation applies to one moving scalar at each minimal level $1,2,3,4$. For the equation fixing level $h+1$ at $8+h$ waves, place its precursor of level $h$ at a bivalent theta site with $h$ other primaries, or at a trivalent branch with $h-1$ primaries. The other branch is empty. In the first case use the integrated eight-site theta with fused charge. In the second case ground the marked branch, place the eight other primaries on the paths and multiply by two for the two branches. All classical vertex factors here equal one. The genus-one vertex, when present on a hexagon, carries the descendant; its two internal momenta have opposite displacements from their common mean. Each five-row determinant in that vertex is affine in this displacement, so the position integral requires only constant and quadratic moments of the circle Green function. The genus-two cubic cut uses the previously determined descendant coefficient on one fewer wave. A cut carrying $h+1$ additional primaries with the descendant uses the seven-field primary vertex. At $h=0$ the starting formula is the longitudinal eight-wave theta minus the fused primary polynomials. Grounded differences at two nonparallel momenta determine the coefficient; degenerate configurations follow by polynomial continuation. [Single-descendant graph solver](two_loop_single_descendant.py) · [Eight-wave seed and four direct theta integrals](two_loop_moving_descendant.py) · [Elliptic insertion and two-vertex tree](two_loop_elliptic_insertions.py)

Differentiating the same action with respect to the level-four scalar gives the specified representative of the rank-five point class $[E_4]_5$. Put the primary at $e^1$ in the polarization $e_2\wedge e_3\wedge e_4\wedge e_5$ and retain the vector at $-2e^1+e^2$ in the polarization $e_1+2e_2$. Both are transverse; their ordered exterior product has unit trace. In the differentiated source the marked momentum $e^1-e^2$ is supplied by the test field. The reduced two-loop coefficient on these ten primary waves is $+19579/12960$. Source differentiation omits the one-edge bracket of the elementary scalar. The full point class is nonzero already at classical order; this evaluation fixes its value in the chosen quantum representative. [Source variation and nonexactness test](operator/README.md#closed-point-classes-from-scalar-source-insertions)

This individual source coefficient can be removed by changing its representative. Use a local primitive with one primary trivector of momentum $k$, one primary vector $v$, and seven labeled scalars. Sum over all permutations of their momenta $s_1,\ldots,s_7$ the product
$$
\det(k,s_1,s_2,s_3,s_4)^3\,T(s_2,s_3,s_4)\,k(v),
$$
where $T$ is the trivector polarization. The primitive has nineteen holomorphic derivatives and the required determinant character. On the ten primary waves above, its classical variation can supply $T$ only by fusing the rank-four polarization with a scalar of momentum $q=\pm e^j$ for $2\le j\le5$. For each $q$ the determinant enforces $s_1=-q$ and one covector on each remaining axis in $s_2,s_3,s_4$. The $3!\,2^3\,3!=288$ assignments yield classical variation $+2304$ after summing the eight channels in trivector-then-vector order. Adding the corresponding quantum-exact variation at genus two with reduced coefficient $-19579/29859840$ cancels this matrix element and preserves the lower-genus representative. [Exact improvement and independent contraction](level_four_source_improvement.py)

### A constant first descendant in the primary point differential

Use the scalar coefficient $A_1$ of the first descendant in Li–Wang (1.10), with the classical vertex $+A_1\operatorname{Tr}(\alpha_0^3)/3!$ for a constant insertion. Hold the primary fields and the normal ordering of the point observable fixed. In this convention the entire genus-$n$ interaction contribution of the primary point differential satisfies
$$
\left.Q_n f(\phi)\right|_{A_1\ {\rm constant,\ other\ fluctuations\ primary}}
=(1-A_1)^{-(4n+1)}
  \left.Q_n f(\phi)\right|_{A_1=0,\ {\rm primary}} .
$$
The equality is a formal series in the constant background, coefficientwise on arbitrary primary jets, for the flat descendant quantization. The classical scalar equation gives the positive sign: with only $\alpha_0,A_1$ the implicit equation for the scalar series is $U=\alpha_0+A_1U$. If the entire field series is literally rewritten as $A(\zeta)=\alpha(-\zeta)$, the same formula in the coefficient of $t^1$ of $\alpha(t)$ has $1+\alpha_1$ in the denominator. The sign is fixed by which variable multiplies the positive classical vertex.

At a local vertex of genus $g$ with $l=2g+3$ primary ends, give each of $m$ labeled extra marks zero momentum and one cotangent-line insertion. The Abel–Jacobi form is pulled back from the $l$ original marks since all extra momenta vanish. Its top power already has top degree on that base. On each smooth base curve, compactify the fiber of the $m$ additional points by stable rational tails. The cotangent class of the last point integrates to $2g-2+l+m-1$: the cotangent line on the moving curve acquires one twist by each of the existing marked points. For the other moving marks their cotangent classes differ from their pullbacks only on the collision divisor with the last mark. The last cotangent class restricts to zero there because the two marks lie on a rigid three-pointed rational tail. The fiber integral therefore iterates, giving
$$
V_{g;(0,\ldots,0,1,\ldots,1)}(P,0,\ldots,0)
=\left[\prod_{j=0}^{m-1}(2g-2+l+j)\right]V_{g;0}(P).
$$
The metric independence and absolute convergence of Li–Wang Lemma 1.6 allow smooth stable cotangent metrics and integration over the fibers. The zero marks have the Taylor divisor $m!$ in (1.10). Summing $m$ thus dresses the vertex by $(1-A_1)^{-(2g-2+l)}$. At genus zero (1.11) independently gives the labeled coefficient $m!$ at a cubic vertex.

For a marked graph of geometric loop number $h$, the fivefold form degree requires $5h+1$ edges and therefore $4h+1$ action vertices. At a vertex of genus $g_v$, removing all zero waves leaves exactly $2g_v+3$ ends by the dimension condition (1.9); no unstable vertex remains. The closed propagator has two primary ends, so all constant descendants are accounted for at the action vertices. The total exponent is
$$
\sum_{v\in\mathrm{action\ vertices}}(2g_v-2+2g_v+3)
=4h+1+4\sum_{v\in\mathrm{action\ vertices}}g_v=4n+1.
$$
Thus the genus expansion in $\lambda^{2n}$ on these primary probes has a common factor $(1-A_1)^{-1}$ and the substitution $\lambda^2\mapsto\lambda^2/(1-A_1)^4$. This is the dilaton equation for the specified **twisted scalar descendant**. In the flat local complex with polynomial gauge parameters, the constant first-descendant direction is exact: $(\bar\partial+t\partial_\Omega)(z_1\partial_{z_1})=t$. This persists along the classical pure constant background. A contribution to the variation of the gauge parameter with $m$ background fields would require a vertex of valence $m+2$ and descendant level $m$, whereas the allowed level there is $m-1$. The rescaling above is measured at fixed primary fields and fixed point coordinates. The polynomial gauge parameter grows at infinity; extending this local statement to a global boundary condition or an untwisted ten-dimensional dilaton requires a separate identification.

On the elementary two-loop probe with scalar momenta $e^1,e^1,e^2,e^3,e^4$ and vector $e^5$ polarized by $-e_4$, the reduced coefficients with zero, one and two **labeled** constant first descendants are
$$
-\frac{1}{31104},\qquad -\frac{9}{31104}=-\frac1{3456},\qquad
-\frac{9\cdot10}{31104}=-\frac5{1728}.
$$
Each carries the same factor $(2\pi i)^{-10}$. The one-loop primary elementary differential remains zero throughout this constant background. For a scalar square on the rank-two Hessian below the one- and two-loop coefficients instead acquire powers $(1-A_1)^{-5}$ and $(1-A_1)^{-9}$. — Codex. [Labeled insertion and two-loop probe replay](constant_descendant.py) · [Li–Wang, (1.9)–(1.11)](https://arxiv.org/abs/2610.21927)

The same calculation includes arbitrary hard descendant fluctuations. The required property of the worldsheet integral is that $\Theta_g(P)^{5g}$ pairs with ordinary cohomology classes on the compactified moduli space. To see this, push the smooth closed form on Li–Wang's real resolution down as a current. For a smooth test form pulled back from the ordinary compactification, Stokes restricts its differential pairing to boundary faces. On a face with $E$ nodes the only factors containing the $E$ node angles in (1.13) are $d_\lambda f_e\wedge d\tau_e$. The relative-length variables $\lambda$ have dimension $E-1$, so the product supplying all $E$ angles vanishes. The smooth test form has no angle component there. This proves the pushed current is closed, also when the test representative is changed by an exact form.

Let $\pi$ forget one extra zero-momentum mark $x$. On the compactification $\psi_i-\pi^*\psi_i$ is the divisor on which $i,x$ lie on a three-pointed rational tail. The restriction of $\psi_x$ to that divisor is zero. Consequently $\psi_x\prod_i\psi_i^{k_i}=\psi_x\,\pi^*(\prod_i\psi_i^{k_i})$ in ordinary cohomology for any hard descendant levels $k_i$. Inserting smooth representatives on the right into the current pairing, the Betti form on the interior also pulls back. Integrating $\psi_x$ over each moving-point fiber gives $2g-2+l$ as before. If the hard levels have total $D_v$, the initial valence is $l=2g+3+D_v$. With total hard level $D$ on the action legs of a marked graph, the genus-$n$ dressing is thus $(1-A_1)^{-(4n+1+D)}$.

Equivalently let $Q_n f$ denote the interacting contribution on a point observable involving only primary scalar jets. In the positive quartic convention of (1.10), holding the fluctuations $\alpha_j$ fixed,
$$
(Q_nf)(\alpha_0,A_1+\alpha_1,\alpha_2,\alpha_3,\alpha_4)
=\frac{1}{(1-A_1)^{4n+1}}
 (Q_nf)\left(\alpha_0,\frac{\alpha_1}{1-A_1},
 \frac{\alpha_2}{(1-A_1)^2},\frac{\alpha_3}{(1-A_1)^3},
 \frac{\alpha_4}{(1-A_1)^4}\right).
$$
The shift on the left adds a constant to the scalar component of the level-one field; all components and derivatives at hard level $j$ acquire the factor $(1-A_1)^{-j}$ on the right. At genus zero (1.11) gives the independent check $(D_v+m)!/\prod_i k_i!$ with $m$ added labeled marks. The formula includes the positive-genus vertices in Li–Wang's announced quantization. — Codex. [Descendant-level and marked-graph replay](constant_descendant.py) · [Boundary formula (1.13)](https://arxiv.org/abs/2610.21927)

The same closed-current argument evaluates constant primary insertions on the elementary scalar. Add one such mark $x$. Write $D_{i,x}$ for the collision divisor of $i$ with $x$. Distinct $D_{i,x}$ are disjoint and $D_{i,x}^2=-\pi^*\psi_i D_{i,x}$. Thus expanding each $(\pi^*\psi_i+D_{i,x})^{k_i}$ leaves the pure pullback and one term $D_{i,x}\,\pi^*(\psi_i^{k_i-1}\prod_{j\ne i}\psi_j^{k_j})$ for each positive $k_i$. The pure pullback vanishes by degree on the base. A smooth divisor representative has fiber integral one. Fubini gives the string equation on stable hard data,
$$
V_{g;(k,0)}(P,0)=\sum_{i:\,k_i>0\ \mathrm{hard\ marks}}V_{g;k-e_i}(P).
$$
Multiple constant primaries lower the hard levels repeatedly. At fixed positive genus a primary one-letter insertion receives just the one-edge tree with the quantum vertex: extra hard descendants do not change the primary internal kernel, the form degree or the bridge null vector. The absorbed leg is primary. There is consequently no genus-zero unstable correction on this observable at $n\ge1$.

With simultaneous constant scalar inputs $A_0,A_1$ at levels zero and one, substitute the following fields into the zero-background elementary $Q_n\phi$ and multiply by $(1-A_1)^{-(4n+1)}$:
$$
\alpha_j\longmapsto
\frac{1}{(1-A_1)^j}
\sum_{r=0}^{4-j}\frac{1}{r!}
\left(\frac{A_0}{1-A_1}\right)^r\alpha_{j+r}.
$$
Here the sum runs over the available higher descendant levels. First remove the constant level-one marks by the dilaton equation. On the rescaled hard fields, removing $m$ constant primary marks lowers hard level $j_i$ by $r_i$ in $m!/\prod_i r_i!$ labeled orders for $\sum_i r_i=m$, which gives the stated substitution. On hard primaries alone the positive-genus elementary response is independent of the constant $A_0$. Replacing any one hard scalar on the displayed two-loop probe by its level-one wave and adjoining one labeled constant primary gives exactly $-1/[31104(2\pi i)^{10}]$. Adjoining one more labeled constant level-one scalar gives $-10/[31104(2\pi i)^{10}]$, since the hard total descendant level before removing the primary mark is one. This also follows by applying the string equation first: nine terms in the dilaton factor come from lowering the hard descendant and one from lowering the added dilaton. — Codex. [Mixed insertions and exact primary coefficient](constant_descendant.py)

### Exact real zero locus of the primary vertex

At genus $g\ge2$ in target dimension $N$, take $n=(N-3)g+3$ conserved real momenta. In the worldsheet normalization,
$$
V_{g;0}(P)>0
\quad\Longleftrightarrow\quad
\operatorname{rank}P=N,\qquad
|T|\le g\operatorname{rank}\{p_i:i\in T\}\quad\text{for every subset of slots }T .
$$
If either test fails the integral vanishes. On the fivefold this says that the waves span five dimensions, every wave is nonzero, and at most $g$ waves lie on any one line. On the fourfold only full rank and absence of zero waves are required. Indeed a subset of rank $r<N$ leaves at least $N-r+1$ waves outside it, since their projections span the quotient and sum to zero. Thus its size is at most $2g-3+r$ on the fivefold and $g-2+r$ on the fourfold. For the elementary scalar one also multiplies by the displayed root contraction.

For a fixed smooth curve $C$, the derivatives in the mark directions have columns $p_i\otimes\operatorname{ev}_{z_i}$, where $\operatorname{ev}_{z_i}$ evaluates its $g$ holomorphic one-forms at the $i$th mark in a local coordinate. A subset violating the inequality gives dependent columns at every point. Conversely, Edmonds' matroid partition theorem gives a partition of the slots into $g$ linearly independent sets of momenta. Temporarily replace each evaluation in the $j$th set by the $j$th vector of a basis of the dual space of holomorphic one-forms. The resulting columns are independent in disjoint summands. A nonzero minor is multilinear in the individual evaluation vectors. Actual curve evaluations span the dual, since a form vanishing everywhere vanishes identically. Hence that minor is nonzero at some tuple of distinct marks, and on a dense open set. A deficient total momentum rank separately forces vanishing by Li–Wang's rank identity.

For sufficiency over the reals, pass to a finite level cover. Record the integral Abel–Jacobi differences $z_i-z_n$ for all $i<n$ in the $(n-1)$-fold power of the universal Jacobian over the Siegel moduli space. This is a mixed Shimura variety of Kuga type, of complex dimension $D=g(n-1)+g(g+1)/2$. The Zariski closure of the marked Torelli image has dimension at most $3g-3+n=gN$. Its very general point is Hodge-generic in the full Kuga variety. In fact the curve monodromy projects onto a finite-index subgroup of $\operatorname{Sp}_{2g}(\mathbb Z)$, and moving each mark $z_i$, $i<n$, around the homology of $C$ gives the full period lattice in the $i$th factor alone. The connected algebraic monodromy is therefore the full derived group $\mathbb G_a^{2g(n-1)}\rtimes\operatorname{Sp}_{2g}$. No proper special subvariety can contain this image. Proper special subvarieties form a countable collection.

In Siegel coordinates $\tau$, write $w_i$ for the lifted Abel differences. A fixed fiber of the real weighted Betti map is cut out upstairs by
$$
\sum_{\substack{i<n\\\mathrm{marks}}} p_i^{\alpha} w_i=x^{\alpha}+\tau y^{\alpha},
\qquad 1\le\alpha\le N ,
$$
where $x^{\alpha},y^{\alpha}$ are the fixed real Betti vectors. The equations are holomorphic and affine in $(w,\tau)$, of complex codimension $gN$ by full momentum rank and conservation. Suppose the maximal real Betti rank were less than $2gN$. Through a very general point of maximal rank with injective mark derivative there would be a positive-dimensional complex fiber germ. Its image in the Kuga domain varies in the base, since the mark derivative excludes any vertical tangent there and Torelli has generically finite fibers. Take an irreducible component of the intersection of the affine equations with the inverse image of the closed marked Torelli image containing this germ. The sum of the dimensions of its upstairs and downstairs Zariski closures is at most $(D-gN)+gN=D$. By mixed Ax–Schanuel its bi-algebraic closure is therefore proper.

At the chosen Hodge-generic point a bi-algebraic closure is an orbit of a connected normal subgroup of the full derived Kuga group. Its projection to the adjoint symplectic group is trivial or surjective by simplicity. In the latter case commutators with the translation representation contain the entire translation group, since the symplectic standard representation has no trivial quotient. The orbit then fills the full domain. In the former case the orbit stays over one base point. Both possibilities contradict a proper closure containing the germ above. The Betti rank is consequently $2gN$ on a nonempty open set. The squared determinant in the integrand makes the integral strictly positive. [Gao, Theorems 3.4–3.5, Definition 3.9 and proof of Proposition 5.1](https://arxiv.org/abs/1806.01408) · [Exact partition and obstruction certificates on rational input](primary_rank.py)

### Generic rank and sharp simultaneous soft limits

For $h\ge2$ conserved real hard momenta of rank $d$ on a genus-$g$ curve, $g\ge2$, the generic complex half-rank of their Betti derivative is
$$
r=\min(gd,3g-3+v),\qquad
v=\min_{T\subseteq\{1,\ldots,h\}}\bigl(h-|T|+g\,\operatorname{rank}p_T\bigr).
$$
Here $p_T$ is the list in the indicated slots. The moving-mark rank $v$ follows from the matroid union rank formula and the multilinear evaluation argument above. The integral Abel differences of all $h$ marks together with the curve have image dimension $3g-3+h$. Indeed the unweighted mark columns are independent for $h\ge2$ at general marks, since two generic canonical evaluations are nonproportional, and Torelli is generically immersive. A regular weighted fiber has complex dimension $3g-3+h-r$, of which $h-v$ dimensions are vertical over curve moduli. If $r<gd$ then the Ax–Schanuel bound above puts its bi-algebraic closure in dimension at most $D-gd+r<D$ in the ambient Kuga variety of dimension $D$. At a very general point a proper closure must be vertical. Hence $r=3g-3+v$ in this case. The two elementary upper bounds give the formula also at full target rank.

In a primary vertex on $\mathbb C^N$ let $k$ marked momenta start at zero with velocities $q_i$, and let the other $h=(N-3)g+3-k\ge2$ hard momenta have values $p_j$ at zero. The first possible power of a real conserved variation with parameter $t$ is $t^{2k}$. Its coefficient is independent of the hard recoil. It is positive exactly when both of the following hold:

1. All $h+k$ hard values and soft velocities together can be partitioned into $g$ linearly independent sets.
2. The soft velocities projected modulo the span of all the hard waves contain $g$ disjoint bases of that quotient. For a full-rank hard frame the quotient has dimension zero.

In the holomorphic derivative factor one power of $t$ from each soft mark column. The remaining determinant at $t=0$ has the hard derivative on $\mathcal M_{g,h}$ and columns $q_i\otimes\operatorname{ev}_{z_i}$. Thus arbitrary hard recoil or higher soft jets cannot enter. Its squared norm integrates coefficientwise by polynomial interpolation and Li–Wang's absolute convergence. For a subset $S$ of the soft labels let $d_S$ denote the rank of its velocities and $r_S$ the hard Betti half-rank after quotienting by their span. The hard rank must be $r_\varnothing=gN-k$. The subspaces of possible soft columns in $S$ span a space of dimension $g d_S-r_\varnothing+r_S$ modulo the hard image. Rado's independent transversal theorem requires this to be at least $|S|$ for every $S$; multilinearity again permits the actual curve evaluations. The empty subset enforces the hard rank through the same inequality with $gN-k$ in place of $r_\varnothing$. Substituting the displayed Betti rank formula gives precisely
$$
|S|+|T|\le g\,\operatorname{rank}(q_S,p_T),\qquad
k-|S|\ge g\bigl(N-\operatorname{rank}(p_{\rm hard},q_S)\bigr)
$$
for all subsets $T$ of hard labels and $S$ of soft labels. The first is matroid partition. The second says that the union of $g$ copies of the soft quotient matroid has rank $g(N-\operatorname{rank}p_{\rm hard})$, giving the disjoint bases. The same computation applies to unequal soft powers $t^{m_i}$, with $q_i$ the first nonzero coefficient of that slot and first possible power $2\sum_{i\in\mathrm{soft\ slots}}m_i$.

With full hard rank on the fivefold the coefficient is positive precisely when every combined wave is nonzero and at most $g$ lie on any one line. Indeed conservation among the hard waves alone puts at least $6-d$ hard labels outside any span of dimension $d<5$. In four dimensions the coefficient is positive for all nonzero velocities and nonzero full-rank hard configurations. For a single soft fivefold wave with the hard slots packable, its normal Hessian has rank five unless a line already contains $g$ hard waves. Then the kernel is that line and the rank is four. Two such lines would leave two hard waves whose projections to their quotient sum to zero, incompatible with total rank five. [Exact rank and quotient-base checks](worldsheet_soft_rank.py) · [Direct three-loop coefficients](soft_scalar_limit.py)

### Transverse order along a collinear fivefold configuration

Let $g+\ell$ nonzero waves lie on a single line, $\ell>0$, and let every other wave be nonzero and off that line, with conserved full rank five. Then $1\le\ell\le g-2$. Write the line momenta as $a_i u$ and denote their velocities in the quotient by the line divided by $a_i$ as their transverse slopes. The primary vertex starts at order at least $2\ell$. Its coefficient at this order is positive precisely when no $g+1$ of those slopes coincide. For $\ell=1$ the transverse Hessian is positive definite of rank $4g$; its kernel in the conserved momentum space is the tangent space of the collinear locus.

To see this, forget the line marks and project the $h=g+3-\ell$ remaining hard waves modulo the line. Their four-dimensional Betti image $T$ has rank $4g-\ell$. The mark columns are independent: a subset of rank $s<4$ contains at most $h-(5-s)$ slots by conservation, hence at most $gs$ for $s>0$. At zero deformation all $g+\ell$ line columns lie in a $g$-dimensional space. In the leading determinant exactly $\ell$ must contribute transverse velocities; all complementary columns use their constant terms. The $4g-\ell$ independent projected hard columns can therefore be removed. The residual columns are the graphs of the slope maps from the $g$-dimensional evaluation space to the $4g$ normal coordinates modulo $T$. Their line components are the canonical evaluations times $a_i$. A common slope supplies one $g$-dimensional graph. Two different slopes span the whole residual space of dimension $g+\ell$: quotient the four normal directions by their nonzero difference. The hard waves still span three dimensions, so their mark rank is at least three and the rank formula gives $3g$ in that quotient. Consequently $T$ plus the evaluation space in the difference direction fills the normal target. Rado's inequalities for the residual columns now reduce to the multiplicity bound above. When $\ell=1$, vanishing of every slope difference is exactly the linearized condition for the $g+1$ waves to stay collinear. [Independent genus-three Hessian calculation](four_collinear_hessian.py)

**Collinear dimensional reduction.** For $g\ge2$ let $g$ labeled momenta be $a_1e^1,\ldots,a_ge^1$, with the other marks remaining stable after they are forgotten. The last four sets of periods pull back from the moduli space of the remaining marks, with sources projected modulo $e^1$. They already fill its top degree. At fixed curve the first set of periods is the Abel–Jacobi map of the $g$ forgotten points with weights $a_j$ and a constant shift. Its holomorphic Jacobian has columns $a_j(\omega_\mu(z_j))_{\mu=1}^g$, where $\omega_\mu$ are normalized holomorphic differentials. Its fiber volume acquires precisely $\prod_j a_j^2$. For unit weights the Abel map from the ordered $g$-fold product of the curve has degree $g!$: a generic degree-$g$ line bundle has one section up to scale by Riemann–Roch, with $g$ distinct zeros. Deleted diagonals and coincidences with retained marks have zero fiber volume. The principal polarization on the Jacobian has normalized top volume one. Fubini for real sources and polynomial continuation give
$$
\left.V_{g;0}\right|_{\mathbb C^5}
(a_1e^1,\ldots,a_ge^1,p_{g+1},\ldots,p_{2g+3})
=g!\left(\prod_{j=1}^g a_j^2\right)
\left.V_{g;0}\right|_{\mathbb C^4}
(\bar p_{g+1},\ldots,\bar p_{2g+3}).
$$
Here $\bar p$ denotes the quotient covector in the remaining four coordinates. This formula uses the worldsheet normalization in each dimension. At $g=4$ its factor is $24\prod_j a_j^2$. Direct fivefold graph evaluations on two- and three-loop scalar probes with independently shifted longitudinal components give strength multipliers $36$ and $900$ in both period channels where present. [Graph checks](worldsheet_primary.py)

### Exact fourfold order along arbitrary real momentum curves

At genus $g\ge2$ take $g+3$ conserved fourfold waves with real formal power-series entries in $t$. Write $m$ for the least order of a four-row momentum minor and $m_i$ for the least order among such minors containing slot $i$. A zero series has infinite order. If the waves have rank below four over the Laurent-series field or one vanishes identically, the vertex vanishes identically. Otherwise its exact order is
$$
\operatorname{ord}_t V_{g;0}=2g m+2\sum_{\substack{i=1,\ldots,g+3\\\mathrm{marked\ slots}}}(m_i-m).
$$
The first nonzero coefficient is positive. Choose a matrix of four rows attaining $m$ and multiply all waves on the right by its inverse. This removes a determinant factor of order $2gm$ from the vertex. Cramer's rule makes all transformed rows integral power series, and the chosen rows are now the coordinate basis. The smallest component order of row $i$ equals $m_i-m$: every minor containing that row has at least its component order, and replacement in the coordinate basis attains it.

Divide each mark tangent column by $t^{m_i-m}$. At $t=0$ the zero-order hard waves span four dimensions and conserve. Each remaining column is the first nonzero velocity times evaluation of a holomorphic form, independently of its power of $t$. The leading determinant is thus the simultaneous soft determinant. For any proper rank-$r$ subset of its combined nonzero directions, the hard waves projected to the quotient span dimension $4-r$ and conserve there. Hence at least $5-r$ slots lie outside the subset. Its size is at most $g-2+r\le gr$ for $1\le r\le3$. The full set has $g+3\le4g$. The simultaneous soft criterion gives strict positivity. For formal series, a finite truncation preserving momentum conservation computes the same first coefficient. In the fivefold sector with $g$ collinear insertions the fiber identity above adds twice the sum of their strength orders to the order of the fourfold quotient, in a fixed unit line frame. [Exact minor algorithm, genus-two and genus-three graph interpolation, and covering checks through genus seven](fourfold_contact_order.py)

The converse use of the fiber identity determines the fourfold worldsheet vertex on arbitrary conserved momenta at genera two and three. At genus two its five inputs give
$$
\left.V_{2;0}\right|_{\mathbb C^4}(q_1,q_2,q_3,q_4,q_5)
=\frac{\det(q_1,q_2,q_3,q_4)^4}{62208}.
$$
The determinant character fixes the polynomial from one basis evaluation. At genus three the three coefficients on the matching, disjoint-edge and triangle polynomials below equal $-1/3!$ times the three rows of the fivefold sextic table. On $e^1,e^1,e^2,e^3,e^4$ and their balancing sixth wave the resulting integral on $\mathcal M_{3,6}$ is
$$
\frac{49188779}{16295634432000}
-\frac{7783}{3143448000}\zeta(3)>0 .
$$
This is an irrational primary period. The replay compares the fourfold formulas on another rational frame with the original fivefold graph evaluator in each period channel. [All three sextic rows](elementary_scalar.json) · [Fourfold evaluator and checks](worldsheet_primary.py)

At genus four, seven conserved fourfold covectors have $40$ symmetric trace symbols of determinant weight eight, $34$ after imposing single softness, and $28$ after imposing quadratic softness. Orbit sums of four squared minors span $26$ dimensions. Two signed products of eight minors complete the basis. Evaluation on 28 integer frames gives determinant $870\bmod 1009$. On three collinear waves of strengths $a,b,c$ and four waves projecting to a conserved transverse simplex, every nonzero minor is the common oriented transverse volume times one of $a,b,c,a+b+c$, up to sign. Full symmetry and quadratic softness leave two polynomials. In a unit transverse frame, two explicit basis restrictions are $24(abc)^2(a^2+b^2+c^2)$ and $4(abc)^2(a+b+c)^2$. Their determinant on strengths $(1,1,1),(1,1,-2)$ is $-20736$. [Character calculation, basis and determinant replay](operator/fourfold_quartic_basis.py)

On this frame the fivefold lift with four unit collinear waves has bracket polynomials $A=12(abc)^2(a^2+b^2+c^2)$, $B=24(abc)^2(a+b+c)^2$, $D=0$. Whenever $5(a+b+c)^2=11(a^2+b^2+c^2)$, the free trace coefficient drops out and the fourfold integral evaluates to
$$
\left.V_{4;0}\right|_{\mathbb C^4}
=-\left(\frac{11K_4}{2880}+\frac{J_4}{5}\right)(abc)^2(a^2+b^2+c^2).
$$
Single real measurements include $(a,b,c)=(1,1,(5\pm\sqrt{22})/3)$. There is no nonzero rational triple on this cone. In primitive integral coordinates reduction modulo $11$ forces $a+b+c=0$ and then $a^2+ab+b^2=0$ modulo $11$. Its discriminant $-3$ is nonsquare, forcing all coordinates divisible by $11$. [Bracket evaluator](operator/four_loop_planar.py) · [Exact period coefficients](quartic_scalar_response.py)

An eight-input fourfold Ward component measures the remaining scalar period. Use line strengths $(1,1,-2)$, five normal waves $(1,0,0),(0,1,0),(0,0,1),(-1,0,1),(0,-1,-2)$, and subtract the second hard-slot potential from the first in the third normal direction. Fusions with a line wave have zero contraction. At each hard-hard fusion the transverse volume is the determinant of the three untouched normal rows. Their signed eighth-power sum against the partner's third component is $252$. In the Li–Wang flux convention the result is
$$
252\left.V_{4;0}(1,1,-2)\right|_{\mathbb C^4}
+\int_{\rm bare}\Phi_H
+\frac{1057669573}{1123411161600000}
-\frac{143249}{181621440000}\zeta(3)=0 ,
$$
where the bare integral sums all purely cubic four-loop graphs with their symmetry factors. At the polygon with one genus-three vertex, the integrated normal factor is $3/280$ and its line factor is $22599229/252194342400-(143249/1945944000)\zeta(3)$. At the two-loop theta with a genus-two vertex, the two normal factors are $-159293/91683900000$ and $-1423/63504000$ for the branch and interior placements, respectively. Their line factor is $1/1296$. The fourfold factors also follow independently from the fivefold factors on adding four equal collinear waves and dividing by $24$.

In the bare term the three line momenta occupy three distinct reduced segments with independent cycle columns. For fixed segment lengths every other current in the missing longitudinal determinant row is affine in the three placement positions. For each unordered triple take the signed three-column minors of its cycle matrix and contract them with the endpoint potentials of the reduced normal cofactor. Write $m_i$ for the mean of the two resulting endpoint values on segment $i$, and $d_i$ for their difference. Integration on the three intervals and summation over all six line-label orders give
$$
(abc)^2(a^2+b^2+c^2)\left[\sum_{i<j\,\text{ in the triple}}(m_i-m_j)^2
 +\frac16\sum_{i\,\text{ in the triple}}d_i^2\right],\qquad a+b+c=0 .
$$
Multiply by the three lengths, sum over the $\binom{14}{3}$ triples, and divide by the square of the total length. This gives the integrated squared cofactor before the four-dimensional graph denominator and paired voltage. The replay also restores all 17 edges for a nonzero triple, checks the six terms by exact cubature, and includes the previously computed lower-vertex factors. The remaining integration is over the 14-segment projective simplex and the hard placements. [Exact Ward and bare integrand replay](operator/fourfold_quartic_ward.py)

Six labeled bare graphs have complete projective integrals. On the two-double core, in the edge order of its evaluator, the eight bivalent marks occupy paths $(0,0,0,2,5,7,7,7)$. Their label order is $(4,0,6,5,2,1,3,7)$, permuting the three line labels $0,1,2$ and holding the five normal labels fixed. With the negative-strength slot at bivalent position $1,4,5$ respectively (numbering from zero), the integrals are
$$
-\frac{1873392799}{321949590346137600000000},\qquad
-\frac{1529408143}{143088706820505600000000},\qquad
-\frac{376313593}{257559672276910080000000}.
$$
Each occurs twice. Their sum is $-80346571/2235761044070400000000$. The calculation uses all 17 original edge lengths, the squared current cofactor and a single voltage numerator over the thirteenth power of the Kirchhoff polynomial. Series and parallel beta integrals on the two arms give the displayed rationals. [Integrated graphs](operator/fourfold_quartic_ward.py)

The three line placements can instead be integrated on the 14-edge graph before enumerating the five hard labels. For each vertex $v$ assign a unit longitudinal source there with its sink at vertex zero. Multiply each of the four oriented auxiliary cofactor components of this dipole by the fourth power of the Kirchhoff polynomial $U$, obtaining a vector of polynomials $D_v$. For a triple of segments let $k$ be the vector of signed minors of its three cycle columns. On segment $i$ with endpoints $l_i,r_i$ put $m_i=k\cdot(D_{l_i}+D_{r_i})$, $d_i=k\cdot(D_{l_i}-D_{r_i})$. For strengths $(1,1,-2)$ the numerator before multiplication by the voltage numerator is
$$
2\sum_{\substack{I\\\text{three-segment subsets}}}
\left(\prod_{i\in I}t_i\right)
\left[3\sum_{\substack{i<j\\i,j\in I}}(m_i-m_j)^2+2\sum_{i\in I}d_i^2\right].
$$
The graph denominator after multiplication by the voltage numerator is $U^{13}$. The resulting density is projective in 14 lengths. The 120 hard-normal labelings give 30 distinct ordered Gale coordinate tuples, with coordinates $(-1,0,1,-1,0)$ up to an affine change. Reordering labels with identical Gale coordinates preserves the squared cofactor by an integral change of transverse frame of determinant $\pm1$. Sum their voltages, then divide by the graph automorphisms. A subgraph of loop number $h$ needs at least $3h+1$ reduced segments: the $3h$ normal cycle rows on its edges already annihilate its positive length vector.

| Cubic core | Hard path multisets | Hard shape orbits | Labeled orbits with trivial stabilizer | Labeled orbits with stabilizer two |
|---|---:|---:|---:|---:|
| Two doubled edges | 69 | 11 | 510 | 15 |
| Three doubled edges | 0 | 0 | 0 | 0 |
| One doubled edge | 171 | 36 | 2550 | 30 |
| Triangular prism | 786 | 80 | 7830 | 60 |
| Complete bipartite | 1134 | 26 | 1890 | 0 |

A nonzero value modulo $10007$ of the fully label-summed density is certified at positive integer lengths on each of the 153 admissible shapes. All eleven shape integrals on the two-double core are strictly negative after including every placement of the three line momenta, all five hard labels and graph automorphisms. Two individual values and the complete topology integral are
$$
\begin{array}{c|c}
\text{five hard paths}&\text{integrated contribution}\\ \hline
(0,0,0,7,7)&-1100123/49289264640000000\\
(0,0,0,7,8)&-51031/9241737120000000\\
\text{all eleven orbits}&-950396813/399243043584000000 .
\end{array}
$$
At hard ordering $(0,4,1,2,3)$ on the first shape only 56 of the 364 triples contribute. Their integral is $-53/70413235200000$ before the hard-label sum. Direct integration at 17 edges and the 14-edge numerator give the same value. The higher-density shapes use linear combinations of endpoint polynomials separated across the two arms before taking their square. A matching on the support of each auxiliary determinant detects structural zeros. Its strongly connected components give diagonal blocks for the nonzero determinants. [Full orbit enumeration and exact integration](operator/fourfold_quartic_orbits.py) · [All shape and label-group values](operator/fourfold_quartic_periods.json)

On the single-double core, the two complete opposite-edge orbits give
$$
\begin{array}{c|cc}
\text{five hard paths}&1&\zeta(3)\\ \hline
(0,0,0,8,8)&-10068993766933/17300531888640000000&184663/381405024000\\
(0,0,8,8,8)&-618691001569/469890989568000000&5106053/4661616960000 .
\end{array}
$$
Their coefficient determinant is $-6752092288289/886246734825330278400000000000$. They span both period channels over the rationals. Separating the bubble variables at fixed degree permits local beta integration before the remaining five paths are compressed. The residual tetrahedron powers sum to $39$ at $q=13$. [Single-double orbit integrals](operator/fourfold_quartic_single.py) · [Tetrahedron period bank](operator/quartic_seven_plane_single.py)

On the triangular prism, the hard-path orbit $(3,3,3,7,7)$ has six zero groups and 24 nonzero groups among its 30 Gale orderings. Its fully integrated contribution to the bare eight-input Ward component is
$$
-\frac{233682040128782429}{80879986579392000000000}
+\frac{85714400153}{35661369744000000}\zeta(3)<0 .
$$
The two marked paths divide into four and three ordered segments. Integrating each path split against its Dirichlet measure and eliminating the two unmarked trivalent stars reduces the numerator to tetrahedron periods. Summing through $k=64$ and bounding the remaining $\zeta(3)$ tail by $1/(2\cdot64^2)$ gives the strict sign. The other 79 prism orbits remain in the Ward sum. [Prism reduction](operator/fourfold_quartic_prism.py) · [Thirty exact group integrals](operator/fourfold_quartic_periods.json)

The same Ward component isolates a primary interaction at every genus $g\ge3$. Take $g-1$ nonzero line strengths $a_i$ summing to zero and the five three-dimensional normal waves above. After a hard fusion the other three normal waves form the transverse frame. Allow the strengths to vary independently by balancing their sum on one hard wave. Every nonzero four-row determinant is linear in the strengths. Quadratic softness in each independent line slot makes their degree-$2g$ product sum divisible by every $a_i^2$. After extracting $\prod_i a_i^2$ the remaining symmetric quadratic on the zero-sum hyperplane is a multiple of $\sum_i a_i^2$. Rescaling the transverse frame multiplies the vertex by its determinant to power $2g$. The signed fusion weights at absolute determinants $0,1,2$ are $3,-4,1$. Thus the coefficient multiplying the genus-$g$ primary interaction in this Ward equation is $4^g-4$. For $g=3,4,5$ it is $60,252,1020$.

On a bare graph the $g-1$ longitudinal jumps occupy distinct edges with independent cycle columns. Taking their signed minors leaves one longitudinal row. At fixed positions its value after extracting $\prod_i a_i$ is linear in the strengths, with one interpolated endpoint potential per mark. At $k=g-1$, summing the square over all $k!$ labelings and integrating the independent positions gives the following factor for the edge subset $I$ of size $k$. Here $m_i$ and $d_i$ are the endpoint mean and difference contracted with the vector of signed minors for $I$.
$$
(k-2)!\left(\prod_{i=1}^k a_i^2\right)\!
\left(\sum_{i\,\mathrm{line\ slots}}a_i^2\right)
\left[\sum_{\substack{i<j\\i,j\in I}}(m_i-m_j)^2
+\frac{k-1}{12}\sum_{i\in I}d_i^2\right].
$$
Indeed the permutation average of two distinct strengths is minus the mean square strength divided by $k-1$. For the full placement sum weight each subset by its edge-length product. In $\operatorname{adj}(C\operatorname{diag}(t)C^{\mathsf T})$, the coefficient of that product is the outer square of the signed minors by Cauchy–Binet. Multiplication of its first and mixed second length derivatives by $t_i$ and $t_it_j$ sums those outer squares over all subsets containing the specified one or two edges. This evaluates the complete placement sum with quadratically many edge pairs at fixed lengths for every genus. Direct subset sums and the original unsuppressed graph determinants at genera three through five reproduce the adjugate calculation. The remaining projective length integration and the lower-genus vertices enter separately. [Exact adjugate and determinant replay](operator/collinear_ward_adjugate.py)

### Reduction of the fourfold period to five and six worldsheet marks

Keep the four hard normal momenta at the unit three-dimensional simplex. Put the longitudinal recoil on its fourth mark, denoted $b$. For moduli with these four retained normal marks write $\Theta_g(\bar P)$ for their three-coordinate Betti form. For two marks $i,j$, $\Theta_g(z_i-z_j)$ is the one-coordinate Betti form with sources $+1,-1$. The fourfold primary with line strengths $a_1,\ldots,a_{g-1}$ then equals
$$
\begin{aligned}
\left(\prod_{i\,\mathrm{line\ slots}}a_i^2\right)\Bigg[&
\left(\sum_{i\,\mathrm{line\ slots}}a_i\right)^2(g-2)!
\int_{\mathcal M_{g,5}}\frac{\Theta_g(z_1-z_b)^2}{2!}
                     \frac{\Theta_g(\bar P)^{3g}}{(3g)!}\\
&+\frac{\sum_{i\,\mathrm{line\ slots}} a_i^2-
         (\sum_{i\,\mathrm{line\ slots}}a_i)^2}{2}(g-3)!
\int_{\mathcal M_{g,6}}\frac{\Theta_g(z_1-z_2)^3}{3!}
                     \frac{\Theta_g(\bar P)^{3g}}{(3g)!}\Bigg],\qquad g\ge3 .
\end{aligned}
$$
The four hard marks and the one or two marks in the displayed Abel difference are the only retained slots. After dividing by the squared line strengths the vertex is a symmetric quadratic. One nonzero strength $1$ and two nonzero strengths $1,-1$ determine its two coefficients. To evaluate either configuration take the first nonzero soft coefficient before forgetting the other marks. A soft mark contributes its curve evaluation column to the longitudinal determinant. At fixed retained moduli the longitudinal derivatives in the other columns specialize to those of the displayed Abel difference. Riemann's bilinear identity integrates the outer square of one evaluation to the identity in a period-orthonormal basis. Integrating $k$ labeled soft columns therefore gives $k!$ times the squared norm of the remaining exterior product by Cauchy–Binet. The two remaining exterior powers have degree two and three, respectively. The hard simplex contributes its pulled-back normal form throughout. Coefficient integration follows by polynomial interpolation from the absolutely convergent unscaled volumes.

At genus four $V(1,1,-2)$ equals twelve times the six-mark integral, while $W=V(1,1,1)$ equals eighteen times the five-mark integral minus three times the six-mark integral. At genus three the two integrals in rational and $\zeta(3)$ coefficient order are
$$
\left(\frac{740466719}{423686495232000},-\frac{469153}{326918592000}\right),
\qquad
\left(\frac{841479311}{211843247616000},-\frac{53359}{16345929600}\right).
$$
[Normalized soft fiber and independent genus-three primary checks](collinear_abel_moments.py)

The auxiliary marks can also be integrated from these two moments as differential forms. In de Jong–van der Lugt's notation $h$ is the curvature of the Arakelov diagonal metric, $e^A$ the curvature of the relative tangent line, $e_1^A$ the fiber integral of $(e^A)^2$, and $\nu$ the two-point fiber integral of $h^3$. Their equation (6.2) gives $\Theta_g(z_x-z_y)=h-(e_x^A+e_y^A)/2$. The normalization agrees on a fixed curve, where both forms have integral $g$ in either moving mark. Integrating in the one or two auxiliary marks gives the forms on $\mathcal M_{g,4}$
$$
\int_x\frac{\Theta_g(z_x-z_b)^2}{2!}
=\frac{e_1^A}{8}-\frac g2 e_b^A,
\qquad
\int_{x,y}\frac{\Theta_g(z_x-z_y)^3}{3!}
=\frac{g-1}{4}e_1^A+\frac{\nu}{6}.
$$
For the first formula use their calculation after (6.2). For the second, the terms in $h^2(e_x^A+e_y^A)$ and $h(e_x^A+e_y^A)^2$ cancel with their binomial coefficients by Lemma 5.3. Integration of $(e_x^A+e_y^A)^3$ gives $6(2-2g)e_1^A$. Thus at $g-1$ unit line strengths the terms in $e_1^A$ cancel at every genus. The unit-strength genus-four period is
$$
W=\int_{\mathcal M_{4,4}}\left(-36e_b^A-\frac{\nu}{2}\right)
                     \frac{\Theta_4(\bar P)^{12}}{12!}.
$$
The difference $\nu-e_1^A=(\pi i)^{-1}\partial\bar\partial\varphi$ contains the Kawazumi–Zhang invariant $\varphi$. Its integral against the normal Betti form requires boundary control on the noncompact moduli space. [de Jong–van der Lugt, equations (1.5), (6.2) and Lemma 5.3](https://arxiv.org/abs/2207.05428)

### Four-hard diagram formula for the remaining elementary-scalar period

Use $D_g$ for the threefold simplex vertex with one cotangent-line insertion and $X_g$ for its five-mark vertex with two insertions on an additional zero wave. For a degree-one function of the graph lengths, the brackets below mean its integral with the four-hard measure $24c^2U^{-3}/|\operatorname{Aut}\Gamma|$ summed over the cubic threefold graphs. Here $c$ is the projective current cofactor and $U$ the Kirchhoff determinant. At genus four the result is
$$
W=-\frac{K_4}{64}+2D_4-\frac12\left(X_4-\langle\epsilon+\varphi\rangle\right).
$$
Both new terms use the original 50 four-hard graphs. The two- and three-loop checks reproduce both coefficients of the five- and six-mark Abel moments above. [Two-descendant Ward integrals](../kodaira-spencer/zero_wave_descendant.py) · [Reduced four-loop numerator](../kodaira-spencer/four_loop_arakelov.py)

Indeed the three normal coordinates contain every Abel difference of the four hard marks. Wedging such a difference form with the normal top form gives zero. Away from collisions its diagonal curvature $h$ is $-dd^c\log G$ for the Arakelov–Green function $G$. Averaging the relation $2h=e_i^A+e_j^A$ against the top form over distinct hard pairs therefore gives
$$
\int e_b^A\,\frac{\Theta_g(\bar P)^{3g}}{(3g)!}
=\frac{K_g}{g(g!)^2}.
$$
In a degeneration the marked Green potential has the admissible graph Green coefficient. For the canonical probability measure $\mu$ that chooses a uniform position on an edge with probability $t_e\partial_e U/(gU)$, the mean graph Green function of two distinct hard slots is $(2R_{\mu\nu}-R_{\mu\mu}-R_{TT})/2$ in the earlier resistance notation. This is the integrand of $K_g/[g(g!)^2]$. [Green metric degeneration, Theorem 1.4](https://arxiv.org/abs/1511.06567) · [Direct resistance normalization](bulk_three.py)

For the other two curvature forms the pairings are
$$
\begin{aligned}
\int e_1^A\,\frac{\Theta_g(\bar P)^{3g}}{(3g)!}
   &=X_g-4D_g-\langle\epsilon\rangle,\\
\int(\nu-e_1^A)\,\frac{\Theta_g(\bar P)^{3g}}{(3g)!}
   &=-\langle\varphi\rangle .
\end{aligned}
$$
Here $\epsilon$ and $\varphi$ in graph brackets are Zhang's metric-graph invariants. The smooth pointed logarithmic $\kappa_1$ pairing is $X_g$ by forgetting the zero wave. The Deligne metric induced by a smooth metric on its cotangent line is continuous on the compact stable base. Removing the four marked sections subtracts $4D_g$ and leaves rational-tail divisors. A graph with a bridge among at least two edges has zero current cofactor. For a single separating edge with at least two hard marks on a rational tail, the other side has at most two hard marks and the node, so its three-dimensional primary form is rank deficient. All these tail terms vanish. The Arakelov Deligne pairing differs from the unmarked $\kappa_1$ by the graph height $-\epsilon$, and $\nu-e_1^A$ is minus $dd^c$ of the archimedean Kawazumi–Zhang invariant in the Green-form convention. [Continuity of Deligne metrics, Theorem A](https://arxiv.org/abs/math/9904070) · [Deligne pairing extension, Theorem 1.4](https://arxiv.org/abs/1511.06567) · [Uniform graph-height asymptotics, Theorems 1.2, 5.7 and Propositions 6.3–6.4](https://arxiv.org/abs/2306.06387)

For the boundary integral exhaust the smooth moduli space by total plumbing height $\rho=-\log\|s_D\|$, where $s_D$ is the canonical section of the boundary divisor in a smooth Hermitian metric. Test a height against $dd^c\chi(\rho/T)$ for a smooth cutoff equal to one below $T$ and zero above $2T$. The total variation against the positive normal top form on that shell is $O(1/T)$ by the chain rule and the bounded level flux. Indeed the level flux of $d^c\rho$ equals the interior integral of the smooth curvature by Stokes. Thus a uniform $o(\rho)$ metric error contributes zero in the limit. For the homogeneous slope the level measures converge to the boundary angular current: Stokes gives the convergence on smooth functions, and the uniformly bounded positive masses give it on continuous functions. The contraction with the angular form of the total boundary then weights the ordinary squared-cofactor measure by the slope. The horizontal term vanishes by $\sum_e\lambda_e\,df_e=0$. Only faces with $E=3h+1$ give the required angle count, and their positive-genus vertex factors vanish at the four-hard simplex. This uses the sign for a removed puncture in the boundary orientation. [Compactification and angular flux, Propositions 1.5 and 1.8](https://arxiv.org/abs/2610.21927)

On a trivalent core put mass $1/2$ at each branch vertex to obtain $\beta$. In terms of the constant resistance potential $\rho$ of $g\mu-\beta$ from the planar calculation and the total edge length $\ell$,
$$
\epsilon=\frac{2}{g}\bigl((g-1)\rho+R_{\beta\beta}\bigr),
\qquad
\varphi=-\frac{\ell}{4}+\frac{3\rho}{2}+\frac{\epsilon}{2}.
$$
Substitute $g_\mu(x,x)=R_{x\mu}-R_{\mu\mu}/2$ in the defining integrals
$\epsilon=\int g_\mu(x,x)((2g-2)\mu+2\beta)$ and
$4\varphi=-\ell+\int g_\mu(x,x)((10g+2)\mu-2\beta)$.
The integrated $\epsilon,\varphi$ at genus two are $1/3240,1/10368$. At genus three their coefficient pairs in $1,\zeta(3)$ order are
$$
\left(-\frac{913195967}{127105948569600},\frac{305801}{49037788800}\right),
\qquad
\left(\frac{2078212919}{317764871424000},-\frac{131371}{24518894400}\right).
$$
[Exact metric-graph integrals](../kodaira-spencer/arakelov_graph_periods.py)

### Collinear leading coefficients from fourfold soft vertices

For $g+\ell$ line marks of nonzero strengths $a_i$ in the fivefold, take $1\le\ell\le g-2$ and keep the four-dimensional quotient at full rank. Let $s_i$ be their transverse slopes, and $\bar P$ the $g+3-\ell$ projected hard momenta. For each matching $M$ of $\ell$ disjoint unordered pairs in the line slots, use their slope differences as $\ell$ soft velocities in the fourfold primary vertex. Balance the soft waves by arbitrary hard recoil there. In the worldsheet conventions,
$$
[t^{2\ell}]V^{(5)}_{g;0}
=(g-\ell)!\prod_{i\in\mathrm{line\ slots}}a_i^2
\sum_{\substack{M\\\ell\ \mathrm{disjoint\ pairs}}}
[t^{2\ell}]V^{(4)}_{g;0}\bigl(\bar P(t),\,\{t(s_i-s_j):\{i,j\}\in M\}\bigr).
$$
The operation $[t^d]$ takes the coefficient of $t^d$. The choice of orientation within each pair disappears upon squaring its evaluation column. Each term is nonnegative; it is positive exactly when all its paired slopes differ.

At a fixed hard configuration, use the metric from the period polarization and choose an orthonormal basis for the $g$ holomorphic evaluations. Riemann's bilinear relation gives identity covariance for the evaluation at each forgotten mark, integrated over the curve. Project the transverse column at slot $i$ orthogonally off the hard image; write $X_i$ for this slope map from the $g$ evaluation components to the remaining $\ell$ components. At the leading order all other columns contribute their projected hard determinant. For independent marks, integrating the remaining squared determinant equals the squarefree coefficient of the determinant of their summed covariance matrices. This follows by expanding both alternating determinants and integrating one column at a time. With scalar indeterminates $x_i$ and $S=\sum_{i\in\mathrm{line\ slots}}x_i$, the covariance sum has upper block $S\,1_g$, lower left block $\sum_i x_iX_i$, and lower right block $\sum_i x_iX_iX_i^*$. Its determinant by the Schur complement is
$$
S^{g-\ell}\det\left(\sum_{\substack{i<j\\\mathrm{line\ slots}}}
x_i x_j(X_i-X_j)(X_i-X_j)^*\right).
$$
A squarefree term selects exactly the disjoint pairs of one matching. The leftover slots occur in $S^{g-\ell}$ with coefficient $(g-\ell)!$. For a fixed matching the mixed determinant is the fourfold soft coefficient at the same hard configuration, with the same hard Gram factor. Integrating over the hard moduli proves the identity. Coincident marks and the locus of deficient hard rank have measure zero. At $\ell=0$ the same formula reproduces collinear dimensional reduction. [Exact mixed determinant checks through genus six](collinear_soft_pairs.py) · [Four-loop integrated trace check](worldsheet_fourfold_positivity.py)

### Exact fivefold order on conserved real momentum curves

On a real analytic or formal momentum curve at $t=0$, the order of the genus-$g\ge2$ primary is twice the minimum total order of $g$ five-row momentum determinants whose labels together cover all $2g+3$ slots. A zero determinant has infinite order. The minimum can be computed using one minimum-order basis and pairwise projective contacts.

Indeed choose a five-row basis of least determinant order $m$ and take coordinates in that basis. If none exists the vertex vanishes identically. By Cramer's rule all new rows are integral power series, and their constant terms span five dimensions and conserve momentum. For each slot factor $t^{s_i}$ from its moving-mark column, where $s_i$ is the minimum order of the five new components. An identically zero slot kills the vertex. Call the resulting row a unit row. The fivefold determinant character and the column factors together contribute $2gm+2\sum_{i\in\mathrm{slots}}s_i$. For each positive integer $j$, group the unit rows by equality of their projective directions modulo $t^j$. Equivalently all two-row minors of two rows in the same group are divisible by $t^j$. If $M_j$ is the size of the largest group, the exact order is
$$
2gm+2\sum_{i\in\mathrm{slots}}s_i
 +2\sum_{\substack{j\ge1\\\mathrm{jet\ orders}}}\max(M_j-g,0).
$$
Infinite total means the primary vanishes to every order. For a finite order its leading real coefficient in $V_{g;0}$ is strictly positive.

With no oversized group at $j=1$ this is the multi-soft theorem. Otherwise the group of size $g+\ell$ at level one is unique. At least five hard labels are outside its line since the hard projections to the normal four-space span and conserve. Thus $\ell\le g-2$. The deeper groups with size greater than $g$ are nested and unique. At the first level where the last oversized group splits into classes of size at most $g$, choose a reference wave in that last group. A regular shear of determinant one moves its projective direction to the coordinate line. After division by longitudinal units each initial line row takes the form $(1,s_i(t))$ with $s_i$ four transverse components. Let $b_i$ be their minimum order, possibly infinity. The sum of the $\ell$ smallest $b_i$ equals the excess sum in the display, by counting, at each level $j$, how many more than $g$ rows agree with the reference. In the holomorphic determinant at most $g$ line components can be used. At least $\ell$ clustered mark columns must supply a transverse component, giving the lower bound even with all hard and soft recoil included.

For strictness, forget the clustered marks and project modulo their line. The other normalized marks and the curve deformation directions have a joint image $T$ of dimension $4g-\ell$. Indeed append $\ell$ arbitrary nonzero soft directions to reach fourfold primary arity and use the sharp fourfold soft result. The original hard projections conserve and span four. Their further projections modulo any one nonzero normal direction span three, so the generic Betti rank there is $3g$ by the rank formula. Consequently evaluations times any one such direction fill the $\ell$-dimensional normal quotient by $T$. Fix a hard base point satisfying this for the finitely many directions used below.

Let $L$ be the first slope coefficient where the last oversized group splits into classes of size at most $g$ (so the classes first differ modulo $t^{L+1}$). All columns with $b_i<L$ must use their transverse leading terms to reach the lower bound. Choose their evaluations independently in the normal quotient and divide by their span. For each remaining mark its columns at the required order are the graph of its slope coefficient at $L$, from the $g$ line-period coordinates into the remaining normal quotient. A single graph spans $g$ dimensions and two different graphs span the residual target, since their difference is a nonzero normal direction. Rado's inequalities hold by the last splitting. The leading determinant uses exactly $g$ longitudinal columns, since every additional transverse choice has positive order. All other columns therefore use their projected constant terms, independent of the forgotten positions. Multilinearity now realizes the nonzero determinant on actual, distinct marked points. If an oversized group persists at every level its rows are exactly proportional and the determinant is identically zero. Squaring and integrating proves the real coefficient claim. For formal arcs a finite jet determines every finite coefficient; the same argument applies by polynomial truncation.

For the covering-determinant formula each basis pays the common order $m$ and the contents of its slots. Within each basis at most one clustered row can use its longitudinal component. Hence any covering pays at least the $\ell$ smallest transverse orders. To attain this, partition the $g+\ell$ line slots into $g$ groups using $\ell$ disjoint pairs. Pair the $A$ early departures ($b_i<L$) with distinct members of the last oversized group. Choose their partners to leave no final slope class larger than $g-A$. This is possible: each original class has size at most $g$, and their total size is at most $2g-A$, so the sum of their excesses over $g-A$ is at most $A$. The remaining $g+\ell-2A$ slots admit $\ell-A$ pairs with different final slopes, since at least $2(\ell-A)$ slots remain and no class exceeds $g-A$. This elementary matching bound follows by repeatedly pairing two largest distinct classes. Each early pair now pays its departure order and each later pair pays $L$. Assign four independent hard normal waves to one singleton group, and one remaining normal slot to each other singleton group. In the paired groups the leading difference is a single nonzero normal vector. Extend every group to a five-row basis using hard normal waves, all of zero content. Every slot is covered and each pair determinant attains its assigned order. With no oversized group the same equality follows from matroid partition of the normalized rows and extension by hard rows to full bases. [Minimum-basis algorithm, independent covering search and integrated three-loop replay](fivefold_contact_order.py)

The absolute leading coefficient also reduces to fourfold soft data. In the minimum-order basis, keep the contents and the central shear above, and divide each clustered unit row by its longitudinal component. For each matching of $\ell$ pairs form the sum of the orders of the transverse projective differences $s_i(t)-s_j(t)$. Keep the matchings of minimum cost only. In the fourfold soft formula insert the first nonzero transverse vector of each such difference, and one further soft velocity for each originally soft slot outside the cluster, using its normalized projected leading vector. The hard fourfold sources are the conserved constant projections of the original hard rows off the line. Sum these fourfold leading coefficients with the factor $(g-\ell)!$, the squared longitudinal strengths of the unit cluster rows, and the $2g$-th power of the leading determinant coefficient of the chosen basis. Indeed at minimum order exactly $g$ clustered columns are longitudinal, so the projected hard Gram factor is independent of every clustered position. The Schur complement calculation above now applies coefficientwise with series slopes. The matching construction attains the central lower bound, and each surviving summand is positive. At four loops, six unit longitudinal waves with transverse coordinates $(t,0,t^2,2t^2,3t^2,4t^2)$ and a seventh unit transverse planar wave give the first term $-150(W+11K_4/576+J_4)t^6$ in the period notation below. Multiplying the first wave by $t$ shifts it to degree eight with the same coefficient. [Three-channel graph replay](worldsheet_fourfold_positivity.py)

This also bounds the full momentum dependence without integrating each graph. Sum $(\det p_{B_1}\cdots\det p_{B_g})^2$ over unordered multisets of five-row bases covering every slot, and call the resulting polynomial $G(P)$. For each fixed genus there are constants $0<c_g\le C_g<\infty$ such that on every conserved real frame
$$
c_g G(P)\le V_{g;0}(P)\le C_g G(P).
$$
Both polynomials are homogeneous of degree $10g$, and their orders agree on every real analytic arc by the covering formula and positivity of the summands. If either bound failed on the unit sphere, semialgebraic curve selection on the graph of $V/G$ (or $G/V$) would give a convergent real Puiseux arc with ratio tending to zero. Reparameterizing it by an integer power gives conflicting orders. At genus two both optimal constants equal $1/622080$ by the integrated seven-field formula. Independent normalization of an individual conserved frame does not enter these bounds. [Exact genus-two coefficient](operator/README.md#first-soft-derivative-of-the-genus-three-bulk-equation)

### Normal Hessian at $g+1$ collinear fivefold momenta

Suppose exactly $g+1$ nonzero waves equal $a_i u$ on a line, the other waves project to nonzero covectors of rank four modulo $u$, and total momentum is conserved. Their projections form a conserved fourfold tuple. Let $H_4(s)$ be the $t^2$ coefficient of the fourfold genus-$g$ vertex on these projected hard waves with one additional soft wave $ts$ and any conserved recoil. Use the quotient volume induced by $u$ and the fivefold volume. For each line wave let $s_i$ be its first transverse derivative divided by $a_i$. The fivefold $t^2$ coefficient is
$$
(g-1)!\left(\prod_{i=1}^{g+1}a_i^2\right)
\sum_{\substack{1\le i<j\le g+1\\\mathrm{collinear\ slots}}}H_4(s_i-s_j).
$$
At first order in the holomorphic determinant, exactly one of the $g+1$ mark columns must vary transversely. The other $g$ span the longitudinal period space in each surviving term. All base derivatives and remaining columns then project to the fourfold quotient at order zero. Extracting $\prod_i a_i$ leaves a linear polynomial in the slopes independent of all strengths and longitudinal lifts. Its integrated squared norm is symmetric in the slopes and vanishes on common slopes, since a common variation moves the line. Polarization therefore expresses it as a sum of one quadratic form on the pair differences. Varying just one slot leaves $g$ fixed collinear marks; fiber reduction gives $g!H_4$ for that slot, fixing the factor $(g-1)!$. Fourfold single-soft positivity makes $H_4$ positive definite. The normal rank is $4g$. An analytic arc first leaving the moving line at order $m$ has exact vertex order $2m$: shear out the common direction through order $m-1$ and apply the same positive quadratic form to the first differing slopes. At genus three, interpolation of the fourfold periods reproduces both graph coefficient channels for unequal strengths, hard longitudinal lifts and independent transverse variations. [Fivefold and fourfold replay](four_collinear_hessian.py)

## The marked scalar square

At loop number $n$, the tight graph has $5n+1$ edges and $4n+1$ interaction vertices in addition to the point operator. For $\phi^2/2$, both scalar slots contract. The cubic traces require total primary rank $5(4n+1)$ and the internal lines supply $4(5n+1)$. The surviving $2n+3$ fields therefore consist of one vector and $2n+2$ scalars. The operator absorbs the negative sum of their momenta.

At one loop let $p_1,\ldots,p_5$ denote the five external momenta and let $p_v$ denote the momentum of the vector with polarization $v$. With normalized Dolbeault delta and the holomorphic volume forms to the left, the bare marked coefficient before the closed-string loop weight is

$$
-\frac{\det(p_1,\ldots,p_5)^2}{24(2\pi i)^5}(p_1+\cdots+p_5+p_v)(v).
$$

The basis covectors $p_j=e^j$, vector momentum $p_v=e^5$ and polarization $v=e_1$ give $-1/[24(2\pi i)^5]$. The six-field bulk contact below and the integrated bilinear have vanishing projections on their respective primary tests. [Marked heat facets, symmetry factor and operator redefinition](operator/README.md)

At two loops take the six **labeled** scalar momenta $e^1,e^1,e^2,e^2,e^3,e^4$, the vector momentum $e^5$, and $v=e_1+e_2+e_3+e_4$. The vector is transverse. The following coefficients include all placements and symmetry factors in the reduced graph convention.

| Theta path lengths | Nonzero labeled path partitions | Integrated coefficient |
|---|---:|---:|
| $(1,5,5)$ | 24 | $1/2430$ |
| $(2,4,5)$ | 96 | $17/17010$ |
| $(3,3,5)$ | 144 | $1/1890$ |
| $(3,4,4)$ | 288 | $1/1890$ |
| Sum | 552 | $1/405$ |

At fixed conductance fractions $w_P$ on the three paths the density times the odd contraction is $Q(L)^2\prod_P w_P^{6-n_P}(\Phi_v-\Phi_{\rm root})(v)$, with $n_P$ the lengths and $Q(L)$ the theta determinant in the repository. Summing all orders puts each labeled tail resistance fraction $s_i$ independently in $[0,1]$. The grounded Green matrix between interior vertices $i,j$ is

$$
s_i s_j+\mathbf 1_{P_i=P_j}\frac{\min(s_i,s_j)-s_i s_j}{w_{P_i}}.
$$

The two repeated covector directions each fill both cycle rows whenever the determinant survives. After removing them its square depends on four tail times. One integration expands the original ten-row determinant in monomials. A second uses the complementary six-row determinant on all 24 orderings of these four times and a simplex cubature exact through degree seven. The two integrations agree graph by graph. All 44 independently formed heat-facet determinants also agree with the theta formula at unequal resistances. [Exact integration](marked_theta.py) · [Independent cubature](theta_cubature.py) · [Literal propagator minors](point_graph.py)

In the same normalized heat convention as at one loop, the two-loop answer on this projection is $1/[405(2\pi i)^{10}]$ before its closed-string factor $\lambda^4$, with $\lambda$ the source string coupling. The closed divergence crosses five antiholomorphic differentials and differentiates the Gaussian, giving two minus signs per edge. The ten gauge-fixing contractions supply $4^{10}$. The normalized Gaussian factor is $(-32i)^{-2}(4\pi)^{-10}$. Gathering nine vertex volumes has sign $(-1)^{25\binom{9}{2}}=+1$. [Kernel sign convention](operator/README.md)

## The three-loop scalar square

At the point insert $\phi^2/2$. Use $n$ labeled scalar jets in each of the directions $e^1,e^2$, one scalar jet in each of $e^3,e^4$, and one vector of momentum $e^5$ polarized by $e_1+e_2+e_3+e_4$. Each scalar jet is that of $\exp(p\cdot z)-1$. At $n=3$ the coefficient including the string coupling is

$$
-\frac{\lambda^6}{(2\pi i)^{15}}\left(
\frac{411065203}{326918592000}
-\frac{3617}{3503500}\zeta(3)\right).
$$

Here $\zeta(3)=\sum_{m\ge1,\ m\in\mathbb Z}m^{-3}$. Each removed set of $n$ identical momentum directions must occupy a cotree, the complement of a spanning tree, by the same column subtractions as in the fan calculation below. The two cotree sums remove two Kirchhoff denominator powers and give $(n!)^2$ from labeled assignments. The remaining three-dimensional graph has four bivalent positions and trivalent ordinary vertices. At its four bivalent positions the retained covectors are the three basis covectors and their negative sum. Every permutation acts by a linear map of determinant $\pm1$, so their squared current cofactor is unchanged.

Let $C,D$ denote the vector and the point, and $R_{ij}$ the effective resistance between two positions. A unit source at $x$ contributes $(R_{CD}+R_{xD}-R_{xC})/2$ to the voltage at $C$ with $D$ grounded. Swapping $C,D$ shows that its average is half the mean pair resistance. This holds for each of the $2n$ inserted sources, and for each of the two retained scalar sources with its label held fixed. Thus the fivefold coefficient in reduced units is $(n!)^2(n+1)$ times the threefold coefficient with a scalar coefficient ring. At one and two loops this gives $2/12=1/6$ and $12(1/10080+29/272160)=1/405$.

At three loops nine catalog graphs have the required valences. For each graph let $U$ be its Kirchhoff polynomial, $c$ the three-dimensional current cofactor and $F_{ij}=U R_{ij}$ the separating two-forest polynomial. Its integrand after summing the 24 assignments is

$$
\frac{4}{|\mathrm{Aut}\,\Gamma|}
\frac{(U^3 c)^2}{U^{10}}
\sum_{\substack{i<j\\\mathrm{bivalent\ positions}}}F_{ij}.
$$

The reduced threefold integral over all nine graphs is $411065203/47076277248000-3617\zeta(3)/504504000$. Multiplying by $(3!)^2(3+1)=144$ gives the displayed fivefold coefficient. The replay also restores both sets of subdivided edges, solves the original five-dimensional currents for all 24 external assignments on each graph and checks the voltage and determinant factors at unequal resistances. [Exact nine-graph calculation](square_three.py) · [Projective period integration](../nonSP/verify/k4_periods/README.md)

For the genus-one and genus-two bulk action chosen below this is also the coefficient of the quantum differential on this fixed scalar square. All genus-one action terms in this choice have at least one external descendant. In a hexagon with the genus-two seven-primary vertex, four other bulk vertices each carry one wave. For a nonzero heat determinant those four singleton directions must be distinct. The five waves at the seven-primary vertex then occupy at most three directions, and its two incident edge momenta add at most one more independent covector by conservation. All five-row minors in its symbol $P$ vanish. Replacing any scalar wave there by its constant part cannot increase this rank. At graph loop number zero only one scalar letter of the point can contract; the remaining scalar vanishes at the point. The irrational reduced coefficient is a matrix element in this normal ordering. The local scalar square has a nonzero classical variation.

### The component invariant under redefinition of the point observable

On the same waves take the transverse vector $v=e_3+e_4$. Swapping the four retained positions in the resistance calculation gives equal mean voltage per unit source. In the four polarization directions the relative weights are $(n,n,1,1)$. At two and three loops the normalized matrix elements on $v$ are therefore

$$
\frac{\lambda^4}{1215(2\pi i)^{10}},\qquad
-\frac{\lambda^6}{(2\pi i)^{15}}
\left(\frac{411065203}{1307674368000}-\frac{3617}{14014000}\zeta(3)\right).
$$

Keep the bulk action above, the output field coordinates and the classical point observable fixed. An arbitrary finite translation- and $SL_5$-equivariant holomorphic redefinition of this marked observable leaves these two homogeneous matrix elements invariant. Indeed, in the one-vector primary output the contributing nonlinear part of a genus-$s$ change has $2+2s$ primary scalar slots and derivative degree $10s$. The exterior rank count around traced vertices leaves total polyvector rank one more than the rank of the mark. The free descendant variation pairs $p_v$ with $v$, which is zero. In the nonlinear tree variation, $v$ can differentiate only the unique $e^3$ or $e^4$ wave. Its momentum fuses with $e^5$. The resulting scalar momenta span four dimensions, and every positive-degree scalar invariant vanishes. The genus-one primary redefinition already vanishes identically since it has only four scalar slots. In a hexagon with a genus-two marked correction, the five active wave momenta must be independent. They consist of the vector and one scalar in each of the four scalar directions. All four scalar spectators at the mark lie in the $e^1,e^2$ plane. Even after including both internal currents, the six scalar arguments of the marked invariant span at most four dimensions. This exhausts the changes through genus three.

The restriction on the polarization matters. For $n=2,3$, take a symmetric local point symbol on $2+2n$ scalars given by the sum of all ordered products of $n$ squared five-row minors whose row sets together cover every scalar slot. Its classical variation on the four basis polarizations is $(4,4,0,0)$ at two loops and $(108,108,0,0)$ at three loops. Thus the coefficient in the $e_1+e_2$ direction can be changed. On this momentum family, the $e_3+e_4$ test gives the surviving component of the two-dimensional space of symmetric transverse responses. The fixed point square has a nonzero classical variation; these assertions concern redefinitions of the insertion in the fixed operator complex. [Counterterm and fusion replay](operator/square_redefinitions.py) · [All nine graphs with the four polarizations separated](square_three.py)

### The fixed scalar square detects the threefold descendant at every genus

In the primary fivefold take $2n$ scalar jets $e^{u_i\cdot z}-1$ with $u_i$ in the $e^1,e^2$ plane, scalar jets $e^{z_3}-1,e^{z_4}-1$, and vector $e^{z_5}(e_3+e_4)$. Here $e^i$ are the coordinate covectors and $e_i$ the dual vectors. Let $D_n$ be the [flat threefold one-descendant coefficient on the unit simplex](../kodaira-spencer/README.md#the-one-descendant-four-field-interaction-on-the-flat-threefold). For the normal-ordered point observable $\phi(0)^2/2$, the full polarized coefficient on these waves is
$$
Q_n\!\left(\frac{\phi(0)^2}{2}\right)
=\frac{(-1)^n}{(2\pi i)^{5n}}\,2D_n\,
  [x_1\cdots x_{2n}]\det\!\left(\sum_{\substack{i\\\mathrm{planar\ inputs}}}
             x_i u_i u_i^{\mathsf T}\right)^n .
$$
The square matrix in this determinant is formed in the indicated plane. The string expansion multiplies the coefficient by $\lambda^{2n}$. To evaluate the bare coefficient, subtract the consecutive current columns at each planar insertion. For an edge $e$ receiving slot $i$, the resulting column is $u_i\otimes C_{\cdot e}$, with $C$ the cycle matrix. Integration over positions followed by the sum over edge choices gives the mixed determinant above times the square of the Kirchhoff polynomial. Indeed it is the coefficient of all $x_i$ in the determinant of $(\sum_i x_i u_i u_i^{\mathsf T})\otimes(C\operatorname{diag}(t)C^{\mathsf T})$, by Cauchy–Binet, with $t$ the retained lengths. The polarization kills every planar contribution to the marked voltage. The remaining four positions give the three-dimensional scalar-square period, exactly $2D_n$. The coefficient of the mixed determinant is $(n!)^2$ for $n$ copies of each plane covector $e^1,e^2$. [Planar matching formula](#six-arbitrary-planar-waves)

No positive-genus primary bulk vertex can enter this matrix element. At positive geometric loop number $\ell$, a nonzero heat cofactor admits at most $2\ell$ planar waves on ordinary cubics, by the independent current-column differences along the chains (at one loop, start the chain at the mark). A genus-$h$ bulk trace admits at most $2h-1$ planar external waves: with $2h$ such waves, its other three momenta sum to zero modulo the plane and cannot supply all three transverse directions. With $m>0$ quantum vertices the resulting bound is $2n-m<2n$. Descendants cannot propagate. On a one-edge tree the uncontracted scalar at the insertion vanishes. Thus the bare expression is the full operation on the chosen fields at every genus of the fixed flat quantization.

The same component survives finite translation- and $SL_5$-equivariant holomorphic changes of the marked observable with the original classical square, bulk action and output coordinates fixed. In this one-vector projection a primary change at genus $s>0$ takes $2s+2$ scalar arguments. Its component at the measured derivative degree has $10s$ derivatives. Every positive-degree invariant on these covectors needs five-dimensional rank, so at most $2s-1$ external planar waves fit there. A derivative-free component has no nonzero external scalar spectators at the mark, since all the scalar jets vanish at the origin. Adding this bound to those at ordinary cubics and quantum bulk vertices gives at most $2n-1$ planar waves. In the exceptional one-edge classical tree the vector can differentiate only the scalar with momentum $e^3$ or $e^4$. After their fusion with $e^5$ the marked arguments have rank at most four. The free descendant change is proportional to $e^5(e_3+e_4)=0$.

Equivalently, apply one planar scalar background with quadratic polynomial $z^{\mathsf T}Hz/2$, using the symmetric Hessian matrix $H$ in the fixed two coordinates. Keep the other three external fields multilinear. Each planar slot supplies exactly two derivatives, and Taylor division by $(2n)!$ gives
$$
\left.Q_n\!\left(\frac{\phi(0)^2}{2}\right)\right|_{H,e^3,e^4,e^5}
=\frac{(-1)^n}{(2\pi i)^{5n}}\,2D_n\,(\det H)^n .
$$
The same background isolates the second derivative of *every* polynomial point observable $f(\phi(0))$. Give each of the two retained scalar jets and the vector its own parameter and take the coefficient linear in all three. A nonzero term with $k$ contracted scalar letters has $k+2n+1$ external fields. Its degree-$10n+1$ momentum symbol consists of $2n$ five-row determinants and one momentum evaluated on $e_3+e_4$. Each determinant must use exactly two planar momenta, and the last contraction annihilates the plane. This requires $4n$ planar derivatives. A quadratic field at the origin supplies exactly two, so there are precisely $2n$ planar inputs and $k=2$. All undifferentiated scalar spectators vanish at the origin. Consequently
$$
\left.Q_n f(\phi(0))\right|_{H,e^3,e^4,e^5}
=\frac{2(-1)^n D_n}{(2\pi i)^{5n}}\,(\det H)^n f''(0),
\qquad n\ge1 .
$$
For $\det H\ne0$ this component vanishes exactly when $f''(0)=0$. The square already varies classically; the calculation measures its quantum differential on cochains. [Direct two-loop polarization and Taylor expansion](operator/square_redefinitions.py)

#### A rank-one quadratic background truncates exactly at two loops

Put the quadratic scalar $H z_1^2/2$ in the first coordinate. Use the three independent scalar sources $e^{z_2}-1,e^{z_3}-1,e^{z_4}-1$ and the transverse vector $e^{z_5}(e_2+e_3+e_4)$. Extract the coefficient linear in all four of these sources. On the normal-ordered scalar observable the entire loop contribution is
$$
-\frac{\lambda^2 H}{8(2\pi i)^5}f''(0)
+\frac{\lambda^4 H^2}{20736(2\pi i)^{10}}f'(0).
$$
At one loop the hexagon gives the momentum sum on the vector polarization, equal to $3$, times $-1/24$. At two loops the seven-field trace on $e^1,e^1,e^2,e^3,e^4,e^5$ and the balancing wave is $-1/31104$ in the reduced graph convention. The root contracts to $-3$; the two identical quadratic inputs have Taylor divisor $2!$. This gives $3/(2\cdot31104)=1/20736$.

More generally put a quadratic scalar on the first $r$ coordinates, $0\le r\le3$, with one scalar source in each remaining direction through $e^4$ and the vector of momentum $e^5$ polarized in those scalar directions. Each of the $2n$ determinants uses exactly $r$ momenta from the quadratic sector, and the final vector contraction uses none. At the origin this requires $rn$ quadratic inputs. Including the $5-r$ remaining fields and equating to the primary arity $k+2n+1$ gives
$$
k=(r-2)n+4-r,\qquad 1\le k\le n+1 .
$$
Rank zero selects only $n=1,k=2$. Rank one selects $(n,k)=(1,2),(2,1)$. Rank two selects $k=2$ at every $n$. Rank three selects $k=n+1$ at every $n$. Thus changing the dimension of the quadratic source distinguishes finite loop corrections, an all-loop second-order term, and the [highest-order term](#nonvanishing-of-the-bare-marked-graph-sum-at-every-loop-order) of the same bulk differential. [Two-loop elementary and marked symbols](longitudinal_theta.py)

## Nonvanishing of the bare marked graph sum at every loop order

Put the point vertex at the origin of the odd polyvector variables. At fixed positive resistances $t_e$, take the product of the four-forms $\iota_{y_e}\prod_{i=1}^5(\theta_{e(0)}^i-\theta_{e(1)}^i)$ on all edges. Its degree is $4(5n+1)$ in $5(4n+1)$ variables. Evaluating this product against one final row defines a vector of cofactors. If the product is nonzero, its common null line is the grounded potential $\Phi$ since $\Phi_{e(0)}-\Phi_{e(1)}=t_e y_e$. Fix its scalar multiple by writing the cofactor vector as $c\Phi$. If all four-forms have full rank in their individual blocks they are decomposable. A nonzero product has the claimed one-dimensional null space. A zero current makes the odd contraction vanish.

At the external vector vertex the odd top coefficient is $c(\Phi_v-\Phi_{\rm root})(v)$. The antiholomorphic heat facet at $e$ uses precisely the same four-forms in the antiholomorphic coordinate differences, with the $e$th contracted volume restored to the full five-form. Wedge with a row dual to $y_e$ to restore this volume. The resulting coefficient is $ct_e$. Four-forms commute, so the product of the two exterior coefficients is

$$
t_e c^2(\Phi_v-\Phi_{\rm root})(v).
$$

Every remaining Gaussian and graph-symmetry weight has a common sign at fixed loop number in the indicated order of forms. Areas where rank drops contribute zero. Pick real momenta positive on $v$ for every scalar and zero on $v$ for the vector wave. The point absorbs their sum. On each connected network the grounded scalar potential is positive at the vector vertex whenever removal of the ground leaves a connected graph, by the discrete minimum principle. A graph separating at the ground has a separate resistance rescaling on each block and zero geometric density. Therefore all nonzero graphs and all their placements contribute with the same sign. Convergence also follows from bounded tree currents in the original unsquared propagator expression and the absolute-integrability argument for the scalar current Jacobian. [Signed current coordinates](../nonSP/README.md#signed-current-coordinates)

For existence begin with the hexagon. Keep the point and vector vertices bivalent. Add a parallel path of five edges to an edge whose two endpoints are bivalent scalar slots, and use an edge between two of the four new slots at the next step. The project chain move multiplies the old nonzero density by a positive resistance power times $\det(y^{\rm red},\lambda_1,\lambda_2,\lambda_3,\lambda_4)$ for the four new sources. The two old scalar endpoint momenta prescribe the two combinations with coefficients $\sigma_j$ and $1-\sigma_j$, where the four tail fractions $\sigma_j$ are distinct. If these old momenta and $y^{\rm red}$ are independent, choose the four new sources spanning a four-plane containing the old pair and excluding $y^{\rm red}$. The two prescribed combinations can be realized by a basis of this plane because their two coefficient rows are independent. The current on each edge of the new path has a nonzero coefficient of $y^{\rm red}$ modulo the new sources, so the next edge has the same independence property. In the initial hexagon generic conserved momenta have this property, as follows directly from the cumulative-current formula.

The resulting nonzero polynomial at fixed times stays nonzero on a dense set of momenta. The constraint $p_v(v)=0$ still admits a nonzero evaluation for fixed $v\ne0$ by a linear change of coordinates taking a vector in the kernel of $p_v$ to $v$. Strict positivity of the projected scalar momenta is a nonempty open condition on this hyperplane. It therefore contains rational inputs with a strictly positive graph integral at every $n$. The replay exhibits nonzero minors directly through eight loops. This proves nonvanishing of the bare marked differential on the scalar square at each loop order. The scalar square already has a nonzero classical variation, so these evaluations do not provide a protected local cohomology class. [Project chain move](../../PLAN.md#31-established) · [Exact graph witnesses](point_graph.py)

## The highest number of scalar letters differentiated

On polynomials in the undifferentiated primary scalar, the bare marked differential at loop number $n$ has order exactly $n+1$ as a differential operator in $\phi$. If $k$ edges meet the marked vertex, deleting it leaves a graph of cycle rank $n-k+1$ whenever the remaining graph is connected. If deletion disconnects it, the geometric density vanishes by the independent rescalings above. Hence $k\le n+1$.

For attainment start with the hexagon and use the same parallel five-edge path replacement, now always on an edge from the point to a bivalent scalar vertex. Each replacement increments the loop number and the root valence by one. The old scalar endpoint momentum prescribes one nonzero linear combination of the four new momenta. Since the point absorbs momentum, its reduced source then follows from conservation. Choose the four new momenta as a basis of a plane containing the old endpoint momentum but excluding the old edge current. The current on the next edge out of the point again has a nonzero component along the old edge current modulo the new four-plane. The resulting determinant is nonzero by induction. Evaluation on $\phi^{n+1}/(n+1)!$ with all $n+1$ slots contracted sets the odd coordinates at every root end to zero. The square and positivity proof applies unchanged. Scalar spectators in lower-order operations can alternatively be removed by taking the $(n+1)$-fold commutator with multiplication by $\phi$.

At two loops take seven scalar momenta $e^1,e^1,e^2,e^2,e^3,e^3,e^4$, a vector of momentum $e^5$ and polarization $e_1+e_2+e_3+e_4$. Put the point at one trivalent endpoint. For each repeated covector its two labels must occupy different paths to fill its two cycle slots. The $e^4,e^5$ pair must likewise split. The remaining four-row minor has square $s_7^2s_8^2$ for their two resistance tails. The mean of its product with the global part of the vector potential is $6/24+1/16=5/16$. Each other scalar on the vector path contributes a Brownian bridge term of mean $1/(120w_P)$. The four symmetry-weighted reduced contributions in the same shape order are
$$
\frac{43}{30240}+\frac{29}{7560}+\frac{1}{480}+\frac{1}{420}=\frac{7}{720}.
$$
Thus the term differentiating all three input scalars is nonzero independently of terms leaving a scalar spectator. Graphs containing a positive-genus bulk counterterm have lower graph loop number at the same total order, so the root-valence bound excludes them from this highest-order part, provided the marked operator and ordinary multiplication are held fixed. [Direct determinants and two-loop symbol](scalar_symbol.py)

## An exact matrix element of the scalar differential at every loop order

At the origin insert $\phi^{n+1}/(n+1)!$. Use $n$ separately labeled scalar waves of each momentum $e^1,e^2,e^3$ and one scalar wave of momentum $e^4$. For a scalar momentum $p$ take the jet of $\exp(p\cdot z)-1$. Take one vector wave of momentum $e^5$ and polarization $e_1+e_2+e_3+e_4$. Its divergence vanishes. Every scalar vanishes at the marked point, so the graph must consume all $n+1$ scalar slots. At a bivalent bulk vertex the constant part gives zero by the free series ratio. The full polarized bare result with normalized propagators, before multiplication by $\lambda^{2n}$, is

$$
\boxed{\displaystyle
\frac{(-1)^n}{(2\pi i)^{5n}}
\frac{(3n+1)n!}{2^n(n+1)(2n+1)!}},\qquad n\ge1 .
$$

The first three reduced coefficients are $1/6,7/720,1/2688$. The first differs from the single-component vector test in the marked-hexagon calculation by the sum of four equal scalar projections. The root has unit Wick derivative in this normalization.

**Removing three momentum directions.** Suppress the bivalent vertices carrying $e^1,e^2,e^3$ while retaining the total resistance of each resulting edge. In the determinant whose edge columns are $(1,C_{\cdot e}\otimes y_e)$, subtract consecutive columns along each subdivided edge. A source $e^k$ then gives $(0,C_{\cdot e}\otimes e^k)$. For each of the first three directions there are exactly $n$ such columns. Their determinant is zero unless the corresponding edges are a cotree (the complement of a spanning tree). In that case its magnitude is one by the reduced incidence identity. Eliminating those three blocks leaves precisely the two-component determinant for the sources $e^4,e^5$. The sum of all resistances is unchanged, so the absolute cofactor $|c|$ is unchanged.

The reduced graph has $n+2$ vertices and the point still has $n+1$ incident half-edges. A loop carries zero reduced current. Two parallel edges support a cycle whose two component rows are dependent, because their currents are parallel. For a nonzero determinant the reduced graph is simple. Its point joins every other vertex, and deleting it leaves a connected tree. Every remaining vertex then has tree degree at most two, making this tree a path. Its endpoints are the scalar $A$ and vector $B$. Thus the reduced graph is the single two-dimensional fan with distinguished endpoints. Conversely all three choices of cotree give the full-rank subdivision determinant whenever the fan determinant is nonzero. The marked fan with its endpoints distinguished has no further automorphisms. Each cotree admits $n!$ assignments of its labeled waves.

If an edge $e$ of total resistance $T_e$ receives $m$ subdivisions, the ordered positions have measure $T_e^m$ times the ordered unit-simplex measure. Summing the orders puts the positions independently in $[0,1]$. Write $G_{Bj}$ for the grounded voltage at vertex $j$ due to a unit source at $B$. Reciprocity identifies it with the response at $B$ from that vertex. It interpolates linearly along any edge before insertion. A unit scalar source uniformly placed on $e=(i,j)$ therefore contributes $(G_{Bi}+G_{Bj})/2$. Each cotree is weighted by $\prod_{e\in S}T_e$. Its total weight is the Kirchhoff polynomial $U$. The three independent sums remove three of the five denominator powers of $U$. The remaining measure is the two-component fan measure $c^2\,dT/U^2$ times $(n!)^3$ and the mean voltage from the wave at $A$ and the three cotrees.

**Fan integrations.** Index the path by $0,\ldots,n$ from $A$ to $B$, and ground the marked vertex. Let $R_i$ be the effective resistance of the prefix through $i$. From $i$ to $i+1$ add a path edge $S_i$ and a direct ground spoke $T_{i+1}$. Put $l_i=R_i/(R_i+S_i)$ and $b_i=T_{i+1}/(R_i+S_i+T_{i+1})$. Set $R_n=1$. Then $G_{Bj}=\prod_{i=j}^{n-1}l_i$, including $G_{Bn}=1$.

At step $i$ the surviving first-component source is $\prod_{j<i}l_j$. The second component of the reduced spoke current is $\prod_{j>i}b_j$. The two-component chain move thus gives the absolute current Jacobian in the split variables as $\prod_i l_i^{n-1-i} b_i^i(1-b_i)$. In these same variables $U=\prod_i R_{i+1}/[b_i(1-b_i)]$, and the projective change of resistances has Jacobian $\prod_i R_{i+1}^2/[b_i^3(1-b_i)^2]$. Consequently
$$
\frac{c^2}{U^2}\,dT=
\prod_{i=0}^{n-1}l_i^{2(n-1-i)} b_i^{2i+1}(1-b_i)^2\,dl_i\,db_i .
$$

For one cotree, let $m_i$ be the probability that the sampled prefix is in the two-component forest state separating $i$ from ground. Begin at $m_n=0$. The series branch at step $i$ is disconnected with probability $h_i=1-b_i(1-m_{i+1})$. Then $m_i=l_i h_i$, the path edge is excluded with probability $(1-l_i)h_i$, and the spoke at $i+1$ is excluded with probability $1-h_i+m_{i+1}$. The first spoke is excluded with probability $m_0$. Multiplying the two edge probabilities by their respective sums of endpoint voltages gives $G_{B,i+1}(1+m_{i+1})-G_{Bi}m_i$. Telescoping, including the first spoke, gives $\sum_{j=1}^n G_{Bj}$. Hence the voltage multiplying the displayed measure is $G_{B0}+\tfrac32\sum_{j=1}^n G_{Bj}$.

The mass of the product measure is
$$
\frac{2^n n!}{(2n)!(n+1)!(2n+1)!}.
$$
The normalized mean of $G_{Bj}$ is $4^{j-n}\binom{2(n-j)}{n-j}$. The sum of the first $n$ central-binomial ratios from order zero to $n-1$ is $2n\binom{2n}n/4^n$ by telescoping. Multiplication by $(n!)^3$ gives the boxed coefficient. The factor $(-1)^n/(2\pi i)^{5n}$ follows from the normalized Gaussian $(-32i)^{-n}(4\pi)^{-5n}$ and $5n$ gauge-fixing factors of four; $25\binom{4n+1}{2}$ is even. The replay checks original five-component graph cofactors, grounded voltages and explicit weighted cotree sums through six loops. The full hexagon and theta integrations above verify the first two results directly. [Fan and five-component replay](scalar_fan.py) · [Independent two-loop symbol](scalar_symbol.py)

This top scalar term is invariant under multiplicative changes of elementary holomorphic jet coordinates that reduce to the identity at $\lambda=0$. For a jet polynomial $f$, let $M_f$ denote multiplication by $f$. On any point monomial each contracted end differentiates one input letter; an incident holomorphic derivative changes its current factor and preserves the heat-form degree. Thus the graph bound of $m+1$ contracted letters at loop number $m$ bounds its differential order on the full holomorphic jet algebra. A positive-genus bulk vertex decreases the graph loop number at fixed total genus. Let $T$ be such a formal algebra automorphism. Conjugating an $(n+1)$-fold graded commutator of $T^{-1}QT$ with $M_\phi$ gives the corresponding commutator of $Q$ with $M_{T\phi}$. Every $Q_m$ with $m<n$ vanishes in that commutator by the order bound; at genus $n$ all $n+1$ multipliers reduce to $M_\phi$. Hence the top symbol of $Q_n$ is unchanged. Its evaluation above is nonzero for every $n$. At every genus to which the bulk quantization extends, ordinary multiplicative field redefinitions consequently preserve the growth of the differential order.

At arbitrary rational sources the highest scalar differential order also has rational reduced coefficients: removing the point leaves a tree, and adding one vertex to a tree gives treewidth at most two. The original graph is series–parallel. The all-integer Schwinger monomial reduction applies to its finitely many extra current factors. [Integer series–parallel periods](../nonSP/verify/k4_periods/README.md#all-holomorphic-three-loop-graphs)

### Arbitrary scalar momenta in a three-plane

For the same marked $\phi^{n+1}/(n+1)!$, let $p_1,\ldots,p_{3n}$ be the momenta of the labeled scalar waves in any common three-plane, $p_A$ the last scalar momentum, and $p_B$ the vector momentum with polarization $v$ satisfying $p_B(v)=0$. The determinant uses the fivefold volume above. The full highest-order matrix element, with the genus weight, is

$$
\frac{(-1)^n\lambda^{2n}}{(2\pi i)^{5n}}\,
\frac{(p_A+p_1+\cdots+p_{3n})(v)}
 {2^n(n+1)(2n+1)!\,n!}
\sum_{\substack{M\ \mathrm{unordered\ partitions\ of}\ \{1,\ldots,3n\}\\
                  \mathrm{into\ triples}}}
\prod_{\{i,j,k\}\in M}\det(p_i,p_j,p_k,p_A,p_B)^2 .
$$

Complete a basis of the three-plane by $p_A,p_B$ on the full-rank locus. In the plane coordinates write $S=\sum_{i\in\mathrm{wave\ labels}}x_i p_i p_i^{\mathsf T}$, with one formal variable per wave. For the reduced fan let $C_e$ be its cycle column on edge $e$ and $L=\sum_{e\in\mathrm{edges}}T_e C_e C_e^{\mathsf T}$. Subtracting adjacent columns at the removed vertices leaves $C_e\otimes p_i$ for each assigned wave. Cauchy–Binet sums their squared determinants with independent edge assignments. The coefficient of $x_1\cdots x_{3n}$ in $\det(L\otimes S)=U^3\det(S)^n$ gives the result at fixed times before the voltage. The coefficient in $\det(S)^n$ counts each unordered triple partition $n!$ times.

Holding the placement of wave $i$ on edge $e$ gives a rank-one determinant variation. Its ratio to the complete coefficient is $T_e C_e^{\mathsf T}L^{-1}C_e/n$, since $x_i\partial_{x_i}$ acts as one on the squarefree coefficient and $\partial_{x_i}\det(S)^n=n\det(S)^n p_i^{\mathsf T}S^{-1}p_i$. The fan voltage identity above then gives mean voltage $\tfrac1{2n}\sum_{j=1,\ j\in\mathrm{path}}^nG_{Bj}$ from each such unit source. Its normalized integral equals that of $G_{B0}$ from the retained scalar. The cotree mass and the integrated voltage therefore multiply the total scalar momentum evaluated on $v$. Restoring the fivefold volume gives the displayed determinants. Polynomiality extends the result to rank-deficient data.

At two loops take plane wave coordinates $(1,0,1),(2,-1,0),(0,1,1),(-1,2,1),(1,2,0),(3,0,-1)$ and retain $p_A=e^4$, $p_B=e^5$. All $5^6$ edge assignments, with the overfull assignments vanishing, give the four transverse components $(41/5,82/15,41/15,41/30)$ in reduced units. The replay checks the separate edge marginals, a restored five-dimensional graph at unequal times, and the repeated-direction specialization through six loops. [Exact calculation](scalar_fan_waves.py)

The momentum factor in this formula depends on the three-plane hypothesis. For $n=2$ take the seven scalar momenta $e^1,e^2,e^3,e^4,(1,1,1,1,1),(1,2,1,0,-1),(2,1,0,1,2)$ and $v=e_1-e_2$. Their sum annihilates $v$. For $p_B=e^5$ the four theta shapes give $-91/8640$ in total reduced units. For $p_B=(-5,-5,-3,-3,-2)$ all eight external momenta sum to zero, yet the term of differential order three on the marked scalar cubic is

$$
\frac{\lambda^4}{(2\pi i)^{10}}
\left(\frac{733}{181440}+\frac{197}{25920}+\frac{5}{1344}+\frac{121}{60480}\right)
=\frac{5\lambda^4}{288(2\pi i)^{10}}.
$$

Thus a holomorphic total derivative of the entire local expression cannot account for the highest scalar term on general jets. Each entry uses the squared fivefold theta determinant and the rooted voltage in [the path integration](scalar_fan_waves.py). Terms of lower differential order also enter the integrated observable below.

### Normal derivatives on every spoke of the marked scalar fan

In the repeated-direction probe above put $d_j$ derivatives along $e_5$ on scalar letter $j$ at the point, $0\le j\le n$. The point insertion is $\prod_j\partial_5^{d_j}\phi/(n+1)!$. Write $D$ for the sum of the derivative orders. In the paper take the two fan endpoint sources with unit dotted contraction and use the displacement component that pairs to one with the last source and zero with the first. Compare the normalized geometric kernels at fixed $n$, summing the Wick assignments of the bulk letters and uniformly permuting the corresponding boundary spokes. Their derivative moments, i.e. derivatives of the displacement exponentials at zero before Taylor factorials, are
$$
\frac{(2n+1)!}{(2n+1+D)!}\prod_{j=0}^{n}(d_j+1)!
\quad\hbox{(flat IIB)},\qquad
\frac{n!}{(n+D)!}\prod_{j=0}^{n}d_j!
\quad\hbox{(boundary fan)}.
$$
Multiply the first by the signed undifferentiated fivefold coefficient above and the second by the paper's geometric coefficient $1/[(n!)^3(n+1)]$ to restore the original weights. The first formula also multiplies the three-plane coefficient above when the differentiation direction annihilates every scalar momentum and pairs to one with the vector momentum. Its validity includes zeros of that coefficient by linearity.

For the ordered calculation write $s_j$ for the fraction of the unit last-endpoint current which enters the mark on spoke $j$. Upon attaching the last path edge, the new spoke receives $1-b_{n-1}$ and each old fraction is multiplied by $b_{n-1}$. The initial one-spoke fraction is one. The three suppressed source directions inject zero projected current, so each segment of a subdivided spoke has the same $s_j$. Holomorphic displacements of the individual marked letters multiply the integrand by the exponential of the spoke currents in the paper's incoming-displacement sign. A derivative acts on the Gaussian; the edge numerators depend only on antiholomorphic displacements. Spectator scalar jets along $e_5$ vanish with the undifferentiated spectators. Thus the derivative insertion multiplies the integrand by $\prod_j s_j^{d_j}$ without changing the contributing diagram shapes.

For $1\le i<n$, use $r_{2i-1}=l_{i-1}$ and $r_{2i}=1-b_{i-1}$ in (1-41) of the paper. In the last triangle use $r_{2n-1}=l_{n-1}(1-b_{n-1})$ and $r_{2n}=b_{n-1}$. Its absolute Jacobian supplies $1-b_{n-1}$. The boundary density on the cube is therefore $\prod_{i=0}^{n-1}l_i^{n-1-i}b_i^i(1-b_i)$. The fivefold current determinant, Kirchhoff denominator and change of variables give the square of this density times $\prod_{i=0}^{n-1}b_i$, the rooted voltage and $(n!)^3$. After cotree averaging the voltage depends only on $l_i$, also for arbitrary retained momenta in the three-plane formula. Normalize each projected integral by its own undifferentiated weight. At the step with $i$ previous spokes, $1\le i\le n$, the forward fraction $b_{i-1}$ consequently has density proportional to $b_{i-1}^{2i-1}(1-b_{i-1})^2$ in the bulk and to $b_{i-1}^{i-1}(1-b_{i-1})$ on the boundary. These fractions are independent.

For fixed ordering the ratio of the two normalized densities is
$$
\frac{(2n+1)!}{n!}\prod_{i=0}^{n-1} b_i^{i+1}(1-b_i)
=\frac{(2n+1)!}{n!}\prod_{j=0}^{n}s_j .
$$
The beta normalizations telescope. Already on the ordered projected kernels, one additional displacement derivative on **each** boundary spoke gives the bulk integral after multiplication by $(2n+1)!/n!$. For the symmetrized moments, start on the simplex with density proportional to $\prod_j s_j^{h-1}$, where $h=2$ for the bulk and $h=1$ for the boundary. Select the last spoke with probability equal to its share. Its size has its marginal density multiplied by that share, hence is beta with parameters $h+1,nh$. Conditional on this share, the old fractions rescaled to sum to one have the original simplex law on $n$ spokes by the simplex Jacobian. Iteration reproduces exactly the ordered fractions just computed. A uniform permutation restores the symmetric simplex law; its monomial integrals give the two factorial formulas.

For one twice-differentiated letter summed over its $n+1$ possible locations, the two ratios to the corresponding zero-jet amplitudes are $3/(2n+3)$ and $2/(n+2)$. Their difference in boundary-minus-bulk order is $n/[(n+2)(2n+3)]>0$ at every positive loop order. An overall scalar normalization at each genus cannot identify even these two geometric jet coefficients. The differential relation just derived holds on this projected fan; comparing the physical supercharges also requires the D3 field and operator map and the boundary color and superfield contractions. [Full current and ordered beta tests](scalar_fan_jets.py) · [Fan displacements and sign, equations (1-39)–(1-42)](../../sources/LLC_main.tex#L2281)

The factors can also be compared in the original heat conventions. Take $\mu=y e^{z^5}e_4$ and $\phi=((z^1)^2+(z^2)^2+(z^3)^2)/2+x(e^{z^4}-1)$ as in the next calculation. The vector now couples only to the retained scalar endpoint, so the voltage is $\prod_i l_i$. Cauchy–Binet on the three repeated transverse inputs gives the same result as extracting ordinary Taylor powers $n,n,n$ on the three coordinate exponentials, including the divisor $(n!)^3$. On the product of the $n+1$ scalar jets above,
$$
(2\pi i)^{5n}[xy]\,Q_n
\left(\frac{\prod_{j=0}^{n}\partial_5^{d_j}\phi(0)}{(n+1)!}\right)
=\left(-\frac12\right)^n
\frac{\prod_{j=0}^{n}(d_j+1)!}{(n!)^2(n+1)(2n+1+D)!}.
$$
The fraction on the right is exactly the symmetrized derivative of the paper's fan with $d_j+1$ derivatives on each spoke, in its unit endpoint contraction. At zero jet orders the first two signed bulk coefficients are $-1/24$ and $1/5760$. This equality holds coefficientwise at every loop order and uses the same quadratic field profile throughout.

There is also a change of the resistance variables for the full two-component current in the maximally contracted graph. Keep the differentiation vectors in the last two coordinate directions, which have zero injection at all subdividing sites. Relative to the boundary measure at $(l,b)$, the bulk density with the single endpoint voltage and the three Taylor divisors included has ratio
$$
\left(\prod_{j=0}^{n}s_j\right)\prod_{i=0}^{n-1}l_i^{n-i}.
$$
Consequently setting $L_i=l_i^2$ gives $l_i^{2(n-i)-1}dl_i=L_i^{n-i-1}dL_i/2$. Integrating any polynomial in the two retained spoke-current components on the bulk fan therefore gives $2^{-n}$ times its boundary integral with the extra factor $\prod_j s_j$ and with the currents evaluated at $(\sqrt L,b)$ in the boundary measure at $(L,b)$. For displacements along $e_5$ alone the currents are independent of $L$, which gives the unmodified displacement derivatives above. This two-component statement is for the maximal differential-order symbol. [Exact chart Jacobian, two-component current replay and both heat coefficients](scalar_fan_jets.py)

#### Exact two-endpoint kernel in displacement coordinates

Keep the fixed quadratic, the unique scalar wave in direction four, and the vector polarized along $e_4$ above. For a fan with $n$ loops, write $a_j,b_j$ for the fractions of unit current from the first and last endpoint arriving at spoke $j$, $0\le j\le n$. Each vector is in the open unit simplex. Differentiating a marked letter along $e_4,e_5$ inserts the corresponding powers of $a_j,b_j$. In this section $b_j$ denotes a *current*; the resistance split in the previous section is denoted by $\beta_i$.

On this field profile the fully contracted coefficient gives the full mixed-jet measurement $[xy]Q_n$ for $n+1$ marked scalar letters. A scalar spectator can survive only if all its derivatives are along $e_4$ and it uses the unique wave of amplitude $x$. After removing it, every remaining input covector annihilates $e_4$. The scalar symbol is linear in the sole vector polarization $e_4$. All its other vector arguments are the directions $e_4,e_5$ of marked holomorphic derivatives. An invariant scalar contracts the polarization with a covector, giving zero, or uses a determinant of five vectors, which vanishes in their two-dimensional span. Products of volume and inverse-volume tensors reduce to the same evaluations. This applies to the flat equivariant holomorphic vertices at every genus. Quadratic scalar jets follow by differentiating plane waves in the first three directions. There can be no other spectator, and consuming $n+1$ letters at this order excludes positive-genus insertions into the bare fan.

Let $A_j,B_j$ be the cumulative sums of $a_k,b_k$ through $k=j$, with $A_{-1}=B_{-1}=0$. The path current from $j$ to $j+1$ is $(1-A_j,-B_j)$. Write
$$
H_j=B_{j-1}a_j+(1-A_{j-1})b_j,\qquad
H_0=b_0,\quad H_n=a_n.
$$
The positive fan is exactly the chamber $a_j/b_j>a_{j+1}/b_{j+1}$. Indeed, up to a common positive multiple, the spoke and path resistances are $1/H_j$ and $(a_j b_{j+1}-a_{j+1}b_j)/(H_jH_{j+1})$. Substitution into the path voltage equation gives its two components. Conversely adjacent voltage drops with positive resistances yield the same ratios, and the strict inequalities follow from positivity of each $H_j$.

The boundary fan with unit endpoint contraction has Lebesgue density one in the $2n$ independent spoke currents on this chamber. To check the Jacobian, set $B_j=\prod_{i=j}^{n-1}\beta_i$ and $L_j=\prod_{i<j}l_i$. The boundary path equation gives $(1-A_j)-\beta_j(1-A_{j+1})=(1-\beta_j)L_{j+1}$. At fixed $\beta$, solving from the right is unit triangular before the factors $1-\beta_j$ and the cumulative $l$ products. The two Jacobians together give exactly $\prod_{i=0}^{n-1}l_i^{n-1-i}\beta_i^i(1-\beta_i)$. Also $H_j=b_jL_j$, including $L_0=1$. Thus the additional **unnormalized** bulk density in these coordinates, before the universal signed heat factor, is $\prod_j H_j$. At fixed spoke ordering this is a polynomial in the endpoint currents. Replacing each current by its displacement derivative therefore gives an exact finite-order differential operator on the paper's ordered kernel for arbitrary jets in both directions.

Uniformly permuting the input letters extends the boundary measure to the uniform probability on the product of two simplices after normalization. On the ordered chamber $H_i=\sum_k\min(a_i b_k,a_k b_i)$ by the ratio inequalities. All chambers consequently have the same symmetric formula. Relative to the normalized boundary measure, the normalized IIB measure is
$$
\frac{2^n(2n+1)!}{n!}\prod_{i=0}^{n}
  \left(\sum_{k=0}^{n}\min(a_i b_k,a_k b_i)\right).
$$
The constant follows from the two undifferentiated fan integrals. Across a generic interior wall $a_i b_k=a_k b_i$ the two factors $H_i,H_k$ have a nonzero jump in their first normal derivatives with the same sign; every other factor is smooth and positive. Therefore no single polynomial in the symmetrized displacement derivatives implements this unintegrated density on the whole product. Equality of every Taylor coefficient would force equality of the compactly supported measures by polynomial density.

For a direct computation, a unit current from the first endpoint moves forward to site $k$ with fraction $L_k-L_{k+1}$, where $L_{n+1}=0$, and returns to a ground spoke in the prefix through $k$ in proportions $b_i/B_k$. In particular
$$
a_i=b_i\sum_{k=i}^{n}\frac{L_k-L_{k+1}}{B_k}.
$$
For the ordered moment with prescribed powers $d_i$ of $a_i$ and $e_i$ of $b_i$, put each of the $d_i$ labeled factors at a turning site $k\ge i$. At site $j$ with $c$ first-endpoint factors already assigned, $\sum_{i\le j}d_i-c$ are eligible. Assigning $m$ gives their binomial choice and the $l_j$-moment with exponents $\sum_i d_i-c-m,m$ on $l_j,1-l_j$. At $j>0$ the $\beta_{j-1}$-moment has exponents $\sum_{i<j}(d_i+e_i)-c,d_j+e_j$. The last site assigns all remaining factors. The independent distributions of $l_i,\beta_i$ are beta with parameters $(h(n-i),1)$ and $(h(i+1),h+1)$ respectively, with $h=1$ for the boundary and $h=2$ for IIB. Storing only $c$ requires $O(n(1+\sum_i d_i)^2)$ rational transitions for an ordered jet. For indistinguishable letters the replay also stores the multiplicities of the remaining types to perform uniform permutation without enumerating all factorially many orders. [Mixed-jet evaluator, inverse resistances and grounded cubature](scalar_fan_joint.py)

**Strict positivity on real endpoint polynomials.** Write $K(a,b)=\prod_i H_i$ with the symmetric $H_i$. Each factor is a positive Gram kernel since
$$
H_i=a_i b_i\sum_k\int_0^\infty
{\bf1}_{t<a_k/a_i}{\bf1}_{t<b_k/b_i}\,dt .
$$
Multiplication gives a sum of integrals of squares in the pairing with $f(a)f(b)$ for any real polynomial $f$. In one term choose the same $k$ for every $i$. The $i=k$ integral is one. Vanishing forces all upper-box integrals of $f(a)\prod_i a_i$ in the coordinates $(a_k/a_i)_{i\ne k}$ to vanish. These coordinates are a bijection from the open simplex to the positive orthant, with smooth inverse. Mixed distributional differentiation recovers the pushed-forward density; hence $f$ vanishes on the simplex. The boundary pairing of $f$ with itself is just its uniform mean squared. Consequently every nonzero real zero-mean jet polynomial in these currents detects a difference between the two symmetrized kernels. This is a statement about the indicated Taylor measurements of the differential.

At one loop put $x=a_0,y=b_0$ in this paragraph. The unnormalized additional density is $\min(x,y)-xy$, the Dirichlet Green function of $-d^2/dx^2$ on $[0,1]$. The bulk probability density is twelve times this function. Two integrations give
$$
\mathbb E_{\rm IIB}[x^p y^q]=\frac{12}{(p+2)(q+2)(p+q+3)}
\qquad(p,q\ge0).
$$
The conditional averaging operator on the marginal density $6x(1-x)$ acts on $y^q$ by $2(1+x+\cdots+x^q)/[(q+1)(q+2)]$. It is self-adjoint and preserves each polynomial degree filtration, so the shifted Jacobi polynomial of degree $q$ for weight $x(1-x)$ has eigenvalue $2/[(q+1)(q+2)]$. At general $n$ the marginal Lebesgue density is $(2n+1)!\prod_i a_i$, and the joint density is $2^n(2n+1)!n!K(a,b)$. On the diagonal $K(a,a)=\prod_i a_i$. The Gram representation gives a Hilbert–Schmidt feature map for conditional averaging: its squared norm integrates the diagonal of the joint density divided by the marginal density with respect to ordinary simplex volume. The resulting positive operator has trace $2^n n!\cdot(1/n!)=2^n$. [Moment checks](scalar_fan_joint.py)

For the derivative on one common input letter, summed over its possible locations, the complete normalized moments are
$$
\mathbb E_{\rm IIB}\sum_{i=0}^{n}a_i b_i^m
=\prod_{k=1}^{n}\frac{4k(k+1)}{(2k+m)(2k+m+2)},\qquad
\mathbb E_{\rm paper}\sum_{i=0}^{n}a_i b_i^m
=\binom{n+m}{m}^{-1}.
$$
To evaluate the first-endpoint mixture, at site $k$ the prefix total $B_k$ is independent of the relative prefix fractions. Removing spokes in size-biased order from the symmetric Dirichlet law leaves the normalized remainder Dirichlet with the same parameter. Write $(v)_m=v(v+1)\cdots(v+m-1)$ with $(v)_0=1$. The expected sum of $(m+1)$st powers in the normalized prefix is $(h+1)_m/(h(k+1)+1)_m$. Multiplication by the independent beta moments of $B_k^m$ simplifies this to $(h+1)_m/(h(n+1))_m$ times $\prod_{i=k+1}^{n+1}i/(i+m/h)$. The expected mixture weight $L_k-L_{k+1}$ is the product of $i/(i+1/h)$ from $i=n-k+1$ to $n$, times $(1/h)/(n-k+1/h)$. Summing the products over $k$ and using the binomial convolution gives
$$
\frac{(h+1)_m}{(h(n+1))_m}\,
\frac{(n+1)!\,(1+(m+1)/h)_n}{(1+m/h)_{n+1}(1+1/h)_n}.
$$
The ratio at successive $n$ yields the two expressions above. For $m=1$ the bulk moment strictly exceeds $1/(n+1)$ for every $n>0$; the successive ratio relative to the paper is $1+1/[(2n+1)(2n+3)]$. For $m=2$ the two moments agree at every genus; for $m>2$ their order reverses. The centered differential that sums the two derivatives on the same letter and subtracts $1/(n+1)$ times their unrestricted double sum cancels on the boundary mixed projection. At two loops it equals $64/175-1/3=17/525$ times the undifferentiated bulk coefficient, giving $17/3024000$ in the signed heat convention. [Independent Kirchhoff cubature and general jet orders](scalar_fan_joint.py)

### The full differential on a quadratic scalar jet

Fix a point insertion $f(\phi(0))$ with ordinary multiplication and the marked normal ordering above. In the first three coordinates let $H_{ij}$ be the constant Hessian of a scalar quadratic form. Use the primary fields

$$
\phi(z)=\tfrac12 H_{ij}z^i z^j+x(e^{z^4}-1),\qquad
\mu(z)=y e^{z^5}e_4,\qquad 1\le i,j\le3 .
$$

Here $x,y$ count the retained scalar and transverse vector. The coefficient linear in each is

$$
[xy]\,Q_n f(\phi(0))=
\frac{(-1)^n}{(2\pi i)^{5n}}\,
\frac{(\det H)^n f^{(n+1)}(0)}{2^n(2n+1)!(n+1)!\,n!} .
$$

For $n=0$ this is $\partial(\phi\mu)$ acting on the insertion by the chain rule. At positive genus let $k$ be the number of contracted scalar letters at the point. Every spectator vanishes since $\phi(0)=0$. The primary arity is $k+2n+1$, including the vector. After the retained scalar there are $k+2n-1$ quadratic fields. Indeed, a scalar coordinate has degree two, a vector coordinate degree one, and the reduced genus-$n$ operation has degree $4n+1$. Its one-vector symbol has momentum degree $10n+1$. Volume invariance gives products of $2n$ five-row determinants and one momentum contracted against $e_4$. On these sources each determinant requires exactly three momenta from the first three directions. Each quadratic input supplies at most two derivatives at the point, so $2(k+2n-1)\ge6n$. The graph bound $k\le n+1$ forces equality. Any positive-genus bulk vertex lowers the geometric loop number and gives a strict upper bound on $k$ here. Consequently the full mixed coefficient comes from the highest-order fan at each genus where the bulk action is defined. Its bare value exists at every genus.

In the plane-wave formula the extra momentum contracted with $e_4$ gives one. Each plane wave enters precisely one squared determinant and supplies its Hessian. For a decomposition $H_{ij}=\sum_{\ell\in\mathrm{waves}} a_\ell(p_\ell)_i(p_\ell)_j$, the sum of one squared three-row minor over ordered triples with weights $a_\ell$ equals $3!\det H$ by Cauchy–Binet. The $(3n)!/((3!)^n n!)$ unordered triple partitions and the Taylor divisor $(3n)!$ on the repeated plane scalar give $(\det H)^n/n!$. This argument extends by polynomiality to arbitrary complex Hessians.

For the single formal observable $e^{\phi(0)}$ the mixed perturbative series sums to the entire function

$$
\sum_{\substack{n\ge0\\\rm loop\ orders}}\lambda^{2n}[xy]\,Q_n e^{\phi(0)}
={}_0F_3\!\left(;1,\tfrac32,2;-\frac{\lambda^2\det H}{8(2\pi i)^5}\right).
$$

The generalized hypergeometric series on the right has successive coefficient ratio $1/[(n+1)^2(n+3/2)(n+2)]$ in its displayed argument and constant term one. For $H$ the identity, the reduced magnitudes at genera one, two and three are $1/24,1/5760,1/5806080$. Independent polarization of an indefinite Hessian with determinant $-4$ checks the coefficients through genus four. The two lower-order genus-two exponential-wave symbols scale with the twelfth power of the plane momenta, beyond the eight or ten derivatives available on the corresponding quadratic jets. [Hessian and graph replay](scalar_quadratic_jet.py)

## The integrated cubic and its classical completion

In full BCOV the constant holomorphic top polyvector $\eta$ with $\Omega(\eta)=1$ gives a closed nonpropagating background $t\eta$, where $t$ here denotes the descendant parameter. Differentiate the classical interaction with respect to that background and restrict to the minimal fields. The free differential preserves the minimal subspace and every Poisson kernel end belongs to it. Differentiation of the classical master equation therefore proves closure. On primary fields the resulting integrated functional is $\int\Omega\,\phi^3/3!$, with the top Dolbeault component understood. Its nonzero derivative-free cubic cannot be a classical boundary in translation-invariant local functionals, since the free differential, interacting differential and total divergences all take a derivative. [Classical descendant completion](operator/README.md) · [Coefficient replay](operator/integrated_descendants.py)

Use the eight labeled holomorphic momenta $e^1,e^1,e^2,e^2,e^3,e^4,e^5,-2e^1-2e^2-e^3-e^4-e^5$. The vector has the last momentum. Its polarization $e_1-2e_3$ is transverse. At a trivalent marking all three fields of the integrated cubic contract. At a bivalent marking one external scalar is a spectator. Descendant letters cannot contract. On this all-primary projection the bare two-loop reduced coefficients are

| Theta | Trivalent marks | Bivalent marks | Sum |
|---|---:|---:|---:|
| $(1,5,5)$ | $41/816480$ | $23/204120$ | $19/116640$ |
| $(2,4,5)$ | $11/102060$ | $389/1632960$ | $113/326592$ |
| $(3,3,5)$ | $1/18144$ | $11/90720$ | $1/5670$ |
| $(3,4,4)$ | $1/22680$ | $17/181440$ | $5/36288$ |
| All | $1/3888$ | $11/19440$ | $1/1215$ |

The external momenta sum to zero, so this term survives integration over the marked point. In the normalized heat convention its coefficient is $1/[1215(2\pi i)^{10}]$ before $\lambda^4$. The sign and all graph automorphisms are as for the point operator. With the far endpoint grounded, summing over markings at vertex $i$ gives the potential difference $(\Phi_v-\Phi_i)(v)$ on each one. The mark at the vector contributes zero. Independent ordered-region cubature reproduces all eight entries in the first two columns. [Full marking sum](integrated_cubic.py) · [Independent integrals](integrated_cubature.py)

**Two probes for redefinitions of the marked observable.** Replace the fourth momentum by $b e^2$, balance the last momentum, and choose vector polarization $e_2-(1+b)e_5$. At $b=1,2$ the two integrated coefficients are $1/1215$ and $2/405$. At order $\lambda^4$, a holomorphic, translation- and $SL_5$-invariant local primitive contributing through the primary classical differential must have seven scalar inputs and twenty derivatives. For seven conserved covectors choose the first six as independent *slots*. Their six five-row determinants are algebraically independent: on the dense rank-five locus they are the signed coordinates of the second linear relation, modulo the constant relation among all seven. More explicitly the $SL_5$ invariants on a five-by-six matrix are generated by its maximal minors, by contraction of the covector slots against copies of the volume tensor. Every five-vector in six dimensions is decomposable. Hence the six minors take all six-tuples of values. Center the seven coefficients of the second relation by subtracting their mean. Permuting the seven scalar slots permutes the centered entries up to the sign of the permutation. A homogeneous quartic is insensitive to the common sign. Symmetric quartics on this sum-zero space are spanned by the square of the sum of squares and the sum of fourth powers.

For the integer evaluation in the replay each centered coordinate is multiplied by seven. To take the classical differential, fuse the vector with each scalar in turn and multiply by the pairing of its polarization with that scalar momentum. The vector itself is transverse. The two boundary columns and the bare marked graph give

| $b$ | Square of sum of squares | Sum of fourth powers | Graph |
|---:|---:|---:|---:|
| 1 | 19208 | 9604 | $1/1215$ |
| 2 | 172872 | 86436 | $2/405$ |
| Second minus nine times first | 0 | 0 | $-1/405$ |

The free differential on a scalar descendant gives the divergence of the vector and vanishes here. An order-$\lambda^2$ primary marked redefinition acted on by the one-loop graph would need five scalar fields and ten derivatives. Five conserved covectors span at most four dimensions, so every positive-degree $SL_5$-invariant polynomial of these covectors vanishes. Descendants in that redefinition would be external. This rules out cancellation of the two-row residual by changing the marked operator alone in the fixed bare holomorphic jet model. The bulk term below contributes to this differential. [Exact two-source comparison](marked_boundaries.py)

## The genus-two bulk action and integrated cubic

Let $q_1,\ldots,q_7$ be seven conserved primary momenta. In the integer centering used above, $c_i=(-1)^{i-1}\det(q_1,\ldots,\widehat q_i,\ldots,q_6)$ for $1\le i\le6$, $c_7=0$, and $a_i=7c_i-\sum_{j\in\mathrm{slots}}c_j$. The hat removes that row. In the trace of seven labeled primary polyvectors the momentum polynomial is

$$
P=-\frac{31}{2987228160}\left(\sum_{i\in\mathrm{slots}}a_i^2\right)^2
  +\frac{1}{71124480}\sum_{i\in\mathrm{slots}}a_i^4 .
$$

The polarized genus-two action term is $-P\,\mathrm{Tr}$, where the trace selects top polyvector and Dolbeault degree. For the completed integrated scalar cubic $F_1$ the polarized seven-scalar correction is $-9P$. Both have the common heat multiplier $\lambda^4/(2\pi i)^{10}$. These are local symbols in the project's normal ordering.

On seven conserved real momenta, $P<0$ precisely at rank five with every momentum nonzero. In the disjoint-pair formula each summand in $-P$ is positive when both complementary cofactor differences are nonzero. A complete multipartite graph on the seven cofactor coordinates has two disjoint edges between unequal values unless six values coincide or all agree. At full rank six equal values mean the remaining momentum vanishes, since the difference of the two independent relations then selects that single momentum. [Disjoint-pair formula](operator/README.md#first-soft-derivative-of-the-genus-three-bulk-equation)

For the soft limit, take six conserved hard momenta $u_i$ spanning five dimensions, with $D=\det(u_1,\ldots,u_5)$. Add a soft momentum $t v$ and shift every hard momentum by $-tv/6$. Write $v=\sum_{i\in\mathrm{hard\ slots}}b_i u_i$ with $\sum_{i\in\mathrm{hard\ slots}}b_i=0$ and let $p_k=\sum_{i\in\mathrm{hard\ slots}}b_i^k$. The exact trace along the entire line is
$$
P(t)=-\frac{D^4t^2}{1244160}
\left(48p_2-24t p_3+t^2(31p_2^2-42p_4)\right).
$$
Indeed, the five-row minor omitting the soft wave and hard slot $i$ is $\pm D(1+t b_i)$ by the rank-one determinant identity. The minor omitting hard slots $i,j$ is $\pm tD(b_i-b_j)$. Inserting these in the disjoint-pair sum gives the three displayed coefficients. In particular the coefficient of $t^2$ is strictly negative for every nonzero real $v$. For $v\ne0$ the discriminant gives $3p_3^2\le p_2(31p_2^2-42p_4)$ with equality precisely when $v$ lies on one of the six hard lines. In fact the shifted waves and $tv$ still span the hard frame, and an additional real zero requires $tv=6u_i$. The quartic coefficient is strictly negative for $v\ne0$ by the seven-coordinate difference formula applied to $(b_1,\ldots,b_6,0)$. [Exact rational recoil replay](seven_scalar_recoil.py)

For eight conserved primary covectors $p_i$ let $P_{ij}$ be $P$ on the seven slots obtained by fusing $p_i,p_j$. Let $\Phi_i$ be the integrated potential at slot $i$, grounded at the second trivalent tip, including all four theta shapes and their symmetries. Put $C_i=\sum_{j\ne i,\ j\in\mathrm{slots}}p_jP_{ij}$ and $E_i=\Phi_i-C_i$. The exact polynomial identity is

$$
p_i\wedge\left(10 E_i-\sum_{s\in\mathrm{slots}}E_s\right)=0
\qquad (1\le i\le8).
$$

For the bulk trace this gives a common translation and a multiple of each leg's own momentum. The common translation drops out of the trace and the latter terms vanish on divergence-free external legs. The bare defect has the negative potential-contraction sign. The classical variation of $-P\,\mathrm{Tr}$ supplies the compensating $C_i$ term on these legs. An independent cubature of a transverse potential difference gives

$$
\frac1{81648}+\frac{37}{1632960}+\frac1{90720}+\frac1{181440}=\frac1{19440}.
$$

The four summands use the shape order of the preceding table. The two trace polynomials have independent boundary images. Three of the four individual shape tensors lie outside their span. The symmetry-weighted sum fixes both displayed coefficients. At this arity a genus-one primary vertex in a wheel would have five slots, and two genus-one tree vertices would have ten slots between them. A primary genus-one term has five more holomorphic derivative covectors than vector slots by scale weight. At least one volume contraction therefore takes five momenta. Five or fewer conserved momenta cannot give a nonzero such determinant. Descendants on a quantum vertex are external since the propagator and bracket have primary ends. Consequently an eight-leg transverse response at this order requires a genus-two action correction in this holomorphic model. [Potential calculation](bulk_eight.py) · [Independent transverse cubature and its source](bulk_cubature.py)

The individual theta sectors admit an exact classification on arbitrary transverse primary waves. Let $S_1,S_2,S_3,S_4$ be their integrated potentials modulo common and longitudinal shifts, in path-length order $(1,5,5),(2,4,5),(3,3,5),(3,4,4)$. On the fourteen unisolvent transverse measurements the first two columns are independent and
$$
3S_3=-S_1+2S_2,\qquad 9S_4=-20S_1+13S_2,\qquad
S_3=\frac{3}{14}(S_1+S_2+S_3+S_4).
$$
There is one independent obstruction to reweighting the four shapes. For weights $c_i$, a seven-slot trace cancels their transverse source precisely when $2c_1+c_2-3c_4=0$. Its coefficient relative to the full trace is $(-9c_1+3c_3+20c_4)/14$. Necessity can be measured directly: on an eight-wave fixture with three components $1,2,-3$ outside a four-plane, take one finite momentum-transfer difference for each of the six pairs joining two inside slots to one outside slot. Every trace boundary has zero sixfold difference. Direct integrations on the 64 configurations give $102820139(2,1,0,-3)/136080$ for the four shapes. The two global column identities give sufficiency. [Integrated shape and momentum-transfer replay](operator/fusion_transfer_box.py) · [Fourteen unisolvent transverse measurements](operator/genus_two_probes.py)

In fact both shape identities hold for each raw slot covector when the second trivalent endpoint absorbs arbitrary net momentum. Before any transverse quotient the eight external momenta are independent. Each squared current determinant contains at least two powers of every bivalent wave by its current jump, and the voltage supplies the last covector. The resulting degree-21 array has determinant weight four. In the character computation, take one distinguished slot and the seventh symmetric power of the independent-slot series of monomials of degree at least two. The multiplicity of $(5,4,4,4,4)$ is exactly 38. Symmetrized products of four five-row minors times one covector give a 38-row evaluation minor $397$ modulo $1009$ at three unconstrained momentum configurations. Direct heat integration gives both shape relations in every raw component there. Thus the $(3,3,5)$ sector contributes $3/14$ of the full grounded potential also for arbitrary momentum entering that endpoint. [Unconstrained character, certificate and graph integration](operator/grounded_theta_shapes.py)

For the cubic, the marked graph covector at the vector leg $i$ is $10\Phi_i-\sum_{s\in\mathrm{slots}}\Phi_s$. The bulk action contributes $\sum_{k<l,\ k,l\ne i,\ k,l\in\mathrm{slots}}(p_k+p_l)P_{kl}$. Adding the classical variation of $-9P$ leaves a multiple of $p_i$ by the displayed identity. On the preceding two transverse sources the cancellations in reduced units are

$$
\frac1{1215}-\frac{109}{77760}+\frac1{1728}=0,\qquad
\frac2{405}-\frac{263}{25920}+\frac1{192}=0.
$$

The separating polynomial tests have ranks fourteen for the bulk and twenty for the joint identity. Their determinants are $96$ and $927$ modulo $1009$. Each rational identity on these tests thus extends to all conserved momenta, including degeneracies. On general fields the remaining longitudinal factors have polynomial coefficients since the five components of each $p_i$ are relatively prime. They are free variations of terms with one descendant. At arity $N>8$ the remaining bulk and cubic defects have descendant level $N-8>0$. The BV consistency identity and exactness of the free complex at positive level give formal local lifts of the bulk action and $F_1$ through genus two. This calculation uses the explicit classical descendant action in the translation-invariant, $SL_5$-equivariant holomorphic jet model. The identification with a ten-dimensional untwisted interaction remains open. [Polynomial test certificates](operator/genus_two_probes.py) · [Free jet complex](operator/README.md#eliminating-higher-descendant-equations-in-the-holomorphic-jet-complex)

For the integrated scalar quartic the seven-field bulk correction enters together with the one-loop action on its six-scalar genus-one correction. A rank-51 polynomial certificate gives a genus-two lift of this observable as well. [Integrated quartic and its four coefficients](operator/README.md#the-integrated-quartic-at-genus-two)

### Closed two-loop potential on all conserved primary waves

The eight covector potentials $\Phi_i$ with the second trivalent endpoint grounded admit an expression in five-row determinants. With $C_i$ as above, the exact result is

$$
\Phi_i=C_i+p_i f_i+\tfrac12\sum_{j\in\mathrm{slots}}p_j f_j .
$$

For $f_i$ delete slot $i$ and order the remaining seven labels. For local labels $j<k$ starting at zero let $[jk]$ be $(-1)^{j+k-1}$ times the determinant omitting rows $j,k$, and $[kj]=-[jk]$. On two- and three-element subsets use $W_{jk}=[jk]^2$ and $W_{jkl}=[jk][kl]+[kl][lj]+[lj][jk]$. Sum over disjoint subsets:

$$
f_i=\frac{1}{4860}\sum_{\substack{\{S,T\}\|S|=|T|=2,\ S\cap T=\varnothing}}W_S W_T
+\frac{1}{10368}\sum_{\substack{S,T\\|S|=2,\ |T|=3,\ S\cap T=\varnothing}}W_S W_T
+\frac{1}{17280}\sum_{\substack{\{S,T\}\|S|=|T|=3,\ S\cap T=\varnothing}}W_S W_T .
$$

Indeed the eight-slot Ward identity gives $E_i-\sum_j E_j/10=p_i f_i$. Polynomial divisibility follows from the relatively prime components of $p_i$. Summing over all eight slots gives the common half-sum above. In any of the seven other slots let $p_\ell=\epsilon q$ and maintain conservation. Each graph cofactor occurs twice, so $\partial_\epsilon\Phi_j|_0=0$. For $j\ne\ell$, $\partial_\epsilon C_j|_0=qP$ since every term retaining the soft slot in the seven-field trace vanishes quadratically. At $j=\ell$ the derivative is $-qP+4qP=3qP$ by conservation and the determinant weight four of $P$. Thus $10C_i-\sum_j C_j$ has vanishing constant and linear terms at every such soft slot. The seven remaining momenta are independent variables and each occurs at least twice in $f_i$.

The $SL_5$ determinant invariant theorem expresses $f_i$ as a linear combination of products of four maximal minors. In each product a label is omitted at most twice, by taking the multihomogeneous components in the seven slots. The four omitted pairs form a multigraph of maximum degree two. Reversing an odd path or an odd cycle kills its permutation sum. The remaining possibilities are a four-cycle, two double edges, a four-edge path, a double edge and a two-edge path, or two two-edge paths. The complementary minors obey the two-row Plücker relation. Under symmetrization it changes the four-cycle into half the two-double-edge sum and the four-edge path into minus half the double-edge and two-edge-path sum. Three independent graph evaluations fix the three displayed coefficients. The replay checks both orbit identities, the rank-three fit and all forty potential components against the four original theta integrations. In particular the bivalent part of $Q_2(\phi^2/2)$ on six scalar waves and one vector wave is $(\Phi_v-\Phi_{\rm root})(v)$, with the point absorbing the balancing momentum. [Closed expression and exact graph replay](longitudinal_theta.py)

For every polynomial $f(\phi)$ at a point this also gives a finite rule for the whole primary component of $Q_2 f(\phi)$ in the chosen bulk action and original point coordinates. Its three coefficients multiply $f'(\phi),f''(\phi),f'''(\phi)$, respectively. With five scalar waves and one vector wave the first coefficient is $p_{\rm root}(v)P$ on those six momenta and their balancing root. With six scalar waves it is the potential difference above for the second coefficient. With seven scalar waves the third coefficient is $2\Phi_v(v)$ computed with the balancing charge at the second trivalent endpoint, allowing arbitrary total momentum on the eight bivalent slots. The marked endpoint removes the endpoint-exchange automorphism. The genus-one bulk choice retains a descendant at each genus-one vertex, while a genus-two primary vertex contributes only through the one-edge tree. Every positive-loop graph with a single end at the mark has a bridge and zero heat cofactor. This exhausts the terms on primary waves. Each reduced coefficient has heat multiplier $\lambda^4/(2\pi i)^{10}$. [Three coefficients on arbitrary rational waves](longitudinal_theta.py)

## A transverse three-loop bulk term

Take three momenta $e^1$, three $e^2$, one each of $e^3,e^4,e^5$, and the balancing covector $-3e^1-3e^2-e^3-e^4-e^5$. Let $I$ be one of the $e^1$ slots and $T$ the $e^3$ slot. Insert a vector $e_2$ at $I$, a vector $e_3$ at the $e^4$ slot, the four-vector $\iota_{e^3}\eta$ at $T$, and scalars elsewhere. Here $\eta$ is the unit holomorphic top polyvector. Every primary is transverse. In the potential order of the eight-leg calculation the surviving difference is

$$
(\Phi_I-\Phi_T)(e_2)=
-\frac{4059837269}{70614415872000}
+\frac{1285913}{27243216000}\zeta(3).
$$

This is the reduced bulk potential pairing at three loops. The $e^3$ pairing with the second vector is one, and its pairing with the first is zero, so the other term in the polarized trace drops out. Both terms above come from the complete nine-graph sum. [Bulk period replay](bulk_three.py)

For the resistance reduction, write $n$ for the loop number, and use the four retained bivalent positions and the two independent cotrees of the scalar-square calculation. Let $\mu$ be the law of one labeled source from a cotree: it chooses an edge with probability $t_e(\partial U/\partial t_e)/(nU)$, then a uniform position on it. Write $\nu$ for the uniform measure on the four retained positions. Write $R_{\alpha\beta}$ for mean effective resistance integrated against two measures, and $R_{TT}$ for the mean between two *distinct* retained positions. Ground at the balancing slot. Each of the $n$ unit $e^2$ sources contributes half the difference of four pair resistances. Averaging over the two labels involved gives the fixed-network factor

$$
\frac n2\left(2R_{\mu\nu}-R_{\mu\mu}-R_{TT}\right).
$$

For the exact period integral let $\beta$ give mass $1/2$ to each trivalent vertex, and let $\rho$ be the mean resistance from a vertex to the signed measure $n\mu-\beta$. It is independent of the chosen vertex. Indeed, subdividing one edge with endpoints $u,v$ at fractional position $s$ gives resistance to a vertex $d$ equal to $(1-s)R_{ud}+sR_{vd}+s(1-s)(t_e-R_{uv})$. The resulting finite expression for $\rho U^2$ agrees as a polynomial for each core vertex in every catalog graph; subdivision leaves the edge measure unchanged along the path since the probability density is $(\partial U/\partial t_e)/U$. The same interpolation shows that the average is constant along an edge: the quadratic term times the unit total mass is canceled by the common-edge Green bridge term using the edge mass $1-R_{uv}/t_e$. Consequently the fixed-network factor equals

$$
\frac{\rho}{2n}+R_{\beta\nu}-\frac{R_{\beta\beta}}{2n}-\frac n2 R_{TT}.
$$

The two-forest numerator of each vertex resistance has denominator $U$. Integrating the squared current polynomial against the path fractions leaves core denominator $U^{2n+5}$. The calculation checks the individual graph formulas at unequal times by a separate average over edge endpoints, including the coincident-edge bridge, and checks both projective cuts for each tetrahedron period. At two loops the same reduction gives $-7/77760$, in agreement with direct integration of all labeled eight-leg theta diagrams. [Period integration](../nonSP/verify/k4_periods/README.md)

With the genus-one bulk solution chosen above, every genus-one vertex retains a descendant. For a genus-two primary vertex in a one-loop hexagon, the five singleton external momenta must have nonzero determinant. Of the $\binom{10}{5}$ choices, 42 meet this condition. In every case the five remaining momenta span at most three dimensions, and the two internal currents supply at most one extra independent covector. The genus-two polynomial $P$ therefore vanishes. A free genus-three divergence vanishes on the transverse inputs. Thus any genus-three bulk solution in these conventions must supply a nine-primary local term with a nonzero $\zeta(3)$ contribution on this projection.

This test admits a local boundary. On nine conserved momenta take the sum over five-slot subsets of the sixth powers of their determinants and multiply by the polarized trace. The fused boundary pairing in the same $C_i-C_T$ convention equals $3(2948-1674-1674+0)=-1200$. There is also a single trace matching the entire transverse tensor on the following momentum family.

Replace the first three waves by $a_i e^1$, the next three by $b_j e^2$, and balance the last wave again. Here $i,j=1,2,3$ label the waves. Each removed momentum column contributes its scalar strength to the geometric cofactor. The square supplies $\prod_i a_i^2b_i^2$. Modulo a common potential and multiples of the individual row momenta, all four retained rows are zero. Indeed their distinct pair distances have a common mean, and their mean distance to any one inserted label is independent of the retained slot. The $e^2$ component on each $a_i e^1$ row is the displayed unit-strength bulk pairing times
$$
\left(\prod_{j=1}^3 a_j^2b_j^2\right)
\frac{b_1+b_2+b_3}{3}.
$$
The $e^1$ component on each $b_j e^2$ row interchanges $a,b$. All other components vanish in this quotient: the sum of the four retained momenta belongs to the common two-plane, and the potential on an inserted row along its own line is longitudinal. This gives the full bulk tensor on the family, including singular strengths by polynomiality. Simultaneous linear transformations give the result for any independent five basis covectors with the usual determinant weight.

On nine conserved slots, sum over unordered triples of five-slot subsets which cover all nine slots and whose three pairwise intersection sizes are $2,2,3$. For each triple take the product of its three squared five-row determinants. The classical boundary of this trace symbol has unit-strength pairing $288$. With the displayed unit-strength bulk period divided by $288$ as its coefficient, its boundary agrees with the complete bulk potential on this whole family. The replay verifies all ten row wedge equations over the integer polynomial ring in the six strengths. The [eleven-primary consistency calculation](#three-loop-bulk-consistency-on-general-conserved-waves) proves existence on arbitrary conserved momenta. Four independent coefficient combinations are evaluated below. The symmetric nine-slot trace ansatz has dimension 50; the unrestricted transverse potential quotient has dimension $4013-493-795=2725$. These dimensions count, respectively, the scalar trace symbols and the covector arrays modulo common and longitudinal shifts, by the same cycle character formula as at two loops. [Polynomial trace identity](operator/three_loop_trace.py) · [Character multiplicities](operator/quartic_unisolvence.py)

The one-loop graph with the genus-two action vertex is known for arbitrary ten conserved waves. Group five labels at that vertex and ground it. If the other five momenta are $p_1,\ldots,p_5$, set $u=\sum_{j\in\mathrm{singletons}}s_jp_j$, where $0\le s_j\le1$ is the fraction of circle resistance from the ground. The other incident momentum is $\sum_{j\in\mathrm{singletons}}(1-s_j)p_j$. Evaluate $-P$ on these two currents and the five grouped momenta. Multiply by the squared determinant of the five singletons and by $\min(s_i,s_j)-s_i s_j$ for the slot-$i$ voltage from $p_j$. Integration over the cube with reflection factor $1/2$, followed by all five-label groups, gives the complete insertion. Independent degree-seven simplex cubature gives $421/447897600$ on one grounded diagonal. On the ten-wave source in the replay a transverse difference is $44224151/1175731200$, so this contribution must be retained for the general bulk equation. [Generic seven-field insertion in the hexagon](bulk_three_hexagon.py)

### Six arbitrary planar waves

The six planar momenta may be arbitrary covectors $u_1,\ldots,u_6$ in $\operatorname{span}(e^1,e^2)$. Keep the other three basis momenta and the balancing tenth wave. Write $[u_i,u_j]$ for the two-row determinant in the indicated plane. For $2n$ planar sources let $Z$ be the coefficient of $\prod_i x_i$ in $\det(\sum_i x_i u_i u_i^{\mathsf T})^n$ in this basis. Expansion of the two-dimensional determinant gives
$$
Z=n!\sum_{\substack{M\\\mathrm{perfect\ matchings\ of\ the\ planar\ labels}}}
          \prod_{\{i,j\}\in M}[u_i,u_j]^2 .
$$
For real momenta $Z>0$ precisely when at least one matching pairs nonparallel nonzero covectors. In the point scalar square with polarization $e_3+e_4$, the full reduced coefficient is $Z$ times the scalar-ring threefold coefficient. At three loops its $\zeta(3)$ part is $-3617 Z/504504000$. The marked redefinition proof above applies to all six planar covectors: a nonzero marked hexagon requires both singleton scalar directions among its four active scalars, leaving only planar spectator momenta at the corrected mark.

In fact this bare scalar-square matrix element is nonzero at every $n\ge1$ for every real planar source family with $Z>0$. The four retained scalar-ring momenta always span three directions. Start with the threefold square and two adjacent scalar vertices whose edge current lies outside the plane of their independent momenta. Add a parallel path with three edges and move the two sources to its new interior vertices. Their effective endpoint charges are the two prescribed old charges: the two columns of interpolation weights $(1-s_i,s_i)$ are independent at distinct interior positions $s_i$. Choose the path and original resistances so their parallel combination gives the old edge resistance. The chain determinant is the volume of the old reduced current and the two new charges, times a nonzero resistance factor. The new middle-edge current has a nonzero component along the old current modulo the plane of the two sources. This repeats at every loop. The four retained momenta still have rank three and sum to zero, so a change of basis realizes the same three coordinate covectors and their balancing covector at every step. On the resulting open resistance region the squared cofactor and mean pair resistance are strictly positive. [Unequal-time ladder minors through eight loops](planar_moments.py)

For the bulk, write $K$ for the reduced unit-strength period above. A second reduced period is required:
$$
J=-\frac{2078212919}{635529742848000}
 +\frac{131371}{49037788800}\zeta(3).
$$
On nine slots denote by $P_{223}$ and $P_{233}$ the sums of three squared five-row minors whose row sets together cover the slots and whose pair-intersection sizes are the three indicated integers. The unordered triples are counted once. The full three-loop bulk potential on six arbitrary planar waves is the transverse boundary of the polarized trace symbol
$$
\frac{K}{864}(P_{233}-P_{223})
+\frac{J}{168}(4P_{223}-P_{233}).
$$
Here the boundary sign is fixed by the convention $\Phi=C$ used above. The proof covers all external polyvector ranks. [Six-planar-wave period integration](planar_bulk.py) · [Two symbolic trace identities](operator/planar_trace.py)

For every $n>1$ the second bare period $J_n$ is strictly negative. At fixed positive edge resistances put $L=\sum_{e\in\mathrm{edges}}t_e C_e C_e^{\mathsf T}$ for the cycle Gram matrix, with $C_e$ the cycle column. Set $A_e=t_eL^{-1/2}C_e C_e^{\mathsf T}L^{-1/2}$ and $\mu_e=\operatorname{tr}(A_e)/n$. Since the $A_e$ sum to the identity, the correlated minus independent edge probabilities form a positive semidefinite Gram matrix with zero row sums,
$$
B_{ef}=\frac1n\operatorname{tr}\!\left[(A_e-\mu_e I)(A_f-\mu_f I)\right],
\qquad B_{ee}=(n-1)\mu_e^2.
$$
The grounded Gaussian field with resistance as its increment variance interpolates between endpoint values by independent Brownian bridges on the edges. With $G_{\rm mid}$ the covariance of the endpoint arithmetic means, two independent uniform positions on edge $e$ give integrated bridge covariance $t_e\int_0^1\int_0^1(\min(s,u)-su)\,ds\,du=t_e/12$. The one-point contributions cancel by the zero row sums. The correlated minus independent mean resistance is therefore
$$
-2\operatorname{tr}(B G_{\rm mid})-\frac{n-1}{6}\sum_{e\in\mathrm{edges}}t_e\mu_e^2<0 .
$$
The fixed-network $J_n$ transfer is $n/2$ times this difference. The squared current cofactor has positive integral on the threefold construction above. For $n$ planar copies of $(1,0)$, $n-1$ copies of $(0,1)$, and the last planar wave $(-n,1-n)$, the total planar momentum vanishes. The transverse $e_1$ projection at a unit $(0,1)$ slot against a retained slot equals $-(n!)^2(n-1)J_n>0$.

At four loops the cubic cores have six vertices and nine edges. A doubled edge joins two bivalent vertices of the simple support. Every connected bridgeless core is obtained by doubling a matching covering precisely those vertices of a biconnected support. For a connected core edge set $S$ of minimum degree two and cycle rank $h(S)$, tightness after inserting $m_e$ bivalent slots on edge $e$ requires $\sum_{e\in S}(1+m_e)\ge 3h(S)+1$. Any partial path adds an acyclic appendage, so these tests suffice. Four indistinguishable slots give the following orbit counts. Every surviving shape has a nonzero thirteen-edge current cofactor at positive integral resistances.

| Cubic core | Automorphism order | Surviving shapes |
|---|---:|---:|
| Two doubled edges | 16 | 3 |
| Three doubled edges | 48 | 0 |
| One doubled edge | 8 | 8 |
| Triangular prism | 12 | 28 |
| $K_{3,3}$ | 72 | 11 |

[Exact census, literal current determinants and independent metric-graph resistance test](quartic_planar_cores.py)

On the two-doubled-edge core each bubble carries two retained slots. The three shapes correspond to path splits $(2,0),(2,0)$, then $(2,0),(1,1)$, then $(1,1),(1,1)$. The graph weights in the normalization of $K,J$ are

| Splits on the two bubbles | $K_4$ | $J_4$ |
|---|---:|---:|
| $(2,0),(2,0)$ | $-1453/136914624000000$ | $-541/8762535936000000$ |
| $(2,0),(1,1)$ | $-4877/410743872000000$ | $-1889/26287607808000000$ |
| $(1,1),(1,1)$ | $-73/22819104000000$ | $-269/13143803904000000$ |

Squaring the linear current minor and integrating each pair of positions at fixed path resistances gives $2488$ and $13824$ monomials in the respective core integrands. Series and parallel beta integrations sum to
$$
K_4^{\rm SP}=-\frac{211}{8214877440000},\qquad
J_4^{\rm SP}=-\frac{1}{6490767360000},\qquad
-1728 J_4>\frac{9}{33806080000}.
$$
The strict bound includes the common sign of the other $47$ shapes. On every four-loop cubic core, exact polynomial division cancels one $U$ in the canonical-resistance numerator and two $U$ in the correlated-minus-independent resistance numerator. Both integrands then have denominator $U^{12}$. The series-parallel term counts reduce to $261$ and $253$.

The first series-parallel shape has a sign-changing $K_4$ density even after integrating the four retained positions. Its reduced numerator equals $-834/49$ at unit core resistances and $757810515778611603$ at $(1,7,7,1,7,7,1,1,7)$. The respective Kirchhoff denominators are $45^{12}$ and $12993^{12}$.
Summing all three shapes and averaging over all sixteen core symmetries still gives $-3408/49$ at unit resistances and $65625600000/49$ at $(1,1,1,10,10,1,10,1,1)$.

On the one-doubled-edge core, contract the parallel pair and suppress its two bivalent endpoints. For powers $\alpha,\beta$ the parallel contraction contributes $\mathrm B(12-\alpha,12-\beta)$ with new power $\alpha+\beta-12$; the series move contributes $\mathrm B(\alpha,\beta)$ with new power $\alpha+\beta$. Here $\mathrm B(i,j)=(i-1)!(j-1)!/(i+j-1)!$. The eight shapes reduce to $1258$ tetrahedron exponent orbits. Their two propagator cuts give
$$
K_4^{\rm one\ double}
=\frac{11144483303}{16895050672500000}
-\frac{186919}{340540200000}\zeta(3),\qquad
J_4^{\rm one\ double}
=\frac{197793221359}{36248733480960000000}
-\frac{507937}{111878807040000}\zeta(3).
$$
The other $39$ shapes have strictly negative $J_4$ contributions, so
$$
-1728 J_4>
\frac{507937}{64744680000}\zeta(3)-
\frac{173064182111}{18355116780000000}
>\frac{9}{33806080000}.
$$
[Exact rational and tetrahedral periods and unequal-time current checks](quartic_planar_periods.py)

On the prism and $K_{3,3}$ the cut is a massless three-loop propagator in dimension $24$. For the nine exponents $\alpha_i$, numbered in the propagator replay order, the projective monomial period equals the cut integral at unit external momentum squared times $\prod_{i\in\mathrm{edges}}(\alpha_i-1)!/11!$. In the prism's first reduced shape the orbit of $(4,4,5,9,6,4,5,5,6)$ has $J_4$ numerator coefficient $-1/14$. Ten dimension shifts and two inequivalent cuts of that monomial give the projective period
$$
-\frac{531372573405073}{2048959660011264000000}
-\frac{80588077}{116154175737600}\zeta(3)
+\frac{15}{14226212}\zeta(5).
$$
On the first $K_{3,3}$ shape the orbit $(4,4,5,4,5,6,5,9,6)$ has $J_4$ numerator coefficient $9/700$. Cutting at two different powers gives the same projective moment
$$
-\frac{186964375676833}{9579551657195520000}
-\frac{53055501437}{1016349037704000}\zeta(3)
+\frac{1693}{21339318}\zeta(5).
$$
These fix the $\zeta(5)$ terms of two scalar integrals. The coefficient in the full primitive graph sum requires its other moments. Fraction-free elimination of the current determinant and ordered path moments preserve the $K_4,J_4$ graph numerators. The weighted orbit file stores the nine positive powers followed by the two rational weights including the cut normalization. [Primitive numerator construction and two-cut scalar controls](quartic_master_weights.py)

The same weights give rational convergent series for both complete primitive periods. Divide each weight by its cut normalization before applying the [projective Gegenbauer expansion](../nonSP/verify/hepp_series/README.md). Its radii on the prism and $K_{3,3}$ are respectively $123/155$ and $4/5$, and its Gegenbauer parameter is $12$. For a numerator monomial, the tropical moment on a proper edge subset $S$ has divisor $\sum_{i\in S}\alpha_i-12 h(S)$, where $h(S)$ is its cycle rank. Each divisor is an integer between $1$ and $48$. Multiplication of the subset moment by $\operatorname{lcm}(1,\ldots,48)^{|S|}$ makes the recurrence integral at every proper subset. The projective moment on all nine edges is the sum over its eight-edge predecessors. The native replay computes the signed and absolute zeroth moments of both numerators using these integers. [Native exact moment](quartic_hepp_moment.cpp) · [Finite partial sums and absolute tails](quartic_master_weights.py)

| Core | Nonzero weight orbits for $K_4,J_4$ | Union | Numerators at nine unit resistances $K_4,J_4$ |
|---|---:|---:|---|
| [Prism](quartic_prism.weights.gz) | $1598521,\ 1597581$ | $1598657$ | $-103417964820487/308700,\ -472498704667/403200$ |
| [$K_{3,3}$](quartic_k33.weights.gz) | $465544,\ 465496$ | $465564$ | $-18390175079148/214375,\ -1675082349/4375$ |

The final column also follows by substituting the nine path lengths before integrating the four positions. Both numerator banks retain their original graph symmetry factors. The negative of the $J_4$ density is a positive measure after averaging each representative over the core automorphism group. With $x=U/\max_I t_I$ for the spanning-tree complements $I$, the integral of $x^{-12}$ against its tropical measure is $-J_4$ on that core. Four rational moments, of orders $0,1,2,3$, give rational Gauss and Lobatto bounds. To evaluate either rule, invert $x^{12}$ modulo the polynomial of its nodes and integrate the resulting polynomial using the four moments. The lower-rule Hermite error has the sign of the product of two squared linear factors. For the upper rule the single factors at $1$ and the tree count have opposite signs on the interior. [Moment coefficients and input digests](quartic_primitive_moments.json) · [Exact rational quadrature bounds](../nonSP/verify/hepp_series/moment_bounds.py)

At order three the prism and $K_{3,3}$ contain $14717581$ and $3989677$ exponent orbits after multiplication by $U^3$. Independent four-position integration also reproduces both numerator banks with each of the nine edges in turn doubled in resistance and the nine values averaged. The four moments give the rational certificates
$$
-J_4^{\rm prism}>\frac{973289}{10^{20}},\qquad
-J_4^{K_{3,3}}>\frac{124891}{10^{20}},\qquad
-1728J_4>\frac{1772146069}{10^{18}}.
$$
The last bound uses the exact eleven-graph subtotal, the unrounded rational quadrature endpoints, and eighty terms of Euler's transform of the alternating series for $\zeta(3)$. [All-core bound replay](quartic_master_weights.py)

#### Integration by parts on the complete numerator

Write $I_q(P)$ for the projective integral of $P/U^q$. For a convergent monomial $t^m/U^{q-1}$ and edge $e$, differentiate $t_e t^m/U^{q-1}$ in a chart fixing a different edge. All positive cotree terms are separately convergent. The two endpoint values vanish, giving
$$
(q-1)I_q(t^m t_e\partial_e U)=(m_e+1) I_q(t^m U).
$$
Writing $U=A+t_eB$, integration of $t_e^{a-1}/(A+t_eB)^q$ gives $\mathrm B(a,q-a) A^{a-q}B^{-a}$. A numerator factor $B^a$ leaves a contracted graph at power $q-a$; a factor $A^{q-a}$ leaves a deleted graph at power $a$. Every edge deletion in the prism and $K_{3,3}$ reduces to a series-parallel or tetrahedral period.

On the $1481$-moment prism control bank at $q=6$, the $1323$ homogeneous identities have rank $852$. Adding $111$ independently integrated face identities gives rank $877$ over the rationals and modulo $1009$. A dense input supported on all $1481$ moments reduces to $604$ unknown integrals and known rational and $\zeta(3)$ terms, reproducing its complete rational, $\zeta(3)$ and $\zeta(5)$ columns. Sampling $256$ orbit representatives of $t^m$ in each actual $q=12$ bank gives $1277$ prism and $1818$ bipartite pivots. The affected prism $K_4,J_4$ rows fall from $13186,13183$ to $11871,11868$; the bipartite rows fall from $15232,15232$ to $13409,13408$. Reversing the row operations recovers every original rational weight. [Exact elimination and closed-face controls](quartic_projective_ibp.py)

In the replay's cycle basis the prism Gram matrix has $L_{14}=0$, $L_{11}=t_1+t_2+t_5$ and $L_{44}=t_3+t_4+t_8$ (zero-based edges). Let $D_S$ be the determinant of the symmetric matrix of entry derivatives on the indicated cycle slots, with half the ordinary derivative on off-diagonal entries. Schur complementation reduces a principal $k$-row derivative to the derivative of a $k$-dimensional symmetric determinant. At the identity the two-row coefficient is $s^2+s/2$ and the three-row coefficient is $s^3+3s^2/2+s/2$. Congruence transports them to each invertible block. Consequently
$$
D_{23}U^{s+1}=(s+1)(s+3/2)L_{11}L_{44}U^s,\qquad
D_{123}U^{s+1}=(s+1)(s+3/2)(s+2)L_{44}U^s .
$$
The product rule for $D_{234}(L_{44}U^s)$ and $\det L_{\{1,4\}}=L_{11}L_{44}$ give
$$
D_{234}D_{123}U^{s+1}=(s+1)(s+3/2)(s+2)^2 D_{23}U^s.
$$
For the monomial $H$ with exponent vector $(5,4,4,4,4,4,5,5,6)$, all endpoint terms vanish and the least proper-subset margins of the two differentiated numerators are $3$ and $2$. Therefore
$$
11550 I_{12}(D_{23}H)=I_{11}(D_{234}D_{123}H).
$$
All four orbit representatives on the left occur in both physical banks on the prism. The two-row identity also reproduces $84$ finite relations in the integrated $q=6$ control. [Literal graph derivatives, endpoint tests and integrated checks](quartic_prism_shift.py)

For $K_{3,3}$ the normal $N$ to the nine cycle squares has $N_{14}=N_{41}=-1$, $N_{23}=N_{32}=1$ and all other entries zero. Thus $N^2=1$, $\det N=1$ and $\operatorname{tr}(NL)=0$. Choose a complex matrix $S$ with $S^TS=N$ and put $X=SLS^T$. This symmetric matrix is traceless and has $\det X=U$. Write $e_j$ for the degree-$j$ elementary symmetric function of a matrix and $\nabla$ for the symmetric entry gradient with its trace removed. For the triples $\alpha,\beta$ in the coefficient table, $e(X)^\alpha$ denotes the corresponding product of $e_2,e_3,e_4$; the same convention applies to $e(\nabla)^\beta$. The nine rational coefficient functions $c_{\alpha\beta}(s)$ satisfy
$$
\sum_{\substack{\alpha,\beta\\ \text{nine table entries}}}
c_{\alpha\beta}(s)e(X)^\alpha e(\nabla)^\beta U^s=U^{s-1}.
$$
To calculate each symbol at a diagonal $X$, expand $\det(1+Z)^s$ through its derivative order. If $T_j$ denotes its degree-$j$ part, logarithmic differentiation gives $nT_n=\sum_{j=1}^{\min(n,4)}((s+1)j-n)e_j(Z)T_{n-j}$. Diagonal rescaling restores the inverse eigenvalues of $X$. Rewrite the result in elementary symmetric polynomials and set $e_3(X^{-1})=0$, since $\operatorname{tr}X=0$. Each column is a polynomial of degree at most eight in $s$. Nine exact Taylor evaluations reconstruct each column, and substitution of the rational solution into every polynomial row gives zero residual. Orthogonal conjugation on the dense open set of symmetric matrices with simple spectrum extends the identity to the traceless hyperplane after clearing denominators. Two evaluations at positive nine-edge metrics use the original graph coordinates and an independent Taylor expansion of $U$. The common denominator is nonzero at every integer $s=-3,\ldots,-11$. [Polynomial proof and graph-coordinate replay](quartic_bipartite_shift.py) · [Exact coefficients in ascending powers of $s$](quartic_bipartite_shift.json)

Boundary terms can be retained before specializing the propagator powers. Let $a$ be their nine-component vector with $\sum_{e\in\mathrm{edges}} a_e=4q$. The cut at unit external momentum squared, with momentum measure $d^{2q}k/\pi^q$ per loop, equals
$$
\frac{\Gamma(q)}{\prod_{e\in\mathrm{edges}}\Gamma(a_e)}
I_q\!\left(\prod_{e\in\mathrm{edges}}t_e^{a_e-1}\right).
$$
Indeed, integrate the cut edge using $U=A+t_eB$ and then integrate the common scale of the eight uncut edges in the Schwinger representation of the propagator. The respective factors are $\mathrm B(a_e,q-a_e)$ and $\Gamma(q-a_e)$. Denote the displayed cut by $F_q(a)$. Expand the graph-coordinate operator in the preceding identity as a sum of $c_{\mu\nu}(1-q)t^\mu\partial^\nu$ with $|\nu|-|\mu|=4$. Integration by parts gives the meromorphic identity
$$
F_q(a)=(q-1)\sum_{\substack{\mu,\nu\\ \text{graph-coordinate operator}}}
(-1)^{|\nu|}c_{\mu\nu}(1-q)
\left(\prod_{e\in\mathrm{edges}}\prod_{j=0}^{\mu_e-1}(a_e+j)\right)
F_{q-1}(a+\mu-\nu).
$$
Each derivative monomial leaves one of the nine variables undifferentiated, since its order is at most eight. Fix that edge to one for its integration. Near $a_e=4q/9$ at sufficiently large real $q$, every original and shifted proper-subset margin is positive. Indeed $4|S|>9h(S)$ for every nonempty proper edge subset $S$ of this graph; the finitely many shifts change the margins by bounded integers. Powers and derivatives also vanish at the requisite coordinate endpoints. Moving the derivatives produces falling factorials of $a_e+\mu_e-1$. Combining these with the shifted gamma normalization gives precisely the displayed rising products. In each ordered Schwinger sector a successive ratio for an edge subset $S$ has exponent $\sum_{e\in S}a_e-qh(S)-1$. The remaining polynomial factor has positive constant term by the greedy spanning tree. Successive one-variable Taylor subtractions at zero meromorphically continue both sides. A generic dimensional regulator with the compensating power on the cut edge keeps all eight uncut indices integral. Nonpositive shifted uncut indices then correspond to ordinary momentum numerators or contracted lines. When shifted terms have poles, take their finite part only after adding their complete Laurent coefficients.

At $q=3$, every positive convergent propagator power is either one or two, since the complement of any edge has cycle rank three. Their sum is twelve. Enumerating all placements of the doubled powers and testing all proper subsets gives six prism, three bipartite and nine wheel orbits. All eighteen periods have been integrated at two cuts in dimension six. For instance the bipartite powers $(1,1,2,1,2,2,1,1,1)$ give $-2\zeta(3)+5\zeta(5)/2$. The full primitive numerators are evaluated below by convergent fixed-dimension relations. [Two-cut integrals](quartic_terminal_periods.py) · [Terminal periods](quartic_terminal_periods.json)

For a finite check retaining endpoint distributions, take $H=t_0t_6t_7t_it_j$ with $i\in\{1,2,5\}$ and $j\in\{3,4,8\}$. Applying $4D_{23}$ to $U^{-3}$ gives $30 L_{11}L_{44}U^{-4}$. On each rung the second derivative integrates as $\int_0^\infty t_e\,\partial_e^2(U^{-3})\,dt_e=U^{-3}|_{t_e=0}$. The mixed derivatives give twice the three bulk periods obtained by removing the corresponding pair of rung factors from $H$. The diagonal derivatives subtract the three wheel periods on the contracted rungs. All bulk and face integrals have positive proper-subset margins. For $i=1,j=3$ the resulting identity in the coefficient order $1,\zeta(3),\zeta(5)$ is
$$
30(1/20,1/45,-1/18)=(-1/3,-11/3,20/3)+(11/6,13/3,-25/3).
$$
The last triple is the signed face subtotal. The nine choices of $i,j$ give two distinct identities, with six different upper prism moment orbits. Both cuts verify every coefficient. [Integrated contact replay](quartic_prism_contacts.py)

#### Period class of the complete four-loop pair

The complete $K_4,J_4$ belong to $\operatorname{span}_{\mathbb Q}\{1,\zeta(3),\zeta(5)\}$. Set $d=4-2\epsilon$ and cut one edge in each absolutely convergent primitive Schwinger monomial. Denote by $B$ the Gram determinant of the external and the three loop momenta at unit external momentum squared. The ten Gaussian shifts convert the $d+20$ scalar cut to the $d$ cut with numerator $B^{10}$, multiplied by
$$
\prod_{j=0}^9\frac{8}{(d+2j-1)(d+2j-2)(d+2j-3)} .
$$
This factor is regular at $d=4$. Expanding $B^{10}$ requires only the pole and constant coefficients of ordinary three-loop massless propagator integrals with integer indices and polynomial numerators. In the $G$-scheme their Laurent coefficients through the constant term are rational combinations of $1,\zeta(3),\zeta(4),\zeta(5)$, with the only irrational pole allowed proportional to $\zeta(3)/\epsilon$. The no-even-zeta theorem removes $\zeta(4)$ in the finite sum. The common change of scheme has unit limit on the finite whole. The projective normalization $\prod_e(a_e-1)!/11!$ and both complete graph weights are rational. These facts fix the space of constants in each physical sum. The complete prism coefficients are evaluated below. [Three-loop propagator structure, Table 1 and §6](https://arxiv.org/pdf/1804.10088#page=13) · [Finite-combination theorem 4 and its three-loop case, §7.2](https://arxiv.org/pdf/1004.1153#page=28)

#### Exact prism periods under conformal transformations

In the dual position integral the two triangular and three rectangular faces are five points. Each edge gives the squared distance between its adjacent faces to the negative of its propagator power. Cut a rung, fix its adjacent points at $0$ and a unit vector, and integrate the three other points in dimension $2q$, with $q=12$ for this bank. The two triangle points have no mutual edge. Add a sixth point at infinity with powers $a_{i\infty}=2q-\sum_{j\ne i,\,j\text{ finite}}a_{ij}$. All six incident sums are $2q$ since the nine original powers sum to $4q$. Fixing any finite separation to one gives the same integral: translate one endpoint to the origin, integrate overall scale against a radial cutoff of unit logarithmic integral and use rotation invariance. Inversion about either unit-separation endpoint interchanges it with infinity with the displayed sixth-point powers. Thus all positive prism charts of one completed weighted graph have the same normalized cut integral. The tabulated physical weights multiply exactly this cut integral.

If one completed point has at most three positive incident powers, send a positive neighbor to infinity and integrate the point with two remaining denominators of powers $b,c$. The incident sum $2q$ and the upper bound $q-1$ on each positive power require at least three positive powers. Write $m$ for the total power of its polynomial numerator in squared distances. Then $b+c=2q-a_{v\infty}+m$. Combine the denominators with a Feynman parameter $t$ and shift to their weighted center. For an even radial moment of order $2j\le2m$, the radial gamma factors contain $\Gamma(q+j)\Gamma(b+c-q-j)$, with $b+c-q-j\ge q-a_{v\infty}>0$. The parameter powers are $t^{q-c+j-1}(1-t)^{q-b+j-1}$ times a polynomial in $t$. Their endpoint indices are positive integers. Angular metric contractions and these beta integrals therefore give a rational linear combination of integer powers of the remaining squared distances. Its combined two-loop integral is finite by the original positive integral.

Continue only the two remaining loop measures to dimension $2q-2\epsilon$ while keeping their rational integrand fixed. In the Gram coordinates of the two momenta and the unit external vector the measure changes by a power $-\epsilon$ of the Gram determinant and an analytic factor. The convergent integral extends to a neighborhood of zero by the strict endpoint margins for its combined rational integrand. The five remaining variable squared distances generate all its loop scalar products. Triangle reduction and bubble integration give two gamma products with rational-function coefficients even when individual scalar terms require continuation. After extracting their common factor of unit limit, their ratio is $1-6\zeta(3)\epsilon^3+O(\epsilon^4)$. A pole beyond order three in its rational multiplier would leave an irrational divergence by Apéry's theorem, since the other coefficient is rational. The combined finite value belongs to $\mathbb Q+\mathbb Q\zeta(3)$. [Two-loop reduction and gamma ratio](../nonSP/verify/k4_periods/README.md)

There is a further identity at a point with four incident powers $0<a_i<q$, $\sum_{i\in\mathrm{neighbors}}a_i=2q$. Its one-loop Schwinger integral has denominator $D(t)^q$ with $D(t)=\sum_{i<j,\,\mathrm{neighbors}} x_{ij}^2t_it_j$ and numerator $\prod_i t_i^{a_i}$ in the logarithmic projective measure. Here $x_{ij}^2$ is the squared distance between neighbors. For complementary pairs $(i,j),(k,l)$ choose positive $d_i$ with $d_i d_j=x_{kl}^2/x_{ij}^2$. These six equations are consistent because opposite right-hand sides are reciprocal. The positive projective substitution $t_i=d_i/u_i$ gives $D(t)=D(u)/\prod_i u_i$. Thus replace the four incident powers by $q-a_i$ and each power joining two of their neighbors by $a_{ij}+a_i+a_j-q$. The original cut integral is the transformed one times
$$
\prod_{i\in\mathrm{four\ neighbors}}\frac{(q-a_i-1)!}{(a_i-1)!}.
$$
Every step is an identity of positive convergent integrals. Permuting all six points and repeating the four-point substitution reaches a finite set of integer weights. A collision of any two finite integration points forces every positive weight to be less than $q$; any pair can be made finite by the preceding change of chart. The six row sums then bound negative weights as well. The replay checks the exact multiplier around every cycle it encounters.

Of the $1{,}598{,}657$ physical rows, $142{,}571$ immediately reduce to the two-loop class. Reassigning the six points leaves $674{,}287$ distinct cuts among the other $1{,}456{,}086$ rows. The four-point substitutions reach a tensor bubble from another $73{,}313$ of these cuts. The rest form $166{,}420$ classes; in one class both physical weights sum to zero. Thus the weight-five coefficient has $166{,}414$ nonzero class weights for $K_4$ and $166{,}394$ for $J_4$, with union $166{,}419$. In the $q=6$ control the six-point permutations give $469$ independent duplicate rows with equal values in all three period columns, and the $159$ direct bubble tests all give zero in the weight-five column. The reflected reduction of two dense control numerators reproduces $-7741325/27216$ and $-150875/6048$ for their fully integrated weight-five coefficients. The measured prism orbit $(4,4,5,9,6,4,5,5,6)$ has sixteen physical representatives after reflection. Their complete contributions to the coefficients multiplying $\zeta(5)$ in the prism parts of $K_4,J_4$ are
$$
-\frac{4363040109}{69708438800},\qquad -\frac{322449769}{1912002892800}.
$$
The full weight-five projection is evaluated below. [Exact controls](quartic_prism_completion.py) · [Weighted orbit reduction](quartic_prism_completion.cpp)

For an integrated vertex $v$, choose a reference neighbor $b$ in a chart with $i$ at infinity. Denote the normalized cut by $F(a)$ and let $e_{jk}$ add one to the power on pair $jk$. Integration of the divergence of $x_v-x_b$ times the full density gives
$$
\begin{aligned}
(a_{vi}-a_{vb})F(a)
={}&\sum_{\substack{j\ne v,b,i\\\mathrm{finite\ neighbors}}}a_{vj}
\big[F(a+e_{vj}+e_{bi}-e_{vb}-e_{ji})\\
&\hspace{38mm}-F(a+e_{vj}+e_{bi}-e_{jb}-e_{vi})\big].
\end{aligned}
$$
Here $2(x_v-x_b)\cdot(x_v-x_j)=x_{vb}^2+x_{vj}^2-x_{jb}^2$; completing the shifted finite weights gives the changes at $i$. Keep rows with all terms convergent. For a collision of a set $S$ of points the radial margin is $q(|S|-1)-\sum_{j<k\in S}a_{jk}$. Complementary sets have equal margins by the six incident sums. Pair and triple inequalities therefore cover all proper collisions after fixing three points. Nested collisions give the same inequalities on successive scales. Positive integral margins kill the boundary flux in this integration by parts.

At $q=6$ enumerate both ordered four-tuples of apex powers with sum $12$, each entry between $-3$ and $5$. The four remaining incident sums on the equator determine its six edges from two integers between $-8$ and $5$. Testing margins and dividing by point permutations gives $5447$ completed charts. Reflections and the weight-five bubble zeros leave $727$ variables. The displayed relations have rank $726$ modulo $1009$. Fix one normalized cut coefficient to its independent value $-1087/2268$, solve the selected system over the rationals, and substitute the nonzero vector into every generated rational relation. This proves exact rank $726$. The solution has $721$ nonzero components and reproduces all $1481$ independently integrated weight-five controls. This establishes a one-dimensional space of weight-five coefficients for these identities at $q=6$. The $q=12$ system uses the same relations with its own integration bounds. [Finite enumerator and rational solution](quartic_prism_position_ibp.py)

At the physical value $q=12$, enumerate the two apex four-tuples of sum $24$ with entries between $-9$ and $11$, and the two independent equator weights between $-20$ and $11$. There are $2{,}959{,}630$ convergent charts up to the six point permutations. Tensor bubbles and four-point reflections leave $271{,}099$ variables, including all physical prism classes. Of $47{,}593{,}846$ nonzero position-space identities modulo $1000000007$, choose one for each attainable highest variable, breaking ties by support size and then the sum of the variable indices. This gives $270{,}524$ triangular pivots. Evaluating the residual equations on blocks of $24$ of the $575$ free coordinate vectors produces another $574$ independent rows. Thus one independently measured normalized cut fixes the entire weight-five projection. For powers $(4,4,5,9,6,4,5,5,6)$ its coefficient is
$$
\frac{15}{14226212}\frac{11!}{\prod_{e\in\mathrm{edges}}(a_e-1)!}
=\frac{1}{41191991672832}.
$$
Six $30$-bit primes reconstruct all $271{,}099$ fractions. Their numerators and denominators have at most $62$ and $83$ bits. Substitution over the rationals verifies the $270{,}524+574$ equations and the normalization. Exactly $267{,}543$ variable coefficients are nonzero. This computes coefficients in the explicit period decomposition: one may equivalently perform the argument modulo $\mathbb Q+\mathbb Q\zeta(3)$ throughout, with the rational and $\zeta(3)$ channels discarded by the bubble test.

Summing the original $1{,}598{,}657$ rationally weighted orbits gives the complete prism contributions
$$
[\zeta(5)]K_4^{\rm prism}=-\frac{13914257373753679}{45754008300000},
\qquad
[\zeta(5)]J_4^{\rm prism}=-\frac{1633443038914372549}{1045387581638400000}.
$$
The $166{,}419$ independently summed four-point-reflection weights reproduce both fractions. All $47{,}726{,}374$ generated identities with retained terms also vanish modulo the separate prime $1000000009$. Independent evaluations in the two inequivalent propagator cuts give $575/92681981263872$ at powers $(5,5,5,5,5,5,6,6,6)$ and $227/27259406254080$ at $(4,6,5,6,5,5,6,6,5)$, equal to the predictions of the fixed linear system. The full system below retains these tensor-bubble terms. [Native exact solver](quartic_prism_physical_ibp.cpp) · [Reproducible physical sum and controls](quartic_prism_physical_ibp.py)

Retaining the tensor bubbles gives $899031$ reflected variables at $q=12$. The same convergence-safe identities give $887415$ direct pivots, $11277$ additional sparse pivots and $336$ further independent equations modulo $1000000007$. The remaining three coordinates can be fixed by elementary or independently integrated graphs. With $b=2q/3$, one scalar bubble with two unit-separated endpoints has value
$$
\frac{\Gamma(q-b)^2\Gamma(2b-q)}{\Gamma(b)^2\Gamma(2q-2b)}.
$$
The completed graph with three points in each part and all nine cross-powers $b$ gives the cube of this value by sending one point to infinity. A single such bubble times the symmetric two-loop $K_4$ cut with all six powers $q/2$ gives the second anchor. Subtract $2b-q$ from the endpoints' existing connecting power before inserting the bubble. The third anchor is the previously integrated prism $(4,4,5,9,6,4,5,5,6)$. In the order $1,\zeta(3),\zeta(5)$ the three $q=12$ rows are
$$
\begin{pmatrix}
1/208215748530929664000000000&0&0\\
-3541028833/9363576310742384640000000000&2273/7224981721251840000000&0\\
-531372573405073/88991429257295623618560000000&
-80588077/5044865604155080704000&1/41191991672832
\end{pmatrix}.
$$
The anchor equations lift the selected matrix to full rank modulo $1000000007$. Ten primes reconstruct the $11616$ direct-triangular free coordinates in all three columns. Substitution over the rationals recovers all other coordinates and checks each selected relation. Clearing common denominators and packing the three integer residuals in radices exceeding their individual bounds verifies all $102596912$ generated identities exactly. The three triangular anchor rows prove the coefficient columns linearly independent, so the exact rank of the entire homogeneous system is $899028$. Thus the actual normalized periods are the displayed linear combinations of the anchor values, independently of any rational-independence conjecture for the zeta values. At $q=6$ all three columns reproduce all $1481$ independently integrated cuts. [Scalar bubble and two-loop reduction](../nonSP/verify/k4_periods/core.py) · [Independent physical anchor](quartic_master_weights.py)

The original $1598657$ physical rows and the regrouped $241543$ reflected class weights independently give
$$
\begin{aligned}
K_4^{\rm prism}={}&
\frac{43896298966379923306171160083}{586747539003225600000000000}
+\frac{768730384185255903297563}{3841799362521120000000}\zeta(3)
-\frac{13914257373753679}{45754008300000}\zeta(5),\\
J_4^{\rm prism}={}&
\frac{476340584507847676272946020919}{1239210802374812467200000000000}
+\frac{758352372011658482516989}{737625477604055040000000}\zeta(3)
-\frac{1633443038914372549}{1045387581638400000}\zeta(5).
\end{aligned}
$$
The two physical cuts $(5,5,5,5,5,5,6,6,6)$ and $(4,6,5,6,5,5,6,6,5)$ reproduce their independent FORCER values in all three columns. Sixty-four terms of the Euler-transformed alternating zeta series give
$$
-\frac{2382}{10^{12}}<K_4^{\rm prism}<-\frac{2381}{10^{12}},
\qquad
-\frac{1043}{10^{14}}<J_4^{\rm prism}<-\frac{1042}{10^{14}}.
$$
For integer $s>1$ the alternating binomial sum $\sum_{\substack{0\le k\le j\\\mathrm{finite\ difference}}}(-1)^k\binom jk(k+1)^{-s}$ is the moment of $(1-t)^j$ in the probability measure $(-\log t)^{s-1}dt/(s-1)!$ on $[0,1]$. This bounds the positive remainder by $2^{-64}/(1-2^{1-s})$. [All three physical columns and controls, run with `--all-periods`](quartic_prism_physical_ibp.py) · [Independent rows of the rational system](quartic_prism_rank_rows.gz)

The following independent calculation isolates the nonplanar master contribution. For the prism, the cut three-loop propagators are planar. Panzer's Theorem 1.1 gives ordinary multiple zeta values for each Laurent coefficient of the planar masters. For the nonplanar core there is one additional master, the unit-power eight-line nonplanar propagator in $d=4-2\epsilon$, whose finite value is $20\zeta(5)$. Its ordinary-master coefficient for every integer-power cut in dimension $24-2\epsilon$ is regular at $\epsilon=0$. [Planar and nonplanar propagators, Theorem 1.1, Figure 2 and (3.28)](https://arxiv.org/pdf/1305.2161) · [Five-master basis and Gram normalization](../nonSP/verify/k33_forcer/README.md)

Let $y_1,\ldots,y_8$ be the squared internal momenta of the $K_{3,3}$ cut with unit external momentum squared. Put $x=-2p_6\!\cdot p_7$ in the repository routing and let $B$ be the Gram determinant of the external momentum and $p_6,p_7,p_8$. Directly solving the nine linear squared-momentum equations gives $B(0,x)=x^2(1-x)^2/16$. In dimension $d$ the Gram measure has power $(d-5)/2$. Ten dimension shifts insert $B^{10}$ with a rational factor regular at $d=4$. Taking the eight contour residues at $y_i=0$ for positive propagator powers $a_i$ selects the coefficient of $\prod_{i=1}^8y_i^{a_i-1}$ in $B^{19/2-\epsilon}$. Every lower sector has a missing denominator and gives zero on this contour. The ratio to the same contour of the unit nonplanar master therefore gives precisely its coefficient in the ordinary basis.

Expand the required coefficient at $y=0$. On the remaining interval $0<x<1$, each term is a polynomial in $x$ times $[x(1-x)]^{19-2k-2\epsilon}$ for an integer $k\ge0$, with a scalar factor holomorphic near $\epsilon=0$. Termwise integration by continuation gives beta functions $\mathrm B(A-2\epsilon,B'-2\epsilon)$ with integers $A,B'$. If both gamma functions in the numerator have a pole then so does $\Gamma(A+B'-4\epsilon)$ in the denominator. Thus each term has at most a simple pole. The unit master contour is a nonzero analytic factor times $\mathrm B(-2\epsilon,-2\epsilon)$, with a simple pole. The common contour phases cancel for integer shifts of the powers. Dividing proves regularity. Hence the unit master contributes only its finite value $20\zeta(5)$ times a rational number. The other four ordinary masters have missing lines and are planar. Their Laurent coefficients, including those needed for spurious poles, belong to the ordinary multiple-zeta algebra. This argument applies to the full rationally weighted banks. The prism evaluation above supplies all of its coefficients.

More explicitly, at $d=4$ the cut coefficients are rational functions of $x$. The unit integrand is $4/[x(1-x)]$. For the dimension raised by $2m$, take the coefficients of $x^{-1}$ and of $(1-x)^{-1}$ at the two endpoints in the required $y$ coefficient of $B^{m-1/2}$. Divide their sum by $8$ and multiply by $\prod_{j=0}^{m-1}8/[(2j+1)(2j+2)(2j+3)]$. The two residues extract the common simple pole before taking the ratio. Twenty-one ordinary graph reductions in dimensions four, six and eight check this formula independently. For instance, in dimension six the uncut powers $(1,1,1,1,1,2,2,1)$ have projection $-1/12$, and in dimension eight $(1,2,1,1,3,1,2,2)$ have projection $1/240$. [Endpoint residue and independent FORCER comparison](quartic_nonplanar_cut.py)

At either endpoint translate the auxiliary variable by a linear form in $y$ and split $16B$ into homogeneous pieces of degrees two, three and four. Translation preserves the formal residue at each order, since its positive-order terms are total derivatives. In the first chart the quadratic is $x^2-4y_6y_7$ and in the second it is $x^2-4y_5y_8$. Give $x$ weight one and the $y_i$ nonnegative integer weights $w_i$. If every monomial in $16B$ has weighted degree at least two, then every term in the expansion of $B^{19/2}$ at that endpoint has weight at least nineteen. A necessary condition to contribute to the $x^{-1}y^{a-1}$ residue is $\sum_{i=1}^8w_i(a_i-1)\ge20$. Checking weights in $\{0,1,2\}^8$ gives four componentwise minimal tests at each endpoint. Applied to all nine cuts these exclude the nonplanar master in at least one cut for $395810$ of the $465564$ positive-power orbits in the complete $K_{3,3}$ numerator bank. The other $69754$ orbits survive these inequalities in every cut.

For the physical control with ordered powers $(9,6,5,6,4,5,4,5,4)$ cut at the first edge, the projection is zero. A direct reduction at $24-2\epsilon$, followed by the exact beta-function ratio for the unit master, checks this zero independently. The full projective $\zeta(5)$ coefficient $1693/21339318$ in this control therefore comes from planar lower sectors in the stated cut. At another physical orbit, $(1,1,3,3,3,9,9,9,10)$ in the displayed order, the nonplanar projection is $-1004423629/95274484236288000$. To compute the residues, scale all variables at the endpoint together. At total $y$ degree $n$ only the extra degree $r=n-2m$ above the quadratic contributes. The cubic and quartic pieces give a two-step power-series recurrence through degree $r$. A finite rectangular coefficient array, modular reconstruction and an explicit absolute integer bound evaluate the resulting coefficient exactly. The complete period is cut independent; its decomposition into master coefficients can depend on the cut. [High-dimensional cross-check and cut-selection replay](quartic_nonplanar_cut.py)

The nonplanar master has also been projected out of both **complete** $K_{3,3}$ banks. Choose one of the certified zero-master cuts wherever available. On the remaining $69754$ orbits put a highest-power line on the cut, taking the lexicographically least uncut representative if necessary. In the order $K_4,J_4$, the exact coefficients multiplying the unit four-dimensional nonplanar propagator are
$$
-\frac{4528302169622311}{2134332979178400000},\qquad
-\frac{592700482300543}{58541704571750400000}.
$$
For the recurrence, set the auxiliary coordinate to one in the homogeneous cubic and quartic and denote them by $C,D$. Write $A=4y_6y_7$ or $4y_5y_8$ according to the endpoint. Expand $(1-A+zC+z^2D)^{19/2}$ in $z$. Its degree-$r$ coefficient, after factoring $(1-A)^{19/2-r}/(2^r r!)$, is a polynomial $S_r$ with $S_0=1$ and
$$
S_r=(21-2r)C\,S_{r-1}+4(r-1)(21-r)(1-A)D\,S_{r-2}.
$$
The second term is absent at $r=1$. Homogeneity selects the degree-$20+r$ physical monomials after restoring $(1-A)^{19/2-r}$. The two endpoint square-root factors contribute $4^{-19}$. The $69754$ rows require degrees $9\le r\le14$ and a downward closure of $6160953$ exponent vectors. Applying the recurrence to the absolute coefficient one-norm bounds each signed integer below $2^{191}$. Thus three $64$-bit limbs give exact signed coefficients throughout, followed by rational weight summation. The two numbers displayed multiply $20\zeta(5)$ in this cut convention. Four lower ordinary masters can contribute further $\zeta(5)$ to the complete period. [Exact native convolution](quartic_nonplanar_bank.cpp) · [Full-bank replay and independent cut-preserving relabeling](quartic_nonplanar_complete.py)
There is an exact relation between the planar six-line master and the endpoints of other cuts of the same vacuum graph. Keep the six propagators other than $u=y_5,v=y_8$ and set their squared momenta to zero. With $w=-2p_6\!\cdot p_7$, direct substitution gives
$$
16B=w\big[(u-v)^2+(w-1)(u+v-w)^2-w(w-1)\big].
$$
For $w>1$, the variables $(u+v-w)/\sqrt{w}$ and $(u-v)/\sqrt{w(w-1)}$ parametrize a disk. Integrating its radial and angular moments gives the six-line projection at generic dimension. If both remaining denominators are present, move all extra propagator derivatives onto the Gram. Reduce its residual powers using $b=B(0,0,w)=w^2(w-1)^2/16$. In the remaining $f(w)B^s/(uv)$ term solve $f=c+g'+s b'g/b$. The constant $c$ gives the unit eight-line master. The other terms are one-contact disk moments. The second six-line contour interchanges the pairs $\{5,8\}$ and $\{6,7\}$. Summing their lower projections gives the ordinary planar master coefficient. Two independent implementations reproduce both rational coefficients in 66 FORCER checks, including missing lines, numerator insertions and a dimension shift. [Six-line projection](quartic_planar_subcut.py)

Write $m=q-2$ for the number of Gram insertions. Let $T$ denote the unit nonplanar coefficient at $\epsilon=0$ and $A$ the $\epsilon^{-2}$ coefficient of the six-line planar master. For a transferred Taylor numerator with $L=\sum_{i=1}^8(a_i-1)-m\ge0$, a monomial $u^i w^k B^{-L-1/2-\epsilon}/v$ in the normalized disk gives the double pole
$$
\frac{16^L}{2}(-1)^{i+k}\binom{i}{2L}\binom{i-2L}{2L-k}.
$$
Only the angular term with index $i-L$ survives at this order. For the pure $w^k/(uv)$ term the lower double pole is $16^L\binom{k-2L-1}{2L}/4$ when $k\ge4L+1$, and zero otherwise. Indeed the Hermite primitive has at most a simple pole; its residue is a multiple of $w(w-1)$. Division by $w(w-1)$ and extraction of $w^{-1}$ at infinity gives the coefficient, with normalization fixed by $k=1,L=0$. Both projector terms have zero triple pole. [Exact Laurent expansion and 24 independent terminal reductions](quartic_planar_double_pole.py) · [Rational projections](quartic_planar_terminal.json)

For positive propagator powers these disk residues simplify further. In each chart and for each missing edge $a$, substitute $y_a\mapsto y_a-1-\sigma-\tau$, $w=-\tau$. With the Taylor branch $\sqrt{16B}=\sigma\tau$ at zero external Taylor variables, the one-contact binomials extract the coefficient of $\sigma^{-1}\tau^{-1}$ in $(16B)^{m-1/2}/(1+\sigma+\tau)$. After summing the two axes of each chart the pure terms contribute $2T$. Moving the $\sigma$ residue to infinity crosses the simple pole $1+\sigma+\tau=0$. Its two endpoint contributions cancel this $2T$. Restore the squared cut momentum as the ninth homogeneous variable. The homogeneous Gram has degree four and the nine Taylor indices sum to $4m-1$. With $x=1+\sigma+\tau$ the residual coefficient in $(16B(y-xe_a,w))^{m-1/2}$ has simultaneous $(x,w)$ degree $-1$. Consequently the $x^0w^{-1}$ term at large $x$ is the zero-endpoint residue at fixed $x=1$ on the cut of edge $a$. Restoring this edge as the external square $y_a-1$ multiplies the other eight residues by $(y_a-1)^{a_a-1}$; its leading coefficient is one. The branch and change of contour supply the sign in the reduction. Thus $A$ is the sum of **four ordinary single-endpoint nonplanar projections with different cut edges**, each carrying the previously defined unit-master factor $1/8$. In the first chart the new cut edges are 5 and 8, with respective zero-endpoint quadratic pairs $\{2,4\}$ and $\{1,3\}$. The exchanged chart gives the other two edges. Exact polynomial identities verify each graph relabeling and its affine auxiliary-momentum change. [Disk, four-vertex and recut comparison](quartic_planar_vertex.py)

In the ordinary basis let $C$ denote the $\epsilon^{-2}$ coefficient of the next, elementary five-line master. When the four lower Laurent coefficients have pole orders at most $(2,2,1,2)$, cancellation of the rational $\epsilon^{-5}$ and $\zeta(3)\epsilon^{-2}$ terms gives $20T+120A+40C$ for the complete cut's $\zeta(5)$ coefficient. All 24 cuts of the three convergent terminal bipartite periods and nine terminal wheel periods satisfy these pole bounds and this formula. Projective normalization at $q=3$ divides by two. A further finite check gives the same $\zeta(5)$ coefficients by summing $20T$ over all nine deleted edges with projective factorial weight, including deletion of the absent line in the wheel. On the independently integrated $q=12$ orbit $(9,6,5,6,4,5,4,5,4)$ the nine-cut sum gives $1693/21339318$.


#### The five-line gamma master at the double pole

Use the ordinary master basis of the three-loop propagator obtained by cutting $K_{3,3}$. The eight-line integral is $M_0$, and $M_2$ deletes lines $5,6$ in the ordinary basis. The five-line units with missing triples $(1,6,8),(2,5,7),(3,7,8),(4,5,6)$ each reduce to $\epsilon M_2/(3\epsilon-2)$. All other five-line units have zero $M_2$ coefficient. The four displayed sectors are permuted by the stabilizer of the cut edge. On each such contour subtract the unit eight-line contour times its maximal-cut coefficient. The sum of the four remainders times $\epsilon/(3\epsilon-2)$ extracts precisely the ordinary $M_2$ coefficient. The other ordinary basis integrals have no denominator support on these cuts. [Ordinary basis](../nonSP/verify/k33_forcer/replay.py) · [Independent five-cut reconstruction and FORM checks](quartic_five_cut.py)

On the first cut use the positive variables $a=p_1^2$, $b=-p_6^2$, $c=-p_8^2$ and the original auxiliary Baikov coordinate $x$. The cut Gram determinant $B$ is symmetric in $a,b,c$. Put $l=x+b+c$, $b=ly$, $c=lz$. As a quadratic in $a$, $16B$ has leading coefficient $l^2$ and roots
$$
\left(\sqrt{(1-l)(1-y)(1-z)}\ \,\pm\sqrt{yz}\right)^2 .
$$
Integrate between the roots and over $0<l,y,z<1$. For Gram exponent $s=-1/2-\epsilon$ the unit five-cut measure after angular integration is a constant times
$$
l^{2s+2}(1-l)^{s+1/2}\,[y(1-y)z(1-z)]^{s+1/2}\,dl\,dy\,dz .
$$
Write $e=s+3/2=1-\epsilon$. The only surviving cut of the ordinary $M_2$ requires one inverse positive variable. By symmetry its normalized value is $\langle b^{-1}\rangle=\langle l^{-1}\rangle\langle y^{-1}\rangle=[(3e-1)/(2e-1)][(2e-1)/(e-1)]=(3\epsilon-2)/\epsilon$, fixing the five-cut normalization. Integration by parts first holds for sufficiently large real $s$, where the Gram factor vanishes at the boundary. Continuation of the beta integrals defines the dotted moments. Transfer the three free-line propagator derivatives onto the Gram polynomial. Each remaining monomial has the form $a^ib^jc^kx^nB^{s+h}$ with $i,j,k\ge-1$, $n\ge0$ and integer $h$. Normalize by the unit five-cut at $s$. If $i\ge0$ and $h=-L\le0$, its coefficient of $\epsilon^{-3}$ is
$$
-\mathbf 1_{j=k=-1}\,16^L(-1)^n
\binom{i}{2L}\binom{i-2L}{2L-n},
$$
with the binomials zero unless $0\le 2L-n\le i-2L$. The coefficient vanishes for $h>0$. To see this, expand the angular moment of $a^i$ in powers $(yz)^m[(1-l)(1-y)(1-z)]^{i-m}$. At $\epsilon=0$ its coefficient can be nonzero only for $m,i-m\ge L$. The independent $y,z$ beta factors then both have poles only at zero, forcing $j=k=-1$, $m=L$, and the constant term of $(1-y-z)^n$. The angular coefficient becomes $16^L\binom{i}{2L}$. The $l$ integral gives the remaining signed binomial. Every angular coefficient is regular, so the three beta integrals also bound the pole order by three. [Exact Laurent evaluator](quartic_five_pole.py)

For the same moment take $-1/4$ times the sum of the three cyclic residues obtained by first taking the residues at $b=c=0$ and then at $r=x+a-1=0$, $x=0$, on the branch $\sqrt{B}=xr/4$ after the first two residues. The other axes use cyclic permutations of the positive variables. At $b=c=0$ one has $16B=x^2r^2$ and $a=1+r-x$. Binomial expansion gives exactly the preceding open-moment formula. Call this sum $D$ for the remainder of the calculation.

It remains to fix the fully polar monomials with $i=j=k=-1$. The maximal eight-line cut has $B_0=x^2(x-1)^2/16$. At fixed $h$, reduction of $x^n$ there uses a primitive $x(x-1)q(x)$. Its coefficient of $x^j$ in $q$, $0\le j<n$, is
$$
\frac{(j+2+2s+2h)_{n-j-1}}{(j+2+4s+4h)_{n-j}},
\qquad (u)_k=u(u+1)\cdots(u+k-1).
$$
It has at most a simple pole. The terms left off the maximal cut contain $(4x-2)B-x(x-1)B_x$, which vanishes at $a=b=c=0$, so they are open moments. Multiplication by $B=B_0+(B-B_0)$ shifts $h$ by one. In the reverse shifts at $h<0$ the leading maximal-cut multiplier is $(2h)_2^2/[16(4h)_4]\ne0$. Forward shifts use regular coefficients. Thus the subtracted fully polar moment has at most a fourth-order pole.

Its fourth-order coefficient vanishes on open moments. Multiplication by $B$ and integration by parts in $x$ therefore identify this coefficient with a linear functional on rational one-forms in $x$ with poles at $0,1,\infty$, modulo derivatives. Indeed at $\epsilon=0$ the forms are polynomial multiples of $B_0^{h-1/2}dx$ and the derivatives have primitives $x(x-1)q(x)B_0^{h-1/2}$. Allowing all integer $h$ gives every such rational form and primitive. Partial fractions leave two classes, represented by $B_0^{-1/2}dx$ and $xB_0^{-1/2}dx$. At $h=0$ the subtracted lower moments for $n=0,1$ are respectively
$$
0,\qquad\frac{-54\epsilon^3+81\epsilon^2-39\epsilon+6}{4\epsilon^3}.
$$
Both fourth-order coefficients vanish on the two generators and hence on all moments. The normalized maximal-cut coefficient on a fully polar monomial is $16^{-h}(2s+1)_{n+2h}(2s+1)_{2h}/(4s+2)_{n+4h}$, with the negative rising factorial defined by $(u)_{-k}=1/(u-k)_k$. This ratio is regular at $\epsilon=0$ by counting its linear factors. The third-order coefficient minus $D$ now descends to the same quotient, since all fourth-order coefficients have vanished and ordinary residues kill exact forms. On these two classes its values are $0-(-3)=3$ and $3/2-0=3/2$. The normalized maximal cut takes values $1$ and $1/2$. Consequently the third-order coefficient on each five-cut is $D$ plus three times the specialized maximal-cut coefficient, for all dots and nonnegative dimension shifts.

For positive propagator powers write $T$ for the ordinary eight-line coefficient at $\epsilon=0$. Write $R_+$ for the sum of the four normalized double residues on axes $1$ through $4$ with the open denominator coordinate shifted by $1+\sigma+\tau$ and $x=-\tau$. Write $R_-$ for the sum on axes $5$ through $8$ with shift $-1-\sigma-\tau$. In each case $16B=\sigma^2\tau^2$ on setting the denominator-shift variables to zero. Higher propagator powers are Taylor coefficients in these variables. This commutes with the iterated residue: equivalently, transfer the holomorphic derivatives onto $B^{s+h}$ before taking it. Under the four cut relabelings the positive axes occur once and the negative axes twice. The residue at the simple open propagator includes $(1+\sigma+\tau)^{-1}$ for either orientation after passing to the positive variable. Multiplying the four third-order coefficients by the leading $-\epsilon/2$ from the five-line normalization gives
$$
[\epsilon^{-3}]M_2\text{ coefficient}=0,\qquad
[\epsilon^{-2}]M_2\text{ coefficient}=-6T+\frac{R_+}{8}+\frac{R_-}{4}.
$$
The rational-function five-cut and its independent Laurent evaluation reproduce all 24 terminal double poles. The double residues also agree on positive dotted controls at four different dimension shifts. [Eight-node replay](quartic_five_vertex.py)

#### Four-loop weight five on the complete bipartite core

For an edge $c$ deleted to form a propagator integral, denote its normalized zero and one endpoint residues by $E_{c,0},E_{c,1}$, with $T_c=E_{c,0}+E_{c,1}$. Their difference obeys a global constraint. With cut square $\lambda$, the positive branch of $\sqrt{16B}$ at undotted internal lines is $w(\lambda-w)$. The coefficient of $w^{-1}$ at large $w$ equals the sum of the two oriented finite residues, hence $E_{c,0}-E_{c,1}$ after normalization. Restoring all nine squared-momentum variables $y$ makes $16B$ homogeneous of degree four. Write $\alpha_i=a_i-1$. Since $\sum_i\alpha_i=4m-1$, the auxiliary residue after extracting eight Taylor coefficients is homogeneous of degree $\alpha_c$ in the remaining cut square. Its value at $\lambda=1$ equals its coefficient of $\lambda^{\alpha_c}$. Thus the residue at infinity can be computed from the same homogeneous polynomial on every cut.

In the reference chart $c=0$, choose on each re-cut edge the endpoint corresponding to $w=0$ at vanishing other squares. The exact change of the Gram coordinate is $w'=s w+l(y)$ with $l$ linear, $l(e_a)=0$, and $s=\pm1$. To verify it, write $L,L'$ for the cubic-in-$w$ coefficients of the two monic quartics. The cubic term gives $l=(sL-L')/4$; substitution verifies the remaining three powers identically. The signs for the nine edges in order are $(+,-,-,-,-,+,+,+,+)$. Translation preserves the formal residue at infinity, and reversal changes its sign. Therefore $E_{a,0}-E_{a,1}=s(E_{0,0}-E_{0,1})$ in these charts.

For a positive open-axis substitution $y_a\mapsto y_a+1+\sigma+\tau$ the simple pole $1+\sigma+\tau=0$ contributes the old zero endpoint and the contour at infinity contributes the re-cut zero endpoint with positive sign. For the negative substitution the re-cut term has negative sign. Indeed at $x=1+\sigma+\tau\to\infty$ the initial square root is $-wx$ and $d\tau=-dw$. The resulting branch is $wx$ on the scale $y_a\pm x$; the standard local square root has leading term $w' (\pm x)$. The factor from $w'=s w+l$ cancels the oriented differential factor $s$. Consequently one normalized positive residue divided by eight gives $E_{0,0}+E_{a,0}$, and one negative residue gives $E_{0,0}-E_{a,0}$. Summing the four edges of each type eliminates the endpoint differences:
$$
\frac{R_+}{8}=2T_0+\frac12\sum_{\substack{a\ \mathrm{incident}\ \mathrm{to}\ 0\\a\ne 0}}T_a,
\qquad
\frac{R_-}{8}=2T_0-\frac12\sum_{\substack{a\ \mathrm{disjoint}\ \mathrm{from}\ 0}}T_a.
$$
The six-line disk formula above has two disjoint edges from each of the old zero and one charts, so their signed endpoint differences cancel as well. The double poles of the six-line and five-line coefficients in the ordinary basis, respectively $A$ and $C$, are therefore
$$
A=\frac12\sum_{\substack{a\ \mathrm{disjoint}\ \mathrm{from}\ 0}}T_a,
\qquad
C=\frac12\sum_{\substack{a\ \mathrm{incident}\ \mathrm{to}\ 0\\a\ne0}}T_a
-\sum_{\substack{a\ \mathrm{disjoint}\ \mathrm{from}\ 0}}T_a .
$$

The remaining two coefficients have orders at most one and two. Otherwise, first consider the deepest pole below $\epsilon^{-5}$. Only the last two masters could contribute; their leading residues are $2\zeta(3)$ and $1/6$. Rational independence of $1$ and $\zeta(3)$ and convergence of the complete input exclude such a pole. This bounds their coefficient orders by four and two. At a pole $\epsilon^{-k}$ with $k\ge3$, the three masters with rational leading residue have only rational contributions since their $\zeta(3)$ terms start at order zero. The leading $2\zeta(3)/\epsilon$ of the remaining master thus bounds its coefficient order by one. This supplies all four bounds for the weight-five extraction above. In the standard multiple-zeta expansion of every convergent positive-power bipartite graph,
$$
[\zeta(5)]\,P_{\rm cut}=20T_0+120A+40C
=20\sum_{\substack{a\ \mathrm{edge}\ \mathrm{of}\ K_{3,3}}}T_a .
$$
The equality holds before the projective factorial. Direct ordinary reductions give this identity on all $35$ convergent bipartite orbits at $q=4$, the three terminal bipartite orbits at $q=3$, and the independent physical $q=12$ control. The nine wheel terminal values also pass the finite test upon including deletion of the absent line. [Affine transformations, disk projections and endpoint signs](quartic_planar_vertex.py) · [Independent finite reductions and complete bank replay](quartic_recut_bank.py)

The bipartite contributions to the coefficients multiplying $\zeta(5)$ in $K_4,J_4$, including the projective weights, are respectively
$$
\frac{249372365272156031}{806626220400000},\qquad
\frac{2209968587039007449}{1393850108851200000}.
$$
For comparison, the ordinary six-line double poles in the previous nonplanar reporting convention are respectively
$$
-\frac{169378632737740819356497}{149403308542488000000},\qquad
-\frac{746978483911572523339}{241054077648384000000}.
$$
On orbits admitting a certified zero top coefficient, the convention fixes the highest-power such cut, then its lexicographically least uncut tuple. Otherwise it uses the highest-power cut with the same tie-break. Of the eighteen single endpoints on each of the $465564$ orbits, the exact weight inequalities leave $3383833$ distinct retained representative rows after combining stabilizer repetitions. The recurrence through homogeneous degree nineteen acts on a downward closure of $42624703$ cells. Its absolute coefficient one-norm bounds every intermediate by $2^{202}$. The signed calculation uses four 64-bit limbs and exact rational final weights. A second run uses the largest endpoint representative on every orbit and both endpoints of every disjoint edge. Its downward closure has $42939230$ cells. Both nine-cut sums and both planar double poles agree exactly. [Bounded-memory native residue convolution](quartic_recut_bank.cpp)

#### Positive angular expansion for the bipartite period

The bipartite period also obeys exact relations at fixed dimension, including its four-wheel boundary. In the two-dimensional cycle frame use the three columns $(1,0),(0,1),(-1,-1)$. The four-dimensional cycle columns of $K_{3,3}$ are the tensor products of the row and column vectors, up to an edge sign. This gives integral flow coordinates: the first two rows and first two columns freely prescribe the four currents, and the remaining currents follow from the zero row and column sums. The cycle Gram determinant is consequently the cotree polynomial $U$. Acting on the row frame by the infinitesimal matrix $\operatorname{diag}(1/2,-1/2)$ changes the three Schwinger parameters in each column by $(t_{1j}+t_{3j},-t_{2j}-t_{3j},0)$ in the Gram matrix. The determinant and the integration volume are invariant. Fix any third-row parameter to one as projective scale. Integration by parts on the positive orthant gives a diagonal coefficient $\sum_{j\,\mathrm{columns}}(a_{1j}-a_{2j})$ on the original normalized cut and coefficients $+a_{3j},-a_{3j}$ on the shifts transferring one power from row one or two to row three. This uses the original normalization $\Gamma(q)/\prod_e\Gamma(a_e)$. At $a_{ij}=1$, the derivative of the normalized half-line power is the delta distribution at zero, with unit weight. The zero-power shifts are therefore the normalized contracted wheel integrals.

For a boundary proof replace every $a_e$ by $a_e+z$ and $q$ by $q+9z/4$. Every proper positive-cycle convergence margin increases since $4|S|>9h(S)$ on proper edge sets of $K_{3,3}$. At positive $z$ all original endpoint derivatives and shifted integrals are ordinary, and their limit on a zero-power edge gives the delta just described. The four-loop bank has $695485$ strictly positive convergent orbits and $3444054$ nontrivial equations. Every shifted term entering these equations individually passes the convergence inequalities before collecting any graph symmetry. [Row generator](quartic_bipartite_orbit_ibp.cpp) · [All individual boundary checks](quartic_bipartite_regulator.cpp)

With the wheel values set to zero, every integer $q\ge3$ admits an explicit homogeneous solution of the individually convergent identities. Write $s_v$ for the sum of the three powers in row or column $v$. Its support consists of matrices having an entry $a_{ij}=q-1$ with $s_i=s_j=2q-1$, where $i$ labels its row and $j$ its column. Two rows cannot attain this sum since the remaining positive row contributes at least three. The same holds for columns. At the supported entry the value is
$$
\frac{(q-3)!((q-2)!)^3}{((2q-5)!)^2}
\frac{\prod_{v\notin\{i,j\}}(2q-2-s_v)!}
{\prod_{e\ {\rm edges}}(a_e-1)!}.
$$
The first factor makes the value one on the matrix with rows $(1,1,1),(1,q-2,q-1),(1,q-1,q-1)$.

Multiply each unknown by $\prod_e(a_e-1)!$. A unit transfer from row $1$ to row $3$ now has coefficient $a_{1j}-1$ in column $j$ and vanishes at a contracted edge. Row $3$ cannot have a saturated entry before the transfer since shifting into it violates $a_{3j}<q$. For row $3$ to acquire a saturated entry through an allowed identity, its old row sum would be $2q-2$, its receiving entry $q-2$, and the column sum $2q-1$. The complementary two-by-two block would then have total $q+1$. A transfer from either donor into row $3$ at another column would lower that total to $q$, violating strict convergence. Thus a term survives only when the source was already supported, with the distinguished row equal to a donor. For distinguished row $1$, only shifts from row $2$ preserve it. They multiply the factorial product in the numerator by $(2q-1-s_2)/(2q-2-s_3)$. Their coefficients sum to $s_2-3=2q-2-s_3$. If $s_2=3$, every entry in row $2$ equals one and the distinguished column forces $a_{3j}=q-1$, excluding this orientation. The surviving sum therefore cancels the diagonal $s_1-s_2=2q-1-s_2$. Interchanging donors and transposing the matrix give the remaining identities. At $q=12$ the displayed solution has exactly $3410$ nonzero orbits. The rank with known wheel boundary is $695484$ on $695485$ variables, and all $465564$ physical orbits lie outside this support. [Exact support and physical weights](quartic_bipartite_physical.py)

**Arbitrary twelve-wave bare graph.** Fix the fractional positions on each subdivided $K_{3,3}$ core path. Let $m_e$ count the twelve insertions on the core edges of resistances $T_e$. Write $c$ for the projective heat cofactor and $\Phi$ for a transverse relative voltage. The projective numerator at denominator $U^{14}$ is
$$
\left(\prod_{e\ {\rm core\ edges}}T_e^{m_e}\right)(U^4 c)^2(U\Phi).
$$
Here $U^4c$ is the degree-fifteen current determinant and $U\Phi$ a degree-five two-forest polynomial. Currents are bounded uniformly in the resistances by the spanning-tree flow formula. For a proper set $S$ of core paths of cycle rank $h>0$, its $5h$ cycle-current rows are supported on $|S|+\sum_{e\in S}m_e$ segments and annihilate their time vector. A nonzero cofactor requires at least $5h+1$ segments there. Scale all paths of $S$ by $\epsilon$. The order of $U$ is $h$, a segment outside $S$ can be deleted to bound $c$, and all relative voltages remain bounded. Every positive-power moment in this numerator therefore satisfies
$$
\sum_{e\in S}a_e\ \ge\ |S|+\sum_{e\in S}m_e+9h\ \ge\ 14h+1.
$$

Distinguish a core edge and multiply its resistance by $R$. Divide the four opposite-square resistances by $R$, keeping its four adjacent resistances fixed. Every cotree uses an opposite edge. Some cotree uses the distinguished edge and exactly one opposite edge. Thus $U$ has Laurent degree zero at $R=\infty$. Let $m$ and $b$ count insertions on the distinguished edge and the opposite square. Removing the distinguished path leaves three cycles on $20-m$ segments; the square carries one cycle on $4+b$ segments. The preceding rank bound gives $m\le4$ and $b\ge2$. Deleting a longest segment of the distinguished path bounds $|c|$ by a constant divided by $R$. Every voltage difference grows at most linearly in $R$. For $m\le3$, the numerator has Laurent degree at most $m-b-1\le0$.

For $m=4$, project the cycle basis onto three cycles avoiding the distinguished edge. Take the opposite square and two further cycles, one using each pair of paths incident to an endpoint of the distinguished edge and closed by paths in the square. On the sixteen surviving segments the latter two pairs give disjoint column combinations. Weight the segments in each pair by their times. Each combination in the fifteen projected rows is a voltage difference across the small square, bounded by a constant divided by $R$; all other projected cycles have zero entries on that pair. At $R=\infty$ the projected rank is consequently at most fourteen. The other five rows can bring the full minor only to rank nineteen. All currents are rational and bounded at infinity. The deleted twenty-row minor therefore has a zero of order at least one in $1/R$, and dividing by the long segment gives $|c|$ bounded by a constant divided by $R^2$. The numerator degree is at most $m-b-3\le-1$.

These bounds hold at arbitrary positive core times and at every interior fractional position. Polynomiality then enforces them separately on each monomial, before and after the exact chamber integrals. As the powers $a_e$ include one from each logarithmic measure, this gives
$$
a_{\rm distinguished}-\sum_{e\in\mathrm{opposite\ square}}a_e\le-3 .
$$
The homogeneous solution above requires $a_{\rm distinguished}=q-1$ and opposite sum $q+1$. Every bare twelve-wave numerator therefore annihilates this direction for arbitrary external momenta and transverse polarization. At $q=14$ completeness of the relative integration-by-parts kernel is open. In three original twenty-one-edge tests over $\mathbb F_{1009}$ the numerator degrees are $-1,0,-5$. At fixed fractions and a root of $U$ the numerators give $404,594,739$, respectively, verifying that the unrestricted fixed-position integrands can require $U^{14}$. [Original five-component determinant and independent polynomial checks](quartic_bipartite_unrestricted.py)

For a strictly convergent orbit put $b_e=q-a_e$. Inverting the nine Schwinger parameters exchanges the four-edge cotree polynomial with the five-edge tree polynomial and replaces $a_e$ by $b_e$. In the six-vertex position graph fix two vertices of the same bipartition at $0,e$, where $e^2=1$. Gaussian integration of the other four positions gives the two-root forest determinant $F$ and exponential $\exp(-\widetilde U/F)$, with $\widetilde U$ the tree polynomial. The scale integral is $\Gamma(q)(F/\widetilde U)^q$, so the cut integral is exactly
$$
P_{\rm cut}=\prod_{e'\in E}\frac{\Gamma(b_{e'})}{\Gamma(a_{e'})}
\int\frac{d^{2q}x}{\pi^q}\prod_{j=1}^3
\int\frac{d^{2q}y_j}{\pi^q}\,
(y_j^2)^{-b_{0j}}\big((y_j-e)^2\big)^{-b_{1j}}
\big((y_j-x)^2\big)^{-b_{2j}}.
$$
For one integrated star with incident Schwinger powers $a_1,a_2,a_3$, suppose $h=q-a_1-a_2-a_3\ge0$. The corresponding position powers are $b_i=q-a_i$. The radial Schwinger integration on that star leaves $(t_1+t_2+t_3)^h$ divided by the distance quadratic to power $q+h$. Expand the numerator with $k_1+k_2+k_3=h$ and invert the three parameters. The quadratic becomes a linear form divided by their product. The three remaining logarithmic-measure exponents are $a_i+h-k_i>0$ and sum to $q+h$. The Dirichlet integral then gives a prism position graph. Its equator edge opposite neighbor $i$ has power $a_i+h-k_i$ and its other six powers are $q-a_e$. For each composition the multiplier, including the Fourier factors above, is
$$
\frac{h!}{k_1!k_2!k_3!}
\prod_{i\in\mathrm{star}}\frac{\Gamma(a_i+h-k_i)}{\Gamma(a_i)}
\prod_{e\notin\mathrm{star}}\frac{\Gamma(q-a_e)}{\Gamma(a_e)}.
$$
Each term is positive before period reduction, so strict convergence of the original integral gives convergence of every prism term. At $q=12$ this reduces $572690$ of the $695485$ positive orbits. At $h=0$, all $237400$ star choices give $216725$ distinct interior orbits, with agreement across choices. The exact interior solution also agrees with the expanded prism sum at every $h=0,\ldots,9$ at $q=12$ and $h=0,\ldots,3$ at $q=6$. [Star expansion and exact checks](quartic_bipartite_stars.py)

The three rows of $b$ refer to the fixed positions $0,e$ and to $x$ and the columns refer to the opposite bipartition. Writing $\lambda=q-1$, the degree-$n$ coefficient of a propagator of power $b$ is a finite sum over $0\le k\le\lambda-b$. Its term is
$$
\frac{(b)_n(b+n)_k(b-\lambda)_k}{(\lambda)_n(n+\lambda+1)_k k!}
\frac{\min(r,\rho)^{n+2k}}{\max(r,\rho)^{n+2k+2b}}.
$$
The full coefficient is positive: its terminating hypergeometric polynomial is the beta mean of $(1-tX)^{\lambda-b}$, with $t=(\min/\max)^2$ and beta parameters $n+b,\lambda+1-b$. Orthogonality in the direction of $y_j$ gives $\lambda/(n+\lambda)$ times the remaining Gegenbauer polynomial of the angle between $x$ and $e$.

For $0<r=|x|<1$, split the radial integral of $y$ at $r,1$. With $d=2q-2b_0$, $u=n+2k_1$, $v=n+2k_2$ its three pieces sum to
$$
\frac{r^{d+u-2b_2}}{d+u+v}
+r^v\frac{1-r^{d+u-v-2b_2}}{d+u-v-2b_2}
+\frac{r^v}{u+v+2b_1+2b_2-d}.
$$
When the middle denominator vanishes its term is $-r^v\log r$. Inversion covers $r>1$ by $b_0\mapsto2q-b_0-b_1-b_2$ and an overall $r^{-2b_2}$. In the outer angular integral the three degrees are $j+k,k+i,i+j$ for nonnegative integers $i,j,k$. With $s=i+j+k$, the mean of the three Gegenbauer polynomials is
$$
\frac{\lambda(\lambda-1)!}{(2\lambda-1)!}
\prod_{d=\lambda+1}^{2\lambda-1}(s+d)
\prod_{t\in(i,j,k)}\binom{\lambda-1+t}{t}.
$$
The last radial integral uses $\int_0^1 r^{p-1}(\log r)^h\,dr=(-1)^hh!/p^{h+1}$. At distinct radii the angular expansion converges absolutely. After angular integration each mode has positive radial integrand, so the remaining sums and integrals commute by monotone convergence. The full result is a positive rational triple series. Its mode denominators have only the forms $j+k+c$, $k+i+c$, $i+j+c$ and $i+j+k+c$, for integer shifts $c$. On all three terminal bipartite graphs the shells $s\le6$ lie strictly below the independently integrated periods in either choice of bipartition. [Exact radial coefficients and shell replay](quartic_bipartite_angular.py)


At fixed shell $s=i+j+k$, the connection coefficient with the radial ratio removed is
$$
\frac{(-1)^k\binom{\lambda-b}{k}(\lambda-1)!}{(b-1)!}
\left(\prod_{t=b+k}^{\lambda-1}(n+t)
\prod_{t=\lambda+1}^{\lambda+k}(n+t)\right)^{-1}.
$$
Its two factors have disjoint simple poles. Each triangle mode times $r^{-n}$ therefore has a finite partial-fraction expansion in $r^{2l}(\log r)^h/(n+a)^p$, with $h\le1$, $1\le p\le3$. For a product of three such terms absorb each directional binomial into the polynomial numerator, ordering the three angular indices as $s-i,s-j,s-k$. The three shifted denominator variables sum to $2s+a+b+c$. Division by this sum and decrement of each power reduces to at most two positive denominator powers. A nonpositive power is removed by the binomial theorem. For two denominators with powers $p,q>0$ and shifts $a,b$, the remaining domain is the square $0\le u,v\le s$ with the lower triangle $u+v<s$ excluded. The square gives a product of finite harmonic sums. On the excluded diagonal $u+v=t$, partial fractions give
$$
\sum_{h=1}^{p}\binom{p+q-h-1}{q-1}
\frac{\sum_{u=0}^t(u+a)^{-h}}{(t+a+b)^{p+q-h}}
+\sum_{h=1}^{q}\binom{p+q-h-1}{p-1}
\frac{\sum_{v=0}^t(v+b)^{-h}}{(t+a+b)^{p+q-h}} .
$$
Summing these diagonals up to $s-1$ expresses every complete shell in rational functions of $s$ times strict harmonic sums $H_s(w)$ of depth at most two. Integer shifts telescope by $H_s(h,w)-H_{s-1}(h,w)=s^{-h}H_{s-1}(w)$. At every tested $q=3,4$ orbit, the apparent half-integer poles cancel in the complete shell. Its outer sum thus has a rational harmonic primitive. At infinity, the same difference equation recursively gives its constant coefficient in powers of $s^{-1}$ and $\log s+\gamma$. The graph is convergent, so this constant gives the integral after subtracting the lower endpoint. Ordinary shuffle, finite-sum stuffle, reversal duality and the multiple-zeta sum formula reduce the result exactly.

Independent cuts with FORM reproduce all three resulting coefficients on the $35$ convergent $q=4$ orbits. Each of the three $q=3$ bipartite periods also agrees in both choices of bipartition. For example, the $q=4$ orbit $(1,2,2,2,2,1,2,2,2)$ gives $1/20+29\zeta(3)/30-7\zeta(5)/6$ before its projective factorial. The fixed-dimension integration-by-parts calculation below evaluates all three columns at $q=12$. [Triangle partial fractions](quartic_bipartite_sums.py) · [Harmonic shell algebra](quartic_harmonic_shell.py) · [Exact period and 35-orbit comparison](quartic_harmonic_period.py)

For a convergent four-wheel boundary period at dimension $2q$, the common rim harmonic degree $n$ has an exact reflection. Put $\lambda=q-1$. In a chamber of ordered radii let $S$ be a proper initial set of vertices, and $b_e$ the power on a crossing rim edge. The radial denominator is half the number of crossings times $n$, plus its zero-mode shift $h_0(S)$ and the sum of the crossing connection indices $k_e$. The identity $h_0(S)+h_0(S^c)=\sum_{e\,\mathrm{crossing}}b_e$ follows by summing the radial powers; both shifts are positive by convergence. Reverse the chamber and send $k_e$ to $\lambda-b_e-k_e$ and $n$ to $-2\lambda-n$. Each of its three radial denominators changes sign. The connection in each rim kernel has roots $b_e+k_e,\ldots,\lambda+k_e$ with $\lambda$ omitted; its sign and the signs of those factors cancel under reflection. The common angular numerator $\prod_{a=1}^{2\lambda-1}(n+a)/(n+\lambda)^3$ is even in $n+\lambda$. Thus the full summand is odd in $n+\lambda$.

For opposite first vertices, adding their two radial orders replaces $(1/x+1/y)/(x+y)$ by $1/(xy)$. Every remaining denominator then has two crossings. Its shift is at most $b_i+b_j-1+(\lambda-b_i)+(\lambda-b_j)=2\lambda-1$. The full summand has pole order at most five at $n=-\lambda$ and six at every other integer pole. The harmonic series for the wheel converges: at unequal radii its angular kernel series converges absolutely, and its degree-$n$ radial coefficient is nonnegative by the Euler integral for the terminating hypergeometric connection (its integrand contains $(1-r^2 u)^{\lambda-b_e}$ for radius ratio $0<r<1$). Integration of the nonnegative terms commutes with summation. Consequently the odd rational summand decays at least cubically. Multiplication by
$$
(n+\lambda)^5\prod_{s=1}^{\lambda-1}\big((n+\lambda)^2-s^2\big)^6
$$
gives a polynomial in $(n+\lambda)^2$ of degree at most $6\lambda-5$. At $q=12$, **62 angular evaluations determine all three wheel coefficients**. Paired partial fractions cancel every even zeta column.

There is also a common denominator for the exact reconstruction. For an individual connection write $l=\lambda-b$ and $a=b+k+j\ne\lambda$. Its simple-pole coefficient is the integer $(-1)^{k+j}(\lambda-a)\binom{l}{k}\binom{l}{j}\binom{\lambda-1}{b-1}$. At $b=\lambda$ the connection is one. Expanding all four connections in simple poles gives at most nine linear denominator factors including the angular double pole and the three radial factors. If at least one connection supplies a pole, at least one factor cancels against the integer angular numerator. At most eight denominator factors remain. Splitting distinct poles introduces one inverse integer difference at each step; a coefficient at order $h$ therefore has denominator dividing $\operatorname{lcm}(1,\ldots,2\lambda-2)^{8-h}$ before the common angular normalization. Its finite harmonic correction costs at most $h$ additional powers. A common denominator for all wheel coefficients is
$$
(\lambda!)^3(2\lambda-1)!\,\operatorname{lcm}(1,\ldots,2\lambda-2)^8 .
$$
On circles of radius $1/2$ about the poles, the remaining integer angular polynomial has magnitude at most $(2\lambda-1)^{2\lambda-2}$ and each remaining inverse factor has magnitude at most two. Summing the absolute integer connection coefficients, the twenty radial chambers and the at most eight principal-part coefficients bounds every numerator with that denominator by $2^{429}$ at $q=12$. The direct 62-node bank and the longer interpolation agree on all eight accessed $q=3$ faces and 77 faces each at $q=6,12$ modulo $999999937$. Fifteen primes now determine every coefficient on all $531038$ wheel faces. The maximum reduced numerator has 103 bits. In column order $1,\zeta(3),\zeta(5)$ the observed least common denominators are $9420623616932216131523502617395200000000000$, $588678304537010550236774400000000$, and $9269940547634960793600000$. Five dimension-24 faces also agree with direct rational partial fractions. The reconstructed four-column table has SHA-256 `7693fc7f7cfbc4718066b2b656d456b43cb4e8ab0084fa537d6a6f9df0404597`. [Odd interpolation and integer bounds](quartic_wheel_odd.py) · [Native radial evaluation](quartic_wheel_bank.cpp) · [Full exact reconstruction](quartic_wheel_exact.py)

To derive the two periods, compress the six planar source columns of the current determinant. If the $i$th source occupies edge $e$, its column is the cycle incidence column $C_e$ tensored with $u_i$. Put $L=\sum_e t_e C_e C_e^{\mathsf T}$, so $U=\det L$. The squared determinant summed over placements with measure $\prod_i t_{e_i}$ is $U^2Z$ by Cauchy–Binet. For a fixed pair $i,j$, let $Z_{\widehat i\widehat j}$ be the same coefficient for the remaining $2n-2$ waves at loop number $n-1$. Put
$$
B_{ij}=\frac{Z-n^2[u_i,u_j]^2 Z_{\widehat i\widehat j}}{n(n-1)}.
$$
To keep the two assigned edges fixed, take the two rank-one variations of $\det(L\otimes\sum_{h\ne i,j}x_hu_hu_h^{\mathsf T})$. The product of the diagonal inverse entries has coefficient $B_{ij}+[u_i,u_j]^2 Z_{\widehat i\widehat j}$, and the squared cross entry has coefficient $-B_{ij}$. Their difference follows from the two-dimensional adjugate determinant identity. Summing all edges gives $Z$ since the cycle projection has trace $n$. In particular every single planar label has the same edge marginal $\mu$ as before.

Besides two independent draws from $\mu$, the cross entry gives an ordered pair of edges $(e,f)$ with probability $t_e t_f(C_e^{\mathsf T}L^{-1}C_f)^2/n$. Put one independent uniform position on each chosen edge, including when the same edge is chosen twice. The two marginals are $\mu$ by the projection identity. The second period is the integral of $n/2$ times the difference between this mean resistance and $R_{\mu\mu}$, multiplied by the squared threefold cofactor measure and the four retained label count. The calculation uses
$$
U^2(C_e^{\mathsf T}L^{-1}C_f)^2
=(\partial_e U)(\partial_f U)-U\,\partial_e\partial_f U,
$$
where $\partial_e$ differentiates the edge resistance. On the two points of the same edge of length $t_e$, linear endpoint interpolation has the integrated Green bridge correction $-t_e/6$ in the pair resistance. The other interpolation terms cancel because the marginals agree. The core denominator is $U^{2n+6}$. At $n=2$ this second period is $-1/20736$ and the complete polarized result agrees with direct eight-leg integration on two independent planar fixtures. [Literal placement determinants and both edge marginals](planar_moments.py)

In a representative with zero potential at each retained slot, the residual at planar slot $i$, modulo its own momentum, is the sum of $K Z\sum_j u_j/108$ and $J\sum_{j\ne i}B_{ij}u_j$ at $n=3$. Fusing two retained slots gives zero in either trace because the remaining momenta span at most four dimensions. For each planar-retained fusion the surviving five planar rows determine the trace: a nonzero five-minor uses two of those rows and three retained rows, or one planar row and all four retained rows. Its square is respectively a squared planar bracket or the squared bracket with the sum of the five planar momenta. This gives the two polynomial boundary identities checked in the replay. The first two independent planar rows can be taken as the coordinate basis; the eight remaining coordinates are polynomial indeterminates. Relabeling then verifies every slot, and polynomiality covers degenerate momenta. Any five independent singletons in a hexagon must use at least three retained positions. The other five waves therefore span at most three dimensions, which kills the genus-two bulk insertion throughout the family.

At every positive genus, the quadratic-soft primary trace is uniquely determined by its transverse cubic variation. The rank is $2398$ at genus four and $1005099$ at genus five. This determines uniqueness of the primary counterterm at those genera, conditional on cancellation of the bulk source. [Joint soft symbols and the length-polynomial proof](operator/README.md#injectivity-of-the-transverse-trace-variation-at-every-positive-genus)

On unrestricted momenta at four loops the constant-bivector Ward identity has rank $77199$ on $213877$ polynomial potential classes. This leaves $136678$ reciprocal classes. At a generic conserved thirteen-wave configuration, the trace on its $2431$ double-fusion partitions has image dimension $1918$. The local consistency sequence is exact, with $1918+780+656=3354$ potential coordinates. These trace values glue on any finite tree of generic momentum splits and on the overlapping thirteen-wave configurations obtained from eleven first fusions with a common anchor among fourteen generic waves. All thirteen anchor fusions leave at most four additional compatibility conditions. [Coupled overlap replay](operator/fusion_overlap.py) A global polynomial trace depending only on the eleven final momenta remains open. [Global Ward rank](operator/README.md#exact-global-rank-of-the-constant-bivector-test) · [Exact fusion matrices and overlap tests](operator/README.md#exactness-of-the-four-loop-fusion-equations-on-a-generic-thirteen-wave-configuration)

### The complete genus-four equation on eight planar momenta

Take eight planar momenta $u_i$, three independent transverse basis covectors $e^3,e^4,e^5$, and their common balancing wave. Every known lower-genus primary interaction drops out of the transverse genus-four source. Indeed, at geometric loop number $\ell>0$ every ordinary planar wave enters a bivalent cubic. Subtract the two incident current columns along its path. The resulting column lies in the $2\ell$-dimensional tensor product of cycle coordinates and the momentum plane. These differences are independent column operations along open chains. On a one-loop polygon the distinguished quantum vertex supplies the cut. Thus a nonzero heat cofactor permits at most $2\ell$ planar waves at ordinary cubics. A genus-$h$ primary trace vertex has $2h+3$ conserved arguments. If $2h$ of its external waves are planar, at most three arguments remain. Their projections modulo the plane sum to zero, so all arguments together span at most four dimensions. Its five-row minors vanish. With $m>0$ positive-genus trace vertices and total genus $n$, the number of possible planar waves is therefore at most $2\ell+\sum_{j\in\mathrm{quantum\ vertices}}(2h_j-1)=2n-m$. The tree sewing of two positive-genus vertices obeys the same bound. At $n=4$ this excludes the genus-three hexagon, the genus-two theta insertion and the tree sewing of the two genus-two vertices. The chosen genus-one vertices have external descendants. The free variation vanishes on transverse primaries.

For $n\ge2$, take $n$ copies of $e^1$, $n-1$ copies of $e^2$ and the single planar momentum $-n e^1-(n-1)e^2$. Their sum is zero. At one $e^2$ slot insert $e_1$, at $e^4$ insert $e_3$, and at $e^3$ insert the contraction of the unit top polyvector with $e^3$. The other slots are scalars. The $K$ term vanishes. Every nonzero matching pairs each $e^1$ with a wave outside its line. The $n!$ matchings have common squared bracket product $(n-1)^2$, giving $Z=(n!)^2(n-1)^2$. For the chosen $e^2$ slot the transfer $B_{ij}$ vanishes against each $e^1$ slot and equals $Z/[n(n-1)]$ against the distinguished planar wave. Thus the transverse bulk pairing is $-(n!)^2(n-1)J_n$. At four loops this is the full known source before the new tree counterterm. Integrating the three series–parallel shapes and using strict negativity on the other 47 gives
$$
-1728J_4>\frac{9}{33806080000}.
$$

A local eleven-primary trace cancels the entire eight-planar source. Following the convention $\Phi=C$, define $P_{3^4,2,1^6}$ and $P_{4,3^3,1^7}$ by summing unordered products of four squared five-row determinants whose eleven label multiplicities are respectively the displayed partitions. The bulk action uses the negative of the following trace:
$$
\frac{J_4}{27}P_{3^4,2,1^6}
+\left(\frac{K_4}{13824}-\frac{J_4}{216}\right)P_{4,3^3,1^7}.
$$
To check its boundary, fusing two retained slots leaves rank below five. In every other fusion seven waves remain in the plane and four project to a conserved transverse simplex. A minor with two planar waves reduces up to sign to their bracket. A minor with one planar wave reduces to its bracket with the sum of all seven. In the first multiplicity pattern, one planar wave occurs in two pairs and all other waves occur in one each. The four missing transverse slots can be assigned in $4!$ ways. In the second pattern, one minor uses one planar wave and the other three use disjoint planar pairs. There are again $4!$ assignments. These restrictions are the first and second scalar symbols in the replay multiplied by 24.

For the transverse component at $u_1=e^1$ put $F=Z\sum_{i\in\mathrm{planar\ slots}}(u_i)_2$ and $G=12\sum_{j\ne1,\ j\in\mathrm{planar\ slots}}B_{1j}(u_j)_2$. After subtraction of the common retained potential the two global trace boundaries are $3F/4+9G/4$ and $6F$. The displayed linear combination gives $K_4F/2304+J_4G/12$ as required. These are polynomial identities over the integers after clearing factors. A second independent planar slot can be fixed to $e^2$ by change of frame. The replay expands all twelve remaining coordinates and checks both identities coefficient by coefficient. Symmetry in the eight slots and polynomial continuation give the result on every planar configuration, in every transverse polarization.

The restricted scalar primitive is ambiguous. Its momentum degree on seven planar slots is sixteen, with at least two per slot. The only individual degree patterns are $(4,2^6)$ and $(3,3,2^5)$. For the first pattern the fourth symmetric power of the plane occurs once in the sixth symmetric power of its quadratic forms. For the second, the symmetric square of cubics contains the sixth and second symmetric powers, each once in the fifth symmetric power of quadratic forms. This bounds the dimension by three. Three independent bracket polynomials attain it. Their transverse boundaries have rank two, and the replay checks the kernel identity on the same twelve variables. Consequently the planar calculation fixes two trace combinations. Both $K_4,J_4$ have complete convergent series with certified tails. Their full exact period reductions are given in the elementary-scalar calculation below. The unrestricted twelve-wave diagram reconstruction remains open. The worldsheet primary vertex above solves the effective equation. [Exact planar boundaries and lower-vertex exclusion](operator/four_loop_planar.py)

### A four-loop measurement on the elementary scalar

Use the original point scalar and take its absorbed momentum to be $e^4$. Put the vector wave at $e^3$ with polarization $e_4$. For each list of seven planar scalar momenta below, append scalar waves at $e^5$ and at minus the sum of all the displayed waves and $e^3,e^4,e^5$. In this order let the four entries of the table specify four measurements of the reduced $Q_4\phi(0)$. The divergence of the one-edge tree supplies $e^4(e_4)=1$. All geometric loops at a one-letter mark have a bridge. Thus each measurement reads exactly the eleven-primary trace at its respective configuration whenever the bulk completion exists.

Put $u=e^1$, $v=e^2$, $w=-4u-3v$. Superscripts in parentheses below give the number of copies. The last three columns evaluate the three scalar bracket symbols in the preceding rank test.

| Planar scalar waves | First | Second | Third |
|---|---:|---:|---:|
| $u^{(3)},v^{(2)},w,u+v$ | 3258 | 0 | −1134 |
| $u^{(3)},v^{(3)},w$ | 3492 | 216 | −1728 |
| $u^{(4)},v^{(2)},w+v$ | 288 | 0 | 0 |
| $u^{(4)},v^{(3)}$ | 36 | 216 | 0 |

The last row reduces the whole seven-planar eleven-field trace to one fourfold worldsheet integral. In the following equation $\left.V_{4;0}\right|_{\mathbb C^4}$ is evaluated at three copies of $e^1$, one copy of $e^2,e^3,e^4$ each, and their balancing seventh momentum, in the quotient by the repeated fivefold direction $u$. The polarized fivefold action on the last row is $24\left.V_{4;0}\right|_{\mathbb C^4}$. Its elementary response has the opposite sign. Solving that row gives the following coefficient of $22A-5B+18D$ in the elementary response:
$$
\frac{3K_4/8+8J_4+24\left.V_{4;0}\right|_{\mathbb C^4}}{288}.
$$
Positivity bounds the new period. On $3u,3v,t(u+v)$ divide the elementary response by $18t^2$ and use $y=t-1$. Its quadratic in $y$ has symmetric matrix
$$
\begin{pmatrix}
11K_4/576+J_4+\left.V_{4;0}\right|_{\mathbb C^4} & K_4/72\\
K_4/72 & K_4/72-8\left.V_{4;0}\right|_{\mathbb C^4}/3
\end{pmatrix}.
$$
The matrix is negative semidefinite by the real Betti-form calculation above. Since $K_4<0$, its off-diagonal entry is nonzero, and both diagonal entries are strictly negative. In particular $0\le\left.V_{4;0}\right|_{\mathbb C^4}<-11K_4/576-J_4$. The determinant gives a further quadratic bound. This reduces the undetermined coefficient to a single scalar integral on $\mathcal M_{4,7}$. [Independent bracket evaluation](worldsheet_primary.py)

More sharply, write $V(a,b,c)$ for the fourfold integral on the three collinear slots above and the unit transverse simplex, and put $W=V(1,1,1)$. The bracket restriction gives
$$
V(1,1,-2)=-\frac{11K_4}{32}-18J_4-22W,\qquad
0<W<-\frac{K_4}{64}-\frac{9J_4}{11}.
$$
Indeed, for $abc\ne0$,
$$
\frac{V(a,b,c)}{(abc)^2}=
\frac{(a+b+c)^2}{9}W+
\frac{a^2+b^2+c^2-(a+b+c)^2/3}{24}V(1,1,-2).
$$
Both values are positive by the fourfold rank test. There is also a direct check independent of Ax–Schanuel. On the smooth model of $y^3=x^5-1$, the five points with $y=0$, three with $x=0$, and infinity have torsion differences from infinity, by the divisors of $x-\alpha$ and $y-\beta$. At the seven marks specified in the replay, logarithmic differentials with residue $+1$ at a finite mark and $-1$ at infinity are $d\log(x-\alpha)/3$ and $d\log(y-\beta)/5$. For each of the four integral source rows, fifteen times its degree-zero divisor is principal, with a function $f_j$. The derivative relative to the flat torsion section is dual, up to a common nonzero factor, to multiplication by $d\log(f_j)/15$ from four copies of $H^0(K)$ into $H^0(K^2(\sum_i z_i))$. Indeed a gluing vector field $\xi$ varies the transitions of the trivializing function by $\xi(f_j)/f_j$. Both vector spaces have dimension sixteen.

If the volume vanished, the Betti determinant would vanish identically. Complexify the local real-analytic expression $dw-(d\tau)(\tau-\bar\tau)^{-1}(w-\bar w)$ in independent conjugate marking variables on this fixed curve. Continuing a conjugate mark round a cycle translates $\bar w$ by a period weighted by its source row. The rank-four real source matrix gives a Zariski-dense group of translations in the four period vectors (already the $A$-cycles span after complexification). Thus the determinant polynomial would vanish for all translations, including arbitrary additions of holomorphic forms to the four logarithmic differentials by the infinitesimal period formula. In the basis $(1/y,1/y^2,x/y^2,x^2/y^2)\,dx$, the shifted multiplication evaluated at sixteen ordinary points has determinant $172$ for $(1,1,1)$ and $129$ for $(1,1,-2)$ modulo $211$. The simple roots used for the marked and evaluation points lift to characteristic zero; all denominators are units at the chosen reduction. Hence each determinant is a nonzero characteristic-zero minor. [Exact multiplication matrices and 125 bracket restrictions](worldsheet_fourfold_positivity.py)

Take the first measurement minus the second minus the third plus the fourth. Fusing the selected $v$ with a planar wave in the transverse $e_1$ component gives four copies of the $u$ fusion and weight $-4$ for the $w$ fusion. The common retained row subtracts the traces with the corresponding planar wave omitted. The selected-slot subtraction drops out because the planar momenta sum to zero. Dividing the resulting bulk pairing by four gives the exact combination of elementary-field variations, with the heat normalization restored:
$$
(Q_4\phi)_1-(Q_4\phi)_2-(Q_4\phi)_3+(Q_4\phi)_4
=\frac{-432J_4}{(2\pi i)^{20}}.
$$
The reduced coefficient is strictly positive. The computed local counterterm realizes this value; every global trace solution would give the same value. Indeed the four-row difference of the table is $(-486,0,594)$, which annihilates the ambiguous polynomial with coefficients $(22,-5,18)$. On the explicit primitive the first restriction has coefficient $24J_4/27$. The measurements use fixed point coordinates and their common power of the string coupling is $\lambda^8$.

Two individual measurements with integer planar momenta determine both periods. Each pair in the first column gives the components along $e^1,e^2$; the three retained momenta and the balancing wave are as above.

| Seven planar waves | $(A,B,D)$ | $(2\pi i)^{20}Q_4\phi(0)$ |
|---|---|---|
| $(1,0),(1,0),(2,0),(0,2),(0,2),(1,2),(-2,-1)$ | $(168192,516096,-62208)$ | $896K_4+92160J_4$ |
| $(2,0),(2,0),(0,-2),(0,1),(1,1),(2,2),(1,-1)$ | $(244224,1299456,62464)$ | $2256K_4+72704J_4$ |

Both rows satisfy $22A-5B+18D=0$. The response $BK_4/576+(8A-B)J_4/9$ gives the displayed coefficients, whose determinant is $-142770176$. Taking $141$ times the first response minus $56$ times the second gives $8923136J_4<0$, so at least one of these two integer measurements is nonzero. To construct the rows, scale each of seven fixed planar directions by an independent nonzero integer. After removing the product of the seven squared scales, the ambiguity is a homogeneous quadratic in those scales by the two degree patterns above. Fixing six scales and testing the last quadratic for a square discriminant produces rational solutions; clearing denominators gives the listed waves. [Exact scalar replay](operator/four_loop_planar.py)

A single complex momentum configuration also eliminates the undetermined trace coefficient. Keep the same vector, root, and final two scalar momenta. Use the seven planar scalars
$$
u,u,u,v,v,\quad t u+v,\quad -(3+t)u-3v,
\qquad 44t^2+45=0 .
$$
Their sum is zero. At arbitrary $t$ the first scalar restriction $A$ is $18(64t^2+117)$, the second $B$ vanishes, and the third $D$ satisfies $22A-5B+18D=576(44t^2+45)$. Hence the entire kernel of the restricted trace variation vanishes on this test. At either root $A=10206/11$. The solved planar interaction, and every global completion by a doubly soft primary trace, therefore give the **same individual elementary-field value**
$$
Q_4\phi(0)=\frac{9072 J_4}{11(2\pi i)^{20}},
\qquad (2\pi i)^{20} Q_4\phi(0)<-\frac{37215067449}{44\cdot10^{18}}.
$$
The inequality uses the complete fifty-shape bound above, multiplied by $9072/(11\cdot1728)=21/44$.

The complete weight-five component of this elementary differential is now fixed. Write $[\zeta(5)]$ for the coefficient in the graph-period reduction to $1,\zeta(3),\zeta(5)$. The doubled-edge cores give zero in this column. Adding the $28$ prism and $11$ bipartite shapes, with their original symmetry weights, gives
$$
\begin{aligned}
[\zeta(5)]K_4
&=-\frac{13914257373753679}{45754008300000}
+\frac{249372365272156031}{806626220400000}
=\frac{109867352441461633}{21778907950800000},\\
[\zeta(5)]J_4
&=-\frac{1633443038914372549}{1045387581638400000}
+\frac{2209968587039007449}{1393850108851200000}
=\frac{13733372208504593}{597364332364800000}.
\end{aligned}
$$
Consequently, at either root of $44t^2+45=0$ on the seven-planar-wave test,
$$
[\zeta(5)]\big((2\pi i)^{20}Q_4\phi(0)\big)
=\frac{9072}{11}\frac{13733372208504593}{597364332364800000}
=\frac{13733372208504593}{724317422400000}.
$$
The four real-wave measurements in the preceding difference have combined reduced $\zeta(5)$ coefficient $-13733372208504593/1382787806400000$. The genus-four coupling is $\lambda^8$. These numbers hold for the solved planar restriction and for every extension by a doubly soft primary trace. [Prism reduction](quartic_prism_physical_ibp.py) · [Independent weight-five bipartite cuts](quartic_recut_bank.py) · [Planar polynomial identities](operator/planar_trace_all_genera.py)

The complete fifty-shape result in the same reduced convention is
$$
\begin{aligned}
K_4={}&-\frac{45974661690227570774909}{37045922424310800000000}
-\frac{63695932815585779293}{19189807005600000000}\zeta(3)
+\frac{109867352441461633}{21778907950800000}\zeta(5),\\
J_4={}&-\frac{6436431208614965917972531}{1138050736874827776000000000}
-\frac{1593985240704793516397}{105375068229150720000000}\zeta(3)
+\frac{13733372208504593}{597364332364800000}\zeta(5).
\end{aligned}
$$
For either root of $44t^2+45=0$ the elementary field therefore has the full coefficient
$$
(2\pi i)^{20}Q_4\phi(0)=
-\frac{6436431208614965917972531}{1379911607762688000000000}
-\frac{1593985240704793516397}{127769593311360000000}\zeta(3)
+\frac{13733372208504593}{724317422400000}\zeta(5).
$$
The three coefficients multiply the same wave configuration and the same power of $\lambda$. In units of $10^{-20}$ the reduced response lies strictly between $-1141009435895$ and $-1141009435894$. All $3444054$ bipartite integration-by-parts equations hold over the integers after clearing denominators in each of the three columns. Fixing the exact wheel faces and star–triangle values gives rank $695485$ on the $695485$ interior orbits. Three dimension-raise sums reproduce the original dimension-six integrals in all three columns. The aggregate covers $3+8+28+11=50$ shapes. Direct evaluation of the unrestricted genus-four graph sum on arbitrary twelve-wave configurations remains open. The worldsheet integral above gives the analytic vertex. [Exact residuals](quartic_bipartite_certificate.cpp) · [Bipartite graph weights](quartic_bipartite_totals.cpp) · [Independent original-weight sum and elementary response](quartic_scalar_response.py)

Even without the star–triangle inputs, the homogeneous interior equations with zero wheel faces have a one-dimensional kernel. Its normalized rational vector satisfies every equation and has support on exactly $3410$ orbits. The support is disjoint from the $465564$ orbits of the two physical numerators, so each contributing bipartite period is individually fixed by four-wheel boundary data. The additional calibration fixes the remaining nonphysical direction. From the project root, `python3 explore/iib-bcov/quartic_bipartite_replay.py` reconstructs these period tables. [Sequential exact replay](quartic_bipartite_replay.py)

These five elementary-field measurements survive finite translation- and $SL_5$-equivariant holomorphic changes of the marked scalar with fixed classical field, bulk action and output coordinates. A genus-$s$ marked change contributing to this one-vector primary component has $2s+1$ scalar arguments and $10s$ derivatives by the arity and rank counts. For $s>0$ each invariant contains a five-row determinant, so at most $2s-2$ external planar waves can enter the mark: the other arguments must supply three independent directions modulo their plane. At geometric loop number $\ell>0$, ordinary cubic vertices can take at most $2\ell$ planar waves. Each positive-genus bulk trace of genus $h$ takes at most $2h-1$. At total genus $n$ this gives at most $2n-2$ planar waves, while each probe has $2n-1$. The tree with a positive-genus bulk trace satisfies the same bound. In the remaining tree with a classical cubic, the vector with polarization $e_4$ annihilates every planar scalar momentum and its own momentum $e^3$. Fusion with a nonplanar scalar leaves all $2n-1$ planar waves on the mark. The free scalar descendant also gives $e^3(e_4)=0$. This proves invariance for each measurement separately in the specified marked complex. The value on the single complex configuration additionally annihilates the bulk ambiguity left by the planar equations.

#### Complete four-loop response on seven planar waves

With the scalar root, transverse vector and balancing wave fixed above, all three bracket polynomials $A,B,D$ now have exact coefficients. The last graph in the four-hard auxiliary integral was prism 16, with $732693$ projective moment orbits. The $3+8+28+11=50$ auxiliary graphs give, in the existing fourfold normalization,
$$
W=\frac{136438245856046429534016581}{7966355158123794432000000000}
+\frac{191114209038169441669}{4172089805452800000000}\zeta(3)
-\frac{827038033199479}{11879404336800000}\zeta(5).
$$
For arbitrary seven planar momenta the elementary response is
$$
(2\pi i)^{20}Q_4\phi(0)=\frac{B}{576}K_4+\frac{8A-B}{9}J_4
+\frac{22A-5B+18D}{288}\left(\frac{3K_4}{8}+8J_4+24W\right).
$$
Four $e^1$ waves and three $e^2$ waves give $-24W$. For arbitrary rational pairs the function `elementary_four_loop(plane)` returns the three exact coefficients of $1,\zeta(3),\zeta(5)$. [All fifty auxiliary amplitudes](../kodaira-spencer/four_loop_arakelov_periods.json) · [Elementary-field evaluation](../kodaira-spencer/arakelov_elementary_measurements.py)

At $a^2+b^2+c^2=3$ the exact period difference satisfies
$$
5872943497<10^{20}\big[V(1,1,-2)-8W\big]<5872943498.
$$
In particular $V(1,1,-1)-V(1,1,1)=[V(1,1,-2)-8W]/9>0$. For fixed absolute strengths the maximum therefore puts the largest one opposite the other two, minimizing $(a+b+c)^2$. A maximum with positive value has $abc\ne0$. At a constrained critical point logarithmic differentiation of $(abc)^2$ times the quadratic form above shows that $1/a,1/b,1/c$ are affine functions of $a,b,c$ with common coefficients. Thus the strengths take at most two distinct values. Up to permutations and a simultaneous sign their maximizing ratio is $(1,1,-t)$, $t\ge1$. After rescaling to the specified sphere the vertex equals
$$
\frac{9t^2\,[4W(2-t)^2+V(1,1,-2)(1+t)^2]}{4(2+t^2)^4}.
$$
Its derivative has the opposite sign to
$$
24W(t-1)(t^2+2)+[V(1,1,-2)-8W](2t^3+3t^2-4t-2).
$$
This cubic increases strictly on $[1,\infty)$ and has its unique zero between $1016629/10^6$ and $1016630/10^6$. The corresponding absolute fivefold elementary response is $24V(a,b,c)$ in the reduced convention. [Certified rational signs and root interval](../kodaira-spencer/arakelov_period_sum.py)

### Two real elementary-field measurements fixing both four-loop periods

Keep the same elementary field, transverse vector and pair of additional scalar waves. For the seven planar scalar momenta take three copies each of $u,v$ and one wave $t(u+v)$. The three planar trace symbols on this family are
$$
A=18t^2(t^2+4),\qquad B=36t^2(t+3)^2,\qquad D=18t^2(2t-1).
$$
The undetermined trace multiplies $22A-5B+18D=72t^2(3t^2-6t-5)$. Both $t=1\pm2\sqrt{6}/3$ give nonzero real waves and eliminate it. The balancing scalar has momentum $-(3+t)(u+v)-e^3-e^4-e^5$. Using the $24$ assignments in the two determinant traces gives the individual reduced elementary-field value
$$
(2\pi i)^{20}Q_4\phi(0)
=\frac{B}{576}K_4+\frac{8A-B}{9}J_4
=\left(\frac{25}{9}+\frac{25t}{6}\right)K_4+(80+96t)J_4 .
$$
$$
\frac{30639247090529969233}{784040686228800000}
\ \pm\ \frac{1517542313707568639}{98005085778600000}\sqrt{6}.
$$
The marked-field redefinition count above applies here: seven planar waves exceed the six allowed by a positive-genus change, and $e_4$ annihilates every planar momentum in the classical cubic fusion. Both real matrix elements are invariant under the stated changes with fixed bulk action and output coordinates. [Polynomial identity on all seven waves](operator/planar_trace_all_genera.py)

### Cancellation on the saturated planar sector at arbitrary genus

For $2n$ planar waves and four retained waves use the same $F,G$ as above, with $G=n(n-1)\sum_{j\ne1,\ j\in\mathrm{planar\ slots}}B_{1j}(u_j)_2$. The bare heat source is $K_n F/[n(n!)^2]+J_nG/[n(n-1)]$. On $2n-1$ planar slots let $A$ be half the sum, over the choice of one repeated slot, of the perfect-matching polynomial of the augmented list. Each unordered matching has coefficient one and weight the product of its squared planar brackets. Let $B$ be the same matching polynomial after adding the sum of the slots as the last wave. Their transverse classical variations obey
$$
n\,n!\,\delta A=(n-1)F+(2n+1)G,\qquad
n!\,\delta B=2(n-1)F\qquad(n\ge2).
$$

Here is a finite polynomial proof uniform in $n$. Fix the tested slot. Polarize all other $2n-1$ waves against one arbitrary finite measure on planar covectors, meaning an independent finite weighted sum in each slot. The coefficient of distinct weights recovers every labeled evaluation by symmetry. A change of frame sets the second-moment matrix of the measure to the identity on its invertible locus over the complex numbers. The mean squared bracket of two independent waves is then $2$. Polynomiality extends the identities across the singular locus.

Distinguish the fusion partner. There are $2n-1$ choices. Write $r=n-1$ so there are $2r$ further waves. In each evaluation of a fused or omitted trace the matching has a selected slot, an appended slot, and $t$ distinguished further waves from the appended argument. For $A$ at most one further wave is selected. Expanding the two copies of the appended sum for $B$ selects at most two. For any fixed coincidence pattern, let $h$ be the number of pairs between these $2+t$ distinguished slots. There remain $2+t-2h$ partners to choose from the undistinguished slots and $r-1-t+h$ ordinary pairs. The ordered choices of the $t$ waves, partner choices and unordered pairings, multiplied by the mean $2$ per ordinary pair, give
$$
\frac{(2r)!}{(r-1-t+h)!}.
$$
Small $r$ with insufficient slots contributes zero. After division by $(2n-1)!/r!$ including the fusion-partner choice, this is a falling factorial of degree $1+t-h$. If $1+t-h=0$ the product is one. The remaining moment contraction is independent of genus. Thus the two averaged left sides, after dividing out their common factorials, have degree at most three in $r$: the extra $n$ in the first identity raises its degree from at most two to at most three. The matching in $F/n!$ has just the fixed slot and the distinguished partner. Its averaged expression has degree at most one with the same divisor. The extra term in $G/n!$ fixes the pair and multiplies by $n$, again giving degree at most one. Four genera therefore suffice. The replay verifies both identities as full polynomials on $4n-4$ free coordinates for $n=2,3,4,5$. The four $F$ polynomials have $9,95,987,9963$ monomials. The four $G$ polynomials have $5,77,907,9613$.

For $n\ge4$ both restrictions extend to symmetric products of squared five-row minors on all $2n+3$ conserved trace slots. For $A$ prescribe $n$ pair minors, one repeated planar wave and all other planar waves occurring once. For $B$ prescribe $n-1$ pair minors and one minor with one planar argument and all four retained arguments. Distribute the $n$ or $n-1$ omissions as evenly as possible over the four retained labels. If $k$ is the omission count at a retained label, prescribe its multiplicity to be $n-k$. Each retained multiplicity is at least three. Symmetrize by summing all unordered minor products with the prescribed multiplicity partition. On this family any planar label of multiplicity at least three would force at least $2n+1$ planar occurrences, while five-row minors allow at most $2n$. The four high-multiplicity labels must therefore be precisely the retained ones. All contributing minors have distinct planar row sets. Dividing by the number of assignments of the omission multiset to the distinct pair minors, summed over its distinct permutations on the four retained slots, gives exactly $A$ or $B$.

Consequently the following restriction of a local trace has the required boundary in the convention $\Phi=C$:
$$
\frac{n!J_n}{(n-1)(2n+1)} A+
\left(\frac{K_n}{2n(n-1)n!}
-\frac{n!J_n}{2n(n-1)(2n+1)}\right)B .
$$
At four loops the two global assignment counts are both $24$. At higher genus the exclusion of lower vertices on this locus applies whenever the lower bulk primaries are traces with the arities above. The result cancels the bare planar source at all such orders; the flat worldsheet vertex supplies the global momentum polynomial at every genus. [All-genus matching identity and four elementary measurements](operator/planar_trace_all_genera.py)

### An elementary-scalar period at each genus

The preceding complex configuration extends to every $n\ge2$. Keep the vector wave at $e^3$ polarized by $e_4$, take absorbed momentum $e^4$, and place two scalar waves at $e^5$ and $-e^3-e^4-e^5$. Use these $2n-1$ further scalar momenta in the plane:
$$
u^{(n-1)},\ v^{(n-2)},\ tu+v,\ -(n-1+t)u-(n-1)v,
\qquad t^2=-\frac{2(n-1)^2(3n-2)}{n^2(3n-1)},\quad u=e^1,\ v=e^2 .
$$
On this configuration the reduced elementary-field value fixed by the planar equation is
$$
(-1)^n(2\pi i)^{5n}Q_n\phi(0)
=\frac{n!(n-1)!(n-1)^2(2n-1)}{3n-1}\,J_n<0 .
$$
The same value is accessible with real momenta. Along this conserved family, at fixed bulk choice, the trace is affine in $t^2$ by the three scalar restrictions below. Multiply its real value at $t=0$ by $1+2(n-1)^2(3n-2)/[n^2(3n-1)]$ and subtract its real value at $t=1$ multiplied by $2(n-1)^2(3n-2)/[n^2(3n-1)]$. The result is exactly the displayed reduced period in the reduced convention. At least one of the two real evaluations is nonzero. Both evaluations are individually invariant under the finite marked redefinitions above.

The sign on the left is the source heat conversion. Directly reconstructed bulk actions give this result at genera two and three. At higher genus the doubly soft primary traces of Li–Wang satisfy the lower-genus hypotheses above. The marked-redefinition count in the four-loop calculation applies verbatim. This is a single elementary field at the mark at every genus. The external waves depend on $n$.

Both resistance periods can be recovered from two **individual real-momentum** elementary-field measurements. Take $n-1$ copies each of $u,v$ and one $t(u+v)$. Use the same root, vector and $e^5$ scalar, and balance the last scalar by subtracting $(n-1+t)(u+v)$ as well as $e^3+e^4+e^5$. In the three trace symbols divide by the common factor $(n-1)(n-1)!t^2$. Then
$$
(A,B,D)\big/\big((n-1)(n-1)!t^2\big)
=\big(t^2+n,\;2(t+n-1)^2,\;2t-1\big).
$$
For $A$, duplicating the oblique wave pairs its two copies once with each axis. Duplicating an axial wave forces the oblique one to pair with that axis. For $B$ the appended wave is $(n-1+t)(u+v)$; it must pair with the axis opposite the first oblique wave's partner. For $D$ the distinguished triple must contain one wave of each direction. Its three choices of center give $t^3,t^3,-t^2$, leaving $(n-2)!$ matchings. These counts give the stated polynomials at every $n\ge2$.

In these units the ambiguous combination $2(3n-1)A-5B+2(2n+1)D$ is $2(n-2)(3t^2-6t+3-2n)$. Put $t_-=1-\sqrt{2n/3}$ and $t_+=1+\sqrt{2n/3}$. Both roots are real and nonzero. At either root the individual reduced elementary coefficient is the preceding explicit linear combination of $A,B$ with $K_n,J_n$; all marked-field redefinition counts apply. To see that the two measurements are independent, set $x=t_-$ and $y=t_+$. Then
$$
(x^2+n)(y+n-1)^2-(y^2+n)(x+n-1)^2
=(y-x)\frac{2n(2n+1)}{3}\ne0.
$$
This follows from $x+y=2$, $xy=1-2n/3$. The two-by-two map from $K_n,J_n$ to the coefficients of $A,B$ is also invertible. Since $J_n<0$, at least one real measurement is nonzero. Both choices have integer momenta when $n=6m^2$ for a positive integer $m$. These exact measurements apply to the flat worldsheet solution at every genus. [Matching and polynomial replay](operator/planar_trace_all_genera.py)

Two individual **integer-momentum** probes separate $K_n,J_n$ at every $n\ge2$. Their $2n-1$ planar scalar waves can all lie on the three positive rays $u,v,u+v$. Write $m=n-1$, take axial strengths $a_i,b_i$ of mean one in each group of $m$ waves, and give the oblique wave strength $t$. Divide the existing bracket symbols $A,B,D$ by their common factor $m\,m!\,t^2(\prod_i a_i b_i)^2$. In these units their values are
$$
t^2+\frac{n}{2m}\sum_{\substack{i=1,\ldots,m\\\mathrm{axial\ slots}}}(a_i^2+b_i^2),\qquad
2(m+t)^2,\qquad 2t-1.
$$
Duplicating an axial slot gives $(m+1)!$ matchings and duplicating the oblique slot gives $2m\,m!$. For $B$, pair the appended wave with an axial slot; its bracket with the oblique slot vanishes at equal axial means. For $D$ the triangle uses one wave of each direction, with $(m-1)!$ matchings left. Its three centers give $t a_i+t b_j-a_i b_j$ after removing the squared strengths.

At $n\ge5$ insert pairs $1+x_j,1-x_j$ into the axial groups and fill each remaining slot with one. There are $2\lfloor m/2\rfloor$ pairs. The ambiguity $2(3n-1)A-5B+2(2n+1)D$ vanishes exactly when
$$
\sum_{\substack{j\\\mathrm{axial\ pairs}}}x_j^2
=\frac{m(n-2)\big(2n-3(t-1)^2\big)}{n(3n-1)}.
$$
For $t=1,2$ the right side is positive and less than the number of pairs. For a right side $p/q$, represent $pq$ by four integer squares and divide the four entries by $q$. Append zeros. Rational orthogonal reflections of this vector can approach the vector with all entries equal and positive on the same sphere. Indeed, take reflection normals tending to the difference of the two vectors. Thus one can choose every $|x_j|<1$. Clearing denominators in all planar waves scales $A,B,D$ by the same positive power; the last transverse scalar is balanced as before.

At $n=2$ all axial strengths are one. At $n=3$ take pair offsets $(1/2,1/2)$ and $(1/2,0)$ for $t=1,2$. At $n=4$ the axial triples and oblique strength after clearing denominators are $(1,16,16),(3,12,18),11$ and $(5,29,32),(11,23,32),44$. In every case $D/B$ has the two distinct values $1/(2n^2)$ and $3/[2(n+1)^2]$. The ambiguity equation then gives distinct $A/B$, so the two integer measurements recover both periods. This provides rational-jet tests even at genera where the equal-strength solutions above are irrational. [Constructive integer waves and direct matching checks](operator/integer_planar_probes.py)

There is a sign constraint on **every** real saturated planar test which eliminates the trace ambiguity. In the expansion on a distinguished triple let $a,b,c$ be its three cyclically oriented planar brackets. Its three choices of center contribute $A=(ab)^2+(bc)^2+(ca)^2$ and $D=-(ab)(bc)-(bc)(ca)-(ca)(ab)$ times the same matching of the other slots. Thus $2(A+D)$ is the sum of the three squared differences of $ab,bc,ca$, while $A-2D=(ab+bc+ca)^2$. All remaining matching weights are nonnegative. On the full test
$$
-A\le D\le A/2,\qquad
2(n-2)A\le5B\le(8n-1)A
$$
when $2(3n-1)A-5B+2(2n+1)D=0$. For $n>2$ a nonzero such test has $0<B<2nA$. Hence both $K_n$ and $J_n$ enter its elementary coefficient with strictly positive rational-polynomial prefactors in reduced units. Extracting $J_n$ alone requires combining real measurements or using complex momenta on this locus.

In particular the **full four-loop weight-five coefficient is positive on every nontrivial real test on this locus** in the stated polarization. The prism and bipartite totals above give $[\zeta(5)]K_4/576-[\zeta(5)]J_4/9>0$. Substituting the two bounds on $B$ gives
$$
\frac{1064888151925289}{41927309424000000}\,A
\ \le\ [\zeta(5)]\big((2\pi i)^{20}Q_4\phi(0)\big)
\ \le\ \frac{923572185515976769}{15680813724576000000}\,A.
$$
For the complete four-loop differential all three columns give $K_4<0$, $J_4<0$ and $K_4/576-J_4/9<0$. The bounds $4A/5\le B\le31A/5$ therefore give, for every real test here with $A>0$,
$$
\left(\frac{31K_4}{2880}+\frac{J_4}{5}\right)A
\ \le\ (2\pi i)^{20}Q_4\phi(0)
\ \le\ \left(\frac{K_4}{720}+\frac{4J_4}{5}\right)A<0.
$$
The exactly summed elementary response has this sign throughout the stated family. Two independent rational enclosures use the Euler transform and the positive ordinary zeta series with integral tail bounds. [Signed intervals](quartic_prism_physical_ibp.py) · [Complete 50-shape sum](quartic_scalar_response.py)

For the third restriction $D$ in the rank test, use the three-slot factor $[p_j,p_k]^2[p_j,p_i][p_k,p_i]$, summed over its center $i$ and the unordered pair $j,k$, times the matching on the other slots. Its boundary satisfies
$$
n\,n!\,\delta D=(n-1)F-(3n-1)G .
$$
In the finite-measure argument above, if the selected slot belongs to the three-slot factor, two further slots are selected. Otherwise all three must be selected along with the partner of the selected slot. After dividing by $(2n-1)!/(n-1)!$, the two counting coefficients have degrees at most one and two in $n-1$, respectively. Multiplication by $n$ gives degree at most three. The four symbolic genus tests verify this third identity as well. The two outputs $F,G$ are independent: the zero-total planar fixture above has $F=0$ with nonzero transfer, and equal counts of $u,v$ give nonzero $F$.

The possible scalar restrictions have individual degrees $(4,2,\ldots,2)$ and $(3,3,2,\ldots,2)$. In the first case the multiplicity is one since the fourth symmetric power of the plane appears once in the $(2n-2)$nd symmetric power of quadratics. In the second, the symmetric square of cubics contains the sixth and second symmetric powers, each once in the $(2n-3)$rd symmetric power of quadratics for $n\ge3$. These multiplicities follow by subtracting adjacent weight counts in the monomial basis of binary quadratics of weights $2,0,-2$. For $n=2$ the sixth power is absent, so $A,B$ already span the space. For $n\ge3$ the dimension is exactly three: on the displayed zero-total family $B=0$, $A$ varies with $t^2$, and $D$ is a nonzero constant as calculated below. Elsewhere $B$ is nonzero. Thus the entire kernel on the restricted symbols is generated by $2(3n-1)A-5B+2(2n+1)D$.

To evaluate on this family write its $n$ waves outside the $u$ line as $(r_i,s_i)$ and put $h_i=r_i/s_i$. Their second coordinates are $1$ with multiplicity $n-1$ and $1-n$ once. In $D$ exactly one of the three distinguished slots must be a copy of $u$: two would make the triangle vanish, and zero would leave too many $u$ slots in the matching. The remaining matching contributes $(n-2)!$ times the squared second coordinates of its non-$u$ waves. After summing the three choices of center the result is $(n-1)!(n-1)^2$ times
$$
\sum_{\substack{i<j\\i,j\,\mathrm{outside\ the\ }u\mathrm{\ line}}}
\big[s_i s_j(h_i-h_j)^2+(s_i-s_j)(h_i-h_j)\big]
=-(n-1)^2-n(n-1).
$$
Here the first sum uses $\sum_i s_i=0$ and $\sum_i s_i h_i=1-n$. For $A$, duplicating $u$ leaves only matchings between the two groups. Duplicating any other wave requires one pair inside the other group. After removing the common factor $(n-1)!(n-1)^2/2$, the first choice contributes $n(n-1)$. The second contributes
$$
\left(\sum_{i\in\mathrm{slots\ outside\ the\ }u\mathrm{\ line}}s_i^2\right)
\sum_{j<k,\ j,k\in\mathrm{slots\ outside\ the\ }u\mathrm{\ line}}(h_j-h_k)^2
+\sum_{i,j\in\mathrm{slots\ outside\ the\ }u\mathrm{\ line}}s_i^2(h_i-h_j)^2 .
$$
All these sums run over the $n$ waves outside the $u$ line. Their slopes are $0$ with multiplicity $n-2$, $t$, and $1+t/(n-1)$. Consequently
$$
\begin{aligned}
A&=(n-1)!(n-1)^2\left[(n-1)(n^2-n+1)+\frac{n^3(n-2)}{2(n-1)}t^2\right],\\
B&=0,\qquad D=-(n-1)!(n-1)^3(2n-1).
\end{aligned}
$$
At the prescribed $t^2$ the kernel vanishes. Substitution in the planar primitive gives the displayed elementary-field period. On the independent integrated seven- and nine-primary vertices the reduced values are $-1/17280$ and $30J_3$. The replay recovers the latter by exact interpolation at $t^2=-7/9$, checking separately the rational and $\zeta(3)$ columns. [Polynomial, graph, and integrated-vertex replay](operator/planar_trace_all_genera.py)

### Comparison with the four-dimensional Yang–Mills supercharge

Both highest-order symbols differentiate $n+1$ input letters. Their field counts differ. The paper's cubic $\mathcal N=4$ fan on $\mathbb C^2$ replaces those letters by two output letters, so a monomial of length $r$ maps to length $r+1-n$. In the fivefold primary sector a graph at geometric loop number $\ell$ has $4\ell+1$ bulk vertices and $5\ell+1$ edges. A primary trace vertex of genus $h$ adds $2h$ half-edges beyond the cubic. At total genus $n=\ell+\sum h$ the output length is therefore
$$
r+3(4\ell+1)+2\sum_{\substack{\mathrm{bulk}\\\mathrm{vertices}}}h-2(5\ell+1)=r+2n+1 .
$$
In particular the Yang–Mills fan annihilates a single letter for every $n\ge1$, while the fivefold $Q_2\phi$ contains five scalar letters and one vector. Even a loop-independent multiplicative substitution of a homogeneous boundary composite of length $L\ge1$ for $\phi$, and length $R\ge0$ for the vector, gives lengths $L-1$ and $5L+R$ on these two sides. It cannot intertwine this nonzero primary component if it preserves the component and its image. [Yang–Mills fan and output fields](../../sources/LLC_main.tex#L209)

The reduced two-dimensional fan itself provides a normalization check. At unit endpoint momentum pairing the paper's zero-displacement fan period is $1/[(n!)^3(n+1)]$. Before cotree labels and the voltage insertion, the fivefold calculation above integrates the squared current measure with mass $2^n n!/[(2n)!(n+1)!(2n+1)!]$. The ratios at $n=1,2,3$ are $1/6,1/90,1/2100$. An overall factor and a fixed rescaling of the loop coupling cancel in their combination $(1/6)(1/2100)/(1/90)^2=9/14$. [Paper's zero-displacement kernel](../../sources/LLC_main.tex#L2303)

At genus three compare the protected elementary-field probe on three $e^1$ waves and two $e^2$ waves with the rational interpolation of the two real elementary probes giving $30J_3$. Their coefficient pairs in $(1,\zeta(3))$ have determinant
$$
\begin{vmatrix}
-49188779/2715939072000 & 7783/523908000\\
-2078212919/21184324761600 & 131371/1634592960
\end{vmatrix}
=\frac{171285139}{96092097118617600000}.
$$
Their ratio is irrational. The paper's fan Taylor coefficients are rational in rational jet and colour conventions after stripping the common phase and heat normalization. Thus no single genus-three normalization can identify both of these elementary matrix elements with rational evaluations of that fan. A comparison allowing quantum mixing of operators would need the corresponding map on the observable complexes. [Protected elementary probe](#the-three-loop-differential-on-the-elementary-scalar) · [Fan Taylor coefficients](../../sources/LLC_main.tex#L2308)

At four loops, the weight-three and weight-five momentum tensors of the elementary scalar are linearly independent on the planar restriction. The full bipartite equations at the prime $999999937$, including the dimension-raise calibration, determine the three coefficients of both periods in this field. The isolated rational coefficients are reduced only at primes that preserve an invertible calibrated minor. Adding the prism and doubled-edge terms gives the following certificate, with columns $K_4,J_4$:
$$
\begin{pmatrix}
450469298&683141421\\
396949588&272572399\\
790650093&875668526
\end{pmatrix}
\pmod {999999937},\qquad
396949588\cdot875668526-272572399\cdot790650093
\equiv323157369\pmod {999999937}.
$$
The rows are the rational, $\zeta(3)$ and $\zeta(5)$ coefficients. Hence the combination $([\zeta(5)]J_4)K_4-([\zeta(5)]K_4)J_4$ has a nonzero $\zeta(3)$ coefficient and is irrational by Apéry's theorem. The two integer elementary measurements span $K_4,J_4$ over the rationals. Adjoining the reduced scalar-fan measurement $13/1209600$ gives exactly the rational span of $1,\zeta(3),\zeta(5)$. A single common rescaling cannot make all three measurements rational. This conclusion uses no independence assumption about $\zeta(5)$. It applies to the solved restriction and to every full bulk trace extending it in the stated complex. [Full field check and physical sum](quartic_scalar_response.py) · [Calibrated bipartite equations](quartic_bipartite_adjoint.cpp)

## The integrated bilinear and simultaneous scalar deformations

The constant top-polyvector translation used for $F_r$ also exists at $r=0$. Write $\eta$ for the constant holomorphic five-vector of unit volume. It is closed, and the primary propagator pairs ranks summing to four. Differentiating the full classical action in direction $\eta$ and restricting to minimal fields gives $F_0$ by exactly the same descendant formula. Its leading term is $\int\Omega\,\phi^2/2$. With only $\phi$ and its level-one scalar descendant $\phi_1$ present the full integrand is $\phi^2/[2(1-\phi_1)]$. The nonzero zero-momentum symbol makes this integrated degree-minus-one functional a nontrivial classical class. [Background coefficient and classical closure](operator/integrated_descendants.py)

In every positive-loop primary graph its bilinear mark has two internal ends and momentum zero. Current conservation identifies the two incident currents up to orientation. Every closed cycle traverses the two edges with that same relative sign. The two columns of the heat determinant, whose entries are cycle coefficient times current component and whose projective row is one, coincide. A univalent marking also gives zero: its incident bridge has zero cycle column and the remaining columns have a positive-time Kirchhoff relation. These arguments apply to holomorphic derivatives of bilinear letters, which multiply the heat density by incident momenta. Thus the bare integrated bilinear has zero primary contribution at every positive graph loop number. [Original graph cofactor checks](operator/integrated_square.py)

At one loop the chosen bulk action retains a descendant in its tree bracket. The first marked equation at five primary slots vanishes. At two loops the signed bulk tree term on the vector leg $i$ is the sum of the six singleton spectator momenta times $P$, the same polynomial of all seven external momenta for every spectator. Conservation gives $-p_iP$. The free variation of a polarized scalar correction with one level-one descendant and symbol $-P$ supplies $+p_iP$. Every higher marked arity has positive descendant level, $N-5$ or $N-7$ at one or two loops, respectively. The free jet complex therefore supplies a formal genus-two lift of the integrated bilinear. [Free jet complex](operator/README.md#eliminating-higher-descendant-equations-in-the-holomorphic-jet-complex)

The lifted scalar backgrounds can be turned on simultaneously in the local master equation. At graph loop number $n$, a local heat graph has $5n+1$ internal edges and $4n+2$ vertices. With $k$ scalar marked vertices and ordinary bulk terms of total polyvector rank five per vertex, the required total rank on the external fields is

$$
5(4n+2-k)-4(5n+1)=6-5k .
$$

For $k\ge2$ this is negative, so every mixed local term vanishes. The free differential lowers descendant level by one and raises the input polyvector rank by one on polarized functionals. Its primitives for the bulk defect of rank six and marked defect of rank one can therefore be taken of ranks five and zero, respectively, at every step of the arity induction. For every $r\ge0$ the completed integrated observable $F_r$ begins with $\int\Omega\,\phi^{r+2}/(r+2)!$. Adding them with independent background parameters of degrees $-(2r+3)$ gives a simultaneous formal local deformation through genus two in the source's unweighted bulk interaction. The primary single-mark solution uses thirteen scalar terms. Its last equation is $419/40-241344/23040=0$. The remaining higher-arity single-mark equations follow from positive-descendant exactness. The rank count above kills the local terms of degree at least two in the background parameters. At a finite renormalization scale the effective interaction also contains induced mixed nonlocal terms. [Complete scalar family, coefficients and certificates](operator/README.md#genus-two-lift-for-every-constant-scalar-background)

## Compactification on $K3$ and $T^4$

In the compactified action the primary scalar is $\gamma\in\Omega^{0,\bullet}(\mathbb C^3)\otimes H^\bullet(Y)[2]$. Write $\mathrm{Tr}_Y$ for projection onto the normalized top class of $Y$. Take $(\mathrm{Tr}_Y\gamma)^r/r!$ at the origin, $r\ge2$. Removing the point cuts $k$ legs and leaves a connected graph of loop number $n-k+1$ on every contributing graph. The ring factors at a cut end insert the unit, dual under the trace to the operation $\mathrm{Tr}_Y$ at the point. Contract a spanning tree. Every remaining line inserts the product of the two elements of a dual basis, summed with graded trace sign. This is $24$ times the top class on $K3$ and zero on $T^4$. The square of the top class vanishes. Thus the bare primary point differential has loop order at most $r$ on $K3$ and $r-1$ on $T^4$. [Compactified action, (2.3.8)–(2.3.16)](https://arxiv.org/pdf/2404.14318#page=14) · [Ring and graph tests](compactified.py)

For the highest differential order, $k=r=n+1$, deletion leaves a tree. Choose one top-class external scalar and all remaining ring labels equal to the unit. In the fan reduction take $n$ scalar waves of momentum $e^1$, one of momentum $e^2$, and one vector wave of momentum $e^3$ and polarization $e_1+e_2$. Scalar jets vanish at the mark as in the fivefold. There is now one cotree color, giving a multiplicity $n!$ and the voltage $G_{B0}+\tfrac12\sum_{j=1}^nG_{Bj}$. The reduced answer for either compactification is
$$
\frac{1}{2^n n!(2n+1)!}.
$$
The first two values are $1/12$ and $1/960$. This attains the torus bound for every $r$. The compactified cubic and the original odd-edge kernels are those in the [threefold bulk calculation](../kodaira-spencer/README.md).

For $K3$ and $n=r$, keep every external ring label equal to the unit. Only $k=r$ can contribute and the ring factor is $24$. In three holomorphic dimensions the construction of the point graph uses a parallel path of three edges, with two new bivalent sources. Grow the point valence to $r$ by replacing an incident scalar edge at each step. At $r=2$ start with the square with two adjacent scalar vertices. At $r>2$ the last inserted path has two adjacent scalar vertices. Use this edge for one additional path replacement away from the point. The two old endpoint momenta prescribe two independent linear combinations of the two new sources. Their span excludes the old edge current, by the same independence induction as for the fivefold. Hence the three-component determinant survives. The voltage-square argument gives nonzero matrix elements with all scalar projections positive and the vector transverse. This attains the $K3$ bound for every $r$.

In particular at $r=n=2$, take unit ring labels, scalar momenta $e^1,e^2$ and vector momentum $e^3$ with polarization $e_1+e_2$. After summing the 24 placements and dividing by the graph symmetries, the scalar-ring values of the two thetas are $1/10080$ and $29/272160$. On these same inputs the reduced two-loop differential is therefore
$$
K3:\quad 24\left(\frac{1}{10080}+\frac{29}{272160}\right)=\frac{2}{405},
\qquad T^4:\quad 0 .
$$
Higher bare loops vanish on this quadratic point operator in both compactifications. For an unmarked bulk graph the remaining ring loop number is $n$ itself, which recovers the bulk vanishing from two loops onward. [Exact compactified replay](compactified.py) · [Independent threefold theta matrices](../kodaira-spencer/theta_periods.py)

### Middle cohomology at the point

For $Y=K3,T^4$ take any polynomial in the linear point fields $\mathrm{Tr}_Y(a\gamma(0))$ with $a\in H^2(Y)$. Cutting $k$ contracted letters inserts their $k$ middle classes into the connected remaining graph. Its $n-k+1$ loops each insert the degree-four handle class. Before including any external colors the trace degree is therefore
$$
2k+4(n-k+1)=4n+4-2k\ \ge\ 2n+2,\qquad k\le n+1 .
$$
The trace has degree four. All bare primary corrections vanish at $n\ge2$, for every polynomial degree. At $n=1$ precisely two letters contract and all three outgoing ring labels have degree zero. The color contraction on the two letters is their intersection pairing $\mathrm{Tr}_Y(ab)$.

For the normalized square $(\mathrm{Tr}_Y(a\gamma(0)))^2/2$, let $p_1,p_2$ be the two scalar momenta and $p_3$ the vector momentum with $p_3(v)=0$. In the preceding reduced normalization the complete box equals
$$
\frac{\det(p_1,p_2,p_3)^2}{24}\,(p_1+p_2)(v)\,\mathrm{Tr}_Y(a^2)
$$
on unit external colors. Indeed, at positions $i,j=1,2,3$ around the root, the integrated grounded Green function is $\min(i,j)(4-\max(i,j))/120$. For either scalar and the vector, summing their six placements gives $1/12$ before the root-preserving reflection factor $1/2$. At coordinate momenta and $v=e_1+e_2$ the box reads $\mathrm{Tr}_Y(a^2)/12$. Polarization contracts the Hessian of an arbitrary polynomial with the same intersection form. The common kernel of all positive-loop bare primary corrections on these polynomials thus consists exactly of the polynomials harmonic for this form. Every Hessian can be tested at an arbitrary constant middle-valued scalar background since the contracted external waves have unit color. In degree two the respective kernel dimensions are $22\cdot23/2-1=252$ and $6\cdot7/2-1=20$. This kernel concerns the loop terms on cochains. [Direct rooted-box and both Frobenius-ring contractions](../kodaira-spencer/middle_point.py)

On $T^4$ allow the colors of the marked linear fields to have any strictly positive degree. The zero handle class forces $k=n+1$. Their combined degree is at least $n+1$, hence all bare primary terms vanish for $n\ge4$. At $n=3$ equality requires four degree-one colors and unit outgoing colors. Take an oriented basis $a_1,\ldots,a_4$ of $H^1(T^4)$ with $\mathrm{Tr}_{T^4}(a_1a_2a_3a_4)=1$. The four odd factors $\mathrm{Tr}_{T^4}(a_j\gamma(0))$ are multiplied in this order. On three labeled scalar waves of momentum $e^1$, one scalar of momentum $e^2$ and the vector of momentum $e^3$ with polarization $e_1+e_2$, their product has reduced top-symbol coefficient
$$
4!\,\mathrm{Tr}_{T^4}(a_1a_2a_3a_4)\frac{1}{2^3\,3!\,7!}
=\frac1{10080}.
$$
Here the sign of permuting four odd marked letters equals the sign of permuting their four exterior-algebra colors, so all $4!$ Wick assignments contribute equally. Reversing the orientation of the four input classes reverses the coefficient. The one-color fan replay restores the full three-dimensional current determinants and integrates the voltage using the cotree sums. Thus the upper bound of three loops on the whole positive-degree torus sector is attained. The corresponding bound on $K3$ is one loop because every positive degree is at least two. [Odd Wick signs, restored fan and integral](positive_torus_point.py)

The same three-loop quartic admits a finite rational evaluator on arbitrary scalar momenta $p_1,\ldots,p_4$ and one vector of momentum $p_5$ with polarization $v$, $p_5(v)=0$. All five outgoing colors are units. An uncontracted marked letter evaluates to $\mathrm{Tr}_{T^4}(a_j)=0$. All four must therefore contract. At geometric loop order $\ell$ a contributing graph contracts at most $\ell+1$ letters. Thus a positive-genus bulk vertex cannot change this particular matrix element at total genus three, with the point product and classical bulk interaction fixed.

Each external field occurs at its own bivalent vertex. A cubic vertex carrying two external fields would attach through a bridge, killing the projective current cofactor. Suppressing the five bivalent vertices leaves the four-valent root $0$ and two trivalent vertices $1,2$. An ordinary self-loop would attach to the rest through a bridge. A root self-loop would disconnect the graph after deleting the mark and give an independent edge-scale relation. The five remaining paths are consequently $01,01,02,02,12$. Each pair of parallel paths must contain at least two of the five bivalent vertices since its three cycle-current rows have a positive-time null vector. Sorting the two paths in each pair and allowing $1,2$ to be exchanged leaves seven path-occupancy shapes. All three-loop integrals are rational. In the numerator use the square of $U^3$ times the current cofactor, multiplied by $U$ times the grounded voltage, over $U^{10}$, where $U$ is the Kirchhoff polynomial. Integrating along each path gives an ordinary simplex moment. Parallel paths of resulting powers $p,q$ contribute $\mathrm B(10-p,10-q)$ and leave power $p+q-10$; the final triangle gives another simplex moment.

Take the five distinct, conserved momenta and the transverse vector
$$
(p_1,p_2,p_3,p_4,p_5)=
((2,1,0),(-1,0,2),(1,-2,-1),(0,1,-1),(-2,0,0)),\qquad v=(0,1,0).
$$
The seven shape contributions on the ordered four-letter product, in lexicographic occupancy order, give
$$
-\frac{3084401}{49612500}+0-\frac{29009957}{198450000}
-\frac{210359}{29767500}+0-\frac{8907809}{127575000}
-\frac{1445041}{51030000}=-\frac{59263}{189000}.
$$
Since the total outgoing momentum vanishes, this component survives the quotient by total holomorphic derivatives. Its collinear specialization gives the independent fan result $1/10080$. The replay compares two series–parallel reductions on every retained exponent tuple and the original currents and voltage on each contributing topology. [Complete polarized evaluator](torus_quartic.py)

The integrated four-letter product already varies classically. Differentiating the cubic $\alpha\mu\partial\gamma$ in the compactified action gives the scalar transport $\mu(\partial\gamma)$ in the primary projection. Let $b_j\in H^3(T^4)$ be the basis dual to $a_i$ under the trace. Take four scalar colors $(b_1,b_2,b_3,1)$ with momenta $(-1,-1,0),0,0,(1,0,0)$, and a vector of color $b_4$, momentum $(0,1,0)$ and polarization $(1,0,0)$. The momenta sum to zero. Because $b_4 b_j=0$, only the unit scalar can differentiate. Its momentum on the polarization is one, and the four successive pairings $\mathrm{Tr}_{T^4}(a_i b_i)$ are all one. The absolute coefficient of the classical variation is one. [Exterior-algebra tree replay](positive_torus_point.py)

More generally, let $F$ be a homogeneous derivative-free scalar density of degree $m\ge2$, and denote its graded polarization by $F(x_1,\ldots,x_m)$ on ring colors $x_j$. Closure of its integral under all primary colored divergence-free vectors forces
$$
F(bx_1,x_2,\ldots,x_m)=F(x_1,bx_2,\ldots,x_m).
$$
Here odd colors carry odd test parameters so the arguments have even total parity; extraction restores all Koszul signs. Give two scalar slots opposite momenta and vary by a constant vector of color $b$ whose polarization pairs to one with the first momentum. Conservation gives precisely the displayed equation. The other slots are interchangeable by graded symmetry. Transferring their colors onto one slot gives $F(x_1,\ldots,x_m)=F(x_1\cdots x_m,1,\ldots,1)$. Conversely any density of this form varies by a total derivative since every component of the vector is divergence-free. By nondegeneracy of the ring trace this is a density $\mathrm{Tr}_Y(a\gamma^m)$ for some $a\in H^\bullet(Y)$. Its ring tensor is a single trace at the root. Contracting a spanning tree now leaves one handle element for each graph loop, so on $T^4$ every bare positive-loop graph vanishes on these integrated primary scalar densities. Four degree-three colors have zero product. Thus no derivative-free scalar completion carrying the displayed $H^3(T^4)$ quartic coefficient can satisfy even the classical transport equation.

The obstruction extends to translation-invariant local holomorphic differential polynomials in all three superfields of the minimal compactified action. On holomorphic tests, set $\alpha=0$ and keep only the term linear in one vector $\mu$. Varying $\alpha\mu\partial\gamma$ with respect to $\gamma$ gives a factor of $\alpha$ or its derivative in $Q_0\alpha$. Varying the action with respect to $\mu$ gives only $\mu^2$ and $\alpha\partial\gamma$ in its interaction part, so $Q_0\mu$ vanishes to first order on this slice. The free Dolbeault variations vanish as well. Terms of the density containing $\alpha$ or $\mu$ therefore make no contribution to the linear-vector test. A scalar quartic of holomorphic jet order $d$ contributes at order $d+1$ under uniform momentum scaling. The order-one term comes entirely from the derivative-free scalar quartic, where the preceding transfer identity excludes the four-three-form coefficient. One may take the coefficient of the top Dolbeault component in a spectator superfield when testing the integrated holomorphic density. [Fernández–Paquette–Williams, (2.3.13)–(2.3.16)](https://arxiv.org/abs/2404.14318)

For every polynomial in point scalar letters paired with strictly positive cohomology colors, the full bare primary calculation on $T^4$ uses just the one-loop box, two two-loop thetas, and the seven three-loop shapes above. On $K3$ only the box survives. Indeed at a nonzero term the connected remainder after cutting the point is always a tree: this follows from the vanishing torus handle and, on $K3$, from the degree bound. For a single even test color $a$ of positive degree, with odd colors multiplied by odd test parameters, the ring and spectator factors for $(\mathrm{Tr}_Y(a\gamma(0)))^r/r!$ at loop order $n$ are
$$
\frac{(\mathrm{Tr}_Y(a\gamma(0)))^{r-n-1}}{(r-n-1)!}
\,\mathrm{Tr}_Y(a^{n+1}x_1\cdots x_{n+2}),\qquad r\ge n+1,
$$
where $x_1,\ldots,x_{n+2}$ are the outgoing ring colors. For $r<n+1$ the coefficient is zero. Polarization gives every product of distinct letters and its graded signs. The remaining geometric coefficient is evaluated on $n+1$ labeled scalar waves and one transverse vector. At two loops the two path lengths are $(1,3,3)$ and $(2,2,3)$. On two scalar momenta $e^1$, one scalar $e^2$ and vector $e^3$ polarized by $e_1+e_2$ their contributions are $7/14400$ and $1/1800$, respectively. One implementation integrates the grounded voltage against the squared six-row current determinant in the three path-position cubes. A second integrates the original seven segment times by series–parallel reductions. They agree separately on both graphs and on arbitrary-wave fixtures.

At zero total momentum of the contracted-wave kernel, the one-loop three-row determinant vanishes. At two loops three independent scalar momenta then determine the vector momentum as minus their sum. Exchanging the three scalar slots permutes a basis, while the fourth power of its determinant is permutation invariant. The covector to be paired with the vector polarization consequently has the same coefficient on all three basis covectors and vanishes on a transverse polarization. The identity extends to rank-deficient waves by polynomiality. The torus quartic is therefore the first fully contracted positive-color primary scalar symbol that survives at zero total momentum. [Exact rational kernels at all contributing orders](positive_color_scalar.py)

## The bulk primary hexagon

For the primary bulk contact at $n$ loops the total external rank is $5(4n+2)-4(5n+1)=6$, whereas the number of primary external legs is $3(4n+2)-2(5n+1)=2n+4$. Thus only one loop can contribute with every external field of vector rank one. On six divergence-free primary waves the polarized one-loop contact vanishes already at fixed Schwinger times, for arbitrary valid polyvector ranks.

Indeed the ungrounded odd integral has one common five-component translation $\theta$ and one residual odd variable $\xi$. Its external factor is $\int d^5\theta\,d\xi\,\prod_i F_i(\theta+\Phi_i\xi)$. The ordinary and odd cycle factors together are proportional to $\det(p_1,\ldots,p_5)^2$, symmetric on all six conserved momenta. The remaining expression is linear in the potentials. Averaging all external placements gives a common diagonal and a common off-diagonal coefficient in the resistance Green matrix. Conservation reduces its action to a multiple of $p_i$ at slot $i$. The final contraction is proportional to $\iota_{p_i}F_i=0$ at each slot. At unit total resistance the integrated diagonal minus off-diagonal response is $1/1440$ with simplex volume included. One transverse six-vector test gives 704 nonzero ordered exterior coefficients among the 720 orders and zero polarized sum. [Raw thirty-variable exterior evaluation and longitudinal control](hexagon.py)

## Three-loop bulk potentials outside the six-planar locus

On a labeled cubic graph push the squared current cofactor forward at fixed total path resistances. On each ordered path, expand the cofactor in polynomials orthogonal for the uniform distribution on the simplex of segment lengths. The joint degree bound on one path is two. This turns the integral of its square into a diagonal sum of core polynomials. To sum all orders on a path, use shifted Legendre polynomials for its independent labeled positions in $[0,1]$. Subtracting consecutive edge columns gives the fixed jumps $C_P\otimes p_i$, up to their permutation. The remaining initial current depends on the resistance mean of the sources, linear in their individual positions. The square therefore has the same position polynomial in every ordering chamber. The grounded Green function adds the term $T_P(\min(s_i,s_j)-s_i s_j)$ within path $P$ of resistance $T_P$ and retains linear interpolation from the branch vertices. For nonnegative integer exponents $k,l$ the cube integration uses
$$
\int_0^1\!\int_0^1 s^k t^l\min(s,t)\,ds\,dt
=\frac{1/(k+2)+1/(l+2)}{k+l+3}.
$$
For the ordered ten-wave tetrahedron with path lengths $(3,3,3,3,2,2)$ the replay's transverse potential difference is $1284691025791082490901/67085096068410900480000000-84144127971551\zeta(3)/5281958307225600000$. This is one labeled graph. [Current determinant](../nonSP/verify/k4_periods/README.md#a-linear-determinant-for-the-current-numerator) · [Path integration](path_moments.py) · [Ordered controls](labeled_bulk.py)

Use three copies of $e^1$, two of $e^2$, two of $e^3$, one of $e^4$, one of $e^5$ and their balancing covector. Suppressing the three $e^1$ slots gives one cotree sum with factor $3!$ and removes one power of $U$ from the denominator. Write $L,A,B$ for the three repeated groups, respectively, and $T$ for the three slots carrying $e^4,e^5$ and the balancing covector. No six momenta of this source lie in one plane.

For subsets of slots, let $R_{AB}$ denote the average resistance over one slot from each indicated group, using distinct slots for $R_{TT}$. Take the vector polarization $e_3$. The potential at an $A$ slot minus that at the $e^4$ slot has average
$$
-R_{AB}+R_{AT}+R_{BT}-R_{TT}.
$$
Indeed the two $B$ sources each return to the balancing slot. The voltage from a unit dipole is half the difference of the four crossed resistances. The squared cofactor is unchanged under all six permutations of $T$ because the corresponding linear changes of source coordinates have determinant $\pm1$. For the potential at an $L$ slot minus an $A$ slot, average also under exchange of $A,B$. This gives
$$
R_{LT}-\tfrac12(R_{LA}+R_{LB})+R_{AB}-\tfrac12(R_{AT}+R_{BT}).
$$
The uniformly sampled $L$ slot has integrated mass $T_P(\partial U/\partial T_P)/(3U)$ on path $P$ and is uniform in path resistance there. This follows by marking one of the three edges in the Kirchhoff cotree. Integrating the Brownian bridge against that point gives $s(1-s)/2$.

The tetrahedron and the cubic double square have automorphism orders $24,16$. Choose two distinct paths for each of $A,B$ and a multiset of three paths for $T$. The $12600$ configurations per core reduce to $316$ and $595$ orbits after core symmetries and the exchange of $A,B$. Each orbit counts its size times $2!2!3!/\prod_{P\in\mathrm{paths}}t_P!$ label assignments, where $t_P$ is the number of retained $T$ slots on path $P$. Multiply by $3!$ for the deleted slots and divide by the core automorphism order. Parallel unsubdivided edges or dependent current jumps give zero before integration. The two-loop version on three parallel paths has $16$ assignment orbits and reproduces the independent eight-wave longitudinal difference $17/77760$ and the cotree-to-retained transverse difference $-7/77760$.

| Contribution in reduced units | $(\Phi_{e^2}-\Phi_{e^4})(e_3)$ | $(\Phi_{e^1}-\Phi_{e^2})(e_3)$ |
|---|---:|---:|
| Tetrahedron | $-9261411221/56491532697600+459733\zeta(3)/3405402000$ | $-85970041/38516954112000+23449\zeta(3)/14859936000$ |
| Double square | $-829523/1975228416000$ | $-23503/219469824000$ |
| Genus-two vertex in the hexagon | $-1/1451520$ | $-1/13063680$ |
| Total | $-23310136147/141228831744000+459733\zeta(3)/3405402000$ | $-37217287/15406781644800+23449\zeta(3)/14859936000$ |

The two forms in the last row are nonzero by irrationality of $\zeta(3)$. The two corresponding trace-boundary rows raise the rank of the six-planar constraints from two to four among the eighteen doubly soft nine-field trace symbols. The other fourteen coefficient combinations are undetermined by these four tests. [Complete core sums](fourfold_bulk.py) · [Genus-two insertion](bulk_three_hexagon.py) · [Exact trace boundaries](operator/three_loop_soft_jet.py)


### Three-loop bulk consistency on general conserved waves

The local obstruction at genus three, evaluated on ten divergence-free primary fields, vanishes quadratically when an external momentum tends to zero, modulo a common shift of the potentials and the longitudinal terms on the individual slots. The bare heat cofactor occurs twice. In the lower-genus hexagon a wave at a singleton belongs to both copies of the five-singleton determinant, and a wave at the seven-field vertex enters at least one squared minor in each term of its trace symbol.

The resulting space of symmetric covector potentials has dimension $770-96-156=518$. The subtractions are the common translations and longitudinal shifts. Cubic BV consistency acts at eleven primary inputs by summing pair fusions using holomorphic divergence on the product of the fused fields. The lower-genus master equations imply that the actual ten-primary graph sum lies in its kernel. An integer minor of size 500 in this map is $27$ modulo $1009$. At a second conserved eleven-wave configuration, 720 independently sampled transverse tests also have rank $500$ modulo $1009$. The eighteen independent variations of the nine-field trace give eighteen kernel elements by the Jacobi identity. Consequently
$$
518-500=18,\qquad
\ker(\text{eleven-primary consistency})=
\operatorname{im}(\text{nine-primary trace variation})
$$
in the quadratic-soft space. The full three-loop transverse primary bulk defect has a unique nine-field counterterm in this trace basis. Positive-descendant exactness supplies the formal higher-arity lift in the flat minimal holomorphic theory. All eighteen coefficients have been integrated in the rational and $\zeta(3)$ channels. The exact matrix has determinant $515\pmod{1009}$, and the action contains the negative of the reconstructed trace. Every integrated scalar background lifts simultaneously through genus three. On the one-vector primary defect, the bivector consistency equation gives a horizontal gradient plus a scalar times the gradient of the field. Two-vector consistency makes the latter scalar a variational derivative. This constructs the required local scalar primitive. The same implication holds at any higher genus through which the bulk master equation is solved. [Exact eighteen coefficients and diagram periods](operator/three_loop_coefficients.json) · [Local primitive and degree count](operator/README.md#genus-three-lift-of-every-integrated-scalar-observable) [Consistency map and signs](operator/soft_consistency.py) · [Soft target and determinant certificates](operator/README.md#the-full-quadratic-soft-three-loop-bulk-target) · [Descendant resolution](operator/README.md#eliminating-higher-descendant-equations-in-the-holomorphic-jet-complex)

All eighteen reduced coefficients lie in $\mathbb Q+\mathbb Q\zeta(3)$. On the projective simplex the tree formula bounds the edge currents, and the projective cofactor is bounded by deleting an edge of time at least $1/E$, where $E$ counts the graph edges. Voltages are bounded sums of edge drops. Thus the squared-cofactor voltage integrals are absolutely convergent by comparison with the absolute single-cofactor kernel. The monomial convergence argument and the series–parallel or tetrahedral reduction apply to the squared numerator and the two-forest polynomials. The seven-field hexagon insertion is rational. A rational full-rank evaluation matrix for the eighteen trace boundaries recovers the counterterm with coefficients in the same two-dimensional period space.

The two evaluated nonplanar rows give the following determinant of the rational and $\zeta(3)$ coefficients, in table order,
$$
\left(-\frac{23310136147}{141228831744000}\right)\frac{23449}{14859936000}
-\left(-\frac{37217287}{15406781644800}\right)\frac{459733}{3405402000}
=\frac{168255455023}{2562455923163136000000}.
$$
Hence the rational and $\zeta(3)$ parts of the required primary counterterm are linearly independent over $\mathbb Q$ as nine-field symbols. [Convergence and period reduction](../nonSP/verify/k4_periods/README.md#schwinger-numerator-and-convergence) · [Trace-boundary injectivity](operator/three_loop_soft_jet.py)

### Eleven evaluated combinations of the nine-field trace

Keep three external covectors equal to $e^1$. Deleting these vertices removes the first current component from the cofactor, including the first component of every retained momentum. A deleted source still contributes to voltages. For polarizations $v_i$ with $p_i(v_i)=0$ and zero total polarization, the reduced voltage contracted with these polarizations is the resistance sum with weights $-[p_i(v_j)+p_j(v_i)]/2$ on unordered distinct slots. These weights sum to zero at each slot. No deleted-deleted distance occurs, since the deleted vectors have transverse polarization. Thus every term involving a deleted point uses its one-point cotree marginal. This also applies when the first components of the retained momenta differ.

In the two ordinary projected directions take the four slots to be $(e^2,2e^2,e^3,e^3)$, then $(e^2,2e^2,e^3,2e^3)$. A third source adds respectively $(1,2,3,4)$ times $e^1$ to those four momenta of the second source. A fourth replaces the last slot of the first source by $e^2+e^3$. The two further basis waves and the balancing wave have projected coordinates interchangeable in the squared cofactor. Summing their six assignments, all path orders, both cubic cores, and the hexagon with the seven-field vertex gives seven independent equations in addition to the four above. As a check, on restored five-dimensional graphs all six permutations of the retained labels and all six placements of the deleted labels reproduce the current minors and the voltage contractions at unequal resistances.

For example, on the first source the potential difference between the $2e^2$ and $e^2$ slots in polarization $e_3$ is
$$
\frac{1945204969}{1681295616000}-\frac{305933}{324324000}\zeta(3).
$$
The eleven-row trace-boundary matrix has a minor equal to $798$ modulo $1009$. The seven remaining measurements can all be taken with polarization $e_1$ on three explicit fixtures. The augmented eighteen-row determinant is $338$ modulo $1009$. The stored rows use the ordered trace basis in the soft-symbol calculation and the boundary convention $\Phi=C$. [All eleven exact rational and zeta coefficients](operator/three_loop_constraints.json) · [Labeled sums and restored controls](scaled_bulk.py) · [Rank certificate](trace_constraints.py) · [Seven remaining measurements](trace_reconstruction.py)

### Four additional exact nine-field trace measurements

Take, in order and with zero-based labels, the ten momenta $e^1,2e^1,e^2,3e^2,e^3,4e^3,e^4,5e^4,e^5$ and their balancing covector. Evaluate the potential at slot 2 minus that at slot 3, 4, 5 or 6 on the polarization $e_1$. These four rows bring the rank of the evaluated system to 15 modulo 1009. For example the first complete reduced potential, including the lower-genus hexagon, is

$$
-\frac{926178603915493}{58845346560}+\frac{142177795049}{11351340}\zeta(3).
$$

Both core numerators were integrated at independent primes until their reconstruction moduli exceeded twice the respective proven numerator bounds. Each reconstructed fraction also passed the complete cubature at the reserved prime 1009. [Four exact measurements with core and hexagon contributions](operator/three_loop_four_pair.json) · [Reconstruction and replay](trace_reconstruction.py)

### The three-loop differential on the elementary scalar

Fix the above bulk action and elementary point coordinate $\phi$. At positive graph-loop number a one-letter insertion has a bridge. Its cycle column is zero; on the other edges the Kirchhoff identity supplies a right null vector for the remaining square current matrix. Thus only the one-edge tree with the quantum bulk vertex survives. Let $p_{\rm root}$ be minus the sum of the external momenta, and let $v$ be the polarization of the primary vector. The nine-primary action enters as $-C\operatorname{Tr}$ in the reduced convention above. Its contribution to $Q_3\phi(0)$ is $p_{\rm root}(v) C$ on the eight external momenta together with the root. Indeed the divergence acts on their product, giving the sum of their momenta evaluated on $v$; its sign combines with the minus in the action. At two loops the same formula uses the seven-field symbol $P$.

Use the vector momentum $e^5$ with $v=-e_4$, and the scalar momenta in the following table. Then $p_{\rm root}(v)=1$. The coefficients use the reduced graph convention.

| Operation | Scalar momenta | Elementary-field matrix element |
|---|---|---|
| $Q_2\phi$ | $e^1,e^1,e^2,e^3,e^4$ | $-1/31104$ |
| $Q_3\phi$ | $e^1,e^1,e^1,e^2,e^2,e^3,e^4$ | $-49188779/2715939072000+(7783/523908000)\zeta(3)$ |
| $Q_3\phi$ | $e^1,e^1,e^2,e^2,e^3,e^3,e^4$ | $-568206164203/14405340837888000+(347115523/11115232128000)\zeta(3)$ |

The first two matrix elements survive any finite translation- and $SL_5$-equivariant holomorphic redefinition of the marked scalar preserving its classical value, with the bulk action and output coordinates fixed. At genus $s$ a contributing nonlinear marked change on these primary fields has $1+2s$ scalar slots and $10s$ derivatives. At genus one it vanishes because three scalar momenta cannot fill a five-row determinant. For the genus-two or genus-three tree change only the $e^4$ scalar can contract with $v$. Replacing it by $e^4+e^5$ gives momenta of rank four, so every positive-degree scalar invariant vanishes. At genus three the remaining hexagon on a genus-two mark requires one active momentum on each coordinate line. All six active label sets on the second row leave two $e^1$ and one $e^2$ as spectators. Even with its two internal currents the five-slot genus-two mark has rank at most four. The free variation of a scalar descendant is killed by $e^5(v)=0$. This fixes a nonzero three-loop correction on the elementary scalar in this marked-observable sense.

For seven arbitrary rational scalar momenta and one transverse vector wave, the complete eighteen coefficients now evaluate the same differential directly. For vector momentum $e^5$, polarization $-e_4$, and scalar momenta $e^1,e^2,e^3,e^4,e^1+e^4+e^5,e^2-e^5,e^3-e^4$, the nine conserved covectors have pairwise distinct directions and $p_{\rm root}(v)=1$. The reduced result is
$$
-\frac{597710892167297}{43216022513664000}
+\frac{36860229443}{3334569638400}\zeta(3).
$$
The full expression has momentum degree $31$, from six five-row minors and the one root contraction. [Rational-wave evaluator and scaling check](elementary_scalar.py) · [All eighteen coefficient pairs](operator/three_loop_coefficients.json)

The full reduced trace $C$ on eight external momenta supported on five independent lines follows from six rational pairs. In a coordinate basis its degree on each line is six by determinant covariance. Every external wave supplies at least two powers of its individual strength. Four momenta on one line force zero. For multiplicities $(3,2,1,1,1)$ or $(2,2,2,1,1)$, pull out the sixth power of each singleton strength and the square of every strength on a repeated line. On each two-wave line with strengths $a,b$ the remainder is a symmetric quadratic. Its basis is $(a+b)^2,ab$. Sum over choices of those lines supplying $ab$, with all other pairs supplying $(a+b)^2$. Ordered by the number chosen, the two or four coefficients are in the linked table. For another five-line frame multiply by its determinant to the sixth power. Six trace rows obtained by reversing selected second strengths all lie in the span of the fifteen integrated boundary measurements; the replay verifies the row identities over the rationals and checks two additional unequal-strength inputs on all eighteen trace coordinates. [Six exact period pairs](elementary_scalar.json) · [Trace extraction and marked redefinitions](elementary_scalar.py)

For real external momenta this also determines the complete zero set on five independent lines. At each paired line write the residual quadratic in the basis $(a+b)^2/4,(a-b)^2/4$. Its coefficients are precisely its values at $(1,1)$ and $(1,-1)$. The three-pair case uses the tensor product of these bases. By permutation symmetry there are six distinct corner values across the two families. Writing each as $r+s\zeta(3)$, the rational replay gives $s>0$ and $1728r+2089s<0$ in all six cases. The series tail after the fourth term is smaller than $\int_4^\infty x^{-3}\,dx=1/32$, giving $\zeta(3)<2089/1728$. All six corner values are strictly negative. The basis weights are nonnegative and at least one product is positive when every wave is nonzero. Hence the reduced $C$ is strictly negative for both $(3,2,1,1,1)$ and $(2,2,2,1,1)$ at every nonzero real strength. It vanishes for $(4,1,1,1,1)$, for any zero wave, or for span below five. This sign refers to the specified reduced trace convention; the elementary differential also carries the factor $p_{\rm root}(v)$.

The sign extends to any nine conserved real waves with three collinear inputs $a_1u,a_2u,a_3u$. Denote the other six by $q_1,\ldots,q_6$. Each of the six determinants in a trace word contains at most one of the three collinear waves. Quadratic softness requires exactly two appearances of each, so every nonzero term factors out $(a_1a_2a_3)^2$. Its remaining factors are the four-row minors of the six momenta projected modulo $u$, using $\det(u,\cdot,\cdot,\cdot,\cdot)$ as volume. Let $d_{ij}$ be the square of the minor omitting slots $i,j$.

The reduced trace divided by $(a_1a_2a_3)^2$ is a linear combination of the following three symmetric polynomials:

$$
\sum_{\substack{\{e,f,h\}\ \mathrm{perfect\ matchings}}}d_e d_f d_h,\qquad
\sum_{\substack{(e,f)\ \mathrm{ordered\ disjoint\ pairs}}}d_e^2d_f,\qquad
\sum_{\substack{i<j<k\\\mathrm{six\ retained\ slots}}}d_{ij}d_{jk}d_{ik}.
$$

The three coefficients are the `three_collinear_sextics` entries in the [exact period table](elementary_scalar.json). To see completeness, the six projected waves sum to zero. At rank four their kernel has dimension two and contains the constant vector. The signed four-row minors are therefore differences $z_i-z_j$ for six numbers defined up to a common shift. Center the numbers to sum zero. Their symmetric degree-six polynomials have dimension four, with power-sum partitions $6,4+2,3+3,2+2+2$. Vanishing when five numbers coincide imposes at least one independent condition, since the sixth power sum is nonzero there away from the origin. The three displayed polynomials have independent values on three momentum fixtures. The trace values there follow from the fifteen exact graph measurements. Two further fixtures verify the polynomial row identity independently for all eighteen basis traces.

Order the six real numbers and take their consecutive gaps. Reflection permits the first gap to exceed or equal the fifth. Subtract the fifth from the first, giving five nonnegative independent variables. After expansion all 205 nonzero coefficients are strictly negative. The replay checks this with the rational upper and lower integral bounds for the tail of $\zeta(3)$ after four terms. The sixth powers of the second, third, fourth and fifth gap variables each have strictly negative coefficient. The sole remaining ray, where five numbers coincide, has zero trace. Thus $C<0$ exactly when $a_1a_2a_3\ne0$, the six projected momenta span four dimensions, and each projected momentum is nonzero. If a projected momentum vanishes, its coordinate unit vector is in the kernel, forcing the other five numbers to coincide. The converse follows from the same two-dimensional kernel. [Exact sign and interpolation replay](elementary_scalar.py)

### Two distinct collinear pairs in the three-loop trace

Suppose two disjoint pairs among the nine conserved momenta lie on distinct lines. At rank five choose one representative of each pair and three unpaired momenta as the frame. In this frame the first eight momenta are $e^1,se^1,e^2,te^2,e^3,e^4,e^5,w$. The ninth balances their sum. Each of the eighteen trace words is the sum over the $9!$ permutations of six signed five-row determinants. Group equal products of the 126 symbolic minors before expansion. The resulting $18$ polynomials have $1676$ distinct monomials in $s,t,w_1,\ldots,w_5$. Their coefficient matrix has rational rank $12$. Every coefficient row lies in the span of the fifteen exact graph measurements: solving on fifteen pivot columns and expanding the other three columns gives polynomial identities over the integers after clearing denominators. Thus both period coefficients are known throughout this locus. A change of frame multiplies the trace by its determinant to the sixth power. [Exact two-pair trace and independent numerical-input integer replay](two_pair_scalar.py)

For real momenta on this locus the reduced trace $C$ is strictly negative at rank five with all nine waves nonzero. It vanishes when a wave is zero or the rank is lower. To verify strictness, project the five unpaired momenta along the two pair directions into three dimensions. The two kernel vectors of their row matrix are the constant vector and the five cofactor coordinates, defined up to an affine change. Place their extrema at $0,1$ on the last two slots and order the other three. Their complement is an independent triple: a dependence supported on those three slots would give an affine function vanishing at both extrema. The five cofactor coordinates in the resulting frame are $-w_3,-w_4,-w_5,1,0$. Reflection lets the first of their four ordered gaps be at least the last. Subtract the last gap from the first, so the four new gaps range over the nonnegative cone after homogenizing their scale.

Pull out the squared strengths of all four waves in the two pairs. In each pair use the half-sum and half-difference of its strengths, and translate the corresponding component of $w$ by the pair sum times the mean of the four largest cofactor coordinates. In these variables the negative trace is quadratic in each resulting triple. The sign replay homogenizes to degree ten in the gaps, by a strictly positive factor. For both rational endpoints $1202056/10^6$ and $1202057/10^6$, which bound $\zeta(3)$ by the integral test after 256 terms, it supplies one positive semidefinite Gram matrix for each of the 286 gap monomials. It checks every principal minor of every nonzero block. Of the 286 matrices, 282 have rank nine and four have rank five at each endpoint. The pure tenth powers of the last three gaps have full rank. If only the first gap survives, the four largest cofactor coordinates coincide. The kernel of its Gram matrix consists of the four products involving only the half-sums and half-differences. On a decomposable input with nonzero pair strengths, this kernel is reached precisely when both translated components vanish. They are then minus the first two components of the ninth momentum, whose other three components already vanish. This gives the stated zero set. Four collinear waves on a single line always give zero: their minimum total degree is eight, and six determinants can use at most six such slots. The elementary-field differential includes the further factor $p_{\rm root}(v)$. [Exact chamber and Gram-matrix replay](two_pair_sign.py)

On the larger locus with a single collinear pair the eighteen-dimensional trace space restricts to dimension seventeen. A kernel element is the permutation sum of the product of six minors with zero-based row sets $01257,01367,01458,02358,15678,23467$. Every pair of the nine labels occurs together in at least one minor. Thus each summand vanishes whenever any two momenta are collinear. Every label occurs at least twice, so this word belongs to the quadratic-soft space. Its first unrestricted integer evaluation is $-98484480$. Seventeen single-pair fixtures give a minor equal to $164$ modulo 1009, proving the matching lower bound in characteristic zero. [Exact kernel coordinates and rank replay](operator/collinear_trace.py)

### Cubic minors and the four-collinear normal variation

The simultaneous real zero set of the eighteen trace words consists exactly of rank below five, a zero momentum, or four momenta on one line. For the converse choose five independent momenta. If three of the four remaining labels lie on a line, none of the chosen five lies on it. Exchange one of the three for a basis vector with nonzero coefficient in its basis expansion. Each line now contains at most two remaining labels, so the four can be split into two independent pairs. Extend each pair to a basis. The three basis minors are nonzero, and their row sets cover the nine labels. Their squared product summed over label permutations is strictly positive and belongs to the eighteen-dimensional quadratic-soft space. The forward implication follows by counting the eight required occurrences of four collinear labels in six minors.

Before symmetrizing, the covering products of three minors span a 5460-dimensional space on conserved momenta. To compute its character at a permutation, choose subsets of its cycles by inclusion–exclusion for the degree-zero slots. For the retained cycles of lengths $L$, take the Schur character of shape $(3,3,3,3,3)$ on their permutation representation minus the trivial representation. Its complete symmetric characters have generating function $(1-t)\prod_L(1-t^L)^{-1}$, with one factor per retained cycle. The five total-momentum relations give the numerator by the soft-ideal resolution used for the genus-three trace. The symmetric invariant Gram matrices have sixteen blocks. Their multiplicities give $68$ symmetric entries in all, with largest block $5$ by $5$. An independent five-variable Weyl extraction gives the same character on all thirty conjugacy classes. [Partition multiplicities and both character calculations](operator/cubic_soft_blocks.py)

At a full-rank configuration with four nonzero collinear waves write them $a_i\ell$ for $i=1,\ldots,4$. The other five waves project to a conserved simplex in the four-dimensional quotient by $\ell$. For a variation of wave $i$, express its projection divided by $a_i$ in this simplex, using the unique five coordinates $z_{ik}$ summing to zero. Let $D$ be the determinant of the frame formed by $\ell$ and four of the other waves. For a linear variation with parameter $t$, the constant and linear coefficients of the reduced trace $C$ vanish. Its quadratic coefficient is
$$
[t^2]C=D^6(a_1a_2a_3a_4)^2
\left(-\frac{740466719}{169474598092800}
+\frac{469153}{130767436800}\zeta(3)\right)
\sum_{\substack{1\le i<j\le4\\\rm collinear\ pairs}}
\sum_{\substack{1\le k\le5\\\rm simplex\ coordinates}}(z_{ik}-z_{jk})^2 .
$$
Each of the six determinants must use one unperturbed line slot, since at least eight appearances of the four waves must be filled. At order two this exhausts their minimum individual degrees and all available variations. The remaining entries therefore use only the projected simplex. After extracting the squared strengths and $D^6$, a common shift of the four slopes comes from moving the entire line, equivalently a shear of the frame. Permutation symmetry makes the quadratic a multiple of the displayed sum: the standard representations on the four centered slopes and the five centered simplex coordinates each have one invariant bilinear form. On $p_1=p_2=e^1,\ p_3=p_4=e^1+te^2$ with remaining frame $e^2,\ldots,e^5$, the sum is $16/5$ and the quadratic row is $768$ times trace coordinate eight minus coordinate two in zero-based order. A second test with unequal strengths, longitudinal shears and two independent variations reproduces all eighteen coordinates. The displayed scalar is strictly negative by $\zeta(3)<2089/1728$. Hence the full-rank, nonzero four-collinear locus has twelve negative normal eigenvalues and the Hessian vanishes on its tangent space. [Exact normal-variation replay](four_collinear_hessian.py)

### Global sign of the three-loop trace on real momenta

For all nine conserved real momenta, the reduced trace $C$ is strictly negative precisely when the span has dimension five, every momentum is nonzero, and no four are collinear. In every other case it vanishes. The elementary observable on the same waves has the additional signed linear factor $p_{\rm root}(v)$.

For a product of three five-row minors using every slot, sum its square over all $9!$ permutations. The sum of the ten such orbit sums, one for each incidence type, is positive off the common zero set proved above. At rank five with no zero wave or four collinear waves there is a covering product with all three minors nonzero. When four are collinear, six determinants can use at most six of their slots, whereas quadratic softness requires eight.

At each endpoint $1202056/10^6,1202057/10^6$ substituted for $\zeta(3)$, the negative trace equals $2^{-48}$ times this ten-orbit sum plus sixteen positive definite rational quadratic forms, each summed over permutation cosets. Each form uses explicitly listed cubic products of minors. In the row-wise numbering of its Young diagram, symmetrize a product over permutations within each row, then antisymmetrize over permutations within each column. Scale each resulting component by its recorded inverse power of two. Relabeling the input permutes the column cosets and multiplies every component on a coset by the same column sign. Its quadratic form summed over cosets is therefore permutation invariant. Every monomial in it carries at least two powers of every wave and exactly six determinants, so it belongs to the same eighteen-dimensional trace space.

The replay evaluates all forms on eighteen conserved integer fixtures where the original trace basis has determinant $999\bmod1009$. There are 68 independent symmetric matrix entries before imposing the 18 polynomial equations. Fifty rational entries are fixed in the certificate. The replay solves the eighteen remaining entries over the rationals at each endpoint, verifies all equations, and checks the positivity of every leading principal minor in the sixteen matrices. The integral bounds for the zeta-series tail after term 256 place $\zeta(3)$ strictly between the two endpoints. Linear interpolation of these identities proves the sign on every real configuration. [Integer and rational replay](operator/three_loop_sign.py) · [Orbit words, rational matrices and fixtures](operator/three_loop_sign_certificate.json)

### Sharp simultaneous soft order of the reduced trace

At total rank five let exactly $k$ labeled momenta vanish, and write each of them along a conserved variation as $t v_s+O(t^2)$. The $9-k$ nonzero hard momenta at $t=0$ span five dimensions and sum to zero, so $1\le k\le3$. The hard momenta may vary arbitrarily subject to conservation. The first possible power of the reduced trace $C$ is $t^{2k}$. Its coefficient is strictly negative exactly when each soft velocity $v_s$ is nonzero and the nine vectors formed by the hard momenta at zero and the soft velocities have no four on the same line. Otherwise that coefficient vanishes.

Indeed, the soft order in a six-minor word is at least the number of occurrences of its soft labels, hence at least $2k$. At order $2k$ every soft label must occur exactly twice. Every hard entry must use its constant value, since every extra hard variation increases the order. Every soft entry must use its first derivative. A zero velocity therefore kills the coefficient. Four collinear vectors in this combined list require at least eight entries on their line and cannot occupy the six determinants.

Conversely split the nine vectors into three independent subsets by the covering construction above, assigning every soft label to just one subset. All hard vectors together span five dimensions, so each subset can be extended to a basis using only hard labels. The corresponding product of three minors has leading order exactly $k$ with nonzero coefficient. Its orbit square supplies a strictly positive coefficient at order $2k$ in $-C$: the unrestricted positive Gram certificate retains $2^{-48}$ times the sum of all ten covering orbit squares, and every other square in that certificate starts at order at least $2k$. This proves strictness independently of the hard variations. For one soft wave the negative of its normal Hessian has rank five when the eight hard waves have no collinear triple. In the presence of such a triple the rank is four and its kernel is that line. Two distinct hard triples are impossible: the remaining two waves project to a conserved pair in the quotient by the two triple directions, giving total rank at most three. [Positive Gram certificate](operator/three_loop_sign.py)

For an exact check take the five coordinate covectors, $3-k$ additional copies of $e^1$, and their balancing hard momentum. Set every soft velocity to $e^1+e^2$. For $k=1,2,3$ the rational and $\zeta(3)$ coefficients at order $2k$ are respectively
$$
\left(-\frac{740466719}{70614415872000},\frac{469153}{54486432000}\right),\quad
\left(-\frac{103647659}{3851695411200},\frac{164077}{7429968000}\right),\quad
\left(-\frac{841479311}{35307207936000},\frac{53359}{2724321600}\right).
$$
Changing the common velocity to $e^1$ gives zero in all three cases. The evaluator selects degree two in each soft label before summing the permutations, with the balancing hard momentum fixed at its zero-parameter value. Interpolation along the conserved curves checks the resulting eighteen-coordinate rows at all three soft orders. [Arbitrary rational soft velocities and exact checks](soft_scalar_limit.py)

The same certificate bounds $-C$ between $2^{-48}$ and $2^{-19}$ times the sum of the ten cubic orbit squares. Indeed each row sum and signed column sum of a cubic has squared norm, summed over column cosets, at most its original orbit square times the column-group order and the square of the row-group order. This follows from Cauchy–Schwarz on the finite row and column sums. Apply Cauchy–Schwarz again to each Gram cross term, and bound each orbit square by the ten-orbit sum. The replay checks the resulting upper bound, including every inverse power of two in the cubic scalings, at both endpoints.

The holomorphic trace also has smooth complex zeros at full rank. In the two-pair frame take $s=t=1$ and $w=xe^1+e^2$, with $e^5$ on the vector and polarization $-e_4$, so $p_{\rm root}(v)=1$. The exact restriction is $C(0)+x^2(C(1)-C(0))$ with the two coefficients
$$
C(0)=-\frac{49188779}{2715939072000}+\frac{7783}{523908000}\zeta(3),
\qquad
C(1)-C(0)=-\frac{7896562031}{211843247616000}+\frac{2500061}{81729648000}\zeta(3).
$$
Both are strictly negative by the preceding upper bound on $\zeta(3)$. At $x=\pm i\sqrt{C(0)/(C(1)-C(0))}$ the derivative $2x(C(1)-C(0))$ is nonzero. Thus the holomorphic zero set is locally a smooth hypersurface. There are nearby full-rank configurations on it with all nine directions distinct: in the space of eight freely chosen complex momenta each specified collinearity has codimension four at nonzero momenta, including pairs involving the balancing wave by permutation symmetry. Their finite union cannot contain the hypersurface. [Symbolic two-pair restriction](two_pair_scalar.py)

Along any conserved real formal momentum curve, let $\operatorname{ord}_t$ give the first power with nonzero coefficient and assign $+\infty$ to the zero series. For a five-element set $I$ of labels, write $\det(p_I)$ for its five-row momentum minor in increasing label order. Then
$$
\operatorname{ord}_t C =
2\min_{\substack{|I|=|J|=|K|=5\\I\cup J\cup K=\{1,\ldots,9\}}}
\left(\operatorname{ord}_t\det(p_I)+\operatorname{ord}_t\det(p_J)
+\operatorname{ord}_t\det(p_K)\right).
$$
Every projected cubic in the sign certificate is a linear combination of these covering products, so its order is at least their minimum. The positive ten-square term contains the square of every such product up to permutation. At any finite minimum its leading coefficient is strictly positive, while the other positive Gram forms have nonnegative coefficients at that order. The first nonzero coefficient of $C$ is therefore strictly negative. If the minimum is infinite all covering products vanish identically, hence so does $C$. The same proof applies to convergent real analytic curves through their power series. Three successive minimizations over unions of label sets compute the exact result from the 126 orders, retaining a minimizing triple. For three hard momenta on $e^1$ and a fourth soft momentum $te^1+t^m e^2$ in the full-rank frame above, the order is $2m$ for $m\ge2$. Multiplying the fourth coordinate by $t$ and the fifth by $t^2$ adds $18$. [Exact polynomial-curve evaluator and minimizing minors](operator/trace_contact_order.py)

### Full three-loop trace modulo 1009

Use the two regular three-pair measurements followed by the coefficient quadratic in the marked soft strength. The last three right-hand sides in the rational and $\zeta(3)$ columns are $(91,492),(392,549),(967,830)$ modulo 1009, including the lower-genus hexagons. The eighteen-row determinant is $515$ in this field. Solving both columns shows that all eighteen trace coefficients have nonzero rational and irrational parts in the specified orbit basis. Each of the 153 two-row minors of the solution is nonzero modulo 1009. Consequently no two coefficients in this basis are rational multiples of one another. The integration denominator bounds and the unit determinant modulo 1009 justify lifting these nonvanishing statements to characteristic zero. [Both complete coefficient columns](operator/three_loop_reserved.json) · [Cubature and linear solve](trace_reconstruction.py)

### Integrating a fixed-time cubic-core numerator

After all ten labeled positions and their path assignments are integrated at fixed core resistances, the squared-cofactor voltage density is a symmetric homogeneous numerator of degree $30$ divided by $U^{12}$. Symmetry here means the core automorphism group, with external wave labels held fixed. For a term of exponent $k_e$ on core edge $e$, set $a_e=1+k_e$. Then $\sum_{e\in\mathrm{edges}}a_e=36$. For every nonempty proper edge subset the sum of its $a_e$ exceeds twelve times the number of independent cycles in that subset by monomial convergence. In particular $1\le a_e\le11$. These inequalities give $2540$ allowed monomial orbits on the tetrahedron and $841$ on the double square. For each orbit the rational and $\zeta(3)$ integration moments follow from the projective period times its orbit size.

At prime $1009$, evaluating the orbit sums on respectively $2540$ and $841$ nonsingular core-time configurations gives invertible matrices. Solving their transposes in the node-by-orbit convention gives the two integration weights at each node. All orbit moments replay. Symmetrizing an independently integrated ordered fivefold graph gives the same residue in both channels on each core. [Modular nodes, both weights and graph tests](core_cubature.py)

For exact reconstruction, the ordered path moment contains the denominator $\prod_{e\in\mathrm{core}} k_e!$. Therefore twice the automorphism order times the common denominator of the individual core moments divided by this factorial product bounds the denominator of each bare channel. The numerator bound uses the $24$-column linear current determinant before eliminating its constant columns. If $H_i$ is the sum of positive external entries in source component $i$, each particular tree flow is bounded in magnitude by $H_i$. Divide the three source rows of that component by $H_i$. At unit complex times the sixteen constant columns have squared norm at most $24$ and the eight variable columns squared norm at most $16^2\,8$. Parseval and the bound $\binom{23}{8}$ on the number of degree-eight monomials bound the squared coefficient one-norm by $\binom{23}{8}24^{16}16^{16}8^8\prod_{i\in\mathrm{components}}H_i^6$. There are at most $\binom{16}{4}$ separating two-forests for each pair. Multiplication by the absolute sum of its resistance weights, the maximal absolute core channel moment, and the number $10!\binom{15}{5}$ of ordered assignments per core gives an explicit reconstruction bound after division by the automorphism order. [Integer bounds for each measurement](period_bounds.py)

### Four-loop bare quadratic variation on the two-double core

Use the five normal momenta $(1,2,0),(3,-1,0),(-2,1,1),(0,2,2),(-2,-4,-3)$ and measure the fifth-component potential at the first normal vertex minus that at the second. The seven planar waves are three copies each of $e^1,e^2$ and one wave $t=t_1e^1+t_2e^2$. Vary the planar lifts $\ell_i$ linearly with $\lambda$ at fixed total. Their centered Gale contraction $Y$ varies as
$$
\partial_\lambda Y
=4\partial_\lambda\ell_0+8\partial_\lambda\ell_1
+14\partial_\lambda\ell_2-7\partial_\lambda\ell_3 .
$$
These coefficients annihilate the three normal rows. On the lift slopes $(1,-1,0,0,0)$ in the $e^1$ direction they give $\partial_\lambda Y=-4e^1$. Adding a normal row to either lift row leaves each signed current minor unchanged.

With two normal slots on each bubble and the fifth on a remaining path, the complete two-double-core contribution to the coefficient of $\lambda^2$ is
$$
\frac{367378814983}{7129340064000000}\,
\frac{36}{16}\left(4t_2^2(\partial_\lambda Y_1)^2+
4t_1^2(\partial_\lambda Y_2)^2-2t_1t_2(\partial_\lambda Y_1)(\partial_\lambda Y_2)\right).
$$
For the planar factor, let $M$ be the sum of the seven matrices $x_i p_i p_i^{\mathsf T}$ for the seven planar momenta $p_i$. Cauchy–Binet gives the coefficient of $x_1\cdots x_7$ in $(\det M)^3(\partial_\lambda Y)^{\mathsf T}\operatorname{adj}(M)(\partial_\lambda Y)/16$. For $t_1t_2\ne0$ the displayed real quadratic form has rank two on the eight-dimensional space of conserved lifts.

After summing the seven planar insertions, put $L=C\operatorname{diag}(T)C^{\mathsf T}$ for the cycle Gram matrix of the fourteen remaining segment times and $U=\det L$. Write $D_r$ for $U^4$ times the projective minor comprising the twelve normal current rows and the lifted current row of cycle $r$, evaluated on the normalized scalar lift $(1,-1,0,0,0)$. Write $V$ for $U$ times the observed potential difference. The remaining density is $V D^{\mathsf T}\operatorname{adj}(L)D/U^{13}$. At fixed exponents of the middle path the numerator matrices in the six or seven local variables on either side admit exact rational factorizations. Series and parallel beta integrals then integrate each side and the final three parallel arms. The seven placement orbits contain $16,8,4,4,8,4,1$ path multisets and $345$ inequivalent normal label orders, representing $5400$ ordered placements before division by the core automorphism order $16$. All seven orbit integrals are positive. The full twelve-wave Ward identity also involves the other occupancies, the remaining cores and the lower-vertex terms. [Exact shape periods and original determinant checks](operator/quartic_seven_plane_double.py)

When three normal slots occupy one bubble and two occupy the other, the same planar tensor has coefficient $-3161360627/475289337600000$. Adding the two possible occupancies gives the complete bare two-double-core coefficient
$$
\frac{367378814983}{7129340064000000}
-\frac{3161360627}{475289337600000}
=\frac{22854171827}{509238576000000}.
$$
For $t=(2,3)$ and $\partial_\lambda Y=-4e^1$ the planar factor is $1296$, giving $22854171827/392931000000$. [Independent integration of the three-plus-two occupancy](operator/quartic_seven_plane_double_three_two.py)

On the one-double core the two parallel paths, in series with their two connectors, contract to a tetrahedron edge. The other edge disjoint from it is the opposite edge. Place all five normal labels on one member of the parallel pair and on this opposite edge. For occupancies $(2,3)$ and $(3,2)$ the respective sums over all $120$ label orders, before the graph factor, are
$$
-\frac{7652526749891}{33713479800000}+\frac{1403408}{7432425}\zeta(3),
\qquad
\frac{14566533501799}{58998589650000}-\frac{4579976}{22297275}\zeta(3).
$$
Each placement has two choices of the parallel member and core automorphism order eight. The two sectors together therefore contribute $4698446757959/943977434400000-(92438/22297275)\zeta(3)$ times the planar tensor. For the five labels in the order $(0,1,2,3,4)$ on the first placement the $3584$ monomials cancel in tetrahedron exponent orbits. The order $(0,2,1,3,4)$ on the same placement gives $350728$ monomials and $185$ surviving exponent orbits. Cutting each tetrahedron along two edges gives the same coefficients. [Exact single-double integrals](operator/quartic_seven_plane_single.py)

At genus $n$ on $2n-1$ planar waves, form the $3\times4$ matrix $P$ of normal components on four independent conserved slots. The remaining current minor transforms by $\det(g)^n$ under $P\mapsto gP$, $g\in GL_3$. It is linear in a tangential lift row $\delta$ and invariant under adding any row of $P$ to $\delta$. Thus its four lift coefficients are proportional to the four signed maximal minors of $P$. These minors have greatest common divisor one, so after dividing by $\det(\delta;P)$ the quotient is a polynomial. It has determinant weight $n-1$. The determinant invariant theorem expresses it as a homogeneous degree-$n-1$ polynomial in the four maximal minors. Every three-vector in four dimensions is decomposable, so the four generators obey no Plücker relation. At most $\binom{n+2}{3}$ independent normal orderings determine the signed current on every label permutation at fixed geometry. The four-loop normal frame above has rank $20$. [Signed determinant interpolation and independent-order checks](operator/quartic_seven_plane_current_basis.py)

The tetrahedron integration after the one-double contraction has denominator power $13$. Its six integral powers sum to $39$, lie between $1$ and $12$, and sum to more than $13$ on each triangle. The $88100$ labeled possibilities give $3872$ tetrahedron automorphism orbits. Both complementary propagator cuts give the same rational and $\zeta(3)$ coefficients in every orbit. Exactly $1088$ orbits have a vertex whose three incident powers sum to at most $13$. Direct star–triangle expansion and rational beta integrals give their complete values. Each of the other $2784$ orbits has a nonzero $\zeta(3)$ coefficient. [Two-cut replay](operator/quartic_seven_plane_single.py) · [Exact moment coefficients](operator/quartic_seven_plane_tetra.weights.gz)

Suppose the one-double graph has a tetrahedron vertex outside the doubled macro-edge with three incident paths containing no normal insertion. Its quadratic lift period is rational for each retained normal ordering. On the three bare star edges set $t_i=T_j T_k/S$ for complementary indices, where $S=T_1+T_2+T_3$. The triangle edge of resistance $T_i$ connects the other two neighbors. The Jacobian is $T_1T_2T_3/S^3$ and eliminating the central potential gives $U_{\rm old}=U_{\rm new}/S$. The old currents have regular limits at a generic point of $S=0$ because each star current is the sum of the two outgoing new triangle currents at its neighbor. All observed potentials are on surviving vertices. Dropping a star column in the projective minor divides a bounded current determinant by $t_i$, so $D_r$ has a pole of order at most three. The voltage numerator has order at most one. The three star cycle columns, oriented outwards, sum to zero, so each entry of $\operatorname{adj}(L)$ has a pole of order at most two. The transformed density therefore vanishes to order at least $13-3-6-1-2=1$ on $S=0$. Its only candidate numerator denominator was a power of $S$. Thus its numerator is polynomial over $U_{\rm new}^{13}$.

The new graph is a triangle with two parallel networks on each side. Each network is series–parallel, and the convergent projective monomials have rational beta integrals. For the ordered path multiset $(0,0,0,4,6)$, summing all $120$ retained label orders gives $-238826266733/22279187700000000$. The replay checks the exact divisions by $S$, its remaining factor in the numerator, and equality of the pure beta integral with the two tetrahedron cuts. [Source-free star integration](operator/quartic_seven_plane_single_star.py)

For the five normal momenta and normalized lift above, all $36$ single-double orbits give the coefficient of the same planar tensor
$$
-\frac{10215037955553251914183}{6684978482582400000000}
+\frac{4684553723887133}{3684401721000000}\zeta(3).
$$
The $171$ path multisets give $20520$ ordered placements before division by $8$. Independently these are $2550$ labeled graph isomorphism classes with trivial automorphism and $30$ with automorphism order two. Exactly the six source-free star orbits have zero $\zeta(3)$ coefficient. Twenty-six orbit amplitudes are positive and ten negative. Adding the two-double result and the zero three-double result, the vanishing of the complete lift variation requires the joint primitive prism and bipartite contribution
$$
\frac{10214737939698010154383}{6684978482582400000000}
-\frac{4684553723887133}{3684401721000000}\zeta(3).
$$
It lies strictly between $-348566265409/10^{15}$ and $-348566265408/10^{15}$. This last number is inferred from the complete Ward identity. [All shape coefficients](operator/quartic_seven_plane_single_totals.json) · [Direct single-double period integration](operator/quartic_seven_plane_single_projection.py)

For the sparse normal frame $q=((1,0,0),(0,1,0),(0,0,-1),(-1,0,1-c),(0,-1,c))$, keep the same lift slopes and normal polarization. The complete single-double readings at $c=2,3$ are, in order,
$$
-\frac{7162011046352686943}{41928185042756812800000000}
+\frac{1642429680839}{11554283797056000000}\zeta(3),
\qquad
-\frac{12785599230834755447}{8385637008551362560000000}
+\frac{117291814373}{92434270376448000}\zeta(3).
$$
Together with the preceding frame these give all three quadratic normal tensor coefficients of the transverse pair measurement. In their degree-seven generator basis the determinant of the two-double rational, single-double rational and single-double $\zeta(3)$ columns is $253\bmod1009$. The required prism-plus-bipartite coefficient vector has signs $(-,+,+)$. Its two exact coefficient columns are recorded with the $108$ direct single-double shape and frame integrations. [Complete coefficients and conditional primitive targets](operator/quartic_seven_plane_single_totals.json) · [Basis reconstruction and sign checks](operator/quartic_seven_plane_single_projection.py)

Write $R(c)$ for the sum of the three reducible-core readings in this frame, with the planar tensor stripped as above. Restricting the three normal generators gives
$$
R(c)=c(c-1)(2c-1)\left[
\frac{R(3)-5R(2)}{120}\left(c-\frac12\right)^2+
\frac{125R(2)-9R(3)}{480}
\right].
$$
Indeed the restrictions contain the common factor $c(c-1)(2c-1)$, with quotients affine in $c(c-1)$. Eight integer substitutions verify these polynomial identities through degree seven. The exact periods satisfy $R(3)-5R(2)>3/10^{10}$ and $125R(2)-9R(3)>1/10^9$. Thus the reducible part has precisely the three real simple zeros $0,1/2,1$ and has the sign of $c(c-1)(2c-1)$ elsewhere. All five normal momenta are nonzero at these three full-rank configurations. The complete Ward identity fixes the opposite signed primitive contribution throughout this curve. [Polynomial restriction and original four-loop replay](operator/quartic_seven_plane_single_projection.py)

On the primitive prism with paths ordered $(01,02,03,12,14,25,34,35,45)$, the complete $120$-ordering sums for $(3,3,3,7,7)$, $(3,5,5,5,7)$ and $(3,3,3,5,7)$ have placement multiplicities $12,6,12$ and core automorphism order $12$. For each placement the stars at $0$ and $4$ contain no normal source. Eliminating both gives a tetrahedron on $(1,2,3,5)$ with pairs of parallel lines on $12,13,35$. For each star sum $S$, the two currents, voltage and Gram adjugate have total pole order at most $3+3+1+1=8$. The transformed numerator therefore has a factor $S^{13-3-8}=S^2$ after including the Jacobian. Every parallel pair of powers $p,q$ contributes $\mathrm B(13-p,13-q)$ and leaves power $p+q-13$. The three integrals use $742,1003,2014$ tetrahedron exponent orbits on the default normal frame. All three periods are strictly negative. Their sum, including the graph factors and stripping the same planar tensor, is
$$
\frac{2095589196653513523367}{13339581119865000000000}
-\frac{20591840995151}{157544183671875}\zeta(3).
$$
It is strictly between $-19487463182/10^{15}$ and $-19487463181/10^{15}$. The three orbits account for $3600$ ordered placements. Subtraction from the joint Ward target leaves $77$ prism and $26$ bipartite orbits whose required total is strictly between $-329078802228/10^{15}$ and $-329078802227/10^{15}$. Independent grounded-current and ordered-simplex evaluations on nine rational fixtures check both star Jacobians on all three orbits. Replay with `python operator/quartic_seven_plane_prism.py --all` and `python operator/quartic_seven_plane_prism_check.py`. — Codex. [Direct star integrals](operator/quartic_seven_plane_prism.py) · [All three normal-frame readings](operator/quartic_seven_plane_prism_totals.json) · [Independent current checks](operator/quartic_seven_plane_prism_check.py)

The other two prism placements admitting the same pair of empty stars are $(3,3,5,5,7)$ and $(3,3,5,7,7)$, with multiplicities $12$ and $6$. Their default $120$-ordering sums before the graph factor are
$$
\frac{1030266064282981614304633}{14818368969724320000000000}
-\frac{6749236524005117}{116672721165000000}\zeta(3)<0,
\qquad
-\frac{17729247232415445359727989}{55568883636466200000000000}
+\frac{3251659963394214649}{12250635722325000000}\zeta(3)>0 .
$$
The exact reductions visit $2654$ and $2550$ tetrahedron moment orbits. Grouping the current products by voltage coefficients before expanding the position polynomials retains all $120$ orderings. Splitting the largest product at fixed degree $m_4+2m_{11}$ in its two indicated segment exponents bounds the memory. The five directly integrated orbits comprise $48$ prism placements and $5760$ assignments. Their graph-weighted total is
$$
\frac{104397999393562651831321499}{1555928741821053600000000000}
-\frac{1368122818212743441}{24501271444650000000}\zeta(3)
\ \in\ \frac1{10^{15}}(-24572272817,-24572272816).
$$
The same four negative and one positive readings occur on both sparse frames. Subtracting the computed default value from the Ward target leaves the required sum on the other $75$ prism and $26$ bipartite orbits in $(-323993992593,-323993992592)/10^{15}$. All five representations use only the rational and weight-three columns of the tetrahedron cuts. [Paired-orbit calculation](operator/quartic_seven_plane_prism_paired.py) · [Fifteen direct frame readings and reconstructed tensor](operator/quartic_seven_plane_prism_totals.json)

The occupied-star obstruction is visible before pairing the new edges. Fix triangle resistances $(2,3,S-5)$, the six other path times $(7,11,13,17,19,23)$, and normal order $(0,2,1,3,4)$. After integrating all five positions, the transformed numerator has order $2$ at $S=0$ for marks $(3,3,3,7,7)$, order $-6$ for $(0,3,3,7,7)$ and order $-10$ for $(0,1,3,7,7)$. Multiplying by $S^{33}$ gives the factor $6(S-5)$ times the degree-$43$ original path numerator at affine times. Exact interpolation at $45$ integer points determines this polynomial and an additional point checks it. The occupied cases leave an additional denominator in the triangle resistance sum. [Ordered-position calculation with `--occupied-star`](operator/quartic_seven_plane_prism_check.py) — Codex.

On the sparse normal frame $q(c)$ above, each of the five symmetrized orbit contributions obeys the same degree-seven restriction already proved for $R(c)$. For the four placements with a negative reading at $c=2$, both $R(3)-5R(2)$ and $125R(2)-9R(3)$ are strictly negative. For $(3,3,5,7,7)$ both are strictly positive. Thus the first four orbits have the opposite sign to the fifth for every real $c$ outside the three common simple zeros $0,1/2,1$. The graph-weighted five-orbit sum has sign opposite to $c(c-1)(2c-1)$. The normalized coefficients of these signs follow by rational bounds on $\zeta(3)$ applied to the two exact sparse readings in each orbit. [Exact sign checks](operator/quartic_seven_plane_prism.py)

### A Ward projection for the mixed genus-four D3 pole

Let $C_g(v)$ denote the leading coefficient $[t^{2g-4}]V^{(4)}_{g;0}$ on the four hard covectors $e^1,\ldots,e^4$, the $g-2$ identical soft covectors $tv$, and the balancing covector. Use eight fourfold waves with zero-based labels in the order
$$
e^1,\ e^2,\ e^3,\ e^4,\ e^4,\ -e^1-e^2-e^3-2e^4-2tv,\ tv,\ tv,
\qquad v=e^3-e^2.
$$
For the physical potential difference $(A_0-A_1)(e_4)$, the fusion weight of partner $j$ is $p_j(e_4)$. Both soft partners drop out. Fusing the balancing partner leaves two identical $e^4$ and rank at most three at $t=0$, so each of the eight four-minors contains a factor of $t$. For each $e^4$ partner the four remaining hard rows form a simplex of absolute determinant one. The coordinates of $v$ are, up to simplex permutation, common shift and overall sign, $(1,-1,0,0,0)$ for fusion with leg zero and $(1,-1,-1,0,0)$ for fusion with leg one. The symmetric quartic interpolation on the five simplex coordinates gives, respectively,
$$
4C_4(e^1)+C_4(e^1+e^2),\qquad
4C_4(e^1)+\frac{11}{3}C_4(e^1+e^2).
$$
Consequently the order-four coefficient of the tree-channel potential, with one classical cubic vertex and one genus-four trace, is $-16 C_4(e^1+e^2)/3$. Direct polynomial evaluation checks it on all $41$ square orbits and both signed trace words. For the interpolation, every four-row minor is affine in $t$ since the varying rows are parallel to $v$. Each retained trace is divisible by $t^4$, so five nonzero values recover its leading coefficient. [Trace basis](operator/fourfold_quartic_basis.py) · [Exact fusion replay](operator/fourfold_mixed_soft_ward.py)

At total genus four the linked pair of genus-two action vertices has zero $[t^4]$. If both soft legs meet one vertex, that determinant vanishes. If the soft legs split, each vertex is a fourth power of a determinant divisible by $t$. The primary genus-one fourfold vertex vanishes on conserved momenta by the rank bound. This leaves the cubically generated four-loop term, the two-loop theta with a genus-two vertex, and the pentagon with a genus-three vertex. In the Li–Wang flux sign convention they obey
$$
-\frac{16}{3} C_4(e^1+e^2)
+[t^4]\!\int_{\mathrm{bare}}\!\Phi_H
+[t^4]\!\int_{\mathrm{theta},\,V_{2;0}}\!\Phi_H
+\frac{740466719}{23832365356800000}
-\frac{469153}{18389170800000}\zeta(3)=0 .
$$

For the pentagon ground its six-point genus-three vertex. The other four vertices have one external wave each. Summing their orders fills the four position intervals with reflection factor $1/2$. If both soft legs are singletons, their determinant vanishes. If both are at the genus-three vertex, its other four hard momenta have rank at most three at $t=0$ by conservation, and all six four-minors at that vertex are divisible by $t$. Thus exactly one soft leg occupies a singleton. The two assignments cancel the reflection factor. The remaining three hard singleton positions are $s_i$; at $t=0$ the current into the grounded vertex is $\sum_i s_i p_i$. For the five hard rows at this vertex take its first three external waves and this current as a basis. Write $D$ for their determinant and $D_i$ for the determinant with row $i$ replaced by $v$. The genus-three soft numerator in units of $C_3(e^1)$ is
$$
\frac{D^4}{4}\left(5\sum_{i\in\mathrm{basis\ rows}}D_i^2-
\left(\sum_{i\in\mathrm{basis\ rows}}D_i\right)^2\right).
$$
Multiplying by the squared singleton determinant with its soft column replaced by $v$, and integrating the potential difference with kernel $\min(s_i,s_j)-s_is_j$, gives $4C_3(e^1)/225$. Independent integration of the full genus-three trace before coefficient extraction gives the same rational and $\zeta(3)$ terms in the Ward equation. [Both pentagon integrations](operator/fourfold_mixed_soft_ward.py)

For the theta term, a soft leg at the genus-two vertex gives a fourth power of $t$. A second soft leg there kills the vertex determinant, and an ordinary second soft leg gives a further factor $t^2$ in the squared graph cofactor. Both soft legs must therefore be ordinary at order four. Project the remaining four geometric events by $q\mapsto(q_1,q_2+q_3,q_4)$, with kernel $\mathbb Z v$ on the momentum lattice. Inserting the two labeled soft legs on segments $i,j$ multiplies the squared projected cofactor by $t^4\det(C_i,C_j)^2$, where $C_i$ is the two-component cycle column of segment $i$. With segment lengths $\ell_i$ the placement sum gives
$$
2\sum_{i<j\,;\,\mathrm{segments}}\ell_i\ell_j\det(C_i,C_j)^2
=2\det\!\left(\sum_{i\,;\,\mathrm{segments}}\ell_i C_i C_i^T\right).
$$
This cancels one Kirchhoff denominator. The genus-two vertex still uses the full four-dimensional current determinant to the fourth power divided by $62208$. There are $\binom62$ choices of hard inputs at a distinguished endpoint, and $\binom63$ at an interior vertex with endpoint reflection factor $1/2$. Integrating the ordered positions and conductance simplex gives, respectively,
$$
\frac{896261}{3802314700800000}
+\frac{1}{23702740992}
=\frac{3170033}{11406944102400000}.
$$
Thus the only unevaluated diagram sum in this Ward component is the bare four-loop integral,
$$
C_4(e^1+e^2)=
\frac{3}{16}[t^4]\!\int_{\mathrm{bare}}\!\Phi_H
+\frac{30426931}{5176678632652800}
-\frac{469153}{98075577600000}\zeta(3).
$$
The replay restores all $21$ distinct segment pairs at two rational positions and strengths for each vertex location. Independently the same integral, with the three-collinear test momenta, recovers the branch and interior normalization separately. [Theta integration and current checks](operator/mixed_soft_theta.py)

For each bare cubic graph delete the two soft bivalent vertices. The remaining six-mark graph has fifteen segment lengths $\ell_i$, cycle columns $C_i\in\mathbb Z^4$ and hard currents $j_i$. At one soft insertion the difference of its two segment columns is exactly $t(0,C_i\otimes v)$. The coefficient of $t^2$ in the signed current cofactor therefore uses the fifteen columns $(1,C_i\otimes j_i)$ evaluated at zero and these two jump columns divided by $t$. It is independent of both positions on the segments. Define the skew matrix $K_{ab}$ by appending $(0,e_a\otimes v),(0,e_b\otimes v)$ to the fifteen hard columns and dividing the determinant by $\sum_i\ell_i$; here $e_a$ labels the coordinate vectors in the cycle space. For $M=\sum_{i\in\mathrm{segments}}\ell_i C_i C_i^T$, the complete square after integrating the two positions is
$$
2\sum_{i<j\,;\,\mathrm{segments}}\ell_i\ell_j(C_i^T K C_j)^2
=-\operatorname{tr}(MKMK).
$$
Equal segments have zero determinant. The hard voltage is evaluated at zero and the Kirchhoff denominator is unchanged by subdivision. The replay checks all $105$ distinct segment pairs on a bipartite core, restores the full seventeen edges at two rational position pairs, and extracts the signed cofactor directly. This leaves a six-hard-mark projective integral with six determinants per metric. [Placement sum and restored rational fixtures](operator/fourfold_mixed_soft_ward.py)

The two identical hard $e^4$ charges can also be integrated before the remaining four hard marks. On retained segments with cycle columns $C_e,C_f$, their jumps occupy two fourth-coordinate rows. Proportional columns give zero. The other two fourth-coordinate rows of the current determinant are affine in both fractional positions, since each charge splits between its endpoints with weights $1-s,s$ or $1-t,t$. The voltage compares two retained vertices and is also affine. The signed square times voltage has total degree at most five, so each position monomial integrates to $1/((i+1)(j+1))$ for powers $i,j$. On the bipartite fixture with retained paths $(0,1,4,8)$, $74$ of $78$ unordered segment pairs survive. A six-by-six charged-graph quadrature independently reproduces each surviving pair. [Four free insertions and the signed fixed-metric value](operator/mixed_free_pair.py)

For all six hard marks, replace the segment columns on each core path by their resistance-weighted mean and the successive hard jumps. The absolute Jacobian of this column operation is one. The six jump columns and the two soft columns are constant. The nine mean columns depend affinely on the six hard positions since a charge at $s$ splits between its endpoints with weights $1-s,s$. Subtracting one mean column leaves an eight-by-eight determinant of affine entries. At voltage insertion $i$, a source on the same path contributes the additional interval Green function $T_e(\min(s_i,s_j)-s_i s_j)$ times its fourth charge. The two required moments are
$$
\int_0^1\!\int_0^1s^k t^l(\min(s,t)-st)\,ds\,dt
=\frac{1}{(k+2)(l+2)(k+l+3)},\qquad
\int_0^1 s^k(s-s^2)\,ds=\frac{1}{(k+2)(k+3)}.
$$
These integrate all eight insertions at fixed core lengths with unrestricted hard positions. At nine unit lengths the complete placements $(0,0,5,5,8,8)$ and $(0,0,0,5,5,8)$ give respectively $-133651127/125547836992825968000$ and $-7130789987/1883217554892389520000$ before the factor $U^{-4}$. The $(0,0,1,1,3,3)$ value is zero at this metric. Direct segmented graphs, ordered simplices and the original ordered numerator check the cube formula. [Unordered position integrals](operator/mixed_unordered_hard.py)

Summing all $54$ bipartite placements at the unit metric covers $2575$ generically nonzero ordered color orbits. The two identical fourth-coordinate charges may be exchanged. Exchanging the two $T$ charges sends the soft target vector to its negative, and exchanging the two $U$ charges is a unimodular target reflection fixing that vector. Their squares agree, so only the voltage needs to be averaged over these exchanges. The complete bare density at this point, including $U^{-4}$ with $U=81$, is
$$
\frac{5033220989519}{20266585166938719172690980000}.
$$
The original ordered six-variable numerators reproduce every placement independently modulo $1009$. [Exact unit-metric sum by placement](operator/mixed_bipartite_unit_metric.json)

The five cubic cores have unit-length Kirchhoff determinants $45,36,56,75,81$ and contribute $797,0,3995,11732,2575$ generically nonzero color orbits. Integrating all eight insertion positions and including each graph's factor $U^{-4}$ gives the following bare density at nine unit lengths:
$$
\frac{9277624090747753060520411797721155138060662396539}
{1795300005626402189032326209293740672000000000000000000000000000}>0.
$$
The ordered numerators reproduce all $377$ placement sums modulo $1009$. The projective length integrals and lower-loop Ward terms still enter the four-loop contact. [Five core values](operator/mixed_core_unit_metric.json)

For the complete placement $(0,0,5,5,8,8)$ fix eight core lengths at one and integrate the length $x$ of path $1$. Here $U=45+36x$. The numerator of the density with denominator $U^{13}$ has degree eight and exactly two positive zeros, one in $(0,1/4)$ and one in $(3/2,2)$. Although the density at $x=1$ is negative, its integral over $0<x<\infty$ is
$$
\frac{8626038015317273}{361577770539338787840000000000000000}.
$$
For a fixed placement with $m$ hard marks on the varying edge, a grounded mean current has an affine numerator over $U$. The squared eight-row cofactor therefore contributes at most degree sixteen over $U^{16}$. The compound soft Gram has degree at most one since the varying edge contributes a rank-one matrix. The voltage contributes at most degree two over $U$, and the hard measure contributes $x^m$. Multiplication by $U^{17}(x+8)^2$ before the final factor $U^{-4}$ gives a polynomial of degree at most $19+m$. Rational interpolation and exact division give the density above. Each term $x^k/(45+36x)^{13}$ integrates to $45^{k-12}36^{-k-1}k!(11-k)!/12!$. [One-edge calculation and root counts](operator/mixed_unordered_hard.py)

Two core-length integrations can also be performed on the complete placement numerator averaged over core automorphisms. Set the other seven lengths to one. If $x,y$ are the two integration variables, the Kirchhoff polynomial is $U=a+bx+cy+dxy$. A numerator monomial $x^iy^j$, with $0\le i,j\le11$, contributes
$$
\int_0^\infty\!\int_0^\infty \frac{x^iy^j\,dx\,dy}{U^{13}}
=\frac{i!(11-i)!}{12!}\int_0^\infty\frac{y^j\,dy}{(a+cy)^{12-i}(b+dy)^{i+1}}.
$$
Partial fractions give a rational number plus a rational multiple of $\log(bc/ad)$. Reversing the integration order reproduces the two coefficients separately. On prism placement $005566$, all $20$ colors and $1728$ nonzero labeled orders combine into $716102$ scalar-power orbits. The unit density is $269424934217/1131510970652103424072265625000$. In the scalar propagator edge order, the two-length integral on pair $(0,1)$ lies in $(1172,1173)10^{-20}$, while pair $(0,6)$ gives $(-478,-477)10^{-20}$. The respective logarithm arguments are $152/77$ and $49/24$. For bipartite matching $005588$ both adjacent and disjoint edge pairs give negative integrals in $(-398,-397)10^{-20}$ and $(-810,-809)10^{-20}$, respectively. Its complete projective period is positive as evaluated below. [Two-length coefficients and exact interval check](operator/mixed_two_edge_marginal.py)

For a connected bridgeless subgraph $H$ of a cubic core, let $h$ be its cycle rank and $c_i$ the restriction of the cycle column on hard path $i$ to the cycles supported on $H$. Use the same restriction for the two soft paths. After subtracting neighboring segment columns at all marks, there are $|E(H)|$ undifferenced columns on $H$. Kirchhoff's voltage law gives one relation with positive path lengths as coefficients on those columns. Rank $4h$ on the supported cycle rows therefore requires
$$
\operatorname{rank}\big(c_i\otimes p_i\mid i\text{ runs over all eight marks}\big)
\ \ge\ 4h+1-|E(H)|,
$$
with both soft momenta replaced by $v$. A single pair of soft paths must pass this test for every $H$. In the three-double core each of its three disjoint double edges needs three insertions, exceeding the eight available.

If $H$ has $b$ boundary vertices and its hard momenta span dimension $r$, its hard currents lie in a space of dimension at most $r+b-1$ for $b>0$. Its $b$ incoming boundary currents sum to the prescribed internal source. If this bound misses two target dimensions, the common soft direction cannot fill them. If it misses one, $v$ must lie outside the hard momentum span, and the two soft columns restricted to $H$ must span all $h$ cycle directions.

The six hard labels give $9^6$ path assignments on each core. The first three counts successively impose the segment, tensor and boundary conditions with a common pair of soft paths. The rightmost counts include all orders along each path, quotient by the oriented core and three label exchanges, and apply the two-edge cut test below.

| Cubic core | Segment count | Tensor rank | Boundary bound | Chamber orbits after cuts | Certified nonzero integrands |
|---|---:|---:|---:|---:|---:|
| Two double edges | 84,320 | 35,072 | 29,760 | 821 | 797 |
| Three double edges | 0 | 0 | 0 | 0 | 0 |
| One double edge | 215,810 | 102,696 | 92,264 | 4,173 | 3,995 |
| Prism | 511,230 | 392,076 | 381,132 | 11,732 | 11,732 |
| $K_{3,3}$ | 519,300 | 464,040 | 449,064 | 2,575 | 2,575 |

Swapping hard labels $0,5$ acts by the determinant-minus-one involution exchanging $e^1$ with the balancing covector and fixing $e^2,e^3,e^4$. Swapping $1,2$ reverses $v$, and the waves at $3,4$ are identical. All three exchanges preserve the squared soft cofactor. Let $R_{ij}$ be the effective resistance between hard sites $i,j$. For the pairs $U=\{0,5\}$, $T=\{1,2\}$ and $Z=\{3,4\}$, the average of $(A_0-A_1)(e_4)$ over these eight permutations is
$$
R_{05}+\frac14\left(
\sum_{i\in T,\,j\in Z}R_{ij}
-\sum_{i\in U,\,j\in Z}R_{ij}
-\sum_{i\in T,\,j\in U}R_{ij}\right).
$$
The first term includes the covariance of the positive root with the opposite balancing source.

Fix the two $Z$ sites. The four other hard momenta $e^1,e^2,e^3,-e^1-e^2-e^3-2e^4$ can be permuted by unimodular maps fixing $e^4$. Each transposition has determinant $-1$, so the squared cofactor is invariant when the two soft insertions transform with the hard momenta. In a fixed hard frame the six assignments therefore use the six soft directions $p_i-p_j$. Write $W_{ij}$ for the preceding mean voltage with $T=\{i,j\}$ and $U=\{k,l\}$ its complement among the four sites. Then
$$
\begin{aligned}
W_{ij}+W_{kl}&=R_{ij}+R_{kl}-\tfrac12(R_{ik}+R_{il}+R_{jk}+R_{jl}),\\
W_{ij}-W_{kl}&=R_{kl}-R_{ij}+\tfrac12\sum_{\substack{z\in Z\\\mathrm{hard\ sites}}}
(R_{iz}+R_{jz}-R_{kz}-R_{lz}).
\end{aligned}
$$
Summing the first line over the three pairings gives $\sum_{\substack{i<j\\\mathrm{four\ hard\ sites}}}W_{ij}=0$ at every positive metric and every hard position. Thus the part of the soft square invariant under all four-site permutations drops out before integration.

Let $D_{ij}$ be the six signed soft cofactor minors at direction $p_i-p_j$, with $\|D_{ij}\|^2$ the integrated soft square from the cycle Gram matrix. For each complementary pairing the contribution is
$$
\tfrac12(W_{ij}+W_{kl})(\|D_{ij}\|^2+\|D_{kl}\|^2)
+\tfrac12(W_{ij}-W_{kl})(\|D_{ij}\|^2-\|D_{kl}\|^2).
$$
Eliminating the six common hard jump columns reduces all directions to the same $45$ maximal minors of a $10\times8$ matrix. At unit $K_{3,3}$ core lengths the two terms, integrated over all marked positions and summed over all $54$ placements, give respectively
$$
-\frac{37649396752847}{941608777446194760000}
+\frac{3181055915459}{62773918496412984000}
=\frac{5033220989519}{470804388723097380000}.
$$
The separately labeled calculation agrees at all $54$ placements. On the six-path placement $(0,1,2,3,4,5)$, doubling the first core length gives $171325529103593/66366466912649939490000$ by both constructions. The metric integration for this placement remains open. [Hard permutation calculation and unit-metric replay](operator/mixed_bipartite_hard_s4.py)

For a two-edge cut, scale all lengths on each separated component independently and vary each of the $k$ connector segment lengths independently. Only the one cycle passing between the components imposes a voltage condition. Let $Q$ be the span of the connector marks and the total hard momentum on the first component. All connector currents lie in $Q$ plus the line of one connector current. The voltage drops internal to the two components agree modulo the actual connector currents by the crossing-cycle equation. If all hard sources on either component belong to $Q$, its voltage drop already belongs to $Q$ plus the connector line. With $\epsilon=1$ if both components have an individual hard source outside $Q$, and zero otherwise, the $k+2$ length variations have at most $\min(4,1+\dim Q+\epsilon)$ independent conditions. When this is at most $k$, there is a second length variation independent of overall scale preserving all currents. The fifteen hard columns then have rank at most fourteen.

For each positive-integer metric the augmented fifteen hard columns must have rank $15$ and their images in the quotient by the four soft cycle directions must have rank $13$. Together these ranks imply a strictly positive integrated soft square. At the first test metric the ranks and nonzero resistance average hold modulo $1000003$ for $19\,099$ chamber orbits, proving a nonzero rational integrand on each. [Rank census, exact witnesses and chamber weights](operator/mixed_bare_ranks.py)

The remaining $24+178$ chambers vanish identically. Scale every length in a selected two-terminal block by a common factor, leaving every other segment independently variable. For an effective branch $e$ oriented from $a$ to $b$, let $s_i$ be the harmonic measure of endpoint $a$ at hard mark $i$. Eliminating the branch deposits $S_e=\sum_{i\in e}s_i p_i$ at $a$ and $\sum_{i\in e}p_i-S_e$ at $b$. Its initial current is the effective conductance times the endpoint voltage difference minus $S_e$. On an ordinary branch all currents are spanned by its voltage difference and its hard momenta. On a branch containing a scaled block, the independent length variations give columns in the span of its endpoint voltage difference, $S_e$, its hard momenta outside the block and the total momentum inside the block. Indeed each outer segment current differs from the initial current by the cumulative source, and the block voltage equals the full voltage minus the drops on those outer segments.

Multiply these columns by the corresponding cycle column of the reduced graph. Express the endpoint potentials by applying the symmetric grounded Green matrix to the deposited sources. Taking all six $s_i$ and all Green entries as independent indeterminates gives a polynomial matrix whose rank bounds the physical rank at every positive metric. On a theta reduction the common Green scalar can be set to one because it rescales only voltage-difference columns. Fraction-free elimination over the integer polynomial ring gives the following ranks. The numbers of length variations include the independent scale of each selected block.

| Original core | Reduced core | Chambers | Length variations | Polynomial rank |
|---|---|---:|---:|---:|
| Two double edges | Theta, one five-path block | 14 | 8 | 6 |
| Two double edges | Theta, both double edges scaled separately | 10 | 8 | 6 |
| One double edge | Theta, one five-path block | 80 | 8 | 6 |
| One double edge | Tetrahedron, one double edge scaled | 98 | 11 | 9 |

In each case the kernel has dimension at least two. Fixing total length leaves a nonzero current-preserving variation, so the fifteen hard columns are dependent. This gives exactly $797,0,3995,11732,2575$ nonzero mixed bare integrands in core order. Their projective integrals remain to be evaluated. [Polynomial certificates](operator/mixed_bare_symbolic.py)

One colored chamber on the two-double core can be integrated completely. Orient the nine paths as $01,01,02,13,23,24,35,45,45$, numbered from zero in this order. Place hard labels $0,1$ in order on path $0$, labels $2,3$ on path $4$, and labels $4,5$ on path $7$. Each bubble requires a third mark by the cycle-rank test, so the two soft marks must occupy opposite bubbles. For soft labels $6,7$ on left and right respectively, the four signed contributions on path pairs $(0,7),(0,8),(1,7),(1,8)$ are
$$
\frac{2741}{24638999261184000000}(4,2,2,1).
$$
With both soft assignments the fixed hard order contributes $2741/1368833292288000000$. Averaging the voltage over the eight hard-label exchanges and restoring the orbit weight $8$ gives the contribution of this chamber orbit to the bare Ward integral,
$$
[t^4]\!\int_{\text{one chamber orbit}}\!\Phi_H
=\frac{858337}{167682078305280000000}.
$$
For paths $p,q$ in opposite bubbles, let $T_p,T_q$ denote their total lengths, $U$ the Kirchhoff polynomial, and $V$ the voltage numerator. Multiplying the signed projective cofactor by $U^4$ gives a degree-fifteen polynomial $D_{pq}$. The position-integrated density is $T_pT_q D_{pq}^2V/U^{13}$. At fixed powers of the three central segments its numerator splits into products of the two bubble polynomials. Series and parallel beta integrals evaluate both bubbles and the remaining three parallel arms. The original grounded currents, all soft placements, the eight voltages, and the elementary projective monomial formula give independent checks. [Exact chamber replay](operator/mixed_bare_double_chamber.py)

For the full two-double family with two hard marks in each bubble and both remaining hard marks on one of the five other paths, the only surviving placement of the last pair is the central path $23$. For example, if the pair is on $02$ or $13$, take the five-path subgraph $23,24,35,45,45$. Its cycle rank is two; the signed current cofactor requires at least $4\cdot 2+1-5=4$ marks on it. It contains only two hard marks and the one soft mark forced onto its bubble. Exchanging the two bubbles proves the exclusion on the other connectors.

On $23$ the three possible splits of the two hard marks between the parallel paths in each bubble give the following signed contributions. The counts include the orders of all six labeled hard marks along their paths. Each integral includes the two soft assignments and division by the core automorphism order $16$.

| Left and right bubble splits, up to interchange | Hard chambers | Positive, negative, zero colored orbits | Bare Ward integral |
|---|---:|---:|---:|
| $(2,0),(2,0)$ | $2880$ | $9,9,9$ | $-12917/226842638400000000$ |
| $(2,0),(1,1)$ | $2880$ | $8,9,13$ | $-258047/4990538044800000000$ |
| $(1,1),(1,1)$ | $720$ | $3,3,6$ | $-29399/3208203028800000000$ |

Together with the four connector placements, the family has $45\cdot6!=32400$ ordered hard chambers in $300$ colored orbits. Its complete fourth-order soft coefficient is
$$
[t^4]\!\int_{\text{two marks on one remaining path}}\!\Phi_H
=-\frac{211663}{1796593696128000000}.
$$
Direct enumeration of the $32400$ labeled chambers fixes every orbit multiplicity. For every representative, summing the eight hard-label assignments as formal coefficients of the pair resistances reproduces the signed voltage used in the integral. [All shape subtotals, sign counts and replay](operator/mixed_bare_double_family.py)

#### Complete two-double contribution to the mixed Ward projection

For the measured voltage $(A_0-A_1)(e_4)$, group the hard chambers by the two bubble populations in decreasing order. Each chamber weight includes the two soft-label assignments, the eight hard-label exchanges and the order-$16$ cubic-core automorphism factor.

| Hard marks on the bubbles | Retained colored chambers | Integrated fourth soft coefficient |
|---|---:|---:|
| $(4,2)$ | $69$ | $-3761867/2096025978816000000$ |
| $(3,3)$ | $68$ | $-7280081/50304623491584000000$ |
| $(3,2)$ | $412$ | $-20278742261/1006092469831680000000$ |
| $(2,2)$ | $272$ | $101611318673/3018277409495040000000$ |

Their sum and the corresponding term in the Ward expression for the D3 bulk coefficient are
$$
[t^4]\!\int_{\text{two-double core}}\!\Phi_H
=\frac{99774853}{8623649741414400000},\qquad
\frac{3}{16}[t^4]\!\int_{\text{two-double core}}\!\Phi_H
=\frac{99774853}{45992798620876800000}.
$$
Inserting this sum, both theta placements, and the pentagon into the mixed pole gives
$$
C_4(e^1+e^2)=\frac{3515607603589}{597906382071398400000}
-\frac{469153}{98075577600000}\zeta(3)
+\frac{3}{16}[t^4]\!\int_{\text{one-double, prism and bipartite cores}}\!\Phi_H .
$$
The segment count leaves $2541$ colored chambers. The six signed cofactor polynomials vanish identically on $1744$. Of the remaining $797$ integrated chamber weights, $456$ are positive and $341$ negative. Proportional skew minors share one integrated square after adding their path-length factors. All $36$ path pairs are checked against the original grounded current cofactor at lengths $1,\ldots,15$ on each retained chamber. An independent integration on hard paths $(0,0,0,1,7,7)$ with color word $(U,U,Z,T,T,Z)$ expands the full fifteen-segment numerator. Its $1218668$ monomials give $8773$ nine-path moments and the same signed chamber value $9409/5589402610176000000$. The one-double and two primitive cores still enter the four-loop Ward projection. [Every signed chamber and multiplicity](operator/mixed_bare_double_core.json) · [All integrals, polynomial zero checks and direct monomial replay](operator/mixed_bare_double_core.py)

#### Tetrahedron periods for the mixed single-double core

The single-double core reduces to the $Q=13$ tetrahedron periods. Each signed soft-pair cofactor is a degree-$15$ numerator divided by $U^4$, where $U$ is the four-loop Kirchhoff polynomial. The two-forest numerator of the signed voltage has degree five and denominator $U$. The two placement lengths have degree $2$, and the physical measure gives a further $U^{-4}$. The complete fifteen-segment numerator therefore has degree $2\cdot15+2+5=37$ over $U^{13}$. At the first segment length $-866351$ and the other fourteen lengths $66184j$ for $j=2,\ldots,15$, $U=0$ but twice the signed numerator of the first chamber below, before its symmetry weight, is $897822\bmod1000003$.

Absolute convergence holds before soft extraction: on the unit projective simplex the voltage and the second current cofactor are bounded, while the first cofactor times $U^{-4}$ is the absolute Jacobian of the bounded current regions. Extracting a polynomial coefficient by finite interpolation and integrating the two positions preserves absolute convergence. On scaling any proper nonempty segment subset, the leading coefficient of the numerator is nonzero on an open cone unless that entire face of its polynomial support is absent. Thus every surviving monomial obeys the strict cycle-rank inequalities. [Signed current coordinates and convergence](../nonSP/verify/k4_periods/README.md#schwinger-numerator-and-convergence)

For a fifteen-segment monomial with exponents $m_i$, sum $1+m_i$ along the two parallel paths to obtain $a,b$, and along the two adjacent connector paths to obtain $c$. The strict inequalities give $a,b<13<a+b$. Local beta integration gives the effective tetrahedron edge power $a+b+c-13$. The other five powers are their respective path sums of $1+m_i$. The six powers sum to $39$, each lies between $1$ and $12$, and each triangle sum exceeds $13$. The $3\,872$ existing convergent tetrahedron orbits therefore give every one of the $3\,995$ single-double chamber periods in $\mathbb Q+\mathbb Q\zeta(3)$. [Local beta integrals and tetrahedron table](operator/quartic_seven_plane_single.py)

Exactly seventeen of the certified chamber orbits place all six hard marks on only two core paths. Their common path multiset is $(0,0,0,8,8,8)$. For the color order $(U,U,Z,T,T,Z)$, the labeled multiplicity divided by the core automorphism order is $8$. The full chamber and the full seventeen-orbit subtotal, respectively, are
$$
-\frac{2113894596943}{4795707439531008000000}
+\frac{989}{2697078384000}\zeta(3)
\ \in\ \frac1{10^{20}}(-294068,-294067),
$$
$$
-\frac{65542158940903}{19982114331379200000000}
+\frac{53663}{19666196550000}\zeta(3)
\ \in\ \frac1{10^{20}}(223747,223748).
$$
For the first chamber the full soft placement sum has one primitive twelve-term determinant. Independent expansion gives $151\,767$ degree-$37$ monomials and $185$ tetrahedron exponent orbits. Integrating this expansion gives the same rational and $\zeta(3)$ terms as the factorized calculation. Opposite kite cuts reproduce the first and last of the $185$ exponent orbits. [Seventeen exact integrations and both checks](operator/mixed_bare_single.py)

If vertex $4$ or $5$ of the original single-double core has no hard mark on any of its three incident paths, the chamber period is rational. Under the star substitution $t_i=T_j T_k/S$, $S=T_1+T_2+T_3$, the old Kirchhoff polynomial becomes $U_{\rm new}/S$ and the Jacobian has denominator $S^3$. Each old star current is a sum of two triangle currents and has a regular limit at generic $S=0$. Omitting a bare star column in the projective cofactor divides a bounded determinant by $t_i$, so each cofactor numerator $U^4D$ has at most a third-order pole. The voltage numerator has at most a simple pole, and the two placement lengths have at most a double pole. The transformed numerator has order at least $13-3-2\cdot3-1-2=1$ at $S=0$. The remaining triangle with parallel networks has only rational beta integrals. [Star beta integral](operator/quartic_seven_plane_single_star.py)

The fourteen hard shapes contain $542$ chamber orbits with total labeled multiplicity $34\,112$ and core automorphism order $8$. Of their exact signed periods, $299$ are positive and $243$ negative. The beta calculation and the independent tetrahedron cuts agree for each chamber. Summing with all label and graph factors gives
$$
[t^4]\!\int_{\text{single-double, empty star}}\!\Phi_H
=-\frac{1151829647875967}{38030295359637504000000000}.
$$
After including the complete two-double core the signed rational contribution to the Ward expression for $C_4(e^1+e^2)$ from these two sectors is
$$
\frac3{16}\left(
\frac{99774853}{8623649741414400000}
-\frac{1151829647875967}{38030295359637504000000000}
\right)
=-\frac{711822546145967}{202828241918066688000000000}.
$$
Other hard shapes on the single-double core and both primitive cores still enter this Ward component. [Fourteen exact subtotals](operator/mixed_bare_single_star.json) · [Classifier and independent beta replay](operator/mixed_bare_single_star.py)

Classify each placement by the two hard populations on the parallel paths and the combined population on the two adjacent connectors. The other marks occupy the five exterior paths. With the two parallel populations unordered, two complete classes give
$$
\begin{aligned}
[t^4]\!\int_{(1,1,1)}\!\Phi_H
&=\frac{8403930341168089853}{300167688374281728000000000}
-\frac{1100062224257}{47267524624320000000}\zeta(3)
\ \in\frac{(2185103968,2185103969)}{10^{20}},\\
[t^4]\!\int_{(3,0,1)}\!\Phi_H
&=-\frac{579577849255781221}{5395170869472384000000000}
+\frac{37948781861}{424789845480000000}\zeta(3)
\ \in\frac{(-3903433512,-3903433511)}{10^{20}}.
\end{aligned}
$$
The first class contains $15$ placements, $369$ nonzero color orbits and $22032$ labeled hard orders. The second contains $5$, $231$ and $14688$, respectively. The $720$ signed hard orders determine the coefficients of all fifteen two-point resistances. The integrator collects proportional signed cofactors and sums their resistance polynomials before taking the parallel-edge beta integral. Clearing the common denominator of the local moments leaves integer contractions on the five exterior paths. The grouped coefficients agree with the separate color integrations on all $51$ orbits of $(0,0,0,2,4,7)$ and all $15$ of $(0,1,2,8,8,8)$. [Grouped integrals](operator/mixed_bare_single_sum.jsonl) · [Replay and overlap check](operator/mixed_bare_single_sum.py)

Among the $3995$ nonzero chamber integrands, the two connector paths carry at most one hard mark in total. Summing the $(1,1,1),(2,0,1),(2,1,1),(3,0,1)$ classes covers every such insertion. Their $40$ placements carry $1527$ nonzero color orbits and $95472$ labeled orders. The resulting signed bare contribution is
$$
[t^4]\!\int_{\text{single-double, occupied connector}}\!\Phi_H
=-\frac{162208752620670243793}{660368914423419801600000000}
+\frac{4250166819953}{20797710834700800000}\zeta(3)
\ \in\frac{(1572661674,1572661675)}{10^{20}}.
$$
For hard color word $(U,U,Z,T,T,Z)$, direct grounded electrical solutions and explicit multiplication of the exterior polynomials reproduce the contractions on $(0,1,2,8,8,8)$, $(0,1,2,4,5,8)$ and $(0,0,0,4,7,8)$. The last test requires seven distinct soft cofactors. Run `mixed_bare_single_sum.py --controls` from the operator directory. [Projection checks](operator/mixed_bare_single_sum.py)

Summing all ten parallel-arm and connector populations gives the complete single-double core on the mixed Ward fixture,
$$
[t^4]\!\int_{\text{single-double core}}\!\Phi_H
=-\frac{772225330124291741}{1310255782586150400000000}
+\frac{16993070109403}{34662851391168000000}\zeta(3)
\ \in\frac{(-7513250331,-7513250330)}{10^{20}}.
$$
The $102$ placements contain $3995$ nonzero color orbits and $248928$ labeled orders. The individually stored periods close $36$ placements. The grouped sums close $68$, with two fully overlapping placements giving exactly the same period in both calculations. On $(0,0,0,2,4,7)$, combining all $51$ color orbits before the exterior moment contraction reduces the number of polynomial products from $735$ to $101$ and their moment rows from $576857$ to $85001$. On the dense $(0,0,1,4,6,8)$ chambers the six signed minors can be retained while the soft squares are integrated individually: each square carries positive path lengths and the same signed voltage. Sparse integer row elimination factorizes its bubble and tetrahedron variables. [Individual periods](operator/mixed_bare_single_core.jsonl) · [Grouped periods](operator/mixed_bare_single_sum.jsonl) · [Exact coverage and bounds](operator/mixed_bare_single_sum.py)

Including both reducible cores, the theta and the pentagon leaves
$$
C_4(e^1+e^2)=-\frac{1462273332514690607}{13976061680918937600000000}
+\frac{16108735470523}{184868540752896000000}\zeta(3)
+\frac{3}{16}[t^4]\!\int_{\text{prism and bipartite cores}}\!\Phi_H .
$$
The fixed term lies in $(11563257517,11563257518)/10^{20}$. The two primitive terms in this equation include the subfamilies evaluated below.

The complete mixed four-loop contact satisfies
$$
C_4(e^1+e^2)\in\operatorname{span}_{\mathbb Q}\{1,\zeta(3),\zeta(5)\}.
$$
For each ordered hard chamber the fifteen-segment numerator has degree $37$ over $U^{13}$. Its Kirchhoff polynomial depends only on the nine core lengths. Integrating the six hard positions on their ordered simplices raises the numerator degree to $43$ with rational coefficients. A resulting monomial of degrees $k_e$ has positive scalar powers $a_e=k_e+1$ summing to $52$. Absolute convergence enforces each strict subgraph inequality: scale the corresponding lengths together and restrict their ratios and the other lengths to an open set on which the lowest homogeneous numerator is nonzero. The $1386282$ positive bipartite power orbits and $2591953$ completed positive prism classes cover every primitive monomial. Both tables have rational columns on $1,\zeta(3),\zeta(5)$. The three completed bipartite placements with determinant $474$ modulo $1009$ below give the reverse inclusion for the rational span of all individual bare chamber periods. Thus the $19099$ color periods span exactly this space. The numerical coefficients of their full weighted sum still require the remaining physical numerators. The axial contact is in the same space, so the symmetric quartic formula gives the same membership for $C_4(v)$ on every rational soft covector $v$. [Prism domain and integration](operator/mixed_prism_scalar_functional.py) · [Bipartite domain](operator/mixed_bipartite_star_bank.json) · [Axial contact and quartic interpolation](../twisted-holography-N4/README.md#complex-tangent-momenta-at-which-the-local-pole-vanishes)

#### Primitive mixed chambers with Kirchhoff power thirteen

On $K_{3,3}$ all six-hard-mark placements allowing two successive empty core star–triangle moves have zero mixed Ward integrand. After the first move its three neighbors have degree four, so both centers belong to the same bipartition class. The union of their incident paths is $K_{2,3}$, with six edges and two independent cycles. It contains zero hard marks. Even if both soft marks lie there, it has at most eight segments. In the four-dimensional current matrix the two supported cycles supply eight rows. The vector of positive segment lengths annihilates those rows by Kirchhoff's voltage law, so their rank is at most seven. This holds for every order of the hard labels. The six hard marks must all lie on the complementary three-edge star. There are $6\binom82-9=159$ such uncolored occupancies and seven core orbits, indexed by partitions of six into at most three parts. Exact grounded-current computations give zero on one fixture in each orbit. [Census and current replay](operator/mixed_bare_ranks.py)

For the nonzero bipartite chamber with hard paths $(0,0,0,5,5,5)$ and colors $(U,U,Z,T,T,Z)$, the empty star on paths $(4,6,8)$ gives a triangular bipyramid. Write $m_i$ for a transformed monomial exponent and $b_i=12-m_i$ for its planar dual prism power. The $105218$ monomials have $0<b_i<13$ and total power $52$. The projective coefficient of each normalized prism cut includes $\prod_{i\in\mathrm{paths}}(b_i-1)!/12!$. The minimum pair or triple collision margin among all their six-point completions is one. Completing the graphs gives $47678$ cuts, of which $47676$ have nonzero net weights. Four-point reflections leave $37391$ nonzero weighted components. Within these, $586$ components admit a scalar star and a positive tetrahedron chart. Their subtotal is
$$
\frac{13780571003319104322144109}{29035118039206639794610176000000000}
-\frac{1842108151731299}{4665481145845648128000000}\zeta(3).
$$
The table retains the weights of all $37391$ components. Its columns are nine prism powers, the chamber weight, a zero second channel, and a tensor-bubble flag. For the $586$ scalar reductions seven more columns give the six tetrahedron powers and the normalized star factor. At two unequal metrics the transformed integrand including its Jacobian agrees modulo $1009$ with direct position integration. At $q=6$ a separate $1481$-graph test reproduces all three independently evaluated period coefficients and checks all $36$ scalar-star reductions. [Ward chamber replay](operator/mixed_bipartite_one_star.py) · [Weighted components](operator/mixed_bipartite_one_star.components.gz) · [Three-column controls](quartic_prism_completion.py)

All $1242$ further tensor-bubble components also integrate exactly. Send a positive neighbor to infinity, combine the two remaining denominators of powers $b,c$ with parameter $t$, and set $M=tA+(1-t)B$ for their two endpoints. For the numerator polynomial $P$, the Gaussian moment of order $j$ is $(\Delta/4)^j P(M)/j!$. For a term of this moment carrying $t^u(1-t)^v$, the scalar factor from the radial and parameter integrals is
$$
\frac{\Gamma(b+c-q-j)}{\Gamma(b)\Gamma(c)}
\mathrm B(q+j-c+u,q+j-b+v)
\bigl(|A-B|^2\bigr)^{q+j-b-c}.
$$
The negative incident powers in the selected charts have total degrees $1,2,3,4,5$, in respectively $689,343,146,53,11$ components. Recompletion leaves convergent tetrahedron or two-denominator integrals. All three choices of infinity agree on every component, giving $3726$ contractions. At $q=6$ all $61$ bubble components of the independent bank agree in every period column. The combined $1828$ scalar and tensor components contribute
$$
\frac{25095295652441839438597}{62374045197006744993792000000000}
-\frac{383542229667277}{1145907649856825856000000}\zeta(3)
\ \in\ \frac{1}{10^{20}}(-17682,-17681).
$$
This leaves $35563$ weighted components after bubble contraction. [Exact contraction and controls](operator/mixed_bipartite_tensor.py)

Of these components, $10940$ have a four-point reflection with two disjoint missing pairs. Sending one of the other two vertices to infinity and fixing the other at the origin leaves a four-wheel with positive powers. Orthogonality sets the four rim angular momenta equal. At $q=13$ the existing angular kernel takes $68$ samples per wheel. Reconstruction over $16$ prime fields gives the subtotal
$$
\frac{612680153354189825651120954254291}{3629389754900829974326272000000000}
+\frac{247723220447053749700079347}{548660582751448219852800000}\zeta(3)
-\frac{2783935140019231}{4056996571992000}\zeta(5)
\ \in\ \frac{1}{10^{20}}(-39162,-39161).
$$
All three coefficients agree with independent propagator integrals on $199$ completed components at $q=6$ and with direct rational partial fractions on twelve wheels at $q=13$. There are $24623$ components beyond these scalar wheel cuts. [Angular sums, chamber weights and independent checks](operator/mixed_bipartite_wheels.py)

A diagonal numerator $|x_0-x_2|^2$ also integrates explicitly. Its radial terms $r_0^2+r_2^2$ preserve the common angular degree. The term $-2r_0r_2\cos\theta$ couples degree $n$ on one half of the rim to $n+1$ on the other. Relative to the scalar angular factor, each of the two assignments has multiplier
$$
-\frac{(n+q-1)(n+2q-2)}{(n+q)^2}.
$$
For rim powers $6,9,11,12$ this integrates ten further chamber components. Their signed subtotal is
$$
-\frac{15469511507765307278569}{146309488733719525294080000000}
-\frac{5575670452211}{19748068342203801600}\zeta(3)
+\frac{579197}{1349265677760}\zeta(5).
$$
Direct partial fractions reproduce all three independent propagator coefficients on $613$ signed charts at $q=6$. [Two angular degrees and exact checks](operator/mixed_bipartite_linear_wheel.py)

Across the four-point reflections, $10583$ further components have a diagonal of power $-1$ and convergent separated radial terms. This includes the ten just displayed. Their total is
$$
-\frac{21349904253811884693303442588528537}{116140472156826559178440704000000000}
-\frac{269760817092032201695082203}{548660582751448219852800000}\zeta(3)
+\frac{2934584497151629337}{3927172681688256000}\zeta(5)
\ \in\ \frac{1}{10^{20}}(-47677,-47676).
$$
For angular degree $n$ the rational summand has poles at $n=-1,\dots,-2(q-1)$ of order at most six. Each radial cut crosses two rim lines. Its shift and the shift for the complementary cut with the two angular degrees interchanged add to $2q-1$ after reversal of the Gegenbauer shifts, so positivity puts every pole within the stated range. The angular factor cancels one pole from the four connections and three radial denominators away from the central degrees. At each central degree only two shifted connections can have poles. Absolute convergence follows from the two scalar radial terms bounding the cosine term and positivity of each projected radial kernel. The summand is proper, so $12(q-1)$ evaluations determine it. For $q=13$ the $144$ evaluations over $19$ prime fields certify all seven partial-fraction channels, including the zeros. For the common denominator, expand each connection into integer simple-pole terms. The angular factor and the three radial denominators leave at most nine inverse factors per term. Splitting the poles and summing their finite harmonic parts uses at most nine powers of $\operatorname{lcm}(1,\ldots,2q-3)$. The numerator bound uses circles of radius $1/2$ around the poles. The $4484$ direct prism representatives agree separately with their reflected values, and the three period columns agree with the independent $q=6$ propagators on $613$ direct charts and $243$ reflected seeds. In this one-star chamber $37391-1828-10940-10583=14040$ components remain. [Integer reconstruction and checks](operator/mixed_bipartite_linear_bank.py)

For a diagonal of power $-2$, expand $|x_0-x_2|^4$. Besides the scalar terms and the neighboring degrees above, the angular degrees on the two sides of the rim may differ by two. Write $\lambda=q-1$ and $n$ for the lower Gegenbauer degree. For the $4\cos^2\theta$ term the relative factors at equal degrees and at degrees $(n,n+2)$ are respectively
$$
\frac{(n+1)(n+2\lambda)}{(n+\lambda)(n+\lambda+1)}
+\frac{n(n+2\lambda-1)}{(n+\lambda)(n+\lambda-1)},
\qquad
\frac{(n+\lambda)(n+2\lambda)(n+2\lambda+1)}
{(n+\lambda+1)(n+\lambda+2)^2}.
$$
They follow by applying the Gegenbauer three-term recurrence twice. The unequal-degree coefficient includes the ratio of the two angular norms. The two extreme radial terms $r_0^4,r_2^4$ pass all pair- and triple-collision inequalities in every selected chart. For unequal radii the angular expansion converges absolutely. Each rim kernel has positive Gegenbauer coefficients by the Euler integral for its terminating connection polynomial. Multiplication by $\cos\theta$ has nonnegative adjacent-degree coefficients. The sum of the absolute angular contributions at fixed radii is therefore bounded by the integral of the positive rim kernels with $|\cos\theta|$ replaced by one and the numerator bounded by $16(r_0^4+r_2^4)$. This proves convergence of the radial integrals and the angular sum.

For a degree jump $j=0,1,2$, reversing the two halves and the connection shifts pairs the radial pole locations of complementary cuts to sum $2\lambda+j$. Positive cut margins place the poles in $1,\ldots,2\lambda+1$. The four connections have simple poles. Including the three radial denominators and the displayed angular factors gives order at most seven. The convergent rational summand is recovered from $7(2q-1)$ consecutive degrees. Expanding each connection into simple poles with integer residues leaves at most thirteen inverse linear factors per term. Thus $\lambda!^3(2\lambda-1)!\operatorname{lcm}(1,\ldots,2\lambda+1)^{13}$ clears all period coefficients. Half-unit circles about the integer poles give the reconstruction height bound used in the replay.

After the scalar and linear sectors have been removed, $5903$ reflected components admit this quadratic cut. Their signed contribution is
$$
\frac{2946614573584234833492017791125139}{58070236078413279589220352000000000}
+\frac{20683963022424458481150037}{152405717430957838848000000}\zeta(3)
-\frac{810035414794129859}{3927172681688256000}\zeta(5)
\ \in\ \frac{(-6146,-6145)}{10^{20}} .
$$
All eight Laurent-sum coefficients were reconstructed at $175$ degrees in $24$ prime fields. The coefficients at weights $1,2,4,6,7$ vanish separately on every selected cut. The three surviving columns agree with independent propagator data on $230$ direct charts and $113$ reflected seeds at $q=6$. Six $q=13$ charts also agree with direct Laurent summation. This reduces the residual count in the chamber to $14040-5903=8137$. [Exact quadratic angular sum and controls](operator/mixed_bipartite_quadratic_bank.py)

For a numerator $|x_0-x_2|^{2d}$, choose $p$ angular factors and $i$ powers of $r_0^2$ from the remaining $d-p$ factors. The radial powers on the two endpoints are $2i+p$ and $2d-2i-p$, with coefficient $\binom dp\binom{d-p}i$. Each multiplication by $-2\cos\theta$ takes a Gegenbauer polynomial of degree $l$ to degrees $l+1,l-1$ with respective coefficients
$$
-\frac{l+1}{l+\lambda},\qquad -\frac{l+2\lambda-1}{l+\lambda}.
$$
The term below degree zero vanishes. After $p$ steps from $n$ to $n+j$, multiply the resulting coefficient by the angular norm ratio
$$
\left(\frac{n+\lambda}{n+j+\lambda}\right)^2
\prod_{i=0}^{j-1}\frac{n+2\lambda+i}{n+1+i}.
$$
For $j>0$ both assignments to the two half-rims enter. Combining terms at each radial degree and spin jump leaves precisely $(d+1)^2$ radial integrals per chart. The same two extreme radial monomials dominate the angular sum. For the height bound the replay expands all $p$-step transfers into products of integer linear factors and cancels matching factors before counting poles. Complementary radial cuts bound the pole locations by $2\lambda+d-1$. The bank evaluates only the intermediate radial sums actually used by the selected charts.

| Diagonal degree $d$ | Additional reflected components | Angular evaluations per prime | Prime fields | Strict bounds for $10^{20}$ times the subtotal |
|---|---:|---:|---:|---:|
| $3$ | $2583$ | $182$ | $22$ | $(-954,-953)$ |
| $4$ | $710$ | $216$ | $24$ | $(-557,-556)$ |
| $5$ | $75$ | $224$ | $26$ | $(-67,-66)$ |

Direct Laurent sums at $q=13$ reproduce a complete three-column period in each row. At $q=6$ the cubic and quartic reductions agree with $57$ and $6$ direct charts and $39$ and $3$ independently evaluated reflected seeds. Including the bubbles and the scalar, linear and quadratic cuts, the disjoint $32622$-component subtotal on this color chamber is
$$
\frac{223034075807229941528734717303}{1009917149189796166769049600000000}
+\frac{1753641788846304650723}{2968942547356321536000000}\zeta(3)
-\frac{73446572862343}{81816097535172000}\zeta(5)
\ \in\ \frac{(-112241,-112240)}{10^{20}}.
$$
The remaining $4769$ components have no positive-rim numerator chart with both separately convergent extreme radial terms anywhere in their four-point reflection orbits. The degree at each completed vertex is $2q$ and the other four incident powers are at most $q-1$. Hence $d\le 2q-4=22$. Every degree $6,\ldots,22$ returns zero further convergent charts. Dropping the separate convergence test leaves $929$ formal numerator charts, first appearing in degrees $1,\ldots,6$ with counts $(87,232,271,218,109,12)$. The other $3840$ components admit no such polynomial diagonal. [All three coefficients and exact replays](operator/mixed_bipartite_numerator_bank.py) · [Radial and Gegenbauer recurrences](quartic_wheel_bank.cpp)

#### A common angular regulator for all remaining polynomial diagonals

The 929 formal charts above have a single failed collision inequality in either extreme radial monomial. Its complementary set consists of the four rim vertices meeting at an arbitrary common nonzero point. Every radial subcluster at the origin or infinity passes the strict inequalities. The unsplit diagram passes all pair and triple tests in its six-point completion.

Smooth the propagator on rim edge zero by the spherical Poisson kernel with parameter $0<t<1$. This multiplies its degree-$m$ harmonic by $t^m$. Write $\rho$ for the ratio of the two radii, $b$ for the rim power and $\theta$ for the angle between the endpoints. Removing the outer radial power leaves $((1-\rho)^2+2\rho(1-\cos\theta))^{-b}$. Its Poisson average is bounded uniformly in $t$ by a constant times this same function. For $\rho\le1/2$ the ratio of the maximum to the minimum is bounded. For larger $\rho$, the angular factor is comparable, with uniform constants, to $((1-\rho)^2+\theta^2)^{-b}$. On a ball of angular radius $R<\theta/2$ this follows from the pointwise bound. For $R\ge\theta/2$ the ball lies in the ball of radius $3R$ about the singularity. Integrating there and dividing by the volume of the original ball gives a bound by a constant times $((1-\rho)^2+R^2)^{-b}$, since $2b\le2q-2<2q-1$, the angular dimension. The positive decreasing Poisson kernel is a weighted average of these ball means. Thus dominated convergence on the unsplit numerator gives the physical period as $t$ increases to one.

For a Gegenbauer path joining lower degree $n$ to $n+j$, the regulated spin $m$ equals $n+j$ when edge zero carries the higher degree and equals $n$ otherwise. At fixed $t$ each individual radial moment converges. There are finitely many connection coefficients and their radial moments grow at most polynomially in $m$, so the damped sum is absolutely convergent. Collect both spin assignments at each $m$ and integrate the first $d$ shells directly. Above them the integrated shell is a rational function of $m-d$. Its Abel sum is finite by the preceding bound. A polynomial part would give a power divergence, and a nonzero simple harmonic tail would give a logarithm. Hence the collected tail decays at least as $1/m^2$ and its ordinary sum gives the period. The translation of each rational denominator by $d$ or $d-j$ increases the maximal pole location by at most $d$. The same lattice bound from the preceding transfer calculation applies with this enlarged range to the tail and to the $d$ finite shells.

The additional counts in degrees one through six are $(87,232,271,218,109,12)$. The three Laurent-sum columns agree with the independently dimension-shifted $q=6$ propagators on all $424$ available reflected seed measurements, including $39$ whose separated monomials diverge. At $q=13$, direct rational partial fractions reproduce one full three-column check in each degree. The disjoint $33551$-component subtotal in the $(0,0,0,5,5,5)$ chamber is
$$
\frac{1875840319482653063081961833881}{8933882473602043013726208000000000}
+\frac{192571615355536381505743}{342912864219655137408000000}\zeta(3)
-\frac{119708316672311}{140256167203152000}\zeta(5)
\ \in\ \frac{(-198448,-198447)}{10^{20}} .
$$
There are $37391-33551=3840$ unresolved components in this chamber. [Collision census and fixed-spin replay](operator/mixed_bipartite_fixed_spin.py) · [Angular shell evaluation](quartic_wheel_bank.cpp)

#### Complete mixed bipartite period for the color word $(U,U,Z,T,T,Z)$

For a six-point completion with edge powers $a_{ij}$ and a missing pair, take $v$ as an integration point, $s$ as reference and $h$ at infinity, with $s,h$ outside the missing pair. Integrate the divergence with vector field $x_v-x_s$. Let $I$ denote its normalized period. In each summand $I_-$ and $I_+$ have the powers on $vj,sh$ raised by one. For $I_-$ lower $vs,jh$ by one. For $I_+$ lower $js,vh$ by one. The identity is
$$
0=(a_{vh}-a_{vs})I+
\sum_{\substack{j\ne v,s,h\\\mathrm{incident\ edges}}}a_{vj}(I_+-I_-).
$$
Indeed $2(x_v-x_s)\cdot(x_v-x_j)=|x_v-x_s|^2+|x_v-x_j|^2-|x_j-x_s|^2$ and the sum of the incident powers at $v$ is $2q$. Only identities whose shifted integrals all obey $a_{ij}<q$ and $a_{ij}+a_{ik}+a_{jk}<2q$ are used. Hold the other points distinct. At a finite neighbor $j$ the flux across a sphere of radius $r$ is bounded by a constant times $r^{2q-1-2a_{vj}}$ for $0<r<1$. The exponent is positive since the powers are integers. The large sphere flux is bounded by a constant times $r^{-2q+2a_{vh}}$, which tends to zero. The absolute convergence of every shifted period justifies the remaining integrations.

Generating the identities on the reflection orbits of the $3840$ unresolved seeds gives $344019$ candidate equations. Each time one unknown period remains in an equation, solve it over the rationals and propagate that value into the other equations. The $4020$ necessary substitutions terminate on $5196$ wheel charts. This finishes all $37391$ components of the $(0,0,0,5,5,5)$ color word $(U,U,Z,T,T,Z)$:
$$
[t^4]\!\int_{(U,U,Z,T,T,Z)}\!\Phi_H
=\frac{213249810833125731948743}{733743238248244224000000000}
+\frac{17962743011421743}{23108567594112000000}\zeta(3)
-\frac{859794979}{727783056000}\zeta(5)
\ \in\ \frac{(-206524,-206523)}{10^{20}} .
$$
The $3/16$ Ward factor gives $-859794979/3881509632000$ in the weight-five column of $C_4(e^1+e^2)$ from this color word. An independent comparison in dimension twelve recovers all three coefficients of the $1481$ previously integrated propagators. A second reduction in dimension twenty-six starts directly with all $37391$ inputs in reverse order. It uses $4030$ necessary substitutions and $34838$ terminal charts, and recovers all three columns of the displayed period. [Local integration-by-parts reduction](quartic_prism_local_ibp.cpp) · [Separate and full seed replays](operator/mixed_bipartite_local_ibp.py)

#### All bipartite hard insertions on at most two paths

The exact rank census leaves one placement orbit with at most two occupied paths, represented by $(0,0,0,5,5,5)$. Its two edges are disjoint. The $18$ core placements carry $18\cdot6!=12960$ labeled orders. Of the $27$ color orbits, $17$ have a nonzero cofactor and cover $7776$ orders. Sum the ordered-position moments over colors before expanding the star on $(4,6,8)$. The two grounded metric checks hold after summation. Retaining the single color $(U,U,Z,T,T,Z)$ reproduces its original $37391$ weighted reflection components exactly, before evaluating any period.

The full signed numerator has $105222$ monomials and $37391$ nonzero reflection components. The local integration-by-parts reduction terminates on $34928$ wheel orbits. Including the hard permutations and core automorphism factor, the complete two-path contribution is
$$
[t^4]\!\int_{\substack{K_{3,3}\\\text{at most two occupied hard paths}}}\!\Phi_H
=\frac{102163573657672606876837}{40763513236013568000000000}
+\frac{43027855856219197}{6419046553920000000}\zeta(3)
-\frac{6178645771}{606485880000}\zeta(5)
\ \in\ \frac{(663316,663317)}{10^{20}}.
$$
Its contribution to the weight-five coefficient of $C_4(e^1+e^2)$ is $-6178645771/3234591360000$. The $2558$ other nonzero bipartite color orbits have at least three occupied hard paths. [Ordered moments, orbit census and complete period](operator/mixed_bipartite_triples.py)

Every core path supports at most three insertions in a nonzero fourfold graph, including the two soft insertions. Indeed, four marks give at least five segments on a path. Their oriented currents lie in four dimensions, so a nonzero real length variation supported on those segments preserves the total voltage drop along the path. Every cycle traverses the entire path or misses it. The cycle equations and all current-conservation equations still hold at the original currents. This length variation and the common rescaling are independent because the graph has another core path. Their span contains a nonzero variation at fixed total length, forcing the projective current determinant to vanish. Thus two occupied hard paths require three hard marks on each, and both soft marks must lie elsewhere.

Across all five cubic cores, exactly five placement orbits survive with at most two hard paths. In core order their numbers of nonzero color orbits are $17,0,17,47,17$. Including the graph automorphisms, both soft labels and all $14256$ nonzero hard orders gives
$$
[t^4]\!\int_{\text{all cores, at most two hard paths}}\!\Phi_H
=-\frac{1523980123087416073741}{215679964211712000000000}
-\frac{363936831873726773}{19257139661760000000}\zeta(3)
+\frac{6651239339}{231567336000}\zeta(5)
\ \in\ \frac{(1620133,1620134)}{10^{20}}.
$$
The weight-five coefficient contributed to $C_4(e^1+e^2)$ is $6651239339/1235025792000$. [Five-core census and integrated sum with `--all-cores`](operator/mixed_bipartite_triples.py)

Two prism placements with two hard marks on each of three paths each have $38$ nonzero color orbits. In the six-digit path order, their fully integrated bare weights satisfy
$$
10^{20}[t^4]\!\int_{001188}\Phi_H\in(-10059596,-10059595),\qquad
10^{20}[t^4]\!\int_{002266}\Phi_H\in(42351717,42351718).
$$
All $72$ individual numerator banks for $001166$ are saved in `.three-path-integrals/prism-three-paths/001166`. Its integrated weight remains open. On the two completed placement banks and the $020211$ color of $001166$, the streaming rational merge and the independent Python merge reproduce the same $1254373$ nonzero coefficients. — Codex. [Exact periods by color](operator/mixed_prism_direct_periods.json) · [Projective weight merge](operator/mixed_prism_metric_placement.py)

Now require the occupied hard paths to form a perfect matching: three disjoint edges covering the six vertices. The five cores carry $44,0,140,108,36$ nonzero color orbits. All eleven placement orbits together give $42720$ labeled hard orders. At nine unit lengths their summed projective density is
$$
\frac{341900178111254670818124904715920523895306522227}
{14362400045011217512258609674349925376000000000000000000000000000}>0.
$$
Integrating every projective core length reverses the sign on this entire sector:
$$
\begin{aligned}
[t^4]\!\int_{\text{all cores, hard paths a perfect matching}}\!\Phi_H
={}&\frac{2806492431874916206862574308477}{5855271041220988907520000000000}
+\frac{5728260260722488736517423}{4468308181639948800000000}\zeta(3)\\
&-\frac{25663685262992162039}{13171883528643840000}\zeta(5)
\ \in\ \frac{(-1899661956,-1899661955)}{10^{20}}.
\end{aligned}
$$
Run `--perfect-matchings --unit` for the position-integrated current determinant and `--perfect-matchings` for the complete projective periods. [Census and two independent sign evaluations](operator/mixed_bipartite_triples.py) · [Four prism matchings](operator/mixed_prism_direct_periods.json) · [Two bipartite matchings](operator/mixed_bipartite_matching_periods.json)

Every surviving two-path placement has disjoint endpoints. Adding the two-path and perfect-matching values therefore exhausts all graphs in which no core vertex meets two distinct occupied hard paths. The $16$ placements contribute $426$ color orbits and $56976$ labeled orders. Their summed unit density is positive and their full period lies in $(-1898041822,-1898041821)/10^{20}$.

#### Mixed bipartite insertions on three paths with an empty vertex

At denominator $U^{13}$ the convergent positive bipartite powers sum to $52$. Of the $1386282$ graph orbits, $247912$ have all six vertex sums greater than $13$. The row and column identities from the positive angular expansion above give $1511084$ convergent equations touching this interior. The restricted matrix has column rank $247912$ modulo both $999999937$ and $1000000007$, hence over the rationals. The $187418$ required prism cuts and $191519$ contracted four-wheel values give the complete boundary. All $247912$ interior periods were reconstructed on $1,\zeta(3),\zeta(5)$. Substitution of the exact rationals satisfies every one of the $1511084$ integer identities in each column. [Exact period table](operator/mixed_bipartite_positive_periods.gz) · [Checksums](operator/mixed_bipartite_star_bank.json) At denominator $U^4$ the same reduction and all $24$ exact identities reproduce the independent harmonic periods for each of the six interior orbits in the rational, $\zeta(3)$ and $\zeta(5)$ columns. A signed contraction over all $35$ eight-dimensional orbits gives $85487/320-(21965/48)\zeta(3)+(6793/24)\zeta(5)$ by both integrations. [Integer system](operator/mixed_bipartite_positive_ibp.cpp) · [Period reduction](operator/mixed_bipartite_star_bank.py)

The complete dimension-$26$ table gives all $1386282$ positive $K_{3,3}$ period triples and $972294$ contracted four-wheel triples. Holding the scalar contractions fixed leaves $5256$ additional positive powers for direct star reduction. Their $29574$ prism charts occupy $23323$ completed-graph orbits. Permuting the six points of a completed scalar wheel gives $591652$ distinct wheel coordinates. Starting from the archived interior and boundary, $1309$ directly evaluated scalar angular values and the $5256$ stars supply $1411213$ exact single-variable substitutions. Direct angular summation also fixes the $1228$ scalar coordinates unused by the resulting equations. Substitution into all $6983701$ strictly convergent row and column identities in the original coordinates gives zero separately in the rational, $\zeta(3)$ and $\zeta(5)$ columns. Exactly $3983873$ identities involve positive powers alone. [Additional $1731727$ period triples](operator/mixed_bipartite_positive_completion.gz) · [Complete-domain reconstruction](operator/mixed_bipartite_star_bank.py)

For positive integer half-dimension $q$, write $C_e$ for the cycle-matrix column of edge $e$ and $x$ for the four loop momenta on the unit sphere. Feynman parametrization and the spherical Gaussian integral give
$$
F_q(a)=\frac{\Gamma(q)}{2\pi^{4q}}
\int_{S^{8q-1}}\prod_{e\in\mathrm{edges}}|C_e\cdot x|^{-2a_e}\,d\sigma(x).
$$
Hölder's inequality gives log-convexity throughout the convergent domain. Averaging the powers over the edge-transitive automorphism group shows that the uniform powers $4q/9$ minimize the continuous period. Distinct cycle columns are nonparallel, so approaching the zero of one edge momentum while keeping all others nonzero makes the equality condition strict for distinct powers. Rescale the loop coordinates by $(CC^{\mathsf T})^{-1/2}$. The angular Jacobian is $\det(CC^{\mathsf T})^{-q}=81^{-q}$ by the spanning-tree count. On the unit sphere the nine squared edge momenta now sum to one. Their arithmetic–geometric mean inequality gives the dimension-uniform bound
$$
F_q(a)\ge F_q(4q/9,\ldots,4q/9)>
\frac{9^{4q}\Gamma(q)}{81^q\Gamma(4q)}
=\frac{81^q\Gamma(q)}{\Gamma(4q)}.
$$
At $q=13$, rational enclosures of $\zeta(3)$ and $\zeta(5)$ certify positivity for every one of the $2358576$ archived entries. The unique minimizing positive orbit has matrix rows $(5,6,6),(6,5,6),(6,6,6)$. Its period times $10^{60}$ lies strictly between $16076410719557382249365245878$ and $16076410719557382249365245879$. [Full-domain sign and extremum check](operator/mixed_bipartite_star_bank.py)

Eight placement orbits with three occupied paths admit an empty star. The six placements with multiplicities $(3,2,1)$ have $54$ nonzero color orbits each, and the two with multiplicities $(2,2,2)$ have $38$ each. Two further placements occupy a perfect matching and have no empty vertex: $(0,0,0,5,5,8)$ has $28$ nonzero color orbits and $(0,0,5,5,8,8)$ has $8$.

For the perfect matching $(0,0,5,5,8,8)$, the colors $(U,U,T,Z,T,Z)$ and $(U,U,T,Z,Z,T)$ have strictly negative full periods. Their respective contributions to the $\zeta(5)$ coefficient of $C_4(e^1+e^2)$ are $-913276178008155211/40979193200225280000$ and $-8841942850010993/5122399150028160000$. In the first color, summing individual current-minor squares before the graph symmetries reproduces all $266734$ weighted integral orbits of the joint contraction. Applying the untransformed calculation to the former two-triple control reproduces all three period columns. [Matching colors and exact sign bounds](operator/mixed_bipartite_matching_periods.json) · [Direct position moments](operator/mixed_bipartite_matching.py)

An independent contraction reconstructs the angular polynomial at each metric monomial from the rational rank decomposition on the marked paths. At each sum of two hard metric powers it squares the angular polynomials first. For two normalized gap lengths $s,t$, the substitution $t=(1-s)u$ gives the exact moment
$$
\int_{s,t\ge0,\,s+t\le1}s^i t^j\,ds\,dt
=\frac{1}{j+1}\sum_{\substack{0\le k\le j+1\\\mathrm{binomial\ expansion}}}
\frac{(-1)^k}{i+k+1}\binom{j+1}{k}.
$$
The two-path color $002112$ reproduces every one of its $21299$ normalized weights and all three period columns. On the three-path matching $005588$, colors $010212$, $010221$, $011220$ and $012021$ reproduce respectively $293056$, $293062$, $293062$ and $293061$ normalized weights. Every row agrees with the separate metric-square integration, as do both grounded currents and all three period columns. [Duffy contraction and direct archive comparison](operator/mixed_bipartite_factor_contract.py)

The full perfect-matching placement $(0,0,0,5,5,8)$ contains $28$ nonzero color orbits and $15552$ labeled orders. Thirteen individual periods are positive and fifteen are negative. The signed sum is
$$
-\frac{53240001034461628617993228343}{278822430534332805120000000000}
-\frac{126777991585527112397599}{248239343424441600000000}\zeta(3)
+\frac{1419973101995997239}{1829428267867200000}\zeta(5)
\ \in\ \frac{(-198781970,-198781969)}{10^{20}}.
$$
Contracting the combined numerator on its $368755$ nonzero core orbits independently gives all three coefficients. Multiplication by $3/16$ gives $1419973101995997239/9756950761958400000$ in the weight-five column of the contact. [Individual and complete contractions](operator/mixed_bipartite_matching_periods.json)

The second complete matching $(0,0,5,5,8,8)$ contains eight nonzero color orbits and $3456$ labeled orders. Its signed sum is
$$
-\frac{4057923452894466488537085131}{185881620356221870080000000000}
-\frac{17850618195005461188359}{305726349270101760000000}\zeta(3)
+\frac{1168858943446553033}{13171883528643840000}\zeta(5)
\ \in\ \frac{(1410544,1410545)}{10^{20}}.
$$
A separate contraction of the combined numerator on $293053$ core orbits gives the sum of the eight individually integrated colors. To form the numerator, group each signed current minor by its monomial in the nine core lengths, and eliminate the last unit gap on each marked path. For $r$ free gaps with monomial exponents $j_1,\dots,j_r$, the simplex integral is $\prod_i j_i!/(r+\sum_i j_i)!$. Apply the resulting bilinear moment to the two hard factors and accumulate by the sum of their core exponents before multiplying by the soft and voltage polynomial. For the first current square of color $(U,T,U,Z,T,Z)$, $1688$ hard metric monomials have $88032$ distinct pair sums. The $501$ soft and voltage metric monomials use $21$ angular monomials. The complete two-path and matching control colors reproduce all $21299$ and $266734$ weights of the original exterior-factor contraction. [Core-monomial calculation](operator/mixed_bipartite_simplex_square.py) · [Exact moments](operator/simplex_gram_moments.cpp)

The integrated bipartite placements $000558,005588,000113$ span the same rational vector space as $1,\zeta(3),\zeta(5)$. Their three rows of coefficients in that order reduce modulo $1009$ to
$$
\begin{pmatrix}883&266&698\\567&692&716\\360&346&777\end{pmatrix},
\qquad \det\equiv474\pmod{1009}.
$$
Thus rational combinations of these three signed periods give both $1$ and $\zeta(3)$. Their irrational ratio rules out any common rescaling making all three periods rational. This conclusion uses the individual placements of this bare projection. [Exact rank check](operator/mixed_bipartite_matching_periods.py)

All eleven bipartite placements with at most three occupied hard paths now give $453$ nonzero color orbits and $254880$ labeled orders. Their complete bare period is
$$
-\frac{35353476500649492067003536192677}{5855271041220988907520000000000}
-\frac{55180481294366150902934197}{3416941550665843200000000}\zeta(3)
+\frac{1616431563128178150493}{65859417643219200000}\zeta(5)
\ \in\ \frac{(-1446194663,-1446194662)}{10^{20}}.
$$
The weight-five coefficient after multiplying by $3/16$ is $1616431563128178150493/351250227430502400000$. The $2122$ remaining nonzero bipartite colors have four, five or six occupied hard paths. [Exact sectors and placement census](operator/mixed_bipartite_low_occupancy.json)

The four-path placement $(0,0,0,1,3,5)$ has $28$ nonzero colors and $15552$ labeled orders. Integrating their sum on $544016$ core weight orbits gives
$$
\frac{111739010051025299444196343656127}{25094018748089952460800000000000}
+\frac{49414656785827664387305831}{4149143311522809600000000}\zeta(3)
-\frac{203705349125002547}{11254172529600000}\zeta(5)
\ \in\ \frac{(-242467916,-242467915)}{10^{20}}.
$$
The individual color contractions reproduce all three columns. [Higher-occupancy periods](operator/mixed_bipartite_higher_periods.json)

The five-path placement $(0,0,1,2,3,4)$ contains $16$ nonzero colors and $6048$ labeled orders. Its combined numerator and the sum of the color integrals give
$$
\frac{2958546615804524563653960277}{2609300820508462080000000000}
+\frac{44039634819609701715534467}{14522001590329833600000000}\zeta(3)
-\frac{75886861087925533861}{16464854410804800000}\zeta(5)
\ \in\ \frac{(-398996081,-398996080)}{10^{20}}.
$$
The soft compound Gram depends only on the sums of the gaps on each core path. Its degree-two monomials shift the hard metric exponents before applying the common voltage. On the first six-distinct-path color of $(0,1,3,6,7,8)$, applying the voltage separately to each of the $24$ squares returns the same $358146$ normalized weights. [Complete and individual periods](operator/mixed_bipartite_higher_periods.json) · [Joint voltage contraction](operator/mixed_metric_voltage.py)

Move the empty star to $(4,6,8)$ using a core automorphism. At that star project the signed sum of voltage-weighted cofactor squares after integrating the ordered hard positions. Two unequal grounded metrics check the transformed numerator before taking the projective periods. For example, on $(0,0,0,1,1,3)$ the $227274$ nonzero monomials reduce to $108581$ weighted reflection components. The $5240$ components without a direct wheel chart are determined by $5745$ strict integration-by-parts substitutions. The full period on its $31104$ nonzero labeled hard orders is
$$
-\frac{111298014348930480796008751091}{975878506870164817920000000000}
-\frac{59790592893287727727}{196010171557200000000}\zeta(3)
+\frac{10998645992657}{23724300600000}\zeta(5)
\ \in\ \frac{(-165903606,-165903605)}{10^{20}} .
$$
Multiplication by $3/16$ gives $10998645992657/126529603200000$ in the weight-five column of $C_4(e^1+e^2)$. [Projected periods by placement](operator/mixed_bipartite_three_paths.json) · [Contraction, integer reconstruction and independent rank census](operator/mixed_bipartite_three_paths.py)

Summing $(0,0,0,1,1,3)$, $(0,0,0,1,1,5)$ and $(0,0,0,1,1,7)$ covers $162$ color orbits and $93312$ nonzero labeled hard orders. Their combined bare period is
$$
-\frac{26132040591312626178621418507}{121984813358770602240000000000}
-\frac{33018431319359824407107}{57626990437816800000000}\zeta(3)
+\frac{1094147768316101}{1256475458700000}\zeta(5)
\ \in\ \frac{(-425130398,-425130397)}{10^{20}}.
$$
The corresponding weight-five coefficient in $C_4(e^1+e^2)$ is $1094147768316101/6701202446400000$. All three placements pass two grounded metric tests and the labeled-order counts agree with the independent rank census. Grouping common exterior polynomials before expanding each power of the star sum reproduces every coefficient in a $132$-slice test. A second full projection of $000113$ in batches of six colors reproduces the entire $108581$-component weight stream byte for byte.

Across all eight three-path placements with an empty vertex the $400$ nonzero color orbits cover $228096$ labeled orders. Their combined bare period is
$$
-\frac{72988848483165369550563741353}{12524643938440617984000000000}
-\frac{43114280350129380234611953}{2766095541015206400000000}\zeta(3)
+\frac{18573085562524649921}{784040686228800000}\zeta(5)
\ \in\ \frac{(-1249486554,-1249486553)}{10^{20}}.
$$
They contribute $18573085562524649921/4181550326553600000$ to the weight-five coefficient of $C_4(e^1+e^2)$. On $(0,0,1,1,3,3)$ the original current basis, a basis adapted to the marked paths, and a reconstruction from six independent bivector coordinates produce the identical $304544$ reflection components. In the last calculation, if $A$ changes the four cycle rows, a pair of unit soft jumps in the new basis gives the two corresponding columns of $A^{-1}$ in the old basis. Their two-row minors give the six coefficients of this signed projection. Integer row and column lattice reductions of each fifteen-row linear matrix have determinant $\pm1$, whose signs are retained. [Exact placement periods](operator/mixed_bipartite_three_paths.json) · [Signed cycle-coordinate replay](operator/mixed_current_basis.py)

Use the lexicographic core edge order returned by the cubic-core generator. Place three hard marks on path $0$ and three on path $6$ of the prism, or on path $5$ of the bipartite core. Their ordered color word is $(U,U,Z,T,T,Z)$ in both cases. Each colored orbit has weight $8$ before the two soft-label assignments. At fixed core lengths $T_e$ integrate the six hard positions, all soft positions and the voltage average. Write the resulting density as $N/U^{13}$, where $U$ is the degree-four core Kirchhoff polynomial.

Over $\mathbb F_{1009}$ take $T=(2,3,5,7,11,13,17,19,23)$ and then replace $T_0$ by the indicated root of $U$. The last column averages the root numerator after permuting the nine core lengths by all graph automorphisms.

| Core | $N(T)$ | Root $T_0$ | $N$ at the root | Automorphism mean at the root |
|---|---:|---:|---:|---:|
| Prism | $190$ | $824$ | $159$ | $917$ |
| Bipartite | $323$ | $337$ | $36$ | $706$ |

The signed current numerators are evaluated as polynomial determinants directly at $U=0$. On each ordered triple the monomial of powers $i,j,k$ integrates to $1/((i+1)(i+j+2)(i+j+k+3))$, and the two triples also pass the independent three-simplex gap integration. All denominators are units in the chosen field. Since $U$ is a primitive integer polynomial and the final column is nonzero, $U$ does not divide either symmetrized rational numerator.

For the meromorphic bipartite dimension shift the coefficient denominator at $s=-12$ is $-5592001781747700000$. It relates $q=13$ cuts to $q=12$ cuts with shifted powers. Divergent or contracted terms require the face or regulator contributions in that identity. On the prism the existing derivative relation gives $16698\,I_{13}(D_{23}H)=I_{12}(D_{234}D_{123}H)$ when the integration-by-parts boundaries vanish. The two pole tests require an additional reduction before the existing period tables give their integrated values. [Position integration, symmetry average, and both dimension checks](operator/mixed_bare_primitive_pole.py)

One primitive prism contribution does reduce to the existing $q=13$ tetrahedron moments. Move the two triples to paths $3$ and $7$, keeping the same color word. The stars at vertices $0$ and $4$ are source-free. Eliminating them and integrating the three resulting parallel pairs gives the symmetry-weighted chamber amplitude
$$
[t^4]\!\int_{\text{prism},\,(3,3,3,7,7,7),\,(U,U,Z,T,T,Z)}\!\Phi_H
=-\frac{130303752962524429}{449417733427049587200000000}
+\frac{39016635509}{161759973158784000000}\zeta(3).
$$
Only the omitted cycle pair $(1,3)$ survives in the chosen cycle basis. Its signed cofactor has $162$ monomials. The transformed polynomial has $140278$ monomials and gives $741$ tetrahedron moment orbits. Two unequal original metrics agree with the star-transformed densities and Jacobians modulo $1009$. The integral-test tail after $512$ terms of $\zeta(3)$ places the amplitude in $(-197,-196)/10^{17}$. [Exact replay with the prism-period option](operator/mixed_bare_primitive_pole.py)

For the complete prism occupancy $(3,3,3,7,7,7)$ and its core-symmetry images, the sum over all six distinct hard labels is positive:
$$
[t^4]\!\int_{\text{prism},\,(3,3,3,7,7,7)\,\text{orbit}}\!\Phi_H
=-\frac{27376932695374751}{3546541456968510000000000}
+\frac{19673982143}{3063635855280000000}\zeta(3)
\ \in\ \frac{(72,73)}{10^{16}}.
$$
The stabilizer of the two triply occupied paths has order $2$. Its $48$ color orbits include $20$ positive periods, $10$ negative periods and $18$ zero periods, covering all $6\cdot720=4320$ labeled chambers across the six core images. The bounds use the integral-test remainder after $1024$ terms of $\zeta(3)$.

Independently, summing the actual signed voltages on all $720$ hard orders at one fixed path placement gives $432$ nonzero determinants. Their squared primitive polynomials give $36$ groups. The factor $6/12$ for path placements and core automorphisms cancels the two soft-label orders. The combined star-transformed polynomial has $140164$ monomials and contracts to the same two coefficients through $741$ tetrahedron moment orbits. [Every signed color-orbit period](operator/mixed_bare_prism_colors.json) · [Individual integrals and independent labeled contraction](operator/mixed_bare_prism_colors.py)

#### Prism hard placements admitting two disjoint empty stars

All six nonadjacent vertex pairs are conjugate under the prism automorphisms. For the representative pair $0,4$, hard marks can occupy precisely the three-edge chain of paths $3,5,7$. Across all six choices this covers $159$ path multisets in $16$ symmetry orbits. The following five hard populations on $(3,5,7)$ give a nonzero integrand; each population is listed once up to core automorphisms.

| Hard populations | Path multisets | Nonzero colored chambers | Nonzero hard color words at one placement |
|---|---:|---:|---:|
| $(3,2,1)$ | $12$ | $54$ | $54$ |
| $(3,1,2)$ | $12$ | $54$ | $54$ |
| $(3,0,3)$ | $6$ | $30$ | $54$ |
| $(2,3,1)$ | $12$ | $54$ | $54$ |
| $(2,2,2)$ | $6$ | $38$ | $72$ |

Grounded-network evaluations certify the $230$ nonzero chamber integrands. The remaining $1000$ colored chambers each have six identically vanishing signed cofactor polynomials. The $230$ chambers represent $21600$ labeled hard orders. For two further complete occupancies the signed periods, including all soft positions, are
$$
\begin{aligned}
[t^4]\!\int_{(2,3,1)}\!\Phi_H
&=\frac{1093549304340531492169}{203308022264617670400000000000}
-\frac{14329143933271}{3201499468767600000000}\zeta(3),\\
[t^4]\!\int_{(3,2,1)}\!\Phi_H
&=-\frac{4259271731973570113}{919945802102342400000000000}
+\frac{27891974351}{7243211467800000000}\zeta(3).
\end{aligned}
$$
The values lie in $(-13376,-13375)/10^{16}$ and $(-10659,-10658)/10^{16}$, respectively. Both occupancies have $27$ positive, $27$ negative and $36$ zero color-orbit periods. Summing the actual voltages of all $720$ hard orders before integration gives the same rational and $\zeta(3)$ coefficients. Collecting proportional squared determinants leaves $54$ groups in each calculation and uses respectively $1440$ and $2269$ tetrahedron moment orbits. [Full classification and grouped replay](operator/mixed_prism_empty_stars.py) · [Individual signed color periods](operator/mixed_prism_empty_stars.json)

For the denser $(3,1,2)$ occupancy, both star-sum divisors depend only on the six empty-star lengths. Dividing after coefficient extraction in the three marked path lengths lets the integration consume one slice at a time. At color word $(U,Z,T,U,T,Z)$ the $177$ slices contain $1023604$ terms in total and give $2550$ tetrahedron moment orbits. The symmetry-weighted contribution is
$$
-\frac{674063179358209292219}{37641028122134928691200000000000}
+\frac{7747407785891}{518642913940351200000000}\zeta(3).
$$
The original grounded density at unequal lengths agrees with the transformed density including its Jacobian. On a smaller $(2,3,1)$ chamber the $35$ slices sum coefficientwise to the independent $53354$-monomial expansion. [Marked-degree integration check](operator/mixed_prism_empty_stars.py)

The remaining two hard populations give
$$
\begin{aligned}
[t^4]\!\int_{(3,1,2)}\!\Phi_H
&=-\frac{61353993738561592803743}{3513162624732593344512000000000}
+\frac{200912305363961}{13830477705076032000000}\zeta(3),\\
[t^4]\!\int_{(2,2,2)}\!\Phi_H
&=\frac{14121568051760445223361}{9758785068701648179200000000000}
-\frac{61659192605033}{51223991500281600000000}\zeta(3).
\end{aligned}
$$
The positive, negative and zero orbit counts are $(28,26,36)$ and $(25,13,10)$, respectively. Summing the original $720$ hard label orders recovers all fifteen signed resistance coefficients for every one of the $90$ color words before the position integration. [All $(3,1,2)$ periods](operator/mixed_prism_empty_stars_312.json) · [All $(2,2,2)$ periods](operator/mixed_prism_empty_stars_222.json)

Summing all five populations, including their core automorphism factors, gives the complete contribution from prism placements with two disjoint empty stars,
$$
[t^4]\!\int_{\text{two-star prism}}\!\Phi_H
=-\frac{7066378837345229827587769}{307401729664101917644800000000000}
+\frac{26443669806000401}{1383047770507603200000000}\zeta(3)
\ \in\ \frac{(-428900592,-428900591)}{10^{20}}.
$$
The $230$ nonzero colored chambers cover $21600$ labeled hard orders. Their contribution to $C_4(e^1+e^2)$ is $3/16$ of this value and lies in $(-80418861,-80418860)/10^{20}$. The other prism placements and the primitive bipartite core still enter the mixed Ward equation. Run `mixed_prism_empty_stars.py --two-star-sum` to check the sign bounds, orbit counts and all five sums. [Exact replay](operator/mixed_prism_empty_stars.py)

#### Mixed prism with three hard insertions on each of paths zero and six

For $(0,0,0,6,6,6)$ the two empty vertices $2,5$ share path $5$. The star at vertex $2$ replaces paths $(1,3,5)$ by a triangle. Paths $0,5$ then have the same endpoints. For positive Schwinger powers $p_0,p_5<13$ with $p_0+p_5>13$, their parallel contraction gives power $p_0+p_5-13$ and multiplies the four-wheel integral by
$$
\frac{(12-p_0)!\,(12-p_5)!}{(25-p_0-p_5)!}.
$$
The three symmetry-related hard placements have $27$ color orbits. Ten have identically vanishing cofactors. The remaining $17$ account for $1296$ labeled orders. Collecting their signed weights gives $17782$ four-wheel exponent orbits and
$$
\begin{aligned}
[t^4]\!\int_{(0,0,0,6,6,6)}\!\Phi_H
&=-\frac{47798932647136967331974177}{4993530371411662080000000000}
-\frac{3943443518419}{154028952000000}\zeta(3)
+\frac{123892431209}{3184050870000}\zeta(5)\\
&\in\frac{(58244,58245)}{10^{20}} .
\end{aligned}
$$
Each color passes two grounded metric checks. On the separately integrated $(3,3,3,7,7,7)$ color word $002112$ the one-star wheel calculation recovers both tetrahedron coefficients and gives zero in the weight-five channel. For the new occupancy its weight-five coefficient in $C_4(e^1+e^2)$ is $123892431209/16981604640000$. The complete census contains $5573$ nonzero prism chambers admitting one star and no disjoint pair. A further $5929$ have no empty vertex. [Verified color periods](operator/mixed_prism_one_star.json) · [Angular replay and census](operator/mixed_prism_one_star.py)

For an arbitrary occupied prism, number the triangular faces $0,1$ and the rectangular faces $2,3,4$. In the native order $23,02,03,13,12,04,24,34,14$ the nine positive propagator powers sum to $4q$. A projective monomial acquires the angular weight $\prod_{e\in\mathrm{edges}}(a_e-1)!/(q-1)!$. Complete to a sixth point with incident sums $2q$ and require pair powers below $q$ and triple sums below $2q$. Regrounding and four-point reflection give a measured prism, an angular wheel, or a convergent six-point integration-by-parts equation. Equal terminal integrals can be evaluated once for every contributing numerator. At $q=13$ the direct $44469$-term calculation of the $(3,3,3,7,7,7)$ color $002112$ reproduces
$$
-\frac{130303752962524429}{449417733427049587200000000}
+\frac{39016635509}{161759973158784000000}\,\zeta(3)
+0\,\zeta(5).
$$
[Scalar contraction](operator/mixed_prism_scalar_functional.py) · [Independent color control](operator/mixed_bare_prism_colors.json)

There are $67990131$ labeled convergent positive prism power tuples at $q=13$, $5677268$ orbits under the twelve finite prism symmetries and $2591953$ classes under six-point regrounding. Sending point $5$ to infinity and integrating the divergence at a finite face $v$ with a distinct rectangular reference face $r$ gives
$$
0=(a_{v5}-a_{vr})F_q(a)+
\sum_{\substack{j\ne r\\\mathrm{finite\ neighbors\ of}\ v}}
a_{vj}\left[F_q(a+e_{vj}-e_{jr})-F_q(a+e_{vj}-e_{vr})\right],
$$
where $e_{ij}$ raises the indicated finite power by one. A triangular $v$ gives at most five terms and a rectangular $v$ at most seven. Retain equations with all shifted finite powers positive and all completed collision margins strictly positive. The independently integrated $q=6$ propagators satisfy $3030$ triangular-face equations in all three columns. At $q=13$ the triangular-face equations on one representative per angular class determine $1353116$ extra triples from $996529$ initial angular classes. Applying both face identities in every finite chart gives another $226705$ triples from an enlarged set of $2362303$ measured or reconstructed classes. Direct angular integration supplies the remaining inputs. All $2591953$ completed triples satisfy all $41235143$ admissible equations exactly on $1,\zeta(3),\zeta(5)$. [Complete scalar table](operator/mixed_prism_positive_periods.gz) · [Checksum and census](operator/mixed_prism_positive_periods.json) · [Exact verification](operator/mixed_prism_positive_faces.py)

For direct contractions, the $5677268$ ordered finite orbit representatives index the $2591953$ angular values. A common denominator in each period column reduces a sparse diagram numerator to three integer dot products in blocks of $32768$ powers. The five combined numerators at $005566,000566,000556,005556,002466$ reproduce the separate angular reductions in every column. On $005566$, the option `--fast --colors` to `mixed_prism_scalar_functional.py` also reproduces every individual period and sign bound, checks the input hashes, and verifies the color sum at each of the $718126$ powers. [Direct integer contraction](operator/mixed_prism_fast_functional.py)

The three rung powers sum to $m=a_{23}+a_{24}+a_{34}<2q$. For triangle lengths $x$ and rung lengths $y$, the nonzero conductance-Laplacian eigenvalues are $3/x,3/x,2/y,3/x+2/y,3/x+2/y$. The matrix-tree formula gives $U=3x^2(2x+3y)^2$. Rescale the four loop coordinates on the momentum sphere by the inverse square root of the cycle Gram matrix. This supplies $U^{-q}$ and makes the nine squared edge momenta, each multiplied by its edge length, sum to one. Weighted arithmetic–geometric mean with weights $a_e/(4q)$ gives
$$
F_q(a)>
\frac{\Gamma(q)(4q)^{4q}}{\Gamma(4q)\prod_{e\in\mathrm{edges}}a_e^{a_e}}\,
\frac{x^{4q-m}y^m}{U^q}.
$$
The equality condition fails on an open set of the momentum sphere. Maximizing over the ratio gives $y/x=2m/[3(2q-m)]$ and
$$
F_q(a)>
\frac{\Gamma(q)(4q)^{2q}(2m/3)^m(2q-m)^{2q-m}}
{3^q\Gamma(4q)\prod_{e\in\mathrm{edges}}a_e^{a_e}}.
$$

Rational enclosures of $\zeta(3)$ and $\zeta(5)$ verify this bound in all $2591953$ angular classes at $q=13$. The unique minimizing class has powers $(4,6,7,6,7,6,5,5,6)$ in the stated face order and satisfies
$$
34854172321832909061656661637
<10^{60}F_{13}(4,6,7,6,7,6,5,5,6)
<34854172321832909061656661638 .
$$

The prism matching $(0,0,5,5,6,6)$ has all six vertices incident to hard paths. The twenty nonzero color orbits give $1728$ labeled hard orders. Their $11211930$ individually weighted powers occupy $718126$ distinct raw powers. Exactly $716102$ survive the color sum. The separate color integrations and the combined numerator both give
$$
\begin{aligned}
[t^4]\!\int_{(0,0,5,5,6,6)}\!\Phi_H
={}&\frac{287228737714978440594078829}{2859717236249567232000000000}
+\frac{445847562497991969402569}{1659657324609123840000000}\zeta(3)\\
&-\frac{480163611507888199}{1176061029343200000}\zeta(5)
\ \in\ \frac{(11129492,11129493)}{10^{20}}.
\end{aligned}
$$
Thirteen color periods are positive and seven negative. The scalar values entering all twenty pass the preceding sphere bound with rational enclosures of $\zeta(3),\zeta(5)$. The contribution to the weight-five coefficient of $C_4(e^1+e^2)$ is $-480163611507888199/6272325489830400000$. [Twenty color periods and verified weights](operator/mixed_prism_direct_periods.json)

The other three surviving prism perfect matchings are $000566,000556,005556$ in the same path notation. Their $28,30,30$ color orbits each cover $2592$ labeled orders per placement. Together the four placements give $108$ color orbits and $9504$ orders, with
$$
\begin{aligned}
[t^4]\!\int_{\mathrm{prism\ perfect\ matchings}}\!\Phi_H
={}&\frac{540314227753639006815730719287}{780702805496131854336000000000}
+\frac{438877769943963393954359}{237093903515589120000000}\zeta(3)\\
&-\frac{945315064816578959}{336017436955200000}\zeta(5)
\ \in\ \frac{(-580484726,-580484725)}{10^{20}}.
\end{aligned}
$$
The weight-five coefficient after the $3/16$ Ward factor is $-945315064816578959/1792092997094400000$. The separate prism placement $002466$ contributes $41928787200489559/2949947320320000$ to the same coefficient from its $23$ colors and $1872$ hard orders.

#### Mixed prism with a single hard insertion on the opposite triangular face

Take three and two insertions on distinct paths of one triangular face and one insertion on the opposite face. The $2\cdot 3!\cdot 3=36$ placements fall into three prism orbits represented by $(0,0,0,1,1,j)$. Each representative has a unique empty vertex, $54$ nonzero color chambers and $5184$ labeled hard orders including its orbit. The other $36$ color words per representative have all six cofactor polynomials identically zero.

| Singleton path $j$ | Empty vertex | Four-wheel powers with nonzero weight | Strict bounds for $10^{20}[t^4]\int\Phi_H$ |
|---|---:|---:|---:|
| $6$ | $5$ | $34404$ | $(-182758672,-182758671)$ |
| $7$ | $4$ | $44709$ | $(-129626601,-129626600)$ |
| $8$ | $3$ | $55078$ | $(-98629885,-98629884)$ |

The sum on these $36$ placements is
$$
\begin{aligned}
[t^4]\!\int_{\text{opposite-face singleton}}\!\Phi_H
={}&-\frac{92034143393698363730707}{44358113948643855360000000000}
-\frac{5278772839967513}{960449840630280000000}\zeta(3)
+\frac{1012998047}{120993933060000}\zeta(5)\\
&\in \frac{(-411015157,-411015156)}{10^{20}} .
\end{aligned}
$$
Its contribution to the $\zeta(5)$ coefficient of $C_4(e^1+e^2)$ is $1012998047/645300976320000$. Each of the $162$ nonzero color chambers passes two grounded metric checks. For $002112$, moving the empty vertex to vertex zero reproduces all $14734$, $17975$ and $26687$ individual wheel weights, respectively. [Exact color periods](operator/mixed_prism_opposite_face.json) · [Geometry, angular reconstruction and relabeling checks](operator/mixed_prism_opposite_face.py)

#### Three hard insertions on a prism rung

Put three insertions on a rung, two on an incident triangular-face edge and one on any edge of the opposite face. The $3\cdot 2\cdot 2\cdot 3=36$ placements reduce to $(0,0,2,2,2,j)$ with $j=6,7,8$. Each has $54$ nonzero color orbits and represents $5184$ labeled hard orders.

The six hard positions are integrated jointly over colors. For a monomial in the gap fractions on a path with $k$ marks and gap exponents $m_0,\ldots,m_k$, the ordered-position moment is $\prod_{i=0}^{k}m_i!/(k+\sum_{\substack{0\le i\le k\\\mathrm{gaps}}}m_i)!$. At fixed total degree on each marked path, the denominator is common to every color. Two copies of a signed cofactor and the voltage multiplied by its soft-placement weight supply three exterior polynomials. Dividing each by its signed integer content permits collecting every coefficient of the same three factors and star-sum power before multiplication. The sum of the three transformed star lengths contains only unmarked variables, so it divides each resulting marked-degree coefficient. Exchanging the two $U$ sources or the two $T$ sources preserves the squared cofactor and its averaged voltage. Choosing the exchange with the sparsest tree flow reduces the peak size of the intermediate determinants.

| Singleton path $j$ | Nonzero moment contractions | Collected products | Four-wheel powers | Strict bounds for $10^{20}[t^4]\int\Phi_H$ |
|---|---:|---:|---:|---:|
| $6$ | $46848$ | $3876$ | $26699$ | $(-126233898,-126233897)$ |
| $7$ | $46828$ | $4152$ | $28922$ | $(-133756976,-133756975)$ |
| $8$ | $355100$ | $78138$ | $38933$ | $(-150756090,-150756089)$ |

The $j=8$ product count adds the nine ordered pairs of colors on the doubly marked path. The wheel powers are combined across all nine. Two grounded metric tests check each collected polynomial. On five nonzero $j=6$ colors all $13332$ wheel weights agree with separate projections in the original coordinates. The $j=7$ integral reproduces both coefficients of the independent two-star $(2,3,1)$ result and gives zero for $\zeta(5)$. The sum over the three positions is
$$
\begin{aligned}
[t^4]\!\int_{\text{triple on rung}}\!\Phi_H
={}&\frac{118235253217539840096529691}{1191548848437319680000000000}
+\frac{2549015695929591935269}{9604498406302800000000}\zeta(3)
-\frac{24401861098597}{60496966530000}\zeta(5)\\
&\in\frac{(-410746963,-410746962)}{10^{20}} .
\end{aligned}
$$
The contribution to the weight-five coefficient of $C_4(e^1+e^2)$ is $-24401861098597/322650488160000$. The two new placements $j=6,8$ account for $108$ nonzero chambers. Including the $230$ two-star chambers and the preceding $17$ and $162$ one-star chambers, the $517$ disjoint prism contributions give $-8818959362641/129060195264000$ for the normalized weight-five coefficient. There remain $5286$ one-star chambers and $5929$ prism chambers without an empty vertex. [Three exact period columns, with overlap removed in the new subtotal](operator/mixed_prism_rung_sum.json) · [Joint color contraction and independent checks](operator/mixed_prism_rung_sum.py)

#### A determinant for all placements of parallel soft legs

On a fixed bare graph of loop number $g$ in complex dimension $d$, remove $r$ labeled bivalent soft marks of common direction $v$. The remaining $1+dg-r$ segment lengths are $\ell_i$. Take $H$ with columns $(1,C_i\otimes j_i)$ for the hard currents at zero soft momentum, and $B$ the map sending a cycle column $c$ to $(0,c\otimes v)$ in the same row order. After subtracting the two adjacent segment columns, a soft insertion of momentum $t v$ supplies $tB C_i$ up to an orientation sign. Repeated insertions on a single segment give proportional columns. At leading order all remaining currents are evaluated at zero. For $r>g$ this order vanishes. For $r\le g$, choose independent rows of $H$ forming a square block $H_0$ and write $H_1$ for the remaining rows. Split $B$ in the same rows as $B_0,B_1$. If $H$ lacks full column rank every leading minor vanishes. Otherwise put $T=B_1-H_1H_0^{-1}B_0$. The complete position integral of the leading **squared current cofactor** is
$$
\frac{r!\,(\det H_0)^2}{(\sum_{i\in\mathrm{segments}}\ell_i)^2}
\,\det(T M T^{\mathsf T}).
$$
Indeed, appending the columns $BC_i$ for a selected set of segments gives determinant $\pm\det H_0$ times the corresponding minor of $TC$. The $r!$ assignments of labels each carry the product of those segment lengths. Cauchy–Binet sums these squared minors with that measure, since $M=C\,\mathrm{diag}(\ell)C^{\mathsf T}$. At $r=2$ this computes the mixed Ward numerator by a two-row Gram determinant; the skew matrix above is decomposable and its Pfaffian vanishes. The direct calculation agrees on all $105$ pairs and on all $364$ triples for a three-soft fixture.

At $r=g$ a stronger equality holds at **finite** strengths $a_1,\ldots,a_g$ in the last coordinate direction. All $g$ rows in that coordinate are already used by the jump columns. The remaining minor is the current cofactor of the hard network projected on the first $d-1$ coordinates. It is independent of the individual positions and of all hard currents in the last coordinate. The restored squared cofactor summed over positions is exactly
$$
g!\left(\prod_{k\in\mathrm{soft\ labels}}a_k^2\right)\det M
\;\times\;(\text{projected hard current cofactor})^2 .
$$
The Kirchhoff denominator $(\det M)^d$ consequently reduces to $(\det M)^{d-1}$. A voltage difference at retained vertices tested in a projected direction is also unchanged. This proves dimension reduction for the reduced bare density after summing the labeled children of **each hard graph**. The four-loop check restores two placements of the four hard and four collinear waves on a seventeen-edge network with soft strengths $(2,-1,3,5)$ and nonzero longitudinal recoil. It compares the literal current determinant and every retained transverse potential. [Exact determinant and finite-strength replay](operator/fourfold_mixed_soft_ward.py)

The same contraction evaluates a Ward voltage **at a specified soft mark**. On segment $i$, let $h_i$ be the mean of the endpoint voltages from the hard currents. At zero soft strength the voltage is affine along that segment. When $\det(TMT^{\mathsf T})\ne0$, multiplying the previous squared-cofactor integral by
$$
\frac1r\,\operatorname{tr}\!\left[
(TMT^{\mathsf T})^{-1}TC\,\operatorname{diag}(\ell_i h_i)\,C^{\mathsf T}T^{\mathsf T}
\right]
$$
gives its integral against that voltage. Differentiate the Cauchy–Binet polynomial after replacing each $\ell_i$ in its measure by $\ell_i(1+s h_i)$, keeping the hard currents fixed, and divide by $r$ for one specified label. If the Gram determinant vanishes, each placement has zero squared cofactor. For $r=g$ the edge probabilities reduce to $\ell_i\,\partial_{\ell_i}\log\det M/g$, the graph's canonical probability measure; the transverse voltage formula then holds at finite collinear strengths. The mixed four-loop replay compares all $105$ pairs and the soft-site voltages on the restored network. [Voltage-weighted calculation](operator/fourfold_mixed_soft_ward.py)

### A scalar background measuring planar moving descendants

In the flat fivefold coordinates take one scalar wave of momentum $e^1$ at level $k_1$, a vector wave of momentum $e^2$ with polarization $e_1$ at level $k_2$, and the primary scalar background $a e^{z^3}+b e^{z^4}+c e^{z^5}$. Extract the coefficient linear in the first two waves in $Q_g\phi(0)$. The absorbed primary root has negative unit pairing with $e_1$. In the minimal theory the scalar level is at most four and the vector level at most three. Take $k_1+k_2=g$. In the ordinary $A_k$ field coordinates, on the plane $a+b+c=0$ this measurement obeys
$$
(2\pi i)^{5g} Q_g\phi(0)
=(-1)^g(abc)^g
\int_{\mathcal M_{g,3}}\frac{\Theta_g^{2g}}{(2g)!}\,\psi_1^{k_1}\psi_2^{k_2}.
$$
Here $\Theta_g$ is the Betti form of the planar unit triangle $(1,0),(0,1),(-1,-1)$ and $\psi_i$ is the cotangent-line class at mark $i$. The displayed coefficient uses the trace and heat conversion of the flat action above.

The fivefold dimension fixes the number of extra primary waves to $3g$. Expand each ordinary cotangent power into pullbacks after forgetting these waves and terms where primary marks fuse into a retained mark. Fix the unabsorbed wave counts on each axis and the number $s_i$ of fusions at each retained mark. A determinant-one shear, adding multiples of the first two momentum coordinates to the last three, removes the transverse momenta at the first two retained sites independently. Charge conservation then fixes the transverse momentum at the root from the unabsorbed waves alone. Thus the pulled-back trace is independent of which axes supplied the $s_i$ fusions. This use of a shear concerns the trace; the original root contraction has the constant value $-1$. The identical-wave Taylor factorials and the labeled fusion counts sum to the factor $(a+b+c)^{s_i}/s_i!$ at each site. Only the term with no fusions survives. An axis containing more than $g$ free marks is deficient in relative rank because its Jacobian has dimension $g$. Each axis therefore carries exactly $g$ marks. The three ordered Abel maps have degree $g!$ apiece, canceled by the three identical-wave factorials.

On the fixed real background $e^{z^3}+e^{z^4}-2e^{z^5}$, the heat-normalized coefficients at genus three for $(k_1,k_2)=(3,0),(2,1)$ are
$$
8\frac{7}{10368000}=\frac{7}{1296000},\qquad
8\frac{17}{7257600}=\frac{17}{907200}.
$$
At genus one the complete fivefold wave polynomial of the action is $abc/24$. At genus two the unit-triangle values at levels $(2,0),(1,1)$ are $1/3456,11/17280$. The fivefold Ward evaluation agrees on the entire plane $a+b+c=0$. For arbitrary $a,b,c$, subtracting twice the $(2,0)$ response from the $(1,1)$ response gives exactly $(abc)^2/17280$; only the Fourier multiplicities $(2,2,2)$ survive this subtraction. [Wave polynomial and graph replay](simplex_scalar_background.py) · [Genus-three triangle integrations](../twisted-holography-N4/README.md#three-point-genus-three-vertex-in-two-target-dimensions)

### Lower-dimensional vertices from prescribed fivefold descendant fields

Take conserved hard momenta projected to $d$ independent coordinate directions on the fivefold and a primary scalar background $B$ with one unit exponential on each of the other $5-d$ axes. The amplitudes on those axes may be arbitrary. Keep a primary root. For a hard scalar wave $h_i$ of level $k_i$, supply its correlated components in the ordinary action fields,
$$
A_{k_i-j}\ \supset\ \frac{(-B)^j}{j!}h_i,\qquad 0\le j\le k_i .
$$
Polarize once in each hard wave, including any primary vector used for the elementary scalar test. Require that the hard levels sum to $3g-3+n-dg$ with $n$ retained marks. The resulting fivefold action coefficient is the $d$-target vertex on the hard projected momenta times the $g$th power of the product of transverse amplitudes. Arbitrary hard longitudinal shifts are allowed if the root balances the waves. For the rooted operation include the original root pairing at the surviving Fourier multiplicities and the heat conversion. This is a prescription for the test fields in one fixed action.

Indeed after lowering the levels by a total of $s$ the fivefold dimension leaves $(5-d)g-s$ free primary waves. The $s$ factors of $B$ attached to the retained fields supply the other $s$ waves. At fixed labeled multiplicities each term has precisely the sign and factorial of fusing its attached marks into its prescribed site. Summing over them gives the independent cotangent pullbacks by the rational-tail inversion. On a fixed curve with its hard marks the remaining map to each transverse Jacobian depends only on its own free mark positions and on a translation. At each original multiplicity the whole signed sum equals the fiber integral of this map over the labeled forgotten marks, divided by their identical-wave factorials. Too many positions on one axis kill the fiber determinant. The unique surviving multiplicity is $g$ on every axis, with vertical degree $(g!)^{5-d}$. This cancels all the factorials. This derivation applies also to nonsimplex hard data.

For the previous two planar genus-three probes, take $B=e^{z^3}+e^{z^4}+e^{z^5}$, keep their hard charges and primary vector, and use the fields just specified at their four scalar sites. Their two ordinary background responses, in $(2\pi i)^{15}Q_3\phi(0)$, are respectively
$$
\frac{4093949}{2514758400}+\frac{4}{121275}\zeta(3),\qquad
\frac{11147}{13608000}.
$$
Multiplication by $(3!)^3$ gives the labeled measurements. On a nonsimplex elliptic four-mark test, direct fivefold evaluation gives $5/12$ times the product of the three transverse amplitudes. On two independent genus-two threefold descendant tests the full fivefold graph evaluator gives $569/77760$ and $739/77760$ times the square of the product of the two transverse amplitudes. The replay checks every possible transverse multiplicity, including those which vanish. [Exact profile expansion and replay](descendant_background.py) · [Planar source momenta](../twisted-holography-N4/README.md#a-genus-three-fivefold-scalar-measurement-from-the-planar-descendant)

#### The transverse Hessian on a general holomorphic background

Allow any holomorphic transverse scalar background $B$, with the same correlated hard waves. In the supplied affine volume coordinates put $r=5-d$, and take the $r$ by $r$ Hessian of $B$ in the transverse coordinates. The transverse factor in the ordinary action coefficient is exactly
$$
\bigl(\det \partial_\perp^2 B(0)\bigr)^g .
$$
For plane-wave hard inputs, the rooted coefficient in $(2\pi i)^{5g}Q_g\phi(0)$ on a primary vector with polarization $v$ multiplies the lower-target vertex by
$$
(-1)^{g+1}\left[
-\left(\sum_{i\,;\,\mathrm{hard\ inputs\ except\ root}}p_i(v)\right)
 (\det\partial_\perp^2B)^g
-\partial_v(\det\partial_\perp^2B)^g\right]_{z=0}.
$$
Indeed the absorbed momentum is minus the sum of all the fixed hard and background momenta. The sums of transverse momenta are generated by differentiating the background factor. At a singular transverse Hessian both terms vanish for $g\ge2$.

For the vertical calculation, first let $B$ be a sum of holomorphic exponential waves with arbitrary transverse covectors. Temporarily label the $rg$ background occurrences before dividing by $(rg)!$. At a fixed marked curve choose a period-orthonormal basis of its $g$ holomorphic differentials. The top vertical Jacobian density is a squared determinant. Each column at a forgotten mark is its transverse covector times the evaluations of these differentials at that mark. In the two determinant expansions the fiber integral pairs two evaluations of differentials to the identity. A surviving assignment therefore splits the labeled columns into $g$ ordered groups of $r$ columns, one group per differential. Its coefficient is the product of the squared $r$ by $r$ determinants of the covectors in the groups. This argument also covers real, nonintegral, and by polynomial continuation complex momenta. The square is the ordinary algebraic square in the continued formula.

The count of ordered groups is $(rg)!/(r!)^g$. For each group, summing the ordered wave choices with their amplitudes gives $r!\det\partial_\perp^2B(0)$ by Cauchy–Binet. The fivefold Taylor factorial cancels the remaining $(rg)!$. Both sides are polynomials of fixed finite jet order at these arities. Arbitrary finite transverse jets are spanned by exponential jets by polynomial interpolation, proving the general holomorphic formula. For the ordinary simplex waves in the preceding section, $B(0)=0$ cancels every fusion without adding any lower-level components to the hard waves, as long as the vector polarization is tangent to the hard plane. This includes the purely quadratic profile
$$
 B(z)=\tfrac12\bigl((z^3)^2+(z^4)^2-2(z^5)^2\bigr)
$$
with the same two genus-three values $7/1296000$ and $17/907200$. On the same plane, using an additional transverse wave of momentum $e^3+e^4+e^5$, the amplitudes $(1,-1,2,-2)$ and $(1,2,-4,1)$ give Hessian determinants $0$ and $-18$. The ordinary action on these backgrounds is respectively $0,-3/4$ at genus one, and $0,33/160$ at genus two on two level-one waves. [Direct fivefold graphs on all four transverse waves](simplex_scalar_background.py) · [Independent correlated-field and root-derivative checks](descendant_background.py)

#### Planar Hamiltonian waves and their cyclic composition

For planar scalar and vector momenta $p,q$, the three retained momenta are $p,q,-p-q$. The linear map with rows $p,q$ sends the planar unit triangle to these charges (up to transpose according to the column convention). The top Betti product multiplies by $[pq]^{2g}$, where $[pq]=p_1q_2-p_2q_1$ in the unit area coordinates. This follows by applying the same two-dimensional map to each pair of period directions. Rank-deficient pairs follow by polynomiality. For the Hamiltonian vector wave with polarization $(q_2,-q_1)$ the scalar momentum evaluates to $[pq]$ and the vector is divergence free. At $B(0)=0$, for the scalar at level $g$ and the primary vector the heat-normalized $Q_g\phi(0)$ symbol at one, two and three loops is respectively
$$
-\frac{[pq]^3}{24}\det\partial_\perp^2 B(0),\qquad
\frac{[pq]^5}{3456}\bigl(\det\partial_\perp^2 B(0)\bigr)^2,\qquad
-\frac{7[pq]^7}{10368000}\bigl(\det\partial_\perp^2 B(0)\bigr)^3 .
$$
At levels $(2,1)$ the seventh-power coefficient is $-17/7257600$. The Fourier symbols determine the corresponding bilinear differential operators on general tangent jets by polarization. The tree with two primaries gives $[pq]$.

For a level-four tangent scalar and primary Hamiltonian vector take $B=((z^3)^2+(z^4)^2+(z^5)^2)/2$. In the coefficient linear in both waves, the genus-$h$ action has $2h+4$ other primary scalar inputs. At the origin each copy of $B$ uses exactly two normal derivatives. The top Betti product supplies exactly $6h$ normal derivatives and the root pairing is tangent. Consequently $2(2h+4)=6h$ is necessary. Only $h=4$ contributes, and the full three-mark computation gives
$$
(2\pi i)^{20}Q_4\phi(0)=\frac{731}{1024192512000}[pq]^9 .
$$
For the $(3,1)$ and $(2,2)$ scalar and vector levels the coefficients of the same ninth power are $1493/438939648000$ and $157/29262643200$. [All genus-four stable graphs and independent weights](../twisted-holography-N4/README.md#three-point-genus-four-vertex-in-two-target-dimensions)

On the same unit transverse quadratic background, take instead a level-four scalar and a level-one Hamiltonian vector. Both levels belong to minimal fivefold BCOV theory. The ordinary matrix element is
$$
(2\pi i)^{25}Q_5\phi(0)=-\frac{14951}{6083703521280000}[pq]^{11},\qquad Q_h\phi(0)=0\quad(h\ne5)
$$
on the part linear in both waves. At genus $h$ the two waves and the absorbed primary root leave $2h+5$ primary quadratics. The required number of normal derivatives is $2(2h+5)$, while the Betti form supplies $6h$. At $h=5$ the period gives $[pq]^{10}$ and the polarized root gives $-[pq]$ in heat units. For scalar and Hamiltonian vector levels $(3,2)$ or $(2,3)$, exchanging the scalar and vector cotangent exponents preserves the triangle period. Both measurements give
$$
(2\pi i)^{25}Q_5\phi(0)=-\frac{463}{86910050304000}[pq]^{11},\qquad Q_h\phi(0)=0\quad(h\ne5).
$$
The $(5,0,0)$ period would require a level-five hard field in this direct bilinear projection, beyond the minimal scalar and vector levels. [Complete three-mark integral](../twisted-holography-N4/README.md#three-point-genus-four-vertex-in-two-target-dimensions)

**A single five-loop matrix element from the three positive cotangent levels.** Use the same transverse quadratic $B=((z^3)^2+(z^4)^2+(z^5)^2)/2$, and take the component linear in each of three hard fields. At level two the hard scalar fields are $e^{z^1}$ and $-z^1e^{-z^1-z^2}$; the hard vector is $e^{z^2}e_1$. For each hard field use its level-one and primary components obtained by multiplication by $-B$ and $B^2/2$, respectively. The primary scalar background is $B$. All vector components are divergence-free. In the field coordinates of the flat worldsheet quantization,
$$
(2\pi i)^{25}Q_5\phi(0)\big|_{\rm trilinear}
=\frac{51787}{965667225600000},
\qquad Q_h\phi(0)\big|_{\rm trilinear}=0\quad(h\ne5).
$$
To see the genus selection, let the total amount of level lowering be $s$. The four retained sites include the primary root. At genus $h$ there are $2h+3+6-s-4=2h+5-s$ additional primary quadratics. The lowered components provide $s$ more. Every quadratic absorbs exactly two normal derivatives at the origin; the top Betti product has exactly $6h$ normal derivatives and the root contraction is tangent. Thus $2(2h+5)=6h$.

Replace the second scalar by $e^{(-1-\epsilon)z^1-z^2}$ and take its first derivative at zero. The planar absorbed momentum is $\epsilon e^1$. The heat conversion and vector contraction give $(-1)^{5+1}\epsilon$, and the transverse Hessian determinant is one. For the level-two vector the same cancellation of rational tails holds as for a level-two scalar: the vertex is the scalar descendant symbol times its one vector polarization; the divergence at the root is tangent and annihilates $B$. At zero root momentum the string equation on levels $(2,2,2,0)$ lowers each hard level once. The three permutations of levels $(2,2,1)$ on the unit momentum triangle are equivalent by planar changes of determinant $+1$ or $-1$; the vertex transforms by the tenth power. The labeled sum is therefore
$$
3\,\frac{51787}{2897001676800000}=\frac{51787}{965667225600000}.
$$
The same soft derivative on the second scalar also measures the levels $(4,1,1)$ on the first scalar, second scalar, and vector, with each lower component $(-B)^j/j!$ times its hard wave. Their total is again six, so the same counting isolates genus five. The three string lowerings give
$$
(2\pi i)^{25}Q_5\phi(0)\big|_{\rm trilinear}
=\frac{5783}{506975293440000}+2\frac{14951}{6083703521280000}
=\frac{49649}{3041851760640000}.
$$

**All minimal levels for this trilinear soft-root measurement.** Give the first scalar momentum $p$ and the vector momentum $q$ with polarization $(q_2,-q_1)$. Give the second scalar momentum $-p-q-\epsilon\sigma$ and differentiate at $\epsilon=0$, where $\sigma$ is any planar covector. Use the transverse quadratic and lower components above. For hard levels $k_1,k_2,k_3$, the number of quadratic backgrounds after counting those supplied by the lower components is $2h-1+k_1+k_2+k_3$ at genus $h$. Twice that number must be $6h$, so only $g=k_1+k_2+k_3-1$ survives. Write $V_{g;k_1,k_2,k_3}$ for the project's three-mark unit-triangle integral with the indicated cotangent powers. The soft root contributes $[\sigma q]$ and the determinant of the map from the unit triangle contributes $[pq]^{2g}$. The string equation gives the whole symbol
$$
 (2\pi i)^{5g}Q_g\phi(0)\big|_{\rm trilinear}
 =(-1)^{g+1}[\sigma q][pq]^{2g}
 \sum_{\substack{1\le j\le3\\\mathrm{hard\ marks\ with}\ k_j>0}}
 V_{g;k_1-\delta_{1j},k_2-\delta_{2j},k_3-\delta_{3j}} .
$$
There is one term per positive mark, including each repeated level separately. Every permutation of the three charges on the unit triangle is induced by a map of determinant $\pm1$. Thus the periods in this sum depend only on the unordered levels. At genus two, for hard levels $(3,0,0),(2,1,0),(1,1,1)$ the three signed unit symbols are $-1/3456,-1/1080,-11/5760$. Both the fivefold Ward graph and interpolation from five nonzero soft momenta reproduce all three. The interpolation also checks $[pq]=3$ and $[\sigma q]=-3$. [Fivefold and nonzero-momentum replay with the soft-root option](descendant_background.py)

There are exactly four linear relations among minimal level assignments that follow from these string lowerings alone, one in each of genera $1,2,3,5$. The following coefficients act on measured symbols in the indicated order:

| Genus | Unordered hard levels | Coefficients with zero sum of weighted symbols |
|---|---|---|
| 1 | $200,110$ | $-2,1$ |
| 2 | $300,210,111$ | $3,-3,1$ |
| 3 | $400,310,220,211$ | $2,-2,1,0$ |
| 5 | $420,411,330,321,222$ | $12,-6,-9,6,-2$ |

At genus five the full signed coefficients of $[\sigma q][pq]^{10}$, in the last row's order and with the factor $(2\pi i)^{25}$ on the differential, are
$$
 \left(
 \frac{15787}{2027901173760000},\
 \frac{49649}{3041851760640000},\
 \frac{463}{43455025152000},\
 \frac{191417}{5530639564800000},\
 \frac{51787}{965667225600000}
 \right).
$$
The stable-graph value $V_{5;3,2,0}=463/86910050304000$ enters the first, third, and fourth entries with multiplicities $1,2,1$. The five fractions satisfy the displayed integer relation exactly.

Here three successive digits are the three hard levels. To verify exhaustiveness, list all partitions of $g+1$ into three nonnegative levels at most four. Up to total eleven, at least one level is at most three and can be assigned to the vector. Larger totals exceed the minimal-field bounds. Lower each positive mark once and collect equal unordered triples in the target. For $g=1,\ldots,10$ the numbers of input triples are $(2,3,4,4,5,4,4,3,2,1)$ and the exact ranks of the ten integer matrices are $(1,2,3,4,4,4,4,3,2,1)$. These ranks treat the three-mark periods as independent data; arithmetic relations among evaluated periods are separate. [Integer matrix replay](descendant_background.py) · [Three-mark graph integrals](../twisted-holography-N4/README.md#three-point-genus-four-vertex-in-two-target-dimensions)

These particular scalar and Hamiltonian projections fail the Jacobi identity when composed as binary operations on planar functions. Take a unit transverse Hessian. For three momenta $(1,0),(0,1),(-2,-3)$ the inner and outer area brackets in the three cyclic terms are $(1,-1),(2,2),(3,-1)$. On composing the displayed symbols, the tree and one-loop cyclic sums vanish. At two loops the two tree–quintic placements and the cubic–cubic placement give respectively
$$
\frac{-120}{3456}+\frac{36}{24^2}=\frac1{36}.
$$
This calculation tests closure after projecting to the specified scalar descendant and primary Hamiltonian input at each order. The full BV differential also acts on the descendant and transverse fields. [Polynomial identity at arbitrary planar brackets and a determinant-five genus-three check](simplex_scalar_background.py)

#### Linearized bulk fields on arbitrary quadratic backgrounds

Let $B$ be a primary scalar of degree at most two in the five affine volume coordinates. Write $H$ for its constant Hessian, $\delta=\det H$, and $L=(\operatorname{adj}H)^{ij}\partial_i\partial_j$. On finite descendant polynomials, the **local Hamiltonian linearization** of the chosen quantum vertices is isomorphic to the tree complex by an invertible triangular differential operator. It preserves polyvector rank and compact support, only lowers descendant levels, and fixes every primary. This gives an isomorphism of filtered linear state complexes, including the minimal fields. Preservation of the BV pairing or the interacting quantum theory would require additional identities.

In the positive cotangent coefficients $\zeta^rY$, put $d=\partial_\Omega$, let $M$ multiply by $B$ and let $K=[d,M]$ contract with its gradient. Transposing the quadratic Euler formula gives the primary Hessian $T_{r+1}$ on a level-$r$ incoming variation, where $T_s=[x^s]T(x)$ and
$$
 T(x)=R(x,L)e^{xM},\qquad R(x,L)=e^{\lambda^2x^3\delta/24}
 \sum_{\substack{k\ge0\\\mathrm{even\ derivative\ pairs}}}
 \frac{\lambda^{2k}x^{2k}L^k}{4^k(2k+1)!}.
$$
The order matters: the even derivatives act on $B^mY$ after transposition. The holomorphic state differential is $D(\zeta^rY)=-\zeta^{r+1}dY+dT_{r+1}Y$. In particular $T_1=M$ at every loop order. The matrix identity $H\operatorname{adj}H=\delta\mathbf1$ gives $[L,K]=2\delta d$; also $[d,L]=0$ and $dK+Kd=K^2=0$. Thus $D^2(\zeta^rY)=K dT_{r+1}Y-dT_{r+2}dY=0$. Indeed $K$ passes through $L^k$ on a divergence since the commutator contains a second $d$. Then $K d(B^mY)=B^m KdY$, giving exactly the term from $B^{m+1}/(m+1)!$ in $T_{r+2}$. Its constant-power term gives $dL^kd=0$.

Let $[\,]_+$ drop strictly negative powers of $\zeta$. The map $Y(\zeta)\mapsto[T(\zeta^{-1})Y(\zeta)]_+$ intertwines $D$ with
$$
 -\zeta d+K+\big[2\delta(\partial_L\log R(\zeta^{-1},L))\,d\,\cdot\,\big]_+ .
$$
Commuting $-\zeta d$ past truncation leaves minus $d$ of the omitted $\zeta^{-1}$ coefficient, canceled by the primary $dT_{r+1}Y$. The remaining identity follows from $[d,T]=xRKe^{xM}$ and $RKR^{-1}=K+2\delta(\partial_L\log R)d$, where $\partial_L$ differentiates the formal series in $L$. At tree level the transformed differential is $-\zeta d+K$.

To pass from that tree operator to the transformed operator at all genera, send $\zeta^rY$ to $u_r(\zeta,L)Y$. Set $u_0=1$ and
$$
 u_{r+1}=\zeta u_r+2\delta\partial_L u_r-
 \big[2\delta(\partial_L\log R)u_r\big]_+ .
$$
Commuting $u_r$ through $K$ gives this recurrence on each input. Its leading term is $\zeta^r$, with only lower-level corrections. Both transformations have finite inverses on each fixed tower and commute with $\bar\partial$. Their composition with the inverse transformations gives the stated isomorphism. If $\operatorname{rank}H\le3$, then $\delta=L=0$ and the original quantum linearization already equals the tree one. [Quadratic Euler Hessian](operator/README.md#all-quadratic-backgrounds-and-the-unweighted-d3-linking-residue) · [Direct square and both intertwiners on two Hessian ranks](quadratic_linearized_complex.py)

**Test on the five-loop states.** Use the rank-three quadratic and planar covectors $p=(1,0)$ and $q=(0,1)$ above. In the tree coordinates $[e^{B/\zeta}Y]_+$ the dressed level-$(4,1)$ inputs are $\zeta^4e^{z^1}$ and $\zeta e^{z^2}\partial_{z^1}$. For $D'=-\zeta d+K$ both have primitives in the minimal complex:
$$
 D'(-\zeta^3 e^{z^1}\partial_{z^1})=\zeta^4e^{z^1},\qquad
 D'(e^{z^2}\partial_{z^1}\wedge\partial_{z^2})=\zeta e^{z^2}\partial_{z^1}.
$$
Multiplication by $e^{-B/\zeta}$ and truncation returns precisely the correlated components of the planar triangle. The primary scalar at the origin annihilates every holomorphic linearized boundary: $B=dB=0$ there and each primary output is the divergence of $B^{r+1}Y/(r+1)!$. Yet the five-loop bilinear Hamiltonian on these exact inputs has primary coordinate $-14951/6083703521280000$ in the heat units above. At one loop the level-one scalar is exact and the primary vector is closed; the corresponding coordinate is $-1/24$. In the elliptic action, varying the Hessian cofactor in this vector wave gives the $(1,1)$ entry $e^{z^2}$. Its contraction with the scalar Hessian and the root divergence gives $e^{z^1+z^2}$ before the signed factor $1/24$. These binary Hamiltonians fail to induce operations on the linearized field cohomology. The two-input Ward equation on a primitive therefore requires compensating quantum terms of exactly the measured magnitude. [Planar charges and periods](#planar-hamiltonian-waves-and-their-cyclic-composition) · [Primitive and elliptic cofactor replay](quadratic_linearized_complex.py)

For a nonquadratic primary background, the local square can already be nonzero. On $B=\rho(z_1w_1+z_2w_2)+\epsilon P(w_3)$ take $Y=z_1w_1^2\partial_{w_1}\wedge\partial_{w_3}/2$. The elliptic vertex with one level-one field is its trace times the Hessian determinant divided by $24$. The all-primary elliptic vertex vanishes. Thus the order-$\lambda^2$ primary term in the square on $Y$ is $-d\partial_i\partial_j((\operatorname{adj}H)^{ij}dY)/24$. Here $dY=z_1w_1\partial_{w_3}$, the two mixed $z_1,w_1$ cofactor entries equal $\rho^3\epsilon\,\partial_{w_3}^2P$, and the contracted derivative of the cofactor is zero by symmetry of the third derivatives. The result is
$$
 -\frac{\rho^3\epsilon}{12}\partial_{w_3}^3P .
$$
At $P=w_3^3/6$ this is nonzero. It computes the square of the local Hamiltonian alone; the flat quantum master equation also has regulated BV contractions. [Full five-dimensional cofactor check](quadratic_linearized_complex.py)

#### A sharp tangent dimension for elementary corrections on primary fields

Take a primary scalar quadratic in the $5-d$ transverse coordinates, and additional primary fields whose coefficients depend only on the other $d$ tangent coordinates. On the elementary bulk fields, **every positive-loop Hamiltonian correction vanishes for $1\le d\le3$**, at every arity. For divergence-free tangent polyvectors the tree operation preserves this sector: the mixed bracket with the transverse quadratic and the free divergence both vanish, and the binary output is the divergence of the product. In particular the planar Hamiltonian currents retain their classical bracket on this primary sector.

Indeed let $k$ be the number of tangent fields. At genus $g\ge1$ the absorbed primary root leaves $2g+2-k$ transverse quadratics. Eliminate the root momentum by conservation. Each of the $2g$ five-dimensional momentum determinants requires $5-d$ transverse components. The quadratics can absorb at most $2(2g+2-k)$ derivatives at any position, hence $k\le2+(d-3)g$. Each determinant also needs $d$ independent tangent momenta, hence $k\ge d$. These inequalities are incompatible in all three dimensions. Fourier polarization proves the assertion on arbitrary tangent jets. The field-independent term is zero on a scalar background. At the elementary mark a positive geometric loop with a bridge has zero current minor, as in the nine-primary calculation above.

In four tangent dimensions the threshold is attained on **closed, nonexact primary inputs**. Set $B=(z^5)^2/2$. Take three scalar waves and one divergence-free vector wave with tangent momenta $p_1,p_2,p_3,q$ and vector polarization $v$, $q(v)=0$. At the transverse origin each quadratic must absorb two derivatives. For $k=4$ this selects $g=2$. The two forgotten normal marks give a fourfold primary vertex times the squared unit transverse Hessian. In the elementary heat convention the scalar component is
$$
 (2\pi i)^{10}Q_2\phi(0)=
 \frac{(p_1+p_2+p_3)(v)}{62208}\det(p_1,p_2,p_3,q)^4.
$$
On the four distinct coordinate covectors with $v=\partial_{z^1}$ this gives $1/62208$. Evaluation of the primary fields at $z^5=0$ kills every linearized boundary because $B=dB=0$ and the Hessian has rank one. Thus all four inputs have nonzero classes, and this scalar output also passes the boundary test. The four normal Fourier sign choices in the fivefold genus-two graph each give $1/31104$ on this fixture. Their mean and the repeated quadratic factor $1/2!$ recover the result. [Fourfold five-point integral](#all-genus-primary-vertex-from-the-worldsheet) · [Fivefold graph evaluator](longitudinal_theta.py)

The quartic term by itself fails the field Ward identity even on such primary waves. Use scalar momenta $e^1,e^2,e^3$ and two transverse vector waves of momenta $q=e^4$, $s=e^1-e^2$, both with polarization $v=\partial_{z^1}+\partial_{z^2}$. Give the waves at $q,s$ two odd parameters on the left in that order. Contract the inner product by its summed momentum and then the outer product by the total momentum. An inner block with both vector waves vanishes by $v\wedge v=0$. If the inner block contains just the $q$ vector, its odd derivative passes its parameter and gives a minus; the $s$ choice reverses the order of parameters. The total momentum evaluates to $2$ on $v$. In units of $1/62208$ the quartic-inside-binary terms sum to $-4$. For the binary-inside-quartic terms a fusion of the $q$ wave with the first or second scalar gives determinant of magnitude one in the outer quartic. Fusing $s$ there instead gives magnitudes two and zero. The third scalar has momentum zero on $v$. Thus the second sum is $2(-1-1+2^4)=28$. The oriented sum of these two nested divergences is
$$
 \frac{-4+28}{62208}=\frac1{2592}.
$$
The compensating BV graphs can be integrated independently. In the fivefold eight-primary theta integral give the two normal momenta each sign of $e^5$, use the same five tangent waves, and balance the eighth momentum at the root. The odd trace selects the difference of the grounded potentials at the $q$ and $s$ slots evaluated on $v$. The other six derivatives either act on scalars or leave the two parallel vectors wedged together. The magnitude of the root pairing is two, canceled by the $2!$ for the repeated normal quadratic. In path-length order $(1,5,5),(2,4,5),(3,3,5),(3,4,4)$ all four Fourier sign choices give the same result for each graph shape. The bare theta defect carries the opposite potential-contraction sign in the eight-primary Ward equation. With the common heat weight removed the complete test reads
$$
 \frac{-4+28}{62208}
 -\frac{131}{1632960}-\frac{67}{408240}
 -\frac{1}{12096}-\frac{1}{17010}=0 .
$$
This is the flat two-loop quantum Ward identity on nonexact primary states. The isolated quartic operation fails the field Ward test by $1/2592$; the four geometric theta shapes restore it exactly. [Independent fivefold contractions after every fusion and original theta integrations](quadratic_linearized_complex.py) · [Grounding and relative sign](#the-genus-two-bulk-action-and-integrated-cubic)

#### Three nonexact primaries on a two-normal scalar background

Take $x=z^3$, $y=z^5$ and $B=y^2/2+h(x)$ with $h(0)=h'(0)=0$. Keep the coefficient linear in each of the two scalar fields $e^{z^1},e^{z^2}$ and the vector $e^{z^4}\partial_{z^1}$. In the flat Li–Wang quantization the elementary scalar $\phi(0)$ measures a nonzero $Q_g$ on the **same three primary states** for every $g\ge2$ when $h=x^2/2+x^4/24$. The three states are closed and nonexact in the regular linearized field complex. Every positive-genus term in that linearization, including with one incoming descendant, has all its momenta in the two normal directions and the one fluctuation direction. Its fifth Betti direction is missing. The tree primary boundary from a level-$r$ field $Y$ is $\partial_\Omega(B^{r+1}Y/(r+1)!)$ and vanishes on $x=y=0$. The three fields have nonzero restrictions there; their free divergences and contractions with $dB$ vanish. This also holds for every smooth Dolbeault boundary in antiholomorphic degree zero.

Use $I_5,I_6$ for the five- and six-mark integrals in the fourfold reduction at genus $g$, with $I_5=1/62208,\ I_6=0$ at genus two. Put $a=h''(0)$, $c=h'''(0)$ and $b=h''''(0)$. For the indicated polarized coefficient the result is
$$
(-1)^g(2\pi i)^{5g}Q_g\phi(0)=
\begin{cases}
 b/62208,&g=2,\\
 a^{g-2}b I_5+
 a^{g-3}c^2\big((g-2) I_5-I_6/2\big),&g\ge3 .
\end{cases}
$$
For arbitrary tangent momenta $p_1,p_2,q$ on the scalars and vector with polarization $v$ and $q(v)=0$, multiply the expression by $(p_1+p_2)(v)\det(p_1,p_2,q)^{2g}$ in the coordinate order $(z^1,z^2,z^4)$. The all-primary genus-one term is zero.

The $2g+3$ action slots include the three tangent fields and the absorbed root. At the origin the $2g$ normal $y$ derivatives require $g$ quadratic $y$ slots. Forgetting those slots gives the fourfold primary on four hard marks and $g-1$ axial marks, with the Abel degree $g!$ canceled by the ordinary Taylor divisor. Write $t_i$ for the strengths on the axial marks, with the recoil on the absorbed hard mark. After dividing by $\prod_i t_i^2$ the remaining polynomial is a symmetric quadratic, by double softness and homogeneity. It is positive whenever all real $t_i$ are nonzero by the fourfold rank criterion. Its matrix is positive semidefinite by continuity, and permutation symmetry gives identical diagonal entries. A zero diagonal would annihilate the matrix. The common diagonal coefficient is therefore positive; the soft-fiber formula gives $(g-2)!\,I_5$. The $t_i t_j$ coefficient is $2(g-2)!\,I_5-(g-3)!\,I_6$ for $g\ge3$. Assign the two remaining derivatives to one fourth jet, or one to each of two third jets. The multiplicities $g-1$ and $\binom{g-1}{2}$ divided by $(g-1)!$ give the displayed answer. The reduced three-loop fivefold graph code returns the two period columns before the heat conversion $(-1)^g$. [Fourfold periods](#reduction-of-the-fourfold-period-to-five-and-six-worldsheet-marks) · [Direct two- and three-loop fivefold checks](primary_normal_jets.py)

This matrix element is invariant under finite translation- and $SL_5$-equivariant local holomorphic redefinitions of the marked scalar preserving its classical value, at fixed bulk action and output coordinates. A contributing genus-$s>0$ change has $2s+1$ scalar arguments at derivative degree $10s$. Its momentum determinants require three arguments independent modulo the normal plane, leaving at most $2s-2$ normal slots. At geometric loop number $\ell>0$ a planar injection on an ordinary cubic produces a column difference in the tensor product of the cycle space with the two-plane. At most $2\ell$ such differences are independent. A genus-$j>0$ bulk trace has $2j+3$ conserved arguments; at least four must lie outside the plane to span its three-dimensional quotient while summing to zero there. Thus the total capacity of the graph is at most $2g-2$ normal slots. The measurement uses $2g-1$. In the exceptional classical cubic tree the vector can only differentiate the first hard scalar, leaving just two tangent momenta after fusion. The free descendant variation is zero by transversality. These identities hold for arbitrary normal Fourier strengths before extracting the quadratic, cubic or quartic jets. [Graph and grading argument](#a-four-loop-measurement-on-the-elementary-scalar)

**A volume-preserving classical comparison and an exact four-loop value.** Choose $a=c=1,\ b=5/3$. The polynomial $h=x^2/2+x^3/6+5x^4/72$ and the analytic germ specified by $x=u-u^2/6,\ h=u^2/2$ give identical ternary measurements at every genus, since only their first four jets enter. For the analytic germ let
$$
X=u\cos(y/3)+y\sin(y/3),\qquad
Y=-u\sin(y/3)+y\cos(y/3).
$$
Then $X^2+Y^2=u^2+y^2$ and $dX\wedge dY=(1-u/3)du\wedge dy=dx\wedge dy$. The map acts trivially on the three tangent coordinates. In the original affine quantization the exactly quadratic profile $(x^2+y^2)/2$ gives zero, whereas this classical pullback of the quadratic gives, in reduced units,
$$
\begin{aligned}
(-1)^g(2\pi i)^{5g}Q_g\phi(0)
 &= (g-\tfrac13)I_5-I_6/2\\
 &= \frac{W_g}{(g-1)!}+\frac23 I_5
 =-\frac16\left(\frac{3g-1}{(g!)^2}K_g+J_g\right)>0 .
\end{aligned}
$$
Here $W_g$ is the fourfold vertex on unit axial strengths and $K_g,J_g$ are the two integrated four-hard resistance periods of the threefold. The first two nonzero reduced values are
$$
\frac5{186624},\qquad
\frac{3399295819}{1271059485696000}
-\frac{1076227}{490377888000}\zeta(3).
$$
At genus four the fifty-graph period sum gives
$$
\frac{100224698200645422317757433}{20484913263746899968000000000}
+\frac{31025870232372083393461}{2370939035155891200000000}\zeta(3)
-\frac{93558855827167151}{4704244117372800000}\zeta(5).
$$
The difference detects the use of the classical substitution with the affine quantum action held fixed. Comparing quantizations transported by the nonlinear volume-preserving map also requires its quantum action on the fields and observables. [Integrated resistance-period identity](../kodaira-spencer/README.md#the-planar-fivefold-period-measures-zhangs-graph-invariant) · [Coordinate Jacobian and three period columns](primary_normal_jets.py)

**General transverse jets.** On a general two-normal germ with zero constant and first jets, $2g-1$ normal inputs must absorb $4g$ derivatives. Each absorbs at least two. The allowed terms therefore contain one fourth derivative or two third derivatives. For a singular Hessian of rank one, all $2g-2$ or $2g-3$ quadratic slots have collinear momenta; more than $g$ on a line gives zero. Thus every such ternary vanishes at $g\ge4$ at that point. For a nonsingular Hessian $H$, determinant covariance reduces the normal quadratic to the unit form, multiplying the result in that frame by $(\det H)^g$. A rank-one fourth jet $u^{\otimes4}$ gives $I_5(u_i(H^{-1})^{ij}u_j)^2$ inside this factor. For two separately labeled cubics $u^{\otimes3},v^{\otimes3}$ there are just two possible contractions. In the following formula a dot means pairing covectors with this inverse Hessian. The coefficients in increasing genus, in the basis $1,\zeta(3)$, are

| Genus | $(u\cdot u)(v\cdot v)(u\cdot v)$ | $(u\cdot v)^3$ |
|---|---|---|
| $2$ | $1/77760$ | $-1/77760$ |
| $3$ | $28153087/21727512576000-(2971/2794176000)\zeta(3)$ | $-500673587/282457663488000+(158827/108972864000)\zeta(3)$ |

Indeed orthogonal invariance pairs all six slots using the inverse Hessian. For a repeated cubic the ordinary Taylor factor is $1/2!$. Fourier signs at the two labeled cubics extract precisely three derivatives at each mark; binomial multiplicities sum the quadratic insertions on both normal axes. Two normal frames fix the entries and a third checks them independently. For example the real harmonic cubic $(x^3-3xy^2)/6$ on the unit Hessian has zero trace and squared third-derivative tensor norm four. With zero fourth jet the reduced two-loop response is $-1/38880$. [Full fivefold jet projection](primary_normal_jets.py)

Both cubic contractions have an all-genus reduction. Keep the four-hard Betti measure $\Theta_g(\bar P)^{3g}/(3g)!$ in every integral below. At the two cubic marks $1,2$ and the balancing mark $b$ write $A=\Theta_g(z_1-z_b)$, $D=\Theta_g(z_2-z_b)$ and $C=A+D-\Theta_g(z_1-z_2)$. The labeled cubic multipliers are $\int A D C$ and $\int C^3/6$. To see the factorials, divide the $2g-3$ quadratic slots on the unit Hessian into groups of $n_1,n_2$ along its two axes. The fiber factors $n_1!n_2!$ cancel the ordinary Taylor divisors. There remain $g-n_i$ period directions on axis $i$, totaling three. The binomial expansion over these two residual powers gives the degree-three divided power of $t^2(u\cdot u)A+ts(u\cdot v)C+s^2(v\cdot v)D$. Extract $t^3s^3$. Any single-axis exterior power beyond $g$ already vanishes. In particular this derivation includes $g=2$.

Use $h_{ij}$ for the Arakelov diagonal curvature, $e_i^A$ for the tangent curvature at mark $i$, and $e_1^A$ on the *unmarked base* for the fiber integral of the squared tangent curvature, as in the earlier period reduction. In the next display use $z$ for the first auxiliary mark to distinguish its point curvature $e_z^A$ from the unmarked class. The diagonal convolution $\int_y h_{zy}h_{by}=h_{zb}$, the diagonal self-pairing $\int_y h_{by}^2=e_b^A$ and $\int_y h_{by}e_y^A=e_b^A$ give
$$
\begin{aligned}
 \int_y DC &=(g-1)A+\tfrac g2(e_z^A-e_b^A),\\
 \int_z A^2 &=e_1^A/4-ge_b^A,\qquad
 \int_z A(e_z^A-e_b^A)=-e_1^A/2,\\
 \int_{z,y}ADC &=-e_1^A/4-g(g-1)e_b^A .
\end{aligned}
$$
Here each single-coordinate Betti form is $h_{ij}-(e_i^A+e_j^A)/2$. Pairing the last line with the four-hard measure gives $-2I_5-gK_g/(g!)^2$. For two aligned cubics the sum of both multipliers is $2(g-2)I_5-I_6$ by the axial formula. The period identity above gives the second multiplier $-4I_5/3+(K_g/(g!)^2-J_g)/3$. [de Jong–van der Lugt, Lemmas 5.3–5.4 and equation (6.2)](https://arxiv.org/pdf/2207.05428) · [Four-hard curvature pairings](#four-hard-diagram-formula-for-the-remaining-elementary-scalar-period)

In a unit normal Hessian frame, write $T$ and $U$ for the third and fourth derivative tensors of $B$ at the origin. Trace and squared norm in this formula mean contraction by the inverse Hessian, with no conjugation. The entire measured coefficient becomes
$$
\begin{aligned}
(-1)^g(2\pi i)^{5g}Q_g\phi(0)
={}&I_5\left(\operatorname{tr}^2U-|\operatorname{tr}T|^2-\tfrac23|T|^2\right)\\
 &-\frac{gK_g}{2(g!)^2}|\operatorname{tr}T|^2
   +\frac16\left(\frac{K_g}{(g!)^2}-J_g\right)|T|^2 .
\end{aligned}
$$
For a real positive unit Hessian the first parenthesis is minus eight times the relative first correction in the ordinary transverse Gaussian integral. Indeed the average of $T(w,w,w)^2$ with unit Gaussian covariance is $6|T|^2+9|\operatorname{tr}T|^2$, and that of $U(w,w,w,w)$ is $3\operatorname{tr}^2U$. The two Taylor coefficients in the exponential are $1/72$ and $-1/24$. A local area-preserving pullback of the quadratic preserves this stationary-phase coefficient, which is zero for the quadratic.

For the polynomial shear $B=(x^2+(y+x^2)^2)/2$ the three contractions $\operatorname{tr}^2U,|\operatorname{tr}T|^2,|T|^2$ are $12,4,12$. Its reduced ternary is $-2(g-1)K_g/(g!)^2-2J_g$, equal to $11/77760$ at two loops. Subtract twelve copies of the preceding analytic pullback measurement to obtain $4gK_g/(g!)^2$. The two changes of normal coordinates thus measure independent linear combinations of both four-hard graph periods on the same three primary states. The comparison fixes the classical backgrounds and the affine bulk quantization. [Genus-two and three fivefold jets, four-loop period contraction and polynomial shear](primary_normal_jets.py)

**Four-input Ward test.** On $h=x^2/2+x^4/24$ retain the two scalar waves and take the two vector momenta $q=e^4,\ s=e^1-e^2$, both polarized along $\partial_{z^1}+\partial_{z^2}$, with left odd parameters ordered $q,s$. Each positive-loop binary on these tangent waves vanishes by the two normal and two tangent momentum directions. Composing two positive genera at this arity gives zero even when the inner result has a normal profile, since its tangent momentum is still the sum at its inputs. The tree unary vanishes on these divergence-free tangent products. The nested ternary and tree binary, in the same orientation as the preceding quartic test, give
$$
-4 I_5+2(4^g-2)I_5=2(4^g-4)I_5 .
$$
The unit determinants in the first ordering and the determinant magnitudes $1,1,2,0$ after tree fusion give the two terms. At $g=2$ this is $1/2592$. At $g=3$ the original genus-two insertion in the hexagon gives $1/2177280$ in the compensating potential projection. Subtracting it from $120I_5$ fixes the combined bare tetrahedron and double-square projection to
$$
\frac{738845099}{3530720793600}-\frac{469153}{2724321600}\zeta(3).
$$
The last pair is inferred for this quartic background. For an independently integrated three-loop comparison, take the repeated-group fixture of the bulk-potentials section, with three equal $e^1$ background waves and vector slots at $e^2,e^4$, polarized by $\partial_{z^3}$. Here the background is the single quadratic in the first coordinate and there are six tangent primaries. Dividing by $3!$ and including the root magnitude two, direct evaluation of the fourfold nested vertices and the three original topology columns gives the same two coefficients:
$$
-\frac{23310136147}{423686495232000}+\frac{459733}{10216206000}\zeta(3)
=\frac13\left(-\frac{9261411221}{56491532697600}
+\frac{459733}{3405402000}\zeta(3)-\frac{829523}{1975228416000}
-\frac1{1451520}\right).
$$
The signed elementary Hamiltonian alone thus fails the Ward test on the ternary probe. Its extension to a physical operator requires the compensating BV graphs. [Three integrated topologies](#three-loop-bulk-potentials-outside-the-six-planar-locus) · [Exact nested operations and same-background hexagon](primary_normal_jets.py)

**Ward identity for a change of normal coordinates.** The change of the ternary under a nonlinear normal coordinate transformation can itself be computed from the original regulated graphs. At unit normal Hessian polarize the cubic tensor as $T=u^{\otimes3}$, and take the infinitesimal normal vector $\xi(w)=l(w)^2 v/2$ with $l(v)=0$ and $w=(x,y)$. Retain the same three hard tangent states. Quadratic divergence-free normal vectors are spanned by these shears: their Hamiltonians are cubic binary forms and cubes of linear forms span that space.

Take the tangent vector before $\xi$ in the order of the two odd sources. At genus $g$ the geometric graph has $2g$ normal sites, including $\xi$ and the cubic. A geometric loop of order $\ell$ admits at most $2\ell$ normal jumps at ordinary cubic vertices. A positive-genus action insertion of genus $j$ admits at most $2j-1$ normal external momenta. Hence every graph of positive geometric loop order with such an insertion vanishes. This excludes the three-loop hexagon separately even before integrating. In a composition of two positive-genus vertices, the inner vertex requires all three hard tangent momenta. The outer vertex then has only their sum as a tangent direction. The classical generator vanishes with its first jet at the origin and fixes the hard waves. The remaining local variation at the measured point is the derivative of the ternary in the direction $\xi(B)$.

In the bare graph the four retained potentials can be set to zero using the integrated planar resistance formula, modulo common and longitudinal shifts. The two-vector contraction with the absorbed root is invariant under those shifts. The tangent root pairing on these hard waves is one. At the normal vector slot the extra voltage momentum must come from the cubic to supply its third derivative; all normal slots already supply two powers in the squared current cofactor. Integrate the $2g-2$ quadratic donors with their ordinary Taylor divisor. The matching mass $Z$ for the remaining two labeled normal slots is
$$
 [ab]\det(\mathbf1+a\,ll^{\mathsf T}+b\,uu^{\mathsf T})^g
 =g^2(l\cdot l)(u\cdot u)-g(l\cdot u)^2 .
$$
The mass after deleting both slots is one. In the normal potential the $K_g$ column multiplies $Z/[g(g!)^2]$ and the $J_g$ column multiplies $(Z-g^2[l,u]^2)/[g(g-1)]$, where $[l,u]$ is the two-row determinant. This gives the signed geometric Ward contribution
$$
u(v)\left[
\frac{K_g}{(g!)^2}\big(g(l\cdot l)(u\cdot u)-(l\cdot u)^2\big)
+J_g(l\cdot u)^2\right].
$$
The classical variation of the cubic is the sum of the three placements of $v^\flat$ in $l\otimes l\otimes v^\flat$, where the unit Hessian identifies $v$ with $v^\flat$. Thus the variation of its trace is $(l\cdot l)v^\flat$. The ordinary Gaussian stationary-phase coefficient is invariant under the divergence-free flow. Differentiating the all-genus ternary formula gives the exact negative of the displayed geometric contribution. If the generator instead has cubic leading jet $l(w)^3v/6$, the extra derivative must land at that slot and contributes $l(v)=0$; higher jets exceed the normal derivative degree. This proves the projected quantum Ward identity for the normal coordinate changes tangent to the identity, for every $g\ge2$. [Original two-period matching transfer](planar_bulk.py) · [Separate graph columns on oblique normal jets through genus six](primary_normal_jets.py)

For the finite shear $B_\tau=(x^2+(y+\tau x^2)^2)/2$ its reduced ternary equals $\tau^2$ times the value at $\tau=1$. In the coefficient linear in $\tau$ the local two-loop derivative and the original four theta shapes give
$$
\frac1{2592}-\frac1{19440}-\frac1{19440}
-\frac{19}{326592}-\frac7{58320}-\frac{11}{181440}-\frac1{22680}=0 .
$$
In the first three entries $\xi=x^2\partial_y$ fuses with the cubic and the two labeled quadratic donors, respectively. The tangent and normal vector potential contractions both contribute in every original theta shape. The finite change of the local two-loop ternary from $\tau=0$ to $\tau=1$ is $11/77760$; integrating the signed geometric Ward contribution gives $-11/77760$. This is a covariance test in the twisted quantum theory. A nonzero coordinate variation of the elementary matrix element here is compatible with the full Ward identity. [Independent original theta shapes](normal_shear_ward.py) · [Literal individual fusion derivatives](primary_normal_jets.py)

On the same normal shear at three loops the genus-two hexagon is zero. Direct integration of the three four-hard graphs with double-square core and the six with tetrahedral core, including the suppressed normal placements, gives respectively
$$
-\frac{347}{6858432000},\qquad
-\frac{370831}{14370048000}+\frac{22219}{1047816000}\zeta(3).
$$
Their sum equals $2K_3/9+4J_3=-7802719/301771008000+(22219/1047816000)\zeta(3)$, exactly the negative of the local shear variation. Both topology integrals were computed from their graph numerators and both tetrahedron cuts. [Nine integrated diagrams](normal_shear_ward.py)

**Finite transport of the two-normal matrix element.** For a nonsingular normal Hessian $H$, contract all normal indices with $H^{-1}$. Under a holomorphic area-preserving change $f$ fixing the critical point, the determinant of $H$ and the quantity
$$
(\det H)^g\left(\operatorname{tr}^2 U-|\operatorname{tr}T|^2-\tfrac23|T|^2\right)
$$
are invariant. For real positive Hessians this is the first stationary-phase coefficient of the scalar integral with its Gaussian determinant restored; polynomial continuation gives the identity for all invertible complex $H$. In the reduced ternary the other two terms constitute the following primitive of the geometric Ward variation:
$$
(\det H)^g\left[-\frac{gK_g}{2(g!)^2}|\operatorname{tr}T|^2
+\frac16\left(\frac{K_g}{(g!)^2}-J_g\right)|T|^2\right].
$$
The signed integrated graph correction from $B$ to $B\circ f$ is its value at $B$ minus its value at $B\circ f$. This follows by integrating the infinitesimal shear identity: the second jet of an area-preserving change tangent to the identity is a divergence-free quadratic vector, the third jet contributes zero after its required quadratic composition, and the linear $SL_2$ part acts by covariance. Equivalently the local ternary has the explicit split just displayed, and the invariant part takes equal values at the endpoints. For two successive changes the intermediate primitive cancels exactly. This constructs transport on this projected interaction; extension to the full complex of local observables is a further problem.

The invariant component with its prescribed values on $T=0$ is unique on the nonsingular locus. In a unit Hessian frame one can kill an arbitrary cubic by a homogeneous cubic Hamiltonian generating a quadratic coordinate change. In null coordinates the Hamiltonian derivation of the quadratic has eigenvalues proportional to $i-j$ on monomials of degree $i+j=3$, all nonzero. The induced fourth jet then has the displayed double trace. This invariant and its complementary primitive both extend polynomially across singular Hessians for $g\ge3$, since at most three inverse Hessians enter. For $g=2$ a regular extension of this *normalized split* is impossible. On $B=(x^2+t y^2)/2+y^3/6$ the invariant without its factor $I_5$ equals $-5/(3t)$ at $g=2$. Its pole cancels in the full local coefficient and in the finite difference under a regular coordinate change.

At a rank-one Hessian choose a frame with $H=\operatorname{diag}(a,0)$, $a\ne0$. Subscripts here denote ordinary derivatives at the critical point. The full two-loop and three-loop reduced ternaries are respectively
$$
\frac{a^2 B_{yyyy}}{62208}
 +\frac{a(B_{xxy}B_{yyy}-B_{xyy}^2)}{77760},\qquad
a^3 B_{yyy}^2(I_5-I_6/2)
$$
with $I_5,I_6$ taken at genus three in the second expression. All higher genera vanish. Set $r=\partial_y^2 f^x(0)$ and $s=\partial_x\partial_y f^x(0)$ for a change tangent to the identity. Area preservation gives $\partial_y^2 f^y(0)=-s$. The shifts of $B_{xxy},B_{xyy},B_{yyy}$ are $2as,ar,0$; the fourth derivative $B_{yyyy}$ shifts by $6(B_{xyy}r-B_{yyy}s)+3ar^2$. Hence the local two-loop variation is
$$
 \frac{11a^2}{155520}\left(B_{xyy}r-B_{yyy}s+\frac{ar^2}{2}\right),
$$
and the local three-loop variation is zero. For example the shear $x\mapsto x+r y^2/2$ on $B=x^2/2+y^3/6$ gives two-loop value $11r^2/311040$ while the three-loop coefficient stays $I_5-I_6/2$. [Polynomial composition and singular-limit checks](normal_coordinate_transport.py)

#### Planar primary brackets and the normal elliptic curve

Let $x,y,z$ be the three normal coordinates to the tangent $(z^1,z^2)$ plane. Take a homogeneous cubic primary scalar background $B(x,y,z)$. In the same fivefold elementary scalar convention write the reduced coefficient $(-1)^g(2\pi i)^{5g}Q_g\phi(0)$ on one scalar $e^{z^1}$ and one vector $e^{z^2}\partial_{z^1}$ as the value of the planar binary at genus $g$. For brevity below $Q_g|_B$ denotes this *original heat-convention* matrix element, including its heat factor. The genus-one primary term is zero. Every positive-genus primary term on divergence-free tangent fields is binary and its coefficient is constant in the normals. In fact $k$ tangent inputs leave $2g+2-k$ background inputs. The $2g$ five-row determinants need $6g$ normal derivatives while the cubic inputs supply at most $3(2g+2-k)$. Rank in the tangent plane requires $k\ge2$. Thus $k=2$ and each of the $2g$ normal slots supplies exactly three derivatives. With only the normal background and one variation even a descendant vertex sees at most four momentum directions, so the local linearized differential is the classical one. On $x=y=z=0$ every regular holomorphic primary boundary from level $r$ is the divergence of $B^{r+1}Y/(r+1)!$, which vanishes. The two stated tangent waves are closed and nonexact.

On the Hesse family $B=(x^3+y^3+z^3)/6+bxyz$ the original seven- and nine-mark fivefold graphs give
$$
\begin{aligned}
 (2\pi i)^{10}Q_2|_B&=\frac{b^4-b}{51840},\\
 -(2\pi i)^{15}Q_3|_B&=(1-20b^3-8b^6)
 \left(\frac{1163706689}{25609494822912000}
       -\frac{74897}{1976041267200}\zeta(3)\right).
\end{aligned}
$$
For $r$ copies of $xyz$ use $(2g-r)/3$ pure cubes on each axis, an integer forced by multi-degree. The Fourier jets $\sinh x$ and $\sinh x\,\sinh y\,\sinh z$ implement $x^3/6$ and $xyz$: the linear jet in an individual slot vanishes by double softness, and the total derivative degree excludes all higher odd jets. Weight the signed waves by their ordinary Taylor divisors, quotient by normal-axis permutations and reflections, and evaluate the original reduced graph polynomials on the orbit representatives. Independent Weierstrass-form evaluations give $(2\pi i)^{10}Q_2=-1/1080$ at $B=y^2z-x^3-xz^2$, and reproduce the sextic on $B=y^2z-x^3-z^3$ in both period columns. In particular the Fermat three-loop value is nonzero by irrationality of $\zeta(3)$. [Orbit sums and independent graph fixtures](planar_cubic_primary.py)

Both measurements survive finite local translation- and $SL_5$-equivariant holomorphic changes of the marked scalar with its classical value, bulk action and output coordinates fixed. The free descendant variation is the divergence of the transverse vector. In a classical fusion the vector differentiates only the hard scalar and leaves a single tangent momentum at the mark. A genus-$s>0$ marked change contributing at this degree has $2s+1$ primary scalar arguments and $10s$ derivatives. At genus one the three arguments cannot form a five-row determinant. In the only remaining three-loop possibility the genus-two mark in a one-loop hexagon has exactly five scalar slots. Its degree-twenty symbol is a multiple of the fourth power of their determinant. At least two of its three spectators are normal cubics, each differentiated four times. Such a graph vanishes before integrating any heat time. [Marked degree and primary graph count](#the-three-loop-differential-on-the-elementary-scalar)

**Elliptic modulus from the two measured loop coefficients.** Divide each measured coefficient by an actual nonzero value of the same genus. Put $s=Q_2|_B/Q_2|_{xyz}$ and $t=Q_3|_B/Q_3|_{(x^3+y^3+z^3)/6}$. On the Hesse family these ratios equal $b^4-b$ and $1-20b^3-8b^6$. For every homogeneous ternary cubic,
$$
 t^2-64s^3=4096\,\operatorname{Res}(\partial_x B,\partial_y B,\partial_z B).
$$
Here the resultant of three homogeneous quadratics is normalized to one on $x^2,y^2,z^2$. Its zero locus means the normal gradient has a common nonzero zero. On the Hesse family $t^2-64s^3=(1+8b^3)^3$. To compute the resultant directly, set $z=1$ in the first two doubled gradients. Away from $b=0$ these have the origin and three solutions with $xy=4b^2$, counted with multiplicity. There is no common point at $z=0$. The product of the last doubled gradient is $(1+8b^3)^3$. Restoring all three factors $1/2$ supplies $2^{-12}$. Under $B\mapsto B\circ M$ for invertible linear $M$, the genus-$g$ graph polynomial acquires $(\det M)^{2g}$ from the momentum determinants. The gradient resultant acquires $(\det M)^{8+4}$. The Hesse family and these changes are dominant: at the Fermat cubic the nine first-order linear changes give the nine monomials $x_i^2x_j$, and varying $b$ adds $xyz$. Polynomiality proves the identity on the full ten-dimensional space of cubics, including singular ones.

For a smooth normal projective cubic $B=0$, the two loop coefficients therefore determine
$$
 j=-\frac{110592s^3}{t^2-64s^3}.
$$
Indeed in Hesse form set $k=-2b$, use the flex $(1:-1:0)$, and put $X=z,\ Y=x-y,\ Z=x+y+kz$. In the chart $Z=1$ the equation is $Y^2=4(k^3-1)X^3/3-3k^2X^2+2kX-1/3$. Translating and rescaling to a monic cubic gives the Weierstrass coefficients $-16s/3,-16t/27$. The usual $j$ of that equation is the displayed ratio. On the smooth cubic locus equality extends from the dense Hesse charts since both sides are regular elliptic invariants. At $j=0$ the first nonzero planar primary term is genus three. At every other smooth $j$ it is genus two. At $j=1728$ the three-loop term vanishes. This concerns the chosen field matrix elements in the twisted theory. The regular D3 brane has an independent calculation: its obstruction to extending a commuting displacement $u$ is $B(u)\mathbf1$, so its projective curve of extendable displacements has this same $j$. [D3 matrix factorization and signed cubic obstruction](../twisted-holography-N4/README.md#the-cubic-normal-obstruction-on-a-regular-d3-brane)

The regular normal D3 disk with ordered boundary insertion $\eta_1\eta_2\eta_3$ evaluates a holomorphic scalar $f$ at the normal origin. For smooth $B$ its residue numerator is $f\det B''/8$ and the gradient map has local degree eight. This returns both inserted loop coefficients. For the displaced factorization at $\tau u$, $B(u)=0$, the identity is killed by exactly $(\tau^2)$ if $dB(u)\ne0$ and survives after inverting $\tau$ if $dB(u)=0$. The loop discriminant therefore detects the appearance of noncontractible generic displacements. A singular finite Hesse member has three such projective directions. At the marked flex $u=(1,-1,0)$, a further displacement by $-\sigma(2b,0,1)$ gives scalar curvature $\sigma^3(1+8b^3)/6$ relative to $B(w)$. In this Hesse frame, $216$ times the cube of this open cubic obstruction coefficient is exactly $t^2-64s^3$. The singular-brane supercharge requires the puncture and open transfer. [Regular disk, flex and endomorphism products](../twisted-holography-N4/README.md#the-cubic-normal-obstruction-on-a-regular-d3-brane)

The same two outputs have an independent test on the unfluxed punctured fivefold. For a formal strength $\rho$ in $W=\rho B$, the wave $e^{z^1+z^2}$ has annihilator exactly $(\rho^3)$ there on a smooth cubic. Both planar outputs carry at least four powers of $\rho$ and are therefore exact after restriction. At a singular finite Hesse member both coefficients remain nonzero and their scalar classes survive on a punctured critical ray. [Explicit Čech primitive and sharpness](../twisted-holography-N4/README.md#the-planar-scalar-on-the-unfluxed-punctured-cubic)

With nonzero D3 flux formal in $N$, a closed continuation of $\rho e^{z^1+z^2}$ exists on every Hesse member. This lift has annihilator exactly $(\rho^2)$ on smooth members. Both scalar loop outputs thus admit exact continuations there. On a singular finite Hesse member, evaluation at a remaining critical point detects the leading scalar of every closed continuation, even with flux. In particular at $b=-1/2$ the reduced two-loop value is $1/92160$ and the three-loop value is $27/8$ times its Fermat value. The elementary local binary on the explicit closed scalar lift and $e^{z^2}\partial_{z^1}$ equals $\rho^{2g}Q_g|_B$ times the same closed lift of the output wave. Arbitrary hard normal momenta are removed by a determinant-one shear and the remaining root divergence gives all its flux terms. Thus the smooth-versus-singular distinction is attained by the actual local Hamiltonian on these inputs. The complete relative Ward operation, including open terms, remains open. [Explicit flux lift, smooth homotopy and singular detector](../twisted-holography-N4/README.md#the-planar-scalar-on-the-unfluxed-punctured-cubic)

There is a short exact evaluator for general coefficients. Write $D=\det(\partial_i\partial_j B)$ for the Hessian cubic in the specified normal coordinates. With parameter $a$ the Hessian pencil obeys the polynomial identity
$$
 \det\partial^2(B+aD)=D+12asB+a^2(-6tB-12sD)+a^3(48s^2 B+2tD).
$$
On the Hesse family $D=-b^2(x^3+y^3+z^3)+(1+2b^3)xyz$. Substitution gives the identity coefficient by coefficient; determinant covariance and dominance give the general result. Extract $s,t$ from any nonzero monomial of $B$ using the linear and quadratic coefficients of the pencil. This also works when the resultant vanishes. On a Weierstrass cubic $B=y^2z-x^3-Pxz^2-Rz^3$ at arbitrary $g$, exactly $g$ marks carry $y^2 z$. If $i$ further marks carry $Pxz^2$ and $j$ carry $Rz^3$, then $2i+3j=g$ by the normal $z$-degree. Thus each higher genus needs only the $P^iR^j$ coefficients of that weight. In particular the equianharmonic locus $P=0$ permits only $3\mid g$, and the locus $R=0$ permits only even $g$. The remaining coefficients starting at genus four require their own period integrals. [Hessian algorithm and coordinate change checks](planar_cubic_primary.py)

In fact the dependence on the normal cubic at every genus is polynomial in the *two measured ratios* $s,t$. In the short Weierstrass frame the Hessian pencil gives $s=-48P$ and $t=-6912R$. Each allowed $P^iR^j$ transforms after substituting these two ratios with determinant weight $4i+6j=2g$, precisely that of the original graph polynomial. The equality extends from smooth Weierstrass forms and their linear changes to all cubics by polynomiality. The smooth forms are covered because a flex can be moved to infinity and its tangent to the line at infinity; alternatively the explicit Hesse transformation above already covers a dense set. Thus the four-loop shape is $s^2$ times a universal period, and the five-loop shape is $st$ times another. The first genus with two independent coefficients is six, on $s^3$ and $t^2$. The coefficients are unique since $P,R$ vary independently. This determines the allowed cubic dependence, including on singular members. [Weierstrass and Hessian checks](planar_cubic_primary.py)

For the remaining four-loop period only nine conserved eleven-mark evaluations are needed. Choose the four normal linear functions $x,y,z,-x-y-z$, use their cubes divided by six, and multiply the fourth by $h$. Substitution in the pencil gives $s=h$, $t=1-6h-3h^2$. Extract the $h^2$ term. Six marks come from the first three axes and each axis can supply at most eight derivatives, so exactly two marks come from each of the four lines. Represent each cubic by its hyperbolic sine. The plus counts on the two marks of each line are $0,1,2$ with weights $1/8,-1/4,1/8$. Permuting the four vectors is a normal linear change of determinant $\pm1$, and reversing all signs has the same invariance. In the orbit order $0000,0001,0002,0011,0012,0022,0111,0112,1111$, the nine weights times $4096$ are
$$
 (2,-16,8,48,-48,6,-64,48,16).
$$
Each eleven-mark vertex retains the two tangent waves and the balancing root. The same projection at genus two returns $1/51840$ for four distinct cubic lines. As an independent normalization check, sum the product of two labeled genus-two vertices over the seventy ordered four-plus-four splits of the eight slots. Applying the nine weights gives exactly $1/51840^2$. The connected four-loop period still requires its own eleven-mark integrals. [Four-line polarization and product control](planar_cubic_primary.py)

**The planar Jacobi test.** For tangent covectors $p,q$ put $[p,q]=p_1q_2-p_2q_1$. Use the Hamiltonian wave $X_p=e^{p\cdot z}(-p_2\partial_{z^1}+p_1\partial_{z^2})$. The genus-$g$ operation on $X_p$ and $e^{q\cdot z}$ is the reduced binary coefficient times $[p,q]^{2g+1}e^{(p+q)\cdot z}$ after removing the heat convention. The same factor multiplies the vector $X_{p+q}$ on two Hamiltonian waves. This follows from the $2g$ tangent determinants and the divergence of the wedge product; the classical bracket has multiplier $[p,q]$. In the cyclic ordering $[[p,q],r]$ let $a=[p,q],b'=[q,r],c=[r,p]$. At the first nonzero loop order the Jacobi polynomial per unit binary coefficient is the cyclic sum of $a^{2g+1}(b'-c)+a(b'-c)^{2g+1}$. For $p=(1,0),q=(0,1),r=(1,2)$ it is $120$ at genus two and $1680$ at genus three. Thus on every smooth cubic the isolated binary fails Jacobi on closed, nonexact tangent waves already at its first nonzero order.

For the same three hard waves the original projective theta integrals on the Hesse family give, in path order $(1,5,5),(2,4,5),(3,3,5),(3,4,4)$ and in the potential projection,
$$
 (b^4-b)\left(\frac{41}{90720},\frac{11}{11340},\frac1{2016},\frac1{2520}\right).
$$
Their sum equals $120(b^4-b)/51840$. With the bare-graph sign in the quantum Ward identity it cancels the full two-loop Jacobiator. For the same projection at the Fermat cubic the full three-loop contributions, including core symmetry divisors, are
$$
\begin{aligned}
\text{double square}&=-\frac{3361}{27433728000},\\
\text{tetrahedron}&=\frac{13904672281}{177843714048000}
-\frac{74897}{1176215040}\zeta(3),\\
\text{one-loop hexagon with genus-two insertion}&=-\frac1{580608}.
\end{aligned}
$$
The double square is integrated directly. On each normal axis the two cubic slots give the difference of their equal- and opposite-sign current moments divided by four, including their repeated-field divisor. With the recoil grounded, reciprocity reduces the projected tangent voltage to $-G(p,q)+4G(p,r)-3G(q,r)$, where $G$ is the grounded graph Green function and its arguments label the three hard sites. After summing path placements, 841 core cubature nodes per prime integrate the full invariant numerator. Nine primes determine the rational period uniquely by the ordered determinant height bound. The genus-three Ward identity fixes the tetrahedron by subtracting the other two entries from $1680$ times the reduced Fermat binary. Independently, 2,540 tetrahedral cubature nodes at the reserved prime $1009$ give $(416,45)$ in the rational and $\zeta(3)$ columns. The double square gives $(782,0)$ there. The two bare contributions have opposite signs, since $\zeta(3)<263/216$ bounds the tetrahedron below by $38711527/59281238016000$. [Original theta and hexagon](planar_cubic_primary.py) · [Fermat period certificate](operator/fermat_jacobi_core_periods.json) · [Rational reconstruction and independent core cubatures](fermat_jacobi_cores.py)

Each of the three topology sums in this projection is homogeneous of degree six in the cubic and transforms by the sixth power of the normal determinant. The separate normal coordinate degrees give $A+C b^3+D b^6$ on the Hesse family. The Fourier matrix $(\omega^{ij})_{0\le i,j\le2}$, with $\omega$ a primitive third root of unity, has determinant to the sixth power $-19683$. It sends the coefficients of $(x^3+y^3+z^3)/6$ and $xyz$ to $3+6b$ and $3-3b$. Equality of constant and linear powers after the transformation gives $C+D=-28A$ and $12A+3C-6D=0$. Thus every entry is its Fermat value multiplied by $t(B)$ on every ternary cubic, by the Hesse dominance argument above. At smooth $j=1728$ all three topology contributions vanish separately. [Fourier character check](fermat_jacobi_cores.py)

**Two fixed nonexact primary states at every genus.** In the same planar geometry take instead the entire primary scalar $B=z^2/2+e^x+e^y-2-x-y$. Its constant and first jets vanish at the normal origin. The rank bound on one variation and the linearized-boundary argument still apply. Exactly $g$ normally differentiated $z^2/2$ slots contribute. Forgetting them cancels their repeated-field factorial. On the remaining $g$ normal scalar slots the constant and linear parts of each axial wave vanish by double softness. With $i$ unit $x$ waves and $g-i$ unit $y$ waves, each of the fourfold vertices with $1\le i\le g-1$ has full rank and all momenta nonzero, including the recoil on the absorbed root. Its contribution is strictly positive with factor $1/(i!(g-i)!)$. Thus on the same two primary states
$$
 (-1)^g(2\pi i)^{5g}Q_g|_B>0\qquad (g\ge2).
$$
The first two reduced values are $1/62208$ and $49188779/16295634432000-(7783/3143448000)\zeta(3)$. The analytic background is fixed as $g$ varies. The earlier Jacobi and elliptic formulas use homogeneous cubic backgrounds. [Independent fivefold graph and fourfold evaluations](planar_cubic_primary.py)

The entire background also gives this same nonzero tower in the linearized punctured bulk with D3 flux. It is invariant under $(x,y,z)\mapsto(x+2\pi i,y-2\pi i,z)$. The translated origin is a critical point off the brane. On the explicit flux-closed lift of $\rho e^{z^1}$ and the vector $e^{z^2}\partial_{z^1}$, the actual local binary has coefficient independent of $z$, the only normal direction of the lift. Its scalar part at the translated point is $\rho^{2g+1}$ times the flat matrix element above, and evaluation there annihilates every linearized boundary. [Flux correction and global critical-point test](../twisted-holography-N4/README.md#the-planar-scalar-on-the-unfluxed-punctured-cubic)

#### A planar binary on a punctured critical line

Take $B=(z^2+y^2h(x))/2$ with $h$ holomorphic and the same two tangent waves. Every point $(x,0,0)$ is critical and has $B=0$. On the critical line exactly $g$ quadratic $z$ marks supply all $2g$ derivatives in that direction. Forget them by the Abel fiber formula. The remaining $g$ slots each supply exactly two $y$ derivatives, and altogether $2g$ derivatives in $x$. The reduced scalar coefficient is therefore a homogeneous differential polynomial of degree $g$ in the derivatives of $h$, with total derivative order $2g$. The same restriction holds on adding any term $y^3H(x,y)$ to $B$. The graph coefficient off the line is still independent of $z$, the only normal-vector direction in the flux lift. Hence the actual elementary binary on the preceding two closed inputs is this scalar graph coefficient times the same closed output lift. [Exact fourfold and fivefold jets](critical_line_binary.py) · [Flux and brane states](../twisted-holography-N4/README.md#regular-d3-displacements-on-the-critical-line)

For $h=1+a e^{px}+b e^{qx}$, with $p,q$ nonzero real frequencies, extract one power of each amplitude. For the reduced binary at unit background strength on the line one obtains
$$
 [ab](-1)^g(2\pi i)^{5g}Q_g\big|_B
 =\frac{6p^2q^2}{(g-2)!}\,C_g(-q e^3+p e^4)\,e^{(p+q)x}\qquad (g\ge2).
$$
Here $C_g(v)$ is the previously defined leading coefficient on four unit hard covectors, $g-2$ equal soft momenta along $v$ and their balancing root. To polarize the two distinguished $y^2/2$ jets, let their normal momenta be $(p,s),(q,t)$. Give each remaining normal $y$ mark momentum $(0,\epsilon)$. Double softness makes the coefficient of $\epsilon^{2g-4}$ exactly that of the individual quadratic jets. In the hard normal basis the soft vector has coordinates $(-q,p)/(pt-qs)$. The hard basis determinant occurs to power $2g$ and $C_g$ has degree $2g-4$. The resulting polynomial is $(pt-qs)^4 C_g(-q e^3+p e^4)$. Its $s^2t^2$ coefficient is $6p^2q^2 C_g(-q e^3+p e^4)$. Only the $g-2$ equal undifferentiated background marks contribute a Taylor factorial. The derivation extends across $pt=qs$ by polynomiality. The sharp fourfold soft theorem gives strict positivity of $C_g$ for this nonzero direction. At $p=1,q=-1$ the result is constant along the critical line. It is $1/10368$ at two loops and $740466719/47076277248000-(469153/36324288000)\zeta(3)$ at three. At four loops it is exactly $3C_4(e^1+e^2)$, one twenty-fourth of the earlier stripped $72C_4(e^1+e^2)$ D3 pole. [Soft rank theorem](#generic-rank-and-sharp-simultaneous-soft-limits) · [Four-loop Ward graphs](#a-ward-projection-for-the-mixed-genus-four-d3-pole)

A fixed ordinary analytic background already suffices for every genus. Set $h=1+e^x$. Its coefficient on the axis is a polynomial in $e^x$ of degree at most $g$. The constant and linear terms vanish by normal rank or by the fourth determinant power in the one-hard-normal limit: with one $x$-dependent mark, its $y$ momentum can be removed by a shear that fixes all pure $y$ slots. The coefficient of $e^{2x}$ equals $3C_g(e^4-e^3)/(g-2)!>0$, since the two exponential marks are identical. Evaluate all genera on the same punctured formal line $(x,y,z)=(\tau,0,0)$ over $\mathbb C((\tau))$. No nonzero polynomial in $e^\tau$ with complex coefficients vanishes there: its distinct exponentials have an invertible Vandermonde matrix of derivatives at zero. Every primary boundary vanishes under this form-degree-zero evaluation, including boundaries from descendants or an explicit flux. Thus both fixed inputs and every $g\ge2$ output are nonexact on this punctured germ. This uses the local Hamiltonian with the explicit flux lift. With $h=e^x-1$ the same $e^{2x}$ coefficient acquires $(-1)^{g-2}$ and the normal Hessian drops rank at the origin.

At two loops the fourfold determinant gives exactly
$$
 (2\pi i)^{10}Q_2\big|_B=\frac{(h'')^2}{20736}.
$$
At three loops the finite jet uses only $h h'' h^{(4)},h(h''')^2,(h')^2h^{(4)},h'h''h''',(h'')^3$. The exact rational and $\zeta(3)$ coefficients are obtained from the symmetric degree-six polynomial on three normal frequencies, by evaluating seven independent rows and three additional mixed-sign tests. For $h(0)=h'(0)=0$, only the last monomial survives at the origin. Its reduced coefficient is $-144754243/30813563289600+(226271/59439744000)\zeta(3)$. At every genus such a double zero selects only $(h''(0))^g$ by derivative count. [Complete polynomial and full graph checks](critical_line_binary.py)

For a normal-coordinate change $X=X(x),\ Y=y/X'(x),\ Z=z$ the Jacobian is one and $h(x)=\widetilde h(X)/(X')^2$. If $h(0)\ne0$ choose $X(0)=0,\ X'=1/\sqrt h$ to give $\widetilde h=1$. The positive-genus planar binary in that affine frame is zero by momentum rank. If $h$ has a double zero, $c=h''(0)/2\ne0$ is invariant. A normal form is $\widetilde h=cX^2$: integrate $d\log X/dx=\sqrt c/\sqrt{h(x)}$ on the branch with residue one. The difference from $1/x$ is regular, so $X/x$ is a holomorphic unit. In particular the preceding nonzero two- and three-loop values at a double zero survive this full reparametrization subgroup.

There is a coordinate-independent contour test in this class of normal splittings. The two Gaussian normal directions at fixed $x$ give the induced one-form $dx/\sqrt{h(x)}$: for positive real $h$, integrating $e^{-(z^2+h y^2)/2}\,dy\,dz$ gives $2\pi/\sqrt h$. Continue the form on a chosen sheet of the square root. It is unchanged by the normal-coordinate transformation just given. In fact any area-preserving change in $(x,y)$ preserving the branch $y=0$ has $X(x,0)=X(x)$ and $\partial_yY(x,0)=1/X'(x)$. The additional terms in $B$ are at least cubic in $y$ and drop out of the critical-line graph at every genus. At $h(0)=h'(0)=0$, $h''(0)=2c\ne0$, choose $\sqrt h/x\to\sqrt c$. The two first nonzero residues of the reduced local binary times this form are
$$
 \operatorname{Res}_0\left((2\pi i)^{10} Q_2|_B\,\frac{dx}{\sqrt h}\right)
 =\frac{c^{3/2}}{5184},\qquad
 \operatorname{Res}_0\left(-(2\pi i)^{15} Q_3|_B\,\frac{dx}{\sqrt h}\right)
 =8c^{5/2}\left(-\frac{144754243}{30813563289600}
     +\frac{226271}{59439744000}\zeta(3)\right).
$$
Here $Q_g|_B$ denotes the axis coefficient before multiplying by the tangent wave, at unit background strength. The value at the double zero depends only on $h''(0)=2c$, fixed by this full subgroup. At a zero of order $m>2$ the differential polynomial vanishes to order at least $g(m-2)$; multiplying by the Gaussian form gives zero residue, also on the ramified double cover for odd $m$. At a simple zero the form itself is regular on that cover. At $h=(1-x^2)^2$ the choice $\sqrt h=1-x^2$ gives residues $+1/648,-1/648$ at $-1,+1$ in genus two. This contour measures the local bulk scalar on the moving regular brane; its comparison with the complete boundary Ward operation is open.

For the infinitesimal volume-preserving generator $k(x)\partial_x-k'(x)y\partial_y$, the background variation is $\delta h=k h'-2k'h$. Subtract transport of the scalar output along the critical line. Its reduced two-loop variation is
$$
 \frac{2h''}{20736}\big((k h'-2k'h)''-k h'''\big)
 =-\frac{h''(3k''h'+2k'''h)}{10368}.
$$
For $h=e^x,k=e^{ux}$ at the origin, the four original theta potentials in the established path order are
$$
 u^2\left(\frac{82u+131}{2177280},\frac{44u+67}{544320},
          \frac{2u+3}{48384},\frac{3u+4}{90720}\right).
$$
Their sum cancels the displayed variation. For arbitrary $h,k$, polarize $h$ in two frequencies and the shear as the first $v$ jet of $e^{ux+vy}(v\partial_x-u\partial_y)$. The three $y$ frequencies have joint degree five; the jet orders $2,2,1$ can be taken by centered differences. After extracting them the three $x$ frequencies also have joint degree five. The integer degree-five simplex at $u=1$ tests the entire identity, with signed points as holdouts. The bare flux-free two-loop Ward graphs thus measure the finite-coordinate issue in this local matrix element. Identifying its transfer to the intrinsic D3 differential still requires the relative diagrams. [Original graph projection by shape](critical_line_ward.py)

#### Mixing the massive coordinate at a double zero

At the crossing take $c\ne0$ and the volume-preserving shear $z\mapsto z+uxy$. The potential becomes $B=((z+uxy)^2+cx^2y^2)/2$. Along the $x$ branch the transverse Hessian in $(y,z)$ has determinant $cx^2$ for every $u$. With the branch of $\sqrt c$ fixed, the Gaussian form is $dx/(x\sqrt c)$. The residue of an elementary holomorphic output against this form is its value at the crossing divided by $\sqrt c$.

For the two tangent primary waves of the preceding calculation, each genus-$g$ normal axis carries $2g$ derivatives. If $k$ background marks are quartic, exactly $k$ are quadratic and $2g-2k$ are cubic. Their coefficient multiplies $u^{2g-2k}(c+u^2)^k$. The original seven-mark graph gives the three reduced coefficients $1/51840,-1/9720,1/5184$. Thus
$$
 (2\pi i)^{10}Q_2\big|_{B,0}
 =\frac{30c^2+44cu^2+17u^4}{155520}.
$$
For the infinitesimal generator $xy\partial_z$, the four original theta shapes sum to $u^3/7776-11u(c+u^2)/19440$. This cancels the derivative with respect to $u$. The direct transverse brane matrix calculation gives the ordered disk on the scalar $f$ as $-f$ on the critical line whenever $cx^2\ne0$.

On this shear family the two-loop term can be killed at $u^2/c=(-22\pm i\sqrt{26})/17$. The two roots are nonreal. Put $r=u^2/c$. Divide each of the two real coefficient polynomials of the reduced genus-three value divided by $c^3$, in the period columns $1,\zeta(3)$, by $17r^2+44r+30$. The determinant of the two linear remainders, with coefficients ordered by descending power of $r$, is $-1290669321259/4721004731437682688000000$. A real linear polynomial vanishing at either root must be identically zero. The determinant therefore excludes simultaneous vanishing of the two elementary terms **within this shear family**.

The same conclusion holds for holomorphic volume-preserving changes of the three normal coordinates fixing both straight critical branches as sets. Assign weights $1,1,2$ to $x,y,z$. The leading jets of such a change preserving the branches individually are $X=\alpha x$, $Y=\beta y$, $Z=\gamma z+u xy$ in weights $1,1,2$ respectively, with $\alpha\beta\gamma=1$. A branch exchange can first be undone by a linear determinant-one signed permutation. Only the weight-four jet of the potential contributes at the crossing: the $2g$ background slots each start at weight four and the graph differentiates to total weight $8g$. Removing the linear change reduces to the preceding shear with parameter $u/(\alpha\beta)$. Along either branch $X=f(x)$, the transverse Jacobian on the $x$ branch is $1/f'(x)$; the Gaussian density pulls back to $f'(x)dx/(\sqrt c f(x))$ with the transported square-root choice. [Two- and three-loop graph checks and two-loop Ward shapes](critical_line_massive_shear.py) · [Sheared D3 disk](../twisted-holography-N4/code/critical_line_brane.py)

For $z\mapsto z+a x^2+bxy+d y^2$ the two critical branches can bend. On the transformed $x$ branch $y=0$, $z=-a x^2$, the transverse Hessian determinant remains $cx^2$. Use $X=b^2/c$ and $Y=(b^2-4ad)/c$. The two-loop scalar in heat units is
$$
 \frac{c^2}{155520}\left(66X+17Y^2-22Y+30\right).
$$
At $c>0$ every real shear leaves this scalar positive: the last three terms are $17(Y-11/17)^2+389/17$. Over the complex numbers both elementary residues can vanish on the same shear. Substitute $X=-(17Y^2-22Y+30)/66$ into the three-loop graph polynomial. Its reduced value divided by $c^3$ is
$$
\begin{aligned}
 &\frac{6740370862184Y^3+7511864728413Y^2-8402913085068Y-10564526354410}
 {633834996867072000}\\
 &\quad-\zeta(3)\frac{1056331747Y^3+1139722653Y^2-1409555850Y-1712240306}
 {122267553408000}.
\end{aligned}
$$
The derivative is positive on the real line; the cubic has opposite signs at $-1/2$ and $-4/9$. These inequalities follow from the integral bounds on the tail after 64 terms of $\zeta(3)$. For its unique real root in that interval, take $a=1$, $b^2=cX$, $d=c(X-Y)/4$. This is a determinant-one holomorphic change of variables. Direct graph evaluation in a regular frame for the rank-one shear $(x+y)^2$ at $c=1$ gives $49/25920$ at genus two and $-13414400183/60526642176000+(511429/2918916000)\zeta(3)$ at genus three, agreeing with the general formula at $Y=0$.

In the simpler $uxy$ family the genus-three hexagon with a genus-two vertex gives $-1/14515200,1/2721600,-1/1451520$ multiplying $u^5,u^3(c+u^2),u(c+u^2)^2$. Replay the massive-shear calculation with `--hexagon`. The two elementary residues depend on the holomorphic field coordinates; the relative quantum operation on the D3 brane still needs its boundary terms. [General quadratic shear and interval certificate](critical_line_quadratic_shear.py) · [Rank-one holdout and Ward graphs](critical_line_massive_shear.py)

#### A bare Ward integral with five hard sites for every repeated-soft direction

Keep the quadratic background $B=(y^2+z^2)/2$ at unit strength, and take two additional scalar inputs
$$
 \tfrac12 y^2 e^{p x+L_1z^1+L_2z^2},\qquad
 \tfrac12 y^2 e^{q x+M_1z^1+M_2z^2}
$$
alongside the lifted scalar $e^{z^1}$ and the vector $e^{z^2}\partial_{z^1}$. Here $p,q$ are normal frequencies, $L,M$ are tangent momenta and $y,z$ are normal coordinates. The reduced form-degree-zero part of this four-input local binary variation, evaluated on the critical axis at zero tangent position, is
$$
 \frac{6p^2q^2(1+L_1+M_1)}{(g-2)!}
 C_g\big((qL_1-pM_1)e^1+(qL_2-pM_2)e^2-qe^3+pe^4\big)\,e^{(p+q)x}.
$$
Indeed the two hard normal momenta $(p,s),(q,t)$ now have tangent components $L,M$. Expressing the $y$ direction in the four hard covectors gives the four displayed components divided by $pt-qs$. The determinant power and $s^2t^2$ extraction are as before. The root divergence on the tangent vector gives $1+L_1+M_1$. Every generic four-component soft direction occurs: for prescribed third and fourth components choose $-q,p$ equal to them, then solve the two linear equations in the tangent momenta. The root factor can be kept nonzero on a dense open set, which determines the whole polynomial at every genus. The five simplex coordinates can all be distinct. Two such directions have been checked against the original fivefold seven- and nine-mark graph polynomials. [Graph replays](critical_line_binary.py)

Both additional inputs have exact closed lifts on the unfluxed or fluxed puncture. For $P=px+L_1z^1+L_2z^2$ with $L_1\ne0$, the vector $e^P(y\partial_y-\partial_{z^1}/L_1)/2$ has divergence zero and contracts $dB$ to $y^2e^P/2$. Its full linearized differential supplies the exact lift, with the same degree-zero scalar at any flux. If $L_1=0$ use a tangent antiderivative in that coordinate. The explicit five-direction rank test excludes quantum linearizations on this quadratic and the purely normal flux coefficients. The nonzero displayed degree-zero term cannot be a linearized boundary since it survives evaluation on the punctured critical axis. This supplies an exact-input test of the four-input elementary operation and fixes the cancellation required in its Ward equation; the two-input operation on a fixed perturbed background was measured separately above. [Exact source primitives with flux](../twisted-holography-N4/code/critical_line_brane.py)

Insert the first primitive into the Ward equation at form degree zero. The hard sites are the tangent vector, the primitive, the tangent scalar, the second mixed scalar and the balancing site. There are $g$ quadratic $z$ donors and $g-1$ quadratic $y$ donors. A positive-genus action insertion of genus $j$ on primary flags has $2j+3$ flags. At least four flags must be outside the donor plane: the three-dimensional quotient must be spanned by conserved projected charges. At most $2j-1$ donors fit. At equality a shear fixing the donor plane removes both donor-coordinate components of three independent hard flags. The fourth is fixed by conservation. Every one of the $4j$ donor-coordinate derivatives must now fall on the $2j-1$ quadratic donor slots, which can supply at most $4j-2$. The trace vanishes at the required Taylor degree even with arbitrary internal currents, since their projected components are independent of the donor-coordinate momenta. Lower projected rank kills each five-row determinant. A geometric loop graph of genus $\ell$ admits at most $2\ell$ ordinary donor jumps in that plane, by the rank of their cycle columns. The total $2g-1$ donors therefore exclude every positive-genus insertion from the positive-geometric-loop terms. A tree composition of two positive-genus actions also vanishes: each side needs at least three external hard flags in the quotient since its sewing momentum is fixed by conservation, whereas there are only five such flags. These counts use the primary flat bulk Ward expansion.

Among classical fusions only the primitive acting on a $y^2/2$ donor survives. In a hard-hard fusion at most four hard flags remain. Quotient rank three lets a shear remove the two plane components of three flags, with the fourth fixed by the donors. The fused trace is then independent of the two mixed hard $y$ momenta; the vector pairing supplies at most one power against the three required Taylor derivatives. A donor joined to the balancing site disappears from the fused trace and the vector pairing cannot supply its two derivatives. A $z$ donor pairs to zero with either vector. The $g-1$ choices of $y$ donor give precisely the Taylor divisor $1/(g-2)!$. Consequently the signed Ward sum of **bare cubic graphs alone** cancels the displayed elementary four-input coefficient at every genus. This gives a second integral for the four-loop mixed coefficient. The earlier six-hard-site projection has nonzero lower-loop terms.

For $p=2,q=-3,L=(1,-1),M=(2,3)$ the sole nonzero two-loop tree fusion gives $1/72$. The original $(1,5,5),(2,4,5),(3,3,5),(3,4,4)$ theta integrals give $-41/15120,-11/1890,-1/336,-1/420$. For the mixed measurement $p=1,q=-1,L=(1,0),M=(-1,0)$ the corresponding tree value is $1/10368$ and the four theta integrals cancel it in the same proportions. The original genus-three hexagon with a two-loop insertion was also integrated separately on both fixtures and gives zero. [All fusions, theta shapes and hexagon](critical_line_ward.py)

**Five-site reduction of the mixed integral.** In the latter measurement label the five hard sites $T,V,f,M,R$ in their preceding order. After forgetting the $g$ quadratic $z$ donors, the remaining $g-1$ quadratic $y$ donors give the leading soft squared cofactor on the fixed hard network. Use the four unit hard covectors and their balancing covector, with soft direction one on $V,M$ and zero on $T,f,R$. Integrate the squared cofactor over all soft positions and divide by $(g-1)!$ for the repeated donors. Multiply this by
$$
 \frac34\left(2\langle R_{ij}\rangle_{\rm cross}
             - R_{VM} -\langle R_{ij}\rangle_{T,f,R}\right).
$$
Here $R_{ij}$ is effective resistance; the first average is over six pairs joining the two groups, the last over the three distinct pairs in $T,f,R$. To derive the factor, ground the Green function $G$ at $R$. For hard $y$ slopes $s,t$ the squared cofactor divided by the canonical one is $(s+t)^2$: the four hard covectors have determinant of magnitude $|s+t|$ and the $g-1$ soft columns each use its inverse. The extra hard-slot voltage depends linearly on $s,t$. Its $s$ coefficient is $-G_{Vf}+G_{VM}$; its $t$ coefficient is $-G_{TV}-G_{Tf}+G_{TM}+G_{VM}$. Extracting $s t^2$ picks the first plus twice the second. Include the factor $1/2$ from the primitive. Use $2G_{ij}=R_{iR}+R_{jR}-R_{ij}$ and average the two permutations of $V,M$ and the six permutations of $T,f,R$. This gives exactly the displayed factor, with the original Kirchhoff measure and diagram automorphisms. At genus four there are fourteen remaining metric segments and three soft jumps. Summing all 364 segment triples and contracting their three-row Gram matrix agree, as do the twelve unsheared hard assignments after the jet projection, on three rational bipartite networks. The resistance factor takes both signs on these networks. [Direct five-site integrand](critical_line_five_sites.py) · [Soft placement identity](#a-ward-projection-for-the-mixed-genus-four-d3-pole)

On each connected bridgeless subgraph with $k$ core paths and $h$ independent cycles, restrict the three projected momentum coordinates to its cycles. Full projected rank requires $3h$ independent current rows there. On each path the differences of successive current columns are its cycle column tensored with the projected momentum of the corresponding hard site. Modulo these jumps there are $k$ columns, with a nontrivial linear relation given by the positive total path lengths (Kirchhoff's voltage law). Thus the jumps must have rank at least $3h+1-k$. Also the local source momenta and the currents entering at the exposed vertices must span three dimensions. Conservation fixes the sum of the entering currents, leaving at most the number of exposed vertices minus one independent new directions, or zero for the whole core. These bounds exclude individual colorings for all metrics.

At a positive rational metric, full and projected ranks $14,13$ imply that the complementary cycle current is nonzero. The cycle Gram adjugate is positive definite there. Testing those ranks and the voltage modulo $1000003$ therefore certifies a nonzero rational integrand. For the sum over colors we also contract the four complementary minors with the Gram adjugate and add the ten excited pairs at the *same* metric before testing. In every geometrically admitted placement this sum has a nonzero residue for at least one integer metric.

| Core | Hard shapes | Path multisets | Surviving color orbits | Identically zero color orbits within the admitted shapes | Surviving labeled hard orders |
|---|---:|---:|---:|---:|---:|
| Two-double | 11 | 69 | 47 | 16 | 6912 |
| Three-double | 0 | 0 | 0 | 0 | 0 |
| Single-double | 36 | 171 | 217 | 46 | 17820 |
| Prism | 80 | 786 | 675 | 16 | 92484 |
| Complete bipartite | 26 | 1134 | 177 | 6 | 132840 |

The eleven two-double shape integrals, including the core symmetry divisor, sum to $-22550149/6387888697344000000$ in the displayed resistance convention. Of the 47 contributing color orbits, one gives a positive period. The other 46 give negative periods. [Exact support and summed-color tests](operator/five_hard_support.py) · [Integrated two-double chambers](operator/five_hard_mixed_double.py)

#### An exact all-loop scalar response on one quadratic background

In the flat minimal fivefold theory set $B=((z^2)^2+(z^3)^2+(z^4)^2+(z^5)^2)/2$. Take the part of $Q_g\phi(0)$ linear in a primary scalar wave $e^{z^1}$ and a constant primary vector $e_1=\partial/\partial z^1$. In addition to these polarized inputs, fix the scalar fields in the ordinary action coordinates to
$$
A_0=B-aB^3/6,\qquad A_1=aB^2/2,\qquad A_2=-aB,\qquad A_3=a,\qquad A_4=0 .
$$
The same profile gives the complete matrix element at every genus, including the tree:
$$
(2\pi i)^{5g}Q_g\phi(0)=(-a)^g\frac{2(3g)!}{6^g g!(2g+2)!}.
$$
In the $a^m$ coefficient, keep the $m$ level-three scalars, the wave, vector and absorbed root. If the scalar levels are lowered by a total of $s$, the fivefold dimension gives $2g+2m-s$ additional primary copies of $B$. The lowered fields supply $s$ further copies. Each homogeneous quadratic requires exactly two normal derivatives at zero. The top Betti form has $g$ factors per normal direction, each quadratic in the corresponding momenta. It therefore supplies $8g$ normal derivatives. The absorbed root pairs tangentially with $e_1$. Thus $2(2g+2m)=8g$, or $m=g$. This also treats genus zero, where the cubic on the two tangent waves gives $+1$. At equality the correlated scalar components give the cotangent pullback above; the transverse Hessian has determinant one and the root pairing is $-1$.

To evaluate the remaining one-target vertex, use the top charge coefficient of the Buryak–Shadrin–Spitz–Zvonkine marked-point formula on charges $(1,-1,0,0,\ldots,0)$. The three primary marks have zero cotangent exponent; the other $g$ marks each have exponent three and zero charge. Write $X$ for the sum of their cotangent variables. In the permutation formula for $P_n$, setting a cotangent variable to zero forces that mark to be last after choosing another initial mark. Its last determinant and denominator give $(e^{bX/2}-e^{-bX/2})/b$ for charge $b$, with value $X$ at $b=0$. The statement extends from generic momenta by regularity. With all charges zero, $P_n=X^{n-2}$: for $n\ge2$ it has degree $n-2$, the specialization gives the identity on each coordinate hyperplane by induction from $P_2=1$, and their product has degree $n$. For $n=1$ the source formula gives $P_1=1/X$.

Scale the two charges by $u$. The numerator after the three primary specializations is $X^{g+1}(\sinh(uX/2)/(uX/2))^2$. In charge degree $2g$ and cotangent degree $3g$ only the leading $1/X$ from $1/(e^{X/2}-e^{-X/2})$ contributes. Since $(\sinh(y/2)/(y/2))^2=2(\cosh y-1)/y^2$, the coefficient on the $g$ labeled level-three marks is $2(3g)!/[6^g(2g+2)!]$. One repeated scalar amplitude divides it by $g!$. [Marked-point formula, Theorem 2](https://arxiv.org/abs/1211.5273) · [Scalar graph identification](../twisted-holography-N4/README.md#deletion-and-contraction-for-every-scalar-double-ramification-graph)

For a general tangent scalar $f(z^1)$ the same matrix element is the following local differential operator. Replacing $e^{z^1}$ by $e^{p z^1}$ scales every axial charge by $p$, supplying $p^{2g}$ in the period and $p$ from the root. At finite jet order this determines the operator:
$$
(2\pi i)^{5g} Q_g\phi(0)\big|_{f,e_1}
=\frac{2(-a)^g(3g)!}{6^g g!(2g+2)!}\,
\frac{\partial^{2g+1} f}{\partial(z^1)^{2g+1}}(0).
$$
For any chosen positive loop order $g$, the single polynomial $f=(z^1)^{2g+1}/(2g+1)!$ has exactly one nonzero response in this family when $a\ne0$.

For the unit exponential and string coupling $\lambda$, put $x=-\lambda^2 a/[6(2\pi i)^5]$. The branch $C(0)=1$ of $C=1+xC^3$ has coefficient $\binom{3g}{g}/(2g+1)$ at $x^g$. Dividing each coefficient by $g+1$ gives the required response. Using $x=(C-1)/C^3$ to integrate $C$ gives
$$
\sum_{\substack{g\ge0\\\mathrm{loop\ orders}}}\lambda^{2g}Q_g\phi(0)
=\frac1x\int_0^x C(s)\,ds
=\frac{C(3-C)}2 .
$$
The series converges for $|x|<4/27$. At $x=4/27$ its value is $9/8$. This is a bulk matrix element of the announced Li–Wang quantization in the stated field coordinates.

The same computation applies to constant hard scalars of levels $k=2,3,4$ with independent amplitudes $a_k$ and components $a_k(-B)^j/j!$ at level $k-j$. For multiplicities $m_k$, put $m=m_2+m_3+m_4$. Only $m_2+2m_3+3m_4=2g$ contributes. The coefficient of $a_2^{m_2}a_3^{m_3}a_4^{m_4}$ in $(2\pi i)^{5g}Q_g\phi(0)$ is
$$
\frac{(-1)^g\,2(2g+m)!}{(2g+2)!\,(2!)^{m_2}(3!)^{m_3}(4!)^{m_4}\,m_2!m_3!m_4!}.
$$
Equivalently use the classical scalar series $U(s)$ of the genus-one formula above with primary value $s$ and positive scalar levels $(0,a_2,a_3,a_4)$. It solves $U=s+a_2U^2/2+a_3U^3/6+a_4U^4/24$. Lagrange inversion gives the whole genus-$g$ tangent response as
$$
(2\pi i)^{5g} Q_g\phi(0)\big|_{f,e_1}
=\frac{(-1)^g}{g+1}\,[s^{2g+1}]U(s)\,
\partial_1^{2g+1}f(0).
$$
Here the coefficient in the classical inverse series is in the primary argument $s$ with all three descendant strengths fixed. Indeed formal inversion gives $[s^{2g+1}]U=(2g+1)^{-1}[y^{2g}](1-a_2y/2-a_3y^2/6-a_4y^3/24)^{-(2g+1)}$, with precisely the multinomial coefficients above.

In particular a sole level-four amplitude selects $g=3r$, with first correction $(2\pi i)^{15}Q_3\phi(0)=-a_4^2/576$ on the unit exponential. The scalar stable graphs through genus four, the general marked-point permutation formula through genus three, and the genus-one fivefold hexagon and tree Ward sum verify the respective coefficients. [Exact replays](quadratic_loop_selection.py)

The fivefold genus-two Ward equation gives an independent check on the level-three background. Polarize the transverse Hessian along the four unit axes and take two identical donor waves on each axis before dividing by $(2!)^4$. Include the $2!$ for the equal level-three marks. The zero-momentum vector is a primary in the scalar trace. The string equation at that slot gives the sum of the two possible level lowerings. The signed heat coefficients at zero through six donor fusions are
$$
\left(\frac{68081}{25920},-\frac{13483}{1944},\frac{14981}{2160},
-\frac{42011}{12960},\frac{10757}{15552},-\frac{211}{4320},0\right).
$$
At zero fusions the ordinary levels after string are $(2,3)$ on zero momenta. Deform the two marked momenta by $\pm s(e^2-e^3)$ and keep all other inputs fixed. The four five-row determinants give degree at most four in $s$, and exchanging the axes makes it even. At precursor levels $(2,2)$, exchange of those axes also equates the two raised vertices. Their Ward-row difference projected on $e_2$ thus gives $2s$ times either vertex. Only tree cuts containing exactly one marked site survive this difference. The three samples $s=1,2,3$ give $107263/972,\ 2238491/4860,\ 2476217/1620$, each independently reproduced by the full row system. Interpolation in $s^2$ with weights $3/2,-3/5,1/10$ gives $68081/1620$. The two string terms and the identical-field factors multiply this by $1/16$. The complete sum is $68081/25920-67361/25920=1/36$. [Fivefold Ward contributions with `--unfused 1 2 3 4 5 6`](quadratic_two_loop_ward_fusions.py) · [Projected and full Ward rows](genus_two_graph_ward.py)

#### Dependence on the choice of the elementary point observable

Keep the bulk action and the input fields fixed. At one loop add $c\lambda^2 A_3\det(\partial_i\partial_j A_0)/(2\pi i)^5$ to the point observable, with $i,j$ the five affine coordinates. On the level-three background the implicit classical scalar satisfies $u=B+f+a(u-B)^3/6$ to first order in the tangent wave $f$. Thus $u=B+f$ at that order. The constant vector gives $Q_0 A_0=\partial_1u$, and $Q_0A_3=0$ on these inputs. The change in the heat-normalized one-loop response is $c\,a\,\partial_1^3f(0)$. Taking $c=1/12$ cancels this particular coefficient.

At genus $g\ge2$ use $g$ copies of the adjugate of the Hessian $\partial_i\partial_j A_0$, contracting their $2g$ upper slots with the $2g$ derivatives of one further $A_0$. Multiply by $\lambda^{2g}A_3^g/(2\pi i)^{5g}$. This affine-volume-invariant point polynomial has $5g+1$ scalar factors, total level $3g$, and degree $2(5g+1)-6g-4g=2$, the degree of $\phi$. On the background its Hessian adjugate at zero is $e_1\otimes e_1$. If the variation acts on one adjugate, at least one unvaried adjugate contracts two tangent derivatives into the purely transverse background. Only variation of the differentiated $A_0$ survives, giving $a^g\partial_1^{2g+1}f(0)$. The coefficient can therefore be shifted at each specified order while keeping all previous orders fixed. With all three scalar strengths, replacing $A_3^g$ by $A_2^{m_2}A_3^{m_3}A_4^{m_4}$ gives the same test for each of the allowed monomials above.

The highest term acting on $g+1$ distinct input letters has a separate invariance under multiplicative field changes that are the identity at $\lambda=0$. Its $(g+1)$-fold commutator with multiplication by $\phi$ annihilates all lower-loop differential operators. [Finite-change and jet checks in sheared coordinates](quadratic_mark_change.py) · [Scalar differential of order $g+1$](#an-exact-matrix-element-of-the-scalar-differential-at-every-loop-order) · [Classical scalar variation](operator/README.md#closed-point-classes-from-scalar-source-insertions)
