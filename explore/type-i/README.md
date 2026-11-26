# Type I Kodaira–Spencer theory and holomorphic Chern–Simons theory on $\mathbb C^5$

### The hexagon

Use the repository's cyclic edge order $e=(0,1),\ldots,(5,0)$, vertex $5$ at the origin, with incoming holomorphic covectors $p_0,\ldots,p_5$ of zero sum. Set $U=\sum_e t_e$ for the positive edge times. The solution of current conservation has $y_e-y_{e-1}=p_e$ and $\sum_e t_e y_e=0$. For the last heat edge, take the independent differences $w_e=z_e-z_{e+1}$, $0\le e<5$. The last difference has $d\bar w_5=-\sum_{e<5}d\bar w_e$. In the exterior product of five propagator numerators and one heat form the coefficient is
$$
-\det(y_0,\ldots,y_4)=-\frac{t_5}{U}\det(p_0,\ldots,p_4).
$$
The first equality follows by selecting the missing differential in each contracted five-form. The second follows by subtracting successive rows and using the coefficient $t_5/U$ of $p_0$ in $y_0$. Moving the heat edge replaces $t_5$ by its marked time. With $t_5=\epsilon$, the five remaining integrations give
$$
\int_{[\epsilon,\infty)^5}\frac{\epsilon\, dt_0\cdots dt_4}{(\epsilon+t_0+\cdots+t_4)^6}
=\frac1{6\cdot5!}.
$$
Thus the repository's reduced kernel is $-\det(p_0,\ldots,p_4)/120$. The direct Gaussian calculation uses the normalized Dolbeault delta form, $\int d^5z\wedge K_t=1$. In this orientation $d^5z\wedge d^5\bar z=-32i\,d^{10}x$. The scalar Gaussian integrations leave $(-32i)^{-1}(4\pi)^{-5} U^{-5}$; the five propagator contractions supply $4^5$. Consequently the ordered heat wheel is $+\det(p_0,\ldots,p_4)/(120(2\pi i)^5)$. Here the five vertex volume forms stand to the left of the edge forms. [Exact exterior-determinant checks on all six facets](hexagon.py)

For one external adjoint field $X_v$ at each vertex, order the three odd letters as external, incoming, outgoing. Reordering from all six external slots followed by the six paired edge slots gives one odd rotation of twelve slots. The resulting color weight is $-\operatorname{Tr}_{\rm ad}(\operatorname{ad}X_0\cdots\operatorname{ad}X_5)$ with reversal of the trace giving the same symmetrized polynomial. Alternation over the external slots cancels the alternating sign of their momentum determinant. The unoriented hexagon has automorphism order $12$. For identical $X_v=X$, the fully polarized coefficient with the stated momenta is therefore
$$
-\frac{6!}{12\cdot5!(2\pi i)^5}\det(p_0,\ldots,p_4)\,\operatorname{Tr}_{\rm ad}((\operatorname{ad}X)^6)
=-\frac{\det(p_0,\ldots,p_4)}{2(2\pi i)^5}\,\operatorname{Tr}_{\rm ad}((\operatorname{ad}X)^6).
$$
Equivalently the invariant characteristic polynomial in Gwilliam–Williams' holomorphic charge-anomaly convention is $-\operatorname{Tr}_{\rm ad}((\operatorname{ad}X)^6)/(2\cdot6!(2\pi i)^5)$. Their directed wheel for an even field has polynomial $\operatorname{Tr}_V(X^6)/(6!(2\pi i)^5)$. The self-paired adjoint wheel has half the oriented multiplicity, and its Dolbeault coefficient space is $\mathfrak g[1]$, with negative supertrace. For a master field $A$, extracting the ordered plane-wave coefficient of $\operatorname{Tr}_{\rm ad}(A(\partial A)^5)$ gives $6!\det(p_0,\ldots,p_4)\operatorname{Tr}_{\rm ad}((\operatorname{ad}X)^6)$: moving $d^5z$ past the six odd coefficients and replacing $\det(p_1,\ldots,p_5)$ by $-\det(p_0,\ldots,p_4)$ gives two minus signs. Costello–Li give this local tensor in (83) and leave the multiplicative constant in the cancellation unspecified. The right side of their (85) has one fewer $\partial A$ than (83) and the degree-six polarization require. [Costello–Li, §III.D](https://arxiv.org/pdf/1905.09269#page=33) · [Gwilliam–Williams, Proposition 3.9 and Appendix B](https://arxiv.org/pdf/1810.06534#page=43) · [Costello–Li, Lemma 7.2.8](https://arxiv.org/pdf/1505.06703#page=82)

### The closed exchange

For the defining orthogonal representation $V$ of dimension $N$, $\mathfrak{so}(V)=\Lambda^2 V$. The eigenvalues of $X\in\mathfrak{so}(V)$ come in opposite pairs. The identity
$$
\operatorname{Tr}_{\Lambda^2V}(e^{tX})=\tfrac12\big[(\operatorname{tr}_V e^{tX})^2-\operatorname{tr}_V e^{2tX}\big]
$$
gives at order six
$$
\operatorname{Tr}_{\rm ad}((\operatorname{ad}X)^6)
=(N-32)\operatorname{tr}_V(X^6)+15\operatorname{tr}_V(X^2)\operatorname{tr}_V(X^4).
$$
The closed BV kernel pairs the $\mathrm{PV}^{1,*}$ and $\mathrm{PV}^{3,*}$ disk vertices by $\partial_\Omega$. Its unpolarized local tensor is $\operatorname{tr}(A\partial A)\operatorname{tr}((\partial A)^4)$. In characteristic-polynomial normalization its required coefficient is $+1/(96(2\pi i)^5)$. The closed kinetic pairing and the two disk couplings enter only through their product in this exchange. After that exchange the residual polynomial is $-(N-32)\operatorname{tr}_V(X^6)/(1440(2\pi i)^5)$. It vanishes at $N=32$. The rank selection concerns these disk couplings and the orthogonal Chan–Paton system. [Costello–Li, (82)–(85)](https://arxiv.org/pdf/1905.09269#page=32)

### Two loops, eight external gauge fields

The four pure holomorphic Chern–Simons theta graphs have three path lengths $(1,5,5),(2,4,5),(3,3,5),(3,4,4)$. Both internal endpoints have zero external momentum. For each ordering of the eight current inputs along the interiors, the integrated reduced kernels are the universal determinant $Q(\bar L)$ of the repository times the signed constants $1/17280,-1/17280,1/17280,-1/25920$, respectively. The weight of the $j$th internal momentum on a path of length $p$ in $\bar L$ is $(p-j)/p$. Each of the two cubic endpoints has color $\operatorname{tr}([A,B]C)=2\operatorname{tr}(ABC)$ on skew matrices. Transport along a path composes the eight adjoint actions in path order. With the vector-trace pairing, a summed generator joins two index lines with the direct and transposed contractions of weights $1/2,-1/2$. [Universal theta integral](../twoloop-universal/twoloop.py) · [Orthogonal contractions and integration](open_theta.py)

The table holds for arbitrary conserved momenta and all $N$. Each row gives the coefficient of a product of traces of consecutive external matrices of the indicated lengths, multiplied by $N$ to the indicated rank power. Other labelings follow by joint signed permutation. Each pair $ij$ in a term represents the five-momentum determinant on the complement of $i,j$ in $0,\ldots,6$, times the permutation sign of that complement followed by $i,j$. The eighth momentum is eliminated by conservation. Straightening products uses the three-term Plücker identity. Each graph includes its automorphism factor and the full unnormalized exterior polarization. These are the reduced coefficients before the common physical two-loop heat factor. [Exact table](open_theta.json)

For a sparse $\mathfrak{so}(32)$ input, write $T_{ij}=E_{ij}-E_{ji}$ using matrix units. Choose $(T_{01},T_{01},T_{23},T_{23},T_{45},T_{56},T_{67},T_{74})$. The first five momenta form the coordinate basis, followed by $(1,-1,2,0,3)$, $(0,2,-2,1,1)$ and the negative sum. Only the $(2,2,4)$ product of traces survives. In path-length order the reduced contributions are
$$
-\frac7{270},\qquad-\frac7{270},\qquad-\frac7{810},\qquad-\frac7{1620},
\qquad\text{sum}=-\frac7{108}.
$$
In the same ordered heat convention as the hexagon the two-loop multiplier is $(2\pi i)^{-10}$. The eleven normalized heat densities and nine Gaussian position integrations leave $(-32i)^{-2}(4\pi)^{-10}U^{-5}$; the ten propagator contractions supply $4^{10}$. The sign from gathering the vertex volume forms is $(-1)^{25\binom92}=+1$. For a marked edge $j$, choose the single missing differential in each of the ten propagator five-forms. Complementary minors of the grounded incidence matrix tensored with the five-dimensional identity convert the surviving $45\times45$ determinant to the cycle-matrix minor on the ten missing coordinates. Moving them to the end contributes $(-1)^{220+5j}$ after combining the interior-contraction signs. Exchanging cycle and coordinate order contributes $(-1)^{\binom22\binom52}=+1$. This gives exactly the project's $(-1)^j$ ordered facet. Thus the full ordered heat-kernel value on this probe is $-7/(108(2\pi i)^{10})$. The replay checks all $44$ marked determinants independently of the analytic argument.

Direct adjoint-basis contractions at $N=8,9,10$ reproduce all four numbers, summing all $8!$ assignments. The trace calculation is polynomial of degree at most two in $N$, so these checks fix the result also at $N=32$. [Literal contraction](check_open_theta.py)

A gauge-only seven-field primitive at this weight has ten derivatives. With six independent momenta, its kinematic factor is in the span of the 21 products of two five-minors. Orthogonal invariance allows the four trace types $(7),(2,5),(3,4),(2,2,3)$; the $32$-index volume tensor cannot enter. Gram contractions of the gauge indices give the listed traces, and $SL_5$ invariance gives the five-minors. The differential of these 84 candidates has rank five. Its entire $(2,6)$-trace component vanishes. The pure-open graph sum has coefficient $(N+13)/43200$ on the $(01)(23)$ minor product in that channel, equal to $1/960$ at $N=32$. [Exact differential and ranks](open_ce.py)

Indeed the pure-open sum has a nonzero Chevalley–Eilenberg differential. Take nine colors in two disjoint directed cycles on defining-representation indices of lengths two and seven; each directed edge is $T_{ij}$. The first eight momenta are the five coordinate covectors followed by three repeats of the first, second and third covectors. The ninth balances their sum. Only the seven adjacent pairs in the long cycle have nonzero brackets, with $[T_{ij},T_{jk}]=T_{ik}$ for distinct indices. Their signed sum at $N=32$ is $1/240$ in reduced units. The full coupled two-loop calculation must include the closed-sector and lower-order insertions. Its explicit counterterm remains open. [Nine-input replay](check_open_theta.py)

For classical disk and sphere vertices with only gauge fields external, the string loop order is $b_1+E_{\rm closed}-V_{\rm sphere}$, where $b_1$ is the ordinary graph loop number. This follows by gluing disks and spheres along open intervals and closed circles in the Euler characteristic count of Costello–Li (73). At string loop order two the one-cycle graph with one closed edge joins the $\mathrm{PV}^{1,*}$ and $\mathrm{PV}^{3,*}$ disks. Its other four vertices are open cubics. Deleting the closed line gives one open chain, so its orthogonal color polynomial is single trace.

For this mixed hexagon take the six vertex momenta as sums of their incident external momenta. In cyclic order $y_e$ equals the cumulative momentum through vertex $e$ minus the Schwinger-weighted mean of the six cumulative momenta. The divergence in the closed BV kernel supplies $y_{\rm closed}$, and the four disk derivatives supply an external momentum or a signed incident open current. Their five-determinant is affine in the common mean: every term with that mean in two columns vanishes. The projective simplex average therefore replaces it by the arithmetic mean of the six cumulative momenta. Multiplication by the reduced hexagon determinant gives the integrated geometric factor of each derivative assignment with no parameter integrals left. This single-trace contribution cannot alter the displayed $(2,6)$ coefficient. [Mixed-edge integral and independent simplex-corner check](mixed_hexagon.py) · [String loop counting](https://arxiv.org/pdf/1905.09269#page=26)

### The eight-gauge-input two-loop Ward equation

Keep the closed edge directed from the four-open disk at the root to the two-open disk. Put $p_0$ on the two-open disk, $p_1,\ldots,p_4$ on the four successive open cubics, and $p_5,p_6,p_7$ on the four-open disk. For fixed cyclic trace of the three latter labels, define the mixed tensor by unnormalized alternation of
$$
\det(p_0,\ldots,p_4)\det(p_0,p_5,p_6,p_7,3p_1+2p_2+p_3)\,
\operatorname{tr}\!\left([X_4,[X_3,[X_2,[X_1,X_0]]]]X_5X_6X_7\right).
$$
Its coefficient in reduced two-loop units is $-1/8640$. To see the relative sign, put $q=y_5$ on the closed edge. The two-open disk gives $2p_0$ modulo $q$. At the four-open disk the four momenta, including the internal open current, sum to $q$. Alternating the three derivatives gives $-24p_5\wedge p_6\wedge p_7$ modulo $q$, times the symmetrized trace. Under the exterior alternation the six orders of these three labels contribute equally. The mean $q$ is $-(3p_1+2p_2+p_3)/6$ modulo $p_0,p_5,p_6,p_7$. The two-disk tensor thus gives $+8$ times the second determinant. In the tree the same polarized disk tensor gives $-48\det(p_0,\ldots,p_4)$, and the heat edge pulls $d^5\bar z$ back with sign $-1$ at the integrated endpoint. This calibrates the disk product by the previously determined $+1/96$. Moving the five open-edge slot pairs into their vertices gives sign $-1$ in the mixed wheel. The remaining factors are the hexagon volume $1/120$ and $1/3!$ for three external fields on one disk. Their product is $-8/(96\cdot120\cdot3!)$. [Signed tree and wheel checks](mixed_hexagon.py)

After adding this mixed contribution at $N=32$, the single-trace sum is the differential of three seven-field cochains. For a seven-field cochain in the following table, take the indicated trace lengths with consecutive slots, multiply the five-minors on the two indicated index sets, then alternate over all seven inputs without division. The table gives **preimages** of the graph sum; the pure-gauge counterterm has their negatives in units of $(2\pi i)^{-10}$.

| Seven-input trace | First minor | Second minor | Coefficient at $N=32$ |
|---|---|---|---:|
| $7$ | $01234$ | $01234$ | $23/46080$ |
| $7$ | $01234$ | $01235$ | $1/23040$ |
| $7$ | $01234$ | $01245$ | $-11/11520$ |
| $3,4$ | $01234$ | $01234$ | $-19/1244160$ |
| $2,2,3$ | $01245$ | $02345$ | $1/41472$ |

The last two rows cover the $(4,4)$ and $(2,2,4)$ traces, respectively. If the single-trace mixed coefficient is initially left free, exact elimination against the full seven-field gauge image gives residual polynomials with monic greatest common divisor $N-32$. At that rank the unique mixed subtraction is $1/8640$, agreeing with the signed graph. [Single-trace calculation](single_trace_descent.py) · [Coefficients](single_trace_preimages.json)

The remaining $(2,6)$ trace passes through the gauge variation of the closed $\mathrm{PV}^{3,*}$ field, whose contraction with $\Omega$ is a two-form. Fix the pair in the length-two trace and put $q=p_0+p_1$, $p_7=-q-p_2-\cdots-p_6$. Write the kinematic coefficient as $L(p_0;q)$ with $p_2,\ldots,p_6$ held fixed. At fixed $q$ it has degree at most two in $p_0$ because it is a product of two determinants. Antisymmetry under $p_0\leftrightarrow q-p_0$ removes the quadratic part, and the soft zero removes the constant. In particular $L(q;q)=0$. For each component $L_k$ homogeneous of degree $k$ in $q$, evaluation on a decomposable two-form is
$$
H_k(a\wedge b;q)=
\frac{\partial_{q,b}L_k(a;q)-\partial_{q,a}L_k(b;q)}{k+1},
\qquad
H_k(p_0\wedge p_1;q)=L_k(p_0;q).
$$
Here $\partial_{q,b}$ differentiates $q$ along $b$ with the first argument of $L_k$ fixed. The second identity follows by differentiating $L_k(q;q)=0$ and using Euler's identity. The saved representation uses labels $7=q,8=a,9=b$. It contains 60 products of five-minors for each of the powers $N^0,N^1$. Substituting $a=p_0$, $b=p_1$ and straightening restores every stored $(2,6)$ coefficient identically. Averaging the six open slots over cyclic rotations and reversal with their exterior signs makes the lift descend to the cyclic trace. [All multiple-trace preimages](multitrace_preimages.json) · [Exact substitution and CE replay](closed_primitive.py)

To restore the closed-field scale, use the coefficient of $\operatorname{tr}((\partial A)^2)$ in the actual variation of the closed two-form. Denote this disk coefficient by $c$; Costello–Li fix its product with the four-open disk coefficient. On two ordered odd inputs the variation is $-2c\,(p_0\wedge p_1)\operatorname{tr}(X_0X_1)$. In the same order the closed-linear counterterm has coefficient $+H/(2c)$ in reduced two-loop units. Together with the five gauge terms this cancels the eight-gauge projection of the theta and mixed-hexagon sum. Terms in the coupled Ward equation with external closed fields still require a separate calculation.

For an independent single-trace test use the same eight momenta as above and $X_i=T_{i,i+1\pmod8}$. Direct adjoint sums for the four theta graphs, direct nested matrix commutators for all $8!$ mixed assignments, and the Chevalley–Eilenberg differential evaluated on the three seven-field cochains give respectively
$$
-\frac{299}{17280}-\frac{-364}{8640}=\frac{143}{5760}.
$$
The direct theta replay interpolates the rank-two color polynomial from $N=8,9,10$ to $N=32$. [Independent eight-cycle test](check_single_trace.py)

### Two external primary closed fields

Let $U_1$ be the cyclic first disk operation. On a divergence-free rank-$k$ polyvector $\mu$ it acts on ordered plane waves by $\mu(p_1,\ldots,p_k)$ times the ordered matrix product, with the first-disk coefficient of (16) understood. The primary ranks allowed by (76) are $1,3,5$. In dimension five a divergence-free top polyvector is constant in the holomorphic directions. For two such external closed fields, the projected diagrams at string loop order two are trees. Open sewing gives the Gerstenhaber bracket of the two first disks and the Hochschild differential of the second disk. Closed sewing to the cubic sphere gives the Schouten bracket of the inputs. [Costello–Li, (16), (31)–(38), (73), (76)](https://arxiv.org/abs/1905.09269)

The replay solves for a cyclic cochain whose ordinary Hochschild differential is
$$
[U_1(\mu),U_1(\nu)]_{G}-U_1([\mu,\nu]_{S}).
$$
The stored coefficients multiply the two first-disk coefficients; the sign of the term added to cancel the displayed defect is minus. In the ordered $k+l-2$ open arguments, each column evaluates $\mu$ and $\nu$ on two increasing subsets of the momentum labels. A derivative on the opposite closed field is included by inserting its momentum. Every open argument is differentiated, and insertion of a polyvector's own momentum vanishes by transversality. Expanding the Hochschild differential on plane waves merges two adjacent momenta. Rotation uses the spectator momentum equal to minus the sum of all remaining open and closed momenta. Reversal includes the transpose of each skew matrix and the permutation of the suspended open slots. The code compares every coefficient of these three polynomial identities before imposing dimension-five determinant relations.

| Closed ranks | Closed momenta | Columns | Rank | Nonzero coefficients in representative |
|---|---|---:|---:|---:|
| $(1,1)$ | both arbitrary | 1 | 1 | 0 |
| $(1,3)$ | both arbitrary | 3 | 3 | 0 |
| $(1,5)$ | vector arbitrary | 4 | 4 | 0 |
| $(3,3)$ | both zero | 12 | 12 | 12 |
| $(3,3)$ | both arbitrary | 42 | 42 | 40 |
| $(3,5)$ | both zero | 60 | 59 | 48 |
| $(3,5)$ | trivector arbitrary | 120 | 118 | 100 |

For constant trivectors, order the three-slot subsets of $0,1,2,3$ as $012,013,023,123$. The coefficient of $\mu(p_I)\nu(p_J)$ is the $(I,J)$ entry of
$$
\frac12\begin{pmatrix}0&1&2&2\\-1&0&1&2\\-2&-1&0&1\\-2&-2&-1&0\end{pmatrix}.
$$
For two rank-five inputs both coefficients are constant multiples of the volume polyvector, and their antisymmetric Gerstenhaber bracket vanishes. For $(3,3)$ and $(3,5)$ the sphere term vanishes in dimension five because $\partial_\Omega(\mu\wedge\nu)=0$ and the two divergences vanish separately. Thus the cyclic primitive cancels the two-first-disk defect for every divergence-free primary pair. The rank-$(3,5)$ rows use a choice of representative in the solution space. [Exact symbol elimination and replay](primary_disk.py) · [Cyclic formality, Proposition 27](https://arxiv.org/abs/0806.4095)

### One external primary closed field

For a graph with $m$ external closed fields, $b_1$ ordinary loops, $E_{\rm closed}$ closed edges and $V_{\rm sphere}$ spheres, the string loop order is $b_1+E_{\rm closed}-V_{\rm sphere}+m$. With a disk present, each sphere cluster uses at least as many closed edges as it has spheres, so $b_1+m\le2$. The local heat insertion on $\mathbb C^5$ has $5(V-1)=4(E-1)+5$ and $E=V+b_1-1$, giving $V=4b_1+2$ and $E=5b_1+1$. A purely closed graph has $E=E_{\rm closed}$ and $V=V_{\rm sphere}$, hence $2b_1-1+m=2$; the cases $(b_1,m)=(0,3),(1,1)$ violate $3V\le2E+m$. For one external closed field the solutions with disks are an all-open-edge hexagon with one first disk and a tree sewing a second disk to a first disk along a closed edge. For two external closed fields the solutions are precisely the open and closed trees above. These counts concern classical vertices in the local heat projection. [Sewing inventory](external_closed.py)

On the hexagon, put the five gauge cubics in successive slots $0,\ldots,4$ and the rank-$k$ disk at the sixth vertex. Integer Schwinger weights $5,4,3,2,1$ multiply the five incoming covectors; the arithmetic mean divides by six. After factoring out this $1/6$, the common ordered heat coefficient, and the polarized first-disk coefficient, the tensor evaluated in the replay is the full signed gauge permutation of
$$
\det(p_0,\ldots,p_4)\,
\mu(p_5,\ldots,p_{k+3},5p_0+4p_1+3p_2+2p_3+p_4)
$$
times the symmetrized defining trace of the two gauge-line endpoints and the $k-1$ other fields on the disk. The two endpoints are sewn with five ordered adjoint actions and the inverse defining trace form. For $k=1$, cyclic rotation of the adjoint trace averages each integer weight to $3$. The vector then evaluates on minus three times its own momentum, giving zero.

For $k=3$, choose the first seven open momenta of the eight-input example above and the closed momentum $q=(-2,-2,-1,-2,-5)$. Take $\mu$ to be the exterior product of $(-5,0,0,0,2)$, $(0,-5,0,0,2)$ and $(0,0,-5,0,1)$; each is annihilated by $q$. The coefficients of the consecutive trace words of lengths $7$, $(2,5)$ and $(3,4)$ in this tensor are respectively
$$
-23300+\frac{1000N}{3},\qquad 11000,\qquad 0.
$$
The replay evaluates all $7!$ assignments independently using the skew-matrix basis at $N=7,8$. With one external trivector, the tree formed from the $(3,3)$ disk has trace lengths $(2,5)$; the $(1,3)$ second disk vanishes in the representative above. The single trace of the hexagon at $N=32$ therefore survives this classical tree. Its cancellation requires the remaining one-closed-field Ward terms, including lower-order insertions or a primitive. The descent for $u^l\mathrm{PV}^{k,*}$ with $l>0$, $k+l$ odd also enters the full coupled complex. [Wheel projection and literal adjoint calculation](external_closed.py) · [Descendant disks, (19), (76)](https://arxiv.org/abs/1905.09269)

At $N=32$ the previously determined $(2,6)$ lift also fixes the single-trace variation with one external trivector. Let $W$ denote the integer-weight hexagon tensor displayed above, and $H$ the stored $(2,6)$ primitive averaged over the twelve cyclic and reflected orders of its six gauge inputs. With the closed polarization held fixed in the ordinary gauge Chevalley–Eilenberg differential,
$$
(d_{\rm CE}H)_{(7)}=-\frac{1}{5760}\,W_{(7)}.
$$
For $q\ne0$ write the closed two-form as $q\wedge\lambda$, so $\mu(p,r,s)=\det(q,\lambda,p,r,s)$. Use $p_0,\ldots,p_4$ as a basis, leave $p_5,p_6,\lambda$ arbitrary and impose $q=-p_0-\cdots-p_6$. All 210 nonzero coefficients in these variables agree. Polynomiality extends the identity across linearly dependent momenta. At $q=0$, hold a decomposable $\mu$ fixed while $q$ tends to zero in its annihilator; linearity gives the result for a general trivector by expansion into decomposable tensors. Independently, the gauge differential on 120 six-input consecutive-trace candidates has rank four after the dimension-five identities. A preimage of $W_{(7)}$ is the unnormalized alternation of $\det(p_0,\ldots,p_4)\operatorname{tr}(X_0\cdots X_5)$ times
$$
-\tfrac{14}{3}\mu(p_0,p_1,p_2)
-\tfrac{23}{3}\mu(p_0,p_1,p_3)
-\tfrac{23}{3}\mu(p_0,p_1,p_4)
-3\mu(p_0,p_2,p_3).
$$
The disk scale times the absolute ordered hexagon factor on $W$ is $24/(2\cdot2!\cdot6\cdot120)=1/120$ times the four-gauge disk coefficient. The factors are the polarized disk, the two orientations, the two open inputs at the disk, the current mean and the ordered hexagon integral. Its product with $c$ is $1/(96\cdot120)=1/11520$ in reduced units, equal to the magnitude supplied by $H/(2c)$. [Complete transverse polynomial identity](one_primary_descent.py)

In the double-trace channel put $p_0,p_1$ on the first disk. Its induced internal trivector acts on three covectors by $\det(p_0,p_1,\cdot,\cdot,\cdot)$ and has momentum $p_0+p_1$. Insert this as the second input of the unique $(3,3)$ disk cochain above, with four successive arguments $p_2,\ldots,p_5$ and spectator $p_6$. For its ordered symbol $U_2$ the complete identities are
$$
W_{(2,5)}=80\,U_2,\qquad W_{(3,4)}=0.
$$
The first identity agrees on all 142 nonzero coefficients. Both are independent of $N$. On restoring disk magnitudes, the first disk-induced vector-trace pair contributes $2c$; each first operation in the $(3,3)$ primitive contributes four times the coefficient of the four-gauge disk. The cyclic action divides by five and the five-slot skew trace has ten signed cyclic/reflected assignments. The resulting magnitude is $2c\cdot(10/5)\cdot 4^2=64c$ times the square of the four-gauge disk coefficient. This equals $80/120$ times that disk coefficient in the one-loop heat units by the calibrated disk product. The single- and double-trace identities are statements about the ordered tensors; the final common sign in the open–closed BV convention and descendant projections require separate checks. [Exact double-trace identity](one_primary_descent.py)

### One external top polyvector

A holomorphic divergence-free rank-five polyvector has zero momentum. Factor out its multiple of the coordinate volume polyvector. The integer-weight hexagon tensor above now has nine gauge inputs with zero momentum sum, and its two exterior factors are five-minors. For any $N$, the $(3,6)$ and $(4,5)$ trace components vanish. The single trace of length nine is an ordinary gauge Chevalley–Eilenberg differential for every $N$. On eight gauge inputs use products of two five-minors selected among seven independent momenta and an unnormalized alternating single trace. The 231 columns of their differential have rank ten. Both the constant and linear powers of $N$ in the nine-input wheel lie in this image. The replay computes rational preimages.

For the $(2,7)$ trace, insert the internal trivector dual to $p_0\wedge p_1$ as the **first** input of the rank-$(3,5)$ second disk, followed by the constant external top polyvector. Use $p_2,\ldots,p_7$ as its six ordered open inputs and $p_8$ as spectator. The two-dimensional freedom in the second disk after the primary tree Ward equation maps injectively into the $(2,7)$ trace. The required representative satisfies
$$
W_{(2,7)}=-48\,U_2.
$$
In the lexicographic column order of the primary disk calculation the two free entries are $81/280$ and $-11/30$. Replacing $-48$ by $+48$ leaves 194 nonzero polynomial coefficients after elimination. The absolute hexagon coefficient here is $6!/(2\cdot4!\cdot6\cdot120)=1/48$ times the six-gauge disk coefficient. The tree product has the same magnitude after using the previously calibrated disk couplings. The replay verifies the complete identity on 254 nonzero monomial coefficients in the entries of the three unconstrained covectors $p_5,p_6,p_7$, with the first five chosen as a basis. [Exact primary top-form descent and disk selection](top_primary_descent.py)

### Arbitrary first-disk descendants

For $u^l\mu$ of polyvector rank $k$ let $m=k+2l$ be the number of ordered open arguments after the spectator. Write $q$ for the incoming closed momentum. The first disk is a finite sum of evaluations of $\mu$ on $k$ covectors. For an ordered subset $I$ of $k$ open arguments the coefficient of $\mu(p_I)$ is $l!/m!$ when every unselected consecutive block in the linear order $0,\ldots,m-1$ has even length, and zero otherwise. For a subset $J$ of size $k-1$ the coefficient of $\mu(q,p_J)$ is
$$
\frac{l!}{(m+1)!}
\sum_{\substack{j\notin J\\J\cup\{j\}\ {\rm leaves\ even\ blocks}}}
(-1)^{\#\{i\in J:i<j\}} .
$$
For $k=0$ only the first formula occurs. In the simplex integral of Costello–Li the curvature is $\omega=\sum_i d\theta_i\wedge d\theta_{i+1}$. The complement of the derivative slots has a perfect matching by neighboring pairs exactly under the even-block condition, and that matching is unique. Hence $\omega^l$ contributes $l!$. The divergence term inserts the one-form $\eta$; its coefficient at position $j$ is the gap to the next marked point. Its integral over the ordered simplex is $1/(m+1)!$. Moving the associated momentum $q$ to the front gives the displayed sign. In the ordered Hochschild differential the resulting operator obeys $bU_1(u^l\mu)=U_1(u^{l+1}\iota_q\mu)$ by the simplex boundary identity. The replay checks the differential and cyclic projection independently for $0\le k\le5$, $0\le l\le3$. Type I retains $k+l$ odd. These first-disk coefficients include the divergence of a general incoming polyvector. [Coefficient generator and exact replay](descendant_disk.py) · [Costello–Li, Definition 3.1.3 and Lemma 3.2.1](https://arxiv.org/abs/1505.06703)

For an external $u\phi$ of scalar momentum $q$ and a divergence-free $u^l\mu$ of rank $k$, the second disk Ward equation in Kodaira–Spencer coordinates has an all-$l$ solution. In the ordinary ordered Hochschild convention above, its defect has the exact preimage
$$
(-1)^k\frac{k+2l-1}{2}\,U_1(u^l\phi\mu).
$$
Here the defect is the Gerstenhaber bracket of the first disks minus the first disk on the closed bracket. For $m=k+2l$, the first scalar disk is multiplication by $\phi/2$. Terms with undifferentiated $\phi$ in its bracket with $U_1(u^l\mu)$ vanish by the linear Ward equation and $\iota_r\mu=0$ for the momentum $r$ of $\mu$. In each remaining term a derivative hits $\phi$. Each allowed even complement at arity $m+1$ has $l+1$ matched pairs whose collapse gives the same term. Relative to the next descendant coefficient this gives
$$
(-1)^k\,(l+1)\frac{l!}{2m!}\frac{(m+1)!}{(l+1)!}
=(-1)^k\frac{m+1}{2}.
$$
The projected closed bracket contributes $(-1)^k$ on the same insertion of $q$; applying the preceding first-disk identity to $\phi\mu$ gives the preimage. For $k=0$ the defect vanishes. The BCOV coordinates of (38) also contain the quadratic shift from $u^{-1}\Sigma^2/2$ when a descendant is present. [Scalar and descendant replay](scalar_descendant.py)

With one external closed field the polynomial linear differential also reduces the descendant tower. Eliminate one open momentum and regard the closed momentum $q$ as independent. Keep $k+l$ fixed. In the dual local complex $u\partial_\Omega$ acts by $q\wedge$ from a $k$-form coefficient to a $(k+1)$-form coefficient as the descendant index decreases. Set $h$ equal to contraction with each coordinate tangent vector followed by differentiation with respect to the corresponding component of $q$, summed over the five coordinates. On a $k$-form homogeneous of polynomial degree $d$ in $q$ the product rule gives
$$
(q\wedge)h+h(q\wedge)=(d+5-k)\,\mathrm{id}.
$$
At fixed odd $k+l=1,3$, the complex terminates at the primary rank $k+l$; at lower ranks the displayed homotopy removes every cycle. At fixed odd $k+l\ge5$ the full five-coordinate sequence is present, and its only cohomology is the top form at $q=0$. Hence beyond the primaries the one-external-field linear complex leaves only constant top modes $u^{2j}\mathrm{PV}^{5,*}$, $j\ge1$. This reduction precedes the nonlinear open and closed Ward differentials. It uses polynomial differentiation, so each eliminated local term has a polynomial primitive. [Exterior-algebra replay](descendant_koszul.py)

### Signed primary BV sewing

Order the odd closed plane-wave coefficient $\chi$ before the gauge coefficients, followed by the pairs of open-edge slots (tail, head), and put $\chi$ first inside its disk. The five preceding gauge cubics give $15$ odd slots. Moving $\chi$ through them reverses the open-only Wick sign, giving $+1$ for $k=1,3,5$. The ordered disk momenta are its external momenta, the arriving current and minus the departing current. Modulo the closed momentum their sum vanishes; the unnormalized alternating disk derivative therefore gives $-(k+1)!$ times the minor on the first $k$ slots. Modulo the disk external momenta and the closed momentum, the arriving current is minus the integer-weight covector in $W$ divided by six. The ordered positive Gaussian heat factor from the hexagon is $1/(120(2\pi i)^5)$. The signed coefficient on $W$ is consequently $+k(k+1)/[1440(2\pi i)^5]$ times the first-disk coefficient in the component convention (16).

Use $Q(A)=+A^2$ from Costello–Li (10). On open coefficient functions this is minus the ordinary Chevalley–Eilenberg differential used by the replay. The odd Leibniz rule then gives $Q(\chi H)_{\rm open}=+\chi\,d_{\rm CE}H$. Thus the rank-three single-trace entries, in units $1/[c(2\pi i)^{10}]$, sum to
$$
\frac{1}{96\cdot120}\,W_{(7)}-\frac{1}{2\cdot5760}\,W_{(7)}=0.
$$
For a pure-open ordinary CE preimage $B$ with $d_{\rm CE}B=W$, its contribution under this explicit coefficient-field convention is $QB=-W$. The earlier description as a negative preimage uses $d_{\rm CE}$ as the differential in the counterterm equation. The Hamiltonian convention (10) fixes the opposite conversion on pure-open polynomials and the sign above on a closed-linear term.

For two odd primary coefficients $\chi,\chi'$, a first-disk operation appears in $Q(A)$ with minus $(k+1)$ times its cyclic action coefficient: taking the left $A$-derivative passes through $\chi$. In $Q^2(A)$ the two first-disk terms give $-\chi\chi'[U_1,U'_1]_G$. Inserting $A^2$ in an even-arity second operation and inserting that operation in $A^2$ give $-\chi\chi' bU_2$. Therefore the cyclic second action coefficient is minus the product of the first operation coefficients divided by the number of open slots, for the positive primitive $bU_2=[U_1,U'_1]_G-U_1([\mu,\nu]_S)$. On the double-trace gauge input the internal closed coefficient varies by $-2c$ times the trivector dual to $p_0\wedge p_1$. With external coefficient first, $Q(\chi\chi')=-\chi Q\chi'$. At rank five the saved disk places the internal trivector first, supplying one additional minus. With each external first disk and its common heat unit factored out, the two signed cancellations are
$$
\frac{80}{120}-\frac{2\cdot 4^2\cdot2}{96}=0,
\qquad
-\frac{48}{48}+\frac{2\cdot4\cdot6\cdot2}{96}=0.
$$
The rank-five single-trace wheel has the ordinary CE primitive computed above; its closed-linear counterterm has the opposite wheel coefficient because $Q(\chi H)_{\rm open}=+\chi d_{\rm CE}H$. [Ordered signs and exact factors](primary_bv_sign.py) · [Costello–Li, (7), (10), (16), (31)](https://arxiv.org/abs/1905.09269)
