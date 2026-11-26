# Bulk Kodaira–Spencer graphs on the Calabi–Yau threefold patch

In the compactified action of Fernández–Paquette–Williams §2.3 every bulk vertex takes the trace of the product of its labels in $R=H^\bullet(Y)$, with $Y=K3$ or $T^4$. Every internal propagator uses the inverse of the same pairing. The holomorphic and Dolbeault derivatives act on the threefold coordinates alone. The trace selects degree four in $R$. Both ends of an internal propagator therefore have degrees adding to four. For any connected graph with $|V|$ vertices and $|I|$ internal lines, the external labels must have total degree
$$
4|V|-4|I|=4(1-n),\qquad n=|I|-|V|+1 .
$$
All degrees in $R$ are nonnegative. At $n\ge2$ every individual internal label assignment vanishes. This kills both two-loop theta diagrams, all sixteen tight three-loop graphs, and every higher connected bulk loop, for either compactification and for every choice of the $\mu^3$ and $\alpha\mu\partial\gamma$ vertices. It holds before any spacetime integration or local collision limit.

At one loop all external labels must be multiples of the unit. Contract a spanning tree with the inverse trace pairing. The remaining edge gives the supertrace of multiplication by the product of the external labels. With units, the $K3$ factor is $1+22+1=24$, and the $T^4$ factor is $1-4+6-4+1=0$. Thus the $K3$ bulk box equals 24 times the corresponding $R=\mathbb C$ box on unit labels. For a general label at each leg multiply by its unit coefficient. This argument uses only the compactified bulk action and its propagator. The defect diagrams with fixed backreaction vertices in Appendix A have different vertices. For the separate theory with $R=\mathbb C$ the trace has degree zero, and the two-loop theta kernels require their own calculation. [Action and propagator, §§2.3 and 6](https://arxiv.org/abs/2404.14318) · [Explicit multiplication, supertrace and theta checks](compactification.py)

For the product of $r$ top-component scalar fields at a point the loop bound is $r$ on $K3$ and $r-1$ on $T^4$, and both bounds are attained. At two loops the quadratic point operator has reduced coefficients $2/405$ and zero on the same unit-labeled transverse waves. [Marked ring contractions and exact values](../iib-bcov/README.md#compactification-on-k3-and-t4)

## Minimal theory with scalar coefficients

For the full minimal cubic theory with $R=\mathbb C$, each connected positive-loop four-field local bulk contact at zero displacement vanishes separately after polarizing its external fields. This includes the box, both two-loop thetas and the three-loop $K_4$ subdivisions. The cancellation takes place at fixed positive Schwinger times. It concerns the local bulk functional on divergence-free polyvectors.

**Edge kernel.** Write $\theta_v$ for the three odd polyvector coordinates at a bulk vertex, and $F_i$ for a field with conserved holomorphic momentum $p_i$, $1\le i\le4$. The Fourier form of $\partial_\Omega F_i=0$ is $\iota_{p_i}F_i=0$. The source's propagator is the divergence applied to the Dolbeault propagator. On the difference of two odd coordinates its extra edge factor is
$$
\iota_{y_e}\,\delta^3(\theta_v-\theta_w)
=\int d\eta_e\,\delta^3(\theta_v-\theta_w+y_e\eta_e)
$$
up to the common choice of Berezin orientation. Here $\eta_e$ is one odd integration variable and the odd delta is the product of its three arguments. Holomorphic derivatives on heat kernels insert the Kirchhoff currents $y_e$. The source's scalar potential $\alpha$ gives a divergence-free degree-two polyvector by holomorphic differentiation, $\gamma$ has degree zero, and $\mu$ degree one. This kernel includes their mixed contractions.

For $n$ loops the scalar contact count gives $E=3n+1$, $V=2n+2$ and $3V-2E=4$ external legs. A bridge carries no cycle current variation, and the independent projective Jacobian then vanishes. In particular a vertex with two external legs cannot contribute. In a nonzero graph put the four external fields at its four distinct bivalent vertices.

**Odd integrations.** Choose a spanning tree and let $C$ be the $n$ by $E$ signed cycle matrix. Integrating the differences of the vertex odd coordinates on the tree leaves the $3n$ equations whose matrix has column $C_{\cdot e}\otimes y_e$. Their right kernel contains the Schwinger vector $t$, by the voltage law. A rank below $3n$ makes their Grassmann product zero. At full rank, change the odd edge variables into these $3n$ linear forms and one remaining odd variable $\xi$ with $\eta_e=t_e\xi$ on the kernel. Direct exterior expansion gives the signed maximal minor omitting any edge $e$, divided by $t_e$. Integrating the tree differences then leaves
$$
\int d^3\theta\,d\xi\;\prod_{i=1}^4 F_i(\theta+\Phi_i\xi),
$$
where $\Phi_i$ is the Kirchhoff potential at the assigned vertex and $\theta$ is the common odd translation. The sign of every tree and edge orientation is fixed before assigning the labels. The four inputs can be taken to have even total parity by adjoining test coefficients of their component parities. Extracting those coefficients restores the graded permutations.

**Dependence on external momenta.** Fix a reference triple of basis momenta at the first three bivalent vertices and balance the fourth. Every other conserved input is the image of this one by a single linear map on the three-dimensional momentum space. The edge equations are linear and use the same resistance in all three components. The odd maximal minor scales by $\det(p_1,p_2,p_3)^n$, because its rows comprise $n$ triples. The ordinary top form from differentiating the $n$ three-component chord currents scales by the same power. The product scales by the even power $2n$. Exchanging any two conserved input momenta therefore leaves this scalar prefactor fixed. The argument extends to dependent momenta by polynomiality.

The remaining integral is linear in the four potentials because $\xi^2=0$. Sum over placements before integrating $t$. The averaged linear map from the four momenta to the four potentials commutes with every permutation of the slots. Grounding the potential at a fixed graph vertex gives a permutation average with equal diagonal entries and equal off-diagonal entries. Conservation then leaves a common scalar times $p_i$ at each slot. Each resulting term differentiates one $F_i$ by $\iota_{p_i}$ and is zero. A change of potential origin shifts the common $\theta$ and integrates to zero as well. This proves the graph-by-graph cancellation at every loop order. Multiplying by the graded-symmetric trace of ring labels gives the same cancellation in the compactified box.

The replay builds each edge's two-form directly on the $3V$ odd coordinates and contracts it with four external polarization factors. For the test momenta the transverse space at each leg has dimension two, giving $\binom84=70$ independent exterior components in total degree four. At fixed rational times 68 components are nonzero in individual assignments for each graph tested. All 70 vanish after summing the 24 assignments in the box, the two thetas, all six $K_4$ subdivisions and all 50 cubic four-loop graphs. The ordinary current Jacobian is recomputed independently for each placement. [Exact replay](check_bulk.py) · [Propagator and compactified action](https://arxiv.org/abs/2404.14318)

**Integrated fixed-label periods.** There is also a compact expression before summing external placements. In the odd integral above substitute for each $\Phi_i$ the $i$th matrix row applied to the four momenta. Multiply this integral by $\det(p_1,p_2,p_3)^{2n}$. With the repository's scalar kernel orientation, raw odd differences $\theta_v-\theta_w$, exterior test factors in input order, and the graph's final vertex grounded, the three matrices are
$$
-\frac{1}{120}
\begin{pmatrix}3&2&1&0\\2&4&2&0\\1&2&3&0\\0&0&0&0\end{pmatrix},
\qquad
\frac1{3024000}
\begin{pmatrix}40&10&16&4\\10&25&4&1\\16&4&40&10\\4&1&10&25\end{pmatrix},
\qquad
\frac1{2721600}
\begin{pmatrix}25&5&12&3\\5&25&12&3\\12&12&48&12\\3&3&12&30\end{pmatrix}.
$$
The graphs are the directed box and the $(1,3,3)$ and $(2,2,3)$ thetas, in that order. The theta paths and their interior vertices have the order of the project kernel. These are reduced labeled kernels before automorphism factors and common BV and heat normalizations.

For the box, the grounded potential between two vertices at ordered cumulative resistance fractions $u\le v$ has Green function $u(1-v)$. Its simplex integral gives the first matrix. On each theta path, the tail resistance fractions $\sigma_i$ are decreasing along the path. Write $w_P$ for its fraction of the total conductance. The Green matrix entries are $\sigma_i\sigma_j$ plus $\min(\sigma_i,\sigma_j)(1-\max(\sigma_i,\sigma_j))/w_P$ when both inputs lie on the same path $P$. Multiplying the scalar chart density by the odd maximal minor gives, up to its checked orientation, $Q(L)^2\prod_P w_P^{4-n_P}$, with $L=\sum_{i\in\mathrm{interior}}\sigma_i p_i$ and the determinant $Q$ of the project theta reduction evaluated on three coordinate momenta and their balancing fourth. The two powers of $Q$ arise from the two independent current minors. Integrating the Green matrix against this polynomial gives the stated values. For the four vector polarizations $(0,1,0),(0,0,1),(1,0,0),(1,-1,0)$ and the three coordinate momenta, the labeled theta kernels are $1/1008000$ and $1/1360800$. [Two integration checks and raw exterior signs](theta_periods.py)

## The one-descendant four-field interaction on the flat threefold

Use scalar coefficients and the flat worldsheet normalization. For four conserved three-dimensional wave momenta let $V_{g;e_i}$ be the polarized genus-$g$ vertex with one cotangent-line class at mark $i$. At every $g\ge1$ it equals a strictly positive constant times $\det(p_1,p_2,p_3)^{2g}$. The constant is independent of $i$. The complete constants through three loops are
$$
\frac1{24},\qquad \frac1{9720},\qquad
\frac{411065203}{94152554496000}-\frac{3617}{1009008000}\zeta(3).
$$
The three-primary positive-genus vertex vanishes by momentum rank. These four-field coefficients are the first interactions with one BCOV descendant. On one elementary scalar, the one-edge tree reads the four-field vertex with the usual root-momentum contraction against a vector polarization; the external states can be two scalars (one a level-one descendant) and one transverse vector. For instance on $p_1=e^1,p_2=e^2,p_3=e^3$ with the vector at the third slot polarized by $e_1$, the absolute root contraction is one. These periods concern the scalar-coefficient flat theory. The compactified Frobenius rings impose the loop truncations above. [Li–Wang, Lemma 1.2 and Proposition 1.8](https://arxiv.org/abs/2610.21927) · [Fivefold marked-square integrals](../iib-bcov/README.md#the-three-loop-scalar-square)

Indeed in (1.18) choose the radial variation $\eta_i=p_i$. A connected component of positive genus on its right side would have three primary inputs with total momentum zero and hence rank at most two. Only ordinary cubics survive. A graph contributing to the full three-dimensional current determinant must have its four marks at four different vertices, since its currents are linear in the total source at each vertex. The bare graph has $3g+1$ edges. Write $U$ for its Kirchhoff polynomial, $c$ for its signed projective current cofactor, $j_e$ for its minimizing current and $\lambda_e$ for its normalized edge lengths. The radial first column of the flux determinant is $D_\eta f_e=|j_e|^2$. Cofactor expansion gives $c^2U^{-3}\sum_e\lambda_e|j_e|^2$ with the boundary orientation of (1.17). The energy is also $-\sum_{i<j}(p_i\cdot p_j)R_{ij}$ with $R_{ij}$ the effective resistance. At the simplex $e^1,e^2,e^3,-e^1-e^2-e^3$ every label permutation has determinant character one after squaring. Summing all 24 assignments of the energy gives $2\sum_i|p_i|^2\sum_{i<j}R_{ij}$. Consequently the coefficient at the simplex is
$$
\sum_{\substack{\Gamma\\\mathrm{four\ unlabeled\ bivalent\ positions}}}
\frac{2}{|\operatorname{Aut}\Gamma|}
\int_{\Delta_\Gamma}\frac{c^2}{U^3}\sum_{\substack{i<j\\\mathrm{marked\ positions}}}R_{ij}.
$$
Every term is nonnegative. The determinant character under any real change of the target basis and permutation symmetry of the simplex give the stated dependence on momenta and on $i$. The integral is half the reduced threefold marked-square period on this probe. Contracting the independent labeled two-loop Green matrices gives $1/20160+29/544320=1/9720$. The three-loop sum uses all nine graph integrals in the linked replay.

### Two cotangent-line insertions on a zero wave

Add a fifth mark $x$ of momentum zero to the three-dimensional simplex and put both cotangent-line insertions there. The coefficients of $V_{g;2e_x}$ at genera one, two and three are
$$
\frac16,\qquad \frac{59}{77760},\qquad
\frac{84137363}{4707627724800}-\frac{3617}{252252000}\zeta(3).
$$
In the Ward equation (1.18) take $k=e_x$. At the hard sites $i,h,k'$ put $e^2,e^1,-e^1-e^2-e^3$, respectively, with $e^3$ at the remaining site. Move $\epsilon e^2$ from $k'$ to $x$ and use the source variation $e_2$ at $x$ and $-e_2$ at $h$. The linear coefficient on the left is $V_{g;2e_x}(0)$. A positive-genus primary component has three conserved waves and vanishes by rank. A positive-genus component carrying the descendant has four waves, one soft, so all its three-row determinants vanish at zero soft momentum. Its product of $2g_v$ determinants, for component genus $g_v>0$, vanishes at least quadratically. Only ordinary cubics and the single quartic with descendant at $x$ enter the derivative.

Removing $x$ leaves the four-hard cubic graph. At a vertex $w$ attach the quartic by adding the source $\epsilon e^2$ there. If two hard waves share a vertex the unperturbed normal current determinant vanishes, so its square has zero derivative. At four distinct hard vertices write $c$ for the reduced current cofactor. With the resistance $R$ on the graph, the differentiated flux at fixed lengths, stripped of $U^{-3}$, is
$$
c(\partial_\epsilon c)\bigl(R_{wk'}-R_{wi}+R_{hi}-R_{hk'}\bigr)
+\tfrac12c^2\bigl(R_{wk'}+R_{hw}-R_{hk'}\bigr).
$$
Sum over all original vertices $w$ and all 24 orders of the hard slots, dividing by the graph automorphisms. At hard attachment $w$ the determinant character gives $\partial_\epsilon c=g c$ for $w=i$ and zero for the other three. The four hard attachments contribute $2(g+1)$ times the one-descendant simplex vertex. The one-loop box has only hard vertices. At genera two and three the replay also integrates the internal attachments and checks the linear coefficient independently with the original currents. [Exact differentiated diagrams](zero_wave_descendant.py)

**Independent electrical check (Codex).** For the genus-four combination $X_4-\langle\epsilon+\varphi\rangle$, fix the nine core lengths to $1,\ldots,9$ in the sorted core-edge order. Solve the full subdivided electrical network over the polynomial ring in the four segment fractions. The cycle Gram matrix is constant on each product of simplices. Form the twelve-row current determinant and differentiate it at each possible attachment using its adjugate. For the metric-graph weight compute Zhang's $\tau=\tfrac14\int_\Gamma(\partial_x R(0,x))^2\,dx$ from the endpoint resistances on each edge. Writing $L_e$ for the length of edge $e=(i,j)$ and $R_{ij}$ for its endpoint effective resistance, the resistance at fractional position $s$ is $(1-s)R_{0i}+sR_{0j}+s(1-s)(L_e-R_{ij})$. Hence
$$
4\tau=\sum_{\substack{e=(i,j)\\\mathrm{core\ edges}}}
\frac{(R_{0j}-R_{0i})^2+(L_e-R_{ij})^2/3}{L_e}.
$$
Exact monomial integration on the ordered paths gives the following four sums of numerators, including the graph symmetry factors.

| Four-loop core | Number of placements | Sum at lengths $1,\ldots,9$ |
|---|---:|---:|
| Two double edges | 3 | $1026264360551841792000000/49$ |
| One double edge | 8 | $1090813247522259129751783968/245$ |
| Triangular prism | 28 | $1505136918935837291404100399955712677593/69457500$ |
| $K_{3,3}$ | 11 | $11633374558566388404968743130846960836/17364375$ |

At genera two and three every trivalent diagram also agrees separately at successive integer lengths. The two genus-two numerators sum to $797/2$ and the nine genus-three numerators sum to $1078033588583880388/2625$. All eleven two-double and one-double polynomials agree with the independent electrical evaluation at the displayed unequal lengths. The orbit-weight check undoes the nine factorial cut factors and averages the nine tests with exactly one edge of length two and the other eight of length one. This average survives folding by graph automorphisms. [Exact electrical determinant and simplex integrations](arakelov_electrical_controls.py) · [Fourfold Arakelov reduction](../iib-bcov/README.md#reduction-of-the-fourfold-period-to-five-and-six-worldsheet-marks)

The [exact auxiliary graph periods](four_loop_arakelov_periods.json) record the integrated coefficients of $1,\zeta(3),\zeta(5)$ for shapes with an independent electrical check. Shape indices are zero-based within each core. All 50 amplitudes now have independent electrical checks. Their three period columns determine the [seven-planar elementary response](../iib-bcov/README.md#complete-four-loop-response-on-seven-planar-waves) and the [axial D3 pole](../twisted-holography-N4/README.md#two-polynomial-insertions-determine-the-four-loop-highest-pole-on-this-family).

For shape 06 on $K_{3,3}$, the $298135$ nonzero folded projective moments in the equal-strength combination $X_4-\langle\epsilon+\varphi\rangle$ give
$$
\frac{159269029990890087561293053}{7506757745155113984000000000}
+\frac{218009830177606408021}{3841799362521120000000}\zeta(3)
-\frac{1252211884146317}{14519271967200000}\zeta(5)
\ \in\ \frac1{10^{16}}(19801,19802).
$$
The determinant is quartic in a single source-component shift. Its derivative at zero is $[8(N(1)-N(-1))-N(2)+N(-2)]/12$. The 18 derivatives and the undifferentiated contribution can be integrated separately. The largest derivative was split by the four possible powers of one left current's core length. Computing the adjugate in a separate process and loading its signed derivatives bounds the resident memory of the subsequent path integrations. The polynomial at lengths $1,\ldots,9$ equals $102487511604459881061448835072/214375$. Averaging over the nine single-edge doublings of unit lengths gives $1028003829689/7717500$. Doubling two incident edges gives average $18318620357510149/2222640000$ and doubling two disjoint edges gives $3468913745361241/158760000$. The grounded current calculation reproduces each value independently. [Partial numerator and gradient cache](four_loop_arakelov.py) · [Exact streaming merge](arakelov_piece_weights.py) · [Calibrated period contraction](bipartite_descendant.py)

For $K_{3,3}$ shape 09 the nine stabilizer orbits of the eighteen current derivatives, numbered by three coordinates at each branch vertex, are
$$
(0,3),(1,5),(2,4),(6,12),(7,14),(8,13),(9,15),(10,17),(11,16).
$$
A graph automorphism fixing the balancing mark permutes the three source coordinates and the branch vertices. The current and its derivative have the same determinant sign, which cancels in their product. The integrated path fractions and the resistance factors transform by the same core-edge permutation. Thus the undifferentiated term has multiplicity one and each representative has multiplicity two in the folded weights. The first three pairs also agree row by row when calculated independently. At distinct lengths the partners require permutation of the nine lengths before summing. This gives the unsymmetrized numerator $15617000920857546638041468255191424/42875$ and reproduces the grounded-network value. The single-edge doubling average is $218895578704121/4939200$. The incident and disjoint double-edge averages are $257569700779887257/204120000$ and $697617511790903159/476280000$.

The $287038$ folded moments give the reduced diagram amplitude
$$
\frac{9587937603777493795567160447}{609924066793853011200000000000}
+\frac{98678601073393341668639}{2346990156012902400000000}\zeta(3)
-\frac{1870418157900253769}{29270852285875200000}\zeta(5)
\ \in\ \frac1{10^{16}}(13041,13042).
$$
For a source shift in one coordinate, precisely four rows of the augmented current determinant vary. Replacing one row by its particular-flow derivative at a time and adding the four determinants computes the derivative at the unshifted source. The centered source shifts for piece 06 crossed $1595328$ KiB. Row replacement finished at $623232$ KiB. The twelve directly prepared derivative polynomials each agree with finite differences of the original signed grounded current minor at two unequal metrics. Use `--prepare-gradients --gradient-rows --piece` with the selected piece and its cache directory. [Cycle-row derivative](four_loop_arakelov.py) · [Augmented determinant](../nonSP/verify/k4_periods/current_minor.py) · [Orbit-weighted merge](arakelov_piece_weights.py)

### The planar fivefold period measures Zhang's graph invariant

At genus $g\ge2$ the planar bulk period satisfies $2J_g=-\langle\varphi\rangle$ diagram by diagram. The brackets integrate against the four-hard current measure. Let $R(x,y)$ be the effective resistance at two points of the metric core, $\mu$ the cotree probability measure, and $\beta$ the measure of mass $1/2$ at each cubic vertex. The signed canonical measure $g\mu-\beta$ has total mass one and constant resistance potential $\rho$. On interiors of oriented edges with arclength coordinates $x,y$, the mixed derivative off the diagonal is $2C_e^{\mathsf T}L^{-1}C_f$, where $C_e,C_f$ are their cycle columns and $L$ is the cycle Gram matrix. Indeed the resistance is the squared norm of a path current minus the squared norm of its cycle projection. The path length is affine on either side of the diagonal. Moving the two endpoints differentiates its projected current with opposite signs.

Cut the graph at an interior $y$. As a function of $x$, $\partial_y R$ is affine on each open segment, satisfies Kirchhoff at the original vertices, and has opposite limits of absolute value one at $y$. In the distributional second derivative on the metric graph, $\partial_x^2 R=2\delta_y-2(g\mu-\beta)$. On each edge this follows from the quadratic resistance interpolation, with vertex coefficient $\deg(v)-2$. Integrating by parts twice on the cut graph, using $R(y,y)=0$, gives
$$
\int R(x,y)\big(\partial_x\partial_y R(x,y)\big)^2\,dx
=1-\int \big(\partial_y R(x,y)\big)^2\,d(g\mu-\beta)(x).
$$
Integrating by parts once with $x$ fixed also gives $\int(\partial_y R(x,y))^2dy=2\rho=4\tau$, with $\tau$ the graph tau constant. Hence, writing $\ell$ for total length and $R_{\mu\mu}$ for mean resistance of two independent draws from $\mu$,
$$
\frac14\iint R(x,y)\big(\partial_x\partial_y R(x,y)\big)^2\,dx\,dy
=\frac{\ell}{4}-\tau,\qquad
gR_{\mu\mu}=\rho+\frac{\epsilon}{2}.
$$
The second identity follows by expanding $g\mu=(g\mu-\beta)+\beta$ in $R_{\mu\mu}$ and using $\epsilon=2((g-1)\rho+R_{\beta\beta})/g$. The planar pair law has density $(C_e^{\mathsf T}L^{-1}C_f)^2/g$. Twice its resistance contribution to $J_g$ is therefore $\ell/4-\tau-gR_{\mu\mu}=\ell/4-3\rho/2-\epsilon/2=-\varphi$. All eight cubic cores at genera two through four satisfy the cleared polynomial identity coefficientwise. [Planar pair law](../iib-bcov/planar_bulk.py) · [Exact normalization replay](arakelov_pole_normalization.py)

On the conserved planar probe for the elementary scalar $\phi$ in the flat fivefold worldsheet theory, its genus-$g$ local quantum differential consequently measures the same graph invariant:
$$
(-1)^g(2\pi i)^{5g}Q_g\phi(0)
=-\frac{g!(g-1)!(g-1)^2(2g-1)}{2(3g-1)}\langle\varphi\rangle<0 .
$$
The reduced two-loop value is $-1/17280$. The four-loop value is $-4536\langle\varphi\rangle/11$. [Momentum configuration and allowed changes of the marked scalar](../iib-bcov/README.md#an-elementary-scalar-period-at-each-genus)

There is a lower bound in terms of the one-descendant vertex $D_g$. Write $\mu_e=(1-R_e/t_e)/g$ for the probability of choosing edge $e$ of length $t_e$ with endpoint resistance $R_e$. The midpoint contribution to $-2J_g$ is a nonnegative contraction of two Gram matrices. The same-edge contribution to $\varphi$ is $g(g-1)\sum_{e\in\mathrm{edges}}t_e\mu_e^2/6$. On a bridgeless core an alternative path has length at most $\ell-t_e$, hence $R_e\le t_e(\ell-t_e)/\ell$ and $\mu_e\ge t_e/(g\ell)$. Convexity on the $3g-3$ edges gives
$$
\varphi\ge\frac{g-1}{6g\ell^2}\sum_{e\in\mathrm{edges}}t_e^3
\ge\frac{\ell}{54g(g-1)}.
$$
Two edge-disjoint paths between any two marked points have total length at most $\ell$, so their parallel resistance is at most $\ell/4$. Summing the six marked-pair resistances gives $D_g\le\langle\ell\rangle/8$. Thus the elementary scalar matrix element and the first threefold descendant obey
$$
(-1)^{g+1}(2\pi i)^{5g}Q_g\phi(0)
\ge\frac{2((g-1)!)^2(g-1)(2g-1)}{27(3g-1)}D_g>0 .
$$
The single ladder below supplies an explicit positive lower bound on $D_g$ at every genus. [Midpoint Gram and same-edge terms](../iib-bcov/quartic_planar_cores.py) · [Exact ladder amplitude](descendant_ladder.py)

For the fivefold D3 contact, the fourfold soft coefficient on the first coordinate covector $e^1$ is the five-mark integral times $(g-2)!$. In the threefold conventions for the one- and two-descendant integrals $D_g,X_g$ and the transverse bulk period $K_g$, this gives
$$
C_g(e^1)=(g-2)!\left(\frac{X_g-4D_g-\langle\epsilon\rangle}{8}
                       -\frac{K_g}{2(g!)^2}\right).
$$
Thus $72C_4(e^1)=18(X_4-4D_4-\langle\epsilon\rangle)-K_4/8$ agrees identically with the seven-plane Ward projection. Independently, the coefficient of $A^2B$ for tangent strengths $A,B$ in the fivefold trace is obtained from its normal shears $x,y$ as $-2[x^2y^2]$ of the signed trace at $A=B=0$. At genera two and three its reduced values are
$$
4C_2(e^1)=\frac1{15552},\qquad
16C_3(e^1)=\frac{740466719}{26480405952000}
             -\frac{469153}{20432412000}\zeta(3).
$$
Both independent period channels agree with the integrated threefold formula and the fourfold primary vertex. [Fivefold trace and D3 normal pole](../twisted-holography-N4/README.md#a-separate-five-mark-period-for-the-second-insertion) · [Fourfold graph-to-worldsheet conversion](../iib-bcov/README.md#four-hard-diagram-formula-for-the-remaining-elementary-scalar-period)

For the five- and six-mark Abel integrals $I_5,I_6$ at genus $g$ above, integration along the moving marks gives
$$
\begin{aligned}
8I_5&=X_g-4D_g-\langle\epsilon\rangle-\frac{4K_g}{(g!)^2},\\
12I_6&=(3g-1)(X_g-4D_g-\langle\epsilon\rangle)+4J_g,\\
3I_6-2(3g-1)I_5&=\frac{3g-1}{(g!)^2}K_g+J_g .
\end{aligned}
$$
At unit line strengths the fourfold vertex is $(g-1)!((g-1)I_5-I_6/2)$ and adjoining $g$ equal fivefold waves multiplies its magnitude by $g!$ on an elementary-scalar probe with unit root contraction. The stripped axial D3 coefficient is $4((g-1)!)^2 I_5$. Subtracting $4/[g(g-1)]$ times the elementary magnitude from the axial coefficient leaves $2((g-1)!)^2I_6/(g-1)$. This is zero at genus two and strictly positive for $g\ge3$, by the fourfold soft criterion on the line strengths $(1,-1,t,\ldots,t)$ with $g-3$ soft marks. In particular the four-loop difference equals $2V(1,1,-2)$. [Line-wave reduction and soft positivity](../iib-bcov/README.md#reduction-of-the-fourfold-period-to-five-and-six-worldsheet-marks)

In increasing loop order the combination $3I_6-2(3g-1)I_5$ equals
$$
\begin{aligned}
g=2:\quad &-\frac5{31104},\\
g=3:\quad &-\frac{3399295819}{211843247616000}
              +\frac{1076227}{81729648000}\zeta(3),\\
g=4:\quad &-\frac{100224698200645422317757433}{3414152210624483328000000000}
           -\frac{31025870232372083393461}{395156505859315200000000}\zeta(3)
           +\frac{93558855827167151}{784040686228800000}\zeta(5).
\end{aligned}
$$
At genus two $I_6=0$ by rank. The genus-three vertices give both moments independently. At genus four $W=18I_5-3I_6$ and $V(1,1,-2)=12I_6$, giving the independent seven-plane identity $V(1,1,-2)+22W=-11K_4/32-18J_4$. The displayed four-loop value lies strictly between $-73152318930388/10^{24}$ and $-73152318930387/10^{24}$. For the separate mixed ray $e^1+e^2$, the graph of hard-mark pairs in its five-mark Betti product has degrees $(3,2,1,1,1)$, whereas the axial product has degrees $(4,1,1,1,1)$. Permuting the five marks preserves these degrees. The symmetric quartic separating the two restrictions evaluates to $0$ and $48$ on the axial and mixed probes. [Five-mark products and separating orbit](../twisted-holography-N4/README.md#a-separate-five-mark-period-for-the-second-insertion) · [Planar period contraction](../iib-bcov/quartic_scalar_response.py)

### Three integer four-loop elementary-scalar measurements span the graph periods

Use the flat IIB scalar at the origin with absorbed momentum $e^4$ and the vector at $e^3$ polarized by $e_4$. Put $u=e^1$, $v=e^2$, $w=-4u-3v$. Use seven planar scalar momenta in each row below, with multiplicities in parentheses. Append the scalar at $e^5$ and the balancing scalar. The rows give the coefficients of $K_4,J_4,W$ in the heat-stripped local operation $(2\pi i)^{20}Q_4\phi(0)$.

| Seven planar momenta | $K_4$ | $J_4$ | $W$ |
|---|---:|---:|---:|
| $u^{(3)},v^{(3)},w$ | $117/2$ | $4320$ | $3720$ |
| $u^{(4)},v^{(2)},w+v$ | $33/4$ | $432$ | $528$ |
| $u^{(4)},v^{(3)}$ | $0$ | $0$ | $-24$ |

The table has determinant $248832$. In the $1,\zeta(3),\zeta(5)$ period columns the determinant of these same three individual measurements is $99665298\bmod 999999937$. Thus their rational linear combinations give each of $1,\zeta(3),\zeta(5)$. Any common nonzero normalization making the three measurements rational would make the ratio $\zeta(3)/1$ rational, contradicting Apéry's theorem. This uses no arithmetic independence hypothesis on $\zeta(5)$. The marked-field invariance here keeps the flat bulk action and output coordinates fixed. [Seven-wave bracket and marked invariance](../iib-bcov/README.md#a-four-loop-measurement-on-the-elementary-scalar) · [Evaluator for any seven rational planar waves](arakelov_elementary_measurements.py)

All fifty individual auxiliary graph periods have positive lower endpoints from the ordinary positive zeta series. Reconstructing $W$, $V(1,1,-2)$ and $72C_4(e^1)$ from the five- and six-mark moments reproduces the Ward values in all three columns. The last elementary row lies in $(-100640984462,-100640984461)/10^{20}$. The centered bare fourfold Ward integral inferred by subtracting its tree, polygon and theta terms lies in $(-9272874925653,-9272874925652)/10^{20}$. [All diagram coefficients](four_loop_arakelov_periods.json) · [Centered fourfold Ward equation](../iib-bcov/README.md#transverse-order-along-a-collinear-fivefold-configuration)

### Complete four-loop descendant

In the same simplex normalization all 50 nonzero diagrams give
$$
D_4=
\frac{192761752316371576853}{26469813789619200000000}
+\frac{39781234295806109}{2042423903520000000}\zeta(3)
-\frac{411294138769}{13894040160000}\zeta(5),
\qquad
\frac{110899269775}{10^{22}}<D_4<\frac{110899269776}{10^{22}}.
$$
The genus-four cores have respectively $0,3,8,28,11$ contributing placements with three double edges, two double edges, one double edge, a triangular prism and $K_{3,3}$. For each placement integrate the four bivalent positions at fixed core lengths and weight the square of the signed current polynomial by the six marked-pair resistances and $2/|\operatorname{Aut}\Gamma|$. An independent Dirichlet integral in the four unsuppressed positions gives the same polynomial at unequal core lengths. The resulting nine-length integrand has denominator $U^{12}$. The scalar propagator period on powers $a_i$ carries the projective factor $\prod_i(a_i-1)!/11!$.

The 28 prism shapes give 6,377,268 exponent rows; folding them to 240,735 reflected classes reproduces all three integrated coefficients independently. The eleven $K_{3,3}$ shapes give 1,690,458 exponent rows. Each avoids the support of the one-dimensional homogeneous integration-by-parts solution. The full three period columns also reproduce the independent dimension lifts and the original fivefold numerator integrals. Each of the fifty physical integrals is positive using rational enclosures of $\zeta(3),\zeta(5)$. In the fivefold [second-order differential](../iib-bcov/README.md#the-fixed-scalar-square-detects-the-threefold-descendant-at-every-genus) the same sum multiplies $2(\det H)^4 f''(0)/(2\pi i)^{20}$. [Graph integrands and coefficient verification](four_loop_descendant.py) · [All fifty exact amplitudes](four_loop_descendant.json) · [Bipartite moment integration](bipartite_descendant.py) · [Prism period equations](../iib-bcov/quartic_prism_physical_ibp.py)

### One integrated ladder at every genus

Take two rows of vertices $a_0,\ldots,a_g$ and $b_0,\ldots,b_g$, with rung lengths $h_i$ and upper and lower rail lengths $r_i,s_i$ between consecutive columns. Put the four simplex momenta at the corners. On arbitrary positive lengths the squared current numerator factorizes:
$$
(U^g c)^2=\prod_{i=0}^{g}h_i^{\,2g-2}.
$$
To see this, let $d_i$ be the vector potential difference between the ends of rung $i$, $A_i$ the upper rail current from column $i$ to $i+1$, and $S$ the total rail current, independent of $i$. With $l_i=r_i+s_i$ the transfer is $d_{i+1}=d_i-l_i A_i+s_i S$, $A_{i+1}=A_i-d_{i+1}/h_{i+1}$ at interior columns. Hence $\det(S,d_i,A_i)$ is constant. Delete the first rung column in the cycle-current matrix and order rows by the elementary rectangular cycles. The matrix is block triangular. Block $i$ has the two rail currents and $d_{i+1}/h_{i+1}$; its absolute determinant is $|\det(S,d_i,A_i)|/h_{i+1}$. Divide by $h_0$ for the projective cofactor. At the left boundary $S=p_{a_0}+p_{b_0}$ and $A_0=p_{a_0}-d_0/h_0$, so the constant determinant in the simplex basis is the transfer resistance of the two end rungs. Let $B$ be a path incidence matrix (one column per interval) and let $M=\operatorname{diag}(1/h)+B\operatorname{diag}(1/l)B^{\mathsf T}$. For a unit right rung source the total rail current is zero and the rung differences solve this tridiagonal system. Thus the transfer resistance is $(M^{-1})_{0g}=(\prod_i l_i^{-1})/\det M=\prod_i h_i/U$. The last equality follows from $U=\det(\operatorname{diag}(l)+B^{\mathsf T}\operatorname{diag}(h)B)$. This proves the factorization and strict positivity of the full descendant at every genus.

The complete ladder contribution, including its symmetry factor, equals
$$
\frac{2}{|\operatorname{Aut}\Gamma|}
\frac{6^{g+1}}{(2g)!(2g+1)!(2g+3)!}
\left[
4+\frac{2g}{g+2}+
\frac{2(2g+1)}{g(g+2)}
\sum_{\substack{0\le i<g\\\mathrm{rail\ intervals}}}\frac1{2i+1}
\right].
$$
The automorphism order is eight for $g=1$ and four otherwise. The first four ladder values are $1/24$, $1/20160$, $373/114307200000$, $173/3943141171200000$. A separate sum over separating two-forests and direct series–parallel beta integrals reproduces all four.

For the integration put $R_0=h_0$, $R_{i+1}=h_{i+1}(R_i+l_i)/(h_{i+1}+R_i+l_i)$ and fix $R_g=1$ by projective homogeneity. The variables $u_i=R_i/(R_i+l_i)$, $v_i=h_{i+1}/(h_{i+1}+R_i+l_i)$ give $l_i=(1-u_i)R_{i+1}/v_i$ and $h_{i+1}=R_{i+1}/(1-v_i)$. Split $l_i$ uniformly between the two rails. The squared-cofactor density times the projective measure in this gauge is
$$
\prod_{i=0}^{g-1}u_i^{2g-2-2i}(1-u_i)v_i^{2i+1}(1-v_i)^3\,du_i\,dv_i .
$$
Indeed the resistance Jacobian before splitting the rails is $\prod_i R_{i+1}^2/(v_i^3(1-v_i)^2)$, the splitting produces $\prod_i l_i$, and $U=\prod_i R_{i+1}/(v_i(1-v_i))$. The product beta integral is $6^{g+1}/[(2g)!(2g+1)!(2g+3)!]$.

Writing $G=M^{-1}$, the uniform rail-split mean of the six corner resistances is
$$
2(G_{00}+G_{gg})+\frac23\sum_i l_i+
\frac13\operatorname{tr}(B^{\mathsf T}GB).
$$
For a unit flow between opposite ends let $x_i=r_i/l_i$. Completing the square in the rail energy leaves $\sum_i l_ix_i(1-x_i)$ plus the rung system energy. Averaging over the four choices of end corners and independent uniform $x_i$, the rung source vector has zero mean and covariance $(e_0e_0^{\mathsf T}+e_ge_g^{\mathsf T})/4+BB^{\mathsf T}/12$. The other two pairs give $G_{00},G_{gg}$. Reflection makes their projective integrals equal, and $G_{gg}=1$ in the chosen gauge. Each diagonal interval entry of $B^{\mathsf T}GB$ equals $l_i-l_i^2\partial_{l_i}\log U$. In the equivalent gauge $h_0=1$, integration by parts in $l_i$ with density $\prod_i l_i\,(\prod_j h_j)^{2g-2}/U^{2g+3}$ replaces this entry by $2g\,l_i/(2g+3)$ under the integral. The endpoint terms vanish since $U$ is affine in $l_i$ with positive constant and linear terms. The normalized product beta mean is
$$
\mathbb E[l_i]=
\frac{2(2g+1)(2g+3)}
{(2i+1)(2i+3)(2g-2i-1)(2g-2i+1)} .
$$
Summing by partial fractions gives the displayed ladder contribution. [Current factorization, independent resistance cubature and exact beta integrals](descendant_ladder.py)

### Two-loop threefold vertices with moving descendants

At genus two take the five labeled momenta $((1,0,0),(0,1,0),(0,0,1),(1,-1,0),(-2,0,-1))$ in three target coordinates. For the cotangent levels $(2,0,0,0,0)$, $(1,1,0,0,0)$ and $(0,0,0,2,0)$ the respective worldsheet vertices are
$$
\frac{569}{77760},\qquad\frac{739}{77760},\qquad\frac{37}{8640}.
$$
All dual graphs at genus two are series-parallel after taking blocks. The mixed Newton-body moments reduce to integrals of polynomials on cubes. On four unit simplex marks with one cotangent class the same solver recovers $1/9720$. For an independent fivefold Feynman calculation adjoin four labeled unit waves $\pm e^4,\pm e^5$. Cotangent pullback subtracts the single fusions and adds the double fusions into the marked descendants. Each transverse ordered pair has Abel degree two. The three signed sums of ordinary genus-two fivefold graph vertices are
$$
\frac{2021}{38880}-\frac{79}{3240}+\frac{13}{7776},\qquad
\frac{151}{1944}-\frac{73}{1620}+\frac7{1296},\qquad
\frac{1121}{38880}-\frac{121}{9720}+\frac{29}{38880},
$$
each equal to four times its threefold vertex. [Stable-graph integral](../twisted-holography-N4/README.md#mixed-newton-moments-on-series-parallel-graphs) · [Independent fivefold graph replay](../iib-bcov/genus_two_mixed_targets.py)
