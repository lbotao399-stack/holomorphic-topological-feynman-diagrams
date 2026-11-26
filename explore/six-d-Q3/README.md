# The local differential on $\mathbb C^{3|2}$

## Three-loop Ward equation

The heat edge and propagator carry the same odd delta function $(\theta_1-\theta'_1)(\theta_2-\theta'_2)$. Around a cycle the differences in each odd coordinate sum to zero. Their exterior product is therefore zero before the bosonic integral, also when one edge is a heat edge. A tadpole has a zero difference on its single edge. The cubic interaction has no odd-coordinate derivative. Consequently every positive-loop bulk collision vanishes, including its insertions into larger graphs. The tree collision is the Jacobi identity for the adjoint bracket. This proves the bulk quantum master equation for the bare interaction in this heat prescription.

For a local insertion, the quadratic identities of the heat kernels group collision faces into two nested insertions and pure-interaction collisions. The latter sum to the bulk master equation. The former give $\mathbf Q^2\mathcal O=0$ on holomorphic jets. Equivalently, the effective insertion differential is conjugate to $\bar\partial+\hbar\Delta_L$ at positive heat scale $L$, with $\Delta_L$ the BV contraction with the heat edge. This argument uses the dimension-independent collision identity with the three-dimensional kernel; the numbers of cubic vertices here are $2n+1$ at loop order $n$. [Wang, quadratic relations](https://arxiv.org/abs/2401.08113) · [Insertion identity and its proof](../../sources/LLC_main.tex)

At order three, deleting the insertion leaves a tree. There are five allowed rooted classes, each with one trivalent color contraction per vertex. The antisymmetry of the structure tensor makes all leg permutations of a fixed Wick graph scalar signs after the fields and left derivatives are put in their prescribed order. Thus the actual three-loop differential belongs to the five-operator span in the project's [operator conventions](VERIFY.md#c6-wick-tensors-q_3--confirmed-signs-through-c7). The chain and bubble reductions regenerate every coefficient of the five densities and all 100 coefficients of their linear displacement maps as polynomials before integration. This covers arbitrary input jets. [Symbolic comparison](verify/v_symbolic_charts.py)

Write $\kappa$ for the five weights of the integrated kernels in the order $\Gamma_1,\ldots,\Gamma_5$, with the prefactor $i/(512\pi^9)$ of the project. The probes evaluate the anticommutators with $\widehat Q_0$ on four $b$ fields in $\mathfrak{so}(5)$. Their $5\times5$ minor has determinant $6$ modulo $101$; all rational denominators are invertible modulo $101$. The five signed right-hand sides from the lower-order operators give

$$
\kappa=\left(-\tfrac12,\tfrac12,-1,\tfrac18,-\tfrac12\right),\qquad
Q_0Q_3+Q_3Q_0+Q_1Q_2+Q_2Q_1=0 .
$$

The bulk equation gives existence for the universal Wick weights. The nonzero minor fixes them uniquely; the weights are independent of the test Lie algebra. The normalization is relative to the stated $Q_1,Q_2$ and the bare jet representatives. [Integer contraction and rank replay](verify/v_rank.py)

## Ordered Wick weights

Each cubic vertex is an ordered product of three odd superfields. Use the source order at a vertex: incoming ends in edge order, its surviving output, and outgoing ends in edge order. Concatenate these words in vertex order, followed by one odd test letter for each spoke in edge order. This test word reads off the coefficient of the unspecialized left functional derivatives. For each edge remove its tail and then its head. Removing a letter in position $j$ (counting from one) contributes $(-1)^{j-1}$. Pair removals are even, so they commute. Finally divide by the action of the printed $n+1$ left derivatives on the test word, namely $(-1)^{n(n+1)/2}$. The resulting sign multiplies the graph's ordered color tensor and the signed S2 integral. Relative signs at a fixed loop order follow directly from the odd Leibniz rule.

In the real orientation $\Omega\wedge d^3\bar x=8\,d^6x_E$, the ordered *scalar* Gaussian with all vertex volumes already at the front is $+\mathcal I_\Gamma/(8\pi^3)^n$. This corrects the scalar sign $(-1)^n$ in the source equation `eq-B-norm`. For example, for the source-oriented box take the three incoming coordinate covectors and $t_e=1$. Expanding the three contracted antiholomorphic forms and one heat form gives $-1/4$ at every mark in vertex-major antiholomorphic order. Recombining the three holomorphic and antiholomorphic vertex volumes contributes $(-1)^3$; the Kirchhoff polynomial is $U=4$. Each marked density is consequently $+1/[256(8\pi^3)]$.

To check the sign at every loop order, order the columns of the grounded incidence matrix by integrated vertices and orient the fundamental cycles by $\sigma_T=\varepsilon_T\det B_T$. At heat mark $j$ (zero based), delete one coordinate from each of the other three-forms. The complementary minor of the incidence matrix tensored with the three-dimensional identity gives the minor with entries $C_{re}y_e^a$, with cycle $r$ first and coordinate $a$ second. The interior-contraction signs cancel the deleted-coordinate labels. The remaining permutation sign is $(-1)^j$. Differentiating the cycle current equation gives $\sum_e C_{T,re} y_e^a\,dt_e=-M_{rs}\,d y_{c_s}^a$ in the unflipped fundamental-cycle basis, with $M=C_TtC_T^{\mathsf T}$ and cotree edges $c_s$. On all $3n$ rows this supplies $(-1)^{3n}U^3\sigma_T^3$ in the oriented basis. The factor $(-1)^j$ is precisely the positive facet orientation of the S2 chart. Real Gaussian position integration divides by $U^3$ and restores the vertex-wise real measures with sign $(-1)^{\binom{2n+1}{2}}=(-1)^n$. These signs multiply to $+1$. Each contraction contributes $2y$, and each surplus heat normalization contributes $1/[8(2\pi)^3]$, giving $(8\pi^3)^{-n}$. The cycle and coordinate orders agree with the S2 cotree order throughout. [Direct real Gaussian and cotree comparison on all 68 marked facets through three loops](verify/marked_volume.py)

In the full Wick product the volume forms initially stand at their separate cubic vertices. Gathering them before the odd cubic words contributes a further $(-1)^{\binom{2n+1}{2}}=(-1)^n$. The pair of left derivatives is even and the single heat form commutes with all even propagator forms; moving it across the whole vertex-volume block contributes the same sign at every loop order because there are $2n+1$ vertices. The tree fixes this common order convention. Consequently the full coefficient is $i^{2n+1}(-1)^n/(8\pi^3)^n=i/(8\pi^3)^n$ multiplying the ordered fibre contraction and S2 kernel. The box coefficient including its stabilizer is $i/(16\pi^3)$.

For positive loop order every vertex has at most one output. At each vertex the $3!$ assignments of the ordered cubic slots cancel its factorial: permuting slots gives the same sign in the odd fields and in the structure tensor. The labelled orbit then leaves the inverse of the mark-preserving tree automorphism order. The rule reproduces all eight independently fixed weights of the box, the two two-loop graphs and the five three-loop graphs, including $+1/8$ for the fourth three-loop graph. At four loops the multiplicities of automorphism orders $1,2,4,8$ are $3,9,1,2$. At five loops the multiplicities of orders $1,2,4,8,16$ are $12,25,7,4,1$. The five-loop table gives each ordered weight and its product with the scalar integral before evaluation on external colors. [Ordered contractions and tree symmetries](wick_weights.py) · [Source kernel normalization](../../code/htk/CONVENTIONS.md) · [Ordered Wick prescription](../../sources/LLC_main.tex)

## A nonzero matrix element of the full differential at every loop

At loop order $n\ge1$, work in $\mathfrak{sl}_{2n+3}$ with its trace pairing. Use the $2n+3$ successive matrix units on a directed cycle, ending with $E_{2n+3,1}$. Contract the first $n$ units against $b$ and the next against $c$, in this order, forming a product of $n+1$ linear traces. Evaluate the output ghost on the remaining $n+2$ units with independent odd coefficients in cycle order. Give the first output momentum $A$, the last $C$, and all $n$ middle outputs $B$. With $\Delta=\det(A,B,C)$ the coefficient of the ordered product of the output odd parameters in $Q_n$ at the origin is
$$
\frac{i}{(8\pi^3)^n}\,
\frac{(-1)^{n(n-1)/2}2^{n+1}\Delta^n}
{(n+2)!(n!)^2((n+1)!)^2}.
$$
The corresponding term of the full differential carries $\hbar^n$. Embedding the same matrix units into $\mathfrak{sl}_N$ proves nonvanishing for $N\ge 2n+3$.

To see which graphs contribute, expand the cubic color tree in traces of nested commutators. A product of the distinct cycle units has nonzero trace only in cycle order. A leaf of the interaction tree has one spoke and one output. Their units must be adjacent. The cycle has exactly two boundaries between input and output units, so the tree has exactly two leaves and is a path. Adjacent marked vertices would give a triangle with the insertion and violate sparsity. Hence the marks alternate, including both ends. The two color assignments differ by path reflection.

Label the path vertices $0,1,\ldots,2n$. In each box $j$ orient the existing spoke from the insertion, then the two path edges to $2j$, in that order. Orient the terminal spoke from $2n$ back to the insertion. On reducing box $j$, let $w_j$ be the conductance fraction of its three-edge path and $s_{1j}>s_{2j}$ its two tail resistance fractions. Let $A_0=A$ and let $A_j$ be the source transported to $2j$. If $Y_j$ is the total current on the effective spoke directed out of the insertion, then
$$
A_j=(1-s_{1j})A_{j-1}+(1-s_{2j})B,\qquad
Y_{j-1}=w_jY_j+A_j-A_{j-1}-B,\qquad
Y_n=-C-A_n.
$$
The cotree currents are $(1-w_j)Y_j$, with the terminal orientation reversed. In the $C$ component their Jacobian in the $w_j$ is triangular, with diagonal magnitudes $\prod_{k>j}w_k$. At fixed $w$, the two remaining components of $Y_j$ depend on $A_j,\ldots,A_n$ with diagonal $-1$. The $(s_{1j},s_{2j})$ derivatives of $A_j$ have determinant $\prod_{k<j}(1-s_{1k})$ in the ordered $(A,B)$ basis. The terminal sign and the tree incidence sign are both negative. The time chart is positive: at fixed effective scale $R$ its three-edge resistance and direct resistance are $R/w_j$ and $R/(1-w_j)$, and the local four-variable Jacobian is $R^3/[w_j^4(1-w_j)^2]$. The resulting positive density for the forward assignment, divided by $\Delta^n$, is
$$
\prod_{j=1}^n w_j^{j-1}(1-w_j)^2(1-s_{1j})^{n-j}.
$$
The two integrals per box are $2(j-1)!/(j+2)!$ and $1/[(n-j+1)(n-j+2)]$.

For forward output order the input colors run in reverse. Each of the $n-1$ interior marks contributes a minus sign when its color is prepended in the commutator expansion. Reversing the input permutation cancels the sign $(-1)^{n(n+1)/2}$ from differentiating $b^n c$. In the ordered fibre word the first spoke removal crosses $6n+3$ letters in total. The two path edges in box $j$ cross $2j$ and $2j+2$ letters, every subsequent interior spoke crosses $6n+4-4j$ letters for its path vertex $2j$, and the terminal spoke crosses $2n+4$. The Wick sign is therefore $(-1)^{1+n(n+1)/2}$. The reflected color term acquires $(-1)^n$, as does $\Delta^n$. Adding the two orientations cancels the order-two automorphism denominator.

At four and five loops, on coordinate covectors with $\Delta=1$, the physical coefficients in units $i/(8\pi^3)^n$ are $1/186624000$ and $1/587865600000$. This detects the differential on local color components. A cohomological lifting or a gauge-invariant anomalous dimension requires a separate calculation. The replay sums every matrix-unit assignment through three loops and checks the signed density by grounded Kirchhoff derivatives through ten loops. [Exact matrix element](matrix_unit_channel.py)

**Zero total momentum.** Change the last middle momentum from $B$ to $-A-(n-1)B-C$. The forward scalar integral is now the preceding one times
$$
(-1)^{n-1}\frac{(n-1)(3n-2)}4 .
$$
Indeed, writing $A_j=\alpha_j A+\beta_j B$ up to $j=n-1$, the change multiplies the original density by
$$
\left[-\frac{s_{2n}(1-w_n)}{w_n}\right]^{n-1}
\left(\frac{\beta_{n-1}}{\alpha_{n-1}}-(n-1)\right).
$$
This follows from the last recurrence for $Y_{n-1}$, whose $C$ component becomes $s_{2n}(1-w_n)$, and from the last determinant $\det(A_{n-1},-A-(n-1)B-C,C)$. In the normalized last-box measure the moments of $[(1-w_n)/w_n]^{n-1}$ and $s_{2n}^{n-1}$ are $n(n+1)/2$ and $2/[n(n+1)]$. At an earlier box with $k=n-j\ge1$ the expectations of $(1-s_{1j})^{-1}$ and $(1-s_{2j})/(1-s_{1j})$ are $(k+2)/k$ and $(k+1)/k$. Iterating the source recurrence therefore gives
$$
\mathbb E\!\left[\frac{\beta_{n-1}}{\alpha_{n-1}}\right]
=n(n+1)\sum_{\substack{1\le k<n\\\mathrm{boxes}}}\frac1{k(k+2)}
=\frac{3n^2-n-2}{4}.
$$
At every $n\ge2$ the output thus survives the quotient by total holomorphic derivatives. At $n=1$ the conserved three-momentum determinant vanishes.

**Conjugation-invariant insertion.** The same matrix element is realized on a polynomial invariant under simultaneous color conjugation. In the input traces replace each fixed unit $E_{j,j+1}$ by $P_j(g_1)g_2 P_{j+1}(g_1)$. Here $P_j$ is the Lagrange polynomial taking value one at the $j$th of $2n+3$ distinct eigenvalues with sum zero and zero at all the others. Evaluate $g_1$ at the constant diagonal matrix with these eigenvalues and $g_2$ at the cyclic shift. The replacements evaluate to the original matrix units.

A bivalent output vertex carrying zero momentum has equal oriented currents on its two incident edges. Every current depends on their two resistances only through their sum. Varying their difference preserves the projective scale and kills the current Jacobian. Thus a graph emitting a constant scalar at any output gives zero, including before multiplying by an input-jet insertion. Every contributing graph instead emits $n+2$ ghosts. In the chosen coefficient of ghost degree $n+2$, the original $c$ must therefore be differentiated. The background sets $b=0$, so every original $b$ must also be differentiated. These exhaust the $n+1$ functional derivatives. All scalar factors in the input traces remain spectators and the color calculation above applies unchanged. This uses an off-shell scalar background and proves nonvanishing on the conjugation-invariant complex of local polynomials, also modulo total derivatives for $n\ge2$.

## Pure-ghost component identities at every loop order

Give the $m$ ordered input superfields independent copies of the paper's odd coordinate, denoted $\theta_j=(\theta_{j1},\theta_{j2})$ in slot $j$. Project the output onto $c$ and its holomorphic jets by setting $g_1,g_2,b$ and their jets to zero. Every graph, and each ordered composition of graph operators, commutes in the graded sense with simultaneous odd translations. Its odd edge kernels depend on differences, its Berezin measures are translation invariant, and the cubic vertex has no odd derivative. The output projection evaluates on superfields constant in the odd directions. The resulting polynomial in the input coordinates therefore depends only on $\theta_j-\theta_m$, $1\le j<m$. This holds separately at fixed colors and holomorphic momenta.

A graph of loop order $n$ with a single insertion takes $n+1$ functional derivatives at one common odd point. Each derivative on a letter of odd-coordinate degree $d$ contributes $2-d$ powers. Saturating the two integrations makes the total input degree exceed the output degree by $2n$, with a difference of $n$ in each flavor. The differences add under composition, including tree operators at $n=0$. For a pure-ghost output the input polynomial at total loop order $n$ thus has degree $(n,n)$ in the two flavors. It vanishes when $m\le n$. In particular each three-loop composition on $bbb$ vanishes on the five-ghost projection (check C8).

At $m=n+1$ there is a single possible odd polynomial, up to its color and momentum coefficient:
$$
\prod_{j=1}^{n}(\theta_{j1}-\theta_{m1})(\theta_{j2}-\theta_{m2}).
$$
Collecting the odd coordinates to the left in the source superfields gives the coefficient $+1$ for $bbbc$ and $-1$ for $g_1g_2bb$ at $n=3$. For the latter the six selected odd-coordinate factors, in product order, have labels $11,42,41,22,31,32$. The two minus signs from selecting the last slot cancel and ordering the six factors takes seven transpositions. The source components give no further sign for these two words. This proves the term-by-term relation in C10 for arbitrary input jets and colors. In general there are $m^2$ nonzero component coefficients of this universal factor, one for each location of the missing first and second flavors. A single $b^{m-1}c$ computation fixes all of them.

For $m>n$, simultaneous determinant-one changes of odd coordinates also preserve the kernels, vertices, measures and pure-ghost evaluation. In the exterior algebra of the $m-1$ difference doublets, weight zero in degree $2n$ has dimension $\binom{m-1}{n}^2$. Weight two has dimension $\binom{m-1}{n+1}\binom{m-1}{n-1}$. Every nontrivial even highest-weight representation contributes one vector to each of these weights, so their difference counts the invariant tensors. It is
$$
\binom{m-1}{n}^2-\binom{m-1}{n+1}\binom{m-1}{n-1}
=\frac{1}{n+1}\binom{m-1}{n}\binom{m}{n},
$$
with out-of-range binomial coefficients zero. This counts allowed flavor tensors before additional color or momentum identities. [Component signs and exact raising-operator ranks](odd_translations.py)

## Rooted graph count at every loop order

Remove the insertion from a positive-loop graph and mark each vertex originally joined to it. An internal cycle vanishes by the odd delta functions. If the resulting forest has $p$ components, each component with $k$ vertices and $s$ marks obeys $2s\le k+1$ by the subgraph bound after reattaching the insertion. Adding over components gives $2(n+p)\le 2n+1+p$, so $p=1$. This tree has $2n+1$ vertices and $n+1$ marks. Every leaf is marked and degree-three vertices are unmarked.

Give a marked vertex weight $+1$ and an unmarked vertex weight $-1$. For a connected subtree the subgraph bound is precisely that its weight is at most one. After cutting any edge the two weights are integers adding to one; hence they are zero and one. Orient the edge toward the side of weight one. Each marked vertex has outdegree zero and each unmarked vertex has outdegree two. Conversely any tree with these outdegrees and the cubic valence bound satisfies every subgraph inequality: its connected subtrees have weight one minus the number of outward boundary arrows. For several components the individual bounds imply the bound on their union. The orientation is unique.

To enumerate the trees, cut at a marked root. On entry along an incoming arrow, a mark has zero or one further incoming branch, and an unmarked vertex has two unordered outgoing branches. On entry along an outgoing arrow, the unmarked vertex has one other outgoing branch and at most one incoming branch. The replay generates these branches directly and retains one canonical rooting. The unlabelled counts for $n=1,\ldots,8$ are $1,2,5,15,49,180,694,2823$.

For the sum of inverse automorphism orders, an unordered pair of branches contributes one half of the square of their weighted series. Let $x$ count unmarked vertices (loops) and $t$ be $x$ times the series for the branch reached by an incoming arrow. The other branch has series $t/(1-t)$. The first branch recursion gives $x=t(1-t)(1-t/2)$. Rooting at one of the $n+1$ marked vertices gives the series $1+t/(1-t)+\tfrac12(t/(1-t))^2$. Lagrange inversion yields

$$
\sum_{\substack{\Gamma\ {\rm rooted\ cubic\ classes}\\{\rm loop\ number}\ n}}
\frac{1}{|\operatorname{Aut}\Gamma|}
=\frac{[t^{\,n-1}](1-t)^{-n-3}(1-t/2)^{-n}}{n(n+1)},\qquad n\ge1 .
$$

Here brackets extract the indicated power-series coefficient. The first six weights are $1/2,1,21/8,8,429/16,96$. An independent enumeration through four loops tests every connected subtree and counts the automorphisms of each decoration. These are counts before possible cancellations of color tensors or external polarizations. [Graph generator and checks](rooted_trees.py)

Every graph in this family is series–parallel: attach the insertion to every bag of the edge decomposition of the remaining tree. The resulting bags have at most three vertices, so the graph has no $K_4$ minor. Deleting the tree leaves first gives the corresponding series/parallel elimination directly: suppress the degree-two leaf with its spoke, then merge a parallel copy if one occurs. Subsequent elimination keeps the leaf's remaining neighbor and the insertion as the only possible endpoints. In the project ratio charts the currents, their Jacobian and all finite jet insertions are polynomials. Each barycentric monomial on a simplex integrates to the product of its exponent factorials divided by the factorial of the dimension plus total exponent. This gives a terminating rational algorithm for every individual kernel and jet. [Schwinger reduction](../../code/htk/spred.py)

All 15 four-loop kernels are nonzero when the output momenta cycle through the three unit basis vectors in generator vertex order. The eleventh graph is the first without a tight chain/bubble reduction. Its ordinary series–parallel chart gives $-37/19906560$ in the project's six-dimensional orientation. The replay checks each Jacobian independently with the grounded Laplacian. [Signed periods, graph by graph](rooted_periods.py)

At five loops all 49 rooted kernel polynomials are nonzero on output momenta, before Wick contractions. The table records one exact witness per graph, with the insertion momentum fixed by conservation and all unlisted vertex momenta zero. Of the 49 graphs, 43 admit tight chain moves and 6 use general series–parallel reductions. Graphs 5 and 23 vanish on their initial probes. Changing only the external momenta gives $23/5225472000$ and $-7033/80621568000$, respectively. Since all 49 kernels are rational polynomials in the seven ordered output momenta, their product is nonzero. Hence a common generic choice of seven momenta makes them all nonzero. The contribution of their signed sum to the five-loop differential also depends on color, odd variables and Wick signs. [Momenta and signed periods](five_loop_periods.jsonl) · [Replay with `rooted_periods.py 5 --jsonl`](rooted_periods.py)

The integrator contracts a simplex only with the polynomial factors that depend on its coordinates. Factors independent of that simplex pass through the integral. For monomial exponents $k_1,\ldots,k_m$ in its $m$ independent barycentric coordinates, contraction multiplies the remaining monomial by $k_1!\cdots k_m!/(m+k_1+\cdots+k_m)!$. It then combines terms with identical remaining exponents. This permits all tight local Jacobians and polynomial insertions to remain factored until a shared integration variable requires their product. [Implementation](../../code/htk/spred.py)

## Four-corner box ladder

Use the ladder's S4 edge order, sources $A,B$ at its left corners and $C,-A-B-C$ at its lower and upper right corners. Set $\Delta=\det(A,B,C)$. Remove the boxes from left to right, including the last one. In box $j$ let $\rho_j$ be its three-edge path resistance divided by the sum of that resistance and the parallel rung resistance. The ordered fractions $1>s_{1j}>s_{2j}>0$ are the path's last two tail resistances divided by its total resistance. Moving the two current sources to the next pair of corners sends $(A,B)$ to $(s_{1j}A+s_{2j}B,(1-s_{1j})A+(1-s_{2j})B)$. The wedge of the two sources gains a factor $s_{1j}-s_{2j}$.

Modulo the span of the left sources, the current on a replacement rung comes entirely from the right boundary. Expanding box $k$ sends the fraction $1-\rho_k$ of that current through its three-edge path. Thus the projected current of replacement rung $j$ is $-C\prod_{k>j}(1-\rho_k)$ in this quotient. The chain Jacobian at each step contributes $\rho_j^2$ times the determinant of the rung current and the two left sources. For the first $n-1$ boxes the path occupies the first three remaining edge slots. Removing its three ratios reverses the projective orientation, canceling the sign of the cotree-current Jacobian. Their projected rung current supplies one minus sign each. The final polygon has the project's positive box orientation. The complete signed density divided by $\Delta^n$ is therefore

$$
(-1)^{n-1}\prod_{j=1}^n
\rho_j^2(1-\rho_j)^{j-1}(s_{1j}-s_{2j})^{n-j}.
$$

The two elementary integrals give $2(j-1)!/(j+2)!$ and $1/[(n-j+1)(n-j+2)]$. Consequently

$$
\frac{\mathcal I}{\Delta^n}
=\frac{(-1)^{n-1}2^{n+1}}{(n+2)!(n!)^2((n+1)!)^2}.
$$

For $\Delta\ne0$, project a downward rung current $y$ to $\det(A,B,y)/\Delta$. The projection on the leftmost rung equals $-\prod_{j=1}^n(1-\rho_j)$, so this insertion divides the uninserted kernel by $-(n+3)!/(6n!)$. All joint projected rung moments follow from the same independent beta integrals. The replay compares the unintegrated density and every rung projection directly with exact Kirchhoff derivatives through eight loops. [Ladder and joint moments](ladder.py)

## Fixed gauge algebra at all loop orders

Take $H=\operatorname{diag}(1,-1)$, $E=E_{12}$ and $F=E_{21}$ in $\mathfrak{sl}_2$ with the trace pairing. Use $n$ input letters $\operatorname{tr}(Hb)$ followed by $\operatorname{tr}(Hc)$. Evaluate the output on $n+2$ independent odd ghost modes with colors $E,H,\ldots,H,F$ and momenta $A,t_1B,\ldots,t_nB,C$ in that order, where $A,B,C$ are the three coordinate covectors and $t_j=n+1-j$. At every loop order $n\ge1$, a linear combination of input jets of total degree $n^2(n+3)/2$ has a strictly positive coefficient in units $i/(8\pi^3)^n$.

In the color tree, any branch off the path from $E$ to $F$ with two or more external leaves contains a bracket $[H,H]=0$. Thus every vertex lies on that path. Both ends carry a spoke and an output, and adjacent spoke vertices would make a forbidden triangle. The $n+1$ spokes consequently alternate along the path in [the matrix-unit calculation](matrix_unit_channel.py). In its ordering the forward color is $\operatorname{tr}(H[E,(\operatorname{ad}H)^{2n-1}[F,H]])=-2^{2n+1}$, where $\operatorname{ad}H$ acts by $[H,-]$. Interchanging $E,F$ reverses this sign.

Let $p_i^A,p_i^B,p_i^C$ be the three momentum components absorbed by spoke $i$, ordered from $E$ to $F$ ($0\le i\le n$). They are the multipliers of the three holomorphic input derivatives. All endpoint loads $p_i^A,p_i^C$ are positive. For unit injection at the left endpoint, the rightward current across a link between successive marks is the sum of the $A$ loads to its right. For unit injection at the right endpoint, the leftward current is the sum of the $C$ loads to its left. Hence the $A$ potential decreases strictly along the path and the $C$ potential increases strictly. At each mark their ratio is $r_i=p_i^A/p_i^C$, since both loads use the same spoke resistance. Thus $r_0>\cdots>r_n$.

Use the following polynomial as a differential operator on the ordered input slots, with $p_i$ replaced by the three derivatives in slot $i$:
$$
\prod_{0\le i<j\le n}(p_i^A p_j^C-p_i^C p_j^A)
\left(\sum_{i=0}^n p_i^A p_i^B\prod_{k\ne i}p_k^C\right)^{n(n-1)/2}.
$$
Its first factor is alternating and strictly positive in forward path order. The second factor is symmetric in the input slots. Therefore summing the signed assignments of the inputs gives $(n+1)!$ times the ordered polynomial.

Give the boundary end of spoke $i$ potential $r_i$. Write $h_j$ for the harmonic potential at the interior output vertex $j$. To compare the two endpoints of any path link, remove the link. Each side has a positive effective conductance to its fixed boundary, with open-circuit potential a convex combination of its own boundary values. Every left boundary value exceeds every right boundary value. Reconnecting the link therefore gives a strictly positive rightward current. This proves $h_1>\cdots>h_n$. Symmetry of the inverse grounded Laplacian gives
$$
\sum_{i=0}^n r_i p_i^B=\sum_{j=1}^n t_j h_j.
$$
The endpoint loads and $h_j$ are independent of the permutation of the $t_j$. Put $d=n(n-1)/2$. In the signed sum of $(\sum_j t_{\sigma(j)}h_j)^d$, the surviving exponents must be distinct and hence are exactly $0,\ldots,n-1$. The two alternating determinants give
$$
\sum_{\sigma\in S_n}\operatorname{sgn}(\sigma)
\left(\sum_{j=1}^n t_{\sigma(j)}h_j\right)^d
=d!\prod_{1\le i<j\le n}(h_i-h_j),
$$
since $\prod_{i<j}(t_i-t_j)=\prod_{k=0}^{n-1} k!$. The ordered ladder density acquires $\prod_j t_j=n!$ relative to all middle momenta equal to $B$: in the box reduction its determinant at step $j$ is linear in the local $t_jB$ and only the $A$ component of the transported earlier sources contributes. Thus the signed middle-output sum has a strictly positive integrand after multiplying the displayed determinants by $(\prod_i p_i^C)^d$.

Reversing the path contributes $(-1)^n$ from the kernel, $-1$ from color, $(-1)^{n(n+1)/2}$ from input reversal and $(-1)^{(n+1)(n+2)/2}$ from output reversal. Their product is $+1$; the two endpoint assignments cancel the factor $1/2$ from graph symmetry. The forward color, ordered Wick sign and left-derivative sign multiply to $2^{2n+1}$. Integration gives a positive rational coefficient because every original current insertion integrates by finite simplex moments. With each input jet specified by its differential component in the $A,B,C$ coordinates, $(\partial_1 b,c)$ gives $2/3$ at one loop, $(\partial_1 b,\partial_2 b,c)$ gives $1/270$ at two loops and $(\partial_1 b,\partial_2 b,\partial_2^2 b,c)$ gives $499/29767500$ at three loops, with the input trace projections understood. The displayed alternating polynomial itself gives $4/3$ at one loop. [Signed graph integrals and harmonic-current checks](verify/sl2_jets.py)

The same nonvanishing on bare local jets holds for every nonabelian complex reductive gauge algebra with nondegenerate invariant pairing. A root $\mathfrak{sl}_2$ inside a simple factor is closed under all tree brackets and its restricted pairing is a nonzero multiple of the trace. For an abelian gauge algebra every cubic color tensor vanishes. These tests quantify over input jet order and probe the differential before passing to cohomology.

### Nonvanishing modulo total derivatives at fixed rank

In the fixed $\mathfrak{sl}_2$ witness above, multiply the polynomial of input jet variables by
$$
\prod_{j=0}^n p_j^A p_j^C .
$$
Here $p_j^A,p_j^C$ are the components of the current absorbed at spoke $j$ from the two endpoint sources. Both are strictly positive inside the time simplex. The product is symmetric in the input slots, so it preserves the sign of the Wick-summed integrand. The total input derivative degree is $n^2(n+3)/2+2(n+1)$.

For each factor $\operatorname{tr}(H\,\partial^k b)$ or $\operatorname{tr}(H\,\partial^k c)$ of this insertion, use $g_1$ in place of $H$ inside the trace, keeping every indicated derivative on $b$ or $c$. This polynomial is invariant under simultaneous conjugation. Evaluate $b=g_2=0$ and $g_1=H\exp(-P\cdot z/(n+1))$, where $H=\operatorname{diag}(1,-1)$ and $P$ is the total momentum of the $n+2$ ghost waves. Each wave has momentum along just one of the three coordinate covectors $A,B,C$. A jet differentiating in both the $A$ and $C$ coordinate directions vanishes on every such ghost wave. In every summand the original $c$ has these two derivatives. A nonzero graph must hit it, and must hit all $n$ original $b$ factors. This exhausts the $n+1$ functional derivatives. Extracting the coefficient of all $n+2$ distinct ghost parameters requires a ghost at every graph output. All original $g_1$ factors consequently remain spectators at the insertion.

The ghost waves and the $n+1$ spectators have total momentum $P-P=0$. Every total holomorphic derivative evaluates to zero on this coefficient. Thus $Q_n$ induces a nonzero operator on conjugation-invariant local polynomials modulo total holomorphic derivatives for fixed $\mathfrak{sl}_2$ at every $n\ge1$. At one loop the integrated coefficient is $1/45$ in units $i/(8\pi^3)$. The statement is at the cochain level.
