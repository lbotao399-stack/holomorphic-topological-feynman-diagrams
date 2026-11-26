# Commuting-pair scalar operators at every rank and degree

For gauge algebra $\mathfrak{sl}_N$, every polynomial invariant of two commuting traceless matrices $x=g_1$, $y=g_2$ lifts to local quantum cohomology on $\mathbb C^{3|2}$. At bidegree $(p,q)$ the lift terminates by $\hbar^{\min(p,q)}$. Constant commuting matrices detect the same number of independent integrated classes for translation-invariant local densities. The result concerns ordinary functions on the reduced commuting-pair quotient; it counts the restrictions of classes to scalar vacua.

## Marked product in three complex dimensions

Take the marked product of two insertions $A,B$, with its auxiliary edge joining the two **body** points. The auxiliary free pair has complementary Dolbeault shifts summing to $3-1$; its normalized single contraction gives $AB$ at tree level. The dressed quadratic relation gives
$$
\mathbf Q(A\star B)=(\mathbf QA)\star B+(-1)^{|A|}A\star(\mathbf QB).
$$
In the collision identity a subgraph containing neither marked insertion contributes the bulk master equation, one containing a single insertion gives the corresponding term on the right, and both together give the left. The auxiliary pair has no interaction vertex. In this theory each pure bulk loop already vanishes by its odd delta functions, and the bulk tree equation is Jacobi. The proof uses the three-dimensional propagator and graph count in this identity. [Marked contraction and quadratic relations](https://arxiv.org/abs/2403.13049) · [Project's dressed-product construction](../../../../sources/LLC_main.tex) · [Six-dimensional bulk equation](../../README.md)

Use insertions involving only scalar and ghost jets. At loop order $m>0$ a binary graph has $2m$ interaction vertices and $3m+1$ edges, including the mark. Deleting the insertions leaves a forest; a cycle among interaction vertices gives zero odd density. For a component of $v$ interaction vertices with $k$ spokes, the subgraph with both insertions and the mark has $v+2$ vertices and $v+k$ edges. Sparsity gives $2k\le v+2$. Summing over all components gives equality by tightness of the full graph. Thus $v=2k-2$ in every component, and $k\ge2$.

All $k$ spokes differentiate at one common odd coordinate. Each scalar supplies one odd coordinate and a ghost supplies two. Consequently $k=2$, $v=2$, and both spokes hit scalar letters. Each of the two vertices has a spoke since a leaf of the full graph violates sparsity. If both spokes end at the same insertion they form a forbidden triangle. Hence each component is a three-edge path from $A$ to $B$. Each cubic vertex on the path outputs a ghost, as the two scalar derivatives already saturate $d^2\theta$. The path uses one $x$ and one $y$, with signs $+,-$ in root order $(x,y),(y,x)$. Graphs made of $m$ such paths and the mark are tight: each path adds three edges and two vertices, and each partially included path obeys the subset bound, with equality for a complete path between the two insertions. Equivalently each path is a tight chain expansion of the marked edge. For a single flavor polarization the two-derivative odd factor squares to zero. Its ordinary and marked products agree, including for jets. [Independent graph enumeration](check.py)

In an iterated marked product, each loop consumes two scalars of opposite flavor and adds two ghosts. Ghost jets produced at a previous multiplication remain spectators thereafter. This proves the claimed termination for homogeneous bidegree. The inverse automorphism weight is $1/m!$ when ordering the paths with both roots fixed.

## Exact kernel and derivative insertions

Number the marked edge first, then the three edges on each path, all oriented from $A$ to $B$. Number the vertices $A$, the two vertices on each successive path, then $B$ as root. Let $p_i,q_i$ be the outgoing ghost momenta at the first and second vertex of path $i$, and $\lambda_A$ the total momentum of spectators at $A$. The dimensionless times on each path are $(1-a_i,a_i-b_i,b_i)$ with $0<b_i<a_i<1$. The conductance fractions $w_0,w_1,\ldots,w_m$ lie on the unit simplex, with 0 the mark. Kirchhoff's equations give
$$
L=\lambda_A+\sum_{i\ {\rm paths}}(a_i p_i+b_iq_i),\qquad
Y_i=w_iL-a_i p_i-b_iq_i .
$$
The marked current is $w_0L$. The three path currents are $Y_i,Y_i+p_i,Y_i+p_i+q_i$. The cotree consists of the first edge on each path. In its current Jacobian the three columns per path $j$ on row block $i$ are
$$
\delta_{ij}L,\qquad (w_i-\delta_{ij})p_j,\qquad(w_i-\delta_{ij})q_j.
$$
The scalar matrix $H_{ij}=w_i-\delta_{ij}$ has determinant $(-1)^m w_0$. Factoring $H$ in three component directions and expanding along the now block-diagonal $p_i,q_i$ columns leaves one factor of $\det H^{-1}$. The determinant is $w_0^2\prod_i\det(L,p_i,q_i)$. For the stated edge order, the spanning-tree incidence determinant and the tree/cotree permutation sign are both positive. Fixing the marked time to one gives path resistances $R_i=w_0/w_i$, whose Jacobian has sign $(-1)^m$. The S2 kernel is therefore
$$
I^{\mathrm{S2}}_m=(-1)^m\frac{2}{(m+2)!}\int_{0<b_i<a_i<1}
\prod_{i=1}^m\det(L,p_i,q_i)\,\prod_{i=1}^m da_i\,db_i .
$$
For input jets retain the $w$ integrals and multiply by the indicated powers of the endpoint currents, with a minus sign for a derivative at $A$ and a plus sign at $B$. For integer exponents $j,k,d_i\ge0$ the necessary moments are
$$
\int_{0<b<a<1} a^j b^k\,da\,db=\frac{1}{(k+1)(j+k+2)},\qquad
\int_{\Delta^m}w_0^{2+d_0}\prod_{i=1}^m w_i^{d_i}\,dw
=\frac{(2+d_0)!\prod_{i=1}^m d_i!}{(m+2+\sum_{i=0}^m d_i)!}.
$$
The replay processes paths sequentially on a polynomial in three components of $L$. A shift by $ap_i+bq_i$ is evaluated by its finite directional Taylor polynomial. Edge powers keep track of their three accumulated $L$ degrees, so increasing the loop order does not introduce additional momentum coordinates. For all $p_i,q_i$ in the same plane the undifferentiated kernel reduces to $(-1)^m 2^{1-m}\prod_i\det(\lambda_A,p_i,q_i)/(m+2)!$. [Kernel evaluator](kernel.py) · [Comparison against the independent project engine](check.py)

## Physical coefficient

For ordered paths contract the two input color derivatives of each path with the ordered adjoint chain $\operatorname{ad}_{c(p_i)}\operatorname{ad}_{c(q_i)}$ using the action pairing. Use alternating flavor as above. The order-$m$ contribution carries
$$
\frac1{m!}\left(-\frac{\hbar}{8\pi^3}\right)^m I^{\mathrm{S2}}_m
$$
multiplying this ordered contraction and the remaining fields. Multiple functional derivatives count assignments to distinct letters. The common coefficient can be fixed independently at every order by the chain identity, using the project's $Q_0=i\widehat Q_0$ and $Q_1=i\widehat Q_{\rm box}/(16\pi^3)$. [Heat normalization](../../VERIFY.md)

To see this, use $m$ independent color algebras and take at the two roots the products of the quadratic $x$ and $y$ traces for all factors. Each root is a strict cocycle. Select the coefficient with one residual scalar of each type per factor. Put three ghost waves in the first factor with the three unit coordinate momenta, and two waves with momenta $\tau e_1,\tau e_2$ in every other factor. All scalar momenta are zero. The color chains in separate factors force one pair per path, canceling the ordered $m!$. Extract the $\tau^{2m-2}$ coefficient. Every path other than the first supplies its two explicit powers of $\tau$. Write $j=1,2,3$ for the position of the ghost carrying $e_3$ in the first factor's ordered triple. Alternating its three odd polarizations supplies the same sign as the momentum determinant.

On $\widehat Q_0$ of the order-$m$ term, the remaining determinant on the first path vanishes to leading order unless varying the residual $x$ injects the first ghost momentum into $\lambda_A$. Varying a ghost or the right scalar gives only higher powers of $\tau$. The cyclic trace moves the first bracket $[c,x]$ to the other side with a minus sign. The other path determinants then contribute the $(m-1)$st power of respectively $1,a_1,b_1$ according to $j$. Their three moment coefficients on the first triangle are
$$
J_1=\tfrac12,\qquad J_2=\frac1{m+1},\qquad J_3=\frac1{m(m+1)} .
$$
On $\widehat Q_{\rm box}$ of the order-$(m-1)$ term, all its paths belong to the other factors. The new box has positions $1>s_1>s_2>s_3>0$. The momentum of the consumed left scalar has third component $s_j$ in one orientation. Exchanging the root ends and reversing the ghost order gives an equal term by $s_i\mapsto1-s_{4-i}$. The resulting simplex moments of $s_j^{m-1}$ are $J_j/(m+2)$ for all three values of $j$. The ratio of the signed cotree simplex prefactors at orders $m-1,m$ is $-(m+2)$, and the sign from moving the first bracket cancels this minus. The two root assignments therefore give $\widehat Q_{\rm box} I_{m-1}/(\widehat Q_0 I_m)=2$ for the contracted probes. Choosing the three $\mathfrak {sl}_3$ root colors in the [mixed quartic](../mixed.md) gives a nonzero trace; the other factors can use $E_{12},E_{21}$ with a noncentral diagonal color. This fixes the common weight recursively from ordinary multiplication. The unrestricted two-loop check also varies the external scalar momenta and retains every ghost and scalar variation. [Ward replay](check_normalization.py)

## Scalar vacua and counting

A commuting complex pair is simultaneously triangularizable. An invariant polynomial takes its value on the diagonal pair by a diagonal conjugation limit. On the diagonal it is symmetric under permutations of the $N$ eigenvalue pairs. Conversely, any symmetric polynomial is a combination of products of power sums $\sum_i x_i^a y_i^b$: express monomial orbit sums as sums on distinct labels and remove label equalities by inclusion-exclusion. Each power sum is the $s^b$ coefficient of $\operatorname{tr}(x+sy)^{a+b}$ divided by $\binom{a+b}{b}$. The Newton identities reduce powers higher than $N$ using the characteristic polynomial of $x+sy$, whose linear trace vanishes.

Choose a homogeneous linear basis of these invariant polynomials and express each basis vector in the polarized traces of degrees $2$ through $N$. Replace each ordinary multiplication by the marked product, with a fixed parenthesization. The traces are strict cocycles; the marked chain identity proves closure of the resulting finite sum. Restriction to constant commuting pairs annihilates every quantum boundary: the commutators vanish classically, and every graph at positive order has a vanishing current Jacobian at constant external momenta. It also kills each positive-loop correction to the lifted representative, which contains ghost fields. The same restriction kills total derivatives of translation-invariant densities. Thus evaluation splits on the chosen basis as a vector-space map; using other polynomial presentations can give lifts differing by classes of zero scalar restriction.

The coefficient counting these independent restrictions is
$$
(1-s)(1-t)[\nu^N]\prod_{a,b\ge0}\frac1{1-\nu s^a t^b}.
$$
The product counts unordered multisets of $N$ eigenvalue monomials; subtracting the two scalar means removes two polynomial generators of weights $(1,0),(0,1)$. On the diagonal weights $(2,2),(3,3),(4,4)$ the dimensions for $N=3$ are $2,3,5$, and for $N=4$ they are $3,6,13$. An independent permutation-character average checks every bidegree through $(7,7)$ at ranks 2 through 7. [Counting replay](hilbert.py)

After inverting the discriminant of $x$, the construction also lifts all Weyl-invariant functions of finite holomorphic scalar jets. The alternating scalar one-jet has a global polynomial completion already before localization. [Regular scalar jets and integrated example](regular_jets.md)
