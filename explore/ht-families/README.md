# Convergence and the power of $\pi$ for series–parallel kernels

Take zero topological source, $T=1$, and even loop number $n$. Every connected graph with $N-1=n(H+1)$, where $N$ is the number of edges and $H$ the holomorphic dimension, has an absolutely integrable projective top form. Multiplication by a polynomial in the edge currents preserves absolute convergence. At zero displacement, every series–parallel reduction of such a graph gives

$$
K_\Gamma(\lambda;0)\in \pi^{n/2+2E}\,\mathbb Q[\lambda].
$$

Here $E$ adds, over the series nodes of the reduction, the integer part of half the number of children of odd loop number. Polynomial current insertions with rational coefficients obey the same statement. This removes the face-integrability hypothesis in c2 of the [earlier verification](VERIFY.md). For $H=1$, every $E>0$ kernel vanishes identically, including its displaced and polynomial-current insertions. The proof and exact checks cover all six previously open six-loop graphs. In two holomorphic dimensions a four-loop graph of the same type instead contributes $-\pi^4/26880$ at the specified sources. [Four-loop counterexample and independent replays](h2_odd_triangle.md) [Odd-loop triangle cancellation](odd_triangle.md) · [Reduction independence of $E$](positive_period.md) [Kernel conventions](../../code/htk/CONVENTIONS.md) · [Exact replay](verify/sector_regularity.py)

Subdividing a cotree relates the kernels in successive holomorphic dimensions by an exact signed sum. This gives nonzero four-loop $E=1$ kernels for every $H\ge2$ at $T=1$ and irrational three-loop $K_4$ periods for every $H\ge3$ at $T=0$. [Dimension-raising identity](dimension_raise.md) · [Signed replay](verify/dimension_raise.py)

At two loops the density pushed forward to the $2H$ independent currents is the signed hyperbolic area of a projective polygon. Every segment gives one linear inequality on the three path times. For $(1,H+1,H+1)$ four endpoint evaluations of an arctangent give that density with arbitrary current insertions and displacement. [Calculation and pointwise replay for the five- and eleven-dimensional theta graphs](current_polygon.md)

## 1. Extension to the boundary of a Schwinger sector

Order the edge times $t_1\ge\cdots\ge t_N>0$, reorder the columns of the cycle matrix $C$ accordingly, and fix the scale $t_1=1$. In this proof the edge subscripts enumerate this order. Use $s_j=\sqrt{t_{j+1}/t_j}$, so $\sqrt{t_e}=s_1\cdots s_{e-1}$. The ordered sector is the open cube in these variables.

The project period matrix is $M=C\,\mathrm{diag}(t)\,C^{\mathsf T}$. The columns of $F=\mathrm{diag}(\sqrt t)\,C^{\mathsf T}$ form a frame with Gram matrix $M$. Orthogonal projection of the ordinary derivative to their span has connection matrix $M^{-1}F^{\mathsf T}dF=\tfrac12M^{-1}dM$. Its curvature $R$ therefore obeys

$$
R=d\big(\tfrac12 M^{-1}dM\big)
   +\tfrac14M^{-1}dM\wedge M^{-1}dM
   =-\tfrac14 M^{-1}dM\wedge M^{-1}dM .
$$

For the project's original form this gives

$$
\rho_n(M)=2^{n/2}\,
\frac{\mathrm{Pf}(MR)}{\sqrt{\det M}}.
$$

Choose the first $n$ independent rows of $C^{\mathsf T}$ by greedy scanning, with row set $I$. The normalized frame $F F_I^{-1}$ has the identity in rows $I$. Its entry in row $e$, in the column indexed by $i\in I$, is

$$
\big[C^{\mathsf T}((C^{\mathsf T})_I)^{-1}\big]_{ei}
\,\sqrt{t_e/t_i}.
$$

For $e<i$ the constant coefficient is zero. Indeed, a rejected earlier row is a linear combination of the already selected rows. For $e\ge i$ the square root is the monomial $s_i\cdots s_{e-1}$. Thus the frame extends polynomially to the closed cube, retaining full rank there. Its Gram inverse, projected connection and curvature are smooth on the closed cube. Under a frame change by $B$, the lowered curvature becomes $B^{\mathsf T}MRB$, so the quotient of its Pfaffian by the square root of its Gram determinant changes by $\mathrm{sgn}(\det B)$. For $B=F_I^{-1}$ this sign is constant on the open sector. This proves smooth extension of $\rho_n$.

For holomorphic currents write $y=y^0+C^{\mathsf T}\ell$ in each component, with any fixed source flow $y^0$. The loop equation is $M\ell=-C\,\mathrm{diag}(t)y^0$. Cramer's rule and Cauchy–Binet express $\ell$ as the weighted average of the solutions $\ell_J$ to $y^0_J+(C^{\mathsf T})_J\ell_J=0$, over all row bases $J$, with weights

$$
\det((C^{\mathsf T})_J)^2\,\prod_{e\in J}t_e .
$$

Indeed, replacing one column of $C^{\mathsf T}$ by $-y^0$ in the second factor of $C\,\mathrm{diag}(t)\,C^{\mathsf T}$ gives the corresponding Cramer numerator. For every prefix of the edge order, the number of its rows in $I$ is the rank of that prefix. Thus if $I=(i_1,\ldots,i_n)$ and $J=(j_1,\ldots,j_n)$ are increasing, $j_k\ge i_k$ for each $k$. Dividing all the weights by $\prod_{i\in I}t_i$ makes them polynomials in the $s_j^2$. Their sum is at least $\det((C^{\mathsf T})_I)^2>0$ on the entire closed cube. All currents and their differentials extend smoothly. The source form obeys $\psi=-M\,d\ell$ in each holomorphic component by differentiation of the loop equation. The factor $U^{-H}$, where $U=\det M$, cancels against the $H$ exterior determinants of $M$.

The projective integrand is therefore a smooth multiple of $\rho_n$ wedged with the current differentials in each closed cube. Its absolute integral is finite there. The finitely many ordered cubes cover the positive projective domain up to measure zero. The same argument applies with a polynomial in the currents or the factor $\exp(-\sum_{e\in\mathrm{edges}}y_e\cdot z_e)$ at fixed displacement.

## 2. Integrable Dirichlet expansion

For a two-terminal subnetwork, $U$ is the spanning-tree complement polynomial and $R$ in this section is its effective resistance. Series gluing gives $U_{\rm S}=\prod_i U_i$ and $R_{\rm S}=\sum_{i\in\mathrm{children}}R_i$. Parallel gluing uses one child spanning tree and terminal-separating two-forests in all other children. The two-forest complement polynomial is $R_iU_i$ by the matrix-tree identity for effective resistance. Hence $U_{\rm P}=(\prod_i R_iU_i)/R_{\rm P}$. Writing $L_i$ for a child's loop number, induction on the reduction tree gives the exponents in $U$ after fixing the root resistance to one:

| Child ratio $b$ | Exponent of $b$ in $U$ |
|---|---:|
| $R_i/R_{\rm S}$ at a series node | $L_i$ |
| $R_{\rm P}/R_i$ at a parallel node | $-1-L_i$ |

In fact $U_i$ has the homogeneous scale factor $R_i^{L_i}$, and each parallel gluing with $k$ children contributes $k-1$ further loops. The edge times are Laurent monomials in the ratios. The Kirchhoff inverse and the Pfaffian in the source formula have only powers of $U$ in their denominators apart from $U^{-1/2}$. The oriented ratio Jacobian has fixed sign in the interior and is Laurent polynomial there. Thus the entire density times $\sqrt U$ has a Laurent-polynomial representative on the product of ratio simplices. The same holds after any polynomial current insertion.

In each block of barycentric coordinates write this density as a monomial with half-integer powers times an ordinary polynomial. Make the polynomial homogeneous in that block by multiplying each lower-degree monomial by the needed power of the sum of the block coordinates, which equals one. Do this in all blocks. If the resulting polynomial vanishes identically the period is zero. Otherwise divide out every barycentric coordinate that is a common factor, transferring it to the monomial. On the relative interior of each codimension-one face the residual polynomial is generically nonzero. Indeed, homogeneity in each block extends vanishing there by independent positive scaling to an open subset of the corresponding coordinate hyperplane, forcing divisibility. Omit one-child blocks whose sole coordinate equals one.

With the remaining powers denoted by $\alpha_b$, the absolute integral on an interior patch of the face $b=0$ is bounded below by a positive constant times $\int_0^\varepsilon b^{\alpha_b}\,db$. Section 1 therefore gives $\alpha_b>-1$ for every $b$. For symbolic sources one may choose generic real sources so that none of the finitely many nonzero face polynomials becomes identically zero. The exponent bound applies coefficientwise. Every monomial in the resulting polynomial times these powers has an absolutely convergent product of Dirichlet integrals.

If $h$ exponents of a block are half-integral and nonintegral, its Dirichlet numerator contributes $(\sqrt\pi)^h$ and its denominator contributes $\sqrt\pi$ exactly when $h$ is odd. The block contributes a rational multiple of $\pi^{\lfloor h/2\rfloor}$. Its $h$ counts odd $L_i$ at a series node and even $L_i$ at a parallel node. To sum these exponents, write $\phi$ for the total from all blocks of a subtree and $L$ for its loop number. At a series node,
$$
\phi-\lceil L/2\rceil
=\sum_{i\in\mathrm{children}}(\phi_i-\lceil L_i/2\rceil)
 +2\lfloor h/2\rfloor .
$$
At a parallel node the same equality holds with the last term zero, since $L=\sum_{i\in\mathrm{children}}L_i+k-1$. An edge has $\phi=L=0$. At the root this gives $\phi=n/2+2E$ for even $n$. Two allowed reductions with different $E$ force a vanishing kernel, since a nonzero integral power of $\pi$ is irrational.

## Raising the loop order at fixed holomorphic dimension

Consider a pure holomorphic graph at zero displacement. Replace its oriented edge $e:u\to w$ by that same edge in parallel with a new path $u\to x_1\to\cdots\to x_{H-1}\to w$. Keep the original edge slots and append the new path in order. The child has one more loop, $H$ more edges, and $H-1$ more vertices. It remains admissible. In a connected subgraph, prune incomplete portions of the new path. If the complete path and the direct edge both occur, deleting the path removes one loop and $H$ edges. If just the complete path occurs, replace it by the direct edge, decreasing the edge count by $H-1$ at fixed loop number. The parent bound $|E_S|\ge Hn_S+1$ applies in each case.

Choose integer vectors $P_1,\ldots,P_{H-1}$ for the new sources $sP_j$ and subtract $s(P_1+\cdots+P_{H-1})$ from the original root source. Let $r$ be the effective time of the edge, $\rho$ its direct current fraction, and $1>\sigma_1>\cdots>\sigma_{H-1}>0$ the tail fractions of the path. The direct time and total path time are $r/\rho$ and $r/(1-\rho)$. In a projective chart fixing another edge, the absolute time Jacobian is $r^H/[\rho^2(1-\rho)^{H+1}]$. If $y_e$ denotes the reduced edge current, the direct current is $\rho y_e$ and the first path current is $(1-\rho)y_e-s\sum_j\sigma_jP_j$. The reduced endpoint sources acquire $s\sum_j\sigma_jP_j$ at $u$ and $s\sum_j(1-\sigma_j)P_j$ at $w$, in addition to the root adjustment.

Choose a parent spanning tree excluding $e$, and extend it by all but one of the newly appended edges. In the child cotree rows replace the two local currents by their sum and the direct current. The sum differs from $y_e-s\sum_j\sigma_jP_j$ by a constant source vector. Subtract $\rho$ times its differential from the direct-current differential. The remaining $H$-component row is
$$
y_e\,d\rho+\rho s\sum_{j=1}^{H-1}P_j\,d\sigma_j .
$$
Its exterior determinant is $\rho^{H-1}s^{H-1}\det(y_e,P_1,\ldots,P_{H-1})$ in the $d\rho,d\sigma$ coordinates. Wedging the other rows gives the original current Jacobian with the shifted sources held fixed. The sign $\varepsilon_e$ is the orientation of these cotree and coordinate changes, independent of times and sources. It can be evaluated as the exact signed ratio at any nonsingular positive rational time. [Original oriented current convention](../../code/htk/CONVENTIONS.md)

Write $K_{\Gamma_e}(s)$ for the child's bare period at the prescribed sources. The coefficient of $s^{H-1}$ uses the unshifted parent kernel. The $\rho$ integral is $1/H$ and the tail simplex has volume $1/(H-1)!$. Hence
$$
\varepsilon_e[s^{H-1}]K_{\Gamma_e}(s)
=\frac{1}{H!}\int D_\Gamma(t,\lambda)\det(y_e(t,\lambda),P_1,\ldots,P_{H-1})\,dt .
$$
All integrals are absolutely convergent, including the current insertion. The displayed coefficient follows also by polynomial evaluation in $s$ at fixed ratios and then integration; the finitely many parent moments are absolutely integrable. For the incidence $B_{ev}$, positive on outgoing edges at $v$, Kirchhoff conservation now gives
$$
H!\sum_{e\,\mathrm{incident\ at}\ v}B_{ev}\varepsilon_e[s^{H-1}]K_{\Gamma_e}(s)
=\det(\lambda_v,P_1,\ldots,P_{H-1})\,K_\Gamma(\lambda;0).
$$
Choose $v$ with nonzero source and coordinate vectors $P_j$ giving a nonzero determinant. Each child is polynomial in $s$ of degree at most $H(n+1)$, where $n$ is the parent loop order. Integer evaluations at $s=0,\ldots,H(n+1)$ recover the coefficient by rational interpolation. Every parent value at integral sources is therefore a rational linear combination of bare child values at integral sources. The zero value requires no condition on $v$.

The replay checks the pointwise identity at three rational points per edge in dimensions two, three and four. For the seven-edge three-loop graph with the [nonzero five-source fixture](../nonSP/README.md), graft on edge $01$ with $P_1=(0,1)$. Nine direct four-loop integrals give
$$
[s]K_{\Gamma_{01}}(s)=-\frac{222317}{50400}+\frac{3021}{1400}\zeta(3)
=-\frac12\int D_\Gamma\,y_{01}^{1}\,dt .
$$
[Pointwise and integrated replay](verify/loop_raise.py)

### First irrational bare period

| Holomorphic dimension | Loop orders at which some integral-source bare kernel is irrational |
|---|---|
| $H=1$ | None. The two-vertex graph consists of parallel edges. |
| $H=2$ | Every $n\ge4$. The sole non-series–parallel graph at three loops has a rational six-degree source polynomial. |
| $H\ge3$ | Every $n\ge3$. The one- and two-loop graphs are series–parallel. |

In the last row the three-dimensional base value is $(39-32\zeta(3))/67200$. Cotree subdivision raises the dimension, preserving the existence of an irrational value, and the polygon construction above raises the loop number. In two dimensions the four-loop base can be taken as $-789/640+59\zeta(3)/80$. More generally, for **every $H\ge2$ and $n\ge4$**, the field generated by all bare geometric graph periods at integral sources contains $\mathbb Q(\zeta(3),\zeta(5))$: the entire four-loop, two-dimensional field embeds by the two rational reconstruction identities. At three loops in every $H\ge3$ the field equals $\mathbb Q(\zeta(3))$ by the all-current upper bound. [Base fields and witnesses](../nonSP/verify/prism_forcer/README.md) · [Three-loop upper bound](../nonSP/verify/k4_periods/README.md) · [Dimension-raising identity](dimension_raise.md)

### Eliminating current insertions by adding loops

Retain a product of current components on the original edge slots of the child. If $k$ of those factors use the direct copy of $e$, their product acquires exactly $\rho^k$ at $s=0$; every other retained factor is its parent current. The coefficient identity for the child with these insertions consequently gives the parent integral with the same product and one additional factor $\det(y_e,P_1,\ldots,P_{H-1})$, divided by
$$
(H+k)(H-1)! .
$$
The sign is the same $\varepsilon_e$. Taking the $P_j$ to be coordinate vectors selects any desired component of $y_e$ up to sign. Rational interpolation and iteration express **every $m$-current period at $n$ loops as a rational linear combination of bare $(n+m)$-loop periods at the same $H\ge2$**, with integral sources when the parent sources are integral. For repeated insertions the extra factor $k$ in the denominator is essential. The same-edge replay gives
$$
[s]\!\int D_{\Gamma_{01}}\,y_{01}^1\,dt
=\frac{4003781}{907200}-\frac{3772}{1575}\zeta(3)
=-\frac13\int D_\Gamma (y_{01}^1)^2\,dt .
$$

At fixed loop order the integrated vertex equation also expresses $\lambda_v^i$ times any $m$-current moment as the signed sum of moments with one additional component on the edges incident at $v$. Taking $\lambda_v^i\ne0$ shows that the rational span of order-$m$ periods is contained in that of order $m+1$. Therefore every irrational bare case in the table has an irrational moment at **each** current order. At $H=2,n=3$ the bare kernel is rational but its first edge derivative can equal $17/9450-\zeta(3)/700$; an irrational moment exists at every $m\ge1$. At $H=1$ and at loop orders one and two all finite current moments are rational by series–parallel reduction. Together these alternatives decide existence of irrational geometric moments for every positive loop order, holomorphic dimension and finite insertion order. [Three-loop derivative and original graph](../nonSP/README.md) · [Repeated-insertion replay](verify/loop_raise.py)
