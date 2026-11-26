# Three-loop $K_4$ subdivisions with arbitrary current insertions

For every admissible subdivision of $K_4$ on $\mathbb C^H$, every coefficient of the displacement expansion at rational conserved sources lies in $\mathbb Q+\mathbb Q\zeta(3)$. In [kernel.py](kernel.py), `evaluate(edges, source, insertions)` returns these two rational coefficients. An insertion `(e,i)` multiplies by component $i$ of the current on the oriented edge $e$. Repeated entries compute higher derivatives with the sign from $-\partial_z$.

Use the integer momenta $((1,0,0),(0,1,0),(0,0,1),(1,1,-1),(1,-1,1),(-1,1,1),(-1,-1,0),(-1,-1,-2))$ in vertex order. Each of the six $\mathbb C^3$ graphs has an irrational bare period by Apéry's theorem.

| Subdivision in project order | $1$ | $\zeta(3)$ |
|---|---:|---:|
| $[1,1,1,1,3,3]$ | $-20383/1209600$ | $83/4200$ |
| $[1,1,1,2,2,3]$ a | $1333/57600$ | $-23/1200$ |
| $[1,1,1,2,2,3]$ b | $-7495/217728$ | $5/168$ |
| $[1,1,1,2,2,3]$ c | $394519/10886400$ | $-307/8400$ |
| $[1,1,2,2,2,2]$ a | $8663/259200$ | $-29/1200$ |
| $[1,1,2,2,2,2]$ b | $-409/226800$ | $11/8400$ |

## Schwinger numerator and convergence

For a holomorphic graph with $n$ loops, $H$ current components and $Hn+1$ edges, write $U=\det(C\,\mathrm{diag}(t)C^{\mathsf T})$ for the Kirchhoff polynomial of its cycle matrix $C$. The spanning-tree formula expresses $Uy_e^i$ as the sum of each cotree monomial times the component of its tree flow. Thus every current stays bounded on positive times. With $D$ the oriented projective density, use the cofactor at edge $j$:

$$
U^{n+H}D=
\frac{(-1)^j\det\big(C_{re}Uy_e^i\big)_{(i,r),\,e\ne j}}
{t_j U^{n(H-1)}}.
$$

Here $r$ enumerates the cycle rows. To see polynomiality, take an irreducible factor of $U$ and a time on which it depends. Over the field of the remaining times, $U$ is affine with a simple root. At that root its cycle matrix $C\,\mathrm{diag}(t)C^{\mathsf T}$ has corank one, since its determinant has a simple zero. In the cycle solution for the current, its adjugate consequently makes all $H$ edge arrays $(Uy_e^i)_e$ proportional there. Subtract multiples of one block from the other $H-1$ blocks of $n$ rows. The determinant has at least $n(H-1)$ factors of the irreducible polynomial; if every block vanishes its multiplicity is larger. The multiaffine $U$ is square-free. Changing the pivot removes the apparent $t_j$ denominator. Homogeneity gives degree $n^2-1$, and using each varied time as pivot at positive infinity bounds its individual degree by $n-1$. For $K_4$ this gives degree eight and degree at most two per time. Newton differences on the downward-closed degree-eight grid recover the polynomial. Every run checks an unequal-time evaluation against an independent grounded Kirchhoff solution.

The density is absolutely integrable by the [signed current coordinates](../../README.md#signed-current-coordinates). Where its Jacobian is nonzero the positive times are recovered uniquely from the maximal minors. Its absolute integral is the volume of bounded current regions. Multiplication by finitely many bounded currents preserves absolute convergence. For a nonempty proper subset of times scaled to zero, a nonzero lowest-degree numerator coefficient is nonzero on an open cone of ratios with the other times positive. In this cone $U$ has scaling degree equal to the cycle rank of the subset. Absolute convergence therefore forces the strict scaling inequality for each numerator monomial, and the ordered Schwinger sectors give its sufficiency. Equivalently, multiplying by a tree-current numerator increases the degree on each subset by at least its cycle rank, which exactly covers the added denominator $U$.

## Reduction to two bubble products

Let $m$ count the inserted currents and put $q=H+3+m$. On each of the six core paths sum its Schwinger times. The tree polynomial of a subdivision becomes exactly the tree polynomial of $K_4$ in these six sums. For an original monomial $t^k$, a path carries the core power $\alpha=\sum_{e\in\mathrm{path}}(k_e+1)$. Its simplex integral contributes $\prod_{e\in\mathrm{path}}k_e!/(\alpha-1)!$. The projective cut at core edge $01$ then gives the two-loop kite at $Q^2=1$ in dimension $2q$, multiplied by $\prod_e k_e!/(q-1)!$. The six core powers sum to $3q$. Convergence gives $1\le\alpha<q$ on every edge, by scaling its complement. The solver verifies a second cut on the opposite edge. [Cut normalization](../k33_forcer/README.md)

Continue the kite to dimension $d=2q-2\epsilon$ with fixed positive propagator powers. For powers $a,b,c,d',e$ of $k^2,(Q-k)^2,l^2,(Q-l)^2,(k-l)^2$ respectively, integrate the divergence with vector $k-l$ in the $k$ loop. It gives the recurrence
$$
\begin{aligned}
(d-2e-a-b)I(a,b,c,d',e)
={}&a\big[I(a+1,b,c,d',e-1)-I(a+1,b,c-1,d',e)\big]\\
&+b\big[I(a,b+1,c,d',e-1)-I(a,b+1,c,d'-1,e)\big].
\end{aligned}
$$
It terminates since $c+d'+e$ decreases at each step. When one of these powers is zero the remaining bubbles give, up to rational functions of $d$, $G(1,1)^2$ and $G(1,1)G(1,2-d/2)$, with $G(a,b)=\Gamma(a+b-d/2)\Gamma(d/2-a)\Gamma(d/2-b)/[\Gamma(a)\Gamma(b)\Gamma(d-a-b)]$. The native reduction keeps both channels over a common integer-polynomial denominator in $d/2$, removes its common divisor after each step, and exchanges loop momenta to minimize the descending exponent sum. [Integer-polynomial recurrence](kite_native.py) · [Independent symbolic recurrence](../../../../deepen/o3/higher_orders/four_loop_reduction/doubled_kite.py)

Extract the common regular factor $\Gamma(1+\epsilon)^2\Gamma(1-\epsilon)^4/\Gamma(1-2\epsilon)^2$, equal to one at zero. Gamma recursion leaves a rational function plus a rational function times
$$
\frac{\Gamma(1+2\epsilon)\Gamma(1-2\epsilon)^2}
{\Gamma(1+\epsilon)^2\Gamma(1-\epsilon)\Gamma(1-3\epsilon)}
=1-6\zeta(3)\epsilon^3+O(\epsilon^4).
$$
The cubic coefficient follows from $((-2)^3+2\,2^3-2(-1)^3-1-3^3)/3=-6$. If the multiplier had a pole of order $s>3$, the coefficient at $\epsilon^{3-s}$ would be rational plus a nonzero rational multiple of $\zeta(3)$. Apéry's theorem and convergence exclude this. The finite part therefore uses only the constant of the sum of the two rational functions and $-6\zeta(3)$ times the cubic-pole coefficient of the second. This proves the stated period space for every finite list of currents. [Exact Laurent evaluation](core.py)

For fixed real sources the Taylor series in the complex displacements converges everywhere. The bounded tree currents dominate the exponential on every compact set of displacements by a constant times $|D|$.

## All holomorphic three-loop graphs

An admissible graph with at least three vertices has no cut vertex. Splitting there into two edge-containing blocks would give the two proper-subgraph bounds with total right-hand side $H(|V|-1)-2$, one less than $(H-1)|E|$ for a tight graph. After suppressing bivalent chains, the sum of the excess valences over two is $2n-2=4$. There are at most four branch vertices. If every pair of branch vertices has at most one joining edge, minimum degree three forces $K_4$. Otherwise a parallel reduction lowers the loop number. At loop number at most two the excess valences allow at most two branch vertices, giving series and parallel moves to a single edge. The two-vertex graph consists of parallel edges directly.

For a convergent monomial with denominator $U^q$, push forward each child measure along its effective resistance $r_i$. If it is $c_i r_i^{s_i-1}\,dr_i$, then a series node uses $r=\sum_i r_i$ and a Dirichlet factor with arguments $(s_i)$. A parallel node uses $1/r=\sum_i 1/r_i$ and $U=r^{-1}\prod_i(r_iU_i)$. In inverse resistances the Dirichlet arguments are $(q-s_i)$. The new radial exponent is $\sum_i s_i$ in series or $\sum_i s_i-(k-1)q$ for $k$ parallel children. Leaf exponents are $1+k_e$. At the root the radial measure is $dr/r$ and the projective chart fixes $r=1$. Every nontrivial simplex factor has positive integer arguments by absolute convergence. Its value is a factorial ratio.

Consequently **every purely holomorphic three-loop kernel, in every admissible dimension and at every finite current order, has coefficients in $\mathbb Q+\mathbb Q\zeta(3)$**. Series–parallel graphs have rational coefficients. At the eight momenta above all ten series–parallel graphs on $\mathbb C^3$ agree with the original ratio-coordinate engine, in addition to the six $K_4$ subdivisions in the table. [Independent catalog and source tests](check_sources.py)

## A linear determinant for the current numerator

For an integral conserved source take a particular flow $h_e$ satisfying the prescribed vertex divergences. Let $C_e$ be the oriented cycle column for edge $e$, and put $D_e=C_e(C_e,h_e)^{\mathsf T}$. In dimension $H$ and loop number $n$, the number of edges is $E=nH+1$. The two blocks of $X=\sum_{e\in\mathrm{edges}}t_eD_e=(M,N)$ have widths $n,H$. Kirchhoff's polynomial is $U=\det M$. The edge currents are $y_e=h_e-N^{\mathsf T}M^{-1}C_e$.

Vectorize an $n$ by $(n+H)$ matrix by putting its left block first in row order and its right block next in component-then-cycle order. Take the $E$ columns $D_e$, followed by the $n^2-1$ columns $A X$, for the matrix units $A$ in row order omitting the final diagonal unit. The determinant of these columns is $(-1)^{(n^2-1)E}U^n c$, where the oriented projective cofactor $c$ is defined by deleting the first edge from the cycle-current matrix and dividing that minor by $t_0$.

Replace the first column by $X$. This multiplies the determinant by $t_0$. Subtract each left block times $M^{-1}N$ from its right block. The right blocks of $X$ and $A X$ vanish. Bringing the $A X$ columns next to $X$ costs $(-1)^{(n^2-1)(E-1)}$. Their left-block determinant together with $X$ is $(-1)^{n^2-1}(\det M)^n$, by right multiplication by $M$ on each row. The remaining right-block determinant is $t_0 c$. Rational elimination of the constant $D_e$ columns now leaves a homogeneous linear determinant of size $n^2-1$. Dependent constant columns give zero immediately.

Orient a subdivided path consistently, so all its cycle columns equal $C_P$. A simultaneous variation of its edge times changes $X$ by $C_P v^{\mathsf T}$ in the span of the path's constant columns. The varied columns $A C_P v^{\mathsf T}$ span at most $n$ dimensions, including $C_P v^{\mathsf T}$. Their span in the quotient by the constant columns therefore has dimension at most $n-1$. Thus the degree in *all* times of any one path together is at most $n-1$. [Determinant and exact division](current_minor.py) · [Independent current interpolation](kernel.py)
