# Two-loop theta kernels as hyperbolic polygon areas

At $T=1$, zero topological source and $p+q+r=2H+3$, fix the initial currents $b,c\in\mathbb R^H$ on paths $B,C$ oriented toward the root. Conservation fixes the initial current on $A$ as $\lambda_0-b-c$. The current $y_e$ on each successive segment is its initial current plus the cumulative sources on its path.

Let $S$ be the three-row matrix of the path indicators in order $A,B,C$. Use the cycle matrix $C$ whose rows are the indicators of $B-A$ and $C-A$. In the following matrices $C$ denotes this matrix. For current component $i$, the voltage rows are $C\,\operatorname{diag}(y^i)$. Stack those rows over $1\le i\le H$ to form $\Psi$. The original positive Schwinger times with the specified currents and path totals $R_A,R_B,R_C$ are determined by

$$
\begin{pmatrix}\Psi\\ S\end{pmatrix}t=
\begin{pmatrix}0\\R_A\\R_B\\R_C\end{pmatrix}.
$$

The matrix is square. On $St=0$, the initial currents drop out of every voltage row. Eliminate the first segment on each path. The remaining matrix equates the three sums of cumulative source columns weighted by the remaining times. Taking consecutive differences of cumulative columns gives the internal source columns. Subtracting one cycle row block from the other puts the matrix into the block form defining $D_\theta$ in the [kernel conventions](../../code/htk/CONVENTIONS.md). Thus its determinant is $\pm D_\theta$ independently of $b,c$. When $D_\theta=0$ the project theta density is zero. Otherwise let $W$ be the last three columns of the inverse in the displayed system. The admissible fibre is the convex polygon

$$
\mathcal P(b,c)=\{R_A+R_B+R_C=1,\ WR\ge0\}.
$$

Only its interior contributes. Its vertices are the feasible intersections of pairs of the segment boundaries $(WR)_e=0$. Since $SW$ is the identity, $WR\ge0$ already implies $R\ge0$.

The period matrix is $M=\left(\begin{smallmatrix}R_A+R_B&R_A\\R_A&R_A+R_C\end{smallmatrix}\right)$. For its entries $a=M_{11},\,u=M_{12},\,v=M_{22}$ and $\Delta=av-u^2$, the two-loop Euler form is

$$
\rho_2=-\frac{v\,da\wedge du-u\,da\wedge dv+a\,du\wedge dv}{2\Delta^{3/2}}
=-\frac{dx\wedge dy}{y^2},\qquad x=\frac ua,\quad y=\frac{\sqrt\Delta}a .
$$

A projective line has equation $\alpha+\beta x+\gamma(x^2+y^2)=0$, a geodesic in the upper half-plane. The $T=1$ source factors are the exterior product of the $2H$ loop-current differentials up to the project orientation. Consequently, for zero displacement or any edge-current insertion, the pushforward to $b,c$ is the signed hyperbolic area of the fibre polygon. More explicitly, let $s=(-1)^{\lceil H/2\rceil+(H+1)(q+1)}$ be the project's theta sign and let $f$ be the polynomial current insertion or displacement exponential. For a nonempty polygon with $k$ sides and interior angles $\theta_j$, counting an ideal vertex as angle zero,

$$
K_\Gamma[f]=\operatorname{sgn}(sD_\theta)
\int_{\mathbb R^{2H}} f(y(b,c))\left((k-2)\pi-
\sum_{j\in\mathrm{vertices}}\theta_j\right)d^H b\,d^H c .
$$

An empty or degenerate polygon contributes zero. The area factor follows by integrating the curvature $-1$ over a geodesic polygon. In particular, the measure at fixed currents generally contains a nonconstant angle, even when its integrated moments are rational multiples of $\pi$.

For $(1,H+1,H+1)$ fix $R_A=1$ and write $u=R_B/R_A$, $v=R_C/R_A$. On either long path $P$, the admissible ratios are exactly those positive numbers $r$ for which $y_A/r$ is in the interior of the simplex spanned by the $H+1$ path currents. This is an interval, obtained by intersecting the ray through $y_A$ with the simplex, or the whole positive interval or the empty set at $y_A=0$. Thus the fibre is $l_B<u<h_B$ and $l_C<v<h_C$. In these coordinates the positive area density is

$$
\frac{du\,dv}{2[(1+u)(1+v)-1]^{3/2}}.
$$

Its integral from $u$ and $v$ to infinity equals

$$
2\arctan\frac1{\sqrt{(1+u)(1+v)-1}}.
$$

The four endpoint evaluations with signs $+,-,-,+$ give the full conditional density in the $(1,3,3)$ and $(1,6,6)$ sectors. Infinite endpoints contribute zero. For an exact example, let $y_A=(1,\ldots,1)$ and the ordered currents on each long path be $0,2H e_1,\ldots,2H e_H$ with $e_i$ the coordinate unit vectors. The corresponding sources are their successive differences. Both ratios range from $1/2$ to infinity. The unsigned conditional density is $2\arctan(2/\sqrt5)=\arccos(1/9)$. The polytope has two ideal vertices and a third with angle cosine $-1/9$. All seven $H=2,5$ theta shapes and the $H=1$ bitriangle were checked pointwise against the original signed Euler–Kirchhoff density, including the projective Jacobian. [Exact polygon and density replay](verify/current_polygon.py)
