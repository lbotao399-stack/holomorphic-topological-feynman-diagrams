# A nonvanishing odd triangle on $\mathbb C^2\times\mathbb R$

At zero displacement and zero topological source, the four-loop reduced kernel in the [project convention](../../code/htk/CONVENTIONS.md) is

$$
K_\Gamma(\lambda;0)=-\frac{p\pi^4}{26880}.
$$

Here $p$ is the second-component source transfer specified below. The graph has $H=2$, $T=1$, $13-1=4(2+1)$ and satisfies the project’s sparsity condition on every vertex subset. It consists of three cyclically arranged two-terminal networks, each with one loop. Using one side in parallel with the series combination of the other two gives $E=1$. Thus $E>0$ allows a nonzero reduced kernel in two holomorphic dimensions.

**Graph and sources.** Orient the paths in the order $0\to1\to2\to0$. The ordered edges, listed one path at a time, are

$$
(03,31;04,41),\quad(15,52;16,62),\quad(27,70;28,89,90).
$$

The vertex sources in vertex order $0,\ldots,9$ are

$$
(-1,p-1),\ (-1,-1-p),\ (-1,-1),\
(1,0),\ (0,1),\ (1,0),\ (0,1),\ (1,0),\ (1,0),\ (-1,1).
$$

Fix the sum of the three effective resistances to one and write the resistances as $x,y,z=1-x-y$ in cyclic order. In subnetwork $i$, the two path resistances are $R_i/w_i$ and $R_i/(1-w_i)$, with $R_i$ its effective resistance. The time fractions on the first path are $(1-a_i,a_i)$. On the second path these are $(1-b_i,b_i)$ for $i=1,2$ and $(1-b_3-c_3,b_3,c_3)$ for $i=3$. All fractions are positive.

**Curvature.** Put $w_i=\sin^2\theta_i$, $0<\theta_i<\pi/2$. For each path take the Euclidean unit vector supported on its edges with entries equal to the square roots of the time fractions. Let $P_i,Q_i$ be the vectors on the first and second paths. In subnetwork $i$ the unit internal cycle and unit terminal flow are

$$
u_i=\cos\theta_i P_i-\sin\theta_i Q_i,\qquad
h_i=\sin\theta_i P_i+\cos\theta_i Q_i.
$$

With $(r_1,r_2,r_3)=(\sqrt x,\sqrt y,\sqrt z)$, the fourth unit cycle is $h=\sum_{i\in\mathrm{subnetworks}}r_i h_i$. In the orthonormal frame $(u_1,u_2,u_3,h)$ the only independent connection entries are $A_{i4}=r_i d\theta_i$. The disjoint path supports and $\langle P_i,dP_i\rangle=\langle Q_i,dQ_i\rangle=0$ give this formula directly. For the curvature $F=dA+A\wedge A$,

$$
F_{ij}=-r_i r_j d\theta_i\wedge d\theta_j,\qquad
F_{i4}=dr_i\wedge d\theta_i,\qquad
\mathrm{Pf}(F)=-d\sqrt{xyz}\wedge d\theta_1\wedge d\theta_2\wedge d\theta_3 .
$$

The project Euler form has factor $2^{4/2}=4$.

**Currents.** Let $A_i,B_i$ denote the time-fraction averages of the successive cumulative vertex sources on the two paths. Their components on the first two subnetworks and on the third are

$$
A_i=(a_i,0),\quad B_i=(0,b_i)\ (i=1,2),\qquad
A_3=(a_3,0),\quad B_3=(b_3,c_3).
$$

The outgoing branch currents are $q+(0,p),q,q$, respectively, where $q$ is a two-component current. Equating voltage drops on the paths and then summing the three branch drops gives

$$
q=-x\big((0,p)+A_1+B_1\big)-y(A_2+B_2)-z(A_3+B_3).
$$

If $Q$ is a branch's outgoing current, its first path carries initial current $w_iQ+U_i$ with $U_i=-(1-w_i)A_i+w_iB_i$. The three first-path initial currents and $q$ are integral independent loop-current coordinates. Indeed, removing the first edge of each first path and the first edge of the third second path leaves a spanning tree, and the sum of the two initial path currents in the third branch is $q$. Write $d^2q=dq^1\wedge dq^2$ and likewise for $U_i$. Wedging with the three $d\theta_i$ eliminates the $dw_i$ terms. Wedging with $d^2q$ then reduces the holomorphic form to $d^2q\wedge d^2U_1\wedge d^2U_2\wedge d^2U_3$ up to the fixed orientation.

For $i=1,2$, $d^2U_i=-w_i(1-w_i)da_i\wedge db_i$ at fixed angles. For the third branch,

$$
d(a_3+b_3)\wedge d^2U_3=w_3\,da_3\wedge db_3\wedge dc_3,\qquad
dc_3\wedge d^2U_3=0 .
$$

Consequently only the first component of $dq$ contributes a remaining internal differential. The second component contributes
$-(p+b_1-c_3)dx-(b_2-c_3)dy$. With $f=\sqrt{xyz}$, the scalar integrand against the positive coordinate measure in $x,y,\theta_1,\theta_2,\theta_3,a_1,b_1,a_2,b_2,a_3,b_3,c_3$ becomes

$$
4z\big[(b_2-c_3)\partial_x f-(p+b_1-c_3)\partial_y f\big]
\,w_1(1-w_1)w_2(1-w_2)w_3.
$$

The overall sign is fixed in the project's edge orientation by evaluating its original Pfaffian and Kirchhoff forms at $p=1$, $(x,y)=(1/5,3/10)$, $w_i=1/2$, and all path coordinates $1/4$. Its density in the $w_i$ coordinates has value $(-1/100)/\sqrt{48/25}=-\sqrt3/240$. The positive change from $d\theta_i$ to $dw_i$ has Jacobian one at this point. [Signed evaluation and eight further checks](verify/h2_odd_triangle.py) · [Independent edge-basis Pfaffian and period](verify/odd_triangle_h2_edgecheck.py)

**Integration.** The function $f$ is zero on the three sides of the central simplex $\Delta$. Integration by parts gives

$$
\int_\Delta z\partial_x f\,dx\,dy
=\int_\Delta z\partial_y f\,dx\,dy
=\int_\Delta f\,dx\,dy
=\frac{\Gamma(3/2)^3}{\Gamma(9/2)}=\frac{2\pi}{105}.
$$

The central integral leaves $b_2-b_1-p$. The first two independent path coordinates $b_1,b_2$ both range from zero to one. Their mean difference is zero. The last path occupies a simplex of area $1/2$. The angle factors give

$$
\int_0^{\pi/2}\sin^2\theta\cos^2\theta\,d\theta=\frac\pi{16},\qquad
\int_0^{\pi/2}\sin^2\theta\,d\theta=\frac\pi4.
$$

All other path intervals have length one, so

$$
K_\Gamma(\lambda;0)=
-4p\,\frac{2\pi}{105}\left(\frac\pi{16}\right)^2\frac12\frac\pi4
=-\frac{p\pi^4}{26880}.
$$

## The same four-loop graph family for every $H\ge2$

Use the zero-based edge set $I=(0,4,8,10)$ in the displayed edge order. Deleting these four edges gives a spanning tree. In the project's positively oriented cycle matrix $C$ its minor is $\det C_I=-1$. For each of $k\ge0$ successive subdivisions, split the current first segment at each of these four edge positions. Keep all first segments in their old positions and append their second segments in increasing edge order. Insert the four new vertices immediately before the old root in the same order. Assign unit source in the new holomorphic direction to each new vertex and total source $-4$ in that direction to the original root. The other source components extend by zero on the new vertices. Each subdivision preserves the loop number and the three odd networks, so $H=2+k$, $n=4$, $E=1$.

For an edge of time $t_e$ write its split times as $(1-s)t_e,st_e$. The new current one-forms on these four independent splitting directions have matrix $C_I\,\mathrm{diag}(t_e)$. The normalizing determinant $U$ remains the same upon subdivision. The cycle basis obtained by duplicating the split columns retains the project orientation since the number of splits is even. Each subdivision consequently contributes $-\prod_{e\in I}t_e/U$ to the density in all old coordinates and the four new $s$ coordinates. After $k$ rounds the first segment of each chosen edge has inherited all preceding factors $1-s$. Integrating the $k$ parameters on one such edge gives $\prod_{j=1}^k\int_0^1(1-s)^{k-j}ds=1/k!$. The resulting kernel is $(-1)^k/(k!)^4$ times the original integral with insertion $(\prod_{e\in I}t_e/U)^k$. [Canonical signs and subdivision checks](verify/dimension_raise.py)

In the original coordinates $U=xyz/\prod_{i=1}^3 w_i(1-w_i)$, so the insertion before taking the $k$th power is

$$
\frac{t_0t_4t_8t_{10}}{U}=
z(1-w_1)(1-w_2)(1-a_1)(1-a_2)(1-a_3)(1-b_3-c_3).
$$

The two central integrations now give

$$
\int_\Delta z^{k+1}\partial_x f\,dx\,dy
=\int_\Delta z^{k+1}\partial_y f\,dx\,dy
=(k+1)\frac{\Gamma(3/2)^2\Gamma(k+3/2)}{\Gamma(k+9/2)}.
$$

The $b_2-b_1$ term still integrates to zero. The three $a_i$ intervals and the last path simplex give $1/[(k+1)^4(k+2)]$. Each of the first two angles contributes $\Gamma(3/2)\Gamma(k+3/2)/(2\Gamma(k+3))$, and the third contributes $\pi/4$. Hence the reduced kernel in holomorphic dimension $2+k$ is

$$
\frac{(-1)^{k+1}p\pi}
{(k!)^4(k+1)^3(k+2)}
\frac{\Gamma(3/2)^2\Gamma(k+3/2)}{\Gamma(k+9/2)}
\left(\frac{\Gamma(3/2)\Gamma(k+3/2)}{2\Gamma(k+3)}\right)^2 .
$$

At $p=1$ the values for $H=2,3,4$ are $-\pi^4/26880$, $\pi^4/3870720$ and $-25\pi^4/39239811072$. Every factor in the absolute value is strictly positive for each integer $k\ge0$. [Exact original and subdivided density checks](verify/h2_odd_triangle.py)
