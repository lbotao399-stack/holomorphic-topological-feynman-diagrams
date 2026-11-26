# Three series–parallel four-loop graphs

Use the edge order of [graphs 0–2](../four-loop-graphs.json). Write $[pq]$ for the determinant of two holomorphic momenta. With external momenta $p_6,p_7,p_8,p_9$ at the four bivalent vertices, summing to zero, their bare geometric kernels are

$$
K_0=-\frac{\pi^2}{42000}[p_6p_7]^2[p_8p_9]^2,\qquad K_1=K_2=0 .
$$

**Network coordinates.** Orient the three arms from vertex 0 to vertex 3. The two compound arms pass through $(1,2)$ and $(5,4)$. Their outer conductance fractions are $a,b$ and the direct arm has $c=1-a-b$. In a compound arm write $u$ for the first parallel-path share, $x$ for the fraction of the total arm resistance in the parallel part, and $s$ for its right lead fraction, $0<s<1-x$. For the other arm use $v,y,t$. Each arm contains momenta $p,q$ in oriented path order, with total $P=p+q$. Its two path offsets $m_1,m_2$ are the sums of the incoming momenta weighted by the resistance fractions following them on that path. Put $M=m_1+m_2$ and $O=xM+sP$. On the first arm the incoming current $Y$ and first inner path current $Z$ are

$$
Y=a(O_A+O_B)-O_A,\qquad Z=u(Y+M)-m_1 .
$$

On the other arm replace $a,u,O_A$ by $b,v,O_B$. The geometric form is the project Euler form times the eight differentials of the independent currents, with its sign fixed by the positive source cycle basis.

**Euler factor.** In the two internal-cycle directions and two outer-cycle directions the Gram matrix becomes block diagonal after the shear $Z\mapsto Z-uY$. The first two entries are $m_A=x/(a u(1-u))$, $m_B=y/(b v(1-v))$. The outer block is $\operatorname{diag}(1/a,1/b)+\mathbf1\mathbf1^\mathsf{T}/c$. In the sheared frame the original metric differential has off-diagonal entries $-m_A\,du$ and $-m_B\,dv$. Substitution into $\operatorname{Pf}(-dM_{\rm Gram}M_{\rm Gram}^{-1}dM_{\rm Gram}/2)/\sqrt{\det M_{\rm Gram}}$ gives, up to the orientation of the cycle basis,

$$
\frac{du\wedge dv}{4\sqrt{u(1-u)v(1-v)xyc}}\wedge
\left[
 c\,dx\wedge dy
 +y\big(b\,da+(b+c)\,db\big)\wedge dx
 -x\big((a+c)\,da+a\,db\big)\wedge dy
 +3xy\,da\wedge db
\right].
$$

**Current factors.** Since the Euler form contains $du\wedge dv$, one may use $dZ-u\,dY=u\,dM-dm_1$ in its product with the holomorphic forms. Its two holomorphic components use both path-tail coordinates of that arm. Denote those coordinates by $m,n$. For lengths $(1,3)$ they satisfy $1>m>n>0$, while for $(2,2)$ they vary independently over $(0,1)$. In both cases $M=mp+nq$. The coefficient of the two tail differentials is $[pq]$ times $u^2$ for $(1,3)$ and $-u(1-u)$ for $(2,2)$.

The two lead fractions $s,t$ occur only in the outer currents. Their differentials contribute the outer determinant $c$ in the direction of $P_A=-P_B$. For the other component put $K=[P_A,M_A]$, $J=[P_A,M_B]$, and $L=xK+yJ$. The square bracket of the Euler factor wedged with
$d(aL-xK)\wedge d(bL-yJ)$ has coefficient $3cxyKJ$ in the order $da,db,dx,dy$. This identity follows by expanding its six minors. In particular $[P_A,M_A]=(n-m)[p,q]$ and $[P_B,M_B]$ has the same formula on its arm. The identity extends to $P_A=0$ by polynomiality.

After integration over $s,t$ the radial integral is

$$
\frac34\int_{a,b>0,\ a+b<1} c^{3/2}\,da\,db
\left(\int_0^1\sqrt{x}(1-x)\,dx\right)^2
=\frac34\frac4{35}\left(\frac4{15}\right)^2.
$$

For $(2,2)$ the signed tail integral is $\int_0^1\int_0^1(m-n)\,dm\,dn=0$. For $(1,3)$ the two remaining integrals are

$$
\int_{0<n<m<1}(m-n)\,dm\,dn=\frac16,\qquad
\int_0^1\frac{u^2\,du}{\sqrt{u(1-u)}}=\frac{3\pi}{8}.
$$

Only graph 0 has two arms of type $(1,3)$. Its absolute coefficient is therefore $(3/4)(4/35)(4/15)^2(\pi/16)^2=\pi^2/42000$. In the source edge order the coefficient is negative. The replay compares the full signed density to a direct Kirchhoff solution, its Pfaffian and the absolute projective coordinate Jacobian, at two rational points on each graph.

**Cyclic color.** Use the cyclic matrix units of the [order-four probe](../order_four.py). On its four external momenta $(1,0),(0,1),(1,1),(-2,-2)$, the signed sum over all 24 labelings, with both momentum determinants squared, gives $-16N^4$ for graph 0 before its geometric factor. Dividing by the automorphism order gives

$$
\frac14\left(-\frac{\pi^2}{42000}\right)(-16N^4)
=\frac{\pi^2 N^4}{10500}.
$$

This coefficient precedes the common four-loop heat-kernel normalization. Cyclic trace contraction and direct gluing of row and column indices on the four matrix units give the same polynomial. [Exact replay](check.py)
