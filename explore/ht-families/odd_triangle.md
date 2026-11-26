# Odd-loop subnetworks on $\mathbb R\times\mathbb C$

Let a connected graph of even loop number $n$ have $N=2n+1$ edges. Suppose it consists of three connected two-terminal subnetworks arranged around a triangle. Their interiors are disjoint and each touches the other two only at its terminals. If all three loop numbers are odd, then at zero topological source

$$
K_\Gamma(\lambda;z)=0
$$

for all holomorphic sources and displacements. Any polynomial in the holomorphic edge currents may be inserted. The subnetworks need not be series–parallel.

**Cycle frame.** Orient terminal current cyclically. Keep the resistance ratios inside each subnetwork fixed. Write $R_i$ for its effective resistance and take $x_i=R_i/(R_1+R_2+R_3)$ for its share of the total, with $i=1,2,3$ labeling subnetworks. The overall scale is immaterial. Varying $x_i$ scales all edge resistances of the $i$th subnetwork together. In the Euclidean edge space take the frame of weighted cycles $\mathrm{diag}(\sqrt t)C^{\mathsf T}$, as in the [sector proof](README.md). The orthonormal internal cycle frames $U_i$ then depend only on internal ratios. Gram–Schmidt on a fixed internal cycle basis supplies them smoothly in the interior.

Let $h_i$ be the weighted unit terminal current of the $i$th subnetwork divided by $\sqrt{R_i}$. The sum of $t_e y_e^2$ for a unit terminal current is $R_i$, and the loop voltage equation makes $h_i$ orthogonal to $U_i$. The one further normalized cycle vector is
$$
h=(\sqrt{x_1}h_1,\sqrt{x_2}h_2,\sqrt{x_3}h_3).
$$
Indeed the common, unweighted unit current through the three subnetworks is conserved at every vertex. The projected connection in the frame $(U_1,U_2,U_3,h)$ has diagonal internal blocks $A_i=U_i^{\mathsf T}dU_i$. Its off-diagonal column incident to $h$ has blocks $b_i=\sqrt{x_i}\,U_i^{\mathsf T}dh_i$. The derivatives of $\sqrt{x_i}$ drop out by orthogonality. The final diagonal entry is zero.

**Pfaffian transgression.** Interpolate the connection from $A_0=\mathrm{diag}(A_1,A_2,A_3,0)$ to the projected connection $A_1^{\rm full}$ by multiplying its off-diagonal column $b$ and row $-b^{\mathsf T}$ by $s\in[0,1]$. Denote this connection by $A_s$ and its curvature by $F_s=dA_s+A_s\wedge A_s$. Its internal block is
$$
(F_s)_{\rm int}=\mathrm{diag}(dA_i+A_i\wedge A_i)-s^2 b\wedge b^{\mathsf T}.
$$
Write $\mathrm{Pf}(\cdot,\ldots,\cdot)$ for the symmetric polarization of the project's Pfaffian, with $\mathrm{Pf}(F,\ldots,F)=\mathrm{Pf}(F)$. Its covariant derivative satisfies the product rule by orthogonal invariance. Since $\partial_sF_s=d_{A_s}(\partial_s A_s)$ and $d_{A_s}F_s=0$, integration of $\partial_s\mathrm{Pf}(F_s)$ gives
$$
\rho_n=d\!\left(2^{n/2}\frac n2\int_0^1
  \mathrm{Pf}(\partial_s A_s,F_s,\ldots,F_s)\,ds\right)
$$
in a compatibly oriented frame. The initial Pfaffian vanishes because of its final zero row. In every term of the primitive, $\partial_s A_s$ consumes the index of $h$ exactly once. The other factors consequently use only $(F_s)_{\rm int}$. Every factor has central differential degree zero, where the two central variables are the independent entries of $x$. On $x_i=0$, $b_i=0$ decouples internal block $i$. That block has odd size and its indices cannot be paired among themselves. The primitive vanishes on all three faces of the central triangle.

**Currents and integration.** For fixed internal ratios and fixed internal sources, prescribing one terminal current determines all edge currents of a subnetwork affinely. Uniform scaling of its resistances preserves these currents. The three junction conservation laws have one free scalar, the circulating current around the triangle. Thus every current, at fixed internal ratios, depends on $x$ through this single scalar. In fact the voltage drop on each oriented subnetwork is its resistance times the prescribed terminal current plus the drop due to internal sources. The latter scales by the same resistance. The three terminal currents differ by constants fixed by conservation. The sum of the three drops is zero, so the circulating current is affine in $x$ when $\sum_{i\in\mathrm{subnetworks}}x_i=1$. All edge currents are smooth on the closed central triangle.

The holomorphic factor is, up to the project orientation, the wedge of the $n$ independent loop-current differentials $d\ell_r$. The edge currents are globally affine in these loop currents at fixed $\lambda$, since $y=y^0+C^{\mathsf T}\ell$. Therefore multiplying the wedge by any polynomial in $y$ and by the displacement exponential gives a closed form. This form has central degree at most one because the central derivatives of all currents have rank at most one. Wedge it with the preceding primitive. In the resulting form the central degree is still at most one. The component of its exterior derivative with two central directions is therefore precisely the central derivative of its component with one central direction.

Use $\sqrt{x_i}$ as coordinates on the positive octant of the unit sphere. For each fixed interior choice of the internal ratios all these expressions are smooth up to its faces, and the primitive vanishes there. Stokes' theorem on this octant gives zero for the central integral. Absolute convergence on the full projective domain, proved in the [sector calculation](README.md), permits integration over the internal ratios.

**Application to c5.** If the cycle matroid is disconnected, separate the edge set into cycle blocks and choose internal cycle bases. Their period matrices are in direct sum. Independent scaling of the times of a block leaves all currents fixed. In that direction contraction into $dM\wedge M^{-1}dM$ vanishes since it gives $M_{\rm block}M^{-1}dM-dM M^{-1}M_{\rm block}=0$ in every block. Bridge times also leave the form fixed with zero contraction. Any such scaling independent of the global scale annihilates the projective top form pointwise.

Consider the remaining graphs with connected cycle matroid. In a series–parallel reduction with $E>0$, some series node contains at least two odd-loop children. Take two successive odd children, absorbing all even children between them into the first subnetwork. This gives two adjoining two-terminal pieces of odd loop number. Their common junction has no exterior attachment, and their other two terminals are distinct. All remaining edges connect those other terminals into one two-terminal subnetwork. Indeed, if no exterior route joins them, the middle junction is a cut vertex separating the two selected pieces. A component attached at just one terminal would also split the cycle matroid. Both possibilities have already been disposed of. With three pieces around the triangle the total loop number is one plus the sum of their loop numbers. The complement has odd loop number. The theorem therefore applies to every $E>0$ reduction at $H=T=1$.

For the six-loop census, the following triples are the terminal vertices in cyclic order. The final column gives loop numbers for the three cyclic pairs in that order.

| Graph in the repository census | Terminals | Loop numbers |
|---|---|---|
| 6 | $1,2,3$ | $3,1,1$ |
| 13 | $3,4,6$ | $3,1,1$ |
| 16 | $0,2,5$ | $1,3,1$ |
| 25 | $0,2,4$ | $1,1,3$ |
| 29 | $0,2,4$ | $1,1,3$ |
| 30 | $0,2,4$ | $1,1,3$ |

[Exact decomposition and current replay](verify/odd_triangle.py)
