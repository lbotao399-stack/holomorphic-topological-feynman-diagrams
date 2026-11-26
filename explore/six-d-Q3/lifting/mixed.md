# Mixed scalar quartics for $\mathfrak{sl}_3$

Write $x=g_1$, $y=g_2$, and $\langle x,y\rangle=\operatorname{tr}(xy)$ in the defining representation. At bidegree $(2,2)$ the classical scalar functions on commuting pairs are spanned by
$$
\operatorname{tr}(x^2)\operatorname{tr}(y^2),\qquad
\operatorname{tr}(xy)^2.
$$
In fact, commuting complex matrices are simultaneously triangularizable: an eigenspace of either is preserved by the other, giving a common eigenvector, and induction works on the quotient. Conjugating a triangular pair by a diagonal one-parameter subgroup takes it to its diagonal entries in the limit. A polynomial invariant restricts to a permutation-invariant polynomial of the diagonal entries. Sums of products with distinct entry labels reduce to products of unrestricted sums by subtracting repeated-label terms. Tracelessness deletes all singleton traces. Bidegree $(2,2)$ then allows only the two products above and $\operatorname{tr}(x^2y^2)$. Polarizing the traceless $3\times3$ identity $\operatorname{tr}(M^4)=\operatorname{tr}(M^2)^2/2$ gives
$$
6\operatorname{tr}(x^2y^2)=\operatorname{tr}(x^2)\operatorname{tr}(y^2)+2\operatorname{tr}(xy)^2
$$
on commuting pairs. The right-hand side is the strict polarized direction. The complementary density
$$
P=\operatorname{tr}(x^2)\operatorname{tr}(y^2)-\operatorname{tr}(xy)^2
$$
vanishes at $x=y$. It takes value $12$ at $x=\operatorname{diag}(1,-1,0)$ and $y=\operatorname{diag}(1,1,-2)$. The polarized direction is nonzero at $x=y=\operatorname{diag}(1,-1,0)$, so the two are independent.

## Heat-kernel calculation

Use the local conventions $Q_0=i\widehat Q_0$ and $Q_1=i\widehat Q_{\rm box}/(16\pi^3)$ for the coefficient of $\hbar$. Write $A_i=\operatorname{ad}_{\partial_i c}$, the adjoint action of the odd ghost jet. For the holomorphic volume tensor of the box and the trace pairing put
$$
U=\epsilon^{ijk}\langle x,A_i A_j\partial_k y\rangle .
$$
The two functional derivatives on $x,y$ multiply by $\theta_2(-\theta_1)$ and in reverse order by $(-\theta_1)\theta_2$. Their difference takes the alternating part of the mixed Hessian of $P$, which at even test colors $V,W$ equals
$$
6\big(\langle x,V\rangle\langle y,W\rangle-
\langle x,W\rangle\langle y,V\rangle\big).
$$
The ordered three-simplex has volume $1/6$. The three consecutive structure tensors give the matrix product $A_i A_j A_k$. Its contraction with $\epsilon^{ijk}$ is skew with respect to the trace pairing, since the transpositions of the three odd ghosts compensate the sign from reversing the three volume indices. It follows that
$$
\widehat Q_{\rm box}P=2\epsilon^{ijk}\langle x,A_i A_j A_k y\rangle=2\widehat Q_0U.
$$
For the last equality, all covariant ghost actions cancel in the invariant pairing. The inhomogeneous term in $\widehat Q_0\partial_k y= [c,\partial_k y]+[\partial_k c,y]$ supplies the last adjoint action.

Only $x,\partial_k y$ can be hit by a second box. The exchanged ends are paired against the symmetric tensor $\epsilon^{ijk} A_i A_j$. For ordered positions $1>s_1>s_2>s_3>0$ the moments of $1-2s_j$ are $(-1/12,0,1/12)$. With $p_1,p_2$ on the two spectator ghosts and $p_3,p_4,p_5$ in box order, the symbol of the five-ghost adjoint trace is
$$
\frac{1}{12}\det(p_1,p_2,p_5-p_3)\det(p_3,p_4,p_5).
$$
The first product of determinants is $\det(p_3,p_4,p_5)\det(p_5,p_1,p_2)$; the subtracted one is $\det(p_1,p_2,p_3)\det(p_3,p_4,p_5)$. Cyclic rotation by two slots preserves the trace of five odd fields and identifies these terms. Hence $\widehat Q_{\rm box}U=0$. Every positive-loop graph with $n\geq2$ differentiates $n+1\geq3$ letters. On $P$ or $U$, each scalar derivative supplies at least one odd coordinate and each ghost derivative two at the common $\theta$. These graph integrands vanish.

The exact local representative of a general bidegree $(2,2)$ scalar density on the commuting-pair quotient is therefore
$$
a\,\operatorname{tr}(x^2)\operatorname{tr}(y^2)+b\,\operatorname{tr}(xy)^2
-\frac{\hbar(2a-b)}{24\pi^3}\,U .
$$
Here $a,b$ are numerical coefficients. In particular $P$ needs $-\hbar U/(8\pi^3)$. Evaluating on constant commuting pairs annihilates $Q_0$ of every primitive, every positive-loop current determinant, and all total derivatives of translation-invariant local densities. It distinguishes both lifted directions locally and after integration. [Heat and graph conventions](../VERIFY.md) · [Root-valence bound](../README.md)

On the displayed diagonal pair give $x$ the plane-wave momentum $(-1,-1,-1)$, keep $y$ constant, and let the three ghost waves have the standard basis momenta, independent odd coefficients, and colors $E_{12},E_{23},E_{31}$, where $E_{ij}$ is the elementary matrix. The determinant and the ghost permutation have the same sign. The six nested brackets paired with $x$ give $6,6,3,0,3,0$. Thus the reduced box of $P$ has coefficient $36$ at zero total momentum. This tests the necessity of the counterterm even after integration. A second check uses the repository's Grassmann jet calculus with the three-color structure tensor: the 42 nonzero terms of the box agree with $2\widehat Q_0U$, and its box on the 54-term $U$ vanishes. [Exact replay](mixed.py)
