# All quadratic local operators

Let $R=\mathbb C[z^1,z^2,z^3]\otimes\Lambda(\theta_1,\theta_2)$ be the coordinate superalgebra of the source and $\mathfrak g$ a finite-dimensional complex simple gauge algebra with the action pairing $\langle\,,\rangle$. The field-length-two classical cohomology of finite-jet local operators and its quantum lift are
$$
H^2(Q_0)=(\Omega^1_R/dR)^\vee_{\mathrm{fin}},
\qquad Q_n\,\ell\langle C,dC\rangle=0\quad(n\geq1).
$$
Here $d$ is the parity-preserving universal derivation, so $d\theta_i$ is odd, and $\vee_{\rm fin}$ takes functionals supported in finitely many coordinate multidegrees. The quotient by $dR$ acts on source-coordinate forms. The local operators themselves retain total holomorphic derivatives.

## Classical differential

Taylor expand $C=c+\theta_1g_1+\theta_2g_2+\theta_1\theta_2b$ with its usual Taylor factorials. Then $\widehat Q_0C=[C,C]/2$ is the Chevalley–Eilenberg differential of the current algebra $\mathfrak g\otimes R$ on polynomial cochains of finite jet order, with the repository's overall sign. Write $Xu=X\otimes u$ for $X\in\mathfrak g$ and a homogeneous source monomial $u$. The bracket is $[Xu,Yv]=[X,Y]uv$. The constant $\mathfrak g$ action on cochains acts trivially on cohomology by the Cartan homotopy with insertion of a constant color. Every homogeneous field length and multidegree is finite-dimensional. Complete reducibility therefore computes the cohomology from invariant cochains. Since the adjoint representation is irreducible and self-dual through the action pairing, every invariant two-cochain has the form
$$
\omega(Xu,Yv)=\langle X,Y\rangle\psi(u,v),
\qquad
\psi(u,v)=-(-1)^{|u||v|}\psi(v,u),
$$
where $|u|$ is the parity of $u$. The three cyclic color contractions in the cocycle equation all equal $\langle[X,Y],Z\rangle$. Canceling this nonzero tensor gives
$$
(-1)^{|u||w|}\psi(uv,w)+(-1)^{|v||u|}\psi(vw,u)
+(-1)^{|w||v|}\psi(wu,v)=0.
$$
Using graded antisymmetry and supercommutativity gives
$$
\psi(uv,w)=\psi(u,vw)+(-1)^{|u||v|}\psi(v,uw).
$$
Putting $u=v=1$ yields $\psi(1,w)=0$. Equivalently,
$$
\psi(u,vw)=\psi(uv,w)+(-1)^{|v||w|}\psi(uw,v).
$$
This is exactly the Leibniz relation for $u\,d(vw)=uv\,dw+(-1)^{|v||w|}uw\,dv$. Thus $\ell(u\,dv)=\psi(u,v)$ descends to $\Omega^1_R$, and $\ell(dw)=0$. Conversely every $\ell$ on $\Omega^1_R/dR$ gives graded antisymmetry by $\ell(d(uv))=0$, and Leibniz gives the cocycle equation. There are no invariant linear cochains because $\mathfrak g$ has no invariant linear form. Hence these representatives have no two-cochain boundary identifications. Polarization of $\ell\langle C,dC\rangle$ gives twice the displayed two-cochain with the usual suspended coefficient signs. This describes the central terms allowed in the graded Jacobi identity for current brackets. Each homogeneous $\ell$ gives a central generator of the parity on which $\ell$ is supported.

## Box and higher orders

The repository functional derivative contracts a field at odd point $\eta$ with a root slot at $\theta$ through
$$
\delta^2(\theta-\eta)=(\theta_1-\eta_1)(\theta_2-\eta_2).
$$
The two contracted fields in $\langle C,dC\rangle$ have the same $\eta$, and both ends of the box have the same $\theta$. A holomorphic source derivative preserves both odd deltas. An odd source derivative hits precisely one. With left odd derivatives,
$$
(\delta^2)^2=0,\qquad
\partial_{\eta_1}\delta^2=-(\theta_2-\eta_2),\qquad
\partial_{\eta_2}\delta^2=\theta_1-\eta_1,\qquad
\delta^2\partial_{\eta_i}\delta^2=0 .
$$
This holds at arbitrary holomorphic root displacements, before applying $\ell$ to source forms, extracting any Taylor coefficient, or integrating. Hence $Q_1$ annihilates every displayed representative. At order $n\ge2$ the $n+1$ root derivatives annihilate any bilinear. Every $Q_n$ increases field length by exactly one, so the same is true of every differential of the $\hbar$-adic spectral sequence. An incoming differential to length two would require a class of length one. But $[\mathfrak g\otimes1,\mathfrak g\otimes R]=\mathfrak g\otimes R$ gives $H^1(Q_0)=0$. Every quadratic class therefore survives at all orders. [Root arities and normalization](../VERIFY.md) · [Component functional derivative](../s2code/salg.py)

## Exact count

Let $t_1,t_2,t_3,u_1,u_2$ record the three holomorphic and two odd-coordinate weights. Each $dz^i$ and $d\theta_j$ contributes one free copy of $R$ shifted by its coordinate weight. The kernel of $d:R\to\Omega^1_R$ consists of constants, since contraction with the Euler vector field on a homogeneous nonconstant monomial gives its positive total degree times that monomial. The ordinary dimension series, counting both parities positively, is therefore
$$
1+\big(t_1+t_2+t_3+u_1+u_2-1\big)
\frac{(1+u_1)(1+u_2)}{(1-t_1)(1-t_2)(1-t_3)}.
$$
At odd weights $(2,1)$ and $(1,2)$ and zero holomorphic weight, the one-dimensional spaces have representatives $\langle g_1,b\rangle$ and $\langle g_2,b\rangle$, respectively. Both vanish on $b=c=0$ while surviving in local quantum cohomology. At odd weight $(1,1)$ and first holomorphic weight along one axis, the dimension is two. The symmetric scalar derivative and the mixed translation current have independent restrictions to diagonal scalar jets. [Translation current](scalar_products/regular_jets.md)

The replay solves the signed cyclic and graded antisymmetry equations at every holomorphic total degree at most three and both odd weights from zero to two. Its 180 multidegrees give 201 classes. It also evaluates the native odd delta, its two first derivatives, and the classical and box variations of the two scalar–$b$ examples. [Exact replay](quadratic.py)

## Integrated quadratic charges

For every complex simple gauge algebra, the translation-invariant quadratic quantum cohomology modulo total holomorphic derivatives has **nine even and eight odd generators**. They are the divergence-free source vector fields with coefficients independent of the three even coordinates.

Let the pairing of two source coefficients be their product integrated in the three even coordinates and in both odd coordinates; even-coordinate integration denotes the quotient by total derivatives. It is nondegenerate on Fourier modes of opposite momenta and complementary odd monomials. An invariant bilinear cochain is therefore uniquely represented, after integrating by parts, by the gauge trace pairing times a source operator $D$: its value on $Xu,Yv$ is $\langle X,Y\rangle\int(Du)v$. Here $D$ is polynomial in the three holomorphic derivatives, with coefficients acting on the four-dimensional odd-coordinate algebra. Invariant linear cochains vanish for the simple algebra. Averaging over the compact real form commutes with total derivatives and computes the absolute cohomology.

For homogeneous $D$, antisymmetry and the classical cocycle equation, respectively, give
$$
\int (Du)v+(-1)^{|D||u|}\int u(Dv)=0,\qquad
D(uv)=(Du)v+(-1)^{|D||u|}u(Dv).
$$
Indeed, insert $Xu,Yv,Zw$ into the cyclic cocycle equation and move $D$ off $w$ by the first identity. The remaining expression is the product of $\langle[X,Y],Z\rangle$ with the pairing of the displayed Leibniz defect against $w$. Varying $w$ proves the second identity. Conversely the two identities make the cyclic sum vanish.

A derivation is fixed by its values on $z^i,\theta_1,\theta_2$. Translation invariance restricts all five coefficients to the odd-coordinate algebra. Skew-adjointness is exactly vanishing divergence. Thus each of the three $\partial_i$ can be multiplied by $1,\theta_1,\theta_2,\theta_1\theta_2$, giving six even and six odd operators. For the two odd directions the allowed fields are
$$
\partial_{\theta_1},\quad \partial_{\theta_2},\qquad
\theta_1\partial_{\theta_2},\quad \theta_2\partial_{\theta_1},\quad
\theta_1\partial_{\theta_1}-\theta_2\partial_{\theta_2}.
$$
The first two are odd, the last three even. The missing trace of the linear transformation has constant divergence, and the two quadratic coefficients have independent linear divergences. In Fourier variables a spatial derivative monomial of order $k\ge2$ would give a mixed term in $(p+q)^k$ in the Leibniz equation; no such operator survives. Exact tests of skew-adjointness and all odd monomial products give parity dimensions $(3,2)$ at order zero, $(2,2)$ per spatial direction at order one, and $(0,0)$ at orders two through six. [Fourier-symbol rank replay](integrated_quadratic.py)

In the quadratic superfield density one factor has no derivative and the other has at most one. The two box contractions contain the top odd delta times itself or one of its first derivatives, hence vanish identically. At higher order there are more than two spoke derivatives. Thus all seventeen integrated charges have strict representatives at finite rank and no invariant linear primitive can generate an incoming boundary.
