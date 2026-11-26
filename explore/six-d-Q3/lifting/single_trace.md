# Trace descendants and marked products on $\mathbb C^{3|2}$

Let $d$ be the odd de Rham differential on the holomorphic source. Then $dz^i$ is odd, $d\eta_j$ is even, and the differential $dC$ of the odd superfield at $(z,\eta)$ is total-even. For any invariant trace in a finite-dimensional representation, every coefficient of $\operatorname{Tr}((dC)^p)$ is a strict quantum cocycle. This also holds after taking holomorphic and odd-coordinate jets of the whole expression. The action on $dC$ in the classical term is an adjoint action, so cyclicity gives $Q_0\operatorname{Tr}((dC)^p)=0$. [Source definition and classical identity](../../../sources/LLC_main.tex)

## Common odd root factor

The odd part of the propagator from a root at $\theta$ to $C(z,\eta)$ is
$$
\delta^2(\theta-\eta)=(\theta_1-\eta_1)(\theta_2-\eta_2).
$$
The holomorphic gradient on a hit $dC$ differentiates the even distribution and retains $\delta^2$. The odd-coordinate gradient on a hit $dC$ retains the undifferentiated even distribution. With the repository's left functional derivative, its form factor is
$$
-d\eta_1\partial_{\eta_1}\delta^2-d\eta_2\partial_{\eta_2}\delta^2
=d\eta_1(\theta_2-\eta_2)-d\eta_2(\theta_1-\eta_1).
$$
Both slots carry this same linear form, independently of their color and body displacements. Its square vanishes because $d\eta_1,d\eta_2$ commute and $\theta_1-\eta_1,\theta_2-\eta_2$ anticommute. Its product with $\delta^2$ and the square of $\delta^2$ also vanish. The remaining $dC$ factors are total-even. Thus two root derivatives already annihilate the full form polynomial for each fixed assignment of the ordered color slots. Extracting distinct $d\eta_1,d\eta_2$ components from the whole trace includes both assignments of these form variables. Every positive-loop graph contains such a pair at the common odd coordinate, hence
$$
Q_n\operatorname{Tr}((dC)^p)=0\qquad(p\ge1,\ n\ge1).
$$
The replay uses independent commuting form variables $d\eta_1,d\eta_2$ and anticommuting $dz^1,dz^2,dz^3$. It checks all 96 ordered contractions on three distinct color factors before tracing, including the undifferentiated and all three first holomorphic jet assignments on each root. [Native replay](single_trace.py)

## Products

For two such coefficients $A,B$, each tree component of a marked graph has at most one spoke at $A$ and one at $B$, by the same two-slot calculation in the odd coordinate of that tree. If the component has $k$ spokes and $v$ interaction vertices, marked-graph sparsity and tightness give $v=2k-2$, $k\ge2$. Consequently every component has two vertices and two spokes, one at each root. All positive-loop kernels for this binary product are parallel three-edge paths and the marked edge. The holomorphic kernel, including arbitrary displacement insertions, is the previously evaluated path kernel. The superfield component derivatives on the two endpoints retain their native Berezin signs. [Graph count and exact path kernels](scalar_products/README.md)

The marked chain identity makes every iterated marked product of these descendants a quantum cocycle. Assign $c,g_1,g_2,b$ their two source weights $(0,0),(1,0),(0,1),(1,1)$. The coefficient of $\hbar^m$ in any binary product has $2m$ cubic vertices and $3m$ physical propagators, besides the marked edge. Each cubic Berezin integral selects weight $(1,1)$ across its incident fields. Each physical propagator removes one $(1,1)$ by pairing complementary components. The total output weight is therefore the sum of the input weights minus $(m,m)$. At total input weight $(\beta_1,\beta_2)$ every iterated marked product terminates by order $\min(\beta_1,\beta_2)$. This constructs lifts of all products. Whether a lifted product is an incoming boundary depends on the adjacent cohomology.

## A cubic sector at finite rank

Fix $\mathfrak{sl}_N$ and a holomorphic direction $i$. At coordinate weight one in $z^i$ and one in each odd coordinate, write $x=g_1$, $y=g_2$. In the cochains relative to constant colors the only possible cubic letters are $x,y,\partial_i c$. Each nonconstant mode costs at least one of the three units of source weight. For $N\ge3$ the invariant color contractions on three traceless matrices are the two cyclic orderings of a single trace. Indeed the invariant permutations of three defining and three dual matrix slots have only the two three-cycles after tracelessness. The two resulting contractions split into
$$
\operatorname{tr}([x,y]\partial_i c),\qquad
\operatorname{tr}((xy+yx)\partial_i c).
$$
Both are covariant traces and closed. All quadratic invariant boundaries use the bracket tensor. In particular $\widehat Q_0\langle b,\partial_i c\rangle=\langle[x,y],\partial_i c\rangle$. The second trace takes the value $-12\varepsilon$ on $x=y=\operatorname{diag}(1,1,-2,0,\ldots)$, $\partial_i c=\varepsilon x$, where $\varepsilon$ is a single odd parameter. All brackets vanish on this commuting probe. For $N=2$, $xy+yx=\operatorname{tr}(xy)\mathbf 1$, so only the boundary remains.

This relative calculation also gives absolute cohomology at the specified weight and length. Filter the split current algebra by the number of constant-color arguments. With $J$ the source augmentation ideal, the spectral sequence has groups $H^r(\mathfrak g,H^s(\mathfrak g\otimes J))$. At $r+s=3$, the $r=1,2$ groups vanish by the two Whitehead lemmas applied at the finite multidegree, and $r=3,s=0$ has weight zero. Each invariant ideal cocycle extends by zero on constant-color arguments, and restriction recovers it, so every surviving relative class is absolute. The symmetric trace is a coefficient of $\operatorname{Tr}((dC)^3)$, with one $dz^i$ and both even odd-coordinate forms at $\eta=0$. Its outgoing quantum differentials vanish. An incoming differential would start at field length two, where every class has a strict representative. The dimension on every page is one for $N\ge3$ and zero for $N=2$. [All quadratic classes](quadratic.md)

## All cubic local classes for simple gauge algebras

For $\mathfrak{sl}_N$ with $N\ge4$, the independent trace contractions through length four compute the stable absolute cohomology at length three. In positive coordinate weight its only cubic primitives are coefficients of $\operatorname{Tr}((dC)^3)$, whose three-color tensor is symmetric. The constant-coordinate class is $\operatorname{tr}(c^3)$. A product of positive field lengths adding to three would require a linear class, and $H^1(\mathfrak{sl}_N\otimes R)=0$. [Stable trace computation](README.md)

For $\mathfrak{sl}_2$, all invariant trilinear forms are multiples of $\langle[X,Y],Z\rangle$. Indeed the adjoint is the three-dimensional orthogonal representation, and the only cubic invariant tensor is its volume. In the cocycle equation on four colors the resulting tensors are
$$
\langle[X,Y],[Z,W]\rangle,\quad
\langle[X,Z],[Y,W]\rangle,\quad
\langle[X,W],[Y,Z]\rangle.
$$
Jacobi gives the respective coefficients $(1,-1,1)$ as a relation. In an orthogonal color frame each channel is the difference of two of the three pairings $\langle X,Y\rangle\langle Z,W\rangle$ and its permutations, times the same nonzero factor. These three pairings are independent, as testing two mutually orthogonal nonisotropic colors in each pair order shows. Thus Jacobi spans the relation space. It follows, for each current coefficient separately, that any cubic $\mathfrak{sl}_2$ cocycle extends to $\mathfrak{sl}_4$ by using the same graded-symmetric coefficient function with the four-dimensional defining trace. Under the diagonal block embedding its restriction is the original cocycle. In positive weight the stable cubic classes restrict to zero because $xy+yx=\operatorname{tr}(xy)\mathbf1$ on the two-dimensional block. Hence every positive-weight cubic cocycle of $\mathfrak{sl}_2$ is exact. At weight zero only the constant-color Cartan class remains.

The same bracket-tensor argument works for any finite-dimensional complex simple $\mathfrak g$. Restrict to a root $\mathfrak{sl}_2$, on which the invariant action pairing is nondegenerate. Every invariant quadratic primitive there has the form $\langle X,Y\rangle\psi(u,v)$ from [the current cochain calculation](quadratic.md), and extends back using the pairing on $\mathfrak g$. Hence a bracket-tensor cubic cocycle which becomes exact on that root is exact on $\mathfrak g$.

For completeness, the trilinear invariant multiplicity is the multiplicity of the adjoint representation in its tensor square. Write $\lambda$ for the highest root and $\rho$ for half the sum of positive roots. Multiplying the adjoint character by the Weyl character numerator of highest weight $\lambda$, the coefficient of the highest-weight numerator of $\lambda$ is
$$
\sum_{\substack{w\in W\\\mathrm{Weyl\ group}}}
(-1)^{\ell(w)}\operatorname{mult}_{\rm ad}
\big(\lambda+\rho-w(\lambda+\rho)\big).
$$
The identity contributes the rank. For a remaining contribution the difference must be a positive root $\alpha$. Preservation of the norm implies $\langle\lambda+\rho,\alpha^\vee\rangle=1$, forcing $\alpha$ to be a simple root orthogonal to $\lambda$. Since $\lambda+\rho$ is regular, only its simple reflection gives this difference. Thus the multiplicity is the number of nonzero Dynkin labels of the highest root. This equals two for $A_r$ with $r\ge2$ and one in the other simple root systems. In type $A$ the two independent trilinear forms are the antisymmetric bracket and the symmetric defining trace. For $\mathfrak{sl}_3$ the latter sector is already covered by stability at $N=$ field length. All other simple types have only the bracket class. Consequently for the full absolute cubic sector the positive-weight dimension vanishes outside type $A_r$, $r\ge2$. In that type its dual is the source two-forms modulo exact two-forms, as in the stable computation. At weight zero every simple factor gives its Cartan class. Every such cubic class survives, since the outgoing representatives are strict and degree two has no outgoing quantum differentials.
