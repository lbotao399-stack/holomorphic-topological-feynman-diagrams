# First-order deformations in the full translation-invariant five-dimensional theory

Consider the commutative five-dimensional Chern–Simons theory of the flagship manuscript with $\kappa(X,Y)=\operatorname{tr}(XY)$ for $\mathfrak{sl}_N,\mathfrak{gl}_N$. Through order $\hbar^2$, the full translation-invariant local BV complex gives the same necessary equation on the holomorphic cubic deformation as the flagship. For $\mathfrak{gl}_N$ it also gives a compatibility equation for the extra quadratic deformations. The argument uses the local renormalized BV notion of equivalence. [Flagship conventions and cubic equations](../flagship/flagship.tex) · [Regulated contact calculation and local homotopy theorem](../../deepen/h1/verify/regulated_contacts.md)

**Projection to the undressed cubic.** In the manuscript $U=\{\tau,\bar 1,\bar 2\}$ as an index set, and $\psi_S$ inserts the constant exterior factor $dx_{U\setminus S}$ for a subset $S$ of those three directions. Every cohomology class in the full translation-invariant complex has a representative which is a sum of $\psi_S(\phi)$ with holomorphic derivatives only. A tree bracket multiplies the two constant exterior factors. The component with no exterior factor can only come from two undressed inputs. Wang–Williams Proposition 3.21 restricts each connected local graph contribution to the diagonal after applying holomorphic derivatives to the test form. The exterior factors therefore multiply there as well. Proposition 3.23 gives zero for a one-loop graph with one canonical first-order insertion, including insertions with nonempty exterior factor. For undressed holomorphic vertices the surviving collapsed graph is itself undressed and uses only holomorphic derivatives. This excludes positive-weight local graph terms in the canonical heat-kernel representative. [Wang–Williams, Propositions 3.21 and 3.23](https://arxiv.org/pdf/2407.08667#page=32)

A general first-order change by the classical differential $Q_0K$ is handled by the local renormalized family over the unit interval. The first-order path is obtained by adding $(Q_0+d_{[0,1]})(sK)$, with $s$ the interval coordinate. Arity by arity, the relative lifting problem has a closed local obstruction in the product with interval forms; the relative de Rham homotopy at one endpoint preserves translation invariance. Consequently an order-two lift at one endpoint gives one at the other. The [local homotopy argument](../../deepen/h1/verify/regulated_contacts.md) includes the passage from the flat torus to flat space for field-dependent local terms. The [published lifting theorem](https://arxiv.org/pdf/0706.1533#page=53) permits renormalization of insertions containing nonholomorphic derivatives.

For every simple algebra in the flagship no-go theorem, $(\mathfrak g^{*\otimes3})^{\mathfrak g}=\mathbb C f$, $\mathfrak g\ne\mathfrak{sl}_2$, the undressed equation in bidegree $(2,2)$ still has no solution, even with arbitrary translation-invariant first-order terms. For $\mathfrak{sl}_{N\ge3}$ the $U(2)$-equivariant undressed coefficient still satisfies $a^2=\varepsilon U N^2/(256\pi^4)$, where $U=\pm1$ here is the manuscript's graph-orientation sign. All additional equivariant quadratic deformations for this simple algebra are allowed at this order. The identification of $a$ with a coupling in untwisted five-dimensional Yang–Mills requires the separate physical map in [the level calculation](J1_and_5d_CS_level.md).

**All $U(2)$-equivariant classes for the matrix algebra.** A component of $\psi_S H^{|S|}$ with topological and antiholomorphic weights $w=-\mathbf1_{U\setminus S}$ must have holomorphic bidegree $(1+w_{\bar1},1+w_{\bar2})$. When $S=U$, this is the cubic $j$. Omitting one antiholomorphic index leaves a two-cochain with one holomorphic derivative. Invariance gives a symmetric invariant bilinear form $B$; the two directions combine to give
$$
G_\epsilon(B)=\int d^2z\wedge
 \big(d\bar z_1\wedge B(\partial_{z_2}A,A)
       -d\bar z_2\wedge B(\partial_{z_1}A,A)\big).
$$
Omitting only $\tau$ would give a two-cochain at bidegree $(1,1)$, which vanishes by input exchange and symmetry of all invariant bilinear forms here. One-cochains have zero momentum modulo total derivatives, leaving only the manuscript's $T_\lambda$ for $\lambda$ a central covector when $S=\{\tau\}$. Degree zero is reduced. Thus $j$, $G_\epsilon(B)$ and $T_\lambda$ exhaust the equivariant first-order classes in this model.

On the traceless and central summands of $\mathfrak{gl}_N$, let $b_s,b_c$ be the scalar actions of $\kappa^{-1}B$. The polarized binary Hamiltonian $B(\partial_{z_i}A,A)$ acts on each plane-wave slot by $2p^i\kappa^{-1}B$ in the manuscript's shuffle ordering. This follows from the two terms $p^i-(-p^i)$. The corresponding constant derivatives commute on both summands, so their antisymmetrization with the exterior factors gives $\{G_\epsilon(B),G_\epsilon(B)\}=0$. The linear term shifts a central field at zero momentum. Its contraction with $G_\epsilon(B)$ integrates a derivative, and its contraction with $j$ gives a vanishing two-momentum determinant. Two copies of its constant exterior two-form wedge to zero. Thus $T_\lambda$ is unrestricted at this order.

For $N\ge3$ write $j=\alpha\pi^*D_{\mathfrak{sl}}+\chi M_1+\gamma S_0$ as in the flagship. Only $M_1$ contains one central and two traceless inputs. In a cubic supported purely on the center or purely on the traceless summand the action of the binary Hamiltonian gives the common scalar times total momentum, which vanishes. On the mixed component the remaining momentum factor is $2(b_c-b_s)p_{\rm center}^i$. In particular choose $X_0=\mathbf1_N$, $X_1=X_2=\mathrm{diag}(1,-1,0,\ldots,0)$, and the momenta $(1,1),(-1,0),(0,-1)$. Before multiplication by the background form and the common BV orientation, each holomorphic component of the mixed bracket has value
$$
2(b_c-b_s)p_0^i\,\chi[p_0p_1]\operatorname{tr}(X_0)\operatorname{tr}(X_1 X_2)
=4N\chi(b_c-b_s).
$$
All three matrices commute. Every Chevalley–Eilenberg boundary therefore vanishes on these waves. The two exterior components belong to distinct summands of translation-invariant cohomology. The additional obstruction vanishes precisely when $\chi(b_c-b_s)=0$. For nonzero $\chi$, the central and traceless fields must have the same coefficient of the derivative in $G_\epsilon$.

With $c_O=-U/(3072\pi^4)$ the three necessary and sufficient equations at $N\ge3$ are
$$
\frac{\alpha^2}{12N^2}-\frac{\chi^2}{48}=-\varepsilon c_O,
\qquad \chi(N\gamma-\chi)=0,\qquad \chi(b_c-b_s)=0.
$$
At $N=2$, for $j=x_0D+x_1 S_1$, the corresponding equations are $x_1(x_0+x_1)=0$ and $(x_0+x_1)(b_c-b_s)=0$. At $N=1$ only the central summand exists. The undressed equation and every mixed bracket vanish there. [Exact formal trace and Cartan replay](mixed_brackets.py)

## All-order lift with the quadratic and linear terms

The three necessary equations above are sufficient to all orders for the canonical equivariant representatives. Use the [formal local cubic lift](../../deepen/o3/higher_orders/center_classes.md#all-higher-orders-with-a-mixed-cubic). The Hamiltonian of $G_\epsilon(B)$ acts on the partial de Rham forms by the constant endomorphism $2\kappa^{-1}B$ times $d\bar z_1\wedge\partial_{z_2}-d\bar z_2\wedge\partial_{z_1}$, with the overall BV orientation of the manuscript. Denote this odd operator by $\delta$. For the free differential $D$ and gauge-fixing operator $D^*$ of the manuscript,
$$
\delta^2=0,\qquad D\delta+\delta D=0,\qquad
D^*\delta+\delta D^*=-8\kappa^{-1}B\big(\partial_{z_1}\partial_{z_2}-\partial_{z_2}\partial_{z_1}\big)=0.
$$
The flat Laplacian commutes with $\delta$. On fields, $K_t$ acts as the scalar heat operator and $P_\varepsilon^L=D^*\int_\varepsilon^L K_t\,dt$. Thus at finite cutoffs
$$
P_\varepsilon^L\delta P_\varepsilon^L=0,\qquad
K_t\delta P_\varepsilon^L+P_\varepsilon^L\delta K_t=0 .
$$
The second sign is the graded composition sign of the two heat-edge markings. Equivalently it results from the graded commutator of $D$ with $P\delta P$ and $[D,P]=K_\varepsilon-K_L$; the cancellation holds separately at each heat time. These identities act on arbitrary external partial forms. [Symbol calculation on every form component](dressed_kinetic.py)

A local contact graph with a bridge and positive loop number fails the Laman inequalities. After deleting the bridge, each component containing an edge must have at least $3L_i+1$ edges at loop number $L_i$. A component consisting of a single vertex has zero loops and no edges. Even with one such component the full edge count exceeds $3L+1$. Hence a quadratic vertex occurs with two internal edges in every surviving positive-loop contact, and a linear vertex cannot occur. With a quadratic vertex the terms whose heat edge is elsewhere factor through $P\delta P$. The two remaining heat markings sum to zero by the second identity; a self-loop contributes zero in the manuscript's gauge. Equivalently the perturbed differential $D+\hbar\delta$ has the same gauge-fixed Laplacian and propagator. [Local contact conditions](../../deepen/h1/verify/regulated_contacts.md)

The surviving tree contacts are BV brackets. Each cubic of the mixed solution is invariant under $\delta$ since $\chi\ne0$ forces $b_s=b_c$ and the common holomorphic derivative sums to zero by momentum conservation. For $\chi=0$ take the decoupled lifts on which the derivative has a common coefficient within each vertex. The Hamiltonian of $T_\lambda$ is a constant central shift with a constant antiholomorphic two-form. Its contraction with the commutator vertex vanishes by centrality. In an order-$k>0$ cubic, setting one momentum to zero kills $[p_0p_1]^k$ by three-point conservation. Its contraction with the quadratic integrates a derivative, and its self-bracket vanishes. Thus adjoining the original linear and quadratic terms to these cubic lifts solves the complete local equations for these representatives.

On a block where $\hbar\delta=t(d\bar z_1\wedge\partial_{z_2}-d\bar z_2\wedge\partial_{z_1})$, the coordinates $z_1+t\bar z_2,z_2-t\bar z_1$ are killed by the deformed antiholomorphic differential. This identifies the quadratic coupling as a constant deformation of the holomorphic directions. A nonzero mixed cubic requires the same deformation on the scalar and traceless fields. This interpretation is internal to the twisted theory.
