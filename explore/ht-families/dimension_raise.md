# Raising the holomorphic dimension by subdividing a cotree

Let $\Gamma$ have $N$ edges and $n$ loops, with $N-1=n(H+T)$ and $T=0$ or $1$. Use the [project edge orientation and cycle matrix](../../code/htk/CONVENTIONS.md). A cotree $I$ is the complement of a spanning tree. Order $I$ by original edge order and let $C_I$ be the corresponding square minor of the cycle matrix $C$. Then $\det C_I=\pm1$.

For each cotree, subdivide its edges once, keeping the first segment in the original edge slot and appending the second segments in increasing $I$ order. The new vertices occur in that same order immediately before the original root. Retain the old sources in the first $H$ directions. In the new holomorphic direction give every inserted vertex source $+1$ and the original root source $-n$. Write $\Gamma_I$ for this graph with these sources. At zero displacement and zero topological source,
$$
K_\Gamma
=(-1)^{n(H+T+1)}
 \sum_{\substack{I\ \mathrm{cotree}\\\mathrm{of}\ \Gamma}}
  \det(C_I)\,K_{\Gamma_I}.
$$
A fixed polynomial in **old** currents can be inserted on both sides, using either segment of each subdivided edge. A displacement in the old directions is retained by placing its factor on either segment.

**Subdivided form.** Set the times of the first and second segment at $e\in I$ to $(1-x_e)t_e$ and $x_e t_e$, respectively. An uncalibrated cycle matrix of the new graph is $(C\mid C_I)$ because the two segments carry the same cycle flow. Its period matrix, its differential and its Pfaffian coincide with the original ones, including for variable $x_e$. In an old holomorphic component, the inserted vertex has zero source and the two edge currents coincide with the original edge current. In the added component Kirchhoff conservation gives $y_{\rm second}-y_{\rm first}=1$. In the source form $\psi_r$ the coefficient of $dx_e$ is therefore $C_{re}t_e$. All the old forms already fill the $N-1$ directions of the original projective domain. The new component contributes only its $n$-vertical part, equal to
$$
\det(C_I)\left(\prod_{e\in I}t_e\right)\bigwedge_{e\in I}dx_e .
$$
The extra factor $U^{-1}$ comes from the added holomorphic direction, where $U=\det(C\,\mathrm{diag}(t)\,C^{\mathsf T})$. The map from the old times followed by the cube coordinates to the new ordered times has positive full Jacobian $\prod_{e\in I}t_e$ and preserves the radial orientation. The vertical exterior product above already accounts for the change of variables.

**Calibration.** For the original tree $T_{\rm tree}$ complementary to $I$, its signed incidence determinant is $\sigma_{T_{\rm tree}}=\det C_I$. In the new graph all original edge slots form a tree, with the new vertices as leaves. Its incidence determinant equals $(-1)^n\sigma_{T_{\rm tree}}$. Indeed, grouping the original tree rows first, the new leaf block is $-1_n$; the inverse row permutation supplies the original tree–cotree sign. The new cotree consists of the appended segments, whose minor in $(C\mid C_I)$ is $\det C_I$. Calibrating the new cycle basis thus requires determinant $(-1)^n$. The full form has $H+1+T$ cycle-orientation factors. Integrating the $n$ unit intervals gives
$$
K_{\Gamma_I}
=(-1)^{n(H+T+1)}\det C_I
\int_\Gamma \frac{\prod_{e\in I}t_e}{U}\,\Omega_\Gamma .
$$
The projective domain in this integral is the parent domain. Cauchy–Binet gives $U=\sum_{I\ \mathrm{cotree}}(\det C_I)^2\prod_{e\in I}t_e$, so the positive weights in this identity sum to one pointwise. This proves the signed sum. Each child integral converges absolutely whenever the parent form is absolutely integrable. [Signed exact checks in both topological dimensions](verify/dimension_raise.py)

**Sparsity.** The child has $N+n$ edges and $n$ loops. To check a connected subgraph containing a cycle, prune its dangling segments and suppress fully present pairs of subdivided edges. The resulting parent subgraph of loop number $l$ contains at least $l$ edges of $I$, since the other edges belong to a tree. The child therefore contains at least $l$ more edges than this parent subgraph. The parent bound $|I_S|\ge(H+T)l+1$ gives the child bound $|I_S'|\ge(H+T+1)l+1$. Trees containing an edge satisfy this bound with $l=0$. This proves preservation of the project’s admissibility condition.

**Consequences.** Subdivision preserves a given series–parallel reduction and the value of $E$ in that reduction, since the newly inserted series nodes have only loop-free children. Starting with the [four-loop $H=2$, $T=1$ graph](h2_odd_triangle.md), the signed sum implies the existence of an admissible four-loop graph with $E=1$ and a nonzero kernel in every holomorphic dimension $H\ge2$. Each such lift still consists of three one-loop, two-terminal subnetworks around a triangle. At $H=1$ these kernels vanish by [transgression](odd_triangle.md).

The same identity applies to a non-series–parallel graph at $T=0$. The [three-loop $K_4$ subdivision](../nonSP/README.md) with rational sources has coefficient $(39-32\zeta(3))/67200$. Since $\zeta(3)$ is irrational, the signed sum must contain at least one irrational child coefficient with rational sources. Iteration gives such an admissible three-loop subdivision of $K_4$ in every $H\ge3$. This conclusion concerns the geometric graph kernel; a theory's color factors and interaction tensors may impose further cancellations.

At four loops start with $K_{3,3}$ at $H=2$, $T=0$, with the six conserved integer momenta $((-3,-1),(-2,2),(0,-1),(1,-1),(3,1),(1,0))$. Its bare geometric period is $-789/640+59\zeta(3)/80$, which is irrational by Apéry's theorem. Subdivision preserves integral momentum since each added coordinate has entries $1$ at the four inserted vertices, $-4$ at the original root, and zero elsewhere. All signed coefficients are integers. If all child periods were rational, so would be the parent. Induction therefore gives an admissible four-loop subdivision of $K_{3,3}$ with an irrational period for every $H\ge2$ at $T=0$. [Exact base period and replay](../nonSP/verify/k33_forcer/README.md#arbitrary-conserved-sources)

For the prism take the four nonzero sources $((1,0),(0,1),(1,0),(-2,-1))$ at vertices $0,1,2,3$. In the four-source bracket expansion $[01]=1$, $[02]=0$, $[12]=-1$. The [prism coefficients](../nonSP/results_nonSP_C2.json) give

$$
K_{\rm prism}=c_{400}-c_{301}+c_{202}-c_{103}
=-\frac{43}{20160}+\frac{\zeta(3)}{560}.
$$

Apéry's theorem and the same signed sum therefore also give an admissible four-loop prism subdivision with an irrational integral-source period in each dimension $H\ge2$, $T=0$.
