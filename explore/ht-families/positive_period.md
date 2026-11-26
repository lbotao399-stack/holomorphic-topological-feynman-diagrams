# Independence of the series excess at even loop number

Let $n>0$ be the even loop number of a series–parallel graph, and $U$ its spanning-tree complement polynomial. For every reduction, $E$ is the sum over its series nodes of half the number of odd-loop children, rounded down at each node. If the cycle matroid is connected, $E$ has the same value for every reduction, including changes of terminal pair. The [period calculation](README.md) supplies the series and parallel exponents of $U$.

**A positive comparison integral.** Let $J$ run over complements of spanning trees. Each $J$ contains $n$ edges and contributes the monomial $\prod_{e\in J}t_e$ to $U$. For a set $S$ of edges, write $r(S)=\max_J|S\cap J|$. Connectedness of the matroid means $r(S)+r(S^c)>n$ for every proper nonempty $S$. Consequently

$$
\min_J |S\cap J|=n-r(S^c)<r(S).
$$

Choose any set $A$ of $n/2$ edges. With $m$ the number of complements of spanning trees, take

$$
p=(n+1)m+\tfrac12,\qquad
\alpha_e=(n+1)\#\{J:e\in J\}+\mathbf1_{e\in A}.
$$

These integer edge exponents satisfy $\sum_{e\in\mathrm{edges}}\alpha_e=pn$. In the gauge $t_1=1$, integrate the positive density

$$
\frac{\prod_{e\in\mathrm{edges}}t_e^{\alpha_e}}{U^p}
\left|\bigwedge_{e\ne1}d\log(t_e/t_1)\right|.
$$

The rational-power factor is scale invariant, and changing the edge used to fix scale has logarithmic Jacobian of absolute value one. On an ordered sector, let $S$ be the first $j$ edges in order of decreasing time and let $x=t_{j+1}/t_j$. The greedy basis has exactly $r(S)$ edges in $S$ for every prefix. Its monomial in $U$ bounds $U$ below. The resulting upper bound for the density uses $x^{p r(S)-\sum_{e\in S}\alpha_e}\,dx/x$ for this ratio. Since at least one basis has fewer than $r(S)$ edges in $S$,

$$
p r(S)-\sum_{e\in S}\alpha_e
=(n+1)\sum_{J\in\mathrm{tree\ complements}}\big(r(S)-|S\cap J|\big)
+\tfrac12r(S)-|A\cap S|
\ \ge\ n+1-n/2>0 .
$$

Thus the integral exists and is strictly positive.

**Evaluation in any reduction.** For a node with child resistances $R_i$, set $b_i=R_i/R$ at a series node and $b_i=R/R_i$ at a parallel node, where $R$ is the parent resistance. In both cases $\sum_{i\in\mathrm{children}}b_i=1$. For $k$ children,

$$
\left|\bigwedge_{i=1}^k d\log R_i\right|
=\left|d\log R\right|
\frac{\left|db_1\cdots db_{k-1}\right|}{\prod_{i=1}^k b_i}.
$$

Indeed, subtract the last log-resistance row from the others. Up to a common sign, the remaining derivative matrix has entries $\delta_{ij}/b_i+1/b_k$ for $i,j<k$, with determinant $1/\prod_i b_i$. Iterating over the tree removes the root scale. In a child with loop number $L_i$, write $\alpha(i)$ for the sum of the edge exponents in the child. Its parameter in the resulting Dirichlet integral is

$$
\alpha(i)-pL_i\quad\text{at a series node},\qquad
p(1+L_i)-\alpha(i)\quad\text{at a parallel node}.
$$

Omit one-child blocks, which have no integration variables. At every other block, absolute convergence of this positive monomial density on the product of open simplices gives strictly positive parameters. Since $p$ is half-integral and nonintegral, the nonintegral parameters occur exactly at the odd $L_i$ of a series node and the even $L_i$ of a parallel node. The product of Dirichlet factors therefore lies in $\pi^{n/2+2E}\mathbb Q_{>0}$. The integral is identical in every reduction. Transcendence of $\pi$ gives equality of $E$.

**Four-loop check.** On the triangle with one triangular ear on each side, take the project's edge order

$$
(01,02,04,05,12,13,14,34,45).
$$

For $p=5/2$ the weights $(2,1,1,1,1,1,1,1,1)$ already give minimum $p r(S)-\sum_{e\in S}\alpha_e=1/2$ over the 510 proper nonempty edge subsets. Terminal pairs $(0,1)$ and $(0,2)$ each give the positive comparison period $\pi^4/6$ and $E=1$. This scalar period separates reduction independence from the cancellation of the physical kernel. [Exact replay](verify/positive_period.py)

**Split cycle matroid.** When $r(S)+r(S^c)=n$ for some proper nonempty $S$, choose a cycle basis adapted to the two column spans so that the period matrix $M$ is block diagonal. Multiplying all times in $S$ by the same factor rescales only its block. The Kirchhoff currents are constant along this change, since the two block equations $M\ell=-C\,\mathrm{diag}(t)y^0$ scale separately. In contraction with the corresponding scaling tangent, the matrix two-form $dM\,M^{-1}\wedge dM$ vanishes because the affected diagonal block gives $M M^{-1}dM-dM M^{-1}M=0$, and all other terms vanish. Thus the physical form also annihilates this scaling tangent. Its projection to projective space is nonzero. Under the degree condition of the kernel a top form annihilating a nonzero tangent is zero pointwise.
