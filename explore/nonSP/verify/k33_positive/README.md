# Four-source period of the complete bipartite graph

Take the project edge order in [graphs.py](../../graphs.py), zero displacements, and four sources at $0,1,2,3$ whose sum vanishes. The bipartition is $\{0,2,4\}\sqcup\{1,3,5\}$. Write $\langle i,j\rangle=\det(\lambda_i,\lambda_j)$ in the project orientation. The four-loop kernel is
$$
K_{K_{3,3}}=c\,\langle0,2\rangle\langle1,3\rangle
              \big(\langle1,3\rangle^2-\langle0,2\rangle^2\big),
\qquad
\frac{6252279056460127}{8925805535599677505536000}<c<
\frac{6252279056460127}{31603654656000}.
$$

**Source dependence.** The determinant formula is homogeneous of degree eight in the sources and transforms by $\det(g)^4$ under $g\in GL_2$. The invariant polynomials on three independent two-vectors are generated freely by their three pairwise determinants. Apply this to $\lambda_0+\lambda_2,\lambda_0-\lambda_2,\lambda_1-\lambda_3$. The three determinants, up to nonzero constant factors, are $\langle0,2\rangle$, $\langle1,3\rangle$, and $\det(\lambda_0-\lambda_2,\lambda_1-\lambda_3)$. The swaps $(0\,2)$, $(1\,3)$ and $(0\,1)(2\,3)(4\,5)$ induce odd permutations of the nine edges. Each changes the kernel's sign. In a degree-four monomial the first two swaps require odd powers of both $\langle0,2\rangle$ and $\langle1,3\rangle$. The third exchanges these determinants and negates the third determinant. Antisymmetrization leaves the displayed polynomial. Switching two vanishing sources on the same side also proves the complementary four-source sector vanishes. [Sign conventions, §2.2](../../../../code/htk/CONVENTIONS.md)

**The nine-variable numerator.** Use $\lambda_0=(1,0)$, $\lambda_1=(0,1)$, $\lambda_2=(1,1)$ and $\lambda_3=(-2,-2)$, so the displayed polynomial gives $K=6c$. In the project notation $M=CtC^{\mathsf T}$, $U=\det M$, and $D(t)$ is the Schwinger density. The spanning-tree expansion of $U$ has 81 square-free degree-four monomials. The numerator of each edge current $U y_e^i$ is the same spanning-tree sum with each term multiplied by the current on the corresponding tree. This follows by Cramer's rule for the grounded Kirchhoff matrix, or by checking that the tree sum satisfies the vertex and cycle equations.

The product $U^6D$ has total degree fifteen and degree at most three in each time. To see polynomiality, vary one time with the others fixed at generic positive values, so $U$ is affine and has a simple root. The pole in $y$ comes from $M^{-1}=\operatorname{adj}(M)/U$. At the root the adjugate has rank one. Consequently the two edge arrays $Uy^1,Uy^2$ there are proportional, allowing one block of four determinant rows to be cleared by the other. If either block vanishes use that block directly. The eight-row current determinant with denominators cleared therefore has a factor $U^4$. The formula for $D$ contributes $U^{-10}$ before this cancellation. Its projective pivot may be chosen different from the varied edge, avoiding a pole at zero. At positive infinity all currents stay bounded by their finite spanning-tree values, and using the varied edge as pivot gives $D=O(t_e^{-3})$. Hence $U^6D$ is separately polynomial of degree at most three. Its total degree is $6\cdot4-9=15$. Integer-grid interpolation on $\{1,2,3,4\}^9$ therefore recovers every coefficient. [Exact tree sums and interpolation](check.py) · [Integer determinant evaluator](grid.cpp)

**Sign of the integral.** Sum $U^6D$ over the 72 graph permutations acting on the time variables. Its 10,980 monomials form 175 orbits. Express the sum in the basis of orbit-averaged monomials. For an exponent vector $h$ of total degree eleven and edge $i$, the exact divergence used in the certificate has numerator
$$
t^h\big((h_i+1)U-5t_i\partial_{t_i}U\big)
$$
over $U^6$. Indeed it is $\partial_{t_i}(t_i t^h/U^5)$. Fix the projective scale with some edge other than $i$. The integral in $t_i$ has zero boundary contribution. At zero the primitive vanishes since $h_i\ge0$. At infinity it vanishes since the recorded fields have $h_i\le2$ and every edge occurs in some spanning-tree complement. For absolute convergence of the expanded terms, the code verifies for every proper nonempty edge subset $S$ and every exponent $a$ used that
$$
\sum_{\substack{e\in S\\ \mathrm{edges}}}(a_e+1)>6\min_{\substack{T\\ \mathrm{spanning\ trees}}}|S\setminus T|.
$$
In each ordered Schwinger sector these are the strict convergence inequalities in the successive time ratios. This also justifies integration in $t_i$ first.

After multiplying the symmetrized numerator by ten and subtracting the 222 orbit-averaged divergences in the certificate, there are 399 strictly positive orbit coefficients. Their least value is 1103 and their sum is 1262466. The integral of this positive numerator over $U^6$ equals $4320c$, with $4320=10\cdot72\cdot6$. [Rational divergence coefficients](fields.json) · [Exact replay](check.py)

**Rational bounds.** In each ordered sector the largest spanning-tree monomial $m$ satisfies $m<U\le81m$ in its positive interior. Replacing $U$ by $m$ in the positive integral for $c$ gives $6252279056460127/31603654656000$. For a monomial with exponents $a$, this sector integral is the product of reciprocals of the positive margins displayed above, taken along the chain of proper initial subsets ordered from smallest to largest time. A greedy spanning tree minimizes $|S\setminus T|$ simultaneously on this chain. The verifier sums the products over all orderings by recursion on the final edge of each subset. Replacing $U$ by $81m$ gives the stated lower bound, strict off the equal-times locus.
