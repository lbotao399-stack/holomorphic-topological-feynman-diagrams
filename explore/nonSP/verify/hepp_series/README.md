# Rational series for a projective Schwinger period

For a connected graph, $U$ is its sum of cotree monomials $t_I=\prod_{e\in I}t_e$, one for each spanning-tree complement $I$. The input is a homogeneous rational polynomial $P$, with scale-invariant projective density $P/U^\alpha$ against the radial contraction of $\prod_e dt_e$. The exponent $\alpha$ can be any positive integer or half-integer. [Implementation](period.py)

**Convergence.** Write $a$ for an exponent vector occurring in $P$, $S$ for a proper nonempty edge subset and $b_1(S)$ for its loop number. Absolute convergence is equivalent to all the inequalities

$$
\sum_{\substack{e\in S\\\mathrm{edges}}}(a_e+1)>\alpha b_1(S).
$$

For necessity fix one time outside $S$ to one, and scale the times in $S$ as $r$ times positive ratios, with one ratio fixed to one. The order of $U$ at $r=0$ is $\min_I |I\cap S|=b_1(S)$ by extending a maximal forest of $S$ to a spanning tree. The lowest-order part of $P$ is a nonzero polynomial in the ratios and outside times. Its inner and outer total degrees are fixed, so the two fixed coordinates preserve nonvanishing. Its absolute value is bounded below on some compact box of positive ratios with nonempty interior. The radial measure contributes $r^{|S|-1}dr$, yielding the strict inequality for the least degree in $S$ and hence for every monomial. Conversely, in each ordering of the edge times use the successive ratios in $(0,1)$ and fix the largest time to one. The dominant cotree product $\widehat U=\max_I t_I$ satisfies $\widehat U\le U\le (\#\mathrm{cotrees})\widehat U$. On each initial edge subset $S$ of the ordering, its cotree $I$ has $|I\cap S|=b_1(S)$. Each ratio integral converges by the displayed inequality. This proves termwise absolute integrability even when coefficients have both signs.

**Monomial moments.** For fixed exponent vector $a$ and denominator exponent $\beta$, with $\sum_e(a_e+1)=\beta b_1(\Gamma)$, set $H(\varnothing)=1$. The moment with denominator $\widehat U^\beta$ is $H(E)$ from the recursion

$$
H(S)=\frac{\sum_{\substack{e\in S\\\mathrm{edges}}} H(S\setminus\{e\})}
{\sum_{\substack{e\in S\\\mathrm{edges}}}(a_e+1)-\beta b_1(S)},\qquad
H(E)=\sum_{\substack{e\in E\\\mathrm{edges}}} H(E\setminus\{e\}).
$$

The first equality holds for proper nonempty subsets. Removing the last edge sums the linear extensions one edge at a time. For the cone where a given cotree $I$ dominates, retain just the subsets with $|I\cap S|=b_1(S)$. Equivalently, every edge on the complementary tree path of a chord must precede that chord. The unrestricted recursion sums all cones at once. When a proper subset splits as two disjoint matroid restrictions with additive loop number, the recursion gives the product of their moments by induction. The implementation caches such factors on their own edges and can combine exponents under graph automorphisms. The symmetry subgroup on a specified cone must stabilize its cotree.

**Series.** Put $x=U/\widehat U$ and choose a rational $s\ge\sqrt{\#\mathrm{cotrees}}$. With $r=(s-1)/(s+1)$ and $C=(s+1)^2/4$, write $x=C(1-2ru+r^2)$, so $|u|\le1$. For $r=0$ only the constant term is required. For $r>0$,

$$
x^{-\alpha}=C^{-\alpha}\sum_{\substack{k\ge0\\\mathrm{Gegenbauer\ series}}}r^kC_k^{(\alpha)}(u),\qquad
(1-2ru+r^2)^{-\alpha}=(1-r e^{i\theta})^{-\alpha}(1-r e^{-i\theta})^{-\alpha},\quad u=\cos\theta.
$$

The binomial coefficients in the two factors are positive, giving $|C_k^{(\alpha)}(u)|\le (2\alpha)_k/k!$; $(2\alpha)_k$ is the rising factorial. Each term $r^k C_k^{(\alpha)}(u)$ is a rational polynomial in $x$. Its $x^j$ moment uses the polynomial $P U^j$ and the preceding recursion with $\beta=\alpha+j$. If $B$ is the sum of the absolute numerator coefficients times their $\widehat U^{-\alpha}$ moments, truncation through degree $n$ has absolute error at most

$$
B\left(1-C^{-\alpha}\sum_{\substack{0\le k\le n\\\mathrm{binomial\ series}}}\frac{(2\alpha)_k}{k!}r^k\right).
$$

Here $C^{-\alpha}=(1-r)^{2\alpha}$, so the bound and partial sums are rational. For every $r<r'<1$ the bound is at most $B C^{-\alpha}(r/r')^{n+1}(1-r')^{-2\alpha}$. The optional direct binomial expansion uses $1-U/(C\widehat U)$ with a rational center $C\ge(1+\#\mathrm{cotrees})/2$ and the same moment recursion.

```python
from fractions import Fraction
from period import HeppSeries

period = HeppSeries(edges, Fraction(11, 2), numerator)
for order, term, partial, error in period.terms(2, permutations='graph'):
    print(order, partial - error, partial + error)
```

The numerator maps full edge-exponent tuples to rational coefficients. Passing a cotree bitmask as the second argument computes its cone only. The default expansion is `gegenbauer`; `method='binomial'` selects the direct expansion. Python-flint accelerates rational arithmetic. [Independent sector and integral checks](check.py) · [Four-loop skeleton support restrictions](skeleton_support.py)

## Five-dimensional four-loop cyclic probe

For the external matrix cycle and momenta $(1,0),(0,1),(1,1),(-2,-2)$ in the [order-four cochain equation](../../../../deepen/o3/higher_orders/README.md), the bare geometric value preceding the common physical normalization is

$$
\frac{\pi^2 N^4}{10500}
+\sum_{\substack{(c,m)\\\mathrm{table\ below}}}N^m
\sum_{\substack{k\ge0\\\mathrm{Gegenbauer\ series}}}
\int \frac{P_{c,m}}{\widehat U_c^{11/2}}\,C_c^{-11/2}r_c^k C_k^{(11/2)}(u_c).
$$

The first term includes the three series–parallel graphs. The remaining $c$ labels refer to the zero-based graph ranges. Each $P_{c,m}$ is the combined degree-thirteen numerator at rank power $N^m$. Its stored coefficients multiply orbit-averaged monomials under the edge permutations preserving $U_c$. Each graph's symmetry factor is already included in the source numerator. [Series–parallel evaluation](../../../../deepen/o3/higher_orders/four_loop_sp/README.md) · [Numerator construction](../../../../deepen/o3/higher_orders/four_loop_reduction/README.md)

| Graphs $c$ | Rank powers $m$ | Orbit counts in rank order | Cotrees | $r_c$ |
|---|---|---|---|---|
| 3–10 | 0, 2, 4 | 428, 435, 247 | 56 | $13/17$ |
| 11–38 | 0, 2, 4 | 2720, 2975, 624 | 75 | $123/155$ |
| 39–49 | 0, 2 | 780, 817 | 81 | $4/5$ |

The [exact replay](cyclic_series.py) computes $k=0,1,2,3$ for all eight series. Each row records the rational coefficients, partial sum and absolute remainder bound. Multiplying a row's bound by $|N|^m$ and adding gives a bound for the full bare probe. The early bounds are too wide to decide the sign of that sum. [Computed rational terms and bounds](cyclic_order3.jsonl) · [Independent check summing all 47 individual graph files](verify_aggregates.py)

```sh
python3 verify_aggregates.py
python3 cyclic_series.py --order 3
```

At the separate geometric source assignment $((2,1),(1,4),(-2,3),(-1,-8))$, graph 3 has 620 strictly negative numerator coefficients. Every monomial satisfies the strict margins on all 510 proper nonempty edge subsets. The graph period is strictly negative. A two-sided check using its first four raw moments gives the integer enclosure $-247169<K_3<-3$. For the positive measure obtained by reversing the sign, interpolate $x^{-11/2}$ at two positive rational-square double nodes for a lower bound, and at both endpoints and one interior rational-square double node for an upper bound. The fourth derivative is positive, so the respective Hermite remainder products have the required signs. Nodes near the Gaussian and Lobatto nodes improve the bounds, and the sign argument applies to any distinct interior nodes. [Graph numerator](../../../../deepen/o3/higher_orders/four_loop_reduction/reconstruct_graph3.json) · [Period replay](skeleton_period.py) · [Rational Hermite bounds](moment_bounds.py)
