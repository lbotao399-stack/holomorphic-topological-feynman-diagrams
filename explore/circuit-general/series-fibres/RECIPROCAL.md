# Reciprocal core reconstruction with a rigorous stopping certificate

**[drv]** The identities and coefficient bounds below are proved
algebraically. **[chk]** `test_reciprocal.py` compares complete coefficient
dictionaries with the original exact Fraction numerator, including every
choice of projective gauge. The higher-loop period run uses these bounds
to certify all reconstructed coefficients; agreement of successive
reconstructions is not the stopping criterion.

## 1. The polynomial to reconstruct

Keep the notation of [DERIVATION.md](DERIVATION.md): $n$ loops, $H$
holomorphic directions, $N=nH+1$ original edges, and $m$ maximal chains
with lengths $L_c$. The original polynomial is

$$P(t)=U(t)^{n+H}D(t),\qquad \deg P=n^2-1.$$

The series core has $v=m-n+1$ vertices. This implementation requires that
it be connected and loopless; parallel core edges are allowed. Integrate
the chain fractions **without** the series-coordinate Jacobian:

$$
\overline P(R)=\int_{\prod_c\Delta^{L_c-1}}
P(t_{c,j}=R_cx_{c,j})\,dx.
$$

The rank bound in the previous note gives total degree at most $n-1$
on every complete chain. Thus

$$\deg\overline P=n^2-1,\qquad\deg_{R_c}\overline P\le n-1.$$

Define the conductance polynomial immediately before using it:

$$
\boxed{Q(g)=\left(\prod_cg_c^{n-1}\right)
\overline P(1/g).}
$$

It has no negative powers, and its total degree is

$$
\boxed{\deg Q=m(n-1)-(n^2-1)=(n-1)(v-2).}
$$

This is a polynomial transformation for coefficient reconstruction. We do
not change the outer period integration variable to $g$, so no inverse
Jacobian is included in $Q$.

For a six-loop $K_5$ core, $n=6$, $m=10$, $v=5$:

$$\deg\overline P=35,\qquad\deg Q=15,\qquad\deg_{g_c}Q\le5.$$

If $\beta$ is an exponent vector of $Q$, its contribution to the original
resistance integral has powers

$$
\boxed{a_c=n-1-\beta_c+L_c,\qquad \sum_ca_c=n(n+H).}
$$

The $L_c$ term includes both the series Jacobian $R_c^{L_c-1}$ and the
shift from a polynomial exponent to a period power $a_c-1$.

## 2. Graphic-rank support after inversion

For an edge subset $T$ of the core, its graphic rank is the number of
vertices minus the number of connected components of the graph with edge
set $T$; isolated vertices are counted. Its loop number is
$h_T=|T|-\operatorname{rank}_{\rm gr}(T)$.

For every proper core edge subset $T$,

$$
\boxed{\deg_TQ\le n\operatorname{rank}_{\rm gr}(T)-|T|.}
$$

**Proof.** Put $S=E_{\rm core}\setminus T$, which is nonempty. The
original cycle-rank bound, applied to the union of the full chains in
$S$, gives

$$\deg_S\overline P\le n\operatorname{rank}C_S-1.$$

For a monomial $R^\alpha$ of $\overline P$,
$\beta_c=n-1-\alpha_c$, and $\sum\beta=(n-1)(v-2)$. Also

$$\operatorname{rank}C_S=|S|+1-\#\operatorname{components}(T).$$

Substitution yields

$$
\begin{aligned}
\sum_{c\in T}\beta_c
&\le (n-1)(v-2)-(n-1)|S|+n\operatorname{rank}C_S-1\\
&=n\operatorname{rank}_{\rm gr}(T)-|T|.
\end{aligned}
$$

A graphic flat is an edge set containing every available core edge
whose endpoints are already connected in that set. Enlarging $T$ to its
flat does not change its graphic rank and only strengthens this bound.
If its closure is the entire graph, the inequality follows from total
degree: for proper $T$, its right side exceeds or equals
$\deg Q$, with difference $m-1-|T|$. Hence it suffices to enumerate proper
flats, together with the individual degree and homogeneity conditions.

For $K_5$ there are **52** flats, including the empty and whole sets.
The implementation retains **25** nontrivial proper-flat inequalities:

| Proper flat | Number | Bound |
|---|---:|---:|
| Triangle | 10 | $\deg_TQ\le9$ |
| Triangle and disjoint edge | 10 | $\deg_TQ\le14$ |
| $K_4$ on four vertices | 5 | $\deg_TQ\le12$ |

The empty and whole-set inequalities are not imposed. The homogeneous
degree is supplied separately.

### Dehomogenization and interpolation

Fix one $g_c=1$. An inequality involving its exponent remains valid after
dropping that nonnegative exponent. All resulting upper bounds are
downward closed, so the existing Newton lower-set interpolation applies
without changing its transform.

For the $K_5$ core, the counts are **[chk]**:

| Reconstruction support | Points |
|---|---:|
| Previous original-graph lower set, degree 35 | 1,788,632,625 |
| Reciprocal core, degree 15 and individual cap 5 | 877,844 |
| Reciprocal core with all flat bounds | 683,834 |
| Homogeneous exponent vectors passing all bounds | 518,124 |

Every one of the ten conductance gauges gives **683,834** points. The
homogeneous support mask is checked after reconstruction: coefficients
outside it must be exactly zero. The larger downward set is necessary
for this ordinary Newton transform; no unsupported deletion of lower
interpolation points is made.

## 3. Twenty-seven exact fibre evaluations at six loops

The chain-total bound also gives a separate ratio-degree bound on each
chain:

$$\deg_{x_c}P(Rx)\le n-1.$$

For the three length-two chains of an admissible six-loop $K_5$
subdivision, the integrand has degree at most five in each interval
coordinate. The three-point Gauss rule on $[0,1]$ is

$$
\int_0^1f(x)\,dx=
\frac5{18}f\left(\frac12-\frac{\sqrt{15}}{10}\right)
+\frac49f\left(\frac12\right)
+\frac5{18}f\left(\frac12+\frac{\sqrt{15}}{10}\right),
\qquad\deg f\le5.
$$

**Exactness.** Write $u=x-1/2$. Odd moments vanish by reflection. The
displayed weights reproduce the three even moments of degrees $0,2,4$:

$$
\frac5{18}+\frac49+\frac5{18}=1,\qquad
2\frac5{18}\frac{15}{100}=\frac1{12},\qquad
2\frac5{18}\left(\frac{15}{100}\right)^2=\frac1{80}.
$$

These are exactly $\int_0^1u^kdx$ for $k=0,2,4$. The product rule
therefore integrates the complete three-interval polynomial in
$3^3=27$ evaluations. Each value still contains the full contracted
determinant and all original sources.

For modular computation, choose a prime for which $15$ is a nonzero
quadratic residue and reduce the displayed algebraic identity in that
finite field. Either square root gives the same rule, because the two
off-centre nodes are exchanged. The weights and evaluations are exact
finite-field elements. There is no floating-point quadrature or numerical
error estimate.

This uses $683834\cdot27=18,463,518$ compressed determinant evaluations
per prime. The sampler reuses the core inverse at the fibre nodes; the
node count is not a count of separate complete circuit inversions.

## 4. Exact denominator of the integrated coefficients

First suppose the source components are integers. The numerator
determinant is integral, and its divisor $t_{e_0}U^{n(H-1)}$ is primitive
over the integers. The polynomiality proof in
[DERIVATION.md](DERIVATION.md) and Gauss's lemma imply

$$P\in\mathbb Z[t].$$

On a length-two chain, every barycentric monomial has degree at most
five. Its integral is

$$\int_0^1x^a(1-x)^b dx=\frac{a!b!}{(a+b+1)!},\qquad a+b\le5.$$

All denominators divide $\operatorname{lcm}(1,2,3,4,5,6)=60$. Hence for
three such chains,

$$\boxed{216000\,Q\in\mathbb Z[g],\qquad216000=60^3.}$$

The certificate computes the corresponding denominator for arbitrary
chain lengths. For a chain simplex of dimension $b$, it takes the least
common multiple of the denominators of all degree-$n-1$ barycentric
moments. Lower-degree monomials are integer linear combinations of these,
after multiplication by $(\sum_jx_j)^{n-1-k}=1$. The product of the chain
denominators is therefore a valid global denominator.

For rational sources, let $d$ be a common source denominator. The
polynomial has source degree $nH$, so the returned coefficient
denominator is multiplied by $d^{nH}$. The height calculation below is
then performed on the integer sources $d\lambda$.

## 5. A rigorous coefficient height and CRT stopping rule

Write $\|F\|_1$ for the sum of absolute coefficients of a polynomial.
The multiplicative Mahler measure is

$$
M(F)=\exp\left(\int_{[0,1]^N}
\log|F(e^{2\pi i\theta_1},\ldots,e^{2\pi i\theta_N})|\,d\theta\right).
$$

The inequalities used here are

$$
\|F\|_1\le2^{\sum_e\deg_{t_e}F}M(F),\qquad
M(F)\le\|F\|_1,\qquad M(FG)=M(F)M(G).
$$

For a nonzero integer polynomial, $M(F)\ge1$.
These are the coefficient, length and product bounds in K. Mahler,
*On some inequalities for polynomials in several variables*,
Journal of the London Mathematical Society **37** (1962), 341–344,
[author's paper reproduced by EMS](https://ems.press/content/book-chapter-files/27423),
equations (3)–(5) and Section 5. **[src]**

**Short proof of the bound needed here.** For a univariate polynomial,
factorization into roots gives
$|a_j|\le\binom dj M(F)$. Apply this coefficient bound one variable at a
time and integrate its logarithm in the remaining circle variables.
Iteration gives
$|a_\alpha|\le M(F)\prod_e\binom{d_e}{\alpha_e}$.
Summation gives the displayed $2^{\sum d_e}$ bound. The inequality
$M(F)\le\|F\|_1$ follows from the triangle inequality on the unit torus;
the product law follows from addition of logarithms. Finally, univariate
Jensen bounds $M(F)$ below by the measure of its leading coefficient in
one variable. Iterating to a nonzero integer coefficient proves
$M(F)\ge1$ for integer $F$.

Let $F$ now be the determinant numerator

$$F=\det[C_{re}(Uy)_e^i]_{e\ne e_0}.$$

Since $F=\pm t_{e_0}U^{n(H-1)}P$ and the other factor is integral,

$$
\|P\|_1\le2^{N(n-1)}M(P)
\le2^{N(n-1)}M(F)
\le2^{N(n-1)}\|F\|_1.
$$

The current-polynomial norms are obtained exactly from spanning trees:

$$
(Uy)_e^i=\sum_{\text{spanning trees }T}
\left(\prod_{f\notin T}t_f\right)(y_T)_e^i,
\qquad
\|(Uy)_e^i\|_1=\sum_T|(y_T)_e^i|.
$$

Different trees have distinct cotree monomials. Form the nonnegative
integer matrix whose entries are

$$|C_{re}|\,\|(Uy)_e^i\|_1,\qquad e\ne e_0.$$

The permanent of this matrix bounds $\|F\|_1$ by the determinant
expansion and the triangle inequality. The code computes the permanent
by subset dynamic programming and chooses the best of the original
gauge edges.

Every chain Dirichlet moment is positive and at most one. Grouping
original monomials into core monomials and integrating them therefore
does not increase the coefficient one-norm. Reciprocation only permutes
the exponent labels. If $D_{\rm mom}$ is the moment denominator, every
integer coefficient of $D_{\rm mom}Q$ is bounded in absolute value by

$$
\boxed{B=D_{\rm mom}\,2^{N(n-1)}
\min_{e_0}\operatorname{per}
\left[|C_{re}|\,\|(Uy)_e^i\|_1\right]_{e\ne e_0}.}
$$

For rational sources use $d\lambda$ in this expression, and reconstruct
the coefficients of $D_{\rm mom}d^{nH}Q$.

### Six-loop fixture certificate

For the complete source assignment in
[`examples.json`](examples.json), `six_loop_K5_H2`, the exact certificate
is

$$
\begin{aligned}
D_{\rm mom}&=216000,\\
\min\operatorname{per}&=3185461923232071622777592101208064,\\
B&=25384925169104535878250163645356935162028057805002375168000.
\end{aligned}
$$

The bound has **195 bits**. Reconstruction is certified when the product
of the distinct CRT primes **exceeds $2B$**. Seven primes of the sizes
used in the implementation suffice. Once that threshold is crossed,
each signed integer coefficient is unique. Independent extra primes and
exact rational fibre evaluations check the implementation; they do not
replace this height argument.

## 6. No singular modular point on the prescribed core grid

The lower-set nodes have integer conductances $1\le g_e\le n$. The
grounded core Laplacian determinant is a positive sum over spanning
trees. If the core has $\tau$ spanning trees, then

$$1\le\det L(g)\le\tau n^{v-1}.$$

For $K_5$, $\tau=125$ and $n=6$, giving

$$\boxed{1\le\det L(g)\le125\cdot6^4=162000.}$$

A prime above this bound cannot divide the determinant at any grid
point. Also $U(1/g)=\det L(g)/\prod_eg_e$, so the cycle denominator is
nonzero modulo the same prime. This proves nonsingularity over the
entire interpolation grid before sampling it.

## 7. Convergence of the individual core periods

For an admissible original graph, each nonempty proper union of chains
has total original edge count $\sum_{e\in T}L_e>Hh_T$.
For a permitted reciprocal monomial the period powers satisfy

$$
\begin{aligned}
\sum_{e\in T}a_e
&=\sum_{e\in T}L_e+(n-1)|T|-\sum_{e\in T}\beta_e\\
&\ge\sum_{e\in T}L_e+n h_T\\
&>(n+H)h_T.
\end{aligned}
$$

Also $a_e\ge L_e>0$ and the full homogeneity is $\sum a_e=n(n+H)$.
These are the strict Schwinger-face convergence inequalities for the
projective core period. Thus the support theorem supplies a
term-by-term convergence certificate when the original graph's
admissibility has been checked. No divergent monomial is assigned a
finite table value by fiat.

## 8. API, checks and limits

```python
plan = reciprocal_plan(reduction)
conductance_points = plan.sample_points(start, stop)
exponents = plan.exponents()              # original chain order
valid = plan.support_mask(exponents)
period_powers = plan.period_powers(exponents)
certificate = integer_certificate(reduction)
```

`reciprocal_plan` works for a connected loopless series core within its
explicit graph-enumeration and grid-allocation budgets. It raises on a
closed-chain core or an exhausted budget. The present planner concerns
the bare $T=0$ numerator; the finite-current extension needs its
correspondingly shifted degrees and moment denominators.

The node count and degree bounds depend on the core, not on a choice to
discard sources. Other cores may have a larger reciprocal degree than
the original degree; inversion is then not automatically the better
route. Likewise a complete numerator does not supply a missing core
period table. The six-loop $K_5$ target is useful because its complete
$q=8$ period table already exists.

Run the independent controls from the repository root:

```sh
python3 explore/circuit-general/series-fibres/test_reciprocal.py
```

The seven test groups passed in the recorded development run. They
compare complete original-numerator and reciprocal coefficient
dictionaries at every gauge, exact fibre values with the full Jacobian,
the integer and rational-source certificates, graph support counts,
Laplacian bounds, and explicit failure behavior at the resource limits.
