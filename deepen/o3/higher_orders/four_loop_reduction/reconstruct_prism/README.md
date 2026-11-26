# Prism core, graphs 11–38

On the cyclic probe the three coefficients, including the graph automorphism factors and preceding the common four-loop normalization, are

| Rank power | Integrated coefficient |
|---|---|
| $N^0$ | $\pi^2(9796575\pi^2-96791482)/12247200$ |
| $N^2$ | $\pi^2(88921725\pi^2-884922356)/734832000$ |
| $N^4$ | $\pi^2(144949856-14698125\pi^2)/653184000$ |

For a monomial with exponents $d_i$ in the nine core edge coordinates let $a_i=d_i+1$. Then $\sum_i a_i=22$. Cutting any edge and setting the external momentum squared to one gives the three-loop massless propagator integral with the remaining eight indices. With loop measure $d^D k/\pi^{D/2}$ the factor converting this integral to the projective monomial period is $\prod_i\Gamma(a_i)/\Gamma(11/2)$. Indeed $U=R_{\rm cut}U'+F'$, where $U',F'$ are the two Schwinger polynomials of the cut graph. Integrating $R_{\rm cut}$ gives a beta factor $B(a_{\rm cut},11/2-a_{\rm cut})$; the remaining projective integral has exactly the propagator exponents. All stored monomials satisfy the strict subset convergence inequalities; the finitely many cut integrals are analytic in a neighborhood of $D=11$.

Put $\epsilon=(4-D)/2$ and $t=\epsilon+7/2$. The three cut prescriptions use (i) edge 0 with lexicographically maximal index tuple in its automorphism orbit, (ii) edge 1 with maximal tuple and (iii) edge 0 with minimal tuple. The cut index comes first in the tuple. For each rank power $F_i(t)$ denotes the resulting sum of normalized cut integrals multiplied by the monomial coefficients and $\prod_e\Gamma(a_e)$. Their three values at $t=0$ agree. FORCER reduces each $F_i$ to one common non-gamma master and four bubble products. If $A_i$ is its rational coefficient, then

$$
\frac{-A_jF_i+A_iF_j}{A_i-A_j}
$$

has the same limiting value. For every pair $i,j$ the quotient $A_i/(A_i-A_j)$ is regular at zero. The non-gamma master cancels identically. With $z=D/2$ the remaining FORCER constants in the raw measure are

$$
G_j=(j+1)\frac{\Gamma(z-1-j\epsilon)\Gamma(z-1)\Gamma(2+j\epsilon-z)}
{\Gamma(1+j\epsilon)\Gamma(2z-2-j\epsilon)} .
$$

The products are $G_0^3,G_0^2G_1,G_0^2G_2,G_0G_1G_2$. Expanding through the finite coefficient and dividing by $\Gamma(11/2)$ gives the table. Every negative Laurent coefficient cancels. Each of the three pairs of prescriptions gives all three entries independently. [Exact replay](evaluate.py) · [Stored rational coefficients](masters) · [Projective numerators](../reconstruct_cyclic/core11-38.json)

Run `python3 evaluate.py` with SymPy and python-flint. To regenerate from the numerator, run `python3 regenerate.py --form /path/to/form --forcer /path/to/forcer`, with NetworkX available. The recorded reductions use FORM 5.0.2 and FORCER commit `68a30bfb2aec3672a2f674ca8666bfc970105e24`. [FORCER](https://arxiv.org/abs/1704.06650)

**Individual cyclic graph values.** For each individual graph and each rank power, take its rational master vector and subtract the stored rank-four prism cut difference multiplied by the ratio of their non-gamma master coefficients. This ratio is regular at $t=0$ for every entry. The two cuts of each graph use different pairs of auxiliary cuts. They give identical results before summing the twenty-eight graphs. Their sums reproduce each entry of the table above. Coefficient pairs in the JSON multiply $(\pi^2,\pi^4)$. [Values by graph](graph_values.json) · [Two-cut replay and scalar-period function](graph_values.py) · [Stored tagged reductions](graphs)

**Arbitrary sources.** Two cuts agree on all 2,989 exponent orbits and yield 2,608 distinct exact period values. Integrating the 28 graphs gives 268 nonzero coefficients among the 420 quartic terms. Every graph has a nonzero kernel. For a prism numerator in the recorded edge order, `scalar_period(record)` reduces its monomials by graph automorphisms and takes a rational scalar product with the stored periods. The optional arguments `form,forcer` evaluate convergent monomials outside the table. The input includes `h` when supplied. The full replay checks each auxiliary multiplier for regularity, compares both cuts and tests the bank against all 84 cyclic graph and color coefficients. [Period bank](monomial_periods.json) · [Integrated kernels](../kernels/core11.json) · [Monomial calculation](monomials.py)

To regenerate the arbitrary-source table, use the source numerator `core11.json` produced by the parent directory’s `source_kernels.py`. In this directory run `python3 monomials.py --source /path/to/source-numerators/core11.json /path/to/scratch --form /path/to/form --forcer /path/to/forcer --save ../kernels/core11.json`. After FORM has run, omitting `--form` and `--forcer` reuses the tagged reductions in scratch. The replay checks all individual projective integrals for strict positivity with rational Machin bounds and checks the non-gamma master through second order against the independent nonplanar cut.
