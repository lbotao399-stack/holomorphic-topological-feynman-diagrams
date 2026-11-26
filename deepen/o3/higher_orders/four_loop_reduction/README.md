# Four-loop cubic graphs on nine Schwinger parameters

For any of the [50 cubic graphs](../four-loop-graphs.json), let $R_e$ be the sum of the times along core edge $e$, and $h_e$ the number of external bivalent vertices on that chain. There are nine core edges and $\sum_{e\in\mathrm{core}}h_e=4$. With $U=\det M$ and $M=C\operatorname{diag}(R)C^{\mathsf T}$ for the oriented four-cycle matrix $C$, integrating the four splitting fractions gives the projective density
$$
\frac{\prod_{e\in\mathrm{core}}R_e^{h_e}}{U^{11/2}}P(R),\qquad
\deg P=9,\qquad \deg_{R_e}P\le2 .
$$
Here $P$ depends on the four conserved momenta, the orientations and the edge order. At integer momenta, $4\cdot8!P$ has integer coefficients. The same formula holds after any linear combination of source assignments. A fourth-order action coefficient also uses the color contraction, graph symmetry factors and the common heat-kernel normalization.

**Integration along a chain.** Use the orientation of the chain to order its segment currents. If $Q_j$ is the cumulative sum of its external momenta before segment $j$, its current is $y_{\rm first}+Q_j$. Put $x_j=t_j/R_e$ and eliminate the first $x_j$. The time-weighted mean of the offsets is $m_e=\sum_{j\in\mathrm{segments}}x_jQ_j$, with $Q_{\rm first}=0$. For a fixed particular first-segment flow $y^0$ obeying vertex conservation, the mean currents on all chains equal
$$
\bar y=\bigl(1-C^{\mathsf T}M^{-1}C\operatorname{diag}(R)\bigr)(y^0+m).
$$
This follows by subtracting the unique cycle flow that sets all voltage drops around cycles to zero. The source one-forms involve exactly
$$
\sum_{j\in\mathrm{segments}} (y_{\rm first}+Q_j)dt_j
=\bar y_e\,dR_e+R_e\,dm_e .
$$
The topological Pfaffian depends on $R$ only. Its coefficient on four ordered horizontal edges $i,j,k,l$ before multiplication by $U^{-1/2}$ is
$$
\frac{\det C_{ijkl}}4
\bigl(B_{ij}B_{kl}+B_{ik}B_{jl}+B_{il}B_{jk}\bigr),
\qquad B=C^{\mathsf T}M^{-1}C .
$$
Expanding the Pfaffian of $-dM M^{-1}dM/2$ gives this sign since exchanging two selected edges changes both the exterior order and the cycle-column determinant. Of the eight holomorphic one-forms, four supply $R_e\,dm_e$ and four supply mean currents times $dR$. The former give $\prod R_e^{h_e}$. The remaining polynomial has degree at most four in the splitting fractions. Each simplex monomial of exponents $a_j$ contributes $\prod_j a_j!/(h_e+\sum_j a_j)!$ in the independent fractions on its chain. The chart in the script keeps the first eight $R$'s, followed by the split fractions, and compares its oriented differential with the positive projective density. [Original thirteen-edge comparison and exact simplex integration](chain_moments.py)

**Pole order.** Fix the other eight times at generic values and vary a single $R_e$. The spanning-tree polynomial $U$ is affine with a simple root there. At the root $M$ has rank three. Over the local ring of the variable at this root, a regular invertible congruence reduces $M$ to a regular invertible three-by-three block and a scalar with a simple zero. In this moving basis transform the one-form matrix $dM$ by congruence as a tensor; differentiation of the basis is unnecessary. Denote its four entries in the last column by $w_i$. Its last row consists of these same one-forms.

In each of the two blocks of source forms, the singular part of row $i$ is a scalar times $w_i/U$. Indeed the singular part of $M^{-1}$ has just its last diagonal entry in this basis, and the singular term of $\bar y$ is a multiple of $C^{\mathsf T}$ applied to this last direction. A regular row operation on the two holomorphic blocks leaves possible poles in only one block. If both scalar residues vanish there are no holomorphic poles.

A singular inverse chosen in one of the two Pfaffian factors uses two of the $w_i$. With $l=1,2$ singular Pfaffian factors, at most $4-2l$ additional singular holomorphic forms can survive exterior multiplication. The resulting pole order is at most $l+(4-2l)\le3$. With $l=0$, every Pfaffian term still contains a $w_i$ since one of its free row or column labels must be the fourth. Thus at most three singular holomorphic forms survive in this case too. The common factors outside these forms are $U^{-2}$ from the two holomorphic directions and $U^{-1/2}$ from the topological direction. Multiplication by $U^{11/2}$ therefore clears the pole. Choose a different projective pivot from the varied edge to avoid any chart denominator at zero. Since this works in each independent variable, the reduced numerator is polynomial; equivalently its primitive denominator divides a power of $U$, whose nonconstant factors are detected by these one-variable tests.

For the degree bound take the varied edge as projective pivot. At positive infinity all currents remain bounded by the finite tree-current values, and $C^{\mathsf T}M^{-1}C$ is bounded by monotonicity of $M$. After factoring out the four vertical $R_e$'s the remaining density is bounded by a constant times $U^{-5/2}/R_e$. As $U$ is affine this gives $\deg_{R_e}P\le2$. The projective density has scale degree $-9$, $U$ has degree four and the vertical factors have degree four, giving $\deg P=22-9-4=9$.

The unreduced numerator in a projective chart has integer coefficients after multiplication by $4\cdot8!$: the Pfaffian contributes $1/4$, and all products of the simplex factorial denominators divide $8!$. The chart pivot and the spanning-tree polynomial are primitive integral polynomials, so their exact cancellations preserve integrality by Gauss's lemma. Homogeneity and the separate degree bound leave 3139 candidate monomials, the central coefficient of $(1+t+t^2)^9$. The values at $R_8=1$ on the $3^8$ grid with coordinates $1,2,3$ uniquely determine $P$.

The 47 non-series–parallel graphs on the cyclic probe reduce to eight convergent rational series at rank powers $N^0,N^2,N^4$. The replay separately checks all graph-orbit sums, evaluates the first four rational coefficients of each series and computes its remainder bound. [Core coefficients and exact integration](../../../../explore/nonSP/verify/hepp_series/README.md#five-dimensional-four-loop-cyclic-probe)

## Integration of the cyclic probe

For a numerator monomial $\prod_e R_e^{a_e-1}$ the nine positive integers $a_e$ sum to $22$. Cut one edge and take Euclidean momentum $Q^2=1$ at its ends. The cut has three loops. In dimension $D=11$, its momentum integral with measure $d^Dk/\pi^{D/2}$ per loop, multiplied by $\prod_e\Gamma(a_e)/\Gamma(11/2)$, equals the original projective period. Indeed, writing $U=R_e U_0+F_0$ gives the beta factor $\Gamma(a_e)\Gamma(11/2-a_e)/\Gamma(11/2)$, and the remaining Schwinger integral has prefactor $\Gamma(11/2-a_e)/\prod_{j\ne e}\Gamma(a_j)$. The strict subset inequalities give absolute convergence at $D=11$ and in a neighborhood for every chosen cut. [Convergence and a separate rational series](../../../../explore/nonSP/verify/hepp_series/README.md#five-dimensional-four-loop-cyclic-probe)

**Doubled edge.** Cutting edge 1 in the order of graphs 3–10 leaves the parallel pair as a bubble, edge 5 carrying $Q$, and a two-loop propagator integral. In order $0,2,3,6,4$ the five propagator momenta are $k,Q-k,l,Q-l,k-l$. Write $F(a,b,c,d,e)$ for their integral with the indicated powers. Applying $\partial_k\cdot(k-l)$ at generic dimension gives

$$
\begin{aligned}
(D-2e-a-b)F(a,b,c,d,e)
={}&a\big[F(a+1,b,c,d,e-1)-F(a+1,b,c-1,d,e)\big]\\
&+b\big[F(a,b+1,c,d,e-1)-F(a,b+1,c,d-1,e)\big].
\end{aligned}
$$

The sum $c+d+e$ decreases by one. When a line disappears, one-loop beta integrals finish the reduction. Reversing the roles of $k,l$ gives the same two rational-function coefficients before taking $D=11$. [Doubled-edge replay](doubled_kite.py)

**Prism and bipartite core.** Permuting a monomial by a symmetry of $U$ and cutting edge zero gives equal integrals at $D=11$. The continuations away from eleven dimensions can differ. FORCER expresses each cut in a fixed basis of master integrals. For the prism there is one master other than gamma products, for the bipartite core there are two. Take two or three cut expansions respectively, and choose rational-function weights summing to one that annihilate the corresponding master columns. In the bipartite case the weights are the three signed two-by-two minors divided by their sum. In both cases the weights have removable singularities at $D=11$, which preserves the target value. The negative Laurent coefficients of the remaining bubble products cancel. Different choices of the auxiliary cuts reproduce each finite part. [Prism replay](reconstruct_prism/README.md) · [Bipartite replay](bipartite_cut.py) · [Gamma expansion](odd_cut.py)

The table gives the bare, alternating cyclic probe divided by $\pi^2$. The graph automorphism orders are already included.

| Graphs | $N^0$ | $N^2$ | $N^4$ |
|---|---:|---:|---:|
| 0–2 | $0$ | $0$ | $1/10500$ |
| 3–10 | $101/945-1241\pi^2/45360$ | $73903/1020600-37\pi^2/3402$ | $34259/20412000-11\pi^2/829440$ |
| 11–38 | $(9796575\pi^2-96791482)/12247200$ | $(88921725\pi^2-884922356)/734832000$ | $(144949856-14698125\pi^2)/653184000$ |
| 39–49 | $(93253402-9461505\pi^2)/12247200$ | $(734213669-74453325\pi^2)/183708000$ | $0$ |

**Heat kernel.** Taking nine vertex integrals and thirteen edges gives the scalar ratio $(8\pi^{5/2})^{-4}$. For the Gaussian contraction choose a time facet. Assign four of its twelve time differentials to each of $\tau,\bar z^1,\bar z^2$. For each assignment the three complementary nine-edge incidence minors give the position-space wedge. The four real moments pair in three ways using $t^{-1}B(B^{\mathsf T}t^{-1}B)^{-1}B^{\mathsf T}t^{-1}/2$; the eight holomorphic moments evaluate at the source-shifted mean. Moving the time differentials to the front fixes the permutation sign. This expansion on all 50 graphs gave a nonzero signed density at each test point and agreed with the cycle-space Pfaffian. The Wick integration exponent $\binom m2+\binom{m-1}2+5(m-1)$ is even for $m=2,6,10$. The internal edge counts $1,7,13$ are odd, so the common heat-edge and external-block conventions give $\epsilon U=1$ in all three cases, as in the manuscript. [Position-space replay](raw_four.py) · [Order-four assembly](assemble.py) · [Wick dictionary](../../../../notes/flagship/flagship.tex)

Run the full coefficient with the project's Python dependencies, FORM and FORCER. FORM uses a temporary directory and runs serially.

```sh
python3 assemble.py --form /path/to/form --forcer /path/to/forcer
python3 raw_four.py --all
```

**Individual bipartite graph periods.** The two entries in each graph and rank cell multiply $\pi^2$ and $\pi^4$. Both choices of cut give the same value for each of the 19 nonempty cells. The 165 arbitrary-source bracket coefficients span 42 independent Schwinger polynomials, each evaluated at two cuts. Here the graph symmetry factor is already present in the cyclic input and is applied only at the later color contraction for the arbitrary-source kernels. [Cyclic coefficients](bipartite_cyclic.json) · [Individual cut replay](bipartite_graphs.py) · [Arbitrary-source kernels](kernels/core39.json) · [Source replay](bipartite_sources.py)

For $t=\epsilon+7/2$, any two convergent cuts of the same period have a difference divisible by $t$. Two such differences eliminate the nonbubble masters. The elimination weights for the nonplanar master are regular, and those for the triangle master are divisible by $t$. Dividing the two identities by $G_0^2G_1$, which has a simple pole, determines the normalized nonplanar master through degree one and the normalized triangle master through degree two. The residue of $G_0^2G_1/\Gamma(11/2)$ is $\pi^2/6994593273507840$. Thus a finite projective period equals this residue times the linear coefficient of the reduced cut divided by $G_0^2G_1$. [Five normalized master expansions and independent calibration cuts](bipartite_masters.py)

## All external momenta

The four files in [kernels](kernels) contain the full geometric kernels of the 50 individually labeled bare graphs at zero displacement. Each coefficient is a pair of rational multipliers of $\pi^2,\pi^4$, before the automorphism factor or color tensor. The four inputs are the bivalent vertices in order $6,7,8,9$. For their first three momenta let the three brackets be $[01],[02],[12]$. The three powers in each column refer to these brackets in this order. Of 750 columns, 409 are nonzero. The identically zero geometric graphs are 1, 2, 9 and 10. Graph 46 has a nonzero geometric quartic and zero commutator color tensor. [Direct lookup and evaluation](four_point.py)

Each reduced density pairs the four rows of the first holomorphic direction with the four corresponding rows of the second. Fix the first two momenta to $(1,0),(0,1)$ by a change of basis and put the third equal to $(-z,y)$. In each pair of rows the $z$ and $y$ variations multiply the same one-form, because both directions use the same linear Kirchhoff operator. Its exterior square vanishes. The total degree in $(y,z)$ is therefore at most four. Fifteen simplex grid points determine the polynomial. Restoring the fourth power of the first determinant gives the homogeneous quartic in the three brackets, and polynomiality extends it to dependent first momenta. For each graph an additional conserved momentum assignment with $[01]=2$ checks the interpolated numerator before outer integration. [Reconstruction and strict convergence](source_kernels.py)

After division by each graph's automorphism order, alternation of its four inputs, and contraction of its ten commutator tensors, add the trees and both kinds of marked theta graphs. The resulting polynomial in the defining traces, the rank and six independent momentum components equals $d i_4$ coefficient by coefficient. The rational part and the $\pi^2$ part in units of $(\epsilon\pi\nu)^2$ each have 300 nonzero trace and momentum coordinates, with their full rank polynomials. The remainder is zero in both parts. This independently checks the primitive for arbitrary traceless colors and arbitrary conserved momenta. [Full cochain replay](verify_four_point.py)
