# Non-series-parallel holomorphic graph periods

At four loops on $\mathbb C^2$ all sixty zero-displacement Ward identities agree coefficient by coefficient for arbitrary conserved momenta. The 218 edge contractions equal the independent integrals of the 84 lower-loop faces with simple quotients. The Ward matrix has rational rank $2344$ on $2359$ coefficients. Local numerator constraints, integration by parts and two angular-mode sums reconstruct all thirteen kernels. The 23 cubic product forests of the author's paper are fixed at zero displacement by lower loops. [Exact rank, radial integrals and physical-spoke derivatives](verify/full_ward/README.md)

## The three-loop graph on $\mathbb C^2$

Subdivide edge $12$ of $K_4$ at vertex $4$. Use ordered edges $(01,02,03,14,42,13,23)$, root $4$, and arbitrary conserved sources $\lambda_0,\ldots,\lambda_4$. Set $A=\lambda_1-\lambda_2$, $B=\lambda_0-\lambda_3$, $C=\lambda_0+\lambda_3$, $D=\lambda_4$ and write $[XY]=\det(X,Y)$. In the project orientation the zero-displacement kernel is

$$
K_\Gamma(\lambda;0)=\frac{[AD]}{2304}
\left(-3[AB]^2+6[BC]^2-[BD]^2\right).
$$

The determinant density has weight three under $GL(2)$, hence is cubic in the brackets of $A,B,C,D$. The transposition of vertices $1,2$ exchanges three pairs of ordered edges; the transposition of $0,3$ exchanges two. The cycle-basis change has positive sign in two holomorphic dimensions. Thus the source polynomial is odd in $A$ and even in $B$. If $D=0$ the two serial resistances enter only through their sum, so their splitting differential drops out and the kernel vanishes. Eliminating $[AC][BD]$ with the Plücker relation leaves eleven cubic monomials with these parities and a factor containing $D$. Eleven original-graph evaluations give an integer interpolation determinant of $69695105787626017259520$. Exactly three coefficients survive. The other $\mathbb C^2$ graphs through three loops are series–parallel, giving rational coefficients at zero displacement. [Exact interpolation, basis and cross-check](verify/three_loop_sources.py)

At $\lambda=((1,0),(0,0),(0,1),(-1,0),(0,-1))$ one has $[AD]=0$. The two components of the current on edge $01$ give

$$
\begin{aligned}
K_\Gamma(\lambda;0)&=0,\\
-\partial_{z_{01}^1}K_\Gamma(\lambda;0)&=\frac{17}{9450}-\frac{\zeta(3)}{700},\\
-\partial_{z_{01}^2}K_\Gamma(\lambda;0)&=-\frac{79}{50400}+\frac{\zeta(3)}{700}.
\end{aligned}
$$

Both derivatives are irrational by Apéry's theorem. Thus the rationality of lower-loop inputs with derivative insertions, considered in [§5.7 of the previous verification](VERIFY.md), fails at one insertion. This addresses a premise of the four-loop quadratic-relations argument. The rank of those relations is a separate question.

In the four-mark signed standard representation the four-loop Ward matrix has exact rank $197$ on $200$ coefficients. The known prism mark placement leaves one free direction affecting $K_{3,3}$. Supplying its independently integrated coefficient conversely fixes the missing prism rung coefficient at $-237/6400-71\zeta(3)/720+3\zeta(5)/20$. A second Ward equation fixes the nonedge placement at $-1/28800+\zeta(3)/9600$. [Quadratic relations and all four prism coefficients in this representation](verify/four_mark_relations.md)

For an independent fixture with all sources nonzero, $\lambda=((1,2),(2,-1),(-1,3),(-3,-4),(1,0))$, the integrals with insertions $1,y_{01}^1,(y_{01}^1)^2$ are

$$
-\frac{71}{12},\qquad \frac{222317}{25200}-\frac{3021}{700}\zeta(3),\qquad
-\frac{4003781}{302400}+\frac{3772}{525}\zeta(3).
$$

Split the chain $1\to4\to2$ at time fraction $s$. Its endpoint sources become $\lambda_1+(1-s)\lambda_4$ and $\lambda_2+s\lambda_4$. For a star centered at $0$ or $3$ use its three normalized conductances $g_i$, the three parallel fractions $r_i$ and $Q=\sum_i g_i(1-r_i)$. The original Kirchhoff currents have polynomial numerators over $Q$. The bordered determinant of $Q$ and the numerators of the six cotree currents gives the oriented density over $Q^7$. Each inserted current adds its numerator and one power of $Q$. Integrating over $s$ termwise, homogenizing in the $g_i$, and integrating the beta series gives the results above. The minimum absolute corner convergence margin in both charts is one. All Kirchhoff equations are checked as polynomial identities and their orientations are checked against the seven-edge density. [Exact graph and moment replay](verify/three_loop_current.py)

## Signed current coordinates

For a $T=0$ graph with $n$ loops and $Hn+1$ edges, use the project's oriented cycle matrix $C$ and a rational source flow $y^0$ for rational sources. Write each current component as $y^i=(y^0)^i+C^{\mathsf T}\ell^i$ with $\ell^i\in\mathbb R^n$. Stack the rows $C\,\operatorname{diag}(y^i)$ over the $H$ components to form $\Psi(\ell)$. Its signed maximal minors are $\delta_e=(-1)^e\det\Psi_{\widehat e}$ in zero-based edge order. The positive projective Schwinger region away from its zero-density locus maps bijectively to the two explicitly semialgebraic current regions

$$
\{\ell:\ \delta_e(\ell)>0\ \text{for every }e\},\qquad
\{\ell:\ \delta_e(\ell)<0\ \text{for every }e\}.
$$

Indeed $\Psi t=0$ is the cycle voltage equation, and the cofactor vector lies in its kernel. At full rank the inverse is $t_e=\delta_e/\sum_f\delta_f$. Conversely, these positive times have the specified currents by uniqueness of the Kirchhoff solution. At smaller rank every top minor and the Schwinger density vanish. Differentiating the voltage equation gives $M\,d\ell^i=-C\,\operatorname{diag}(y^i)\,dt$, with $M=C\,\operatorname{diag}(t)\,C^{\mathsf T}$. Its exterior determinant cancels $(\det M)^{-H}$ in the density. The sign is that of $\delta_e$, since on the normalized section $D(t)=\sum_e\delta_e/(\det M)^H$. The bare integral is the difference of the two ordinary volumes, and a polynomial current insertion is integrated as the same polynomial on the two regions. This supplies constant-density coordinates for the E6 graph set, including non-series-parallel graphs.

For the paired unit sources in the preceding calculation the seven-edge density is $-6/52521875$ at $(1,1,1,1,1,1,2)$ and $+6/52521875$ at $(1,1,1,1,1,2,1)$. Hence the two regions have equal positive volumes and an irrational difference of first moments. A finite signed rational-polytope decomposition would give a rational integral of every rational linear current, contradicting the evaluated derivative. [All eight E6 graph inverses and signed-density replay](verify/current_chambers.py)

## Higher-loop values

For every three-loop holomorphic graph, every finite current insertion at rational sources has period in $\mathbb Q+\mathbb Q\zeta(3)$. A single eight-source assignment gives a nonzero $\zeta(3)$ coefficient on each of the six $\mathbb C^3$ graphs, including the two whose four-mark benchmark is rational. The exact evaluator accepts arbitrary rational momenta and any list of oriented edge-current components. [All six values, proof and evaluator](verify/k4_periods/README.md)

The three-loop $\mathbb C^3$ kernels at the source assignment in [graphs.py](graphs.py) and [§3.2 of the verification record](VERIFY.md) have the following coefficients.

| $K_4$ subdivision | Coefficient |
|---|---:|
| $[1,1,1,1,3,3]$ | $(39-32\zeta(3))/67200$ |
| $[1,1,1,2,2,3]$ b | $(1296\zeta(3)-1597)/5443200$ |
| $[1,1,1,2,2,3]$ c | $(34-27\zeta(3))/453600$ |
| $[1,1,2,2,2,2]$ b | $(39-32\zeta(3))/67200$ |

The exact replay verifies the Kirchhoff equations on the original graph as polynomial identities, fixes the orientation against the direct density, and integrates the current Jacobian over the chain simplices. A multinomial series and beta moments then integrate the reduced density. Both star centers give the same coefficient for each graph. [Graph replay](verify/graph_k4_series.py) · [Series and convergence checks](verify/reduced_k4_series.py) · [Eight exact results](verify/out_graph_k4_series_exact.log)

Run the following in this directory with Python 3, SymPy, NumPy, NetworkX and python-flint installed.

~~~sh
python3 verify/graph_k4_series.py 0:0 0:1 2:0 2:1 3:0 3:1 5:0 5:1
~~~

Each argument gives the zero-based graph and star-center indices. The output coefficients with keys 0 and 3 multiply $1$ and $\zeta(3)$, respectively.

At four loops on $\mathbb C^2$, the fields generated by all bare graph values at rational sources are $\mathbb Q$ on the five series–parallel graphs, $\mathbb Q(\zeta(3))$ on four graphs with $K_4$ core, and $\mathbb Q(\zeta(3),\zeta(5))$ on each of the two wheel topologies, the prism and the bipartite graph. All eight non-series–parallel graphs admit explicit irrational periods. [Thirteen exact kernels and field witnesses](verify/prism_forcer/README.md)

The complete four-loop prism kernel on $\mathbb C^2$ has 40 signed bracket coefficients for six conserved sources. Both inequivalent cuts agree on every one of its 1,481 projective monomial integrals. At a reflection-symmetric source the bare graph vanishes, while a rational linear combination of first displacement derivatives is $588751/907200-5123\zeta(3)/9450$. [Kernel, exact evaluator and derivative](verify/prism_forcer/README.md)

The nine-edge four-loop prism on $\mathbb C^2$ has been replayed from the original graph at zero displacement. With sources at vertices $0,1,2,3$, all fifteen bracket coefficients agree with the [project table](results_nonSP_C2.json). With $\lambda_2=-\lambda_0-\lambda_1$ and the remaining three vertices source-free, the coefficient of $\langle\lambda_0,\lambda_1\rangle^4$ is $1/34560$ in both star charts. The verifier checks the nine-edge current equations, computes the current Jacobian, integrates the three parallel splitting fractions, and evaluates the resulting beta series. [Exact prism calculation](verify/prism_series.py) · [Coefficient replay](verify/out_prism_series_exact.json)

~~~sh
python3 verify/prism_series.py four 0
python3 verify/prism_series.py three 0
python3 verify/prism_series.py three 3
~~~

For $K_{3,3}$ with sources at $0,1,2,3$, the graph automorphisms force the kernel to be a single scalar times $\langle0,2\rangle\langle1,3\rangle(\langle1,3\rangle^2-\langle0,2\rangle^2)$. The scalar is $(2160\zeta(5)-1414\zeta(3)-539)/14400$, strictly positive in the project orientation. Two choices of cut reduce the Schwinger period to three-loop propagator integrals and give the same finite part. [Exact coefficient and replay](verify/k33_forcer/README.md). The proof subtracts exact divergences from the symmetrized Schwinger numerator and gives rational bounds by ordered sectors. [Proof and bounds](verify/k33_positive/README.md) · [Exact replay](verify/k33_positive/check.py)

~~~sh
python3 verify/k33_positive/check.py
~~~

For all six conserved sources at zero displacement, the full $K_{3,3}$ kernel has five independent signed bracket orbit sums. Four coefficients contain $\zeta(3)$ and $\zeta(5)$, and the coefficient of the signed orbit of $[01][02][12][34]$ is $1/3456$. [Full coefficient table](verify/k33_forcer/README.md#arbitrary-conserved-sources) · [Seven direct graph evaluations](verify/k33_forcer/all_sources.py)

At the integer source assignment $((-3,-1),(-2,2),(0,-1),(1,-1),(3,1),(1,0))$, the $K_{3,3}$ period reduces to $-789/640+59\zeta(3)/80$. Apéry’s theorem gives an irrational rational-source kernel. The signed cotree subdivision identity propagates its existence to all holomorphic dimensions $H\ge2$ at four loops. [Subdivision identity](../ht-families/dimension_raise.md)
