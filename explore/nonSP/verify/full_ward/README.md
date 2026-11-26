# Four-loop Ward identities on $\mathbb C^2$

At zero displacement the thirteen bare kernels are quartics in the brackets $[ij]=\det(\lambda_i,\lambda_j)$ of conserved vertex momenta. Eliminate the last vertex momentum. In a standard bracket monomial the four increasing pairs have weakly increasing first and second entries. There are 490 such monomials on five independent momenta and 1,764 on six. All thirteen kernels have been reconstructed in this basis. The three entries at each stored monomial multiply $1,\zeta(3),\zeta(5)$. [Full polynomial table](polynomials.json) · [Graph integration and reconstruction](polynomials.py)

For the eleven graphs with a series or parallel reduction, set the first Schwinger time to one. The numerator over $U^6$, with $U$ the Kirchhoff polynomial, has total degree fifteen and degree at most three in each time. Transpose interpolation on $\{1,2,3,4\}^8$ to apply the exact integrals of the time monomials directly to the determinant samples. The source interpolation uses 490 independent integer assignments and two more for comparison with a separate monomial integration. For the prism and $K_{3,3}$ the input is the independently cut full kernel. [Dual integration weights](cubature.py) · [Integer determinant](cubature.cpp) · [Primitive cuts](../prism_forcer/README.md)

## Rank and determining integrals

For each graph, sum the actions of its automorphisms on the 490 monomials, weighted by the sign of the ordered edge permutation. This operator squares to its product with the order of the automorphism group. Its trace divided by that order gives the dimension. In catalog order the thirteen dimensions are

$$
(22,266,92,254,490,254,106,86,254,266,224,40,5).
$$

Contraction of an admissible edge merges the two source momenta. Its zero-based edge position gives the sign $(-1)^e$, followed by the sign restoring the quotient edge order. Straightening $[il][jk]=[ik][jl]-[ij][kl]$ for $i<j<k<l$ produces the rational matrix for all sixty seven-vertex sliding graphs. Its 218 edge terms act on 2,359 independent source coefficients. Modular elimination gives rank at least 2,344. Fifteen lifted vectors annihilate every integer row and restrict to an invertible diagonal minor. Thus the rational rank is exactly 2,344. [Contractions](contractions.py) · [Integer null vectors and rank replay](homogeneous.py)

The ranks of restriction from this kernel to each individual graph are

$$
(0,0,0,0,5,0,4,5,10,5,15,10,4).
$$

The five series–parallel kernels are fixed by lower-loop identities alone. The four $K_4$ reductions together detect five free directions. The prism and $K_{3,3}$ together detect ten. Wheel 10 alone detects all fifteen. Fifteen scalar evaluations of that wheel suffice to fix all thirteen kernels once the lower-loop faces are supplied. The absolute determinant on the integer kernel basis used by the replay is $98826674411378747484712419961633505280$. Applying its inverse recovers every $\zeta(5)$ coefficient of all thirteen graphs from the fifteen wheel values. This statement uses the explicit rational coefficient triples throughout. [Fifteen inputs and values](wheel_inputs.json) · [Reconstruction](wheel_inputs.py)

## Analytic evaluation of the two remaining wheel periods

Use the wheel with four spokes. Its eight projective times $t_i$ have Kirchhoff polynomial $U$. At denominator $U^q$ let the propagator powers be $\alpha_i$, so the numerator is $\prod_i t_i^{\alpha_i-1}$ and $\sum_{\substack{i\ \mathrm{edges}}}\alpha_i=4q$. The two radial integrals used here have $q=2$ and $q=3$. Give each spoke of the original wheel power $q-1$ and each rim edge power one. In the planar loop-momentum variables the four rim propagators then have power $q-1$ and each spoke has power one. Cutting a spoke of the dual wheel fixes one of its rim momenta at unit length. The normalized loop measure is $d^{2q}x/\pi^q$. Integrating the cut Schwinger time supplies the projective conversion factor $\prod_i\Gamma(\alpha_i)/\Gamma(q)$. [Cut normalization](../k33_forcer/README.md)

Write $\lambda=q-1$ and let $m$ be the angular harmonic degree. The rim propagator between radii $r,s$ expands in the Gegenbauer polynomials $C_m^\lambda$ with radial coefficient $\max(r,s)^{-2\lambda}\exp(-m|\log r-\log s|)$. With unit angular measure the degree-$m$ projector kernel is $(m+\lambda)C_m^\lambda/\lambda$. Indeed its value on the diagonal is the dimension of degree-$m$ harmonics, obtained by subtracting traced symmetric tensors, and $C_m^\lambda(1)=\binom{m+2\lambda-1}{2\lambda-1}$. The three angular integrations leave $C_m^\lambda(1)(\lambda/(m+\lambda))^3$. The expansion is absolutely convergent at distinct neighboring radii; after angular integration the spectral terms are positive.

In logarithmic radii the radial Jacobians cancel the pure radial powers from the four rim edges and the three uncut spokes. The remaining integral is the fourfold convolution of $\exp(-(m+\lambda)|x|)$ at zero. Fourier transformation gives $5/(2(m+\lambda)^3)$. Each radial measure contributes $2/\Gamma(q)$. Thus the projective period is

$$
\frac{20\lambda^3\,\Gamma(\lambda)^4}{\Gamma(q)^4}
\sum_{\substack{m\ge0\\\mathrm{angular\ harmonics}}}
\frac{\binom{m+2\lambda-1}{2\lambda-1}}{(m+\lambda)^6}
=
\begin{cases}
20\zeta(5),&q=2,\\[2pt]
\dfrac53\bigl(\zeta(3)-\zeta(5)\bigr),&q=3 .
\end{cases}
$$

For the reduction to these periods, differentiation of $t^\beta/U^5$ in $t_i$ with $\beta_i=0$ gives 751 convergent boundary identities. At infinity the boundary vanishes; at zero the wheel contracts to a series–parallel graph or a $K_4$ reduction. On a triangle with resistances $x,y,z$, substitute star resistances $r_1,r_2,r_3$, with $x=S/r_1$, $y=S/r_2$, $z=S/r_3$, $S=r_1r_2+r_2r_3+r_3r_1$. The Kirchhoff polynomials satisfy $U_\triangle=(x+y+z)U_{\rm star}$ and the logarithmic Jacobian has determinant one. The integrand gains $S^{\alpha_x+\alpha_y+\alpha_z-12}$ and each triangle power becomes $6-\alpha_i$. A nonnegative exponent of $S$ gives a finite sum of lower-loop integrals. Triangle reduction and planar inversion give 548 such scalar inputs.

In the 269-dimensional space allowed by the homogeneous wheel integration-by-parts equations these 1,299 conditions have rank 217. Multiplying the two radial integrands by the required equal powers of $U$ in numerator and denominator gives rank 219. Both remaining wheel source evaluations lie in this row span; the replay verifies their exact rational recombination on all 269 columns. The full Ward reconstruction consequently derives all 3,892 nonzero coefficient triples with no four-loop period measurements as input. [Scalar reductions](wheel_faces.py) · [Radial evaluations and rational row spans](radial_wheel.py) · [Full reconstruction](bootstrap.py)

## Lower-loop face integrals

The [lower-face evaluator](lower_faces.py) accepts seven conserved rational two-component sources in the ordered sliding graph. For an inner cut of loop number $n$ and outer loop number $m=4-n$, it reconstructs the exact inner bare kernel on the standard bracket monomials of degree $n$. Their dimensions for $n=1,2,3$ are $1,6,50$. Two further evaluations test each reconstruction. At an inner vertex $v$ substitute its source plus incoming exterior currents minus outgoing exterior currents. The outer Kirchhoff polynomial $U$ clears each current denominator, so the product of the outer numerator with the homogenized inner polynomial has denominator $U^{m+2+2n}$. Series and parallel beta integrals, or the two-loop propagator cut for an outer $K_4$ subdivision, integrate its monomials. Moving the face to the right-hand side multiplies by $-(-1)^{mn}$ and the sign of the split edge ordering.

At sources $((1,2),(-2,1),(3,-1),(-1,-2),(2,3),(-4,1),(1,-4))$, the [60 exact graph fixtures](lower_faces.json) retain the edge orders, every contributing cut and both rational coefficients in the basis $1,\zeta(3)$. Of the 225 proper tight subsets with three to five vertices, 141 collapse to parallel edges. The other 84 comprise 41 triangle, 21 two-loop and 22 three-loop inner cuts. The 60 signed edge-contraction sums agree with the independently integrated right-hand sides. Forty-five values are nonzero, including seven with a nonzero $\zeta(3)$ coefficient. For ordered edges $01,05,06,12,14,23,25,34,45,56$, the only cut is $056$ and the right-hand side is $-55211/25200+3243\zeta(3)/700$.

The all-source check uses the tree expansion of each outer current. In a face with inner loop number $n$, the inner kernel has degree $2n$ in the spanning-tree weights after clearing its current denominators. For a monomial with tree multiplicities $\alpha$, its coefficient is the order-$\alpha$ forward difference at zero divided by the product of the multiplicity factorials. For every such monomial the replay integrates the product with the outer Kirchhoff numerator. Interpolation in each component of the outer source uses the integer simplex of total degree $m$. Its homogeneous cardinal polynomial at a composition $\alpha$ of $m$ is

$$
\prod_j\frac{1}{\alpha_j!}
\prod_{r=0}^{\alpha_j-1}
\left(x_j-\frac r m\sum_{\substack{k\ \mathrm{independent\ outer\ vertices}}}x_k\right),
$$

where $x_j$ is that component at independent outer vertex $j$. These two finite interpolations evaluate the lower-loop terms at arbitrary sources. Independently, the upper side merges the momenta across each contracted edge in the thirteen integrated graph polynomials. Both sides agree in every coefficient channel at 1,764 integer assignments with full bracket-evaluation rank modulo 65,521. They are therefore equal as polynomials. The edge-contraction sums have a nonzero rational channel on 45 graphs, a nonzero $\zeta(3)$ channel on seven, and zero $\zeta(5)$ channel on all sixty. [Tree moment integral](tree_moments.py) · [All 60 polynomial comparisons](verify_identities.py)

Restriction of the homogeneous solution space to at most three, four or five nonzero vertex momenta has rank $1,11,15$, respectively, jointly over the entire catalog. Wheel 10 attains each of these ranks by itself. In particular four directions vanish on every four-mark experiment. [All individual restriction ranks](wheel_inputs.json)

## The cubic operator product in four dimensions

Keep the ordered marked edge $(A,B)$ of the author's product. Deleting its endpoints leaves a forest on four cubic interaction vertices. Each tree has one vertex of full degree two, which carries its surviving momentum. Every degree-three interaction uses all its fields and carries zero momentum. Enumerating placements up to graph automorphisms fixing the ordered pair gives 23 forests. Restricting the exact kernels to the two insertions and the forest roots gives twenty nonzero polynomials. All twenty are rational, including the contributions from graphs whose unrestricted kernels involve $\zeta(3)$ and $\zeta(5)$. The table uses the catalog edge order, before the interaction-specific Wick and color factors. [Product-forest theorem in the author's paper](../../../../sources/LLC_main.tex) · [All 23 source polynomials](product_forests.json) · [Exact restriction](product_forests.py)

The restriction of the fifteen homogeneous Ward solutions to all 23 forests has rank one. On wheel 10 the placement $(A,B)=(2,3)$ has root 1 and gives $-[\lambda_2,\lambda_1]^4/34560$, where $[\lambda_2,\lambda_1]=\det(\lambda_2,\lambda_1)$. The placement $(2,4)$ vanishes. On graph 9 the placement $(0,1)$ has root 4, gives $[\lambda_0,\lambda_4]^4/34560$ and already detects the single homogeneous direction in the product.

## Current insertions in the cubic product

For a marked product forest, the least total derivative order of a nonzero geometric moment on the physical insertion spokes has the following distribution. An odd exchange of two identical nonempty branches kills each word untouched by that exchange. Remaining zeros are checked by exact polynomial integration: fix the first two independent momenta as a basis and retain symbolic components for every other root. This suffices on a dense open set. A single nonzero rational evaluation certifies survival. The oriented current densities and all edge currents agree with a direct Kirchhoff solve and the full resistance-to-ratio Jacobian through five loops. [Enumeration and first-degree replay](first_degree.py) · [Pointwise density checks](forest_chart_check.py)

| Loops | Degree zero | First at degree one | First at degree two | First at degree three |
|---:|---:|---:|---:|---:|
| 1 | 1 | 0 | 0 | 0 |
| 2 | 3 | 0 | 0 | 0 |
| 3 | 7 | 1 | 0 | 0 |
| 4 | 20 | 3 | 0 | 0 |
| 5 | 51 | 14 | 0 | 0 |
| 6 | 149 | 44 | 1 | 0 |
| 7 | 414 | 158 | 4 | 0 |
| 8 | 1216 | 515 | 22 | 0 |
| 9 | 3514 | 1761 | 93 | 1 |

The three identically vanishing bare forest placements acquire nonzero first displacement derivatives. For graph 6 at $(A,B)=(2,3)$ and sources $((0,0),(0,1),(1,0),(-3,-2),(0,0),(2,1))$, inserting the second component of the $02$ current gives $-1/259200$. For wheel 10 at $(A,B)=(2,4)$ and sources $((0,0),(0,1),(1,0),(0,0),(-1,-1),(0,0))$, inserting the first component of the $23$ current gives $-1/414720$. Reversing the marked placement preserves this nonzero physical source assignment. Both numbers agree with independent propagator reductions of the original graph. The ratio evaluator accepts any list of current insertions and integrates each interval or simplex with its dependent factors. [Polynomial current evaluator](forest_jets.py) · [Independent cuts](forest_jet_cuts.py)

For every positive integer $k$, join $k$ three-loop forks to the same marked edge $AB$. In each fork a degree-two root $u$ joins two source-free vertices $p,q$, each incident on $A$ and $B$. Order each block as $(pu,pA,pB,qu,qA,qB)$ and place $AB$ last. Set $\lambda_A=(1,0)$, each $\lambda_u=(0,1)$, and $\lambda_B=(-1,-k)$. Interchanging the leaves of one untouched fork permutes three pairs of edges. Hence every displacement coefficient of total degree less than $k$ vanishes. Inserting the first component of the $pu$ current once in every fork attains degree $k$.

Indeed write $l=y_{uA}^{\rm red}$, $h=y_{uB}^{\rm red}$ and $y=y_{AB}^{\rm red}$ in the paired star formula of the paper. The inserted current is $-r_1 l-r_2 h$. Integrating it against $r_1-r_2$ over the unit square gives $(h-l)/12$. In the $j$th fork the untouched part of the marked current has fraction $1-r_3-r_4$. Each earlier fork contributes two factors of that fraction through $[l,y][h,y]$. Its triangle integral is $1/((2j)(2j+1))$. The product of these denominators is $(2k+1)!$.

In the remaining parallel network let $s_j$ be the fraction of root source sent to $A$, let $r_j$ be the current fraction on the $j$th series branch, and let $r_0=1-\sum_{\substack{1\le j\le k\\\mathrm{branches}}}r_j$. The total reduced current is $L=(1,\sum_{\substack{1\le j\le k\\\mathrm{branches}}}s_j)$. With $P=(0,1)$ the two branch currents are $l_j=s_jP-r_jL$ and $h_j=(1-s_j)P+r_jL$, while $y=r_0L$. In the paper's edge ordering the reduced density is $r_0[P,L]^k$: successive series reductions give the marked-edge fractions, equivalently, the chord-current Jacobian contains the identity matrix with the column of branch fractions subtracted from each column, whose determinant is $r_0$. Thus $(h_j-l_j)^1=2r_j$, $[l_j,y][h_j,y]=r_0^2s_j(1-s_j)[P,L]^2$ and $[P,L]=-1$. The ordered moment in the paper's $\mathcal I$ convention is

$$
\frac{(-1)^k}{6^k(2k+1)!}
\left(\int_0^1 s(1-s)\,ds\right)^k
\int_{\substack{r_j>0\\ r_0>0}} r_0^{2k+1}\prod_{j=1}^k r_j\,d^k r
=\frac{(-1)^k}{36^k(4k+1)!}.
$$

On the spoke $pA$ ending at the physical insertion $A$, the current is $r_1 l-r_3 y$. Its integral over the unit square against $r_1-r_2$ gives $l/12$ and $l^1=-r_j$. Replacing all $pu$ insertions by $pA$ thus gives $1/(72^k(4k+1)!)$. The first nonzero derivative order is therefore unbounded even for jets acting on the physical insertion. The original ratio evaluations at three, six, nine and twelve loops reproduce the formula. [Exact replay](fork_jets.py) · [Paired star and product-forest formulas](../../../../sources/LLC_main.tex)

## A nonzero physical-spoke moment for every cubic product forest

At loop number $n$, set $\lambda_A=(1,0)$, set every tree root momentum to $(0,1)$, and fix $\lambda_B$ by conservation. In each reduced network the first component is a unit current from $A$ to $B$. The voltage difference $\Phi_A-\Phi_B$ therefore has strictly positive first component. At an internal vertex $v$ take the scalar voltage $[\Phi_v-\Phi_B,\Phi_A-\Phi_B]$. Both boundary values vanish. In each remaining tree its only source is the nonzero determinant of the root source with that voltage difference. The inverse Dirichlet Laplacian has strictly positive entries within the connected tree: the discrete minimum principle first gives nonnegativity for a positive source, and any interior zero would propagate to the source. The scalar voltage is thus nonzero at every internal vertex. A prior series reduction transfers only second-component root sources to the boundaries.

It follows that both determinants of an internal-to-boundary current with the marked current are nonzero. Every single star factor and every series factor is nonzero in the positive ratio domain. Each paired star has only the additional interior zero $r_1=r_2$. Thus the original density $D(t)$ is nonzero on an open set for every forest in this source assignment.

Let $e$ be the marked $AB$ edge and $\delta_e$ the signed maximal minor of the current matrix $\Psi$ with its $e$th column removed. Writing $U$ for the Kirchhoff polynomial gives

$$
\delta_e=t_e U^2D(t),\qquad
\int_{\sum t=1}\delta_e D(t)\,dt
=\int_{\sum t=1}t_eU^2D(t)^2\,dt>0.
$$

Here $dt$ is the positive simplex measure used in the graph kernel. Boundedness of the currents and absolute integrability of $D$ ensure convergence. Each term of this minor contains one current component from each of the $2n$ physical edges. Cutting an internal tree edge on the side away from its root leaves only source-free vertices. Kirchhoff conservation expresses that edge current as a signed sum of the boundary spokes of the cut subtree. Consequently some moment of degree $2n$ on spokes incident on $A$ or $B$ is nonzero. This is a statement about the geometric kernel before color and Wick contractions. [Projective density and signed current variables](../../README.md#signed-current-coordinates) · [Original-graph cofactor and internal-edge elimination through five loops](forest_cofactor.py)

## Local numerator restrictions

The joint order-two zero derived below supplies five independent equations. The ten wheel measurements previously sufficient are reduced to two by the wheel integration-by-parts relations below. All 23 physical product coefficients already follow from lower-loop data and the local numerator zero. [Reconstruction](bootstrap.py)

At a bivalent vertex with zero source, the two incident resistances enter the currents through their sum. The current Jacobian therefore vanishes. More generally, for a two-terminal block with terminals $u,v$ and internal source sum $P$, interchange the terminals on the block. Negate its internal sources and add $P$ to each terminal source. The potentials inside can be replaced by the sum of the two terminal potentials minus the old potential. All currents inside change sign in the inherited edge directions. The two-component current form has even cycle-basis parity, leaving only the sign from the catalog edge ordering. This gives the Whitney flip relation also between catalog graphs 3 and 5.

For source $p$ at a bivalent vertex, the ratio of its two resistances to their sum parametrizes transfer of $p$ between its two neighbors. Splitting this unit interval at its midpoint gives the sum of two kernels with bivalent source $p/2$ and the other half at one or the other neighbor. The oriented current coordinates on each half differ from the corresponding kernel's coordinates by constants in the resistances. The 23 bivalent zero-source tests, every two-terminal block flip and all 23 subdivision identities hold separately on the physical source coefficients and on all 15 homogeneous Ward directions. [Exact restrictions](soft_sources.py)

Graph 9 has a stronger simultaneous source zero. Its subnetwork from 0 to 2 consists of the parallel paths $03$ and $043$, followed by $32$. At $\lambda_3=\lambda_4=0$, fix its effective resistance and the division of current between $03$ and $043$. Every edge current is independent of the subdivision ratio on $043$ and of the series ratio between the parallel subnetwork and $32$. These are independent positive resistance coordinates. When the two internal sources scale by $s$, both Jacobian columns are divisible by $s$; the density has a zero of order at least two. Its integrated kernel inherits this polynomial source identity by coefficientwise integration. The 63 standard bracket monomials of combined $\lambda_3,\lambda_4$ degree less than two have zero coefficients. They give five independent conditions on the fifteen homogeneous directions. On the surviving ten directions only wheels 8 and 10, the prism and $K_{3,3}$ have nonzero restriction ranks, respectively $10,10,10,4$. Checking all 23 ordered product forests gives rank zero. Wheel 10 reaches ranks $0,6,10$ with at most three, four and five marked vertices. [Exact quotient and product-forest restrictions](joint_soft.py)

For each wheel graph fix the eight core resistances and integrate the bivalent fraction alone. Its source polynomial satisfies the midpoint identity at these fixed core values. Average over graph automorphisms with their edge-order signs. The dimensions permitted by signed symmetry and midpoint subdivision are 70 and 67 for wheels 8 and 10. For the prism and $K_{3,3}$ signed symmetry permits dimensions 40 and 5. All four upper bounds are attained by the corresponding graph-local density polynomials. To see this, the Kirchhoff numerator has individual time degrees at most three. At fixed series sum its two substituted times make a polynomial of degree at most six in the bivalent fraction. The cofactor of the tree-current numerators, divided by its omitted time, differs from this numerator by $U^4$, constant along the series fraction. The seven interior nodes $1/8,\ldots,7/8$ therefore integrate it exactly with the first seven ordinary moments. For the two unsubdivided graphs use the cofactor directly. On 88 integer source probes, finite-field minors at positive rational resistance choices reach ranks 70, 67, 40 and 5. Consequently the remaining ten directions satisfy every graph-local linear source identity that holds at each choice of the residual resistances, even after the wheel fraction integral. Relations between periods after integration over those remaining resistances impose additional information. [Numerator degree](../k4_periods/README.md) · [Exact span certificates](source_span.py)

The two wheel graphs also admit a common source-space test. Suppressing their bivalent vertices gives the same four-spoke wheel $W_4$. For either subdivision let $R$ be the resistance of the collapsed edge. The coordinate change from the two original times to $R$ and their fraction contributes the factor $R$. The resulting core numerator has degree sixteen, with common projective denominator $U^6$. Average the two numerators together over the eight core automorphisms. Planar self-duality gives an additional common average: if $e^*$ is the dual edge, then $t'_e=1/t_{e^*}$ and $U(t)=(\prod_e t_e)U(t')$. The projective inversion Jacobian and $U^{-6}$ send a core numerator $P(t)$ to $(\prod_e t_e)^4P(t')$.

At 148 shared resistance assignments in $\{1,2,3,6\}^8$, the replay integrates the bivalent fractions, includes both averages and computes integer cofactors. A 136-row minor is nonzero modulo 1009. Solving this minor over the rationals expresses all ten surviving pairs of wheel source polynomials as linear combinations of the same 136 fiber rows. Evaluation at all 88 source probes in each graph verifies these coefficients exactly; the preceding separate rank certificates make these probes separating on the two allowed polynomial spaces. Thus a single common moment functional respecting wheel automorphisms and planar inversion is compatible with every residual homogeneous direction. [Common wheel certificate](joint_wheels.py)


## Wheel integration by parts

With the eight wheel resistances $t_i$ and $U$ their Kirchhoff polynomial, fix one resistance other than $t_i$ to one. For every degree-twelve monomial $t^h$ used in the calculation,

$$
\frac{\partial}{\partial t_i}\frac{t_i t^h}{U^5}
=\frac{t^h}{U^6}\big((h_i+1)U-5t_i\partial_{t_i}U\big).
$$

The primitive vanishes at zero and infinity since $0\le h_i\le3$ and $U$ is affine with positive leading coefficient in $t_i$. For each proper nonempty edge subset $S$, the replay checks $\sum_{\substack{e\in S\\\mathrm{edges}}}(h_e+1)>5\min_T|S\setminus T|$, where $T$ ranges over spanning trees. The expanded degree-sixteen terms satisfy the corresponding strict inequalities for $U^6$. These inequalities ensure absolute integrability in every ordered resistance sector. After wheel symmetry and planar inversion the 4,328 divergences on 1,292 moment orbits have exact rational rank 1,023. The 269 independent integer null vectors and a modular minor certify both bounds.

Integrating the subdivided edge fraction contributes $h_a!h_b!/(h_a+h_b+1)!$ for its two powers. The resulting moment response on 88 source assignments has rank four on either wheel and jointly. Appending the ten Ward directions gives joint rank twelve, leaving two free directions. On each wheel, the prism and $K_{3,3}$ these two directions restrict to ranks zero, one and two at three, four and five marked vertices. [Exact moment and intersection replay](wheel_ibp.py)

Eight integer relations on thirteen assignments of wheel 10 realize the reduction. The square reconstruction uses 2,344 Ward evaluations, five numerator zeros, the eight relations and wheel measurements 10 and 11 in the zero-based input table. It recovers all 3,892 stored coefficient triples and reproduces the 23 physical product restrictions with both wheel inputs set to zero. [Eight relations](wheel_ibp.json) · [Two input values](wheel_inputs.json) · [Exact square solve](bootstrap.py)

## Arbitrary integer powers on an $n$-spoke wheel

Let $n\ge3$. For positive integer Schwinger powers $\alpha_e$ with sum $nq$, integrate $\prod_e t_e^{\alpha_e-1}/U^q$ on the projective wheel simplex, with integer $q$. Assume the strict edge-subset convergence inequalities above. The complement of one edge gives $1\le\alpha_e<q$. For each rim propagator put $\beta=q-\alpha_e$, $\lambda=q-1$, and $d=\lambda-\beta$. The same rule $\beta=q-\alpha_e$ gives the position powers of the spokes. Gaussian integration of the $n$ rim vertices in $2q$ dimensions gives the inverse $q$th power of the spanning-tree polynomial. Inverting all Schwinger variables converts it to the cotree polynomial and sends $\beta_e$ to $q-\beta_e$. The normalized position measure is $d^{2q}x/\pi^q$ per vertex, multiplied by $\prod_e\Gamma(\beta_e)$. Fix one squared radius to one to remove the common logarithmic scale.

For spherical harmonic degree $m\ge0$, put $s=m+\lambda$. Its multiplicity is

$$
\binom{m+2\lambda+1}{2\lambda+1}-\binom{m+2\lambda-1}{2\lambda+1}
=\frac{2s^2}{(2\lambda)!}\prod_{j=1}^{\lambda-1}(s^2-j^2).
$$

After removing the symmetric homogeneous powers of the two radii, the rim eigenvalue in their ratio $r<1$ is $\sum_{k=0}^d E_{d,k}(s)r^{s+2k-d}$, where

$$
E_{d,k}(s)=\frac{\lambda!}{(\beta-1)!}
\frac{(-1)^k\binom dk}{\prod_{j=0}^{d}(s+k-d+j)}.
$$

Indeed expansion of $(1-2rz+r^2)^{-\beta}$ in Gegenbauer polynomials gives the degree-$m$ eigenvalue $r^m(\beta)_m/(\lambda+1)_m$ times the terminating polynomial $\,{}_2F_1(m+\beta,-d;m+\lambda+1;r^2)$, with $(a)_m=a(a+1)\cdots(a+m-1)$. The symmetric radial factor supplies $r^\beta$. The unsplit eigenvalues are nonnegative. The positive Laplace representation of $(1-2rz+r^2)^{-\beta}$ is a positive mixture of $\exp(uz)$ with $u>0$; its power series consists of positive-semidefinite tensor-power kernels. Therefore the cyclic harmonic traces may be summed under the radial integral.

Write $v_i=\log|x_i|$ for the rim radii in cyclic order, $b_i$ for the spoke position power, and $\beta_i$ for the rim power on $i,i+1$. Select one $k_i$ on each rim. The radial integrand at $v_1=0$ is

$$
\exp\left(\sum_{\substack{i\\\mathrm{rim\ vertices}}}(2q-2b_i-\beta_{i-1}-\beta_i)v_i
       -\sum_{\substack{i\\\mathrm{rim\ edges}}}(s+2k_i-d_i)|v_i-v_{i+1}|\right).
$$

In each ordering of the $v_i$, the $n-1$ gaps integrate independently. For the initial set $S$ below one gap the inverse factor is $\sum_{i\in S}(2q-2b_i-\beta_{i-1}-\beta_i)+\sum_{i\in\partial S}(s+2k_i-d_i)$, with $\partial S$ the crossing rim edges. Positivity at $m=k_i=0$ follows by restricting the convergent integral to separated radial annuli and generic open angular sets. Increasing $m$ or $k_i$ preserves positivity.

The radial integral has only integer poles in $s$. Write the drift coefficient at $i$ as $c_i-c_{i-1}$, which is possible since their sum vanishes, and put $p_i=s+2k_i-d_i$. Fourier transformation on the differences $v_i-v_{i+1}$ expresses their convolution at zero as a contour integral of

$$
\prod_i\frac{2p_i}{(z-c_i+p_i)(c_i+p_i-z)}.
$$

The common contour strip follows by applying gap positivity to every contiguous segment and its complement. Poles within either half-plane translate together as $s$ varies. A left pole can meet a right pole only when $2s=c_i-c_j-(2k_i-d_i)-(2k_j-d_j)$. The recurrence for $c_i$ gives $c_i-c_j\equiv\beta_i-\beta_j\pmod2$, so every such $s$ is integral. The ordered-gap formula bounds each pole order by $n-1$.

Under $s\mapsto-s$, replace each $k_i$ by $d_i-k_i$. The $p_i$ all negate. Simultaneous negation of $p_i,c_i$ multiplies each ordered-gap term by $(-1)^{n-1}$. Reversing all $v_i$ makes the integrated expression even in the common sign of the $c_i$. Since $E_{d,d-k}(-s)=-E_{d,k}(s)$ and the angular multiplicity is even, the full summand is odd rational. Each eigenvalue adds at most one pole at any integer, all at $|s|\le\lambda-1$. The multiplicity vanishes there and has a double zero at zero. Thus every odd pole order is at most $2n-3$. Summing partial fractions for integers $s\ge\lambda$ gives finite rational harmonic corrections. Coefficients at opposite poles cancel all even zeta values. The harmonic-series coefficient vanishes by convergence. Every convergent wheel monomial therefore belongs to

$$
\mathbb Q+\sum_{\substack{2\le j\le n-1\\\mathrm{weights}\ 2j-1}}
\mathbb Q\,\zeta(2j-1).
$$

For spoke powers $q-1$ and rim powers $1$, all $d_i$ and all drifts vanish. The radial factor is the convolution of $n$ copies of $e^{-s|v|}$ at zero, with Fourier integral $\int_{\mathbb R}(2s/(s^2+y^2))^n\,dy/(2\pi)$. Including the squared-radius normalization gives

$$
\binom{2n-2}{n-1}
\sum_{s=1}^{\infty}
\frac{2s^2\prod_{j=1}^{q-2}(s^2-j^2)}{(2q-2)!\,s^{2n-1}},
\qquad 2\le q<n .
$$

The additional terms $s<q-1$ vanish. Its positive tail diverges at $q\ge n$; $q\ge2$ is required by positive spoke powers. Increasing $q$ from two to $n-1$ successively introduces the weights $2n-3,2n-5,\ldots,3$ with nonzero coefficient. These scalar periods span exactly the rational span of those odd zeta values. At five loops they are $70\zeta(7)$, $35(\zeta(5)-\zeta(7))/6$, and $7(\zeta(3)-5\zeta(5)+4\zeta(7))/36$.

For a purely holomorphic kernel with wheel core, integrate each subdivided path at fixed total resistance. Its monomial powers combine as $\alpha=\sum_{e\in\mathrm{path}}(1+h_e)$ and its simplex factor is $\prod_{e\in\mathrm{path}}h_e!/(\alpha-1)!$. Finite current insertions increase $q$ and the numerator degree while preserving strict convergence. Thus every displacement coefficient at rational sources of every such graph lies in the displayed odd-zeta space, at every loop order and in every admissible holomorphic dimension. This applies in particular to both once-subdivided wheel topologies on $\mathbb C^2$ for $n\ge4$. [Numerator and convergence](../k4_periods/README.md#schwinger-numerator-and-convergence) · [Exact algorithm and independent cut comparisons](wheel_angular.py)

### Exact reduction to one univariate rational function

Put $w=z-c_i$ temporarily on rim edge $i$. Sum its complete angular eigenvalue before taking radial residues. With $p=s+2k-d$ the link identity is

$$
\sum_{k=0}^d\frac{(-1)^k\binom dk}{\prod_{j=0}^d(s+k-d+j)}
                 \frac{2p}{p^2-w^2}
=\frac{2\,4^d d!}{\prod_{k=0}^d((s+2k-d)^2-w^2)}.
$$

The residues at $w=-p$ agree since $\prod_{j\ne k}((s+2j-d)^2-p^2)=4^d(-1)^kk!(d-k)!\prod_{j=0}^d(s+k-d+j)/p$. Both sides are even in $w$ and decay at infinity. This removes every independent angular summation index.

For each rim vertex $i$, form the integers $c_i-d_i,c_i-d_i+2,\ldots,c_i+d_i$. Across all vertices they have one parity. Subtract their common minimum and divide by two. Let $N_j$ count the occurrences of the resulting integer $j$ and put $P(w)=\prod_j(w-j)^{-N_j}$. The contour strip is nonempty precisely when the largest reduced integer is less than $q-1$. In addition, $\sum_{i\in\mathrm{rim}}\alpha_i>q$ is required for the coincidence of all rim positions away from the hub. Any proper cluster away from the hub has acyclic internal edges and is integrable since each position power is below $q$. The preceding radial-sector inequalities cover hub and infinity clusters.

Expand $P(w)=\sum_{j,k} A_{j,k}/(w-j)^k$ by Taylor expansion at each distinct integer. Let $R=\sum_jN_j$. Changing the contour variable above to half its displacement from the first left pole cancels every power of two from the radial measure and the link identities. In addition to the even spherical multiplicity, the harmonic summand becomes

$$
\bigg(\prod_i (b_i-1)!d_i!\bigg)
\sum_{j,k,l,t}
\frac{(-1)^{R+t}\binom{k+t-2}{k-1}A_{j,k}A_{l,t}}
     {(s+l-j)^{k+t-1}}.
$$

Here $j,l$ range over the occupied integers and $k,t$ run from one through $N_j,N_l$. Indeed the two pole factors after rescaling are $(w-j)^{-k}$ and $(s+l-w)^{-t}$, and the left residue is the coefficient of $(w-j)^{k-1}$ in the second factor. Multiplying by the spherical multiplicity, removing the cancelling polynomial part and summing the remaining integer partial fractions gives the rational and all odd-zeta coefficients exactly. Reflecting the occupied integers around their midpoint leaves the result invariant. Thus monomials with the same multiplicities up to reflection require only one evaluation after the displayed factorial has been removed.

Both the ordered-sector algorithm and this residue algorithm agree with all 1,421 independently cut four-loop wheel periods in every coefficient. The residue algorithm also reproduces the radial formula through ten loops. For the five-loop monomial with cyclic spoke powers $(5,5,4,5,5)$, rim powers $(3,2,2,2,2)$ and $q=7$, both algorithms give
$$
-\frac{45001}{43545600}-\frac{116849}{29937600}\zeta(3)
+\frac{50837}{2993760}\zeta(5)-\frac{7}{594}\zeta(7).
$$
[Exact residue and weighted-numerator evaluator](wheel_residues.py) · [Independent sector integral](wheel_angular.py) · [Highest-weight finite-field check](wheel_highest.cpp)

The same univariate link follows directly from the position propagator. For $\alpha=q-\beta$ let $\Delta$ be the flat Laplacian in $2q$ dimensions. On positive radius, $-\Delta |x|^{-2\beta}=4\beta(\alpha-1)|x|^{-2(\beta+1)}$. The flux of $|x|^{2-2q}$ through a sphere gives $-\Delta |x|^{2-2q}=4\pi^q\delta_0/(q-2)!$. Therefore
$$
(-\Delta)^\alpha |x|^{-2\beta}
=\frac{4^\alpha(\alpha-1)!\pi^q}{(\beta-1)!}\delta_0.
$$
For a spherical harmonic $Y_m$, one has $-\Delta(r^uY_m)=[s^2-(u+\lambda)^2]r^{u-2}Y_m$. Convolve $|x-y|^{-2\beta}$ with $|y|^{z+\beta-2q}Y_m(\widehat y)$ in measure $d^{2q}y/\pi^q$. Factoring $(|x||y|)^{-\beta}$ out of the propagator gives $2|x|^{z-\beta}Y_m(\widehat x)/\Gamma(q)$ times the bilateral Laplace transform of its radial eigenvalue. Applying the preceding distributional identity fixes that transform to
$$
\frac{2\,4^d d!\,\Gamma(q)}{(\beta-1)!}
\prod_{j=0}^d\frac{1}{s^2-(z+d-2j)^2}
=
\frac{2\,4^d d!\,\Gamma(q)}{(\beta-1)!}
\prod_{j=0}^d\frac{1}{(s+2j-d)^2-z^2}.
$$
The integral converges for $|\Re z|<\beta$ and determines the rational transform on its full harmonic strip.

The scalar wheel span is attained even at fixed $q=n+2$. Choose Schwinger spoke powers $(1,n+1,\ldots,n+1)$ and rim powers all two. The occupied multiplicities are $(1,2,\ldots,2,1)$ on $0,\ldots,n$. Before the spherical multiplicity the pole orders are at most three. The double zero at the origin and single zeros at the other possible poles leave no odd zeta coefficient of order at least three. This convergent positive monomial therefore gives a nonzero rational number. For each radial master with $2\le q'<n$, multiply numerator and denominator by $U^{n+2-q'}$. The resulting monomials are individually convergent since every cotree has at least the cycle rank on each subset. Their integrals are at $q=n+2$ and sum to the master. The preceding triangular radial formula then gives the reverse inclusion. The scalar period span is exactly the rational span of $1,\zeta(3),\ldots,\zeta(2n-3)$, without an arithmetic independence assumption. Further multiplication by $U$ gives the same equality at every integer $q\ge n+2$.

## An exact five-loop holomorphic kernel

Use the five-spoke wheel with hub 0 and cyclic rim $1,\ldots,5$. Subdivide $01$ at vertex 6. Use lexicographically sorted edges and ordered vertex momenta $((-1,2),(2,-1),(-2,-1),(1,2),(1,-2),(-2,1),(1,-1))$. The zero-displacement kernel in the project orientation is

$$
\frac{38942526833}{279936000}+\frac{158670023}{453600}\zeta(3)
-\frac{2789753}{9450}\zeta(5)-\frac{1257}{5}\zeta(7).
$$

Fix the first resistance to one. The numerator has degree at most four in each remaining resistance. Its values on $\{1,2,3,4,5\}^{10}$ have magnitude at most $32\cdot5^{10}(187\cdot5^5)^5$. Indeed the cofactor has 32 terms, each tree-current component has magnitude at most five, and the Kirchhoff polynomial is bounded by $187\cdot5^5$. The inverse one-variable interpolation matrix has maximum absolute row sum 56. Every beta denominator divides 2520. Hence each orbit coefficient multiplied by 2520 is an integer of magnitude at most $2520\cdot5^{10}\cdot56^{10}\cdot32\cdot5^{10}(187\cdot5^5)^5$. The product of eleven distinct replay primes exceeds twice this bound. Two further resistance points per prime use independent grounded solves. The exact 28,162 nonzero orbit coefficients give the displayed four channels by the angular integral. The coefficient $-1257/5$ is nonzero in this reduction. Rational independence of the four displayed periods is open. [Exact numerator and lift](five_wheel.py) · [Angular integration](wheel_residues.py)
