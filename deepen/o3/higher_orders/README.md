# Higher orders of five-dimensional commutative Chern–Simons theory

Use the finite-jet momentum cochains of the [five-dimensional manuscript](../../../notes/flagship/flagship.tex), with zero total momentum, the defining trace on $\mathfrak{sl}_N$, and its $U(2)$-equivariant holomorphic vertices. The order-$\hbar^w$ local obstruction in this sector has Chevalley–Eilenberg degree four and momentum bidegree $(w,w)$. Write $C^0_{\mathfrak g}$ for the manuscript's momentum-conserving complex. In the formal trace algebra, equivalently for $N\ge5$ in arities through five, the complete degree-four cohomology at total momentum degree $2w$, $w\ge2$, is
$$
H^4(C^0_{\mathfrak{sl}_N})_{\deg=2w}
=\big\{[p_0P][p_2P]h(P)\ \text{in the }(01)(23)\text{ trace channel},\quad P=p_0+p_1,\quad h\text{ homogeneous of degree }2w-4\big\}.
$$
Here $[pq]$ denotes the manuscript's holomorphic momentum determinant. Permuting the four inputs with their exterior signs supplies the other pairings. Charge conjugation $X\mapsto-X^{\mathsf T}$ fixes this entire space. As an $SL_2$ representation it is the space of binary forms of degree $2w-4$. In particular the middle weight at $w=4$ belongs to spin two, and the invariant part vanishes for every $w\ge3$.

## One trace

The coefficient of the word $\operatorname{tr}(X_0\cdots X_j)$ specifies a single-trace cochain of arity $j+1$. Alternation makes it cyclic with sign $(-1)^j$. Represent the incoming momenta by sides $p_i=x_{i+1}-x_i$ of a closed polygon, with $x_{j+1}=x_0$. This identifies coefficients with polynomials in the vertices invariant under simultaneous translation. In the next cyclic word only adjacent input letters can arise from a commutator. Merging these letters adds their momenta, deleting the intervening vertex. Thus the differential, up to an overall degree sign, is the alternating sum of vertex deletions restricted to cyclic cochains.

Before imposing cyclicity, take the same deletion differential on all simultaneous-translation-invariant polynomials. In the $j$ independent successive increments its formula is the ordinary additive bar differential, with internal faces adding consecutive increments and end faces dropping the first or last. In fixed polynomial degree $m$, normalization removes the degenerate cochains. The normalized complex is the graded dual of the reduced bar complex of the two-generator polynomial algebra: dualizing the coproduct $f(p+q)$ gives multiplication in the graded dual symmetric algebra (over characteristic zero the factorial basis gives this identification). The Koszul resolution therefore puts its degree-$j$ cohomology in internal degree $m=j$, with dimension $\binom2j$. In particular all ordinary bar groups in internal degree $m>2$ vanish.

The cyclic exact sequence applies at each $m$. Indeed rotation, insertion of a zero side, and every face preserve momentum degree and the closed-polygon condition. The finite cyclic norm identifies coinvariants with invariants in characteristic zero. In the complex omitting the last deletion, a contracting homotopy appends a duplicate final vertex with sign $(-1)^{j-1}$ on degree $j>0$. Thus the usual cyclic bicomplex yields
$$
HC^{j-2}_{m}\longrightarrow HC^j_m\longrightarrow HH^j_m .
$$
Here $HC$ denotes the signed cyclic subcomplex and $HH$ the full deletion complex just computed. Starting in negative cochain degree with zero, this proves $HC^j_m=0$ for every $m>2$. The cyclic complex includes all arities during this argument. Singletons absent for traceless matrices have no effect on the arity-three or arity-four cohomology because their adjacent differentials only use arities two through five.

## Two traces

A four-letter invariant with two traces on traceless matrices has word type $(01)(23)$. Write $r=p_0$, $q=p_2$; the four momenta are $r,P-r,q,-P-q$. Let $F(r,q,P)$ be the coefficient of this word. It is odd under exchange within either pair and even under exchange of the pairs. In degree five, fix the two-letter trace and the oriented three-letter trace created by $d$. On the first three momenta $r,s,t$ with $r+s+t=P$, its coefficient, up to a common sign, is
$$
F(r,q,P)+F(s,q,P)+F(t,q,P)
=F(r,q,P)+F(s,q,P)-F(r+s,q,P).
$$
There are no further independent closure conditions except the same identity on the other pair. A polynomial additive in $r$ is linear in $r$. Exchanging $r$ and $P-r$ forces this linear functional to kill $P$. If its coefficients are $A,B$, then $P_1A+P_2B=0$ in the polynomial ring. Coprimality of $P_1,P_2$ gives $(A,B)=(P_2,-P_1)h_1$ with polynomial $h_1$. Applying the second additivity and exchange to $h_1$ gives
$$
F(r,q,P)=[rP][qP]h(P),\qquad h(-P)=h(P).
$$
Every such coefficient is closed. There is no three-letter two-trace word on traceless matrices, and the classical differential preserves the number of traces. This proves the stated cohomology. Stable trace independence through arity five follows already at $N=5$ by testing a word against off-diagonal matrix units arranged along its permutation cycles on distinct indices.

## Local quantum master equation

Once $J_1=aF_D$ satisfies the manuscript's $a^2=\varepsilon U N^2/(256\pi^4)$, every subsequent equivariant local obstruction is exact for $w\ge3$. The induction uses closure of the order-$w$ obstruction after the lower equations have been solved, the heat-kernel $U(2)$ character, and the local full-wedge holomorphic contact theorem for differentiated vertices. They place the obstruction in the invariant part of the group above. In degree three the single-trace calculation gives $H^3_{\deg=2w}=0$ for $w\ge2$. The invariant three-cochain space at bidegree $(w,w)$ is one-dimensional: conservation leaves two independent momenta, with sole determinant $[p_0p_1]$. Multiplying $[p_0p_1]^w$ by the commutator trace for even $w$, or by the symmetrized cubic trace for odd $w$, gives the manuscript's alternating cubic. Hence there is a unique invariant cubic primitive at each of these orders in the formal algebra. This determines existence of a formal local QME lift in the manuscript's scheme at every order. Numerical values of higher counterterm coefficients still require graph integrals. [Local contact and scale hypotheses](../../h1/verify/regulated_contacts.md)

Specialization is performed on the entire local identity. The trace propagator for $\mathfrak{sl}_N$ contracts by the matrix completeness relation including $-1/N$. Each graph at fixed order has coefficients in Laurent polynomials in $N$ with scalar graph periods. Choose either solution with $a/N$ constant. In the formal trace basis the differential has integer coefficients independent of $N$. A rational linear inverse on its image introduces only constant denominators. Induction thus produces a primitive over the same Laurent coefficient ring. At each finite $N\ge3$ evaluation on traceless matrices commutes with contractions and $d$, giving the local QME lift at that rank. Intrinsic cohomology in the presence of small-rank trace relations is a separate question. The decoupled extension to $\mathfrak{gl}_N$ with a free center follows by projection of every vertex to its traceless part.

For $N\ge4$ the primitive's coefficient can be extracted from a single ordered evaluation. Use the standard matrix units $E_{12},E_{23},E_{34},E_{41}$ and the successive momenta $(1,0),(0,1),(1,1),(-2,-2)$. With the manuscript's differential, the differential of the unscaled cubic on the cycle of four matrix units is
$$
-[p_2p_3]^w+[p_1p_2]^w+[p_3p_0]^w-[p_0p_1]^w
=2^w+(-1)^w-1 .
$$
Thus dividing the right-hand side of the order-$w$ cochain equation on this input by $2^w+(-1)^w-1$ recovers the exact primitive coefficient.

## Finite-rank check

The replay alternates products of momentum determinants and trace words. At arity four their dependence on the three determinants $[01],[02],[12]$ is homogeneous of degree $w$. At $[01]=1$ the grid $([02],[12])=(y,z)$, $y,z\ge0$, $y+z\le w$, determines each coefficient by triangular Newton interpolation. Below the stable rank the replay uses the exact Gram form obtained by applying the traceless projector in every matrix slot. It proves lower rank bounds for $d:C^4\to C^5$ by exact minors on integer momenta. The upper bounds follow from $d^2=0$ and, at $w=2$, the double-trace cocycle on two orthogonal commuting diagonal matrices in $\mathfrak{sl}_{N\ge3}$. For that test use $p_0=(1,0),p_1=(0,1),p_2=(1,0),p_3=(-2,-1)$ and colors $X_0=X_1=\operatorname{diag}(1,-1,0)$, $X_2=X_3=\operatorname{diag}(1,1,-2)$, padded with zeros. [Exact ranks](invariant_ranks.py) · [Output, columns $w,N,\dim C^3,\dim C^4,\operatorname{rank}d_3,\operatorname{rank}d_4,\dim H^4$](ranks.txt)

[All-order matrix lift, including mixed and rank-one sectors](center_classes.md#all-higher-orders-with-a-mixed-cubic) · [Independent cubic coefficient check](probe.py)

## Computing the order-four coefficient

For the cyclic matrix-unit probe above, the tree terms have values $\{D,D_3\}=10$ and $\{[p_0p_1]^2f,[p_0p_1]^2f\}/2=-4$. On the theta graphs the marked-vertex coefficients, divided by $\pi\nu$ and with the vertex coupling removed, are

| theta | one $[p_0p_1]^2 f$ | two distinct projected $D$ |
|---|---:|---:|
| $(1,3,3)$ | $-499N^2/2835-256/105$ | $431N^2/11340+718/2835$ |
| $(2,2,3)$ | $-34/27$ | $5/63$ |

Each $D$ is projected on all three color slots before contraction by replacing its tensor with $D-2S_1/N+4S_0/N^2$. With $i_2,i_3$ from the manuscript, trees and decorated two-loop graphs together contribute
$$
(\varepsilon\pi\nu)^2\,
\frac{21221N^4+85304N^2+192480}{3628800}
$$
to the order-four primitive. The bare four-loop contact on the same input is divided by $16$ in this cochain equation. The replay checks the derivative brackets, the two kinds of marked graph and the unprojected combination against the full order-two Moyal graph on this probe. [Calculation](order_four.py)

A bare tight trivalent four-loop graph has ten vertices and thirteen edges. Its four external fields occur at four bivalent vertices: a vertex of internal degree one violates tightness on the remaining nine vertices. Tightness also excludes a cut vertex, since summing the bounds on two subgraphs meeting at one vertex gives a strict bound for the whole graph. Suppressing the bivalent vertices leaves a biconnected six-vertex cubic multigraph. Its underlying simple graph has degrees two and three, and the edges to be doubled form a perfect matching on the degree-two vertices. A loop would disconnect the rest upon deleting its incident vertex, and triple parallel edges would form a disconnected component. Distributing four subdivisions over these nine edges, reducing by the core automorphisms and testing every induced subset gives 50 graphs. Their four surviving core types contribute $3,8,28,11$ isomorphism classes, respectively. [Enumeration](four_loop_graphs.py) · [Oriented chains and automorphism orders](four-loop-graphs.json)

The complete sum of the 50 bare graphs, before the common normalization and divided by $\pi^2$, is

$$
-\frac{172}{945}
+N^2\left(\frac{139229}{48600}-\frac{35701\pi^2}{120960}\right)
+N^4\left(\frac{760981}{3402000}-\frac{130727\pi^2}{5806080}\right).
$$

For the manuscript's Wick dictionary $\epsilon U=1$ and $\nu=-U/(64\pi^5)$. The four-loop factor in the cochain equation is $\epsilon U/(8\pi^{5/2})^4=(\epsilon\pi\nu)^2/\pi^2$. Adding the trees and marked theta graphs and dividing the bare contact by $16$ gives the unique equivariant Taylor cochain in the heat-kernel scheme:

$$
\frac{i_4}{(\epsilon\pi\nu)^2}=
\left[
\frac1{24}
+N^2\left(\frac{1102559}{5443200}-\frac{35701\pi^2}{1935360}\right)
+N^4\left(\frac{4216}{212625}-\frac{130727\pi^2}{92897280}\right)
\right][p_0p_1]^4 f .
$$

Here $f$ is the defining trace of one input against the commutator of the other two, and $[p_0p_1]$ is the determinant of their conserved holomorphic momenta. The $N^0\pi^2$ coefficient in the bracket cancels. The $N^2$ and $N^4$ coefficients in the bracket are strictly positive by $\pi^2<10$, and their $\pi^2$ terms are nonzero with the same sign. Thus at every integer $N\ge3$ the normalized coefficient is positive and transcendental. Both choices of the order-one coupling $a$ give the same $i_4$. The order-four obstruction is a boundary, and its cancellation requires this nonzero cubic in the specified scheme. [Graph sums, signs and replays](four_loop_reduction/README.md#integration-of-the-cyclic-probe) · [Series–parallel integrations](four_loop_sp/README.md)
