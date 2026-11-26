# Four-loop quadratic relations with four marked sources

Use the source's sliding graphs with seven vertices and ten ordered edges. A proper Laman subset $S$ contributes the inner graph on $S$ acting on the quotient, with $\partial_v$ the sum of outgoing displacement derivatives minus the incoming ones. The sign is the signature of the reordering that puts the inner edges first. With the project kernel $K$ the additional relative sign for a split into $m$ and $n$ loops is $(-1)^{mn}$. It follows from $I^{\rm S1}_n=(-1)^{n(n+1)/2}K_n$. [Budzik et al., (3.1), (3.4), (3.11)–(3.13)](https://arxiv.org/abs/2207.14321) · [Project signs](../../../code/htk/CONVENTIONS.md)

For four vertices carrying the conserved source vectors, write $[ij]=\det(\lambda_i,\lambda_j)$. Use the three unordered partitions into pairs, ordering the endpoints of each pair and the pairs themselves lexicographically. The corresponding polynomials are

$$
[ij][kl]\big([kl]^2-[ij]^2\big).
$$

Each changes sign under reversal of one pair and under interchange of the pairs. The three polynomials are independent. Their characters on $1,(12),(12)(34),(123),(1234)$ are $(3,-1,-1,0,1)$, the signed standard representation of $S_4$. The full four-source polynomials at this order have character $(15,3,3,0,1)$: they are quartic polynomials in the three independent brackets of three momenta, which transform in the signed standard representation. The character inner product is one. Thus the coefficient of a fixed displayed polynomial is extracted by the signed mean over its eight pair-preserving permutations, divided by that polynomial at the test momenta. The replay uses $(1,0),(0,1),(2,0),(-3,-1)$, for which the $(12)|(34)$ test value in this ordered four-vector list is $-6$.

If an edge contraction merges two marks, the three-source quartic is a multiple of the fourth power of the bracket of the two unmerged sources. These possibilities belong to the permutation representation on unordered pairs, which has no signed standard component. Thus the projected equations close on four distinct marks. Their terms containing an unknown four-loop period are precisely the contractions of a single edge. A six-vertex Laman subset inside a seven-vertex sliding graph would leave a leaf, giving the tautological edge–graph identity. Leaf-free sliding graphs suffice.

Enumerate the thirteen four-loop graph types and sixty leaf-free sliding types. There are $200$ source coefficients compatible with graph automorphisms and the sign from the ordered edges. For each sliding type and each paired choice of four marks, contract every admissible edge $e$ in zero-based edge order with sign $(-1)^e$. Include the signs from reordering the quotient edges and the two source pairs. Discard zero rows and identify primitive integer rows up to sign. The resulting matrix has $1050$ rows and rational rank $197$. Its three-dimensional homogeneous solution space restricts to the four prism coefficients and the $K_{3,3}$ coefficient as follows, for the ordered project graphs in [graphs.py](../graphs.py).

| Graph and pairing | Row on a basis of the homogeneous kernel |
|---|---:|
| prism $(01)\vert(23)$ | $(-4,0,0)$ |
| prism $(02)\vert(13)$ | $(2,1,-1)$ |
| prism $(01)\vert(24)$ | $(-4,0,0)$ |
| prism $(05)\vert(12)$ | $(0,-4,0)$ |
| $K_{3,3}$ $(02)\vert(13)$ | $(0,0,-4)$ |

The first two prism rows are already fixed by the graph integral with marks $0,1,2,3$. They leave one homogeneous direction, proportional to $(0,1,1)$ in the table. It changes the fourth prism coefficient and the $K_{3,3}$ coefficient equally. This remains true with all lower-loop right-hand sides specified. [Exact enumeration and rational nullspace replay](four_mark_relations.py)

The first two prism coefficients are $-13/28800+\zeta(3)/3200$ and $-13/57600+\zeta(3)/6400$ in table order. Denote the four prism coefficients in the table by $p_1,\ldots,p_4$ and the last coefficient by $k$. One projected Ward row for $p_3-p_1$ and a rational combination of sixteen rows for $p_4-k+4p_2+2p_1$ involve only a triangle acting on the three-loop seven-edge $K_4$ subdivision. For an original external edge current oriented from $i$ to $j$, the induced triangle source at $v$ is its prescribed source plus $(\delta_{vj}-\delta_{vi})y_{ij}$, summed over the external edges. The ordered triangle contributes $-[\lambda_i+\partial_i,\lambda_j+\partial_j]/2$ after this substitution. Its three-loop graph numerator therefore has common denominator $Q^9$ in either star chart. The signed eight-point projections comprise 16 evaluations for the first Ward combination and 144 for the second. Integrating their summed numerators gives

$$
\begin{aligned}
p_3-p_1 &= \frac1{2400}-\frac{\zeta(3)}{4800},\\
p_4-k+4p_2+2p_1 &= -\frac9{6400}+\frac{\zeta(3)}{1200}.
\end{aligned}
$$

Both star centers give the same two results with absolute corner convergence margin one. The independently integrated value $k=(2160\zeta(5)-1414\zeta(3)-539)/14400$ then fixes both remaining prism coefficients:

$$
p_3=-\frac1{28800}+\frac{\zeta(3)}{9600},\qquad
p_4=-\frac{237}{6400}-\frac{71}{720}\zeta(3)+\frac3{20}\zeta(5).
$$

These are all prism coefficients in this four-mark representation. The $\zeta(5)$ term in $p_4$ comes from $k$; the two lower-loop remainders have been evaluated explicitly. No rational-independence conjecture for odd zeta values enters. [Ward combinations, original three-loop integrands and both exact replays](prism_rung.py) · [Independent $K_{3,3}$ integral](k33_forcer/README.md)
