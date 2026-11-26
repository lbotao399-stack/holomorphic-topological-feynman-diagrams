# Twistor BF local operators

Use the three cubic vertices $b[a,a]/2$, $g\{h,h\}/2$ and $b\{h,a\}$ on the flat twistor patch. A propagator joins $b$ to $a$ or $g$ to $h$. The first two momentum coordinates in the replay are the Poisson fiber coordinates. The following coefficients use the repository's reduced holomorphic kernel, its ordered edge differences, and the rooted graph weight $1/|\mathrm{Aut}\,\Gamma|$. The external fields are polarized by the unnormalized alternating sum. [Action, §§2–3](https://arxiv.org/abs/2208.12701) · [Reduced kernel conventions](../../../code/htk/CONVENTIONS.md)

## Graphs in every degree

Write $n$ for the loop number and $k$ for the number of contracted letters of the local insertion. A tight graph on $\mathbb C^3$ has $3n+1$ edges and $2n+1$ interaction vertices. The number of surviving interaction fields is
$$
3(2n+1)-2(3n+1)+k=k+1 .
$$
For each connected component remaining after deletion of the insertion, reattach the insertion with its incident edges. Each resulting graph obeys $2|E'|\leq3(|V'|-1)-1$. Summing this bound over the components and using equality for the full graph gives at most one component. The two-vertex bound permits just one edge between a pair of vertices.

Orient each propagator from its multiplier. Every cubic vertex has exactly one multiplier. If $m$ surviving fields are multipliers, the insertion supplies $(3n+1)-(2n+1-m)=n+m$ multipliers. The connected graph after deletion has cycle rank $n-k+1$. Consequently:

| Contracted letters | Interior | Orientations |
|---|---|---|
| $k=n$ | One cycle with attached trees | All root letters are multipliers. Every tree points toward the cycle, with two choices of cycle direction. |
| $k=n+1$ | Tree | Either one root letter is ordinary, or one surviving field is a multiplier. All interior edges point toward the vertex carrying this slot. |

The connected interior has maximum degree three. Every leaf gets a spoke, a bivalent vertex can get one, and a trivalent vertex gets none. Keep precisely the resulting graphs satisfying $2|E'|\leq3(|V'|-1)-1$ on each subset of at least two vertices. This construction enumerates all admissible graphs. The replay removes isomorphic marked trees by their canonical rooted codes and isomorphic marked cycles by rotations and reversal of their cyclic list of rooted trees. An independent enumeration starts from all unmarked tight graphs through $n=3$. [Enumerator](rooted_graphs.py) · [Checks](check_operators.py)

| $n$ | Tree | One internal cycle |
|---:|---:|---:|
| 1 | 1 | 0 |
| 2 | 2 | 3 |
| 3 | 5 | 18 |
| 4 | 15 | 95 |

At $n=1$ the graph is the square rooted at a vertex. At $n=2$ suppressing bivalent vertices gives the theta graph. The allowed triples of path lengths are $(1,3,3)$ and $(2,2,3)$. Rooting at a trivalent endpoint gives the two trees. Rooting in a path of length three in the first graph or in a path of length two or three in the second gives the three internal cycles.

For an operator of length $l$, all $n>l$ vanish. For an integrated bilinear at positive loop order both root letters must contract, leaving three output momenta whose sum is zero. They span at most two holomorphic directions. Every cotree current lies in this span, so the $3n$-dimensional current Jacobian vanishes before Schwinger integration. This proves vanishing modulo total holomorphic derivatives for arbitrary bilinear jets in these cubic BF theories. Undifferentiated or differentiated spectator fields in longer operators can carry further momentum.

## Gauge bilinear

Take the invariant contraction $\operatorname{tr}(b\partial_i b)$, with $i$ a holomorphic derivative direction, and three outgoing $a$ fields of momenta $p_1,p_2,p_3$. Write $f$ for the Lie bracket paired with the third color by the trace in the action. In a checked orthonormal three-color basis the covariant structure tensor has $f_{012}=1$. For a simple Lie algebra write $C_A$ for the Casimir with two brackets contracted using the action pairing, $f_{imn}f_j{}^{mn}=C_Ag_{ij}$. This gives $C_A=2$ in the three-color check.

Each graph contracts five copies of the structure tensor. Join the two legs at the insertion using the trace pairing. The resulting color graph has five trivalent vertices, six internal edges and three external legs. A self-loop vanishes. A double edge gives $C_A$ times the pairing and a triangle gives $\pm C_A f/2$ by antisymmetrizing the adjoint trace. A simple triangle-free graph on five vertices with six edges must be $K_{2,3}$: a degree-three vertex has three independent neighbors, and the last vertex must join all three. Exchanging the two degree-three vertices reverses the three remaining tensors and kills the contraction. The same contractions on the residual one- and three-vertex graphs prove that every contribution is a universal multiple of $C_A^2f$.

There are no holomorphic derivatives in the gauge vertices. Passing from three independent basis momenta to $p_1,p_2,p_3$ multiplies the two-loop current Jacobian by $\det(p_1,p_2,p_3)^2$. The differentiated root contributes one linear combination of the momenta. Its three coefficients are fixed by the three derivative directions on the ordered basis. The replay gives the following values after division by $f\det(p_1,p_2,p_3)^2(p_1+p_2+p_3)_i$ at $C_A=2$.

| Graph, path containing the root | Root automorphisms | Coefficient |
|---|---:|---:|
| $(1,3,3)$, length three | 1 | $-1/180$ |
| $(2,2,3)$, length two | 2 | $7/2160$ |
| $(2,2,3)$, length three | 2 | $0$ |

Thus the complete reduced two-loop symbol is
$$
-\frac{C_A^2}{1728}\,
f\,\det(p_1,p_2,p_3)^2(p_1+p_2+p_3)_i .
$$
The zero at one loop holds already locally for this color contraction: joining the insertion legs leaves the triangle color tensor, and all three linear derivative coefficients on basis momenta vanish after the external permutations. At two loops the sum is a nonzero local total derivative. [Full-color contraction and exact simplex moments](kernels.py) · [Polarization and rank checks](check_operators.py)

## Ordered-jet checks

The replay takes the multiplier and ordinary BV superfields as odd. At each interaction it orders flags by incident edge number followed by the surviving field. The Wick sign moves the flags into ordered edge pairs followed by outputs. A root derivative contributes minus the current at the tail or plus the current at the head of the edge. At each Poisson or mixed vertex the numerator is the fiber determinant of the first two incoming momenta. Conservation makes the three cyclic pairs equal. The gauge tensor is totally alternating, and a mixed vertex pairs its two gauge colors. Root-edge assignments and output permutations are summed before dividing by the rooted automorphism order.

Take $p_1=(1,0,1)$, $p_2=(0,1,2)$, $p_3=(1,2,0)$ and differentiate the second root letter in the first holomorphic direction. Superscripts $(0),(1),(2)$ below label the three test colors. On the one-loop square the polarized coefficients are
$$
\begin{array}{c|c|r}
\text{root letters}&\text{output fields in momentum order}&\text{coefficient}\\ \hline
b^{(0)},\partial_1 b^{(0)} & b^{(0)},a^{(1)},a^{(2)} &-5/12\\
g,\partial_1 g & g,h,h &1/12\\
b^{(0)},\partial_1g & b^{(1)},a^{(2)},h &2/3
\end{array}
$$
For the two tree interiors at $n=2$ take the root letters $b^{(0)},\partial_1 b^{(1)},b^{(2)}$, the fourth momentum $p_4=(2,1,-1)$ and outputs $b^{(0)},a^{(0)},a^{(1)},a^{(1)}$. For $(1,3,3)$ and $(2,2,3)$ respectively the two graph coefficients are $6835/864$ and $355/54$.

For $g\partial_1g$ with three outgoing $h$ fields, the pure Poisson coefficients of the three internally cyclic graphs, in the order of the table, are
$$
\frac{8411}{1190700},\qquad
\frac{110939}{57153600},\qquad
\frac{72283}{57153600},\qquad
\text{sum }\frac{559}{54432}.
$$
Adding the gauge sector multiplies each graph by $1+\dim\mathfrak g$. With gravitational external fields, the gauge flags of every mixed vertex have one incoming and one outgoing edge, and any pure gauge vertex would require an additional incoming branch with an external gauge leg. Equivalently, counting gauge heads minus tails over the entire graph gives the number of pure gauge vertices, which must vanish. Thus gauge lines are disjoint closed cycles. On the unicyclic interior only its unique cycle can be colored, and its Kronecker contractions give $\dim\mathfrak g$. The same argument applies for any jets and any loop number in the $k=n$ channel with gravitational external and root fields. [Species-resolved sum and exact replays](check_operators.py)

These are local BV graph symbols before the common position-space normalization. The $\eta$ interactions in the anomaly-cancelling extension have additional propagators. Quantum cohomology of that extended theory requires its own differential.
