# A source-symmetry bound for zero holomorphic derivatives

For any local classical class at zero holomorphic derivative degree, let $(p,q)$ be its two odd-coordinate weights. On the perturbative spectral sequence an outgoing differential $d_r$ can be nonzero only if
$$
1\le r\le \max\{0,\min(p,q)-3\}.
$$
This applies at every gauge rank and field length, including antifields and absolute constant-ghost factors. In particular every class with one odd-coordinate weight at most three has a finite quantum lift. There are no incoming positive-order pages at zero holomorphic degree.

## Two mixed source transformations

The source vector fields $\theta_j\partial_i$ preserve every edge kernel since the odd delta kills the shear of its body displacement. The opposite odd vector fields $z^i\partial_{\theta_j}$ preserve the classical action and BV pairing. Their action on perturbative pages follows from their integrated quadratic Hamiltonian
$$
\frac12\int\operatorname{Tr}\big(C\,z^i\partial_{\theta_j} C\big)
$$
in the action pairing and the same volume as the kinetic term. The divergence of the vector field is zero. Graded integration by parts gives its Hamiltonian variation on the field, up to the common ordering sign. The vector field commutes in the graded sense with the kinetic differential and is a derivation of the cubic vertex. These are the classical Ward conditions. The calculation uses the polynomial Taylor coefficient $z^i$ at the insertion and compactly supported test fields for integration by parts. [Quadratic Hamiltonian and component pairing](quadratic.md#integrated-quadratic-charges) · [Direct shear](type_d.md)

Add this Hamiltonian once to the interaction with a formal parameter of square zero. In any connected bulk graph with $V$ vertices and $L\ge1$ loops there are $V-1+L$ edge delta factors in each of the two odd coordinates. In the coordinate different from the one differentiated by the Hamiltonian, all these are untouched linear differences in at most $V-1$ independent odd generators. Their product is zero. This includes heat edges and subgraphs with external fields, and applies whether the differentiated leg hits a delta or a surviving field. Hence the deformed bulk master equation holds through first order in the parameter. The local insertion Ward identity then makes its first-order variation a filtered chain operator. Its tree part has a single edge between the Hamiltonian and the local insertion, giving the stated infinitesimal source action on jets. Additional interaction-dressed terms define its action at positive perturbative orders. The coefficient $z^i$ only requires the linear Taylor term at a body point in these local graphs. [Insertion Ward identity and bulk cancellation](../README.md)

The induced operator on page one is the classical generator. Filtered naturality induces the same relations on every subsequent page by passage to cycles modulo images. This needs no identification of the chain-level commutator of two dressed Hamiltonians. Work on the absolute complex, which already contains the relative classes by the constant-color splitting. The linear source transformations fix the origin, and their page-one actions preserve the split subspaces. The bracket of opposite mixed generators acts by the sum of their chosen spatial and odd-coordinate weights; the two changes at loop order $r$ cancel in this sum. [Constant-color splitting](relative.md) · [48 native mixed-generator and Ward tests](source_shears.py)

## The three spatial directions

Choose the scalar flavor of weight $q$. Denote by $E_i$ the odd generator taking the first ghost derivative in the $i$th spatial direction to this scalar, and by $F_i$ the reverse generator. Fix their signs so that $E_iF_i+F_iE_i$ is the sum of these two weight operators. For $i\ne j$ their anticommutator takes the $i$th spatial direction to the $j$th. It commutes with each $F_l$ for $l\ne i$. These statements follow by applying the two linear vector fields to the five source coordinates.

Every zero-holomorphic class is killed by every $E_i$ on page one, since applying $z^i\partial_{\theta_j}$ at the origin requires a spatial derivative on an input jet. The holomorphic $SL_3$ acts trivially on the class. Thus on any page reached by this class its image $v=d_r[F]$ is killed by all $E_i$ and off-diagonal spatial rotations. Its spatial weights are $(r,r,r)$ and its chosen odd weight is $q-r$.

For $k=1,2,3$ apply $F_k$ successively in distinct directions. Before $F_k$ acts, $E_k$ kills $F_{k-1}\cdots F_1v$: commuting it through the previous factors produces only off-diagonal rotations which commute through the other factors and kill $v$. The weight identity gives
$$
E_kF_k F_{k-1}\cdots F_1v
=(q-r-k+1+r)F_{k-1}\cdots F_1v
=(q-k+1)F_{k-1}\cdots F_1v .
$$
If $v\ne0$, the first $\min(q,3)$ successive factors must therefore give a nonzero vector. Its odd-coordinate weight $q-r-\min(q,3)$ is nonnegative, since every cochain has nonnegative source weights. Applying the same argument to the other flavor proves the bound. With no outgoing pages the usual lift recursion stops by order $\min(p,q)$, when a further term would have a negative odd-coordinate weight. [Exact graph weight shifts](relative.md)
