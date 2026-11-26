# Permanent scalar-jet classes on $\mathbb C^{3|2}$

Fix the gauge algebra $\mathfrak {sl}_N$ and a nonzero polarization $u^1,u^2$ of the two even scalar letters $g_1,g_2$. For $2\leq j\leq N$ take the trace in the defining representation,
$$
P_j=\operatorname{tr}\big((u^1g_1+u^2g_2)^j\big).
$$
Every polynomial in the $P_j$ and their holomorphic derivatives is a strict all-order cocycle. These representatives are linearly independent over $\mathbb C[[\hbar]]$ in local cohomology when expressed in distinct monomials of the invariant jets. The same map for translation-invariant integrated local functionals is injective after taking total derivatives on its domain. These statements concern cocycles and the additive cohomology. [Component conventions](../verify/vops.py) · [Graph support](../README.md)

## Differential and detection

In the project's heat normalization $Q_0g_i=i[c,g_i]$, so $Q_0P_j=0$ by cyclicity, and $Q_0$ commutes with holomorphic derivatives. In each connected graph for $Q_n$, $n\geq1$, the root has $n+1$ contracted fields and the interaction vertices form one tree after deletion of the root. All surviving superfields and all input functional derivatives therefore have the same odd coordinate $\theta$. Acting on any jet of $u^1g_1+u^2g_2$, the component derivative in each contracted slot carries $u^1\theta_2-u^2\theta_1$. Thus every graph contains
$$
(u^1\theta_2-u^2\theta_1)^{n+1}=0.
$$
This is an exact zero, independently of the heat prefactor, displacement polynomial and color tensor. Hence every spectral-sequence differential out of these representatives vanishes without counterterms.

Choose constants $v_1,v_2$ with $u^1v_1+u^2v_2=1$. On an arbitrary traceless matrix-valued formal holomorphic arc $H(z)$ set $g_i=v_iH(z)$, $c=b=0$. The four classical variations restrict to zero, including $[g_1,g_2]=0$ in $Q_0 b$. In every positive-loop graph, each of the $n+2$ surviving interaction fields restricts to $(v_1\theta_1+v_2\theta_2)H$ or a holomorphic derivative of it. The square of the common odd factor vanishes. Thus evaluating any $\mathbf Q$-boundary on these arcs gives zero coefficient by coefficient in $\hbar$.

On a diagonal arc put $h_N=-\sum_{i=1}^{N-1}h_i$ for the last entry. For $j=2,\ldots,N$ the Jacobian of the invariant values has entries
$$
\frac{\partial P_j}{\partial h_i}
=j(h_i^{j-1}-h_N^{j-1}),\qquad
\left|\det\frac{\partial(P_2,\ldots,P_N)}{\partial(h_1,\ldots,h_{N-1})}\right|
=N!\prod_{i<l}|h_i-h_l|.
$$
This follows by subtracting the last column of the ordinary Vandermonde matrix from the others and expanding its constant row. At a regular diagonal value this Jacobian is invertible. For every nonzero derivative multi-index, differentiating $P_j$ gives the same Jacobian on the highest-order derivative of $h$ plus terms at lower orders. In increasing total derivative order the finite jet map is block triangular with invertible diagonal blocks, for jets in all three holomorphic directions. Its image contains a neighborhood. A polynomial in the invariant jets which vanishes on every arc hence vanishes identically. Evaluating a putative incoming spectral-sequence differential gives the same contradiction. All these classes survive with their displayed representatives.

## Integrated functionals

For a polynomial density $F$ in the invariant jets, let the variational derivative in $P_j$ be
$$
\frac{\delta F}{\delta P_j}
=\sum_{\alpha\ {\rm jet\ multi\text{-}indices}}
(-\partial)^\alpha\frac{\partial F}{\partial(\partial^\alpha P_j)} .
$$
If $F$ is $\mathbf Q$-exact modulo total derivatives, its restriction to every diagonal arc is a total derivative. Varying $h_i$ and integrating by parts gives
$$
0=\frac{\delta F(P(h))}{\delta h_i}
=\sum_{\substack{2\leq j\leq N\\\mathrm{invariant\ generators}}}\frac{\partial P_j}{\partial h_i}
\left.\frac{\delta F}{\delta P_j}\right|_{P(h)} .
$$
Invertibility of the finite jet maps gives $\delta F/\delta P_j=0$ identically. Scaling all invariant jets by $s$ and integrating $dF(sP)/ds$ from zero to one gives, modulo total derivatives,
$$
F(P)-F(0)=\int_0^1 \sum_{\substack{2\leq j\leq N\\\mathrm{invariant\ generators}}}P_j\left.\frac{\delta F}{\delta P_j}\right|_{sP}\,ds=0.
$$
Evaluation on the zero arc sets the constant $F(0)$ of a boundary to zero. This proves injectivity on the integrated quotient with translation-invariant local primitives.

For an explicit differentiated class take $P_2=\operatorname{tr}(g_1^2)$ in $\mathfrak{sl}_2$. On $g_1=\operatorname{diag}(1,-1)(w+w^{-1})$, with $w=e^{z^1}$, the derivative $\partial_1$ acts as $w\partial_w$. Then
$$
P_2=2(w+w^{-1})^2,\qquad
(\partial_1P_2)^2=16(w^4-2+w^{-4}),\qquad
[w^0](\partial_1P_2)^2=-32 .
$$
Here $[w^0]$ extracts the constant Laurent coefficient. It annihilates $w\partial_w$, and derivatives in the other two directions vanish on this arc. Thus $(\partial_1P_2)^2$ survives in integrated cohomology as well as locally. The replay checks polarized input cancellations and output projection in all eight available graph operators through three loops, the Vandermonde determinant for $2\le N\le9$, and an explicit four-block mixed jet Jacobian. [Replay](check.py)

At scalar bidegree $(2,2)$ for $\mathfrak{sl}_3$, the independent mixed trace density requires one explicit counterterm and then closes to all orders. [Full sector and replay](mixed.md)

For every rank and bidegree, marked products lift all polynomial functions of commuting scalar pairs. The heat kernels have one explicit graph shape per loop and terminate by the smaller scalar degree. [Physical normalization, kernels and counting](scalar_products/README.md)

All field-length-two local cohomology classes of a complex simple gauge algebra survive every quantum differential, including both scalar–$b$ bilinears. [Current-algebra computation and dimension series](quadratic.md)

Every component of the full trace descendants and every iterated marked product has a finite all-order lift. At first holomorphic weight and odd weights $(1,1)$ the cubic cohomology has dimension one for $N\ge3$ and zero for $N=2$. [Root contraction, graph reduction and cubic sector](single_trace.md)

## Stable rank at fixed field length

For local finite-jet polynomials of field length $p$, the absolute Chevalley–Eilenberg complex includes the undifferentiated ghost. For $\mathfrak{gl}_N$ and $\mathfrak{sl}_N$ the perturbative spectral sequence satisfies
$$
E_1^p=E_\infty^p,\qquad N\ge p .
$$
The equality is of associated graded spaces. The marked products give finite lifts. The fixed-rank calculation of nonzero $Q_n$ concerns the action before passing to this cohomology.

**Trace complex.** At fixed length and source multidegree the cochains are finite-dimensional. Constant color conjugation acts trivially on classical cohomology by Cartan's insertion homotopy, so taking invariants computes absolute cohomology. The invariant multilinear color tensors on $p$ matrices are spanned by the traces indexed by permutations: each cycle gives one trace. They are independent for $N\ge p$ by the permutation action on $(\mathbb C^N)^{\otimes p}$, tested on $p$ distinct basis vectors. Graded antisymmetrization acts on the source slots as well. The classical differential replaces one letter by the product of two letters in its cyclic word with the commutator signs. Thus the formal trace algebra is the free graded-commutative algebra on the finite dual of the cyclic chain complex shifted by one. Its differential preserves each trace. In characteristic zero its cohomology is the free graded-commutative algebra on the shifted cyclic homology dual. This argument compares classical complexes in finite degrees and needs no quantum maps between different ranks. [Loday–Quillen, matrix Lie algebra homology](https://doi.org/10.1007/BF02566367)

For $N=p$, the first relation, at length $p+1$, is the alternating sum of permutation contractions. Indeed the kernel of the symmetric-group action consists of the matrix blocks with more than $N$ rows, here just the one-dimensional sign representation. The trace relation is symmetric under simultaneous relabelling of its color slots, which conjugates the permutations. At fixed homogeneous source inputs, complete symmetrization in the color slots kills every term in a classical coboundary, since each term contains a color bracket of two of the slots. Thus no nonzero formal coboundary can become this first relation. The kernel in degree $p$ and the image from degree $p-1$ coincide with the formal answers also at $N=p$.

**Primitives.** For $R=\mathbb C[z^1,z^2,z^3]\otimes\Lambda(\theta_1,\theta_2)$, the cyclic homology of degree $p-1$ at positive source weight is $\Omega_R^{p-1}/d\Omega_R^{p-2}$. To see this, resolve the diagonal over $R\otimes R$ one coordinate at a time. An even difference gives a two-term Koszul resolution. An odd difference gives the complex with a divided power in every nonnegative resolution degree; multiplication by the odd difference has image equal to its kernel. Tensoring with the diagonal makes these differentials zero. On Hochschild chains the map $u_0\otimes\cdots\otimes u_k\mapsto u_0\,du_1\cdots du_k/k!$ induces the resulting forms isomorphism: graded antisymmetrization of the coordinate differentials is its inverse on homology, including the divided powers for odd coordinates. It takes the Hochschild boundary to zero and Connes' cyclic operator $B$ to $d$. In the cyclic complex the lower form-degree summands of a cycle therefore give de Rham cohomology. At strictly positive weight contraction with the Euler vector field, divided by that weight, kills all those summands. At weight zero one copy of $\mathbb C$ survives in every even cyclic degree.

At positive weight $d$ identifies the primitive quotient with the image of $d$ in $\Omega_R^p$. Under the cyclic trace map its dual is generated by taking coefficients of $\operatorname{tr}((dC)^p)$: differentiating $\operatorname{tr}(C(dC)^{p-1})$ gives exactly this expression. Each such coefficient is killed by every graph differential. In fact two functional derivatives on $dC$ at the common odd point multiply either two identical odd linear factors or an odd linear factor by the top odd delta, both zero. The weight-zero classes are $\operatorname{tr}(c^p)$ for odd $p$; two derivatives on the undifferentiated superfields give the square of the odd delta. Their tree variations vanish by cyclicity. Taking marked products lifts all classical classes of length $p$ and of length $p-1$. Each $Q_n$ increases field length by one, so these lifts rule out both outgoing and incoming differentials. Finally $\mathfrak{gl}_N=\mathfrak{sl}_N\oplus\mathbb C\,1$ is orthogonal for the trace. All interacting vertices vanish on a central color, leaving the center as the free algebra of its jets. Restricting to $\mathfrak{sl}_N$ removes precisely the length-one primitive generators. [Strict representatives and marked products](single_trace.md)

**Count.** Fix a nonzero multidegree $\mu$ in the ordered coordinates $(z^1,z^2,z^3,\theta_1,\theta_2)$. Let $s$ count its nonzero entries, and set $r=\max(\mu_4-1,0)+\max(\mu_5-1,0)$. The polynomial counting forms by form degree $t$ is $t^r(1+t)^s$: an occupied even coordinate can contribute either no differential or one, and an occupied odd coordinate of weight $a$ contributes $a-1$ or $a$ commuting differentials. Exactness in positive weight divides this polynomial by $1+t$ to count the ranks of the outgoing de Rham maps. The number of independent primitive classes at length $p$ is therefore
$$
\binom{s-1}{p-1-r}
$$
with value zero outside the binomial range. Their parity as local polynomials is $p+\mu_4+\mu_5$ modulo two. At weight zero there is one primitive per odd length in the range above. The replay computes the graded cyclic boundaries themselves, checks their squares, and verifies 112 dimensions, including odd weights greater than one. [Exact cyclic ranks](stable.py)

For $\mathfrak{sl}_2$ the entire sector at zero holomorphic derivative degree survives, with a closed dimension series for all odd-coordinate weights and field lengths. [Koszul computation and exact replay](sl2_zero_jets.md)

For a simple gauge algebra the full integrated field-length-two sector consists of nine even and eight odd source-symmetry charges at every rank. [Integrated quadratic classification](quadratic.md#integrated-quadratic-charges)

The full cubic sector for every complex simple gauge algebra is computed by its trilinear invariant tensors. For $\mathfrak{sl}_2$, all field lengths and odd-coordinate weights at total holomorphic degree one are also determined. [All cubic classes](single_trace.md#all-cubic-local-classes-for-simple-gauge-algebras) · [First holomorphic derivative and exact ranks](one_jet/README.md)

The routine `local_dimension` multiplies the symmetric-algebra factors for even primitives and the exterior-algebra factors for odd primitives, truncated to the requested field length and five source weights. For $\mathfrak{sl}_N$ with $N\ge4$, the field-length-four multiplicities at weights $(0,0,0,2,2)$ and $(1,1,1,1,1)$ are $3$ and $24$, respectively.

## Integrated cohomology in the same rank range

For $N\ge p$, every classical class of field length $p$ modulo total holomorphic derivatives survives to all perturbative orders, for both matrix gauge algebras above. Closure modulo derivatives requires keeping the classical descent equations when constructing the quantum representatives.

Use the source differential $d$ and $\widehat Q_0=Q_0/i$, with $\widehat Q_0C=C^2$. For each integer $L\ge1$ take the homogeneous form components of the ordinary Chern–Simons transgression
$$
L\int_0^1\operatorname{tr}\!\left(C\,[t\,dC+t(1-t)C^2]^{L-1}\right)dt .
$$
Applying $d+\widehat Q_0$ gives $\operatorname{tr}((dC)^L)$. Indeed $t\,dC+t(1-t)C^2=(d+\widehat Q_0)(tC)-(tC)^2$, and the covariant differential $d+\widehat Q_0-[tC,-]$ kills this curvature. Under the trace the differentiated integrand is the $t$ derivative of the $L$th power of the curvature. Thus on the component of form degree $k$, the classical differential is $-i\,d$ of the component of degree $k-1$. Every monomial in these components uses only $C$ and the same $dC$ at a common superpoint. The products of odd delta factors in the [single-trace calculation](single_trace.md) therefore kill every $Q_n$, $n>0$, on these components individually.

Take one finite-jet dual copy of $\Omega_R^k$ for each $0\le k<L$, in field length $2L-1-k$, with differential $-i\,d^\vee$ to the preceding copy, using the graded transpose on coefficient functionals. Map a functional to its value on the corresponding transgression component. On the free graded-commutative algebra generated by these copies, send monomials to the graded symmetrization of left-associated marked products. The binary marked product obeys the chain Leibniz identity even on nonclosed inputs. This gives a chain map and it commutes with all three holomorphic translations; the Schwinger displacements depend only on differences. It preserves field length: an order-$m$ product adds $2m$ cubic vertices and contracts $3m$ physical propagators. The map reduces to the ordinary product at $\hbar=0$.

In the classical formal trace complex this map is a quasi-isomorphism in every length. At positive weight the truncated forms complex has cohomology only in its top copy $k=L-1$, where the representative is the indicated functional on $\operatorname{tr}(C(dC)^{L-1})$. At weight zero only its bottom copy contributes, a nonzero beta-integral multiple of $\operatorname{tr}(c^{2L-1})$. These are exactly the cyclic generators computed above. Length zero splits off as constants.

At positive length both complexes before taking cohomology are free modules over the polynomial algebra in the three translations. For each ordered product of factors their spatial jets are polynomial in the individual momenta. Change coordinates to total momentum and differences from the mean. Permutations fix the total momentum and act on the differences, so graded symmetrization and the trace identifications preserve this free factor. Quotient by total derivatives sets total momentum to zero. Resolving zero momentum by the three-variable Koszul complex consequently computes the ordinary quotient on both sides. Applying that bounded resolution to the classical quasi-isomorphism gives a quasi-isomorphism after the quotient as well.

Finally specialize the formal *classical* trace complex to rank $N$. In degrees through $p$ the matrix invariants are independent for $N\ge p$, also with coefficients modulo total momentum. At $N=p$ the only relation in degree $p+1$ still has an all-singleton trace component with coefficient one; a classical coboundary at this degree has no such component because inserting a bracket preserves the number of trace cycles. Hence the integrated cohomology at field lengths $p$ and $p-1$ agrees with its stable value. Evaluating the marked chain map directly at this finite rank lifts every model cocycle modulo derivatives. This excludes both directions of spectral-sequence differential. Removing the central fields gives the traceless case by the quantum chain retraction.

For $\mathfrak{sl}_N$, there are no length-one generators. At lengths two and three there are consequently no products of shorter generators. In the integrated primitive complex the quotient retains constant spatial Taylor coefficients of the forms; its differential acts only on the two odd coordinates. At length two there are 20 coefficients of one-forms and the incoming de Rham map has rank three, leaving 17 charges. At length three the two-form space has dimension 48, the one-form differential has rank 14, and one constant Cartan class comes from $\operatorname{tr}(c^3)$. The integrated cubic quantum cohomology thus has **16 even and 19 odd generators for $N\ge3$**. The verifier independently imposes zero total momentum on cyclic chains and recomputes 40 multigraded dimensions. [Noncommutative transgression identities and exact ranks](stable_integrated.py)

For a simple algebra outside type $A_{r\ge2}$ the invariant cubic color tensor is the bracket. On a root $\mathfrak{sl}_2$ two independent contractions of four colors remain independent even with momentum coefficients modulo total momentum. Thus each closed integrated bracket cochain on $\mathfrak{sl}_2$ extends to $\mathfrak{sl}_4$, and its quadratic primitives restrict. Conversely a quadratic primitive on the root algebra extends using the invariant pairing. In the matrix calculation the 34 two-form classes come from $\operatorname{tr}(C(dC)^2)$ with fully symmetric cubic color (the two $dC$ slots are total-even). They restrict to zero in this case. The unique integrated cubic class is the odd Cartan representative $\operatorname{tr}(c^3)$ in a representation with nonzero invariant pairing. [Finite-color channel calculation](single_trace.md)

Every local quantum spectral page splits additively by the constant-ghost cohomology. At every matrix rank the maximal relative length at each multidegree is the full Weyl-invariant algebra of diagonal scalar and first ghost modes, and all its classes survive. [Relative graph argument and exact ranks](relative.md)

The full sector with zero holomorphic derivatives, including every scalar, antifield and constant-ghost degree, is computed at every rank in types $A,B,C,D$, and for $G_2$. All reductive types have an exact classical Cartan description. [Koszul duality, Euler identity and dimension series](zero_jet_ranks.md)

The dressed binary bracket supplies the extra type-$D$ Pfaffian generators at zero holomorphic order and at maximal relative length. The finite-rank $\mathfrak{sl}_2$ calculation also certifies the full quantum cohomology in every multidegree of total weight at most ten. [Pfaffian graph operation](type_d.md) · [Finite-rank exact certificates](singular/README.md)

For every complex simple color algebra the box kills an undifferentiated scalar polynomial exactly when its mixed Hessian is symmetric. In the exceptional $F_4$ scalar sector $(4,4)$ all eight classes lift. Three are strict, four more have marked-product representatives, and the remaining Hessian has a classically exact bare box. [Scalar box, exact Weyl dimensions and Hessian representative](f4.md)

The mixed source Ward identities exclude every zero-holomorphic outgoing page above $\max(0,\min(p,q)-3)$ for the two odd-coordinate weights $(p,q)$. This proves global permanence at all ranks when either weight is at most three, including exceptional antifield sectors. [Opposite shears and exact selection rule](source_shears.md)
