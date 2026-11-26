# Order-$g$ obstruction in the current complex

In the $g=0$ divergence-free quotient, $X$ is the odd vector field and $D_a$ is the holomorphic total derivative. The scalar $\beta'=\beta+X^a\gamma_a$ already appears in T5. The strict differential gives

$$
Q\beta'=-X^aD_a\beta',\qquad
Q\gamma_a=D_a\beta'-X^bD_b\gamma_a-(D_aX^b)\gamma_b .
$$

At $w=1$ integrate by parts all derivatives of the unique $\gamma$ or $\beta'$ letter. This gives a unique density $\gamma_aJ^a+\beta' K$ modulo total derivatives. Independence follows by taking the variational derivatives in $\gamma,\beta'$. Write $s=Q+X^aD_a$ on vector jet polynomials. Since $D_aX^a=0$, the induced differential on the two coefficients is

$$
(J^a,K)\longmapsto
\left(sJ^a-(D_bX^a)J^b,\ -D_aJ^a-sK\right).
$$

In particular the two-loop graph $X_1$ lies in the current term, has $K=0$ and satisfies $D_aJ^a=0$ by the longitudinal Ward identity. The ten interior momentum columns in $D_\theta$ give at least one derivative on each external field before integration by parts. Moving derivatives from $\gamma$ preserves that property for all nine vector fields. [Strict differential](../cohomology_11d.md) · [Momentum determinant](../../VERIFY.md) · [Longitudinal identity](../gamma_ward.md) · [Exact differential on 72 jet inputs](check_differential.py)

**Removing constant and linear jets.** Filter by the number of undifferentiated $X$ letters. The operator $s$ kills each undifferentiated $X$. On differentiated jets the undifferentiated terms in $Q$ cancel $X^aD_a$. The diagonal differentials are the Chevalley–Eilenberg differentials of the formal divergence-free vector fields vanishing at the origin, with coefficients in the tangent space or the scalars. On the associated graded, divergence differentiates only nonconstant jets. The $k$ undifferentiated letters contribute a spectator tensor $\Lambda^k\mathbb C^5$. The action of $\mathfrak{sl}_5$ on either diagonal cohomology is zero by Cartan's insertion identity for the linear vector fields. Therefore the invariant cohomology of the associated graded vanishes for $1\le k\le4$, by the exact sequence between the two diagonal complexes. For $k=5$, the spectator determinant has weight $+1$. At total degree $N=8,9$ the remaining tensors must have determinant weight $-4$ with at most four vector jet arguments. On decomposable jets each vector contributes only one distinct momentum to a determinant. At most one determinant can use the single covector from $\gamma$; the other determinants need five vector momenta each. The invariant cochain group itself is zero. Here the invariant tensor theorem allows opposite volume tensors to be eliminated in pairs. Consequently inclusion from $k=0$ induces an isomorphism on invariant degree-nine cohomology at the target weight.

The origin stabilizer is the semidirect product of $\mathfrak{sl}_5$ and the ideal of divergence-free fields vanishing through first order. This ideal acts trivially on either output fiber. At fixed jet weight its graded cochain spaces are finite-dimensional. Filtering the Chevalley–Eilenberg bicomplex by the number of arguments in the ideal gives $H(\mathfrak{sl}_5)$ tensored with the invariant relative cohomology. Indeed the invariant subcomplex is a direct summand as a complex of $\mathfrak{sl}_5$-modules, and all nontrivial isotypic components have zero $\mathfrak{sl}_5$ cohomology. The primitive exterior degrees of $H(\mathfrak{sl}_5)$ are $3,5,7,9$.

For relative degree $q$ count invariant cochains on $q$ positive jets, with one undifferentiated $\gamma$ (current) or one $\beta'$ (scalar). The sum of the positive jet degrees, Taylor order minus one, is $14$ or $15$ respectively. The determinant weight is $-3$ in both columns.

| $q$ | Current | Scalar |
|---:|---:|---:|
| 0 through 4 | 0 | 0 |
| 5 | 52 | 13 |
| 6 | 1,853 | 568 |
| 7 | 21,051 | 8,331 |
| 8 | 109,050 | 55,495 |
| 9 | 302,166 | 197,563 |
| 10 | 482,720 | 409,734 |
| 11 | 455,662 | 516,471 |
| 12 | 249,062 | 398,798 |
| 13 | 72,301 | 182,956 |
| 14 | 8,552 | 45,369 |
| 15 | 0 | 4,632 |

The Weyl projection in the replay uses unordered distinct 32nd roots of unity. After multiplication by the inverse determinant character and the Weyl factor, every individual Laurent exponent has absolute value below 32. Finite averaging gives the exact constant term. Reconstruction uses pairwise coprime primes whose product exceeds the unreduced dimension at the identity. Direct exterior monomial enumeration checks the same method in dimensions two, three and five at lower weights. [Character calculation](relative_dimensions.py) · [Independent check](check_dimensions.py)

**First differential.** For five transverse vector inputs with polarizations $u_i$ and momenta $p_i$, let $D$ be the determinant of the five ordered momenta. At relative degree five a scalar basis consists of $D^3$ times the orbit sums of products of five pairings, one pairing $p_{f(i)}(u_i)$ per input, with $f(i)\ne i$. Here $f$ labels a map of the five slots, up to simultaneous relabeling. For currents mark a slot $j$ and replace one $D$ by the determinant with the test covector in slot $j$. The remaining two determinants equal $D$. These give 13 unmarked and 52 marked candidates. A finite evaluation proves simultaneously that they are bases and that both degree-five differentials are injective. The six vector test inputs are $u_i(\exp(p_i(z))-1-p_i(z))$, so their values and first derivatives vanish at the origin. In each bracket term the code retains powers at least two in each of the two merged momenta. Integer minors of the evaluation matrices have residues

$$
\det_{13}=736898\pmod{1000003},\qquad
\det_{52}=746682\pmod{1000003}.
$$

The dimensions in the table give the matching upper bounds. Thus both relative degree-five cohomology groups vanish, and at least 52 current cochains at relative degree six are boundaries. [Determinant certificates](first_differential.py)

**Second differential.** On six arguments let $D_i$ be the oriented determinant omitting momentum $i$. For current outputs let $E_{ij}$ be the oriented determinant omitting $i,j$ and appending the test covector. Signed orbits of six transverse pairings times $D_iD_jD_k$ or $D_iD_jE_{kl}$ give scalar or current relative cochains. Retain terms with at least two momentum powers in each argument. Each seven-input row evaluates the differential on homogeneous vector fields of Taylor orders at least two. The seven orders sum to 22 for scalars and 21 for currents. In a fused bracket the polynomial uses momentum $s p_i+t p_j$ and polarization $t p_j(u_i)u_j-s p_i(u_j)u_i$; the two Taylor orders specify the powers of $s,t$ to extract.

| Output | Signed orbits | Test rows through last pivot | Rank modulo $251$ | Independent interpolation checks |
|---|---:|---:|---:|---:|
| Current | 5,210 | 3,882 | 1,801 | 20 |
| Scalar | 1,097 | 1,555 | 555 | 15 |

Adding signed columns in groups and evaluating on test chains can only decrease rank. A transverse input over the finite field lifts over the rationals by solving for one polarization coordinate using a nonzero momentum coordinate modulo $251$. The computed ranks thus bound the rational ranks from below. The first differential and the identity $d^2=0$ give matching upper bounds $1853-52=1801$ and $568-13=555$. Hence both diagonal relative cohomology groups vanish in degrees five and six. [Exact elimination with independent 25-point coefficient extraction](second_differential.cpp)

```sh
c++ -O3 -std=c++17 second_differential.cpp -o /tmp/sugra-relative-second
/tmp/sugra-relative-second current
/tmp/sugra-relative-second scalar
```

**Remaining complex.** At degree eight the only possible extra linear-jet term has $\mathfrak{sl}_5$ degree three and relative degree five. At degree nine it has $\mathfrak{sl}_5$ degree three and relative degree six. Both vanish in each diagonal complex. Thus the inclusion of relative cochains induces isomorphisms on the current and scalar cohomology in both degrees eight and nine. Comparing the long exact sequences of the two current complexes gives an isomorphism on the full degree-nine invariant obstruction group. In particular every class has a representative in which each vector has Taylor order at least two.

There are $23$ derivatives in total. For nine vector arguments in the $\gamma$ term, the excess over two derivatives per argument is $23-18=5$, with seven unordered profiles. For eight vector arguments in the $\beta'$ term it is $23-16=7$, with fifteen unordered profiles. The profile $2,2,2,2,2,2,2,9$ has zero invariant dimension, leaving fourteen. The character replay gives each individual dimension. At $w=1$ the remaining bound is

$$
\dim H^9(V)^{SL_5}_{w=1,\,t=(-3)^5}
\le 302166+55495=357661 .
$$

The current and scalar relative Euler characteristics at this weight are $5$ and $-2$. Since the scalar term is shifted by one, the relative complex has Euler characteristic $7$. The cohomology in degree nine and the projection of the two-loop graph onto it are open.

### Exact elimination of divergence at degree nine

For one positive vector jet, a polynomial in its covector momentum $p$ is linear in its vector polarization and has momentum degree at least two. The transverse relation is $p(u)=0$. Over $S=\mathbb C[p_1,\ldots,p_5]$ the module of these polynomials is the degree-at-least-two submodule of $S^5/Sp$. It embeds into a free module: send the coefficient vector $v$ to $p\wedge v$. The kernel on $S^5$ is $Sp$, since $p_i v_j=p_j v_i$ and the indeterminates $p_i,p_j$ are coprime. This also works for the degree-at-least-one jets of a vector field vanishing at the origin.

On $q$ slots the total derivatives $D_i$ multiply by the sums $P_i$ of their respective momenta. If one tensor factor is free on its slot, the entire tensor product is free over $\mathbb C[P_1,\ldots,P_5]$: replace the independent momenta of that slot with $P$ minus the sum on the other slots. Algebraically the binomial expansion of this translation defines mutually inverse polynomial module maps even when the other slots have relations. Embed each of the $q$ slot modules in a free one as above. The resulting $q$ short exact sequences, tensored over $\mathbb C$, shift positive $\operatorname{Tor}$ degree by one per slot. Explicitly the degree-$i$ Tor of the tensor product equals degree $i+q$ for the tensor product of the $q$ cokernels, for every $i>0$. The Koszul resolution of the residue field of $\mathbb C[P_1,\ldots,P_5]$ has length five. For $q\ge5$ all positive Koszul homology vanishes; antisymmetrization is a commuting idempotent, so the same holds on alternating jets. Restriction to the determinant isotype preserves exactness. In particular the kernel of divergence is exactly the image of contraction by $P$ on antisymmetric two-tensor coefficients.

For the scalar Euler degree $15$, the Koszul term in exterior degree $j$ consists of $q$ jets of Euler degree $15-j$ tensored with $\Lambda^j(\mathbb C^5)^*$. Its dimension in the determinant isotype is obtained by the same Weyl integral, inserting the $j$-th elementary symmetric polynomial in the inverse torus eigenvalues. Increasing $j$ lowers the jet degree by $j$ and subtracts at most one from each coordinate weight. The 32-point grid thus retains the earlier exact constant-term bound. Integer reconstruction uses the full exterior dimension at the identity times $\binom5j$. The scalar cokernel dimension is the alternating sum of the six multiplicities. [Exact character replay](divergence.py) · [Independent Koszul matrices at four small weights](check_divergence.py)

| vector jets | scalar | current | exterior degrees 2, 3, 4, 5 | rank of divergence | kernel | cokernel |
|---:|---:|---:|---|---:|---:|---:|
| 7 | 8331 | 21051 | 21735, 11246, 2860, 280 | 7982 | 13069 | 349 |
| 8 | 55495 | 109050 | 86438, 33637, 6234, 427 | 50442 | 58608 | 5053 |
| 9 | 197563 | 302166 | 181993, 51948, 6624, 271 | 165768 | 136398 | 31795 |
| 10 | 409734 | 482720 | 214923, 42586, 3345, 59 | 307097 | 175623 | 102637 |

Filter the relative current/scalar complex by its number of vector arguments. At fixed arity the first differential is precisely divergence, so total degrees $8,9,10$ have dimensions $58957,141451,207418$ after this differential. In particular
$$
\dim H^9(V)^{SL_5}_{w=1,\,t=(-3)^5}\ \le\ 136398+5053=141451.
$$
The scalar cokernel at arity eight can be evaluated at $P=0$. Indeed, embedding all eight slots into free modules one at a time gives cokernel filtrands with seven embedded-or-free slots and one quotient slot. The preceding Tor shift in the other slots kills their degree-one Tor with the residue field. Restriction modulo $P$ therefore injects the scalar cokernel into the ambient free polynomial quotient. The original graph is divergence-free before removal of linear jets. Determining its final cohomology class still requires the CE differential and the induced map from the scalar cokernel to the current kernel in the filtered complex.

The scalar cokernel at arity six has dimension nine and its relative CE differential has rank nine. For six balanced momenta the six oriented five-minors coincide up to the orientation sign. Before imposing balance, choose three minors omitting distinct slots. Their product has at least two momentum powers in every slot. Multiply by one transverse pairing per polarization and alternate. After imposing balance this is the common determinant cubed times a symmetric orbit of pairings; there are 40 orbits of maps without fixed slots. Overall nonzero orbit multiplicities do not affect rank.

On balanced transverse inputs the scalar cokernel differential is obtained by merging two plane waves with the usual commutator. Indeed a discarded constant test vector gives total momentum of the other inputs contracted with its polarization, zero by balance and transversality. A linear test vector acts on the other momenta and polarizations by the traceless matrix $u_i p_i$, zero by invariance. A term in a split bracket cannot have Taylor order at most one in both new slots: the original slot has order at least two and its bracket adds one derivative. These exhaust the truncation corrections. A $9\times9$ evaluation minor on seven balanced waves has determinant $198$ modulo $251$. The rational lift keeps the last momentum balanced and solves one transverse polarization coordinate on each wave. This proves injectivity over the rationals. [Balanced-wave certificate](balanced_six.py)

### The current profile with eight quadratic fields and one seventh-order field

Write $R_j$ for the divergence-free fields of Taylor order $j+1$, as in the jet character calculation. A test chain on the profile $(2,2,2,2,2,2,2,2,7)$ takes values in the determinant component of $\Lambda^8R_1\otimes R_6\otimes (\mathbb C^5)^*$. This space has dimension $175$. Bracketing two quadratics or bracketing a quadratic with the seventh-order field gives degree-eight cochain dimensions $321$ and $8$, respectively. [Profile characters](extreme_character.py)

To construct the chains as intertwiners into $\Lambda^8R_1$, dualize the final two factors and tensor with $\det^{-3}$. The high jet has highest weight $(1,0,0,0,-7)$, from the trace-free part of $\mathrm{Sym}^7(\mathbb C^5)^*\otimes\mathbb C^5$. Tensoring it with the test covector, dualizing and shifting by $\det^{-3}$ gives three highest weights with multiplicity one. Their multiplicities in $\Lambda^8R_1$ are $14,145,16$. At output arity seven the corresponding three weights increase the first coordinate by one and have multiplicities $0,8,0$. In the character replay the eighth and seventh exterior powers give totals $175$ and $8$, already bounded below the finite-field prime by the character calculation above. Every component multiplicity is nonnegative, so these lifts are unique.

| Highest weight | Test chains | Independent quadratic-bracket images modulo $251$ |
|---|---:|---:|
| $(5,-3,-3,-3,-4)$ | 14 | 13 |
| $(4,-2,-3,-3,-4)$ | 145 | 144 |
| $(4,-3,-3,-3,-3)$ | 16 | 16 |

Use the basis of divided-power monomial vector fields, subtracting a coordinate-5 monomial of the same divergence whenever necessary. The verifier generates highest vectors by the ordered product over positive roots. In the factor for $i<j$, its term of order $n$ on weight $\lambda$ is
$$
\frac{(-1)^n E_{ji}^n E_{ij}^n}{n!\prod_{r=1}^n(\lambda_i-\lambda_j+j-i+r)} .
$$
It checks all four simple raising maps on every sampled highest vector. The characteristic-zero rank bound also follows from integral base change. Tensoring $R_1$ by $\det^2$ gives a polynomial $GL_5$ module of degree $9$, and its eighth exterior power has degree $72<251$. The degree-$72$ Schur algebra in this characteristic is semisimple by Schur–Weyl duality and Maschke's theorem. Its Weyl characters and highest multiplicities agree with characteristic zero. Thus the integral highest-weight kernels reduce onto the sampled spaces. The dual-bracket coefficients are computed over $\mathbb Q$ by decomposing the action on the high input into words in the four lowering generators. The combined quadratic and high brackets have rank at least $174$. [Exterior calculation](extreme_profile.cpp) · [Rational dual bracket](extreme_fusion.py)

The remaining vector has highest weight $(5,-3,-3,-3,-4)$. Its $32{,}676$ nonzero exterior coefficients lift from the primes $251,257$ to integers of maximum absolute value $27{,}840$. The least exterior word has coefficient $1{,}568$ in this certificate. Exact integer evaluation kills every simple raising operator and all $51{,}474$ occupied coordinates of the quadratic bracket. The mixed bracket gives zero on every coordinate of its three highest-weight inputs, whose weight-space dimensions are $200,2870,3680$. These are all irreducible summands of the input to the intertwiner. Hence the two degree-eight boundary projections together have **rank $174$ over $\mathbb Q$**. The shared replay reconstructs the vector in a temporary directory. The optional two-column output records eight increasing zero-based basis labels packed from left to right into seven-bit words, followed by the signed integral coefficient. [Integer check](extreme_cycle.cpp) · [Reconstruction](replay_extreme.py)

**Pairing with the integrated graph.** Replace the high vector and the covector by a single scalar-test vector of Taylor order eight, retaining the same highest vector on the eight quadratics. Its dual irreducible representation has highest weight $(5,-3,-3,-3,-4)$ after the determinant shift. The quadratic bracket still vanishes. The other bracket would map it to a scalar chain on seven quadratics and one ninth-order vector. The required multiplicity in $\Lambda^7R_1$ is zero, the first entry of $(0,8,0)$. This scalar chain is therefore closed.

In the relative two-column complex the boundary of a pure scalar chain is the negative dual divergence in the current column plus the negative scalar CE boundary. Dual divergence lowers a Taylor order by one and appends the derivative covector. All eight quadratics then have zero contribution, since an order-one field drops out of the positive-jet complex. On the high slot the vector $(z^1)^8(e_5)^a/8!$ gives $(z^1)^7(e_5)^a(e^1)_b/7!$. This map is nonzero on the highest test, hence its invariant completion is exactly the one-dimensional current cycle up to an overall rational factor. Its scalar component is zero. Thus the surviving diagonal current cycle is a boundary in the full relative two-column complex.

Every relative closed representative pairs to **zero** on that cycle. For the raw integrated order-$g$ graph, the degree-nine isomorphism above gives a relative closed representative differing by a full coboundary. In the current column that difference is the diagonal CE differential. Evaluation of this differential on the positive test fields uses only their brackets, because their values and first derivatives vanish at the origin. The same integer cycle annihilates it. This establishes zero for the raw integrated graph pairing without projecting its pointwise Ward identity through the linear jets.

Consequently a degree-eight current primitive removes the entire extreme profile of the integrated graph. Splitting the seventh-order slot of the first primitive gives the other profiles $(2^7,3,6)$, $(2^6,3,3,5)$, $(2^6,3,4,4)$. Splitting the eighth-order slot of the second gives $(2^7,3,6)$ and $(2^7,4,5)$ in addition to the extreme component. Here a superscript on a Taylor order counts repeated inputs. A primitive also changes the scalar cochain by its divergence. The six other degree-nine current profiles and the scalar column remain in the obstruction calculation.

```sh
python3 extreme_character.py
python3 replay_extreme.py --rank
```

### Scalar cycles on the remaining current profiles

For scalar tests on the first two profiles $(2^8,8)$ and $(2^7,3,7)$, the source dimension is $14+495=509$. Their brackets have profiles $(2^5,3,3,7)$, $(2^6,4,7)$ and $(2^6,3,8)$ of dimensions $344,76,22$. The remaining possible output $(2^7,9)$ has dimension zero. On the middle $76$-dimensional space the quadratic bracket has rank at least $53$. Therefore the preceding scalar chain space contains at least
$$
509-(344+76+22)+53=120
$$
cycles. Exactly one cycle lies on $(2^8,8)$. Lowering the high field on $(2^7,3,7)$ injects its scalar tests into the $(2^7,3,6)$ current tests, since the degree-six slot is unique and the gradient embeds a positive homogeneous jet. Lowering the cubic contributes at $(2^8,7)$. Once the extreme current has been removed, the next $4567$ coefficients thus obey at least **$119$ independent homogeneous equations**.

More generally select a set of scalar profiles and take the union of their possible bracket outputs. The rank of the scalar boundary is bounded by the dimension of this union minus any certified rank for its next boundary. On profiles with no cubic field, bracketing two quadratics creates a unique cubic. Distinct such profiles have independent output projections. In the determinant component the ranks on $(2^6,4,7)$, $(2^6,5,6)$ and $(2^5,4,4,5)$ are at least $53,98,1076$. The calculation projects $\Lambda^qR_1$ onto its highest weights, applies the integral quadratic bracket, and multiplies each rank by the multiplicity in the untouched higher jets. A signed projection of the output coordinates reduces the size of the elimination. [Character and support calculation](profile_boundaries.py) · [Highest-weight rank checks](quadratic_rank.cpp)

| Scalar profiles admitted cumulatively | Independent scalar cycles, lower bound |
|---|---:|
| $(2^8,8)$ | 1 |
| add $(2^7,3,7)$ | 120 |
| add $(2^7,4,6)$ | 120 |
| add $(2^6,3,3,6)$ | 3033 |
| add $(2^6,4,4,4)$ | 3033 |
| add $(2^5,3,3,3,5)$ | 11975 |

Dual divergence is injective on this cumulative space. On each new profile its high derivative gives, respectively, the current profiles $(2^8,7)$, $(2^7,3,6)$, $(2^7,4,5)$, $(2^6,3,3,5)$, $(2^6,3,4,4)$ and $(2^5,3,3,3,4)$. All other derivatives land on earlier current profiles or vanish. On three equivalent quartic inputs the derivative inserts one cubic and retains the other two quartics. Composing with a left inverse of the one-field derivative and wedging the three quartics gives three times the identity. Thus every scalar cycle counted in the table gives an independent total boundary, on which the integrated current pairs to zero. The numbers count equations jointly across the indicated current profiles.

For the $8105$-dimensional current profile $(2^7,4,5)$ the quadratic-bracket verifier tests all $998$ directions with first highest-weight coordinate at least four and obtains rank at least $867$. Its primitives have profile $(2^5,3,4,5)$. Splitting their cubic contributes only to $(2^7,4,5)$, and splitting either higher field contributes to later profiles in the table. These $867$ eliminations preserve both earlier current profiles. Their overlap with the scalar-cycle equations is uncomputed.

### Seven balanced scalar inputs

The scalar differential after quotienting by divergence has ranks $9$ and $340$ on six and seven vector jets, respectively. On seven slots take three ordered five-momentum determinants and one contraction of each transverse vector with another slot's momentum. If slot $i$ occurs as a target $s_i$ times and is omitted from $r_i$ determinants, its Taylor order is $3+s_i-r_i$. Requiring this to be at least two puts every constructed cochain in the positive-jet complex before imposing balance.

For each determinant attach the sign of its five labels followed by its two omitted labels in increasing order. Extending the resulting function antisymmetrically in those two labels, a permutation of all seven slots supplies its sign times the permuted function. Three determinant factors therefore supply the full exterior sign. The remaining evaluation is an ordinary permutation sum over the three ordered pairs of omissions and the directed contractions. The 100 orbits of maps without fixed inputs permit the same products of seven contractions to be reused for many omission terms.

Exact evaluation on balanced transverse plane waves gives rank $349$ before the differential. For eight inputs merge each pair with the plane-wave Lie bracket and its ordinary Chevalley–Eilenberg sign. The differential has rank at least $340$ at each of the primes $251$ and $257$; the respective products of nonzero pivots are $2$ and $32$. The polarizations are signed four-by-four minors orthogonal to four test momenta. Consequently every finite-field test lifts to rational balanced transverse data by lifting the first seven momenta, balancing the last, and taking the same minors. The upper bound follows from the rank-nine incoming differential and $d^2=0$. Thus
$$
349-9-340=0,\qquad
\dim H^8(\mathrm{scalar\ cokernel})\le 5053-340=4713 .
$$
At the next page of the divergence filtration this gives the order-$g$ degree-nine upper bound $136398+4713=141111$ before the outgoing differentials. [Two-prime evaluation certificate](balanced_seven.cpp)

### Counterterm on block-diagonal vector backgrounds

Restrict the nine vector inputs of the order-$g$ defect to divergence-free fields on $\mathbb C^2\times\mathbb C^3$ of the form $X_{(2)}(z^1,z^2)+X_{(3)}(z^3,z^4,z^5)$. Allow arbitrary one-form and scalar test fields. The integrated graph is exact on this restriction. By the degree-nine relative reduction above it suffices to contract its relative representative. The restriction commutes with both the CE differential and divergence: mixed brackets vanish, and a derivative of either vector block in a coordinate of the other block vanishes.

The weight $\det_5^{-3}$ restricts to $\det_2^{-3}\det_3^{-3}$. In either the scalar column or the current with output in the three-dimensional block, the planar positive jets have total Euler degree six. For a current with planar output this degree is five, with the planar fundamental representation in the coefficient. In coordinates $x,y$ on the plane, a divergence-free vector field of Taylor order $h+1$ is the Hamiltonian field of a homogeneous binary form of degree $h+2$. Choose its action on a polynomial $G$ as $\partial_x F\,\partial_y G-\partial_y F\,\partial_xG$. For binary forms $F,G$ set $(F,G)_r$ equal to the evaluation of $r$ copies of $\partial_x\otimes\partial_y-\partial_y\otimes\partial_x$ on $F\otimes G$, followed by ordinary multiplication. Thus the Hamiltonian bracket is $(F,G)_1$.

In the planar scalar complex only degrees two and three occur, each of dimension one. The degree-two generator $S$ on two quintics is $(F,G)_5$. Its differential on $x^3,y^4,x^3y^2$ equals $17280$. The contraction sends $dS$ to $S$. In the planar current complex the dimensions in degrees two, three and four are $1,2,1$. With arguments in ascending degree, use
$$
P(F_4,G_5)=(F_4,G_5)_4,\qquad
A(F_3,G_3,H_5)=((F_3,G_3)_1,H_5)_4,\qquad
B(F_3,G_4,H_4)=(F_3,(G_4,H_4)_3)_2 .
$$
The subscripts indicate polynomial degrees; extend by alternation. The planar symplectic form identifies the linear-polynomial outputs with planar vectors. Define $T=dB$; direct Poisson brackets on the full monomial bases give
$$
dP=-A+7B,\qquad dA=7T,\qquad
T(y^3,xy^2,x^2y,x^3y)=-2592y .
$$
The homotopy is $h(A)=-P$, $h(B)=0$, $h(T)=B$, $h(P)=0$. It obeys $dh+hd=1$. The dimensions follow by taking the weight-zero minus weight-two multiplicities for scalar outputs and weight-one minus weight-three multiplicities for current outputs in each product of exterior powers of binary forms. [Dimension and monomial certificates](block_homotopy.py)

Tensor the appropriate planar homotopy with the identity on the three-dimensional block. Its anticommutator with the three-dimensional CE differential vanishes by the graded tensor-product sign. For a closed relative pair $(J,K)$ write $h_C,h_S$ for these current and scalar contractions. Then a primitive of the restricted graph is the pair
$$
\left(h_CJ,\ -h_S\left(K+\operatorname{div}(h_CJ)\right)\right).
$$
Indeed $d_C h_CJ=J$ and $d_S(K+\operatorname{div}(h_CJ))=0$, so applying the two-column differential reproduces $(J,K)$. The formula works for every linear change of the $2+3$ splitting. A nonzero full obstruction therefore requires vector-field tests that mix the coordinate blocks. The conclusion concerns the class: individual graph coefficients on block backgrounds can be coboundaries.

### One mixed planar vector

Let $x$ be the two planar coordinates and $y$ the three transverse coordinates. Add one horizontal perturbation $\epsilon f(y)Z(x)$ to the previous block background, where $\epsilon^2=0$, $Z$ is a planar Hamiltonian vector and $f$ is arbitrary. Take linear combinations of such perturbations. This is a Lie module: the planar bracket acts on $Z$ and a transverse vector $Y$ acts by $Y(f)$. It tests an independent mixed input in the order-$g$ cocycle.

After subtracting a primitive of the separated part, filter the linear-in-$\epsilon$ cochains by the number $s$ of derivatives on $f$. Use filtration degree $s+1$ for a current with transverse output and degree $s$ for the scalar and planar output. A transverse divergence differentiating $f$ preserves this degree; differentiating a transverse vector lowers it by one. Inserting $Y(f)$ for a positive transverse field $Y$ of Taylor order $a+1$ replaces $s$ by $s-a$. The ordinary vertical bracket preserves the degree. Thus cochains of bounded filtration degree form subcomplexes. Possible terms from the marked field acting on the unmarked output use at most its value and first derivative; the separated pair has already been made zero. At filtration degree at least three all marked jets have $s\ge2$, so every planar Taylor order of the mixed field, including constant and linear, is allowed in the relative five-dimensional jets.

At fixed transverse slots the planar test chain has one marked binary form, with every unmarked planar Hessian or higher jet alternating as before. If $h$ denotes planar Taylor order minus one, marked $h$ starts at $-1$ and unmarked $h$ at $1$. The Euler totals remain six for scalars and transverse currents, five for planar currents. Bracketing two binary monomials $x_1^i x_2^{h+2-i}$ and $x_1^j x_2^{t+2-j}$ gives coefficient $i(t+2)-j(h+2)$; a bracket containing the mark retains it. The following are the full planar homology multiplicities.

| Planar output | Arity 3 | Arity 4 | All other arities |
|---|---:|---:|---:|
| Scalar | 0 | 1 | 0 |
| Planar vector | 1 | 1 | 0 |

Dual planar divergence maps the unique scalar four-class onto the unique planar-current four-class. The exact calculation imposes the scalar boundary and the $\mathfrak{sl}_2$ raising operator on weight-zero chains, obtaining nine invariant cycles before boundaries. Apply the two coordinate derivatives to the Hamiltonian slots and take the quotient by current boundaries. The induced image has dimension one. Derivatives taking ordinary quadratics out of the positive-jet complex drop out; the marked slot includes all planar orders. [Exact planar differential and induced divergence](mixed_planar.py)

In the horizontal current–scalar cone the arity-four planar-current class consequently removes the scalar class. If transverse divergence also hits that scalar, the planar-current class supplies its unique compensating coefficient in cohomology. The remaining horizontal classes are the planar-current line at arity three and the transverse-current copy of the scalar line at arity four. In total degree nine these need, respectively, six or five ordinary transverse vectors. Their vertical Euler degrees are $9-s$ and $8-s$. Since each ordinary transverse vector uses at least one Euler degree, filtration degrees above four have no candidate.

For the two remaining high filtration levels the ordinary three-dimensional CE calculation is as follows. The symmetric power of the transverse momentum of the marked function is included in the highest-weight condition. The triples are highest weights for the transverse vector representation on test chains. The ranks in the table refer to the chain boundary out of and into the indicated arity.

| Filtration degree | Output | Transverse Euler degree and arity | Highest weight | Dimension | Outgoing | Incoming |
|---:|---|---|---|---:|---:|---:|
| 4 | transverse | $(5,5)$ | $(1,-3,-3)$ | 3 | 3 | 0 |
| 4 | transverse | $(5,5)$ | $(0,-2,-3)$ | 6 | 6 | 0 |
| 3 | planar | $(6,6)$ | $(0,-3,-3)$ | 7 | 7 | 0 |
| 3 | transverse | $(6,5)$ | $(0,-3,-3)$ | 21 | 14 | 7 |
| 3 | transverse | $(6,5)$ | $(-1,-2,-3)$ | 24 | 19 | 5 |

At fixed filtration degree the bracket of $Y$ with $f$ has lower degree, so the vertical differential induced on these horizontal groups is precisely the ordinary transverse jet differential. For a transverse output the weights follow from $\det^{-3}\operatorname{Sym}^s(\mathbb C^3)\otimes\mathbb C^3$ and for a planar output from $\det^{-3}\operatorname{Sym}^s(\mathbb C^3)$. The replay evaluates the monomial brackets on each Weyl weight at primes $251$ and $257$. To compare multiplicities over the rationals, tensor the chain spaces with $\det^{12}$. At Euler degrees five and six the result is polynomial of degree $36-5$ or $36-6$, below either prime, since each term has at most twelve covariant slots. Tensor powers in these degrees are semisimple by the symmetric-group idempotents, whose denominators divide the degree factorial. The differential preserves the isotypic projections. Finite-field ranks on each multiplicity space are therefore lower bounds in characteristic zero; outgoing and incoming ranks already add to the dimension. [Two-prime transverse calculation](mixed_vertical.py)

All degree-nine associated-graded cohomology at filtration degree at least three therefore vanishes. Successive local counterterms on this **one-mixed-vector restriction** remove those leading terms of the integrated defect. A representative requires at most two derivatives of $f$ in the scalar or planar current and at most one in the transverse current. The class at the remaining low transverse orders is undetermined.

### The antisymmetric transverse-current class

Project the one-mark complex to the transverse current and take the nonconstant part of the marked function. On tests at the origin this means $f(0)=0$. Every such marked field vanishes at the origin, including $f$ linear and $Z$ constant. That last field is a linear shear in the full five-dimensional algebra. Keep it as a test in the restricted complex relative to the block linear algebra. An ambient relative cochain extends to this test by its actual value, which is zero for an ambient relative representative. Restricted counterterms can take nonzero values on the shear. The transverse projection is a cochain map here: the marked vector has zero transverse component, its value at the origin vanishes, and the block linear actions are separate. The planar marked module therefore includes $h=-1$ at every $s\ge1$. The transverse differential is the tensor sum of the planar marked differential and the three-dimensional differential with the marked function acted upon by $Y(f)$. Both the constants and the functions vanishing at the origin are preserved by the positive transverse vector fields.

At filtration degree two the planar-current component has $s=2$ and the transverse current has $s=1$. For transverse Euler degree seven the diagonal computations give

| Output | Vertical arity | Highest weight | Dimension | Outgoing chain rank | Incoming chain rank | Homology |
|---|---:|---|---:|---:|---:|---:|
| Planar | 6 | $(-1,-3,-3)$ | 28 | 22 | 6 | 0 |
| Transverse, symmetric part | 5 | $(-1,-3,-3)$ | 51 | 29 | 22 | 0 |
| Transverse, alternating part | 5 | $(-2,-2,-3)$ | 33 | 17 | 15 | 1 |

For the alternating weight the only nonzero chain spaces, in arities $3,4,5,6,7$, have dimensions $(4,21,33,18,3)$ and outgoing ranks $(0,4,17,15,3)$ modulo $251$ and $257$. Tensoring with $\det^{14}$ gives polynomial degree $35$, below either prime. The ranks force all homology outside arity five to vanish over the rationals. The Euler characteristic $-1$ then fixes the rational arity-five multiplicity to one. The two transverse indices contract against the derivative of $f$ and the test one-form. The divergence term containing two derivatives of $f$ vanishes on their alternating part. [Diagonal rank certificates](mixed_vertical.py)

The coupling between different $s$ in the projected transverse current can be computed over the rationals. Take a primitive on four ordinary transverse vectors and $1\le s\le3$. Its three spaces have dimensions $57,49,26$. On five vectors the corresponding dimensions are $84,45,9$. Alternate products of three momentum determinants and one contraction per ordinary vector, inserting the test one-form once. Keep at least two momentum powers on each vector and exactly $s$ on the function. Evaluation gives bases in all six spaces. On a fused pair of ordinary wave vectors use their ordinary Lie bracket, subtracting the constant and first Taylor coefficient in each of the two input momenta. The two subtractions have disjoint monomial support since the combined momentum degree is at least three. For the bracket acting on the function only the Taylor orders zero and one of the vector input need subtraction.

Two hundred integer five-vector evaluations have rank $138$ on the target, verified modulo $251$ and $257$. Fraction-free elimination of the integer evaluations of the full four-to-five differential gives
$$
\operatorname{rank}_{\mathbb Q}d=88
=(29+17)+(14+19)+(3+6).
$$
The right side is the sum of the three diagonal ranks. Consequently the filtered image contributes no further boundary on the alternating line at $s=1$. The four-vector space at $s=4$ has dimension five, with zero diagonal cohomology, so one can subtract a preceding boundary to remove the leading term of any primitive there. Higher $s$ have no four-vector space. A five-vector cochain with one derivative on $f$ cannot create a term with undifferentiated $f$ through $Y(f)$. Tensoring with the one-dimensional arity-four planar marked scalar cohomology therefore gives **exactly one degree-nine transverse-current class with nonconstant marked function** in this restricted projection. Its coefficient for the actual graph, and its extension through the scalar and planar currents, remain to be determined. [Exact filtered differential](mixed_function.cpp)

For a graph test, the planar dual class has a representative on a single Taylor profile. The three ordinary planar fields have orders $2,2,3$ and the marked planar field has order $3$. Imposing the planar boundary and raising operator at this profile produces a cycle outside the full arity-five image. Its primitive integer representative has $24$ terms and maximum absolute coefficient $126$. The function `graph_test()` returns the ordered monomials and their coefficients. All four planar fields in this probe lie in the ambient positive jets. [Planar graph test](mixed_planar.py)

```sh
c++ -O3 -std=c++17 -I/opt/homebrew/opt/gmp/include -L/opt/homebrew/opt/gmp/lib mixed_function.cpp -lgmpxx -lgmp -o /tmp/sugra-mixed-function
/tmp/sugra-mixed-function 251
/tmp/sugra-mixed-function 257
```

### Evaluation of the graph on the mixed-current class

The planar test reduces to two products of powers of linear forms. On unrestricted homogeneous Hamiltonians of degrees $3,3,4,4$, alternating the two cubics, the invariant space has dimension two. Indeed $\Lambda^2\mathrm{Sym}^3=\mathrm{Sym}^4\oplus\mathbf1$, and both summands occur once in $\mathrm{Sym}^4\otimes\mathrm{Sym}^4$. Evaluate the two invariants $(F,G)_3(H,I)_4$ and $((F,G)_1,(H,I)_2)_4$ on the test. Their ordered values are $(-1036800,-8709120)$. On $(x^3,y^3,x^4,y^4)$ they are $(20736,124416)$, and on $(x^3,y^3,(x+y)^4,(x-y)^4)$ they are $(331776,-995328)$. Therefore every planar invariant evaluates on the 24-term cycle as $-190/3$ times the first product plus $5/6$ times the second.

There is a single transverse vector-momentum list. In transverse coordinate order it is
$$
(1,0,0),\ (1,1,0),\ (0,1,0),\ (0,1,1),\ (0,0,1).
$$
Each vector test is its polarization times $\exp(p\cdot y)-1-p\cdot y$. Use $\exp(k\cdot y)-1$ for the marked function. The $89$ rows of the weight table combine polarizations, test covectors, and $k$ equal to the first coordinate covector, the second, or twice the first. In the order of the table the columns are the primitive integer weight, the five-bit polarization selector, the function selector $(0,1,3)$, and the zero-based transverse covector axis. The replay constructs each polarization orthogonal to its momentum. The weighted evaluation kills all $132$ four-vector primitive columns over the integers. At $s=1$, the outgoing five-to-six cochain map has rank $37=22+15$. On its kernel the test has nonzero evaluations $128$ and $199$ modulo $251$ and $257$ respectively, using the computed kernel basis at each prime. Since the degree-five and degree-six rational ranks are attained modulo these primes, the five-cochain kernel lifts locally at each prime. This proves that the integer test detects the unique rational line. All test vectors have zero value and first derivative at the origin, including the marked vector, so the raw integrated graph may be evaluated on the test. [Integer weights](mixed_graph_probe.tsv) · [Boundary and nonzero-pairing certificates](mixed_graph_probe.cpp)

In each planar term put scales $x_1,x_2,x_3$ on the three pure planar momenta, $a$ on the planar momentum of the mark, $b$ on its transverse function momentum and $t_1,\ldots,t_5$ on the transverse vectors. The one-form momentum is minus the sum of the nine vector momenta. In the ten-column graph determinant precisely four columns enter the two planar directions in the two loop rows. Thus either the mark supplies a planar column and the one-form a transverse column, or conversely. After the eight pure columns are factored out the determinant is
$$
\frac{D_\theta}{x_1x_2x_3t_1t_2t_3t_4t_5}
=b\sum_{i=1}^{3}c_i x_i+
 a\sum_{i=1}^{5}c_{3+i}t_i+ab\,c_9 .
$$
Each integer coefficient is a product or sum of two products of a $4\times4$ planar determinant and a $6\times6$ transverse determinant. All terms have planar degree ten after including the Wick tensor. Set $a=1$ and extract $x_1x_2x_3^2$ in the residual polynomial to recover the planar orders $2,2,3,3$. Multiplication by $3!3!4!4!=20736$ converts these exponential coefficients to the Hamiltonian monomials used above. On the factored residual polynomial take the difference of its values at one and zero in each $t_i$ to remove linear transverse jets. Take the same difference in $b$ to remove the constant function. Thus the positive-jet calculation requires $64$ evaluations in the three-variable ring modulo $(x_1^2,x_2^2,x_3^3)$ at each Schwinger point. External odd signs and the ordering sign of the ten determinant columns must both be restored when inserting a test into a graph.

| Graph | Labeled path partitions | Surviving partitions | Tested vertex order | Projected numerator modulo $1009$ |
|---|---:|---:|---|---:|
| $166$ | 252 | 132 | `0146723589` | 554 |
| $256$ | 1260 | 660 | `0146723589` | 985 |
| $346$ | 2520 | 1320 | `0146723589` | 350 |
| $355$ | 3150 | 2244 | `6901472358` | 902 |
| $445$ | 4200 | 3192 | `0691472358` | 412 |

Digits $0,1,2,3$ are the three ordinary planar inputs and the mark, $4,5,6,7,8$ are the five transverse inputs and $9$ is the one-form. On each path of length $n_P$ choose the edge resistance fraction $2j/[n_P(n_P+1)]$ at edge $j$, starting at one. The three conductance fractions are $(1,2,3)/6$. Each entry includes the full integer detector, the planar Hamiltonian factors, the odd external sign, $D_\theta$ and all four Wick orientations at the indicated order. These nonzero residues prove that each shape has a nonzero ordered rational integrand after projection. The graph symmetry weights, sum over all orders and integration over the Schwinger parameters still enter the class coefficient. The determinant replay compares against direct ten-dimensional integer determinants. The positive-jet Wick replay checks its planar coefficient independently by scalar interpolation through degree seven and quartic polarization. [Determinant and path partitions](mixed_graph_kernel.py) · [Projected Wick integrands](mixed_graph_point.cpp)

The resistance integrals and the sum over orders on one path have a common subset recursion. For interior label set $S$ let $Q_S$ be the sum of its momenta. For any ordering let $Q_e$ be the prefix before edge $e$, and let $x_e$ be the ordered edge's resistance fraction. The source uses $\hat m_P=\sum_e x_eQ_e$, with the $x_e$ uniform on the simplex. For an auxiliary vector $u$ paired with these momenta, the generating function of its ordered path moments is
$$
\prod_{e=1}^{n_P}(1-u\cdot Q_e)^{-1}.
$$
The coefficient of $u^\alpha$ is the moment of $\hat m_P^\alpha$ multiplied by $(n_P-1+|\alpha|)!/[\alpha!(n_P-1)!]$. This follows by expanding each geometric series and integrating each simplex monomial.

Hold the initial current $\ell_P$ formal. Let $T_j(\ell_P+Q_{S\setminus\{j\}})$ be the local Wick matrix on appending label $j$, with its incoming and outgoing orientations fixed. Starting with the identity matrix on the empty set, the ordered sum times its moment generator is
$$
M_S=(1-u\cdot Q_S)^{-1}\sum_{j\in S}
  M_{S\setminus\{j\}}\,T_j(\ell_P+Q_{S\setminus\{j\}}).
$$
A prescribed one-form position restricts the final-label sum according to the prefix size. The four Wick signs depend on that position and the shape. Permutations of the vector labels inside each fixed set supply the same combined exterior and determinant sign. The recurrence therefore needs at most $m2^{m-1}$ matrix transitions for $m=n_P-1$ labels. To assemble the paths, retain the polynomial in their three initial currents and substitute $\ell_P=w_P\sum_Q\hat m_Q-\hat m_P$ before applying the three moment functionals. The remaining conductance moments are the Dirichlet moments above. For both one-form orientations with a fixed and an unrestricted position, direct nested ordered-time integration matches all $4\cdot13\cdot25=1300$ integer matrix coordinates through moment degree twelve. [Subset recurrence and independent ordered integration](mixed_graph_paths.py)

### Highest conductance degree of the complete $166$ shape

Write $M=\hat m_B+\hat m_C$ for this shape. The conductance-dependent part of a current on path $P$ is $w_PM$. In each Wick monomial there is one volume multivector, the test covector, nine vector polarizations and thirteen momentum covectors. The volume needs at least four momenta if it contains the test covector and five otherwise. At most one of these can be $M$. Thus the conductance degree is at most ten. At degree ten every vector polarization must pair with $M$ and the test covector must lie in the volume. The same count bounds the degree in a single path mean by ten.

The degree-ten coefficient can therefore be extracted with every vector polarization set equal to one vector $v$ with $M(v)=1$. Project the external momenta and the one-form along $M$ into the annihilator of $v$; the volume already contains $M$. A forward transfer at that degree uses $v\otimes b-w_P\Pi$, where $b$ is the transverse part of the edge current and $\Pi=1-v\otimes M$. A segment of $s$ consecutive forward matrices is $(-w_P)^{s-1}$ times its first matrix: $\Pi v=0$, $b\Pi=b$, $b(v)=0$. A reverse segment keeps its last matrix by transposition. The unprojected currents also contain a component of $b$ along $M$ before the projection; using it in place of $w_PM$ lowers the conductance degree and does not change the stated coefficient.

Place the one-form on path $B$. The other long path is exchanged with it by the graph automorphism. Let $Q$ be the prefix before the one-form, $L$ the total $B$ momentum and $p_\gamma$ the one-form momentum. The surviving endpoint and segment momenta, modulo $M$, lie in the span of $Q,L,p_\gamma,\hat m_B$. The cubic vertex at an endpoint leaves an ordinary differentiated exterior one-form at its bivalent source; with the cubic vertex in the interior its two adjacent momenta wedge to one current wedged with $p_\gamma$. Thus the universal degree-ten coefficient is a linear combination of the three orientations of $\det(\eta,M,p_\gamma,Q,L)$, $\det(\eta,M,p_\gamma,Q,\hat m_B)$ and $\det(\eta,M,p_\gamma,L,\hat m_B)$, multiplied by the nine pairings with $M$. For a one-form at the first or last vertex, respectively $Q=0$ or $Q=L-p_\gamma$, leaving one independent determinant.

In the leading homogeneous term put $w_A=-w_B-w_C$. The eleven independent exterior fixtures across the five positions and eleven values of $w_B/w_C$ suffice by interpolation. Exact integer Wick contractions give zero for all $121$ summed evaluations, with $170$ nonzero individual orientation evaluations. The conductance degree is therefore at most nine. For path degree $d\le10$ the path denominator divides $(5+d)!$. For conductance degree at most nine the odd denominator factors end at $29$. The complete $166$ contribution divided by $\pi$ cannot have a prime $31$ in its denominator. This conclusion includes all orders, the factor $1/4$ from graph automorphisms and the theta normalization $1/158400$. The value at the remaining primes is undetermined. [Exact top-symbol cancellation](mixed_residue_angular.cpp)

### The highest dyadic pole of the symmetry-weighted $166$ shape

Remove the source factor $\pi/158400$ and sum one representative per graph automorphism orbit. For integral external momenta and polarizations the two independent path integrals give this normalized shape a denominator with at most nine factors of two. The following rules out the a priori tenth factor after summing orders.

At degree $d$ in a single path mean the moment of the component monomial of multi-index $\alpha$ is $5!\alpha!/(5+d)!$ times an integral coefficient of the ordered path moment generator. In five components the maximal denominator exponents at $d=0,\ldots,10$ are
$$
(0,1,1,4,4,5,4,6,5,6,5).
$$
The joint degree in the two path means is at most eleven: at least two of the momentum covectors in the volume are independent of their span. The only possible tenth-factor bidegrees are $(3,7),(4,7),(5,5)$ and their transposes. At the smaller selected degree $d=3,4,5$, only squarefree component monomials attain the indicated maximum. All conductance moments are odd units, with residue one. Hence $w_A=w_B=w_C=1$ in the leading residue, and the initial currents on $B,C$ are $\hat m_C,\hat m_B$, respectively.

Keep the larger-degree path ordered for now and sum the smaller-degree path over all orders. Introduce distinct commuting square-zero label variables for its five waves and square-zero variables for the $d$ selected observer components. A pass-through, with entering current $k$ and external transverse wave $(p,u)$, acts on the directed transverse fiber by a first-order Lie derivative. In characteristic two its matrix in the forward direction is $(k\cdot u)1+u p^{\mathsf T}$, and in the dual direction it is $(k\cdot u)1+p u^{\mathsf T}$. For example the difference between the first matrix and the original forward tensor is $u(k+p)^{\mathsf T}$, annihilated by the next directed vertex or by the transverse endpoint. The shifted argument after the vertex is $k+p$. The ordinary generators thus have scalar first-order symbol on both fibers. The one-form source changes fibers by a zeroth-order matrix, and the cubic one-form sink does so by a first-order matrix. The latter contains $\epsilon(\eta,k,p_\gamma)$ in its two remaining slots.

The summed moment is the coefficient of all selected labels in the $(5+d)$-th power of the sum of these generators and the diagonal observer. To see the insertion count, insert each observer in any of the six gaps: its symbol in that gap is the paired prefix momentum. Replacing the prefix by the entering current adds the central scalar obtained by pairing the observer with the initial current. Its square vanishes on the squarefree observer component. In position variables an ordinary generator sum has the form $X^a\partial_a+A$ with scalar principal symbol. Its square is first order: cross derivatives cancel in characteristic two, and $(X^a)^2=0$ in the square-zero label algebra. The square again has scalar principal symbol. The part linear in the distinguished one-form is a commutator with its matrix differential operator; its order stays at most one. Iterating gives order at most one for powers $2,4,8$. The central shift can contribute once, with the same order count. Thus for $5+d=8,9,10$ the resulting transfer has degree at most $1,2,2$ in the opposite path mean. Independent subset evaluation on all four fibers verifies these bounds over the field with sixteen elements. [Marked and ordinary transfer checks](mixed_binary_path.cpp)

The two endpoints add at most two powers of that mean. A cubic endpoint wedges two currents, so it can use the same mean at most once. All transfers on the larger-degree path see the smaller mean as their initial current. The largest degrees available to feed the other path are therefore $3,4,4$, strictly below the required $7,7,5$. Every possible leading bidegree vanishes.

The automorphism quotient can be taken before this unrestricted sum. Put the one-form on $B$ by the long-path exchange. Fix the simultaneous reversal by choosing the lexicographically smaller ordering on the path held fixed. The groups of leading bidegrees with path valuation pairs $(4,6),(6,4),(5,5)$ are stable under reversal. Indeed translation $\hat m_P\mapsto L_P-\hat m_P$ sends a monomial of multi-index $\alpha$ into indices $\beta\le\alpha$. The ratio of the translated moment's factorial prefactor, including its binomial coefficient, to the original prefactor is the integer $(5+|\alpha|)!/[(5+|\beta|)!(\alpha-\beta)!]$. No path can acquire a worse denominator under reversal. This permits choosing the reversal representative separately for each leading valuation group.

For the specified mixed detector both planar wave coefficients, after multiplying by $20736$, are divisible by $2^7$ in the rationals localized at two. Since $158400=2^6\cdot2475$, the complete symmetry-weighted $166$ detector coefficient divided by $\pi$ has at most **eight** factors of two in its denominator.

### Integrated graph-symmetry orbit on the mixed-current detector

For the $166$ order `0146723589` the normalized Schwinger mean, including the full $89$-row detector, both planar tensors, the nine odd external signs and all four Wick orientations, is

$$
\frac{7852756428219172854406518899931430703856}{165629555}.
$$

This order has four distinct graph-symmetry images. Exchanging the long paths sends the first five determinant rows to the last five minus the first five. Its determinant sign is $(-1)^5=-1$ and moving the one-form by five slots has the same sign. Reversing both paths sends the one-form from slot nine to five and preserves both signs. Resistance fractions reverse on each path, $\hat m_P\mapsto L_P-\hat m_P$ for total momentum $L_P$, and each oriented edge current reverses. Signed incoming momenta at the vertex are thereby preserved. Relabeling the internal components of all twelve primitive Wick tensors matches the four orientations individually with their total Wick signs for each automorphism. The two long paths have equal Dirichlet parameters. The four equal integrated terms cancel the order-four graph-symmetry divisor. Consequently this orbit contributes

$$
\frac{14872644750415100103042649431688315727}{49688866500}\,\pi
$$

with the source factor $\pi/158400$. The sign of the sum over the remaining orbits and the value of the complete class are open.

The polynomial contraction uses the uncontracted initial currents on $B,C$ as its ten variables. At an ordinary forward vertex with momentum $p$, polarization $u$, and initial edge current $k$, the original matrix $u k^{\mathsf T}-(k\cdot u)1$ equals $-u p^{\mathsf T}-(k\cdot u)1$ on the incoming transverse fiber of momentum $k+p$. The reverse matrix on the dual fiber is $-p u^{\mathsf T}+(k\cdot u)1$. Every directed chain starts at a cubic or bracket output transverse to its edge. The one-form and endpoint vertices retain their full primitive tensors. The three transverse waves on $B$ select eight polarization masks. For each mask the code first sums the weighted $C$-path matrices including the one-form tests, then multiplies the two path tensors.

For a monomial in the two initial currents substitute $\ell_B=(w_B-1)\hat m_B+w_B\hat m_C$ and $\ell_C=w_C\hat m_B+(w_C-1)\hat m_C$. Integrate each conductance monomial by the Dirichlet formula and the two resistance monomials by the prefix recurrence. All these linear functionals are multiplied by a common integer denominator before contraction. The reduced determinant is linear in the marked transverse scale or in a pure transverse scale. After removing the three pure planar scales, the integrand is homogeneous of planar degree seven. The desired exponents on the four planar scales are $(1,1,2,3)$. The centered differences at the integral nodes $2j-d$ in each direction of order $d$ have combined denominator $2^7\,1!1!2!3!=1536$. Homogeneity pairs the global sign reversals, leaving $24$ nodes. At each node the six transverse differences retain precisely the required positive jets. An independent finite-field extraction on $2\times3\times4$ roots of unity gives the normalized integral $735\pmod{1009}$. The fixed Schwinger-point evaluation from the primitive word algorithm gives $554\pmod{1009}$ with either Taylor extraction.

[Exact mixed-detector integral](mixed_integrated_order.cpp) · [Sixteen primitive symmetry identities](mixed_order_orbit.cpp) · [Independent exact scalar control and full scalar shape](mixed_integration/README.md)

### Integrating zero-momentum insertions in the path-order sum

When a positive-jet difference evaluates a pure wave at zero momentum, its full Lie matrix is the scalar $\pm(\ell_P+Q)(u)$ times the identity, where $Q$ is the prefix momentum at its position. After summing path orders let the five labeled tail times be independent Lebesgue coordinates on $[0,1]^5$, with common prefactor $5!$. At time $t$ the prefix is $\sum_{i:\,t_i>t}p_i$. The zero-momentum insertion appears in no other prefix or path mean. Its conditional time integral is

$$
\int_0^1 \left(\ell_P+\sum_{i:\,t_i>t}p_i\right)(u)\,dt
=(\ell_P+\hat m_P)(u)=w_P M(u).
$$

Here $M=\hat m_B+\hat m_C$. Moving this scalar across the one-form switches its orientation sign and the odd-test sign $(-1)^{9-\mathrm{slot}}$ together. Its sign can thus be read with every removed insertion placed at the last position of its path. For $k$ zero-momentum vectors the shortened normalized path acquires the factor $5!/(5-k)!$ and the product of their $w_PM(u)$ insertions. An independent enumeration with labeled insertions, integer simplex gap lengths and all choices of marked path matches the conditional integral in $368$ cases. [Conditional integration and signs](mixed_order_orbit.cpp)

### Which transverse derivatives of the marked vector detect the class

Restrict the five-vector cochains by the Taylor degree $s$ of the marked transverse function. On degrees $s=1,2$, the four-to-five differential has rank $88$ and adjoining a rational closed representative of the line raises its rank to $89$. On degrees $s=1,3$ the same representative agrees exactly with a boundary. The representative is supported in $s=1$. Hence detection requires comparing the linear and quadratic marked responses. Cubic marked data can be omitted.

For each of $110$ integer inputs the evaluation divided by the function momentum scale is quadratic in that scale. Its values at $1,2,3$ recover the three homogeneous coefficients by two finite differences. The target evaluations have ranks $84,45,9$ on the respective degrees. The closed representative is an integer vector in the kernel of the $65$ outgoing evaluations of rank $37$. On $s=1,3$ an $87$-column rational primitive solves every coordinate exactly. On $s=1,2$ the augmented rank is attained modulo $251$; the full rational differential already has rank $88$.

One single transverse function momentum suffices for a complete detector. The five ordinary transverse momenta in order are $(1,0,0),(1,1,1),(0,1,0),(1,1,-1),(0,0,1)$. The marked function momentum is $(1,0,-1)$. The polarization selector has the same definition as in the original table. The $94$ rows list the integer weight, mask and one-form axis. All $132$ boundary pairings vanish exactly. On the degree-nine cohomology line the original detector divided by the single-momentum detector has the value

$$
-\frac{1029229447699740945723408289}{114415887442730778653770948594609901093939178844554683860158369901482612}.
$$

This normalization applies to the fully summed closed graph. [Single-momentum weights](mixed_single_weights.tsv) · [Jet matrices](mixed_probe_jets.cpp) · [Exact primitive and normalization comparison](mixed_probe_jets.py)

### Separate current degrees and the prime eleven for the $166$ integral

Place the one-form on $B$ and let $B,E$ be that path's initial and final covectors. Write $C,F$ for the unmarked path's covectors and $A=-B-C$ for the short edge, so $F=C+B-E$. At an ordinary wave $(p,u)$ with entering current $k$, put $\alpha=k(u)$. Removing its scalar gives the forward transfer $-g$ and the reverse transfer $g^{-\mathsf T}$, where

$$
g=1+u p^{\mathsf T}/\alpha,\qquad \det g=1,\qquad g^{\mathsf T}k=k+p.
$$

The identities hold as rational functions at $\alpha\ne0$ and the original Wick numerator is polynomial. Let $L,R$ be the ordered transports to the left and right of the one-form, $G=LR$, and $H$ the ordered transport on the unmarked path. Set $J=H^{-\mathsf T}$, $\eta'=L^{-\mathsf T}\eta$ and $p'=G^{-\mathsf T}E-B$. After the nine scalar factors have been removed, the two possible matrices on the marked path are

$$
S=(\eta' p'^{\mathsf T}-p'\eta'^{\mathsf T})G,
\qquad T=\epsilon(\eta',B,p')G^{-\mathsf T}.
$$

In the second matrix $\epsilon(\eta',B,p')^{ij}=\epsilon^{ijklm}\eta'_k B_l p'_m$. Write $T(x,y)=x_iT^{ij}y_j$. Then $T(B,y)=T(x,E)=0$ and $H^{\mathsf T}C=F$. For a cubic on the marked path the two endpoint contractions, before their relative signs $-,+$, are

$$
\begin{aligned}
&-T(H^{\mathsf T}A,F)+T(F,H^{\mathsf T}A)
 +\operatorname{tr}(H)T(A,F)-T(A,H^{\mathsf T}F),\\
&-T(C,JA)+\operatorname{tr}(J)T(C,A)+T(JA,C)-T(JC,A).
\end{aligned}
$$

For the two endpoint cubics the signs are $+,+$. Writing $U^{ijk}(x,y)=\epsilon^{ijklm}x_l y_m$, their two contractions are

$$
\begin{aligned}
&-U^{ijk}(B,C)\,[S_{ji}(JE)_k-E_i(SJ^{\mathsf T})_{jk}],\\
&\phantom{-}U^{ijk}(E,F)\,[S_{ij}(H^{\mathsf T}B)_k-B_i(S^{\mathsf T}H)_{jk}].
\end{aligned}
$$

The local forward signs and exterior one-form position sign give this same four-sign list at all five positions. Direct evaluation agrees separately with all $140$ signed primitive orientations on $35$ conserved inputs. [Transport and original-word comparison](mixed_transport.cpp)

First set $C=tM$ with $B,E$ fixed and put $Q=B-E$. At infinity $H=1+t^{-1}K$ modulo $t^{-2}$, with $K^{\mathsf T}M=Q$ and $\operatorname{tr}K=0$. The endpoint-cubic pair has degree at most one in $t$; its degree-one terms cancel because the remaining two covectors in the volume are $(E,B)$ and $(B,E)$. In the interior-cubic pair the two degree-two coefficients before their relative signs are both $-4T(M,M)$. After inserting $F=tM+Q$ and $A=-tM-B$, the sum of the signed degree-one coefficients is

$$
2\big[T(M,Q)+T(B,M)-T(M,B)\big]=0 .
$$

The remaining transfer gains at most five powers from its five unmarked scalars. Thus the full Wick polynomial has degree at most five in $C$.

Now set $B=tM$ with $C,F$ fixed. Both $G$ and $L$ equal the identity modulo $t^{-1}$. To leading degree the source is the two-form $\eta\wedge p_\gamma$, and the coefficient of $t$ in the sink is the antisymmetric matrix $V=\epsilon(\eta,M,p_\gamma)$. The reduced numerator has degree at most two: each endpoint cubic uses both the varying marked and the fixed unmarked current, while each interior-cubic term has one bracket differentiated on $A$ and the other on $C$ or $F$. At degree two the four contractions in order are

$$
2 V(JM,C),\qquad 2 V(H^{\mathsf T}M,F),\qquad
-2 V(H^{\mathsf T}M,F),\qquad -2 V(JM,C).
$$

Their sum vanishes. The marked path has four ordinary scalar factors, giving degree at most five in $B$ as well.

In each orientation the short edge supplies $w_A$: an endpoint cubic contains $B\wedge C$ or $E\wedge F$, and momentum conservation writes each as a wedge with $A=-B-C=-E-F$. With an interior cubic, one of the bracket endpoints differentiates its short edge. Since $A=w_A(\hat m_B+\hat m_C)$, absorb $w_A$ into the conductance measure. Its parameters change from $(11/2,1/2,1/2)$ to $(13/2,1/2,1/2)$ and the scalar prefactor becomes $\pi/187200$.

Put $u=1-w_B$ and $z=w_C/u$. Under the shifted measure these are independent beta variables with parameters $(7,1/2)$ and $(1/2,13/2)$ respectively. The initial currents are $B=-u\hat m_B+(1-u)\hat m_C$ and $C=uz\hat m_B+(uz-1)\hat m_C$. A term of degree $d$ in $\hat m_B$ therefore contains $u^d$ before division by $w_A=u(1-z)$. The reduced polynomial has $z$-degree at most four by the degree-five result and $u$-degree at most eight by the conductance-degree-nine result. For $k=5,6,7,8$, the moment of $u^k$ is $(7)_k/(15/2)_k$, divisible by eleven in the rationals localized at eleven; $(a)_k$ denotes the product of $a,a+1,\ldots,a+k-1$. The four possible positive $z$ moments have denominators prime to eleven. Hence each term with $d\ge6$ acquires a factor eleven from the conductance integral. Interchanging $B,C$ gives the same statement for the other resistance integral. A resistance pole at eleven requires $d\ge6$, has order at most one, and cannot occur on both paths since their total degree is at most eleven. **Every ordered $166$ coefficient, including its scalar kernel and after removing $\pi$, has denominator prime to eleven** for integral external data.

### Sharp degree of the $166$ Wick polynomial

The separate degrees in $B,C$ are at most four. Both bounds are attained in the independent literal-word checks. In particular all normalized resistance and conductance integrations may use the joint degree-eight rule.

Let $W=H^{\mathsf T}$, $N=W-1$, $R=W^{-1}-1$, $Q=F-C=B-E$. Use $T(B,\cdot)=T(\cdot,E)=0$ and $F=WC$ in the two signed interior-cubic contractions above. Their exact sum reduces to
$$
T(C,(N+R)E)-T((N+R)B,F)+T(Q,NB)+T(RE,Q)
+(\operatorname{tr}W-\operatorname{tr}W^{-1})T(C,F).
$$
For $C=tM$ both $N,R$ have a factor $t^{-1}$ in their expansions at infinity. The identity $N+R=N^2W^{-1}$ gives a factor $t^{-2}$. Since $\det W=1$, the degree-one trace term in $N$ is zero and $2\operatorname{tr}N-\operatorname{tr}(N^2)$ has a factor $t^{-3}$, by expanding $\det(1+N)=1$ through second order. Thus $\operatorname{tr}W-\operatorname{tr}W^{-1}$ has a factor $t^{-3}$. The entire interior-cubic sum has a factor $t^{-1}$.

Write $D(S,X,Y,Z)=\epsilon^{ijkab}S_{ij}X_kY_aZ_b$. The endpoint-cubic sum equals
$$
D(S,W^{-1}E,B,C)+D(SW^{-\mathsf T},E,B,C)
+D(S,WB,E,F)+D(WS,B,E,F).
$$
Its terms with $W$ replaced by the identity add to $2D(S,B,E,Q)=0$. Modulo terms with a factor $t^{-1}$, the remaining sum is
$$
-D(S,NE,B,C)-D(S,E,NB,C)-D(NS+SN^{\mathsf T},E,B,C).
$$
The infinitesimal determinant identity rewrites this expression as
$D(S,E,B,NC)-(\operatorname{tr}N)D(S,E,B,C)$. Here $NC=Q=B-E$, and $\operatorname{tr}N$ has a factor $t^{-2}$. The endpoint-cubic sum also has a factor $t^{-1}$. The five ordinary scalar factors give unmarked degree at most four.

For $B=tM$ the finite transport $W$ is arbitrary. Its diagonalizable elements are dense, so choose a determinant-one eigenbasis and let the five eigenvalues be free Laurent variables subject to product one. Write $P=C-F$, so $E=B+P$. If $G^{\mathsf T}=1+t^{-1}K$ modulo $t^{-2}$, then $\operatorname{tr}K=0$ and the transported one-form momentum at infinity is $p=P-KM$. The leading degree-two coefficient above vanishes identically for arbitrary one-form test and momentum, so their subleading changes do not enter the degree-one coefficient. Insert
$$
S=(\eta p^{\mathsf T}-p\eta^{\mathsf T})(1+t^{-1}K^{\mathsf T}),\qquad
T(x,y)=t\,\epsilon(x,y,\eta,M,p)-\epsilon(x,Ky,\eta,M,p)
$$
to the displayed precision. Expand the exact endpoint-cubic and interior-cubic sums in the free entries of $C,M,\eta,K$ and the four independent eigenvalues. Their degree-one coefficients each have $14280$ nonzero integer Laurent monomials. The coefficients cancel individually. The reduced sum has degree at most zero in $t$, and the four marked-path scalars supply at most four powers. [Exact Laurent-polynomial identity](mixed_transport_symbol.py) · [Literal finite-field degree checks](mixed_transport.cpp)

The conductance degree is now at most eight, making every $166$ coefficient divided by $\pi$ integral at $29$. At total path degree at most eight the only bidegrees that could produce nine factors of two in the normalized resistance denominator are $(3,5)$ and $(5,3)$. The degree-three path moment must use squarefree components. Its eight-letter path power in characteristic two contributes at most one power of the opposite mean and the endpoints at most two; that falls short of five. The reversal-orbit argument above applies to the two distinct valuation groups $(4,5),(5,4)$. The normalized symmetry-weighted average has at most eight factors of two in its denominator. Including $20736$ times the planar coefficients and $1/158400$ leaves at most **seven** factors of two in the coefficient of $\pi$ for the original mixed detector.

### Sharp current and conductance degrees for all five order-$g$ theta shapes

The four displayed endpoint contractions apply to every theta shape. Cyclically label the marked path $B$ and write $G_A,G_C$ for the ordinary transports on the other paths. A change of frame at the right endpoint by $G_A^{-1}$ turns the first path transport into the identity. In the previous formulas replace $G$ by $LR G_A^{-1}$ and $H$ by $G_C G_A^{-1}$, and transport all right endpoint covectors by $G_A^{-\mathsf T}$. Thus $F=H^{\mathsf T}C$ and $E=B+C-F$ in the common frame. The primitive Wick sign, external sign and ordering of the two bracket inputs give the same four-orientation signs $(+,+,-,+)$ times a shape sign.

| Shape | Shape sign | Primitive orientations checked | Observed joint degree | Observed collinear degree |
|---|---:|---:|---:|---:|
| $166$ | $+1$ | 200 | 8 | 7 |
| $256$ | $-1$ | 200 | 8 | 7 |
| $346$ | $+1$ | 200 | 8 | 7 |
| $355$ | $-1$ | 200 | 8 | 7 |
| $445$ | $+1$ | 200 | 8 | 7 |

The sign and each contraction are checked separately at all ten marked positions, including the shape with one empty interior path. [Common finite transport and original Wick words](mixed_transport.cpp)

For the upper bounds write $W=H^{\mathsf T}$, $U=G^{\mathsf T}$ and set $N=W-1$, $R=W^{-1}-1$, $Z=N+R$, $Q=F-C$. The exact sum of the endpoint-cubic pair, using the infinitesimal determinant identity in the preceding section, is
$$
(\operatorname{tr}N)D(S,Q,B,C)+D(S,ZE,B,C)-D(SZ^{\mathsf T},Q,B,C)+D(S,NB,B,Q).
$$
The exact interior-cubic pair is
$$
T(C,ZE)-T(ZB,F)+T(Q,NB)+T(RE,Q)
 +(\operatorname{tr} W-\operatorname{tr}W^{-1})T(C,F).
$$
Each term contains at least two powers of $W-1$ near the identity, using $Q=NC$, $Z=N^2W^{-1}$, and $\det W=1$.

Scale both independent initial currents by $t$ and divide the initial and final endpoint currents by $t$ in these two identities. Their new values are denoted by $B,C$, $F=WC$, $E=B+C-F$. The marked covector momentum divided by $t$ becomes $p'=U^{-1}E-B$. The change of scale leaves an overall $t^4$ multiplying the two pairs. Including their nine ordinary scalar factors gives total degree at most thirteen before the following cancellations. Both relative transports approach the identity in the formal variable $r=1/t$. Their local coordinates can be taken to be the logarithms. It suffices to compute the homogeneous Taylor coefficients along $W=\exp(rL), U=\exp(rK)$ with arbitrary traceless matrices $L,K$: these straight lines in logarithm coordinates determine each homogeneous polynomial there.

The $r^0,r^1,r^2$ coefficients vanish because $p'$ contributes a further factor of $r$. The matrix $L$ can be diagonalized on a dense set. In that basis a diagonal basis change sets all components of a generic $B$ to one, preserving diagonal $L$ and equivariance up to the common volume weight. Linearity and coordinate permutations reduce the independent marked test covector to one coordinate covector. Keep the other initial covector $C$ arbitrary. The exact cubic and quartic Taylor coefficients then have $414$ and $3060$ nonzero monomials in each pair, respectively. Every signed sum vanishes. This proves divisibility by $r^5$, and total initial-current degree at most $13-5=8$.

For collinear initial directions set $C$ equal to an arbitrary scalar multiple of $B$ in the same calculation. At order $r^5$ both endpoint pairs vanish individually. Therefore the collinear current degree is at most $13-6=7$. Bounded offsets in the initial currents preserve this conclusion: the degree-eight homogeneous polynomial vanishes on collinear leading currents. Each conductance variation varies the currents in multiples of the single covector $\sum_P\hat m_P$, and scaling any single path mean likewise gives collinear leading currents. The bounds thus apply to every individual resistance mean and to the entire conductance polynomial. All five pairs of upper bounds are attained by the independent primitive evaluations above. [Formal Taylor identities in arbitrary traceless matrices](mixed_theta_degree.py)

### Lifting the mixed transverse class through the full restricted differential

On the one-mark block algebra let $C$ be a planar current four-cochain with one marked planar Hamiltonian, relative to the block linear algebra. Choose $C$ closed and put $P=\operatorname{div}_x C$. The exact construction in [mixed_current_lift.py](mixed_current_lift.py) imposes all incoming chain boundaries and the planar raising operator. In 224 weight-zero coordinates these 289 equations have joint rank 215. The resulting primitive integer $C$ has 114 nonzero coordinates. Its divergence has 148 and pairs with the integer planar four-cycle as

$$
P(\text{planar four-cycle})=2365440.
$$

The replay verifies closure and invariance separately for $C$ and $P$. The planar dual derivative omits derivatives of ordinary quadratics that would leave the relative complex. The marked planar slot includes the affine jets as in the transverse projection.

Take the closed transverse five-cochain $T$ representing the alternating line above, with one marked function $f$ vanishing at the origin. Use its representative supported on the first derivative of $f$. Its transverse divergence $\operatorname{div}_yT$ is a scalar cochain with the same five vector arguments. It is closed: divergence intertwines the vector-field action with the scalar action because for every divergence-free transverse vector $Y$ one has $\operatorname{div}[Y,T]=Y(\operatorname{div}T)$. This identity also commutes with the action on the marked argument. In the divergence, the coefficient of the second derivative of $f$ vanishes by antisymmetry of its two indices.

On the mixed field $f(y)Z(x)$ evaluate the planar cochain on $Z$ and the transverse cochain on $f$. Put all planar slots before the transverse slots. Choosing zero cochains on unmarked inputs, the following pair is closed on the one-mark restriction:

$$
J_{\rm transverse}=P\otimes T,\qquad
J_{\rm planar}=-C\otimes\operatorname{div}_yT,\qquad K=0 .
$$

Indeed both tensor products are closed for the diagonal differential. The marked action on an unmarked output vanishes with that output. Differentiating on the product of the two coordinate blocks gives

$$
\operatorname{div}_x J_{\rm planar}
+\operatorname{div}_y J_{\rm transverse}
=-(\operatorname{div}_xC)\otimes\operatorname{div}_yT
+P\otimes\operatorname{div}_yT=0 .
$$

Its transverse projection pairs nontrivially with the product cycle and the transverse projection is a cochain map. Thus **the unique projected transverse line lifts to a nonzero class in the full one-mark restricted current–scalar complex**. The value of the sum of the five Feynman graphs on that line and extension from the one-mark algebra to the entire five-dimensional formal vector-field algebra are separate calculations.

The planar current in this lift is necessary even when a scalar cochain is allowed. On five transverse vector inputs, the current cochain dimensions indexed by the degree of the marked function are $84,45,9$. A scalar primitive has four transverse vector inputs. For marked degree $s=1,2,3,4,5$ its multiplicity is that of highest weight $(s-3,-3,-3)$ in four positive vector jets of total Euler degree $9-s$. Weyl alternation gives $31,36,22,10,1$, respectively. These bases are realized by products of three volume tensors and contractions of all vector polarizations.

For each of 172 transverse wave tests, evaluate the current differential on six vectors. On five vectors evaluate the divergence by setting the test covector equal to the sum of the five vector momenta and the marked function momentum. On the same five vectors evaluate the scalar differential of each primitive. All vector Taylor orders are at least two, and the marked function vanishes at the origin. The integer detector on the 138 current coordinates is an exact rational combination of 110 closure and divergence equations. With all 100 scalar primitive coordinates adjoined, the detector extended by zero is an exact rational combination of 123 equations of the form $dT=0$ and $\operatorname{div}_y T+dK=0$. The replay solves on a nonsingular minor and verifies the resulting integer identity on every one of the 238 columns.

Thus transverse divergence induces a nonzero map from the unique transverse class into the scalar cohomology. Pairing the planar slots of a full cochain with the planar four-cycle kills the planar differential on the scalar term. Consequently every representative detecting this line requires a planar current; a transverse current and a scalar alone have zero pairing. [Integer wave evaluations](mixed_divergence_probe.cpp) · [Exact row identities and scalar multiplicities](mixed_divergence_test.py)

At arities two through eight the marked transverse scalar cochain dimensions are $(3,30,100,150,105,31,3)$. These follow from the same Weyl calculation at each arity. The six successive differentials have rank at least $(3,27,73,76,28,3)$ over the rationals from nonzero minors modulo 251. At arity five the divergence class proved above is nonzero, so the incoming and outgoing ranks sum to at most 149. At the other arities the square-zero identities and the dimensions force equality in all six lower bounds. The scalar cohomology has dimension one in arities five and six and vanishes elsewhere. In particular transverse divergence is an isomorphism from the transverse-current five-class onto the scalar five-class. [Middle ranks](mixed_divergence_test.py) · [Outer ranks](mixed_scalar_tail.cpp)

Consider the one-mark subcomplex with zero unmarked component and at least one transverse derivative on the marked function, relative to the block linear algebra. Its current and scalar diagonals split into the indicated planar and transverse tensor factors. In every total degree the planar divergence from planar current arity four onto planar scalar arity four is an isomorphism tensor the identity on transverse scalar cohomology. This makes divergence onto scalar cohomology surjective. Its kernel in total degree nine has two lines. The first is planar current arity three tensor transverse scalar arity six. Its planar divergence goes to planar scalar arity three, whose cohomology vanishes. The second is the diagonal cancellation of planar current arity four tensor transverse scalar arity five against planar scalar arity four tensor transverse current arity five. Each of the latter divergence maps is an isomorphism onto the same scalar line. The exact sequence of the current–scalar complex therefore gives **dimension two** for degree nine on this specified subcomplex. The graph detector calculated here measures the second line.

The scalar-divergence test admits either the first and second or the first and third marked transverse Taylor coefficients. At each choice the 94-row detector is an exact combination of 123 closure and selected scalar equations. There is also an exact current and scalar primitive for which the first and fourth scalar equations vanish while the detector is nonzero. The scalar target bases were selected by independent evaluations on these same waves, so zero on these evaluations enforces the corresponding homogeneous scalar identity. The current closure rows attain rank $138-88-1=49$, proving closure of the witness. [Four-grade replay via the optional `--graded` flag](mixed_divergence_probe.cpp)

### Planar current with six transverse vector fields

At planar scalar Euler degree four, the relative chain dimensions at arities two through five are $(1,5,4,1)$ and the outgoing boundary ranks are $(0,1,3,1)$. Choose the closed scalar three-cochain from the 36 monomial coordinates in the replay. For scalar coefficient $F$ its Hamiltonian vector field has components $(-\partial_yF,\partial_xF)$, in the convention $X_F(G)=(F,G)_1$. This current has 120 nonzero monomial coordinates. Its divergence vanishes identically and its planar three-class is nonzero.

In fact there are integral three-cycles on the planar Hamiltonian degree profiles $(3,6,2)$ and $(3,5,3)$, with the mark last. The cycles have 24 and 30 terms and their pairings with the Hamiltonian current are $-2520$ and $2520$. On the first profile the sextic is the top representation in the tensor product of the cubic, the marked quadratic and the test covector. Its multiplicity is one. The same current on the fixture consisting of $x^3,y^6,x^2$ and test covector $dx$ is $360$. Thus every planar invariant evaluates on the first cycle as $-7$ times this single fixture. [Integral cycles and closed Hamiltonian current](mixed_planar_stream.py)

Tensoring this planar Hamiltonian current with the closed transverse scalar six-class gives a strictly divergence-free representative of the first degree-nine line. Both restricted lines in fact admit a closed antisymmetric bivector potential. For this first line take the planar area bivector times the planar stream scalar times the transverse scalar cocycle, with the sign of the Hamiltonian convention. For the second line take mixed components $C^iT^j$ with the planar slot first and extend by antisymmetry. Differentiation in the planar slot gives $P T$, while differentiation in the transverse slot gives $-C\operatorname{div}_y T$. Thus both restricted couplings can be written using the exterior derivative of the external one-form after integration by parts.

For the first line use six positive transverse waves with momenta $(1,0,0),(1,1,1),(0,1,0),(1,1,-1),(0,0,1),(1,0,1)$. The marked function is an exponential minus its constant, with momentum $(1,0,-1)$ or $(0,1,1)$. The 77 weighted polarizations cancel all 150 incoming scalar differential columns exactly. In the 105-dimensional transverse scalar six-cochain basis, the outgoing differential has rank 28. Appending the evaluation raises that rank to 29 modulo 251. Hence the integer test detects the scalar six-class over the rationals. [Wave evaluations](mixed_six_probe.cpp) · [Integer weights](mixed_six_weights.tsv) · [Rank and boundary checks](mixed_six_detector.py)

On this fixture the four Schwinger determinant columns with planar components must use the two pure Hamiltonians, the marked Hamiltonian and the one-form. Write the three planar wave momenta as $s\,dx,t\,dy,a\,dx$. Only $-t\,dy$ contributes from the one-form column, so the determinant has exactly the factor $s t^2 a$. The diagonal planar weights force the full momentum multidegree to have exactly three $dx$ factors and five $dy$ factors, in addition to the test one-form. Setting the marked planar scale to zero in the Wick tensors thus extracts the required term $s^2 t^5 a$ at $s=t=1$. The normalization for the integer cycle is $-7\cdot 3!6!2!=-60480$. Horizontality of the sum over the five shapes removes every linear ordinary vector, including the constant-function term where the marked planar vector itself is linear. The summed detector therefore uses exactly two Wick geometries. The separate shape contributions in this reduced calculation acquire their class interpretation after all five are added. [Two-geometry integration](mixed_integration/six_vector_nodes.py)

The planar current three-class has sharp order three in the total external momentum on this restricted subcomplex. Take the two independent planar momenta as the coordinate covectors and let $(a,b)$ be the sum of all three planar momenta. This fixes the marked planar momentum to $(a-1,b-1)$. On these Hamiltonian exponentials the stream scalar above is exactly
$$
-\frac{a^2b^2-a^2b+a^2-ab^2-ab+b^2}{24}.
$$
The two current components are $-b$ and $a$ times this expression. A change of planar frame covers every independent pair, and polynomiality extends the vanishing through second order to dependent pairs. The two planar current coboundaries have rank two on the coefficients through cubic order; appending the displayed current raises the rank to three. Thus no planar primitive removes its cubic coefficient. [Exact polynomial and ranks](mixed_planar_soft.py)

For the transverse scalar six-class one can also impose exact momentum conservation. Set the marked function momentum to minus the sum of the six vector momenta at each of 34 wave fixtures. The resulting integer test kills all 150 incoming columns, and raises the rank of the outgoing equations from 28 to 29 modulo 251. Pairing with this transverse cycle kills transverse coboundaries on the planar three-vector component, leaving just the two planar boundaries in the preceding calculation. The representative constructed above therefore has the maximum possible order three even among the representatives of its full restricted class with zero unmarked component. [Balanced replay via the `--balanced` option on both the generator and checker](mixed_six_detector.py)

### Order of the mixed four-plus-five current at zero one-form momentum

In the four-planar, five-transverse sector take the closed planar current $C$ and transverse current $T$ in the factorization above. Let $p_x$ and $p_y$ be the sums of the external vector momenta in their respective coordinate blocks, so the one-form carries momentum $-(p_x,p_y)$. Both factors can be chosen to vanish on their own balanced momenta.

For $C$, substitute the marked planar momentum $p_x$ minus the three ordinary planar momenta. Expand the marked Hamiltonian exponential and alternate the other three. The 86 constant-coefficient equations in $p_x$ raise the rank of the 289 closure and invariance equations from 215 to 219 on 224 coordinates. A primitive integer solution has 116 nonzero coordinates. Its planar divergence has 142 and pairs to $73920$ with the same integer four-cycle. Both first partial derivatives of $C$ have 48 nonzero alternating monomial coefficients. Adding all first-derivative equations gives rank 222 and makes the divergence pairing zero. On the three ordinary planar momenta $(1,2),(2,1),(3,-2)$ the coefficient of $(p_x)_1^2$ in $\operatorname{div}_x C$ is $70312/9$. [Exact balanced planar polynomial](mixed_balanced_four.py)

For $T$, substitute the marked function momentum $p_y$ minus the sum of the five transverse vector momenta. At $p_y=0$ the result belongs to the unmarked five-vector current at Euler degree eight. Its multiplicity at highest weight $(-2,-3,-3)$ is
$$
19360-17582-15717+13090+11150-10226=75.
$$
Every marked term has positive function degree. Thus its image in these 75 coordinates vanishes when the five vector momenta separately sum to zero. Evaluation on those further balanced waves has rank at least 16. Evaluation of the original marked restriction has rank at least 59. The two maps compose to zero and $16+59=75$, proving that 59 independent rows impose the complete balanced restriction. The modular minors use integer momenta and transverse polarizations. Together with the 49 independent closure equations the 169 integer tests of each condition admit a primitive current on 106 of the 138 coordinates with nonzero 94-row pairing. Imposing also the first normal derivatives makes the detector an exact rational combination of 132 constraints. Thus the chosen $T$ has a simple transverse zero. [Integer wave tests](mixed_balanced_five.cpp) · [Rational witness and first-derivative obstruction](mixed_balanced_five.py)

The representative
$$
J_{\rm transverse}=(\operatorname{div}_x C)\otimes T,\qquad
J_{\rm planar}=-C\otimes\operatorname{div}_y T
$$
has at least two planar and one transverse total-momentum factors in its transverse component, and at least one planar and two transverse factors in its planar component. The displayed planar coefficient and the simple transverse zero give an actual term of total degree three. The marked transverse factor contains only positive jets of $f$, so the whole representative vanishes identically at $f=1$. The supported planar Euler orders of the mark in $C$ and its divergence are $0,1,2$. A horizontal linear shear $y^j\partial/\partial x^i$ has constant marked planar vector, of Euler order $-1$, so insertion also vanishes on every such shear. Its cubic zero is compatible with the double zero of the graph and uses a zero unmarked block. The following product test also bounds the order for general representatives.

On the planar scalar four-class take three independent wave slots with planar momenta $(1,0),(0,1),(s,t)$, and set the marked momentum to $(-1-s+a,-1-t+b)$. Use the five pairs $(s,t)=(1,1),(1,-1),(2,-1),(-2,-1),(1,2)$. In the coefficients through first order in $(a,b)$ the planar boundary matrix and its augmentation by $P=\operatorname{div}_x C$ both have rank four. Through second order they have ranks five and six. Therefore a rational combination of these Taylor evaluations kills every planar scalar boundary and evaluates nontrivially on $P$. [Exact planar four-scalar calculation](mixed_planar_soft.py)

The transverse current test above supplies a boundary-annihilating evaluation of order at most one on $T$. Its product with the planar test is a cycle, since each factor kills its incoming boundaries. On the zero-unmarked differentiated-function subcomplex the projection to the transverse current is a cochain map. The product therefore pairs nontrivially with every representative of the four-plus-five line, including representatives with unfactorized coefficients, using at most three derivatives in the total one-form momentum. The three-plus-six line has zero transverse projection. Together with its independent balanced test above, this shows that a nonzero degree-nine restricted class can attain soft order three and cannot attain order four.

### Extension of the planar three-class to triangular fields

Consider origin-preserving divergence-free fields on the fixed splitting with planar component $X(x,y)$ and transverse component $Y(y)$, with each component divergence-free in its own variables. The horizontal fields whose coefficients lie in the square of the transverse maximal ideal form a Lie ideal: $Y(0)=0$ preserves this square and the horizontal bracket differentiates only planar variables. In the quotient the horizontal coefficients are affine in the transverse variables and the bracket of two transverse-linear horizontal fields vanishes.

The scalar six-class admits a representative using only the first derivative of its marked function. The first-derivative sector has dimension 70, with closure rank 22 on the wave tests. A rational closed cochain supported on those 70 coordinates pairs nontrivially with the six-vector test. The integer replay verifies every closure equation on the full basis; those equations attain the complete outgoing rank 28. [First-derivative replay with `--linear` on both programs](mixed_six_detector.py)

In the horizontal first-jet quotient take the pure planar bivector given by the planar area tensor, times the stream scalar three-cochain, times this transverse scalar six-cochain. Give it zero value on two or more transverse-linear horizontal arguments and on a purely separated input. On a single transverse-linear argument closure is the known tensor-product identity. On two such arguments their bracket vanishes in the quotient and their derivative at the origin annihilates the planar output plane. On three or more such arguments every term of the differential vanishes. Insertion of a horizontal shear $y^j\partial_{x^i}$ vanishes because the stream has no constant planar marked vector in its support. Cartan's identity gives invariance under these shears. Pullback therefore gives a closed bivector on the full triangular origin stabilizer, relative to its two block linear algebras and horizontal shears. Its divergence is a planar current with zero scalar component. If the pulled-back current were exact, a primitive restricted to zero marks would be closed. The block homotopy above makes that primitive exact on zero marks; extend its preimage by zero on mixed arguments relative to the block linear algebra and subtract its differential. The resulting primitive on zero and one marks has zero unmarked component and contradicts the nonzero three-plus-six pairing. This also proves nontriviality after forgetting the shear-relative condition. Planar differentiation preserves the transverse square ideal throughout.

### Two horizontal marks in the mixed current

In the horizontal first-jet quotient consider the four-plus-five transverse current $P T$ with $T^\alpha(f)=S^{\alpha\beta}\partial_\beta f$ in its alternating first-derivative representative. Keep the ordinary transverse inputs fixed in this tensor. A planar correction on two marked fields $fZ,gW$ with transverse coefficient $S^{\alpha\beta}(\partial_\alpha f)(\partial_\beta g)$ is symmetric in the two marked planar Hamiltonians. It still alternates in the two ordinary planar Hamiltonians. Impose zero insertion of every constant planar marked vector, which is a horizontal shear at this transverse degree.

On three ordinary planar fields and two marked ones, the ordinary bracket preserves each mark, while the two marked fields commute in the quotient. The action of a marked shear contributes its planar vector times the scalar $P$ on the other mark. With $P$ allowed to vary over the invariant closed one-mark scalar four-cochains, the exact system for this planar correction has 602 coordinates and rank 600. A solution has 148 nonzero scalar coordinates and 222 in each planar component. Its pairing on the original planar four-cycle is $-665280$.

Requiring this correction to cancel divergence on two marks adds the terms $P(U,V,W,Z_{\rm marked})+P(U,V,Z,W_{\rm marked})$ on the two ordinary planar inputs $U,V$. They come from transverse differentiation of the formerly unmarked factor and have the same sign because $S$ is antisymmetric. The joint matrix now has rank 601, and the detector vanishes on its entire kernel over the rationals. Thus a planar two-mark current alone, carrying this alternating transverse coefficient, cannot give a strictly conserved completion of the four-plus-five ansatz. Other two-mark tensor components and further transverse jets enter the extension problem. [Rational shear and divergence equations](mixed_two_horizontal.py)

### Number of transverse momentum directions in a four-plus-five detector

Suppose the five ordinary transverse wave momenta lie on at most three complex lines. For three independent lines choose them as the coordinate axes. Up to permutation their multiplicities are $2+2+1$ or $3+1+1$. Give every wave its own independent scale and keep all three components of the marked function momentum free. A multilinear current is determined by the two coordinate polarizations perpendicular to each wave and the three choices of test covector. In each of these 96 evaluations expand the alternating basis products into monomials in the eight free variables.

For $3+1+1$ every coefficient of all 138 current symbols vanishes. For $2+2+1$ the full coefficient matrix has rational rank 18. A 138-dimensional evaluation minor of the unrestrained current basis is invertible modulo 251. Inverting it on the four-to-five differential and restricting to the axial polynomial coefficients gives rank 18 modulo 251. Hence the 132 primitive columns span the full rational image of the restriction. The same primitive works after arbitrary changes of transverse frame and permutations of slots by equivariance. Dependent triples and coinciding lines follow by polynomial specialization, keeping each polarization transverse as its line varies. Thus every boundary-annihilating detector for the transverse class must use at least one configuration with four distinct ordinary transverse momentum lines. [Symbolic wave coefficients](mixed_transverse_directions.cpp) · [Rational dimension and boundary minor](mixed_transverse_directions.py)

Four ordinary lines attain the lower bound. Choose directions $(e_1,e_1,e_2,e_3,e_1+e_2+e_3)$ in the transverse momentum frame, with independent scales. The five scale and marked-momentum tuples in the replay give 408 polarization evaluations. On them the incoming boundary matrix has 86 independent columns modulo 251. A primitive integer linear combination of 87 evaluations kills every one of the 132 incoming columns exactly. Appending its functional to the 49 independent closure equations raises the rank to 50 modulo 251. Thus exactly four distinct ordinary transverse directions suffice for a rational detector of this line.

For the six-vector two-geometry expression the complete $166$ contribution, divided by $\pi$, is
$$
\frac{25824931468885850344108541896448563183380640499855602823157604309082276610559548863876486659268483}{152108775}.
$$
The 28 partitions at both marked momenta were evaluated over 20 prime fields. The reconstruction modulus has 451 bits and the numerator bound at common denominator $293772917487329280000$ has 442 bits. A separate modulo-1009 calculation gives $472$. This is one term in the Ward-reduced five-shape sum. [Exact reconstruction](mixed_integration/six_vector_166_exact.py) · [Prime residues](mixed_integration/six_vector_shapes/166_crt.tsv)

### Minimum number of ordinary transverse momentum lines for the scalar six-class

Suppose three ordinary momenta lie on one line. The three volume determinants use at most three covectors from that line, and only the other three polarizations can contract a covector on it. Each of the three positive vector slots must therefore have degree exactly two. Their dependence on the three wave scales factors as their squares. What remains is alternating in three polarizations from a two-dimensional plane, hence zero. Four or more on one line already violate the degree count.

For at most four ordinary lines the surviving multiplicities are therefore $2+2+2$ and $2+2+1+1$. For the first use the three coordinate axes, and for the second the axes and their sum, with the first two axes repeated. Expand in the six independent scales and three marked momentum coordinates for every transverse polarization. The rational coefficient matrices on the 105 scalar symbols have ranks $8,50$. The 150 primitive columns restricted to the same coefficients have ranks at least $8,50$ modulo 251, using an invertible 105-row evaluation minor to recover the differential modulo this prime. Hence a single primitive agrees with a given six-cochain on the full four-line chart, including all scales and polarizations. Equivariance extends it to every projective frame. Dependent frames and the three-line pattern follow by polynomial specialization with the polarizations kept perpendicular. [Symbolic restrictions](mixed_six_directions.cpp) · [Rank replay](mixed_six_directions.py)

Five lines suffice at one balanced geometry. In the transverse coordinate frame use ordinary momenta
$$
(-1,0,0),\ (1,0,0),\ (0,-1,0),\ (0,0,-1),\ (1,1,1),\ (-1,1,-2),
$$
and marked momentum $(0,-1,2)$. The 64 polarization evaluations admit an exact combination on 34 masks annihilating all 150 incoming columns, with a nonzero pairing on a closed scalar six-cochain. This establishes five as the minimum number of distinct ordinary transverse momentum lines in any detector. [Direct differential](mixed_six_probe.cpp) · [Exact integer test](mixed_six_conserving.py)

For integration use the distinct ordinary momenta
$$
(1,0,0),\ (0,1,0),\ (0,0,-1),\ (0,1,1),\ (-1,-1,0),\ (-1,0,-1).
$$
The marked momentum is their negative sum, $(1,-1,1)$. Sixty polarizations carry integer weights of absolute value at most 16. Their combination annihilates all 150 incoming columns. Its pairing on the entire 77-dimensional cocycle kernel, multiplied by
$$
-2217817151632288990651635154940307443560736538601897342476896593773383364911585160693312,
$$
equals the original 77-polarization detector. The two cycles therefore give the same calibrated evaluation on the closed sum of the five graph shapes. This equality need only hold on the sum. The one-form carries purely planar momentum in every Wick geometry. The numbers of supported path partitions are $28,280,560,524,768$. On the first two shapes the small-cycle coefficients after dividing by $\pi$ are $822,669$ modulo 1009. For $166$, factoring the original integer planar multiplier $-60480$ out of the universal denominator bound gives denominator $4857356439936000$ and a 125-bit common-numerator bound. [Wick geometries](mixed_integration/six_conserving_nodes.py) · [Reconstruction with `--conserving`](mixed_integration/six_vector_crt.py)

The marked Hamiltonian in this planar cycle is quadratic. Its vector vanishes at the planar origin. Since the marked function also vanishes at the transverse origin, the mixed field and its first derivative at the origin both vanish. Hence the mixed action on an unmarked output contributes zero to this pairing, even before subtracting the block-diagonal primitive. The planar cycle kills the planar differential. Pairing any closed full current with it therefore gives a closed transverse scalar six-cochain, to which the calibration on all 77 cocycles applies.

The six large-prime residues for the $166$ small cycle have a 136-bit product. With the 125-bit bound they give the exact contribution
$$
\frac{3568529267}{12168702}\,\pi .
$$
Its coefficient reduces to $822$ modulo 1009 in the independent full partition sum. [Six prime residues](mixed_integration/six_conserving_shapes/166_crt.tsv)

### Nonzero integrated order-$g$ obstruction on the three-plus-six cycle

On the balanced 60-polarization cycle the complete contributions, with the common factor $\pi$ removed, are

| Theta graph | Exact rational coefficient | Modulo $1009$ |
|---|---:|---:|
| $166$ | $3568529267/12168702$ | $822$ |
| $256$ | $390374345071/1338557220$ | $669$ |
| $346$ | $130687286941/514829700$ | $704$ |
| $355$ | $1673839750871/20078358300$ | $50$ |
| $445$ | $344930551643/4015671660$ | $781$ |
| Sum | $1008$ | $1008$ |

The order-$g$ unrooted local defect and the rooted $Q_2\mathcal I$ pairing are therefore
$$
-\frac{1008g}{4096\pi^{10}}=-\frac{63g}{256\pi^{10}},\qquad
12\left(-\frac{63g}{256\pi^{10}}\right)=-\frac{189g}{64\pi^{10}}.
$$
Before the planar Taylor and cycle multiplier $-7\cdot3!6!2!=-60480$, the same graph sum is $-\pi/60$. For $445$ the denominator divides $2023898516640000$. The four large-field residues are $2583851,4822054,3321786,2615945$ at $6000011,6000023,6000041,6000047$. Including the independent $1009$ field gives modulus $1307690589355868259432719918779$, larger than twice the certified common-numerator bound $73530548029193358613398258230$. The replay checks each saved large-field sum against every partition key, reconstructs all five rationals, and checks the independently calibrated 77-polarization modular total. [Full ledger replay with `--verify-exact`](mixed_integration/six_conserving_nodes.py)

On the 60-polarization conserving cycle above, the five complete symmetry-weighted shape contributions divided by $\pi$ have residues
$$
(822,669,704,50,781)\pmod{1009},\qquad
822+669+704+50+781=3026=3\cdot1009-1.
$$
All $28+280+560+524+768=2160$ supported path partitions enter this sum. The common denominator of the five rationals divides $12143391099840000$, which is coprime to $1009$. The rational graph coefficient is therefore nonzero. For the second shape the six large-prime residues give $390374345071\pi/1338557220$, independently reducing to $669\pi$ in the small-prime control. [Partition ledgers](mixed_integration/six_conserving_shapes) · [Complete replay with `--verify`](mixed_integration/six_conserving_nodes.py)

The 24-term planar cycle evaluates invariant cochains as $-7$ times the single Hamiltonian fixture on this profile. Its tensor product with the transverse integer cycle annihilates every current coboundary. On the marked slot the mixed action on an unmarked output vanishes as above. The scalar component of a primitive only contributes to the scalar output. Hence the nonzero pairing obstructs a local two-loop counterterm for the bare order-$g$ defect. Every potential order-zero-in-$g$ primitive entering through the $g$ part of the classical differential would need $w=-1$, since that part raises $w$ by two. The field algebra has $w\ge0$. In the source's unrooted heat convention the physical coefficient multiplies this rational number by $-g/(4096\pi^{10})$ at two loops. [Physical prefactor and contact check](../five_theta/README.md#ordered-zero-time-contact)

### The ordinary transverse scalar six-class

For scalar cochains on positive divergence-free vector jets in three dimensions, of weight $\det^{-3}$ and Euler degree nine, the dimensions at arities three through nine are $2,16,53,74,44,11,1$. The six successive differentials have ranks
$$
2,14,39,34,10,1 .
$$
Indeed these six ranks are attained as lower bounds over the rationals by nonzero minors modulo 251. An additional degree-six class follows from the next exact test. The square-zero identities then force equality in every degree. In particular the ordinary scalar cohomology is one-dimensional at arity six and vanishes in all other arities at that weight. [Dimension and middle-rank replay](mixed_unmarked_scalar.py) · [Outer ranks with `--unmarked`](mixed_scalar_tail.cpp)

For the first-derivative representative, the marked scalar depends on the transverse function through $S^a\partial_a f$, where $S^a$ is a closed tangent-valued six-cochain. Under the divergence-free action, $\operatorname{div}[Y,S]=Y(\operatorname{div}S)$. Thus $-\partial_a S^a$ is a closed ordinary scalar six-cochain. Setting the function momentum to minus the sum of the six ordinary momenta evaluates precisely this divergence with its minus sign. On the conserving geometry above a 40-mask integer combination annihilates all 53 ordinary scalar primitives and all 150 marked scalar primitives. It pairs nontrivially with a closed first-derivative representative. Consequently the transverse factor of the nonzero anomaly, after integration by parts against $f$, generates this ordinary scalar cohomology. [Integer test on both primitive spaces](mixed_unmarked_scalar.cpp) · [Rational reconstruction of the test](mixed_unmarked_scalar.py)

### Finite one-loop counterterms on the divergence-free detector

Use full-wedge holomorphic local counterterms, formal integer powers of $g$, and the graph weights. At one loop the coefficient of $g^m$ has $w=2m$ and determinant character $\det^{-(1+m)}$. Each tree vertex supplies one incoming letter and two more per cubic one-form vertex, while each sewing consumes one. The holomorphic kernel contributes one inverse determinant per loop and the cubic one-form tensor one per power of $g$. In $\tfrac12\{I_1,I_1\}$ the coefficient of $g$ thus pairs the vector-only sector with the sector containing two one-form or scalar-antifield letters. A BV sewing removes one unit of $N=\#X+\#\beta$. An even full-wedge counterterm has even $N$. Their $N$-degrees sum to ten for the nine-vector detector.

Set $\nu=0$ and take divergence-free vectors. The first sector must have $N=6$ and the second $N=4$ to contribute. With at most four vectors, even the first variation in $\nu$ has at most five wave momenta summing to zero, which cannot supply a determinant. If the second sector has $N\le2$, two one-forms and at most two vectors have at most three independent momenta. Their two one-form polarizations cannot complete both five-determinants. The same bound is stronger for terms containing $\beta$. At $N=4$ the $\beta'\gamma$ term has five wave momenta and one one-form polarization; $\beta'^2$ has four momenta. Both densities vanish by this count. Here $\beta'=\beta+X^a\gamma_a$. Since $Q_0\gamma_a$ contains $\partial_a\beta'$, closure in the $\beta'\gamma$ channel requires invariance under adding an exact one-form to either one-form argument.

For the vector-only density the six conserved momenta have a single independent five-determinant. It changes sign under each permutation of the six slots. The remaining factor is symmetric in those slots and uses one pairing $p_j(u_i)$ per vector polarization $u_i$. For row $i$ its diagonal entry vanishes and its row sum vanishes. On a generic momentum frame the other four coordinates of each row vary freely. For a permutation with cycle lengths $r$, the trace on the product of row spaces is the product, over its cycles, of the number of labels fixed by its $r$-th power minus two. Averaging gives dimension nine. Evaluation on seven conserved transverse waves gives outgoing Chevalley–Eilenberg rank at least six modulo 1009.

The three closed tensors use the polarizations of $(\operatorname{tr}J^2)^3$, $(\operatorname{tr}J^3)^2$ and $\operatorname{tr}J^2\operatorname{tr}J^4$, where $J^a{}_b=\partial_bX^a$. The corresponding densities evaluate the polarized polynomial on $J,\partial J,\ldots,\partial J$, with wedge product of the five holomorphic differentials. Indeed $-J$ satisfies the Maurer–Cartan identity for the vector-field Chevalley–Eilenberg differential acting by Lie derivative on the form coefficients. Adding the holomorphic exterior derivative gives curvature $-\partial J$. Its sixth characteristic form vanishes by holomorphic degree. The component of the transgression in holomorphic degree five is the displayed density, so its differential is a total derivative. The three densities have independent evaluations modulo 1009. At five vector arguments no conserved momentum determinant exists. This gives exactly three scalar classes on the divergence-free restriction. [Scalar rank replay](one_loop_vector_counterterms.py)

For the order-$g$ lift use four vector and two one-form inputs. In a generic six-wave frame the coefficients of each one-form in the six momenta are defined modulo a common constant. The longitudinal Ward identity also removes the coefficient in the one-form's own momentum direction. Four differences remain. Multiplication by the five-momentum determinant expresses each difference polynomially by row replacement. Its square clears both one-form coordinates. Thus all two-form symbols arise from four independent coordinates in each of the six slots, alternated in the four vector slots and symmetrized in the two one-form slots. Conversely all these polynomials have the required determinant character and twelve derivatives. The same cycle-trace formula, averaged over $S_4\times S_2$ with the sign on $S_4$, gives dimension $83$. On generic frames this parametrization is injective.

The seven-wave equation is $Q_0 I_1^{(1)}=-Q_g I_1^{(0)}$, where $Q_g$ denotes the coefficient of $g$ in the tree differential. On the left the brackets are the vector bracket and the Lie derivative of a one-form. On the right insert the vector dual, via the five-volume, to the wedge of the two exterior derivatives of the one-forms. With the determinant and density conventions above, the ranks are

| Map | Rank modulo 1009 |
|---|---:|
| Three-vector, two-one-form primitives into the 83-dimensional slot | 3 |
| Differential from that slot | 80 |
| Same differential augmented by the three scalar sources | 83 |

The three incoming coboundaries give the rational upper bound $83-3=80$ since the differential squares to zero. Each modular minor lifts to characteristic zero by lifting the conserved momenta and solving one pivot equation per transverse polarization. The scalar and two-form replays give the same ranks at the independent prime 251 with `TWIST_COUNTERTERM_PRIME=251`; both also use 1009 by default. Thus all three scalar sources are independent modulo the rational image. No nonzero scalar combination lifts. Every closed two-one-form symbol in this slot is a boundary on the restriction. [Complete lift differential and rank certificates](one_loop_lift.py)

Consequently an order-$\hbar$ cocycle has vanishing vector-only restriction in the only contributing degree. Its vector variation annihilates all divergence-free variations, and the one-form variation of the other factor is divergence-free by the longitudinal identity. The $\nu,\beta$ sewing vanishes because the $\beta'\gamma$ symbol vanishes. Therefore the order-$g$ self-bracket of such counterterms is zero on the divergence-free nine-vector, one-form detector. This is a statement about the homogeneous full-wedge jet complex.

### Periodic local counterterms on a flat compact quotient

The detected translation-invariant local cocycle descends as a differential polynomial to the product of a circle with a flat complex five-torus, with the constant holomorphic volume form. It is nonexact there even among local two-loop primitives with arbitrary smooth periodic coefficients. Indeed, if the flat classical BV differential of a periodic local functional were this cocycle, average the functional over the compact group of all translations with Haar mass one. The differential commutes with translations, and the cocycle is fixed, so the average would be a translation-invariant primitive. Averaging leaves only the constant mode of each coefficient of its finite-jet density. Varying smooth fields supported in a single coordinate chart then gives the same constant-coefficient variational identity on flat jets. Its restriction to the three-plus-six cycle contradicts the integral pairing. The assertion concerns the descended local class and fixed lower-genus action. Finite-size effective vertices and globally holomorphic on-shell amplitudes require their own computations.

### Transverse multivectors and divergence at determinant weight $-3$

For outputs of exterior degrees $0,1,2,3$, the ordinary transverse Euler degrees are $9,8,7,6$. In increasing Chevalley–Eilenberg arity beginning at three, the cochain dimensions and the successive differential ranks are
$$
\begin{array}{c|l|l}
\text{exterior degree}&\text{dimensions}&\text{ranks}\\ \hline
0&2,16,53,74,44,11,1&2,14,39,34,10,1\\
1&4,31,75,70,25,3&4,27,47,22,3\\
2&4,21,33,18,3&4,17,15,3\\
3&2,5,4,1&2,3,1 .
\end{array}
$$
For degree one, contracting the output with the first derivative of a marked function identifies the complex with the homogeneous first-derivative sector. The mod-251 minors give the five displayed lower bounds. The nonzero divergence class at arity six and the square-zero identities give their upper bounds. For degree two the alternating first-derivative sector is the bivector calculation above. Degree three is scalar weight $(-2,-2,-2)$ after using the volume. [Tangent rank replay with `--current`](mixed_scalar_tail.cpp) · [Bivector weights](mixed_vertical.py) · [Trivector weight check](mixed_unmarked_scalar.py)

At arity five the bivector generator has divergence equal to the tangent generator. Indeed contracting the bivector with $df$ gives the marked five-current; its divergence pairs nontrivially with the marked scalar five-cycle above. The Hessian of $f$ drops out by antisymmetry. At arity six the tangent generator maps to the scalar generator by the 40-mask test. These are the only four cohomology lines; the two induced divergence maps are isomorphisms. The total complex of all four transverse multivector columns with the divergence differential is therefore acyclic at this weight. In the physical three-plus-six coupling the planar current pairs with the scalar six-class in the transverse block.

This also computes scalar Gelfand–Fuchs cohomology for the full formal divergence-free transverse algebra at weight $(-2,-2,-2)$. In the complex relative to $\mathfrak{sl}_3$, filter by the number $k$ of constant vector inputs. The other inputs are positive jets. At fixed $k$, the coefficient is an alternating $k$-form in the translations, identified by the volume with exterior degree $3-k$ in the preceding table. Its positive-jet Euler degree is $6+k$. The differential preserving $k$ is the positive-jet differential. The next differential inserts $[\partial_i,Y]=\partial_iY$ with the exterior translation sign. It is the exterior derivative on the translation slots, hence divergence under the volume identification up to a common degree sign. The two isomorphisms just computed make the second page zero. Thus every relative scalar cohomology group vanishes at this weight. For absolute cohomology the Hochschild–Serre filtration for the reductive subalgebra $\mathfrak{sl}_3$ has the relative cohomology tensored with $H^\bullet(\mathfrak{sl}_3)$ on its second page, so the same vanishing holds in every arity.

The old 77-polarization detector has now completed both marked geometries on every supported partition. Its five residues after division by $\pi$ are $(472,106,290,904,865)$, with sum $619$ modulo 1009. The integer calibration gives a factor $390$ modulo 1009, and $390(-1)=619$ in that field. The two detector calculations share the native angular integrator and use different momentum lists, different polarization cycles and different numbers of marked geometries. [Two-geometry partition records](mixed_integration/six_vector_shapes) · [One-geometry records](mixed_integration/six_conserving_shapes)

For one-loop graphs with one full-wedge holomorphic genus-one insertion, reflect all relative topological coordinates across the root time. Each normalized heat form acquires one minus sign from its topological volume factor. The gauge-fixing operator commutes with the reflection, so each propagator also acquires one minus sign. Holomorphic vertex derivatives and full-wedge vertices introduce no further relative time factor. The pullback sign from the $E$ edges divided by the relative orientation sign of the $V-1$ integrations is $(-1)^{E-V+1}=-1$ at one loop. This applies at every positive cutoff to each placement of the heat edge. Non-self zero-time contacts inherit the same parity by taking the limit at fixed other cutoffs; contracting that edge preserves $E-V+1$. For a heat self-edge the pullback of the topological difference one-form to its own diagonal is zero. Thus the one-loop sewing term with this genus-one insertion also vanishes on the detector.

For the three-plus-six planar soft test, one can include the full separated planar-current contribution before taking the Taylor coefficients. The two unmarked invariant three-current columns are $A,B$ in the block homotopy above. In the positive-function part a primitive with unmarked transverse output contributes through $(\partial_a f)Z$ acting on that output. Its two ordinary planar scalar slots are the two quintic Hamiltonians paired by $(F,G)_5$. On the frame of the soft test this adds one current column $(1-b,a-1)$, the constant marked vector. Alongside the two planar coboundaries these give five possible additional columns. The coefficient matrices of these columns and their augmentation by the Hamiltonian current have ranks $(4,4)$ through quadratic order and $(5,6)$ through cubic order. Thus the cubic Taylor pairing survives both the separated contribution and the mixed action from an unmarked primitive on this relative three-plus-six test. [Exact Taylor ranks](mixed_planar_soft.py)

The automorphism action on a fully polarized assignment is free: the ten interior positions have distinct labels, and a graph automorphism fixing every position fixes the endpoints and the three paths. All assignments in an orbit have the same signed contraction. Thus the sum over assignments divided by the graph symmetry equals the sum over orbit representatives. In the native kernel the explicit reversal factor cancels the final division by two. For a *complete* shape one may consequently use the ordered component denominator at symmetry one, then remove the common factor shared with $60480$. In shape order these denominators are
$$
1214339109984000,\quad4047797033280000,\quad2023898516640000,\quad4047797033280000,\quad2023898516640000.
$$
Their least common multiple is $12143391099840000$. The positive-resistance operator bound applied before summing the path orders bounds the common numerators at $98,99,98,98,96$ bits respectively. The rational $346$ and $355$ contributions on the same cycle are $130687286941\pi/514829700$ and $1673839750871\pi/20078358300$. [Resistance norm bound](mixed_integration/theta_spectral_bound.py) · [Rational shape replay](mixed_integration/six_vector_crt.py)

### Translation descent for the field-strength cochain

Hahner, Paquette and Raghavendran give the strict field-strength map in §2.2.4: it takes $\gamma$ to $F=\partial\gamma$, kills $\beta$, and fixes $\mu,\nu$. They identify its flat minimal model with $E(5|10)$. In the tree differential the ideal generated by $\nu$, $\operatorname{div}X$ and their derivatives is preserved at every $g$. In particular the extra vector in $QX$ is dual by the volume to $gF\wedge F/2$ and has zero divergence since $\partial F=0$. On the quotient $QF=-\mathcal L_XF$. Thus the strict field-strength differential is the binary Chevalley–Eilenberg differential of the divergence-free vectors and closed two-forms with the source odd-odd bracket. All terms of the projected two-loop defect have ordinary cochain arity $(11-2m)+(2m-1)=10$. The individual longitudinal Ward identities express their densities in these field strengths. [Strict field-strength model and minimal algebra, §§2.2.4 and 4.2.1](https://arxiv.org/html/2506.02246v1) · [Ward identities](../gamma_ward.md)

The passage to densities modulo the five holomorphic total derivatives is exact in the arities used here. For one momentum $p$, the divergence-free vector symbol module is $\mathbb C[p]^5/\mathbb C[p]p$. It embeds in the exterior square by $v\mapsto p\wedge v$. The two-form symbol module is the span of $p_i\gamma_j-p_j\gamma_i$ in the free potential module. Its quotient there is the ideal $(p_1,\ldots,p_5)$, by exactness of the polynomial Koszul complex on the momentum coordinates. Each of these modules is torsion-free and free locally wherever one momentum coordinate is invertible.

In polarized arity $n$ take their tensor product over the separate momentum slots. Off the joint origin one slot is locally free. The five components of the total momentum form a regular sequence there: before localization they are independent polynomial coordinates in that slot over the tensor product of the other slots. At the joint origin the depth is at least $n$, because one momentum coordinate from each slot gives a regular sequence on the tensor product. For $n\ge5$, at each of the first five total momentum components, the origin has positive remaining depth after dividing by the preceding components. It cannot be an associated prime. Together with regularity off the origin this proves regularity everywhere. Polarization symmetries are direct summands in characteristic zero. The same argument applies with one slot replaced by the ideal $(p_1,\ldots,p_5)$. The long exact Koszul sequence therefore keeps the inclusion from field strengths to potentials injective after imposing total momentum zero, one slot at a time. In particular a zero of the integrated field-strength differential can be checked in the potential symbols.

On the branch with vanishing one-loop defect, BV consistency gives closure of the two-loop defect. Its integrated strict field-strength class is nonzero at each fixed $g\ne0$: in the one-two-form sector a primitive would have to produce the graph by the vector bracket and Lie derivative alone, since the odd-odd part raises the number of two-form inputs by two. Pulling such a primitive back to potentials contradicts the three-plus-six cycle. This uses only a finite number of coupling sectors at arity ten.

Resolve the integrated strict cochains by the Koszul complex of the five total derivatives. Its total degree is ordinary cochain arity minus Koszul exterior degree. The regularity just proved identifies it with the integrated complex in degree ten and its neighboring degrees. Cartan insertion of each constant vector field gives a nullhomotopy for its corresponding total derivative on absolute cochains. The five insertions anticommute and commute with translations. Conjugation by the finite exponential of their sum tensored with contraction of the matching Koszul generator removes the Koszul differential. Thus in total degree ten the absolute summands have cochain arity $10+j$ tensored with exterior degree $j$ of the translations, for $0\le j\le5$. The inner $\mathfrak{sl}_5$ action on absolute cohomology vanishes by Cartan insertion. Taking $SL_5$ invariants leaves precisely $j=0,5$. Consequently the detected class has a nonzero image in
$$
H^{10}_{\rm GF}(E(5|10))\;\oplus\; H^{15}_{\rm GF}(E(5|10))
$$
with ordinary super Lie algebra cochain arity. The [arity-fifteen test](#nonzero-absolute-arity-fifteen-projection) below distinguishes the projections. Relating the physical obstruction degree to the twelfth group discussed in §4.3.2 of Hahner, Paquette and Raghavendran requires their grading comparison.

### Commuting momentum eigenmodes

Every two-loop sector at $\nu=0$ contains $11-2m\ge1$ vector inputs. Suppose one nonzero vector wave $u\exp(q\cdot z)$ commutes with every other vector wave and every closed two-form wave. Divergence freedom gives $q(u)=0$. For another vector wave $v\exp(p\cdot z)$, its commutator has coefficient $p(u)v-q(v)u$. If $u,v$ are independent this gives $p(u)=0$. For proportional polarizations this follows directly from $p(v)=0$. For a nonzero two-form wave $F\exp(p\cdot z)$ the vanishing Lie derivative gives
$$
p(u) F+q\wedge\iota_uF=0,\qquad p(u)\,\iota_uF=0,
$$
where the second equality is the contraction of the first with $u$. If $p(u)\ne0$, the first equality now forces $F=0$. Hence every external momentum annihilates $u$. Both five-component blocks of all ten interior columns of $D_\theta$ then lie in the four-dimensional annihilator of $u$, so their combined rank is at most eight. All five graph integrands vanish in every coupling sector before Wick or angular integration. This excludes detection on pairwise commuting divergence-free momentum eigenmodes. It also suffices for a single vector input to centralize the other nine inputs.

More generally, for any $r$-dimensional span of vector polarizations, each column whose momentum annihilates that span lies in a common space of dimension $2(5-r)$. An invertible $D_\theta$ admits at most $2(5-r)$ such columns. For each vector input this requires two momenta with nonzero contraction with its polarization; these two must include different paths since each path has only one circuit direction. Each of these waves has a nonzero bracket with that vector. Any mutually commuting collection of vector waves with polarizations spanning four dimensions already forces a zero: its at least four columns lie in a space of dimension two.

### Nonzero absolute arity-fifteen projection

The ordinary arity-ten summand in translation descent maps to integrated cochains by taking the density of an absolute cocycle modulo total derivatives. Indeed the Cartan conjugation acts identically on exterior degree zero. Consider this map on the component with exactly one two-form. There the two-form–two-form bracket has no source, so the component of every closed absolute ten-cochain is closed for the $g=0$ differential. Pull it back to the potential complex by $F=\partial\gamma$. In the variables of the strict tree differential,
$$
Q\beta'=-X^aD_a\beta',\qquad Q\gamma_a=D_a\beta'-\mathcal L_X\gamma_a.
$$
Use ordinary cochain degree one for $X,\gamma$; the displayed differential gives degree two for $\beta'$. At degree ten and $w=1$ the scalar subcomplex has eight vector arguments and one $\beta'$, and the quotient has nine vector arguments and one $\gamma$. On each diagonal the marked field transforms by Lie derivative.

Filter either diagonal by ordinary cochain degree minus the number of derivatives of the marked field. For its marked momentum $p$, the filtration-preserving differential is multiplication by $-X^i(0)p_i$, with the five odd translation coordinates $X^i(0)$. Indeed only the undifferentiated vector in the transport term raises the marked derivative order by one. The bracket on the vector jets and every other Lie derivative term preserve or lower that order. At fixed torus weight and arity this is a bounded filtration. The marked potential jet module is a free polynomial module in $p$, with either its covector fiber or its scalar fiber. The five-variable Koszul cohomology therefore leaves just $X^1(0)\cdots X^5(0)$, the undifferentiated marked field, and the remaining differentiated vector jets.

The translation product has determinant weight $+1$. On the component of the graph at $w=1$ the target weight is $\det^{-3}$, so the remaining tensor must have weight $\det^{-4}$. In the one-form diagonal it has four vector arguments. On decomposable vector jets these supply at most four distinct momentum covectors. Each covariant five-determinant must therefore use the marked one-form, whose degree is one. The required four excess covariant determinants cannot occur. In the scalar diagonal there are three vector arguments and no marked covector. Polarization over characteristic zero spans all jets, and taking invariants commutes with the filtration. Hence both diagonal absolute cohomologies, and therefore the whole absolute potential cohomology at degree ten in this component, vanish.

An elementary character test gives the same vanishing. The four vector-jet momenta annihilate a common line. Scale that line by $t$ and fix a complement. Each momentum stays fixed, every vector polarization is polynomial in $t$, and the single marked covector has at most a simple pole. A tensor of character $\det^{-4}$ must scale by $t^{-4}$, forcing its value to vanish. With three vector jets and a marked scalar there is no pole. The divergence-free pure jets used in this test span each constrained jet space: a functional on degree-$k$ pure jets has the form $F_i(p)u^i$ for a homogeneous polynomial covector $F_i(p)$, where $p$ is the momentum and $u^i$ the polarization. If it annihilates every $u^i$ with $p_i u^i=0$, then $p_iF_j=p_jF_i$. Coprimality of the distinct momentum coordinates gives $F_i=p_iG$ for a polynomial $G$ of degree $k-1$, exactly a functional through the divergence.

If the projected two-loop class came entirely from ordinary absolute arity ten, its one-two-form component at the graph weight would consequently be exact even before integrating after pullback to potentials. The integrated three-plus-six pairing is nonzero and kills every such primitive. The absolute arity-fifteen component must therefore be nonzero:
$$
H^{15}_{\rm GF}(E(5|10))\ne0
$$
for continuous polynomial-jet super Lie algebra cochains in the ordinary arity grading. This argument determines the arity-fifteen projection by its nonzero pairing after descent. It leaves open the dimension of that group, a comparison with the published physical degree twelve, and whether there is an additional arity-ten contribution in the strict field-strength model. [Translation descent](#translation-descent-for-the-field-strength-cochain) · [Integrated witness](#nonzero-integrated-order-g-obstruction-on-the-three-plus-six-cycle)

### Holomorphic position-dependent two-loop counterterms

The potential complex itself admits the same translation descent at $w=1$. Its marked slot is a free polynomial jet module. Cartan insertion of a constant vector acts only on $X$ in the displayed $\beta',\gamma$ variables and has commutator with $Q$ equal to minus the corresponding total derivative. At degree ten its invariant absolute arity-ten cohomology vanishes at the graph weight by the preceding calculation. Its nonzero integrated class therefore has nonzero exterior-degree-five component also after pullback to potentials.

Allow polynomial coefficients in the holomorphic position $z$ in the full-wedge local functional. The Koszul differential now uses the sum of the symbol total momentum and the derivative on each position coefficient. Filtering by polynomial degree in $z$ recovers the regular symbol Koszul complex, since differentiating the coefficient lowers its degree. This resolves the quotient by all five total derivatives in the same arities. The Cartan conjugation removes the symbol momentum, leaving the position derivatives alone. By volume duality their homological Koszul complex is the polynomial de Rham complex in complementary degree. Contraction with the Euler vector field divided by polynomial degree plus form degree contracts every positive total degree. Thus its homology is exactly the constants in exterior degree five, and the inclusion from constant coefficients acts identically on this summand. The nonzero exterior-degree-five component therefore obstructs every polynomial-coefficient primitive on a flat chart.

A finite-jet primitive with holomorphic position coefficients would give a polynomial one. Take its Taylor expansion at the chart origin and project to the scalar-dilation weight of the order-$g$ equation. Each fixed jet monomial admits at most one homogeneous Taylor degree at this weight. There are finitely many monomials at the specified field arity and finite jet order, and the classical order-zero differential and total derivative relation respect the dilation. The resulting finite polynomial coefficient identity contradicts the preceding calculation. Thus the obstruction persists for arbitrary holomorphic position dependence of a full-wedge local two-loop counterterm on a flat chart, at the prescribed lower-genus action.

### Smooth component counterterms on a flat chart at fixed lower genus

At order $g$ work in the divergence-free potential quotient and retain the unique marked input $\gamma$ or $\beta'$. Restore all component forms of the superfields in $\Omega^\bullet(\mathbb R)\otimes\Omega^{0,\bullet}(\mathbb C^5)$, with free differential $d_\tau+\bar\partial$. A component with form index set $I$ has ordinary cochain degree equal to the bottom degree minus $|I|$. The full-wedge component of the ten-input graph has degree $10-6=4$. All external form indices on a theta cut combine by the exterior product: the twelve propagators and the heat kernel already have $12\cdot5+6=66=6(12-1)$ relative form indices and only contain coordinate differences. Every additional relative form therefore vanishes. Terms with additional topological or antiholomorphic derivatives can be separated as below.

Resolve component local functionals by the Koszul complex for all eleven translations. On a polarized symbol the marked potential slot is free in all eleven momentum variables. Filter by derivative order on that slot, assigning one additional unit per Koszul exterior generator. The leading Koszul differential multiplies by those eleven independent variables; it is exact in positive exterior degrees. The remaining terms differentiate the position coefficient or another slot. Eliminating the highest filtered piece and repeating terminates for every finite-jet symbol. Thus the resolution also works with arbitrary smooth coefficient functions and arbitrary component contractions at $w=1$.

For each of the six directions of $d_\tau+\bar\partial$, contraction of the superfield by its coordinate vector gives a derivation of component observables: acting on a bottom field extracts the corresponding one-form component; in higher degree it adds that index with its exterior sign. Its commutator with the free differential is the total field derivative in that direction by the de Rham Cartan identity. It commutes in the graded sense with the tree differential on the superfields, because their brackets use wedge products and holomorphic field derivatives. These six contractions anticommute with each other and with the five insertions of constant holomorphic vectors. Hence the eleven Cartan conjugations leave the absolute component differential and the eleven derivatives of the smooth position coefficients. Via the constant volume, the latter Koszul complex is the smooth de Rham complex in complementary degree. On a ball its homology consists of constants in exterior degree eleven. Evaluation of bottom fields with every positive form component and every $\tau,\bar z$ derivative set to zero is a chain retraction of absolute component cochains onto the holomorphic bottom cochains.

Embed the holomorphic Koszul representative in the bottom fields, append the six Koszul generators for $\tau,\bar z$, then undo their six Cartan conjugations. In exterior degree zero in those six generators this extracts, up to the common orientation sign, exactly the six-form component of the polynomial in the superfields, by the exterior Leibniz rule. Thus the full-wedge lift sends the holomorphic exterior-degree-five contribution to exterior degree eleven, where the bottom-field retraction recovers its nonzero absolute arity-fifteen class. For any additional terms, count the sum of the number of $\tau,\bar z$ field derivatives, the number of component form indices, and the number of corresponding Koszul generators. The translation resolutions and both types of Cartan conjugation preserve this count at constant coefficients. A term of external form degree six and positive derivative order has count greater than six, and therefore has zero bottom-field retraction in exterior degree eleven. This leaves the computed holomorphic graph component.

It follows that the order-$g$ defect at fixed lower-genus action has no finite-jet local two-loop primitive on a ball, even with smooth position coefficients and arbitrary component contractions. This statement concerns the two-loop primitive. Allowing a different one-loop solution also introduces its self-bracket and one-loop sewing terms.

### Smooth representatives of one-loop counterterms at the source weights

For this comparison restore $\rho=1-\nu$, $X=\mu/\rho$ and $\beta'=\beta+X^a\gamma_a$ before taking the divergence-free quotient. The tree differential gives $Q\rho=-\partial_a(\rho X^a)$, $Q\gamma=\partial\beta'-\mathcal L_X\gamma$. Its terms at $g=0$ on $X,\beta'$ are $-X^a\partial_aX$ and $-X^a\partial_a\beta'$. At nonzero $g$ the additional term in $QX$ is the vector dual to $gF\wedge F/(2\rho)$ and that in $Q\beta'$ is its contraction with $\gamma$. Both are independent of $X$. Insertion of a constant vector into $X$ therefore has commutator minus the holomorphic total derivative on all four fields, including $\rho$. The other six component Cartan operations act as above.

Before imposing constraints every polarized nonconstant field symbol is free in the eleven momentum variables of any chosen slot. The same resolution thus computes the local complex in all positive arities; the graded polarizations are direct summands. After the eleven Cartan conjugations, smooth position coefficients on a ball retain only exterior degree eleven. At the level of absolute component cochains the bottom-field inclusion is a quasi-isomorphism. At fixed polarized arity the one-letter $d_\tau+\bar\partial$ complex is the six-variable polynomial Koszul complex of the individual nonholomorphic momenta and form components. Its homology retains the bottom component at zero nonholomorphic momentum. The linear holomorphic differential commutes with that contraction. Tensor the contractions in separate slots, then take the graded polarization symmetry. All remaining tree terms increase arity, so recursively removing their higher-arity remainders gives the same comparison on field-degree completions, with finitely many operations in each arity.

It follows that every smooth finite-jet one-loop counterterm class on a ball, at the source genus and coupling weights, has a representative obtained by appending the eleven translation generators to an absolute holomorphic bottom cocycle and undoing their Cartan conjugations. In exterior degree zero the five holomorphic insertions give a constant-coefficient holomorphic jet polynomial; the six component contractions take its full wedge. The construction respects $SL_5$, and the simultaneous rotation of the translation generators and their contractions preserves the source character. In the cotangent grading each holomorphic insertion removes only an $X$, while component contraction preserves field type. The resulting representative is in the full-wedge complex used for the genus-one rank test.

Replacing a closed order-$\hbar$ deformation by a cohomologous one is an order-$\hbar$ BV canonical transformation generated by its primitive. At finite cutoff the quantum BV differential transforms the exponential of the action by this canonical transformation. The induced change in the first nonzero, order-$\hbar^2$ defect is a classical BV coboundary after the corresponding shift of the order-$\hbar^2$ local action. Equivalently, integrating the canonical homotopy includes its order-$\hbar^2$ self-bracket and BV-Laplacian terms. The nonzero obstruction class is therefore common to smooth representatives of each one-loop branch covered by the source weights and the finite rank calculation. Field-independent constants have zero sewing and bracket. [One-loop restriction and ranks](#finite-one-loop-counterterms-on-the-divergence-free-detector)

### Arity and the loop power in the quantum BV exponential

The source assigns only total parity to the physical minimal twist. Ordinary absolute cochain arity and the degree of the obstruction acting on the quantum exponential can nevertheless be compared directly. On full-wedge densities the tree interaction has ordinary degree $3-6=-3$. Sewing an $X,\gamma$ pair removes two ordinary field degrees and one set of six form indices, so the BV bracket has degree $-2+6=4$. The $\nu,\beta'$ pair also has total degree two. The classical differential consequently has degree $4-3=1$. Keeping $Q+\hbar\Delta$ homogeneous in this arity count assigns the loop power the forced bookkeeping shift $1-4=-3$.

At loop number $n$, a trivalent defect graph of the source has $5n+2$ vertices and $6n+1$ edges. Its external arity, full-wedge degree, and degree of the top absolute descent are respectively
$$
3(5n+2)-2(6n+1)=3n+4,\qquad 3n+4-6=3n-2,\qquad 3n-2+11=3n+9.
$$
The order-$\hbar^n$ master-equation defect enters the quantum differential of $\exp(S/\hbar)$ with the power $\hbar^{n-1}$. Including that power gives $3n+9-3(n-1)=12$ after top descent. At two loops the calculation is $15-3=12$. This distinguishes the degree-twelve quantum obstruction with its loop parameter from the ordinary arity-fifteen cochain computed above. The single sentence in §4.3.2 of Hahner, Paquette and Raghavendran supplies no convention identifying its stated group with this loop-graded group. [Published discussion, §§1.2 and 4.3.2](https://arxiv.org/abs/2506.02246)

The absolute class also has a definite principal weight, measured on its inputs. In the standard depth-two grading of $E(5|10)$ a vector of Taylor order $k$ has degree $2k-2$ and a two-form of Taylor order $k$ has degree $2k-1$. These degrees add under each of the vector, mixed and two-form brackets. In the $g^m$ graph the total Taylor order after extracting the $2m-1$ individual field strengths is $(22+m)-(2m-1)=23-m$. Its ten inputs therefore have principal weight
$$
2(23-m)-2(11-2m)-(2m-1)=25.
$$
Each holomorphic total derivative increases the input weight by two, so its Koszul generator has weight two. The absolute arity-fifteen component in exterior degree five therefore has input weight $25-5\cdot2=15$. Every $g^m$ contribution contains an odd number of two-form inputs. This specifies the nonzero absolute component by ordinary arity fifteen, input principal weight fifteen, and odd internal parity.

On the leading one-two-form component the ordinary absolute fifteen-cochain has fourteen vector inputs and one closed two-form. Removing the five Koszul directions changes its torus character from $\det^{-3}$ to $\det^{-2}$. Writing $S_5$ for the formal divergence-free vector fields and taking the continuous dual of the closed formal two-forms gives
$$
H^{14}\!\left(S_5;(\Omega^2_{\rm cl})^\vee\right)_{\det^{-2}}\ne0.
$$
Its pullback to the dual of one-form potentials is also nonzero. Indeed setting $\beta'=0$ is a cochain map at $g=0$ and the three-plus-six detector kills all current coboundaries already in that quotient. Its ordinary absolute arity-ten potential cohomology at $\det^{-3}$ vanishes by the marked-momentum Koszul calculation. The same descent argument therefore detects the image in $H^{14}(S_5;(\Omega^1)^\vee)_{\det^{-2}}$. [Component Cartan sign and grading replay](component_descent.py)

On the quotient $\beta'=0$ this absolute representative is explicit. In the local volume frame multiply the integrated current density $\gamma_aJ^a$ by $X^1(0)\cdots X^5(0)$, the five translation components in that order, and evaluate at the origin. The current closure equation at the beginning of this file gives $Q(\gamma_aJ^a)=-X^bD_b(\gamma_aJ^a)$. The differential of the product of translation components vanishes by $D_aX^a=0$ and its product with any $X^b(0)$ is zero. This proves absolute closure. Successive insertions of the five ordered constant frame vectors strip the translation product with unit coefficient: all nine vector slots in $J^a$ have positive derivative order. Thus the top Koszul descent recovers the measured integrated density in this orientation.

More generally multiplication by the five translation components identifies the cochains of the origin stabilizer on an undifferentiated marked one-form with the surviving row of the marked-momentum Koszul filtration. It intertwines the differentials: the transport of every positive jet contributes a repeated translation component, the Jacobian acts on the marked fiber, and the translation product has zero differential. The inclusion induces an isomorphism on that first page, hence on cohomology. The absolute degree-fourteen coefficient class consequently has origin-stabilizer degree nine and residual input character $\det^{-3}$. In its Levi decomposition the only surviving relative degree is nine. For Levi degrees $5,7,8,9$ the relative degree is at most four and the cochain spaces vanish. For Levi degree three the [second differential](second_differential.cpp) gives zero relative degree-six current cohomology.

In a trivalent $n$-loop sector with $m$ cubic one-form vertices and no external $\nu$, the same count gives $4n+3-2m$ external vectors and $2m-n+1$ one-forms. The degree-$5n$ holomorphic circuit factor and $5n+2+m$ vertex derivatives give $10n+2+m$ momentum powers. Extracting the field strengths leaves Taylor order $11n+1-m$. Hence the input principal weight before descent is
$$
2(11n+1-m)-2(4n+3-2m)-(2m-n+1)=15n-5.
$$
The top absolute component has input principal weight $15(n-1)$ and arity $3n+9$. In particular five times its ordinary arity minus its input principal weight equals $60$ in every such sector. This is a selection rule for possible graph classes.

For reduction along a flat complex torus of dimension $r$, all external zero-mode momenta annihilate its tangent space. The ten columns of $D_\theta$ then have rank at most $2(5-r)$ on every graph. Thus the flat eleven-dimensional two-loop local defect restricts to zero on the entire holomorphic zero-mode sector for every $r\ge1$. A nonzero contribution with Kaluza–Klein momenta needs at least $2r$ external inputs carrying compact momentum. For $r=1$ at least two distinct paths must carry such an input. These statements evaluate the eleven-dimensional local vertex on external modes; computing the effective action after integrating out massive compact modes is a separate operation. [Circuit determinant](../../VERIFY.md)

On the original 77-polarization test the independently reconstructed $256$ coefficient after division by $\pi$ is
$$
-\frac{89550553293577907425050546913660403371042796950660845529966824566311694175152257027627295466951638}{98423325},
$$
which reduces to $106$ modulo $1009$. Its exact five-shape sum is $1008$ times the integral calibration above. The original $445$ can therefore be recovered by subtracting its first four exact coefficients and checking the independent direct $445$ partition sum modulo $1009$. [Two-geometry replay with `--calibrated-table`](mixed_integration/six_vector_crt.py)

### Constant two-form background and the further twist

At $\nu=0$ each interior two-form carries at least one derivative after passing from $\gamma$ to $F=\partial\gamma$. Its vertex supplies the individual exterior derivative and its distinct column in $D_\theta$ supplies one additional momentum. The five two-loop graph densities are therefore invariant under a constant shift of $F$. In particular every insertion of the constant decomposable background $dz_1\wedge dz_2$ vanishes in the strict graph. This is precisely the field strength of the further-twist potential $\gamma_{nm}=(z_1dz_2-z_2dz_1)/2$ in the published flat model. [Raghavendran–Saberi–Williams, §4.2, equation (45)](https://arxiv.org/abs/2111.03049)

The differential of the further-twisted formal field-strength Lie algebra takes a vector $X$ to its Lie derivative of $dz_1\wedge dz_2$ up to sign, and takes a closed two-form $F$ to the vector dual to $dz_1\wedge dz_2\wedge F$. For a closed vector, $\iota_X(dz_1\wedge dz_2)$ is a closed planar one-form independent of the other three coordinates. Its planar component is a Hamiltonian vector field in $z_1,z_2$. The remaining vector is divergence-free in the transverse three variables. Its contraction with the transverse volume is transversely closed; choose a transverse one-form primitive by the formal Poincaré homotopy. The full exterior derivative of that primitive has the prescribed transverse two-form component, proving that every such transverse vector is a boundary. For a closed odd element $F$, its two-transverse component is zero. Take a formal primitive of $F$ and subtract an exact one-form to remove its transversely closed transverse component. The resulting planar one-form is $\iota_X(dz_1\wedge dz_2)$ for a planar vector $X$. Complete $X$ to a divergence-free vector by integrating its negative planar divergence in one transverse coordinate. This writes $F$ as a boundary. Thus the planar Hamiltonian vectors embed as a strict quasi-isomorphism on the formal field-strength algebra.

The shifted two-loop density evaluates identically to zero on these residual representatives: the constant background cannot saturate an $F$ slot, and every original graph requires at least one such slot. Equivalently their two-dimensional holomorphic momenta give circuit rank at most four. This computes the restriction of the minimal-theory vertex on the further-twist representatives.

The absolute fifteen-class also survives the potential-field central extension. The central cocycle $\kappa(X,X';F)=F(X,X')|_0$ has ordinary arity three on $E(5|10)$. Its pullback to potentials is $Q\beta(0)$ in the source's polarized convention. Thus the pullback of every product of $\kappa$ with a closed ordinary twelve-cochain is exact. The top-descent class just detected has nonzero potential pullback already in its one-two-form sector. Consequently its ordinary arity-fifteen class has nonzero image in
$$
H^{15}_{\rm GF}(E(5|10))\big/\bigl(\kappa\,H^{12}_{\rm GF}(E(5|10))\bigr).
$$
The analogous leading component in $H^{14}(S_5;(\Omega^2_{\rm cl})^\vee)$ survives quotient by products of the same cocycle, now counted in two vector arguments, with scalar degree-twelve classes of $S_5$. [Central cocycle in the source model](../cohomology_11d.md)

### Action of the obstruction on the diffeomorphism ghost

In the source convention $Q_0\mu^a=\delta_L W/\delta\gamma_a$. At order $g$, on the strict quotient the obstruction has the variational form $\int\gamma_aJ^a$ with nine vector inputs. Its derivative with respect to the even one-form is exactly $J^a$. For the quantum BV differential of the effective action, the BV derivation identity, the odd Jacobi identity and the square-zero BV Laplacian give
$$
\left(Q+\{I,-\}+\hbar\Delta\right)^2
=\left\{QI+\tfrac12\{I,I\}+\hbar\Delta I,-\right\} .
$$
Thus the order-$g\hbar^2$ failure of nilpotence on the holomorphic diffeomorphism ghost has the current $J^a$ as its local coefficient. Pairing the output with the specified planar one-form and the nine-vector product cycle gives $-63g\hbar^2/(256\pi^{10})$ in the unrooted convention. The divergence $D_aJ^a$ vanishes identically on the restriction by the longitudinal Ward identity, so its pairing with any exact one-form is zero. This specifies a component of the local quantum gauge-symmetry failure directly on the ghost. [Source tree differential](../cohomology_11d.md) · [Exact normalized pairing](#nonzero-integrated-order-g-obstruction-on-the-three-plus-six-cycle)

### Exact defect in a constant scalar background

Let $\rho=1-\nu$ be constant and invertible. Holding the original nine $\mu$ inputs fixed, all eleven $b$ vertices in the order-$g$ theta sum contribute $\rho^{-1}$ each. Holding the divergence-free vectors $X=\mu/\rho$ fixed instead contributes $\rho^9$ from those inputs. Thus on the same balanced cycle the unrooted pairings at fixed $\mu$ and at fixed $X$ are respectively
$$
-\frac{63g}{256\pi^{10}\rho^{11}},\qquad -\frac{63g}{256\pi^{10}\rho^2}.
$$
At $\nu=0$ the corresponding pairings with $K$ additional labeled zero-momentum $\nu$ inputs are the zero-background value times $(K+10)!/10!$ and $(K+1)!$, respectively. These statements include the graph-symmetry factors and the complete five-shape sum. [Scalar attachment generating function](../constant_nu.md)

For a trivalent $n$-loop sector with $m$ cubic one-form vertices, the $5n+2-m$ remaining vertices each contribute $\rho^{-1}$ and its $4n+3-2m$ external vectors each contribute $\rho$ at fixed $X$. Its coupling dependence is therefore
$$
\rho^{4n+3-2m-(5n+2-m)}g^m=\rho^{1-n}(g/\rho)^m.
$$
This also follows on the strict constant-density fields from the tree action: its $\gamma DX$ and $X^2\partial\gamma$ terms acquire $\rho$, while $g\gamma(\partial\gamma)^2$ acquires none. The quantum exponential depends on $S(g/\rho)/(\hbar/\rho)$ there.
