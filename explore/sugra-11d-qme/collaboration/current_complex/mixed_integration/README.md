# Exact Schwinger integration on the 166 shape

## The prime seven after summing path orders

Fix the five labels on each long path. The determinant in label order is constant on this partition. Ordinary vector permutations have the same combined exterior and determinant sign. The remaining position sign for the one-form is $(-1)^{9-g}$ at its zero-based slot $g$. The two graph symmetries and the scalar denominator $158400$ are units at seven.

For five waves on a path with six resistance fractions, let $Q$ be the prefix momentum in the current gap, $s$ an observer direction, and $e_i$ the commuting square-zero wave labels. In the coefficient with component multiplicities $\alpha$ and $d=|\alpha|$, the normalized resistance moment of its summed transfer is

$$
\frac{5!\,\alpha!}{(5+d)!}[e_1\cdots e_5 s^\alpha]T^{5+d}.
$$

A wave letter in $T$ inserts its transfer, an observer letter multiplies by $s\cdot Q$. Expanding the power places observations in the six gaps. A specified gap with multiplicities $\beta$ has $|\beta|!/\beta!$ orders, exactly the factor obtained by expanding the simplex moment.

Every ordinary directed segment starts at a transverse output. The trivalent bracket sends its two vector inputs $z_1,z_2$ into $k(z_2)z_1-k(z_1)z_2$, annihilated by its outgoing momentum $k$ for arbitrary inputs. The cubic tensor is separately transverse in each output. Hence each ordinary matrix can be replaced successively by its vector-field Lie derivative, or its dual when the path is read backwards. This also applies to a circuit through a trivalent bracket. The displayed endpoint tensors can be retained throughout the contraction.

Multiplying the ordinary transfers on one side of the marked letter by minus one absorbs its position sign up to a constant. The marked source acts by contraction with $d\gamma$. The marked sink is the volume dual of $d\alpha\wedge d\gamma$ on a one-form $\alpha$. Both maps commute with the divergence-free Lie action in the sense $[\mathcal L_X,G_\gamma]=G_{\mathcal L_X\gamma}$. The projected observer adds the central subtraction $-s\cdot\ell$ to the constant vector-field derivative, where $\ell$ is the first-edge momentum.

Over the integers localized at seven, write $T^7=U+V+7R$. For each positive wave subset, the reduction of $U$ is its coefficient in the characteristic-seven restricted Lie derivative. Mixed seventh-power terms of ordinary Lie derivatives are iterated Lie brackets. Equivalently the effective vector is the positive wave part of the sevenfold vector-field derivation applied to the coordinate functions. The term containing the marked letter is $G_{\mathcal L_Z^6\gamma}$ for the ordinary vector field $Z$. This follows from $\sum_{j=0}^6 X^j G X^{6-j}=(\mathrm{ad}\,X)^6G$ in characteristic seven. These expressions have integral lifts, with initial-current degree one for an ordinary transfer or marked sink, and zero for a marked source. At degree one the ordinary symbol is a scalar identity. The empty wave component is $V=\sum_i s_i^7 Q_i^7$ modulo seven. Every such term pays a factor of seven through $\alpha!$. The block $R$ with $r$ waves has initial-current degree at most $r\le5$ since its original seven-letter word contains $7-r$ observers of degree zero. The source saves one further degree when present. At zero waves the degree is zero.

Put $d=d_B+d_C$. The two $5+d_P$ powers contain $k=\lfloor(5+d_B)/7\rfloor+\lfloor(5+d_C)/7\rfloor$ full seven-letter blocks. This is the total denominator valuation from the two factorials for $d_P\le10$, $d\le11$. Every $V$ pays one factor through the factorial of its observed component and every $R$ pays an explicit seven. Each remaining $U$ reduces the initial-current degree bound by six. Counting both path lengths, the two endpoints and the marked letter therefore gives a total bound $d+12-6q$ for a term with $q$ good blocks $U$. The two bracket endpoints contribute two degrees for an interior cubic sink. A cubic endpoint and a bracket contribute three, and the source on the path removes one. The substitution $\ell_P=w_P(\hat m_B+\hat m_C)-\hat m_P$ requires total initial-current degree $d$ for this moment grade.

Three good blocks fall below degree $d$. Two good blocks with an additional $V$ or $R$ also fall below it, since this additional block has degree at most five for its seven slots. For precisely two good blocks, equality requires every remaining wave transfer to take its highest symbol, every observer to occur in a fused block, and every endpoint momentum to use its initial part. All ordinary symbols are scalar. The only possible covectors entering the volume are then $\ell_B,\ell_C$, and the momentum and polarization of the bare or fused one-form. Since $\ell_A=-\ell_B-\ell_C$, their fifth exterior power vanishes. Every surviving term has at most one unpaid factor of seven.

Consequently the normalized **order-summed resistance integral has a denominator with at most one factor of seven**. The subsequent conductance moment of degree at most ten has denominator dividing $13\cdot15\cdots31$, with one factor of seven. The normalized full shape, and the mixed detector on it after restoring its stated normalization, have at most two factors of seven in their denominators. The argument applies to each fixed partition, before the determinant, physical detector weights, or orientation sum.

`resistance_check.cpp` independently expands the original cubic Wick words over $\mathbb Z/49$, sums $120$ orders on each long path, inserts two observer directions per path, and extracts the total grade $d_B+d_C$ in the two first currents. Six conserved fixtures give 1008 checks divisible by seven, 48 divisible by 49, and 533 lower-grade controls nonzero modulo seven. `frobenius.py` and `gamma_fusion.py` compare the full and fused transfers on the ordinary, source and sink fibers.

## Ordered scalar fixture

Take vertex order `0146723589`, branch zero, polarization mask zero, function selector zero and one-form axis zero of the mixed-current test. Set its planar wave scales to $(1,2,3,1)$, and all transverse and function scales to one. For the signed determinant $18$, `exact_scalar.py` contracts the four integral tensors, expands the resulting 760 current monomials and integrates both resistance simplices and the conductance simplex. Its normalized Wick mean and determinant-weighted mean are respectively

$$
\frac{4170228116417}{25758708393600},\qquad
\frac{4170228116417}{1431039355200}.
$$

`cubature.py` and `scalar166.cpp` independently give $750$ and $383$ modulo 1009 using 6188 resistance pairs and 42 conductance nodes. `integer_rule.py` gives $7^4$ times these two integrals as $35,42$ modulo 49. The ordered fixture therefore has three factors of seven in its denominator. `pole_cubature.py` gives the same leading pole on its separately specified degree range. These numbers are before positive-jet extraction, the full order sum, the symmetry factor $1/4$ and the scalar $\pi/158400$.

```sh
python3 exact_scalar.py
python3 frobenius.py
python3 gamma_fusion.py
c++ -O3 -std=c++17 resistance_check.cpp -o /tmp/twist-resistance-check
/tmp/twist-resistance-check
```

## All orders of the scalar fixture

For each path let $Q_j$ be its five noninitial prefix momenta and $\alpha$ the observed component powers. Removing the common Dirichlet factorial gives an integer path moment

$$
\alpha!\,[s^\alpha]\prod_{j=1}^5(1-s\cdot Q_j)^{-1}.
$$

`partition166.py` exports the coefficient of each first-current monomial in each of the two endpoint-indexed path matrices, retaining its wave order. `partition_sum.cpp` first sums its product with this path moment over the $120$ orders. In each term of the product of the two paths let $a,b,c,d$ count the choices $(w_B-1)\hat m_B,w_B\hat m_C,(w_C-1)\hat m_C,w_C\hat m_B$ from the first currents. After the integer contractions, the only factor to restore in this bin is

$$
\frac{(5!)^2}{(5+a+d)!\,(5+b+c)!}\,\mathbb E[(w_B-1)^a w_B^b(w_C-1)^cw_C^d].
$$

The default C++ integrator checks all intermediate 128-bit products and sums for overflow. Its final arithmetic uses exact rationals. The optional prime argument reduces the summed path tensors at that prime. The full Lie-path replacement reproduces the original integer current polynomial on twelve shuffled fixture orders. Both implementations give the same $14400$-order answer for the partition $01467\mid23589$. The first two orders on each path give the signed determinant-weighted four-term control $1602408298997/804959637300$ by both direct and grouped integration.

For $01467\mid23589$ the order-summed signed determinant-weighted coefficient is $91175657551/141968190$. Fixing the marked one-form on the second long path selects $126$ unordered partitions. Exactly $34$ have nonzero determinant for these waves. Summing their $34\cdot(5!)^2$ labeled terms from `scalar_partitions.tsv`, then dividing by the remaining reversal and restoring the scalar kernel, gives

$$
-\frac{3761204910364}{212952285}\frac{\pi}{2\cdot158400}
=-\frac{85481929781}{1533256452000}\pi.
$$

The coefficient of $\pi$ is $82$ modulo $1009$. Multiplying it by $49$ gives $1$ modulo seven. The seven-adic bound is attained by the complete scalar fixture. This evaluates the $166$ graph before projection to the mixed-current cohomology class.

`detector_factors.py` exactly reconstructs the $89$ integer detector weights as a sum of outer products under each of the $126$ path partitions, separately for the three function momenta. There are at most $8,4,3$ outer products for selectors $0,1,3$ respectively. `detector_node.py` accepts the planar branch and four scales, the six-bit positive-jet switch mask, and the partition order. It uses the determinant with the eight pure wave scalars factored out. On the ordered fixture the $89$-row value is $368$ modulo $1009$ from both the factored contraction and the sum of $89$ individually integrated rational expressions.

```sh
c++ -O3 -std=c++17 -I/opt/homebrew/opt/gmp/include -L/opt/homebrew/opt/gmp/lib partition_sum.cpp -lgmpxx -lgmp -o /tmp/twist-partition
python3 partition166.py | /tmp/twist-partition
python3 shape166.py /tmp/twist-partition scalar_partitions.tsv
python3 detector_factors.py
python3 detector_node.py /tmp/twist-partition --limit 1 --check
```

For the $89$-row calibration at branch zero, scales $(1,2,3,1)$ and all switches enabled, the same partition summed over $14400$ orders gives $758$ modulo $1009$. The exact single-order value is $3634438258175223115275162886058862305319133/77276125180800$. `summed_node` accepts another integer weight table and a field factory that returns its momenta, polarizations and residual determinant columns. Compile `partition_sum.cpp` with `-DARBITRARY_PRECISION` for unrestricted integral weights and pass `--prime 0` to `detector_node.py` for rational output.

## Short-edge factor

The first and last momentum on the length-one path $A$ is $\ell_A=w_A(\hat m_B+\hat m_C)$. When the cubic interaction sits at an endpoint, its two differentiated momenta give $k_B\wedge k_C=\pm k_A\wedge k_B$ by conservation there. When it sits at the interior one-form, exactly one bracket endpoint differentiates the length-one path. Each orientation thus contains $w_A$ before summation. Absorbing it shifts the Dirichlet parameters from $(11/2,1/2,1/2)$ to $(13/2,1/2,1/2)$ and gives

$$
\frac{\pi}{158400}\,\frac{11}{13}=\frac{\pi}{187200}.
$$

At moment grade $d_B+d_C\le11$ at most one of $(5+d_B)!,(5+d_C)!$ contains eleven. The shifted conductance moment has no factor eleven in its denominator. Hence the normalized ordered Schwinger mean is integral at eleven and its coefficient after restoring the source has at most one denominator factor eleven.

The full $89$-row rational calibration on $01467\mid23589$ at that scalar node is

$$
\frac{77256275400207643868076571812881391403352}{638856855}.
$$

## Single-geometry positive-jet projection

With the $94$ integer weights of `../mixed_single_weights.tsv`, take $02345\mid16789$. Only branch one with the marked transverse momentum switched on has a nonzero residual determinant. Eighteen of the twenty-four centered planar nodes survive. The $576$ surviving cases, together with the $5!^2$ path orders, give the projected normalized partition

$$
302\pmod{1009},\qquad \frac{302}{2\cdot158400}=66\pmod{1009}.
$$

The second number includes the remaining reversal and gives this partition's coefficient of $\pi$. At centered planar scales $(1,-1,-2,-3)$ and switch mask $63$, one order gives $60$ and its whole partition gives $450$ before the centered-difference and switch weights. Both values agree with the independent path contraction. `single_partition_02345.tsv` contains the normalized residues at each nonzero-determinant case.

To integrate a path without enumerating the inactive waves, average a zero-momentum vector at its own time. Its scalar insertion becomes $\pm u(\ell_P+\hat m_P)$. With $k$ zeros the remaining path has $5-k$ active waves, and the sum over orders gives $5!/(5-k)!$ times its normalized simplex mean. After inserting this factor the moment denominator of active degree $d$ is $(5-k+d)!/5!$. The native program expands the product in the initial current and the own mean. The full Wick polynomial has degree at most four in each of the two initial currents. Fix the active ordering and average over the zero insertions. Their common scalar product has leading initial-current degree $k$, so its active tensor has degree at most $4-k$. Each surviving monomial has initial-current degree plus explicit own-mean degree at most four on each path. Forty-one selected full detector nodes reproduced the uncut Wick result with this reduction.

The dense product of the path-order coefficient table and its moment table uses double precision on integer residues. Each dot product has at most $120$ nonnegative terms, and the build enforces $120(p-1)^2<2^{53}$. Thus every product and partial sum is an exactly representable integer. `single_projection.py` checkpoints completed cases. `exact_projection.py` sums certified residues until the modulus exceeds twice an absolute numerator bound and checks against the held-out modulus $1009$.

```sh
c++ -O3 -std=c++17 -DDIRECT_GLUE -I/opt/homebrew/opt/gmp/include -L/opt/homebrew/opt/gmp/lib native_partition.cpp -lgmpxx -lgmp -framework Accelerate -o /tmp/twist-native-partition
python3 single_projection.py /tmp/twist-native-partition --order 0234516789 --checkpoint single_partition_02345.tsv
python3 exact_projection.py --bound-only
```

The complete $02345\mid16789$ partition has now been reconstructed over the rationals. Its normalized coefficient is

$$
-\frac{106846132799324898264173673602557534332699721113304189295076620711006615636337335463232}{4732273}.
$$

After the reversal factor and the source kernel its contribution is

$$
-\frac{151770074999041048670701240912723770358948467490488905248688381691770760847070078783}{2129522850}\,\pi.
$$

The nineteen recorded prime residues have a 428-bit product, exceeding twice the proved numerator bound. The independent modulus $1009$ gives $302$ before and $66$ after the reversal and scalar-kernel factors. [Prime residues](single_exact_0234516789.tsv) · [Exact reconstruction](exact_projection.py)

## Planar projection of the complete path sum

Group the 126 partitions with the marked one-form on the second path by their subset of transverse labels $4,5,6,7,8$ on the first path. Each group is closed under every relabeling of the four planar vectors. The order sum, including the signed determinant, is alternating in these inputs. At zero transverse momentum of the marked function all four become ordinary planar vectors. The extracted Hamiltonians have degrees $3,3,4,4$. Under planar $SL_2$ the two exterior squares decompose as

$$
\Lambda^2\mathrm{Sym}^3=\mathrm{Sym}^4\oplus\mathbf1,
\qquad
\Lambda^2\mathrm{Sym}^4=\mathrm{Sym}^6\oplus\mathrm{Sym}^2 .
$$

They share no irreducible summand, so their product contains no invariant. This holds for each transverse switch mask on each partition group after planar coefficient extraction.

Write the planar scales as $(x_1,x_2,x_3,a)$ in the established order and factor out the eight pure external momentum scales in the determinant. The residual planar polynomial has degree seven. On the first planar branch its degrees in $(x_1,x_3)$ and $(x_2,a)$ are three and four. Indeed the two planar axes separately scale all three volume tensors, and the four planar polarizations occupy the two axes twice each. Before dividing by the pure determinant columns each planar coordinate has degree five. The removed columns have coordinate counts two and one. Exchanging the first and third pure waves negates the residual polynomial. Consequently its dependence on $(x_1,x_3)$ is a sum of multiples of $x_1^3-x_3^3$ and $x_1x_3^2-x_1^2x_3$. At $(1,\omega)$ for a primitive cubic root of unity the first vanishes and the second equals $\omega^2-\omega$. A four-point Fourier sum in $a$ extracts $a^3$ at $x_2=1$.

On the second branch reflect the planar coordinates by interchanging the two axes and swap the first two pure labels. The last planar direction changes sign. Reflection contributes $-1$ through the volume tensors, the normalized momentum columns and four planar polarizations supply a further $-1$ together, and the exchange of the two odd labels supplies $-1$. Thus the residual polynomial at $(x_2,x_1,x_3,-a)$ is minus its value at $(x_1,x_2,x_3,a)$. A three-point Fourier sum in $x_3$ and a four-point sum in $a$ isolate powers $2,3$ exactly. Their next nonnegative alias has planar degree at least eight. The remaining quadratic in $(x_1,x_2)$ is symmetric. At $(1,i)$, $i^2=-1$, its pure squares cancel and the mixed term has factor $i$.

Over $\mathbb F_{1009}$ these two rules take four and twelve samples per nonconstant transverse mask. Their coefficient is the same $[x_1x_2x_3^2a^3]$ as the centered rule. The planar weights remain $-1313280$ and $17280$. Exactly 16 of the 31 transverse label subsets require evaluations, with 31350 nonzero input nodes after determinant and rank tests. The old centered rule has 130328 nonzero-determinant nodes before rank tests. Individual partition contributions in the two rules can differ. The exact projection identities apply to each complete group.

`projective_shape166.py --check` compares the two coefficient functionals on every relevant degree-seven monomial after imposing its reflection symmetry. `projective_shape166.py SOLVER --all` stores individual finite-field evaluations in `single_shape_projective/`. The independent centered calculation `single_shape166.py SOLVER` stores the actual projected value of each partition in `single_shape_166/partitions.tsv`.

## Denominators of ordered integrals on all five one-form shapes

For integer external data the normalized Schwinger mean of every ordered one-form graph is integral at eleven. Put $m_P=n_P-1$ for the number of interior vertices on each path. The conductance Dirichlet parameters are $(13/2-n_P)_P$ with sum $13/2$. Modulo eleven these are $(-m_P)_P$ with sum one. The angular denominator is a unit through degree eight. Thus a conductance monomial integrates to zero modulo eleven whenever the exponent on some path $Q$ exceeds $m_Q$.

A resistance denominator at eleven requires mean degree $d_P\ge11-m_P$. It occurs at most once since total degree is at most eight and $m_P\le5$. In the coefficient of $\hat m_P$ each first current contributes $w_Q-\delta_{PQ}$ on path $Q$. This lies in the ideal generated by all $w_R$ with $R\ne P$, using $w_P-1=-\sum_{R\ne P} w_R$. Consequently every conductance monomial of this coefficient has total degree at least $d_P>10-m_P$ on the other paths after eliminating $w_P$. Their total capacity is exactly $10-m_P$, so angular averaging cancels the resistance pole. The scalar kernels contribute at most one factor of eleven to the physical denominator. The four long-path shapes attain this one-factor bound on the previously computed integer fixture. The stronger eleven-integrality proof for the physical $166$ kernel uses its short-edge factor.

A complete common denominator can also be obtained from the sharp degree bounds. At current degree eight the numerator vanishes on collinear independent currents. Its homogeneous part is therefore an integer sum of a $2\times2$ current minor times a degree-six monomial. To see integrality, substitute $C=tB$. Monomials with the same combined exponent and same number of $C$ factors have identical images. Exchanges $B_i C_j\leftrightarrow B_j C_i$ connect every pair in one such group. Substituting the means in a minor subtracts the two corresponding conductance products. For a distribution of $d_P$ mean factors among the $m_P$ nonzero prefix gaps with multiplicities $t_j$, the normalized resistance monomial gives $m_P!\prod_j t_j!/(m_P+d_P)!$. The integrator enumerates these factorial divisors and the exact angular moments, and takes the gcd over the three choices of dependent first current. For $166$ it also imposes degree four on each long initial current.

In shape order $166,256,346,355,445$, denominators after the scalar kernel, graph automorphism factor and division by $\pi$ divide respectively

$$
293772917487329280000,\quad
489621529145548800000,\quad
244810764572774400000,\quad
979243058291097600000,\quad
489621529145548800000 .
$$

These are universal divisors before any additional order and detector cancellations. [Exact rational enumeration](theta_denominators.py)

## Summing assignments of undifferentiated wave insertions

At a fixed finite-difference node let $Z$ contain the ordinary transverse waves with momentum switched off. Fix all active path assignments. Suppose $k$ waves of $Z$ occupy path $B$. By the zero-wave time integral, their combined scalar insertion in an orientation is

$$
(\pm w_B)^k(\pm w_C)^{|Z|-k}\prod_{j\in Z} M(u_j),
\qquad M=\hat m_B+\hat m_C .
$$

The marked-path position sign is included in the second sign. Each sign depends only on the oriented path. The two active simplex dimensions, the factors $5!/(5-k)!$ and $5!/(5-|Z|+k)!$, both endpoint tensors and every active prefix are fixed. Consequently every placement has the same integrated Wick tensor for each polarization.

For the $166$ graph the signed residual determinants now sum to one determinant in external label order. Keep $(q_j,0)$ for active waves assigned to $B$ and $(0,q_j)$ for those assigned to $C$. Insert $(q_j,q_j)$ for $j\in Z$, using the removed transverse momentum directions $q_j$ in these columns. Multilinearity expands this determinant over all assignments, and its two five-row blocks force exactly $k$ choices of $B$. The inactive columns span a diagonal copy of a three-dimensional space. Four such columns therefore give zero. The argument holds at every Schwinger point after integrating the insertion times of the inactive waves. For a theta shape with three occupied paths, replace the common diagonal embedding by the sum of the three loop-column embeddings with one formal variable per path. Extracting the coefficient with the specified inactive path counts gives the same rank-three vanishing, since the summed embedding is common to every inactive column. In the unfactored graph, setting an ordinary scale to zero after dividing by its pure momentum column extracts the linear Taylor field $u_j(q_j\cdot y)$. Thus each fully labeled shape has zero evaluation whenever four of these transverse inputs are linear fields. More generally, a linearly dependent set of their momentum directions already forces this cancellation. The argument is pointwise in the active Schwinger data.

| Number of inactive transverse waves | Original nonzero evaluations | Evaluations after summing assignments |
|---|---:|---:|
| 0 | 1083 | 1083 |
| 1 | 5274 | 5274 |
| 2 | 10263 | 7641 |
| 3 | 9946 | 3640 |
| 4 | 4784 | 0 |
| Total | 31350 | 17638 |

The table is for the single-momentum detector modulo $1009$. The checker compares 325 sums against direct ten-dimensional determinants and 12 full path-order evaluations across interchangeable assignments. On the complete planar orbit with transverse labels $4,5$ on $B$, the normalized sum is $677$ modulo $1009$. The full-shape driver reuses projective evaluations and verifies their Wick values whenever multiple earlier evaluations occur in the same group.

```sh
python3 zero_shape166.py --count
python3 zero_shape166.py SOLVER --verify
python3 zero_shape166.py SOLVER
```

Three independent transverse linear fields can survive. On planar branch zero set the four planar scales to $(1,1,2,1)$, choose transverse polarization mask zero, and take the test one-form along the second transverse coordinate. Use ordinary waves at transverse labels $4,5$, the marked function momentum $(1,0,-1)$ and the linear fields $u_j(q_j\cdot y)$ at labels $6,7,8$. The five transverse directions, in label order, are those of the single-momentum detector. With these data, 44 of the 126 signed residual determinants are nonzero and 24 groups survive their assignment sums. Restoring the removed planar factors, the path reversal and the scalar kernel gives the complete $166$ coefficient

$$
\frac{2}{2}\frac{\pi}{158400}\frac{58316291}{1277713710}
=\frac{5301481\pi}{18399077424000}.
$$

The normalized numerator bound at common denominator $765033639289920000$ is $568481908240729630243160064000000000000$. The reconstruction modulus is $46657835164917879513235769216268016187723$. A separate evaluation of all 44 nonzero determinants gives the normalized sum $994$ modulo $1009$. This fixture measures the graph on linear external jets. [Exact all-order replay](linear_shape166.py)

For arbitrary linear transverse vectors, let $A_i{}^a{}_b$ be the traceless linear map of input $i$, with $a,b$ belonging to the transverse vector space. At each active Schwinger point the entire dependence on the chosen $k$ linear fields factors through the exterior product of the covectors $M_a A_i{}^a{}_b$. For a rank-one field this is precisely $M(u_i)q_i$. Multilinearity extends it to every traceless matrix, since the off-diagonal matrix units and the differences of diagonal units are spanned by traceless rank-one matrices. The exterior product takes values in symmetric degree $k$ in $M$ and exterior degree $k$ in the lower slots. Exact coordinate ranks are

| Linear inputs | Exterior domain dimension | Rank of the universal tensor | Kernel dimension |
|---|---:|---:|---:|
| 1 | 8 | 8 | 0 |
| 2 | 28 | 18 | 10 |
| 3 | 56 | 10 | 46 |
| 4 | 70 | 0 | 70 |

For three inputs the target is $\mathrm{Sym}^3\mathbb C^3\otimes\Lambda^3(\mathbb C^3)^*$, of dimension ten. The exact graph value above shows that the $166$ shape uses this component. [Exact ranks in a traceless matrix basis](linear_transverse.py)

In particular the invariant transverse Cartan three-cycle pairs to zero with each shape, for arbitrary remaining inputs on this restriction. Raising the three slots of the closed form $(A,B,C)\mapsto\mathrm{tr}(A[B,C])$ by the trace pairing gives a cycle with six nonzero exterior coefficients in the traceless basis. Its pairing with the same closed form is $-8$. The replay computes both Chevalley–Eilenberg differentials and its zero image in the ten-dimensional graph tensor.

The same factorization applies to traceless linear vector fields on the entire holomorphic space $V=\mathbb C^5$. For $2\le k\le5$ every target monomial in $\mathrm{Sym}^k V\otimes\Lambda^k V^*$ is attained by a wedge of traceless matrices. To construct it, match the multiset of upper labels to the distinct lower labels while avoiding equal pairs. A set of upper labels containing two distinct values admits every lower label as a neighbor. A set containing copies of only one value admits at least $k-1$ neighbors. Thus the only failed matching has all $k$ upper labels equal to one of the lower labels $i$. In that case use $E_{ii}-E_{jj}$ at the lower label $i$, with $j$ any other lower label. The $E_{jj}$ term repeats that lower label in the exterior product. All other factors are off-diagonal matrix units. The replay constructs all $150,350,350,126$ target monomials at $k=2,3,4,5$ individually.

In particular the graph kills every invariant exterior chain of positive degree of $\mathfrak{sl}_5$ among its linear inputs. For $k=1$ the traceless adjoint has no invariant vector. For $2\le k\le5$, an invariant tensor would define an $SL_5$ map $\Lambda^k V\to\mathrm{Sym}^k V$. The diagonal torus forces the image of the wedge of the first $k$ basis vectors to be a multiple of their symmetric product. Adding the first vector to the second fixes the wedge and changes that product, so the multiple is zero. Permuting the coordinate lines gives zero on every exterior basis vector. These identities hold shape by shape with the remaining wave data fixed.

Momentum conservation removes the entire fifth-linear-input component of the actual graph. Work over the fraction field of the three path variables. Project the two-dimensional loop space onto its quotient by the common line used for inactive columns. The inactive columns project to zero. Each active column projects to a scalar multiple of its own physical momentum. A nonzero ten-dimensional determinant requires these projected columns to span five dimensions. With five or more linear fields there are at most five remaining momenta, including the one-form momentum, and their sum is zero. Hence every labeled shape vanishes on five linear vector fields even when their directions span all of $V$. In the single-momentum calculation, all $6492$ rank-four native nodes among the original $31350$ cancel through this rank condition.

Four independent linear fields attain the ambient bound. On the preceding three-linear-field fixture, also replace the third pure planar wave by its linear field. The four planar scales in the factored coefficient are now $(1,1,0,1)$. All other data remain as in that fixture. The shape sum from the 24 surviving assignment groups gives

$$
\frac{\pi}{2\cdot158400}\left(-\frac{2915}{4483206}\right)
=-\frac{53\pi}{25823266560}.
$$

The same reconstruction modulus exceeds twice the numerator bound $87069075336229823623987200000000000000$. Direct visits to all 44 nonzero determinant partitions give the normalized sum $687$ modulo $1009$. [Exact replay with the four-linear fixture selected by the --four option](linear_shape166.py)

## Integrated ordered detector on all five shapes

Use the 94-row single-momentum detector and the indicated vertex orders. For each order the three-path calculation includes all four orientations, all 3072 centered planar and transverse cases, all resistance integrals and the conductance integral. The entries use the normalized Schwinger mean, including the signed residual determinant and positive-jet weights. The scalar kernel and the sum over other orders enter the full graph.

| Shape | Vertex order | Signed scalar fixture over the rationals | Full detector modulo 1009 |
|---|---|---:|---:|
| $166$ | `0146723589` | $4170228116417/1431039355200$ | $920$ |
| $256$ | `0146723589` | $-1360113875759/1192532796000$ | $10$ |
| $346$ | `0146723589` | $26094891211/56787276000$ | $594$ |
| $355$ | `6901472358` | $-5763274313/13449618000$ | $573$ |
| $445$ | `0691472358` | $11269210233/37858184000$ | $296$ |

The scalar column uses the existing `fixture()` in `exact_scalar.py`. Independent modular cubature applied to the primitive Wick matrices agrees with all five rational integrals. It uses 36 conductance points and the tensor products of simplex rules through total resistance degree eight. On the full 94-row input at centered case 35, direct rational summation gives signed residues $546,430,411,704,120$, in shape order. [Rational three-path integral](ordered_shapes.py) · [Primitive cubature](check_ordered_shapes.py) · [All-polarization rational control](ordered_node_check.py) · [Native contraction](ordered_detector.cpp) · [Projection and saved-node replay](ordered_detector.py)

## The complete five-shape sum on one linear vector field

Choose one ordinary vector field with polarization $u$ and linear dependence $q\cdot y$. After factoring its momentum column, set its external momentum to zero. Fix the paths, order and Schwinger data of the other nine waves. Write $m_P$ for their number on path $P$. If the linear field is inserted on $P$, the three original edge counts are $n_R=1+m_R+\delta_{RP}$. Integrating its time and summing the positions multiplies its orientation sign by $(m_P+1)w_P M(u)$. The other path matrices and both endpoint tensors use only the active prefixes and means.

The Dirichlet density had parameter $13/2-n_R$ on path $R$. Inserting $w_P$ changes this to $11/2-m_R$ on every path and multiplies its normalized mean by $(13-2n_P)/13$. The scalar kernel divided by $\pi$ has the five reciprocal denominators in the source shape table. For fixed $(m_A,m_B,m_C)$, its product with
$$
\frac{(m_P+1)(13-2n_P)}{13}
$$
is independent of the chosen feasible $P$. Equivalently the scalar kernel contains the shape dependence $\prod_R\Gamma(13/2-n_R)/(n_R-1)!$. The position count cancels the extra factorial and the Dirichlet shift advances the gamma function. A parent path with seven edges has six columns in one five-dimensional fiber, so its determinant is zero; its formal placement can be included with coefficient zero. This uses only the regular shifted density.

Use ordered paths and put the common factor $1/12$ in front of their sum. Each underlying shape has $6/\prod_r m_r!$ orderings, where here $m_r$ counts paths of length $r$, and its automorphism group has size $2\prod_r m_r!$. In the source determinant the three path columns of a wave $q$ are $(q,0),(q,q),(0,q)$. Interchanging paths preserves the combined signed quantity: a bracket endpoint acquires the sign of the permutation of its two vector legs, and the two alternating factors of a cubic endpoint compensate. The loop-basis determinant sign, the positions of the odd contractions and the one-form position sign are included in the replay.

At each position the zero-wave scalar from an ordinary Lie vertex has sign $-$ in the forward direction and $+$ in reverse. Multiply by the sign of the permutation of the thirteen edge tails and by $(-1)^{9-g}$ for the zero-based one-form slot $g$. For each fixed active order and orientation, the resulting three signs for inserting the wave on $A,B,C$ are $(+,-,+)$ times one common sign, separately in every insertion gap. There are only 25 ordered triples of nonnegative active path counts summing to nine with each at most five, nine possible marked slots and four orientations. The exhaustive integer replay verifies all 10,800 gap signs, all 60 feasible parent weights and all 1,200 path permutation signs.

The three determinants consequently sum by multilinearity in the single column
$$
(q,0)-(q,q)+(0,q)=0.
$$
Thus the full symmetry-weighted order-$g$ graph vanishes whenever any ordinary external vector field is linear and traceless. Constants already vanish by their momentum column. This is a statement about the integrated graph itself before addition of a counterterm. On the four-linear fixture the other four shapes sum to $+53\pi/25823266560$, exactly canceling $166$. On the three-linear transverse fixture they sum to $-5301481\pi/18399077424000$.

The separate cubic-word replay integrates the insertion time at a fixed active geometry by summing the actual gap widths, divides each path value by its $w_P$, and evaluates the original two Wick vertices. Of 300 fixed geometries, 292 individual path choices are nonzero; every three-path sum is zero modulo $1009$. [Exhaustive scalar and Koszul identities](horizontal_theta.py) · [Independent cubic-word contraction](horizontal_words.cpp)

For the complete positive-jet projection on the mixed detector, every mask with an ordinary transverse switch off therefore drops out after adding the five shapes. The term with constant marked function drops out by the planar exterior-square calculation. Only the all-active mask $63$ is required for the total class coefficient. The isolated $166$-shape calculation retains its full switch sum.

For $g^m$ with $1\le m\le5$ and no external scalar waves, mark an arbitrary set of $2m-1$ external one-forms. Each mark reverses the arrow on its path. The sign for external labels in the fixed-label determinant is the parity of the full label permutation times the parity of its ordinary-vector subsequence. The scalar normalization and the signed column identity above still apply. The exhaustive checks by coupling order are

| $m$ | Insertion-gap signs | Path-permutation signs | Nonzero individual path choices in 300 independent cubic-word tests |
|---|---:|---:|---:|
| 1 | 10800 | 1200 | 292 |
| 2 | 100800 | 14400 | 296 |
| 3 | 151200 | 30240 | 273 |
| 4 | 43200 | 14400 | 307 |
| 5 | 1200 | 1200 | 270 |

Every three-path sum vanishes. Thus at zero external $\nu$, **all five fully summed two-loop coupling sectors vanish on any constant or traceless linear ordinary vector field**. The rightmost column uses the original cubic vertices and exact finite-field arithmetic. This cancellation uses no cohomology closure assumption.

On a constant scalar background the same zeros are multiplied by $(1-\nu)^{-(12-m)}$, so each coefficient in the constant $\nu$ expansion also vanishes. [Constant-wave combinatorics](../../constant_nu.md)

The independent three-path moment contractions reproduce the linear-insertion cancellation after all possible positions of one specified linear field. In particular, after extracting $\pi$ the three ordered path subtotals in rational fixture 12 are $24514148863/85003737698880000$, $24514148863/10625467212360000$, and $-24514148863/9444859744320000$. Their sum is zero. Two further rational controls compare different parent shapes. Tests modulo 1009 cover all five shapes with 120 active-word fixtures, 70 nontrivial cancellations and 60 checks under path relabeling. [Direct integrated replay](horizontal_integrals.py)

## Two geometries for the complete positive-jet detector

For the sum of all five shapes, horizontality at each pure vector eliminates every transverse switch with one of the five ordinary transverse momenta set to zero after its determinant column is factored out. The planar four-cycle eliminates the constant marked function, as proved above. Only switch mask 63 remains. Let $F(a)$ denote the 94-row residual amplitude summed over all five shapes at the second planar branch, scales $(-38,-38,29,a)$ and this mask, with the same scalar and graph weights in every evaluation. The original class detector on this sum is exactly

$$
-\frac{80}{10469}\bigl(F(2)-2F(1)\bigr).
$$

Indeed the degree-seven residual polynomial is divisible by $x_1x_2x_3$ by the three planar linear Ward identities. After extracting $a^3$ the remaining factor beyond $x_1x_2x_3$ is linear in the three pure scales. Its three terms place the ordinary quartic Hamiltonian on the first, second, and third wave. On the two invariant transvectants their ordered vectors, including the signs for moving that quartic to the third position, are $( -1,-78),(-1,-78),(16,-48)$ in units of $20736$. The corresponding weighted vector is $(540,4536)$, which becomes the required $(-50,-420)$ on multiplication by $-5/54$. Reflecting the planar axes and exchanging the first two equal scales sends $a$ to $-a$ and negates the residual polynomial. Its degree in $a$ is at most four, so $[a^3]F=(F(2)-2F(1))/6$. Restoring $20736$, the factor $x_1x_2x_3=(-38)^2\cdot29$ and $-5/54$ gives the displayed coefficient. The rational identity was checked on all 174 products of seven planar volume brackets with the required order on every pure wave. It uses the five-shape Ward identity. [Exact projection and 166 partition replay](total_projection.py)

Summing all $126$ gamma-on-second-path partitions at each of the two geometries, the $166$ order sums are $251,880$ modulo $1009$. The symmetry-weighted coefficients of $\pi$ are $262,199$. The contribution of this shape to the two-geometry expression for the complete class detector is $231$ modulo $1009$. [All 252 order sums](total_projection_166.tsv)

## Complete partitions of an arbitrary one-form theta shape

For the two all-active scalar geometries, fix the labeled path partition. On path $C$, sum the coefficient of each initial-current monomial times the complete normalized path moment over every order, separately by matrix element and external polarizations. Contract the 94 polarization weights before joining the other two paths. For a monomial in the $B,C$ initial currents distribute each factor into the three resistance means using $\ell_P=w_P(\hat m_A+\hat m_B+\hat m_C)-\hat m_P$. These choices fix the conductance moment and the moment multi-index on each path. The $A$ order fixes its moment. For each remaining monomial the program multiplies the whole array of $C$ responses by the array of $B$ order moments, then contracts against the two endpoint tensors and remaining path fibers.

The total initial-current degree bound is eight after the four Wick orientations are summed. For each retained monomial the allowed $C$ observer multi-indices are subindices of its total five-component exponent. At degree eight there are at most $3^3 2^2=108$ subindices. Each integer dot product in the two dense matrix operations therefore has length at most $120$ and residues at most $1008$. Its exact integer value is less than $2^{53}$. The reversal representative fixes the order of two distinguishable labels on one active path. The routine restores its factor of two before applying the graph weight. Path-permutation representatives leave a common division by two in the outer sum. The native program subdivides its contraction tables by fiber element and the driver reaps each native process after twelve completed partition entries.

On the first nonzero partition at the first class geometry the five complete normalized sums in shape order are $(312,967,120,271,707)$ modulo $1009$. Keeping at most two orders on each path gives $(39,551,640,884,926)$ by independent evaluations of each ordered graph. All 72 partitions of each $166$ computation give $262,199$ after the scalar kernel and graph symmetry factors. [General path-order integrator](general_partition.cpp) · [Partition driver](complete_nodes.py) · [Class geometries](total_projection.py)

```sh
c++ -O3 -std=c++17 -I/opt/homebrew/opt/gmp/include -L/opt/homebrew/opt/gmp/lib general_partition.cpp -lgmpxx -lgmp -framework Accelerate -o /tmp/twist-general-bounded
python3 complete_nodes.py /tmp/twist-general-bounded 256 1
```

For each path-A order, the native kernel now multiplies the two endpoint tensors through path A first. This gives one polynomial from each path-B matrix entry to each path-C matrix entry. It then accumulates the coefficients of every path-B order on each common monomial before joining the path-C moment. The third-path contraction is shared by the two flows with the same arrow direction, with the original flow sign restored. Polynomial entries and monomial indices are stored contiguously. This ordering reproduced the earlier kernel on 50 partially summed cases across the two detectors and the five full-order residues $(478,417,240,240,337)$ at partition index 21, scale two, in shape order. The full two-geometry replay validates every determinant-supported partition key before forming the class coefficient.

```sh
python3 total_projection.py --complete-shapes complete_shapes
```

## One planar wave cycle for the complete four-plus-five class

Use three unmarked planar Hamiltonian exponentials of momenta $(1,0),(0,1),(1,1)$ and one marked Hamiltonian exponential of momentum $(-1,-1)$. Their Hamiltonian vector polarizations are $(-p_2,p_1)$. Keep the original transverse momenta, marked function and 94 polarization weights. On the scalar planar relative complex evaluate each ordered Hamiltonian monomial by its exponential coefficient, alternating the first three slots. On all six invariant scalar three-cochains the values of their differentials are zero. On the closed scalar four-cochain $P$ from the original lift the value is $32/9$. Its value on the original integral cycle is $2365440$. Since the regular scalar four-class is one-dimensional, the original cycle and $665280$ times the new wave cycle are homologous. These six boundary tests and the nonzero pairing are exact rational equalities. [Planar evaluation and native input](single_planar_nodes.py)

The complete graph is horizontal on each ordinary planar and transverse field. At constant marked function, all four planar fields are ordinary and only positive planar jets can contribute. For a transverse output the planar coefficient is scalar, of Euler degree six. The unmarked positive-jet scalar multiplicities at this degree are one each in arities two and three and zero in arity four. Hence the constant marked function drops out of this transverse component on every planar Taylor profile. Evaluation of the two cycles on the complete closed graph gives the same projected coefficient. The normalization relating them requires the sum over shapes.

The closure used here follows from the BV consistency identity. At a fixed heat scale $L$, let $I[L]$ be the effective interaction and let $A[L]$ be its QME defect. The square of $Q_{\rm BV}+\hbar\Delta_L$ is zero. Conjugation by $e^{I[L]/\hbar}$ and application of the square to the unit give
$$
0=(Q_{\rm BV}+\hbar\Delta_L+\{I[L],-\}_L)A[L]+\hbar^{-1}A[L]^2
=(Q_{\rm BV}+\hbar\Delta_L+\{I[L],-\}_L)A[L].
$$
The last equality uses the odd parity of the defect. On the local UV collision faces the same identity is the sum over an inner defect and its outer insertion, with additive loop orders. At two loops the three choices of outer and inner loop orders are $(2,0),(1,1),(0,2)$. The tree defect vanishes by the classical master equation. The one-loop inner weight vanishes since its real heat factor is the Pfaffian of an odd skew matrix. The remaining term is the tree differential on the full two-loop defect. At $g=0$ all twelve vertices supply just twelve incoming $\gamma$ or $\beta$ letters, while thirteen propagators require thirteen. Thus the coefficient of $g$ in the consistency identity is $Q_{g=0}X_1=0$. Projection to the divergence-free, zero-$\nu$ quotient commutes with the strict differential. A common nonzero scalar multiplying the two-loop weights leaves this identity valid. [Paper, proof of the QME lemma](../../../../../sources/LLC_main.tex) · [Two-loop letter count and normalization](../../../VERIFY.md)

Reflection exchanging the two planar axes swaps ordinary labels zero and one and fixes the other momenta. Each of the four Hamiltonian polarizations gains the determinant of this reflection. The ten-dimensional momentum determinant gains the square of the frame determinant and a minus sign from permuting its columns. The unique five-dimensional volume tensor supplies a second minus sign. Equivalently, equivariance in planar determinant weight $-3$ and alternation of the two ordinary labels give two minus signs. The complete order sum on each partition therefore equals that on its reflected partition, including the relabeling of equal-length paths. The five supported partition counts are $(50,500,1000,886,1272)$; their reflection orbit counts are $(25,250,500,443,636)$. The two members of the first orbit of each shape independently gave $(827,383,829,772,513)$ modulo 1009 with the factor $665280$ included and before shape normalization.

```sh
python3 single_planar_nodes.py 166 --verify-planar --solver /tmp/twist-general-bounded
```

In the connected cutoff expansion the unrooted symmetry factor also follows directly. For a specified heat edge the bracket $\tfrac12\{I,I\}$ supplies two unordered endpoint vertices; the ten other interaction insertions carry $1/10!$. The ordinary Wick contractions come from the exponential of the propagator contraction, so their edge orderings cancel their exponential factorial. Orbit–stabilizer gives $1/|\mathrm{Aut}(\Gamma,e)|$ for the graph with its heat edge $e$. Summing over heat-edge orbits equals $1/|\mathrm{Aut}(\Gamma)|$ times the sum over all heat-edge choices on a representative. That edge sum already belongs to the UV boundary kernel. Equivalently the $12!$ in the twelve-entry labeled bracket cancels the QME Taylor coefficient $1/12!$. To pass to the integrated density, use any one vertex as the simultaneous-translation coordinate; subtracting its position from all others has unit absolute Jacobian. Choosing the representative with root $N'$ adds no multiplicity.

The edge-form reordering can be checked separately. An ordinary propagator form has degree five, a heat form degree six, and each edge contraction is odd. With the heat edge in position $j$ of the thirteen ordered edges, moving all contractions rightward gives the sign $(-1)^{78-(j-1)}$. The exterior derivative of the ordered product of thirteen propagator forms gives the cut sign $(-1)^{j-1}$. Their ratio is $(-1)^{78}=1$, independent of the cut. This fixes the edge interchange in the common kernel convention used above.

### Positive-resistance bound for the complete current coefficients

Take the real momenta and polarizations of either integral detector. Drive unit current between the two endpoints of the chosen edge. The potential attains its extrema at these vertices, and their effective resistance is at most the edge resistance. By reciprocity the edge transfer coefficients are those vertex potentials divided by the edge resistance. Subtracting their common minimum leaves coefficients in $[0,1]$. By total momentum conservation, the current on any edge belongs to the convex hull of the subset sums of the external momenta. Thus its Euclidean norm is bounded by the largest norm $R$ of the $2^{10}$ subset sums, and its evaluation on an ordinary polarization $u$ is bounded by $\tfrac12\sum_{i\ {\rm external}}|p_i(u)|$.

For an ordinary wave with momentum $p$ and $p(u)=0$, the transverse transport matrix is $\alpha\,1\pm u p^{\mathsf T}$ or its transpose, where $\alpha$ is the entering edge current evaluated on $u$. Restricting its singular-value calculation to the plane of the orthogonal vectors $u,p$ gives
$$
\|\alpha\,1\pm u p^{\mathsf T}\|_2
=\tfrac12\bigl(\sqrt{4\alpha^2+\|p\|^2\|u\|^2}+\|p\|\,\|u\|\bigr).
$$
A bracket endpoint has Frobenius norm at most $\sqrt{8}R$ and a cubic endpoint has norm at most $\sqrt{6}R^2$, by contracting their two Kronecker terms or their alternating three-index tensor. A one-form source has spectral norm $\|\eta\wedge p_\gamma\|$. A one-form sink is the dual of the wedge of this two-form with an edge current, of norm at most $R\|\eta\wedge p_\gamma\|$. Contract the two endpoints against the tensor product of the three transports. The two cubic-endpoint flows and two interior-cubic flows together have absolute value bounded by
$$
2(\sqrt{48}+8)R^3\|\eta\wedge p_\gamma\|
\prod_{i\ {\rm ordinary}}\|\alpha_i\,1\pm u_i p_i^{\mathsf T}\|_2 ,
$$
using the scalar bounds for $|\alpha_i|$ in the product. All resistance and conductance measures are positive normalized averages. Every radical in the replay is rounded upward using integer square roots. [Exact arithmetic bound](theta_spectral_bound.py)

On the single planar four-cycle the first three complete coefficients divided by $\pi$ are $125,907,415$ modulo 1009 in shape order. The first coefficient is
$$
-\frac{38829246720120915439256439531003532212060110643364504843308360883391655159085332721}{21729825}.
$$
The independent two-geometry calculation at $256$ gives $895,943$ before forming its contribution to the complete class. In the full labelled sum the graph automorphism acts freely on the external orders. Its symmetry divisor cancels on each complete orbit. Using this cancellation for all five shapes, the total coefficient divided by $\pi$, multiplied by $1103944645440000$, is an integer multiple of $114415887442730778653770948594609901093939178844554683860158369901482612$. The absolute value of the remaining integer is at most $27360787075254838090685647461259407$. [Shape and total reconstruction](single_planar_crt.py) · [Independent order sums](complete_shapes)

For the three-path contractions, the two initial five-component currents have total degree at most eight, hence exactly $\binom{18}{10}=43758$ possible monomials including degree zero. For each exponent in the first current the allowed exponents in the second form an initial segment of the degree-ordered table. The integrator groups the path weights by this exact integer index and then sums repeated third-path vectors separately within each monomial. Blocks split before contributing their result when they reach the storage bound. Agreement with the earlier contraction covered all eight complete benchmark partitions at 1009 and 6000011, and 25 independently polarized capped-order checks across all five shapes. [Three-path contraction](general_partition.cpp)

### Three small transverse geometries for the four-plus-five coefficient

Use the same planar cycle with five transverse momenta. The first four are $(e_1,e_2,e_3,-e_1-e_2-e_3)$. In the three geometries the fifth momentum and the marked function momentum are respectively
$$
(-e_1-e_2+e_3,e_1),\qquad
(-e_1-e_2+e_3,e_1-e_3),\qquad
(e_1+e_2-e_3,e_2-e_3).
$$
The 201 integer polarization weights have maximum absolute value 179 and sum of absolute values 4460. Weighted evaluation annihilates all 132 incoming primitives. On all 89 independent closed transverse currents the ratio of the original single-momentum test to this test is
$$
-\frac{1830654199083692458460335177513758417503026861512874941762533918423721792}{3557}.
$$
At least one closed current has a nonzero pairing. The exact check evaluates the two cycles in the same 138-coordinate basis. [Weights](four_small_weights.tsv) · [Integer differential evaluations](four_small_cycle.cpp) · [All 89 rational comparisons and spectral bound](four_small_bound.py)

Write $A$ for the original complete coefficient divided by $\pi$ and $r=114415887442730778653770948594609901093939178844554683860158369901482612$ for its preceding divisibility factor. The ratio just computed is $-16r/3557$. Both integer cycles have denominator dividing $1103944645440000$. Since $\gcd(3557,16)=1$, applying the second cycle proves
$$
\frac{1103944645440000}{16r}\,A\in\mathbb Z,\qquad
\left|\frac{1103944645440000}{16r}\,A\right|
\le460823791161032393567293574956 .
$$
The inequality is the positive-resistance bound on the three displayed geometries. The five prime fields $1009,6000011,6000023,6000041,6000047$ have product $1307690589355868259432719918779$, exceeding twice this bound. [Reconstruction](single_planar_crt.py)

The $355$ shape gives $158$ modulo $1009$ after integrating all 443 reflection orbits and dividing by twice the scalar denominator $798336$. The first four shapes therefore sum to $596$ in this field. The path-B fiber polynomial depends only on the B current. At fixed C exponent, a degree-eight B monomial has at most $3^3 2^2=108$ factorizations across the endpoint link and the B fiber. Across 25 matrix entries the direct accumulation uses at most 2700 products before reduction. The native contractor accumulates these products by their integer monomial addresses and stores the nonzero coefficients of each path order. It reproduced the 16 full-partition benchmark values at primes $1009,6000011$ and 25 independently polarized capped-order values. [Partition ledger](single_planar_1009/355.tsv) · [Contraction](general_partition.cpp)
