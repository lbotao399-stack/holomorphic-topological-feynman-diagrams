# Boundary fan with a single scalar spectator

Apply the paper's $Q_n$ at $n\ge2$ to the undifferentiated gauge-invariant word $\operatorname{tr}(\beta_1\beta_2\beta_3 b^{n-2}\gamma^I)$. Only $\gamma^I$ can be the spectator. The three $\beta$ derivatives supply $-\theta_1\theta_2\theta_3$ together in the paper's ordering, each $b$ derivative supplies $-ig$, and the two output superfields supply differentiated ghosts. Contracting the scalar instead needs at least four Grassmann coordinates. All $(n+1)!$ assignments of the remaining slots have the same Grassmann factor. The spinor contraction of the two ghost jets has exchange parity $(-1)^{n+1}$. [Component differential](../../sources/LLC_main.tex#L4788) · [Two-loop normalization](NOTE.md#4-loop-orders-drv)

Take traceless external color matrices $D,E$ on the ghost jets and $F$ on the scalar spectator. The full color dependence after factoring the common generator-pairing and commutator conventions is a polynomial in $N$ multiplying $\operatorname{tr}(DEF)+(-1)^{n+1}\operatorname{tr}(DFE)$. Transposition exchanges the cyclic orders and reverses the sign of each of the $n+1$ commutators. A disconnected external trace contains a single traceless matrix and vanishes. $U(N)$ completeness gives the same contractions for $SU(N)$ because a central generator in any spoke has zero commutator. Set $m=n+1$, the number of spokes. The coefficient of $N^{n-2h}$ has ribbon genus $h$. The planar coefficient is one. On this same physical cochain the exact ratio of the one-handle term to the planar term is
$$
\frac{\binom{n+2}{4}+n(n-1)+\binom n3+2+(-1)^n(n-1)}{N^2}.
$$
At $n=2,3,4,5$ the complete color polynomials are
$$
N^2+6,\quad N^3+12N,\quad N^4+36N^2+120,\quad N^5+63N^3+360N .
$$
For odd $n$ the symmetric color trace requires a nonzero cubic invariant. For even $n$ the antisymmetric trace tests the commutator. Both are available for $SU(N)$ with $N\ge3$. These are matrix elements on cochains. The holographic comparison to one-loop BCOV on the punctured fivefold is open. [Color replay](code/fan_trace_handles.py) · [Two-loop direct color check](code/two_loop_trace.py)

With the paper's $\operatorname{tr}(t_A t_B)=\tau\delta_{AB}$, the full coefficient of the two cyclic traces times the color polynomial and the spectator field is
$$
\frac{i^{n-2}g^n\tau^n}{2(32\pi^2)^n(n!)^3(n+1)} .
$$
The two ghost jets carry $n$ lower and $n$ raised dotted derivative indices contracted in order. The chain of $n+1$ structure constants converts to its nested commutator with factor $i^{-(n+1)}/\tau$, and the $n+1$ spoke contractions supply $\tau^{n+1}$. The Grassmann and ghost factors give $g^2(-ig)^{n-2}$. Multiplying these by the paper's fan coefficient gives the phase above. [Undifferentiated input formula](../../sources/LLC_main.tex#L298)

## Isolating a ribbon genus of the boundary supercharge by derivatives on the input

Keep the cochain $\operatorname{tr}(\beta_1\beta_2\beta_3 b^{n-2}\gamma^I)$ for $n\ge2$, with the $m=n+1$ active letters numbered from 1 in trace order. The scalar spectator is last. Holomorphic derivatives preserve the Grassmann count, so the same two differentiated ghost endpoints and constant spectator give the entire projected output. Dotted spinor derivatives always have the conventions of the paper. Fix their two outgoing holomorphic momenta $p,q$. On the projected fan, a derivative on the letter attached at spoke $i$ measures the momentum $a_i p+b_i q$. Here $a_i,b_i$ are the endpoint-current fractions, with joint boundary distribution uniform on the part of two unit simplices where $a_i/b_i$ decreases in path order. The topology at fixed $n$ is measured by the power $N^{n-2h}$, with $h$ the number of handles. [Endpoint chamber and inverse resistance calculation](../iib-bcov/README.md#exact-two-endpoint-kernel-in-displacement-coordinates) · [Supercharge on derivative insertions](../../sources/LLC_main.tex#L3991)

For any two tangent directions $u,v$, differentiate the active slots by
$$
\sum_{i=1}^{m}\partial_{i,u}\partial_{i,v}
-\frac1m\left(\sum_{i=1}^{m}\partial_{i,u}\right)
           \left(\sum_{j=1}^{m}\partial_{j,v}\right).
$$
In every color assignment the polynomial is symmetric in the spokes, so its mean equals the uniform product-simplex mean. It multiplies the original two-ghost amplitude by $n[p(u)p(v)+q(u)q(v)]/[(n+1)(n+2)]$. The mixed coefficient cancels for every $N$. On the fixed flat fivefold scalar profile the same momentum test instead gives $17/525$ times the original heat coefficient at two loops when $p(u)=q(v)=1$, $p(v)=q(u)=0$. This compares prescribed tests on the two complexes; an identification of observables under the D3 transfer remains open.

**An exact two-loop cochain with only one handle.** Let $D_i$ differentiate only beta slot $i$ along a fixed tangent vector $u$ before coincidence. Apply
$$
2D_2^2-2D_1D_2-2D_2D_3+D_1D_3
$$
to the four-letter input. In the two-ghost amplitude the complete factor $N^2+6$ of the undifferentiated trace is replaced by
$$
\frac{p(u)^2+q(u)^2}{2}-\frac43p(u)q(u).
$$
Both planar cyclic traces cancel as polynomials in $p,q$. In path order $123,132,213,231,312,321$, the coefficients in columns $p(u)^2,p(u)q(u),q(u)^2$ are respectively $(0,0,0),(0,-1/3,1/4),(1/4,-1/3,0),(1/4,-1/3,0),(0,-1/3,1/4),(0,0,0)$. Only orders $123,321$ contribute a planar trace. Each order has one handle with color $\operatorname{tr}([D,E]F)$, where $D,E,F$ denote the three external color matrices. At $p(u)=1,q(u)=0$ the original three Levi-Civita tensors and four Pauli matrices give $-1/16$ for one ordered ghost component with color indices $(1,2,3)$ in the paper's $\tau=1/2$ convention. The other ghost assignment restores the path reversal factor when extracting the coefficient of two independent odd waves. This probes the cochain differential. For example take $I=2$, set $\beta_2$ and its jets to zero, keep $\beta_1,\beta_3$ independent constant odd multiples of $t_1$, set $\gamma^2=t_3$, $\gamma^1=t_2$ and let $\gamma^3$ have only its second directional jet $t_1$ at the point. With ghosts zero the classical variation of $2D_2^2\beta_2$ contributes the nonzero color trace $2\operatorname{tr}(t_1[t_1,t_2]t_1t_3)=-i/4$ before the classical coupling and fermion-order factors. [Two-loop sign and independent Pauli contraction](code/fan_jet_handles.py)

**Four derivatives that isolate the maximum genus through six loops.** In the table a pair $ij$ means to apply $(\partial_{i,+\dot a}\partial_{j,+}{}^{\dot a})^2$ to the indicated slots. The square uses two independent dotted contractions. Let $[pq]=p_{+\dot a}q_+{}^{\dot a}$. Multiply the last column by $N^{n-2h}[pq]^2$, with $h=\lfloor n/2\rfloor$, to replace the entire undifferentiated color polynomial in the amplitude at the top of this file. The two cyclic traces retain their relative sign $(-1)^{n+1}$, and the remaining phase, heat factors, trace normalization and original $n$ dotted contractions are unchanged.

| $n$ | Input derivative combination | Surviving $h$ | Coefficient |
|---|---|---|---|
| 2 | $-3(12)+(13)$ | 1 | $-1/2$ |
| 3 | $-2(13)+(14)$ | 1 | $-9/50$ |
| 4 | $-12(12)-13(14)+9(15)$ | 2 | $-64/5$ |
| 5 | $17(14)-12(15)+3(25)$ | 2 | $480/49$ |
| 6 | $62(16)-81(17)-9(24)+133(25)$ | 3 | $2025/2$ |

For fixed assignment, the dotted jet on two path positions $i<j$ gives $[pq](a_i b_j-a_j b_i)$. Write $d=j-i$. The moment of its $r$th power after stripping $[pq]^r$ is
$$
\frac{(2)_r^2\,(d)_r}{(m)_r^2(m+1)_r},
\qquad (v)_r=v(v+1)\cdots(v+r-1).
$$
Indeed take two sets of $m$ independent unit-rate exponentials, one per spoke per endpoint. Normalizing each by its sum gives independent uniform simplex points. Both sums are independent of the fractions and their ratio order. For each spoke the sum of its two exponentials is gamma of shape 2, independent of the uniform fraction belonging to the first endpoint. The ratio order is the order of these uniform fractions. The original minor is the product of the two spoke sums times the difference of the two ordered uniform fractions, divided by both total endpoint sums. The uniform difference of rank $d$ has beta parameters $d,m+1-d$ by integrating the uniform spacings. Taking moments gives the displayed expression. At $r=2$ the only angular weight needed is $36d(d+1)/[m^2(m+1)^3(m+2)]$.

Multiply this weight by each original nested-commutator trace with the slot pair in its actual path positions. Fundamental completeness splits or joins the trace at each repeated generator. Sorting the empty index loops by their number produces the separate powers of $N$. Both cyclic traces in all lower powers cancel exactly for every row in the table. An independent calculation for the even loop orders evaluates the original Pauli products and their adjoint Levi-Civita chains at rank two, with the angular factors from iterated integrals of ordered uniform points. At odd $n$ a nonzero symmetric cubic trace requires, for example, $SU(N)$ with $N\ge3$. [Exact ordered moments, all color columns and finite-rank checks](code/fan_jet_handles.py)

### Existence of a finite local cochain selecting any ribbon genus

At every $n\ge2$ and $0\le h\le\lfloor n/2\rfloor$, integer combinations of input jets on this single trace give a nonzero two-ghost coefficient supported on just $N^{n-2h}$. The jet order in the general construction may exceed four. The proof also works when the external color test has only the antisymmetric trace.

First treat the weights of spoke permutations as independent. With one active slot, the two color traces contract to $\tau\operatorname{tr}([D,E]F)$. Replace that slot in the input matrix word by two adjacent slots. Take weights $+1,-1$ for their forward and reverse adjoint path orders. Their commutator and completeness replace this block by $N\tau$ times the original slot. For example $\sum_{A,B}t_A t_B\otimes[t_A,t_B]$, paired with matrices $X,Y$ in the two factors, gives $\tau^2(\operatorname{tr}X\operatorname{tr}Y-N\operatorname{tr}(XY))$. In the nested adjoint contraction the central term vanishes and the path-order sign gives $+N\tau$. For three adjacent input slots $A,B,C$, keep the two adjoint path orders $ACB$ and $BCA$ with unit weights. Stripping $\tau^3$, the first contraction gives
$$
\operatorname{tr}([D,E]F)+
\operatorname{tr}D\operatorname{tr}E\operatorname{tr}F
-N\operatorname{tr}(DE)\operatorname{tr}F .
$$
The reverse gives the same commutator trace and the opposite last two terms. Their sum is exactly $2\tau^2$ times the single-slot tensor, including its contraction inside any larger word. These identities follow by the same trace splitting on arbitrary $D,E,F$. The center on internal lines has zero adjoint action, so $U(N)$ completeness applies also to the $SU(N)$ graph. Inserting $h$ triples and $n-2h$ pairs multiplies the single-slot color by $2^h\tau^n N^{n-2h}$. The weights have path-reversal parity $(-1)^n$. At six loops the explicit nested contractions give $N^6,2N^4,4N^2,8$ in increasing handle number after stripping common trace factors. [Block expansions and six-loop checks](code/fan_jet_handles.py)

It remains to realize these weights by actual local derivatives, without an external-momentum-dependent choice of tangent direction. For each pair of letter labels form the signed minor $a_i b_j-a_j b_i$. An input dotted contraction realizes it times the common factor $[pq]$. For each permutation $\sigma$ let $L_\sigma$ be the sum of these minors over the pairs in the order $\sigma$. Let $V$ be the product of the squares of all distinct pair minors. On the chamber of $\sigma$ the magnitude of $L_\sigma$ is the sum of all absolute minors. The reverse chamber has the same magnitude and opposite sign. On any other chamber with $V>0$ both signs occur in that sum and the triangle inequality is strict. The common permutation-invariant majorant $V(\sum_{i<j}|a_i b_j-a_j b_i|)^2$ has a positive maximum on the compact product of simplices. The supremum on the target chamber reaches that maximum by relabeling, and can be approached on sets of positive volume. On the closure of each chamber outside the reversal pair, $V L_\sigma^2$ has strictly smaller maximum; on a shared ordering wall $V=0$.

For even $n$ use $(V L_\sigma^2)^k$; for odd $n$ multiply by one more $L_\sigma$. All trial polynomials at fixed $k$ have the same homogeneous degree and the required reversal parity. Their integrals on the target reversal pair dominate the sum of absolute integrals on all other chambers for sufficiently large finite $k$, by the strict gap in maxima. For the odd choice $L_\sigma$ has a constant sign on the target and is nonzero at every positive maximizing point. On a set of one representative per reversal pair, the resulting rational moment matrix is strictly diagonally dominant after row normalization and hence invertible. Substituting the dotted input derivatives in the rational linear combination for the desired weights gives the same color selection for arbitrary $p,q$ times a common power of $[pq]$. Clearing denominators gives integer coefficients. The common degree has parity $n$, so the base spinor bracket together with these extra brackets has even exchange parity, as required for the antisymmetric trace on the two odd ghost modes. A nonzero $SU(2)$ commutator tests the result at every $n$. This constructs matrix elements of the gauge-theory supercharge at prescribed ribbon topology. Their comparison with quantum IIB diagrams requires the interacting D3 map on these composite operators.

### Which two-ghost measurements survive the boundary differential?

Take the finite jets of two ghosts of colors $D,E$ and one scalar $\gamma^I$ of color $F$. Write $p,q,s$ for their tangent momenta and $T=p+q+s$ for the total. Relative cochains contain a derivative on each ghost. The two cyclic traces split into $\operatorname{tr}(F[D,E])$ and $\operatorname{tr}(F\{D,E\})$, where braces denote the matrix anticommutator. Their polarized momentum polynomials are respectively symmetric and antisymmetric in $p,q$. The Lie algebra automorphism $X\mapsto-X^{\mathsf T}$ separates these two sectors. Every positive $Q_n$ vanishes on both and on a bilinear ghost–scalar primitive, since their odd-coordinate degree is one and $Q_n$ lowers it by $3n$. [Classical fields](../../sources/LLC_main.tex#L174) · [Loop grading](../../sources/LLC_main.tex#L4827)

**Commutator trace.** Let $f_T(p,q)$ denote the symmetric polynomial at fixed total momentum, so the third input has momentum $T-p-q$. Applying $Q_0$ and using invariance of the trace gives the coefficient $f_T(p,q)+f_T(p+q,r)$ of $[[X,Y],Z]$ in the cyclic Jacobi expression on three ghosts of momenta $p,q,r$. Choosing $X=Y=t_1$, $Z=t_2$ in the Pauli basis and permuting the slots shows that the only relation among the three color tensors is Jacobi. Closure is therefore equivalent to
$$
 f_T(p,q)+f_T(p+q,r)=f_T(q,r)+f_T(p,q+r).
$$
Differentiate in $r$ at zero. If $D_2$ denotes the directional derivative in the second argument, this gives $D_2 f_T(p,q)=D_2 f_T(p+q,0)-D_2 f_T(q,0)$. Symmetry in $p,q$ and a mixed derivative show that $D_2 f_T(v,0)$ is a closed one-form in $v$. Its polynomial primitive is $H_T(v)=\int_0^1 D_2 f_T(tv,0)[v]\,dt$. Since $f_T(p,0)=0$,
$$
 f_T(p,q)=H_T(p+q)-H_T(p)-H_T(q).
$$
This is exact for the full boundary differential. In component normalization, a primitive obtained by applying the symbol $H(p,s)$ to $c^A(\gamma^I)_A$ has polarized coefficient
$$
 g\,C_{DEF}\big[H(p,q+s)+H(q,p+s)-H(p+q,s)\big].
$$
Here $C_{DEF}$ is evaluated on the two ordered ghost colors and the scalar. The first two terms come from varying the scalar and moving the odd differential past the first ghost. The third comes from $Q_0c^A=-g C^A{}_{BC}c^Bc^C/2$. Thus $H(p,s)=-H_{p+s}(p)/g$ gives $C_{DEF} f_T(p,q)$. All formulas retain any overall graph coefficient. For the loop measurements the primitive is of nonnegative order in $g$.

For the actual directional two-loop handle write $K(p,q)=[pq]^2((p(u)^2+q(u)^2)/2-4p(u)q(u)/3)$. The unpolarized output has symbol $C_{DEF}K(p,q)$, multiplied by $ig^2\tau^3/[48(32\pi^2)^2]$. In $SU(2)$ put $u=e_1$, use odd ghost waves of colors $(1,1,2)$ with momenta $(e^1,2e^1,e^2)$, and a constant scalar of color $2$. Subtract one from each exponential so the ghost vanishes at the insertion. The two product-derivative terms in the variation are $-26/3,-2/3$; the ordered variation with both ghost assignments is $2g(-26/3+2/3)=-16g$. With $\tau=1/2$ the result is $-ig^3/[24(32\pi^2)^2]$. Pure exponentials and scalar momentum $-3e^1-e^2$ give the same variation with zero total momentum. Taking $H(p,s)=K(p,s)/g$ instead gives an exact three-letter cochain whose polarized coefficient at $s=0$ is $2C_{DEF}K(p,q)$, precisely the measured two-ghost coefficient. All changes differentiate the scalar. This holds for any symmetric $K$ vanishing on the axes, including the all-genus commutator selectors. The raw diagram is nonclosed. Its complete set of constant-scalar measurements can be reproduced by a quantum boundary.

**Symmetric trace.** For $SU(N)$, $N\ge3$, denote the antisymmetric momentum coefficient by $h_T(p,q)$. Pair with a traceless scalar and use the traceless parts of $\{[X,Y],Z\}$ and its cyclic permutations. These three tensors are independent. Taking $X=Y=\mathrm{diag}(1,-1,0)$ and $Z=E_{13}$, where $E_{ij}$ is the matrix unit, forces the coefficients of the last two to agree. Cyclic permutations make all three equal. Taking $X=\mathrm{diag}(1,-1,0)$, $Y=E_{12}$, $Z=E_{21}$ makes all three equal to the traceless part of $2\,\mathrm{diag}(1,1,0)$, which has trace pairing $4$ with $\mathrm{diag}(1,1,-2)$. The same traceless tests embed in all higher ranks. Expanding $[\{X,Y\},Z]=\{[Y,Z],X\}-\{[Z,X],Y\}$ gives
$$
 h_T(p+q,r)=h_T(p,r)+h_T(q,r).
$$
Antisymmetry gives additivity in the other slot. Polynomial additivity in two tangent dimensions leaves exactly $h_T(p,q)=B(T)[pq]$, with $B$ an arbitrary polynomial representing total tangent derivatives. In the paper's conventions the $dz_1 dz_2\,d\theta_I$ component of $\operatorname{Tr}(dC)^3$ gives $3g^2[pq]\operatorname{tr}(F\{D,E\})$. Every total derivative is quantum closed. A nonzero $B$ has a nonvanishing Cartan evaluation at $\mathrm{diag}(1,1,-2)$; alternatively the only bilinear primitive has commutator color. Thus these classes are nonzero at leading formal order. The symmetric trace vanishes identically for $SU(2)$. [Graviton representatives](../../sources/LLC_main.tex#L4101)

This classification also gives a direct test of an incoming diagram. In a three-letter output of any single $Q_n$ fan at least one of the two differentiated endpoints is a ghost. Its two endpoint momenta enter with their area bracket to power $n$. At fixed $T$, take the part linear in each of $p,q$ by Taylor expansion at $p=q=0$. For $n\ge2$ it vanishes, including when the scalar occupies one endpoint. Thus any individual higher fan whose three-letter output is closed has zero coefficient of the cubic graviton and all its tangent translations in this sector. Lower-order terms of a transferred differential require a separate check. [Original endpoint tensor](../../sources/LLC_main.tex#L3993)

At three loops the four-derivative handle selector above has symmetric-trace kernel proportional to $[pq]^5$. For $s=0$ choose $(p,q)=((1,0),(0,1))$ or $((2,0),(-1,1))$. Both totals are $(1,1)$. The required values of $B(1,1)$ after removing the common graph coefficient would be $1$ and $16$. Hence even adding arbitrary derivatives of the scalar cannot complete these measured data to a closed three-letter cochain. The same fixed-total test excludes every odd power greater than one. [Exact component variation, conserved test and both color sectors](code/fan_jet_cohomology.py)

The regular classical D3 disk gives a control in the surviving channel. The normal Koszul dictionary sends the regular closed polyvector $\partial_{w_2}\wedge\partial_{w_3}\wedge\partial_{z_1}\wedge\partial_{z_2}$ to $\theta_2\theta_3\partial_{z_1}\wedge\partial_{z_2}$ on the brane. For the undifferentiated scalar of flavor 1 and two ghost waves vanishing at the point, the two tangent derivatives give $[pq]$ and the two odd ghost orderings give $\operatorname{tr}(F\{D,E\})$. The two ordered disk positions have configuration weight $\int_{0<t_1<t_2<1}dt_1dt_2=1/2$. This identifies the regular disk tensor. The normal Laurent flux is singular on the brane and its quantum transfer requires the punctured and contact diagrams. [Costello–Li, §7.2](https://arxiv.org/pdf/1606.00365#page=36) · [Cyclic disk configuration weights, Definition 3.1.3](https://arxiv.org/pdf/1505.06703#page=34)

#### Regular D3 product and integrated disk coefficient

For the regular source above, measure its strength by its actual tangent-coordinate commutator $[z_1,z_2]_\star=\epsilon\theta_2\theta_3$. In the integer grading $\epsilon$ has degree $-2$. Costello–Li's tangent cup representative is $\partial_1 a\,\partial_2 b$; its skew part differs by a Hochschild coboundary, the second-derivative product rule. In the cyclic representative the deformed product on $\Omega^{0,*}(\mathbb C^2)\otimes\Lambda(\theta_1,\theta_2,\theta_3)$ is
$$
 a\star b=ab+\frac{\epsilon\theta_2\theta_3}{2}
 (\partial_1a\,\partial_2b-\partial_2a\,\partial_1b).
$$
It is associative to all powers of this source: the order-one associator cancels by Leibniz and every order-two term contains $(\theta_2\theta_3)^2=0$. Tangent integration by parts with compact smearing gives the usual closed trace for the product at zero total momentum. Accordingly the linear change of the cubic term $\int\operatorname{tr}(C\star C\star C)/3$ in trace units has coefficient
$$
 \frac16\int\Omega\,d^3\theta\,\theta_2\theta_3
 \operatorname{tr}\big(C(\partial_1C\,\partial_2C-\partial_2C\,\partial_1C)\big)
$$
per unit of the measured commutator. The factor $1/3$ belongs to the cyclic action and $1/2$ to the skew binary product. The ordered disk simplex area alone does not specify a diagonal action coefficient. To use the paper's component action with $\operatorname{tr}(t_A t_B)=\tau\delta_{AB}$, divide this trace-unit action by $\tau$. The physical normalization of an external singular bulk wave relative to the coordinate commutator remains open. [D3 product and tangential cup representative](https://arxiv.org/pdf/1606.00365#page=38) · [Cyclic trace normalization](https://arxiv.org/pdf/1505.06703#page=25) · [Component action](../../sources/LLC_main.tex#L4784)

Evaluate the cubic density on two ordered odd ghost waves $ig\,\eta_1D,ig\,\eta_2 E$ and the scalar $\theta_1 F$, with tangent momenta $p,q,s$ summing to zero. Use $dz_1,dz_2$ on the left and $\int d^3\theta\,\theta_1\theta_2\theta_3=1$. The sum of the three scalar placements gives $[pq]+[qs]+[sp]=3[pq]$. The disk and the $dz_1dz_2d\theta_1$ component of $\operatorname{Tr}(dC)^3$ have signed symbols, in trace units,
$$
 -\frac{g^2}{2}[pq]\operatorname{tr}(F\{D,E\}),\qquad
 3g^2[pq]\operatorname{tr}(F\{D,E\}).
$$
For the disk the two ghosts supply $i^2=-1$. In the form coefficient moving a ghost coefficient through a tangent differential gives a second minus sign. At unit area and $D=E=F=H=\mathrm{diag}(1,1,-2)$ the symbols are $6g^2,-36g^2$. This fixes the relative sign and the ratio $-1/6$ per measured coordinate commutator in the displayed matrix action. A constant scalar in a conserved three-wave measurement would force $[pq]=0$; the nonzero measurement includes its momentum and all its differentiated placements. [Exact exterior and matrix replay](code/regular_cubic_disk.py)

The cubic regular disk has no higher-surface correction on this projection among vertices obeying the fivefold open–closed degree axiom. The regular rank-four closed input has degree two; the two ghosts and one scalar have degrees $-1,-1,0$ and the tangent Dolbeault smearing supplies two. With surface genus $g_{\rm s}$, $h$ boundaries carrying marks, and $f$ empty boundaries, the functional degree is $2(2g_{\rm s}-2+h+f)$. Thus $2g_{\rm s}+h+f=1$. The marked open fields require $h\ge1$, leaving only the disk. On this disk a separated tree with a sphere vertex would give at most two closed flags per sphere on average, counting its single external closed flag and its required connection to a disk. Sphere vertices have at least three. With no sphere, the sourced disk has three open flags; adding any ordinary cubic open vertex in a tree increases the number of external open flags. Thus the single sourced disk gives the full three-field coefficient within these vertices. A descendant of positive level adds twice its level to the input degree and is excluded on the same test. These degree statements use the integer-graded flat model; the fixed-flux punctured transfer has additional source and boundary data. [Fivefold degree axiom, Definition 4.3.1(4)](https://arxiv.org/pdf/1505.06703#page=42) · [Grading on the flux background](https://arxiv.org/pdf/1606.00365#page=49)

At rank two, for $t_i=\sigma_i/2$, $\{t_i,t_j\}=\delta_{ij}\mathbf1/2$. The deformation of a bracket of traceless currents therefore takes values in the central **matrix** direction. For example $[z_1t_3,z_2t_3]_\star=\epsilon\theta_2\theta_3\mathbf1/4$. Pairing the disk with three traceless fields vanishes, and pairing it with $(t_1,t_1,\mathbf1)$ at unit area gives $-g^2/2$ in trace units. At rank three,
$$
 [z_1 H,z_2 H]_\star=\epsilon\theta_2\theta_3(2\mathbf1-H).
$$
Its traceless Cartan part is nonzero. A linear change of representatives on the commuting currents cannot remove that part: ordinary commutators with $H$ have zero diagonal. At rank two each ordinary cubic open vertex contains a commutator and annihilates a central matrix leg. Any connected open tree with one sourced vertex and traceless external fields therefore feeds only traceless Lie words into the three slots of the sourced vertex, giving zero. With two or more such vertices and componentwise open propagators the odd diagonal kernels identify all their $\theta$'s, so the square-zero coefficient kills the tree. This explains the decoupling of this regular source from pure traceless rank-two open trees, including higher external valences. The linear bulk Laurent mode and the full $U(2)$ source both still exist.

#### D3 states and amplitudes on the regular quadratic superpotential

The regular scalar background used in the bulk measurement has another exact test when treated as a B-model superpotential. Put $W=\rho(z_1w_1+z_2w_2+c w_3^2/2)$ with $\rho$ the amplitude of the regular background. The brane is $w_1=w_2=w_3=0$. Use the exterior algebra of its three normal directions with even vacuum. Let $E_i$ denote exterior multiplication by the $i$th generator and $H_i$ the dual contraction, so $H_iE_j+E_jH_i=\delta_{ij}$. The Koszul matrix factorization with the convention $D^2=W$ is
$$
 D=w_1 H_1+w_2 H_2+w_3 H_3+
 \rho(z_1E_1+z_2E_2+c w_3 E_3/2).
$$
The open differential is the graded commutator with $D$. It sends the three $E_i$ to $w_i$, and the $H_i$ to $\rho z_1,\rho z_2,\rho c w_3/2$, respectively. In particular the closed odd generator $\eta=H_3-\rho c E_3/2$ satisfies $\eta^2=-\rho c/2$. Ordering the $E_i$ on the left gives the ordinary Koszul resolution of $(w_1,w_2,w_3)$ in the underlying complex. After removing it the remaining differential sends $H_1,H_2$ to $\rho z_1,\rho z_2$ and kills $\eta$. This proves the following description without setting $\rho$ to a nonzero value prematurely.

Over the ring $S=\mathbb C[c,z_1,z_2][[\rho]]$, the even coefficient sector is $S/(\rho z_1,\rho z_2)$. Its additional odd generator $z_1H_2-z_2H_1$ has coefficient module $S/(\rho)$ and squares to zero. Indeed the syzygies of $(z_1,z_2)$ are generated by $(-z_2,z_1)$ and the boundary of $H_1H_2$ multiplies it by $\rho$. Tensor this algebra with the one-generator Clifford algebra of $\eta$ over the coefficient ring, with the graded sign. At fixed $\rho\ne0$ only the Clifford factor remains and the open states are supported at $z_1=z_2=0$. At $\rho=0$ in the original complex all three exterior $H_i$ and both tangent coordinates survive. This specialization is different from setting $\rho=0$ after computing the cohomology over $S$.

The primary cubic regular tensor has a first-order obstruction on this background. Contracting $dW$ with $\partial_{w_2}\wedge\partial_{w_3}\wedge\partial_{z_1}\wedge\partial_{z_2}$ gives four distinct trivectors with coefficients $\rho z_2,-\rho c w_3,\rho w_1,-\rho w_2$. The first coefficient is nonzero even for $c=0$. In the open algebra the skew disk symbol on $(z_1,z_2)$ is $H_2\eta$ at zeroth order. Its variation gives $\rho z_2\eta$. The antisymmetrized Hochschild coboundary of a unary cochain vanishes on these two central even coordinates at $\rho=0$, so the displayed obstruction cannot be canceled there. In the primary closed Landau–Ginzburg complex on the punctured fivefold, contraction by $dW$ sends $(L/w_1)\partial_{z_1}$ to $\rho L$ for any tangent-independent scalar Laurent mode $L$. This vector has zero holomorphic divergence. Thus the cubic Laurent mode is nonzero at $\rho=0$ and killed by multiplication by $\rho$ in this formal primary complex. For $c\rho\ne0$ the critical point is at the origin on the excised brane. On the patch $w_1\ne0$ contraction has homotopy given by wedging $(\rho w_1)^{-1}\partial_{z_1}$; the $w_2$ patch uses $(\rho w_2)^{-1}\partial_{z_2}$ and the $w_3$ patch uses $(\rho c w_3)^{-1}\partial_{w_3}$. Hence the primary Koszul complex is acyclic off the brane in this massive regular model.

For $c\rho\ne0$ the bulk Jacobi ring, the functions modulo the five gradients of $W$, is one-dimensional. The Kapustin–Li disk formula with its unit bulk insertion puts the boundary operator on the right of $(\partial D)^{\wedge5}$. In the coordinate order $(z_1,z_2,w_1,w_2,w_3)$ its normalized antisymmetrized supertrace on $\eta$ is $\rho^3 c$. To compute it, write $N_i=E_iH_i$. The ordered product of five derivatives followed by $\eta$ is $-\rho^2(\rho c/2)N_1N_2(2N_3-1)$, whose supertrace is $\rho^3c$. Swapping two derivatives contributes only anticommutators times lower Clifford degree with zero supertrace. The five antisymmetrized slots therefore supply the same value after division by $5!$. The Hessian determinant of $W$ is $\rho^5 c$. With $\operatorname{STr}$ even minus odd this gives
$$
 \langle1;\eta\rangle_{\rm disk}=\frac{1}{\rho^2},\qquad
 \langle1;\eta A,\eta B,\eta C\rangle_{\rm disk}
 =-\frac{c}{2\rho}\operatorname{tr}(ABC),
$$
where $A,B,C$ here are three ordered Chan–Paton matrices. The disk with the boundary identity vanishes by parity. Sewing $h$ handles onto the single boundary with insertion $\eta$ multiplies the amplitude by $(\rho^5 c)^h$. An additional empty boundary gives zero since its boundary state is zero in the Jacobi ring. These are primary matter correlators of the regular Landau–Ginzburg model at a nondegenerate critical point. The negative powers of $\rho$ prevent Taylor extraction at the flat vacuum. The Euler-source measurement and the punctured string diagrams use the separate formal expansion at that vacuum. [Kapustin–Li, Eq. (4) and §3.2](https://arxiv.org/pdf/hep-th/0305136#page=13) · [Exact matrix factorization and all 120 signed terms](code/quadratic_d3_factorization.py)

#### The cubic normal obstruction on a regular D3 brane

Let $w=(x,y,z)$ be the three D3 normal coordinates and $B=(x^3+y^3+z^3)/6+bxyz$. For a constant displacement $u$ of the brane, use the preceding exterior and contraction operators to set
$$
D_u=\sum_{i\in\mathrm{normal\ directions}}
\left((w_i-u_i)H_i+E_i\int_0^1(\partial_iB)(u+v(w-u))\,dv\right).
$$
The Clifford identities and the fundamental theorem of calculus give $D_u^2=(B(w)-B(u))\mathbf1$. Thus a point displaced to $u$ has curvature $-B(u)$ relative to the fixed superpotential. To use deformation generators with positive contraction $H_i$, expand $D_{-t u}=D_0+t\eta(u)+t^2 A_2(u)$. The components are
$$
\eta(u)=\sum_{i\in\mathrm{normal\ directions}}u_i\eta_i,\qquad
\eta_i=H_i-\frac16\sum_{j\in\mathrm{normal\ directions}}(\partial_i\partial_jB)(w)E_j,\qquad
A_2(u)=\frac13\sum_{i\in\mathrm{normal\ directions}}(\partial_i B)(u) E_i .
$$
For the graded commutator $d=[D_0,-]$ this gives $d\eta(u)=0$, $dA_2(u)=-\eta(u)^2$, $A_2(u)^2=0$ and the exact cubic obstruction
$$
\eta(u)A_2(u)+A_2(u)\eta(u)
=\left(\frac{u_1^3+u_2^3+u_3^3}{6}+b u_1u_2u_3\right)\mathbf1 .
$$
The identity is nonexact since $D_0|_{w=0}=0$. Also $dE_i=w_i\mathbf1$ and $\{\eta_i,\eta_j\}=-(\partial_i\partial_j B)(w)\mathbf1/3$. Ordering the $E_i$ before the distinct $\eta_j$ gives the Koszul complex of the normal coordinates for each ordered subset of the $\eta_j$. Thus its cohomology is the exterior algebra of the three $\eta_j$ over tangent functions. Changing $A_2$ by a closed odd endomorphism changes the cubic coefficient by an anticommutator vanishing in that algebra. The single-generator cubic coefficient is $1/6$ and the coefficient polarized once in each of the three directions is $b$.

Every direction with $B(u)=0$ extends to the exact family $D_{tu}$. Every direction with $B(u)\ne0$ is obstructed at order three. Their projectivized zero locus is the Hesse cubic $x^3+y^3+z^3-3\lambda xyz=0$ with $\lambda=-2b$. In coordinates $X=z$, $Y=x-y$, $Z=x+y+\lambda z$, the equation is $ZY^2=4(\lambda^3-1)X^3/3-3\lambda^2X^2Z+2\lambda XZ^2-Z^3/3$. On the smooth locus $1+8b^3\ne0$ it gives
$$
j=\frac{110592\,b^3(1-b^3)^3}{(1+8b^3)^3}.
$$
This recovers the elliptic modulus from integrable commuting brane displacements in the regular open matter theory. Comparing an ordered disk operation or a singular D3 quantum Ward identity requires the respective open–closed maps. [Eight-dimensional matrices and signed Weierstrass reduction](code/cubic_d3_obstruction.py)

On a direction with $B(u)=0$, the cohomology of the displaced brane over $\mathbb C[[t]]$ detects an additional obstruction to persistence. Write $A_i$ for the integral in $D_{tu}$ and put $L_{ij}=\int_0^1(1-v)(\partial_i\partial_jB)(tu+v(w-tu))\,dv$. Subtracting the value of $A_i$ at $tu$ gives
$$
\left\{D_{tu},H_i-\sum_{j\in\mathrm{normal\ directions}}L_{ij}E_j\right\}
=t^2(\partial_iB)(u)\,\mathbf1,\qquad
D_{tu}|_{w=tu}=t^2\sum_{j\in\mathrm{normal\ directions}}(\partial_j B)(u)E_j .
$$
If $dB(u)\ne0$, the identity class has annihilator exactly $(t^2)$: the first equation kills it by $t^2$, while evaluating any boundary at $w=tu$ makes every entry divisible by $t^2$. After inverting $t$ the brane is contractible. For $dB(u)=0$, evaluation at $w=tu$ kills every boundary even over $\mathbb C((t))$, and the identity survives. Thus on a smooth projective cubic every nonzero direction on the curve gives a curvature-free formal family with contractible generic fiber. For a singular Hesse member $\lambda^3=1$, the three projective directions with noncontractible generic fiber are $(\omega^k:\lambda^2\omega^{-k}:1)$ for $k=0,1,2$ and a primitive cube root $\omega$. These solve $x^2=\lambda yz$, $y^2=\lambda xz$, $z^2=\lambda xy$. [Displaced matrix homotopies](code/cubic_d3_obstruction.py)

**Endomorphisms along a ray and the nine flexes.** Work over $\mathbb C[[t]]$ on $B(u)=0$. For a constant normal vector $v$, let $\eta_t(v)$ be its contraction with the odd operators on the left of the preceding homotopy, so $d\eta_t(v)=t^2dB(u)(v)\mathbf1$. The ordered products of the three $E_i$ followed by distinct shifted contraction operators still form a basis of the endomorphism complex over normal functions. Indeed reordering the change of generators is triangular in the number of contraction operators. Its differential on that basis is the Koszul differential of $w-tu$ followed by that of $t^2dB(u)$. The three independent entries of $w-tu$ form a regular sequence. If $dB(u)\ne0$, choose one vector on which the gradient is one. Its residual complex is multiplication by $t^2$ on $\mathbb C[[t]]$. The remaining two directions span the tangent plane $\ker dB(u)$. Their ordered monomials therefore give a free basis of rank four over $\mathbb C[[t]]/(t^2)$. If $dB(u)=0$, all three odd generators survive and their eight ordered monomials give a free basis over $\mathbb C[[t]]$.

For closed generators associated to $v,v'$ the product in cohomology obeys
$$
 \big\{[\eta_t(v)],[\eta_t(v')]\big\}
 =-t\,B''(u)(v,v')\,[\mathbf1].
$$
Indeed $L_{ij}=\partial_i\partial_j B(w)/6+t\partial_i\partial_j B(u)/3$ and the anticommutator before taking cohomology is $-2L_{ij}$. The difference from its value at $w=tu$ is exact by $dE_i=w_i-tu_i$. These relations and the ordered basis determine the whole graded algebra. At the singular Hesse node $u=(1,1,1)$, $b=-1/2$, the Hessian has diagonal entries one and off-diagonal entries $-1/2$. Its kernel is the radial direction and its transverse rank is two.

On a smooth member take $v$ tangent and independent of $u$. Euler's identity gives $B''(u)(u,v)=2dB(u)(v)=0$. For a vector $l$ with $dB(u)(l)\ne0$, the Hessian determinant in the frame $(u,v,l)$ is $-4dB(u)(l)^2B''(u)(v,v)$. Thus all odd anticommutators vanish precisely where the projective tangent line has contact of order three, a flex. In Hesse coordinates $\det B''+6b^2 B=(1+8b^3)xyz$. On a smooth member this gives the nine distinct points $B=xyz=0$.

For the second displacement from $tu$ to $tu-sv$, the scalar curvature relative to $B(w)$ is exactly
$$
 -B(tu-sv)=-\tfrac12t s^2 B''(u)(v,v)+s^3 B(v).
$$
At an ordinary point the first obstruction is quadratic in $s$, with a nonzero multiple of $t[\mathbf1]$. At a flex the first obstruction is cubic, with a nonzero unit coefficient $B(v)$ because a smooth cubic contains no line. Changing the second-order primitive by a closed odd endomorphism adds an anticommutator with $[\eta_t(v)]$, zero at a flex. At $b=-5/6$, $u=(1,1,2)$, $v=(1,-1,0)$ the quadratic coefficient is $-8t/3$. At the marked flex $u=(1,-1,0)$ choose $v=(2b,0,1)$. The cubic coefficient is $(1+8b^3)/6$. Its cube times $216$ is the bulk discriminant $t^2-64s^3$ in the notation of the planar loop ratios; the letters $s,t$ in that bulk expression are the two measured ratios. This last normalization uses the displayed Hesse frame and displacement. [Exact matrix anticommutators and both tangent lines](code/cubic_d3_obstruction.py) · [Bulk ratios](../iib-bcov/README.md#planar-primary-brackets-and-the-normal-elliptic-curve)

**Normal disk projection of the scalar loop output.** Let $W=\rho B$ for any homogeneous ternary cubic with smooth projective zero locus and $\rho\ne0$. Use $W$ in $D_0,\eta_i$ and take the ordered boundary insertion $\eta_1\eta_2\eta_3$ on the right. In the normal coordinate order $(x,y,z)$ the Kapustin–Li numerator is
$$
\frac1{3!}\operatorname{STr}\big((\partial D_0)^{\wedge3}\eta_1\eta_2\eta_3\big)
=\frac{\det(\partial^2 W)}8\,dx\wedge dy\wedge dz .
$$
Every product of fewer than six linear Clifford factors has zero supertrace on the exterior algebra of three generators. In its top exterior symbol each $\partial_iD_0$ can therefore be replaced by $\partial_iD_0-\eta_i=\sum_{j\in\mathrm{normal\ directions}}(\partial_i\partial_j W)E_j/2$. The sign follows from $E_1E_2E_3H_1H_2H_3=-(E_1H_1)(E_2H_2)(E_3H_3)$, with supertrace $+1$. The three quadratic gradient equations have no common point at infinity, so a regular value of the gradient map has $2^3$ simple preimages by Bézout. The residue of its Jacobian is therefore eight. With a holomorphic scalar insertion $f$, all these preimages tend to the origin as the value tends to zero. This gives the tangent-function-valued normal disk
$$
\langle f;\eta_1\eta_2\eta_3\rangle_{\rm disk,normal}
=\operatorname{Res}_0
\frac{f\,\det(\partial^2 W)\,dx\wedge dy\wedge dz}
{8(\partial_x W)(\partial_y W)(\partial_z W)}
= f|_{w=0}.
$$
Both numerator and denominator carry $\rho^3$. At the Fermat cubic the integrand is exactly $f\,dx\wedge dy\wedge dz/(xyz)$. At $\rho=0$ and at a singular projective cubic the isolated normal residue formula does not apply. In the fivefold the two tangent derivatives of $W$ vanish. For the product theory, pairing the resulting tangent function with the compact radial top-form test $\tau_z$ of unit integral used below gives its value at the tangent origin by angular averaging. The planar bulk scalar on $e^{z_1}$ and $e^{z_2}\partial_{z_1}$ is its constant normal coefficient times $e^{z_1+z_2}$. Hence this disk test returns exactly the two- and three-loop bulk scalar coefficients, with the original factors $\rho^4$ and $\rho^6$ respectively when the background is scaled. It evaluates an inserted bulk scalar in the regular matter theory. An intrinsic singular-brane differential would also involve puncture transfer and defect terms. [Planar bulk coefficients](../iib-bcov/README.md#planar-primary-brackets-and-the-normal-elliptic-curve) · [Generic symmetric Hessian and Fermat check](code/cubic_d3_obstruction.py) · [Kapustin–Li, Eq. (4)](https://arxiv.org/pdf/hep-th/0305136#page=13)

#### The planar scalar on the unfluxed punctured cubic

Remove $w=0$, set the D3 flux to zero, and work over $\mathbb C[[\rho]]$ at $W=\rho B$. Let $f=e^{z_1+z_2}$, the output wave of the planar binary. For a smooth projective normal cubic,
$$
 \operatorname{Ann}_{\mathbb C[[\rho]]}[f]=(\rho^3)
$$
in the local linearized punctured descendant theory. Consequently its measured loop outputs $\rho^4 Q_2|_B f$ and $\rho^6 Q_3|_B f$ are both exact there, including at formal $\rho$. On a singular projective member the same wave has zero annihilator. In particular both outputs survive on every singular finite Hesse member.

To construct a primitive write $a_i=\partial_i B$ in the three normal coordinates. The sets $a_i\ne0$ cover the puncture in the smooth case. On these sets let $Y_i=(\partial/\partial w_i)/a_i$. For any normal multivector $Y$ write its divergence-free lift as
$$
 fY-\partial_{z_1}\wedge f\partial_\Omega Y .
$$
Here $\partial_\Omega$ is divergence in the coordinate volume and $\partial_{z_1}f=f$. The internal primary differential on this lift is contraction with $dW$, since the free descendant divergence vanishes and the tree variation is $\partial_\Omega(WY)$ on a general primary $Y$. Contraction with $dW$ commutes with the lift. Each positive-genus linearization has only the three normal momentum directions and a single incoming momentum in its five-row determinants, hence vanishes. In Čech degrees zero, one and two use the lifts, respectively, of
$$
 \rho^2Y_i,\qquad \rho Y_i\wedge Y_j,\qquad -Y_1\wedge Y_2\wedge Y_3.
$$
With alternating Čech differential plus $(-1)^p$ times the internal differential in degree $p$, their total differential is $\rho^3f$.

For sharpness, send a level-$r$ field $Y$ to the pure normal component of $W^rY/r!$, retaining its normal and tangent coefficient function. This is a chain map to the normal Koszul differential $\iota_{dW}$ and is the identity on primary scalars. Indeed the original free term at level $r+1$ is $-\partial_\Omega Y$ and the primary term is $\partial_\Omega(W^{r+1}Y/(r+1)!)$. Their images sum to the pure normal projection of $W^r\iota_{dW}Y/r!$. It therefore suffices to rule out a primitive of $\rho^2f$ for the normal Koszul complex. The gradient charts and their intersections are Stein, so holomorphic Čech cochains compute the complex. In a hypothetical primitive let $V_i$ be the coefficient of $\rho$ in its normal vector on chart $i$, and let $C_{ij}$ be the constant coefficient of its bivector on the overlap. The degree-zero, one and two equations imply
$$
 \iota_{dB} V_i=f,\qquad \delta V=\iota_{dB} C,\qquad \delta C=0.
$$
The punctured three-space has $H^1(\mathcal O)=0$. On the coordinate puncture cover this follows by splitting Laurent series according to the coordinates with negative powers; only the empty set and the full set can give cohomology, in degrees zero and two. The same calculation holds with holomorphic tangent parameters. Thus $C=\delta M$ for normal bivectors $M_i$. The vectors $V_i-\iota_{dB}M_i$ glue and extend across the missing plane by Hartogs. At $w=0$ their contraction with $dB$ vanishes while $f$ is nonzero. This excludes $\rho^2f$ and all lower powers. On a singular member, evaluating the chain map on any remaining nonzero critical ray instead kills every normal Koszul boundary for every power of $\rho$. At $b^3=-1/8$ the bulk ratios are $-9b/8,27/8$, both nonzero. This comparison uses the unfluxed linearized bulk complex; relative quantum graphs and open terms are needed for the intrinsic D3 supercharge. [Seven-chart replay](code/cubic_puncture.py) · [Local planar loop outputs](../iib-bcov/README.md#planar-primary-brackets-and-the-normal-elliptic-curve)

**Turning on D3 flux.** On every finite Hesse member, let $U$ denote the normalized normal angular two-form and $b_x,b_y$ its one-form primitives for $xU,yU$. Write $F=U\partial_{z_1}\wedge\partial_{z_2}$ and let $N$ be its strength. For the tangent wave $f$ set $X_f=\partial_\Omega(f\partial_{z_1}\wedge\partial_{z_2})$. One power of $\rho$ is necessary and sufficient for $\rho f$ to admit a closed lift from $N=0$ in the local linearized complex regular in $N,\rho$. Necessity follows at $\rho=0$ from the first-order flux variation $U X_f$. Its primary normal two-form has nonzero Laurent residue, whereas a free descendant divergence raises the level. For sufficiency put
$$
 T=\partial_x^2 U,\qquad
 k=\tfrac12\partial_x^2(x b_x)-2\partial_x b_x+b\partial_x^2(z b_y),
 \qquad \bar\partial k=(x^2/2+b yz)T-U .
$$
The last equality follows from the product rule and $\bar\partial b_i=w_iU$. Take the normal two-form valued bivector
$$
 V=T\partial_x\wedge X_f-(\partial_x T)f\partial_{z_1}\wedge\partial_{z_2}.
$$
Its divergence vanishes and its contraction with $dW$ gives $\rho(\partial_x B)T X_f$. The field $kX_f$ has zero divergence and annihilates $W$. The full flux bracket on $V$ vanishes in normal form degree four. Its bracket on $kX_f$ vanishes because $X_f$ preserves the tangent bivector and all angular coefficients are independent of tangent position. Consequently the exact closed lift is
$$
 \rho f-N V+N\rho k X_f .
$$
All local positive-genus linearizations still vanish by the five-direction rank test, since the background flux coefficients depend only on the normals. At every singular finite Hesse member, multiplying this lift by $\rho^{2g-1}Q_g|_B$ gives a nonexact closed extension of the measured scalar output for $g=2,3$. Evaluation of the form-degree-zero primary at a remaining critical point kills every boundary, including terms with an explicit flux. This constructs representatives in the linearized punctured complex. Identifying the intrinsic quantum operation on lifted inputs requires the relative Ward terms. [Angular component identities](code/cubic_puncture.py) · [Normalized angular representatives](#a-cyclic-homotopy-on-smooth-polynomial-angular-modes)

For smooth $B$ the displayed closed lift of $\rho f$ has annihilator exactly $(\rho^2)$ over $\mathbb C[[\rho]]$ with flux formal in $N$. In particular both measured scalar outputs have exact lifts at that flux. To see this, suppress the common tangent wave and project a divergence-free primary to its terms without $\partial_{z_1}$. This projection is invertible onto the exterior algebra on the three normal vectors and $\partial_{z_2}$. The missing $\partial_{z_1}$ term is uniquely fixed by divergence, since $\partial_{z_1} f=f$. The projected internal differential is $\rho\iota_{dB}$ together with exterior multiplication by $NU\partial_{z_2}$ (with the form grading). Indeed $\partial_\Omega(FY)$ can contribute without $\partial_{z_1}$ only when its divergence removes that tangent vector. The whole flux term is zero on the summand already containing $\partial_{z_2}$.

Choose a smooth normal vector $H$ on the puncture with $dB(H)=1$, for instance the complex conjugates of the three gradient components divided by their squared norm in the coordinate Hermitian metric. Exterior multiplication by $H$ contracts the normal Koszul differential, with its form-degree sign. Construct a primitive of $\rho^3$ first at $N=0$ by successively removing the lowest normal form degree with this contraction and dividing by $\rho$. Its components in degrees zero through two are divisible by $\rho^2,\rho,1$ respectively, and their normal factors are $H$, $H\wedge\bar\partial H$, $H\wedge(\bar\partial H)^2$. The final Dolbeault derivative is zero because $\bar\partial H$ takes values in the rank-two kernel of $dB$. Subtract the full differential of this primitive from $\rho^2$ times the closed lift above. The difference lies entirely in the $\partial_{z_2}$ summand and is closed. In normal form degrees one, two and three its coefficients are divisible by $\rho^3,\rho^2,\rho$, respectively. In degree one this follows from the coefficient of $kX_f$. Degrees two and three also receive the flux differential of the degree-zero and degree-one primitives. In succession remove each remaining lowest degree by wedging with $H$ and dividing by $\rho$. Its Dolbeault remainder in the next degree loses at most one power; degree four vanishes identically. All primitives have polyvector rank at most four. Lifting them by the inverse divergence-free projection proves exactness in the minimal descendant complex, to all powers of the background flux in this local linearization. One lower power would give a primitive of $\rho^2f$ after setting $N=0$, contradicting the sharp unfluxed test. On singular finite Hesse members the critical-ray evaluation above instead survives every power of $\rho$.

**Value of the local binary on closed inputs.** Write $\widehat{\rho f}$ for the explicit closed lift above. Use $\widehat{\rho e^{z_1}}$ and $e^{z_2}\partial_{z_1}$ as the two inputs; the vector is divergence free, preserves the flux and annihilates $W$. For the elementary local genus-$g$ binary in the heat convention, on every finite Hesse member and at every $g\ge2$ its complete value is
$$
 \rho^{2g} Q_g|_B\,\widehat{\rho e^{z_1+z_2}} .
$$
Here $Q_g|_B$ is the original flat scalar matrix element at unit background strength. To evaluate the scalar trace kernel, give each of the two tangent waves an arbitrary additional normal momentum. A determinant-one shear fixing every normal momentum removes both shifts, since their tangent parts are the two coordinate covectors. Every five-row determinant, including those with the balancing root, is unchanged. The kernel therefore differentiates neither hard normal coefficient. All $2g$ cubic slots are differentiated three times, giving the same constant $\rho^{2g} Q_g|_B$ as the flat scalar measurement. If a normal scalar background slot is replaced by $NF$, the flux tangent bivector wedges to zero with $\partial_{z_1}$ before the root divergence, including on every differentiated angular coefficient.

It remains to apply the root divergence to the product of the actual inputs. Write $f=e^{z_1+z_2}$ and $P=\partial_{z_1}\wedge\partial_{z_2}$. In form-first order their product is
$$
 \rho f\partial_{z_1}+NTf\partial_x\wedge P-N\rho k f P .
$$
Its divergence is $\rho f-NV+N\rho kX_f=\widehat{\rho f}$, including the sign from the one-form $k$. At $b=-1/2$ the actual two-loop binary on these inputs has leading scalar $\rho^5 f/[92160(2\pi i)^{10}]$. Its restriction at a nonzero critical point detects a nonexact output of the local linearized bulk theory with D3 flux. On smooth members the same operation on these inputs is exact by the preceding homotopy. Establishing the full boundary supercharge still requires relative Ward diagrams and open terms. [Exterior product, divergence and flux differential](code/cubic_puncture.py) · [Independent normal shears of the two- and three-loop symbols](../iib-bcov/planar_cubic_primary.py)

**A single global punctured background with nonexact local outputs at every genus.** Use the entire fivefold and take $B=z^2/2+e^x+e^y-2-x-y$ in the normal coordinates of the disk computation. Each point $(2\pi i m,-2\pi i m,0)$ for a nonzero integer $m$ lies off the brane and satisfies $B=dB=0$. Translation by this vector preserves the entire background exactly. Now $\partial_z B=z$. In the closed scalar lift replace $T,k,\partial_x$ by $-\partial_z U,-\partial_z b_z,\partial_z$, respectively. The defining identity is $\bar\partial(-\partial_z b_z)=z(-\partial_z U)-U$.

The same normal-momentum shear removes both hard coefficient derivatives from the scalar trace kernel. Each of its quadratic $z$ slots requires two derivatives by double softness and supplies at most two. Hence exactly $g$ such slots occur. The other $g$ slots give only $x,y$ exponentials. The resulting scalar coefficient depends on $x,y$ alone, so the remaining root divergence on the lifted product cannot differentiate it: its only normal vector is $\partial_z$. Thus the elementary local binary on the same lifted scalar $\widehat{\rho e^{z_1}}$ and tangent vector $e^{z_2}\partial_{z_1}$ is the scalar graph coefficient times the closed output lift. At any of the displayed punctured points the scalar coefficient is exactly the flat value at the normal origin. With $\rho=1$ and tangent coordinates zero, the primary form-degree-zero measurement therefore obeys
$$
 (-1)^g(2\pi i)^{5g}Q_g>0\qquad(g\ge2),\qquad
 (2\pi i)^{10}Q_2=\frac1{62208},\qquad
 -(2\pi i)^{15}Q_3=\frac{49188779}{16295634432000}-\frac{7783}{3143448000}\zeta(3).
$$
The strict inequality follows from the full-rank fourfold periods with at least one mark on each of the $x,y$ axes. The two fixed inputs and every measured output have nonzero primary evaluations at the same critical point. Evaluation kills all linearized boundaries there, including those from descendants and those with explicit flux. This is a result on the entire punctured background for the elementary local bulk Hamiltonians with flux. The proof uses critical points outside the germ at the excised origin. [All-genus positivity and lower-genus graph evaluations](../iib-bcov/README.md#planar-primary-brackets-and-the-normal-elliptic-curve) · [Flux lift with a nonconstant scalar coefficient](code/cubic_puncture.py)

#### Regular D3 displacements on the critical line

On $B=(z^2+y^2h(x))/2$ use the closed punctured scalar lift for $\partial_z B=z$ above. At a moving point $(u,0,0)$ of the critical line the regular normal D3 factorization uses the same exterior multiplications $E_i$ and contractions $H_i$ on the three-dimensional normal exterior algebra:
$$
 D_u=(x-u)H_x+yH_y+zH_z+\tfrac12yh(x)E_y+\tfrac12zE_z,
 \qquad D_u^2=B\mathbf1 .
$$
Its identity survives over $\mathbb C[[u]]$ and after inverting $u$, since $D_u$ vanishes at the moving point. The three closed generators $H_x,H_y-h(x)E_y/2,H_z-E_z/2$ anticommute in distinct directions. Their squares in cohomology are $0,-h(u)/2,-1/2$. The ordered monomials give a free basis of rank eight: order the $E_i$ first and resolve the regular sequence $x-u,y,z$ by their differential. A holomorphic scalar acts on this family by its value at $(u,0,0)$ times the identity, since its difference from that value belongs to $(x-u,y,z)$. In particular the measured planar scalar part of the local bulk output acts nontrivially for every genus on the same regular family at $h=1+e^x$ as a formal germ in $u$. At $h=1+a e^x+b e^{-x}$ its $ab$ coefficient is constant on the axis and equals $6C_g(e^1+e^2)/(g-2)!$ in reduced units at unit background strength and zero tangent position. This is the scalar multiplication operator in regular matter. The full open-string supercharge requires relative quantum transfer. [Matrix computation](code/critical_line_brane.py) · [Bulk calculation, coordinate dependence and theta Ward graphs](../iib-bcov/README.md#a-planar-binary-on-a-punctured-critical-line)

At fixed $u$ with $h(u)\ne0$, restrict to the massive $y,z$ factors and use the ordered boundary insertion $\eta_y\eta_z$, where $\eta_y=H_y-h(u)E_y/2$ and $\eta_z=H_z-E_z/2$. With exterior forms ordered $dy\wedge dz$, the antisymmetrized Kapustin–Li numerator is $-h(u)$. Dividing by $(h(u)y)z$ and taking the transverse residue gives $-f(u,0,0)$ on a scalar insertion $f$. This checks the sign of the fiberwise disk independently of the endomorphism-class argument, which also applies on degenerate fibers. On the quadratic background the two extra transverse-Hessian sources have exact flux lifts. Their form-degree-zero local four-input coefficient probes the full repeated-soft polynomial at every genus, with the required Ward cancellation measured by a five-site bare graph integral. [Fiberwise trace and source boundaries](code/critical_line_brane.py) · [Five-site bulk Ward calculation](../iib-bcov/README.md#a-bare-ward-integral-with-five-hard-sites-for-every-repeated-soft-direction)

#### The flux class on the mixed quadratic in the descendant theory

Work to first order in an external flux on $\rho B$ with $B=z_1w_1+z_2w_2+cw_3^2/2$. Use the state differential in the positive-$\zeta$ vertex convention, with free part $\bar\partial-\zeta\partial_\Omega$. For a primary variation $Y$, its nonlinear primary part at tree level is $\partial_\Omega(\rho B Y)$. On tangent-independent primary tangent polyvectors the free divergence is zero, and the primary part contracts $\partial_{z_i}$ to $\rho w_i$. This restricted differential follows from the actual BCOV vertex. At every positive genus the local two-primary Hamiltonian linearization on a quadratic vanishes: the level-zero formula for $E_0(B;Y)$ would require $1-2g-j\ge0$ with $j\ge0$. In particular $\rho L$ is a boundary for **every** tangent-independent Laurent scalar $L$ on the triple overlap, with primary primitive $(L/w_1)\partial_{z_1}$. This includes the three cubic modes $-\partial_{w_I}f$. These statements concern the local bulk Hamiltonians before punctured quantum transfer. [Quadratic Euler derivative](../iib-bcov/operator/README.md#all-quadratic-backgrounds-and-the-unweighted-d3-linking-residue)

On the open patches $U_i=\{w_i\ne0\}$ use the alternating Čech coboundary plus $(-1)^p$ times the internal state differential in Čech degree $p$. Put $u_i=\partial_{z_i}$. For each positive integer $n$ the tangent-bivector mode has a closed extension with the following nonzero entries:
$$
 U_{123}:\ \frac{u_1\wedge u_2}{w_1w_2w_3^n},\qquad
 U_{i3}:\ -\frac{\rho u_i}{w_iw_3^n}\quad (i=1,2),\qquad
 U_3:\ -\frac{\rho^2}{w_3^n}.
$$
The internal differential of the first entry is $\rho(u_2/w_2-u_1/w_1)/w_3^n$. The alternating coboundary of the pair of middle entries cancels it. On each double overlap their signed internal differential is $+\rho^2/w_3^n$, canceled by the last entry. At $n=1$ the triple entry has unit normal residue. Any primary Hamiltonian boundary has tangent-bivector coefficient given by normal divergences of a rank-three polyvector, each with zero Laurent residue. A Čech coboundary from the double overlaps also has zero residue. The completed unit-flux state therefore remains nonexact for this linearized classical descendant theory, even allowing all descendant inputs.

There is a second-order effect in varying its representative. Multiply the whole extension by a function $P(w_3)$ such that $q=P(w_3)/w_3^n$ is regular. Subtract the total differential of the cochain with entries $q u_1\wedge u_2/(w_1w_2)$ on $U_{12}$ and $\rho q u_i/w_i$ on $U_i$, $i=1,2$. The remaining cochain is precisely the **global primary scalar** $-\rho^2 q$, on every chart. For $n=1$, $P=w_3^2/2$, this gives $-\rho^2 w_3/2$. On the mixed-only quadratic ($c=0$) it is nonexact. At $z_1=z_2=w_1=w_2=0$, $w_3\ne0$, the background and its gradient vanish, so every primary output $\partial_\Omega((\rho B)^{k+1}Y/(k+1)!)$ from a regular level-$k$ primitive vanishes there. The scalar has nonzero value there. For invertible $c$ it instead has the divergence-free primary primitive $-\rho\partial_{w_3}/(2c)$. Thus multiplication of the chosen flux lift and transport of a free-exact variation differ by a specific regular scalar. An exact lift on the mixed background adds $+\rho^2w_3/2$. The separated puncture and BV heat contributions to a Ward identity on the compact sources are still independent questions. [Signed Čech differential on all seven overlaps and exact primary primitives](code/quadratic_d3_factorization.py) · [Two compact source measurements](#two-compact-linking-sources-with-the-flux-varied-on-its-own-leg)

#### A two-normal bulk interaction on three closed puncture states

Take the scalar background $B=(w_1^2+w_2^2)/2+w_1^4/24$ and the D3 flux $NF=NU\,\partial_{z_1}\wedge\partial_{z_2}$ on the punctured fivefold. At genus $g\ge2$ take the primary scalar parts and the tangent vector
$$
 w_3\frac{z_1^{2g}}{(2g)!},\qquad
 w_3\frac{z_2^{2g-1}}{(2g-1)!},\qquad
 \frac{w_3^{2g-2}(z_1+z_2)^2}{2(2g-2)!}
 (\partial_{z_1}-\partial_{z_2}).
$$
The scalar $w_3H$ lifts by $-N b_3 X_H$, with $X_H=\partial_\Omega(H\,\partial_{z_1}\wedge\partial_{z_2})$ and $\bar\partial b_3=w_3U$. The correcting vector and the displayed third input have zero tangent divergence, preserve $F$, and annihilate $B$. This gives three closed states of the classical linearization. They are nonexact already on evaluating their degree-zero primary components at points of $w_1=w_2=0,\ w_3\ne0$: the Dolbeault differential and terms containing $F$ raise form degree, and $\partial_\Omega(B^{r+1}Y/(r+1)!)$ from every level-$r$ primitive $Y$ vanishes there since $B=dB=0$. The free divergence raises descendant level.

Evaluate the primary scalar at $z_1=z_2=w_1=w_2=0,\ w_3=1$. This point observable kills every classical linearized boundary by the same argument. Only the three displayed components enter the degree-zero output of the local bulk vertex. For its two scalar waves and vector wave take momenta $(\alpha,0,s),(0,\beta,t),(\gamma,\gamma,u)$ in the hard directions $(z_1,z_2,w_3)$. The root divergence is $\alpha-\beta$ and their determinant is $\alpha\beta u-\beta\gamma s-\alpha\gamma t$. On the critical axis the quartic normal jet measures the five-mark period $I_5>0$. The required mixed Taylor coefficient of the stripped vertex is
$$
\begin{aligned}
[s t u^{2g-2}](\alpha-\beta)
(\alpha\beta u-\beta\gamma s-\alpha\gamma t)^{2g}
 &=2g(2g-1)(\alpha-\beta)(\alpha\beta)^{2g-1}\gamma^2,\\
(-1)^g(2\pi i)^{5g}Q_g\phi\big|_{w_1=w_2=0}
 &=2g(2g-1)I_5,\qquad
 (2\pi i)^{10}Q_2\phi=\frac1{5184}
 \quad\text{on this axis}.
\end{aligned}
$$
The second line takes $[\alpha^{2g}\beta^{2g-1}\gamma^2]$. The determinants and divergence use all $4g+1$ tangent derivatives and all $2g$ derivatives in $w_3$ on these polynomial inputs, so no undifferentiated factor survives in these three coordinates. This follows also from the determinant character in each coordinate: the background has no momenta there, while the last divergence differentiates tangentially. Hence the local scalar coefficient throughout the punctured space depends only on $w_1,w_2$. It commutes with the flux and defines a nonzero classical linearized output class. The direct seven-field genus-two trace, including the repeated quadratic divisor $2!$, gives $1/5184$ on extracting the jets. [Two-normal bulk formula](../iib-bcov/README.md#three-nonexact-primaries-on-a-two-normal-scalar-background) · [Direct trace and brane replay](code/two_normal_puncture.py)

For a separate regular D3 test, use the same $B$ as a B-model superpotential. With the exterior and contraction matrices $E_i,H_i$ used above, set $D=w_1H_1+w_2H_2+w_3H_3+(w_1/2+w_1^3/24)E_1+w_2E_2/2$. Then $D^2=B$ and $\{D,H(z)E_3\}=w_3H(z)\mathbf1$. Both incoming scalar multiplications are boundaries on this brane. Since $D|_{w=0}=0$, every open boundary evaluates to zero there, whereas the outgoing scalar multiplication evaluates to $2g(2g-1)I_5\mathbf1$ in reduced units. Thus the local contact cannot by itself define a trilinear operation on the images of these scalar states in regular brane cohomology. The transferred operation with puncture and defect graphs, and its comparison with the Yang–Mills supercharge, remain open.

### Protection and count of all symmetric-color cubic operators

Work in the paper's relative polynomial superfield jet complex for $SU(N)$, $N\ge3$. The symmetric-color cubic classical cohomology consists precisely of the nonzero independent tangent derivatives of the components of $\operatorname{Tr}(dC)^3$, including the components in the expansion in $\theta$. Every one of these classes survives after inverting $\hbar$.

To see completeness, evaluate a cubic cochain on three exponentials in the two even tangent coordinates and three odd $\theta_I$ coordinates. Their momenta have two even and three independent odd components; pairing each odd momentum with its odd coordinate makes each exponential even. Use independent odd amplitudes on the three waves. The coefficient $A(p,q,s)$ of the symmetric cubic trace is then alternating in its three full arguments. The components recover all original jets by taking coefficients in the odd momenta. Extend a relative cochain by zero if any input is a constant current. Invariance under constant color rotations ensures that its differential also vanishes on such an input, so one can use unsubtracted exponentials.

At fixed total supermomentum $T$, put $A_T(p,q)=A(p,q,T-p-q)$. Apply the classical current bracket to four waves. The three brackets involving the last wave act on the symmetric product by the adjoint action, since the color pairing is invariant. The coefficient of $\{[X,Y],Z\}$ after the same trace pairing as above is, up to the overall classical coupling and orientation,
$$
A(p+q,r,s)-A(p,r,q+s)-A(q,r,p+s),\qquad T=p+q+r+s .
$$
Its two cyclic permutations give the other coefficients. The independent $SU(3)$ matrices above still test these identities with arbitrary formal odd momenta. Thus $A_T$ is additive in both arguments. Additivity of a polynomial on the even and odd vector space implies linearity: the substitution $p\mapsto 2p$ removes every degree except one. Hence $A_T(p,q)=b(T)(p,q)$ for a polynomial alternating two-form $b$ on supermomentum space. Cyclicity now reads $b(T)(q,T-p-q)=b(T)(p,q)$, or $i_E b=0$ for the Euler vector field $E$ whose value at $T$ is $T$. If the coefficients of $b$ have homogeneous degree $k$, the de Rham differential $d$ in $T$ and contraction $i_E$ give
$$
 (i_Ed+di_E)b=(k+2)b,\qquad
 A(p,q,s)=\frac{(db)(T)(p,q,s)}{k+2},\qquad T=p+q+s .
$$
Indeed evaluating the contraction on $p,q$ gives $(db)(T)(p+q+s,p,q)=(db)(T)(p,q,s)$ on even test supervectors. Odd-coordinate differentials may repeat. This writes the cubic as one superfield derivative on each letter, antisymmetrized in the slots, with a polynomial in total momentum supplying the remaining derivatives. These are exactly the components and total derivatives of $\operatorname{Tr}(dC)^3$. Conversely each satisfies the identity. At length three every invariant trace on traceless matrices is a combination of the two cyclic traces: invariant fundamental indices contract by permutations and every factor of length one vanishes. The quadratic invariant giving a possible classical primitive has commutator color. There is therefore no classical boundary in the symmetric sector.

For the quantum statement, the paper's fan annihilates the full superfield expression $\operatorname{Tr}(dC)^3$ at all positive orders before component extraction. Translation invariance gives the same for total derivatives. An incoming first transferred differential $d_1$ preserves word length and color-automorphism parity. Its possible sources are exactly the classically closed symmetric cubics just classified, on which $Q_1$ vanishes. For an incoming $d_r$ with $r\ge2$, take any highest weight of the dotted $SL_2$ acting on the two tangent derivatives. With $m$ total tangent derivatives and $j=0,1,2$ tangent form indices in $dC$, the representation is contained in
$$
\operatorname{Sym}^m(\mathbb C^2)\otimes\bigwedge^j(\mathbb C^2).
$$
For $j=0$ the highest weight has zero minus derivatives. For $j=2$ it has one. For $j=1$ the symmetric component has zero, and contracting one pair gives the other component with one. The source of $d_r$ requires $r$ fewer minus derivatives. Thus there is no such source at a highest weight. Equivariance excludes its whole irreducible representation, including all multiplicity spaces. This establishes protection of every nonzero class in this sector at finite rank. [Grading of transferred differentials](../../sources/LLC_main.tex#L4886) · [All-order single-graviton identity](../../sources/LLC_main.tex#L4119)

There is a simple count by the total number of tangent and odd derivatives on the three superfield jets. Let $t$ count one such derivative. The coefficient polynomials in the two even and three odd components of $T$ have series $(1+t)^3/(1-t)^2$. The free one- and two-form generators contribute respectively $5t$ and $13t^2$; the latter counts $1$ even–even, $6$ mixed, and $6$ odd–odd pairs. Contraction with $E$ gives an exact complex at positive total degree by the displayed homotopy. Its image on one-forms is the positive-degree polynomials. Taking the alternating difference of the dimensions of two-forms, one-forms and polynomials therefore gives the series of horizontal two-forms
$$
 \frac{(1+t)^3}{(1-t)^2}(13t^2-5t+1)-1
 =\frac{t^3(25+34t+13t^2)}{(1-t)^2}.
$$
Direct signed four-current color equations on the cubic jets have nullities $25,84,156$ at degrees $3,4,5$. At fixed source degrees $(n_+,n_-;\theta_1,\theta_2,\theta_3)=(k,0;2,1,1)$ for $k=0,1,2,3$, they give respectively $2,3,3,3$ closed cochains among $3,10,22,38$. The exterior sign check includes the three repeated odd form directions of $\operatorname{tr}((\gamma^1)^3)$. [Independent relative-jet enumeration](code/cubic_graviton_jets.py)

**The remaining cubic color tensor.** In the commutator trace the three-superwave symbol $A(p,q,s)$ is symmetric. Put $f_T(p,q)=A(p,q,T-p-q)$. Pairing the fourth current with the other three by the invariant metric gives precisely the symmetric cocycle identity for $f_T$ above, now for full supermomenta. This follows by expanding the coefficient of each nested commutator and using Jacobi. The same Pauli test shows necessity already in $SU(2)$.

There is an explicit relative quadratic primitive. Write $f_{T,m}$ for the part of degree $m$ in $p,q$ at fixed $T$. Relativity and symmetry give $m\ge2$. Applying the cocycle identity to four inputs $(p,p,q,q)$ in the two orders gives
$$
 f_T(2p,2q)-2f_T(p,q)
 =f_T(p+q,p+q)-f_T(p,p)-f_T(q,q).
$$
Thus $H_T(v)=\sum_{m\ge2\,;\,\mathrm{homogeneous\ degrees}} f_{T,m}(v,v)/(2^m-2)$ satisfies $f_T(p,q)=H_T(p+q)-H_T(p)-H_T(q)$. The relative condition in the **third** slot now gives $H_T(p)+H_T(T-p)=H_T(T)$. To obtain a relative bilinear, subtract from $H_T(v)$ the linear function of $v$ obtained as follows: for each positive homogeneous component of $H_T(T)$ of degree $j$ in $T$, take its derivative along $v$ and divide by $j$. At $v=T$ this returns $H_T(T)$ by Euler's identity, including the odd variables. There is no constant term. This subtraction preserves $f_T$ and gives $H_T(0)=H_T(T)=0$ and $H_T(T-p)=-H_T(p)$ for the corrected polynomial. Pair its two slots with the invariant color metric. The cyclic sum of its three bracket coefficients is
$$
 H_T(p+q)+H_T(q+s)+H_T(s+p)
 =H_T(p+q)-H_T(p)-H_T(q)=A(p,q,s).
$$
Restoring the classical differential's sign and coupling gives a relative bilinear primitive. Its full differential gives a quantum-exact completion $Q_0$ of the primitive plus $\hbar Q_1$ of the same primitive; all higher fans require more than two input letters. On the four weights $(k,0;2,1,1)$ for $k=0,1,2,3$, the closed commutator dimensions $1,3,5,7$ are independently reproduced by the ranks of the relative bilinear differentials. Hence the preceding counting series is the **entire associated-graded cubic cohomology** for $SU(N)$, $N\ge3$. It vanishes for $SU(2)$, whose only cubic invariant is the commutator trace. The word length refers to the leading term of a perturbative state. [Direct calculation in both color sectors](code/cubic_graviton_jets.py)

The paper's component Konishi relation illustrates why the leading word length matters. On the primitive $(\gamma^1)_A(\beta_1)^A$, the classical gauge terms cancel and the flavor-one Euler derivative of the superpotential gives $-\bar W(\gamma,0)$. The one-loop color contraction has adjoint Casimir $2N\tau$ in the paper's trace convention. Thus
$$
 \mathbf Q\big((\gamma^1)_A(\beta_1)^A\big)
 =-\bar W(\gamma,0)+\hbar\,\frac{g^2 N\tau}{32\pi^2}\,
 \delta_{DE}(\partial_{+\dot a}c^D)(\partial_+{}^{\dot a}c^E).
$$
No higher fan acts on this bilinear. The ghost bilinear is nonzero on two first jets with the same Cartan color. At its derivative counts $(1,1)$ only $d_1$ can enter. Its source has two undifferentiated letters and one copy of each odd-coordinate weight, so its invariant cochains are the three flavor-diagonal scalar–beta bilinears. All have the same $Q_0$ image $-\bar W$ and the same $Q_1$ image. Hence $Q_1$ annihilates every closed source. The ghost bilinear survives after inverting $\hbar$. In particular the unchanged scalar cubic $\bar W$ represents $\hbar$ times that nonzero quadratic class, even though its classical cubic class is zero. [Paper's one-loop coefficient and physical identification](../../sources/LLC_main.tex#L4067) · [Superpotential](../../sources/LLC_main.tex#L4761)

This example generates the entire possible quantum class of a classically exact **unchanged** cubic in the commutator channel. For a relative quadratic primitive $B$ write its alternating metric symbol at fixed $T$ as $H_T(p)$. Its classical closure condition is exactly $H_T(p+q)=H_T(p)+H_T(q)$. Thus $H_T$ is linear in $p$, and relativity gives $H_T(T)=0$. The Euler contraction argument applied to a one-form writes every such quadratic as total superderivatives of components of $\operatorname{Tr}(dC)^2$. Each is annihilated by all positive $Q_n$. There are no relative invariant linear primitives for $SU(N)$. The same first-page and dotted highest-weight argument as at length three excludes every incoming quantum differential. Consequently all nonzero closed quadratics survive over $\mathbb C((\hbar))$.

By nilpotence $Q_1 Q_0 B=-Q_0 Q_1B$. If it vanishes, $Q_1B$ is a closed quadratic and its higher differentials all vanish. Then $Q_0B$ is quantum closed and $[Q_0B]=-\hbar[Q_1B]$. This quadratic is independent of the choice of primitive. In the one-loop fan both endpoints have the same odd coordinate and carry one contracted pair of tangent derivatives. On two polarized superwaves the symbol of $Q_1B$ has the factor $[pq]$ in its tangent momenta, and its odd momenta enter only through their sum. At fixed total $T$, $[pq]=[p,T_{\rm tan}]$ is already linear in the split momentum $p$, where $T_{\rm tan}$ is the tangent part of $T$. It is a nonzero divisor even with exterior coefficients. For the full quadratic symbol to be linear in $p$, the remaining polynomial must depend only on $T$. Conversely every polynomial of total superderivatives times this area bracket is realized by differentiating the displayed Konishi relation in $z,\theta$; the kernel commutes with these supertranslations. Thus the possible quantum classes of unchanged commutator cubics are exactly $\hbar$ times the total superderivatives of the $dz_1 dz_2$ component of $\operatorname{Tr}(dC)^2$.

For a failure of closure take the same flavor-one primitive and put two tangent derivatives in direction $u$ on the scalar alone. Integrating the two spoke assignments of the paper's one-loop triangle multiplies the quadratic area bracket by $(p(u)^2+q(u)^2)/3+p(u)q(u)/2$. At fixed $p(u)+q(u)=1$ the values for $p(u)=0,2$ are $1/3,2/3$. Its unchanged classical cubic therefore fails the quantum closure test. [One-loop triangle](../../sources/LLC_main.tex#L2697) · [Exact split-derivative integral](code/fan_jet_cohomology.py)

**Linear bulk spectrum.** In the project's punctured-fivefold dictionary the Laurent normalizable modes map to polynomial forms on $\mathbb C^{2|3}$ and the divergence maps to $d$. A regular divergence-free source brackets with the flux by taking the divergence of its product with the flux. That product, if nonzero, maps to a zero-form: multiplication by the flux kills positive normal powers and already fills both tangent polyvector slots. Its divergence cannot reach form degree three. The surviving bulk Laurent modes in degree three are therefore closed three-forms, equivalently exact three-forms by Euler contraction at positive weight. Sending the horizontal $b$ above to $(db)/(k+2)$ identifies their weights and multiplicities with the protected cubic boundary states for every $SU(N)$, $N\ge3$. Direct independent matrices for the bulk Laurent divergence and the boundary four-current equation agree at all $231$ weights of total degrees $3,4,5$. The same linear bulk sector is present at nonzero flux $N=2$, whereas the cubic boundary sector vanishes there by the trace identity. This is a finite-rank requirement on a quantum holographic dictionary. [Signed linear state dictionary](NOTE.md#2-bulk-single-particle-states-drv-chk-at-t0) · [Two independent matrices at each weight](code/cubic_graviton_jets.py)

## Counting the faces

After summing the spoke permutations, introduce a matrix $X$ for the legs in the commutators. The terms of $\operatorname{tr}((\operatorname{ad}_X)^mD\,E)$ are $\operatorname{tr}(X^{m-k}DX^kE)$ with coefficient $(-1)^k\binom mk$. Every $X$ contracts with a different matrix at the input trace. Number the second-vertex gaps $0,\ldots,m-1$ in cyclic order $\alpha$. The first-vertex order transported across the strips is an $m$-cycle $\beta$. Its $m$ possible positions of $F$ recover all $m!$ assignments. Following both ends of a strip makes the face permutation $\beta\alpha$. With $c$ cycles there are $c-1$ empty color faces and $m-c=2h$ when $D,E,F$ occupy one face.

The gap of $E$ is $m-1$. At $k=0$ the number of placements giving the cyclic trace $DEF$ is its face length. At $k=m$ the order is $DFE$. At $0<k<m$ the gap of $D$ is $m-1-k$. If it belongs to the same face as $E$, the number of placements giving $DEF$ is its forward distance from $E$ in that face. Each step traverses precisely one first-vertex gap. The planar permutation is the identity and only $k=0$ contributes to $DEF$.

At genus one, $\beta\alpha$ is the descending cycle on three chosen gaps or two transpositions joining opposite pairs of four cyclically ordered gaps. These are the unique choices on their supports for which $\beta$ is an $m$-cycle. At $k=0$ their combined weight is $\binom{m+1}{4}+2\binom{m-1}{2}+\binom{m-1}{3}$. For $0<k<m$ the three-cycle contribution is $m+k-3$ and the two-transposition contribution is $(k-1)(m-1-k)$. Their sum is $k(m+1-k)-2$. Its $m$th finite difference vanishes for $m\ge3$. Subtracting the endpoint terms from that difference gives the displayed one-handle coefficient.

## Fixed-genus computation and the sum over all genera

Let $P_{m,h}$ denote the coefficient of $N^{m-1-2h}$ above. A face permutation of defect $2h$ moves at most $4h$ gaps. Deleting fixed gaps from it deletes the same gaps from the cyclic order and preserves the single-cycle requirement on the first vertex. Conversely choose the moved support of size $s$ and a permutation on this ordered support that moves every gap, has defect $2h$, and satisfies that requirement. Write $\ell$ for the face length at its last gap, and $d_j$ for the distance from the last gap to gap $s-2-j$, with zero for a different face. Its $k=0$ contribution summed over embeddings is $\binom ms+(\ell-1)\binom{m-1}{s-1}$. At $0<k<m$, prescribing $j$ moved gaps between $D$ and $E$ gives the multiplicity
$$
\binom{k-1}{j}\binom{m-1-k}{s-2-j}.
$$
As a polynomial in $k$ this has degree $s-2<m$. Its alternating binomial sum on internal $k$ is the negative of its values at $0,m$ with their endpoint signs. Thus the entire contribution of this support permutation to $P_{m,h}$ is
$$
\binom ms+(\ell-1)\binom{m-1}{s-1}
-\sum_{\substack{0\le j\le s-2\\\mathrm{gaps\ between}\ D,E}}
d_j\left[(-1)^j\binom{m-1}{s-2-j}
+(-1)^{m+s-2-j}\binom{m-1}{j}\right].
$$
At fixed $h$ this counts only permutations on $2h+1\le s\le\min(m,4h)$ for $h>0$. The direct face enumeration and the reduced algorithm agree through $n=9$.

For the limit, write $B_{m,h}$ for the unmarked number of transported cyclic orders at defect $2h$. The first binomial in the preceding formula sums to $B_{m,h}$. Every other term has degree at most $4h-1$ separately on even and odd $m$. The unmarked generating polynomial is
$$
\sum_{\substack{h\ge0\\\mathrm{ribbon\ genera}}} B_{m,h}N^{m-2h}
=\frac{\prod_{j=0}^{m}(N+j)-\prod_{j=0}^{m}(N-j)}{m(m+1)}.
$$
Indeed $N^{\#\mathrm{cycles}}$ is the permutation character on the $m$th tensor power of the fundamental representation. In its Schur–Weyl decomposition the $m$-cycle class sum acts by $(m-1)!$ times the character at a cycle divided by the representation dimension. Taking the trace against the second $m$-cycle squares that character. Only the hooks $(m-r,1^r)$ contribute, each with squared character one. The hook-content dimension formula gives the summand $\prod_{j=0}^{m-1}(N-r+j)/m$ for $0\le r<m$. The consecutive products telescope to the displayed difference.

Consequently $B_{m,h}$ is the elementary symmetric polynomial of degree $2h+1$ in $1,\ldots,m$ divided by $\binom{m+1}{2}$. Its leading coefficient as a polynomial of degree $4h$ in $m$ is $1/[2^{2h}(2h+1)!]$. Indeed in the $(2h+1)$st power of $\sum_{j=1}^m j$ each repeated index loses at least one power of $m$. This gives the same leading coefficient for $P_{m,h}$.

The limit can be taken in the full finite-$N$ color polynomial. Since $d_j,\ell\le s$ and $\binom{m-1}{j}\le2^s\binom ms$ for $0\le j\le s$, the absolute value of each support contribution is bounded by $[1+4h+2(4h)^2 2^{4h}]\binom ms$. Positivity of the unmarked counts gives
$$
\frac{|P_{m,h}|}{N^{2h}}
\le [1+4h+2(4h)^2 2^{4h}]
\frac{1}{(2h+1)!}\left(\frac{m(m+1)}{2N}\right)^{2h}.
$$
This bound is summable uniformly for $m^2/N$ in any bounded interval. If $n,N$ tend to infinity with $n^2/N\to x$, the exact limit of the full color coefficient divided by its planar value is
$$
\lim\sum_{\substack{0\le h\le\lfloor n/2\rfloor\\\mathrm{ribbon\ genera}}}
\frac{P_{n+1,h}}{N^{2h}}
=\sum_{\substack{h\ge0\\\mathrm{ribbon\ genera}}}
\frac{x^{2h}}{2^{2h}(2h+1)!}
=\frac{\sinh(x/2)}{x/2},
$$
with value one at $x=0$. The ratio compares contributions to the same boundary $Q_n$ with all Grassmann, spinor and propagator factors held fixed.

The unmarked long-trace two-point function has this same scaling function. The signed fan first differs at the next order. For fixed $h>0$ only permutations with $s=4h$ contribute degree $4h-1$ to $P_{m,h}-B_{m,h}$. Each such permutation consists of $2h$ disjoint transpositions, hence $\ell-1=1$. Its leading marked correction is $\binom{m-1}{4h-1}$, compared with the unmarked $\binom m{4h}$. The finite-difference terms have degree at most $4h-2$. Consequently the coefficient of $m^{4h-1}$ in the difference is $4h/[2^{2h}(2h+1)!]$. To exchange limit and sum after multiplying the difference by $m$, use $m\binom{m-1}{r}=(r+1)\binom m{r+1}\le s2^s\binom ms$ for $r<s$. The preceding summable bound then acquires at most an extra factor $4h$. For the same sequences as above,
$$
\lim m\sum_{\substack{0\le h\le\lfloor n/2\rfloor\\\mathrm{ribbon\ genera}}}
\frac{P_{m,h}-B_{m,h}}{N^{2h}}
=2\left(\cosh(x/2)-\frac{\sinh(x/2)}{x/2}\right).
$$
Thus the commutator insertions in this supercharge matrix element preserve the leading long-trace scaling and change the subleading coefficient. [Published unmarked scaling, equation (2.8)](https://arxiv.org/pdf/hep-th/0205033#page=7)

## One flux insertion in the local IIB hexagon

On $X=\mathbb C^2_z\times(\mathbb C^3_w\setminus0)$ use the project's unit normal Čech representative
$$
F=(w_1w_2w_3)^{-1}\partial_{z_1}\wedge\partial_{z_2}.
$$
Pull back the chosen flat, one-descendant genus-one local vertex, $\frac1{24}\int\operatorname{Tr}\alpha_1\det(\partial_i\partial_j\alpha_0)$ in reduced heat units, and substitute $\alpha_0\mapsto NF+\alpha_0$. Here $\alpha_0$ is primary and $\alpha_1$ is the level-one field. A holomorphic derivative of the smooth Dolbeault representative of $F$ still has normal antiholomorphic form degree two. Two such factors at one local vertex wedge to zero in three normal complex directions. Consequently the translated genus-one vertex starts with five fluctuating fields; the genus-two and genus-three vertices of arities at least seven and nine start with at least six and eight. The quadratic variation at $NF$ of these pulled-back positive-genus local interactions vanishes. [Flat local vertex and choice](../iib-bcov/operator/README.md#the-genus-one-bulk-term-on-the-one-descendant-sector) · [D3 background, §9](https://arxiv.org/pdf/1606.00365)

Polarize the four primary scalar waves by the quadratic jets $z_1^2/2,z_2^2/2,w_2^2/2,w_3^2/2$. In the Hessian determinant they fix four diagonal entries, leaving only $\partial_{w_1}^2$ on the flux. The coefficient pairing with $\alpha_1$ is
$$
\frac N{24}\partial_{w_1}^2F
=\frac{N}{12w_1^3w_2w_3}\,\partial_{z_1}\wedge\partial_{z_2}
\quad\longleftrightarrow\quad
\frac N{24}(d\theta_1)^2 .
$$
The arrow is the existing linear dictionary for the Laurent modes. Pair with the regular normal trivector $w_1^2\partial_{w_1}\wedge\partial_{w_2}\wedge\partial_{w_3}$ and take the normal residue (coefficient of $(w_1w_2w_3)^{-1}$). The result is $N/12$. The same test with primaries $z_1w_2,z_1w_3,z_2w_2,z_2w_3$ gives $-N/6$: for coefficients $s_1,\ldots,s_4$ the complementary four-row minor is $(s_1s_4-s_2s_3)^2$.

These four mixed inputs lift to closed linearized source modes at $NF$. Write the smooth flux as $F=f\,\partial_{z_1}\wedge\partial_{z_2}$ with $f$ a normal $(0,2)$-form. The class of $w_j f$ vanishes because its Laurent representative has exponent zero in direction $j$, so write $w_j f=\bar\partial b_j$ with $b_j$ a purely normal $(0,1)$-form. If $\{F,z_iw_j\}$ is the signed product of $w_j f$ with the other tangent vector, subtract the same signed product of $N b_j$ and that vector from the source. Its divergence and its bracket with $F$ vanish by tangential translation invariance. Its Dolbeault variation cancels the flux action. In the one-flux contact every tangential-vector addition to the four inputs wedges to zero with the tangent bivector already present. The value $-N/6$ persists for these four closed inputs.

Replacing the last two diagonal quadratics by $(p\cdot w)^2/2,(q\cdot w)^2/2$ selects a squared derivative along their complementary normal direction, using the chosen normal volume. The three coordinate directions and their three distinct pair sums give six linearly independent responses, with absolute coefficient determinant eight in the monomial basis $(d\theta_i)(d\theta_j)$, $i\le j$. These are the six directions of the scalar bilinear graviton in the linear spectrum. This calculation evaluates an off-shell local vertex on the punctured background. $\alpha_1$ is a descendant source, and the ordinary normal-trivector residue test has nonzero divergence. A quantum matrix element on physical cohomology would also need the transfer to external states and diagrams with separated vertices. [Laurent mode dictionary](NOTE.md#2-bulk-single-particle-states-drv-chk-at-t0) · [Polarization and residue](code/flux_hexagon.py)

## A primary vector from the genus-two and genus-three IIB vertices at D3 flux

At genus two take one flux $NF$ and five labeled regular scalar sources. Their generating jets are $e^{a w_1},e^{b w_1},e^{c w_2},e^{d z_2+s w_3},e^{e z_1+t w_3}$, with $a,b,c,d,e,s,t$ scalar parameters. For a normal momentum $q$ write $q_3$ for its $dw_3$ coefficient. The balancing momentum is the negative of the six inputs' sum. The signed five-row cofactors of the six inputs are, up to a common sign, $c d e q_3(b,-a,0,0,0,0)$; the seventh cofactor in the centering convention is zero. The seven-field action is minus the trace symbol $P$. Substitution in its two quartic invariants gives exactly
$$
P=-\frac{a^2b^2c^4d^4e^4}{31104}\,q_3^4 .
$$
The Euler derivative of the action with respect to the primary field is $-P$ times the remaining field product. The primary differential is minus its holomorphic divergence. The flux carries the ordered tangent bivector $\partial_{z_1}\wedge\partial_{z_2}$, so its divergence on this product supplies $e\partial_{z_2}-d\partial_{z_1}$. Select the coefficient of $a^2b^2c^4d^4e^5st$. The five inputs are therefore
$$
w_1^2/2!,\quad w_1^2/2!,\quad w_2^4/4!,\quad
w_3 z_2^4/4!,\quad w_3 z_1^5/5! .
$$
In reduced genus-two units the normal Laurent class in the elementary primary vector variation is
$$
-\frac{N\,\partial_{z_2}}{1296\,w_1w_2w_3^3}.
$$
Indeed $\partial_{w_3}^4(w_3^{-1})=24 w_3^{-5}$, and the two scalar factors $w_3$ leave a third-order pole. At the same parameter degree the term along $\partial_{z_1}$ has an extra power of $d$ and vanishes. The physical heat multiplier at genus two is $\lambda^4/(2\pi i)^{10}$. The rank-two Euler derivative here gives a rank-one elementary-field response; the hexagon in the preceding calculation differentiated the action with respect to a descendant. [Two-loop trace and its normalization](../iib-bcov/README.md#the-genus-two-bulk-action-and-integrated-cubic) · [Primary BV sign](../iib-bcov/operator/README.md#closed-point-classes-from-scalar-source-insertions)

The five scalar inputs lift to linearized closed modes. For a tangential monomial $w_3 H(z_i)$ use its derivative $H'(z_i)$ in the $b_3$ correction above. The additional tangent vector has zero divergence and brackets to zero with $F$: differentiating $H'$ with the flux produces a second copy of that same tangent vector. Its product with the explicit flux inside the trace vanishes. The primary output is closed for $\bar\partial+t\partial_\Omega+N\{F,-\}$. For exactness reduce modulo the descendant parameter $t$. A rank-one normal $(0,2)$ output that is a boundary is, modulo a Dolbeault boundary, the bracket of $NF$ with a globally holomorphic scalar. Indeed the form-degree-one part of a potential primitive's differential forces its scalar of form degree zero to be holomorphic. That scalar extends across the missing codimension-three plane. Bracketing it with $F$ differentiates only the scalar in tangent directions and gives at most a simple pole in each normal variable in the Laurent module. The displayed triple pole therefore represents a nonzero class of the **linearized** target. The linear dictionary takes $F$ to one and intertwines divergence with $d$. Since $\partial_\Omega(z_1 F)=f\partial_{z_2}$, this output maps to $-N\,dz_1(d\theta_3)^2/2592$. This tests the local quintic contribution to the variation. The complete operation transferred to physical cohomology also involves other vertices, propagators and choices of representatives.

The same local genus-two vertex can put **all four normal derivatives on the flux** with no loss of pole order. Use the five regular scalar waves
$$
e^{w_1},\quad e^{w_2},\quad w_1e^{z_1},\quad w_2e^{z_2},\quad
(w_1+w_2)e^{z_1+z_2}.
$$
In the ordered coordinates $(w_1,w_2,z_1,z_2)$ introduce their generating covectors $(a,0,0,0),(0,b,0,0),(s,0,d,0),(0,t,0,e),(u,u,h,h)$ and differentiate in $s,t,u$ at zero. Every nonzero five-row minor contains the normal flux momentum or its balancing recoil and is linear in that flux momentum. Since these five scalar covectors lie in one four-plane, the only surviving flux component is $q_3$. At $a=b=d=e=h=1$ the seven signed cofactors are $q_3(s-u,t-u,-1,-1,1,0,0)$ up to an overall sign. More generally the fully differentiated coefficient of $P$ is
$$
[stu]P=-\frac{q_3^4}{38880}\,a^2b^2d^2e^2h^2(a+b)(de-dh-eh).
$$
Any derivative falling on the scalar exponentials instead of $P$ leaves a factor $w_1$ or $w_2$ and hence vanishes in normal Laurent cohomology. For unit strengths, $[stu]P=q_3^4/19440$. Applying the positive divergence to $P$ in the primary differential and using $\partial_{w_3}^4 f=24/(w_1w_2w_3^5)$ gives
$$
\frac{N e^{2z_1+2z_2}(\partial_{z_2}-\partial_{z_1})}{405\,w_1w_2w_3^5}
\quad\longleftrightarrow\quad
\frac{N}{19440}d(e^{2z_1+2z_2})(d\theta_3)^4 .
$$
The four determinants in $P$ bound the number of derivatives on the flux by four. The polynomial input obtained by selecting $a^3 b^2 d^4 e^3 h^2$ consists of $w_1^3/3!,w_2^2/2!,w_1z_1^4/4!,w_2z_2^3/3!,(w_1+w_2)(z_1+z_2)^2/2!$. Its output is $-N\partial_{z_2}/(1620w_1w_2w_3^5)$, corresponding to $-N\,dz_1(d\theta_3)^4/38880$.

For a scalar input $w_j H(z)$ the flux bracket is $w_j f$ times the tangent Hamiltonian vector of $H$ with its bracket sign. Use the opposite signed $N b_j$ times that vector in the lift. A tangent Hamiltonian vector preserves the tangent bivector, has zero tangent divergence and does not differentiate the normal coefficient. The corrected inputs are closed, and all tangent-vector corrections wedge to zero against the explicit tangent bivector at this one-flux contact. The fifth-order pole survives the simple-pole primary boundary test above. Scaling $d=e=h$ by one generating parameter extracts gradients of every power of $z_1+z_2$. Tangent $SL_2$ polarization spans the gradients of all tangent polynomials. The normal $SL_3$ preserves the unit residue class, and its pure fourth powers span the fifteen symmetric quartic normal forms. Thus the contact image contains $d\mathbb C[z_1,z_2]\otimes\operatorname{Sym}^4\langle d\theta_1,d\theta_2,d\theta_3\rangle$ in the linearized target. This is the image of the specified local quintic Taylor coefficient, before transfer of the full quantum operation. [Independent five-dimensional replay and polynomial witness](code/flux_quartic_pole.py)

The same mechanism evaluates the saturated planar trace at arbitrary genus. Use the project's $2n-1$ zero-total planar covectors in the $dw_1,dw_2$ plane for which the reduced elementary scalar period is
$$
\frac{n!(n-1)!(n-1)^2(2n-1)}{3n-1}\,J_n .
$$
Here $J_n<0$ is the second planar resistance period and $n\ge2$. Give each of these covectors its scalar exponential jet, and take the last two regular scalar inputs to be $w_3 e^{z_2},w_3 e^{z_1}$. Exponentials stand for coefficientwise generating jets. Together with the flux and the balancing slot the four retained momenta are a frame complementary to the planar pair, with one relation among the four. More explicitly apply to the original five-dimensional coordinate covectors the linear map with ordered images $dw_1,dw_2,dz_2+s\,dw_3,dz_1+t\,dw_3,q$. Its determinant equals $q_3$ up to a fixed orientation sign, and all $2n$ determinant factors in the trace acquire that same multiplier. Permuting the four retained labels leaves the trace invariant. The primary symbol is therefore precisely the displayed period times $q_3^{2n}$. The derivatives in $s,t$ attach the two factors of $w_3$ to the input waves. Their product with the other normal-plane exponentials is $w_3^2$ because the planar momenta sum to zero. This gives the reduced flux-linear contact
$$
N(2n)!\,\frac{n!(n-1)!(n-1)^2(2n-1)}{3n-1}J_n\,
\frac{e^{z_1+z_2}(\partial_{z_2}-\partial_{z_1})}{w_1w_2w_3^{2n-1}} .
$$
At genera two and three the independently integrated full local vertices give the coefficients $-N/720$ and $21600N J_3$, respectively. At higher genera the Li–Wang flat worldsheet solution supplies the doubly soft primary traces and the lower-genus hypotheses. The pole argument applies at every $n\ge2$. In the linear state dictionary the indicated vector Laurent mode without its prefactor maps to $d(e^{z_1+z_2})(d\theta_3)^{2n-2}/(2n-2)!$, an exact form of degree $2n-1$ in the single-graviton spectrum. Its geometric degree specifies the output channel. Matching a boundary matrix element requires the transferred quantum operation. [Planar sources, proof and complete low-genus checks](../iib-bcov/README.md#an-elementary-scalar-period-at-each-genus)

At genus four, summing all fifty flat graphs fixes the entire flux-linear local contact on this nine-scalar source family. Its reduced coefficient multiplying $N e^{z_1+z_2}(\partial_{z_2}-\partial_{z_1})/(w_1w_2w_3^7)$ is
$$
8!\frac{9072}{11}J_4
=-\frac{6436431208614965917972531}{34223998208400000000}
-\frac{1593985240704793516397}{3168888723000000}\zeta(3)
+\frac{96133605459532151}{125749552500}\zeta(5)<0 .
$$
The $8!$ comes from differentiating $w_3^{-1}$ eight times before the two lifted source monomials supply $w_3^2$. This is a local consequence of the solved planar restriction. The unrestricted genus-four flat vertex is given by the worldsheet integral. Its full diagram evaluation and the transfer on the punctured background remain open. [Four-loop period and marked scalar](../iib-bcov/README.md#a-four-loop-measurement-on-the-elementary-scalar)

Multiplication by the BCOV descendant parameter makes this vector exact. The tangent bivector with the same scalar Laurent coefficient is closed under $\bar\partial+N\{F,-\}$, and its divergence gives the displayed tangent vector. Its full linearized differential is therefore $t$ times the vector. The nonzero class at $t^0$ survives precisely this first test and dies upon inverting $t$. [Seven- and nine-slot replay on shifted momenta](code/flux_primary.py)

At genus three one can also attain the highest normal pole of this single-flux vertex. Use two copies of each scalar $e^{w_1},e^{w_2}$ and the three further scalars $w_1 e^{z_1},w_2 e^{z_2},(w_1+w_2)e^{z_1+z_2}$. Introduce slopes $s,t,u$ on their three linear normal factors by differentiating scalar exponentials as above. These three inputs admit flux-closed lifts, since the coefficient of the flux action is a multiple of $w_1 f$ or $w_2 f$ and the correcting tangent vector is perpendicular to the tangential exponential momentum. The corrections wedge to zero at the flux insertion. The positive normal powers from undifferentiated exponentials in $w_1,w_2$ vanish in the normal Laurent quotient. The surviving coefficient is the $stu$ coefficient of the trace. Each five-row minor contains exactly one power of $q_3$ if nonzero. Its four-row complement contains at most three mixed inputs and requires two tangent rows, hence is affine in $s,t,u$. The integrated trace contains six such minors. Polarize its degree-three part by inclusion-exclusion on the seven nonempty subsets of the three slopes and extract the odd cubic on each ray. Two interpolation stencils and independent normal momenta give the coefficient
$$
q_3^6\left(\frac{20401023829}{105921623808000}
-\frac{6463627}{40864824000}\zeta(3)\right).
$$
Acting with six derivatives on the unit flux supplies $6!/(w_1w_2w_3^7)$. Taking the divergence supplies $2e^{2z_1+2z_2}(\partial_{z_2}-\partial_{z_1})$. Thus the reduced output in the linear dictionary is the coefficient in parentheses times $N\,d(e^{2z_1+2z_2})(d\theta_3)^6$. Its seventh-order pole is maximal because the flux momentum can enter each determinant only linearly even when that determinant contains the balancing slot. This component is nonzero by the irrationality of $\zeta(3)$. [Complete eighteen-coefficient evaluation](code/flux_primary.py)

### A ninth-order normal pole at four loops

Let $z_1,z_2$ be the brane coordinates and $w_1,w_2,w_3$ its normal coordinates. For nine labeled scalar inputs take $w_1^2e^{2z_1}/2!$, six copies of $w_1e^{z_1}$ and two copies of $w_2^4e^{z_2}/4!$. Each scalar vanishes on the brane. For a scalar $w_j^k H(z)$ its flux bracket has the Dolbeault primitive built from $w_j^{k-1}b_j$ times the tangent Hamiltonian vector of $H$. The correcting tangent vector wedges to zero with the explicit tangent bivector. The eleven-primary local trace gives, in reduced units,
$$
N\,8!\left(\frac{2625}{2}K_4+50400J_4\right)
\frac{e^{8z_1+2z_2}(8\partial_{z_2}-2\partial_{z_1})}{w_1w_2w_3^9}.
$$
Both periods are strictly negative. The ninth-order pole is maximal because each of the eight momentum determinants is at most linear in the flux momentum after conservation. Modulo the descendant parameter, the flux action on a holomorphic scalar can only give simple normal poles. Thus this contact reaches a nonzero class of the linearized target. Its coefficient in the linear dictionary on $N\,d(e^{8z_1+2z_2})(d\theta_3)^8$ is $2625 K_4/2+50400J_4$. The flat worldsheet interaction solves the unrestricted bulk equation. Its transfer with separated vertices remains open.

Rescale tangent momenta in the first seven slots by $\alpha$ and the other two by $\beta$. The symbol acquires $(\alpha\beta)^8$, and the tangent differential becomes $d(e^{8\alpha z_1+2\beta z_2})$. Coefficient extraction gives $d$ of every tangent monomial. Normal $SL_3$ rotations of $(d\theta_3)^8$ span the $45$ independent symmetric eighth powers. Thus the linear span of these local outputs contains $d\mathbb C[z_1,z_2]\otimes\operatorname{Sym}^8\langle d\theta_1,d\theta_2,d\theta_3\rangle$ in the linear state dictionary. [Independent eleven-slot determinant sums](code/flux_maximal_minors.py)

The normal component of its variation under the lifted physical special conformal generator $K_2^N$ is
$$
8N\,8!\left(\frac{2625}{2}K_4+50400J_4\right)
\frac{e^{8z_1+2z_2}\,\partial_{w_3}}{w_1w_2w_3^8}.
$$
Indeed the bracket of $K_2=z_2(E_z-E_w)$ with a tangent vector of second component $V^2$ has normal component $V^2E_w$. Only $w_3\partial_{w_3}$ survives in the displayed Laurent mode. Transforming a scalar argument at zero flux still yields a tangent output for the shifted trace; the unshifted eleven-field trace has no operation on nine inputs. A strictly equivariant transferred operation must compensate this coefficient. [Lifted conformal generators](#the-two-physical-special-conformal-generators)

At genus $n\ge2$, use $2n-2$ copies of $w_1e^{z_1}$ and two copies of $w_2^n e^{z_2}/n!$, keeping the first input. Give the first $2n-1$ waves generating covectors $s_i\,dw_1+h_i\,dz_1$, with $h_0=2$ and the remaining $h_i=1$. Let the other two covectors be $u\,dw_2+dz_2,v\,dw_2+dz_2$. The linear map from the fivefold planar frame has determinant $(u-v)q_3$ in coordinate order $(w_1,z_1,w_2,z_2,w_3)$, where $q_3$ is the $dw_3$ flux momentum. This multiplies the trace by $(u-v)^{2n}q_3^{2n}$. Extract $s_0^2\prod_{i=1}^{2n-2}s_i$ and $u^nv^n$. Terms differentiating the exponentials in place of the symbol leave $w_1$ or $w_2$ and vanish in the normal Laurent quotient.

For the planar coefficients put $L=2n-2$ and divide by $(-2)^{n-2}(L-3)!!$, using $(-1)!!=1$. Each ordinary matched pair contributes $-2$ from its squared bracket. For $A,D$ the slot with two normal derivatives must lie in the distinguished triple. If its two other tangent strengths are $h_j,h_k$, the three possible centers give respectively $2(2h_0^2-h_j^2-h_k^2)$ and $h_0(h_j+h_k)-2h_jh_k$ after removing $h_jh_k$. All other strengths are one. The $\binom L2$ triples therefore give $A=6L(L-1)$ and $D=L(L-1)$ in these units. For $B$, pair the appended total with the distinguished slot or one of the $L$ other slots. The first case gives $-2(L-1)L^2$. In the second case the factor from the appended bracket is $-2(L+1)$, and its other normal derivative may hit the distinguished slot, contributing $-4(L-1)$, or one of the $L-1$ remaining slots, contributing $L-1$. Together these give $B=2L(L-1)(2L+3)$. Hence $3B=(4n-1)A$, and the undetermined trace $2(3n-1)A-5B+2(2n+1)D$ vanishes on this coefficient.

Multiplication by $[u^nv^n](u-v)^{2n}=(-1)^n\binom{2n}{n}$ gives the reduced trace coefficient of $q_3^{2n}$,
$$
\binom{2n}{n}\frac{(2n-2)!}{n!}
\left(\frac{4n-1}{n!}K_n+n!J_n\right).
$$
The normal derivative contributes $(2n)!$, and tangent divergence acts on $e^{2n z_1+2z_2}$. At four loops the planar coefficients before the $u,v$ extraction are $(2160,10800,360)$. At two and three loops the independent full seven- and nine-field vertices reproduce the formula on both period channels where present. Above genus three it evaluates the saturated planar restriction of the Li–Wang flat solution. [Exact coefficient extraction and five-dimensional checks](code/flux_maximal_planar.py)

### Three tangent-dependent scalars suffice for the highest pole at every genus

At genus $g\ge2$, let $w_1,w_2,w_3$ be normal brane coordinates and $z_1,z_2$ tangent coordinates. Use $g-1$ scalar sources $e^{w_1}$, $g-1$ sources $e^{w_2}$, and the three further regular sources
$$
\frac{w_1^2}{2!}e^{z_1},\qquad w_2e^{z_2},\qquad w_2e^{A z_1+B z_2}.
$$
Here $A,B$ specify the tangent momentum of the last source. Every tangent-dependent input vanishes on the brane and has a flux-closed lift as above. In the fourfold let $e^1,\ldots,e^4$ be the coordinate covectors, take them as four hard waves with a fifth balancing wave, and add $g-2$ soft waves of momentum $\epsilon(e^1-Ae^2-Be^3)$, with the balancing wave recoiling. In the reduced graph convention, the coefficient of the highest normal pole in the shifted fivefold primary trace, before the flux derivatives and tangent divergence, is
$$
\frac{4A^2B((g-1)!)^2}{(g-2)!}\,
[\epsilon^{2g-4}]V^{(4)}_{g;0}.
$$
For $g=2$ the indicated fourfold vertex has just the five hard waves. For $g>2$ its coefficient is strictly positive by the fourfold soft theorem for every real $A,B$. Thus this trace coefficient has the sign of $B$ whenever $AB\ne0$. In the primary vector variation multiply it by
$$
N(2g)!\,
\frac{e^{(1+A)z_1+(1+B)z_2}
\big((1+A)\partial_{z_2}-(1+B)\partial_{z_1}\big)}
     {w_1w_2w_3^{\,2g+1}} .
$$
Here $N$ is the D3 flux. At $(A,B)=(1,1)$ the reduced trace coefficients are $1/15552$ at genus two and $740466719/7565830272000-(469153/5837832000)\zeta(3)$ at genus three. Both were evaluated from the full fivefold graph traces and from the independent fourfold formulas. The source-dependent heat factor is $(-1)^g\lambda^{2g}/(2\pi i)^{5g}$. [Two independent vertex evaluations and shears](code/flux_three_tangent_sources.py) · [Fourfold primary normalization](../iib-bcov/README.md#all-genus-primary-vertex-from-the-worldsheet)

To derive the coefficient, put $s\,dw_1$, $t\,dw_2$, $u\,dw_2$ on the three generating waves and take $[s^2tu]$. A shear fixing all normal covectors removes the first two normal shifts and leaves $x\,dw_1+y\,dw_2$ on the third wave, with $x=-As$, $y=u-Bt$. All scalar momenta lie in $dw_1,dw_2,dz_1,dz_2$. Thus the normal flux momentum enters only by $q_3^{2g}$, where $q_3$ is its $dw_3$ component. Set $q_3=1$. The positive worldsheet trace $V^{(5)}_{g;0}(x,y)$ has separate degree at most two in $x,y$. For $x=1/\epsilon$ multiply the first normal coordinate of every momentum by $\epsilon$. The determinant character supplies $\epsilon^{2g}$, whereas each of the $g-1$ pure first-coordinate waves supplies at least two soft powers. Their difference gives the bound two; use the second normal coordinate for $y$. Hence $[s^2tu]$ of the negative worldsheet trace is $2A^2B[x^2y^2]V^{(5)}_{g;0}$. Normal powers left on the input exponentials annihilate the simple $w_1,w_2$ poles. The tangent-vector corrections in the source lifts wedge to zero with the explicit tangent flux bivector.

Now take $x=y=1/\epsilon$ and scale both normal coordinates by $\epsilon$. After extracting both sets of $g-1$ soft columns, the hard momenta consist of $dz_1,dz_2,dw_1+dw_2+A\,dz_1+B\,dz_2,dw_3$ and their balancing wave. Both normal hard period rows are the Abel difference between the third hard mark and the balancing mark. Fix the five hard marks and curve. In a period-orthonormal basis the covariance of a holomorphic evaluation integrated over a forgotten mark is the identity. Integrating the $r$ labeled evaluation columns of one soft direction in the leading squared determinant gives $r!$ times the remaining hard polarization power $\Theta^{g-r}/(g-r)!$ in that direction. Indeed the covariance of the two alternating $r$-row minors is $r!\,\delta_{IJ}$ for increasing index sets $I,J$. The three other target directions contribute their common top powers. In five dimensions the two sets with $r=g-1$ therefore give $((g-1)!)^2\Theta^2$ of the same hard Abel difference. In four dimensions a single set of $r=g-2$ gives $(g-2)!\Theta^2/2!$. In the fourfold basis of the hard simplex the soft direction becomes $e^1-Ae^2-Be^3$. Taking the ratio and integrating on $\mathcal M_{g,5}$ proves the formula. Polynomial interpolation and absolute convergence justify coefficientwise integration; diagonal coincidences have zero fiber measure.

Three is the minimum number of tangent-dependent regular scalar sources for a pole of order $2g+1$ on flux-closed lifts. If the two tangent momenta are independent, a determinant-one shear removes all normal momenta from both mixed waves, fixing the pure normal waves and the flux. The pure-trace symbol and its tangent divergence are therefore independent of both removed normal shifts. If the tangent span has dimension at most one, the momentum span has dimension at most four and the trace vanishes. A flux-closed regular scalar restricts to a constant on the brane, since its tangent Hamiltonian derivative times the simple normal residue must be Dolbeault-exact. After dropping that constant from the differentiated vertex, each mixed input supplies a factor in the ideal $(w_1,w_2,w_3)$. With two mixed inputs these factors remain on the output, and its normal pole degree drops by at least two. This argument works on Taylor coefficients of polynomial sources by polarization. For $A=B=1$, independent rescaling of the two tangent momenta and normal $SL_3$ rotation give the full span $d\mathbb C[z_1,z_2]\otimes\operatorname{Sym}^{2g}\langle d\theta_1,d\theta_2,d\theta_3\rangle$ in the linearized state dictionary. Multiplication by the descendant parameter gives a boundary, by divergence of the tangent bivector with the same scalar Laurent coefficient. The full operation on punctured-space cohomology also contains separated vertices.

These local measurements recover the repeated-soft fourfold polynomial from a bulk contact. The five hard simplex waves sum to zero. Expand a soft direction in them; its five real coefficients are defined up to a common shift. The coefficient $[\epsilon^{2g-4}]V^{(4)}_{g;0}$ is symmetric in these coefficients and homogeneous of degree $2g-4$. The three-source tests above give it whenever two of the five coefficients agree: permute them to the final two positions, shift those to zero, then use homogeneity on $(1,-A,-B,0,0)$ and polynomial continuation from $AB\ne0$. A symmetric polynomial vanishing on this union of hyperplanes is divisible by every squared pairwise difference, by reflection in each hyperplane, hence by their product of degree twenty. Consequently these tests determine the complete repeated-soft fourfold polynomial for $2\le g\le11$. At higher genus they determine it modulo this discriminant. Numerical or period evaluation of the higher-genus measurements is a separate integration problem.

#### Two polynomial insertions determine the four-loop highest pole on this family

At genus four, extract $[A^2B]$ in the third mixed source, giving $w_2z_1^2z_2/2!$. Write $C_4(v)$ for the coefficient of $\epsilon^4$ with two soft waves in direction $v$ on the fourfold simplex specified above. In the reduced trace the insertion gives
$$
72C_4(e^1)=-\frac{11}{16}K_4-36J_4-36W,\qquad W=V(1,1,1).
$$
Indeed the soft waves give the collinear strengths $(1,\epsilon,\epsilon)$. The coefficient of $\epsilon^4$ in this fourfold vertex is $W/9+V(1,1,-2)/36$. The centered fourfold relation and $0<W<-K_4/64-9J_4/11$ give the strict bounds
$$
0<-\frac{K_4}{8}-\frac{72J_4}{11}
<72 C_4(e^1)<-\frac{11K_4}{16}-36J_4.
$$
The common vector factor for this insertion is $N8!e^{z_1+z_2}(\partial_{z_2}-\partial_{z_1})/(w_1w_2w_3^9)$. A second insertion takes $A=-1$ and extracts $[B]$, giving $w_2z_2 e^{-z_1}$. Its stripped trace is $72C_4(e^1+e^2)$ and its vector factor is $-N8!e^{z_2}\partial_{z_1}/(w_1w_2w_3^9)$. Both fourfold coefficients are strictly positive. [Collinear periods and their bounds](../iib-bcov/README.md#a-four-loop-measurement-on-the-elementary-scalar)

The complete fifty-graph sum fixes the axial ninth-order pole,
$$
72C_4(e^1)=
\frac{1141633397798354170399423}{2593214569701756000000000}
+\frac{141363079741091006137}{120056230078785000000}\zeta(3)
-\frac{38974345636001537}{21778907950800000}\zeta(5).
$$
It lies strictly between $112386871456/10^{20}$ and $112386871457/10^{20}$. The mixed insertion still requires its own bare four-loop graph sum. The strict comparison below gives $10^{20}\,72C_4(e^1+e^2)>18731145242$. [Exact integrals and axial conversion](../kodaira-spencer/four_loop_arakelov_periods.json) · [Mixed Ward equation](../iib-bcov/README.md#a-ward-projection-for-the-mixed-genus-four-d3-pole)

**A boundary test of the mixed ninth-order output.** Restrict every component of the boundary superfield to a formal Cartan-valued profile depending on one common tangent coordinate. Write $v_{\dot a}$ for its tangent derivative covector in the paper's dotted spinor notation. The classical variation of every jet vanishes by the Cartan commutator. At positive loop order $n$, the two output superfields in every Taylor term of the boundary fan have $n$ explicitly paired dotted derivatives. All their additional displaced derivatives act on the same profiles. The paired factor after restriction is
$$
(v_{\dot a}v^{\dot a})^n=0 .
$$
This annihilates the image of each $Q_n$ on **every** boundary cochain, including arbitrary derivatives on the inputs and products of traces. With independent endpoint waves of tangent momenta $(p_1,\tau p_2)$ and $(q_1,\tau q_2)$ the contracted endpoint factor is exactly $\tau^n(p_1q_2-p_2q_1)^n$ in the orientation of the replay. Thus the vanishing can also be tested to order $n$ transverse to the one-coordinate profiles. [Full Taylor kernel](../../sources/LLC_main.tex#L256) · [Classical and component conventions](../../sources/LLC_main.tex#L4771)

The linear dictionary sends the mixed local contact to $72N C_4(e^1+e^2)d(e^{z_2})(d\theta_3)^8$ in reduced units. At $z=0$ the corresponding boundary component of $\operatorname{Tr}(dC)^9$ is $9ig\,\operatorname{Tr}((\partial_{z_2}c)(\gamma^3)^8)$. Set $\partial_{z_2}c=pH$, $\gamma^3=xH$, with $p$ odd, $x$ even, and $H=\operatorname{diag}(1,1,-2,0,\ldots,0)$ in the defining representation of $SU(N)$, $N\ge3$. These are one-coordinate Cartan profiles, and the boundary component evaluates to
$$
9ig\,p x^8\operatorname{Tr}(H^9)=-4590\,ig\,p x^8 .
$$
Therefore this particular local output, with the linear state dictionary, cannot itself be an image of any ordinary boundary $Q_n$ or their sum. At $SU(2)$ the same cochain vanishes identically by Cayley–Hamilton. At $SU(3)$ it is a nonzero product of lower trace components: for traceless matrices $A,B$ with $B$ even and $p_j=\operatorname{Tr}B^j$,
$$
\operatorname{Tr}(A B^8)=\tfrac14 p_2^2 p_3\operatorname{Tr}(AB)
+\left(\tfrac18p_2^3+\tfrac19p_3^2\right)\operatorname{Tr}(AB^2).
$$
This follows by polarizing $\operatorname{Tr} B^9=3p_2^3p_3/8+p_3^3/9$. For $SU(N)$ the invariant $\operatorname{Tr}B^9$ is an independent generator precisely at $N\ge9$. The same Cartan test on $dz_2(d\theta_3)^{2g}$ gives $(2g+1)ig[2-2^{2g+1}]p x^{2g}$ for every positive genus $g$. A proposed intertwining map for the ordinary differential must account for this nonzero evaluation through the remaining punctured-bulk or dictionary terms. [Exact endpoint contractions and polarized trace](code/finite_shear_states.py) · [One-particle dictionary](../../code/holography/derham.py#L86)

The axial insertion gives a fully evaluated test of the required correction. At the origin the linear state dictionary takes the numerator $8!(\partial_{z_2}-\partial_{z_1})$ over $w_1w_2w_3^9$ to $(dz_1+dz_2)(d\theta_3)^8$. On the same $z_2$-dependent Cartan fields, its $dz_1$ term vanishes. If a completed operation is an image of the paper’s boundary differential, its remaining transfer and dictionary terms together must give $+4590\,igNpx^8\,72 C_4(e^1)$ on this test. In the reduced heat units, division by $igNpx^8$ here means extracting its coefficient. That exact coefficient is
$$
\frac{1141633397798354170399423}{564970494488400000000}
+\frac{141363079741091006137}{26156041411500000}\zeta(3)
-\frac{38974345636001537}{4744860120000}\zeta(5).
$$
This is a conditional matching equation for the omitted terms. [Check of the dictionary signs and the three coefficients](code/finite_shear_states.py)

For a direction $v$, subtract the mean of its five simplex coefficients and denote their power sums of orders two and four by $p_2,p_4$. A symmetric degree-four polynomial on this centered hyperplane is spanned by $p_2^2,p_4$. Evaluating on the two directions of the insertions gives the full answer
$$
C_4(v)=\frac{30p_4-7p_2^2}{8}\,C_4(e^1)
       +\frac{13p_2^2-20p_4}{12}\,C_4(e^1+e^2).
$$
Both multipliers are nonnegative for every real $v$. Fix $p_2=1$. At a constrained extremum each of the five centered coefficients solves the same depressed cubic. With two distinct roots the multiplicities $(1,4),(2,3)$ give $p_4=13/20,7/30$. With three distinct roots their unweighted sum is zero. Their multiplicities $(3,1,1)$ force the tripled root to vanish, giving $p_4=1/2$. Multiplicities $(2,2,1)$ force the lone root to vanish, giving $p_4=1/4$. This proves the bounds and gives the exact relation $C_4(e^1-e^2)=4C_4(e^1)+C_4(e^1+e^2)$. [Exact interpolation](code/flux_three_tangent_sources.py)

#### Complex tangent momenta at which the local pole vanishes

Analytically continue the three regular scalar sources in $A,B$. Their fourfold soft direction has five simplex coefficients $(1,-A,-B,0,0)$. The soft coefficient is a polynomial in these momenta, so each real polynomial identity above extends to complex momenta. The continuation acts on the momentum in both factors of the squared determinant.

Take a primitive cube root of unity $\omega$ and choose $(A,B)=(-\omega,-\omega^2)$. A cyclic permutation of the first three simplex coefficients multiplies their tuple by $\omega$. Permutation symmetry and homogeneity of the leading soft coefficient give
$$
[\epsilon^{2g-4}]V^{(4)}_{g;0}
=\omega^{2g-4}[\epsilon^{2g-4}]V^{(4)}_{g;0}.
$$
Thus the local D3 contact has zero coefficient at pole order $2g+1$ whenever $g\not\equiv2\pmod 3$. Its tangent divergence has both coefficients $1-\omega,1-\omega^2$ nonzero. At genus two the stripped pole coefficient equals $4A^2B/62208=-\omega/15552$. The number of pure normal scalar inputs is $2g-2$ at each genus. At $g=5,8,\ldots$ this symmetry allows a coefficient whose value requires the corresponding period.

At genus four the two independent period multipliers vanish together at exactly eight ordered pairs $(A,B)$ over $\mathbb C$. Write $p_2,p_4$ for the power sums of the five centered coefficients. The determinant of their two coefficient rows in the basis $(p_4,p_2^2)$ is $125/48$. Thus the conditions are $p_2=p_4=0$. The five centered coefficients sum to zero and two of them have the common value $r=-(1-A-B)/5$. Newton's identities reduce their monic root polynomial to $X^5-e_3X^2-e_5$, where $e_3,e_5$ are their third and fifth elementary symmetric functions.

* If $r=0$, the three nonzero entries are $1,\omega,\omega^2$ up to exchanging the last two.
* If $r\ne0$, the repeated root forces $e_3=5r^3/2$ and $e_5=-3r^5/2$. On substituting $X=r(1+z)$ the polynomial becomes $r^5z^2(2z^3+10z^2+20z+15)/2$. For each ordering $(u,v,w)$ of the three roots of this cubic, the remaining pair is $(A,B)=(-v/u,-w/u)$. The cubic has distinct, nonzero roots, with sum $-5$, so their six normalized orderings are distinct.

The genus-three leading vertex on these same eight tangent momenta is zero in each of the rational and $\zeta(3)$ channels. Interpolating the fourfold six-wave trace at six integer directions gives $C_3(v)=(5p_2/4)C_3(e^1)$ in both channels, with an extra direction as a control. All eight have nonzero genus-two trace $A^2B/15552$. The quartic integral can have additional zeros depending on the periods.

On the tangent-momentum line $B=1-A$, the three nonzero simplex coordinates sum to zero, so $p_2=2(A^2-A+1)$ and $2p_4=p_2^2$. The local pole before tangent divergence and flux differentiation in genera two through four is therefore
$$
\begin{aligned}
g=2:\quad &\frac{A^2(1-A)}{15552},\\
g=3:\quad &40 A^2(1-A)(A^2-A+1)\,C_3(e^1),\\
g=4:\quad &72 A^2(1-A)(A^2-A+1)^2
            \bigl(4C_4(e^1)+C_4(e^1+e^2)\bigr).
\end{aligned}
$$
The positive real periods make the zeros at $A=-\omega,-\omega^2$ exactly simple at genus three and exactly double at genus four. For real probes the ratio of these stripped poles at $A=1/3$ and $A=1/2$ is $(16/27)(28/27)^{g-2}$ for $g=2,3,4$. More generally the symmetric polynomial on this line belongs to the ring generated by the elementary symmetric functions of degrees two and three of the triple with zero sum. At degree $2g-4$, each monomial contains at least $(g-2)\bmod 3$ powers of $A^2-A+1$: the power of the degree-three generator must be even. This gives the corresponding lower multiplicity at all genera. [Symbolic elimination and graph checks](code/flux_three_tangent_sources.py)

#### Strict comparison of the two soft directions

Let $E$ be the real hyperplane of five centered simplex coefficients. At fixed curve and seven marks the coefficient of $t^2$ in the signed Betti determinant for two soft waves $tv,tv$ is a quadratic in $v\in E$. The two soft marking columns each supply one $t$; all other columns take their hard value. Polarize this quadratic before integrating its squared norm. The result is a positive semidefinite Hermitian form on $\operatorname{Sym}^2(E_{\mathbb C})$ whose value on $v\otimes v$ is $C_4(v)$. Permuting the five hard marks preserves it, since the induced change of their four-dimensional momentum frame has determinant $\pm1$.

Under $S_5$ the symmetric square splits without multiplicities into the trivial representation, the four-dimensional standard representation, and the five-dimensional representation of partition $(3,2)$. Indeed the ten symmetric pairs of distinct labels form a basis for the symmetric square of the standard representation: if $p_i$ are the five simplex vectors, use $p_i^2=-\sum_{j\ne i}p_i p_j$ from $\sum_i p_i=0$. An equivariant matrix on the ten pairs has three free entries, according as two pairs coincide, meet in one label, or are disjoint. Removing the trivial and standard summands therefore leaves one irreducible summand. In the Euclidean metric of the centered hyperplane, write $P=I-\mathbf1\mathbf1^{\mathsf T}/5$ for its orthogonal projection. The trivial tensors are multiples of $P$ and the standard tensors are $P\operatorname{diag}(c)P$ for centered $c$. Their squared norms on projecting $v\otimes v$ are
$$
\frac{p_2^2}{4},\qquad
\frac53\left(p_4-\frac{p_2^2}{5}\right),\qquad
\frac{13p_2^2-20p_4}{12}
$$
in the three summands, respectively. The middle factor follows from $\operatorname{tr}((P\operatorname{diag}(c)P)^2)=3\lVert c\rVert^2/5$. Thus the projection norms squared at $e^1$ are $(4/25,12/25,0)$ and at $e^1+e^2$ they are $(9/25,2/25,1)$. With $\lambda_1,\lambda_4,\lambda_5$ the eigenvalues of the integrated Hermitian form in this metric,
$$
C_4(e^1+e^2)-\frac16 C_4(e^1)=\frac{\lambda_1}{3}+\lambda_5>0 .
$$

For strictness use the [torsion-curve multiplication determinant](../iib-bcov/README.md#a-four-loop-measurement-on-the-elementary-scalar) on $y^3=x^5-1$. Take the seven marks and sixteen evaluation points of its replay in the same order. Give the first five marks the hard momenta $p_i=e^i$ for $i\le4$ and $p_5=-p_1-p_2-p_3-p_4$. The two soft marks carry $tv,tv$ and the fifth hard mark recoils by $-2tv$. Vanishing of an eigenvalue would force the corresponding polarized determinant component to vanish at every marked curve by positivity and continuity. Complexify the conjugate coordinates of the first four hard marks and continue each around cycles at fixed curve. The Betti derivative $dw-(d\tau)(\tau-\bar\tau)^{-1}(w-\bar w)$ then adds a lattice of holomorphic differentials to each hard row. The four hard momenta there are the same basis for every $v$ and every $t$. Thus the shifts in the test are simultaneous in every compared determinant. Their lattice is Zariski dense in the four copies of $H^0(K)$, the space of holomorphic one-forms, so an identically vanishing component would vanish at the replay's common algebraic shift.

Write $D(v)$ for the signed coefficient at that shift, extracted in the same evaluation frame before taking any norm. The four columns corresponding to the changed momentum direction bound the determinant degree by four in $t$. Over $\mathbb F_{211}$ the five $D(p_i)$ are $21,91,144,205,120$. Their sum is $159$ and the first difference is $141$. They detect the trivial and standard components, respectively. In the ordered pair directions $(1,2),(3,4),(1,3),(2,4)$, $D(p_i+p_j)$ gives $172,77,180,155$. The sum of the first two minus the last two is $125$. Any quadratic in the trivial and standard components takes the form $\sum_i c_i v_i^2$ on centered coordinates; its values on these four pairs cancel in this difference. This detects the five-dimensional component. All reductions use good denominators and evaluate the same characteristic-zero determinant polynomial. Consequently all three eigenvalues are strictly positive. The abstract $S_5$ average of the squared rank-two determinant $(v_1-v_2)(v_1+v_2)$ is supported entirely on the standard summand and attains the factor $1/6$; strictness here uses the curve calculation. — Codex. [Projection norms, common shift and finite-field determinants](code/soft_quartic_gram.py)

#### A separate five-mark period for the second insertion

Substitution of the four-hard diagram formula for $W$ gives the first stripped pole directly in the integrated period channels:
$$
72C_4(e^1)=-\frac{K_4}{8}-36J_4-72D_4
             +18\bigl(X_4-\langle\epsilon+\varphi\rangle\bigr).
$$
Here the periods and the graph averages use the conventions of the [four-hard diagram formula](../iib-bcov/README.md#four-hard-diagram-formula-for-the-remaining-elementary-scalar-period).

To separate the second insertion, label the seven momentum rows from zero and denote their four-row determinants by brackets. Sum the following monomial over its 1,260 distinct relabelings:
$$
[0123]^2[0456]^2[0124]^4 .
$$
Each row occurs at least twice. The first two determinants together use all seven rows, so any three collinear rows force one of them to vanish. The polynomial is thus invisible on the entire triple-collinear restriction fixing $W$. On two equal soft waves $tv,tv$ the degree-four coefficient uses only the 180 relabelings in which each soft row occurs once before squaring. At $v=e^1$ it is zero. At $v=e^1+e^2$ precisely 48 terms give one and the rest vanish. The balancing hard recoil changes only higher orders. This symmetric polynomial has the required determinant weight and quadratic softness. A relation fixing the mixed measurement requires information beyond these constraints and triple-collinear data.

There is an explicit five-mark integral for that additional datum. Label the hard simplex by $1,\ldots,5$, with the fifth wave balancing, and write $\Theta_4(z_i-z_j)$ for the one-coordinate genus-four Betti form with charges $+1,-1$ at these two marks. For the mixed ray take projected normal hard momenta $q,-q,r,s,-r-s$ in an integral unit frame and longitudinal hard charges $0,1,0,0,-1$. These follow from replacing the first two target coordinates of the hard momenta by their difference and the second coordinate; $e^1+e^2$ becomes the unit longitudinal direction. At order $t^4$ both soft marks enter by their respective curve evaluation columns. Integrating their two outer squares by Riemann's bilinear identity leaves two longitudinal exterior factors and a multiplicity $2!$. Hence
$$
C_4(e^1+e^2)=2!\int_{\mathcal M_{4,5}}
\frac{\Theta_4(z_2-z_5)^2}{2!}
\frac{\Theta_4(z_1-z_2)^4}{4!}
\frac{\Theta_4(z_3-z_5)^4}{4!}
\frac{\Theta_4(z_4-z_5)^4}{4!}.
$$
The axial formula, after interchanging labels 1 and 2, has $\Theta_4(z_1-z_5)^4$ in place of $\Theta_4(z_1-z_2)^4$. The two real soft coordinates of each new mark force exactly one factor of $t$ in each of the holomorphic determinant and its conjugate. Thus the longitudinal recoil at the fifth hard mark does not affect the displayed coefficient. Absolute convergence of the unscaled volumes gives coefficient integration by finite polynomial interpolation.

This leaves only two retained marks after integration along the three independent curve fibers. Let $\pi$ forget the moving mark $x$ and put $\pi_*(\Theta_4(x-z)^4/4!)$ for integration over the smooth genus-four curve with fixed mark $z$. The forms extend across coincident points of a fixed smooth curve. Over $\mathcal M_{4,2}$ the mixed integrand is
$$
2!\,\frac{\Theta_4(z_2-z_5)^2}{2!}
\ \pi_*\!\left(\frac{\Theta_4(x-z_2)^4}{4!}\right)
\left[\pi_*\!\left(\frac{\Theta_4(x-z_5)^4}{4!}\right)\right]^2 .
$$
For the axial integral all three pushed factors are based at $z_5$. This reduction specifies a separate integration target for the second pole. The exact replay checks the orbit evaluations and the factor $2!$ against both original sixteen-row determinants. — Codex. [Determinant replay](code/flux_three_tangent_sources.py) · [Fourfold trace basis](../iib-bcov/operator/fourfold_quartic_basis.py)

#### A rational upper bound on the mixed pole coefficient

The curve pushforward of $\Theta_4(x-z)^4/4!$ is a strongly positive $(3,3)$-form, as follows by contracting its squared holomorphic four-form with the vertical tangent and integrating. Positivity and exchange of the two retained marks alone give no scalar Hölder upper bound for its mixed product. In a nine-dimensional quotient take the positive coordinate volume forms on three disjoint three-planes. The sum of the first two and the third each have zero cube, whereas the square of the sum times the third is positive. Symmetrizing over the two assignments preserves this discrepancy.

On the five-mark space write $\Theta_{ij}=\Theta_4(z_i-z_j)$. For a positive complementary product $S$ of Betti forms and smooth closed positive $(1,1)$-forms from the compact stable moduli space, the Green identity gives
$$
\int\Theta_{ij}\,S\ \le\ \frac12\int(\psi_i+\psi_j)S .
$$
Indeed $\Theta_{ij}=(\psi_i+\psi_j)/2+d\beta$. On a resolved graph face with normalized lengths $\lambda_e$ and node angles $\tau_e$, contraction with $\sum_e\lambda_e\partial_{\tau_e}$ annihilates $S$. This follows from $\sum_e\lambda_e\,df_e=0$ for the minimizing currents of each Betti factor, with $f_e$ half the squared current; smooth forms from the stable space have no node angles. The boundary restriction of $\beta$ therefore contributes $-(\sum_e\lambda_e f_e)(\sum_e d\tau_e)\wedge S$, using the currents for the selected pair. This is a nonpositive measure. The total-angle measure is the outward limit of $d^c\rho\wedge S$ on positive-height levels, where $\rho$ is the negative log norm of the total boundary section: positivity on each level follows from $d\rho\wedge d^c\rho\wedge S\ge0$, and Stokes gives convergence on smooth tests on the compact resolution. The cutoffs have bounded level mass as in the four-hard calculation. [Green identity and boundary values, (1.7) and Proposition 1.5](https://arxiv.org/abs/2610.21927)

Every cotangent line is nef. Add $\delta>0$ times a fixed ample class to choose a positive smooth representative of each shifted cotangent class. Replacing a Betti factor by its shifted cotangent average increases the integral. This permits iteration with positive remaining factors and then $\delta\downarrow0$ by polynomiality of the compact intersections. If all the remaining factors in $S$ come from the compact stable space, only a single separating node can contribute: more than one node would require additional angles, and the minimizing current on a single self-loop is zero. The separating edge carries current of absolute value one exactly when the two marks lie on opposite sides. Thus in this last integration the exact divisor to use is $(\psi_i+\psi_j-\sum\delta_{h:I})/2$. The sum runs over stable separating divisors, oriented by $i\in I$, $j\notin I$ with genus $h$ on the side of $i$. [Nefness of the cotangent classes, Corollary 4.2](https://arxiv.org/abs/2304.08568)

Keep one $\Theta_{12}$ to the last step. On $\overline{\mathcal M}_{4,5}$ this gives
$$
\begin{aligned}
C_4(e^1+e^2)
&\le
\frac{
\displaystyle\int
(\psi_2+\psi_5)^2(\psi_1+\psi_2)^3
(\psi_3+\psi_5)^4(\psi_4+\psi_5)^4
\left(\psi_1+\psi_2-\sum_{\substack{h,I\\1\in I,\ 2\notin I\\\mathrm{stable\ separating\ divisors}}}\delta_{h:I}\right)
}{2^{14}(4!)^3}
\\
&=\frac{379101134941}{39450992600678400}.
\end{aligned}
$$
With every factor replaced by its cotangent average the value is $1535473758853/157803970402713600$. For each separating divisor the marked cotangent monomials restrict to the two stable factors, with cotangent exponent zero at each new node. Applying the descendant recursion with either the smallest or largest positive exponent after the string and dilaton reductions gives both rational values. — Codex. [Exact compact intersection replay](code/soft_quartic_psi_bound.py) · [Descendant recursion](https://arxiv.org/abs/0708.0565)

#### Exact intersection with two Betti factors

Write $D_{ij}=(\psi_i+\psi_j-\sum\delta_{h:I})/2$, summing the stable separating divisors that separate $i,j$. For a product $S$ of compact cotangent classes of complementary degree, the two-factor integral is
$$
\int_{\mathcal M_{g,n}}\Theta_{ij}\Theta_{kl}S
=\int_{\overline{\mathcal M}_{g,n}}D_{ij}D_{kl}S
-\frac1{12}\sum_{\substack{\Gamma\,\mathrm{two\ parallel\ edges}\\i,j\,\mathrm{and}\ k,l\ \mathrm{separated}}}
\frac{1}{|\operatorname{Aut}\Gamma|}\int_{\overline{\mathcal M}_\Gamma}\xi_\Gamma^* S .
$$
Here $\xi_\Gamma$ glues the two components along their two pairs of branches. The automorphism exchanging the two edges has order two. The formula counts that symmetry once.

To see the correction, intersect first with general ample divisors to obtain a transversal surface. Polarization and multilinearity recover the compact cotangent tests. This argument can be made in local stack charts, with intersections divided by the chart degree. The Chern–Weil formula for the Abel section with its semipositive biextension metric identifies the two-form product with the limit of intersections of its divisor extensions on modifications. For different pairs polarize the square. The metric has toroidal singularities and the input pairs have integral charges; the normalization of the single-node trace is the preceding Green identity. On the surface only two-node degenerations require further modifications. If every current is fixed by conservation, its energy is linear in the node lengths and the single-node traces already account for it. For two parallel edges of positive lengths $x,y$ separating the marks, minimizing the energy of the two currents summing to one gives $xy/[2(x+y)]$. This is the only nonlinear case. [Biextension singularities and Chern–Weil, Theorems 5.20 and 6.3](https://arxiv.org/abs/2112.09007) · [Graph energy and boundary normalization, Proposition 1.5](https://arxiv.org/html/2610.21927v1)

For adjacent primitive rays $(a,b),(c,d)$ with $|ad-bc|=1$, put $u=a+b$, $v=c+d$. The change in the exceptional coefficient on inserting their sum has magnitude $1/[2uv(u+v)]$. Each blowup subtracts the product of the two jumps from the intersection. On the unit segment the normalized profile is $s(1-s)/2$. Across a subinterval of length $h$ the squared derivative error of its chord integrates to $h^3/12$. The inserted mediant splits a length $1/(uv)$ into $1/[u(u+v)]$ and $1/[v(u+v)]$; the decrease in the total error is $1/[4u^2v^2(u+v)^2]$. The mesh tends to zero, giving total correction $1/12$ when both pairs cross. Bounded remainders in the toroidal metric have zero divisorial jump. As a normalization check, on $\overline{\mathcal M}_{1,2}$ the compact square $\int D_{12}^2=1/24$ is canceled exactly by the graph correction $1/(12\cdot2)$, as required by the rank-one Betti form.

In the mixed genus-four integral keep two of the four factors $\Theta_{12}$ and replace the other twelve factors by their cotangent averages. The sharpest of the ten choices of retained factors gives
$$
C_4(e^1+e^2)\le\frac{19566009209}{2191721811148800}.
$$
The two distinct separating divisors in the compact product are evaluated on ordered three-component trees. A repeated node contributes the negative sum of its two branch cotangent classes. Choosing the first or last positive exponent in descendant recursion gives the same result for all ten choices. For target dimension one at genus two or dimension two at genus one, (1.9) itself has just two Betti factors; the identity then gives the complete descendant vertex. For conserved momenta $P$, substitute $\Theta_g(P)=-\sum_{i<j\,;\,\mathrm{marks}}(p_i\cdot p_j)\Theta_{ij}$ and divide its squared integral by $2!$ as in (1.9). For two opposite unit momenta at genus two the descendant levels $(3,0),(2,1)$ give $1/1920,1/640$. [Stable-graph sum and replay](code/soft_quartic_two_betti.py)

#### Scalar graph heights and compact cotangent tests

On a stable graph with edge lengths $\ell_e$, let $p$ be its conserved scalar vertex charges and $j$ any current with boundary $p$. The height of the Abel section is $h(\ell)=\min_{\partial j=p}\sum_{e\,;\,\mathrm{edges}}\ell_e j_e^2/2$. Its upper convex Newton body consists of all $a$ with $a_e\ge j_e^2/2$ for some such $j$. In the coordinates $\rho_e=\sqrt{2a_e}$ the circulation criterion gives one inequality per bond (a minimal separating edge cut) $C$,
$$
\sum_{e\in C\,;\,\mathrm{cut\ edges}}\rho_e\ge |p(C)|,
$$
where $p(C)$ is the charge on one side. This cut description uses scalar currents. For vector charges the energy constraint is the squared norm of the full edge current.

For positive edge variables $X_e$ the Laplace integral over the Newton body is
$$
F_\Gamma(X)=\int_N \prod_{e\,;\,\mathrm{edges}} X_e e^{-X_ea_e}\,da_e .
$$
For a monomial ideal on a regular toric chart, Aluffi's integral over the complementary Newton region computes the pushforward of the exceptional Segre class. In degree $d$ its kernel coefficient exceeds the corresponding Laplace coefficient by $d!$ (integrate the exponential kernel against $e^{-s}ds$ after $X\mapsto sX$). On the good region the resulting class is the pushforward of $\exp(-E)$ for exceptional divisor $E$. For one edge it is $\exp(-X_ep(C)^2/2)$. Powers from its full edge support are pushed down along the gluing map after removing one $X_e$ per edge and substituting the normal class $X_e=-\psi_h-\psi_{h'}$. Multiply by the marked-leg exponential $\exp(\sum_{i\,;\,\mathrm{marks}}p_i^2\psi_i/2)$ and divide by the graph automorphism order. [Monomial Segre formula, Theorem 1.1](https://arxiv.org/abs/1307.0830) · [Pushforward on the cone stack, Propositions 70 and 76](https://arxiv.org/abs/2410.22546)

To compute the coefficient with branch powers $k_e$ after this substitution, apply inclusion–exclusion over failed cuts. Only sets of cuts covering every edge survive. Each contributes its sign times the integral of $\prod_e \rho_e^{2k_e+1}/(2^{k_e}k_e!)$ on a bounded rational polytope. For a bundle of parallel edges separating charge $s$, the simplex integral gives $-|s|^{2d}\prod_e (2k_e+1)!/[(2d)!\prod_e 2^{k_e}k_e!]$, where $d=\sum_e(k_e+1)$. The general solver triangulates the cut polytope over the rationals and integrates barycentric monomials. [Capacity integrals](code/betti_flow_polytope.py) · [Decorated stable graphs and cotangent integrals](code/betti_power_graphs.py)

For the compact identification, choose a smooth metric on the relative curve at the marked sections. Write $A$ for the marked degree-zero divisor with coefficients $p_i$. Its regularized self-height uses a Green function with singularity $-\log|w_i|$. In Li–Wang's Green function the singularity is $+p_i\log|w_i|$. Subtracting both diagonal singularities in the same coordinate gives their compact cotangent potential $F_h=-\langle A,A\rangle'_\mu/2$ up to the bounded sum of the logarithms of smooth marked norms. The self-height theorem on the full versal marked family gives
$$
F_h=-h_\Gamma((-\log|t_e|)_e)+O(1)
$$
uniformly in the plumbing parameters $t_e$. Here $O(1)$ denotes a bounded function on the entire punctured chart. The Schur complement in the graph polynomial is the minimum current energy, including when the smoothing rates differ. Thus the compact reference $\sum_i p_i^2\psi_i/2$ carries a semipositive singular rational line metric with the required toroidal singularities. Tensor by smooth ample metrics and polarize its top Chern–Weil identity to insert compact divisor classes, expressing each cotangent class as a difference of ample classes. The non-pluripolar product assigns zero mass to the boundary. [Uniform height asymptotic, (2.3) and Theorem 6.5](https://arxiv.org/abs/1512.04862) · [Green potential and Betti curvature, (1.6)–(1.8)](https://arxiv.org/abs/2610.21927) · [Toroidal metrics and Chern–Weil, Proposition 3.11 and Theorem 5.20](https://arxiv.org/abs/2112.09007)

In a chart the multiplier ideal at integer level $m$ contains precisely the monomials with exponent $u$ such that $u+\mathbf1\in\operatorname{int}(mN)$. This follows by radial integration of the monomial against the toroidal metric; a supporting direction at equality gives divergence, while strict containment gives exponential decay. Contracting edges projects $N$ to the body of the contracted graph: solve the remaining vertex balances on a spanning tree of each contracted component. Localizing the ideal therefore gives the projected ideal, by projection of the interior and rounding the free exponents upward. Node permutations preserve it. This constructs the global monomial ideal on the Artin fan at each $m$, including strata with monodromy. Compute the finite-level exponential pushforward in equivariant Chow on a cone using the toric chart and its free finite-dimensional approximations. Cone restrictions determine an operational class on the Artin fan; pullback and proper pushforward give the decorated strata with the branch normal classes and graph symmetries used above. [Artin-fan Chow, §2.1 (8)–(10), Propositions 70 and 76](https://arxiv.org/abs/2410.22546)

Let $N_m$ be the ideal's Newton polyhedron divided by $m$. Rounding upward gives $N\subset N_m$ and $N_m+\mathbf1/m\subset N$. Minimal exponents have a uniform coordinate bound after rescaling: removing directed cycles leaves flows of magnitude at most the total positive charge. For independent coordinates apply this to each current separately. Coordinate projections obey the same bound. The alternating subtraction of all projections thus has support in one fixed box, since above the bound a coordinate constraint is redundant. The full-support polynomial moments converge by dominated convergence on that box. The nef approximants in the proof of Theorem 5.18 are the multiplier-ideal divisors along a dyadic subsequence plus a fixed compact divisor divided by $m$. The latter terms vanish by multilinearity and boundedness of the lower-degree moments. The same computation applies on a smooth projective alteration of the stable stack resolving the inverse image of the boundary. A pulled-back node is a monomial times a unit, so the Lipschitz potential pulls back with bounded error. The pulled-back base ideals and the upstairs multiplier ideals have the same limiting Newton body, with an error bounded by a chart-dependent constant divided by $m$, by the two inclusions above and the nonnegative exponent matrix. Both have bounded projected moments. Their limiting degrees agree; at each finite level the pulled-back base ideals satisfy the projection formula. Division by the degree of the alteration recovers the stack integral. This argument uses the limiting moments when the cover ramifies. [Approximating divisors, (4.2), Corollary 5.16, Theorems 5.18 and 5.20](https://arxiv.org/abs/2112.09007)

In particular every scalar-target descendant integral at integral charges is rational, at every genus. Only finitely many graphs with at most $g$ edges enter the top power $\Theta(p)^g/g!$; the capacities, triangulation, and ordinary cotangent intersections are all rational. The deletion–contraction proof below identifies the full scalar vertex at every positive genus with $[t^{2g}]\,\int \mathrm{DR}_g(tp)\prod_i\psi_i^{k_i}$ for every allowed tuple of marked levels. The scalar graph tests reproduce every cotangent distribution for $(1,-1),(2,-3,1)$ at genera two, three and four, and $(2,-3,0,1)$ at genera two and three. On $\overline{\mathcal M}_{4,5}$ retaining $\Theta_{12}^4$ or $\Theta_{35}^4$ gives respectively $1651193/611529523200$ and $1584203/684913065984$ before the fifth-factor improvements below. [Independent descendant formula, Theorems 1–2](https://arxiv.org/abs/1211.5273)

#### Scalar double ramification on cycles of parallel bundles

In Pixton's regularization let $r$ be the residue modulus, $u$ an edge residue and $f_r(u)=u(r-u)$ for $0\le u<r$, extended periodically. Put $m_e=k_e+1$, $d=\sum_e m_e$. In the part homogeneous of charge degree $2d$ at $r=0$, replace the normalized residue sum by the normalized integral over the circle-valued flows. Indeed, dilate $r$ and the charges by the same positive integer. Cutting the cycle coordinates at multiples of $r$ gives a finite disjoint collection up to lower-dimensional faces. On each piece the scaled sum is a Riemann sum of a polynomial. Equivalently the integral is the term of highest joint degree $2d$ in the normalized weighted lattice sum. Lower joint degrees cannot contribute charge degree $2d$ at $r=0$. On these graphs the residue sums on each ordering chamber are obtained by successive Faulhaber summation. [Regularized graph formula](https://arxiv.org/abs/1602.04705)

For a single cycle shift the cumulative charges $P_i$ so that $\min P_i=0$ and take $r>\max P_i$. Then $f_r((y-P_i)\bmod r)=-(y-P_i)^2+r|y-P_i|$ for $0\le y\le r$. At order $r^0$ in the normalized cycle integral, the term with no absolute value gives $\prod_i P_i^{2m_i}=0$ up to sign. In the linear-in-$r$ term, changing the sign of $y-P_j$ below $P_j$ supplies the only further constant. The result is
$$
(-1)^{d-1}\sum_{j\,;\,\mathrm{cycle\ edges}}2m_j
\int_0^{P_j}(P_j-y)^{2m_j-1}\prod_{i\ne j}(P_i-y)^{2m_i}\,dy .
$$
The capacity test on this cycle asks whether the intervals $[P_i-\rho_i,P_i+\rho_i]$ have a common point. In the Laplace integral the capacities have independent tails $\Pr(\rho_i\ge t)=e^{-X_it^2/2}$. Condition on which interval has the smallest right endpoint $y$. Integrating $X_j(y-P_j)\exp[-\sum_i X_i(y-P_i)^2/2]$ over $y\ge P_j$ and summing over $j$ computes $F_\Gamma$. Its full-support coefficient is exactly the displayed integral divided by $2^d\prod_i m_i!$. Normal-class substitution gives the remaining $(-1)^{d-|E|}$, the sign of the Pixton edge expansions.

For a bundle of $b$ parallel edges the continuum residue kernel is their normalized circular convolution. On $0<s<r$ it is a homogeneous polynomial of degree $2d$ in $s,r$. Differentiate $2m_e$ times on each factor inside the convolution. The derivative of $f_r^{m_e}$ has constant density $(-1)^{m_e}(2m_e)!$ off zero, a delta mass at zero canceling its mean, and derivatives of delta masses. Under normalized convolution the mean-zero distribution $1-r\delta_0$ satisfies $(1-r\delta_0)^{*b}=(-1)^{b-1}(1-r\delta_0)$. Thus the highest $s$ coefficient of the kernel is $(-1)^{d+b-1}\prod_e(2m_e)!/(2d)!$. Reflection $s\mapsto r-s$ fixes the term linear in $r$ as well. Modulo $r^2$, the entire bundle can therefore be replaced by $f_r(s)^d$ times $(-1)^{b-1}\prod_e(2m_e)!/(2d)!$. On an outer cycle the only possible denominator is $r$ and all integration endpoints are polynomial in $r$, so the discarded terms have zero constant. In the capacity integral the parallel capacities add. The change to their total has the same factorial factor by the Dirichlet simplex integral after the Pixton edge factorials and normal signs are included. Bridges factor by conservation. Any self-loop has no highest charge contribution.

A simple graph with two independent cycles needs at least five edges. Consequently every scalar graph with at most four edges is covered after merging parallel bundles, including all branch and leg decorations in degrees up to four. This proves equality of the two graph prescriptions through that degree for arbitrary conserved scalar charges. Direct discrete residues, with an independent interpolation first in the modulus and then in the marked charges, give $217/120$ on the triangle with charges $(1,2,-3)$, $3653/13440$ when its first edge is doubled, and $-9529/960$ on the four-cycle with charges $(3,-1,2,-4)$. [Discrete checks including a five-edge graph](code/pixton_flow_cycles.py)

#### Deletion and contraction for every scalar double-ramification graph

On a connected graph orient the edges and prescribe real vertex charges of total zero. Put $m_e=k_e+1$ for the edge with branch cotangent power $k_e$ and $d=\sum_{e\,;\,\mathrm{edges}}m_e$. Compare the full-support Newton moment, with the divisors $\prod_e k_e!$ and the branch signs already included as in the scalar evaluator, with the charge-degree-$2d$ part of Pixton's regularized edge coefficient. At modulus $r$ the corresponding edge function on the circle is $q_m=(-1)^{m-1}f_r^m/(2^m m!)$, with $f_r(u)=u(r-u)$ on $0\le u<r$, extended periodically.

For integral charges, Pixton's fixed-graph Proposition $3''$ states that the unnormalized discrete sum of any half-edge polynomial is eventually a polynomial in $r$ divisible by $r^{h^1(\Gamma)}$. This applies to each edge weight: the product of the two opposite half-edge residues equals $f_r$, including at zero. Zero-charge additional legs permit arbitrary vertex valences in the graph identity. The normalized sum has degree at most $2d$. To extract its top charge part, dilate the charges by an integer $t$ and evaluate the polynomial at $r=tR$ at $2d+1$ distinct fixed integers $R$ larger than the total positive charge. In this range the affine flow chambers are fixed as $r$ varies: choose each edge residue between $0$ and $r$, and lift the divergence equalities by integer multiples of $r$. The possible lifts form a finite set. The incidence matrix is totally unimodular. At each vertex of a flow polytope the free forest currents are cut sums of the charges plus integer multiples of $r$, with every other current at $0$ or $r$. Thus the signs defining its tangent cones are fixed throughout this range. Integer vertex positions and fixed tangent cones give the same weighted lattice polynomial there by the vertex generating function, including for zero cut sums. Positive edge powers vanish at both ends, so coincident closed residue-box boundaries add zero weight. [Fixed-graph polynomiality, §1.4 and Appendix A.1–A.3](https://arxiv.org/abs/1602.04705)

After division by $t^{2d}$ the values converge, by Riemann summation in independent cycle coordinates, to the normalized continuous circle-flow integrals at $R$. Polynomial interpolation at the fixed nodes is linear. Thus the continuous integrals in the large-modulus chamber are themselves values of a polynomial of degree at most $2d$ in $R$, and its constant equals the coefficient of $t^{2d}$ in the discrete constant term. Indeed interpolation for the dilated polynomial in $R$ uses the same evaluations at zero for every $t$. Polynomial dependence on positive integral $t$ also follows from weighted Ehrhart dilation at each node. This proves regularity at modulus zero separately for every graph and every set of nonzero edge powers. Constant factors on a nonbridge edge correspond to its contraction. The continuous constant depends continuously on the charges, since the same interpolation nodes work throughout a neighborhood and the periodic integrands are continuous. The construction on rational rays extends to real charges by homogeneity and continuity.

Shift the charges at the endpoints of a nonbridge edge $e$ oppositely by $s$. Translate the particular flow along that edge alone. For $m_e=1$, the second derivative on the periodic circle is $q_1''=-1+r\delta_0$. The constant term $-1$ removes the edge weight, giving the flow on the contracted graph. The normalized delta fixes the edge current to zero, giving the flow on the deleted graph: projection to that edge is surjective with uniform circle measure. Thus the second charge derivative of the constant term is the deleted-graph coefficient minus the contracted-graph coefficient. At $m\ge2$ there is no point mass; differentiating on the open edge interval gives
$$
q_m''=(2m-1)q_{m-1}+\frac{r^2}{4}q_{m-2},\qquad q_0=-1.
$$
The second term has zero constant by the fixed-graph regularity above. These distributional differentiations commute with extraction of the constant: use common large-modulus interpolation nodes on a neighborhood for all the polynomials involved.

In the Newton integral fix the other capacities. The allowed transfers between the endpoints of $e$ using only those edges form an interval $[L,R]$, possibly empty, by convexity of the scalar flow box. Freeing an edge in the other Möbius subtractions can give an unbounded interval. A finite edge of capacity $\rho$ allows transfers within distance $\rho$ of the interval. Integrating the finite-edge indicator minus the freed-edge indicator against $\rho^{2m-1}d\rho/[2^{m-1}(m-1)!]$ therefore gives zero for an empty interval and otherwise
$$
-\frac{\operatorname{dist}(s,[L,R])^{2m}}{2^m m!}.
$$
At $m=1$ its second derivative is the feasible-transfer indicator for deletion minus the indicator that any transfer is possible, which is the contracted-graph test. At $m\ge2$ its second derivative is $2m-1$ times the integrated edge of power $m-1$. Integrating the remaining full-support differences gives the same recurrences as on the circle. The remaining capacities may be truncated to a fixed box uniformly on every compact set of charges, because cycle removal bounds every necessary current by the total positive charge. Thus the recurrences also hold across charge walls as distributional identities.

A bridge carrying fixed current $s$ factors in both prescriptions with coefficient $-s^{2m}/(2^m m!)$. A self-loop has zero full-support moment, and its circle integral is a constant times $r^{2m}$ multiplying a regular polynomial. The edgeless connected graph has coefficient one. Induct in $d$. For a graph with no bridges or self-loops, all second derivatives of the difference along the edge directions of any spanning tree vanish by the recurrences. These edge directions form a basis of the conserved vertex-charge space. The difference is affine separately in each basis coordinate and has degree at most $|V|-1$. Both coefficients are homogeneous of degree $2d>|V|-1$, so the difference is zero. Bridge factorization covers the remaining graphs. This identifies the scalar Newton moment with the top homogeneous Pixton coefficient for arbitrary graphs and all branch powers. Multiplication by the common compact marked cotangent factors and the graph symmetry factors gives, for a single target at every positive genus,
$$
\int_{\mathcal M_{g,n}}\frac{\Theta_g(p)^g}{g!}\prod_i\psi_i^{k_i}
=[t^{2g}]\int_{\overline{\mathcal M}_{g,n}}\mathrm{DR}_g(tp)\prod_i\psi_i^{k_i},
\qquad \sum_i k_i=2g-3+n .
$$
The identification of the left side uses the height intersection formula above. On the first simple two-cycle graph, the five-edge diamond with path lengths $(1,2,2)$, direct summation over both free residues and an independent capacity-polytope integration give $-113/604800$ at charges $(1,-1,-1,1)$ and $-8599331/75600$ at $(1,2,3,-6)$, in the order source, two intermediate vertices, sink. [Diamond replay including higher edge powers](code/diamond_flow.py)

In the polynomial solver choose the current on each edge of a spanning tree as a coordinate of the conserved vertex charges. To recover a monomial of total degree $2d$, choose its first tree coordinate $j$ with exponent $\alpha_j\ge2$. The coefficient equals the known coefficient at exponent $\alpha-2e_j$ in the edge's second derivative divided by $\alpha_j(\alpha_j-1)$. Contracted charges add at the merged vertex. For a bridge the deleted graph contributes zero and the same recurrence applies; a self-loop contributes zero. The difference equation therefore integrates by division of individual rational coefficients, with no undetermined constants. This replaces high-dimensional capacity triangulation in the scalar compact tests at degree five and above. For complete graphs on four and five vertices with charges $(1,2,3,-6)$ and $(1,2,3,4,-10)$, all undecorated edges give respectively
$$
-\frac{5277073217}{79833600},\qquad
\frac{14913154107446216480777}{608225502044160000}.
$$
The first also agrees with six-dimensional capacity triangulation. In genera five and six, on charge lists $(1,-1)$ and $(2,-3,1)$, the compact graph sum reproduces all 179 complementary cotangent distributions of the independent marked-point formula. [Polynomial recurrence and checks](code/scalar_flow_recursion.py)

#### Mixed Newton moments on one cycle

On a cycle choose oriented vector-valued currents $x-P_i$ on its edges. The differences of successive offsets $P_i$ are fixed by charge conservation. Write $A_i(x)=|x-P_i|^2/2$ in the supplied target metric. The Newton body of the minimum quadratic energy is the union, over $x$, of the orthants $a_i\ge A_i(x)$. It is convex by convexity of the squared norm. This calculates the cycle symbol associated with the height.

Use nonnegative integer powers $k_i$ and select edge 1. The full edge-support moment is the constant term in $R$ of the integral of $\prod_i a_i^{k_i}$ over the body inside $[0,R]^E$, where $E$ is the number of edges. Indeed, when any constraints are dropped their free coordinates integrate to positive powers of $R$. The alternating sum over the projections is compactly supported: a nearest feasible $x$ may be projected into the convex hull of all offsets, decreasing every distance. The full moment is therefore independent of the cutoff for sufficiently large $R$ after the subtractions.

Apply the divergence theorem with vector component $a_1^{k_1+1}\prod_{i\ne1}a_i^{k_i}/(k_1+1)$ in direction 1. The upper face has zero constant term. On the lower boundary let $J$ be the constraints active with positive multiplier. The unique minimizing $x$ is their convex combination. Conversely a positive convex combination minimizes the corresponding strictly convex weighted energy, hence gives a unique boundary parameter. An affinely dependent active set has dimension less than $E-1$ and has zero surface measure. Only $J$ containing edge 1 contribute. Write $x=P_1+\sum_{j\in J\setminus\{1\}} t_j(P_j-P_1)$ on the unit simplex $t_j\ge0,\ \sum t_j\le1$. The projected area Jacobian is $\det G_J(1-\sum t_j)$, where $G_J$ is the Gram matrix of these differences: the derivative of $(A_j)_{j\in J\setminus\{1\}}$ is $(\mathbf 1t^{\mathsf T}-I)G_J$. Integrating every inactive $a_i$ from $A_i(x)$ to $R$ now gives the full moment
$$
-\sum_{\substack{J\ni1\\\mathrm{affinely\ independent}}}\det G_J
\int_{\mathrm{simplex}}(1-\textstyle\sum t_j)
\frac{A_1(x)^{k_1+1}}{k_1+1}
\prod_{j\in J\setminus\{1\}} A_j(x)^{k_j}
\prod_{i\notin J}\frac{-A_i(x)^{k_i+1}}{k_i+1}\,dt .
$$
All integrands are polynomials. For three edges let $a,b,c$ be the squared distances of the three pairs of offsets, equivalently the squared vertex charges. At $k_i=0$ the simplex integrations give
$$
-\frac{(a+b+c)^3}{960}+\frac{(a+b+c)(ab+bc+ca)}{120}-\frac{7abc}{240}.
$$
The value for three unit squared charges is $17/960$. With charges $+1,-1$ at the terminals of two charge-free paths of lengths $r,s$, the squared half-capacity of each path is the minimum of the independent exponential edge energies. Its rate is the sum of the edge rates. Expanding the two-edge failed-cut integral then independently gives $(-1)^{r+s-1}rs(2r-1)!(2s-1)!/[2^{r+s-2}(2r+2s)!]$ for the full square-free coefficient. [Mixed simplex moments and independent cuts](code/cycle_height_moments.py)

Parallel lines between the same vertices are equivalent to one capacity $\rho=\sum_i\rho_i$, where $\rho_i=\sqrt{2a_i}$: a Euclidean flow of any norm up to this sum can be split proportionally among the lines. Freeing any one line removes the whole constraint. For $m$ parallel lines the fiber integral of $\prod_i\rho_i^{2k_i+1}d\rho_i/2^{k_i}$ therefore replaces them by one line with power $K=\sum_i k_i+m-1$ in the moment, with multiplier
$$
\frac{2^{m-1}\prod_i(2k_i+1)!}{(2K+1)!}.
$$
A bridge carrying fixed current $p$ contributes $-(|p|^2/2)^{k+1}/(k+1)$ to the signed moment. Self-loops are absent from every cut. After parallel lines are combined, a connected graph of at most four edges is a tree or has a single cycle. The projection formula and the cycle simplex integral consequently evaluate all mixed-target graph symbols through four Betti factors. The Taylor divisors $\prod_i k_i!$ refer to the original edges. For a triangle with one doubled edge and three unit squared charges the square-free four-edge coefficient is $5/8064$. A unit source and sink at the ends of the doubled edge instead give $1/6720$, also obtained directly as $\int_0^1 r^3(1-r)^4dr/24$.

On the mixed four-loop D3 integral there are 19 admissible choices of three retained factors and 31 choices of four. Expanding the other factors as their cotangent averages, the best values at each degree retain only $\Theta_{35}$ or only $\Theta_{45}$ (one-based marked-point labels). They give
$$
C_4(e^1+e^2)\le\frac{66530581}{10146860236800},
\qquad
C_4(e^1+e^2)\le\frac{1584203}{684913065984}
$$
respectively, in the height intersection formula. Independently the full elliptic action in target dimension four reproduces all 110 cotangent distributions on the five-, six- and seven-mark fixtures. [Polarization, exhaustive bounds and elliptic comparison](code/betti_mixed_cubic.py) · [Stable-graph evaluation](code/betti_power_graphs.py)

Retaining the separating trace of one more factor improves the upper bound. In the toroidal height intersection, the difference between the pullback of $D_{ij}$ and the true Betti divisor is effective: bridge currents are fixed by conservation, and the rest of the minimizing energy is nonnegative. Pair this difference with four true nef Betti factors and nine nef cotangent averages. On a separating node, the Jacobian and the energy split over the two components. The opposite currents at the new marks balance their charges. For component genera $h,4-h$, only the powers $h,4-h$ of the marked-pair Betti forms can contribute by the rank of the Jacobians. The restricted divided fourth power is the product of the two divided powers. With both charged marks on the same side the other component must have genus zero. The fifth pair orients the separating node, with no residual edge exchange.

For retained $\Theta_{35}^4$, the fifth choices $D_{12},D_{25},D_{45}$ give, respectively,
$$
\frac{18671669}{8561413324800},\qquad
\frac{1172347}{570760888320},\qquad
\frac{17739487}{8561413324800}.
$$
Testing all nine choices of the repeated pair and the fifth pair gives the smallest bound $C_4(e^1+e^2)\le1172347/570760888320$ for $\Theta_{35}^4D_{25}$ and $\Theta_{45}^4D_{25}$. The replay compares the degree-zero gluing calculation with the compact divisor and 1323 degree-two cuts with the independent two-Betti intersection. [Exact separating contributions](code/separating_betti_bound.py)

#### Mixed Newton moments on series-parallel graphs

Use conserved Euclidean vector charges $p$ and the Newton body $N_\Gamma(p)$ of energy coordinates $a_e\ge |j_e|^2/2$ for feasible currents $\partial j=p$. Its signed full-support moment of powers $k_e$ is the constant term of its integral over $[0,R]^E$ against $\prod_e a_e^{k_e}\,da_e$. Projection subtraction gives the same constant, since a freed coordinate contributes a strictly positive power of $R$. Each feasible coordinate flow can have its oriented cycles removed componentwise, decreasing all norms and bounding the remaining currents by the individual total positive charges. The projection subtraction is therefore supported in a fixed box, uniformly for charges in compact sets.

Choose a root edge $r$. Apply the divergence theorem in its coordinate, with primitive $a_r^{k_r+1}/(k_r+1)$. On a regular lower face let $J$ be the support of its inward conormal. Only $r\in J$ contributes. Contract all inactive edges and, on the quotient, solve the minimum-energy flow at positive lengths on $J$, defined up to common scale. The current on $J$ is uniquely determined by this minimization; set $A_e=|j_e|^2/2$ there. The remaining current has boundary $p-\partial j_J$ on each connected component $H$ of the inactive edges. The fiber of the face over its active energies is exactly the product of the bodies $N_H$ on these residual charges. Indeed all feasible completions with those active currents have the same minimum energy at the chosen conormal. Conversely the strictly positive active lengths fix all active currents for any completion on this face.

With coordinates $u$ on the positive projective lengths of $J$, the face contributes
$$
-\int \frac{A_r^{k_r+1}}{k_r+1}
\prod_{j\in J\setminus\{r\}}A_j^{k_j}
\left|\det\frac{\partial(A_j)_{j\in J\setminus\{r\}}}{\partial u}\right|
\prod_{H\,;\,\mathrm{inactive\ components}} M_H\,du ,
$$
where $M_H$ is the raw signed full-support moment on $H$. An isolated vertex contributes one. The Jacobian is the area projected away from the root coordinate. A smooth boundary point has a unique projective conormal. The singular boundary has zero projected area; a deficient Jacobian gives zero as well. An active self-loop carries zero current. If the active quotient has two distinct blocks, rescaling all lengths in just one block fixes every minimizing current, so its Jacobian vanishes. All positive powers of $R$ disappear after full-support subtraction. In detail, expand each inactive fiber integral as a polynomial in $R$ by projection subtraction on every coordinate subset. Its coefficients stabilize uniformly since the residual charges are uniformly bounded. Only the product of its constant terms contributes. The root upper face has no constant term for the same reason. This proves the recursion by decreasing edge count. Divide its final answer by $\prod_e k_e!$ for the normal-class coefficient used in the stable graph sum.

For a block with no $K_4$ minor take the endpoints of its root edge as terminals. A biconnected series-parallel graph has a two-terminal decomposition for such a terminal edge. In each binary series composition use $u=R_1/(R_1+R_2)$ for the effective resistances $R_1,R_2$ of the children. In a parallel composition use $u=R/R_1$ for their parallel resistance $R$. Together these fractions give a bijection of positive projective edge lengths with a unit cube: starting at unit effective resistance, recursively restore the two positive child resistances. Write the terminal voltage drop as $R(I-b)$ with input current $I$ and effective source shift $b$; a single edge has $b=0$. In series the child inputs are $I,I+P$, with $P$ the charge at the common vertex plus all interior charges of the first child. In parallel they are $b_1+u(I-b_1-b_2)$ and its complement in $I$. The parent shift is $ub_1+(1-u)(b_2-P)$ in series and $b_1+b_2$ in parallel. Induction gives all currents as polynomials in the cube variables, linear in the original charges.

The determinant in the face integral has fixed sign. In the chart $\ell_r=1$, it is the Hessian determinant of the concave minimum energy in the remaining $|J|-1$ lengths, hence has sign $(-1)^{|J|-1}$ wherever nonzero. The change to the cube is a diffeomorphism, whose orientation is calculated at its center. Contraction and edge deletion preserve the exclusion of $K_4$ minors. The recursion thus integrates polynomials on cubes with rational coefficients for every such graph, any target dimension, and rational charges. In particular all genus-two descendant vertices in the flat worldsheet formula are rational at rational charges: their dual graphs have at most two independent cycles, while $K_4$ has three. At genus three a nonrational part requires a stratum whose graph has a $K_4$ minor. [Exact cube solver, root changes and lower-degree tests](code/betti_series_parallel.py)

The five-edge diamond has a direct independent integral. Let its two terminals be joined by one direct edge and two paths of length two. Use source charge $p$, intermediate charges $q_1,q_2$, and conserve charge at the sink. On a face visible from the direct edge every path must retain an active edge; otherwise the direct current vanishes. Put $s_i$ equal to the second edge's fraction of its path length ($0$ or $1$ if only one segment is active). Let $a_0,a_1,a_2$ be the three positive fractions of the total conductance, summing to one. With $W=p+s_1q_1+s_2q_2$ the direct current is $a_0W$. Path $i$ carries $a_iW-s_iq_i$ then $a_iW+(1-s_i)q_i$. The nine active path choices therefore integrate over a triangle and up to two unit intervals. Inactive segments have fixed current. For source and intermediate charges the three orthogonal unit vectors, the square-free coefficient is $-43/16200$. The cube recursion reproduces it with all five choices of root, including choices whose inactive fiber has a cycle. It also reproduces independent scalar deletion–contraction and polarized first-variation coefficients with higher branch powers. On the theta graph with three paths of length two, source $(1,1,0)$ and intermediate charges $(0,-1,1),(2,-1,0),(0,1,1)$ give $13416305641/4191264000$ on all six roots. [Independent diamond integral](code/mixed_diamond_moment.py)

#### Three-point genus-three vertex in two target dimensions

For ordered charges $(1,0),(0,1),(-1,-1)$, equation (1.9) gives the following compact cotangent integrals. The nonnegative exponents sum to $3$. Permutations of the three marks give all ordered levels.

| Cotangent exponents $(k_1,k_2,k_3)$ | $\displaystyle\int_{\mathcal M_{3,3}}\frac{\Theta_3^6}{6!}\,\psi_1^{k_1}\psi_2^{k_2}\psi_3^{k_3}$ |
|---|---:|
| $(3,0,0)$ | $7/10368000$ |
| $(2,1,0)$ | $17/7257600$ |
| $(1,1,1)$ | $11/2177280$ |

Replacing the first two charges by the columns of a real two-by-two matrix multiplies all three integrals by its determinant to the sixth power, by the top exterior power of the two-target Betti form at genus three. Every tetrahedral distribution of three marks leaves a vertex with zero charge. Each of the $296$ monomials multiplying $\zeta(3)$ in the full six-edge interior uses at least one charge pairing incident to each vertex. Thus the tetrahedral $\zeta(3)$ coefficient vanishes. For four marks with charges $(1,0),(0,1),(1,2),(-2,-3)$ and exponents $(1,1,1,1)$ the complete vertex is
$$
\frac{545359}{33530112}+\frac{2}{1617}\zeta(3).
$$
An independent sum over labeled placements and explicit vertex genera reproduces the three rational integrals using the opposite order in the cotangent recursion. An added zero-charge mark reproduces their string and dilaton equations. [Full graph sum and independent checks](code/mixed_sextic_vertex.py) · [Tetrahedron period](code/mixed_k4_interior.py)

#### Three-point genus-four vertex in two target dimensions

For the same three charges $(1,0),(0,1),(-1,-1)$, the four ordered representatives give the complete vertices

| Cotangent exponents | $\displaystyle\int_{\mathcal M_{4,3}}(\Theta_4^8/8!)\,\psi_1^{k_1}\psi_2^{k_2}\psi_3^{k_3}$ |
|---|---:|
| $(4,0,0)$ | $731/1024192512000$ |
| $(3,1,0)$ | $1493/438939648000$ |
| $(2,2,0)$ | $157/29262643200$ |
| $(2,1,1)$ | $433/37623398400$ |

At charges $p,q,-p-q$ the multiplier is $(p_1q_2-p_2q_1)^8$. Each of the four period directions supplies a squared planar determinant. Permuting the three charges preserves this eighth power. The first three channels admit a primary absorbed root in the fivefold differential. The fourth measures the closed action with positive level at every mark.

Expand the eight Betti factors into leg and normal powers on connected stable graphs. Self-loops have zero full-support moment. At an empty trivalent vertex of star resistances $s_i>0$ replace the star by triangle resistances $R_{ij}=s_i s_j(1/s_1+1/s_2+1/s_3)$. Its inverse is $s_i=R_{ij}R_{ik}/(R_{12}+R_{13}+R_{23})$ for distinct $i,j,k$. At neighbor $i$ the original current into the star is the sum of the two outward triangle currents. For a charged bivalent vertex with load $p$ and successive resistances $s_1,s_2$, put $t=s_2/(s_1+s_2)$. Move $tp,(1-t)p$ to the two endpoints. With reduced current $I$ from the first endpoint to the second, the original consecutive currents are $I-tp,I+(1-t)p$. Parallel division adds its own independent interval. Restoring the positive resistances in reverse order proves that these fractions cover the whole projective cone. All currents are polynomial at each successful elimination. The concave minimum-energy Hessian fixes the sign of each projected Jacobian. The signed cube integral therefore reduces to monomial integrals. In particular all divided moments on a tetrahedron with one empty vertex are rational at rational charges: all its proper connected face quotients and proper inactive components are series-parallel.

By edge count $0,\ldots,8$ there are $(1,1,2,5,12,33,103,333,1183)$ loopless connected graph shapes, and $(1,24,240,1548,5296,10204,9107,3721,562)$ distinct graph, load and normal-power inputs after stability and dimension. All $84+245+89=418$ inputs requiring more than series-parallel reduction agree with a second root projection. At six edges the independent tetrahedron integral is $-4657/217728000$ on the three unit-triangle loads and an empty fourth vertex. A global leg and half-edge expansion with first-KdV intersections independently reproduces all $61731$ channel weights. Grounded Laplacians on $244$ rational networks reproduce every transformed edge current. Run the linked planar evaluator with `--genus 4 --check-roots` and its assembly check without arguments. [Complete graph sum](code/planar_four_loop.py) · [Independent weights and current checks](code/planar_four_loop_assembly.py) · [Star charts](code/betti_star_charts.py) · [Independent tetrahedron integral](code/mixed_sextic_vertex.py)

At genus five, delete edges $01,02$ from the complete graph on $0,\ldots,4$ and place the three charges at $0,1,2$. The only initial move eliminates the charged bivalent vertex $0$. The resulting parallel reduction leaves a tetrahedron with all four vertices charged. This eight-edge full face occurs at levels $(3,1,1)$ with genus one at $0$ and one cotangent at each end of $34$. The divided branch coefficient has stable multiplier $1/24$. At lexicographic conductances $1,\ldots,8$ every root energy projection has nonzero Jacobian. A cut integrating both current fractions and conductance $c_{12}$ leaves four positive metric variables. With $c_{34}=1$ and $(c_{13},c_{14},c_{23},c_{24})=(1,2,3,1)$ the divided lower-face density integrates on this cut to
$$
-\frac{465489528164329449978660994907}{22686383745706107326428233439728280322801141376}.
$$
Integrating also $c_{13}$ leaves three positive metric variables. At $(c_{14},c_{23},c_{24},c_{34})=(2,3,1,1)$, the positive projected density on this second cut is
$$
-\frac{71161938221472988962540094215797989902554359889}{13598660191230190000557188561633690815707322368000000000000000}
+\frac{23062673138156617081}{5980085648358972073353442597384800}\log(225/56).
$$
Here $15c_{13}+14$ and $4c_{13}+15$ are the deletion and contraction tree factors. On the full positive projective cone the five nonroot energy gradients have six complementary minors. After multiplication by $U^3/\ell_{34}$ they are polynomials, where $U$ is the tetrahedron Kirchhoff polynomial and $\ell_{34}$ is the root resistance. The numerator has degree at most twelve in the splitting fraction. The resistance density has denominator $U^{16}$. Its $40489$ convergent monomials give $5648$ tetrahedral symmetry orbits. Each orbit agrees with its reciprocal-length dual. The positive projected area is
$$
\frac{99133977499}{26125187194397491200000}-\frac{\zeta(3)}{3692000632320}.
$$
The proper Newton faces sum to $17914642896173/15342609934164344832000$. Subtracting the projected area gives the complete divided coefficient of this decorated graph,
$$
\frac{3273677772720991}{2812811821263463219200000}+\frac{\zeta(3)}{3692000632320}.
$$
Projection from edge $12$ gives different proper and full rational parts and reproduces the same divided coefficient after addition. This single term carries the stable multiplier $1/24$.

For the specified marking and levels $(3,1,1)$, stability and degree force the genus-one component to lie at $0$. Vertices $3,4$ each carry one cotangent insertion, with four half-edges available at each vertex. The sum over all $16$ placements, including the genus-one intersection $1/24$ and the order-two graph symmetry, gives the complete contribution of this marked graph,
$$
\frac{3963229030865953661}{7088285789583927312384000000}-\frac{343723}{34557125918515200}\zeta(3).
$$
For levels $(4,1,0)$ the same graph has marks $1,2,0$ in charge order. Genus one is forced at vertex $1$, with the same $16$ flag placements at $3,4$. Root projections $34$ and $12$ both give the stable contribution
$$
\frac{145775286798138013}{26059874226411497472000000}-\frac{9817}{863928147962880}\zeta(3).
$$
Across all graphs with at most eight edges, these two stable contributions give the entire $\zeta(3)$ coefficient in their respective channels. [Full-face rank and stability](code/genus_five_star_obstruction.py) · [Two conductance cuts](code/genus_five_k4_cut.py) · [Full graph, dual periods and eight-edge Jacobian checks](code/genus_five_k4_period.py)

**Genus-five three-mark census.** In the order $(5,0,0),(4,1,0),(3,2,0),(3,1,1),(2,2,1)$, the stable marked orbits remaining after the polynomial-chart and rank tests are counted below

| Edges | Non-series-parallel shapes | Stable three-occupied-vertex orbits | Unreduced orbits in the five channels |
|---:|---:|---:|---|
| 6 | 1 | 1 | $(0,0,0,0,0)$ |
| 7 | 3 | 33 | $(0,0,0,0,0)$ |
| 8 | 20 | 304 | $(0,1,0,1,0)$ |
| 9 | 118 | 743 | $(0,3,0,3,0)$ |
| 10 | 714 | 714 | $(0,6,0,6,0)$ |

The two eight-edge orbits are evaluated above. At most five edges every shape is series-parallel. Coincident marks have scalar conserved loads and rational flow-polytope moments. For distinct marks, contract each inactive connected piece. If the quotient has at most two nonzero loads, all active currents are collinear. With $E$ active edges and $V>2$ quotient vertices its energy differential then has rank at most $E-V+1<E-1$. For two independent target currents the rank bound is $2(E-V+1)$, since each current varies only along cycle flows. In particular a ten-edge full face on seven vertices has rank at most eight and its nine-dimensional projected Jacobian vanishes. For a surviving face, a root with series-parallel vertex deletions at both endpoints ensures that each inactive piece is series-parallel because the endpoints lie in distinct pieces. For the sole remaining ten-edge graph the recursive check retains the possibly nonzero loads at each active-edge endpoint. All ten root projections admit polynomial charts on every contributing piece.

Thus the complete $(5,0,0),(3,2,0),(2,2,1)$ genus-five vertices are rational. Every remaining nine- or ten-edge obstruction in $(4,1,0)$ or $(3,1,1)$ has exactly nine active edges. The twelve ten-edge orbits project onto the same six marked nine-edge quotients as the full nine-edge faces. These give two unmarked shapes. Positive series, parallel and empty-star moves take both to a charged tetrahedron. All lifted currents are polynomial in the independent splitting fractions and in tetrahedron tree-current numerators divided by $U$. After integrating the fractions the projective denominator is a single power of $U$. On each Schwinger sector absolute integrability forces strict convergence of each numerator monomial. Scaling any proper subset of lengths to zero leaves a lowest-degree polynomial nonzero on an open set of positive ratios. The tetrahedron cut therefore expresses every period in $\mathbb Q+\mathbb Q\zeta(3)$.

Both remaining $\zeta(3)$ columns cancel. Write $J_0,\ldots,J_3$ for the four chord currents of the lexicographically first spanning tree on each quotient, ordered by edge. Combining nine-edge full faces and ten-edge proper faces gives the three insertions
$$
-\frac{|J_2+J_3|^2}{96},\qquad
-\frac{|J_2+J_3|^2+2|J_0+J_2|^2}{96},\qquad
-\frac{|2J_0+J_2+J_3|^2+2|J_2|^2+2|J_3|^2}{96}.
$$
Each multiplies the positive projected area and the root energy. Grounded currents at unequal positive resistances give a strictly negative value in each case. The first term in each period sum below is the complete eight-edge marked graph. The other three terms are the periods of the displayed insertions.
$$
\begin{aligned}
(4,1,0):&\quad -\frac{9817}{863928147962880}+\frac{25457}{17278562959257600}
-\frac{341063}{86392814796288000}+\frac{597739}{43196407398144000}=0,\\
(3,1,1):&\quad -\frac{343723}{34557125918515200}+\frac{1877}{533288980224000}
+\frac{93371}{345571259185152000}+\frac{2127563}{345571259185152000}=0.
\end{aligned}
$$
On the empty-star graph with triangle resistances $X,Y,Z$ and $S=X+Y+Z$, the original three star resistances are $XY/S,XZ/S,YZ/S$. The star polynomial is $U_{\Delta}/S$ and the absolute change of measure is $XYZ/S^3$. Exact division cancels all factors of $S$ before the parallel integrations. The resulting tetrahedra have denominator $U^{14}$. Every convergent monomial agrees with its reciprocal cut. Thus **all five complete genus-five three-mark vertices are rational**. The eight-edge marked graph and the three grouped nine-active-edge terms contribute rational subtotals $352407233223796349/63175452670088478720000000$ and $1172053333004275717/2147965390783008276480000000$ in channels $(4,1,0)$ and $(3,1,1)$. [Summed graph insertions](code/genus_five_nine_faces.py) · [Successive series and parallel integrations](code/genus_five_nine_series.py) · [Empty-star period and complete cancellation](code/genus_five_nine_star.py)

Adding the rational charts at every edge count to these exceptional terms gives the complete unit-triangle values
$$
\int_{\mathcal M_{5,3}}\frac{\Theta_5^{10}}{10!}\,\psi_1^4\psi_2
=\frac{14951}{6083703521280000},\qquad
\int_{\mathcal M_{5,3}}\frac{\Theta_5^{10}}{10!}\,\psi_1^3\psi_2\psi_3
=\frac{5783}{506975293440000}.
$$
For these two channels the complete polynomial chart sums at nine edges are $-11951111019345795501914309/567062863166714184990720000000$ and $-482890555953316000326163/23627619298613091041280000000$. At ten edges they are $4528990077096507699697/272626376522458742784000000$ and $895146285295273962968201/70882857895839273123840000000$. The independent half-edge expansion checks all stable weights. Replay with `genus_five_rational.py --channels 1 3 --check-weights`.

The other two unit-triangle values are
$$
\int_{\mathcal M_{5,3}}\frac{\Theta_5^{10}}{10!}\,\psi_1^5
=\frac{179}{442451165184000},\qquad
\int_{\mathcal M_{5,3}}\frac{\Theta_5^{10}}{10!}\,\psi_1^3\psi_2^2
=\frac{463}{86910050304000}.
$$
Exchange the first and third marks while retaining the ordered charges $(1,0),(0,1),(-1,-1)$. Their squared lengths are $1,1,2$, so this changes the individual leg weights and edge-count sums. The complete integral is invariant because the induced target map has determinant $-1$ and the top Betti product has weight ten. For levels $(0,0,5),(0,2,3)$ the nine-edge contributions are $100248615079129507927/7364452768398885519360000$ and $668074668662909828471/245481758946629517312000000$. At ten edges they are $-677933885186781521/323002314403459891200000$ and $112278650612830526521/30685219868328689664000000$. With the original marks the nine-edge terms instead equal $-5195601237791028319477/1472890553679777103872000000$ and $-1280692740894716008217/49096351789325903462400000$. The two full sums agree exactly after all eleven edge counts. Replay both evaluations with `genus_five_rational.py --channels 0 2 --check-weights`, adding `--swap-marks` for the second marking.

For the ordered unit-triangle charges the complete $(2,2,1)$ vertex is
$$
\int_{\mathcal M_{5,3}}\frac{\Theta_5^{10}}{10!}\,\psi_1^2\psi_2^2\psi_3
=\frac{51787}{2897001676800000}.
$$
The $0$–$8$ edge subtotal is $24942350470871679106273/3221948086174512414720000000$. The nine- and ten-edge contributions are $-153185788543248886281611/6443896172349024829440000000$ and $2295582284289643992637/143197692718867218432000000$. At charges $p,q,-p-q$ multiply the complete result by $(p_1q_2-p_2q_1)^{10}$. For each of the five period directions the top exterior power of the two-target Betti form supplies the squared determinant.

Parallel capacities add by convolution. An uncharged bivalent vertex takes their minimum, giving additive exponential rates. Articulation blocks factor with their attached charges added at the cut vertex. For two target currents, Cauchy–Binet factors the energy Jacobian through the chord currents whenever their dimension equals the face dimension or exceeds it by one. FLINT integrates each cube antiderivative at its upper endpoint before advancing to the next coordinate. On the ten-edge prism with the additional charged diagonal the complete moment $-7820529758891/53645489280294912000000$ took 11 seconds in the factored evaluator and 53 seconds in the direct replay. The capacity reductions agree with all original graph moments through five edges and the independent full six-edge sum. A separate global half-edge expansion with the first KdV recursion reproduces all $302666$ stable weights. [Complete rational integral](code/genus_five_rational.py) · [Independent high-edge moments](code/genus_five_rational_check.py) · [Stable-orbit and face replay](code/genus_five_graph_census.py) · [Two-root period and charge-permutation replay](code/genus_five_k4_period.py)

#### One further independent Betti form from the first variation

Give the second target current with charges $q$ energy weight $t\ge0$, keeping weight one on the scalar charges $p$. The Newton body is the Minkowski sum $N(p)+tN(q)$. On an open facet of $N(p)$ for a nonzero bond charge $p(C)$, an inward conormal in the energy variables is $\ell_e=1/\rho_e$ on $C$ and zero elsewhere. Contract the two connected sides. Minimizing the second current energy at this conormal gives
$$
h_q(\ell)=\frac{q(C)^2}{2\sum_{e\in C}\rho_e}
=\frac{q(C)^2}{2|p(C)|}.
$$
Consequently the weighted first variation is obtained by differentiating the scalar cut bound with this velocity. Indeed differentiating a moving level surface gives the facet density times conormal displacement divided by the conormal length. The facet intersections have zero surface measure. This first-variation formula follows also by taking the weighted boundary-layer integral on compact subsets of the regular facets, then exhausting the boundary of the convex body; bounded minimal flow energies give uniformly bounded displacement. On coordinate facets only a bridge can carry nonzero second energy by itself. A bridge of zero primary charge translates in its energy coordinate by $tq(C)^2/2$. This contributes $-q(C)^2/2$ to its linear term when its branch power is zero. For multiple nonzero primary cuts differentiate the inclusion–exclusion integrals on the corresponding capacity polytopes. Dropping any participating $\rho_e$ with unit coefficient gives the projected facet measure. The leg exponential supplies in addition $\sum_iq_i^2\psi_i/2$ times the previous pure graph sum.

Thus $\Theta(p)^{d-1}\Theta(q)/(d-1)!$ has a rational stable-graph algorithm in every degree $d$. For the mixed D3 integral the best of the five inequivalent choices with four repeated forms and one independent form keeps $\Theta_{12}^4\Theta_{25}$. Replacing the remaining nine forms by cotangent averages gives
$$
0<C_4(e^1+e^2)\le\frac{171793}{152882380800}.
$$
The three fifth-pair choices $12,25,45$ after four retained $35$ give $3226883/1712282664960$, $12053471/8561413324800$ and $2102297/1426902220800$. The other inequivalent choice $\Theta_{12}^4\Theta_{35}$ gives $3629/1651507200$. The independent mixed-vector solver reproduces the local derivatives on all 34 four-edge-or-decorated controls and both full mixed quartic stencils tested. [First variation and exact five-factor graph sums](code/betti_first_variation.py)

Retain one further separating trace $D_{ij}$ with four copies of $\Theta_{12}$ and one true $\Theta_{25}$. On a separating node of genera $h,4-h$ the four copies split as $h,4-h$ by rank. The additional $\Theta_{25}$ restricts once to each side in the sum, and eight remaining cotangent averages pull back to the corresponding components. A component with zero $\Theta_{25}$ charge contributes zero when assigned that form. The effective height difference between $D_{ij}$ and the sixth actual Betti class pairs nonnegatively with the preceding five nef classes and the eight averages. For $D_{25}$ the resulting upper bound is $231473/214035333120$. Either $D_{35}$ or $D_{45}$ gives
$$
0<C_4(e^1+e^2)\le\frac{113203}{107017666560}.
$$
This is the minimum among the fifteen tests across all five inequivalent choices of the true factors. [Component polarization, five-factor limit and all separating cuts](code/separating_mixed_betti_bound.py)

#### Genus two in one target dimension and the double ramification cycle

For integral charges $p_j$ summing to zero, let $\mathrm{DR}_2(tp)$ be the genus-two double ramification cycle with all charges multiplied by $t$. For nonnegative descendant levels $k_j$ summing to $n+1$ on $n$ labeled marks, the one-dimensional BCOV vertex satisfies
$$
V_{2;\boldsymbol k}(p)=[t^4]\int_{\overline{\mathcal M}_{2,n}}\mathrm{DR}_2(tp)\prod_{j\,;\,\mathrm{marks}}\psi_j^{k_j}.
$$
Polynomiality extends this equality to arbitrary charges. Conservation turns the weighted sum of the preceding divisors into $\tfrac12(\sum_i p_i^2\psi_i-\sum_I p_I^2\delta_I)$, where each stable separating divisor occurs once and $p_I$ is the total charge on one side. Its half square gives the terms on compact-type graphs. In Pixton's genus-two formula the prefactor is $1/4$. With no edges the leg exponentials give $(\sum_i p_i^2\psi_i)^2/8$. At one separating edge the constant term in the regularization modulus $r$ of the edge weight is $-p_I^2$. The linear leg terms give the cross term; the next term of the edge factor gives $-p_I^4(\psi_h+\psi_{h'})/8$ at its two branches. On two distinct separating edges the weights multiply. These are precisely the intersections of the compact divisors, including the normal bundle of a repeated divisor.

A self-loop has a free residue and its averaged edge factor is independent of the marked charges. With one other edge or one extra degree of decoration the remaining charge degree is at most two. The remaining connected two-edge graph has two parallel edges. Choose the sides so that their total charge $s$ lies between zero and $r$, taking $r$ large. On the cyclic group put $f(u)=u(r-u)$ on residues $0,\ldots,r-1$. Splitting at $u=s$ and summing powers gives
$$
\frac1r\sum_{u=0,\ldots,r-1\,;\,\mathrm{residues}} f(u)f((s-u)\bmod r)
=\frac{r^4-1}{30}-\frac{s(r-s)[s(r-s)+2]}6 .
$$
The charge-quartic part at $r=0$ is $-s^4/6$. Multiplication by $1/4$ gives $-s^4/24$ before the edge-exchange symmetry. This equals the correction to the half square from the two-Betti integral. The graph replay compares 131 distributions of marked charges and cotangent powers, including zero charges. [Pixton formula, §0.4 and Theorem 1](https://arxiv.org/abs/1602.04705) · [Exact graph replay](code/soft_quartic_two_betti.py)

With all the cotangent power at one mark $i$, Theorem 1 of Buryak–Shadrin–Spitz–Zvonkine gives the following result after extraction of the charge quartic:
$$
V_{2;(n+1)e_i}(p)
=\sum_{j\ne i\,;\,\mathrm{marks}}\frac{p_j^4}{1920}
+\sum_{j<l,\ j,l\ne i\,;\,\mathrm{marks}}\frac{p_j^2p_l^2}{576}.
$$
Indeed $\sinh(z/2)/(z/2)=1+z^2/24+z^4/1920+\cdots$, and the charge-independent denominator contributes only below charge degree four. Their Theorem 2 computes the corresponding coefficient for any distribution of descendants. For charges $(2,-3,1)$, levels $(4,0,0),(3,1,0),(2,1,1),(0,2,2)$ give respectively $7/120,1/6,83/240,133/576$. Both the regularized graph and the independent marked-point formulas include the labeled-mark multiplicities. [Published descendant formulas](https://arxiv.org/abs/1211.5273) · [Independent polynomial replay](code/one_target_dr.py)

For any $g\ge2$, the same two-Betti calculation with one marked cotangent power gives
$$
\frac12\int\Theta_g(p)^2\psi_i^{3g-5+n}
=\frac{1}{24^{g-2}(g-2)!}
\left(\sum_{j\ne i\,;\,\mathrm{marks}}\frac{p_j^4}{1920}
+\sum_{j<l,\ j,l\ne i\,;\,\mathrm{marks}}\frac{p_j^2p_l^2}{576}\right).
$$
To verify it, polarize the two Betti factors with independent charge lists. Each pair of marked differences uses at most four labels beyond $i$; the string equation forgets all other labels. Divide by $1/(24^g g!)$. In the compact divisor product the sum of the powers away from $i$ is at most two. String, dilaton and Virasoro reduce it to a polynomial of degree at most two in $g$, with the only new two-point ratio $(12g^2-12g+5)/5$ for levels $(2,3g-3)$. On a separating tree each component away from the distinguished one has at least two flags and at most one cotangent power, so its contribution requires genus zero. In the parallel-edge term the distinguished genus drops by one, giving ratio $24g$. Symbolic independent charges at $g=2,3,4$ therefore determine every coefficient. At $g>2$ a full one-dimensional BCOV vertex uses $g$ Betti factors. [Polarized polynomial and five-genus replay](code/one_target_dr.py)

### The highest-pole coefficient under changes of field coordinates

At genus $n=2,3$ take the normal residue after multiplication by $w_3^{2n}$ and evaluate the tangent vector at $z=0$. In the exponential probes above replace each pure normal wave $e^{w_j}$ by $e^{w_j}-1-w_j$ and each mixed wave $L(w)e^{H(z)}$ by $L(w)(e^{H(z)}-1)$, where $L,H$ are the linear functions in the probe. The primary trace vanishes quadratically in each regular scalar momentum with its balancing recoil. These replacements therefore preserve its residue. Every new scalar has zero value and zero first derivatives at the origin. The genus-two polynomial probe already has this property. Its projected coefficient along $\partial_{z_2}$ is $-N/1620$.

Fix all lower-genus coordinates. An order-$\lambda^{2n}$ holomorphic local change of the flat fields adds its commutator with the classical differential. Consider constant-coefficient jets equivariant under the fivefold volume group, or just under the block group $S(GL_2\times GL_3)$ preserving the tangent plane, normal plane and their product volume. The latter group preserves the tangent bivector times the normal residue. The part of the coordinate change that contributes through a classical interaction has $10n$ holomorphic derivatives because the classical vertex supplies the remaining divergence. On the binary term its $2n+1$ primary inputs preserve total polyvector rank. A descendant output in a higher classical vertex leaves additional scalar spectators by the classical descendant-level rule, with the same rank balance and derivative count.

For the block group, in an invariant contraction subtract the number of normal input vector indices from the total number of normal covector indices, including covectors testing output vectors. The analogous tangent difference and the normal difference sum to $10n$ by rank balance. Invariance under the relative dilation $(z,w)\mapsto(c^3z,c^{-2}w)$ sets three times the tangent difference equal to twice the normal difference. The differences are $4n,6n$. Pairwise contractions and the two volume tensors generate the invariant contractions. In all terms with the flux inside the coordinate change there are no normal input vector indices, so the normal covectors form exactly $2n$ alternating triples. Every triple contains the same normal slot momentum at most once.

If the classical interaction is outside the coordinate change and contains the flux directly, it differentiates the flux at most once. Otherwise its scalar spectators receive at most one derivative in total. The coordinate change supplies at most $2n$ normal flux derivatives by the alternating-triple bound. To contribute to the specified residue requires at least $2n$ derivatives in $w_3$ alone. The one outer derivative can then differentiate a spectator once, or supply one extra flux derivative whose normal pole absorbs one power from that spectator. Either way, its vanishing value and first jet kill the residue at $z=0$.

If the classical interaction is inside the change, all its inputs are primary, so only its binary term occurs. The scalar-scalar term vanishes. Fusing the flux of normal momentum $q$ with scalar momentum $p$ gives the tangent vector obtained by contracting $p$ into the flux bivector, at momentum $q+p$. A pure normal scalar gives zero. For a mixed scalar, attaining $q_3^{2n}$ uses the normal part $q$ of the fused slot in every normal triple. All further derivatives on the fused scalar are tangent. They leave the factor $L(w)$, which vanishes in the Laurent quotient with simple $w_1,w_2$ poles. Finally the Dolbeault differential commutes with holomorphic jets, and the free descendant term vanishes on primary scalars and the divergence-free flux after projection to level zero. At first order in the flux, using one of the tangent-vector source lifts at a vertex with no explicit flux gives normal form degree one. Holomorphic jets preserve this degree, and the Dolbeault differential commutes with the change. Such terms give zero in the normal form-degree-two projection.

The same bound holds for holomorphic coefficients on the punctured background preserving tangent translations, the block group and the derivative character. They depend only on $w$ and extend across $w=0$ by Hartogs. In Taylor degree $k$ their $k$ coordinate vectors $w$ raise the required derivative degree to $10n+k$. The normal covector excess after subtracting these $k$ vectors is still $6n$. A normal determinant made of three coordinate vectors vanishes because all three are $w$. Hence the normal contraction consists of $2n$ covector triples and $k$ direct pairings with $w$. It contains at most $2n+k$ flux derivatives. The degree-$k$ coefficient itself consumes $k$ powers of normal pole under the residue. In an inner fusion every triple and every one of the $k$ pairings must therefore take a flux momentum; the remaining mixed scalar still vanishes on the brane. In an outer fusion the single divergence can also differentiate the coefficient, removing one $w$ and using up its one available derivative. The spectator's zero first jet again kills the term. Thus the two measured residues survive these same-genus coordinate changes. Contributions of separated vertices to a transferred operation require a separate calculation.

For the three tangent-dependent sources above the highest-pole residue is invariant under the same class of genus-$g$ field changes at every $g\ge2$, with all lower coordinates fixed. At the origin subtract the constant and linear jets from the pure normal exponentials, using quadratic softness. From each $w_2$ mixed input subtract its first jet $w_2$. These latter subtractions also kill the maximal residue by the two-mixed-source shear test. Every input now has zero first jet. In an outer classical variation the spectator test above therefore applies. In an inner fusion all $2g$ alternating normal triples must differentiate the flux; any $w$-dependent coefficient consumes its own extra derivatives as above. The mixed scalar's positive normal factor survives into the simple $w_1,w_2$ poles and kills this term. Thus the displayed nonzero residue persists under a same-genus coordinate change at every genus. The argument fixes the lower-genus coordinates and concerns the local vertex.

### A normal Dolbeault homotopy edge at the puncture

Use the smooth unit-residue form $U=c\,\iota_{\bar w}(d^3\bar w)/s^3$ of the flux, where $s=|w|^2$, $\bar w$ is the antiholomorphic Euler vector, and $c$ is the fixed unit-residue normalization. For the $j$th coordinate antiholomorphic vector $e_j$, set
$$
b_j=\frac{c}{2s^2}\iota_{\bar w}\iota_{e_j}(d^3\bar w),
\qquad \bar\partial b_j=w_j U .
$$
Indeed $\bar\partial\iota_{\bar w}\iota_{e_j}(d^3\bar w)=2\iota_{e_j}(d^3\bar w)$ and $\bar\partial s\wedge\iota_{\bar w}\iota_{e_j}(d^3\bar w)=s\iota_{e_j}(d^3\bar w)-w_j\iota_{\bar w}(d^3\bar w)$.

Choose a smooth radial function $\chi(s)$ increasing from zero to one across a compact normal annulus. Pair normal forms by integrating their product with $d^3w$. Divide each $\bar w_i\bar\partial\chi$ by the same integral of $\bar w_1\bar\partial\chi\wedge w_1U\wedge d^3w$ and denote the resulting one-form by $l_i$. For the functional $\ell_i(\omega)=\int l_i\wedge\omega\wedge d^3w$,
$$
\ell_i(w_j U)=\delta^i{}_j .
$$
The off-diagonal integrals vanish under relative coordinate phases of determinant one. The diagonal integrals agree by signed permutations and are nonzero because $\bar\partial\chi\wedge U=c\,\chi'(s)d^3\bar w/s^2$. The map $h\omega=\sum_{i\in\{1,2,3\},\,\mathrm{normal\ coordinates}} b_i\ell_i(\omega)$ on degree two, zero in other degrees, contracts precisely these three exact pairs. In fact $p=\bar\partial h+h\bar\partial$ is a chain projection, $p^2=p$. Its kernel contains representatives of all normal Laurent classes. A contraction onto those classes can therefore split off this finite block. For the cyclic normal Dolbeault pairing on test forms one can include the dual pairs $(l_i,\bar\partial l_i)$, sending $\bar\partial l_i$ to $l_i$. Within the first block the pairing vanishes because $b_i$ and $U$ are horizontal to $\bar w$. Within the dual block a product $l_i\wedge\bar\partial l_j$ repeats $\bar\partial\chi$. The two blocks pair perfectly by the displayed identity and integration by parts.

Tensor this normal edge with the identity kernel along $\mathbb C^2_z$. The tangent factor preserves holomorphic tangent translations and $SL_2$. For its normal endpoints $x,y$, the original block has kernel $\sum_i b_i(x)l_i(y)$, with the same coordinate sum. The dual block adds its graded transpose. Take $x$ outside the cutoff annulus and $y$ inside at a point with $\chi'(|y|^2)\ne0$, so the transpose vanishes. The part of the contraction by the common relative antiholomorphic Euler vector that hits the source one-form contains
$$
-\chi'(|y|^2)\left(|y|^2-\sum_{j\in\{1,2,3\},\,\mathrm{normal\ coordinates}}y_j\bar r_j\right)
\iota_{\bar x}\iota_{\bar y}(d^3\bar x),
$$
up to the nonzero common factor $c/(2|x|^4)$ divided by the normalization integral. Here $r$ is the normal position of the root. For $x=(2,0,0)$, $y=(0,1,0)$ and $r=(0,0,1)$ the displayed expression is $2\chi'(1)d\bar x_3$; choose the cutoff annulus inside $0<s<4$ with $\chi'(1)>0$. The contraction hitting the target one-form has a different distribution of degrees at the endpoints. Thus the common radial null vector of the translation-invariant flat propagators need not persist on a genuine Dolbeault contraction of exact normal modes. The tangent identity is supported on the tangent diagonal. Extending a Dolbeault contraction to a full BCOV propagator and computing a renormalized loop integral require additional data.

On the two-loop theta with path lengths $(3,3,4)$ the saturated Dolbeault form with nine flat edges and one such puncture edge is a nonzero current. Take the three vertex paths $(0,1,2,3)$, $(0,4,5,3)$ and $(0,6,7,8,3)$. Removing edges $(4,5)$ and $(6,7)$ gives an eight-edge spanning tree. Use its five-component complex edge differences as the forty relative coordinates. Eight flat angular forms on the tree each omit their own displacement line. On those eight lines evaluate the two tangent diagonal covectors of edge $(4,5)$ and its target and source normal covectors, then the five difference covectors of $(6,7)$. Append one column consisting of four zeros followed by the displacement of $(6,7)$ to form its flat angular contraction. At the integer positions in the replay the absolute nine-row determinant is $1116$. The full forty-row exterior product independently gives $1116$ in absolute value after removing the same nonzero scalar kernel factors. The tangent delta is evaluated on its smooth support, all normal positions are nonzero, and all ordinary edge displacements are nonzero. This tests the antiholomorphic form of one edge assignment. The full propagator, polyvector contractions, graph sum and integration still enter the physical amplitude. [Integer coordinates and both exterior-product evaluations](code/puncture_theta_form.py)

### Sum of the separated one-loop pentagons

Take one flux, one test with normal trivector and antiholomorphic form $d\bar z_1\wedge d\bar z_2\wedge d\bar w_1$, and three primary scalar legs. A connected one-loop graph of cubic vertices with five distinct external positions has five vertices and one external leg at each; it is a pentagon. Use the preceding normal contraction on one line, including its cyclic transpose, and ordinary flat propagators on the other four. In $(z_1,z_2,w_1,w_2,w_3)$ order put the five positions at
$$
(0,0,1,1,1),\quad(1,0,1,1,2),\quad(0,1,1,2,1),\quad
(-1,0,1,1,2),\quad(0,0,0,1,1).
$$
Use vertex labels $0,\ldots,4$ in this order, with the flux at $1$ and the test at $2$. Only $0,4$ have coincident tangent coordinates. Their squared normal radii are $3,2$; choose the radial cutoff derivative nonzero at $2$ and supported away from $3$. Only the block with target $0$ and source $4$ contributes to the principal tangent delta derivative. All six unoriented labeled pentagons containing this edge are obtained by ordering $1,2,3$ along the path from $0$ to $4$. Each has trivial labeled automorphism group.

Remove the common flux radial factor and the scalar normalization of this puncture block. Each of the six paths has squared edge distances $(2,4,4,3)$. In fixed vertex coordinates each antiholomorphic determinant equals $-2$. The odd determinants of the two tangent delta derivatives are $0,2$ on every path. Thus the second derivative direction has combined angular coefficient $6(-2)(2)=-24$. The four Schwinger integrations also agree path by path: at squared distance $D>0$ the scalar part of a differentiated flat propagator uses $\int_0^\infty t^{-7}\exp(-D/(4t))\,dt=4^6\,5!/D^6$. Removing its universal scalar prefactor on each line gives the two coefficients
$$
0,\qquad -\frac{24}{(2\cdot4\cdot4\cdot3)^6}=-\frac{1}{32614907904}.
$$
This proves noncancellation in the separated one-loop graph sum for the specified puncture block, including all six orderings and all four proper-time integrations. The remaining bulk vertex integrations against closed sources and the other blocks of a full punctured propagator enter the physical matrix element. [Twelve pairs of absolute-coordinate exterior determinants](code/puncture_pentagon.py)

### Integrating the flux position

Keep external positions $0,2,3,4$ of the pentagon calculation and integrate vertex $1$ over the punctured fivefold with the actual flux profile $\iota_{\bar w}d^3\bar w/|w|^6$. Retain the preceding common kernel factors, extracting the flux constant $cN$ and the conversion from complex volume to real Lebesgue measure. Let $I_1,I_2$ denote the resulting coefficients of the two holomorphic tangent delta derivatives on edge $(0,4)$. They satisfy
$$
I_2=-I_1,\qquad 1<10^{12}I_1/\pi^5<11/10 .
$$
In particular the sum survives integration over the full background insertion. Its direction changes from the value at the single flux position above.

For each oriented path from $0$ to $4$, let $A,B$ be the neighbors of the flux vertex and $D$ the product of the squared lengths of the two edges disjoint from it. Write the inverse widths of the two incident flat Gaussians as $r(1-t),rt$, and that of the normal flux Gaussian as $rh/(1-h)$, with $0<t,h<1$. The angular numerator depends antiholomorphically on the flux position. Its Gaussian integral therefore evaluates it at the conjugate mean. For these real external positions the mean has tangent part $(1-t)A_z+tB_z$ and normal part $(1-h)((1-t)A_w+tB_w)$. Write $P(t,h)$ for the product of the two exterior determinants at this mean. The scale integral is $\int_0^\infty r^9e^{-rV}\,dr=9!/V^{10}$ with
$$
V=t(1-t)|A-B|^2+h|(1-t)A_w+tB_w|^2.
$$
The contribution divided by $\pi^5$ is exactly
$$
\frac{63}{5D^6}\int_0^1\!\!\int_0^1
\frac{t^5(1-t)^5 h^2 P(t,h)}{(1-h)V^{10}}\,dh\,dt .
$$
Here $63/5=9!/(2(5!)^2)$ includes the two flat denominator powers six and the flux denominator power three. At finite flat heat cutoffs the positional Gaussian integration is convergent. In the resulting parameters each determinant product is bounded by a constant times $(1-h)h(h+t(1-t))$ and $V$ is bounded below by a positive constant times $h+t(1-t)$. Integrating the bound $t^5(1-t)^5h^3/(h+t(1-t))^9$ first in $h$ proves convergence as the flat cutoffs are removed. These estimates apply after the Gaussian angular integrations.

The $h$ integration gives a rational function of $t$. After summing the six paths, Hermite reduction gives $5I_1/(63\pi^5)$ as the sum over these three rows of $R+S\int_0^1 dt/(t^2+a t+b)$:

| $a$ | $b$ | $R$ | $S$ |
|---:|---:|---:|---:|
| $-4$ | $-3$ | $-33125311/513995447802066370560$ | $-374735/1269124562474237952$ |
| $-1$ | $-3$ | $98394703/376567746524778332160$ | $189727/232449226249863168$ |
| $3$ | $-6$ | $5971348121/12162470186183892664320$ | $2461843/1583654972159361024$ |

For $I_2$ every $R,S$ changes sign. Writing $\Delta=a^2-4b$, each last integral is $-2\operatorname{arctanh}(\sqrt{\Delta}/(-2b-a))/\sqrt{\Delta}$. Its first 24 odd Taylor terms and a rational geometric tail give the stated strict bounds. The replay computes the angular numerators from the two exterior products, verifies their signs against the original absolute-coordinate determinants, and differentiates every Hermite primitive back to its rational integrand. [Flux-position integration](code/pentagon_flux_integral.py)

This coefficient can be isolated within the full cubic graph sum linear in this puncture block. On any cycle of length $l$ the $l$ polyvector edge factors have degree $4l$ in $5(l-1)$ independent odd differences. For $l<5$ they vanish identically, including when the flux shares an interaction vertex with another external leg. With five cubic vertices the only surviving topology is the pentagon. At positive heat cutoffs a block incident to the integrated flux position gives a smooth tangent kernel for the remaining external legs, since its delta fixes the free tangent position of the flux. A block on any other external pair misses the chosen tangent diagonal.

This separation persists for the first delta derivative after removing the cutoffs. Test the tangent difference $d=z_0-z_4$ with $d_k\rho(|d|^2/\eta^2)$, where $\rho$ is a smooth compactly supported function equal to one near zero, and send $\eta$ to zero. A delta derivative returns minus its $k$th coefficient. With the flux on the $b_i$ end of the special line its angular form vanishes since $U\wedge b_i=0$ and holomorphic normal derivatives of $b_i$ are scalar multiples of $b_i$. With the flux on the cutoff end, only its other neighboring flat edge can become singular. In the normal coordinate difference $u$ from this neighbor its leading scalar denominator is $(|d|^2+|u|^2)^6$; its two numerator factors are antiholomorphic differences. The leading normal integration kills every factor containing $\bar u$ by phase symmetry. For tangent indices $i,j$ it gives
$$
\int_{\mathbb C^3}\frac{\bar d_i\bar d_j\,d^6u}{(|d|^2+|u|^2)^6}
=\frac{\pi^3}{60}\frac{\bar d_i\bar d_j}{|d|^6}\qquad(d\ne0).
$$
Integration by parts on the special delta can take one holomorphic tangent derivative of this expression. This leading term has rotation character $-3$ in $d$ at every radial heat cutoff. Its pairing with the character-$1$ test is zero. Taylor remainders of the smooth factors gain at least one power of $(d,u)$; their pushforwards have scaling degree at most four in the four real tangent coordinates and vanish on the shrinking linear test. Thus these other placements cannot cancel the measured derivative. This determines the distributional coefficient for one specified puncture block before the three scalar sources and the test are transferred to physical states. [Exact angular and radial moments](code/pentagon_flux_integral.py)

### Constraint on a componentwise completion of the propagator

Suppose a correction to the ambient Dolbeault homotopy acts identically on all polyvector components and commutes with holomorphic BV divergence (with the total-complex signs). Suppose its coefficients, after testing in the tangent variables, are smooth at distinct nonzero normal positions $x,y$. If the correction tends to zero as $|y|\to\infty$ at every fixed normal difference $x-y\ne0$, it vanishes at **all distinct normal positions**.

Indeed, applying the commutator to a single normal vector polarization extracts its normal coefficient derivative. Integration by parts on compactly supported input fields gives $(\partial_{x_i}+\partial_{y_i})H(x,y)=0$ for every coefficient $H$ and each of the three normal directions. Fix $d=x-y\ne0$. Then $a\mapsto H(a+d,a)$ is antiholomorphic on $\mathbb C^3\setminus\{0,-d\}$. It extends through both points. In a polydisk containing one missing point, the Cauchy integral on a circle in the first conjugate coordinate extends the function across the center. It agrees with the old function for every nonzero transverse pair of coordinates, and hence on the full punctured overlap. The resulting entire antiholomorphic function tends to zero at infinity and is zero by Liouville on each complex line. The same argument applies after testing a tangent delta or its derivatives.

In particular the smooth exact-pair block used in the pentagon sum cannot occur with a nonzero **net** off-diagonal correction in a completion satisfying all these hypotheses. Its separated coefficient would cancel there. The argument leaves open kernels supported at coincident normal coordinates, additional singular supports, and propagators acting on the polyvector components. The alternating Čech contraction below commutes with coefficient derivatives on Laurent monomials; its chart cochains have a different kernel domain. This constraint concerns completions with the specified smoothness and decay. It gives no vanishing theorem for the full punctured IIB amplitude.

### A cyclic homotopy on smooth polynomial angular modes

Use the unit Čech form $U=2\sum_{\rm cyclic} r_i\,\bar\partial r_j\wedge\bar\partial r_k$ and $b_i=r_j\bar\partial r_k-r_k\bar\partial r_j$, with $r_i=\bar w_i/|w|^2$. Work on forms generated by polynomials in $w,r$ and by the three $\bar\partial r_i$. Their relations include $w_i r_i=1$, $r_i b_i=0$ and $w_i\bar\partial r_i=0$ (normal coordinate contractions). Two-forms are multiples of $U$, one-forms are combinations of the $b_i$, and
$$
\bar\partial b_i=w_i U,\qquad
\bar\partial r_i\wedge b_j=(w_j r_i-\delta_{ij})U/2,\qquad
b_i\wedge b_j=\varepsilon_{ijk}r_k U/2.
$$
Here the coordinate Kronecker symbol denotes the identity on this fixed three-chart coordinate frame.

For polynomial bookkeeping take $w,r$ independent before restricting $w_i r_i=1$. Decompose every polynomial into terms homogeneous of degrees $(k,m)$ in $w,r$ annihilated by the trace $\sum_i\partial_{w_i}\partial_{r_i}$. In unconstrained variables the decomposition multiplies these terms by powers of $w_i r_i$; on the normal space this factor is one. Existence and uniqueness follow by taking successive traces using $\Delta((w_i r_i)P)=(3+k+m)P+(w_i r_i)\Delta P$ on a bihomogeneous $P$, where $\Delta$ is the polynomial trace. For each trace-free term $P$ set
$$
h(PU)=\frac{2}{k(m+2)}\sum_{i\in\mathrm{normal\ coordinates}}(\partial_{w_i}P)b_i\quad(k>0),
\qquad h(PU)=0\quad(k=0).
$$
The displayed partial derivative holds $r$ fixed. The degree-two projection keeps exactly $k=0$, the normal Laurent modes. The degree-zero projection keeps exactly $m=0$, the regular holomorphic modes. For $k>0$, applying $\bar\partial$ to the numerator in $h$ gives $[k(m+2)P-\Delta P]U/2$, hence $\bar\partial h(PU)=PU$.

This $h$ commutes with the actual holomorphic coefficient derivatives. On degree $(k,m)$ their action on $PU$ is $(\partial_{w_i}P-(m+3)r_iP)U$; on $P b_j$ it uses $m+2$ instead. For trace-free $P$ the trace-free part of $r_iP$ is $r_iP-(w_jr_j)\partial_{w_i}P/(k+m+2)$. Substitution and $r_jb_j=0$ give $\partial_{w_i}h(PU)=h\partial_{w_i}(PU)$, including $k=0,1$.

For a one-form $A_j b_j$ take the polynomial $\varepsilon_{ijk}r_k\partial_{w_i}A_j$. For each of its trace-free components of bidegree $(k,m)$ with $m>0$, divide by $-m(k+2)$ and add the results to obtain $h(A_j b_j)$. Set $h=0$ on functions. This formula is the graded transpose of $h$ on two-forms for the linking residue. Indeed the unit residue of $w^\alpha r^\beta U$ is zero unless $\alpha=\beta$, in which case it is $2\alpha!/(|\alpha|+2)!$ by the spherical monomial integral. It pairs a trace-free component only with its reversed bidegree: every index must contract with one belonging to the other polynomial. Also $\varepsilon_{ijk}r_k\partial_{w_i}$ integrates to zero, by this same factorial identity. Integration by parts in the two-form formula now gives the stated one-form formula. The residue of a total holomorphic coefficient derivative vanishes by the factorial identity, so commutation with these derivatives extends by transposition.

The second application of $h$ vanishes because the curl of each polynomial gradient vanishes. On functions the contraction identity follows by transposition from the one on top forms. For a one-form $\omega$, the remainder $\omega-\bar\partial h\omega-h\bar\partial\omega$ is closed and killed by $h$. Every closed one-form in this model has a polynomial primitive. On $w_i\ne0$ integrate it along the affine fiber $w_jr_j=1$ from $r=e_i/w_i$. The differences of the three primitives are regular holomorphic Čech one-cocycles and glue after subtraction of chart functions by $H^1=0$. A regular function on all three charts is polynomial on the affine hypersurface: expanding a large power of $1=(w_jr_j)$ clears the chart denominators. The function contraction identity now forces the remainder to vanish. This proves the cyclic contraction on the whole polynomial angular complex. Its field domain and linking residue are explicit; radial bulk propagation is still required for a loop on the fivefold.

For the unit flux this model gives an independent finite normal trace. With $M_i$ multiplication by $w_i$, $\iota\pi$ the above projection, and $\operatorname{Alt}$ the signed sum over the three coordinate assignments, the operator
$$
\operatorname{Alt}\,[\iota\pi,M_1][h,M_2][h,M_3]\,M_U
$$
has rank one and supertrace $-1$. Its image consists of regular functions. Equivariance under unitary normal coordinate changes allows only complete traces of the $r$ indices to map a bidegree $(k,m)$ to a regular polynomial of degree $k-m$, so every trace-free component with $m>0$ maps to zero. On a regular homogeneous $H$ of degree $k>0$, the coefficients of $b_i$ after the first commutator with $M_l$ are $A_i=(\delta_{il}H-w_l\partial_{w_i}H/k)/(k+1)$. The next commutator gives $r_i V_i$ with
$$
V=\frac{-H e_j\times e_l+(w_l e_j-w_j e_l)\times\nabla H/k}{(k+1)(k+2)}
$$
in the fixed normal frame. The last commutator on this function gives $V_i/(k+3)-w_i\operatorname{div} V/[(k+2)(k+3)]$. The alternating sums of $V_i$ and $w_i\operatorname{div}V$ are $-2H/[(k+1)(k+2)]$ and $-2H/(k+1)$ respectively, hence cancel. On $H=1$ a single ordered contraction sends $U$ to $b_l$, then to $-\varepsilon_{jli}r_i/2$, then to $-\varepsilon_{jli}/6$ for the selected final coordinate (no sum in the last component). The alternating sum sends $1$ to $-1$. This reproduces the earlier cutoff residue with the ordinary smooth wedge product and a holomorphic-derivative-compatible contraction. [Exact differential, residue pairing and trace replay](code/normal_flux_homotopy.py)

## Sewing the local primary vertices

For the same theta assignment include the polyvector divergence on the puncture block. Its tangent derivative differentiates the tangent delta. Factor out the nine nonzero flat scalar kernel weights, the flux radial coefficient and the common puncture-block scalar. On the other nine edges both angular factors contract with their flat displacement. On the extra edge the coefficient of a derivative of the delta in coordinate direction $j$ contracts the odd five-form with the $j$th basis direction. In the order $(z_1,z_2,w_1,w_2,w_3)$ the five relative odd determinants are
$$
(1584,-1026,558,1116,-1674).
$$
Their pairing with the flat displacement on the extra edge vanishes. Only the first two entries multiply derivatives of the tangent delta. Put the flux at vertex 1 and a test at vertex 2 with normal trivector and antiholomorphic form $d\bar z_1\wedge d\bar z_2\wedge d\bar w_1$. The two external polyvectors fill the common Berezin volume. Including all five common antiholomorphic directions gives the angular determinant $-2232$. Thus the first tangent delta derivative has a coefficient of absolute magnitude $2232\cdot1584=3535488$ in the stated kernel factors. Normal derivatives of the smooth puncture block have lower delta order and cannot cancel this principal coefficient in this edge assignment. This is a nonzero polarized graph current for a specified exact-pair block. A cyclic propagator for the full theory, its full sum over edges, and its value on the five prescribed source modes remain to be determined. [Independent rows, including both angular factors](code/puncture_theta_form.py)

At nonzero $N$, include the unary flux bracket in the free differential. At most one normal $(0,2)$ flux may enter a single local holomorphic vertex. The minimum nonlinear rooted arities at genera $0,1,2,3$ are consequently $2,4,5,7$ for the chosen interaction. A rooted tree on $r$ inputs obeys $r=1+\sum_{v\in\mathrm{vertices}}(k_v-1)$ for rooted arities $k_v$. At genus two and $r=5$ only a single shifted genus-two vertex fits. At genus three and $r=7$ a genus-one plus genus-two pair would need eight inputs and three genus-one vertices would need ten. Only the shifted genus-three vertex fits. These counts hold for every genus-independent linear homotopy on the punctured fields. On these primary slots the shifted vertex commutes in graded degree with the linearized differential. Its inputs containing descendants vanish, the divergence of its output vanishes by $\partial_\Omega^2=0$, and the flux bracket on either side would wedge two normal $(0,2)$ flux jets. The measured class is therefore independent of representatives in this contribution from trees sewn out of the fixed local vertices.

At first order in $N$ this tree exclusion extends to every $g$ on $2g+1$ regular sources. The rooted arity of an unshifted genus-$h$ vertex is at least $2h+2$. One flux at a nonlinear vertex reduces the required arity by one. Hence a sewn tree with $V$ nonlinear vertices requires at least $1+2g+V-1=2g+V$ sources. Putting the flux on a lifted external leg or on a propagator preserves the unshifted lower bound $2g+V+1$. At four loops and nine sources, only the single shifted local vertex occurs among these trees.

#### The genus-three insertion on a separated pentagon in the mixed ninth-order pole

Take the nine scalars of the mixed four-loop measurement, three copies each of $e^{w_1},e^{w_2}$ and the three mixed inputs $w_1^2e^{z_1}/2,w_2e^{z_2},w_2z_2e^{-z_1}$. In a pentagon with one genus-three primary vertex let it absorb seven regular scalars and have two internal ends. The other four sites carry the remaining two scalars, the D3 flux, and the test. If one internal momentum is $q$, its other end at this vertex balances $q$ and the seven scalar momenta. All seven lie in the hyperplane with zero $dw_3$ component. Each five-row minor in the nine-field trace therefore has a factor $q_3$, the $dw_3$ component of $q$, with quotient independent of $q$. The six-minor trace gives precisely $q_3^6$ times the input coefficient. Thus the inserted scalar momentum polynomial applies six $w_3$ derivatives on one of the two internal lines. The pointwise coefficient multiplying those derivatives depends only on $w_1,w_2,z_1,z_2$.

Exactly nine of the 36 labeled seven-scalar sets yield the zero vertex. Three sets contain six pure waves and only one mixed wave, so even including $q$ their momenta span at most four dimensions. Eighteen sets contain five pure and two mixed waves. The three collinear pure waves saturate one normal direction and reduce the vertex to $3!W_3$ times the sixth power of the two tangent momenta's determinant and $q_3^6$. Here $W_3$ is the fourfold primary on two equal axial waves and four other marks projecting to the unit transverse simplex. Independent tangent waves permit a shear removing their normal components; polynomiality extends the formula to dependent tangent waves. Generate the third mixed wave by differentiating $w_2 e^{-z_1+Bz_2}$ once in $B$ at zero. In its six pairings with the first mixed wave the tangent determinant contributes $B^6$, so these six vanish. The twelve pairings of the other two types have unit squared tangent determinant.

The remaining fifteen sets contain all three mixed waves. In their position polynomials keep the highest normal monomial $w_1^2w_2^2 z_2/2$ multiplying the product of the exponentials and $q_3^6$. This term selects zero normal momenta at each mixed slot and $B=0$ in the trace. At $q_3=1$ all five directions span, every momentum is nonzero, and at most three waves lie on one line. The worldsheet primary is strictly positive. In the period basis $1,\zeta(3)$ its two coefficients for pure multiplicities $(3,1)$ or $(1,3)$ are $841479311/35307207936000,-53359/2724321600$. For $(2,2)$ they are $763894480043/14405340837888000,-156171991/3705077376000$. These 27 nonzero insertions are local factors inside separated graphs. Their four remaining vertex integrations and the full punctured propagator still enter the physical residue. [Exact nine-field trace and all 126 minors](code/separated_mixed_pentagon.py)

#### Integrated highest-pole test with one puncture block

Use exactly one copy of the cyclic finite-rank puncture block $b_i(x)l_i(y)$ plus its graded transpose. The $l_i$ are supported in a fixed normal annulus. Take every other line to be the flat propagator with proper times in $[\epsilon,L]$, $0<\epsilon<L<\infty$. Integrate every bulk position except the point testing the output, with the actual holomorphic scalar waves of the mixed ninth-pole measurement. For the resulting primary $(0,2)$ vector coefficient on the punctured normal space the limiting integral of $w_3^8 d^3w$ against its restriction to the oriented sphere $|w|=R$ is **zero**. This holds separately for every pentagon, every placement of all nine scalar inputs, and both orientations of the special line.

In fact the same argument applies at every genus $g\ge2$ with $2g+1$ scalar inputs, one flux and a positive number of geometric loops. Only primary vertices enter by the degree count below. If removing the special edge disconnected the graph, the other $4(E-1)=5(V-2)+1$ odd differences would lie in a space of dimension $5(V-2)$, giving zero. The remaining ordinary lines consequently contain a spanning tree. At fixed test point their Gaussian decay controls every relative position and every polynomial times holomorphic exponential from these sources, uniformly for test coordinates in a compact set. The tangent delta and its derivatives at the puncture line remove a tangent integral or act on a smooth factor. Internally use derivatives of the locally integrable $U,b_i$ as distributions, with their derivatives moved onto the other smooth factors. At distinct normal integration variables these are tensor products. At one vertex the wedge of $U$ with $b_i$ and their holomorphic derivatives vanishes on the punctured domain: all are horizontal to the same antiholomorphic normal Euler vector. This test extends that zero product by zero. These prescriptions give smooth functions of the remaining test coordinates except for singular factors incident directly to the test vertex.

At genus $h$ of that vertex, the rooted primary momentum polynomial uses $2h$ five-row determinants. After eliminating the test momentum each determinant is affine in each other slot, since a single momentum variation has rank one. There are at most $2h$ derivatives on the singular slot. The flux form has coefficient scaling $R^{-5}$ and $b_i$ has coefficient scaling $R^{-3}$, with each further holomorphic normal derivative costing at most $R^{-1}$. The tangent flux bivector and the test fill the common odd trace, so the final divergence differentiates only tangentially. A test at $l_i$ vanishes near the brane. With a flux at the test the worst power is $R^{-5-2h}$; with an angular $b_i$ at the test the possible BV derivative on its own line gives $R^{-4-2h}$. All other singular factors have already been integrated. Pair against $w_3^{2g}d^3w$ on the link, supplying $R^{2g+5}$ including the two antiholomorphic differentials. Since $h\le g-1$, the absolute value is bounded by a fixed constant times $R^2$ or $R^3$, respectively, for $0<R<1$. In contrast the actual local contact containing $\partial_{w_3}^{2g}U$ has normalized moment $(2g)!$.

The shrinking-link limit here is taken at fixed proper-time cutoffs. If all ordinary vertex coincidences are excluded by a smooth cutoff at positive separation, Gaussian decay away from the diagonal also permits $\epsilon\downarrow0$ before the linking test there. For the unrestricted renormalized operation, simultaneous approach of ordinary vertices to the brane and the remaining sectors of the punctured propagator still require evaluation. [All 79,200 ordered four-loop assignments, angular identities and local residue control](code/pentagon_linking_residue.py)

#### Differentiated angular puncture modes on the normal link

The same finite-cutoff test includes an arbitrary finite number of differentiated exact pairs. Let $f(w)$ be a homogeneous holomorphic test of degree $m$, let $\psi$ be a smooth normal $(0,1)$-form, and take any multi-index $\kappa$ of holomorphic derivatives, of total order $k$. In the normal coordinates of the preceding calculation,
$$
\partial_w^\kappa b_i=(-1)^k(k+1)!\,\frac{\bar w^\kappa}{s^k}\,b_i .
$$
On $|w|=R$, a Taylor coefficient of $\psi$ of bidegree $(p,q)$ gives $f\,d^3w\wedge\partial_w^\kappa b_i\wedge\psi$ the radial power $m+2-k+p+q$. Under the common phase of the three coordinates its weight is $m-k+p-q$. Sphere integration kills every nonzero weight. Every surviving term therefore has power $2+2q$. Taylor expansion through degree $\max(0,k-m)$ also bounds the smooth remainder by a constant times $R^3$. In particular the integral is $O(R^2)$ for **every** $k$, including $k>m$. The constant in this estimate can depend on the mode.

With the sphere integral normalized by the unit residue of $d^3w\,U$, the result for $f=w^\kappa$ and $\psi=d\bar w_i$ is exactly
$$
(-1)^k\frac{\kappa_1!\kappa_2!\kappa_3!\,(k+2-\kappa_i)}{(k+2)(k+3)}R^2 .
$$
Indeed the restriction of $b_i\wedge d\bar w_i$ divided by $U$ is $(s-|w_i|^2)/2$. For the actual test $w_3^8$, the eighth derivative of $b_3$ gives $(8064/11)R^2$ with this smooth one-form. The eighth derivative of $U$ gives the radius-independent $8!$.

There are genuine exact pairs at arbitrarily high order already in the ray $\partial_{w_3}^k b_1$. Their Dolbeault derivatives equal $w_1\partial_{w_3}^k U$. Pair these images with $\bar w_1w_3^j\bar\partial\chi$ on the original annulus. Relative to the seed normalization, the pairing vanishes unless $j=k$, in which case it is $(-1)^k3k!/(k+3)$. Separate coordinate phases also kill the pairing on every $\partial_w^\kappa U$. Dividing by the diagonal value gives duals to any finite selection of the ray. Horizontal wedges and the repeated $\bar\partial\chi$ give the same cyclic transpose construction as at order zero. At a test endpoint on the annulus the link still misses the support. At an angular endpoint the preceding bound absorbs both vertex derivatives and the BV derivative. A flux at a separate test vertex retains the original $R^2$ bound. The highest residue thus vanishes at finite cutoffs for one insertion of any such finite block. The estimate is termwise. Summing infinitely many modes and extending across ultraviolet collisions require additional bounds. [495 exact link coefficients and 169 annular pairings](code/pentagon_linking_residue.py)

#### Summing the differentiated puncture sector on the physical waves

At fixed $0<\epsilon<L<\infty$ the ordinary heat edges and the nine prescribed waves are entire functions of the normal coordinates and their independent complexified conjugates. This permits summation of **all** holomorphic derivatives of the three exact seed pairs in the integrated linking measurement.

The corresponding exact image at order $k$ consists of $H(w,r)U$, where $r=\bar w/|w|^2$ and $H$ is trace-free of bidegree $(1,k)$. Its primitive is $2(\partial_{w_i}H)(w,r)b_i/(k+2)$ with the coordinate contraction, taking the polynomial derivative before substituting $r$. Normalize $H$ on the unit sphere by the unit-residue measure. Its annular dual is $\overline{H(y,\bar y)}\bar\partial\chi(y)$ divided by $\int |y|^2\bar\partial\chi\wedge U\wedge d^3y$. Indeed the exact image has radial power $1-k$ relative to $U$ and the dual polynomial has degree $k+1$. Sphere orthogonality gives the identity on this harmonic space and kills the Laurent representatives. The unitary reproducing kernel for this angular space, with independent target variables $(w,r)$, is
$$
\frac{(k+1)(k+2)(k+3)}2
\left[(w\cdot\bar y)(r\cdot y)^k
-\frac{k}{k+2}(w\cdot r)|y|^2(r\cdot y)^{k-1}\right].
$$
At $k=0$ omit the second term. Harmonicity fixes the subtraction, and its unit diagonal equals the dimension $(k+1)(k+3)$. Equivalently the displayed kernel acts as the identity on the harmonic projections of $w_i\bar w^\kappa$, as the unit-sphere factorial moments verify directly. Differentiating in $w$ to take the primitive kills the trace subtraction by $r_i b_i=0$. Hence the angular-to-annular kernel at degree $k$, before its common annular normalization, is
$$
(k+1)(k+3)(r(x)\cdot y)^k
\left(\sum_{i\in\mathrm{normal\ coordinates}}\bar y_i b_i(x)\right)\bar\partial\chi(y).
$$

Let $A$ bound $|y|$ in the cutoff annulus and fix $B>\max(A,1)$. At an integrated angular endpoint replace $(r(x)\cdot y)^k b_i$ by $(-1)^k(y\cdot\partial_x)^k b_i/(k+1)!$ and integrate by parts in the normal variable $x$. Apply the one-variable Cauchy estimate on the disk of normal shifts of radius $B$ in direction $y/|y|$. All remaining factors depending on $x$ are entire there. On this disk each finite-scale Gaussian is bounded by a cutoff-dependent constant times a Gaussian with doubled variance on the original real coordinates. The remaining vertex integral converges with that bound as in the finite-block test. The flux at another vertex is a locally integrable distribution of fixed derivative order. A flux at the angular endpoint has the previous zero horizontal product. Thus the integrated order-$k$ contribution is bounded by a constant times a fixed polynomial in $k$ times $(A/B)^k$. Fixed vertex and edge derivatives at either endpoint only enlarge that polynomial.

At an angular test endpoint, $d$ further holomorphic normal derivatives give coefficient degree $-3-k-d$ and central weight $-2-k-d$ including the differential. A smooth one-form Taylor term can survive only if $p-q=k+d-m$. Cauchy estimates of the entire smooth coefficients on a polydisk of radius $B$ then bound the signed angular integral by $R^2$ times a fixed polynomial in $k$ times $(A/B)^k$, uniformly for sufficiently small $R$. The finitely many low modes satisfy the preceding quadratic estimate separately. With the flux at a separate test endpoint the internal bound multiplies the original flux bound $O(R^2)$. An annular test endpoint misses the link. These estimates include the cyclic transpose. Consequently the sum over this whole sector and the shrinking-link limit commute for the actual pentagon waves at fixed cutoffs, and its highest residue is zero.

For example, take the third unit axis in the angular polynomial's $y$ slot, use $w_3^m\exp(pw_3)$ as test, and pair with $d\bar w_3$. Dividing out the annular factor, the sum after sphere integration is exactly
$$
2R^2\int_0^1 t^{m+1}\exp(pt)\,dt .
$$
The coefficient of $p^n$ comes solely from $k=m+n$ and equals $2R^2/[n!(m+n+2)]$. Before angular testing the kernel series need not converge when $|r(x)\cdot y|\ge1$.

More generally fix any positive holomorphic degree $a$. Take an orthonormal basis of trace-free polynomials $H$ of bidegree $(a,k)$ on the unit sphere. The primitive is $2(\partial_{w_i}H)(w,r)b_i/[a(k+2)]$, as follows from the polynomial flux contraction. Its dual uses $\overline{H(y,\bar y)}|y|^{2-2a}\bar\partial\chi$ with the same radial normalization. Since $H(w,r)$ has radial power $a-k$, their paired radial weight is again $|y|^2$. The angular endpoint now has an additional $a-1$ in both its radial power and its phase weight. With $d$ holomorphic normal derivatives the nonzero Taylor terms consequently obey $p-q=k+d+1-a-m$, and still have radial power $2+2q$.

For completeness the bounds are summable in $k$ for every fixed $a$. The number of modes at $k$ is bounded by $\binom{a+2}2\binom{k+2}2$. The supremum of any unit-norm mode on the sphere is bounded by the square root of this number, by unitary invariance of the sum of squared basis values. On a unit great circle its restriction is a trigonometric polynomial of degree at most $a+k$. Termwise differentiation of its finite Fourier series costs a polynomial in $a+k$ for any fixed derivative order. Together with homogeneity for radial derivatives this bounds all ambient derivatives up to order $a$ polynomially in $k$ on the sphere. Thus in $H(w,\bar w)=\sum_{|\alpha|=a} w^\alpha P_\alpha(\bar w)$ each $P_\alpha$ has polynomial supremum, by taking $\alpha$ holomorphic derivatives. On the unit coordinate polydisk the homogeneous degree-$k$ polynomial $P_\alpha$ is bounded by $3^{k/2}$ times this supremum. Cauchy's inequality bounds every monomial coefficient by the same quantity. There are polynomially many monomials. At an integrated angular endpoint their $r^\kappa b_i$ again give derivatives divided by $(k+1)!$, now acting on a smooth factor times a polynomial of degree at most $a-1$ in $x$. Its Gaussian integral remains finite. Choosing the Cauchy radius $B>\sqrt3 A$ absorbs the extra $3^{k/2}$ and proves absolute convergence after testing. Derivatives on the dual produce polynomial factors on the fixed annulus. The same bound applies to the Taylor selection at an external angular endpoint.

For $H(w,r)=w_1^a r_3^k$, the sphere integral of its primitive against $w_3^k\bar w_1^{a-1}d\bar w_1$ times $d^3w$ is exactly $2(a-1)!k!R^{2a}/(a+k+2)!$. Thus the finite-cutoff maximal residue vanishes separately at every holomorphic harmonic degree, including its entire sum over the antiholomorphic degree. [Polynomial flux contraction](code/normal_flux_homotopy.py) · [Harmonic projection, 91 higher-degree tests and pentagon census](code/pentagon_linking_residue.py)

The two degrees can also be summed **simultaneously at an angular test endpoint**. On the unit sphere $\partial_{w_i}H$ is a spherical harmonic of total degree $a+k-1$, since the ambient derivative commutes with the Laplacian. Multiplication by a coordinate changes spherical harmonic degree by one up or down. The same holds for differentiating a radial power times a spherical harmonic. Indeed differentiating its homogeneous harmonic polynomial lowers its degree by one, while differentiating the radial power inserts one coordinate. The coefficient of $b_i$ adds one antiholomorphic coordinate, and restriction to the link adds at most one more coordinate by the normal contraction. After $d$ extra holomorphic derivatives the angular coefficient has no harmonic degree below $a+k-3-d$. Thus, in addition to the central phase rule, a smooth Taylor term on the degree-$m$ test can contribute only for
$$
p+q\ge a+k-3-d-m .
$$
Let $A_0>0$ and $A$ be inner and outer radius bounds for the fixed annulus. Its dual factor has magnitude bounded by $A_0^{-a} A^k$ times a polynomial in $a+k$ and a common constant. This follows directly by evaluating the unit sphere harmonics on $y/|y|$. Normal derivatives of fixed order on either endpoint cost a polynomial in $a+k$, by the great-circle and radial derivative bounds. For $R<1$ the phase-selected power $R^{2+2q}$ is at most $R^2$. Entire Taylor coefficients on a normal polydisk of radius $B>\max(1,A_0^{-1},A)$ give a geometrically summable tail in $p+q$. By the degree inequality the signed sphere integral is bounded by
$$
C R^2(1+a+k)^D (A_0 B)^{-a}(A/B)^k ,
$$
where $C,D$ are independent of $a,k,R$. The polynomially many modes at each bidegree can be absorbed into $D$.

The same absolute summation after angular testing holds for an integrated angular vertex restricted to a bounded normal ball $|x|\le M$, with fixed $M\ge1$. Its three holomorphic volume differentials and its smooth complementary $(0,2)$-form impose $p-q=k+d+2-a$. Including the radial volume its surviving power is $R^{3+2q}\,dR$. This controls $R<1$. For $1\le R\le M$ the radial prefactor before the smooth Taylor term is bounded by $M^{a+1}$, and the entire Taylor bound uses $(M/B)^{p+q}$. Taking $B$ larger than both $M^2/A_0$ and $MA$ gives geometric summability there as well. The radial integrals with holomorphic derivatives agree with the distributional prescription: a contact term at $x=0$ with coefficient homogeneity $a-k-d-4$ would require a delta derivative of order $k+d-a-2$, whereas its coefficient has phase $a-k-d-2$ of absolute weight four larger whenever the order is nonnegative. This excludes a contact of the required homogeneity and phase. A smooth compactly supported radial localization gives the same estimate.

Consequently the exact angular modes of the polynomial flux complex $H(w,r)U$, summed over both degrees and including their transpose, give zero maximal link residue at finite heat cutoffs with a bounded normal localization on the angular vertex. If the angular vertex is the output itself, the same conclusion holds without any extra localization. The following path calculation controls the unrestricted integrated vertex of the pentagon. The other propagator sectors and the joint ultraviolet collision limit still require a calculation. The replay includes $1112$ exact moments killed by spherical harmonic degree despite satisfying central phase, alongside $4973$ nonzero controls. [Two-degree tests](code/pentagon_linking_residue.py)

#### Removing the angular-vertex localization for the separated pentagon

For this pentagon the bound is uniform in the normal integration radius $M$. Cut the unique puncture edge. The remaining ordinary graph is a four-edge path from the angular endpoint $x$ to the annular endpoint $y$. If the test is at $y$ it misses the annulus on a sufficiently small link. The angular test at $x$ was covered by the preceding two-degree estimate. Otherwise the fixed test normal coordinate $w$ separates $x$ and $y$ on this path. The tangent delta on the special line has no normal component. Thus integrating the intermediate vertices creates no Gaussian cross term between $x$ and $y$.

At fixed proper times the smooth normal factor on the $x$ arm, term by term after vertex differentiation, is a polynomial of bounded degree times
$$
\exp[-c|x|^2+u\cdot x+v\cdot\bar x],
\qquad c>0,\qquad v=c\lambda w,\quad 0\le\lambda\le1 .
$$
The nine scalar sources are polynomials times holomorphic exponentials, so their slopes change only $u$. Uniformly on $[\epsilon,L]$ the precision $c$ is bounded above and below by positive constants, and $u$ is bounded for $|w|\le1$. Without an intervening flux the Gaussian is centered at $w$. For an internal flux first move its finite holomorphic derivatives to smooth neighboring factors and write $|f|^{-6}=\int_0^\infty s^2e^{-s|f|^2}\,ds/2$ at its normal position $f$. If the flux lies between $x$ and $w$, let $\alpha,\beta$ be the inverse of four times the summed proper times from $x$ to $f$ and $f$ to $w$. Completing the square in $f$ gives
$$
c=\frac{\alpha(\beta+s)}{\alpha+\beta+s},
\qquad \lambda=\frac{\beta}{\beta+s}.
$$
The determinant contributes $(\alpha+\beta+s)^{-3}$. The numerator $\bar f$ of $U$ is either a Gaussian mean or belongs to a covariance with a polynomial factor. Both carry a further $(\alpha+\beta+s)^{-1}$ times a bounded-degree polynomial in the endpoints, leaving an integrable bound proportional to $s^2/(1+s)^4$. The same bound applies on the other arm when the flux lies there. A flux at the annular endpoint is away from zero. A flux at the test retains its previously bounded external derivatives, and a flux at the angular endpoint has zero horizontal wedge.

To integrate $x$ over its entire normal space, take a monomial of the smooth polynomial factor of bidegree $(j,s')$. Taylor-expand $\exp(u\cdot x)$ and $\exp(v\cdot\bar x)$ at orders $p,q$. For a normalized mode $H$ of bidegree $(a,k)$ and $d$ further holomorphic angular derivatives, central phase requires
$$
p+j-q-s'=k+d+2-a .
$$
The spherical harmonic bound gives $p+j+q+s'\ge a+k-3-d$. Hence $q\ge a-D$ and $p\ge k-D$ with $D$ uniform in both harmonic degrees. Including $d^6x$, the remaining radial power is $3+2(q+s')$. Up to the fixed surface volume its integrated magnitude, including both Taylor coefficients, is bounded by
$$
\frac{\Gamma(q+s'+2)}{2q!\,c^{s'+2}}\,
\frac{|u|^p}{p!}\left(\frac{|v|}{c}\right)^q .
$$
The gamma ratio is polynomial in $q$. This bound also holds with any radial cutoff valued between zero and one. Choose $0<R_*<\min(1,A_0)$ with $A_0$ the inner annular radius. Then $|w|\le R_*$ gives $(|v|/c)^q\le R_*^{a-D}$, canceling the annular growth $A_0^{-a}$ with a geometric factor $(R_*/A_0)^a$. The annular growth in $k$ is bounded by $A^k$ for the outer radius $A$, and $p\ge k-D$ makes its sum converge by the remaining $p!$. The angular norms, degeneracy, and the finitely many endpoint derivatives cost only polynomials in $a+k+p$. This is a uniform bound after angular integration, including the entire radial integral at infinity. It justifies removing $M$ and summing both degrees before the linking limit.

An exact normal-moment check, with the two endpoint slopes in distinct coordinate directions, is
$$
\int_{\mathbb C^3}
\frac{x_1^{a-1}\bar x_3^{k+1}}{|x|^{2k+4}}
e^{-c|x|^2+u x_3+v\bar x_1}\,d^6x
=\frac{\pi^3a!}{(a+k+2)!}\frac{u^{k+1}v^{a-1}}{c^{a+1}} .
$$
The angular integral selects exactly $p=k+1$, $q=a-1$. Independently, introducing a Schwinger parameter $t$ for the denominator gives the beta integral of $[t/(c+t)]^{k+1}[c/(c+t)]^a$ in the variable $t/(c+t)$. Both evaluations use the preceding distributional prescription at $x=0$.

The integrated angular sector is thus uniformly bounded in the remaining test coordinate. If the flux sits at the test, its original $O(R^2)$ highest-pole bound applies. The other integrated placements are smooth or bounded there. Including the angular test and the cyclic transpose, all exact angular modes of the polynomial flux complex, with **no normal localization**, contribute zero maximal residue on the nine actual sources at finite proper-time cutoffs. This covers the pentagon linear in this block. The remaining propagator sectors and the collision limit $\epsilon\downarrow0$ still enter the renormalized boundary comparison. [72 path reductions, 210 independent unrestricted Gaussian moments, and 79,200 assignments](code/pentagon_linking_residue.py)

#### A regulated simultaneous collision in the axial ninth-pole test

Take each of the four ordinary proper times in $[\tau,2\tau]$, then send $\tau$ and the normal link radius $R$ to zero. The axial source replaces the third mixed scalar by $w_2z_1^2z_2/2$. This is a fixed polynomial times a holomorphic exponential, so the preceding cut-path Gaussian calculation also applies to this source.

For an interior test the path to the annular endpoint has total proper time at most $8\tau$. Its normal Gaussian supplies $\exp[-(A_0-R)^2/(32\tau)]$. Derivatives and flux integration lose at most powers of $\tau$ and allow replacement of $32$ by $64$ in this bound: the positive radial flux mass only increases the minimized Gaussian action; its determinant and numerator give the integrable $s$-bound above after rescaling $s$ by $\tau$. On the angular arm $|u|\le C+R/(4\tau)$, with $C$ uniform because holomorphic source shifts are convex combinations of finitely many slopes. The sum of the $(a,k)$ polynomial bounds against $A^k |u|^p/p!$ for $p\ge k-D$ costs at most a polynomial in $\tau^{-1}$ times $\exp[C'(1+R/\tau)]$. The earlier geometric bound in $a$ still holds for $R<A_0/2$. The term $C'R/\tau$ is absorbed by the annular suppression for sufficiently small $R$. After the exterior linking factors the result on this simultaneous-collision region is therefore bounded by $C''\tau^{-D'}e^{-\delta/\tau}$ for some $\delta>0$, uniformly as the link shrinks. A flux at the test retains the original additional $R^2$ bound. Integrating the proper times over the displayed shell preserves the exponential bound, as does summing dyadically to zero. This excludes a primitive simultaneous collision for the entire angular sector on the $47\,520$ interior-test assignments. Another $15\,840$ assignments have annular test and miss the link by support. For any **fixed finite** angular block, the remaining angular-test placements are also suppressed: Taylor expansion to the fixed mode order leaves a finite number of heat derivatives across the normal annulus, each bounded by a power of $\tau^{-1}$ times that same exponential. This last finite-mode argument does not pass directly to the infinite tower.

To test the failure of a uniform estimate at an angular test, keep the genuine $(1,k)$ puncture block, the annular flux input $y_j U$, and the normalized **scalar normal heat kernel** $(c/\pi)^3e^{-c|x-y|^2}$. Here $c$ is inverse heat width. This collision model removes the pentagon's tensor numerators and tangent integrations. Let $d\mu(\rho)$ be the radial measure induced by $\bar\partial\chi\wedge U\wedge d^3y$ on the annulus and put $\mu_j=\int\rho^j d\mu(\rho)$. The block normalization is $\mu_2$. The output is $F(|x|^2)b_j(x)$ by unitary covariance and $\bar x_i b_i=0$. To compute the transverse coefficient take $x=R e_3$, $j=1$. In the annular integral $|y_1|^2$ multiplies the angular series $(k+1)(k+3)(y_3/R)^k$. With $L=c\rho^2$ and $\theta=cR^2$, its direction-averaged value including the unnormalized heat exponential is
$$
D=2\rho^2e^{-L-\theta}
\sum_{k,j\ge0}
\frac{(k+1)(k+3)L^{k+j}\theta^j}{j!(k+j+3)!},
\qquad F(R^2)=\frac{c^3}{\pi^3\mu_2}\int D\,d\mu(\rho).
$$
Indeed the $j$ holomorphic powers of the heat exponential require $k+j$ antiholomorphic powers, and the unit-sphere moment is $2(k+j)!/(k+j+3)!$. Equivalently let $X,J$ be independent Poisson variables with means $L,\theta$. Then
$$
D=\frac{2\rho^2}{L^3}\,
\mathbb E\big[(X-J)(X-J-2)\,\mathbf1_{X-J\ge3}\big].
$$
Without the indicator the expectation is $L^2-(2\theta+1)L+\theta^2+3\theta$. At fixed $\theta>0$ the discarded tail and its first eight $\theta$ derivatives are exponentially small uniformly on the annulus. For example the exponential weight $e^{J-X}$ bounds the lower-tail probability by $e^2\exp[L(e^{-1}-1)+\theta(e-1)]$; polynomial moments and fixed Poisson differences preserve exponential decay. Thus $\lim_{c\to\infty}cD=2$ with $R^2=\theta/c$, whereas each fixed $k$ gives zero even after multiplying by any fixed power of $c$.

The axial eighth-link calibration of this model pairs $w_3^8d^3w$ against $\partial_{w_3}^8(F b_3)\wedge d\bar w_3$. In unit-residue normalization it is exactly $s^{11}(d/ds)^8(F(s)/s^2)/495$ at $s=|w|^2$. Denote this moment by $M(c,\theta)$ at $s=\theta/c$. Substitution of the Poisson moments gives the exact limits
$$
\lim_{c\to\infty}\frac{M(c,\theta)}c
=\frac{2\cdot9!\,\theta\,\mu_0}{495\pi^3\mu_2},
\qquad
\lim_{c\to\infty}\left[M(c,\theta)-
\frac{2\cdot9!\,\theta\,\mu_0}{495\pi^3\mu_2}c\right]
=-\frac{(2\cdot9!\,\theta+4\cdot8!\,\theta^2)\mu_{-2}}{495\pi^3\mu_2}.
$$
The second limit varies with the annular profile and the ratio $\theta$. It supplies a counterexample to removing the cutoff mode by mode at an angular test. The eight derivatives here calibrate the target residue; the physical genus-three insertion and the tensor edges distribute their derivatives according to the pentagon graph. The actual remaining axial Cartan coefficient is still constrained by $+4590\,igNpx^8\,72C_4(e^1)$ if the transferred operation is a boundary differential image. This scalar model does not determine its value. [260 independent sphere/Poisson coefficients and the eighth derivative](code/pentagon_linking_residue.py) · [Axial boundary matching test](#two-polynomial-insertions-determine-the-four-loop-highest-pole-on-this-family)

#### Exterior tensors at the angular test on the axial sources

The actual pentagon eliminates $13\,248$ of the $15\,840$ angular-test assignments. Label the cut path from the angular test $x$ to the annulus $y$ by $0,1,2,3,4$. The tangent delta sets $z_x=z_y$. For the four ordinary lines write $D_l$ for their conditional mean antiholomorphic displacement, in the order $(z_1,z_2,w_1,w_2,w_3)$. Holomorphic vertex derivatives of ordinary kernels insert antiholomorphic polynomials. Introduce the source plane waves before their polynomial coefficient extractions. An interior flux can be handled by its radial Schwinger mass. Conditional Gaussian averages of all remaining antiholomorphic polynomial factors are obtained by substituting their means, because the antiholomorphic–antiholomorphic covariance is zero. The flux mass acts only in the normal directions. The two tangent components of $D_l$ divided by their proper times consequently have successive jumps given, up to common sign and heat units, by the total tangent source momenta at the three interior vertices. Their unweighted displacements sum to zero. A holomorphic tangent derivative of the special delta differentiates the smooth test against it. It cannot differentiate the antiholomorphic tangent endpoint in these means.

The ordinary Dolbeault four-forms on the path each leave one displacement direction free. On these four directions the two tangent delta covectors give the two rows of tangent components of the $D_l$. Every normal form of the puncture, flux and linking test leaves these rows unchanged. Thus at least two independent tangent momenta must be injected at **distinct interior vertices**. The only mixed marks have tangent momenta $(1,0),(0,1),(A,B)$, with $[A^2B]$ required on the axial source. Their interior rank is tested exactly by the three evaluations $(A,B)=(0,0),(1,0),(0,1)$ since every minor is affine in $A,B$. First $6048$ assignments have flux at the test, where $U$ wedges to zero with the angular one-form and its holomorphic derivatives. An additional $6744$ have deficient tangent rank. In another $240$ the genus-three vertex has fewer than two marks among the three mixed scalars, flux and test. The other marked waves span only the two pure normal directions, and momentum conservation leaves at most one further independent internal momentum. Every five-row minor then vanishes.

Another $216$ assignments have flux at the annulus and the genus-three vertex in the interior. The rank tests force two mixed waves at that vertex and the third at a different interior vertex, with a pure normal wave at the remaining one. Its local genus-three trace is proportional to the sixth power of the third normal current and is independent of all normal shifts of the mixed waves: removing the two shifts by a shear preserves the pure normal waves, and the two mixed tangent momenta must both enter the reduced four-row minor. This is a polynomial identity and extends to dependent tangent momenta. The current rows span the endpoint difference $\Delta=\bar y-\bar x$ in the normal directions and the three interior momenta. Every nonzero four-row minor therefore has exactly one factor linear in $\Delta$. Denote its cofactor in the five-dimensional odd tensor by $C$. Its normal components use both mixed tangent directions and the pure normal momentum. Its tangent components use that pure normal momentum and at most one mixed normal shift. The Dolbeault determinant with annular flux uses only the **normal** components of $C$, since its two tangent rows are supplied by the special delta.

Under common phase rotation of the normal coordinates the product of the angular and dual mode coefficients has weight $-2$ before their differential factors; this holds at every bidegree. The annular radial form wedged with the flux pairs $y$ with $\bar y$ at weight zero. Projecting onto the link at $x$ adds a normal covector with coefficient $x$ of weight $+1$. The two four-row minors each supply weight $-1$ from $\Delta$, and the genus-three current supplies $-6$. The test $w_3^8$ then leaves net weight $-1$ for a tangent special BV derivative and $-2$ for a normal derivative. In the tangent channel at most one of the four requested mixed normal shifts can hit a tensor coefficient; in the normal channel none can. Thus at least three or four holomorphic powers of endpoint coordinates, respectively, must come from the external exponentials. Their common phase cannot balance in either channel. This eliminates all $216$ after angular integration at every heat scale.

The $2592$ assignments left by these tests have $720$ for each of the three interior flux sites and $432$ with annular flux. In $360$ the interior tangent rank drops to the first tangent coordinate at $B=0$. Every normal odd cofactor also needs both tangent rows and therefore supplies a second $B$ beyond the one in the Dolbeault tensor. The same holds for the odd cofactor completed by the first tangent direction. Only the **second tangent** special BV derivative can survive the axial $[B]$ projection on these $360$.

There are also admissible placements with the opposite selection. Take flux at $y$ and put the genus-three vertex at $x$ with the test, three first-direction pure waves, two second-direction pure waves and the first mixed wave. The three interior source momenta in path order are $dz_2+t\,dw_2$, $dw_2$, and $A\,dz_1+B\,dz_2+u\,dw_2$. Replace the four displacement rows by $\Delta$ and these three jumps, suppressing the nonzero change-of-basis determinant. In this row order their signed odd cofactor is exactly
$$
C=(0,0,-A\Delta_3,0,A\Delta_1).
$$
Here $\Delta_i$ denotes the $i$th normal component. In the $(1,k)$ angular sector the root one-form has coefficient $\bar y\times\bar x$ times its scalar radial and harmonic factors. After wedging the annular radial covector with the flux, the Dolbeault determinant has columns given by the two tangent coordinates, $(\bar y\times\bar x)\cdot D_{l,\mathrm{normal}}$, and $x\cdot D_{l,\mathrm{normal}}$. The normal BV derivative at $x$ adds an odd contraction with $\bar x$, since the holomorphic derivative of the radial kernel is proportional to $r_i=\bar x_i/|x|^2$. Multiplying that Dolbeault determinant by $C_{\mathrm{normal}}\cdot\bar x$ gives, with the same row order,
$$
A^2(x\cdot\Delta)
(\bar x_1\Delta_3-\bar x_3\Delta_1)^2 .
$$
Both tangent BV components vanish identically. The displayed normal tensor is nonzero away from alignment and has a quadratic transverse zero when $\bar y$ and $\bar x$ are parallel. At $x=R e_3$ its transverse factor is $R^2\bar y_1^2$; the $|y_1|^2$ angular test of the scalar collision model samples a different component. This is a tensor-level survivor. The differentiated genus-three vertex, remaining angular and proper-time integrations, and sum of placements still determine its axial Cartan residue. [Exact census, four-row cofactors and Dolbeault determinant](code/pentagon_linking_residue.py)

#### Integrated cancellation of the annular-flux pentagons on the axial source

The angular integration and the genus-three vertex eliminate **all $432$ remaining annular-flux assignments** above. They also eliminate $432$ with interior flux. The mechanism fixes the full normal Taylor degree, so it applies at arbitrary positive proper times and link radius.

Rotate the three normal coordinates independently by $g=\operatorname{diag}(g_1,g_2,g_3)$ with $|g_i|=1$, and transform every incoming normal momentum by $g^{-1}$. After conditioning on the normal endpoints, the normal displacement columns and the normal flux vector transform by $g^{-1}$. The odd determinant of the four current rows completed by the special-edge BV derivative has character $(\det g)^{-1}$, also when that derivative acts on the radial profile. The rooted Dolbeault determinant has the same character. To check it, transform the three common normal column vectors by $g^{-1}$. The angular covector with coefficient $\bar y\times\bar x$ contributes $(\det g)^{-1}$ on these columns. So does insertion of the flux vector in the normal volume form. The link covector with coefficient $x$ and the annular radial covector with coefficient $y$ contract invariantly. Returning the three common columns to their fixed basis multiplies the determinant by $\det g$. This counting applies to the full projector at every harmonic bidegree: its scalar contractions are unitary invariant, and the primitive's remaining normal covector is along $\bar y\times\bar x$. Indeed the unitary scalar projector with two vector arguments depends on their pairwise contractions; differentiating its target holomorphic variable leaves only $\bar y_i b_i$, since $r_i b_i=0$. Differentiating a tensor is equivariant with its additional covector index. The same character calculation therefore covers derivatives distributed over the tensor and its scalar factors.

The genus-three local trace contains six five-row momentum determinants and contributes $(\det g)^{-6}$. The linking test $w_3^8$ supplies $g_3^8$. Consequently the angularly integrated graph as a function of the incoming scalar normal momenta has degree $(8,8,0)$ in the three normal directions. This uses a Fourier coefficient on the compact phase group and does not require a radial limit. The marked mixed polynomials use two derivatives in each of the first two normal directions. The remaining **six pure scalar waves have total degrees exactly $(6,6)$**. These degrees count both vertex momenta and all Taylor terms of propagated exponentials.

Whenever the high vertex carries five pure scalars and one mixed scalar, its six fixed momenta are either three copies of $dw_1$ and two of $dw_2$, or the reverse, followed by the mixed momentum $M$. Call two exceptional momenta $q_1,q_2$. The third is eliminated by conservation. All nine momenta are then linear combinations of the five columns $dw_1,dw_2,M,q_1,q_2$ with constant coefficients on the six prescribed unit waves. Equivalently, allowing independent strengths for the repeated pure waves, every nonzero five-row minor contains one pure factor from **each** of the first two normal axes: expand the eliminated row first, then use that the fixed marks span three directions. Its common determinant is
$$
\det(dw_1,dw_2,M,q_1,q_2)
=(q_2)_3[M_T,(q_1)_T]-(q_1)_3[M_T,(q_2)_T]
$$
in the corresponding ordered frame, where $[\, ,\,]$ is the tangent two-row determinant and $T$ labels tangent projection. The normal shift of $M$ drops out. Thus the six determinant factors spend all six allowed pure powers on each axis. They differentiate exceptional normal legs only in the third direction. If the high vertex instead carries five pure scalars and two mixed scalars, the same count follows from its four fixed directions and one independent internal momentum. In both cases the surviving link uses the constant Taylor term of the remaining pure wave and of all pure-wave exponentials outside the differentiated vertex factors.

If the remaining pure scalar is alone at an ordinary interior vertex away from the flux, its current jump is therefore zero. Two ordinary odd edge rows become proportional and the complete odd determinant vanishes. This kills the $432$ annular-flux assignments, including the preceding nonzero *unintegrated* normal cofactor. It also kills $144$ assignments at each of the three interior flux sites. These are zeros of the physical integrated subfamilies at each positive time, on separated as well as approaching configurations. Their regulated simultaneous and hierarchical collision limits remain zero. The other $1728$ candidate assignments have interior flux: at each of the three sites, $216$ have a constant pure insertion at the annular endpoint and $360$ have four pure marks at the high vertex. [Exact census, normal determinant character and 504 five-row minor comparisons](code/pentagon_linking_residue.py)

There is an additional endpoint-specific zero for the $648$ constant annular insertions. At that cubic vertex the ordinary BCOV propagator is a polyvector divergence $\partial_\Omega P$. If the special BV derivative acts at the annular end as well, its kernel there has the form $\partial_\Omega h$. The full local vertex integral is then proportional to
$$
\int\operatorname{Tr}\big((\partial_\Omega P)(\partial_\Omega h)\big)
=\pm\int\operatorname{Tr}\big(P\partial_\Omega^2 h\big)=0 .
$$
Integration by parts in the holomorphic and odd variables has no normal boundary on the support of the annular factor. Tangentially the pairing is distributional with the compactly supported relative delta. This proves the cancellation of that derivative placement at all cutoffs. If the special derivative is applied at the opposite endpoint, obtaining the same conclusion requires divergence closure at the annular end. This identity alone does not establish that closure for a separately truncated puncture block. The full physical axial coefficient still requires the other candidate assignments, propagator sectors and renormalized transfer.

#### Flux scale in the simultaneous angular-test collision

All $1728$ remaining axial assignments have interior flux. Put the four positive proper times in $[\tau,2\tau]$ and the angular test radius in $\theta_0\tau\le R^2\le\theta_1\tau$ with fixed $0<\theta_0<\theta_1$. At the interior flux position $f$ use $|f|^{-6}=\tfrac12\int_0^\infty s^2e^{-s|f|^2}\,ds$. For any fixed $0<b\le1$, the entire contribution with $s\ge\tau^{-b}$ vanishes faster than every power of $\tau$, including both angular harmonic degrees and the actual tensor numerators.

Let $\alpha,\beta$ be the inverses of four times the summed proper times on the two arms from $f$ to the angular point $x$ and annular point $y$, respectively. Completing the square gives the normal exponential
$$
\exp[-A R^2-B\rho^2+2C\operatorname{Re}(\bar x\cdot y)],\qquad
C=\frac{\alpha\beta}{\alpha+\beta+s},\quad
A-C=\frac{\alpha s}{\alpha+\beta+s},\quad B-C=\frac{\beta s}{\alpha+\beta+s},
$$
where $\rho=|y|$ lies in the fixed compact annulus away from zero. The normal degree $(8,8,0)$ proved above allows each normal source exponential to be truncated before taking the integrated axial coefficient. Move derivatives off the flux onto the ordinary kernels. The remaining conditional Gaussian moments are polynomials of bounded degree in the endpoints, with coefficients bounded by powers of $\tau^{-1}$. The factor $\bar f$ of the flux supplies an extra $(\alpha+\beta+s)^{-1}$ from its Gaussian mean or its contraction against a holomorphic insertion. Including the normal Gaussian determinant, the $s$ measure is bounded by a power of $\tau^{-1}$ times $s^2(\alpha+\beta+s)^{-4}\,ds$.

The interior flux integration itself can be done at arbitrary positive times. Put $h=\alpha+\beta$, take independent Gaussian sources $J,\bar J$ for $\bar f,f$, and put $d=J\cdot\bar J/h$. The substitution $u=h/(h+s)$ gives the exact generating kernel
$$
\begin{aligned}
\int_{\mathbb C^3}e^{-h|f|^2+J\cdot\bar f+\bar J\cdot f}
\frac{\bar f_i}{|f|^6}\,d^6f
&=\frac{\pi^3\bar J_i}{2h}\int_0^1(1-u)^2e^{ud}\,du\\
&=\frac{\pi^3\bar J_i}{h}\sum_{j\ge0}\frac{d^j}{(j+3)!}.
\end{aligned}
$$
Differentiate in the two sources for all polynomial insertions, evaluate at $J=\alpha x+\beta y$ and its conjugate, and restore $e^{-\alpha R^2-\beta\rho^2}$. In particular with $m=(\alpha x+\beta y)/h$ the normalized Gaussian bridge average of the flux coefficient is
$$
\frac{\bar m_i}{|m|^6}
\left[1-e^{-h|m|^2}\left(1+h|m|^2+\frac{h^2|m|^4}{2}\right)\right].
$$
Its value at $m=0$ is zero. This integration uses the normal part of the actual flux at all three interior sites.

At angular bidegree $(a,k)$ the radial dependence of the paired mode contributes $(R/\rho)^{a-k}$ times fixed powers of $R,\rho$. On the unit annular sphere the dual is harmonic of degree $a+k$ before the finitely many derivatives and polynomial factors. Its product decomposes into spherical degrees $l$ with $|l-a-k|\le D$, for one fixed $D$. The angular norms, projection on each degree and mode multiplicity cost at most a polynomial in $a+k$. With normalized sphere measure the eigenvalue of the remaining exponential on degree $l$ is
$$
2\sum_{j\ge0}\frac{(CR\rho)^{l+2j}}{j!(l+j+2)!}.
$$
Rotation covariance reduces this identity to the harmonic $y_1^l$ on the unit sphere at the first real axis. The $j$th term uses the moment $\int_{S^5}|y_1|^{2(l+j)}=2/[(l+j+1)(l+j+2)]$.

Put $X=CR^2$ and $Y=C\rho^2$. The numerator including the radial weight is exactly
$$
(R/\rho)^{2a}(\rho/R)^{a+k-l}\,Y^{l+j} X^j .
$$
The second factor costs at most $(\rho/R)^D$ for $R<\rho$. At any fixed shift $l-a-k$ and fixed $l$, summing the first factor over $a\ge a_0$ costs at most $(R/\rho)^{2a_0}/[1-(R/\rho)^2]$. Summing $l+j$ and $j$ with their factorials gives at most $e^{Y+X}$ times a fixed polynomial in $1+X+Y$ by the exponential moment identity. Thus the integrated tail has a bound, for constants $K,d,\delta>0$ independent of $s,a_0,\tau$, of
$$
K\tau^{-d}(R/\rho)^{2a_0}
\frac{s^2\,ds}{(\alpha+\beta+s)^4}
\exp[-(B-C)\rho^2-(A-C)R^2].
$$
Here the displayed bound is prior to radial annulus integration. Since $(B-C)\rho^2\ge\delta\min(s,\tau^{-1})$, integration gives the claimed suppression. At each fixed power of $\tau$ only finitely many holomorphic degrees $a$ can consequently contribute, even at bounded $s$.

For any fixed $0<\eta<1/2$, the terms with $|k-C\rho^2|>\tau^{-1/2-\eta}$ also vanish faster than every power after integration. Indeed the majorant for $l+j,j$ divided by $e^{Y+X}$ is the product of Poisson weights with means $Y,X$, multiplied by a bounded-degree polynomial. Here $X$ is bounded and $Y$ is at most a constant times $\tau^{-1}$. If $|(l+j)-Y|\le\tau^{-1/2-\eta}/2$ in this tail, then $a+j\ge\tau^{-1/2-\eta}/2-D$. The geometric factor in $a$ and the bounded-mean Poisson tail in $j$ eliminate those terms. In the complementary region the Poisson bound at deviation $v$ is $2\exp[-v^2/(2(Y+v))]$, which absorbs every polynomial cost.

Both estimates allow an arbitrary hierarchy of individual proper times **within either arm**, provided each arm sum stays between fixed positive multiples of $\tau$. A normalized ordinary kernel has two antiholomorphic displacement factors divided by its edge time, and each additional holomorphic derivative can insert another such factor. On an arm with total time $T$, label its vertices in path order, let $T_i$ be the elapsed time at vertex $i$, and write $Z_i$ for any one of its complex coordinates. Condition first on the two arm endpoints. The conditional mean of $(\bar Z_e-\bar Z_{e-1})/t_e$ is the full conjugate endpoint difference divided by $T$. Its conditional covariance with $Z_i$ is exactly
$$
4\left(\mathbf1_{e\le i}-T_i/T\right).
$$
Its covariance with every antiholomorphic variable vanishes. Thus Wick evaluation of the tensor numerators and holomorphic waves introduces powers of the inverse arm totals alone. The normalized undifferentiated heat kernels on each arm convolve to the kernel at its total time. Tangent directions use the full path or the same two-arm calculation. This gives estimates uniform down to individual zero-time faces and allows integration over the proper-time simplices at fixed arm totals. Writing $T_x,T_y$ for the angular and annular arm totals, the annular suppression exponent at arbitrary arm ratio is exactly $sT_x\rho^2/(T_x+T_y+4sT_xT_y)$. The bounded-mass region in the surviving harmonic window and degenerations of the arm ratio still enter the physical residue. [Exact bridge matrices, conditional covariances and spherical coefficients](code/pentagon_flux_scale.py)

#### Genus-one pentagons of the third source descendant

Test the rank-two component of $E_1$ with the level-one normal trivector source, on four labeled copies of $B=z_1w_1+z_2w_2$ and one explicit flux $NU\,\partial_{z_1}\wedge\partial_{z_2}$. Evaluate at $z_1=z_2=0$ on the normal link with weight $w_3^2$. A pentagon has a single classical quartic vertex: its one level-one input is the distinguished source. It therefore contains exactly one other external field. With a $B$ at this root the vertex has the undifferentiated factor $B(0,w)=0$. A tangent delta derivative on an adjacent puncture edge can be evaluated in the other endpoint, leaving this multiplication at the root. With the flux at the root, an angular root gives the wedge of $U$ with a one-form horizontal to the same normal antiholomorphic radial vector, which is zero. A root in the cutoff annulus misses the shrinking link. In a labeled cyclic ordering with fixed root and one oriented puncture line these alternatives account for $960,48,48$ assignments, respectively, leaving $144$ interior roots with flux. These counts retain both cyclic orders.

Cut the puncture line for one of the remaining terms, and denote by $T_x,T_y$ the ordinary time totals from the fixed interior root to the angular and annular ends. Write $R,\rho$ for the root and annular radii. The four scalar polynomials have bounded degree. At fixed endpoints, Gaussian integration on the ordinary arms gives bounded-degree polynomials multiplying the separate normal Gaussians. The displacement covariances above make their coefficients polynomially bounded in the two inverse arm times, uniformly at short individual-edge faces. Move any special-end holomorphic derivatives at the integrated angular endpoint onto smooth factors. This uses the distributional angular prescription above.

The normalized Gaussian on the angular arm has precision $\alpha=1/(4T_x)$. For a homogeneous spherical harmonic $H_l$ of degree $l$, $\nu\ge1$ and $l+2-\nu\ge0$, its exact action is
$$
\left(\frac{\alpha}{\pi}\right)^3\int_{\mathbb C^3}
e^{-\alpha|x-w|^2}\frac{H_l(x)}{|x|^{2\nu}}\,d^6x
=\frac{\alpha^\nu H_l(w)}{(\nu-1)!}
\int_0^1 u^{\nu-1}(1-u)^{l+2-\nu}e^{-\alpha R^2u}\,du .
$$
First take $2\nu<l+6$, where the ordinary integral at the origin exists. The Schwinger integral rescales the Gaussian center by $1-u$, and convolution of a harmonic polynomial with a centered Gaussian preserves that polynomial. Spherical testing supplies another $l$ powers of radius in the first nonzero Taylor term. The identity thus extends analytically to $\nu<l+3$ with the indicated angular prescription. A contact distribution of the same homogeneity would have derivative order $2\nu-l-6<l$ and could not carry angular degree $l$. The absolute value is consequently bounded by $\|H_l\|_\infty R^{l-2\nu}\Pr\{P\ge\nu\}$, with unit-sphere norm and a Poisson variable $P$ of mean $\alpha R^2$.

For the full exact-pair sector of bidegree $(a,k)$ the angular coefficient has numerator a harmonic derivative of degree $a+k-1$ times one coordinate, and denominator $|x|^{2(k+2)}$. Multiplication by a polynomial of homogeneous degree $d$ followed by harmonic decomposition gives $l=a+k+d-2j$, $\nu=k+2-j$ for $j\le d+1$. Thus $l+2-\nu=a+d-j\ge0$ because $a\ge1$. The bounded shifts cover all tensor insertions. At $\nu\le0$ there are only bounded $k$ and the Gaussian moments are polynomial. The annular harmonic test has the independent Poisson majorant of the preceding calculation, now with means $Y=\rho^2/(4T_y)$ and $X=R^2/(4T_y)$. Up to fixed polynomial costs and bounded shifts its order $k$ equals the difference of the two Poisson indices minus $a$, and its holomorphic-degree factor is $(R/\rho)^{2a}$.

Fix $\delta>0$. In the region
$$
\frac{R^2}{4T_x}\le (1-\delta)\frac{\rho^2}{4T_y},\qquad R/\rho\longrightarrow0 ,
$$
the full integrated angular sum is bounded by a polynomial in $R^{-1},T_x^{-1},T_y^{-1}$ times $\exp[-c_\delta\rho^2/(4T_y)]$ for some $c_\delta>0$. For $k\ge(1-\delta/2)Y$ the angular-arm Poisson upper tail gives this bound. Below that threshold the annular Poisson index of mean $Y$ must deviate by a positive fraction of $Y$, or $a$ or the other annular index must exceed a positive fraction of $Y$. The latter alternatives are suppressed by $(R/\rho)^{2a}$ and $X/Y=R^2/\rho^2$, respectively. Polynomial moment factors preserve the exponential bound. This proves zero contribution to every power for algebraic simultaneous limits with $T_y\to0$ in the displayed region, including hierarchies of individual ordinary times. The much shorter angular arm with $T_x$ of order $R^2T_y/\rho^2$ or smaller, and a nonshrinking annular time, still require evaluation. [Root census and independent sphere, radial and Schwinger moments](code/pentagon_flux_scale.py)

#### Difference of two genus-one linking measurements

At tangent origin use the same rank-three source test at level two and take the unweighted link. Keep the term linear in $c$ on four mixed inputs $B=z_1w_1+z_2w_2$ and one axial input $c w_3^2/2$, still at first order in the flux $N$. In a pentagon the root is now the classical quintic trace, with two external partners besides the source. A mixed partner multiplies the root by zero. The sole remaining pair consists of the flux and the axial scalar. Removing the axial leaf leaves the preceding level-one quartic graph, with the same kernels, ordered exterior factors and source orientation. Its integrand acquires precisely $c w_3^2/2$ at the root.

The polarized classical descendant weights are $2!/2!=1$ at the quintic and $1!/1!=1$ at the quartic. The axial field is even; the shift between source levels changes parity by zero. With four labeled mixed fields, fourteen of the fifteen root pairs vanish, accounting for $3360$ of $3600$ assignments with one oriented puncture line. The other $240$ map bijectively to the previous flux-at-root assignments, including each of its $144$ interior roots. Taylor evaluation of the four repeated mixed fields divides both coefficients by $4!$. Consequently the unweighted graph sum linear in $c$ is exactly $c/2$ times the weighted level-one graph sum. This equality holds for each set of regulated proper times, also in the short angular-arm corner. It uses only the local multiplication at the root and the common regulator. It applies with additional componentwise puncture lines whenever their regulated products exist. A cycle shorter than five vanishes because its $4l$ relative odd factors occupy at most $5(l-1)$ dimensions for length $l<5$.

In worldsheet action units the local genus-one primary Euler derivatives are $E_1^{(1)}=\det(\partial^2 A_0)/24$ and $E_2^{(1)}=A_0E_1^{(1)}$. Take $A_0=B+c w_3^2/2+Nf\,\partial_{z_1}\wedge\partial_{z_2}$ with $f=(w_1w_2w_3)^{-1}$. The mixed Hessian has transverse determinant one. Adding the axial Hessian gives fivefold determinant $c$; the flux dependence of $E_1^{(1)}$ is $N\partial_{w_3}^2 f/24$ times the tangent bivector. The two local flux coefficients to compare, with this bivector removed, are
$$
\frac{Nc}{24}\left(f+\frac{w_3^2}{2}\partial_{w_3}^2 f\right)
\quad\hbox{and}\quad
\frac{Nc}{24}\frac{w_3^2}{2}\partial_{w_3}^2 f .
$$
The first is the axial coefficient of the level-two descendant. The second is $c/2$ times the weighted level-one descendant. Their difference has unit-normalized Laurent residue $Nc/24$, using $\iota_{\partial_{z_1}\wedge\partial_{z_2}}\Omega=dw_1\wedge dw_2\wedge dw_3$. All matched pentagons cancel in this difference before removing a common regulator. An independent change of source representatives or of defect contact terms would be additional input to a boundary interpretation. [Exact labeled-graph bijection and Hessian calculation](code/pentagon_flux_scale.py)

#### Two compact linking sources with the flux varied on its own leg

Use four mixed scalar insertions from $B_0=z_1w_1+z_2w_2$, a holomorphic normal polynomial $P$, and the tangent bivector with normal Laurent coefficient $NF_n$, where $F_n=1/(w_1w_2w_3^n)$ and $n\ge1$. In the second experiment multiply the same smooth representative of the bivector by $P$ before contracting it. Use the constant annular weight in both experiments. In the preceding Euler convention, writing $\mathcal M_r$ for the resulting genus-one, rank-two linking measurement at source level $r$, the inherited graphs give
$$
\mathcal M_2(B_0^4,P;NF_n)-\mathcal M_1(B_0^4;NPF_n)
=\frac{N}{12}\operatorname{Res}\!\left(F_n\partial_{w_3}^2P\right)
=\frac{Nn(n+1)}{12}[w_3^{n+1}]P(0,0,w_3).
$$
Here four identical mixed fields carry their Taylor factorial, and the extra scalar is polarized once. Both experiments use $D(({-t})^r\chi\tau_z\eta_w)$ with the annular orientation fixed above, $t=-\zeta$, and $D=\bar\partial+t\partial_\Omega$. Using $t^r$ itself in both potentials changes the displayed subtraction to addition. The level-two rank-three input is external to the minimal field complex. All scalar and flux inputs are free closed. The flux corrections to the mixed scalar lifts exist since $w_1 F_n,w_2 F_n$ vanish in normal Laurent cohomology. They can be multiplied by $P$ for the second experiment.

At a classical pentagon root a mixed scalar factor vanishes at tangent origin. A smooth compact rotation-invariant tangent top form of integral one selects exactly this value: translating every integrated tangent coordinate by the root makes the remaining dependence polynomial and holomorphic in the root. Each nonzero quintic root carries $P$ and the flux together. It equals the quartic root with the composite leg $PF_n$ at the same position and with no derivative on either factor. The equality is before integration at a common regulator, including graphs with componentwise puncture lines. For the rank-three source with explicit flux there are $V=4\ell+1$ vertices and $E=V+\ell-1$ edges at geometric loop number $\ell$. Using a lifted scalar or the rank-two term of the compact source each lowers $V$ by one. Every such extra graph has more than $5(V-1)$ odd differences in its $4E$ edge factors, or has no vertex at $\ell=0$. Independent defect vertices and kernels mixing polyvector components are outside this identity.

The contact is the difference of the coefficients from $A_0\det(\partial^2 A_0)/24$ and $\det(\partial^2 A_0)/24$. The two tangent rows of the mixed Hessian select the $w_1,w_2$ columns. The only remaining normal derivative is $\partial_{w_3}^2$. Stripping $N$, the contact difference before taking its residue is
$$
\frac{F_n\partial_{w_3}^2P+P\partial_{w_3}^2F_n-\partial_{w_3}^2(PF_n)}{24}
=-\frac{(\partial_{w_3}P)(\partial_{w_3}F_n)}{12}.
$$
For $n=1$ and $P=w_3^2/2$ the two contact residues are $N/12$ and zero. Their puncture graphs can have a nonzero common value. The variation of the free flux representative has a normal Dolbeault primitive $(w_3/2)b_3\,\partial_{z_1}\wedge\partial_{z_2}$, where $\bar\partial b_3=w_3 U$ and $U$ represents the unit flux. Moving its differential across a cycle opens a nonseparating edge with a BV heat kernel or residual projection. Its degree permits a contribution. [Root bijection and independent full Hessians at three pole orders](code/pentagon_flux_scale.py)

The annular experiments also admit compact external sources closed about the **classical scalar background** to the orders measured. Take $\rho B_0+\epsilon P$ with formal scalar strengths $\rho,\epsilon$, and initially take $\tau_z$ to be the tangent delta top form at the origin. The genus-zero derivative of the vertex with a level-$r$ state has primary Hamiltonian coefficient $(\rho B_0+\epsilon P)^{r+1}/(r+1)!$ in positive cotangent coefficients. This follows by taking two derivatives of the classical potential, with $r+1$ scalar slots left. Its divergence on the rank-three potential $\chi\tau_z\eta_w$ differentiates only normal coordinates. Since $z_i\tau_z=0$, that extra term starts at $\epsilon^{r+1}$. Thus applying the true linearized differential to the cutoff potential gives the free annular source modulo $\epsilon^{r+1}$. It is closed there by the classical master equation. In the first experiment $r=2$ and only one $P$ is used. In the second $r=1$ and $\epsilon=0$.

For a smooth representative, let $H$ be the translation-invariant tangent Dolbeault homotopy kernel multiplied by a radial cutoff equal to one near zero, and put $S=1-[\bar\partial,H]$, with graded commutators. The fundamental-solution identity makes $S$ convolution with a smooth compactly supported kernel; the cutoff derivative is supported away from its singularity. Both operate only tangentially. If $D_B$ is the classical linearized differential, $1-[D_B,H]$ commutes with $D_B$. Its difference from $S$ lowers tangent antiholomorphic degree by one and polyvector rank by one: the free holomorphic divergence graded-commutes with $H$, and every remaining Hamiltonian term applies one divergence. The cube of $1-[D_B,H]$ therefore smooths the tangent delta. A term with three degree lowerings vanishes in tangent dimension two, and every other term contains the smoothing convolution. The resulting source is smooth, compact and annular to the stated orders. Its leading rank-three component uses the smooth radial smearing $S^3\tau_z$, of integral one. Every additional component has rank at most two. With at most one explicit flux bivector and scalar backgrounds it cannot fill the five common odd directions in any inherited componentwise graph. This proves the same difference for the classically closed sources. The argument does not impose the quantum Ward identity on an individual puncture measurement.

#### Quantum bulk completion of the annular potential on a quadratic

Use $B=\rho(z_1w_1+z_2w_2)+c w_3^2/2$ and the tangent delta top form $\tau_z$. For the level-$r$ potential $\zeta^r\chi\tau_z\eta_w$ the genus-$h$ primary Hamiltonian applies $\partial_\Omega$ to the transposed bulk Hessian. On the collar $\chi=1$, $B^m\tau_z=(cw_3^2/2)^m\tau_z$, as a distribution before differentiation. In the cofactor operator only $\rho^4\partial_{w_3}^2$ can act on this product: every mixed derivative also differentiates in $w_1$ or $w_2$. Put $k=h-j$ and $m=r+1-2h-j$ in the quadratic Euler formula. For $m\ge k$ the coefficient before divergence is $\tau_z\eta_w$ times
$$
 \frac{\rho^{4h}c^{r+1-2h}(2m)!\,w_3^{2(r+1-3h)}}
 {(2m-2k)!\,2^m m!\,24^j j!\,4^k(2k+1)!}.
$$
It vanishes for $m<k$, including a missing term of negative $m$. Take the unweighted source at $r=3g-1$. For $h>g$ the collar coefficient vanishes. At $h=g$ it is constant in the normal variables, with value
$$
 \frac{\rho^{4g}c^g}{24^g g!}\int_0^1(1+3s^2)^g\,ds .
$$
Indeed $(2k)!/(2k+1)!=1/(2k+1)$ integrates the $s^{2k}$ term. Its divergence is zero. At $h<g$ every term contains $c^{3g-2h}$, hence at least $c^{g+2}$. The full quantum Hamiltonian differential of the cutoff potential is therefore annular through order $c^{g+1}$. For the level-one comparison at $c=0$ the collar also vanishes. These differentials are closed because the local linearized Hamiltonian squares to zero on **every** field about a quadratic. The tangent smoothing by $1-[D,H]$ described above applies to the quantum operator as well. Every correction still lowers polyvector rank by at least one and has finite differential order at a fixed level. Thus the two $N/12$ measurements have smooth compact completions closed in the actual linearization of all flat bulk Hamiltonians. Their puncture Ward identity still includes additional diagrams. [Bulk differential and triangular equivalence](../iib-bcov/README.md#linearized-bulk-fields-on-arbitrary-quadratic-backgrounds) · [Direct collar derivatives through genus nine](../iib-bcov/quadratic_source_residues.py)

For a general extra axial polynomial, write the background as $\rho B_0+\epsilon P(w_3)$. The level-two, genus-one *state* Hessian on the potential, to first order in $\epsilon$, is
$$
 \frac{\rho^4\epsilon}{24}\tau_z\eta_w
 \left(\chi\,\partial_{w_3}^2P+\partial_{w_3}^2(P\chi)\right).
$$
This includes both the determinant and the transposed cofactor in the variation of $B\det(\partial^2B)/24$. The primary differential on the collar is $+\rho^4\epsilon\,\tau_z(\partial_{w_3}^3P)\partial_{w_1}\wedge\partial_{w_2}/12$. The dual coordinate Hamiltonian has the opposite sign. A regular change of the potential starting at this same order in $\epsilon$ cannot cancel a nonzero axial germ by the classical primary differential on the mixed background. Integrate tangentially against the constant holomorphic test and restrict $w_1=w_2=0$. The tangent divergences integrate to zero; the remaining $\partial_{w_3}$ derivative acts on a multiple of $B_0^{k+1}$ for every descendant level $k$, still in the ideal $(w_1,w_2)$. The primary Dolbeault correction integrates to zero in the tangent top degree without a normal antiholomorphic index. For $P=w_3^2/2$ the germ vanishes, as required. For $P=w_3^3/6$ this calculation detects it. A change already present at $\epsilon=0$ or an independent defect term would change the problem.

#### The nonseparating pentagon with an annular residual projection

Move the normal Dolbeault differential from the smooth unit-flux variation primitive $(w_3/2)b_3\,\partial_{z_1}\wedge\partial_{z_2}$ at the level-one root. Consider a resulting **incident residual edge** in the componentwise angular Dolbeault model, with tangent identity delta. Use an annulus for the root strictly inside the annuli dual to the normal modes. A Laurent-output half-edge at the root vanishes since every holomorphic derivative of $U$ is horizontal to the same antiholomorphic radial vector as $b_3$. An annular dual at the root misses its support. For the regular output, summing the Taylor coefficients dual to $x^\alpha$ gives the kernel $U(p)\wedge\bar\partial\chi(p)/(1-x\cdot r(p))^3$, since $(-1)^{|\alpha|}\partial^\alpha U/\alpha!=(|\alpha|+2)!r^\alpha U/(2\alpha!)$. Here $x$ is the root and $r(p)=\bar p/|p|^2$.

With four ordinary lines the regular-output term has zero shrinking-link limit. Their chain between $x$ and the annular point $p$ has total time $T$ and heat density $(4\pi T)^{-5}\exp(-|x-p|^2/(4T))$ after integrating the three intermediate positions. A conjugate edge increment divided by its time has conditional mean equal to the conjugate endpoint displacement divided by $T$. Its covariance with any holomorphic position is bounded by four. Two conjugate increments have zero covariance. Thus the polynomial insertions cost only a polynomial in $1/T$ uniformly over all subdivisions. The separated endpoint radii give exponential decay as $T$ tends to zero. The four time variables at infinity give the integrable bound $T^3 T^{-5}dT$. The root normal measure scales as $R^6$, its annular differential as $R^{-1}$, and the primitive coefficient as $R^{-2}$. Even assigning two further holomorphic derivatives to the latter still leaves a factor of $R$. The regular kernel and its derivatives are uniformly bounded on the shrinking root annulus.

If at least three of the five lines carry tangent identity deltas, the tangent Dolbeault product vanishes. For exactly three, the two ordinary four-forms must each leave one tangent antiholomorphic differential. They contract their tangent two-forms with opposite conjugate displacements after the three delta identifications, so their wedge vanishes. Four identifications leave no relative tangent differential for an ordinary four-form. Five delta volumes already exceed the relative tangent dimension. Holomorphic delta derivatives preserve these antiholomorphic equalities. This leaves a sector with one additional puncture propagator. If its annular end meets the degree-three regular dual, the normal wedge vanishes. When the additional propagator is nonadjacent to the root, each finite harmonic bank at fixed positive heat cutoffs gives zero on the shrinking link by ordinary Gaussian smoothing. Simultaneous ultraviolet limits in these nonadjacent chambers remain open.

Take the adjacent ordering with the regular-output and angular-primitive ends at the root. Denote their two remote annular ends by $p,q$ respectively, and put $T=t_1+t_2+t_3$ for the three ordinary times on the path between them. Both special tangent deltas must be differentiated, one on each attached mixed scalar. This gives $[p,q]=p_1q_2-p_2q_1$. In the directed cycle from the root through $p$, the two intermediate nodes and $q$, the odd trace is $-[p,q]\det(\Delta\bar w_1,\Delta\bar w_2,\Delta\bar w_3)$ before the three inverse edge times. The pulled ordinary Dolbeault tensor is $[\bar z_2,\bar z_3]\,q\cdot(\bar q-\bar p)$ before its inverse edge times. Here $z_2,z_3$ are the two intermediate tangent positions. For the Gaussian convention above, the two-node covariance in a single complex coordinate has determinant $d=16t_1t_2t_3/T$. The tangent contraction against the two intermediate mixed scalars gives $d[w_2,w_3]$. Contracting this normal area with the three conjugate normal increments gives $2d(\bar q_3-\bar p_3)$. The resulting scalar, with the sign of the ordered odd trace factored outside, is
$$
\frac{512}{T^2}(4\pi T)^{-5}e^{-|p-q|^2/(4T)}
[p,q]\,q\cdot(\bar q-\bar p)(\bar q_3-\bar p_3).
$$
The edge normalizations, signs of the two delta derivatives, radial cutoff densities, and other cuts remain separate from this ordered tensor.

For a trace-free harmonic $H$ of bidegree $(a,k)$ the root pairs $w^\alpha(w_3/2)b_3$ with $2(\partial_{w_i}H)b_i/[a(k+2)]$. Angular integration by parts transfers the indicated derivative to $w^\alpha$ and gives the negative of the pairing with
$$
\frac{w_3}{2a(k+2)}
\left(r_2\partial_{w_1}-r_1\partial_{w_2}\right)w^\alpha .
$$
The numerator has zero polynomial trace because the two mixed derivatives commute. Hence $a=1$ and $k=|\alpha|$ are required. At unit radii, $\alpha=(1,0,0)$ and $H=w_2r_3$ give root residue $-1/72$. Their normalized remote angular mean in the constant Taylor coefficient of the Gaussian is $1$. Exchanging directions one and two gives $+1/72$ and $-1$. The angular total in this coefficient is $-1/36$.

In the shell order $|q|<|p|$, the harmonic sum followed by the Taylor sum gives
$$
-\frac{q_3\overline{[p,q]}}
{2|p|^2|q|^2(1-q\cdot\bar p/|p|^2)^3}.
$$
This formula uses unit-sphere harmonic duals. The unintegrated annular dual convention above has the further common factor $|q|^2/\int |y|^2\bar\partial\chi\wedge U\wedge d^3y$. The sum follows from $\sum_{k\ge1}k(k+1)s^{k-1}=2/(1-s)^3$ for $|s|<1$. Set $A=|p|^2$, $B=|q|^2$, and $C=q\cdot\bar p$. After multiplying the scalar polynomial in the preceding path formula, its mean over simultaneous $U(3)$ rotations is
$$
-\frac{(AB-|C|^2)(B-C)^2}{24AB(1-C/A)^3}.
$$
Indeed $\overline{p\times q}$ is Hermitian orthogonal to $p$ and $q$, where the cross product uses the stated normal volume. The fourth moment of a unit vector gives $(AB-|C|^2)(B-C)/12$ for the mean of $|[p,q]|^2 q_3(\bar q_3-\bar p_3)$. At equal radii take the Abel limit from the convergent harmonic sum. The variable $s=C/A$ then has disk density $2(1-|s|^2)d^2s/\pi$. The remaining angular factor is $-A^2(1-|s|^2)/[24(1-s)]$. Its real part is strictly negative in the open disk after multiplication by the positive angular Gaussian at every $T>0$. The equal-shell Abel limit is integrable: with Taylor regulator $b<1$ the pre-limit factor at fixed shells contains $(1-s)^2/(1-bs)^3$, bounded for $b\ge1/2$ by $8/|1-s|$. This supplies a full ordered angular survivor in the split model. The physical Ward coefficient still requires radial and small-time integration, the opposite endpoint order, and the remaining residual and defect terms. [Exact exterior determinants, Gaussian contractions and harmonic tests](code/pentagon_residual_edge.py)

**The full time integral and its radial contact.** Keep the ordered tensor and the same stripped edge conventions. Write $I_\epsilon(A,B)$ for its angular average integrated over the three positive ordinary times with $T\ge\epsilon$. The time simplex supplies $T^2dT/2$, leaving $256\,T^{-5}dT/(4\pi)^5$ times the angular mean just derived. For $0<B<A$ the resulting integral is **exactly zero** at $\epsilon=0$. Indeed set $t=\sqrt{B/A}$ and $s=C/\sqrt{AB}$. Integrating the Gaussian in $T$ gives $6\cdot4^4/[A(1+t^2-2t\operatorname{Re}s)]^4$. The $S^5$ measure has marginal $2(1-|s|^2)d^2s/\pi$. Its product with the exterior tensor supplies a second $1-|s|^2$. The resulting disk weight is a constant multiple of the first complex coordinate marginal $3(1-|s|^2)^2d^2s/\pi$ of the unit sphere $S^7$. The eight-real-dimensional ball Poisson kernel is $(1-t^2)/(1+t^2-2t\operatorname{Re}s)^4$. Its integral against the holomorphic function
$$
f(s)=\frac{(B-\sqrt{AB}\,s)^2}{(1-ts)^3}
$$
is $f(t)=0$. Absolute integration follows from $|p-q|\ge\sqrt A-\sqrt B>0$ and $t<1$. This proves the full cancellation for smooth profiles supported on separated shells in that order.

At $A=B>0$ the fixed-shell and integrated radial limits differ. Taking the equal-shell Abel value above, write $x=1-\operatorname{Re}s$ and $y=\operatorname{Im}s$. The unintegrated angular mean has the factor $-A^2/(12\pi)$ times the Laplace integral in $x\in[0,2]$ with weight $e^{-Ax/(2T)}$ and slice
$$
\int_{-\sqrt{2x-x^2}}^{\sqrt{2x-x^2}}
\frac{x(2x-x^2-y^2)^2}{x^2+y^2}\,dy .
$$
The slice divided by $x^2$ tends to $4\pi$: substitute $y=xv$ and dominate by $4/(1+v^2)$. The angular mean divided by $T^3$ therefore tends to $-16/(3A)$. Consequently
$$
\lim_{\epsilon\downarrow0}\epsilon I_\epsilon(A,A)
=-\frac{4096}{3A(4\pi)^5}.
$$
For a smooth compactly supported radial weight $W(A,B)$ away from zero, include both radial Jacobians and any common dual-profile factors in $W$. On the single ordered side the distributional limit is
$$
\lim_{\epsilon\downarrow0}\sqrt{\epsilon}
\int_{0<B<A}W(A,B)I_\epsilon(A,B)\,dA\,dB
=-\frac{32768}{15\sqrt{\pi}(4\pi)^5}
\int_0^\infty W(A,A)\frac{dA}{\sqrt A}.
$$
To extract it at fixed $A$, set $h=\sqrt{T/A}$, $t=1-hu$ and $s=1-h^2v+ihw$. The domain tends to $u\ge0,\ v\ge w^2/2$. The Gaussian becomes $\exp(-u^2/4-tv/2)$. With the shell Jacobian, the angular coefficient and the time density give $-256T^{-3/2}/[6\pi(4\pi)^5\sqrt A]$ times the limiting integral
$$
\int_{u\ge0,\ v\ge w^2/2}
(2v-w^2)^2\frac{(u+iw)^2}{(u-iw)^3}
e^{-u^2/4-v/2}\,du\,dw\,dv
=64\cdot\frac{2\sqrt\pi}{5}.
$$
The $v$ integral gives $64e^{-w^2/4}$; in the $(u,w)$ half-plane the angular integral is $\int_{-\pi/2}^{\pi/2}\cos(5\theta)d\theta=2/5$. For $t$ bounded below by $t_0>0$, the identity $|1-ts|^2-|t-s|^2=(1-t^2)(1-|s|^2)$ bounds the rescaled rational factor by $1/\sqrt{u^2+t_0^2w^2}$. The Gaussian and the bound $2v-w^2-h^2v^2\le2v$ give integrable domination; the compact radial support supplies $t_0$. Integrating $T^{-3/2}$ at the lower limit supplies the factor two. The coefficient is a one-sided ultraviolet contact for this total-time cutoff. Its cancellation or finite remainder depends on the other orderings, subtractions and intrinsic defect terms. [Moment, time and slice checks](code/pentagon_residual_edge.py)

#### Ordinary ambient lines

For an ordinary finite-scale Feynman graph built from the ambient propagator, write $\ell$ for its cycle number, $n$ for the total genus, $M$ for its external valence, and $s_v$ for the valence of vertex $v$ minus $2g_v+3$, where $g_v$ is its interaction genus. All $s_v\ge0$ in the constructed genera. For $V$ vertices and $E=V+\ell-1$ edges, the half-edge and genus counts give
$$
4E-5(V-1)=2n+3-M+\sum_{v\in\mathrm{vertices}}s_v .
$$
The antiholomorphic part of every ordinary propagator on the edge $(i,j)$ is a scalar times the contraction of the five antiholomorphic differences with $\bar x_i-\bar x_j$. Each is annihilated by simultaneous radial dilation of the $V-1$ relative antiholomorphic coordinates. The product vanishes for $4E\ge5(V-1)$ and $V>1$: equality gives a horizontal top form. Holomorphic vertex derivatives leave the common null vector intact. This is a pointwise identity at finite positive heat cutoffs. For the one-flux primary coefficient $M=2n+3$ every graph with an ordinary ambient loop therefore vanishes. Taking the flux from an external lift gives $M=2n+2$, also zero. More generally if $m$ fluxes are inserted and $u$ of the scalar legs are replaced by their vector lifts, then $M=2n+2+m$. A primary test for the rank-two Euler derivative has Dolbeault degree three. The external degree condition $2m+u+3=5V-4E$ becomes $m+u+\sum_v s_v=1$, so the fixed component cannot receive higher powers of flux from these graphs. Heat-edge Ward graphs have a separate five-form kernel. [Flat kernel factors](../iib-bcov/operator/README.md#signs-of-the-closed-heat-edges-at-two-loops) · [Radial subgraph test in the paper](../../sources/LLC_main.tex#L2043)

For the nonseparating cut of the preceding genus-one pentagon, the **odd** kernel eliminates the term in which all five edges are ordinary, including the BV heat edge. At fixed positions $x_v$, each edge contributes the odd four-form obtained by contracting the five odd differences with the metric covector of $\bar x_u-\bar x_v$. Replacing an ordinary propagator by its heat endpoint changes the Dolbeault degree and preserves this odd contraction. Modulo common odd translations, simultaneous variation of the odd vertex coordinates in their metric-dual antiholomorphic position directions annihilates every edge factor, since contraction twice by an edge displacement is zero. Away from the diagonal this is a nonzero relative direction, whereas the five edges have degree twenty in the twenty relative odd directions. Their product vanishes identically at every positive choice of cutoff times. Equivalently the five covectors are dependent since the spatial differences sum to zero around the pentagon. Holomorphic derivatives of ordinary Gaussian lines preserve the odd null direction. The identical argument applies on any fully ordinary subgraph whose odd degree fills its relative odd space. For four ordinary pentagon edges and one componentwise puncture edge, the remaining odd contraction depends only on the fifth covector modulo the span of the four ordinary displacements. A puncture projection or a singular contact requires its own covector and integral.

Ambient heat smoothing also measures its own failure on puncture cohomology. For multi-index $\alpha$ of normal derivatives, $\bar\partial(\partial_w^\alpha U)$ is the same derivative of the unit normal delta current on ambient space. Pair the smoothed output on a linking sphere of radius $R$ with $w^\alpha$, using normal orientation divided by its full-space residue. For proper time $\tau$ in the heat kernel $\exp(-|w|^2/(4\tau))$ and $x=R^2/(4\tau)$ the fraction of that residue enclosed is
$$
1-e^{-x}\sum_{\substack{0\le j\le|\alpha|+2\\\mathrm{Taylor\ terms}}}\frac{x^j}{j!}.
$$
Indeed Stokes gives the integral of $w^\alpha\partial_w^\alpha K_\tau$ over the normal ball. The holomorphic derivatives supply $(-\bar w/(4\tau))^\alpha$ and the Gaussian radial density acquires the factor $x^{|\alpha|}$. Division by the whole-space moment yields the integral of $v^{|\alpha|+2}e^{-v}/(|\alpha|+2)!$ from $0$ to $x$. Distinct normal multi-indices pair to zero by coordinate phase rotations. For positive $R,\tau$ the enclosed moment varies with $R$. The original Laurent class has the same residue on every linking sphere. Thus a heat transfer preserving Laurent states must account for the normal boundary current. The preceding vanishing determines the contribution with ordinary ambient lines; loops with a puncture-sensitive line remain to be evaluated.

For a homotopy acting identically on polyvector components, each internal edge also contributes four odd differences at its endpoints, including a puncture-sensitive edge with distributional coefficients. On $2g+1$ regular sources at first order in $N$, taking the flux from a source lift gives $4E=5(V-1)+1+\sum_v s_v$. There are only $5(V-1)$ independent differences, so that entire class vanishes before integration. With explicit flux the Laurent output degree requires $\sum_v s_v=0$. Now $4E=5(V-1)$ fills the relative odd volume. The common odd volume must come from the flux tangent bivector and the test. Consequently these graphs give a tangent bivector before the output divergence. At genus four, the possible geometric loop numbers with an explicit flux are $0,1,2,4$: geometric loop number three would require exactly one genus-one primary vertex, which is identically zero. For geometric loop number one, every cycle has length at least five by the degree of the relative odd kernel on that cycle. The five vertices must therefore form the pentagon used above.

### Laurent inputs in a rooted primary bulk operation

Work at descendant level $t^0$. Use $l$ external normalizable Laurent states of Dolbeault degree two, $p$ regular sources of degree zero, and extract order $N^b$ in the D3 flux. For a rooted connected graph let $n$ be its total genus and $s_v$ the excess valence over $2g_v+3$ at each local vertex, as above. Its output Dolbeault degree is
$$
q_{\rm out}=l-p+2n+2+b+\sum_{v\in\mathrm{vertices}}s_v .
$$
Indeed for $b$ explicit flux vertices the integrated test has $M=l+p+b+1$ external slots. The $E$ BV kernels supply Dolbeault degree $4E$ and integration away from the output removes $5(V-1)$, giving $2l+2b+4E-5(V-1)$ before testing. The preceding half-edge identity yields the result. Moving a flux from the diagram onto a lifted external state adds one to its Dolbeault degree: the flux bracket adds two and its homotopy removes one. This preserves the formula. It holds also for puncture kernels of the same total Dolbeault degree, including their tangent delta factors.

The Dolbeault cohomology of the punctured fivefold has only rows zero and two. The projection to a Laurent output consequently requires $p=l+2n+b+\sum_v s_v$; a regular output requires two further sources. In particular **every nonlinear operation with only normalizable external states vanishes after this primary projection**, for any number of flux insertions, in this rooted graph sector. A state of degree greater than two cannot reappear at $t^0$ by perturbing the projection using the flux bracket, since each further bracket with its homotopy increases the degree by one. This count applies to the Taylor operations of the bulk action. A point observable with several marked legs has a different edge count.

#### A four-input tree with two normalizable single-graviton states

At tree level the first Laurent response to two Laurent states requires two regular sources. Turn on the commuting shears $w_1T,w_2T$, where $T=\partial_{z_1}$, with exterior parameters $\epsilon_1,\epsilon_2$ on the left. Take the two states $\partial_{w_3}^a U\,X_H$ and $\partial_{w_3}^b U\,X_K$ for integers $a,b\ge1$, where $X_H=(\partial_{z_1}H)\partial_{z_2}-(\partial_{z_2}H)\partial_{z_1}$ and $U$ is the unit normal flux form. The dictionary sends them to $dH(d\theta_3)^a,dK(d\theta_3)^b$. All tangent Hamiltonian vectors here preserve the flux and have zero divergence.

In the polynomial angular contraction the lift of the first input through both shears is
$$
\partial_{w_3}^a
\left(U X_H+(\epsilon_1 b_1+\epsilon_2 b_2)X_{TH}
+\epsilon_1\epsilon_2 r_3 X_{T^2H}\right).
$$
This follows from $\bar\partial b_i=w_i U$ and $w_1 b_2-w_2 b_1=-\bar\partial r_3$; $\partial_{w_3}$ commutes with both shears. The zeroth and first coefficients of the bracket of the two lifts vanish since $U\wedge U=U\wedge b_i=0$. Thus there are no lower output terms to propagate. The bracket of the shears vanishes as well. The coefficient of $\epsilon_1\epsilon_2$ sums the remaining rooted cubic trees and is
$$
\frac{\partial_{w_3}^{a+b+1}U}{(a+b+3)!}
\left(-a!(b+2)![X_{T^2H},X_K]
      -(a+2)!b![X_H,X_{T^2K}]
      +2(a+1)!(b+1)![X_{TH},X_{TK}]\right).
$$
In the two terms where the shears act on different states, moving one exterior parameter past a one-form gives a minus. Their ordered sum uses $2b_1\wedge b_2=r_3U$. Differentiating $U,b_i,r_3$ successively in $w_3$ gives the factors $(-1)^j(j+2)!/2$, $(-1)^j(j+1)!$, $(-1)^j j!$ times $r_3^j$ times the respective original factor. This fixes all three coefficients.

For $a=b=1$, $H=z_1^3/6$, $K=z_2^2/2$, only the first bracket survives. The coefficient and its boundary state are
$$
\frac{1}{20}\partial_{w_3}^3 U\,T
\quad\longleftrightarrow\quad
-\frac{1}{20}\,dz_2(d\theta_3)^3 .
$$
This represents a nonzero Laurent class in the ordinary linearized flux complex: flux boundaries of regular scalars give at most the simple pole, whereas $\partial_{w_3}^3 f=-6/(w_1w_2w_3^4)$. There is a different result when both shear parameters are included in the differential. Put $Y=z_1\partial_{z_1}-z_2\partial_{z_2}$, so $[T,Y]=T$. Add $\epsilon_1{\cal L}_{w_1T}+\epsilon_2{\cal L}_{w_2T}$ to the linearized differential. A primitive of the full $\epsilon_1\epsilon_2$ output is
$$
\frac{1}{20}\left(
-\epsilon_2\,\partial_{w_1}\partial_{w_3}^3U\,Y
+\epsilon_1\epsilon_2\,\partial_{w_1}\partial_{w_3}^3 b_1\,T\right).
$$
The two terms proportional to $w_1\partial_{w_1}\partial_{w_3}^3 U$ cancel; differentiating the $w_1$ in $\bar\partial b_1$ gives the measured coefficient. In particular this mixed four-input operation is tree-exact in the primary graph family. With $l=p=2$, positive genus would require $q_{\rm out}\ge4$.

One quadratic correction suffices simultaneously for every finite collection of these states. Write the displayed general pair coefficient as $\partial_{w_3}^m U\,X_L$, with $m=a+b+1$ and $L$ its polynomial Hamiltonian. Take the polynomial $J$ with $TJ=L$ and $J(0,z_2)=0$. A primitive is
$$
-\epsilon_2\,\partial_{w_1}\partial_{w_3}^m U\,X_J
+\epsilon_1\epsilon_2\,\partial_{w_1}\partial_{w_3}^m b_1\,X_L .
$$
Its bracket with every lifted input is identically zero. A term with a repeated exterior parameter vanishes, and each remaining term has normal antiholomorphic degree at least three in the two-dimensional annihilator of the normal Euler vector. Holomorphic differentiation preserves that annihilator. Two primitives both contain $\epsilon_2$ and their bracket vanishes as well. Taking the negative primitive as the bilinear correction therefore solves all the classical bracket identities on this family, with no higher terms. It uses a Laurent derivative transverse to the original $w_3$ direction. All terms are tangent Hamiltonian vectors, so $t\partial_\Omega$ and $N\{F,-\}$ vanish on them. With at most these two exterior shear sources and at least one Laurent input, the degree bound $p\ge l+2n+b$ also rules out every positive-loop primary operation in the stated rooted graph family. [Signed tree, coherent quadratic correction and primitive replay](code/two_laurent_transfer.py)

#### Three commuting shears detect the bracket of two Laurent states

Take the bulk shears $g_i T$ with $T=\partial_{z_1}$ and $(g_1,g_2,g_3)=(w_1,w_2,w_3^k)$, $k\ge2$. They commute and preserve the D3 flux at any $N$. In the regular dictionary the first two are super-bivectors and the third has polyvector degree $k+1$. All three belong to the extended bulk symmetry algebra outside the finite superconformal algebra. Pull back the unit angular form $U$ and its primitives $b_i,r_i$ along $y_i=g_i(w)$, and denote the resulting forms by $U_y,b_i(y),r_i(y)$. The Čech class of $U_y$ is $1/(g_1g_2g_3)$. The physical flux representative remains $F$.

Take the vector state $U_y(z_1^2/2)\partial_{z_2}$ and scalar state $z_2U_y$. Both define nonzero ordinary linearized classes. A scalar has no incoming flux bracket at $t=0$. For a vector boundary of form degree two the regular scalar primitive must be holomorphic, and its flux bracket has at most a simple normal pole in each variable. The vector here has pole order $k$ in $w_3$.

Write $D_\epsilon$ for the linearized differential with the three exterior shear parameters on the left. The vector has the closed lift
$$
\left(\tfrac12 z_1^2 U_y
+z_1\sum_{i\in\mathrm{shears}}\epsilon_i b_i(y)
+\sum_{\mathrm{cyclic}\ (i,j,h)}\epsilon_i\epsilon_j r_h(y)\right)\partial_{z_2}.
$$
The scalar is strictly closed. In the vector lift the divergence and flux bracket vanish termwise, and the remaining terms cancel by the angular identities above. The flux bracket with the scalar contains four normal antiholomorphic differentials in three normal dimensions. In the vector-first bracket the coefficients of exterior degrees zero and one vanish by $U_y\wedge U_y=b_i(y)\wedge U_y=0$. The oriented quadratic coefficient at $(i,j)$ is $C_{ij}=r_h(y)U_y$ for cyclic $(i,j,h)$. Since $\partial_{y_h}U=-3r_h U$ before pullback,
$$
[C_{23}]=\frac{1}{3g_1^2g_2g_3},\qquad
[C_{31}]=\frac{1}{3g_1g_2^2g_3},\qquad
[C_{12}]=\frac{1}{3g_1g_2g_3^2}.
$$
Multiplication by the complementary shear and addition gives an identity of forms and a unit normal residue:
$$
g_1 C_{23}+g_2 C_{31}+g_3 C_{12}=U_y,\qquad
\operatorname{Res}_{w=0}\frac{w_3^{k-1}\,dw_1\wedge dw_2\wedge dw_3}{g_1g_2g_3}=1 .
$$
For a possible $D_\epsilon$-boundary reduce modulo $t$ and project to scalars. Every flux bracket raises polyvector rank. The exterior-degree-zero component of its scalar primitive is a closed $(0,3)$ form, which can be removed by a full coboundary since $H^3(X,\mathcal O)=0$. The linear coefficients of the primitive are then closed $(0,2)$ forms $P_i$. The quadratic coefficients of its image have Laurent classes $[g_iTP_j-g_jTP_i]$. Their complementary weighted sum vanishes identically. The unit residue therefore proves nonexactness in the full three-shear complex, also with formal descendants. At $k=2$ the linear Laurent dictionary takes the detected scalar to $dz_1dz_2\,d\theta_3$. This is a bracket of the specified equivariant lifts at genus zero. The bare states alone do not specify their exterior-positive continuations.

In fact the bracket with this fixed vector detects **every finite scalar Laurent class**. For a monomial $1/\prod_i w_i^{h_i}$ with $h_i\ge1$, put $(m_1,m_2,m_3)=(1,1,k)$ and $\alpha_i=\lfloor(h_i-1)/m_i\rfloor$. Use $z_2$ times the following normal two-form $V$ as the second state:
$$
V=\frac{\prod_i w_i^{m_i(\alpha_i+1)-h_i}}{(-1)^{|\alpha|}\prod_i\alpha_i!}\,
    (\partial_y^\alpha U)(g(w)).
$$
Each holomorphic derivative of $U$ is a scalar multiple of $U$ pointwise. The same calculation therefore gives $C_{ij}=r_h(y)V$. From $r_h\partial_y^\alpha U=-\partial_y^{\alpha+e_h}U/(|\alpha|+3)$, where $e_h$ increments the derivative in the complementary coordinate,
$$
[C_{ij}]=\frac{\alpha_h+1}{|\alpha|+3}\frac{[V]}{g_h},\qquad
g_1 C_{23}+g_2 C_{31}+g_3 C_{12}=V .
$$
Division in this expression denotes the displayed Laurent monomial with that additional pole. Pairing against $\prod_i w_i^{h_i-1}$ measures its unit residue after the weighted sum. Different Laurent coefficients can thus be recovered separately.

Every two-shear restriction is exact. Drop the common holomorphic prefactor in $V$ and compute on the $y$ variables before pullback. For the oriented pair $(i,j)$ and remaining coordinate $h$, put $\beta=\alpha+e_h$. A scalar primitive is
$$
\frac{\epsilon_j z_1\,\partial_y^{\beta+e_i}U
      -\epsilon_i\epsilon_j\,\partial_y^{\beta+e_i} b_i}
     {(|\alpha|+3)(\beta_i+1)} .
$$
Its quadratic image follows by differentiating $\bar\partial b_i=y_iU$. Its flux bracket vanishes by normal form degree on the first term and tangent independence on the second. Restoring the prefactor after pullback commutes with its differential. Thus the individual restrictions vanish while their joint class recovers the full Laurent coefficient. With two Laurent inputs and at most three regular sources, the degree count above gives $q_{\rm out}\ge3$ at positive genus for the primary rooted action. The scalar result here uses no positive-genus vertex. — Codex. [Angular identities, pair primitives and 1,372 exact residues](code/three_shear_two_states.py) · [Regular symmetry dictionary](../../code/holography/symmetry.py)

The same detection is possible with **three linear shears and two ordinary nonzero states**. Now set $g_i=w_i$. Use the scalar $z_1^2U/2$ and its lift given by the parenthesis of the preceding vector lift with $y=w$. Its linear coefficient has zero flux bracket by $b_i\wedge U=0$ and the quadratic coefficient $r_h$ is tangent-independent. For a normal derivative $\partial_w^\alpha U$ of any total order $m=|\alpha|\ge1$, the second input is the vector
$$
(\partial_w^\alpha U)
\left(z_1\partial_{z_1}+(m-1)z_2\partial_{z_2}
      +w_1\partial_{w_1}+w_2\partial_{w_2}+w_3\partial_{w_3}\right).
$$
The vector in parentheses commutes with every $w_iT$. Its divergence on the displayed normal two-form is $1+(m-1)+3-(m+3)=0$ by holomorphic normal homogeneity. Its flux bracket has four normal antiholomorphic differentials. It is therefore strictly closed with all three shears present. In ordinary cohomology its normal vector component along $\partial_{w_s}$ is $-\alpha_s\partial_w^{\alpha-e_s}f$ whenever $\alpha_s>0$, where $f=(w_1w_2w_3)^{-1}$. A flux bracket on a regular scalar has only tangent components. The scalar input likewise survives because it has no incoming flux bracket.

With the vector first, its bracket with the lifted scalar vanishes at exterior degrees zero and one. All positive normal derivatives of a horizontal angular form remain horizontal for the antiholomorphic radial vector, so the degree-three terms vanish in its two-dimensional annihilator. At the pair $(i,j)$ the normal Euler vector acts on $r_h$ with weight $-1$. Hence
$$
C_{ij}=-r_h\partial_w^\alpha U,\qquad
w_1 C_{23}+w_2 C_{31}+w_3 C_{12}=-\partial_w^\alpha U .
$$
The negative of the preceding pair primitive now works on each restriction. The full three-shear class is nonzero for every $\alpha$ with $|\alpha|\ge1$ and the detected classes are linearly independent. For $\alpha=(0,0,1)$ the normal residue after multiplication by $w_3$ is $+1$ and the linear dictionary sends the detected $-\partial_{w_3}f$ to $dz_1dz_2\,d\theta_3$. The replay verifies all nineteen derivative multi-indices of orders one through three, including the full flux differential, divergence, shear action and bracket. — Codex. [Linear shears and signed full-complex calculation](code/three_shear_two_states.py)

**Finite-rank boundary test.** The linear dictionary sends $-\partial_w^\alpha f$ to $dz_1dz_2(d\theta)^\alpha$. Put $m=|\alpha|\ge1$ and take the $\theta=0$ coefficient of this form monomial in $\operatorname{Tr}(dC)^{m+2}$. Polarize the two first ghost jets and the three scalars along one traceless Cartan matrix $H$, writing $p_1,p_2$ for the odd jet coefficients and $x_I$ for the even scalar coefficients. In the source ordering the coefficient is
$$
g^2\frac{(m+2)!}{\alpha_1!\alpha_2!\alpha_3!}
\operatorname{Tr}(H^{m+2})\,p_1p_2\,x_1^{\alpha_1}x_2^{\alpha_2}x_3^{\alpha_3}.
$$
Indeed the ghost component of $C$ is $igc$ and moving $p_1$ across $dz_2$ cancels $i^2$. For $SU(2)$ take $H=\operatorname{diag}(1,-1)/2$. Its trace power is $2^{-m-1}$ for even $m$ and zero for odd $m$. At odd $m$ the entire cochain vanishes: the entries of $dC$ are even and commute, so the traceless $2\times2$ Cayley–Hamilton identity gives $(dC)^2=\operatorname{Tr}((dC)^2)\mathbf1/2$. For $SU(3)$ the test matrix $H=\operatorname{diag}(1,1,-2)$ gives $2+(-2)^{m+2}\ne0$ for every $m\ge1$. Each nonzero evaluation excludes a $Q_0$-boundary, since the classical variation and its jets vanish on commuting fields. The paper's single-graviton closure then excludes a quantum primitive over $\mathbb C[[\hbar]]$ by its constant coefficient. To test $\hbar$-torsion one must also evaluate the possible incoming quantum images. This comparison uses the linear single-particle dictionary. — Codex. [Superfield normalization](../../sources/LLC_main.tex#L4771) · [Single-graviton closure](../../sources/LLC_main.tex#L4107) · [Signed Laurent dictionary](../../code/holography/derham.py#L86)


The possible incoming differential can also be excluded. In the target the two derivative counts are both one. An incoming $d_1$ must come from the zero-derivative sector. Its letter count is $m+2$ and its total odd-coordinate degree is $m+3$. The relative complex there contains no undifferentiated ghost, so the source contains exactly one $\beta$ and $m+1$ scalars. Write it as $\sum_{I\in\mathrm{flavors}}\operatorname{Tr}(\beta_I F_I(X))$, with $X_I=\gamma^I$ and $F$ equivariant under conjugation. Its classical closure is $dW(F)=0$ for $W=\operatorname{Tr}X_1[X_2,X_3]$. Rotate the nonzero flavor vector $(x_1,x_2,x_3)$ to $(1,0,0)$ by $SL_3$, so the evaluation point becomes $(H,0,0)$ and $dW=0$ there. For a root $\rho$ with $\rho(H)\ne0$, torus equivariance gives flavor matrices $A_+,A_-$ for the derivative of $F$ on the opposite root lines. Infinitesimal conjugation of $(H,0,0)$ gives
$$
(A_+)_{11}=(A_-)_{11}=\frac{\rho(F_1)}{\rho(H)} .
$$
On opposite root vectors with invariant trace pairing one, the Hessian of $W$ on flavors $2,3$ is $\rho(H)J$, where $J_{23}=1$, $J_{32}=-1$ and the diagonal vanishes. Differentiate $dW(F)=0$ in these two root directions. Writing $A_{\pm,\perp}$ for the $2,3$ blocks gives
$$
\rho(H)(A_{+,\perp}^{T}J+JA_{-,\perp})+\rho(F_1)J=0,
\qquad \operatorname{tr} A_+ + \operatorname{tr} A_- =0 .
$$
The last identity follows by multiplying the first by $J^{-1}$, taking the trace and restoring the two flavor-$1$ entries. In the two-ghost projection of the published one-loop fan the two structure constants act by the adjoint action of $H$ twice and the mixed $\gamma^I,\beta_I$ derivatives trace over flavor. The contribution of each pair is therefore a common factor times $\rho(H)^2(\operatorname{tr}A_++\operatorname{tr}A_-)=0$. Zero roots and roots with $\rho(H)=0$ give zero directly. The chosen Cartan test annihilates $Q_1$ of every possible classical cocycle and all $Q_0$-boundaries. For $r\ge2$ the input to $d_r$ would need derivative counts $(1-r,1-r)$. Thus the nonzero $SU(2)$ and $SU(3)$ evaluations above certify survival after inverting $\hbar$ as well. In $U(N)$, choosing the central test $H=\mathbf1$ gives trace $N$ and kills all quantum images directly: the classical brackets vanish on central jets and every positive-loop fan has a central output color in a structure constant. — Codex. [Component one-loop fan](../../sources/LLC_main.tex#L4067) · [All-order weight shift](../../sources/LLC_main.tex#L443) · [Hessian and 679 polarizations](code/finite_shear_states.py)

##### Finite-rank count for invariant polynomials in the boundary superfield

For $SU(N)$ and degree $L$, take any homogeneous conjugation-invariant polynomial $P$ on traceless matrices and extract the $\theta=0$ form coefficients of $P(dC)$. The $\gamma,c$ sector above has ordinary multiplication and all its positive-loop differentials vanish on cochains. The single-trace factors are classically closed; therefore every coefficient of $P(dC)$ is closed. On the same Cartan fields, for a set $S\subseteq\{1,2\}$ of tangent differentials, $d=|S|$ and $|\alpha|=L-d$, the coefficient in increasing $dz$ order is
$$
(ig)^d(-1)^{d(d-1)/2}
\frac{L!}{\alpha_1!\alpha_2!\alpha_3!}\,
P(H)\left(\prod_{i\in S}^{\mathrm{increasing}}p_i\right)
x_1^{\alpha_1}x_2^{\alpha_2}x_3^{\alpha_3}.
$$
Indeed on these fields $dC$ is $H$ times a scalar even form. The two first ghost jets and the three scalar coefficients are independent. For $d=0,1$ at least one derivative count is zero and every incoming positive-loop differential would have a negative count. For $d=2$ the opposite-root cancellation just proved applies to every possible one-loop preimage, including multi-traces; higher incoming differentials have negative derivative counts. This argument includes $L=2,d=2$ with zero scalar spectators. Hence the Cartan evaluations distinguish these classes also over $\mathbb C((\hbar))$.

Restriction to a diagonal identifies invariant polynomials with symmetric polynomials in $N$ eigenvalues whose sum is zero. The dense set of diagonalizable matrices gives injectivity, and Newton's identities give free generators $\operatorname{Tr}H^j$, $2\le j\le N$. For fixed $d\le L$, the exact dimension of this family in quantum cohomology is consequently
$$
\binom{2}{d}\binom{L-d+2}{2}
[x^L]\prod_{j=2}^{N}(1-x^j)^{-1}.
$$
Here $[x^L]$ extracts the degree-$L$ coefficient. Summing the three possible values of $d$ gives $(2L^2+2L+1)[x^L]\prod_{j=2}^N(1-x^j)^{-1}$. This counts only the displayed $\theta=0$ invariant-polynomial family. For $SU(3)$ at $L=6$ the two color polynomials are $(\operatorname{Tr}H^2)^3$ and $(\operatorname{Tr}H^3)^2$. The diagonals $(1,-1,0)$ and $(1,1,-2)$ give the evaluation matrix with rows $(8,0)$ and $(216,36)$, of determinant $288$. Thus the $56$ zero-ghost, $84$ one-ghost and $30$ two-ghost components survive independently. — Codex. [Exact trace and multiplicity checks](code/finite_shear_states.py)

For the S-dual gauge algebras $\mathfrak{so}(2N+1)$ and $\mathfrak{sp}(N)$ the opposite-root argument uses the same invariant trace pairing on each pair of complex root lines. In the defining representation the Cartan eigenvalues occur in pairs $h_i,-h_i$, with one extra zero in the odd orthogonal case. The signed permutations act independently on the $h_i$. Symmetric polynomials in the $h_i^2$ are freely generated by their first $N$ power sums, which are half the even defining traces in either algebra. The entire protected family has the same degree-$L$ multiplicity on both sides,
$$
(2L^2+2L+1)[x^L]\prod_{j=1}^{N}(1-x^{2j})^{-1}.
$$
At $N=3,L=6$ the color polynomials are the degree-six trace, the product of the degree-four and degree-two traces, and the cube of the degree-two trace. Their $\theta=0$ components give $255$ independent classes in each theory, including $45$ with both differentiated ghosts. This excludes an incoming perturbative lifting pair for every state in this family even in the $SO(7)$ theory, whose other multi-graviton sectors have known incoming differentials. — Codex. [Scope of the paper's finite-rank lifting statement](../../sources/LLC_main.tex#L4087)

### Changing a source by a Laurent state

The five polynomial scalar sources of the third-order-pole example belong to an exactly compatible classical family. Write their sum with independent couplings as $A(w)+w_3H(z)$ and put $X_H=\partial_\Omega(H\,\xi_{z_1}\xi_{z_2})$. The strictly closed lift adds $-N b_3 X_H$, with $\bar\partial b_3=w_3U$. Two lifts for $H$ and $K$ also have zero classical bracket: their scalar–vector cross terms are proportional to $w_3 b_3[X_H(K)+X_K(H)]=0$, and the vector–vector terms have $b_3\wedge b_3=0$. The normal scalar $A$ differentiates neither tangent vector. The lifts are divergence-free, so all free descendant equations hold. With only primary fields the classical Hamiltonian has just its binary term. These sources can therefore be turned on together to arbitrary classical order at fixed flux.

The polynomial five-source coefficient $-N/1296$ is a nonzero class with the input states held fixed. It can be changed if source continuation at order $\lambda^4$ may introduce independent singular bulk states. In reduced units, put $u=w_2^4/4!$, one of these five sources, and write $\xi_j=\partial_j$ with the coordinate named in the subscript. The Laurent bivector
$$
Y=-\frac{N}{216}\partial_\Omega
\left(\frac{z_1\,\xi_{z_1}\xi_{z_2}\xi_{w_2}}{w_1w_2^4w_3^3}\right)
$$
is divergence-free and closed under the linearized flux differential. Its bracket with $u$ contracts the bivector with $w_2^3dw_2/6$, giving $+N\xi_{z_2}/(1296w_1w_2w_3^3)$. Polarize the five sources by independent even square-zero parameters. In the coordinate convention $Q\alpha_0=-\partial_\Omega E_0$, add $-Y$ times the product of the other four parameters at genus two. Its linear variation vanishes, and the classical coordinate variation $-\{u,-Y\}$ cancels the measured coefficient in this polarized family. This change injects an independent normal Laurent mode at fourth order in the sources. Its existence concerns the freedom in continuing states, beyond the local equivariant coordinate changes tested above.

For any normal scalar $u(w)$ with nonzero derivative in direction $w_j$, every tangent vector $\partial_\Omega(H(z)L\,\xi_{z_1}\xi_{z_2})$ with a finite Laurent class $L$ occurs in its bracket image on closed Laurent states. Choose $L'$ with $(\partial_{w_j}u)L'=L$ and take $Y=-\partial_\Omega(HL'\xi_{z_1}\xi_{z_2}\xi_{w_j})$. To see that $L'$ exists, pair Laurent classes with normal Taylor jets by residue. If $L$ measures jets of degree at most $l$ and $\partial_{w_j}u$ starts at degree $h$, multiplication by $\partial_{w_j}u$ injects jets modulo degree $l+1$ into jets modulo degree $l+h+1$. The lowest homogeneous components multiply in a polynomial domain. The dual map is surjective. This proves the bracket assertion without making a canonical choice of singular continuation. [Mode and sign checks](code/flux_state_bracket.py)

## Simultaneous tangent shears of a Laurent scalar

The three infinitesimal shears $v_i=w_i\partial_{z_1}$, $1\le i\le3$, commute, preserve the flux and have zero divergence. Write $U$ for the normal degree-two Dolbeault representative of the unit Čech class $f=(w_1w_2w_3)^{-1}$. On $H(z_1)U$ each shear acts by the exact form $w_i U H'$. Simultaneous equivariant completion means closure for the linearized flux differential plus $\sum_{i\in\mathrm{shears}}\epsilon_i{\cal L}_{v_i}$, with exterior parameters $\epsilon_i$ placed to the left. At nonzero $N$ such a completion exists **exactly when $H'''=0$**.

Set $r_i=\bar w_i/|w|^2$. For cyclic $(i,j,k)$ use $b_i=r_j\bar\partial r_k-r_k\bar\partial r_j$ and $U=2\sum_{\mathrm{cyclic}\ (i,j,k)}r_i\bar\partial r_j\wedge\bar\partial r_k$. The factor two is fixed by the Čech unit. In the cover $w_i\ne0$ put $A_i=b_i/w_i$ and $B_{ij}=-\varepsilon_{ijk}r_k/(w_iw_j)$ for distinct indices, where $\varepsilon_{123}=1$. Then $\bar\partial A=U$, $\delta A=\bar\partial B$, and $-\delta B=1/(w_1w_2w_3)$ with total differential $\delta+(-1)^p\bar\partial$ in Čech degree $p$. Also
$$
\bar\partial b_i=w_iU,\qquad
w_i b_j-w_j b_i=-\varepsilon_{ijk}\bar\partial r_k,\qquad
U\wedge b_i=0,\qquad
\partial_{w_i} U=-3r_iU .
$$
Here and below in each cyclic expression $k$ is the remaining normal coordinate. The first two identities follow from $\sum_i w_i r_i=1$ and its Dolbeault derivative. The last follows directly from the cubic denominator in $|w|^2$.

Use the partial lift
$$
H U+H'\sum_{i\in\mathrm{shears}}\epsilon_i b_i
+H''\sum_{\mathrm{cyclic}\ (i,j,k)}\epsilon_i\epsilon_j r_k .
$$
The entire residual consists of the following two terms in exterior degrees two and three:
$$
-\frac{N H'''}{3}\sum_{\mathrm{cyclic}\ (i,j,k)}
\epsilon_i\epsilon_j\,(\partial_{w_k} U)\,\xi_{z_2}
+H'''\epsilon_1\epsilon_2\epsilon_3 .
$$
The first term is the flux bracket on the last scalar correction, since $\partial_\Omega(HF)=H'U\xi_{z_2}$. The cubic follows from $\sum_i w_i r_i=1$. All displayed fields have zero divergence, so the $t\partial_\Omega$ terms vanish.

The quadratic residual is intrinsic modulo first-order choices. Reduce modulo $t$. The linearized differential raises Dolbeault degree minus polyvector rank by one, and the shears preserve it. A closed first-order ambiguity has this difference equal to one. Its scalar $(0,1)$ component is Dolbeault-exact because $H^1(X,\mathcal O)=0$. Subtract its linearized primitive. The next component is now a closed Laurent vector. Changing the scalar primitive adds the flux bracket of a holomorphic scalar, whose $\xi_{z_2}$ coefficient has only the simple class $f$. The Lie derivatives of normal vector polarizations also give terms along $\xi_{z_1}$. Thus in the $\xi_{z_2}$ projection, modulo multiples of $f$ by tangent functions, the first-order choices change the coefficient at $(i,j)$ only by $w_i\partial_{z_1}Y_j-w_j\partial_{z_1}Y_i$ for Laurent coefficients $Y_i$.

Multiply the three oriented pair coefficients by the complementary $w_k$ and add. This kills those differences and all simple $f$ ambiguities. In the Laurent quotient $w_k\partial_{w_k} f=-f$ for each normal coordinate, so the detector on the residual equals $NH'''f$. For $H=z_1^3/3!$ its value is $Nf$. Hence no choice of the first two homotopies removes the simultaneous obstruction, even in the full descendant complex. For $H'''=0$ the displayed partial lift already closes. Requiring equivariance for all six $w_i\partial_{z_a}$ with $a=1,2$ gives precisely the polynomials of total tangent degree at most two. Necessity follows by applying the transverse-vector detector in each tangent direction (a second tangent coordinate in $H$ also gives a longitudinal component, removed by this projection). For sufficiency replace each $H'$ by its derivative in the corresponding shear direction, and each $H''$ by the symmetric pair of tangent derivatives. All residual terms then contain three tangent derivatives. This is a classical higher action on the linearized bulk state complex. [Unit normalization, descent and detector](code/shear_descent.py)

For a general primary Laurent scalar with finite normal expansion, write it as $\sum_{\alpha\in\mathrm{normal\ derivatives}}H_\alpha(z)\partial_w^\alpha f$. In scalar rank zero at $t^0$ no flux bracket can be an incoming boundary. The six first-order equations require $w_i\partial_{z_a}$ of the Laurent class to vanish for each $i,a$. Since $w_i\partial_w^\alpha f=-\alpha_i\partial_w^{\alpha-e_i}f$ in that module, every $H_\alpha$ with $|\alpha|>0$ must be constant. All such constant-coefficient normal derivatives are strictly invariant and closed. Together with the six tangential polynomials of degrees zero through two multiplying $f$, they exhaust the jointly equivariant primary Laurent scalars.

There is a graph selection rule for the rank-one Laurent projection of this obstruction. Take one scalar input, $p$ external shear vectors, $k$ flux bivectors, and a testing leg of rank $4-r$ for a Hamiltonian output of rank $r$ after its divergence. On the chosen pure-primary action every interaction vertex takes total polyvector rank five and every internal primary BV line carries rank four across its ends. This holds also for a puncture Dolbeault contraction preserving polyvector rank before the single BV divergence. With $M=p+k+2$ total external slots one has
$$
p+2k+4-r=5V-4E=M+2-2n-\sum_{v\in\mathrm{vertices}}s_v,
\qquad
r=2n+k+\sum_{v\in\mathrm{vertices}}s_v .
$$
At genus zero the measured rank-one term uses one flux. No positive-genus graph in this primary trace and propagator family can directly correct it. The explicit first- and second-order source lifts in the calculation are scalar and obey the same rank count. At $t^0$ the obstruction has Dolbeault degree minus polyvector rank equal to one; its only possible Laurent cohomology component is a rank-one $(0,2)$-vector. This explains why the nonzero secondary action is visible already in the classical punctured theory. The separate five-scalar loop measurement has four more scalar inputs and obeys a different rank count. [Primary fields and kernel ranks](../iib-bcov/operator/README.md)

For the componentwise Dolbeault homotopy there is a second bound, valid for any number of the six shears. Represent the polyvector slots at a vertex by five odd generators. Each ordinary BV edge carries a degree-four polynomial in the differences of the generators at its endpoints, since one holomorphic divergence acts on the degree-five diagonal polyvector kernel. A puncture kernel on the Dolbeault coefficients has the same property when it acts identically on the polyvector components. Write the three normal odd generators at each vertex as their values at a fixed root plus differences from that root. The product of the vertex traces must contain all three common normal generators. Edge factors contain none of them. The initial Laurent scalar and its scalar Čech lifts, every tangent shear, and every flux bivector also contain none. For the chosen bulk vertices the momentum polynomial multiplies the wedge trace. Thus only the testing leg can supply the three common normal generators, including when the graph is a single local vertex. It follows that
$$
3\le 4-r,\qquad r=2n+k+\sum_{v\in\mathrm{vertices}}s_v,\qquad n=0 .
$$
If a scalar Čech lift has absorbed some shears already, it still supplies a single rank-zero input and the same calculation applies with the remaining shears. At $t^0$ all positive-genus connected primary Hamiltonian graphs of this family therefore vanish on the shear action on one Laurent scalar. This covers the two-shear detector and mixed directions among the six shears, with any number of flux insertions. The five-scalar flux contact instead has $r=2n+k-4$ at minimum valence, allowing $r=1$ at genus two. The proof uses the wedge-trace action and componentwise polyvector kernel; marked composite operators or additional defect vertices require their own contractions.

For a general primary graph on the D3 plane, let $s$ count scalar inputs, $p$ tangent vectors and $k$ tangent bivectors, including the explicit flux insertions. Their coefficient functions may depend on all five coordinates. For the cubic and integrated primary trace vertices of valence $2g_v+3$, the same calculation gives
$$
r=1+2n+k-s,\qquad 0\le r\le1.
$$
Hence the rank-one response requires $s=2n+k$, and the scalar response requires $s=2n+k+1$. In particular a response linear in an explicit flux first admits a vector with five scalars at two loops and with seven scalars at three loops. The calculated local contacts attain both lower bounds. At rank one the testing leg has exactly three polyvector slots, all normal. Its complementary Euler derivative is therefore a tangent bivector, and its divergence is a tangent Hamiltonian vector. If a regular scalar source contributes its tangent-vector lift, count that leg as a vector and include the lift's external power of $N$. For $s_0$ such original scalar sources the total flux power in the rank-one response is exactly $N^{s_0-2n}$. These are selection rules for the specified primary vertices and propagator, independent of the Schwinger integrals and the number of external shear vectors.

## A finite normal trace at the puncture

Use alternating Čech cochains on the three opens $w_i\ne0$. A monomial $w^\alpha$ occurs on an intersection indexed by $I$ exactly when $I$ contains $S=\{i:\alpha_i<0\}$. In positive cochain degree set
$$
(hc)_I=\frac{1}{3-|S|}\sum_{\substack{j\notin S\\\mathrm{normal\ charts}}}c_{jI}
$$
when $|S|<3$, with a repeated index giving zero. Put $h=0$ on degree zero and on the fully negative sector. Write $\iota\pi$ for the projection onto the diagonal average of the three vertex values for regular monomials and onto the degree-two cochain for fully negative monomials. It vanishes on the other components. Inserting each $j$ in the alternating coboundary gives $\delta h+h\delta=1-\iota\pi$. Antisymmetry gives $h^2=0$. Each nonzero holomorphic derivative preserves $S$, so the homotopy commutes with coefficient derivatives and extends componentwise to polyvectors with the total-complex signs.

Pair monomials with exponents $\alpha,\beta$ when $\alpha+\beta=(-1,-1,-1)$ by integrating their Whitney forms on the oriented normal two-simplex. This is the linking residue pairing. A vertex pairs with the top face with weight $1/3$. Two distinct oriented edges pair with weight $\pm1/6$. For a degree-$p$ cochain $c$ the resulting identity is $\langle hc,d\rangle=(-1)^p\langle c,hd\rangle$. The projection is self-adjoint. The eight sign sectors exhaust the proof by the finite matrices in the replay.

For three normal polynomials $g_1,g_2,g_3$, write $M_{g_i}$ for ordinary multiplication on these cochains. For a finite normal Laurent class $f$ in degree two, $M_f$ is the Whitney product with its top cochain. It multiplies the average of vertex values by $f$ and kills positive cochain degrees. The supertrace, the alternating trace in Čech degree, is finite and satisfies
$$
\operatorname{str}\!\left([\iota\pi,M_{g_1}][h,M_{g_2}][h,M_{g_3}]M_f\right)
=-\frac16\operatorname{Res}_{w=0}\!\left(f\,dg_1\wedge dg_2\wedge dg_3\right).
$$
The residue extracts the coefficient of $(w_1w_2w_3)^{-1}dw_1\wedge dw_2\wedge dw_3$. For monomials $g_j=w^{u_j}$, a commutator changes the sign sector only when at least one exponent crosses from negative to nonnegative. Three nonzero commutators therefore require precisely one distinct crossing apiece. This bounds all three incoming exponents. Only the term of $f$ of exponent $-u_1-u_2-u_3$ can contribute to the diagonal. If the $j$th multiplication crosses coordinate $\sigma(j)$, the number of incoming monomials is $\prod_j (u_j)_{\sigma(j)}$. For each, the top cochain successively contracts to the opposite edge, to a vertex with factor $1/2$, and to the regular diagonal with factor $1/3$. The remaining oriented coefficient is $-\operatorname{sgn}(\sigma)$. Whitney multiplication by the top reads off the vertex average, which pairs to one on the regular diagonal. Thus the sum is $-\det(u_1,u_2,u_3)/6$ times that coefficient of $f$. Linearity gives the displayed trace. For the unit D3 Laurent class $1/(w_1w_2w_3)$ it measures only the normal Jacobian at the brane. For $g_i=w_i$ its value is $-1/6$. The tangential and polyvector contractions, physical graph weights and loop integration are separate inputs to the bulk amplitude. [Exact contraction and 614 monomial traces](code/cech_normal.py)

There is an independent evaluation at the cutoff. Let $\Pi_L$ retain all Laurent exponents in $[-L,L]^3$ and let $\operatorname{Alt}$ be the signed sum of the six permutations of the $g_i$, with no factorial. For every box containing the support of the finite diagonal trace,
$$
\operatorname{Alt}\operatorname{str}\!\left(
 \iota\pi M_{g_2}h M_{g_3}h M_f [M_{g_1},\Pi_L]\right)
=-\operatorname{Res}_{w=0}\!\left(f\,dg_1\wedge dg_2\wedge dg_3\right).
$$
Insert $\Pi_L$ to the left of the preceding trace and alternate. Expand its three commutators. Adjacent regular multiplications cancel under alternation, and $h^2=\iota\pi h=0$. Of the two surviving terms, rotate one using cyclicity of the trace with the finite-rank cutoff. The regular multiplication commutes with $M_f$; its commutator with $\Pi_L$ is exactly the displayed outgoing boundary term. For $g_i=w_i$ and flux coefficient $N/(w_1w_2w_3)$ its value is $-N$ at every $L\ge1$. This recovers the normal three-dimensional residue cocycle on one Laurent and three regular arguments, with its coefficient fixed by the ordered contractions. Identifying its coefficient in an IIB observable requires the physical vertex and field traces. [Boundary-mode replay](code/cech_normal.py) · [Faonte–Hennion–Kapranov, equation (2.3.1)](https://arxiv.org/pdf/1701.01368#page=35)

The alternated residue is independent of how the available charts are weighted in $h$, provided the weights in each sign sector sum to one and the regular projection uses the degree-zero weights. Indeed fix the coordinates crossed in order from left to right. The last projection supplies the regular weight of the first coordinate, the intervening homotopy the weight of the second among the remaining two, and the initial homotopy has one available chart. Alternating the three monomials gives the signed exponent determinant times these weights. Summing over the orders of the coordinates multiplies the determinant by one. The replay checks each of the six choices of a strict chart priority separately.

The same contraction evaluates three commuting shears $w_i^{m_i}\partial_{z_1}$ for arbitrary positive integers $m_i$. A source $H(z_1)/\prod_i w_i^{b_i}$ with $1\le b_i\le m_i$ admits individual primitives on the three complementary edges. After two shears the oriented scalar lift at the complementary vertex $k$ has coefficient $H''w_i^{m_i}w_j^{m_j}/\prod_i w_i^{b_i}$. The third residual is the diagonal regular scalar
$$
H'''(z_1)\prod_{i=1}^3 w_i^{m_i-b_i}
\pmod{(w_1^{m_1},w_2^{m_2},w_3^{m_3})}.
$$
In the primary scalar projection first-order changes of primitives are exact since $H^1=0$. Second-order changes then contribute only elements of the displayed ideal. The $m_1m_2m_3$ monomials form its quotient basis, so different Laurent inputs have independent cubic residuals. This scalar projection applies at $t=0$ for every flux $N$, since the flux bracket raises polyvector rank.

For sufficiency in the smooth flux complex, take $u_i=w_i^{m_i}$ and pull back $r_i,b_i,U$ of the shear descent along $w\mapsto u$. A unit Čech representative of the D3 flux coefficient is $\prod_i w_i^{m_i-1} U(u)$, and the source uses $\prod_i w_i^{m_i-b_i} U(u)$. Their alternating Čech classes are respectively $1/\prod_i w_i$ and $1/\prod_i w_i^{b_i}$. The pulled-back forms obey $U(u)\wedge b_i(u)=0$. The first and second scalar primitives of the earlier descent now use $b_i(u)$ and $r_k(u)$ with the same source prefactor. Every residual contains $H'''$. Changing between these flux representatives adds a Dolbeault-exact coefficient times the tangent bivector. Its primitive times the same bivector has vanishing mutual bracket with both fluxes, commutes with these shears and has zero divergence, so the corresponding linearized complexes are conjugate. Thus simultaneous completion holds precisely for $H'''=0$. For all six shears with either tangent direction, polarization gives the same criterion with all third tangent derivatives. On a finite Laurent sum, the first-order equations already force every tangentially varying coefficient to lie in the indicated box of pole orders. [Exact shear lifts on 216 inputs](code/cech_normal.py)

## Shear Ward test for a one-loop trilinear operation on scalar states

Consider the part linear in $N$ of a trilinear operation on polynomial regular scalar classes at fixed nonzero flux. Project its primary output at $t^0$ to the Laurent tangent vector. For the componentwise propagators above the Euler derivative is a tangent bivector, so the output after divergence is a tangent Hamiltonian vector. Its part with coefficient $f$ is the boundary $\{NF,H(z)\}$ after dividing by $N$. Keep the parts with at least one normal derivative of $f$.

Under a block-unitary linear coordinate change $g\in U(2)\times U(3)$ each BV line has character $(\det g)^{-2}$, one inverse determinant from the antiholomorphic kernel and one from the polyvector volume. The trace on pulled-back fields has the same character. A graph with equally many lines and vertices therefore defines an invariant integrated functional. Its Euler derivative evaluated on pulled-back inputs equals $(\det g)^2$ times the pullback of the Euler derivative. The single background transforms as $g^*F=(\det g)^{-1}F$. If the three scalar inputs have homogeneous normal degrees $n_i$, tangent degrees $d_i$, and the output potential is a tangent polynomial of degree $m$ multiplying $\partial_w^\alpha F$, the two circle weights imply
$$
|\alpha|+n_1+n_2+n_3=6,\qquad m+4=d_1+d_2+d_3.
$$
Polarize three symmetric normal input jets by covectors $p_1,p_2,p_3$. Write normal derivatives of the flux as a polynomial in a fourth variable $p_4$, so multiplication by $w_j$ acts on it as $-\partial_{(p_4)_j}$. The $SL_3$ invariant tensors in six normal slots are products of two alternating triples. Thus each of $n_1,n_2,n_3,|\alpha|$ is at most two.

For the shear $v=w_j\partial_{z_a}$ the ordinary linear action preserves the regular scalar classes. Indeed, for a scalar $A$ of positive normal degree $n$, put $A_i=(\partial_{w_i}A)/n$ so $A=w_iA_i$ with the normal Euler contraction. Its closed lift is $\widehat A=A-N\sum_i b_iX_{A_i}$, where $X_G$ is the tangent Hamiltonian vector of $G$. The difference between $L_v\widehat A$ and $\widehat{vA}$ is
$$
\frac{N}{n+1}\sum_{i\in\mathrm{normal\ coordinates}}(w_i b_j-w_j b_i)X_{\partial_{z_a}A_i}.
$$
The parenthesis is a Dolbeault derivative of $-c\epsilon_{ijk}r_k/2$ in the angular convention above. The corresponding primitive times the tangent Hamiltonian vector has vanishing divergence and flux bracket, so the difference is exact for the full linearized differential.

Suppose the trilinear map commutes on states with these linear shears and tangent translations. On tangent exponentials let $k_i$ be the three momenta and $k_4=k_1+k_2+k_3$. On polynomial inputs translation and the tangent degree imply that the vector symbol is polynomial of momentum degree five: its value at the origin uses only tangent jets of total degree five. The divergence forces its two components to be the transverse tangent vector of $k_4$ times a scalar polynomial of degree four, since the two components of the generic $k_4$ are relatively prime.

If all three inputs vanish on the brane, all four normal degrees in the projected symbol are positive. Their only partition is $(2,2,1,1)$. With $[p_i,p_j,p_k]$ the normal volume determinant, the six possible normal tensors have doubled pair $i,j$ and take the form
$$
[p_i,p_j,p_k]\,[p_i,p_j,p_l],\qquad \{i,j,k,l\}=\{1,2,3,4\}.
$$
A shear differentiates all four normal variables in one common direction, weighted by the corresponding $k_i(\partial_{z_a})$. Differentiating a singleton drops out: the regular scalar slots vanish on the brane and the simple-flux output is a boundary. Fix the one label $j$ that stays doubled. The derivatives of the three possible other doubled labels have exactly one relation, their sum zero. In the normal plane modulo $p_j$, write $\omega(u,v)=[p_j,u,v]$. For the other three labels write $A=\omega(v,p_1)\omega(p_2,p_3)$ and $B,C$ for its two cyclic permutations, temporarily numbering those labels $1,2,3$. The plane identity gives $A+B+C=0$. The three derivatives are $C-B,A-C,B-A$ and span a plane. The three coefficients of the original tensors, multiplied by their respective tangent covectors $k_i$, must therefore agree. For generic $k_1,k_2,k_3$ these covectors lie on different lines, forcing all six coefficients to vanish. Polynomiality extends this to all tangent momenta. If an input has zero normal degree it is a constant regular class. The only potential normal tensor is then $[p_i,p_j,p_4]^2$. Its three shear derivatives have distinct multidegrees and again force zero. [Independent exact rank calculation](code/pentagon_shear_ward.py)

Thus a strict shear-equivariant trilinear state map with the displayed weights has zero Laurent tangent-vector response. The isolated puncture line fails the required Ward identity. For a tangent wave of momentum $k$, its commutator with $w_i\partial_{z_a}$ on $e^{kz}w_jU$ is $k_a e^{kz}w_i b_j$: $h(w_jU)=b_j$, whereas $h(w_iw_jU)=0$ by normal phases. The cyclic transpose pairs to zero on both by radial horizontality. Consequently a physical application of the vanishing requires the full punctured transfer, or a Ward proof that its nonlinear transformations on the external states give no terms at this genus and arity. The value of an individual separated kernel tests a different functional.

The scalar inclusion itself has a nonzero second shear homotopy on regular modes. Use $v_i=w_i\partial_{z_1}$, put $A=w_sH(z_1)$, and write $Y_j(A)$ for the preceding primitive of $L_{v_j}\widehat A-\widehat{v_j A}$. The same formula on the quadratic normal source $v_iA$ gives
$$
Y_j(A)=-\frac{Nc}{4}\epsilon_{sjk}r_k X_{H'},\qquad
Y_j(v_i A)=-\frac{Nc}{12}(\epsilon_{ijk}w_s+\epsilon_{sjk}w_i)r_k X_{H''}.
$$
Since $\sum_{k\in\mathrm{normal\ coordinates}} w_kr_k=1$, the second residual is
$$
L_{v_i}Y_j(A)-L_{v_j}Y_i(A)-Y_j(v_iA)+Y_i(v_jA)
=\frac{Nc}{6}\epsilon_{sij} X_{H''}.
$$
For $H=z_1^3/6$ the output is the regular constant vector $\partial_{z_2}$. To test choices, alternate over $s,i,j$, evaluate the regular output at $w=0$, and select its $\partial_{z_2}$ component at $t^0$. A closed change of a first primitive has Dolbeault-degree-minus-polyvector-rank $-1$, so its lowest component is a holomorphic vector; it extends across the missing plane. Acting with $v_i$ on its coefficient vanishes there, and the action on its normal-vector polarization points along $\partial_{z_1}$. Changes evaluated on $v_i A=w_iw_sH'$ cancel by antisymmetry. A linearized boundary cannot give a regular vector at form degree zero. Thus the alternating residual is independent of all such changes. This supplies a classical higher transformation on the regular source sector; the physical one-loop Ward equation still requires its appropriate contact and transfer terms.

The published finite twisted superisometry is $\mathfrak{psl}(3|3)$, with even part $\mathfrak{sl}_3\oplus\mathfrak{sl}_3$. A normal $\mathfrak{sl}_3$ with nontrivial commuting tangent $\mathfrak{sl}_2$ lies in one simple factor: a diagonal embedding has zero centralizer. The even adjoint therefore contains no representation simultaneously a tangent doublet and a normal fundamental or dual. The six coordinate shears have precisely these two labels and ordinary even flow parameters. A Ward identity for them requires the larger bulk coordinate symmetry with its boundary conditions and transfer, beyond the finite physical superisometry. [Twisted isometry, section 9](https://arxiv.org/abs/1606.00365)

### The two physical special conformal generators

In the repository coordinates the two published special conformal vector fields are $K_a=z_a(E_z-E_w)$, $a=1,2$, with $E_z,E_w$ the tangent and normal Euler vector fields. Write $P=\xi_{z_1}\wedge\xi_{z_2}$, so $F=UP$. Since $E_wU=-3U$ and $[K_a,P]=-3z_aP+E_w\wedge X_{z_a}$, one obtains
$$
[K_a,F]=U E_w\wedge X_{z_a},\qquad
K_a^N=K_a-N\sum_{i\in\mathrm{normal\ coordinates}} b_i\,\xi_{w_i}\wedge X_{z_a}.
$$
Here $[K_a,F]$ is the ordinary Lie derivative. The project's suspended BV bracket gives $\{F,K_a\}=[K_a,F]$. The added term has Dolbeault derivative $-N[K_a,F]$ and zero divergence by $\partial_{w_i}b_i=-2r_i b_i$ and $\sum_i r_i b_i=0$. Its bracket with $F$ vanishes because $X_{z_a}\wedge P=0$ and $b_i$ is independent of $z$. Thus $K_a^N$ is closed under the flux differential. [Costello–Li, §9.2, Proposition 9.2.2](https://arxiv.org/html/1606.00365#S9.SS2) · [Project BV sign](../../code/holography/bulk_x.py#L128)

The two lifted special conformal generators commute in cohomology. Their bracket is $2N\sum_{\mathrm{cyclic}\ (i,j,k)}\bar\partial r_k\,\xi_{w_i}\wedge\xi_{w_j}$, the full differential of the same sum with $r_k$ in place of $\bar\partial r_k$. The divergence and the flux action on this primitive each contract the alternating normal tensor with two copies of $r$. [Exact commutator](code/projective_flux.py)

**The seven flux-square brackets on the finite symmetry algebra.** The projective vector $z_b\partial_{\theta_I}$ in the supercoordinate dictionary corresponds to $-z_b w_I$ in the regular bulk fields. Its smooth lift is $-z_b w_I+N b_I X_{z_b}$. Take the suspended BV bracket with $K_a^N$ first and use $\epsilon_{12}=1$ for the tangent coordinates. For $a\ne b$ the result is
$$
\{K_a^N,-z_b w_I+N b_I X_{z_b}\}
=\bar\partial\left(\epsilon_{ab}N
\sum_{\substack{i,k=1,2,3\\\mathrm{normal\ coordinates}}}
\epsilon_{iIk}r_k\,\xi_{w_i}\right).
$$
The six primitives have zero divergence and zero flux bracket: $\partial_{w_i}r_k=-r_i r_k$ and $\partial_{w_i}U=-3r_iU$ each contract a symmetric product with $\epsilon_{iIk}$. The same contractions annihilate the preceding two-normal-vector primitive for $\{K_1^N,K_2^N\}$. Thus all seven brackets are full linearized boundaries and, on these smooth lifts, their coefficients of $N^2$ already vanish as cochains. This gives zero on each weight-allowed pair for the transferred $\ell_4(F,F,-,-)$ in $\operatorname{coker}\varphi_*$ at genus zero. [Seven signed polynomial identities](code/projective_flux.py) · [Seven target weights](symmetry/S1.log)

There is also a representation-theoretic check on the transfer. Insertion of a normal $\mathfrak{sl}_3$ generator in this Laurent-valued two-cocycle vanishes by tangent weights: the sum of the two input weights must be positive in both tangent coordinates. The Chevalley–Eilenberg cocycle identity with this insertion then implies normal $\mathfrak{sl}_3$ equivariance. Write $V$ for the three coordinates $\theta_I$ with their normal action; $d\theta_I$ span the same representation and commute. At zero tangent label the Laurent cokernel in these weights consists of closed forms of degree at least two. The pair $(K_1,K_2)$ is a normal singlet, whereas the target at total normal degree three lies in $(V\otimes\operatorname{Sym}^2 V)\oplus\operatorname{Sym}^3 V$, with no singlet. For each mixed triplet the source is $V^*$ and the target lies in $\operatorname{Sym}^2 V$, which contains no copy of $V^*$. Hence even the cocycle in the weight-preserving cokernel vanishes, for any equivariant transfer realizing the project's classical bracket. The computation here fixes the binary symmetry extension; quantum operations on additional external states require their own transfer. [Symmetry and Laurent conventions](../../code/holography/symmetry.py) · [Costello–Li, Proposition 9.2.2](https://arxiv.org/abs/1606.00365)

Work modulo $t$ and $N^2$ to extract the primary term linear in flux, assuming the operation descends modulo $t$ and is regular in $N$. For each lifted polynomial scalar $\widehat A$ above, $[K_a^N,\widehat A]-\widehat{K_a A}$ is $N$ times a closed $(0,1)$-vector. The tangent polynomial $K_a A$ has the same normal degree as $A$. Since $H^1(X,\mathcal O)=0$, a vector primitive removes the discrepancy modulo $N^2$. At zero flux the three-source graph count gives a scalar of form degree one. For a closed output, subtract the flux differential of a scalar primitive of this zero-flux part. Its contribution at order $N$ is tangent, so the resulting representative is $N$ times a closed tangent vector of form degree two.

Let its Laurent coefficient be $V=V^b\xi_{z_b}$, with $b$ a tangent component index. The **normal** coefficient of $[K_a^N,NV]$ is $N V^a E_w$. Strict equivariance on three states forces this coefficient to vanish in Laurent cohomology: each transformed scalar input again gives a tangent output. At this bidegree every boundary has normal coefficient Dolbeault-exact since the flux bracket of a scalar has only tangent components. In the Laurent module
$$
w_j\partial_w^\alpha f=-\alpha_j\partial_w^{\alpha-e_j}f
\quad(\alpha_j>0),\qquad w_j\partial_w^\alpha f=0\quad(\alpha_j=0).
$$
The joint kernel of the three multiplications is therefore exactly $\mathbb C f$ in the normal variables. Applying the two $K_a$ leaves only a simple $f$ in each component of $V$. Its tangent Hamiltonian part $Nf X_H$ is the flux boundary of $H(z)$. Hence a primary trilinear map in this graph sector satisfying the two **strict physical** Ward identities has zero tangent-vector class to first order in $N$. A nonzero completed response would require an additional term in those identities or a failure of the specified primary projection or graph assumptions. The published Lie action on classical cohomology supplies the two generators; the quantum identities on three transferred sources remain to be computed.

This gives an explicit special conformal test on the five polynomial inputs of the genus-two vertex above. Its reduced one-flux output has $V^2=-1/(1296w_1w_2w_3^3)$ and $V^1=0$. Acting with $K_2^N$ and retaining the normal Laurent component gives
$$
-\frac{N\,\xi_{w_3}}{1296w_1w_2w_3^2}.
$$
The other two components vanish by Laurent multiplication. Acting on any of the five scalar inputs instead gives a tangent output for this isolated shifted vertex. The unshifted seven-field vertex has no five-input operation. Thus the displayed normal vector is a nonzero Ward discrepancy of the reduced local quintic on these closed inputs. Any strict special conformal identity for the full five-source response must cancel this precise Laurent coefficient. It can be detected by taking the normal residue of $w_3$ times the $\xi_{w_3}$ coefficient.

### A special conformal transformation mixes scalar and Laurent states at fixed flux

For $A=w_s H(z_2)$, let $D$ be the full linearized flux differential. The exact smooth action of $K_1^N$ in the angular normalization $U=2\sum_{\rm cyclic} r_i\bar\partial r_j\wedge\bar\partial r_k$ above is
$$
\{K_1^N,\widehat A\}-\widehat{K_1 A}
=D\!\left(-N H'\sum_{i,k\in\mathrm{normal\ coordinates}}
\epsilon_{isk}r_k\xi_{w_i}\right)
-\frac{N^2}{2} H''\sum_{i,k\in\mathrm{normal\ coordinates}}
\epsilon_{isk}\,\partial_{w_k} U\,\xi_{w_i}\wedge\xi_{z_1}.
$$
Indeed $K_1 A=w_s z_1(z_2H'-H)$. The order-$N$ normal component is $NH'\sum_i(w_i b_s-w_s b_i)\xi_{w_i}$, the Dolbeault derivative of the displayed primitive. At order $N^2$, the bracket of the two corrections gives $H'' b_i\wedge b_s$ in polarization $\xi_{w_i}\wedge\xi_{z_1}$. Subtracting the flux bracket of the primitive adds $H''\epsilon_{isk}r_k U$. The identities $b_i\wedge b_s=\epsilon_{isk}r_k U/2$ and $\partial_{w_k}U=-3r_k U$ give the formula. The pure tangent component vanishes by $\sum_i r_i b_i=0$ and antisymmetry. Both the primitive and the residual have zero divergence. For a primary boundary of rank two and form degree two, its regular vector primitive is holomorphic. The mixed normal–tangent part of its flux bracket has at most a simple $f$; differentiating $U$ only occurs in the pure tangent bivector part. Hence the displayed normal derivative of $f$ survives in cohomology whenever $H''\ne0$. [Signed exterior-algebra replay](code/projective_flux.py)

This mixing obstructs a conformally invariant linear choice of lifts of the entire normal-linear polynomial scalar sector. Take $H=z_2^2/2$. For tangent degree $d$ the physical dilation $E_z-2E_w/3$ acts on $w_s H$ with weight $d-2/3$. A mixed bivector with a first normal derivative of $f$ must then have tangent coefficient degree $d-3$, and a pure tangent bivector with two normal derivatives of $f$ must have degree $d-2$. Only these two types of Laurent changes of the input lift can reach the measured mixed coefficient under $K_1$. At $d=2$ the first degree is negative, and the second would require an $SL_2$ map from tangent homogeneous quadratics to the scalar tangent volume. There is no such map. A change of lift of $K_1 A$, of degree three, would require an $SL_2$ map from homogeneous cubics to a constant tangent vector, also zero. The spaces of binary forms of degrees $0,1,2,3$ are irreducible tangent $SL_2$ modules of dimensions $1,2,3,4$, respectively. Thus the mixed class persists for every lift compatible with tangent rotations and physical dilations. For $s=3$ its $\xi_{w_1}\wedge\xi_{z_1}$ component is $-N^2/(2w_1w_2^2w_3)$. The physical generators $w_s$ and $w_s z_a$ in this scalar family have vanishing second derivative. The extension detected by the quadratic source belongs to the classical flux complex.

For $s=3$ and $H=z_2^2/2$, the linear single-particle dictionary sends the full residual to the exact two-form
$$
\frac{N^2}{2}dz_2(\theta_1 d\theta_2-\theta_2 d\theta_1)
=\frac{N^2}{2}d(dz_2\,\theta_1\theta_2).
$$
The order-$N^2$ mixing therefore reaches a definite single-graviton channel with a fixed coefficient in the unit Čech convention. [Linear dictionary](../../code/holography/derham.py#L86)

## Ordinary boundary products in four component sectors

For arbitrary gauge-invariant local operators $A,B$ constructed from $c,\gamma^1,\gamma^2,\gamma^3$ and their holomorphic jets, the paper's canonical marked product and supercharge satisfy at every rank
$$
A\star B=AB,\qquad Q_n A=0\quad(n\ge1).
$$
Use only operators independent of the undifferentiated ghost in the relative complex. The same identities hold on each of the three sets $\{c,\beta_i,\gamma^j:j\ne i\}$, with the flavor index $i$ fixed. Classical gauge variation preserves these algebras; the superpotential variation of $\beta_i$ uses the two remaining scalars. In particular any number of ordinary products of closed operators from one such sector are already quantum cocycles. The first sector also follows from the paper's bound by the number of $b,\beta$ letters. [Cubic marked forests](../../sources/LLC_main.tex#L1049) · [Original termination bound](../../sources/LLC_main.tex#L3052)

In the component functional derivative each $c$ slot supplies $\theta_1\theta_2\theta_3$, each $\gamma^i$ the two complementary odd coordinates, each $\beta_i$ the single $\theta_i$, and $b$ a constant. At positive loop order a fan differentiates at least two input letters at the same odd coordinate. Removing both marked insertions from a marked product graph leaves trees meeting both insertions. The odd delta functions on each internal tree identify all its odd coordinates, so its two attachments already multiply two derivative factors. In each listed sector that product vanishes identically, before the color contraction or any Schwinger integral. Arbitrary spacetime derivatives preserve the odd factors. This argument also gives $A\star B=AB$ for two different sets whenever every derivative factor from one annihilates every derivative factor from the other. [Component derivative and propagator](../../sources/LLC_main.tex#L4792)

These are all maximal coordinate-letter sets with this pairwise vanishing property. More generally, classify flavor-graded linear spans by their derivative polynomials in the exterior algebra on the three odd directions. A nonzero constant term precludes universal pairwise vanishing. With no linear terms the maximal space consists of all quadratic and cubic terms. Otherwise the linear terms must span a single line generated by a form $\eta$, since independent one-forms have nonzero wedge. The quadratics annihilating $\eta$ are exactly those divisible by $\eta$, and all cubic terms are allowed. These maximal spaces are indexed by the projective plane of choices of $\eta$. Fixing the coordinate basis recovers the three sectors containing a $\beta$.

For a general monomial let $m$ count all $b$ letters, including their jets. The largest number of non-$b$ factors with pairwise disjoint derivative supports is three if all three distinct $\beta$ species occur. Otherwise it is two if two distinct $\beta$ species occur or some $\beta_i$ and $\gamma^i$ both occur, one if any non-$b$ letter remains, and zero otherwise. Write this number as $r$. A genus-$n$ fan is then zero whenever $n+1>m+r$. This cutoff is sharp for the odd-coordinate integral; kinematic and color factors enter afterwards. [Direct exterior algebra on unreduced fan and tree coordinates](code/boundary_odd_sectors.py)

In the D3 dictionary the two inputs of the [four-input tree](#a-four-input-tree-with-two-normalizable-single-graviton-states) have forms $d(z_1^3/6)d\theta_3$ and $d(z_2^2/2)d\theta_3$. Their boundary representatives and the output component $dz_2(d\theta_3)^3$ belong to the first sector: at $\theta=0$ only the ghost and scalar fields occur in $dC$. For the output, the coefficient of $dz_2(d\theta_3)^3$ in $\operatorname{Tr}(dC)^4$ is $4ig\operatorname{Tr}((\partial_{z_2}c)(\gamma^3)^3)$ in the paper's superfield normalization. For $SU(2)$ setting both colors to $\sigma_3/2$ gives $\operatorname{Tr}((\sigma_3/2)^4)=1/8$, so this graviton survives restriction to commuting colors. The mixed bulk tree probes the response to shear sources. The undeformed boundary product of its two graviton inputs has zero positive-loop correction. [Boundary superfield](../../sources/LLC_main.tex#L4772)

## The three-trace joining symbol and the four-point bulk tree

In [Request 5](../../AGENT_DISCUSSION.md#L122), take the three odd factors in cyclic order, $\operatorname{tr}(\beta_1\gamma^2)$, $\operatorname{tr}(\beta_2\gamma^3)$, $\operatorname{tr}(\beta_3\gamma^1)$. For cyclic $(i,j,k)$, the coefficient of $\theta_k$ in $\operatorname{tr}((\partial_{\theta_j}C)^2)$ is $-2\operatorname{tr}(\beta_i\gamma^j)$. They are single-graviton components. Their classical variations vanish by $\operatorname{tr}([\gamma^j,\gamma^k]\gamma^j)=0$, and each restricts nontrivially to the Cartan. Denote their ordered product by $\mathcal O_0$. In the ternary symbol $\Phi^3_{Q_2}$, every term on a proper subproduct vanishes because three $\beta$ derivatives are required. Thus the measured output is $Q_2\mathcal O_0$.

Use $SU(2)$ with $t_A=\sigma_A/2$ and specialize the three spectator matrices to $t_3$. Write $u^A_{ab}=\partial_{+\dot a}\partial_{+\dot b}c^A$ and $p_+^A=\partial_{+\dot+}c^A$ for the ghost jets. The color factor is six times the cube of $(\epsilon_{3DR})$, hence $-6\epsilon_{3DE}$. The two-loop fan gives
$$
Q_2\mathcal O_0\big|_{\gamma^1=\gamma^2=\gamma^3=t_3}
=\frac{ig^2}{32768\pi^4}
\left(u^1_{++}u^2_{--}-2u^1_{+-}u^2_{+-}+u^1_{--}u^2_{++}\right).
$$
In the relative complex $Q_0 c^A=-g\epsilon^A{}_{BC}c^Bc^C/2$. Differentiate before setting $c=0$. The ordered monomial $p_+^2p_+^3u^2_{--}$ arises only from $Q_0u^1_{++}=-2g p_+^2p_+^3$ in the first term. Its coefficient in $Q_0Q_2\mathcal O_0$ is $-ig^3/(16384\pi^4)$. Thus even on these three closed inputs the raw joining symbol is a nonclosed output. Nilpotence fixes $Q_1^2\mathcal O_0=-Q_0Q_2\mathcal O_0$. The paper's closed quantum product of the three gravitons supplies $\mathcal O_1,\mathcal O_2$ with $Q_2\mathcal O_0+Q_1\mathcal O_1=-Q_0\mathcal O_2$; the transferred second differential of this product vanishes. [Component fan](../../sources/LLC_main.tex#L4080) · [Multi-graviton theorem](../../sources/LLC_main.tex#L4125)

The iterated one-loop fan independently gives the cancelling term. Set $c=0$ and the differentiated scalar jets to zero after both operations, retaining all terms with two first ghost jets, one second ghost jet and the three undifferentiated scalar spectators. Only a first $\beta_i,\beta_j$ contraction followed by a contraction of its differentiated $\gamma^k$ with the remaining $\beta_k$ survives. Let $p_a^D=\partial_{+\dot a}c^D$. With $A,B$ the incoming adjoint colors, the second fan's component is
$$
Q_1\big((\partial_{+\dot b}\gamma^k)_A\,\beta_k^B\big)\big|_{cc}
=\frac{g^2}{64\pi^2} C_{ARD}C_{BRE}
\left(\frac23 u^D_{ba}p^{Ea}+\frac13 p_a^D u_b^{E\,a}\right).
$$
Indeed the averages of $1-r_1$ and $r_2$ on the triangle $r_1,r_2\ge0$, $r_1+r_2\le1$, are $2/3$ and $1/3$; exchanging the two spokes gives the same expression. For all three spectators equal to $t_3$, the first color sum doubles the color-two ghost term. The second color contraction and commuting the odd $Q_1$ past that spectator each give a minus. The coefficient of $p_+^2p_+^3u^2_{--}$ in $Q_1^2\mathcal O_0$ is consequently
$$
3\left(\frac12\right)^3\frac{ig}{32\pi^2}\frac{g^2}{64\pi^2}\frac13
=\frac{ig^3}{16384\pi^4}.
$$
The three factors of $1/2$ are the trace norms and the factor three counts cyclic choices. Polarizing the three independent scalar colors gives $27$ tests. In their full exterior ghost basis all $168$ nonzero coefficients of $Q_1^2\mathcal O_0$ cancel the corresponding coefficients of $Q_0 Q_2\mathcal O_0$ before cohomology. This checks the boundary identity for arbitrary $SU(2)$ scalar spectators in this jet projection. [Displaced one-loop fan](../../sources/LLC_main.tex#L2697) · [Exact two-fan replay](code/ternary_tree.py)

The cancellation extends to every color Lie algebra with an invariant nondegenerate metric. Polarize the two first-jet colors and the second-jet color separately. Both contractions are trivalent color trees on these three colors and the three scalar spectators. Root each tree on the second-jet color. At every vertex cyclic invariance of $\langle[A,B],C\rangle$ puts the parent edge last; the two children give a Lie bracket. Expand these brackets as commutators in the free associative algebra on the other five external colors. With the two first jets in fixed exterior order there are twelve spinor choices. Eight give zero in each diagram sum. Each of the other four gives 48 nonzero ordered words, with opposite coefficients in the iterated one-loop graph and the classical variation of the two-loop graph. Equality of the free Lie polynomials implies the identity for every metric Lie algebra, including coincident external colors by substitution. [Universal integer word expansion](code/ternary_color_identity.py)

For the primary bulk comparison each of the three inputs has tangent jet weights $(0,0)$. Both dictionaries give tangent bulk weights $(-1,-1)$, for a Laurent operator and for a candidate regular source. An exponent in $z_a$ is nonnegative and each $\xi_{z_a}$ occurs at most once. Both $\xi_{z_1},\xi_{z_2}$ must consequently occur in every input, with no $z$ powers. The two copies of the background $F$ have the same tangent bivector. A primary cubic tree with four physical legs and two explicit flux legs has four cubic vertices and three internal edges. Root it at the testing leg. At any furthest vertex the two incoming leaves each carry the tangent bivector, so their product at the wedge-trace vertex is zero. Their Schouten bracket also vanishes before taking cohomology. This kills every tree in the stated sector for any propagator on the internal edges. It also holds for the diagonal flavor weights. A bulk comparison with the nonzero ordinary-product cochain therefore requires a cochain-level map of composite operators or different external kernels. [Bulk weights and both dictionaries](NOTE.md#6-results-of-the-parallel-checks-claude-2026-11-21) · [Exact boundary coefficient and all 4225 bulk pairs](code/ternary_tree.py)

### The complete marked product for the three cyclic gravitons in SU(2)

Use the ordered three trace factors and the fundamental trace of the preceding calculation. Write the coefficients of the two bracketings as $\mathcal O_j$ for the coefficient of $\hbar^j$. At constant Cartan fields all $Q_0$ commutators vanish and every $Q_n$ with $n\ge1$ has differentiated output. Evaluation on that locus sends both products to
$$
\tfrac18(\beta_1)^3(\beta_2)^3(\beta_3)^3
(\gamma^1)_3(\gamma^2)_3(\gamma^3)_3\ne0 ,
$$
where the $3$ on each component is the Cartan color. The corresponding quantum class is free of $\hbar$-torsion. Its ordinary-product $Q_2\mathcal O_0$ has 63 nonzero color and spinor monomials.

In the left and right bracketing, respectively, the one-loop triangle gives $486,198$ nonzero monomials in $\mathcal O_1$, and the two-loop expression gives $1086,396$ in $\mathcal O_2$. Exact exterior-polynomial evaluation in each case gives
$$
Q_0\mathcal O_1=-Q_1\mathcal O_0,\qquad
Q_0\mathcal O_2=-Q_1\mathcal O_1-Q_2\mathcal O_0,\qquad
Q_1\mathcal O_2=0 .
$$
There are at most two $\beta$ letters in $\mathcal O_1$ and one in $\mathcal O_2$. Each derivative with respect to $\beta,\gamma,c$ contributes respectively one, two, three odd coordinates to the common Berezin integral. Hence $Q_2\mathcal O_1=Q_2\mathcal O_2=0$ and $Q_n=0$ on all three coefficients for $n\ge3$. Both quadratic truncations already solve the all-order Ward equation. Using the two Hessians $u^A_{++},u^A_{--}$ defined above, the coefficients of $u^1_{++}u^2_{--}(\gamma^3)_1(\beta_3)^2$ in $\mathcal O_2$ are
$$
\frac{g^2}{98304\pi^4}\quad\text{on the left},\qquad
\frac{g^2}{196608\pi^4}\quad\text{on the right}.
$$
The upper $2$ on the last field is its adjoint color. The difference of the two bracketings is a quantum boundary. An explicit primitive in the relative gauge complex has $72,315,6$ component monomials at orders one, two, three. At the first two orders the replay solves for four quartic and sixteen cubic spin and color contractions. At order three it uses the quadratic ghost primitive below. [Exact cochain primitive](code/ternary_associator.py) [Marked product](../../sources/LLC_main.tex#L3075) · [Rational replay](code/ternary_product.py)

The degree loss $\deg_\theta\mathcal O_n=9-3n$ terminates the marked products after $\hbar^3$. At that order three ghost jets remain. An interaction–interaction edge reduces the number of interaction output letters by one. With bare input factors the three-loop term would therefore require three interactions each joined to both insertions. The insertion consisting of one graviton has only two letters. For the iterated term consisting of a two-loop inner product followed by a triangle, a derivative of the inner ghost-only result already consumes all three odd coordinates. Thus the sole contributing graph is the direct two-loop fan attached to the inner one-loop triangle and the remaining graviton. The triangle leaves a first ghost jet and a differentiated spectator. For the fan, a monomial $r_1^i r_2^j r_3^k r_4^l$ integrates, including its symmetry factor, to
$$
\frac{j!\,l!}{2(i+1)(k+1)(j+l+3)!}.
$$
The sums over the displaced holomorphic derivatives give
$$
\begin{pmatrix}\mathcal O_{3,\mathrm{left}}\\ \mathcal O_{3,\mathrm{right}}\end{pmatrix}
=\frac{g^3}{(32\pi^2)^3}\begin{pmatrix}1/192\\1/576\end{pmatrix}
C_{ABC}(\partial_{+\dot\alpha}c^A)
(\partial_{+\dot\beta}\partial_{+\dot\gamma}c^B)
(\partial_+{}^{\dot\alpha}\partial_+{}^{\dot\beta}\partial_+{}^{\dot\gamma}c^C).
$$
All $36$ color and spinor components agree. The $SU(2)$ structure constants are $C_{123}=1$ and its antisymmetric permutations. With the adjoint metric $\delta_{AB}$, direct differentiation of $Q_0c^A=-g C^A{}_{BC}c^Bc^C/2$ gives
$$
Q_0\!\left[(\partial_{+\dot\alpha}\partial_{+\dot\beta}\partial_{+\dot\gamma}c^A)
(\partial_+{}^{\dot\alpha}\partial_+{}^{\dot\beta}\partial_+{}^{\dot\gamma}c_A)\right]
=-6g C_{ABC}(\partial_{+\dot\alpha}c^A)(\partial_{+\dot\beta}\partial_{+\dot\gamma}c^B)
(\partial_+{}^{\dot\alpha}\partial_+{}^{\dot\beta}\partial_+{}^{\dot\gamma}c^C).
$$
Every higher fan on the two-ghost primitive vanishes by Berezin degree. The cubic correction in either marked representative is therefore quantum exact. [Displaced fan and primitive](code/ternary_product_cubic.py)

### Exact criterion for an irrational genus-three planar descendant

Take four conserved real planar vertex momenta. Use the existing planar bracket $[p_i p_j]=\det(p_i,p_j)$ in the oriented coordinate basis. Multiply the coefficient of $\zeta(3)$ in the tetrahedron period by $485100$. Its value is
$$
\begin{aligned}
&\sum_{\substack{i\in\{1,2,3,4\}\\\mathrm{vertices}}}
 \prod_{j\ne i}[p_i p_j]^2\\
&\quad+\left(\sum_{\substack{1\le i<j\le4\\\mathrm{edges}}}[p_i p_j]^2\right)
\left([p_1p_2]^2[p_3p_4]^2+[p_1p_3]^2[p_2p_4]^2+[p_1p_4]^2[p_2p_3]^2\right).
\end{aligned}
$$
Substitute $p_4=-p_1-p_2-p_3$ into all 296 monomials of the integrated tetrahedron's zeta column. Subtraction of this bracket expression, with the common denominator cleared, gives the zero polynomial in the six independent momentum components. All proper faces have rational coefficients. For real momenta the displayed expression is nonnegative. A zero vertex charge or rank at most one kills every monomial. If all four charges are nonzero and span the plane, conservation excludes three collinear charges. The charges then admit a matching into two nonparallel pairs, and its term is strictly positive.

At genus three on a two-dimensional target, the Betti power is six. A graph with a tetrahedron minor needs at least six edges, with equality only for the tetrahedron. Its six edges already exhaust the power, and its three independent cycles exhaust the genus. Give the retained marks levels $k_i$ summing to their number. At each of the four genus-zero vertices the three flags have level zero. A set $S$ of marks at that vertex contributes only if the sum of its levels equals $|S|$, in which case its cotangent integral is $|S|!/\prod_{i\in S}k_i!$. Sum over partitions of the labeled marks into four nonempty unordered groups. Each contribution is the product of the four cotangent integrals and the displayed coefficient divided by $485100$, evaluated at the four sums of momenta. Partitions with an empty group have zero zeta coefficient. This proves that the full vertex at rational momenta is irrational exactly when at least one allowed partition has four nonzero group momenta spanning the plane, since every other term is rational and $\zeta(3)$ is irrational. Three marks always give a rational vertex. Four level-one marks with momenta $(1,0),(-1,0),(0,1),(0,-1)$ give
$$
\frac{464971}{2514758400}+\frac{2}{121275}\zeta(3).
$$
[Exact polynomial identity and descendant replay](code/planar_sextic_period.py) · [Tetrahedron coefficients and their Schwinger integrals](code/mixed_k4_interior.py)

### A genus-three fivefold scalar measurement from the planar descendant

Use the flat fivefold scalar $\phi(0)$ and the repository's $A_k$ descendant convention. In the fivefold coordinate covector basis $e^1,\ldots,e^5$, take the first four retained scalar momenta to be
$$
e^1-4e^3-4e^4-3e^5,\quad -e^1,\quad e^2,\quad -e^2
$$
at respective levels $(2,2,1,1)$. The absorbed primary root has momentum $e^3$. The primary vector has momentum $e^4$ and polarization $e_3$, giving root contraction one and vector divergence zero. Add three labeled primary scalars on each of the axes $e^3,e^4,e^5$, each of unit momentum on its own axis. The total momentum is zero.

In the signed sum each extra scalar remains a separate mark or fuses into one of the four retained descendants, with one minus sign and one level removed per fusion. The root and vector stay fixed. This gives the pullback of the retained cotangent classes. The three ordered, degree-three Abel maps each have degree $3!$, so integration of the nine extra marks multiplies the planar vertex by $216$. The root and vector become zero-momentum primary marks on the plane. Forgetting them twice by the string equation lowers two of the four cotangent levels in order. The rational and zeta coefficients of the planar integral are obtained from
$$
\begin{aligned}
&2\left(\frac{464971}{2514758400}+\frac{2}{121275}\zeta(3)\right)
+2\frac{143}{1451520}+8\frac{4909}{43545600}+2\frac{289}{3628800}\\
&\hspace{12mm}=\frac{4093949}{2514758400}+\frac{4}{121275}\zeta(3).
\end{aligned}
$$
Thus the signed sum of $(2\pi i)^{15}Q_3\phi(0)$ matrix elements is
$$
\frac{4093949}{11642400}+\frac{96}{13475}\zeta(3).
$$
The computation uses the flat worldsheet action and fixes its field coordinates. The rational term is obtained from complete four-mark stable-graph sums. The zeta coefficient also follows from the two allowed nonzero four-group partitions before forgetting the root and vector. [Fivefold pullback and heat normalization](../iib-bcov/README.md#an-all-genus-scalar-target-measurement-of-the-fivefold-differential) · [Complete replay](code/planar_sextic_period.py)

Keep all momenta, the polarization and the heat convention fixed. Move one cotangent power from the fourth retained scalar to the third, giving levels $(2,2,2,0)$. There are now just three positive-level marks. Four nonempty tetrahedron groups cannot each have mean level one. The same two applications of the string equation give the signed fivefold response
$$
216\left(2\frac{1661}{36288000}+2\frac{143}{1451520}
+4\frac{4909}{43545600}+\frac{289}{3628800}\right)=\frac{11147}{63000}.
$$
At levels $(2,2,1,1)$ and $(2,2,2,0)$ the sums contain respectively $54721$ and $28072$ labeled fusion assignments. To count them, choose a subset of the nine labeled primary marks for each retained scalar, with disjoint subsets and cardinalities bounded by its level. Since the second response is rational and nonzero and the first has a nonzero rational multiple of $\zeta(3)$, their ratio is irrational. Their rational span in these fixed bulk coordinates contains both $1$ and $\zeta(3)$. A common rescaling cannot make every elementary matrix element in the two sums rational. [Independent all-label and vertex-genus replay with --labeled](code/planar_sextic_period.py)

### A single moving fivefold descendant from the three-loop Feynman Ward equation

In the $(2,2,1,1)$ probe, fuse two $e^5$ labels into the first hard scalar, two $e^4$ into the second, and the remaining $e^4$ into the fourth. The third hard scalar at $e^2$ retains level one. The three free $e^3$ labels and the absorbed root give four identical primary momenta $u=e^3$. The vector remains at $e^4$ with polarization $e_3$. Direct evaluation gives the individual ordinary coefficient
$$
(2\pi i)^{15}Q_3\phi(0)
=\frac{11180467393}{23538138624000}-\frac{220697}{567567000}\zeta(3).
$$

In the ten-primary Ward equation, the four bivalent current jumps on the bare three-loop graph lie in the three-dimensional cycle space tensored with $u$. Their linear dependence has coefficient sum zero and is independent of the positive-resistance Kirchhoff relation. Every maximal current minor vanishes. In the hexagon the determinant of its five singleton waves uses at most one $u$. The genus-two primary at the remaining vertex then contains at least three $u$'s. Its four five-row determinants cannot provide two occurrences to each of these three slots. Thus only the cubic tree with the genus-three primary vertex survives.

Let $V$ be the ordinary positive-convention nine-primary trace after fusing one $u$ into an off-axis momentum $p_i$. The three remaining $u$'s each occur twice across the six determinants, so each determinant uses one $u$. Its other columns can be reduced modulo $u$. Hence every off-axis choice gives the same $V$. Fusing two $u$'s gives $4V$, and fusing two off-axis waves gives zero. The sum of the six off-axis momenta is $-4u$. The tree potentials at an off-axis and an axial marking, respectively, are
$$
4V(p_i+u),\qquad 3(2u)(4V)+(6u-4u)V=26Vu.
$$
Subtracting the common potential $4Vu$ leaves the longitudinal coefficients $4V$ and $22V$. For the displayed probe the level-one wave is off-axis and the root contracts to one. Both period components agree from the full eighteen-term primary trace and from fifteen integrated Feynman measurements. The replay checks all 45 pair fusions, both cubic graph cores at six source placements, and the genus-two hexagon. [Fivefold graph and Ward replay](code/fivefold_fused_primary.py) · [Primary Feynman coefficients](../iib-bcov/operator/three_loop_coefficients.json)

The six-fusion layer is entirely primary. A unit-determinant shear fixing $e^3,e^4,e^5$ moves the transverse momentum within each opposite planar pair onto its positive leg. Only the totals sent into the first and second pairs matter. There are 45 distributions of the three sets of donor labels, with four labels entering the first pair, two entering the second, and three remaining free. In the five-fusion layer the same Ward calculation applies whenever the root or vector axis retains all three of its scalar labels. The heat-stripped subtotals below use the ordinary coefficients before applying the fusion signs.

| Levels | Six fusions | Five fusions with four equal axial marks |
|---|---:|---:|
| $(2,2,1,1)$ | $9024154034239/1212570777600-(56134219613/9356256000)\zeta(3)$, 15120 assignments | $174245248393/65383718400-(17206396/7882875)\zeta(3)$, 2160 assignments |
| $(2,2,2,0)$ | Half of the first row, 7560 assignments | Half of the first row, 1080 assignments |

Consequently the first signed response minus twice the second cancels both displayed layers. Its value $-129911/58212000+(96/13475)\zeta(3)$ comes entirely from the remaining moving descendants. [Fusion multiplicities and exact subtotals](code/fivefold_fused_primary.py)

### Isolating higher descendants in the difference of the two planar probes

Fix all labels sent into each of the two opposite planar pairs. A linear change fixing $e^3,e^4,e^5$ moves transverse covectors between the two legs of each pair while preserving the pair sum. Another exchanges those legs. Both preserve every scalar vertex: a determinant of $-1$ contributes its sixth power in the genus-three trace. The absorbed root evaluates to one in both probes. If two donors enter the last pair, its remaining levels are $(0,0)$, with two labeled splits from $(1,1)$ and one from $(2,0)$. If one donor enters, its two placements from $(1,1)$ give the same vertex as the single placement from $(2,0)$. Thus the first signed response minus twice the second receives contributions only from donor fusions into the first pair. The complete five- and six-fusion layers cancel. With zero through four fusions the remaining counts are $1,18,144,504,756$. Interchanging $e^3,e^4$ in the symmetric scalar trace leaves 24 pairs of ordinary vertices.

The cotangent pullback reduces this further. Temporarily treat the root and vector as ordinary primary marks of the symmetric scalar vertex. Adding the root to the three forgettable $e^3$ marks gives an alternating pullback integral of zero: four vertical complex marking directions enter only the three-dimensional Jacobian of that axis. Splitting this zero sum according to whether the root is fused forces one $e^3$ fusion. Applying the same argument to the vector forces one $e^4$ fusion. The two extra signs multiply to one. The exchange of the opposite planar legs still cancels every term with a fusion into the last pair, including fusions of the two distinguished axial marks. At least one primary on each of the root and vector axes remains available to recover the same rooted contraction.

For each row below, fuse the indicated numbers of axial marks into the first pair. Keep their totals on the positive planar leg by the above shear. The row contributes its integer weight times the ordinary vertex with last-pair levels $(1,1)$ minus twice the ordinary vertex with last-pair levels $(2,0)$. The second column gives the residual first-pair levels, in either order.

| Fused counts along $(e^3,e^4,e^5)$ | First-pair levels | Weight |
|---|---|---:|
| $(1,1,0)$ | $(0,2)$ | $2$ |
| $(1,1,0)$ | $(1,1)$ | $2$ |
| $(1,1,1)$ | $(0,1)$ | $-18$ |
| $(1,1,2)$ | $(0,0)$ | $18$ |
| $(1,2,0)$ | $(0,1)$ | $-36$ |
| $(1,2,1)$ | $(0,0)$ | $108$ |
| $(1,3,0)$ | $(0,0)$ | $36$ |
| $(2,2,0)$ | $(0,0)$ | $54$ |

These eight differences sum to $-129911/58212000+(96/13475)\zeta(3)$. At least one of these eight differences of moving descendants is irrational in heat units. As an independent graph check of the fusion and rank reduction, the genus-one version has one donor per axis and first-pair levels $(1,1)$. Full labeled polarization and the reduced hexagon-plus-tree Ward replay give the two individual sums $5/3$ and $5/6$. [Labeled fusion recurrence and both graph controls](code/planar_fusion_difference.py)

### Longitudinal Ward trees at maximal collinear multiplicity with moving descendants

Take $m$ primary momenta $a_i u$ on one line. Let $k_j$ be the levels on the remaining hard momenta $p_j$, with $K=\sum_j k_j$. For $g=2,3$ the flat fivefold vertex vanishes when $m>g+K$. At $m=g+K$ let $V_{g;0}$ denote the primary vertex with $g$ unit momenta $u$ and the same projected hard momenta modulo $u$. Adjust any one longitudinal hard component to conserve momentum. The saturated primary is independent of this choice: its $2g$ five-row determinants use every one of the $g$ line momenta twice. The ordinary vertex with the stated levels is
$$
\frac{K!}{\prod_j k_j!}\left(\sum_{\substack{S\subseteq\{1,\ldots,m\}\\ |S|=g\text{ primary line marks}}}\prod_{i\in S}a_i^2\right)V_{g;0}.
$$

For a direct graph proof lower one positive hard level by one in the Ward source, giving levels $l_j$ of total degree $K-1$. In a graph of positive geometric loop number $\ell$, a bridge kills the heat cofactor. An ordinary vertex with $b\ge2$ internal ends and level sum $L$ contains $L+3-b$ external slots. It can carry at most $L$ of the line marks, except at a bivalent vertex containing only a line source. Each such exception contributes a current jump in the $\ell$-dimensional cycle space times $u$. More than $\ell$ exceptions kill every maximal minor. A positive-genus vertex of genus $h$ carries at most $h+L$ line marks by induction on genus. Summing these bounds gives at most $g+K-1$ line marks. A one-edge tree with two positive-genus vertices obeys the same bound. The genus-one bound starts from the zero primary and the hexagon with this same current-jump argument.

Ground the remaining tree at its genus-$g$ end. For the potential at hard mark $j$, the classical end must contain that mark. If it contains $b$ hard marks, it can absorb $L+2-b$ line marks. Applying the lower-arity bound at the genus-$g$ end forces $b\le1$ in the saturated case. For $m>g+K$ every hard potential is zero. At equality the classical end contains $p_j$ and a subset $J$ of $l_j+1$ axial labels. Its classical multinomial is one. In the remaining primary symbol the $g$ surviving line labels form a subset $S$ disjoint from $J$, and their strengths contribute $\prod_{i\in S}a_i^2$. For fixed $S$ the other $K$ labels admit $\binom K{l_j+1}$ choices of $J$, or $\binom{K-1}{l_j}$ if one additional prescribed label must belong to $J$. After including the lower-arity multinomial, the coefficient of $p_j$ is the displayed expression with $k_j=l_j+1$ and all other levels unraised. The term in $u$ is independent of $j$: it is the sum of $(K-1)!V_{g;0}(\prod_{i\in S}a_i^2)(\sum_{i\notin S}a_i)/\prod_j l_j!$ over $S$. Removing this common potential gives the claimed longitudinal Ward coefficient. The argument first uses two nonparallel hard momenta to fix the common potential, then extends by polynomiality.

For a fivefold elementary-scalar probe take five unit line momenta $u=e^3$, including the absorbed root. The six hard momenta are
$$
e^1-4e^3-4e^4-e^5,\quad-e^1+2e^4,\quad e^2-e^3,\quad -e^2+e^4,\quad e^5,\quad e^4.
$$
The final wave is a primary vector polarized by $e_3$. With level two on the third hard scalar, the ordinary heat-stripped $Q_3\phi(0)$ equals
$$
\frac{11180467393}{9415255449600}-\frac{220697}{227026800}\zeta(3).
$$
With level one on each of the first and third hard scalars it is exactly twice this value. The multipliers of the saturated nine-primary trace are $\binom53=10$ and $2\binom53=20$. Ten independent genus-two graph checks cover moving levels $(2)$, $(1,1)$, $(2,1)$, $(1,1,1)$ and $(3)$, at both unit and unequal signed axial strengths. [Consecutive Ward cuts and direct genus-two graphs](code/saturated_descendant_ward.py)

### A fixed elementary-scalar polarization whose first response is at three loops

In the preceding eight-row table the remaining axial multiplicities, in decreasing order, determine the ordinary scalar trace. The first pair carries all transverse recoil; changing the basis permutes its three axial components along with the primary marks. Thus there are five entries. In order they have axial multiplicities and residual first-pair levels
$$
((3,3,3),(0,2)),\quad ((3,3,3),(1,1)),\quad
((3,3,2),(0,1)),\quad ((3,3,1),(0,0)),\quad ((3,2,2),(0,0)).
$$
Every entry takes the difference of last-pair levels $(1,1)$ and twice $(2,0)$ as above. Forgetting all nine axial primaries in the first two entries gives two independent rows on this ordered list:

| Row on the five entries | Identical planar four-mark evaluation times $216$ |
|---|---:|
| $(1,0,-9,9,27)$ | $253/168000$ |
| $(0,1,-18,18,54)$ | $-7631/2910600+(48/13475)\zeta(3)$ |

Indeed one donor can enter any of nine positions, and two from the same or distinct axes have respectively $9$ or $27$ choices. With first-pair levels $(1,1)$ these multiplicities double. The ordinary second entry minus twice the first is therefore $-327949/58212000+(48/13475)\zeta(3)$. Three combinations of the five entries remain undetermined by these equations.

There is a fixed field profile measuring this ordinary combination. Give four independently labeled scalar waves momenta
$$
e^1-3e^3-3e^4-3e^5,\quad -e^1,\quad e^2,\quad -e^2 .
$$
On each of the first two waves set the primary, first- and second-descendant amplitudes to $(1,t,-t^2)$ times that wave's formal coefficient. Use $(1,s,-s^2)$ on each of the last two. Add the primary scalar $e^{z^3}+e^{z^4}+e^{z^5}$ and one primary vector $e^{z^4}e_3$. In $Q_g\phi(0)$ extract $t^2s^2$ and the part linear in the vector and in each of the four labeled scalar coefficients. In the flat fivefold action its values satisfy
$$
Q_0\phi(0)=Q_1\phi(0)=Q_2\phi(0)=0,\qquad
(2\pi i)^{15}Q_3\phi(0)=-\frac{327949}{1397088000}+\frac{2}{13475}\zeta(3)
\quad\text{on this polarization}.
$$

The pair coefficient is $(1,1)-(2,0)-(0,2)$ in descendant levels. For any fixed nonempty set of donors fused into either pair, transferring a donor between its legs preserves the transverse pair sum. A shear transfers the momentum while fixing all other legs. For one fused donor the two level-one survivors cancel on either leg; for two fused donors the weights are $2-1-1=0$. More than two fusions exceed the pair's level. Thus the ordinary pair polarization already equals the pulled-back cotangent measurement when all the other primary marks, including the root and vector, are forgotten in the scalar trace. At genus $g$ there are $2g+3$ such marks in the three-dimensional complementary target. Their vertical period rank is at most $3g$. This proves the lower-genus vanishings. At genus zero the classical multinomial gives the same cancellation.

At genus three the fiber coefficient for any nine transverse covectors is the sum over their ordered partitions into three triples of the product of the three squared determinants. To derive it, expand the vertical Abel determinant using the three normalized holomorphic one-forms on the curve. Integration at each labeled point identifies its holomorphic and antiholomorphic form labels by the bilinear relations. Each form label occurs three times. The resulting three determinants per side give the stated sum; it is also the squarefree coefficient of the cube of the determinant of the weighted sum of the nine outer products. For three unit waves per axis this is $216$. The degree is unchanged by transverse recoil at the four hard marks and supplies the full change for unequal transverse momenta while the planar momenta stay fixed.

If the primary scalar multiplicities are $(a,b,c)$ on the three axes at genus three, $a+b+c=7$. The vector adds one fourth-coordinate wave. The balancing root has transverse components $(3-a,2-b,3-c)$. Only $(2,2,3),(3,1,3),(3,2,2)$ admit three independent triples. The root evaluates to $1,0,0$ on the vector polarization in this order. Dividing the sole surviving labeled answer by $2!2!3!$ gives the fixed-profile coefficient. All ten one-loop multiplicities were also tested in the hexagon/tree Ward equation. On the two-loop multiplicities $(2,1,2)$ the four ordinary Ward vertices, in the order $(0,2;0,2),(0,2;1,1),(1,1;0,2),(1,1;1,1)$, give $42253/12960,21407/3240,1069/162,43321/3240$. Their combination with weights $(4,-2,-2,1)$ is zero. [Planar labeled graphs, determinant identity and Fourier replay with --two-loop](code/planar_fusion_difference.py)
