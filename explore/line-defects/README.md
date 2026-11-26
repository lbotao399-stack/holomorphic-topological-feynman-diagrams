# Schwinger currents with fixed defect vertices

The ordered twistor defect OPE coefficient at arbitrary loop and jet order is given by the [ordered interval integral and exact recurrence](higher_f/README.md).

For each complex coordinate, split the graph incidence matrix $B$ into columns $B_I,B_D$ for integrated and fixed vertex coordinates. On a straight complex line the normal coordinates of each defect vertex are fixed at zero. Its tangent coordinate can be either fixed or integrated as prescribed by the observable. Choose $t_e>0$ so the scalar Gaussian on edge $e$ is $\exp(-|Z_{e(0)}-Z_{e(1)}|^2/t_e)$, and write $L=B^{\mathsf T}t^{-1}B$. Each integrated connected component must meet a fixed vertex, or one translation mode must first be removed. At fixed boundary positions $Z_D$, completing the square gives

$$
\langle Z_I\rangle=-L_{II}^{-1}L_{ID}Z_D,\qquad
\langle (Z_I-\langle Z_I\rangle)(\overline Z_I-\langle\overline Z_I\rangle)^{\mathsf T}\rangle=L_{II}^{-1}.
$$

The boundary Gaussian has matrix $L_{DD}-L_{DI}L_{II}^{-1}L_{ID}$. The mean currents are $t^{-1}(B_D Z_D+B_I\langle Z_I\rangle)$, with conjugation when the numerator uses barred coordinates. Each integrated complex coordinate contributes the inverse determinant of its grounded Laplacian and the scalar measure factor of the source. Pull back every propagator differential before doing the Gaussian moments. These formulas apply to the finite Schwinger integrand before any small-time extension. [Grounded network and exact moments](chain.py)

**One-loop boxes in $\mathbb C^3$.** The two bulk vertices $X,Y$ form a chain between fixed endpoints $(w,0,0),(z,0,0)$. Divide its three resistances by their sum to get $u,v,s$ with $u+v+s=1$. The mean tangent coordinates are $(v+s)w+uz$ and $sw+(u+v)z$. In either normal direction the covariance divided by the sum of resistances has diagonal entries $u(v+s),s(u+v)$ and off-diagonal entry $us$. Its determinant on the simplex equals $uvs$.

In Fernández (7.11) the numerator contains $[\bar x,\bar y]x^2y^1$ in the normal directions. The two Wick pairings give minus this determinant. After the scale integral and the three Schwinger factorials, the remaining integral is

$$
-\frac{(-2\pi i)^6}{(2!)^3|w-z|^2}\int_{u,v>0,\ u+v<1}\big((1-u)w+uz\big)\,du\,dv
=-\frac{(-2\pi i)^6(2w+z)}{48|w-z|^2}.
$$

Restoring the differential and the propagator factors gives $\mathcal M_{0000}=(2w+z)d\bar z_0/[96\pi^2(w-z)]$, $z_0=(w+z)/2$. For any transverse holomorphic polynomial $F$, contracting the two indicated barred variables acts by $uvs(\partial_{x^1}\partial_{y^2}-\partial_{x^2}\partial_{y^1})$ on $F$, times the square of the total resistance. The remaining Gaussian pairings give ordinary monomials in the three covariance entries. In a term with $i,j,k$ diagonal-$X$, cross, and diagonal-$Y$ pairs, the remaining endpoint powers after removing $uvs$ are $u^{i+j-k_1-k_2}s^{j+k-l_1-l_2}$. Both exponents are nonnegative because each endpoint barred variable takes one pair incident there. The scale integral at transverse balance has exponent one. Thus arbitrary ordered transverse derivatives and arbitrary holomorphic tangent polynomials at both vertices reduce to finite simplex factorials. The replay reproduces 3,696 individual coefficients of (C.10) through total transverse order eight, and separately checks 25 choices of tangent powers at both vertices. The leading OPE coefficients are $1/(32\pi^2)$ and $h^\vee/(48\pi^2)$ with the paper's Killing form. [Fernández, §§7 and C](https://arxiv.org/abs/2302.14292) · [Coefficient replay](check_gauge.py)

With Bittleston's vertex factors the same network gives, in the order of (C.11), (C.14) and the first line of (C.15),

$$
-\frac{22w+13z}{35\pi^2(w-z)},\qquad
-\frac{2(4w+3z)}{35\pi^2(w-z)},\qquad
-\frac{3}{16\pi^2(w-z)}.
$$

The first two values give the simple-pole coefficient $2/(5\pi^2)$ in the gravity OPE of §5. [Bittleston, §§5 and C](https://arxiv.org/abs/2211.06417) · [Independent Schwinger calculation](chain.py)

### The commuting-current coefficient in (7.3.5)

Use the pole order $J_X[1,0](0)J_X[0,1](z)$ of Costello–Paquette (7.3.5), with the same Cartan element $X$ on each of the three incoming currents and on the output pairing. Take normalized residues in $z,w$, weighting the third current $J_X[0,0](w)$ by $w$. The level-zero residues of the two composites agree because their contracted color tensor is symmetric. For one of them the $w$ contour about the displaced current, followed by the normal-product contour, gives

$$
K^{ij}\operatorname{tr}\!\left(
X\big[[X,t_i],[[X,t_j],X]\big]\right)
=-K^{ij}\operatorname{tr}\!\left((\operatorname{ad}_X^2t_i)(\operatorname{ad}_X^2t_j)\right)
=-\operatorname{Tr}_{\rm ad}(X^4).
$$

Here $K^{ij}$ inverts the trace pairing on the color basis $t_i$. Each tree gauge residue enters twice, so reversing the common sign of the tree gauge pole leaves this contraction fixed. The $F[0,0]$ intermediate state contributes $-\widehat\lambda_{\mathfrak g}^{\,2}\operatorname{tr}(X^2)^2$ in the first contour order. Each of the two alternative contours gives $+\widehat\lambda_{\mathfrak g}^{\,2}\operatorname{tr}(X^2)^2$ via $E[1,0]$ or $E[0,1]$. This uses the full axion OPE in §4, with $\widehat\lambda_{\mathfrak g}$ multiplying the double pole as well. Equating the contours yields

$$
-2 C\,\operatorname{Tr}_{\rm ad}(X^4)
=3\widehat\lambda_{\mathfrak g}^{\,2}\operatorname{tr}(X^2)^2,\qquad
C=-\frac{3\widehat\lambda_{\mathfrak g}^{\,2}}{2\lambda_{\mathfrak g}^{\,2}}
=-\frac{1}{8(2\pi i)^3}.
$$

The last equality uses the paper's quartic trace identity and its axion coupling. The residue fixes the common $C$ in both lines of (7.3.5). At commuting external colors, reversing normal-product order can only add a derivative contracted with an antisymmetric structure constant, which vanishes. Fernández's ordered pole has coefficient $(2\pi i)C=1/(32\pi^2)$. His two cubic vertices and three propagators contribute a net $2\pi i$ when the BF action carries $1/(2\pi i)$, while the two defect couplings occur on both sides of the residue equation. In the later associativity paper the displayed $JJ$ and $J\widetilde J$ formulas put different prefactors in front of its parameter $C$. The normalized Cartan residues here use the OPEs of the earlier paper directly. [Costello–Paquette, §§4, 7.3](https://arxiv.org/abs/2201.02595) · [Later contour computation](https://arxiv.org/abs/2204.05301) · [Fernández](https://arxiv.org/abs/2302.14292) · [Two exact color probes](cartan_residue.py)

### The one-loop three-point form factor in self-dual supersymmetric Yang–Mills theory

For the completed insertion $\mathcal O_{B^2}(x)=\frac1{12}\int d^8\tilde\theta\,\log\det(\bar\partial_{L_x}+\mathcal A|_{L_x})$, the incidence relation $\eta_i=w_i^\alpha\lambda_\alpha$ identifies the chiral supertranslation at fixed $x$ with a translation of an integrated odd coordinate. Its Berezin integral vanishes. In a supersymmetry-preserving prescription apply $Q_{3,\alpha}$ to $\chi_a^{13}[1](z_1)K_b^2[1](z_2)\chi_c^{34}(z_3)$. The middle variation is zero. The minus sign of $Q_3\chi^{34}$ cancels the minus sign for crossing the odd middle field. At separated points the two components give

$$
\begin{pmatrix}
z_1^{-1/2}&z_3^{-1/2}\\
z_1^{1/2}&z_3^{1/2}
\end{pmatrix}
\begin{pmatrix}
\langle\mathcal O_{B^2}|K_a^1[1](z_1)K_b^2[1](z_2)\chi_c^{34}(z_3)\rangle\\
\langle\mathcal O_{B^2}|\chi_a^{13}[1](z_1)K_b^2[1](z_2)K_c^4(z_3)\rangle
\end{pmatrix}=0,\qquad
\det=\frac{z_3-z_1}{\sqrt{z_1}\sqrt{z_3}} .
$$

Thus the coefficient of $\epsilon^{ijkl}f_{abc}[12]$ in §5.2.2 is zero. The Ward identity applies to the sum of the three component terms of the insertion in §5. [Source, §§2, 4, 5](https://arxiv.org/abs/2604.21015) · [Exact replay](sdsym_ward.py)

### Five-dimensional mixed defect

In Oh–Zhou Fig. 12 both end vertices have $t_1=t_3=z_1=z_3=0$. Keep the differential of the endpoint $w$ coordinate on the right. From their propagator (88) the parts selected by the two endpoint integrals are

$$
\left.P_{12}\right|_{d\bar w_1}=\frac{(t_2\,d\bar z_2-\bar z_2\,dt_2)\wedge d\bar w_1}{d_{12}^5},
\qquad
\left.P_{23}\right|_{d\bar w_3}=-\frac{(t_2\,d\bar z_2-\bar z_2\,dt_2)\wedge d\bar w_3}{d_{23}^5}.
$$

The holomorphic vertex derivatives act only on the scalar coefficients. Their product contains the square of the same bulk one-form and vanishes pointwise. In (204), moving $dt_2$ past $d\bar w_1$ gives a minus sign. This changes the sign inside the parenthesis in (132). The separated graph with the paper's $|t_2|\ge\epsilon$ cutoff therefore gives zero. Multiplying each propagator numerator by a smooth radial heat regulator preserves this identity on the diagonals as well. The displayed one-loop diagram does not establish the coefficient in (140). A different local collision term or a pair of defect planes at different times requires its own specification. [Oh–Zhou, §§5.3, B.3](https://doi.org/10.21468/SciPostPhys.10.2.029) · [Exterior algebra replay](oh_zhou.py)

Independently, the scalar integrand printed in (207) can itself be integrated on its stated circular contour $|w_2|=R$. Its ratio to the first summand, using the same measure for the two summands, is

$$
\frac{(1-q)^2(1+5q)}{(1+q)^3},\qquad q=\frac{\epsilon}{\sqrt{\epsilon^2+R^2}}.
$$

Indeed $x=t_2^2+|z_2|^2$ reduces the positive-time radial measure with the factor $t_2|z_2|^2$ to $\pi(x-\epsilon^2)^2dx/4$. The first kernel is $2/x^4$. For the second kernel $(5R^2+2x)/[x^{5/2}(x+R^2)^{5/2}]$ substitute $s=\sqrt{x/(x+R^2)}$. Its integral after division by the first gives one minus the displayed ratio. For (208)'s normalization the first summand is $2\pi^3/(225\epsilon^2)$ on positive times. This evaluates (207) as written. The forms derived from (88) give the pointwise zero above. [Exact radial integration](oh_zhou.py)

**Costello's central parameter at second order.** In §11 the zeroth mode $E_{\alpha\beta}$ maps to the outgoing–incoming marked disk. Its charge-$K$ image is $\hbar^{-1}\sum_{\substack{1\le i\le K\\\mathrm{gauge\ basis}}}\rho_{i\alpha}\frac{\partial}{\partial\rho_{i\beta}}$, whereas the unmarked disk maps to $K$. Taking the gauge trace of (13.0.1), its two adjoint sums cancel and each of the $NK$ framing oscillators contributes one commutator $\hbar$. Uniform-in-$K$ injectivity then gives
$$
(c-N\hbar)D(\varnothing)=\hbar\sum_{\substack{1\le\alpha\le N\\\mathrm{flavor\ basis}}}D(\alpha\Uparrow,\alpha\Downarrow).
$$
In particular $D(\varnothing)$ vanishes in the classical fiber at $c\ne0$, as also follows from Proposition 10.0.5. The marked trace lifts the central identity. A term $\hbar^2D(\varnothing)$ starts at cubic order in this formal algebra. Equation (2.28) in the later DDCA presentation gives the closed zeroth mode as rank divided by $\epsilon_1$. [Costello, §§10–13,16](https://arxiv.org/abs/1705.02500) · [Gaiotto–Rapčák–Zhou, §2.3](https://doi.org/10.1007/s00220-025-05488-y)

At fixed $c$, duality (2.30), followed by the dilation of Lemma 10.0.3, multiplies both $\hbar$ and the marked trace by $c/(N\hbar-c)$. The anti-involution (3.58) fixes both deformation parameters. These transformations permit an arbitrary coefficient for the marked trace times $\hbar^2$. After fixing the bulk constant, inversion of the map in §16 multiplies that second-order coefficient by $-(32\pi^2)^{-3}$. The five-flavor projector for its cubic Poincaré–Birkhoff–Witt symbol annihilates all 120 oriented ribbon trees on its five marked legs. [Central-coordinate and inverse-map check](costello_center.py) · [Ribbon-tree and tadpole check](costello_two_trace.py)
