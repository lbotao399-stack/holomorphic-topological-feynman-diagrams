# $J_1$, the 5d Chern–Simons level, the Vogel factor, and the $e\to0$ limit of twisted M-theory

2026-11-21. Status tags: [chk] = checked by `cs_restriction_check.py` (exact Grassmann algebra); [drv] = derived here, not independently checked; [src] = read in the source named; [open].

## 1. Restriction of $\mathrm{CS}_5$ to the holomorphic–topological fields [drv; chk for (1.3)]

$\mathcal A=\mathcal A^aT_a$, $\mathcal A^a\in\Omega^{0,*}(\mathbb C^2)\hat\otimes\Omega^*(\mathbb R)$ of total degree 1 (odd), no $dz_i$ components; $d=\partial+\bar\partial+d_\tau$, $\partial=dz_i\partial_{z_i}$.

$$\mathrm{CS}_5(\mathcal A)=\mathrm{Tr}\big(\mathcal A\,d\mathcal A\,d\mathcal A+\tfrac32\mathcal A^3d\mathcal A+\tfrac35\mathcal A^5\big),\qquad d\,\mathrm{CS}_5=\mathrm{Tr}\,F^3 .\tag{1.1}$$

Top degree on $\mathbb R\times\mathbb C^2$ needs $dz_1dz_2$, which only $\partial$ supplies; $\mathcal A^3d\mathcal A$ and $\mathcal A^5$ contain at most one $\partial$:

$$\int_{\mathbb R\times\mathbb C^2}\mathrm{CS}_5(\mathcal A)=\int\mathrm{Tr}(\mathcal A\,\partial\mathcal A\,\partial\mathcal A).\tag{1.2}$$

$$\mathrm{Tr}(\mathcal A\,\partial\mathcal A\,\partial\mathcal A)=-dz_1dz_2\,d_{abc}\,\mathcal A^a\partial_{z_1}\mathcal A^b\partial_{z_2}\mathcal A^c,\qquad d_{abc}=\mathrm{tr}(T_a\{T_b,T_c\}).\tag{1.3}$$

With $F_D=\tfrac13\int\Omega\,d_{abc}\,\partial_{z_1}A^a\partial_{z_2}A^bA^c$ (flagship Prop. 6.2, $J_1=aF_D$) and $d_{abc}\partial_1\mathcal A^a\partial_2\mathcal A^b\mathcal A^c=d_{abc}\mathcal A^a\partial_1\mathcal A^b\partial_2\mathcal A^c$ [chk]:

$$\int_{\mathbb R\times\mathbb C^2}\mathrm{CS}_5(\mathcal A)=-3F_D .\tag{1.4}$$

For the bosonic level term $S_\kappa=\frac{\kappa}{24\pi^2}\int\mathrm{CS}_5(A)$, (1.4) gives the first equality below. The second is a hypothesis about the twist of the full supersymmetric level deformation in the BV exponent $e^{S/\hbar}$ [open]:

$$\left.S_\kappa\right|_{A=\mathcal A}=-\frac{\kappa}{8\pi^2}F_D,\qquad J_1=c_{\rm tw}\left.S_\kappa\right|_{A=\mathcal A}\quad\Longrightarrow\quad a=-\frac{c_{\rm tw}\kappa}{8\pi^2}.\tag{1.5}$$

Here $c_{\rm tw}$ is an overall coefficient of the restricted functional. A rescaling of $\mathcal A$ enters cubically in (1.4). Identification with the order-$\hbar$ action term requires the field and volume normalizations, the phase of the action, and the relation between $\hbar$ and the Yang–Mills coupling $g^2$ [open].

$d_{abc}\neq0$ iff $(S^3\mathfrak g^*)^{\mathfrak g}\neq0$ iff $\mathfrak g\supset\mathfrak{sl}_{N\ge3}$ factor (simple case): the algebras admitting a 5d CS level are exactly class (C).

Assuming the twist identification in (1.5), the order-$\hbar^2$ result (flagship Thm 5.10) gives $a^2=\varepsilon U\,N^2/(256\pi^4)$ $\Rightarrow$

$$c_{\rm tw}^2\,\kappa^2=\frac{\varepsilon U}{4}\,N^2 .\tag{1.6}$$

## 2. The pure five-dimensional twist

- For the rank-one minimal twist of pure 5d $\mathcal N=1$ SYM on $\mathbb C^2\times\mathbb R$, set the matter representation to zero in ESW Theorem 9.3. The classical BV theory is perturbatively equivalent to $T^*[-1]\operatorname{Map}(\mathbb C^2\times\mathbb R_{\rm dR},B\mathfrak g)$. Equivalently, the generalized Chern–Simons target in Theorem 9.2 uses $\mathfrak g\ltimes\mathfrak g^*$, with the coadjoint action and the canonical cross-pairing. The flagship action uses one $\mathfrak g$-valued field and its invariant bilinear pairing. [ESW, §9.1.1, Theorems 9.2–9.3](https://arxiv.org/pdf/2002.10517#page=67) · [Example 1.5](https://arxiv.org/pdf/2002.10517#page=19)
- In the split pure BF theory the field $A$ and its cotangent field $B$ have action $\int\langle B,(\bar\partial+d_\tau)A+\tfrac12[A,A]\rangle$. Each cubic vertex supplies one $B$ and each internal $AB$ propagator consumes one. With $V$ vertices and $E$ internal edges, a connected graph has $L=E-V+1$ loops and $N_B^{\rm ext}=V-E=1-L$ external $B$ legs. Since $N_B^{\rm ext}\ge0$, the bare BF vertices generate no connected graph at $L\ge2$. The flagship two-loop obstruction therefore imposes no constraint on the bare pure-SYM twist. A separately added level deformation requires its own twisted vertices. [ESW, §1.6.1, Definition 1.37](https://arxiv.org/pdf/2002.10517#page=28)
- ESW treats the classical theory obtained by reduction from six dimensions. Relating the flagship obstruction to a physical five-dimensional level would require a map from the supersymmetrically completed level deformation to $\hbar J_1$ in the flagship theory. Equation (1.6) gives a conditional matching equation for $SU(N)$, with the map and coefficient still open. [ESW, §9.1.1, Theorem 9.1](https://arxiv.org/pdf/2002.10517#page=66) · [Definition 1.20](https://arxiv.org/pdf/2002.10517#page=23)
- Costello 1610.04144 §1.5 [src]: the 5d theory arises from the twist of the D6-brane theory (7d) on $\mathbb R^3\times\mathbb C^2$ in the $\Omega$-background, gauge algebra $\mathfrak{gl}_{k+1}$; his existence theorem is for $|\delta|\ll|e|$ (§1.3).
- For the flagship theory itself, counterterms outside $V^\eta$ are addressed separately by the $S\ne U$ summands of flagship Thm 3.1 and referee-2 R1.

### Cotangent-field counting and the kinetic pairing

Keep the BF cross-propagator and add the restricted level term (1.4) as a vertex containing only $A$. If a connected graph has $v$ such vertices, each of its other vertices supplies one $B$. Each edge consumes one $B$, so its loop number is $1-v-N_B^{\rm ext}$. A single level insertion therefore occurs only at tree level. Its BV self-bracket is zero. The free BV differential gives a total $(\bar\partial+d_\tau)$ derivative. For the commutator term take the differential with $Q\mathcal A=\mathcal A^2$; changing its overall sign preserves closure. Then $Q\partial\mathcal A=[\mathcal A,\partial\mathcal A]$ in total degree and
$$
Q\operatorname{tr}\big(\mathcal A(\partial\mathcal A)^2\big)
=-\operatorname{tr}\big(\mathcal A^2(\partial\mathcal A)^2\big)
=-\tfrac12\partial\operatorname{tr}\big(\mathcal A^3\partial\mathcal A\big).
$$
The second equality uses the vanishing graded trace of the square of the odd element $\mathcal A\partial\mathcal A$. Compact support kills both total derivatives. Thus adding this specific level vertex to BF gives no two-loop level-matching equation. For any more general completion, two loops with the same propagator require $\sum_{\text{vertices}}(\#B-1)=1+N_B^{\rm ext}$.

There is also an obstruction to obtaining the single-field kinetic pairing by a Lie algebra section of the pure BF target. For a simple algebra, use the Killing form $K$ to write such a section of $\mathfrak g\ltimes\mathfrak g^*\to\mathfrak g$ as $X\mapsto(X,K(DX,-))$. Preservation of the bracket is equivalent to $D[X,Y]=[DX,Y]+[X,DY]$. The pulled-back BF kinetic pairing is

$$
K(DX,Y)+K(X,DY)
=\operatorname{Tr}_{\rm ad}\!\left([D,\operatorname{ad}_X]\operatorname{ad}_Y
+\operatorname{ad}_X[D,\operatorname{ad}_Y]\right)
=\operatorname{Tr}_{\rm ad}[D,\operatorname{ad}_X\operatorname{ad}_Y]=0 .
$$

For a reductive algebra the cocycle equation evaluated on a central argument also kills the mixed derived–central pairing. Only the center can carry a nonzero pullback. For $\mathfrak{sl}_2,\mathfrak{sl}_3,\mathfrak{gl}_2$ the exact dimensions of the space of sections at this linear level are $3,8,4$ and the dimensions of their symmetric images are $0,0,1$. Thus even a field truncation by such a section loses the flagship's kinetic term on the simple factor. This restricts (1.5) independently of any choice of $c_{\rm tw}$. [Coadjoint cocycle replay](cotangent_sections.py)

## 3. The $e\to0$ limit of twisted M-theory [drv]

Costello 1610.04144 §1.2 [src]: $S=\frac1{2\delta}\int\mathrm{Tr}\,A\,dA\,dz_1dz_2+\frac1{3\delta}\int\mathrm{Tr}\,A*_eA*_eA\,dz_1dz_2$; first Moyal order $=$ $j_M=D/2$ (flagship Prop. 5.8). Flagship $\mathfrak{gl}_N$ family $j=s\,j_M+t\,\pi^*D_{\mathfrak{sl}}$, $e=s\hbar$:

$$t(s+t)=-12N^2\varepsilon c_O=\frac{\varepsilon U N^2}{256\pi^4}\quad\Rightarrow\quad t=\frac{\varepsilon UN^2}{256\pi^4\,s}+O(s^{-3}),\qquad \hbar\,t\,\pi^*D_{\mathfrak{sl}}=\frac{\varepsilon UN^2}{256\pi^4}\,\frac{\hbar^2}{e}\,\pi^*D_{\mathfrak{sl}}+\dots\tag{3.1}$$

$e^{-1}\delta^2$ is inside the regime $e^{-r}\delta^n$, $r\le n$, of Costello Thms 9.0.2–9.0.3. At $e=0$ the obstruction $[O_2]\neq0$ is not cancelled by $j_M$ ($\tfrac12\{j_M,j_M\}=d\,i_M$). In the equivariant sector the possible order-$e^{-2}\delta^2$ commutator trace extends to the closed direction $e^{-2}s_e$. Modulo this direction and $e^{-1}\partial_e s_e$, (3.1) is the leading coefficient for zero first-order quantum deformation in $\delta$. [Laurent cochain calculation](../../deepen/moyal/README.md)

## 4. The Vogel factor [drv from verified identities]

$n_4=d(C_4-\tfrac18)$, $C_4=\frac{3t\,t_2-t_3}{16t^3}$ (Okubo quartic index, MSV 1105.0115 (8)), $\Lambda^2\mathfrak g=\mathfrak g\oplus X_2$, $\dim X_2=\tfrac{d(d-3)}2$:

$$8n_4-3d=4d\,(2C_4-1)=\tfrac32\,\dim X_2\,\frac{\alpha\beta\gamma}{t^3},\qquad \mu(F_{223})=2d(8n_4-3d)=3d\,\dim X_2\,\frac{\alpha\beta\gamma}{t^3}.\tag{4.1}$$

For $\mathfrak{sl}_2$, $\dim X_2=0$ and $C_4=\tfrac12$. In the Killing pairing let $P_{X_2}$ project onto the kernel of the Lie bracket in $\Lambda^2\mathfrak g$, and let $F_{223}^{uv}$ be the $uv$ coefficient in the paper's four-input color tensor. In the $(01)|(23)$ channel, crossing the two middle slots of an alternating cochain gives the trace identity for the detector $\mu$ of (4.1), written $\varpi$ in the flagship draft:
$$
\mu(F_{223})=-2d\,\operatorname{tr}(P_{X_2}F_{223}^{uv}),\qquad
\operatorname{tr}(P_{X_2}F_{223}^{uv})=3d-8n_4
=-\frac{3\alpha\beta\gamma}{2t^3}\dim X_2 .
$$
The singlet trace of the same tensor equals $-6$. An identification with an untwisted effective coupling remains open. [Projector, crossing calculation and 22 contraction tables](x2_channel_check.py)
