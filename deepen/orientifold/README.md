# M4. Five-dimensional Moyal brackets and the orientifold

### The coincident D6/O6 projection

Sen's coincident D6/O6 system has orientifold operation $(-1)^{F_L}\Omega I_3$, with $I_3$ reflecting the three transverse coordinates. Its orthogonal M-theory realization identifies $(\vec r,x^4)$ with $(-\vec r,-x^4)$ at the asymptotic circle. Costello's first equivariant parameter $e$ multiplies the rotation $V$ of this circle, and his three-form background is $-e V^\flat\wedge\omega^{0,2}_{\mathbb C^2}$. The generator transforms to $-V$ under the identification while the D6 tangent plane $\mathbb C^2$ is fixed. The specified background therefore requires $e=-e$, hence $e=0$, at a fixed complex parameter. In the IIA description the tangential NS two-form is odd under $\Omega$ and even under the transverse reflection. This excludes a smooth nonzero constant tangential field for the coincident-plane projection. [Sen, equations (12), (15)](https://arxiv.org/pdf/hep-th/9707123) · [Costello, Proposition 6.0.2 and §§7–8](https://arxiv.org/pdf/1610.04144) · [Bonora et al., §3](https://arxiv.org/pdf/hep-th/0006091)

### Rigidity on the ordinary gauge fields

Use the manuscript's determinant $[pq]$ of two holomorphic momenta, with the defining trace pairing. For a matrix $X$ write $X^\sharp=X^{\mathsf T}$ in the orthogonal case, or $X^\sharp=J^{-1}X^{\mathsf T}J$ for the defining symplectic form $J$. The ordinary Lie algebra has $X^\sharp=-X$. The Moyal commutator of two such waves has matrix coefficient
$$
\cosh(e[pq]/2)[X,Y]+\sinh(e[pq]/2)(XY+YX).
$$
The second summand is self-adjoint for $\sharp$. Projecting it out leaves the first summand.

For any simple algebra of class (B) in the manuscript, the only invariant three-tensor is the commutator tensor. Thus a translation- and $SL_2$-invariant full-wedge cubic on these fields gives $K([pq])[X,Y]$ on waves, where $K$ is even by input exchange. Polynomiality at each perturbative order and the two-momentum invariant ring give this form to every order. Take $p=(s,0)$, $q=(0,1)$, $r=(s,1)$ and colors $X,Y,X$. The Jacobi sum is
$$
K(s)\big(K(0)-K(2s)\big)[[X,Y],X].
$$
Take $X,Y$ in a root $\mathfrak{sl}_2$ with $[[X,Y],X]\ne0$. The classical bracket has invertible constant coefficient, so Jacobi forces $K(2s)=K(0)$ identically; linearization forces the same conclusion for a closed first-order cubic deformation. For the projected Moyal bracket the factor is $\cosh(es/2)(1-\cosh(es))$, whose $e^2$ coefficient is $-s^2/2$.

Over $\mathbb C((e))$ even a Laurent choice therefore only rescales the classical cubic and its first-order closed perturbation in this sector. A nonzero rescaling multiplies the two-loop graph by its sixth power and the classical differential by its first power. The bracket of two multiples of the original cubic vanishes. The manuscript's nonzero class (B) two-loop obstruction survives. This argument applies at the stated ranks, including $\mathfrak{so}_5,\mathfrak{sp}_4$. The manuscript treats low-rank isomorphisms and the cubic invariant of $\mathfrak{so}_6=\mathfrak{sl}_4$ separately. [Class (B) and the two-loop obstruction](../../notes/flagship/flagship.tex)

### The family with $e$ reversed

Bonora et al. impose $A(e)=-A(-e)^\sharp$ and then discuss imposing a Seiberg–Witten map to express higher components in terms of the ordinary gauge potential. Work first with the unconstrained higher coefficients of their parity condition. The even powers of $e$ have skew-adjoint coefficients and the odd powers self-adjoint coefficients. Put $u=e^2$. On the basis consisting of skew-adjoint $X$ and $eY$ for self-adjoint $Y$, the trace metric has blocks $\operatorname{tr}(XX')$ and $u\operatorname{tr}(YY')$ and vanishing cross terms. The second inverse block is $u^{-1}$ times the inverse self-adjoint metric.

Over $\mathbb C((u))$ extension to $\mathbb C((e))$ identifies the fields with all matrices via this basis. The parity action respects the Moyal bracket, because transpose reverses matrix multiplication and $e\mapsto-e$ reverses the scalar Moyal product. The quantum cubic at power $e^n$ has momentum exponent $n+2$; its color tensor has charge conjugation parity $(-1)^{n+2}$. It therefore descends termwise to the invariant family. Extending scalars commutes with every propagator contraction and the CE differential, and the field extension is faithful. The [full two-loop matrix counterterm](../moyal/README.md) solves the quantum equation on these localized invariant fields. Its $e^{-1}$ numerator vanishes on three skew-adjoint inputs. With one $eY$ input the pole has finite value, and exchange of the self-adjoint field contributes at net power $e^0$.

At $u=0$ the self-adjoint metric block degenerates. In the reduction of the lattice modulo $u$, the generators $eY$ still occur, whereas evaluation of the matrices at $e=0$ discards them. These operations give different classical field content. The graph cancellation just obtained uses independent $eY$ fields. A Seiberg–Witten-constrained theory has its own BV action and propagator. The string construction in Bonora et al. uses a sign-changing transverse $B$ profile. Bars, Sheikh-Jabbari and Vasiliev instead combine matrix transpose with a reflection of a noncommuting coordinate, keeping the deformation parameter fixed and changing the tangent projection. Their fixed-point fields obey parity conditions on spacetime jets. An application to a specific twisted intersecting-brane system requires its projected BV pairing and fixed-locus interactions. [Bonora et al., equations (2.4)–(2.13) and §3](https://arxiv.org/pdf/hep-th/0006091) · [Bars et al., equations (21)–(22), (52)–(54)](https://arxiv.org/pdf/hep-th/0103209)

### Exact contact on commuting fields

Use momenta $(1,0),(0,1),(1,0),(-2,-1)$ and colors $X,X,Y,Y$. In $\mathfrak{so}_5$ take $X=E_{12}-E_{21}$, $Y=E_{34}-E_{43}$. In $\mathfrak{sp}_4$ take $X=\operatorname{diag}(1,0,-1,0)$, $Y=\operatorname{diag}(0,1,0,-1)$. Here $E_{ij}$ are defining matrix units. In both cases $XY=0$ and $\operatorname{tr}(X^2)\operatorname{tr}(Y^2)=4$. Every classical CE boundary vanishes on these waves. The [replay](check.py) contracts the six commutator tensors with the seven inverse trace pairings and alternates all four slots. Values in units of $\pi\nu$ are

| external colors | internal propagator | $(1,3,3)$ | $(2,2,3)$ | sum |
|---|---|---:|---:|---:|
| $\mathfrak{so}_5$ | skew-adjoint | $3/5$ | $2/5$ | $1$ |
| $\mathfrak{so}_5$ | full matrix | $12$ | $8$ | $20$ |
| $\mathfrak{sp}_4$ | skew-adjoint | $24/5$ | $16/5$ | $8$ |
| $\mathfrak{sp}_4$ | full matrix | $48/5$ | $32/5$ | $16$ |

With $N$ the matrix size, the exact pole bracket of $j_M$ with $\pi^*D_{\mathfrak{sl}}$ is $16/N$. The gradients of both cubics on two skew-adjoint inputs are self-adjoint, so only the self-adjoint inverse block contributes. The graph and tree contact in the descended theory give $4N-(N^2/4)(16/N)=0$. The nonzero first and third rows use only the ordinary gauge fields.
