# Four-flavor test for the unit-singleton projection

Costello's Proposition 16.5.2 fixes the first-order generator on two linear Weyl symbols $z_1,z_2$, with $z_1*z_2-z_2*z_1=c$. The unit-output test of a field $V$ against $E_{ij}$ is the constant Weyl coefficient of $\mathrm{Tr}(E_{ij}V)$. On a symmetric product of two fields take the sum of the two pairings. Use four distinct flavor coordinates.

| Test | First input | Second input | Dual output matrices |
|---|---|---|---|
| $A$ | $E_{12}z_1$ | $E_{23}z_2$ | $E_{34},E_{41}$ |
| $B$ | $E_{12}z_2$ | $E_{23}z_1$ | $E_{34},E_{41}$ |
| $C$ | $E_{12}z_1$ | $E_{34}z_2$ | $E_{23},E_{41}$ |

Each test selects one length-four color cycle. A primitive for this component may be restricted to the invariant single-trace one-cochains with two PBW outputs. Such a primitive on an ordinary marked matrix has the form
$$
H(E_{ij}u)=\ell(u)\sum_{\substack{k\in\{1,\ldots,N\}\\\mathrm{flavor\ basis}}}
(E_{ik}f)(E_{kj}g)
$$
and linear combinations thereof. Here $\ell$ extracts a Weyl coefficient and $f,g$ are Weyl monomials. Interchanging $f,g$ allows the other ordering. In a relative cochain $H(E_{13}1)=0$. For a hypothetical absolute primitive of the displayed deformation one can average under the compact color group. Equivariance then gives $dH(X,Y)=H([X,Y])$ for two constant traceless matrices, and the first-order deformation vanishes on those inputs. Taking $X=E_{12},Y=E_{23}$ gives the same normalization.

The input-bracket terms in $(A-B)(dH)$ containing the quadratic symbol cancel. Their residual argument is $cE_{13}1$. For $C$, the two source matrices have zero product in both orders. In the module-action terms the spectator output must be constant, while the other output must be linear to give a constant after multiplication by a linear input. The only surviving coefficient pairs are
$$
\big((A-B)(dH),C(dH)\big)=(c/2,c/2)
$$
for $(\ell,f,g)=(z_1^*,1,z_1),(z_1^*,z_1,1),(z_2^*,1,z_2),(z_2^*,z_2,1)$. Here $z_i^*$ selects the coefficient of $z_i$. This proves $(A-B-C)(dH)=0$. The first-order bracket itself gives
$$
(A,B,C)(\mu_1)=(-1,1,1),\qquad (A-B-C)(\mu_1)=-3 .
$$
On the stable cyclic complex, tensor this functional with the zero-input, one-output functional extracting the constant Weyl coefficient in a second trace. Its differential vanishes because the CE differential raises input number on each trace factor separately. A fifth distinct flavor $E_{55}$ selects that singleton. The resulting test annihilates coboundaries and evaluates to $-3$ on $\kappa\mu_1$. [Matrix-unit calculation](first_order_cycle.py) · [Costello, §§16.5 and A](https://arxiv.org/abs/1705.02500)

All three terms in this detector have total input Weyl degree two and degree-zero PBW outputs, including the singleton. Both $c$ and $\hbar$ have rotation weight two. At second order the only compatible coefficient is a multiple of $\hbar^2/c$. Thus contributions built solely from classical Moyal vertices, $c$-independent free homotopies and Wick contractions vanish on this detector by rotation. The possible Laurent bulk counterterms and first-order comparison maps are the remaining cases treated in the [principal-symbol argument](README.md).
