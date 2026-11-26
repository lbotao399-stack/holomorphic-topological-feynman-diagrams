# A classical bulk source residue with full flux closure

26 November 2026.

This computation supplies an independent classical bulk control for the
[quantum-supercharge reconstruction](README.md). It does not reproduce a
positive-order boundary supercharge coefficient.

## 1. A nonzero source/normalizable-state residue

Use the repository's tangent coordinates $z_1,z_2$, normal coordinates $w_1,w_2,w_3$, and $P=\partial_{z_1}\wedge\partial_{z_2}$. Work in the three-chart Čech model of the punctured normal space, with $U_i=\{w_i\ne0\}$ and the convention

$$
(\delta K)_{123}=K_{23}-K_{13}+K_{12}.
$$

The flux class is the already-used representative

$$
F=\frac{P}{w_1w_2w_3}.
$$

The normalizable mode in the protected $\operatorname{Tr}(\beta_1\gamma^2)$ sector is

$$
P\wedge\left(
\frac{2\partial_{w_3}}{w_1w_2^3w_3}
-\frac{\partial_{w_2}}{w_1w_2^2w_3^2}
\right).
$$

Take the regular, divergence-free vector fields

$$
A=w_2\partial_{w_1},\qquad
B=w_2\partial_{w_2}-w_3\partial_{w_3}.
$$

Their flux brackets are obtained directly, with the Schouten convention $\{A,F\}=L_AF$:

$$
\begin{aligned}
\{F,A\}
&=-w_2\partial_{w_1}\left(\frac1{w_1w_2w_3}\right)P
=\frac{P}{w_1^2w_3},\\
\{F,B\}
&=-\left(w_2\partial_{w_2}-w_3\partial_{w_3}\right)
\left(\frac1{w_1w_2w_3}\right)P=0.
\end{aligned}
$$

The first bracket is Čech-exact. Set

$$
K_{13}=\frac{P}{w_1^2w_3},\qquad K_{12}=K_{23}=0.
$$

Then

$$
\delta K=-\{F,A\},\qquad
\partial_\Omega A=\partial_\Omega B=\partial_\Omega K=0,
\qquad \{F,K\}=0.
$$

Consequently the full linear differential

$$
D=\delta+t\partial_\Omega+N\{F,-\}
$$

satisfies

$$
D(A+NK)=0,\qquad DB=0.
$$

This construction uses only the already-established linear Čech flux model. It does not assume a complete source/operator pairing or a quantum holographic map.

The exterior product of the normalizable mode with $A$ and $B$, ordered as written, is

$$
\begin{aligned}
&P\wedge\left(
\frac{2\partial_{w_3}}{w_1w_2^3w_3}
-\frac{\partial_{w_2}}{w_1w_2^2w_3^2}\right)
\wedge(w_2\partial_{w_1})
\wedge(w_2\partial_{w_2}-w_3\partial_{w_3})\\
&\quad=
\left(\frac{2}{w_1w_2w_3}-\frac{1}{w_1w_2w_3}\right)
P\wedge\partial_{w_1}\wedge\partial_{w_2}\wedge\partial_{w_3}\\
&\quad=
\frac{P\wedge\partial_{w_1}\wedge\partial_{w_2}\wedge\partial_{w_3}}
{w_1w_2w_3}.
\end{aligned}
$$

Replacing $A$ by its full flux-closed lift $A+NK$ leaves this expression unchanged: the additional term contains two copies of $P$ and vanishes.

Contract against the holomorphic five-volume, and take the transverse Grothendieck residue, normalized by the positively oriented product circles. Its value is exactly

$$
\boxed{
\frac1{(2\pi i)^3}
\oint\!\oint\!\oint
\frac{dw_1\,dw_2\,dw_3}{w_1w_2w_3}=1.
}
$$

Equivalently, for arbitrary regular normal vectors $A=\sum A^i\partial_{w_i}$ and $B=\sum B^i\partial_{w_i}$, the top exterior coefficient is

$$
\frac{2(A^1B^2-A^2B^1)}{w_1w_2^3w_3}
+\frac{A^1B^3-A^3B^1}{w_1w_2^2w_3^2}.
$$

Its transverse residue is

$$
\left.
\left[
\partial_{w_2}^2(A^1B^2-A^2B^1)
+\partial_{w_2}\partial_{w_3}(A^1B^3-A^3B^1)
\right]
\right|_{w=0}.
$$

For the specified pair this is $2-1=1$. The sign and coefficient were checked independently with exact symbolic differentiation.

### A relative normalization fixed by the Cartan weights

More generally replace $B$ by the diagonal divergence-free source

$$
H=h_1w_1\partial_{w_1}+h_2w_2\partial_{w_2}+h_3w_3\partial_{w_3},
\qquad h_1+h_2+h_3=0.
$$

It is strictly flux-closed. The same product has residue

$$
2h_2+h_3=h_2-h_1,
$$

and the ordinary vector-field bracket is

$$
[H,A]=(h_2-h_1)A.
$$

Thus the residue is exactly the Cartan root weight. For the two fixed integral Cartan generators

$$
w_2\partial_{w_2}-w_3\partial_{w_3},\qquad
w_1\partial_{w_1}-w_2\partial_{w_2},
$$

the residues are respectively $1$ and $-2$. Their ratio is $-2$ and does not change under a common rescaling of the normalizable mode or of $A$. This is a classical relative normalization check. It does not fix the overall bulk-to-boundary normalization needed for a quantum $Q_n$ coefficient.

The normalizable mode has the opposite Cartan weight, $h_1-h_2$. Accordingly this calculation has the natural interpretation of a classical $\mathfrak{sl}_3$ source/current weight coupling. It detects no positive-order $Q_n$ by itself. An independent agent audited the Čech lift, exterior signs, root-weight formula, and this scope statement and found them consistent.

The same proof works for every ordered pair of distinct normal directions. For a permutation $(i,j,k)$ of $(1,2,3)$, multiply the state

$$
P\wedge\left(
\frac{2\partial_{w_k}}{w_iw_j^3w_k}
-\frac{\partial_{w_j}}{w_iw_j^2w_k^2}
\right)
$$

by the sign of that permutation, and pair it with $w_j\partial_{w_i}$ and $H$. The residue is $h_j-h_i$. Its source flux bracket has coefficient $1/(w_i^2w_k)$; the Čech primitive lies on $U_i\cap U_k$, with sign chosen so its coboundary is the negative flux bracket. No further correction is needed: its divergence and its bracket with $F$ are both zero.

The standalone standard-library verifier `source_residue.py` checks the full five-direction exterior products, divergence and Schouten brackets, all six Čech source lifts, and 486 exact Cartan residues. Its output is `source_residue.json` (PASS). The finite checks validate the algebra and signs; the displayed formulas prove the unrestricted Cartan identity.

### Scope of this computation

This is a nonzero cyclic residue of one normalizable mode and two flux-closed regular-source lifts in the linear bulk model. It is a possible transverse factor in a bulk boundary contact calculation. It is **not** a full integrated fivefold amplitude: the tangent Dolbeault degrees, radial boundary action, external kernels, and quantum composite renormalization have not been supplied. In particular, it is not equal by derivation to any coefficient of $Q_1$ or $Q_2$.

It establishes a concrete distinction relevant to the user's question: diagrams with normalizable representatives on every leg vanish for an elementary exterior-algebra reason, while mixed source/state couplings need not vanish. We must construct the latter rather than interpreting the former as a failure of holography.


## Sources and reproduction

The linear flux complex and Laurent dictionary are those of the existing
[bulk state calculation](../NOTE.md),
with the normalizable mode from the
[all-order input-state calculation](../boundary_bulk_all_orders/bulk/REPORT.md).
The punctured fivefold background is proposed by Costello–Li,
[*Twisted supergravity and its quantization*, Section 9.1](https://arxiv.org/html/1606.00365v1).

From the repository root:

```bash
python explore/twisted-holography-N4/independent_bulk_reconstruction/source_residue.py
```

The proof fixes the general Cartan coefficient. The independent exact
calculation verifies the full exterior products, source completions, and
486 numerical Cartan choices over the rational numbers.
