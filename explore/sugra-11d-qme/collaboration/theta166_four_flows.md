# The four oriented contractions at fixed external labels on $\theta_{166}$

Fix the order-$g$ sector with $\nu=0$. Use the momenta and polarizations of the [single-flow calculation](theta166_flow.md). The ten external fields occupy prescribed vertices, with $\gamma$ at $B_1$ and $\mu$ at each other interior vertex. This is an off-shell evaluation. In particular the $\mu$ at $C_1$ has momentum $-e^1$ and polarization $e_1$.

All four admissible orientations share the integrand factor

$$
(1+(\ell_B)_1)^3(1-(\ell_C)_1)^4
\det\begin{pmatrix}1&(\ell_B)_3&(\ell_C)_3\\1&(\ell_B)_4&(\ell_C)_4\\1&(\ell_B)_5&(\ell_C)_5\end{pmatrix}.
$$

In the table $x=(\ell_B)_1$ and $z=(\ell_C)_1$ are the first components of the two path currents. Each row is specified by its field orientation on the single edge $A$ and the position of the $\gamma\partial\gamma\partial\gamma$ vertex. The last column is $44464437108000$ times the signed Dirichlet mean of the full numerator.

| $A$ | Cubic $\gamma$ vertex | Remaining signed factor | Scaled mean |
|---|---|---:|---:|
| $N\to N'$ | $B_1$ | $-2x^2-4xz+x-4z^2+3z$ | $-65558839$ |
| $N\to N'$ | $N'$ | $-(x+2z)(2x+1)$ | $-11982421$ |
| $N'\to N$ | $B_1$ | $2x^2+4xz+x+4z^2-z$ | $85509601$ |
| $N'\to N$ | $N$ | $(2x+1)(x+2z-1)$ | $-27353656$ |

The two rows with a cubic $\gamma$ vertex at $B_1$ sum to $2x+2z$. The other two rows sum to $-2x-1$. The complete signed factor reduces to $2z-1$. Since $D_{166}=1$, the ordered contraction, including the source theta measure with $g$ removed, is

$$
\frac{\pi}{158400}\frac{-65558839-11982421+85509601-27353656}{44464437108000}
=-\frac{3877063\pi}{1408633367581440000}.
$$

Putting the five negative $C$ momenta in the same order as the $B$ momenta makes all four integrals vanish. Antisymmetrization over external odd fields, the other four theta topologies and the physical normalization enter the bulk obstruction separately. [Exact four-flow check](check_theta166_four_flows.py)

## Exact cancellation after antisymmetrizing the nine external $\mu$ fields

A coordinate permutation acts twice on the holomorphic five-volume determinants of $D_{166}$ and once on the volume multivector of the unique cubic $\gamma$ vertex. All other tensor contractions pair a vector with a covector. Thus an orthogonal signed coordinate permutation $T$ multiplies the kernel with all external data transformed by $(\det T)^3$. This calculation uses the explicit tensors in the fixed volume convention.

For the polarizations above, interchange axes $3$ and $4$. The external $\gamma$ momentum $e^1$, its polarization $e^2$, and all nine $\mu$ polarizations $e_1$ are fixed. The induced permutation on the nine odd legs is

$$
(+e^3\;\; +e^4)(-e^3\;\; -e^4),
\qquad
\operatorname{sgn}=+1,\qquad(\det T)^3=-1 .
$$

The antisymmetrized kernel therefore equals its negative and vanishes. This pairs terms even with $\gamma$ kept at $B_1$.

For the transverse assignment, give the $\mu$ at $B_2$ (momentum $+e^2$) polarization $e_3$ and the $\mu$ at $C_1$ (momentum $-e^1$) polarization $e_2$. The other seven $\mu$'s retain $e_1$. The ordinary coordinate permutations preserving this data fix axes $1,2,3$, leaving precisely the identity and $(4\;5)$. The latter exchanges $(+e^4\;\; +e^5)(-e^4\;\; -e^5)$ and leaves every polarization fixed, so this antisymmetrized kernel also vanishes.

More explicitly, write a signed permutation as $Te^i=s_i e^{\pi(i)}$, with $s_i=\pm1$ and $\pi$ an axis permutation. For the transverse fixture, preservation of the leg data up to multilinear polarization signs requires $\pi\in\{1,(4\;5)\}$, $s_1=s_2=1$, and allows every choice of $s_3,s_4,s_5$. Each sign reversal exchanges the two odd legs with opposite momentum on that axis. The lone polarization $e_3$ at $B_2$ supplies an additional $s_3$. The three signs are

$$
\operatorname{sgn}(\text{odd legs})=s_3s_4s_5,\qquad
\operatorname{sgn}(\text{polarizations})=s_3,\qquad
(\det T)^3=\operatorname{sgn}(\pi)s_3s_4s_5 .
$$

Eight of these sixteen transformations force zero, precisely those with $\operatorname{sgn}(\pi)s_3=-1$. One has determinant $+1$, namely $\operatorname{diag}(1,1,-1,-1,1)$. For both fixtures, the nonzero fixed-label values are compatible with exact cancellation of the alternating component evaluated on those data.
