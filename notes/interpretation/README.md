# Five-dimensional Chern–Simons theory

## I2. Trace content of the two-loop obstruction

Use $\kappa(X,Y)=\operatorname{tr}(XY)$ in the defining representation. For the four external matrices $X_i$, $p_i$ are holomorphic momenta with $\sum_i p_i=0$, $[p_ip_j]$ is their two-dimensional determinant, and $\operatorname{alt}$ is the signed sum over all 24 input permutations. In the manuscript $\beta=d([q_0q_1]^2 f)$ is the Chevalley–Eilenberg boundary of the three-cochain with $f(X,Y,Z)=\operatorname{tr}([X,Y]Z)$. [Manuscript, section on the classes for $\mathfrak{sl}_N$ and $\mathfrak{gl}_N$](../flagship/flagship.tex)

For $\mathfrak{sl}_N$, both graph cochains satisfy the following identities in the formal trace algebra:
$$
\begin{aligned}
F_{133}+F_{223}&=2N^2\beta,\\
F_{223}&=-24\beta-6N\operatorname{alt}\big([p_0p_1][p_2p_3]\operatorname{tr}(X_0X_1)\operatorname{tr}(X_2X_3)\big).
\end{aligned}
$$

For $\mathfrak{gl}_N$, every external leg enters a commutator vertex, so replacing $X_i$ by $X_i-\operatorname{tr}(X_i)\mathbf1/N$ preserves the graph. This sends each two-letter trace to $\operatorname{tr}(X_iX_j)-\operatorname{tr}(X_i)\operatorname{tr}(X_j)/N$. The term with four one-letter traces vanishes by $[p_0p_1][p_2p_3]-[p_0p_2][p_1p_3]+[p_0p_3][p_1p_2]=0$. The two triple-trace terms coincide after alternation. Hence the additional term in $F_{223}$ is
$$
12\operatorname{alt}\big([p_0p_1][p_2p_3]\operatorname{tr}(X_0X_1)\operatorname{tr}(X_2)\operatorname{tr}(X_3)\big).
$$

In $[O_2]=-U\hbar^2[F_{223}]/(3072\pi^4)$, $U$ is the common orientation factor. With $\hbar=\lambda/N$, the coefficients of the two trace structures above, in the order double trace and triple trace, are
$$
\frac{U\lambda^2}{512\pi^4 N},\qquad -\frac{U\lambda^2}{256\pi^4 N^2}.
$$
The triple trace vanishes on traceless inputs. These are coefficients of unnormalized traces. The manuscript proves the resulting class is nonzero for $N\ge3$ and zero for $N\le2$ where the gauge algebra is defined.

Both graphs have six trivalent vertices and seven internal edges. A ribbon contraction of genus $h$, with $b$ trace factors carrying external matrices and $f$ closed index loops, obeys $6-7+f=2-2h-b$. Each closed index loop contributes $N$. The single-trace terms in $(-3F_{133}+2F_{223})/240$ have coefficient $-(N^2/40+1/2)\beta$ and are exact. The surviving $N$ double trace and constant triple trace are both planar, with two and three external boundaries.

The verifier contracts both graphs and compares each formal trace and momentum coefficient, checks closure of the displayed multi-trace cochains, and tests the ribbon powers. [Exact replay](planar_trace_check.py)

For the five-dimensional theory with traceless matrices, the double-trace channel exhausts the higher four-point cohomology. Its degree-$2w$ part consists of binary forms of degree $2w-4$, with an equivariant obstruction only at $w=2$. [Proof, all-order local lift and finite-rank replays](../../deepen/o3/higher_orders/README.md) · [Matrix center](../../deepen/o3/higher_orders/center_classes.md)
