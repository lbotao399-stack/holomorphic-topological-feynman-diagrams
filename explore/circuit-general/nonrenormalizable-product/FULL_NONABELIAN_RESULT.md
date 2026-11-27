# Complete nonabelian quartic source polynomial

Take four adjoint matter flavors of $SU(2)$, the invariant adjoint metric, and
\[
\bar W=(\gamma_1^D\gamma_{3,D})(\gamma_2^E\gamma_{4,E}),\qquad
A=\beta^{1D}(\partial_{+\dot a}\beta_D^2)\gamma_5,\qquad
B=b_E\beta^{3E}.
\]
The fifth field is neutral, free, and absent from the interaction. These inputs are color-singlet cochains, with no assumption of quantum closure. The sources $\lambda_{\gamma_4},\lambda_c,\lambda_{\gamma_5}$ act on the indicated surviving fields. The complete output color tensor contracts $\gamma_4^D c_D$.

## Complete three-loop result

The raw `htk` current expression, with all twenty-eight matchings and all their color coefficients summed, is

\[
\frac1{2160}\Bigg\{
\lambda_{\gamma_4,+\dot a}\Bigg(
\begin{aligned}
&24\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}\\
&-40\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle\\
&+72\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&+56\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\\
&-33\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-48\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}\\
&+126\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&+12\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}
\Bigg)
+
\lambda_{c,+\dot a}\Bigg(
\begin{aligned}
&15\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}\\
&-68\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle\\
&+25\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&+9\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\\
&+8\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-75\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+3\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}\\
&-36\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&+54\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+51\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}
\Bigg)
+
\lambda_{\gamma_5,+\dot a}\Bigg(
\begin{aligned}
&6\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}\\
&-72\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&+144\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}
\Bigg)
\Bigg\}.
\]

For the coefficient of $\hbar^3$ in the physical product, multiply this expression by $g^2/(8\pi^2)^3$ and replace each source by $-\partial_+/2$ acting only on its indicated surviving field. The integer polynomial in physical derivatives therefore has common coefficient $-g^2/[276480(8\pi^2)^3]$ and acts on $\gamma_4^D c_D\gamma_5$.

Every term contains either a ghost source in a bracket or the free spinor $\lambda_{c,+\dot a}$, so the expression vanishes at $\lambda_c=0$. Setting $\gamma_5=1$ gives
\[
\frac{\langle\lambda_{\gamma_4},\lambda_c\rangle^3}{720}
\left(8\lambda_{\gamma_4,+\dot a}+5\lambda_{c,+\dot a}\right).
\]

## Complete cancellation of the period

The twenty matchings with two matter vertices contribute the following entire coefficient multiplying $\zeta(3)$:

\[
\frac1{1050}\Bigg\{
\lambda_{\gamma_4,+\dot a}\Bigg(
\begin{aligned}
&6\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle\\
&-3\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&+2\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\\
&+4\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&+3\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+2\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}\\
&-2\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&+6\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}
\Bigg)
+
\lambda_{c,+\dot a}\Bigg(
\begin{aligned}
&3\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle\\
&-6\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-3\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\\
&-4\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-2\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-6\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&+2\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-2\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}
\Bigg)
\Bigg\}.
\]

The eight matchings with one matter vertex and one gauge vertex contribute exactly the negative of this displayed polynomial. Therefore the complete coefficient of $\zeta(3)$ is zero for every retained source spinor, including dependent spinors. No relation or independence assumption among different zeta values is used. With $\gamma_5=1$, each displayed intermediate period contribution is already zero.

The full rational and period contributions of both classes, every matching and color contraction, the exact interpolation proof, and five independent complete-source replays are recorded in [nonabelian_wick_certificate.json](nonabelian_wick_certificate.json). [nonabelian_wick.py](nonabelian_wick.py) recomputes the entire deterministic certificate. Normalization and the proof that all higher-loop coefficients vanish for these particular inputs are in [README.md](README.md).
