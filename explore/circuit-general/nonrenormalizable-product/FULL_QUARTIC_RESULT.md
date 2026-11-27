# Complete quartic source polynomial

This is the complete three-loop current polynomial for the anomaly-free abelian theory with charges $(1,-1,0,0,0)$, interaction $\bar W=\gamma_1\gamma_2\gamma_3\gamma_4$, and inputs $A=\beta^1(\partial_{+\dot a}\beta^2)\gamma_5$, $B=b\beta^3$. The fifth field is neutral and absent from the interaction. The sources $\lambda_{\gamma_4},\lambda_c,\lambda_{\gamma_5}$ act on the indicated surviving fields. All contractions and all period coefficients have already been summed.

The raw `htk` current expression, before the common physical normalization, is

\[
\frac1{4320}\Bigg\{
\lambda_{\gamma_4,+\dot a}\Bigg(
\begin{aligned}
&12\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}\\
&+16\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle\\
&+36\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-20\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\\
&+96\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-12\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-24\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}\\
&+54\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-18\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}
\Bigg)
+
\lambda_{c,+\dot a}\Bigg(
\begin{aligned}
&-19\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}\\
&-28\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle\\
&-46\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-18\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\\
&-20\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-48\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+6\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}\\
&-48\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&+36\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-24\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}
\Bigg)
+
\lambda_{\gamma_5,+\dot a}\Bigg(
\begin{aligned}
&12\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}\\
&-36\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&+18\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}
\Bigg)
\Bigg\}.
\]

For the coefficient of $\hbar^3$ in the physical product, multiply this expression by $-g^2/(8\pi^2)^3$ and replace each $\lambda$ by $-\partial_+/2$ acting only on its indicated surviving field. Equivalently, the displayed integer polynomial in physical derivatives has common coefficient $g^2/[552960(8\pi^2)^3]$ and acts on $\gamma_4c\gamma_5$.

Every term contains either a derivative source on $c$ in a bracket or the free source $\lambda_{c,+\dot a}$. The complete expression therefore vanishes at $\lambda_c=0$.

Setting $\gamma_5=1$ gives the complete compact raw expression

\[
\frac{\langle\lambda_{\gamma_4},\lambda_c\rangle^3}{4320}
\left(12\lambda_{\gamma_4,+\dot a}-19\lambda_{c,+\dot a}\right).
\]

The exact all-source reconstruction and every Wick matching are in [quartic_wick_certificate.json](quartic_wick_certificate.json); [quartic_wick.py](quartic_wick.py) recomputes them. This is a complete specified cochain product coefficient, with no claim that the inputs are quantum cocycles or that every nonrenormalizable product is rational.
