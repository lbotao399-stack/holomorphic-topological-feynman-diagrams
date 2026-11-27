# Complete open-color, full-source three-loop tensor

Generated exactly by `render_tensor.py` from the full tensor certificate's computation.

The coefficient tensor of the polarized superpotential is defined by

\[\bar W(\gamma)=\bar W^{i_1i_2i_3i_4}\gamma_{1i_1}\gamma_{2i_2}\gamma_{3i_3}\gamma_{4i_4}.\]

The inputs are \(\beta^{1i_1}(\partial_{+\dot a}\beta^{2i_2})\gamma_5\) and \(b_A\beta^{3i_3}\). All six input/output color indices and the free spinor index remain open. The field \(\gamma_5\) is a separate free gauge singlet, absent from \(\bar W\).

Sources and currents have the manuscript's lower index \(+\dot a\). The angle bracket is its antisymmetric spinor pairing, with

\[\langle\lambda_{\gamma_4},\lambda_c\rangle=\epsilon^{\dot b\dot c}\lambda_{\gamma_4,+\dot b}\lambda_{c,+\dot c},\qquad\epsilon^{\dot+\dot-}=1=-\epsilon^{\dot-\dot+}.\]

The same pairing defines the brackets involving the other original source spinor.

Each of the ten color tensors below multiplies its displayed complete source spinor polynomial.
Their sum acts on \(\gamma_{4i_4}\), \(c^B\), and \(\gamma_5\). Replacing each source by the derivative on its own surviving field and multiplying by

\[\frac{g^2}{128(8\pi^2)^3}\]

gives the complete \(g^2\bar W\) coefficient of \(\hbar^3\). The factor from \(\lambda=-\partial/2\) is already included in this multiplier; it must not be included a second time. See [DERIVATION.md](DERIVATION.md) for all signs, normalization, and the exact source-degree proof.

## The complete rational tensor

### 1

\[(t_A)^{i_1}{}_{j_1}(t_B)^{i_2}{}_{j_2}\bar W^{j_1 j_2 i_3 i_4}\]

Full coefficient:

\[
\begin{aligned}
&\frac{\lambda_{\gamma_4,+\dot a}}{2160}\left[\begin{aligned}
&\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}-2\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle+3\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&+3\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}-6\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+3\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-2\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+6\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-3\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+\dot a}}{4320}\left[\begin{aligned}
&11\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}-18\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle+29\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&+12\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}-42\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+21\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-4\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+12\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-12\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-3\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_5,+\dot a}}{2160}\left[\begin{aligned}
&\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-3\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+3\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]
\end{aligned}
\]

### 2

\[(t_A)^{i_1}{}_{j_1}(t_B)^{i_3}{}_{j_3}\bar W^{j_1 i_2 j_3 i_4}\]

Full coefficient:

\[
\begin{aligned}
&\frac{\lambda_{\gamma_4,+\dot a}}{2160}\left[\begin{aligned}
&4\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}-12\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+12\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+3\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-3\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+\dot a}}{4320}\left[\begin{aligned}
&6\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}-16\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+12\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+3\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-14\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+24\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-18\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_5,+\dot a}}{1080}\left[\begin{aligned}
&-\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+3\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-3\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]
\end{aligned}
\]

### 3

\[(t_A)^{i_1}{}_{j_1}(t_B)^{j_1}{}_{k_1}\bar W^{k_1 i_2 i_3 i_4}\]

Full coefficient:

\[
\begin{aligned}
&\frac{\lambda_{\gamma_4,+\dot a}}{2160}\left[\begin{aligned}
&6\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}-2\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle+18\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-12\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}+24\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-12\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+9\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-15\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+3\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+\dot a}}{4320}\left[\begin{aligned}
&\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}-3\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle+4\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-9\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}+28\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-24\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-2\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+26\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-30\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+18\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_5,+\dot a}}{1080}\left[\begin{aligned}
&-\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+3\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+3\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]
\end{aligned}
\]

### 4

\[(t_B)^{i_1}{}_{j_1}(t_A)^{i_2}{}_{j_2}\bar W^{j_1 j_2 i_3 i_4}\]

Full coefficient:

\[
\begin{aligned}
&\frac{\lambda_{\gamma_4,+\dot a}}{4320}\left[\begin{aligned}
&-2\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}+9\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle-6\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-16\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}+27\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-6\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+14\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-42\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+21\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+\dot a}}{4320}\left[\begin{aligned}
&4\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}-7\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle+16\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&+3\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}-28\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+24\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+18\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-48\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+18\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_5,+\dot a}}{4320}\left[\begin{aligned}
&-7\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+21\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-21\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]
\end{aligned}
\]

### 5

\[(t_B)^{i_1}{}_{j_1}(t_A)^{i_4}{}_{j_4}\bar W^{j_1 i_2 i_3 j_4}\]

Full coefficient:

\[
\begin{aligned}
&\frac{\lambda_{\gamma_4,+\dot a}}{2160}\left[\begin{aligned}
&2\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}-8\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle+6\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&+18\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}-16\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+3\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-22\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+36\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-8\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+\dot a}}{4320}\left[\begin{aligned}
&-\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}+6\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle-4\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-3\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}+48\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-14\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-48\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+78\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}-8\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_5,+\dot a}}{120}\left[\begin{aligned}
&-\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]
\end{aligned}
\]

### 6

\[(t_B)^{i_1}{}_{j_1}(t_A)^{j_1}{}_{k_1}\bar W^{k_1 i_2 i_3 i_4}\]

Full coefficient:

\[
\begin{aligned}
&\frac{\lambda_{\gamma_4,+\dot a}}{2160}\left[\begin{aligned}
&-2\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}+9\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle-6\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-12\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}+15\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+6\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-2\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-24\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+3\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+\dot a}}{4320}\left[\begin{aligned}
&-2\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}-9\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle-8\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-3\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}-52\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}\\
&-8\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-72\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}-12\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_5,+\dot a}}{360}\left[\begin{aligned}
&\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-3\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-2\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]
\end{aligned}
\]

### 7

\[(t_A)^{i_2}{}_{j_2}(t_B)^{i_3}{}_{j_3}\bar W^{i_1 j_2 j_3 i_4}\]

Full coefficient:

\[
\begin{aligned}
&\frac{\lambda_{\gamma_4,+\dot a}}{432}\left[\begin{aligned}
&\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-3\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+3\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+\dot a}}{432}\left[\begin{aligned}
&\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-3\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}+3\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_5,+\dot a}}{864}\left[\begin{aligned}
&\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-3\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+3\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]
\end{aligned}
\]

### 8

\[(t_A)^{i_2}{}_{j_2}(t_B)^{j_2}{}_{k_2}\bar W^{i_1 k_2 i_3 i_4}\]

Full coefficient:

\[
\begin{aligned}
&\frac{\lambda_{\gamma_4,+\dot a}}{432}\left[\begin{aligned}
&\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle+3\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-3\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}\\
&+6\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-3\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+\dot a}}{864}\left[\begin{aligned}
&-2\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}-\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle-2\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&+3\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}-4\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+6\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-2\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+12\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-6\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_5,+\dot a}}{864}\left[\begin{aligned}
&\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-3\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-3\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]
\end{aligned}
\]

### 9

\[(t_B)^{i_2}{}_{j_2}(t_A)^{i_4}{}_{j_4}\bar W^{i_1 j_2 i_3 j_4}\]

Full coefficient:

\[
\begin{aligned}
&\frac{\lambda_{\gamma_4,+\dot a}}{4320}\left[\begin{aligned}
&-2\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}+8\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle-6\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-18\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}+16\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-3\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+22\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-36\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+8\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+\dot a}}{2160}\left[\begin{aligned}
&-6\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}+21\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle-14\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-18\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}+33\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-4\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-3\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+3\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}+2\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_5,+\dot a}}{240}\left[\begin{aligned}
&\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]
\end{aligned}
\]

### 10

\[(t_B)^{i_2}{}_{j_2}(t_A)^{j_2}{}_{k_2}\bar W^{i_1 k_2 i_3 i_4}\]

Full coefficient:

\[
\begin{aligned}
&\frac{\lambda_{\gamma_4,+\dot a}}{4320}\left[\begin{aligned}
&4\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}-3\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle+12\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&+18\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}+3\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+2\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}\\
&+42\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+15\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+\dot a}}{4320}\left[\begin{aligned}
&7\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}-36\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle+13\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-6\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}-46\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-9\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+7\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-26\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+18\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+15\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_5,+\dot a}}{720}\left[\begin{aligned}
&-\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+3\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+5\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]
\end{aligned}
\]

## Complete compact specialization with the free spectator set to one

Set \(\gamma_5=1\), hence \(\lambda_{\gamma_5}=0\). The table retains all ten color words, including the two that become zero. Factor out

\[\frac{\langle\lambda_{\gamma_4},\lambda_c\rangle^3}{4320}.\]

| Color tensor | Remaining full spinor coefficient |
| --- | --- |
| \((t_A)^{i_1}{}_{j_1}(t_B)^{i_2}{}_{j_2}\bar W^{j_1 j_2 i_3 i_4}\) | \(2\lambda_{\gamma_4,+\dot a}+11\lambda_{c,+\dot a}\) |
| \((t_A)^{i_1}{}_{j_1}(t_B)^{i_3}{}_{j_3}\bar W^{j_1 i_2 j_3 i_4}\) | \(0\) |
| \((t_A)^{i_1}{}_{j_1}(t_B)^{j_1}{}_{k_1}\bar W^{k_1 i_2 i_3 i_4}\) | \(12\lambda_{\gamma_4,+\dot a}+\lambda_{c,+\dot a}\) |
| \((t_B)^{i_1}{}_{j_1}(t_A)^{i_2}{}_{j_2}\bar W^{j_1 j_2 i_3 i_4}\) | \(-2\lambda_{\gamma_4,+\dot a}+4\lambda_{c,+\dot a}\) |
| \((t_B)^{i_1}{}_{j_1}(t_A)^{i_4}{}_{j_4}\bar W^{j_1 i_2 i_3 j_4}\) | \(4\lambda_{\gamma_4,+\dot a}-\lambda_{c,+\dot a}\) |
| \((t_B)^{i_1}{}_{j_1}(t_A)^{j_1}{}_{k_1}\bar W^{k_1 i_2 i_3 i_4}\) | \(-4\lambda_{\gamma_4,+\dot a}-2\lambda_{c,+\dot a}\) |
| \((t_A)^{i_2}{}_{j_2}(t_B)^{i_3}{}_{j_3}\bar W^{i_1 j_2 j_3 i_4}\) | \(0\) |
| \((t_A)^{i_2}{}_{j_2}(t_B)^{j_2}{}_{k_2}\bar W^{i_1 k_2 i_3 i_4}\) | \(-10\lambda_{c,+\dot a}\) |
| \((t_B)^{i_2}{}_{j_2}(t_A)^{i_4}{}_{j_4}\bar W^{i_1 j_2 i_3 j_4}\) | \(-2\lambda_{\gamma_4,+\dot a}-12\lambda_{c,+\dot a}\) |
| \((t_B)^{i_2}{}_{j_2}(t_A)^{j_2}{}_{k_2}\bar W^{i_1 k_2 i_3 i_4}\) | \(4\lambda_{\gamma_4,+\dot a}+7\lambda_{c,+\dot a}\) |

In the physical derivative formula the corresponding common multiplier, including this compact denominator, is \(g^2/[552960(8\pi^2)^3]\).

## Full coefficients of the two representation defects

Before applying \([t_A,t_B]=iC^C{}_{AB}t_C\), the complete \(\zeta(3)\) summand consists of the following two displayed spinors multiplying their displayed tensor defects. These spinors retain all three original source spinors. The tensor defects vanish exactly.

### Defect on superpotential slot 1

\[\bigl([t_A,t_B]-iC^C{}_{AB}t_C\bigr)^{i_1}{}_{j_1}\bar W^{j_1 i_2 i_3 i_4}\]

Complete spinor coefficient:

\[
\begin{aligned}
&\frac{\lambda_{\gamma_4,+\dot a}}{2100}\left[\begin{aligned}
&3\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle+3\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+19\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-3\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}-3\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+19\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle
\end{aligned}\right]\\
&+\frac{\lambda_{c,+\dot a}}{2100}\left[\begin{aligned}
&6\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle-3\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-6\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\\
&+15\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-2\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}-6\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&+21\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}-5\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_5,+\dot a}}{700}\left[\begin{aligned}
&\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]
\end{aligned}
\]

### Defect on superpotential slot 2

\[\bigl([t_A,t_B]-iC^C{}_{AB}t_C\bigr)^{i_2}{}_{j_2}\bar W^{i_1 j_2 i_3 i_4}\]

Complete spinor coefficient:

\[
\begin{aligned}
&\frac{\lambda_{\gamma_4,+\dot a}}{2100}\left[\begin{aligned}
&-3\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle+6\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-2\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\\
&+15\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-6\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}-5\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}\\
&+21\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-6\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+\dot a}}{2100}\left[\begin{aligned}
&3\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle+3\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-3\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\\
&+19\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+19\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}-3\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_5,+\dot a}}{700}\left[\begin{aligned}
&\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]
\end{aligned}
\]

The complete unreduced twelve-word tensor, including all rational terms before using the representation relation, is retained in [quartic_tensor_certificate.json](quartic_tensor_certificate.json). No color contraction or special source assignment is required for either cancellation.
