# Complete order-2 input-jet tensor

Generated and checked coefficient by coefficient by [render_jets.py](render_jets.py).

The polarized superpotential is \(\bar W(\gamma)=\bar W^{i_1i_2i_3i_4}\gamma_{1i_1}\gamma_{2i_2}\gamma_{3i_3}\gamma_{4i_4}\). The inputs are

\[\beta^{1i_1}(\partial_{+\dot a}\partial_{+\dot b}\beta^{2i_2})\gamma_5,\qquad b_A\beta^{3i_3}.\]

The output fields are \(\gamma_{4i_4},c^B,\gamma_5\), with every input/output color index and every displayed spinor index retained. The spectator \(\gamma_5\) is free, gauge neutral, and absent from \(\bar W\).

Sources carry the manuscript's lower indices. Before using angle brackets, fix their convention:

\[\langle\lambda_{\gamma_4},\lambda_c\rangle=\epsilon^{\dot b\dot c}\lambda_{\gamma_4,+\dot b}\lambda_{c,+\dot c},\qquad \epsilon^{\dot+\dot-}=1=-\epsilon^{\dot-\dot+}.\]

The same antisymmetric pairing defines the other brackets. Parentheses on the free spinor indices mean the average over every permutation of those indices. These are the actual symmetric derivative indices; all their values are retained. For example,

\[\lambda_{\gamma_4,+(\dot a}\lambda_{c,+\dot b)}=\frac{\lambda_{\gamma_4,+\dot a}\lambda_{c,+\dot b}+\lambda_{\gamma_4,+\dot b}\lambda_{c,+\dot a}}{2}.\]

For three indices the average contains all six index permutations. Each color tensor below multiplies its entire displayed spinor tensor. Replace each source spinor by the derivative acting on its own surviving field, and multiply their complete sum by

\[\frac{g^2}{128(8\pi^2)^3}.\]

This is the complete \(g^2\bar W\) coefficient at \(\hbar^3\). The factors from contracted input derivatives and from \(\lambda=-\partial/2\) have already been included in this common multiplier. See [DERIVATION.md](DERIVATION.md).

## Complete rational tensor after the representation identity

### Color tensor 1

\[(t_A)^{i_1}{}_{j_1}(t_B)^{i_2}{}_{j_2}\bar W^{j_1 j_2 i_3 i_4}\]

Full coefficient:

\[
\begin{aligned}
&\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_4,+\dot b)}}{8640}\left[\begin{aligned}
&\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}-2\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle+3\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&+3\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}-6\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+3\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-3\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+6\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-3\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{c,+\dot b)}}{4320}\left[\begin{aligned}
&3\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}-5\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle+8\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&+6\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}-12\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+6\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-3\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+9\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-3\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_5,+\dot b)}}{8640}\left[\begin{aligned}
&2\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-6\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+3\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{c,+\dot b)}}{8640}\left[\begin{aligned}
&17\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}-24\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle+41\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&+12\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}-48\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+24\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-3\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+9\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-9\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-5\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{\gamma_5,+\dot b)}}{8640}\left[\begin{aligned}
&2\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-6\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+6\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_5,+(\dot a}\lambda_{\gamma_5,+\dot b)}}{8640}\left[\begin{aligned}
&-\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+3\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-3\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]
\end{aligned}
\]

### Color tensor 2

\[(t_A)^{i_1}{}_{j_1}(t_B)^{i_3}{}_{j_3}\bar W^{j_1 i_2 j_3 i_4}\]

Full coefficient:

\[
\begin{aligned}
&\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_4,+\dot b)}}{2880}\left[\begin{aligned}
&2\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}-6\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+6\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+3\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-3\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{c,+\dot b)}}{2880}\left[\begin{aligned}
&3\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}-8\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+6\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-2\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+6\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}-6\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{c,+\dot b)}}{14400}\left[\begin{aligned}
&6\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}-15\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+10\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+3\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-15\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+25\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-15\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{\gamma_5,+\dot b)}}{2880}\left[\begin{aligned}
&-\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+4\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-6\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+4\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_5,+(\dot a}\lambda_{\gamma_5,+\dot b)}}{2880}\left[\begin{aligned}
&\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-3\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+3\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]
\end{aligned}
\]

### Color tensor 3

\[(t_A)^{i_1}{}_{j_1}(t_B)^{j_1}{}_{k_1}\bar W^{k_1 i_2 i_3 i_4}\]

Full coefficient:

\[
\begin{aligned}
&\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_4,+\dot b)}}{2880}\left[\begin{aligned}
&4\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}-2\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle+12\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-6\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}+12\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-6\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+5\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-9\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+3\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{c,+\dot b)}}{2880}\left[\begin{aligned}
&\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}-2\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle+4\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-6\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}+16\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-12\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+2\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-6\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}+6\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_5,+\dot b)}}{1440}\left[\begin{aligned}
&-2\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+3\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{c,+\dot b)}}{14400}\left[\begin{aligned}
&2\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}-4\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle+10\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-14\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}+25\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-10\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+5\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-35\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+25\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{\gamma_5,+\dot b)}}{2880}\left[\begin{aligned}
&\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-8\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+6\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-4\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_5,+(\dot a}\lambda_{\gamma_5,+\dot b)}}{2880}\left[\begin{aligned}
&\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-3\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-3\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]
\end{aligned}
\]

### Color tensor 4

\[(t_B)^{i_1}{}_{j_1}(t_A)^{i_2}{}_{j_2}\bar W^{j_1 j_2 i_3 i_4}\]

Full coefficient:

\[
\begin{aligned}
&\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_4,+\dot b)}}{25920}\left[\begin{aligned}
&-4\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}+12\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle-12\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-35\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}+36\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-12\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+54\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-93\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+24\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{c,+\dot b)}}{25920}\left[\begin{aligned}
&3\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}-28\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle+12\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&+45\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}-112\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+18\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-24\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+204\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-192\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+20\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_5,+\dot b)}}{25920}\left[\begin{aligned}
&-46\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+138\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-69\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{c,+\dot b)}}{129600}\left[\begin{aligned}
&-52\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}+93\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle-260\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-48\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}+465\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-520\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+12\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-180\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+870\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-500\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{\gamma_5,+\dot b)}}{25920}\left[\begin{aligned}
&8\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-78\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+186\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-55\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_5,+(\dot a}\lambda_{\gamma_5,+\dot b)}}{25920}\left[\begin{aligned}
&23\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-69\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+69\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]
\end{aligned}
\]

### Color tensor 5

\[(t_B)^{i_1}{}_{j_1}(t_A)^{i_4}{}_{j_4}\bar W^{j_1 i_2 i_3 j_4}\]

Full coefficient:

\[
\begin{aligned}
&\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_4,+\dot b)}}{14400}\left[\begin{aligned}
&4\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}-15\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle+12\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&+30\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}-30\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+6\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-60\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+60\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-15\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{c,+\dot b)}}{2880}\left[\begin{aligned}
&-6\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}+8\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-3\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+4\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-36\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+16\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-2\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_5,+\dot b)}}{1440}\left[\begin{aligned}
&7\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-12\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+3\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{c,+\dot b)}}{21600}\left[\begin{aligned}
&4\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}-30\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle+20\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&+18\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}-105\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+25\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+45\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-210\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}+30\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{\gamma_5,+\dot b)}}{1440}\left[\begin{aligned}
&9\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-15\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}+\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_5,+(\dot a}\lambda_{\gamma_5,+\dot b)}}{240}\left[\begin{aligned}
&\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]
\end{aligned}
\]

### Color tensor 6

\[(t_B)^{i_1}{}_{j_1}(t_A)^{j_1}{}_{k_1}\bar W^{k_1 i_2 i_3 i_4}\]

Full coefficient:

\[
\begin{aligned}
&\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_4,+\dot b)}}{2880}\left[\begin{aligned}
&-\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}+4\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle-3\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-5\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}+6\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+3\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+6\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-9\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{c,+\dot b)}}{2880}\left[\begin{aligned}
&-\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}-4\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+6\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\\
&-8\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+2\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+28\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-12\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}-4\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_5,+\dot b)}}{2880}\left[\begin{aligned}
&2\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+12\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-3\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{c,+\dot b)}}{14400}\left[\begin{aligned}
&\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}+14\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle+5\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&+4\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}+55\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+20\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-2\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+20\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+120\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+10\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{\gamma_5,+\dot b)}}{2880}\left[\begin{aligned}
&-\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+6\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+24\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+3\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_5,+(\dot a}\lambda_{\gamma_5,+\dot b)}}{1440}\left[\begin{aligned}
&-2\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+6\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+3\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]
\end{aligned}
\]

### Color tensor 7

\[(t_A)^{i_2}{}_{j_2}(t_B)^{i_3}{}_{j_3}\bar W^{i_1 j_2 j_3 i_4}\]

Full coefficient:

\[
\begin{aligned}
&\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_4,+\dot b)}}{12960}\left[\begin{aligned}
&-2\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}+6\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-6\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-3\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+9\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-9\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{c,+\dot b)}}{12960}\left[\begin{aligned}
&-3\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}+8\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-6\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-3\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+6\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-6\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_5,+\dot b)}}{6480}\left[\begin{aligned}
&-7\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+21\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-21\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{c,+\dot b)}}{64800}\left[\begin{aligned}
&-6\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}+15\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-10\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-3\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+15\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}-15\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{\gamma_5,+\dot b)}}{25920}\left[\begin{aligned}
&-\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-24\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+78\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-80\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_5,+(\dot a}\lambda_{\gamma_5,+\dot b)}}{25920}\left[\begin{aligned}
&-13\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+39\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-39\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]
\end{aligned}
\]

### Color tensor 8

\[(t_A)^{i_2}{}_{j_2}(t_B)^{j_2}{}_{k_2}\bar W^{i_1 k_2 i_3 i_4}\]

Full coefficient:

\[
\begin{aligned}
&\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_4,+\dot b)}}{12960}\left[\begin{aligned}
&-4\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}-4\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle-12\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-7\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}-30\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+6\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+7\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-60\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+9\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{c,+\dot b)}}{12960}\left[\begin{aligned}
&-4\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}+16\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle-16\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&+18\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}-20\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-6\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-21\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+30\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-96\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+6\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_5,+\dot b)}}{12960}\left[\begin{aligned}
&40\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-81\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+42\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{c,+\dot b)}}{129600}\left[\begin{aligned}
&-139\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}-81\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle-95\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&+204\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}-195\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+170\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-42\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+240\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}-380\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{\gamma_5,+\dot b)}}{25920}\left[\begin{aligned}
&17\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+12\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-162\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+80\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_5,+(\dot a}\lambda_{\gamma_5,+\dot b)}}{25920}\left[\begin{aligned}
&-13\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+39\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+39\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]
\end{aligned}
\]

### Color tensor 9

\[(t_B)^{i_2}{}_{j_2}(t_A)^{i_4}{}_{j_4}\bar W^{i_1 j_2 i_3 j_4}\]

Full coefficient:

\[
\begin{aligned}
&\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_4,+\dot b)}}{43200}\left[\begin{aligned}
&-4\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}+15\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle-12\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-30\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}+30\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-6\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+60\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-60\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+15\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{c,+\dot b)}}{8640}\left[\begin{aligned}
&-5\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}+16\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle-12\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-30\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}+24\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-3\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+24\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-36\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+2\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_5,+\dot b)}}{4320}\left[\begin{aligned}
&-7\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+12\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-3\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{c,+\dot b)}}{4320}\left[\begin{aligned}
&-8\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}+24\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle-16\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-15\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}+33\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-3\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+2\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{\gamma_5,+\dot b)}}{4320}\left[\begin{aligned}
&3\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-3\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}-\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_5,+(\dot a}\lambda_{\gamma_5,+\dot b)}}{720}\left[\begin{aligned}
&-\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]
\end{aligned}
\]

### Color tensor 10

\[(t_B)^{i_2}{}_{j_2}(t_A)^{j_2}{}_{k_2}\bar W^{i_1 k_2 i_3 i_4}\]

Full coefficient:

\[
\begin{aligned}
&\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_4,+\dot b)}}{25920}\left[\begin{aligned}
&4\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}-16\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle+12\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-3\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}-30\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-6\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-38\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-27\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-12\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{c,+\dot b)}}{25920}\left[\begin{aligned}
&23\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}-12\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle+44\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&+78\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}+12\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-24\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+18\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+84\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-4\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_5,+\dot b)}}{8640}\left[\begin{aligned}
&-2\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-32\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-15\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{c,+\dot b)}}{43200}\left[\begin{aligned}
&39\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}-182\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle+55\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-42\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}-175\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-30\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+24\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-90\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+100\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+40\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{\gamma_5,+\dot b)}}{8640}\left[\begin{aligned}
&-5\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+18\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-4\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-13\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_5,+(\dot a}\lambda_{\gamma_5,+\dot b)}}{4320}\left[\begin{aligned}
&2\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-6\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-13\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]
\end{aligned}
\]

## Complete coefficients of the representation defects

The complete unreduced twelve-word tensor is retained in the certificate. Its entire \(\zeta(3)\) contribution consists of the following two full spinor tensors multiplying the displayed exact representation defects. No independence of period values is used.

### Superpotential slot 1

\[\zeta(3)\bigl([t_A,t_B]-iC^C{}_{AB}t_C\bigr)^{i_1}{}_{j_1}\bar W^{j_1 i_2 i_3 i_4}\]

Full coefficient:

\[
\begin{aligned}
&\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_4,+\dot b)}}{2100}\left[\begin{aligned}
&3\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle+\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+3\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\\
&+7\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-3\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}+3\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}\\
&-3\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{c,+\dot b)}}{1050}\left[\begin{aligned}
&3\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle-3\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+3\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-5\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}+3\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-6\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-2\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_5,+\dot b)}}{2100}\left[\begin{aligned}
&3\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-7\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-3\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{c,+\dot b)}}{700}\left[\begin{aligned}
&-\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-2\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}+3\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-3\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}-2\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{\gamma_5,+\dot b)}}{2100}\left[\begin{aligned}
&9\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-17\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}+\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_5,+(\dot a}\lambda_{\gamma_5,+\dot b)}}{2100}\left[\begin{aligned}
&-\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-3\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]
\end{aligned}
\]

### Superpotential slot 2

\[\zeta(3)\bigl([t_A,t_B]-iC^C{}_{AB}t_C\bigr)^{i_2}{}_{j_2}\bar W^{i_1 j_2 i_3 i_4}\]

Full coefficient:

\[
\begin{aligned}
&\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_4,+\dot b)}}{700}\left[\begin{aligned}
&-\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle-2\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}-2\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}\\
&-3\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+3\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{c,+\dot b)}}{1050}\left[\begin{aligned}
&-3\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle+3\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-5\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\\
&+3\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-2\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-6\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&+3\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_5,+\dot b)}}{2100}\left[\begin{aligned}
&\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-17\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+9\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{c,+\dot b)}}{2100}\left[\begin{aligned}
&\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle+3\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-3\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\\
&+7\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+3\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}-3\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&+3\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{\gamma_5,+\dot b)}}{2100}\left[\begin{aligned}
&-3\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-7\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}+3\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_5,+(\dot a}\lambda_{\gamma_5,+\dot b)}}{2100}\left[\begin{aligned}
&-3\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]
\end{aligned}
\]

## Exact completeness

Every rational and period coefficient before and after the representation relation, every original geometric tensor, every source interpolation value, and all further full-source replays are retained in [jet_order_2_certificate.json](jet_order_2_certificate.json).

The all-finite-jet rational boundary evaluator independently reproduces every entry of this tensor at every proven-degree source interpolation node. Its proof and the strict scope are in [DERIVATION.md](DERIVATION.md) and [the source-simplex derivation](../forest-jet-identities/DERIVATION.md).
