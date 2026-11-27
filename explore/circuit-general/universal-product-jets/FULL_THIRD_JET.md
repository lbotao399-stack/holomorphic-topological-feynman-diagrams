# Complete order-3 input-jet tensor

Generated and checked coefficient by coefficient by [render_jets.py](render_jets.py).

The polarized superpotential is \(\bar W(\gamma)=\bar W^{i_1i_2i_3i_4}\gamma_{1i_1}\gamma_{2i_2}\gamma_{3i_3}\gamma_{4i_4}\). The inputs are

\[\beta^{1i_1}(\partial_{+\dot a}\partial_{+\dot b}\partial_{+\dot c}\beta^{2i_2})\gamma_5,\qquad b_A\beta^{3i_3}.\]

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
&\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_4,+\dot b}\lambda_{\gamma_4,+\dot c)}}{25200}\left[\begin{aligned}
&\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}-2\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle+3\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&+3\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}-6\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+3\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-3\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+6\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-3\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_4,+\dot b}\lambda_{c,+\dot c)}}{100800}\left[\begin{aligned}
&27\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}-46\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle+73\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&+57\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}-114\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+57\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-48\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+87\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-30\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-9\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_4,+\dot b}\lambda_{\gamma_5,+\dot c)}}{8400}\left[\begin{aligned}
&\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-2\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{c,+\dot b}\lambda_{c,+\dot c)}}{50400}\left[\begin{aligned}
&41\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}-59\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle+100\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&+60\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}-120\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+60\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-24\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+72\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-12\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-16\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{c,+\dot b}\lambda_{\gamma_5,+\dot c)}}{33600}\left[\begin{aligned}
&8\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-24\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+9\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+2\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_5,+\dot b}\lambda_{\gamma_5,+\dot c)}}{25200}\left[\begin{aligned}
&-2\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+6\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-3\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{c,+\dot b}\lambda_{c,+\dot c)}}{100800}\left[\begin{aligned}
&159\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}-200\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle+359\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&+80\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}-360\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+180\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-16\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+48\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-48\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-44\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{c,+\dot b}\lambda_{\gamma_5,+\dot c)}}{100800}\left[\begin{aligned}
&12\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-36\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+36\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+13\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{\gamma_5,+\dot b}\lambda_{\gamma_5,+\dot c)}}{100800}\left[\begin{aligned}
&-8\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+24\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-24\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-3\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_5,+(\dot a}\lambda_{\gamma_5,+\dot b}\lambda_{\gamma_5,+\dot c)}}{25200}\left[\begin{aligned}
&\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-3\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+3\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]
\end{aligned}
\]

### Color tensor 2

\[(t_A)^{i_1}{}_{j_1}(t_B)^{i_3}{}_{j_3}\bar W^{j_1 i_2 j_3 i_4}\]

Full coefficient:

\[
\begin{aligned}
&\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_4,+\dot b}\lambda_{\gamma_4,+\dot c)}}{25200}\left[\begin{aligned}
&8\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}-24\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+24\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-5\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+15\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-15\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_4,+\dot b}\lambda_{c,+\dot c)}}{16800}\left[\begin{aligned}
&12\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}-32\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+24\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-3\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+2\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+12\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-18\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_4,+\dot b}\lambda_{\gamma_5,+\dot c)}}{12600}\left[\begin{aligned}
&\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-3\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+3\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{c,+\dot b}\lambda_{c,+\dot c)}}{42000}\left[\begin{aligned}
&24\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}-60\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+40\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+3\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-30\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+65\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-45\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{c,+\dot b}\lambda_{\gamma_5,+\dot c)}}{25200}\left[\begin{aligned}
&-3\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+16\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-30\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+24\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_5,+\dot b}\lambda_{\gamma_5,+\dot c)}}{25200}\left[\begin{aligned}
&\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-3\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+3\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{c,+\dot b}\lambda_{c,+\dot c)}}{252000}\left[\begin{aligned}
&40\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}-96\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+60\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+20\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-102\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+165\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-90\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{c,+\dot b}\lambda_{\gamma_5,+\dot c)}}{126000}\left[\begin{aligned}
&-18\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+75\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-110\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+60\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{\gamma_5,+\dot b}\lambda_{\gamma_5,+\dot c)}}{50400}\left[\begin{aligned}
&9\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-34\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+48\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-30\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_5,+(\dot a}\lambda_{\gamma_5,+\dot b}\lambda_{\gamma_5,+\dot c)}}{6300}\left[\begin{aligned}
&-\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+3\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-3\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]
\end{aligned}
\]

### Color tensor 3

\[(t_A)^{i_1}{}_{j_1}(t_B)^{j_1}{}_{k_1}\bar W^{k_1 i_2 i_3 i_4}\]

Full coefficient:

\[
\begin{aligned}
&\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_4,+\dot b}\lambda_{\gamma_4,+\dot c)}}{25200}\left[\begin{aligned}
&20\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}-12\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle+60\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-24\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}+48\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-24\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+21\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-39\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+15\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_4,+\dot b}\lambda_{c,+\dot c)}}{16800}\left[\begin{aligned}
&6\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}-11\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle+24\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-27\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}+68\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-48\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+9\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-14\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-6\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+18\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_4,+\dot b}\lambda_{\gamma_5,+\dot c)}}{4200}\left[\begin{aligned}
&-3\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+5\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{c,+\dot b}\lambda_{c,+\dot c)}}{42000}\left[\begin{aligned}
&9\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}-22\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle+45\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-50\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}+100\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-50\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+17\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+20\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-105\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+75\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{c,+\dot b}\lambda_{\gamma_5,+\dot c)}}{16800}\left[\begin{aligned}
&5\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-12\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+16\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-16\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_5,+\dot b}\lambda_{\gamma_5,+\dot c)}}{8400}\left[\begin{aligned}
&5\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-7\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{c,+\dot b}\lambda_{c,+\dot c)}}{252000}\left[\begin{aligned}
&5\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}-25\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle+30\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-65\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}+186\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-135\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-15\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+192\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-180\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+70\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{c,+\dot b}\lambda_{\gamma_5,+\dot c)}}{84000}\left[\begin{aligned}
&2\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+15\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+100\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-60\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{\gamma_5,+\dot b}\lambda_{\gamma_5,+\dot c)}}{8400}\left[\begin{aligned}
&-2\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+13\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-7\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+5\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_5,+(\dot a}\lambda_{\gamma_5,+\dot b}\lambda_{\gamma_5,+\dot c)}}{6300}\left[\begin{aligned}
&-\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+3\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+3\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]
\end{aligned}
\]

### Color tensor 4

\[(t_B)^{i_1}{}_{j_1}(t_A)^{i_2}{}_{j_2}\bar W^{j_1 j_2 i_3 i_4}\]

Full coefficient:

\[
\begin{aligned}
&\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_4,+\dot b}\lambda_{\gamma_4,+\dot c)}}{302400}\left[\begin{aligned}
&-16\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}+60\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle-48\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-118\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}+180\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-48\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+265\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-306\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+132\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_4,+\dot b}\lambda_{c,+\dot c)}}{302400}\left[\begin{aligned}
&18\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}-47\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle+72\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&+297\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}-188\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+108\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-501\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+1332\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-426\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+120\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_4,+\dot b}\lambda_{\gamma_5,+\dot c)}}{100800}\left[\begin{aligned}
&-128\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+221\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-58\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{c,+\dot b}\lambda_{c,+\dot c)}}{1512000}\left[\begin{aligned}
&-86\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}+1233\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle-430\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-2106\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}+6165\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-860\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+960\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-9810\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+11610\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-620\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{c,+\dot b}\lambda_{\gamma_5,+\dot c)}}{302400}\left[\begin{aligned}
&282\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-1896\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+1551\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-122\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_5,+\dot b}\lambda_{\gamma_5,+\dot c)}}{302400}\left[\begin{aligned}
&326\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-978\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+489\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{c,+\dot b}\lambda_{c,+\dot c)}}{1512000}\left[\begin{aligned}
&299\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}-501\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle+1794\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&+150\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}-3006\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+4485\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-30\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+1140\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-7755\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+6060\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{c,+\dot b}\lambda_{\gamma_5,+\dot c)}}{100800}\left[\begin{aligned}
&-16\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+174\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-664\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+293\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{\gamma_5,+\dot b}\lambda_{\gamma_5,+\dot c)}}{151200}\left[\begin{aligned}
&-47\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+351\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-771\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+188\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_5,+(\dot a}\lambda_{\gamma_5,+\dot b}\lambda_{\gamma_5,+\dot c)}}{302400}\left[\begin{aligned}
&-163\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+489\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-489\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]
\end{aligned}
\]

### Color tensor 5

\[(t_B)^{i_1}{}_{j_1}(t_A)^{i_4}{}_{j_4}\bar W^{j_1 i_2 i_3 j_4}\]

Full coefficient:

\[
\begin{aligned}
&\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_4,+\dot b}\lambda_{\gamma_4,+\dot c)}}{189000}\left[\begin{aligned}
&20\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}-72\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle+60\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&+135\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}-144\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+30\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-240\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+270\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-72\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_4,+\dot b}\lambda_{c,+\dot c)}}{28000}\left[\begin{aligned}
&\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}-5\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle+4\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-10\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}+20\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-10\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+60\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-100\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+55\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-8\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_4,+\dot b}\lambda_{\gamma_5,+\dot c)}}{12600}\left[\begin{aligned}
&30\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-32\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+9\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{c,+\dot b}\lambda_{c,+\dot c)}}{25200}\left[\begin{aligned}
&3\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}-20\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle+15\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&+54\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}-52\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+12\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-30\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+216\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-92\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+15\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{c,+\dot b}\lambda_{\gamma_5,+\dot c)}}{50400}\left[\begin{aligned}
&-75\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+456\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-178\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+18\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_5,+\dot b}\lambda_{\gamma_5,+\dot c)}}{6300}\left[\begin{aligned}
&-17\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+30\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-8\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{c,+\dot b}\lambda_{c,+\dot c)}}{151200}\left[\begin{aligned}
&-5\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}+45\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle-30\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-15\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}+414\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-123\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-198\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+1017\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}-144\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{c,+\dot b}\lambda_{\gamma_5,+\dot c)}}{50400}\left[\begin{aligned}
&-129\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+444\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}-38\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{\gamma_5,+\dot b}\lambda_{\gamma_5,+\dot c)}}{12600}\left[\begin{aligned}
&-52\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+87\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}-4\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_5,+(\dot a}\lambda_{\gamma_5,+\dot b}\lambda_{\gamma_5,+\dot c)}}{420}\left[\begin{aligned}
&-\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]
\end{aligned}
\]

### Color tensor 6

\[(t_B)^{i_1}{}_{j_1}(t_A)^{j_1}{}_{k_1}\bar W^{k_1 i_2 i_3 i_4}\]

Full coefficient:

\[
\begin{aligned}
&\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_4,+\dot b}\lambda_{\gamma_4,+\dot c)}}{25200}\left[\begin{aligned}
&-4\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}+15\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle-12\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-18\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}+21\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+12\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+21\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-30\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-3\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_4,+\dot b}\lambda_{c,+\dot c)}}{16800}\left[\begin{aligned}
&-4\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}+5\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle-16\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&+9\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}-12\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-27\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}\\
&+56\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-24\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}-12\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_4,+\dot b}\lambda_{\gamma_5,+\dot c)}}{8400}\left[\begin{aligned}
&-8\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+13\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-2\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{c,+\dot b}\lambda_{c,+\dot c)}}{42000}\left[\begin{aligned}
&\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}+17\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle+5\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-50\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}+25\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+50\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-22\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-235\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+70\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+20\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{c,+\dot b}\lambda_{\gamma_5,+\dot c)}}{16800}\left[\begin{aligned}
&-11\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-100\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+46\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+12\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_5,+\dot b}\lambda_{\gamma_5,+\dot c)}}{8400}\left[\begin{aligned}
&-4\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-16\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+5\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{c,+\dot b}\lambda_{c,+\dot c)}}{252000}\left[\begin{aligned}
&-15\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}-65\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle-90\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-25\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}-486\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-165\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+5\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-162\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-1140\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-280\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{c,+\dot b}\lambda_{\gamma_5,+\dot c)}}{84000}\left[\begin{aligned}
&18\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-145\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-520\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-50\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{\gamma_5,+\dot b}\lambda_{\gamma_5,+\dot c)}}{8400}\left[\begin{aligned}
&3\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-16\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-40\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-4\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_5,+(\dot a}\lambda_{\gamma_5,+\dot b}\lambda_{\gamma_5,+\dot c)}}{6300}\left[\begin{aligned}
&5\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-15\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-6\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]
\end{aligned}
\]

### Color tensor 7

\[(t_A)^{i_2}{}_{j_2}(t_B)^{i_3}{}_{j_3}\bar W^{i_1 j_2 j_3 i_4}\]

Full coefficient:

\[
\begin{aligned}
&\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_4,+\dot b}\lambda_{\gamma_4,+\dot c)}}{5400}\left[\begin{aligned}
&\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-3\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+3\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_4,+\dot b}\lambda_{c,+\dot c)}}{3600}\left[\begin{aligned}
&\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-2\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+2\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_4,+\dot b}\lambda_{\gamma_5,+\dot c)}}{3600}\left[\begin{aligned}
&\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-3\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+3\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{c,+\dot b}\lambda_{c,+\dot c)}}{9000}\left[\begin{aligned}
&\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-5\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}+5\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{c,+\dot b}\lambda_{\gamma_5,+\dot c)}}{3600}\left[\begin{aligned}
&\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-2\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+2\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_5,+\dot b}\lambda_{\gamma_5,+\dot c)}}{21600}\left[\begin{aligned}
&13\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-39\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+39\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{c,+\dot b}\lambda_{c,+\dot c)}}{54000}\left[\begin{aligned}
&6\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-15\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}+10\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{c,+\dot b}\lambda_{\gamma_5,+\dot c)}}{18000}\left[\begin{aligned}
&\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-5\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}+5\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{\gamma_5,+\dot b}\lambda_{\gamma_5,+\dot c)}}{21600}\left[\begin{aligned}
&\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+9\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-33\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+35\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_5,+(\dot a}\lambda_{\gamma_5,+\dot b}\lambda_{\gamma_5,+\dot c)}}{43200}\left[\begin{aligned}
&11\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-33\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+33\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]
\end{aligned}
\]

### Color tensor 8

\[(t_A)^{i_2}{}_{j_2}(t_B)^{j_2}{}_{k_2}\bar W^{i_1 k_2 i_3 i_4}\]

Full coefficient:

\[
\begin{aligned}
&\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_4,+\dot b}\lambda_{\gamma_4,+\dot c)}}{21600}\left[\begin{aligned}
&8\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle+8\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}+24\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+60\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-12\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_4,+\dot b}\lambda_{c,+\dot c)}}{21600}\left[\begin{aligned}
&-6\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}-\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle-4\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-14\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}+16\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+24\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-24\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+33\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+72\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-12\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_4,+\dot b}\lambda_{\gamma_5,+\dot c)}}{21600}\left[\begin{aligned}
&-14\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+75\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-18\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{c,+\dot b}\lambda_{c,+\dot c)}}{216000}\left[\begin{aligned}
&-60\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}+125\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle-100\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&+135\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}-175\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+200\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-264\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+45\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+150\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+120\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{c,+\dot b}\lambda_{\gamma_5,+\dot c)}}{7200}\left[\begin{aligned}
&5\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-8\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+34\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-4\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_5,+\dot b}\lambda_{\gamma_5,+\dot c)}}{21600}\left[\begin{aligned}
&-35\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+72\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-39\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{c,+\dot b}\lambda_{c,+\dot c)}}{216000}\left[\begin{aligned}
&-122\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}-87\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle-37\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&+150\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}-117\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+145\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-30\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-84\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+240\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-150\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{c,+\dot b}\lambda_{\gamma_5,+\dot c)}}{108000}\left[\begin{aligned}
&24\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-45\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-75\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+155\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{\gamma_5,+\dot b}\lambda_{\gamma_5,+\dot c)}}{21600}\left[\begin{aligned}
&-8\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-3\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+72\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-35\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_5,+(\dot a}\lambda_{\gamma_5,+\dot b}\lambda_{\gamma_5,+\dot c)}}{43200}\left[\begin{aligned}
&11\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-33\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-33\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]
\end{aligned}
\]

### Color tensor 9

\[(t_B)^{i_2}{}_{j_2}(t_A)^{i_4}{}_{j_4}\bar W^{i_1 j_2 i_3 j_4}\]

Full coefficient:

\[
\begin{aligned}
&\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_4,+\dot b}\lambda_{\gamma_4,+\dot c)}}{756000}\left[\begin{aligned}
&-20\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}+72\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle-60\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-135\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}+144\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-30\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+240\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-270\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+72\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_4,+\dot b}\lambda_{c,+\dot c)}}{504000}\left[\begin{aligned}
&-92\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}+285\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle-228\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-480\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}+435\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-60\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+780\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-600\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+15\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+36\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_4,+\dot b}\lambda_{\gamma_5,+\dot c)}}{50400}\left[\begin{aligned}
&-30\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+32\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-9\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{c,+\dot b}\lambda_{c,+\dot c)}}{100800}\left[\begin{aligned}
&-59\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}+160\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle-120\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-264\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}+192\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-12\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+156\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-216\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-48\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+20\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{c,+\dot b}\lambda_{\gamma_5,+\dot c)}}{100800}\left[\begin{aligned}
&-120\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+192\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-16\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-9\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_5,+\dot b}\lambda_{\gamma_5,+\dot c)}}{25200}\left[\begin{aligned}
&17\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-30\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+8\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{c,+\dot b}\lambda_{c,+\dot c)}}{75600}\left[\begin{aligned}
&-100\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}+270\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle-180\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-132\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}+342\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-24\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+9\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-9\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}+18\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{c,+\dot b}\lambda_{\gamma_5,+\dot c)}}{25200}\left[\begin{aligned}
&3\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-3\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}-4\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{\gamma_5,+\dot b}\lambda_{\gamma_5,+\dot c)}}{25200}\left[\begin{aligned}
&-9\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+9\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}+2\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_5,+(\dot a}\lambda_{\gamma_5,+\dot b}\lambda_{\gamma_5,+\dot c)}}{1680}\left[\begin{aligned}
&\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]
\end{aligned}
\]

### Color tensor 10

\[(t_B)^{i_2}{}_{j_2}(t_A)^{j_2}{}_{k_2}\bar W^{i_1 k_2 i_3 i_4}\]

Full coefficient:

\[
\begin{aligned}
&\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_4,+\dot b}\lambda_{\gamma_4,+\dot c)}}{302400}\left[\begin{aligned}
&24\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}-20\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle+72\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&+104\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}+12\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+57\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}\\
&+240\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+84\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_4,+\dot b}\lambda_{c,+\dot c)}}{302400}\left[\begin{aligned}
&78\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}-220\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle+172\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-48\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}-176\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-60\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-426\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-177\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+288\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+84\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_4,+\dot b}\lambda_{\gamma_5,+\dot c)}}{302400}\left[\begin{aligned}
&202\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+243\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+168\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{c,+\dot b}\lambda_{c,+\dot c)}}{1512000}\left[\begin{aligned}
&1116\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}-453\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle+1555\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&+3090\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}+375\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-1100\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+828\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+1665\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-480\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+570\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{c,+\dot b}\lambda_{\gamma_5,+\dot c)}}{151200}\left[\begin{aligned}
&-54\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-270\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+9\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+28\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_5,+\dot b}\lambda_{\gamma_5,+\dot c)}}{100800}\left[\begin{aligned}
&12\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+178\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+105\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{c,+\dot b}\lambda_{c,+\dot c)}}{1512000}\left[\begin{aligned}
&850\langle\lambda_{\gamma_4},\lambda_c\rangle^{3}-3690\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle+1005\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-960\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}-2886\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-465\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+360\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-1332\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+1860\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+520\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{c,+\dot b}\lambda_{\gamma_5,+\dot c)}}{252000}\left[\begin{aligned}
&-66\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+240\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-160\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-125\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{\gamma_5,+\dot b}\lambda_{\gamma_5,+\dot c)}}{100800}\left[\begin{aligned}
&26\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-92\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-10\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+77\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_5,+(\dot a}\lambda_{\gamma_5,+\dot b}\lambda_{\gamma_5,+\dot c)}}{50400}\left[\begin{aligned}
&-10\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+30\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+77\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
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
&\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_4,+\dot b}\lambda_{\gamma_4,+\dot c)}}{4900}\left[\begin{aligned}
&5\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle+4\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}+4\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&-5\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}-8\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-4\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_4,+\dot b}\lambda_{c,+\dot c)}}{4900}\left[\begin{aligned}
&10\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle-15\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-2\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\\
&-3\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-20\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}+2\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}\\
&-21\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-18\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}-3\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_4,+\dot b}\lambda_{\gamma_5,+\dot c)}}{1225}\left[\begin{aligned}
&-\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-2\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{c,+\dot b}\lambda_{c,+\dot c)}}{4900}\left[\begin{aligned}
&-10\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-10\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-9\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-6\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-21\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}+6\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{c,+\dot b}\lambda_{\gamma_5,+\dot c)}}{4900}\left[\begin{aligned}
&-10\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-11\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+2\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&+7\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_5,+\dot b}\lambda_{\gamma_5,+\dot c)}}{4900}\left[\begin{aligned}
&-5\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+4\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+4\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{c,+\dot b}\lambda_{c,+\dot c)}}{2450}\left[\begin{aligned}
&2\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}-6\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}+7\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{c,+\dot b}\lambda_{\gamma_5,+\dot c)}}{2450}\left[\begin{aligned}
&-10\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}+8\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{\gamma_5,+\dot b}\lambda_{\gamma_5,+\dot c)}}{4900}\left[\begin{aligned}
&-20\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+24\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}+\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_5,+(\dot a}\lambda_{\gamma_5,+\dot b}\lambda_{\gamma_5,+\dot c)}}{980}\left[\begin{aligned}
&\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]
\end{aligned}
\]

### Superpotential slot 2

\[\zeta(3)\bigl([t_A,t_B]-iC^C{}_{AB}t_C\bigr)^{i_2}{}_{j_2}\bar W^{i_1 j_2 i_3 i_4}\]

Full coefficient:

\[
\begin{aligned}
&\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_4,+\dot b}\lambda_{\gamma_4,+\dot c)}}{2450}\left[\begin{aligned}
&2\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}+7\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-6\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_4,+\dot b}\lambda_{c,+\dot c)}}{4900}\left[\begin{aligned}
&-10\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle-9\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}-10\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&+6\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}-21\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-6\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_4,+\dot b}\lambda_{\gamma_5,+\dot c)}}{2450}\left[\begin{aligned}
&8\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-10\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{c,+\dot b}\lambda_{c,+\dot c)}}{4900}\left[\begin{aligned}
&-15\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_5},\lambda_c\rangle+10\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-20\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\\
&-3\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-2\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}-3\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}\\
&-18\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-21\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}+2\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{c,+\dot b}\lambda_{\gamma_5,+\dot c)}}{4900}\left[\begin{aligned}
&7\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+2\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-11\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}\\
&-10\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_4,+(\dot a}\lambda_{\gamma_5,+\dot b}\lambda_{\gamma_5,+\dot c)}}{4900}\left[\begin{aligned}
&\langle\lambda_{\gamma_5},\lambda_c\rangle^{3}+24\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-20\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{c,+\dot b}\lambda_{c,+\dot c)}}{4900}\left[\begin{aligned}
&5\langle\lambda_{\gamma_4},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-5\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}+4\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\\
&+4\langle\lambda_{\gamma_4},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}-4\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle-8\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{c,+\dot b}\lambda_{\gamma_5,+\dot c)}}{1225}\left[\begin{aligned}
&-2\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}-\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{c,+(\dot a}\lambda_{\gamma_5,+\dot b}\lambda_{\gamma_5,+\dot c)}}{4900}\left[\begin{aligned}
&4\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle+4\langle\lambda_{\gamma_5},\lambda_c\rangle\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{2}-5\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle^{3}
\end{aligned}\right]\\
&+\frac{\lambda_{\gamma_5,+(\dot a}\lambda_{\gamma_5,+\dot b}\lambda_{\gamma_5,+\dot c)}}{980}\left[\begin{aligned}
&\langle\lambda_{\gamma_5},\lambda_c\rangle^{2}\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle
\end{aligned}\right]
\end{aligned}
\]

## Exact completeness

Every rational and period coefficient before and after the representation relation, every original geometric tensor, every source interpolation value, and all further full-source replays are retained in [jet_order_3_certificate.json](jet_order_3_certificate.json).

The all-finite-jet rational boundary evaluator independently reproduces every entry of this tensor at every proven-degree source interpolation node. Its proof and the strict scope are in [DERIVATION.md](DERIVATION.md) and [the source-simplex derivation](../forest-jet-identities/DERIVATION.md).
