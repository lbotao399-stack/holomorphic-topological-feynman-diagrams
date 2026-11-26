# An independent Grassmann proof of the determinant color polynomial

This derives the same finite-rank color polynomial as the cofactor proof in
[REPORT.md](REPORT.md), using an independent algebraic representation.
It supplies no independent bulk quantum amplitude.

Use $N$ odd variables $\psi_i$ and $N$ odd variables $\bar\psi_i$, with
integration orientation fixed by

$$
\int d^N\bar\psi\,d^N\psi\,e^{\bar\psi M\psi}=\det M.
$$

The even matrix with entries $\psi_i\bar\psi_j$ obeys

$$
(\psi\bar\psi)^2=(\bar\psi\psi)(\psi\bar\psi).
$$

Here $\bar\psi\psi=\sum_i\bar\psi_i\psi_i$ is an even scalar in the
Grassmann algebra. Left and right multiplication by this even matrix
commute. Expanding the iterated commutator and using its quadratic identity
therefore gives, for every $k\ge2$ and any constant even matrix $E$,

$$
\begin{aligned}
\operatorname{ad}_{\psi\bar\psi}^{k}(E)
={}&(\bar\psi\psi)^{k-1}
\left((\psi\bar\psi)E+(-1)^kE(\psi\bar\psi)\right)\\
&-\left(1+(-1)^k\right)(\bar\psi\psi)^{k-2}
(\psi\bar\psi)E(\psi\bar\psi).
\end{aligned}
$$

In this equation $\operatorname{ad}_{\psi\bar\psi}(E)$ means
$(\psi\bar\psi)E-E(\psi\bar\psi)$. The intermediate binomial sum is

$$
\sum_{j=1}^{k-1}(-1)^j\binom{k}{j}=-1-(-1)^k,
$$

which fixes the coefficient of the final term without a finite-rank
assumption.

## Trace-dual completeness and the signs

With $\operatorname{Tr}(t_A t_B)=\tau\delta_{AB}$, completeness gives

$$
\sum_A t_A(\bar\psi t_A\psi)
=-\tau\psi\bar\psi-\frac\tau N(\bar\psi\psi)\mathbf1.
$$

The central term drops from every commutator. A chain of $k$ insertions
therefore has the sign $(-1)^k$ after dividing out $\tau^k$.
Set $E=E_{ab}$ and test against $E_{ba}$, where $a\ne b$. The two end
terms give the diagonal entries $(\psi\bar\psi)_{aa}$ and
$(\psi\bar\psi)_{bb}$. The middle term gives their product.
For diagonal $M$,

$$
\begin{aligned}
\int e^{\bar\psi M\psi}(\psi\bar\psi)_{aa}
&=-\partial_{\mu_a}\det M,\\
\int e^{\bar\psi M\psi}
(\psi\bar\psi)_{aa}(\psi\bar\psi)_{bb}
&=\partial_{\mu_a}\partial_{\mu_b}\det M.
\end{aligned}
$$

Every factor $\bar\psi\psi$ differentiates the common scalar shift of
$M$. Explicitly this derivative is $\sum_i\partial_{\mu_i}$. Thus the
color polynomial of the main report is exactly

$$
\begin{aligned}
\mathcal C_k=(-1)^k\Bigg[&-
\left(\sum_i\partial_{\mu_i}\right)^{k-1}
\left(\partial_{\mu_a}+(-1)^k\partial_{\mu_b}\right)\prod_i\mu_i\\
&-\left(1+(-1)^k\right)
\left(\sum_i\partial_{\mu_i}\right)^{k-2}
\partial_{\mu_a}\partial_{\mu_b}\prod_i\mu_i\Bigg].
\end{aligned}
$$

Let $e_j$ have the meaning fixed in the main report: the elementary
symmetric polynomial in the entries other than $\mu_a,\mu_b$, zero outside
its degree range. Performing every derivative yields

$$
\boxed{
\mathcal C_k=(k-1)!\left((-1)^{k-1}\mu_b-\mu_a\right)e_{N-k-1}
-\left(1+(-1)^k\right)k(k-2)!e_{N-k}.
}
$$

All expressions are polynomial. Singular matrices require no limiting
prescription or inverse matrix.

## Independent implementation

`determinant_grassmann_check.py` implements the full exterior algebra on
$2N$ generators using sparse integer coefficients. It constructs the
matrix $\psi\bar\psi$, multiplies its entries in that algebra, and
iterates the actual commutator. It then integrates the resulting
polynomial against the complete Gaussian. It imports neither the
cofactor evaluator nor the boundary fan engine.

The committed certificate verifies 60 complete matrix identities,
162 Gaussian integrals at ranks 2 through 7, including zero eigenvalues,
and 78 exact rank-threshold checks through $n=40$. The proof above applies
to every $N$ and every $k\ge2$.
