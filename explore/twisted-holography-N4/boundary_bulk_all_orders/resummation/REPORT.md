# Summing every boundary correction in one protected-product measurement

## 1. Fixed conventions and the observable

Use the manuscript's fields, coupling $g$, spinor derivative
$\partial_{+\dot a}$, and fundamental trace
$\operatorname{Tr}(t_A t_B)=\tau\delta_{AB}$.
The four single-graviton components used in the
[boundary calculation](../boundary/REPORT.md) are

$$
\begin{aligned}
G_1&=\operatorname{Tr}(\beta_1\gamma^2),&
G_2&=\operatorname{Tr}(\beta_2\gamma^3),\\
G_3&=\operatorname{Tr}(\beta_3\gamma^1),&
G_4&=\operatorname{Tr}(b\gamma^1-ig\beta_2\beta_3).
\end{aligned}
$$

Each is annihilated by the complete boundary supercharge.
The first three are odd and $G_4$ is even. For $n\ge2$, define the
ordinary product

$$
\mathcal O^{(n)}_0=G_1G_2G_3G_4^{n-2}.
$$

The superscript labels the family and the subscript zero specifies the
ordinary, uncorrected product. The repeated $G_4$ factors have no
divided-power normalization.

The measurement sets $\beta_I=b=0$ and
$\gamma^1=\gamma^2=\gamma^3=X$, where $X$ is a constant traceless
diagonal matrix. Choose two opposite matrix units $E_{12},E_{21}$,
and write $\Delta$ for the difference between the first two eigenvalues
of $X$. The two ghost probes have independent odd amplitudes
$\varepsilon_+,\varepsilon_-$, colors $E_{12},E_{21}$, and
$\partial_{+\dot a}$ eigenvalues $p_{+\dot a},q_{+\dot a}$.
Their contracted spinor area is

$$
[p,q]=p_{+\dot a}q_+{}^{\dot a}.
$$

All displayed measurements take the ordered
$\varepsilon_+\varepsilon_-$ coefficient at the spacetime origin.
Only the scalar background is diagonal; the output ghosts remain
off diagonal. This is compatible with the manuscript's vanishing
of quantum images when every color is restricted to the Cartan.

The all-order boundary result is

$$
\left.Q_k\mathcal O^{(n)}_0\right|_{\varepsilon_+\varepsilon_-}
=0\quad(k\ne n),\qquad
\left.Q_n\mathcal O^{(n)}_0\right|_{\varepsilon_+\varepsilon_-}
=-\frac{\Delta^{n+1}}{(n!)^2}
\left(\frac{ig\tau[p,q]}{32\pi^2}\right)^n.
$$

The equality for $k=0$ also follows from classical closure of every
factor. The vertical evaluation notation always includes the
specified field specialization and the ordered ghost coefficient.
The full cochain $Q_k\mathcal O^{(n)}_0$ is not asserted to vanish
away from this measurement.

For completeness, a separate pair of protected constituents detects
$Q_1$:

$$
\mathcal O^{(1)}_0=
\operatorname{Tr}((\gamma^1)^2)\operatorname{Tr}(\beta_1\gamma^2),
\qquad
\left.Q_1\mathcal O^{(1)}_0\right|_{\varepsilon_+\varepsilon_-}
=\frac{2\tau g^2\Delta^2[p,q]}{32\pi^2}.
$$

Every $Q_k$ with $k\ge2$ vanishes on this pair by its total odd-coordinate
degree. Thus the family includes every positive boundary order.

## 2. Exact generating function

Let $\hbar$ record boundary perturbative order and let $u$ count the
even constituent $G_4$. Both are formal variables. The ordinary
generating observable is

$$
\frac{G_1G_2G_3}{1-uG_4}
=\sum_{n\ge2}u^{n-2}\mathcal O^{(n)}_0.
$$

The loop-isolation result makes its complete supercharge measurement
diagonal in the two gradings:

$$
\boxed{
\left.
\mathbf Q\!\left(\frac{G_1G_2G_3}{1-uG_4}\right)
\right|_{\varepsilon_+\varepsilon_-}
=-\Delta\sum_{n\ge2}
\frac{1}{(n!)^2}
\left(\frac{ig\tau\hbar\Delta[p,q]}{32\pi^2}\right)^n
u^{n-2}.
}
$$

The standard modified Bessel function of order zero is specified here
by its entire series

$$
I_0(2\sqrt{x})=\sum_{n\ge0}\frac{x^n}{(n!)^2}.
$$

The expression is independent of the choice of square root.
This is [the standard definition](https://dlmf.nist.gov/10.25.E2)
specialized to order zero. Consequently the exact sum is

$$
\boxed{
-\frac{\Delta}{u^2}
\left[
I_0\!\left(2\sqrt{\frac{ig\tau\hbar u\Delta[p,q]}{32\pi^2}}\right)
-1-\frac{ig\tau\hbar u\Delta[p,q]}{32\pi^2}
\right].
}
$$

The apparent $u^{-2}$ has no pole in the formal result: both lower
powers have been subtracted. This formula is a resummation of a
boundary cochain measurement. It is not an independently evaluated
bulk partition function.

For the standard two-by-two block
$X=\operatorname{diag}(1/2,-1/2)$, $\tau=1/2$, and $[p,q]=1$,
the exact first coefficients are

| Boundary order | Measured coefficient |
|---:|---|
| $1$ | $g^2/(32\pi^2)$ on $\mathcal O^{(1)}_0$ |
| $2$ | $g^2/(16384\pi^4)$ |
| $3$ | $ig^3/(9437184\pi^6)$ |
| $4$ | $-g^4/(9663676416\pi^8)$ |
| $5$ | $-ig^5/(15461882265600\pi^{10})$ |
| $n\ge2$ | $-(ig)^n/((64\pi^2)^n(n!)^2)$ |

All orders $n\ge2$ arise in the same representation and with the same
root gap and spinor area. No large-rank limit is needed for the
nonvanishing statement.

## 3. What the quantum composite dictionary must reproduce

Replace the ordinary product by a fixed bracketing of the manuscript's
marked products,

$$
G_1\star G_2\star G_3
\underbrace{\star G_4\star\cdots\star G_4}_{n-2\ \text{factors}}
=\sum_{j\ge0}\hbar^j\mathcal O^{(n)}_j.
$$

For every fixed $n$ this is a finite expansion, by the inherited
multi-graviton completion theorem. Its zeroth term is the ordinary
product above. Every factor is quantum closed, so the full product
is quantum closed. At total order $n$ this yields the cochain identity

$$
Q_n\mathcal O^{(n)}_0+
\sum_{j=1}^{n}Q_{n-j}\mathcal O^{(n)}_j=0.
$$

The new explicit consequence for every $n\ge2$ is

$$
\boxed{
\left.
\sum_{j=1}^{n}Q_{n-j}\mathcal O^{(n)}_j
\right|_{\varepsilon_+\varepsilon_-}
=\frac{\Delta^{n+1}}{(n!)^2}
\left(\frac{ig\tau[p,q]}{32\pi^2}\right)^n.
}
$$

Thus some composite-dictionary correction must contribute to this
measurement at every order. The individual terms need not be
nonzero separately, and the identity does not determine them.
Summing these necessary correction terms against
$\hbar^n u^{n-2}$ gives the negative of the Bessel expression in
Section 2.

The [bulk report](../bulk/REPORT.md) identifies the four protected
one-particle states and proves that local primary trees on their
linear images give zero. The explicit nonzero sum above therefore
quantifies the contribution that must be supplied by the complete
composite-observable map or external interaction kernels in that
comparison. It does not imply a lifting of the protected products.
Their nonexactness has a separate evaluation: put every field on a
constant Cartan direction, with $\gamma^I=X$, $\beta_I=\chi_I X$,
$b=vX$, and $c=0$, where the $\chi_I$ are independent odd amplitudes.
This evaluation annihilates every differential image. Its value on
the zeroth product is
$(\operatorname{Tr}X^2)^{n+1}\chi_1\chi_2\chi_3v^{n-2}\ne0$.
Thus the closed marked completions survive, including after
inverting $\hbar$. The boundary report proves this independently of
the off-diagonal ghost measurement.

## 4. Independent certificate

Run:

```bash
python explore/twisted-holography-N4/boundary_bulk_all_orders/resummation/all_order_response.py
```

The script imports no sibling evaluator. It multiplies the original
fan, Berezin, trace-conversion, and assignment factors using Gaussian
rational arithmetic, then evaluates both ordered root colors by
literal matrix commutators. It checks 299 phase, root-gap, and
spinor-area identities through boundary order 24.

Independently, the equation
$x\,d^2/dx^2+d/dx-1$ determines the entire series from its constant
coefficient one. Writing $a_n$ for its scalar series coefficient,
the recursion is $(n+1)^2a_{n+1}=a_n$.
The script constructs these coefficients by the
recursion and verifies the factorial formula through degree 64.
This check is independent of the finite list of fan coefficients;
the all-order equality follows analytically from the displayed series.

The one-loop normalization is supplied by the separate independent
calculation in the boundary certificate. The present script records
it as an input and does not claim to rederive it.

## Provenance

The fan and component conventions are inherited from
[`sources/LLC_main.tex`](https://github.com/lbotao399-stack/holomorphic-topological-feynman-diagrams/blob/29ef4be3c92b8055ac51067357524c1048e6cda0/sources/LLC_main.tex),
equations `eq-Qnzero` and `eq-A-16`. The all-order isolation, normalization,
and protected constituent test are proved in the sibling boundary
report. The generating function and composite-dictionary sum rule
are the additional results here.
