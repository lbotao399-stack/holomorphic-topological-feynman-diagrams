# Sharp compact-momentum support in the eleven-dimensional theta kernel

This calculation sharpens the necessary Kaluza–Klein support bound in
[the existing eleven-dimensional notes](../../explore/sugra-11d-qme/collaboration/current_complex/README.md).
It concerns the geometric theta kernel before its supergravity tensor
contractions and the sum over graph shapes.

## Exact rank with prescribed routing

For the three oriented paths $A,B,C$, the two-loop holomorphic
determinant has columns
$$
 (p_i,0)\quad(i\in A),\qquad
 (p_i,p_i)\quad(i\in B),\qquad
 (0,p_i)\quad(i\in C).
$$
There are two copies of the five-dimensional holomorphic cotangent
space. For compactification along a complex torus of dimension $r$,
evaluate each momentum on its compact tangent directions; in the
following rank formula $p_i$ denotes these compact covectors.

The rank of the resulting compact columns is exactly
$$
\begin{aligned}
&\dim\operatorname{span}\{p_i:i\in A\}
+\dim\operatorname{span}\{p_i:i\in B\}
+\dim\operatorname{span}\{p_i:i\in C\}\\
&\hspace{8mm}
-\dim\Bigl(
\operatorname{span}\{p_i:i\in A\}
\cap\operatorname{span}\{p_i:i\in B\}
\cap\operatorname{span}\{p_i:i\in C\}\Bigr).
\end{aligned}
$$

To prove this, map the direct sum of the three indicated spans to
the two compact cotangent copies by
$$
 (a,b,c)\longmapsto(a+b,b+c).
$$
Its kernel consists exactly of $(-b,b,-b)$ with $b$ in the
triple intersection. Rank–nullity gives the formula.

The noncompact columns span at most $2(5-r)$ dimensions.
Consequently a nonzero full determinant requires compact rank $2r$,
and therefore at least $2r$ distinct external slots with nonzero
compact momentum. At this minimum, the compact covectors within
each path must be independent and the three path spans must have
zero common intersection. Merely having $2r$ compact excitations
is not sufficient if their routing violates this condition.
For $r=1$, this reproduces the requirement that two nonzero
compact momenta occur on different paths.

## The lower bound is attained for every compact dimension

Take the shape $(1,6,6)$. Let $e_1,\ldots,e_5$ be the ordered
unit holomorphic covectors. Put
$$
 (p_i)_{i\in B}=(e_1,e_2,e_3,e_4,e_5),\qquad
 (p_i)_{i\in C}=(-e_1,-e_2,-e_3,-e_4,-e_5).
$$
The path $A$ has no interior input. Momentum conservation holds.
With the displayed ordering,
$$
 D_\theta
 =\det\begin{pmatrix}\mathbf1_5&0\\
                     \mathbf1_5&-\mathbf1_5\end{pmatrix}
 =-1.
$$
Declare any $r$ of these five directions compact. Exactly $2r$
external slots then carry compact momentum, the compact rank is
$2r$, and the full determinant remains $-1$.

Let $\Gamma$ be the Euler gamma function, as in the original
theta-integral calculation. In the original convention whose
path density starts with $D_\theta/2$, the geometric integral is
$$
\begin{aligned}
K_{166}
&=\frac{D_\theta}{2}
 \frac{\Gamma(11/2)\Gamma(1/2)^2}
      {\Gamma(13/2)(5!)^2}\\
&=-\frac12
 \frac{(945/32)\sqrt\pi\,\pi}
      {(10395/64)\sqrt\pi\,(120)^2}
=-\frac{\pi}{158400}.
\end{aligned}
$$
| Compact complex dimension | Nonzero compact external slots | Compact rank | Full determinant | $K_{166}/\pi$ |
|---:|---:|---:|---:|---:|
| 1 | 2 | 2 | −1 | $-1/158400$ |
| 2 | 4 | 4 | −1 | $-1/158400$ |
| 3 | 6 | 6 | −1 | $-1/158400$ |
| 4 | 8 | 8 | −1 | $-1/158400$ |
| 5 | 10 | 10 | −1 | $-1/158400$ |

Thus the geometric support bound is sharp, including momentum
conservation. The statement does not determine the complete
supergravity coefficient at these inputs: polarizations, the
vertex contractions, and all five theta shapes can impose
additional cancellations.

For actual smooth torus Fourier modes, compact holomorphic
derivatives supply the corresponding complex Fourier momenta.
Multiplying all displayed momenta by the common Fourier
normalization multiplies this degree-ten geometric kernel by
its tenth power. The rank and exact support statements are
unchanged. This is a test of the local vertex on external modes;
it does not integrate out massive compact modes.

## Certificate and provenance

`checks/compact_momentum_rank.py` checks the rank formula in
120 rational examples with two independent linear-algebra
constructions, and verifies all five sharp-support examples
and the gamma-function integral exactly.

The inherited determinant and kernel conventions are in
[VERIFY.md](../../explore/sugra-11d-qme/VERIFY.md), section b2,
and the support bound is in the current-complex notes linked above.
The new contributions here are the exact routing criterion and
the simultaneous sharpness certificate in every compact dimension.
