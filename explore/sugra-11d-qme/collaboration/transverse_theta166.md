# Four orientations with divergence-free external vectors

Use the same momenta and fixed external one-form as in the ordered $\theta_{166}$ calculation. In the coordinate basis and its dual,

$$
p^B=(e^1,e^2,e^3,e^4,e^5),\qquad
p^C=-(e^1,e^2,e^4,e^3,e^5),\qquad
\lambda_N=\lambda_{N'}=0,\quad D_{166}=1,\quad \gamma_{B_1}=e^2 .
$$

Here $p^B,p^C$ are lists of the five interior momenta along each path. The $\mu$ polarizations at $C_1,\ldots,C_5$ are $(e_2,e_1,e_1,e_1,e_1)$. At each of the nine vector insertions the polarization lies in the kernel of its own momentum, $p_a\mu^a=0$. The four admissible edge orientations consist of one with the cubic $\gamma$ vertex at $N$, one at $N'$, and two at $B_1$. [Edge and vertex conventions](theta166_flow.md)

| $\mu$ polarizations at $B_2,B_3,B_4,B_5$ | Sum of the four integrated kernels, divided by $\pi$ |
|---|---:|
| $(e_3,e_1,e_1,e_1)$ | $-25457/278891278848000$ |
| $(e_3,e_4,e_5,e_2)$ | $11120243/48381574026240000$ |

Both values strip the coupling $g$ and use the project's $\mathcal K_\theta$ normalization. The calculation contracts the four oriented vertex tensors over the ten first-edge current components and integrates their sum with the source Dirichlet measure. Each sum contains 72 monomials. Substitution of currents from a Kirchhoff solve at two positive edge metrics reproduces both polynomials' values. [Symbolic contraction and integral](check_transverse_theta166.py)

For the first row, exchanging directions $4$ and $5$ leaves all polarizations fixed. The four affected vector insertions are permuted by two transpositions. The holomorphic volume tensor changes sign, so the fully alternating evaluation is zero. [Permutation proof](theta166_four_flows.md)

For the second row, any linear symmetry of the momentum set preserving the $\gamma$ leg permutes the five coordinate lines with signs, and fixes its momentum $e^1$. At each of the other four directions the vector polarizations on the positive and negative momenta differ: all four negative momenta carry $e_1$, and none of the positive ones does. A signed permutation preserving polarization lines cannot exchange the two signs. The $\gamma$ polarization fixes line $2$. The vector at $+e^2$ then fixes line $3$, that at $+e^3$ fixes line $4$, and that at $+e^4$ fixes line $5$. Only the identity remains. The alternating evaluation on this input is still open.

An independent polynomial integration in path ratios and the two successive conductance fractions reproduces both values. Its enumeration over signed coordinate permutations gives stabilizer orders $16$ and $1$, respectively. [Independent verifier](verify_transverse_theta166.py)

For the antisymmetrization with $\gamma$ fixed at $B_1$, $D_{166}$ factors as the two five-column path determinants. Each path must receive one momentum on each coordinate line. On $B$ the first line is fixed; choosing the signs of the other four lines and ordering their momenta gives $2^4\,4!$ assignments. The five remaining momenta on $C$ have $5!$ orders. Thus $2^4\,4!\,5!=46080$ of the $9!=362880$ labelled assignments pass the determinant test. Further zeros can come from vertex contractions.
