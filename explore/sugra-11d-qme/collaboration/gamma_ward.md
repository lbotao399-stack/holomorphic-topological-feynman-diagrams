# Holomorphic longitudinal $\gamma$ decoupling in the 11d two-loop theta integrand

In the twisted 11d supergravity model on $\mathbb R\times\mathbb C^5$, set $K=0$ (no external $\nu$). On each theta graph, in each admissible flow, every external holomorphic one-form $\gamma$ occurs through its own $\lambda_v\wedge\gamma$, where $\lambda_v$ is its holomorphic momentum. Replacing any one polarization by $\gamma=\lambda_v$ gives zero pointwise in the Schwinger parameters. The statement concerns the Wick contraction summed over letter assignments at each vertex. The integrated obstruction in the quantum master equation remains undetermined.

## The two possible vertices

Write $m$ for the number of cubic $\gamma$ vertices $c$ and $g$ for their coupling. The other $12-m$ vertices are of type $b_0$. In the source action their terms are, with $a,b,c,d,e$ abstract holomorphic tangent or cotangent slots on $\mathbb C^5$ and $\varepsilon^{abcde}$ the source volume multivector,

$$
\mu^a\mu^b\partial_a\gamma_b,\qquad
\frac g6\,\varepsilon^{abcde}\gamma_a\partial_b\gamma_c\,\partial_d\gamma_e.
$$

The odd fields $\mu$ contract with $\gamma$ along 13 edges. A flow orients each edge from its $\mu$ letter to its $\gamma$ letter. Thus the external multiplicities are

$$
\#\mu=2(12-m)-13=11-2m,\qquad
\#\gamma=(12-m)+3m-13=2m-1,\qquad 1\le m\le5.
$$

Each trivalent end consumes all three of its letters. Each of the ten bivalent interior vertices has exactly one external letter. The source defines $\lambda_v$ as the total momentum of all external letters at $v$. Here a vertex carrying an external $\gamma$ has precisely that one letter, so its individual momentum equals $\lambda_v$. Such a vertex is either a $b_0$ with two outgoing edges or a $c$ with two incoming edges.

At a $b_0$ vertex, the ordered left derivatives with respect to the two contracted $\mu$'s in the project's edge order give the two-covector

$$
-\bigl((\lambda_v)_a\gamma_b-(\lambda_v)_b\gamma_a\bigr)
=-(\lambda_v\wedge\gamma)_{ab}.
$$

At an interior $c$ vertex, let $\gamma_1,\gamma_2$ be the covectors inserted in its two contracted slots in edge order. The derivative momenta $J_1,J_2$ in those slots equal $+y_e$ at $e(1)$ and $-y_e$ at $e(0)$ in the source's fixed circuit orientation. The two signed currents and the external momentum obey

$$
J_1+J_2+\lambda_v=0 .
$$

Use $\varepsilon(\gamma_1,\gamma_2,\gamma,\cdot,\cdot)$ for evaluation of the same volume multivector on its five covector arguments. In the $c$ vertex the two differentiated $\gamma$'s contribute equally under exchange of their two pairs of indices. Moving the two derivative slots past the polarization slots requires three transpositions. The six Wick terms therefore give

$$
\begin{aligned}
&-\frac g3\bigl[
\varepsilon(\gamma_1,\gamma_2,\gamma,J_2,\lambda_v)
+\varepsilon(\gamma_1,\gamma_2,\gamma,\lambda_v,J_1)
+\varepsilon(\gamma_1,\gamma_2,\gamma,J_1,J_2)\bigr]\\
&=-\frac g3\,\varepsilon(\gamma_1,\gamma_2,\gamma,J_2-J_1,\lambda_v-J_1)\\
&=-g\,\varepsilon(\gamma_1,\gamma_2,\gamma,J_1,J_2)
=g\,\varepsilon(\gamma_1,\gamma_2,\gamma,J_1,\lambda_v).
\end{aligned}
$$

Both tensors vanish for $\gamma=\lambda_v$. All other vertex factors, the flow sign and the scalar theta kernel multiply these tensors. Each of the $2m-1$ external $\gamma$'s can therefore be differentiated explicitly. After the displayed substitution the number of internal derivative-current slots is

$$
(12-m)+2m-(2m-1)=13-m .
$$

## External $\nu$ with arbitrary momentum

The same longitudinal zero holds with any number of external $\nu$ fields. The exclusion of the bivalent vertex $\beta\,\partial_a\mu^a$ leaves every $\nu$ external. At a $b_k$ vertex write $p$ for the momentum of the individual external $\gamma$. The two contracted $\mu$ slots give

$$
-\nu^k(p_a\gamma_b-p_b\gamma_a).
$$

The derivative in $\nu^k\mu^a\mu^b\partial_a\gamma_b$ acts on $\gamma$ alone. The circuit uses the **total** external momentum $\lambda_v$ at this vertex, including all $\nu$ momenta. Multiplying by that circuit polynomial preserves the zero at $\gamma=p$. A $c$ vertex carries no $\nu$, so its sole external $\gamma$ still has momentum $\lambda_v$ and the preceding calculation applies outright. This proof holds for arbitrary $k$. The enumerations below test $K=0$. [External $\nu$ and the vertex selection rule](constant_nu.md)

## Exact replay

Run `python3 check_gamma_ward.py` from this directory. For each of $\theta_{166},\theta_{256},\theta_{346},\theta_{355},\theta_{445}$ the script enumerates all $2^{13}$ edge orientations. The numbers of admissible flows at $m=1,2,3,4,5$ are $(40,480,1008,480,40)$. Across the five graphs there are 10,240 admissible flows and 51,200 external $\gamma$ occurrences. At one exact rational configuration per graph with vanishing endpoint momenta and nonzero scalar determinant $D_\theta$, the literal Wick-word evaluator agrees in all 2,500 distinct site, vertex-type and component comparisons. The same 2,500 entries vanish for longitudinal polarization, and 1,896 are nonzero with the sampled unrestricted polarizations.

The action, edge signs and Kirchhoff equation are in the [11d note](../../../../研究项目_非Git/Twisted_11d_sugra_loops_2026-11-18/main.tex), equations labelled `R-action`, `eq-Delta` and the holomorphic Kirchhoff equation following `eq-extraction`. The existing [sector classification](../VERIFY.md) supplies the five topologies and admissible vertices. The [Wick-word evaluator](../verify/v_integrand.py) supplies the independent component calculation.
