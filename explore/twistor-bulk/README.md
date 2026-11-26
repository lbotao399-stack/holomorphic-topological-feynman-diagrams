# Twistor bulk boxes

Use Bittleston–Sharma–Skinner's flat twistor patch, the Poisson bracket on its two fiber coordinates, and their normalized heat kernel. In a fixed cyclic order of the four conserved holomorphic momenta, the current on an edge is the cumulative momentum through its tail minus the Schwinger-time weighted average of all four cumulative momenta. This solves conservation and makes the weighted sum of loop currents zero. The normalized times $t_e$ lie on the unit three-simplex. The signed scalar box measure supplies the three-dimensional determinant of the first three momenta. Alternation of the four fields cancels the sign of this determinant under reordering.

At a gauge insertion take the adjoint color matrix. At a gravitational insertion of momentum $p_i$ take the two-dimensional Poisson determinant $[p_i,y_{i,i+1}]$ on the fiber; replacing the outgoing by the incoming current leaves it unchanged. Each such factor is linear in the four simplex coordinates. A monomial of powers $m_e$ integrates to $\prod_{e}m_e!/(3+\sum_e m_e)!$. The replay performs this integration in the three independent fiber brackets after total momentum conservation. There are six cyclic orderings with the first input fixed.

| Physical ghost and three fields | Sum of the cyclic box symbols, after the common determinant |
|---|---|
| Four gauge legs | $1/6$ times the sum of the six adjoint traces |
| One gravitational leg | Zero separately in each of the two cyclic gauge color orders |
| Two gravitational legs $i,j$ | $-2[p_i p_j][p_j p_i]/24$ times the adjoint gauge trace |
| Three gravitational legs | Zero by the one-generator adjoint trace |
| Four gravitational legs | $1/120$ times the sum of the six cyclic products $[p_i p_j][p_j p_k][p_k p_l][p_l p_i]$ |

For a plane wave the connection $s$ in the source is the negative fiber Hessian of its Hamiltonian with one index raised by the supplied Poisson tensor. The products in the table are precisely its two- and four-matrix traces. In the mixed channel two permutations polarize the two like physical fields. The pure gravity and gauge channels each have six. Thus in the source's common heat normalization $(i/(2\pi))^3$ the coefficient of $\operatorname{tr}_{\rm ad}(c(\partial a)^3)$ is $1/3!$. The gauge and gravity variations of the mixed term each have coefficient $-1/24$. The coefficient of $\operatorname{tr}(\psi(\partial s)^3)$ is $(1+\dim\mathfrak g)/5!$, since each adjoint component and the gravitational system give the same four-Hamiltonian loop. These are bulk graphs with external bulk fields. [Bittleston–Sharma–Skinner, §§3–4 and Appendix B](https://arxiv.org/abs/2208.12701) · [Exact replay](box.py)
