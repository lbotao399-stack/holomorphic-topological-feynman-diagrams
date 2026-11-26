# Three-loop $K_4$ periods in higher holomorphic dimension

At zero displacement the three-loop $T=0$ graph on $\mathbb C^4$ with the ordered oriented edges

$$
((0,4),(0,10),(1,2),(1,7),(1,10),(2,3),(2,8),(3,6),(4,9),(5,6),(7,5),(8,4),(9,5))
$$

and sources $\lambda_0=e_1$, $\lambda_3=e_2$, $\lambda_6=e_3$, $\lambda_7=\lambda_8=\lambda_9=e_4$, $\lambda_{10}=-e_1-e_2-e_3-3e_4$ has

$$
K_\Gamma(\lambda;0)=\frac{2592\zeta(3)-3119}{83980800}.
$$

Here $e_i$ are the coordinate unit vectors, and all unlisted sources vanish. The graph satisfies the $H=4$ subgraph inequality $3|I_S|\le 4|S|-5$ for each of its 1,874 edge-containing induced vertex sets $S$. Fourteen sets attain equality. Irrationality follows from Apéry's theorem.

**Cotree and sign.** Take [case 0](k4exact.py), the $H=3$ graph with skeleton chain lengths $[1,1,1,1,3,3]$. Its oriented zero-based edge list is

$$
((0,4),(0,7),(1,2),(1,5),(1,7),(2,3),(2,4),(3,6),(4,5),(5,6)).
$$

Select the three single-edge chains $I=(3,6,8)$ in that list. Their complement is a spanning tree. The positive cycle matrix $C$ in the [project orientation](../../code/htk/CONVENTIONS.md) has $\det C_I=-1$. The [cotree subdivision formula](../ht-families/dimension_raise.md) inserts vertices $7,8,9$ immediately before the old root, now $10$, and puts second segments at the end of the list. It gives the parent integral multiplied by the signed weight

$$
-\frac{t_3t_6t_8}{U},\qquad U=\det(C\,\mathrm{diag}(t)\,C^{\mathsf T}).
$$

**Star coordinates.** Use the parent chart centered at vertex $1$ with outer triangle $(A,B,C)=(2,4,5)$. Its arm conductances are $g_A,g_B,g_C$, summing to one, and its triangle side opposite $i$ has conductance $g_jg_k r_i/(1-r_i)$, with $\{i,j,k\}=\{A,B,C\}$. Set $Q=1-g_A r_A-g_B r_B-g_C r_C$. Suppressing degree-two vertices preserves $U$ when chain resistances are summed. Since no edge in $I$ is subdivided in the parent, its unsigned weight is the product of its complement's skeleton conductances divided by the skeleton spanning-tree polynomial. Eliminating the star center gives the three effective triangle conductances $g_jg_k/(1-r_i)$, and hence

$$
\begin{aligned}
\sum_{\substack{T\ \mathrm{spanning}\ \mathrm{tree}}}
\prod_{e\in T}\mathrm{conductance}_e
 &=\frac{g_Ag_Bg_C Q}{(1-r_A)(1-r_B)(1-r_C)},\\
\prod_{\substack{e\notin I\\\mathrm{skeleton~edges}}}\mathrm{conductance}_e
 &=g_Ag_B\,\frac{g_Ag_Cr_B}{1-r_B},\\
\frac{t_3t_6t_8}{U}&=\frac{g_A r_B(1-r_A)(1-r_C)}{Q}.
\end{aligned}
$$

In the chart centered at vertex $2$ with outer triangle $(1,4,5)$, the weight is $g_A r_C(1-r_A)(1-r_B)/Q$. Both weights are independent of the chain fractions. Thus the parent chain-integrated numerator $P$ can be reused directly in each chart. Its derivation from the ten Kirchhoff edge currents is [replayed here](verify/graph_k4_series.py).

**Iteration.** After each subdivision, select its three just-appended second segments as the next cotree. At stage $j=0,1,\ldots$ the parent has $H=3+j$. With the orientation calibrated, the new cotree minor is the old minor multiplied by $(-1)^3$. The sign at stage $j$ is consequently $(-1)^{3(4+j)}(-1)^{j+1}=-1$. The period matrix remains fixed. At depth $j$ each selected edge time is its original time multiplied by all its preceding second-segment fractions. For $m$ subdivisions, the fractions from step $k$ therefore each occur to power $m-k$. Integration over the three unit intervals at each step gives

$$
\prod_{k=1}^m\left(\int_0^1 x^{m-k}\,dx\right)^3=(m!)^{-3}.
$$

With $P/Q^{10}$ the oriented parent density after the four chain integrations, the exact remaining period at $H=3+m$ is

$$
\frac{(-1)^m}{(m!)^3}
\int_{\substack{g_A,g_B>0\\ g_A+g_B<1}}dg_A\,dg_B
\int_{(0,1)^3}dr_A\,dr_B\,dr_C\;
\frac{P\,[g_A r_B(1-r_A)(1-r_C)]^m}{Q^{10+m}}.
$$

The unsigned weight is a single spanning-tree probability and lies in $[0,1]$, so absolute convergence follows from the parent integral. The [exact beta-series replay](verify/raised_k4_series.py) and its [eight evaluations](verify/out_raised_k4_series_exact.log) record the coefficients below. The replay expands the new numerator, integrates the two star charts independently, and checks densities at three rational points per chart against the unreduced graph in every dimension shown. All eight evaluations have strict corner convergence margin one.

| $H$ | coefficient of $1$ | coefficient of $\zeta(3)$ |
|---:|---:|---:|
| 4 | $-3119/83980800$ | $1/32400$ |
| 5 | $-1127773/724250419200$ | $181/139708800$ |
| 6 | $102037909/4519322615808000$ | $-131/6974263296$ |
| 7 | $-55302943/1220217106268160000$ | $71/1883051089920$ |

~~~sh
python3 verify/raised_k4_series.py 1 2 3 4
~~~
