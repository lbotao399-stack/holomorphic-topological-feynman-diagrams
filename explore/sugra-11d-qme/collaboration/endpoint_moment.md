# The two-current insertion at the trivalent endpoint

Orient the three paths $A,B,C$ of $\theta_{pqr}$ from $N$ to $N'$. At an empty endpoint $N$, the injected holomorphic momentum is $\lambda_N=0$. The incoming symbols of a cubic $\gamma\partial\gamma\partial\gamma$ vertex at $N$ are $-\ell_A,-\ell_B,-\ell_C$, the negatives of the three initial path currents. Their sum vanishes. After contracting the volume tensor with the three polarizations, the two remaining derivative slots select the exterior product $\ell_A\wedge\ell_B$. [Cubic vertex](gamma_ward.md) · [Source, series reduction and theta insertion](../../../../研究项目_非Git/Twisted_11d_sugra_loops_2026-11-18/main.tex)

The resistance fractions $x^P_j$ on each path sum to one, as do the three conductance fractions $w_P$. The path lengths are $n_A=p,n_B=q,n_C=r$. Write $Q^P_j$ for the sum of the interior momenta preceding edge $j$ on path $P$. The source's current relations give

$$
\hat m_P=\sum_{\substack{1\le j\le n_P\\\text{path edges}}}x^P_j Q^P_j,\qquad
\ell_P=w_P(\hat m_A+\hat m_B+\hat m_C)-\hat m_P .
$$

Exterior multiplication and $1-w_A-w_B=w_C$ give the pointwise identity

$$
\ell_A\wedge\ell_B
=w_A\hat m_B\wedge\hat m_C+w_B\hat m_C\wedge\hat m_A+w_C\hat m_A\wedge\hat m_B .
$$

For any barycentric monomial $F=\prod_P w_P^{s_P}\prod_{P,j}(x^P_j)^{k^P_j}$, the theta measure multiplied by $F$ factors over the four simplices. Write $\langle\cdot\rangle_F$ for its mean. The Dirichlet integral gives

$$
\langle w_P\rangle_F
=\frac{13/2-n_P+s_P}{13/2+s_A+s_B+s_C},\qquad
\langle\hat m_P\rangle_F
=\frac{\displaystyle\sum_{\substack{1\le j\le n_P\\\text{path edges}}}
(1+k^P_j)Q^P_j}
{\displaystyle n_P+\sum_{\substack{1\le j\le n_P\\\text{path edges}}}k^P_j}.
$$

In the source normalization $\mathcal K_\theta[F]$ includes $D_\theta/2$, where $D_\theta$ is the determinant of the ten interior momentum columns. Since the two factors $\hat m$ in every term belong to distinct paths,

$$
\mathcal K_\theta[F\,\ell_A\wedge\ell_B]
=\mathcal K_\theta[F]\!
\sum_{\substack{(P,Q,R)=(A,B,C),(B,C,A),(C,A,B)\\\text{cyclic path order}}}
\langle w_P\rangle_F\,\langle\hat m_Q\rangle_F\wedge\langle\hat m_R\rangle_F .
$$

Linearity permits any polynomial numerator to multiply the endpoint insertion, with its own exponents used for each monomial.

For $\theta_{166}$ place the five coordinate covectors $e^1,\ldots,e^5$ in order on path $B$ and $-e^5,\ldots,-e^1$ on $C$. Both endpoint momenta are zero. The source determinant is $D_\theta=-1$. At $F=1$, in this coordinate basis,

$$
\hat m_A=0,\qquad
\langle w_A\rangle_1=\frac{11}{13},\qquad
\langle\hat m_B\rangle_1=\frac{(5,4,3,2,1)}6,\qquad
\langle\hat m_C\rangle_1=-\frac{(1,2,3,4,5)}6 .
$$

The component on the first two coordinate vectors is

$$
\mathcal K_{\theta_{166}}[(\ell_A\wedge\ell_B)_{12}]
=-\frac{\pi}{158400}\,\frac{11}{13}\,\frac{-10+4}{36}
=\frac{\pi}{1123200}.
$$

This value integrates the endpoint insertion with unit remaining numerator. A complete supergravity diagram has eleven additional vertices and their contractions. Direct expansion of the two currents and termwise simplex integration reproduce all ten exterior components on each of the five graphs, for both the unit numerator and a second monomial. [Exact replay](check_endpoint_moment.py)
