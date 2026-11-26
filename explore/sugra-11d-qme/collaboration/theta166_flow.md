# One integrated orientation of the 11d two-loop graph $\theta_{166}$

Keep the project’s edge order $(A_1,B_1,\ldots,B_6,C_1,\ldots,C_6)$ from $N$ to $N'$. The cubic $\gamma\partial\gamma\partial\gamma$ vertex sits at $N$. The external $\gamma$ sits at $B_1$ with covector polarization $e^2$. Each of the remaining nine interior vertices carries an external $\mu$ with vector polarization $e_1$. Here $e_i$ and $e^i$ are the standard vector and dual bases. Take the momentum sequences

$$
\lambda_N=\lambda_{N'}=0,\qquad
\lambda^B=(e^1,e^2,e^3,e^4,e^5),\qquad
\lambda^C=-(e^1,e^2,e^4,e^3,e^5),\qquad D_{166}=1.
$$

Orient each contraction from $\mu$ to $\gamma$. Path $A$ runs to $N$, path $C$ runs to $N$, and the two segments of $B$ run from $B_1$ to opposite endpoints. In increasing edge order, the tail-vertex list in the project’s numbering is $(11,1,1,2,3,4,5,6,7,8,9,10,11)$, with eleven inversions. The Wick sign is $-1$. [Contractions and signed flow numerator](../VERIFY.md)

Write $(\ell_B)_i,(\ell_C)_i$ for the first-edge current components. At a pass-through vertex with external $u=e_1$ and signed incident momentum $J$, the matrix in path order is $uJ^{\mathsf T}-(J^{\mathsf T}u)\mathbf1$, or its transpose for reverse flow. Multiplication of two forward matrices leaves the first times $-J_{\rm second}^{\mathsf T}u$. Multiplication of two reverse matrices leaves the second times $-J_{\rm first}^{\mathsf T}u$. The four $B$ matrices contribute the scalar $-(1+(\ell_B)_1)^3$ and the five $C$ matrices contribute $(1-(\ell_C)_1)^4$. Contracting the remaining matrices with the $B_1$ source and endpoint tensors gives the complete numerator with $g$ removed,

$$
\begin{aligned}
&(1+(\ell_B)_1)^3(1-(\ell_C)_1)^4(2(\ell_B)_1+1)((\ell_B)_1+2(\ell_C)_1-1)\\[-2pt]
&\hspace{23mm}\times\det\begin{pmatrix}
1&(\ell_B)_3&(\ell_C)_3\\1&(\ell_B)_4&(\ell_C)_4\\1&(\ell_B)_5&(\ell_C)_5
\end{pmatrix}.
\end{aligned}
$$

For either six-edge path, let $t_j=x_{j+1}+\cdots+x_6$ be its five tail fractions. The uniform path simplex gives density $5!$ on $1>t_1>\cdots>t_5>0$. Sequential integration gives

$$
\frac{\int_{\Delta^5}t_1^{k_1}\cdots t_5^{k_5}\,dx}{\int_{\Delta^5}dx}
=\frac{5!}{(1+k_5)(2+k_4+k_5)\cdots(5+k_1+\cdots+k_5)}.
$$

The two paths are independent in the normalized theta measure. The Dirichlet parameters of $(w_A,w_B,w_C)$ are $(11/2,1/2,1/2)$. Expand the first-edge currents from $\ell_P=w_P(\hat m_B+\hat m_C)-\hat m_P$ and integrate the resulting monomials. The signed mean above and the value in the project’s kernel normalization are

$$
-\frac{3419207}{5558054638500},\qquad
\frac{\pi}{158400}\left(-\frac{3419207}{5558054638500}\right)
=-\frac{310837\pi}{80035986794400000}.
$$

Putting the $C$ momenta in the same coordinate order as the $B$ momenta, with the opposite signs, gives an exact zero integral. Both evaluations concern one orientation at fixed external labels. Odd-field antisymmetrization, graph multiplicities and the common physical normalization enter the assembled obstruction separately. [Theta integral with insertions](../../../../研究项目_非Git/Twisted_11d_sugra_loops_2026-11-18/main.tex) · [Exact evaluator and pointwise controls](check_theta166_flow.py)

## Independent integration conditioned on the path maxima

Condition on the first tail fractions $s=t^B_1$ and $t=t^C_1$. Their densities are $5s^4$ and $5t^4$. The other four fractions on each path are the decreasing order statistics of four independent uniforms on $[0,s]$ or $[0,t]$. Thus their five conditional means in path order are $(5,4,3,2,1)s/5$ and $(5,4,3,2,1)t/5$. Write the five tails of each path in **coordinate** order, so that the third and fourth $C$ entries are exchanged in the stated fixture. Since $\hat m_B$ is the positive $B$ tail vector and $\hat m_C$ the negative $C$ tail vector, their positive tail vectors map to $(\ell_B,\ell_C)$ by the coefficient matrix

$$
\begin{pmatrix}w_B-1&-w_B\\w_C&1-w_C\end{pmatrix},
\qquad \det=-w_A.
$$

The determinant with rows $3,4,5$ in the numerator is linear in each path's tail vector and invariant under additive constants in either column. Conditional independence gives its exact conditional mean

$$
-w_A\left(\frac{2s}{5}\frac{2t}{5}-\frac{s}{5}\frac{t}{5}\right)
=-\frac{3w_Ast}{25}.
$$

Identical coordinate orders make the same determinant zero after conditional averaging. All remaining factors of the numerator depend only on the first components $(\ell_B)_1=(w_B-1)s-w_Bt$ and $(\ell_C)_1=w_Cs+(1-w_C)t$.

To reproduce the scalar integrations, write $h=w_B+w_C$ and $v=w_B/h$. The Dirichlet density factors into independent beta densities with parameters $(1,11/2)$ and $(1/2,1/2)$ for $h$ and $v$. In a monomial $s^it^jh^kv^\ell$ from the factors multiplying the determinant, the conditional determinant supplies $-3(1-h)st/25$. Its normalized integral is therefore

$$
-\frac{3}{(6+i)(6+j)}
\frac{(1)_k(11/2)}{(13/2)_{k+1}}
\frac{(1/2)_\ell}{(1)_\ell},
\qquad (x)_n=x(x+1)\cdots(x+n-1),\quad (x)_0=1.
$$

The independent verifier contracts the literal Wick words over the polynomial ring in all ten first-edge current components. For the four flows compatible with the specified external letters it checks the complete signed polynomials, both coordinate orders, and 24 additional positive rational resistance configurations. It reproduces the single-flow value above and gives $-3877063\pi/1408633367581440000$ after adding all four flows at these fixed labels. [Independent symbolic and integral checks](verify_theta166_flow.py)
