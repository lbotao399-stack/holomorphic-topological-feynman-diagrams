# Fixed-line graph for the ordered OPE coefficient

For the coefficient $f^{(m)}_{(r,t)}(k_1,\ldots,k_{m+1})$ of Fernández–Paquette use their ordered color tensor and their multi-indices $r,t,k_i$. Write $|r|=r^1+r^2$, $r!=r^1!r^2!$. Each coefficient vanishes unless $r+t=(m,m)+\sum_i k_i$. Introduce two covectors $P,Q$ on the transverse two-plane to generate the input modes, with $[P,Q]=P_1Q_2-P_2Q_1$. The polynomial whose $P^rQ^t$ coefficient is
$$
\binom{|r|}{r^1}\binom{|t|}{t^1}(16\pi^2)^m f^{(m)}_{(r,t)}(k_1,\ldots,k_{m+1})
$$
is the following ordered integral:
$$
-[P,Q]^m\int_{1>u_1>\cdots>u_{m+1}>0}
\prod_{i=1}^{m+1}\frac{(|k_i|+1)!}{k_i!}\,
\big(u_iP+(1-u_i)Q\big)^{k_i}\,d u_1\cdots d u_{m+1}.
$$
Here $P^r=P_1^{r^1}P_2^{r^2}$, also for each linear combination in the integrand. Integrate successively from the last variable. This computes all loop orders by polynomial multiplication and one-variable antiderivatives. [Published OPE, Fig. 1, §§3–4 and Appendix C](https://arxiv.org/abs/2412.17168) · [Coefficient evaluator](interval.py)

## Reduction of the fixed-line graph

The scalar $f$ multiplies the ordered chain color tensor. In its bulk diagram the two external fields sit at the endpoints of a spine of $m+1$ cubic vertices. Each vertex has one spoke ending on the defect. At fixed positive Schwinger times the longitudinal coordinate sees a tree: the $m+1$ defect attachments have independent longitudinal positions. After removing one common translation there are $2m+1$ complex longitudinal differences, one per edge. The incidence matrix of this tree has determinant of magnitude one. The external antiholomorphic tangent form supplies the differential of the common translation. Every edge must then supply its remaining longitudinal differential to saturate relative tangent degree. In edge-difference coordinates each longitudinal Gaussian has unit mass after its heat normalization. A holomorphic longitudinal coordinate from the test ghost averages to the common translation. This gives the tangential reduction for the universal local collision in Schwinger space; the associativity check below fixes the coefficient in the published normal-product prescription.

Both transverse coordinates vanish at every defect attachment. They therefore see the $m$-loop fan of the local paper. Inject the two external transverse momenta at the ends of its spine. Let $\alpha_i,\beta_i$ be the currents flowing through spoke $i$ for respective *unit* injections at the first and last vertex. They are positive, with $\sum_i\alpha_i=\sum_i\beta_i=1$. The path with the final spoke is a spanning tree. Its cotree currents are $\alpha_iP+\beta_iQ$ for $1\le i\le m$. Their Jacobian is $[P,Q]^m$ times ordinary volume in the first $m$ coordinates of each probability simplex. [Paper, fan and current conventions](../../../sources/LLC_main.tex)

In these coordinates the entire positive-resistance region is
$$
\frac{\alpha_1}{\beta_1}>\frac{\alpha_2}{\beta_2}>\cdots>\frac{\alpha_{m+1}}{\beta_{m+1}} .
$$
Indeed, let $A_i,B_i$ be the partial sums of $\alpha,\beta$ through $i$. If $s_i$ is the resistance of spoke $i$ and $d_i$ that of path edge $(i,i+1)$, the two voltage drops are $s_i\alpha_i-s_{i+1}\alpha_{i+1}=d_i(1-A_i)$ and $s_i\beta_i-s_{i+1}\beta_{i+1}=-d_i B_i$. Thus
$$
\frac{s_{i+1}}{s_i}=\frac{\alpha_iB_i+\beta_i(1-A_i)}{\alpha_{i+1}B_i+\beta_{i+1}(1-A_i)},\qquad
\frac{d_i}{s_i}=\frac{\alpha_i\beta_{i+1}-\alpha_{i+1}\beta_i}{\alpha_{i+1}B_i+\beta_{i+1}(1-A_i)} .
$$
These solve the inverse map uniquely up to $s_1$. The region has volume $1/[(m!)^2(m+1)!]$ by simultaneously permuting the two simplex labels.

The scalar fan differentiates the two external transverse monomials by $(\partial_{X^1}\partial_{Y^2}-\partial_{X^2}\partial_{Y^1})^m$ and then substitutes $X=\sum_i\alpha_i z_i$, $Y=\sum_i\beta_i z_i$, where $z_i$ is the transverse displacement of the $i$th defect attachment. Extract $\prod_i z_i^{k_i}$ with the factorials of the source's defect jets.

For nonnegative exponents $a_i,b_i$, the integral of $\prod_i\alpha_i^{a_i}\beta_i^{b_i}$ on the ordered region is
$$
\frac{\prod_i a_i!b_i!}{(m+\sum_i a_i)!(m+\sum_i b_i)!}
\int_{u_1>\cdots>u_{m+1}}\prod_i
\frac{(a_i+b_i+1)!}{a_i!b_i!}u_i^{a_i}(1-u_i)^{b_i}\,du_i
$$
with $0<u_i<1$. To verify it, introduce two independent positive scales for the simplexes, multiply by the two normalized gamma densities of shapes $m+1+\sum a_i$ and $m+1+\sum b_i$, then change each pair of scaled coordinates $(x_i,y_i)$ to $(x_i+y_i,u_i=x_i/(x_i+y_i))$. The common scales preserve ratio order. The $x_i+y_i$ integrals give the displayed factorials. Expanding the two external monomials now gives the ordered integral above: in each spoke factor $a_i$ counts powers assigned to $P$ and $b_i$ those assigned to $Q$. The $m$ differentiations give the determinant power. Its common sign and physical loop factor are fixed by the one-loop seed and the following exact recurrence check.

## Identification with the published recurrence

Write $B(u)=uP+(1-u)Q$ and $A_L=(|L|+1)!/L!$ for a mode $L$. Let $G(k_1,\ldots,k_{m+1})$ denote the displayed generating polynomial. For the last two output modes $K,L$, let $e_1=(1,0)$, $e_2=(0,1)$, and a negative mode give zero. With $W=(k_1,\ldots,k_{m-1})$, integration in the last variable uses
$$
A_L B^L-\sum_{j=1}^2 Q_j A_{L-e_j}B^{L-e_j}
=\frac{A_L}{|L|+1}\frac{d}{du}\left(uB(u)^L\right).
$$
At its upper boundary $u=u_m$, $Q_1B_2-Q_2B_1=-u[P,Q]$. Therefore
$$
G(W,K,L)=\sum_{j=1}^2 Q_jG(W,K,L-e_j)
-Q_1\frac{A_KA_L}{(|L|+1)A_{K+L+e_2}}G(W,K+L+e_2)
+Q_2\frac{A_KA_L}{(|L|+1)A_{K+L+e_1}}G(W,K+L+e_1).
$$
For the auxiliary one-letter term on the right use $G(H)=-A_H\int_0^1B^Hdu$. Its $P^rQ^t$ coefficient for $r+t=H$ is $-\binom{|r|}{r^1}\binom{|t|}{t^1}$ by the beta integral.

To recover exactly this recursion from the source, take its unit step $s=e_j$ in (74), multiply by $t^j/|t|$, and sum over $j$. In the generating polynomial the identity $\binom{|t|}{t^1}t^j/|t|=\binom{|t|-1}{(t-e_j)^1}$ shifts the last input. The term shifting the last output gives the first sum. In the first product term the intermediate mode is $H=K+L+e_{3-j}$. The special one-loop formula gives
$$
(16\pi^2) f^{(1)}_{(H,e_j)}(K,L)=(-1)^{j+1}\frac{A_K A_L}{(|L|+1)A_H}.
$$
This gives both fused terms. In the remaining product the intermediate mode $h=t-(1,1)-L$ is independent of $j$. Whenever $h$ is nonnegative, both steps exist and
$$
\binom{|t|-1}{(t-e_j)^1}(16\pi^2) f^{(1)}_{(t-e_j,e_j)}(h,L)
=(-1)^{j+1}\frac{|h|+1}{|t|}\binom{|h|}{h^1}\binom{|L|}{L^1}.
$$
The two terms cancel. At $m=1$, the four-term source recurrence gives the same shifted term, the same cancellation with $h=K-r$, and the one-letter beta integral in place of the lower-loop product. The sign in $G(H)$ accounts for the positive $f^{(1)}_{(r+t-e_j,e_j)}(K,L)$ on that right-hand side. The formula has no terms of input degree below $m$. Induction on $m$ and then on $|t|$ determines every conserved coefficient. In particular its seed is $f^{(1)}_{(0,1),(1,0)}(0,0)=1/(32\pi^2)$. [Source recurrence replay](coefficients.py)

For arbitrary nonnegative $a_1,\ldots,a_{m+1}$ take $r=(0,m+\sum_i a_i)$, $t=(m,0)$, $k_i=(0,a_i)$. The two evaluations simplify to
$$
(16\pi^2)^m f^{(m)}_{(r,t)}(k_1,\ldots,k_{m+1})
=(-1)^{m+1}\prod_{j=1}^{m}\frac{1+a_j}{\sum_{i=j}^{m+1}(1+a_i)} .
$$
The ordered integral has independent beta densities $(1+a_i)u_i^{a_i}$. Integrating from the last endpoint gives the product. Independently, the source recurrence merges the final two exponents into $a_m+a_{m+1}+1$ with factor $-(1+a_m)/(2+a_m+a_{m+1})$ at each step beyond the one-loop seed. At zero output modes the value is $(-1)^{m+1}/[(16\pi^2)^m(m+1)!]$.

At two loops, $r=(1,3)$, $t=(2,1)$, and ordered outputs $((0,1),(1,1),(0,0))$ give $1/[60(16\pi^2)^2]$. Reversing the outputs gives $-1/[40(16\pi^2)^2]$. The replay compares the published recurrence in both unit-step orders, the ordered integral and the independent Schwinger fan coordinates. [Fan transfer](fan.py) · [Exact comparisons](check.py)

The longitudinal check grounds each of four different vertices at every loop order through nine, and verifies $\det(B^{\mathsf T}t^{-1}B)\prod_e t_e=1$ for the full tree incidence matrix $B$ with the grounded column removed. It also recovers 45 path resistances from the transverse current ratios and checks 840 individual ordered one-loop coefficients against the unreduced fixed-line residue. [Direct heat check](check_heat.py)
