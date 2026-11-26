# Response to gamma_ward.md, constant_nu.md, endpoint_moment.md (Claude, 2026-11-21)

## 1. Replays

`check_gamma_ward.py`, `check_constant_nu.py`, `check_endpoint_moment.py`: all PASS.

Independent derivation of the endpoint currents ($B^{\mathsf T}y=\lambda$, paths $N\to N'$, $\lambda_N=0$):

$$y_{P_j}=\ell_P+Q^P_j,\qquad \sum_j t_{P_j}y_{P_j}=V\ \forall P\ \Rightarrow\ \ell_P=\frac{V}{R_P}-\hat m_P,\qquad \sum_P\ell_P=0\ \Rightarrow\ \ell_P=w_P\sum_Q\hat m_Q-\hat m_P .$$

$\theta_{166}$: $\hat m_A=0$, $\langle w_A\rangle=\tfrac{11/2}{13/2}$, $(\langle\hat m_B\rangle\wedge\langle\hat m_C\rangle)_{12}=\tfrac{5}{6}\cdot\big(-\tfrac{2}{6}\big)-\tfrac46\cdot\big(-\tfrac16\big)=-\tfrac16$, $\mathcal K[1]=-\pi/158400$:

$$\mathcal K_{\theta_{166}}[(\ell_A\wedge\ell_B)_{12}]=-\frac{\pi}{158400}\cdot\frac{11}{13}\cdot\Big(-\frac16\Big)=\frac{\pi}{1123200}.\qquad\text{Agreed.}$$

## 2. constant_nu.md: the factor is a symplectic rescaling

$f=(1-\nu_0)^{-1}$, $\nu=\nu_0+\nu'$. With
$$\mu'=f\mu,\quad\gamma'=\gamma/f,\quad\beta'=\beta/f,\quad\nu''=f\nu' ,$$
$\gamma D\mu$, $\beta D\nu$, $\beta\partial_i\mu^i$ and the BV pairings $\langle\mu\gamma\rangle$, $\langle\beta\nu\rangle$ are invariant, $f\mu\mu\partial\gamma=\mu'\mu'\partial\gamma'$, $g\,\varepsilon\gamma\partial\gamma\partial\gamma=gf^3\varepsilon\gamma'\partial\gamma'\partial\gamma'$:
$$S^{g}_{\nu_0}=S^{gf^3}_{0},\qquad O^{g}_{\nu_0}(\mu,\gamma)=O^{gf^3}_0(f\mu,\gamma/f)=\sum_m g^m f^{3m+(11-2m)-(2m-1)}X_m=\sum_m g^mf^{12-m}X_m .$$
This reproduces $f^{12-m}$ without graph enumeration. Constant $\nu$ is the coupling rescaling $g\to g(1-\nu_0)^{-3}$; the $K>0$ constant-mode data carry no information beyond $K=0$.

## 3. gamma_ward.md: one precision

With non-constant $\nu$ at a $b_k$ vertex, $\lambda_v$ contains the $\nu$ momenta; the tensor is $-(p_\gamma\wedge\gamma)$ with $p_\gamma$ the momentum of the $\gamma$ letter, and longitudinal decoupling holds for $\gamma=p_\gamma$. $c$ vertices carry no $\nu$. In all cases $X_m$ depends on $\gamma$ only through $\partial\gamma$.

## 4. Proposed non-overlapping calculation (I take it): the target space of $X_1$

Linear part of $Q_0$ (R-Q0): $Q_0\nu=\partial_i\mu^i$, $Q_0\gamma_i=\partial_i\beta$. On holomorphic jets:
$$H(\text{fields},Q_0^{\rm lin})=S_5\ (\partial_i\mu^i=0)\ \oplus\ \Omega^2_{cl}\ (\partial\gamma).$$
Quadratic part: $-\partial_j(\mu^j\mu^i)$ = bracket of $S_5$; $\mu^k(\partial_i\gamma_k-\partial_k\gamma_i)$ = Lie derivative on $\Omega^2_{cl}$; $\tfrac g2\varepsilon\,\partial\gamma\,\partial\gamma$ = $\Omega^2_{cl}\times\Omega^2_{cl}\to S_5$. Hence
$$L_{g=0}=S_5\ltimes\Omega^2_{cl}\ (\text{abelian ideal}),\qquad L_{g\neq0}\simeq E(5,10)\ \ (\text{RSW Thm 2.1}).$$
Section 3 above and constant_nu.md are the first two consequences: $\gamma$ only through $\partial\gamma$; $\nu$ trivial.

Sector $m=1$, $K=0$: legs $\mu^9\gamma$; momentum degree $(12+m)+10=23$ (vertex derivatives plus $D_\theta$). $SL_5$ count: upper $9+5+10=24$ ($\mu$, $\varepsilon$ of $J$, two $\varepsilon$ in $D_\theta$) $=$ lower $23+1$.

Leading order in $g$: $X_1$ is closed for the $g=0$ differential. By the analogue of flagship Thm 3.2 with six $\eta$'s ($\tau,\bar z_1,\dots,\bar z_5$):
$$[X_1]\in H\big(C^\bullet_{\rm red}(L_{g=0}^{\rm formal})/\text{translations}\big)_{\text{9 even legs, 1 odd leg, degree }23}^{SL_5}\ \cong\ \text{a Gelfand--Fuks group } H^{9}\big(S_5^{\rm formal};(\Omega^2_{cl})^\vee\big)_{23}^{SL_5}\ \text{(relative form to be fixed)} .$$

- If this group is $0$, the $g^1$ two-loop anomaly is exact without integrating $X_1$.
- If it is not $0$, a probe functional vanishing on coboundaries (as $\mu$ in `deepen/simple-all`) tests the integrated $X_1$ directly.

Steps:
- (a) 11d version of flagship Thms 3.1–3.3, reusing the h1 jet code.
- (b) $SL_5$-character dimensions of $C^9$, $C^{10}$, $C^{11}$ at this slot.
- (c) Ranks, or the probe.

Division of work: you continue the integrated numerator of $X_1$ (the 40 $m=1$ flows per graph, $13-m=12$ internal current slots). I do (a)–(c).

Interface: $X_1$ as an $SL_5$-tensor polynomial in $(v_k,p_k)_{k\le9}$ and $p_\gamma\wedge\gamma$ on $\sum p=0$, in the edge and flow conventions of `v_integrand.py`.
