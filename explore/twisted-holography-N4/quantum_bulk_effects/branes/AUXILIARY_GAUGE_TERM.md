# A nonzero one-loop cross-contraction with auxiliary gauge transformations

**Scope.** This calculation evaluates the paper's displayed superfield fan on an insertion containing the undifferentiated ghost. It gives a concrete obstruction to tensoring that unreduced fan with the classical auxiliary gauge action. The physical comparison in Appendix B is relative to constant gauge transformations and excludes this ghost. Therefore the computation does not, by itself, identify a unique full absolute quantum differential or prove that its missing auxiliary corrections exist. No corrected brane differential is asserted.

Source: Li–Liu–Chang, [arXiv:2611.24160v2](https://arxiv.org/pdf/2611.24160v2), equations (8.4), (8.7), (A.3), (A.4), and Proposition B.1. All derivatives below retain the paper's spinor convention. The undifferentiated ghost is the letter denoted $c$ in (A.2).

## 1. The ghost–field-strength contraction

For an insertion without input derivatives, equation (8.4) at one loop is

$$
Q_1=-\frac{1}{128\pi^2}\int d^3\theta\,
C_{ARD}C_{BRE}
(\partial_{+\dot a}C^D)(\partial_+{}^{\dot a}C^E)
\frac{\delta}{\delta C^A(\theta)}
\frac{\delta}{\delta C^B(\theta)}.
$$

The two relevant terms of (A.4) are

$$
\frac{\delta}{\delta C^A(\theta)}
=\frac{\theta_1\theta_2\theta_3}{ig}
\frac{\partial_L}{\partial c^A}
-ig\frac{\partial}{\partial b_A}
 +\text{the scalar and fermion derivatives}.
$$

Each assignment of the two derivatives to $c^A b^B$ contributes

$$
\frac{\theta_1\theta_2\theta_3}{ig}(-ig)
=-\theta_1\theta_2\theta_3.
$$

The two output superfields consequently contribute only their ghost terms, whose coefficient is $(ig)^2=-g^2$. Their spinor contraction is symmetric in $D,E$. Thus

$$
\boxed{
Q_1(c^A b^B)
=-\frac{g^2}{64\pi^2}
C_{ARD}C_{BRE}
(\partial_{+\dot a}c^D)(\partial_+{}^{\dot a}c^E).
}
$$

This is the negative of equation (8.7) for the ordered scalar–fermion insertion with the same free color indices. Its sign follows directly from $(-\theta_2\theta_3)(-\theta_1)=+\theta_1\theta_2\theta_3$ in that insertion, compared with the minus sign above.

## 2. The auxiliary anticommutator

Let $\psi$ be an odd auxiliary column transforming in the fundamental representation. Write $R_0$ for the ordinary auxiliary gauge action, acting only on auxiliary variables:

$$
R_0\psi=igc^At_A\psi,
\qquad
R_0\bar\psi=ig\bar\psi c^At_A.
$$

These signs use odd $\bar\psi,\psi$ and a left-acting odd differential. Together with $Q_0(igc)=(igc)^2$, they give $(Q_0+R_0)^2=0$ classically. On $b^B\psi$, the physical one-loop fan has insufficient inputs and vanishes. Consequently

$$
\begin{aligned}
\{Q_1,R_0\}(b^B\psi)
&=Q_1(igc^A b^B t_A\psi)\\
&=\boxed{-\frac{ig^3}{64\pi^2}
C_{ARD}C_{BRE}
(\partial_{+\dot a}c^D)(\partial_+{}^{\dot a}c^E)t_A\psi}.
\end{aligned}
$$

It is nonzero already for $SU(2)$. Set $t_A=\sigma_A/2$, $C_{ABC}=\epsilon_{ABC}$, $B=3$, and choose two independent odd ghost amplitudes $\epsilon_+,\epsilon_-$ in color $1$, with derivative spinors satisfying $p_{+\dot a}q_+{}^{\dot a}=1$. Then

$$
\sum_R C_{AR1}C_{3R1}=\delta_{A3},
\qquad
[\epsilon_+\epsilon_-]
(\partial_{+\dot a}c^1)(\partial_+{}^{\dot a}c^1)=2,
$$

and hence

$$
\boxed{
[\epsilon_+\epsilon_-]\{Q_1,R_0\}(b^3\psi)
=-\frac{ig^3}{32\pi^2}t_3\psi\ne0.
}
$$

If a supplied absolute physical differential begins with precisely this $Q_1$, a corrected auxiliary extension $R_0+\hbar R_1+\cdots$ must satisfy

$$
\{Q_0+R_0,R_1\}=-\{Q_1,R_0\}.
$$

The calculation determines a nonzero term on the right. It neither solves this equation nor excludes additional absolute terms that are invisible after relative restriction. In particular, it justifies avoiding a naive joint relative auxiliary complex.

## 3. Independent one-loop determinant normalization

Take the standard $SU(N)$ normalization $\operatorname{Tr}(t_At_B)=\tfrac12\delta_{AB}$. Contracting equation (8.7) with the trace metric and using $C_{ARD}C_{ARE}=N\delta_{DE}$ gives

$$
Q_1\operatorname{Tr}(\gamma^1\beta_1)
=\frac{g^2N}{64\pi^2}
\operatorname{Tr}\big((\partial_{+\dot a}c)(\partial_+{}^{\dot a}c)\big).
$$

At $N=2$, tracelessness gives the exact identity

$$
\operatorname{Tr}\big(\operatorname{adj}(m\mathbf1+\gamma^1)\beta_1\big)
=-\operatorname{Tr}(\gamma^1\beta_1).
$$

For the ordered root ghosts $E_{12},E_{21}$ and spinor contraction $1$, the coefficient in the trace on the right is $2$. Therefore

$$
\boxed{
[\epsilon_+\epsilon_-]Q_1
\operatorname{Tr}\big(\operatorname{adj}(m\mathbf1+\gamma^1)\beta_1\big)
=-\frac{g^2}{16\pi^2}.
}
$$

The determinant calculation gives the same value: its two-spoke color polynomial is $-4$ at $N=2$, so $g^2(\tfrac12)(-4)/(32\pi^2)=-g^2/(16\pi^2)$. For a literal odd source $\eta_1$, the coefficient of $\eta_1$ in $\mathbf Q\det(m\mathbf1+\gamma^1+\eta_1\beta_1)$ is the negative of this operator variation, because $\mathbf Q$ passes the odd source.

The exact verifier `auxiliary_gauge_check.py` checks the three Grassmann derivative signs, the ordered ghost factor, the $SU(2)$ color contraction, and this determinant normalization. Its result is recorded in `auxiliary_gauge_check.json`.
