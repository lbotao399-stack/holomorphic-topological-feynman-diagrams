# The Chevalley–Eilenberg boundary on the ten paired Fourier modes

**Proposition.** In the ordinary algebraic Chevalley–Eilenberg chain complex of $S_5\ltimes\Omega^2_{\mathrm{cl}}$, take the even divergence-free holomorphic vector fields and the abelian odd space of closed holomorphic two-forms. At each of the ten momenta $+e^1,-e^1,\ldots,+e^5,-e^5$ put one Fourier mode, with nine vector fields and one two-form. The chain boundary is injective on the span of these chains. This includes combinations in which the two-form occupies different momenta. All ten-input chains considered here have zero total momentum.

Here $e_i$ and $e^i$ are the coordinate tangent and dual covector bases on $\mathbb C^5$, with axis labels $i=1,\ldots,5$. At covector momentum $p$, a vector mode is $e^{p(z)}u$ with $p(u)=0$. A closed two-form mode at $p\ne0$ has polarization $p\wedge\gamma$, with $\gamma$ defined modulo multiples of $p$. Each of the ten polarization spaces has dimension four. Finite sums of these exponential modes form a momentum-graded subalgebra upon allowing all momenta in the lattice they generate. The boundary sums brackets of pairs of slots with their Chevalley–Eilenberg signs.

**The particular input.** At momenta $e^2,e^3$ the transverse fixture uses the vector polarizations $e_3,e_4$. Their bracket is

$$
[e^{z^2}e_3,e^{z^3}e_4]=e^{z^2+z^3}e_4 .
$$

Only this pair produces the momentum $e^2+e^3$. The other seven vector modes and the two-form stay at distinct original momenta. This component of the boundary is nonzero.

**Arbitrary polarizations.** First fix the two-form slot. On each axis $i$ take one of the vector slots at momentum $s_i e^i$, with $s_i\in\{+1,-1\}$. Such a choice exists even on the axis carrying the two-form, because its opposite momentum carries a vector. Fix basis polarizations at the remaining four vector slots and at the two-form slot; they are spectators. Contributions fusing chosen slots $i\ne j$ have the unique new momentum $s_i e^i+s_j e^j$. Among these ten input momenta it has exactly one decomposition as a sum of a pair. The spectator polarizations consequently remain independent in this component of the boundary.

For two vector polarizations $u,v$, with $u^j=e^j(u)$ and $v^i=e^i(v)$, the merged vector polarization is

$$
s_j u^j v-s_i v^i u .
$$

The Chevalley–Eilenberg reordering sign is constant for this pair and these spectator slots. For fixed polarizations on the three other chosen axes, let $C(x,y)$ denote the coefficient with polarizations $e_x,e_y$ at $i,j$. Transversality gives $x\ne i$ and $y\ne j$. The coordinate-$i$ component of the bracket can only come from $(x,y)=(j,i)$, and the coordinate-$k$ component for $k\notin\{i,j\}$ has exactly two sources. A cycle therefore satisfies

$$
s_j C(j,i)=0,\qquad
s_j C(j,k)-s_i C(k,i)=0 .
$$

For each basis coefficient of the five selected vectors draw $i\to j$ when the vector at $s_i e^i$ has polarization $e_j$. Every axis has one outgoing arrow and no arrow is a loop. Following the arrows gives a directed cycle. If it has length two, its coefficient vanishes by the first equation. Otherwise choose three consecutive distinct axes $i\to j\to k$ on that cycle. The second equation expresses its coefficient as $s_i/s_j$ times the coefficient with arrows $i\to k$ and $j\to i$. These new arrows shorten the cycle by one. Induction on the shortest cycle length kills every coefficient, for every choice of spectator polarizations.

**Varying the two-form slot.** In all target components just used, the fused two-axis momentum carries a vector and the two-form still has its original one-axis momentum. An input with the two-form at a different original momentum could reach that component only by fusing a vector with the two-form. Such a bracket places the two-form at the fused momentum instead. Thus the proofs for the ten choices of two-form slot do not interfere.

**Integral certificate.** On five chosen axes there are $4^5=1024$ basis coefficients. Ordering them by shortest cycle length gives 580, 300, 120 and 24 columns at lengths two through five. For each column the proof selects a boundary row with diagonal entry $\pm1$; any second entry has smaller cycle length. For all 32 sign choices a direct integer replay checks distinct rows and a triangular minor with determinant $\pm1$. For a fixed two-form slot the other $4^5$ spectator choices give independent copies. A homology pairing for the two-loop symbol must use chains outside this prescribed Fourier support. [Exact verifier](verify_fixed_support_cycles.py)
