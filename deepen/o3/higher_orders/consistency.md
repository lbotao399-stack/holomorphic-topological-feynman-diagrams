# Consistency one order above the first local obstruction

In the manuscript's heat-kernel gauge on $\mathbb R\times\mathbb C^2$, let $S=I_0+\hbar J_1+\cdots$ be any formal interaction built from translation-invariant full-wedge vertices with holomorphic derivatives. The classical vertex $I_0$ is cubic. Write $O(S)$ for its local QME curvature, the sum of the free differential $QS$ and all connected total-collapse contact terms, with the single-edge term $\tfrac12\{S,S\}$. If its coefficients below order $w$ vanish, then

$$
(Q+\{I_0,-\})O_{w+1}+\{J_1,O_w\}=0.
$$

Here $O_{w+1}$ is evaluated on the same trial interaction as $O_w$. Its existence as a local cochain does not require solving the order-$w$ equation.

**Finite scale.** Denote the sum of connected effective graphs at scale $L>0$ by $S[L]$. Differentiating the graph sum with respect to a local vertex inserts one marked vertex, with the differentiated symmetry factors. The effective curvature is

$$
O[L]=QS[L]+\hbar\Delta_L S[L]+\tfrac12\{S[L],S[L]\}_L,
\qquad
(Q+\hbar\Delta_L+\{S[L],-\}_L)O[L]=0 .
$$

The second identity follows by applying $(Q+\hbar\Delta_L)^2=0$ to $\exp(S[L]/\hbar)$. The first variation of the effective interaction in the direction $O(S)$ equals $O[L]$. Indeed the compactified Schwinger Stokes formula leaves the free differential and the ultraviolet faces after the $L$-faces combine into the displayed BV terms. A connected collapsing subgraph produces its local contact at the contracted vertex; the remaining macroscopic edges give the effective graph with this single marked vertex. Disconnected collapsing subgraphs have zero contact. Summation over marked subgraphs gives precisely the first-variation multiplicities, by orbit–stabilizer. This argument applies before imposing any order of the QME. It also applies after adjoining an odd square-zero coefficient to an interaction; ordered Grassmann coefficients carry the signs of the marked vertices. [Wang–Williams, equations (69)–(72), Propositions 3.21–3.23](https://arxiv.org/html/2407.08667v3#S3.SS3) · [Costello, equation (12.0.1) and Proposition 12.0.3](https://arxiv.org/pdf/0706.1533)

View the effective interaction as a formal map of graded spaces of vertices, with derivative $D(S\mapsto S[L])$. The preceding equality says that the local curvature vector field and the effective curvature vector field are related by this derivative. Differentiating again in the odd curvature direction gives

$$
0=D^2(S\mapsto S[L])(O(S),O(S))+
D(S\mapsto S[L])\big(DO(S)[O(S)]\big).
$$

The Hessian is graded symmetric and $O(S)$ is odd, so its term vanishes. The remaining derivative is injective on the completed vertex space. At the first nonzero $\hbar$ order and then the least nonzero field arity of an insertion, the uncontracted one-vertex graph returns the insertion itself. A tree with $r>0$ additional classical cubic vertices adds $r$ external legs, and any positive-loop graph or positive-order additional vertex raises $\hbar$ order. Hence $DO(S)[O(S)]=0$.

**Coefficient at $w+1$.** A marked graph has order equal to the insertion order plus its loop number plus the orders of its ordinary vertices. The marked tree contact is the BV bracket: for a tree with $v$ vertices and $v-1$ edges, the contact condition $3v=2(v-1)+4$ forces $v=2$. At order $w+1$ the free differential contributes $QO_{w+1}$; tree contacts contribute $\{I_0,O_{w+1}\}+\{J_1,O_w\}$. All remaining contacts have one loop and one $O_w$ insertion. The marked vertex is a local full-wedge vertex with only holomorphic derivatives, by the same total-collapse theorem. The odd-loop reflection therefore kills each of these contacts. A self-contraction has zero pullback of the heat form or propagator to the diagonal in this gauge, including after holomorphic differentiations. This proves the formula. [Marked contact hypotheses](../../h1/verify/regulated_contacts.md)

**Central class.** After a local counterterm makes the actual $O_w$ purely central, its bracket with $J_1$ is a classical boundary. In the manuscript's finite-jet complex the cohomology class of $\{j,O_w\}$ therefore vanishes in Chevalley–Eilenberg degree five. For this application one needs a detection argument for its class to infer vanishing of the central obstruction. The consistency identity alone does not evaluate the underlying graph coefficients.
