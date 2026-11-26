# Constant external $\nu$ in the two-loop bulk graphs

In the twelve-interaction bracket $\{\mathcal I^{\otimes12}\}_0$, fix one of the five theta graphs and let $m$ count its cubic $\gamma\partial\gamma\partial\gamma$ vertices. The field contractions pair $\mu$ with $\gamma$ and $\nu$ with $\beta$. A vertex $\beta\partial_i\mu^i$ has only two letters. At a trivalent vertex three contractions are required. At an interior vertex both letters would be contracted, giving zero output momentum in a column of $D_\theta$. Thus $\beta$ is absent and every $\nu$ is external. The other $12-m$ vertices are of type $b_k=\nu^k\mu^i\mu^j\partial_i\gamma_j$ with $k\ge0$. The thirteen edges consume thirteen $\mu$ and thirteen $\gamma$, leaving $11-2m$ external $\mu$ and $2m-1$ external $\gamma$. Hence $1\le m\le5$. [Interaction, contractions and theta kernel](../../../../研究项目_非Git/Twisted_11d_sugra_loops_2026-11-18/main.tex) · [Graph checks a2–a4](../VERIFY.md)

Restrict the external even field $\nu$ to its constant mode. Every holomorphic derivative hitting $\nu$ vanishes. Assigning $k$ copies to one $b$ vertex multiplies its integrand by $\nu^k$ and contributes zero additional momentum. Summing its allowed vertices gives $f(\nu)=(1-\nu)^{-1}$ as a formal series. At fixed $m$, each bare integrand is therefore multiplied by

$$
f(\nu)^{12-m}=(1-\nu)^{-(12-m)},\qquad
[\nu^K] f(\nu)^{12-m}=\binom{K+11-m}{11-m},\qquad K\ge0.
$$

Here $[\nu^K]$ extracts the ordinary monomial coefficient. The factor counts the distributions of $K$ equal insertions among the $12-m$ vertex positions. Differentiating $K$ times at zero, as for $K$ separately labelled constant insertions, gives $K!\binom{K+11-m}{11-m}$. These factors apply before Schwinger integration, including when the bare integrand vanishes. Terms containing derivatives of $\nu$ require external momenta away from this specialization.

With $k\le r$ at each $b_k$, the factor becomes $(1+\nu+\cdots+\nu^r)^{12-m}$ and the largest allowed $K$ is $r(12-m)$. For $m=1$, the unrestricted first three coefficients are $1,11,66$. Keeping $k\le1$ gives $1,11,55$. The difference of $11$ counts the single $b_2$ insertion. The action permits arbitrarily large $K$. [Exact checker](check_constant_nu.py) · [Transfer-matrix enumeration](../verify/v_structure.log)

## Arbitrarily many external scalar waves

Fix a graph and a field orientation of sector $g^m$. There are $12-m$ vertices carrying $f(\nu)=(1-\nu)^{-1}$. Label the external scalar waves before taking equal momenta. Placing $n_v$ labeled waves at vertex $v$ gives the vertex factor $n_v!$. All undifferentiated scalar factors have unit polarization here. Allow any labeled set of scalar momenta; a single further zero-momentum scalar multiplies its polarized two-loop coefficient by $12-m+K$, where $K$ is the number already present. Indeed, insertion at $v$ multiplies the local weight by $n_v+1$, and the sum of these factors is $12-m+K$. The currents and the determinant are unchanged in this insertion. This identity holds term by term in the graph and orientation after summing the vertex positions of the soft scalar.

For repeated **nonzero** momenta, let $K_i$ be the number of scalar labels carrying momentum $k_i$, and $K$ their total. To conserve momentum for all multiplicities, fix one of the external vector waves and subtract $\sum_{i\in\mathrm{modes}}K_i k_i$ from its momentum, keeping its polarization fixed. A divergence-free family can use a vector polarization annihilated by its initial momentum and by every $k_i$. Denote the polarized graph coefficient by $A_{K_1,\ldots,K_s}$ for $s$ distinct specified momenta. Then
$$
\frac{A_{K_1,\ldots,K_s}}{(12-m)(13-m)\cdots(11-m+K)}
\quad\text{is a polynomial in }(K_1,\ldots,K_s)
\quad\text{of total degree at most }23-m .
$$
The empty product is one. The displayed denominator removes only the known number of labeled attachments.

To see the degree bound, fix the counts $n_{iv}$ of each mode at each vertex. The momentum in the Kirchhoff equation at $v$ is shifted by $\sum_{i\in\mathrm{modes}}n_{iv}k_i$, in addition to the compensating vector shift. Each edge current is linear in these loads. There are ten momentum factors in $D_\theta$ and $12+m$ in the differentiated vertices. Each external one-form supplies one factor of its own fixed momentum, at each of the $2m-1$ one-form sites, by the [one-form Ward calculation](gamma_ward.md). Cubic one-form vertices never carry $\nu$; at a vertex with two outgoing vector edges the derivative acts only on its one-form. Hence at most $10+12+m-(2m-1)=23-m$ factors depend on the scalar multiplicities. The Schwinger measure is independent of all momenta.

Normalize the weight of an individual assignment of the scalar labels by dividing $\prod_{v\in f\text{-vertices}} n_v!$ by the product in the displayed equation, with $n_v=\sum_{i\in\mathrm{modes}}n_{iv}$. This is a probability distribution since the next label has the $n_v+1$ insertion weights computed above. For a nonnegative integer $l$, write $u^{\underline l}=u(u-1)\cdots(u-l+1)$. For integers $r_{iv}\ge0$, write $r_i$, $r_v$ and $r$ for their row sums, column sums and total. Under the assignment distribution the exact factorial moment is
$$
\mathbb E\!\left[\prod_{i,v}n_{iv}^{\underline{r_{iv}}}\right]
=\frac{\displaystyle\prod_{i\in\mathrm{modes}}K_i^{\underline{r_i}}\,
              \prod_{v\in f\text{-vertices}}r_v!}
             {(12-m)(13-m)\cdots(11-m+r)} .
$$
There are $\prod_i K_i^{\underline{r_i}}$ ways to distinguish and order the selected labels. Any prescribed placement of those labels has the indicated probability, independently of subsequent placements. Expanding a polynomial in falling powers preserves its degree. The expectation of the graph polynomial has the asserted degree, also after summing labels of the vector and one-form inputs, orientations, and topologies. In $s$ specified momentum modes it suffices to know the samples with total multiplicity at most $23-m$. The replay implements exact multivariate Newton interpolation; a component-specific lower degree can be used when proved.

For the $\theta_{166}$ endpoint component of the [unit-numerator calculation](endpoint_moment.md), put the cubic one-form vertex at the first trivalent endpoint, insert $K$ scalar waves of momentum $e^1+2e^2$ and compensate the vector at interior vertex 3. The two path momentum determinants are each affine in the loads by the rank-one determinant identity. The two-current wedge is affine since the shared shift covector wedges with itself to zero. The degree bound here is three, and the value divided by $\pi (11\cdot12\cdots(10+K))$ is
$$
\frac1{1123200}
-\frac{\binom K1}{1235520}
+\frac{\binom K2}{6739200}
+\frac{\binom K3}{2141568}.
$$
This is the endpoint insertion with its remaining numerator set to one. Direct expansion checks all $12{,}376$ weak compositions at $0\le K\le6$, with the circuit determinant and two-current expectation separately recomputed on nonzero loads. [Exact multiplicity algorithm](nu_multiplicities.py) · [Replay](check_nu_multiplicities.py)
