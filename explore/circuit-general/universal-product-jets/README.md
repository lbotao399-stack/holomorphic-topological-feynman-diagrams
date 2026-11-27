# Complete product tensors with arbitrary finite input jets

The three-loop product in the published [universal quartic calculation](../universal-quartic-product/README.md) now has an **exact rational evaluator for arbitrary simultaneous finite jets on all four contracted input fields**. It retains all original sources, all input/output color indices, and every free spinor index. The two difficult graph combinations are integrated through their complete source and time boundaries; no individual complete-four-vertex period is needed.

For the polarized pure-matter quartic interaction, the inputs have field contents

\[
\beta^{1i_1}\beta^{2i_2}\gamma_5,
\qquad b_A\beta^{3i_3},
\]

with arbitrary finite holomorphic derivatives separately on each displayed field. The spectator \(\gamma_5\) is free, gauge neutral, and absent from the interaction. The full output at the computed \(g^2\bar W\hbar^3\) order is supported on \(\gamma_{4i_4}c^B\gamma_5\), with the complete derivative tensor. No gauge group, representation dimension, invariant superpotential tensor, or input color contraction is selected.

## Completed results

| Deliverable | Exact scope |
| --- | --- |
| [Rational product evaluator](rational_product.py) | Every finite polynomial word in the four original input currents, with all ten color words after the representation identity |
| [All-degree geometric proof](../forest-jet-identities/DERIVATION.md) | Two complete graph-pair functionals expressed as finite rational boundary integrals |
| [Original-current routing certificate](routing_certificate.json) | All sixteen original commutator matchings, all four differentiated input fields, both gauge-ghost orderings |
| [Complete second-jet tensor](FULL_SECOND_JET.md) | Every source and both free spinor indices; all ten rational color tensors and both full representation-defect coefficients |
| [Complete third-jet tensor](FULL_THIRD_JET.md) | Every source and all three free spinor indices; all ten rational color tensors and both full representation-defect coefficients |
| [Rank-two certificate](jet_order_2_certificate.json) | Original six-graph integrations at all 21 source nodes and five additional complete source assignments |
| [Rank-three certificate](jet_order_3_certificate.json) | Original six-graph integrations at all 28 source nodes and five additional complete source assignments |

The source-simplex proof handles every finite input derivative order, not just the two printed tensors. For arbitrary pure \(\bar W(\gamma)\) of degree at least three, the existing complete polarization classification extends the rationality theorem to this same fixed field-content family. The new API directly evaluates the polarized quartic coefficient. The scope and the higher-degree argument are in [DERIVATION.md](DERIVATION.md).

## Exact API

Install the shared dependencies in [requirements.txt](requirements.txt), then run from the repository root:

```bash
python -m pip install -r explore/circuit-general/universal-product-jets/requirements.txt
python explore/circuit-general/universal-product-jets/run_checks.py
```

For a complete tensor with two derivatives on the original second beta field:

```python
import sys
sys.path.insert(0, "explore/circuit-general/universal-product-jets")
from rational_product import RationalJetProduct

with RationalJetProduct(((1, 2), (3, -1), (-2, 4))) as product:
    tensor = product.tensor({"beta2": 2})
```

The three source rows are, in order, the full spinors of \(\gamma_4,c,\gamma_5\). Sources must be exact rationals. `tensor()` returns every entry separately symmetric within each input's derivative indices, and every one of the ten open-color words, including zero coefficients. Entry arrays store the complete tensor; they do not select a spinor direction or contract free indices. The complete covariant printed tensors avoid entry notation altogether.

`tensor({"beta1": 1, "beta2": 2, "beta3": 1, "b": 1})` retains all four groups of free derivative indices. `moment({("beta1", 0): 1, ("b", 1): 2})` is the exact storage interface for one specified entry. `pair_moment(polynomial)` evaluates both entire geometric graph-pair functionals on a polynomial in the three common full currents; six exponents are required for each monomial.

The physical common multiplier, after replacing each displayed source by the derivative acting on its own surviving field, is

\[
\frac{g^2}{128(8\pi^2)^3}.
\]

The conversion \(\lambda=-\partial/2\) is already included. A `gamma5` jet means the **physical derivative on the free spectator applied after this conversion**. It has no propagator-head factor. The corresponding source factor in the returned polynomial records that physical spectator derivative; it is not an additional raw current insertion.

`max_terms` bounds polynomial expansions. `max_entries` bounds complete `tensor()` output and the number of cached color-word values. A single `moment()` always returns its fixed ten words and can run with a smaller cache budget without caching. A failed complete-tensor or polynomial budget raises an exception; no truncated answer is returned. These are coefficient-count budgets, not bounds on rational numerator bit sizes or total machine memory. There is no fixed derivative-order cutoff.

## Reproducible independent checks

The default runner verifies original-source routing, complete covariant restoration, exact Markdown coefficient rendering, all 62 time-face certificates, and the rational evaluator against every node and every entry of the two direct-period tensor certificates. It also compares all sixteen entries of each of the two independently integrated mixed-input rank-four graph-pair tensors. During the new evaluator checks, the individual-period worker is disabled.

The additional exact identity runner checks all 720 ordered time sectors, all 62 leading faces on both source triangles, both full outer-current differential corrections, and forty fresh original-graph pair entries at another complete source assignment:

```bash
python explore/circuit-general/forest-jet-identities/run_checks.py
```

To regenerate all original six-graph integrals for the two full-source tensors before comparing their complete certificates:

```bash
python explore/circuit-general/universal-product-jets/run_checks.py --recompute
```

The [check log](CHECKS.log) records the completed replay. No period table is built or modified.

## Limits and credit

The original Wick ledger, exact geometric periods, representation convention, pure-superpotential classification, and field-count termination proof come from the published universal quartic calculation and the manuscript. The new work supplies complete higher-jet tensors, verifies all simultaneous original-current routings, and replaces the hard graph combinations by an executable rational boundary integration.

The theorem concerns this fixed input field content and pure \(\bar W(\gamma)\), with a free neutral spectator. It does not assert rationality of every product forest, rationality of the individual three-loop graphs, or a general result for derivative superpotentials or \(\bar W(\gamma,w)\). It uses no assumption about independence of Riemann-zeta values. The all-order termination statement for these fixed inputs is distinct from displaying every lower-order product coefficient.
