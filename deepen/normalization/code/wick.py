"""BV Wick sign of a labelled graph (superfield rule), independent implementation.

Model (derivation in REPORT): every vertex field is an odd superfield Phi (total degree), the contraction of two
fields is <Phi^a(x) Phi^b(y)> = hbar kappa^{ab} P(x - y) delta^k(theta - theta') with P(x-y) delta even
(P an even form in all cases used here: (H,T,k) = (2,1,0), (3,0,2), (2,1,2)); the raw graph integral is built with
P_e(x_{e(0)} - x_{e(1)}).  A cubic vertex contributes f_{h1 h2 h3} Phi_h1 Phi_h2 Phi_h3 for any order (h1,h2,h3) of its
half-edges (f totally antisymmetric, Phi odd).  A field contracted with a letter of the operator O at the root is
replaced in place by the odd functional derivative delta/delta C (left derivative), with kernel P(x_v - 0); the raw
integral uses P_e(x_{e(0)} - x_{e(1)}), which is -P(x_v - 0) for a spoke oriented O -> v.

sigma = sgn( source -> target ) * (-1)^{# spokes with tail O},
source = concatenation over the interaction vertices (vertex order) of their half-edges in local order,
target = [ (e,0), (e,1) for each internal edge e in edge order ] + [ outputs in the given order ] + [ spoke fields in delta order ].
Half-edges: ('e', e, end) with end 0 = tail e(0), 1 = head e(1);  ('x', v) = output / external leg at v.
"""


def perm_sign(seq):
    seq = list(seq)
    s = 1
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] > seq[j]:
                s = -s
    return s


def wick_sign(vertices, edges, local, outputs, spokes_delta=(), O='O'):
    """vertices: interaction vertices in vertex order; edges: list of (tail, head) (O may appear);
    local: dict v -> list of half-edges; outputs: list of ('x', v); spokes_delta: list of edge indices (delta order)."""
    source = [h for v in vertices for h in local[v]]
    target = []
    for e, (a, b) in enumerate(edges):
        if a != O and b != O:
            target += [('e', e, 0), ('e', e, 1)]
    target += list(outputs)
    nsp_out = 0
    for e in spokes_delta:
        a, b = edges[e]
        assert (a == O) != (b == O)
        target.append(('e', e, 1) if a == O else ('e', e, 0))
        nsp_out += (a == O)
    assert sorted(map(str, source)) == sorted(map(str, target)), (source, target)
    pos = {str(h): i for i, h in enumerate(target)}
    return perm_sign([pos[str(h)] for h in source]) * (-1) ** nsp_out


def pattern_local(vertices, edges, outputs_at, O='O'):
    """S2 leg pattern: in-edges in edge order, output, out-edges in edge order."""
    loc = {}
    for v in vertices:
        ins = [('e', e, 1) for e, (a, b) in enumerate(edges) if b == v]
        outs = [('e', e, 0) for e, (a, b) in enumerate(edges) if a == v]
        loc[v] = ins + ([('x', v)] if v in outputs_at else []) + outs
    return loc


def explorer_LLp_sign(nv, edges, ext):
    """sgn(L -> L') of anomaly-structure cochain.colour_weight, re-implemented from its docstring:
    L = (x_{v_1}..x_{v_k}, (e0,out),(e0,in), (e1,out),(e1,in), ...),
    L' = concatenation over v = 0..nv-1 of [x_v if external] + [half-edges at v sorted by (edge, end)]."""
    L = [('x', v) for v in ext]
    for e in range(len(edges)):
        L += [('e', e, 0), ('e', e, 1)]
    Lp = []
    for v in range(nv):
        Lp += [('x', v)] if v in ext else []
        Lp += sorted([('e', e, s) for e, ab in enumerate(edges) for s in (0, 1) if ab[s] == v], key=lambda h: (h[1], h[2]))
    pos = {h: i for i, h in enumerate(L)}
    return perm_sign([pos[h] for h in Lp]), Lp


def explorer_local(nv, edges, ext):
    loc = {}
    for v in range(nv):
        loc[v] = ([('x', v)] if v in ext else []) + sorted(
            [('e', e, s) for e, ab in enumerate(edges) for s in (0, 1) if ab[s] == v], key=lambda h: (h[1], h[2]))
    return loc
