"""Four-loop ordered products of cubic interactions in the author's C² twist."""
from fractions import Fraction as F
from itertools import permutations
import json
from pathlib import Path

from flint import fmpz_mat
import networkx as nx
import numpy as np

from contractions import basis, exact_product, record, scaled_coefficients, substitution
from homogeneous import exact_kernel


def placements(g):
    maps=list(nx.isomorphism.GraphMatcher(g,g).isomorphisms_iter())
    seen=set()
    for a,b in permutations(range(6),2):
        if not g.has_edge(a,b) or (a,b) in seen:
            continue
        seen.update((iso[a],iso[b]) for iso in maps)
        forest=g.subgraph(set(g)-{a,b})
        if not nx.is_forest(forest) or any(g.degree(v)>3 for v in forest):
            continue
        roots=sorted(v for v in forest if g.degree(v)==2)
        if all(len(set(roots)&c)==1 for c in nx.connected_components(forest)):
            yield a,b,roots


def certificate():
    data,graphs=record()
    den,kernels=scaled_coefficients(data)
    invariant,_,_,lift=exact_kernel()
    offsets=np.cumsum([0]+[m.shape[1] for m in invariant])
    rows=[]
    result=[]
    for i,g in enumerate(graphs):
        for a,b,roots in placements(g):
            # Keep momenta in the order A, roots, and eliminate the B source.
            labels=[a]+roots
            images=tuple(tuple(int(v==j)-int(v==b) for j in labels) for v in range(5))
            restriction=substitution(images)
            physical=exact_product(restriction,kernels[i])
            assert not np.any(physical[:,1:])
            block=restriction @ (invariant[i] @ lift[offsets[i]:offsets[i+1]])
            rows.append(block)
            table=[[s,str(F(int(c),den))] for s,c in zip(basis(len(labels)),physical[:,0]) if c]
            result.append({'catalogue_index':i,'A':a,'B':b,'roots':roots,
                           'coefficients_in_catalogue_edge_order':table,
                           'homogeneous_rank':fmpz_mat(block.tolist()).rank()})
    assert len(result)==23
    rank=fmpz_mat(np.concatenate(rows).tolist()).rank()
    print('product forests',len(result),'nonzero',sum(bool(r['coefficients_in_catalogue_edge_order'])
          for r in result),'homogeneous rank',rank,flush=True)
    return {'joint_homogeneous_rank':rank,'forests':result}


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--save',action='store_true')
    args=parser.parse_args()
    data=json.dumps(certificate(),separators=(',',':'))+'\n'
    path=Path(__file__).with_name('product_forests.json')
    if args.save:path.write_text(data)
    else:assert path.read_text()==data
