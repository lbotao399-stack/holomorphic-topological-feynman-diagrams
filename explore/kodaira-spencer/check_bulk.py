"""Raw odd delta-form determinants for the four-field minimal KS contact."""
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations, permutations
import json
from pathlib import Path
import sys
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'twoloop-universal'))
import twoloop


def sign(order):
    return (-1)**sum(x > y for i,x in enumerate(order) for y in order[i+1:])


def edge_rows(vertices, edges, y):
    """Factor i_y (dtheta_v-dtheta_w)^3 as a wedge of two rows."""
    rows = []
    for (v,w), moment in zip(edges,y):
        pivot = next((k for k in range(3) if moment[k]),None)
        if pivot is None:
            return None
        i,j = (pivot+1)%3,(pivot+2)%3
        for k,scale in ((i,moment[pivot]),(j,1)):
            row = [F(0)]*(3*vertices)
            for endpoint,s in ((v,scale),(w,-scale)):
                row[3*endpoint+k] = s
                row[3*endpoint+pivot] = -s*moment[k]/moment[pivot]
            rows.append(row)
    return rows


def complement(rows):
    """Eliminate the edge rows only. External fields still enter by a literal 4 by 4 determinant."""
    pivots = []
    weight = F(1)
    for j in range(len(rows[0])):
        h = len(pivots)
        pivot = next((i for i in range(h,len(rows)) if rows[i][j]),None)
        if pivot is None:
            continue
        if pivot != h:
            rows[pivot],rows[h] = rows[h],rows[pivot]
            weight *= -1
        pivots.append(j)
        scale = rows[h][j]
        weight *= scale
        rows[h] = [v/scale for v in rows[h]]
        for i,row in enumerate(rows):
            if i != h and row[j]:
                rows[i] = [a-row[j]*b for a,b in zip(row,rows[h])]
    if len(pivots) < len(rows):
        return None
    free = [j for j in range(len(rows[0])) if j not in pivots]
    assert len(free) == 4
    inverse = {p:i for i,p in enumerate(pivots)}
    residual = [tuple(-rows[inverse[p]][j] for j in free) if p in inverse
                else tuple(int(p==j) for j in free) for p in range(len(rows[0]))]
    return weight*sign(pivots+free),residual


def transverse(p):
    i = next(i for i,x in enumerate(p) if x)
    out = []
    for j in range(3):
        if i != j:
            u = [0]*3
            u[j],u[i] = p[i],-p[j]
            out.append(u)
    return out


def det4(vectors):
    return twoloop.det(vectors)


def loop_response(vertices, edges, times):
    incidence, unused = twoloop.incidence(vertices,edges,vertices-1)
    cycles = sp.Matrix.hstack(*sp.Matrix(incidence).T.nullspace()).T
    return -(cycles*sp.diag(*times)*cycles.T).inv()*cycles


def geometric(response, currents):
    # Differentiating loop voltage conservation gives the independent
    # cycle-current variations. Fix the final Schwinger time.
    mat = [[response[i,j]*currents[j][k] for j in range(response.cols-1)]
           for i in range(response.rows) for k in range(3)]
    return twoloop.det(mat)


def check(edges, at):
    vertices = max(max(e) for e in edges)+1
    loops = len(edges)-vertices+1
    times = [F(1+((i+5)**2%31), i+2) for i in range(len(edges))]
    momenta = [(1,0,1),(0,2,1),(1,-1,1),(-2,-1,-3)]
    basis = [transverse(p) for p in momenta]
    response = loop_response(vertices,edges,times)
    base = None
    bins = defaultdict(F)
    nonzero = set()
    for order in permutations(range(4)):
        src = [[0,0,0] for _ in range(vertices)]
        for label,location in enumerate(order):
            src[at[location]] = momenta[label]
        y = twoloop.currents(vertices,edges,vertices-1,times,src)
        raw = geometric(response,y)
        if base is None:
            base = raw
        assert raw == sign(order)**loops * base
        if not raw:
            continue
        rows = edge_rows(vertices,edges,y)
        reduced = complement(rows) if rows is not None else None
        if reduced is None:
            continue
        factor,residual = reduced
        # Preorder of external rows is by field label, including mixed polyvector degrees.
        forms = []
        for label,location in enumerate(order):
            block = residual[3*at[location]:3*at[location]+3]
            for u in basis[label]:
                forms.append([sum(a*b for a,b in zip(u,col)) for col in zip(*block)])
        for modes in combinations(range(8),4):
            value = factor*det4([forms[k] for k in modes])
            if value:
                nonzero.add(modes)
                bins[modes] += sign(order)**loops * value
    assert all(v == 0 for v in bins.values()), (loops,len(nonzero),len([v for v in bins.values() if v]))
    print(loops, 'loops;',len(nonzero),'exterior components nonzero before permutation; all 70 vanish afterwards')
    return len(nonzero)


if __name__ == '__main__':
    assert check([(0,1),(1,2),(2,3),(3,0)],list(range(4)))
    for shape in ((1,3,3),(2,2,3)):
        vertices,edges = twoloop.theta_graph(*shape)
        at = [i for i in range(vertices) if sum(i in e for e in edges) == 2]
        assert check(edges,at)

    sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'nonSP'))
    from graphs import K4SUB
    for name,edges in K4SUB:
        at = [i for i in range(8) if sum(i in e for e in edges)==2]
        check(edges,at)
    if '--four' in sys.argv:
        catalog = Path(__file__).resolve().parents[2]/'deepen/o3/higher_orders/four-loop-graphs.json'
        for graph in json.loads(catalog.read_text()):
            edges = [e for chain in graph['chains'] for e in zip(chain,chain[1:])]
            check(edges,list(range(6,10)))
