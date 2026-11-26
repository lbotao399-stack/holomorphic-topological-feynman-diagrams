"""Projective current numerator from a homogeneous linear matrix of size n^2-1."""
from math import lcm
from fractions import Fraction
from flint import fmpq_mat, fmpz_mat, fmpz_mpoly_ctx
import networkx as nx


def bareiss_det(matrix,ctx):
    """Fraction-free elimination with row and column pivots on the sparse linear minor."""
    rows = [list(row) for row in matrix]
    previous = ctx.constant(1)
    sign = 1
    while len(rows)>1:
        row_counts = [sum(bool(entry) for entry in row) for row in rows]
        col_counts = [sum(bool(entry) for entry in col) for col in zip(*rows)]
        nonzero = ((len(entry),row_counts[i]*col_counts[j],i,j)
                   for i,row in enumerate(rows) for j,entry in enumerate(row) if entry)
        chosen = min(nonzero,default=None)
        if chosen is None:
            return ctx.constant(0)
        _,_,i,j = chosen
        if i:
            rows[i],rows[0] = rows[0],rows[i]
            sign = -sign
        if j:
            for row in rows:
                row[j],row[0] = row[0],row[j]
            sign = -sign
        pivot = rows[0][0]
        rows = [[(pivot*x-row[0]*y)//previous for x,y in zip(row[1:],rows[0][1:])]
                for row in rows[1:]]
        previous = pivot
    return sign*rows[0][0]


def polynomial_det(matrix,ctx):
    size = len(matrix)
    if size>8:
        support = nx.Graph()
        support.add_nodes_from(range(2*size))
        support.add_edges_from((i,size+j) for i,row in enumerate(matrix)
                               for j,v in enumerate(row) if v)
        matching = nx.bipartite.maximum_matching(support,range(size))
        if len(matching) != 2*size:
            return ctx.constant(0)
        # Match each row to a diagonal column. Strong components are the
        # diagonal blocks of a simultaneous block-triangular permutation.
        permutation = [matching[i]-size for i in range(size)]
        graph = nx.DiGraph()
        graph.add_nodes_from(range(size))
        graph.add_edges_from((i,matching[j]) for i in range(size) for j in support[i])
        blocks = list(nx.strongly_connected_components(graph))
        if len(blocks)>1:
            sign = (-1)**sum(x>y for i,x in enumerate(permutation) for y in permutation[i+1:])
            result = sign*ctx.constant(1)
            for block in blocks:
                ids = sorted(block)
                result *= polynomial_det([[matrix[i][permutation[j]] for j in ids] for i in ids],ctx)
                if not result:
                    break
            return result
        return bareiss_det(matrix,ctx)
    terms = {0:ctx.constant(1)}
    for row in matrix:
        next_terms = {}
        for mask,minor in terms.items():
            for j,value in enumerate(row):
                if mask>>j&1 or not value:
                    continue
                new = mask | (1<<j)
                signed = value*minor
                if (mask>>(j+1)).bit_count()%2:
                    signed = -signed
                next_terms[new] = next_terms.get(new,0)+signed
        terms = next_terms
    return terms.get((1<<size)-1,ctx.constant(0))


def numerator(graph, source, flow, omitted_rows=(), *, source_row=None, edge_times=None, determinant=None):
    """Return U^n times the oriented projective cofactor for a particular tree flow.

    An omitted (component,cycle) row is filled by a unit column at zero time.
    The sign follows the augmented matrix with unit columns appended in input order.
    A source_row=(component,cycle,values) replaces that single row's particular
    edge currents in both the constant and length-dependent blocks. Pass edge_times
    to specialize the segment lengths in a common polynomial ring before elimination.
    """
    n,ne,h = graph.n,graph.N,len(source[0])
    if ne+len(omitted_rows)!=n*h+1:
        raise ValueError('Expected one more edge than cycle-current components.')
    if edge_times is None:
        ctx = fmpz_mpoly_ctx.get(tuple(f't{i}' for i in range(ne)))
        times = ctx.gens()
    else:
        if len(edge_times) != ne:
            raise ValueError('Expected one specialized polynomial per segment.')
        times = edge_times
        ctx = times[0].context()
    # The first n^2 rows use matrix order; all source rows use (component,cycle).
    def slot(r,j):
        return n*r+j if j<n else n*n+(j-n)*n+r
    extended = [list(c)+list(y) for c,y in zip(zip(*graph.C),flow)]
    dim = n*(n+h)
    columns = [[0]*dim for _ in range(ne)]
    for e,b in enumerate(extended):
        for r in range(n):
            for j,v in enumerate(b):
                columns[e][slot(r,j)] = b[r]*v
    if source_row is not None:
        component,cycle,values = source_row
        if omitted_rows or not (0<=component<h and 0<=cycle<n and len(values)==ne):
            raise ValueError('Expected one full source row and an edge current per edge.')
        replace = slot(cycle,n+component)
        for e,b in enumerate(extended):
            columns[e][replace] = b[cycle]*values[e]
    if len(set(omitted_rows))!=len(omitted_rows) or any(
            not (0<=component<h and 0<=r<n) for component,r in omitted_rows):
        raise ValueError('Invalid omitted cycle-current rows.')
    for component,r in omitted_rows:
        unit = [0]*dim
        unit[slot(r,n+component)] = 1
        columns.append(unit)
    echelon,rank = fmpq_mat(columns).rref()
    if rank!=len(columns):
        return ctx.constant(0)
    pivot = [next(j for j in range(dim) if echelon[i,j]) for i in range(len(columns))]
    complement = [j for j in range(dim) if j not in pivot]
    det = int(fmpz_mat([[row[j] for j in pivot] for row in columns]).det())
    sign = (-1)**(sum(j-i for i,j in enumerate(pivot))+(n*n-1)*len(columns))
    denominator = 1
    matrix = []
    for j in complement:
        weights = [Fraction(str(echelon[i,j])) for i in range(len(columns))]
        scale = lcm(*(v.denominator for v in weights))
        denominator *= scale
        null = [0]*dim
        null[j] = scale
        for pos,v in zip(pivot,weights):
            null[pos] = -int(v*scale)
        row = []
        for r in range(n):
            contractions = [sum(null[slot(r,k)]*v for k,v in enumerate(b)) for b in extended]
            if source_row is not None and r==cycle:
                contractions = [c+null[replace]*(v-b[n+component])
                                for c,v,b in zip(contractions,values,extended)]
            for s in range(n):
                if (r,s)!=(n-1,n-1):
                    row.append(sum((b[s]*v)*t for b,v,t in zip(extended,contractions,times)))
        matrix.append(row)
    top = sign*det*(determinant or polynomial_det)(matrix,ctx)
    result = top//denominator
    assert result*denominator==top
    return result
