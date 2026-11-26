"""Fourfold equal-strength graph period from a threefold two-descendant and metric-graph heights."""
from argparse import ArgumentParser
from collections import Counter,defaultdict
from fractions import Fraction as F
from gzip import open as gzip_open
from hashlib import sha256
from itertools import combinations,permutations
import json
from math import comb
from pathlib import Path
from tempfile import NamedTemporaryFile
import networkx as nx
from flint import fmpq,fmpq_mat,fmpz_mat,fmpz_mpoly_ctx

from zero_wave_descendant import change,exact_numerator,particular_flow,tree_map
from current_variation import gradients
from current_minor import numerator
from four_loop_descendant import integrate
from bulk_three import canonical_resistance,forest_data
from path_moments import PathMoments
from point_graph import solve_graph
from quartic_planar_cores import cubic_cores,graph_shapes
from square_three import G,skeleton_chains


def source_derivative(graph,tree,sources,w,pivot,i):
    # One momentum component changes four cycle rows. The determinant is
    # quartic in its strength, so this centered derivative is exact.
    def shifted(t):
        return exact_numerator(graph,tree,change(sources,w,pivot,i,t))
    delta = 8*(shifted(1)-shifted(-1)) - shifted(2)+shifted(-2)
    diff = delta//12
    assert 12*diff == delta
    return diff


def row_derivative(graph,tree,sources,w,pivot,i):
    """Apply determinant multilinearity to the cycle rows of one varied source coordinate."""
    flow = particular_flow(graph,tree,sources)
    zeros = [(0,0,0)]*len(sources)
    shift = particular_flow(graph,tree,change(zeros,w,pivot,i,1))
    values = [v[i] for v in shift]
    return sum(numerator(graph,sources,flow,source_row=(i,r,values)) for r in range(graph.n))


def check_row_derivative(graph,sources,w,pivot,i,polynomial):
    """Exact centered differences of the signed grounded-current minor at two metrics."""
    n = graph.n
    half = (n+1)//2
    for times in [list(range(1,graph.N+1)),[1+(j*j+3*j)%17 for j in range(graph.N)]]:
        gram = fmpz_mat([[sum(int(a)*int(b)*t for a,b,t in zip(r,s,times))
                          for s in graph.C] for r in graph.C])
        scale = fmpq(gram.det()**n,times[0])
        def signed(step):
            flows,_ = solve_graph(graph.edges,change(sources,w,pivot,i,step),times,0)
            rows = [[str(c[e]*flows[e][k]) for e in range(1,graph.N)]
                    for k in range(3) for c in graph.C]
            return scale*fmpq_mat(rows).det()
        slope = sum(fmpq((-1)**(h+1)*comb(2*half,half-h),h*comb(2*half,half))*
                    (signed(h)-signed(-h)) for h in range(1,half+1))
        assert slope == polynomial(*times)


def cached_derivative(directory,graph,sources,w,pivot,i,polynomial=None):
    signature = json.dumps([graph.edges,graph.C,sources,w,pivot,i],separators=(',',':'))
    path = Path(directory)/(sha256(signature.encode()).hexdigest()+'.gradient.gz')
    if polynomial is not None:
        path.parent.mkdir(parents=True,exist_ok=True)
        with NamedTemporaryFile(dir=path.parent,delete=False) as f:
            pending=Path(f.name)
        try:
            with gzip_open(pending,'wt') as out:
                print(signature,file=out)
                for j in range(len(polynomial)):
                    powers=tuple(map(int,polynomial.monomial(j)))
                    assert max(powers)<16
                    print(''.join(format(k,'x') for k in powers),polynomial.coefficient(j),file=out)
            pending.replace(path)
        finally:
            pending.unlink(missing_ok=True)
        return
    with gzip_open(path,'rt') as file:
        if file.readline().strip()!=signature:
            raise ValueError('Gradient cache does not match the oriented graph and sources.')
        terms={}
        for line in file:
            ex,c=line.split()
            powers=tuple(int(k,16) for k in ex)
            assert len(powers)==graph.N
            terms[powers]=int(c)
    return fmpz_mpoly_ctx.get(tuple(f't{j}' for j in range(graph.N))).from_dict(terms)


def period_polynomial(edges, piece=None, left_slice=None, gradient_cache=None, prepare=False, gradient_rows=False):
    net = nx.Graph(edges)
    graph = G(len(net),edges)
    g = graph.n
    branch,chains = skeleton_chains(edges,len(net))
    ordered = []
    for vertices,path in chains:
        i,j = branch.index(vertices[0]),branch.index(vertices[-1])
        if i>j:
            i,j,vertices,path = j,i,vertices[::-1],path[::-1]
        ordered.append(((i,j),vertices,path))
    ordered.sort(key=lambda entry:(entry[0],tuple(entry[2])))
    core,vertices,paths = map(tuple,zip(*ordered))
    ctx,u,forest = forest_data(core,len(branch))
    slots = [v for v,d in net.degree if d==2]
    sources = [(0,0,0)]*len(net)
    for i,p in zip(slots,((1,0,0),(0,1,0),(0,0,1),(-1,-1,-1))):
        sources[i] = p
    tree = tree_map(graph)
    directions = 3*len(branch)
    if piece is not None and (g != 4 or not 0 <= piece <= directions):
        raise ValueError('Select a genus-four current derivative or its final undifferentiated term.')
    if prepare:
        if gradient_cache is None or piece == directions:
            raise ValueError('Gradient preparation requires a cache directory and a derivative term.')
        if piece is not None:
            w,i = branch[piece//3],piece%3
            derivative = row_derivative if gradient_rows else source_derivative
            diff = derivative(graph,tree,sources,w,slots[-1],i)
            if gradient_rows:
                check_row_derivative(graph,sources,w,slots[-1],i,diff)
            cached_derivative(gradient_cache,graph,sources,w,slots[-1],i,diff)
            return core,None
    base = exact_numerator(graph,tree,sources)
    if not base:
        return core,ctx.constant(0)
    if piece is None:
        zero = [(0,0,0)]*len(net)
        variations = [particular_flow(graph,tree,change(zero,w,slots[-1],i,1))
                      for w in branch for i in range(3)]
        regenerated,derivatives = gradients(graph,particular_flow(graph,tree,sources),variations)
        assert regenerated==base
        if prepare:
            for j,polynomial in enumerate(derivatives):
                cached_derivative(gradient_cache,graph,sources,branch[j//3],slots[-1],j%3,polynomial)
            return core,None
        derivatives = iter(derivatives)
    at = {v:(i,j) for i,chain in enumerate(vertices) for j,v in enumerate(chain[1:-1],1)}
    space = PathMoments(base,paths,g-1,ctx=ctx)
    complete = space
    if left_slice is not None:
        from copy import copy
        axis,degree = left_slice
        if piece is None or not 0<=axis<len(paths) or not 0<=degree<g:
            raise ValueError('A left current slice requires one piece, core path and valid degree.')
        space = copy(complete)
        space.base = {}
        space.positions = {}
        for mode,poly in complete.base.items():
            terms = {ex:c for ex,c in poly.terms() if ex[axis]==degree}
            if terms:
                space.base[mode] = ctx.from_dict(terms)
    def resistance_factor(left,right,weights):
        # The forest numerator of each pair is quadratic in the ordered path fractions.
        endpoints = {v:(branch.index(v),)*2 for v in branch}
        endpoints.update({v:core[path] for v,(path,_) in at.items()})
        const = ctx.constant(0)
        linear,quadratic = defaultdict(lambda:ctx.constant(0)),defaultdict(lambda:ctx.constant(0))
        for (v,w),count in weights.items():
            if v==w or not count:
                continue
            i,j = endpoints[v]
            k,l = endpoints[w]
            c = fmpq(str(count))
            const += c*forest[i][k]
            if v in at:
                arm = at[v][0]
                curve = ctx.gen(arm)*u-forest[i][j]
                linear[v] += c*(forest[j][k]-forest[i][k]+curve)
                quadratic[v,v] -= c*curve
            if w in at:
                arm = at[w][0]
                curve = ctx.gen(arm)*u-forest[k][l]
                linear[w] += c*(forest[i][l]-forest[i][k]+curve)
                quadratic[w,w] -= c*curve
            if v in at and w in at:
                quadratic[v,w] += c*(forest[j][l]-forest[i][l]-forest[j][k]+forest[i][k])
                if at[v][0]==at[w][0]:
                    bridge = 2*c*ctx.gen(at[v][0])*u
                    linear[min((v,w),key=lambda h:at[h][1])] -= bridge
                    quadratic[v,w] += bridge
        result = space.inner(left.base,right.base)*const
        for v,coef in linear.items():
            if coef:
                result += space.inner(left.position(*at[v]),right.base)*coef
        for (v,w),coef in quadratic.items():
            if coef:
                result += space.inner(left.position(*at[v]),right.position(*at[w]))*coef
        return space.volume*result

    def add(row,v,w,value):
        row[tuple(sorted((v,w)))] += value

    total = ctx.constant(0)
    base_pairs = Counter()
    for w in [*slots,*branch]:
        derivative_pairs = [Counter() for _ in slots]
        for i,h,k in permutations(range(4),3):
            iw,hw,kw = slots[i],slots[h],slots[k]
            for v,z,c in ((w,kw,1),(w,iw,-1),(hw,iw,1),(hw,kw,-1)):
                add(derivative_pairs[i],v,z,c)
            for v,z,c in ((w,kw,1),(hw,w,1),(hw,kw,-1)):
                add(base_pairs,v,z,F(c,2))
        if w in slots:
            for key,coef in derivative_pairs[slots.index(w)].items():
                base_pairs[key] += g*coef
            continue
        # d_4=g N-d_1-d_2-d_3 by the determinant character.
        for key,coef in derivative_pairs[3].items():
            base_pairs[key] += g*coef
        for i in range(3):
            if piece is not None and 3*branch.index(w)+i != piece:
                continue
            pairs = derivative_pairs[i].copy()
            pairs.subtract(derivative_pairs[3])
            if piece is None:
                diff = next(derivatives)
            elif gradient_cache is not None:
                diff = cached_derivative(gradient_cache,graph,sources,w,slots[-1],i)
            else:
                diff = source_derivative(graph,tree,sources,w,slots[-1],i)
            if diff:
                second = PathMoments(diff,paths,g-1,ctx=ctx)
                total += resistance_factor(space,second,pairs)
    if piece is None or piece == directions:
        total += resistance_factor(space,complete,base_pairs)
        rho,rem = divmod(canonical_resistance(core,ctx,u,forest,0),u)
        assert not rem
        epsilon = (2*(g-1)*rho+sum(forest[i][j] for i,j in combinations(range(len(branch)),2)))/g
        phi = -sum(ctx.gens())*u/4+3*rho/2+epsilon/2
        total -= 24*space.volume*space.inner(space.base,complete.base)*(epsilon+phi)
    aut = sum(1 for _ in nx.isomorphism.GraphMatcher(net,net).isomorphisms_iter())
    return core,total/aut


def derivative_images(edges,core,piece):
    """Reorder the core lengths for each equivalent signed derivative with fixed recoil."""
    if len(set(core))!=len(core):
        raise ValueError('The core edge permutation requires distinct endpoint pairs.')
    net = nx.Graph(edges)
    branch,_ = skeleton_chains(edges,len(net))
    slots = [v for v,d in net.degree if d==2]
    images = {}
    for mapping in nx.isomorphism.GraphMatcher(net,net).isomorphisms_iter():
        if mapping[slots[-1]] != slots[-1]:
            continue
        vertex = branch.index(mapping[branch[piece//3]])
        coordinate = slots.index(mapping[slots[piece%3]])
        target = 3*vertex+coordinate
        moving = [branch.index(mapping[v]) for v in branch]
        order = tuple(core.index(tuple(sorted((moving[i],moving[j])))) for i,j in core)
        if target not in images or order < images[target]:
            images[target] = order
    return images


def diagram(name,index,piece=None,left_slice=None,gradient_cache=None,prepare=False,gradient_rows=False):
    distinct = {'two-double':7,'one-double':8,'prism':9,'bipartite':9}[name]
    core, = (c for c in cubic_cores() if len(set(c))==distinct
             and (distinct<9 or nx.is_bipartite(nx.Graph(c))==(name=='bipartite')))
    edges = list(graph_shapes(core))[index]
    ordered,polynomial = period_polynomial(edges,piece,left_slice,gradient_cache,prepare,gradient_rows)
    if not prepare:
        print('Core-length control',polynomial(*range(1,len(ordered)+1)),flush=True)
        if name in ('bipartite','prism') and piece is not None and piece<18:
            for target,order in sorted(derivative_images(edges,ordered,piece).items()):
                if target!=piece:
                    print('Relabeled core-length control',target,
                          polynomial(*(i+1 for i in order)),flush=True)
    return ordered,polynomial


if __name__=='__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('core',choices=('two-double','one-double','prism','bipartite'))
    parser.add_argument('shape',type=int)
    parser.add_argument('--weights',type=Path)
    parser.add_argument('--piece',type=int,help='One numerator derivative (0–17) or the remaining term (18).')
    parser.add_argument('--left-slice',nargs=2,type=int,metavar=('PATH','DEGREE'))
    parser.add_argument('--gradient-cache',type=Path,help='Directory of graph- and source-specific current derivatives.')
    parser.add_argument('--prepare-gradients',action='store_true',help='Write the adjugate derivatives, or the selected centered derivative with --piece, and exit.')
    parser.add_argument('--gradient-rows',action='store_true',help='Use separate cycle-row determinants for the selected derivative.')
    args=parser.parse_args()
    if args.prepare_gradients:
        if args.gradient_rows and args.piece is None:
            parser.error('--gradient-rows requires a single --piece.')
        diagram(args.core,args.shape,args.piece,args.left_slice,args.gradient_cache,True,args.gradient_rows)
        parser.exit()
    if args.gradient_rows:
        parser.error('--gradient-rows is used with --prepare-gradients.')
    result=integrate(args.core,args.shape,args.weights,
                     builder=lambda name,index:diagram(name,index,args.piece,args.left_slice,args.gradient_cache))
    if result is not None:
        print('Integrated',*result)
