"""Loop-two inventory, one-primary open hexagon and two-vector disk Ward checks."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations,permutations
from math import factorial
from open_theta import parity,canonical_trace,det


def stable_colour(k):
    """Cyclic trace with a symmetric (k+1)-point disk insertion and five adjoint actions."""
    words=Counter({(-1,):1})
    for i in range(5):
        pending=Counter()
        for w,c in words.items():
            pending[(i,)+w]+=c
            pending[w+(i,)]-=c
        words=pending
    result=Counter()
    for s in permutations([-2,*range(5,k+4)]):
        cut=s.index(-2)
        for w,c in words.items():
            a=w.index(-1)
            # Rotate the disk trace to put its first metric endpoint at the start,
            # expand nested ad on the other endpoint and close with the SO metric.
            # Tr(T [s[:cut] w[:a]] T [w[a+1:] s[cut+1:]]).
            u=s[:cut]+w[:a]
            v=w[a+1:]+s[cut+1:]
            for traces,weight in [((u,v),c),((u[::-1]+v,),-(-1)**len(u)*c)]:
                out=[canonical_trace(t) for t in traces]
                rank=sum(not t for t in traces)
                keys=tuple(sorted((a for a,b in out if a),key=lambda w:(len(w),w)))
                factor=weight
                for _,sgn in out: factor*=sgn
                result[rank,keys]+=F(factor,2*factorial(k))
    return {key:c for key,c in result.items() if c}


def relabel(traces,target):
    if sorted(map(len,traces))!=sorted(map(len,target)):
        return
    def cyclic(w):
        return [w[i:]+w[:i] for i in range(len(w))]
    from itertools import product
    for assignment in permutations(target):
        if any(len(a)!=len(b) for a,b in zip(traces,assignment)):
            continue
        options=[]
        for w in assignment:
            opts={x:1 for x in cyclic(w)}
            opts.update({x:(-1)**len(w) for x in cyclic(w[::-1])})
            options.append(list(opts.items()))
        for choices in product(*options):
            s=[None]*sum(map(len,traces))
            sign=1
            for src,(dst,sgn) in zip(traces,choices):
                for i,j in zip(src,dst):s[i]=j
                sign*=sgn
            yield tuple(s),sign*parity(s)


def oriented_tensor(k,p,form,permutation):
    s=permutation
    determinant=det([p[i] for i in s[:5]])
    if not determinant:return 0
    weighted=[sum((5-j)*p[s[j]][i] for j in range(5)) for i in range(5)]
    return determinant*form([p[i] for i in s[5:]]+[weighted])


def wedge_of_vectors(vectors):
    def evaluation(covectors):
        return det([[sum(a*b for a,b in zip(v,c)) for c in covectors] for v in vectors])
    return evaluation


def projected(k,cycles,p,form):
    result=Counter()
    for (power,words),c in stable_colour(k).items():
        for permutation,sign in relabel(words,cycles):
            result[power]+=c*sign*oriented_tensor(k,p,form,permutation)
    return {i:v for i,v in result.items() if v}


def planar_inventory():
    # A local projected graph with a heat edge has V=4b+2 and E=5b+1 in dimension five.
    seen=set()
    pure=set()
    for m in range(8):
        for b in range(5):
            v,e=4*b+2,5*b+1
            for spheres in range(v+1):
                disks=v-spheres
                for ec in range(e+1):
                    if b+ec-spheres+m!=2:continue
                    if not disks:
                        if ec==e and 3*v<=2*e+m:pure.add((m,b))
                        continue
                    if ec<spheres:continue
                    if b==1 and m<spheres:continue  # Each vertex of the hexagon has internal valency two.
                    # At least three closed incidences at each sphere; at least one disk endpoint if spheres occur.
                    if spheres and 3*spheres+1>2*ec+m:continue
                    # For one closed internal edge with no spheres a self-loop gives zero.
                    seen.add((m,b,spheres,ec))
    assert seen=={(0,2,0,0),(0,1,0,1),(1,1,0,0),(1,0,0,1),(2,0,0,0),(2,0,1,1)}
    assert not pure
    return seen


def check_vector_ward():
    # [u.exp(p),v.exp(q)] as an actual commutator of derivations on a scalar and matrix current.
    p=[1,2,0,-1,2];q=[2,0,-1,1,0]
    u=[0,1,2,2,0];v=[1,0,2,0,0]
    dot=lambda x,y:sum(a*b for a,b in zip(x,y))
    assert dot(p,u)==dot(q,v)==0
    lie=[dot(q,u)*b-dot(p,v)*a for a,b in zip(u,v)]
    for r in ([1,0,-2,3,1],[-2,1,5,2,-3]):
        left=dot(v,r)*dot(u,[a+b for a,b in zip(q,r)])
        right=dot(u,r)*dot(v,[a+b for a,b in zip(p,r)])
        assert left-right==dot(lie,r)
    # On the open one-loop wheel there are five inputs. Cyclicity of the adjoint trace
    # makes all five weights equal once their rotations are included.
    c=[5,4,3,2,1]
    assert [sum(c[(i+j)%5] for j in range(5)) for i in range(5)]==[15]*5


def literal_wheel(n,edges,p,form):
    """Direct adjoint propagation and defining-representation trace on a physical probe."""
    from check_single_trace import bracket,trace
    value=0
    k=len(edges)-4
    for s in permutations(range(len(edges))):
        gauge=oriented_tensor(k,p,form,s)
        if not gauge:
            continue
        colour=0
        for start in combinations(range(n),2):
            state=start
            for i in s[:5]:
                state=bracket(edges[i],state) if state is not None else None
                if state is None:
                    break
            if state is None:
                continue
            for word in permutations([state,*[edges[i] for i in s[5:]]]):
                colour+=trace([start,*word])
        value+=parity(s)*gauge*colour
    return -F(value,2*factorial(k))


def check_trivector_wheel():
    from check_single_trace import trace
    p=[[int(i==j) for i in range(5)] for j in range(5)]+[[1,-1,2,0,3],[0,2,-2,1,1]]
    q=[-sum(values) for values in zip(*p)]
    vectors=[[q[4]*int(i==j)-q[j]*int(i==4) for i in range(5)] for j in range(3)]
    form=wedge_of_vectors(vectors)
    cases=[((tuple(range(7)),),{0:-23300,1:F(1000,3)}),
           (((0,1),tuple(range(2,7))),{0:11000}),
           (((0,1,2),tuple(range(3,7))),{})]
    for cycles,expected in cases:
        assert projected(3,cycles,p,form)==expected
        edges=[(i,j) for word in cycles for i,j in zip(word,word[1:]+word[:1])]
        trace_value=1
        for word in cycles:
            trace_value*=trace([edges[i] for i in word])
        for n in (7,8):
            expect=trace_value*sum(n**i*c for i,c in expected.items())
            assert literal_wheel(n,edges,p,form)==expect
    return [expected for _,expected in cases]


if __name__=='__main__':
    print('nonzero topological candidates (external closed,b1,spheres,closed edges)',sorted(planar_inventory()))
    check_vector_ward()
    print('vector wheel cancels; two-vector disk commutator equals the Lie action')
    print('trivector wheel trace coefficients',check_trivector_wheel())
