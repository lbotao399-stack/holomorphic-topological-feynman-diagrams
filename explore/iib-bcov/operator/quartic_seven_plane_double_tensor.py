"""All three quadratic tangential coefficients from two-double bare four-loop diagrams."""
from collections import defaultdict
from fractions import Fraction as F
from quartic_seven_plane_double import (
    graph_data, image, orders, polynomials, shapes, symmetries, tree_flow,
)
from quartic_seven_plane_double_three_two import MARKS, integrate_polynomials
from quartic_seven_plane_lift import NORMAL
from quartic_seven_plane_gale import canonical_sources, cubic_nodes, gale as first_four_minors, monomials
from quartic_seven_plane_space import as_matrix, gale, normal_row
from planar_normal_tensors import tensor_row
from current_minor import numerator


def frame(c):
    return ((1,0,0),(0,1,0),(0,0,-1),(-1,0,1-c),(0,-1,c))


def all_shapes():
    return list(shapes())+[(m,len({tuple(sorted(g[h][0] for h in m)) for g in symmetries()})) for m in MARKS]


def cubic_coefficients(data):
    # Indexed by a cubic in the four Gale coordinates, in increasing monomial order.
    graph,ctx=data[:2]
    nodes=cubic_nodes()
    inv=(6*as_matrix([monomials(m) for m in nodes]).inv()).tolist()
    sampled=[]
    for m in nodes:
        src,scale=canonical_sources(m)
        flow=tree_flow(graph,src)
        sampled.append([scale*(-1)**r*numerator(graph,src,flow,[(0,j) for j in range(4) if j!=r])
                        for r in range(4)])
    return [[sum((int(v)*row[r] for v,row in zip(coeff,sampled) if v),ctx.constant(0))
             for r in range(4)] for coeff in inv]


def buckets(marks,q):
    # Put each order into the least Gale tuple along the geometric stabilizer.
    # Equal tuples have identical squared current numerators even when the
    # incidence of the measured labels in the network differs.
    z=gale(q)
    group=[a for a in symmetries() if sorted(a[h][0] for h in marks)==list(marks)]
    buckets={}
    for order,mult in orders(marks):
        config=tuple(tuple(order[i] for i,h in enumerate(marks) if h==j) for j in range(9))
        choices=(sum(image(config,a),()) for a in group)
        x=min(choices,key=lambda g:(tuple(z[h] for h in g),g))
        key=tuple(z[h] for h in x)
        _,resistance=buckets.setdefault(key,(x,defaultdict(int)))
        for i,sign in ((x.index(0),-1),(x.index(1),1)):
            for j,label in enumerate(x):
                if i!=j and q[label][2]:
                    resistance[tuple(sorted((i,j)))]+=mult*sign*q[label][2]
    return buckets.values()


def one_shape(marks,size,q):
    data=graph_data(marks)
    _,ctx,_,_,forests,_,_=data
    coeff=cubic_coefficients(data)
    ans=F()
    for order,resistance in buckets(marks,q):
        if not any(resistance.values()):
            continue
        twice=sum((c*forests[6+i][6+j] for (i,j),c in resistance.items() if c),0)
        if not twice:
            continue
        assert all((v/2).denominator==1 for _,v in twice.terms())
        voltage=ctx.from_dict({k:int(v/2) for k,v in twice.terms()})
        m=first_four_minors([q[i] for i in order])
        slopes=[(1,-1,0,0,0)[j] for j in order]
        d=sum(a*b for a,b in zip(m,slopes))
        values=monomials(m)
        ds=[d*sum((v*row[r] for v,row in zip(values,coeff) if v),ctx.constant(0))//6
            for r in range(4)]
        if order==tuple(range(5)):
            literal,_=polynomials(data,order,normal=q)
            assert ds==literal
        ans+=integrate_polynomials(data,ds,voltage)
    return size*ans/len(symmetries())


ORIGINAL=F(22854171827,509238576000000)
SPARSE={2:F(2693233,354882705408000000),3:F(51018137,638788869734400000),
        4:F(1860953,4655895552000000)}


def recover():
    inputs=[frame(2),frame(3),NORMAL]
    raw=[SPARSE[2],SPARSE[3],ORIGINAL]
    z=[gale(q) for q in inputs]
    # q_0 cross q_1 is parallel to the measured third vector on these fixtures.
    # Matchings on the seven waves and e1 is four times the planar factor of this integral.
    rhs=as_matrix([[v*(q[0][0]*q[1][1]-q[0][1]*q[1][0])/(4*(a[0]-a[1])**2)] for q,a,v in zip(inputs,z,raw)])
    matrix=as_matrix([tensor_row(a,(0,1),7) for a in z])
    assert F(str(matrix.det()))==F(-64012032,5)
    coefs=tuple(F(str(x)) for x in (matrix.inv()*rhs).entries())
    assert coefs==(F(-129013,5110310957875200000),F(-777251,1022062191575040000),F(-19109,69686058516480000))
    direct=as_matrix([normal_row(a,(0,1),2) for a in z])
    assert direct.det()==3110984755200
    assert [F(str(x)) for x in (direct.inv()*rhs).entries()]==[
        F(-5201573,459927986208768000000),F(-5113441,76654664368128000000),F(-502991,98555997044736000000)]
    return coefs


if __name__=='__main__':
    from argparse import ArgumentParser
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('--frame',type=int,choices=(2,3,4))
    parser.add_argument('--shape',type=int,choices=range(11))
    args=parser.parse_args()
    coeff=recover()
    print('quadratic tensor by normal generators',*coeff,flush=True)
    if args.frame is not None:
        q=frame(args.frame)
        total=F()
        for i,(m,size) in enumerate(all_shapes()):
            if args.shape is None or args.shape==i:
                value=one_shape(m,size,q)
                print('shape',i,m,value,flush=True)
                total+=value
        if args.shape is None:
            z=gale(q)
            expected=4*(z[0]-z[1])**2*sum(v*c for v,c in zip(tensor_row(z,(0,1),7),coeff))
            assert total==expected==SPARSE[args.frame]
            print('complete direct frame',args.frame,total,flush=True)
