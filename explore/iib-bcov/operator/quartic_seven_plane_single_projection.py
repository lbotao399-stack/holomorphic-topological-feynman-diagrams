"""Integrate the single-double bubble before expanding a polynomial on the tetrahedron."""
from collections import defaultdict
from fractions import Fraction as F
from functools import cache
from itertools import combinations_with_replacement, permutations
from math import factorial, prod
from pathlib import Path
import json

from flint import fmpq, fmpq_mat, fmpq_mpoly_ctx, fmpz_mpoly_ctx
from quartic_seven_plane_current_basis import expand, integer_poly, polarized_currents
from quartic_seven_plane_double import Q, SCALE, local_moment, primitive
from quartic_seven_plane_single import CORE, cuts, graph_data
from quartic_planar_cores import constraints
from fourfold_bulk import symmetries
from square_three import orbit


def split_variables(poly,nleft):
    """Separate bubble and tetrahedron variables at each total bubble degree."""
    slices=defaultdict(dict)
    for mon,c in poly.terms():
        m=tuple(map(int,mon))
        slices[sum(m[:nleft])][m[:nleft],m[nleft:]]=int(c)
    for degree,entries in slices.items():
        left=sorted({a for a,b in entries})
        right=sorted({b for a,b in entries})
        matrix=fmpq_mat([[entries.get((a,b),0) for a in left] for b in right])
        rows,rank=matrix.rref()
        for i in range(rank):
            pivot=next(j for j in range(len(left)) if rows[i,j])
            l,c=primitive([(m,F(str(rows[i,k]))) for k,m in enumerate(left) if rows[i,k]])
            r,d=primitive([(m,F(str(matrix[k,pivot]))) for k,m in enumerate(right) if matrix[k,pivot]])
            yield l,r,degree,c*d


class Bubble:
    def __init__(self,n,split):
        names=tuple(f'x{i}' for i in range(n))
        self.ctx=fmpz_mpoly_ctx.get(names)
        self.rational=fmpq_mpoly_ctx.get(names)
        self.split=split
        self.ids={}
        self.polys=[]
        self.pairs={}
        self.triples={}
        self.singular_basis=[]

    def register(self,terms):
        if terms not in self.ids:
            self.ids[terms]=len(self.polys)
            self.polys.append(self.ctx.from_dict(dict(terms)))
        return self.ids[terms]

    def project(self,i,j,k):
        """Exact local integral and coordinates of every formally divergent term."""
        key=(min(i,j),max(i,j),k)
        if key not in self.triples:
            if key[:2] not in self.pairs:
                self.pairs[key[:2]]=self.polys[i]*self.polys[j]
            polynomial=self.pairs[key[:2]]*self.polys[k]
            excluded={}
            value=0
            n,m=self.split
            for ex,c in polynomial.terms():
                ex=tuple(map(int,ex))
                a,b=sum(ex[:n])+n,sum(ex[n:n+m])+m
                if min(Q-a,Q-b,a+b-Q)<=0:
                    excluded[ex]=int(c)
                else:
                    value+=int(c)*local_moment(self.split,ex)
            coordinates=expand(self.singular_basis,self.rational.from_dict(excluded)) if excluded else ()
            self.triples[key]=(value,coordinates)
        return self.triples[key]


@cache
def orbital(powers):
    return orbit(powers)


@cache
def factorials(powers):
    return prod(map(factorial,powers))


def integrated_weights(marks,orders=None,*,fold=True,**kinematics):
    data=graph_data(marks)
    bases,voltages=polarized_currents(data,permutations(range(5)) if orders is None else orders,**kinematics)
    return project_currents(data,bases,voltages,fold=fold)


def project_currents(data,bases,voltages,*,fold=True):
    paths,split=data[5:]
    nleft=sum(len(p) for p in paths[:4])
    bubble=Bubble(nleft,split)
    ctx=fmpq_mpoly_ctx.get(tuple(f'x{i}' for i in range(14-nleft)))
    right_factors={}
    def parts(polynomial):
        answer=[]
        for l,r,degree,c in split_variables(polynomial,nleft):
            if r not in right_factors:
                right_factors[r]=ctx.from_dict(dict(r))
            answer.append((bubble.register(l),right_factors[r],degree,c))
        return answer

    pieces={key:parts(p) for key,(p,_) in bases.items()}
    projected={}
    excluded={}
    for (i,j,k,h),v in voltages.items():
        if not v:
            continue
        v,content=integer_poly(v,data[1])
        scalar=content*bases[i,j][1]*bases[k,h][1]
        third=parts(v*data[3][i][k])
        for xl,xr,xg,xc in pieces[i,j]:
            for yl,yr,yg,yc in pieces[k,h]:
                square=xr*yr
                finite={}
                singular={}
                for zl,zr,zg,zc in third:
                    value,bad=bubble.project(xl,yl,zl)
                    weight=scalar*xc*yc*zc
                    rational=fmpq(weight.numerator,weight.denominator)
                    for pos,c in bad:
                        singular[pos]=singular.get(pos,0)+zr*(rational*c)
                    if value:
                        finite[zg]=finite.get(zg,0)+zr*(rational*value/SCALE)
                for pos,kernel in singular.items():
                    excluded[pos]=excluded.get(pos,0)+square*kernel
                for z,kernel in finite.items():
                    degree=xg+yg+z
                    projected[degree]=projected.get(degree,0)+square*kernel
    assert not any(excluded.values())
    weights=defaultdict(F)
    for degree,polynomial in projected.items():
        sL=degree+nleft-Q
        for ex,c in polynomial.terms():
            m=tuple(map(int,ex))
            rest=tuple(sum(m[i-nleft]+1 for i in path) for path in paths[4:])
            assert sL+sum(rest)==3*Q
            # Ordinary path beta denominators cancel the tetrahedron moment prefactor.
            moment=F(factorial(sL-1)*factorials(m),factorial(Q-1))
            key=orbital((sL,*rest)) if fold else (sL,*rest)
            weights[key]+=F(int(c.numerator),int(c.denominator))*moment
    return {p:c for p,c in weights.items() if c}


def shape_orbits():
    actions=symmetries(CORE)
    tests=constraints(CORE)
    unseen=set(combinations_with_replacement(range(9),5))
    while unseen:
        marks=min(unseen)
        group={tuple(sorted(action[h] for h in marks)) for action in actions}
        unseen-=group
        if all(sum(marks.count(k) for k in subset)>=n for subset,n in tests):
            yield marks,len(group)


def recorded_totals():
    dataset=json.loads(Path(__file__).with_name('quartic_seven_plane_single_totals.json').read_text())
    assert (dataset['core_automorphisms'],dataset['normal_orders'])==(8,120)
    assert [(tuple(m),s) for m,s in dataset['shapes']]==list(shape_orbits())
    for frame_data in dataset['frames'].values():
        channels=[[F(a),F(b)] for a,b in frame_data['channels']]
        assert len(channels)==36
        values=[sum(s*row[j]/8 for (_,s),row in zip(dataset['shapes'],channels)) for j in range(2)]
        assert values==list(map(F,frame_data['total']))
    if 'single_tensor' in dataset:
        from quartic_seven_plane_double_tensor import SPARSE, recover
        from quartic_scalar_response import interval
        coefficients=quadratic_tensor({k:v['total'] for k,v in dataset['frames'].items()})
        double=recover()
        target=tuple((-a-b,-z) for (a,z),b in zip(coefficients,double))
        assert coefficients==tuple(tuple(map(F,row)) for row in dataset['single_tensor'])
        assert target==tuple(tuple(map(F,row)) for row in dataset['primitive_target'])
        determinant=F(str(fmpq_mat([[str(x) for x in (b,a,z)] for (a,z),b in zip(coefficients,double)]).det()))
        assert determinant.numerator*pow(determinant.denominator,-1,1009)%1009==dataset['channel_rank_residue_1009']
        assert dataset['channel_rank_residue_1009']==253
        bounds=[interval(row,24) for row in target]
        assert bounds[0][1]<0 and all(lo>0 for lo,_ in bounds[1:])
        r2,z2=map(F,dataset['frames']['2']['total'])
        r3,z3=map(F,dataset['frames']['3']['total'])
        r2+=SPARSE[2]
        r3+=SPARSE[3]
        assert interval((r3-5*r2,z3-5*z2),24)[0]>F(3,10**10)
        assert interval((125*r2-9*r3,125*z2-9*z3),24)[0]>F(1,10**9)
    return dataset


def quadratic_tensor(totals):
    from quartic_seven_plane_double import NORMAL
    from quartic_seven_plane_double_tensor import frame
    from quartic_seven_plane_space import as_matrix, gale
    from planar_normal_tensors import tensor_row
    normals=(frame(2),frame(3),NORMAL)
    z=[gale(n) for n in normals]
    matrix=as_matrix([tensor_row(a,(0,1),7) for a in z])
    assert F(str(matrix.det()))==-F(64012032,5)
    observations=[]
    for key,q,a in zip(('2','3','default'),normals,z):
        conversion=F(q[0][0]*q[1][1]-q[0][1]*q[1][0],4*(a[0]-a[1])**2)
        observations.append([conversion*F(v) for v in totals[key]])
    coefficients=matrix.solve(as_matrix(observations))
    assert matrix*coefficients==as_matrix(observations)
    return tuple(tuple(F(str(v)) for v in row) for row in coefficients.tolist())


def sparse_curve(c,second,third):
    """The full scalar on frame(c) from its two independent readings."""
    u=c*(c-1)
    return u*(2*c-1)*((F(third,120)-F(second,24))*u+F(second,4)-F(third,60))


def verify_sparse_curve():
    from quartic_seven_plane_double_tensor import frame, SPARSE
    from quartic_seven_plane_space import gale
    from planar_normal_tensors import tensor_row
    def generators(c):
        return tensor_row(gale(frame(c)),(0,1),7)
    # Every restricted generator has degree <= 7. Eight values verify the polynomial identity.
    for c in range(8):
        assert generators(c)==[sparse_curve(c,a,b) for a,b in zip(generators(2),generators(3))]
    assert sparse_curve(4,SPARSE[2],SPARSE[3])==SPARSE[4]
    return True


def verify_fourth_frame(dataset):
    from quartic_seven_plane_double_tensor import frame
    from quartic_seven_plane_single import shape_sum
    # Rebuild 120 complete fourteen-variable numerators on the opposite-edge shape.
    i=21
    marks=tuple(dataset['shapes'][i][0])
    second=map(F,dataset['frames']['2']['channels'][i])
    third=map(F,dataset['frames']['3']['channels'][i])
    predicted=tuple(sparse_curve(4,a,b) for a,b in zip(second,third))
    assert predicted==(-F(12928598167,5034546316800000),F(2371,1109908800))
    assert shape_sum(marks,normal=frame(4))==predicted
    return predicted


if __name__=='__main__':
    from argparse import ArgumentParser
    from quartic_seven_plane_double_tensor import frame
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('indices',nargs='*',type=int)
    parser.add_argument('--all',action='store_true')
    parser.add_argument('--fourth-check',action='store_true')
    parser.add_argument('--frame',type=int,choices=(2,3))
    args=parser.parse_args()
    shapes=list(shape_orbits())
    assert (len(shapes),sum(n for _,n in shapes))==(36,171)
    dataset=recorded_totals()
    assert verify_sparse_curve()
    if args.fourth_check:
        print('fourth-frame direct orbit',*verify_fourth_frame(dataset),flush=True)
    key='default' if args.frame is None else str(args.frame)
    source={} if args.frame is None else {'normal':frame(args.frame)}
    for i in range(36) if args.all else args.indices or [3]:
        marks,size=shapes[i]
        weights=integrated_weights(marks,**source)
        value=cuts(weights)
        if key in dataset['frames']:
            assert value==tuple(map(F,dataset['frames'][key]['channels'][i]))
        print(i,marks,size,len(weights),*value,flush=True)
        # Release per-shape caches when replaying the entire topology.
        graph_data.cache_clear()
        factorials.cache_clear()
        local_moment.cache_clear()
