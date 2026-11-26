"""Invariant dimension lowering for the determinant on traceless symmetric 4 x 4 matrices."""
from argparse import ArgumentParser
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations, permutations
from math import factorial, prod
from pathlib import Path
import json

from flint import fmpq, fmpq_mat, fmpz_mpoly_ctx
import sympy as sp


def minor(matrix,indices):
    return sum((-1)**sum(a>b for k,a in enumerate(p) for b in p[k+1:])
               *prod(matrix[i][j] for i,j in zip(indices,p)) for p in permutations(indices))


def elementary(matrix):
    return [sum(minor(matrix,ids) for ids in combinations(range(4),n)) for n in range(1,5)]


def partitions(degree):
    for c in range(degree//4+1):
        for b in range((degree-4*c)//3+1):
            a,rem=divmod(degree-4*c-3*b,2)
            if not rem: yield a,b,c


def symmetric_coefficients(polynomial):
    ctx=polynomial.context()
    x=ctx.gens()
    elementary_x=[sum(prod(items) for items in combinations(x,k)) for k in range(1,5)]
    result={}
    while polynomial:
        lead=tuple(map(int,polynomial.monomial(0)))
        exponent=tuple(lead[i]-lead[i+1] for i in range(3))+(lead[3],)
        assert min(exponent)>=0
        c=int(polynomial.coefficient(0))
        result[exponent]=c
        polynomial-=c*prod(p**v for p,v in zip(elementary_x,exponent))
    return result


class SymmetricDeterminant:
    def __init__(self,strength):
        self.strength=strength
        self.slots=list(combinations(range(4),2))+[(i,i) for i in range(4)]
        self.ctx=fmpz_mpoly_ctx.get(tuple(f'z{i}{j}' for i,j in self.slots))
        self.diagonal=fmpz_mpoly_ctx.get(tuple(f'x{i}' for i in range(4)))
        z=[[None]*4 for _ in range(4)]
        partial=[[None]*4 for _ in range(4)]
        for k,(i,j) in enumerate(self.slots):
            z[i][j]=z[j][i]=self.ctx.gen(k)
        trace=sum(z[i][i] for i in range(4))
        for i,j in self.slots:
            # Four times the actual projected derivative. Its symbol of order
            # r is divided by 4**r before solving the differential identity.
            partial[i][j]=partial[j][i]=4*z[i][j]-trace if i==j else 2*z[i][j]
        self.symbols=elementary(partial)[1:]
        self.parts=elementary(z)
        self.taylor=[self.ctx.constant(1)]
        self.images={}

    def image(self,symbol):
        if symbol in self.images:return self.images[symbol]
        degree=sum(k*p for k,p in zip((2,3,4),symbol))
        while len(self.taylor)<=degree:
            n=len(self.taylor)
            numerator=sum(((self.strength+1)*j-n)*self.parts[j-1]*self.taylor[n-j]
                          for j in range(1,min(n,4)+1))
            quotient=numerator//n
            assert n*quotient==numerator
            self.taylor.append(quotient)
        operator=prod(v**i for v,i in zip(self.symbols,symbol))
        homogeneous=self.taylor[degree]
        result=defaultdict(int)
        for ex,c in operator.terms():
            ex=tuple(map(int,ex))
            incidence=[0]*4
            for n,(i,j) in zip(ex,self.slots):
                incidence[i]+=n
                incidence[j]+=n
            assert all(not v%2 for v in incidence)
            powers=tuple(v//2 for v in incidence)
            result[powers]+=int(c)*int(homogeneous[ex])*prod(factorial(v) for v in ex)
        coeffs=symmetric_coefficients(self.diagonal.from_dict(result))
        self.images[symbol]={ex:F(c,4**degree) for ex,c in coeffs.items() if not ex[2]}
        return self.images[symbol]


def search(strength,last):
    engine=SymmetricDeterminant(strength)
    labels=[]
    columns=[]
    target={(0,0,1):F(1)}
    for degree in range(4,last+1):
        for multipliers in partitions(degree-4):
            a,b,c=multipliers
            for symbol in partitions(degree):
                col=engine.image(symbol)
                columns.append({(i+b,j+a,l-a-b-c):value for (i,j,k,l),value in col.items()})
                labels.append((multipliers,symbol))
        if not list(partitions(degree-4)):continue
        keys=sorted(set(target).union(*(set(col) for col in columns)))
        matrix=sp.Matrix([[col.get(key,0) for col in [*columns,target]] for key in keys])
        echelon,pivots=matrix.rref()
        print('highest derivative',degree,'symbols',len(columns),'rank',len(pivots),flush=True)
        if len(columns) in pivots:continue
        solution={labels[j]:F(echelon[i,-1]) for i,j in enumerate(pivots) if echelon[i,-1]}
        actual=defaultdict(F)
        for label,col in zip(labels,columns):
            for key,v in col.items():
                actual[key]+=solution.get(label,0)*v
        assert {k:v for k,v in actual.items() if v}==target
        print('traceless symmetric determinant shift',strength,solution,flush=True)
        return solution
    return None


def graph_symbols():
    from quartic_planar_periods import core_data
    from four_loop import K33
    from square_three import G
    cycle=G(6,K33).C
    ctx,u,_,_=core_data(tuple(K33))
    x=ctx.gens()
    cotrees=[tuple(map(int,ex)) for ex,_ in u.terms()]
    for size in range(1,9):
        for edges in combinations(range(9),size):
            rank=min(sum(tree[j] for j in edges) for tree in cotrees)
            assert 4*size>9*rank
    slots=[(i,j) for i in range(4) for j in range(i,4) if (i,j)!=(0,3)]
    incidence=fmpq_mat([[cycle[i][e]*cycle[j][e]*(1+(i!=j)) for i,j in slots] for e in range(9)])
    inverse=incidence.inv()
    normal=((0,0,0,-1),(0,0,1,0),(0,1,0,0),(-1,0,0,0))
    assert fmpq_mat(normal).det()==1
    matrix=[[sum(cycle[i][e]*cycle[j][e]*x[e] for e in range(9)) for j in range(4)] for i in range(4)]
    derivative=[[ctx.constant(0) for _ in range(4)] for _ in range(4)]
    for k,(i,j) in enumerate(slots):
        derivative[i][j]=derivative[j][i]=sum(inverse[k,e]*x[e] for e in range(9))
    def act_normal(matrix):
        return [[sum(normal[i][v]*matrix[v][j] for v in range(4)) for j in range(4)] for i in range(4)]
    l=act_normal(matrix)
    assert not sum(l[i][i] for i in range(4)) and minor(l,range(4))==u
    diff=act_normal(derivative)
    trace=sum(diff[i][i] for i in range(4))/4
    for i in range(4):diff[i][i]-=trace
    invariants=elementary(l)[1:]
    symbols=elementary(diff)[1:]
    return u,invariants,symbols


def graph_check(strength,solution):
    u,invariants,symbols=graph_symbols()
    ctx=u.context()
    x=ctx.gens()
    terms=[(prod(p**i for p,i in zip(invariants,a)),prod(p**i for p,i in zip(symbols,b)),fmpq(str(v)))
           for (a,b),v in solution.items()]
    for times in ([1]*9,[2,3,5,7,11,13,17,19,23]):
        value=u(*times)
        shifted=u.compose(*(x[i]+t for i,t in enumerate(times)))
        parts=[{} for _ in range(5)]
        for ex,c in shifted.terms():
            powers=tuple(map(int,ex))
            parts[sum(powers)][powers]=c/value
        parts=[ctx.from_dict(row) for row in parts]
        normalized=[ctx.constant(1)]
        for n in range(1,9):
            normalized.append(sum(((strength+1)*j-n)*parts[j]*normalized[n-j]
                                  for j in range(1,min(4,n)+1))/n)
        total=0
        for multiplier,operator,weight in terms:
            jet=normalized[int(operator.total_degree())]
            coeff=sum(c*jet[tuple(map(int,r))]*prod(factorial(int(k)) for k in r) for r,c in operator.terms())
            total+=multiplier(*times)*weight*coeff
        assert total*value==1
    print('Independent nine-edge Cayley-Taylor check in the original positive graph metric.',flush=True)


def universal():
    # At order n every normalized Taylor coefficient is a polynomial in the
    # strength of degree at most n. Exact interpolation reconstructs that
    # polynomial before solving the operator identity over Q(s).
    s=sp.Symbol('s')
    fixed=search(-11,8)
    labels=list(fixed)
    nodes=range(9)
    images=[]
    for n in nodes:
        engine=SymmetricDeterminant(n)
        rows=[]
        for (a,b,c),symbol in labels:
            rows.append({(i+b,j+a,l-a-b-c):v for (i,j,k,l),v in engine.image(symbol).items()})
        images.append(rows)
    keys=sorted({(0,0,1)}.union(*(set(col) for rows in images for col in rows)))
    def interpolate(values):
        result=0
        differences=list(map(sp.Rational,values))
        polynomial=1
        for n in nodes:
            result+=differences[0]*polynomial
            differences=[b-a for a,b in zip(differences,differences[1:])]
            polynomial*=s-n
            polynomial/=n+1
        return sp.Poly(result,s,domain=sp.QQ).as_expr()
    matrix=sp.Matrix([[interpolate([rows[j].get(key,0) for rows in images]) for j in range(len(labels))]
                      for key in keys])
    target=sp.Matrix([int(key==(0,0,1)) for key in keys])
    sample=matrix.subs(s,-11)
    selected=sample.T.rref()[1]
    assert len(selected)==len(labels)
    square=matrix.extract(selected,range(len(labels)))
    rhs=target.extract(selected,[0])
    system,vector=square.to_DM().unify(rhs.to_DM())
    solved,den=system.solve_den(vector)
    coefficients=[sp.factor(x/system.domain.to_sympy(den)) for x in solved.to_Matrix()]
    common=sp.factor(sp.lcm([sp.denom(c) for c in coefficients]))
    assert all(common.subs(s,n) for n in range(-11,-2))
    top=[sp.factor(c*common) for c in coefficients]
    assert all(sp.expand(v)==0 for v in (matrix*sp.Matrix(top)-common*target))
    print('polynomial strength denominator',common,flush=True)
    for label,c in zip(labels,top):
        print(label,c,flush=True)
    return s,dict(zip(labels,coefficients)),common


def coefficient_table(variable,form,denominator):
    def coefficients(expr):
        return list(map(int,reversed(sp.Poly(sp.expand(expr),variable,domain=sp.ZZ).all_coeffs())))
    return {'powers':'ascending in s', 'common_denominator':coefficients(denominator),
            'terms':[{'matrix':list(a),'derivatives':list(b),
                      'numerator':coefficients(sp.cancel(value*denominator))}
                     for (a,b),value in form.items()]}


if __name__=='__main__':
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('--strength',type=int,default=-11)
    parser.add_argument('--through',type=int,default=12)
    parser.add_argument('--universal',action='store_true')
    parser.add_argument('--save',action='store_true')
    args=parser.parse_args()
    if args.universal:
        variable,form,den=universal()
        table=json.dumps(coefficient_table(variable,form,den),indent=2)+'\n'
        path=Path(__file__).with_suffix('.json')
        if args.save:path.write_text(table)
        else:assert path.read_text()==table
        solution={k:F(v.subs(variable,args.strength)) for k,v in form.items()}
    else:
        solution=search(args.strength,args.through)
    if solution:
        graph_check(args.strength,solution)
