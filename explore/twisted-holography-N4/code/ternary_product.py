"""Marked products of the cyclic three single-graviton components, in SU(2).

The trace factor tau=1/2 is restored after extraction. Each triangle strips
i/(32 pi^2). The component of C at theta=0 is the source's ig c.
"""
from collections import defaultdict
from fractions import Fraction as F
from functools import cache
from itertools import permutations, product

SPIN=((0,1,1),(1,0,-1))
COLORS=tuple((*p,(-1)**sum(p[i]>p[j] for i in range(3) for j in range(i+1,3)))
             for p in permutations(range(3)))
UNIT={():1}


def letter(kind,flavor,color,jet=()):
    return kind,flavor,color,tuple(sorted(jet))


def theta(i):
    return letter('0theta',i,0)


def odd(x):
    return x[0] not in ('gamma','top')


def signed_word(word):
    if len({x for x in word if odd(x)})!=sum(map(odd,word)):
        return None,0
    sign=(-1)**sum(odd(x) and odd(y) and x>y for i,x in enumerate(word) for y in word[i+1:])
    return tuple(sorted(word)),sign


def add(out,word,number):
    w,s=signed_word(word)
    if s:
        out[w]+=s*number


def mul(*factors):
    result=UNIT
    for f in factors:
        out=defaultdict(int)
        for w,k in result.items():
            for v,l in f.items():
                add(out,w+v,k*l)
        result={w:k for w,k in out.items() if k}
    return result


def plus(*terms):
    out=defaultdict(int)
    for poly,factor in terms:
        for w,k in poly.items():
            out[w]+=factor*k
    return {w:k for w,k in out.items() if k}


def derivative(poly,d):
    out=defaultdict(int)
    for word,c in poly.items():
        for i,(kind,flavor,color,jet) in enumerate(word):
            if kind!='0theta':
                add(out,word[:i]+(letter(kind,flavor,color,(*jet,d)),)+word[i+1:],c)
    return {w:k for w,k in out.items() if k}


def graviton(i):
    return {(letter('beta',i,c),letter('gamma',(i+1)%3,c)):1 for c in range(3)}


@cache
def superfield(color):
    poly={(letter('ghost',0,color),):1}
    for i in range(3):
        poly[theta(i),letter('gamma',i,color)]=1
        w,s=signed_word((theta((i+1)%3),theta((i+2)%3),letter('beta',i,color)))
        poly[w]=-s
    poly[theta(0),theta(1),theta(2),letter('top',0,color)]=1
    return poly


def variational(poly,color,copy=0):
    """Even delta/delta C. Its key records physical derivatives on the contracted letter."""
    out=defaultdict(lambda:defaultdict(int))
    for w,c in poly.items():
        for i,(kind,f,k,jet) in enumerate(w):
            if k!=color or kind=='0theta':
                continue
            if kind=='beta':
                t,s=(theta(f),),-1
            elif kind=='gamma':
                t,s=(theta((f+1)%3),theta((f+2)%3)),-1
            elif kind=='top':
                t,s=(),-1
            else:
                t,s=tuple(map(theta,range(3))),1
            t=tuple(theta(h[1]+3*copy) for h in t)
            add(out[jet],w[:i]+t+w[i+1:],s*c)
    return {j:{w:k for w,k in row.items() if k} for j,row in out.items()}


def integrate_theta(polynomial):
    tt=tuple(map(theta,range(3)))
    return {w[3:]:v for w,v in polynomial.items() if w[:3]==tt}


def triangle(A,B):
    """Equation (2-34) with at most one derivative on the pair of contracted inputs.

    A derivative on the left spoke transfers with coefficients (r1,-r2)
    to the new C and the left spectators. For a right spoke they are (1-r1,r2).
    """
    av=[variational(A,y) for y in range(3)]
    bv=[variational(B,z) for z in range(3)]
    total=defaultdict(F)
    for x,y,z,col in COLORS:
        for jl,aa in av[y].items():
            for jr,bb in bv[z].items():
                if len(jl)+len(jr)>1:
                    raise ValueError('This triangle handles at most one displaced derivative.')
                for a,b,eps in SPIN:
                    field=derivative(superfield(x),a)
                    left=derivative(aa,b)
                    if jl or jr:
                        d=(jl or jr)[0]
                        terms=[(derivative(field,d),left,F(1 if jl else 2,6)),
                               (field,derivative(left,d),F(-1 if jl else 1,6))]
                    else:
                        terms=[(field,left,F(1,2))]
                    for f,l,moment in terms:
                        for word,c in integrate_theta(mul(f,l,bb)).items():
                            total[word]+=col*eps*moment*c
    return {w:k for w,k in total.items() if k}


def value(poly,word):
    w,s=signed_word(word)
    return s*poly.get(w,0)


def pair_control():
    # The published scalar-fermion example has two gamma^1 factors at the first mark.
    A={(letter('gamma',0,c),)*2:1 for c in range(3)}
    example=triangle(A,graviton(0))
    w=(letter('ghost',0,0,[0]),letter('gamma',0,1,[1]),letter('gamma',1,2))
    # Multiplying by i/(32 pi^2) and ig converts -1 into +g/(32 pi^2).
    assert value(example,w)==-1


def first_order():
    g1,g2,g3=map(graviton,range(3))
    b12=triangle(g1,g2)
    left=plus((mul(b12,g3),1),(triangle(mul(g1,g2),g3),1))
    right=plus((mul(g1,triangle(g2,g3)),1),(triangle(g1,mul(g2,g3)),1))
    return left,right


def two_ghost_fan(A,B):
    """Only the component with two Hessians of the ghost and two undifferentiated spectators.

    Both interactions must emit a ghost. The asymmetric graphs have one
    surviving interaction letter. Terms differentiating the left spectators
    have at most one Hessian, so only r1*r3 survives in the direct fan.
    """
    av={(y,v):variational(variational(A,v,1)[()],y).get((),{}) for y,v in product(range(3),repeat=2)}
    bv={(z,w):variational(variational(B,w,1)[()],z).get((),{}) for z,w in product(range(3),repeat=2)}
    thetas=tuple(map(theta,range(6)))
    out=defaultdict(F)
    for x,y,z,col in COLORS:
        for u,v,w,other in COLORS:
            rest=mul(av[y,v],bv[z,w])
            if not rest:
                continue
            for a,b,e in SPIN:
                for d,h,f in SPIN:
                    ghosts=(letter('ghost',0,x,[a,h]),letter('ghost',0,u,[d,b]))
                    for word,coef in mul({ghosts:1},rest).items():
                        if word[:6]==thetas:
                            # Right Berezin integral moves past three left odd coordinates.
                            add(out,word[6:],-F(col*other*e*f*coef,48))
    return {w:k for w,k in out.items() if k}


def interaction_field(x,copy):
    return {tuple(theta(v[1]+3*copy) if v[0]=='0theta' else v for v in w):c
            for w,c in superfield(x).items()}


def two_loop(A,B):
    """The direct fan and both asymmetric graphs of (2-37), with zero input jets.

    Units are (32 pi^2)^(-2); the two internal ghost components still denote ig c.
    Indices a,b are lower and raised spin choices are summed using SPIN.
    """
    def pairs(poly,copy):
        return {(y,v):variational(variational(poly,v,copy)[()],y).get((),{})
                for y,v in product(range(3),repeat=2)}
    av,bv=pairs(A,1),pairs(B,1)
    aa,bb=pairs(A,0),pairs(B,0)
    a1,b1=[variational(A,k)[()] for k in range(3)],[variational(B,k)[()] for k in range(3)]
    thetas=tuple(map(theta,range(6)))
    total=defaultdict(F)
    for x,y,z,col in COLORS:
        for u,v,w,other in COLORS:
            for a,b,e in SPIN:
                for d,h,f in SPIN:
                    measure=F(col*other*e*f,12)
                    for first,second in product(range(2),repeat=2):
                        field=derivative(interaction_field(x,0),a)
                        otherfield=derivative(interaction_field(u,1),d)
                        left=av[y,v]
                        if first:
                            otherfield=derivative(otherfield,b)
                        else:
                            left=derivative(left,b)
                        if second:
                            field=derivative(field,h)
                        else:
                            left=derivative(left,h)
                        raw=mul(field,otherfield,left,bv[z,w])
                        for word,c in raw.items():
                            if word[:6]==thetas:
                                total[word[6:]]-=measure*c/2**(first+second)
                    if z!=u:
                        continue
                    output=derivative(derivative(superfield(x),a),d)
                    leftA=derivative(derivative(aa[y,v],b),h)
                    leftB=derivative(derivative(a1[v],b),h)
                    raw=plus((mul(output,leftA,b1[w]),1),(mul(output,leftB,bb[y,w]),-1))
                    for word,c in integrate_theta(raw).items():
                        total[word]+=F(col*other*e*f*c,24)
    return {w:c for w,c in total.items() if c}


def measure(shifts):
    from math import factorial
    powers={(0,0):1}
    for terms in shifts:
        new=defaultdict(int)
        for (a,b),c in powers.items():
            for h,k,w in terms:
                new[a+h,b+k]+=w*c
        powers=new
    return sum(F(c*factorial(a)*factorial(b),factorial(a+b+2)) for (a,b),c in powers.items())


def first_fan_displaced(poly):
    """One-loop differential (2-16), with physical jet derivatives on both spokes."""
    from functools import lru_cache
    av=[variational(poly,k) for k in range(3)]
    one_r=(((0,0,1),(1,0,-1)),((0,1,1),))
    two_r=(((1,0,1),),((0,0,1),(0,1,-1)))
    @lru_cache(None)
    def outputs(d,e,j1,j2):
        return mul(with_derivatives(superfield(d),j1),with_derivatives(superfield(e),j2))
    out=defaultdict(F)
    for y,h,d,s in COLORS:
        for z,r,e,t in COLORS:
            if h!=r:
                continue
            for jr,first in av[z].items():
                for jl,dd in variational(first,y).items():
                    if not dd:continue
                    shifts=(one_r,)*len(jl)+(two_r,)*len(jr)
                    for choices in product(range(2),repeat=len(shifts)):
                        weight=measure([v[k] for v,k in zip(shifts,choices)])
                        j1=tuple(j for j,k in zip(jl+jr,choices) if k==0)
                        j2=tuple(j for j,k in zip(jl+jr,choices) if k==1)
                        for a,b,spin in SPIN:
                            raw=mul(outputs(d,e,(a,)+j1,(b,)+j2),dd)
                            for w,c in integrate_theta(raw).items():
                                out[w]+=s*t*spin*weight*c
    return {w:c for w,c in out.items() if c}


def with_derivatives(poly,jet):
    for d in jet:
        poly=derivative(poly,d)
    return poly


def second_fan_on_inputs():
    """Equation (2-28), after extracting i/[48(32 pi^2)^2] and writing h=ig c."""
    out=defaultdict(int)
    for i,d,r,si in COLORS:
        for j,s,t,sj in COLORS:
            if s!=r:
                continue
            for k,v,e,sk in COLORS:
                if v!=t:
                    continue
                for v1,v2,v3 in permutations((i,j,k)):
                    gammas=(letter('gamma',1,v1),letter('gamma',2,v2),letter('gamma',0,v3))
                    for a,b,eps1 in SPIN:
                        for y,z,eps2 in SPIN:
                            ghosts=(letter('ghost',0,d,[a,y]),letter('ghost',0,e,[b,z]))
                            add(out,ghosts+gammas,si*sj*sk*eps1*eps2)
    return {w:c for w,c in out.items() if c}


def second_order():
    g1,g2,g3=map(graviton,range(3))
    left=plus((mul(two_loop(g1,g2),g3),1),(two_loop(mul(g1,g2),g3),1),
              (triangle(triangle(g1,g2),g3),-1))
    right=plus((mul(g1,two_loop(g2,g3)),1),(two_loop(g1,mul(g2,g3)),1),
               (triangle(g1,triangle(g2,g3)),-1))
    return left,right


def classical_letter(x):
    kind,f,color,jet=x
    pol=plus(*[(mul(superfield(y),superfield(z)),F(col,2)) for a,y,z,col in COLORS if a==color])
    if kind=='ghost':
        head,sign=(),1
    elif kind=='gamma':
        head,sign=(theta(f),),-1
    elif kind=='top':
        head,sign=tuple(map(theta,range(3))),-1
    else:
        head,sign=signed_word((theta((f+1)%3),theta((f+2)%3)))
        sign=-sign
    value={w[len(head):]:c*sign for w,c in pol.items() if tuple(k for k in w if k[0]=='0theta')==head}
    for d in jet:
        value=derivative(value,d)
    return value


@cache
def classical_items(x):
    return tuple(classical_letter(x).items())


def classical(poly):
    out=defaultdict(int)
    for word,c in poly.items():
        sign=1
        for i,k in enumerate(word):
            for w,coef in classical_items(k):
                add(out,word[:i]+w+word[i+1:],sign*c*coef)
            if odd(k):
                sign=-sign
    return {w:c for w,c in out.items() if c}


def first_fan(poly):
    """Equation (2-16) on undifferentiated input, in units -1/(128 pi^2)."""
    colors=range(3)
    out=defaultdict(int)
    for y,h,d,s1 in COLORS:
        for z,r,e,s2 in COLORS:
            if h!=r:
                continue
            dd=variational(variational(poly,z)[()],y).get((),{})
            for a,b,spin in SPIN:
                terms=integrate_theta(mul(derivative(superfield(d),a),derivative(superfield(e),b),dd))
                for w,c in terms.items():
                    out[w]+=s1*s2*spin*c
    return {w:c for w,c in out.items() if c}


def second_order_channel():
    g1,g2,g3=map(graviton,range(3))
    w=(letter('ghost',0,0,[0,0]),letter('ghost',0,1,[1,1]),
       letter('gamma',2,0),letter('beta',2,1))
    # On the left bracketing the direct fan consumes G3, each asymmetric graph
    # emits only one C, and the inner direct two-loop term retains gamma^1.
    nested=triangle(triangle(g1,g2),g3)
    assert value(nested,w)==F(1,12)
    right_nested=triangle(g1,triangle(g2,g3))
    fan=two_ghost_fan(g1,mul(g2,g3))
    assert value(right_nested,w)==0
    # The printed direct fan has +1/[2(32 pi^2)^2], in contrast to the
    # squared imaginary triangle prefactor. Converting both ig*c factors
    # into physical ghosts gives an additional minus on the direct fan.
    assert value(fan,w)==-F(1,24)
    return w,nested,fan


if __name__=='__main__':
    pair_control()
    left,right=first_order()
    w=(letter('ghost',0,0,[0]),letter('beta',0,1,[1]),
       letter('gamma',2,2),letter('beta',2,2),letter('gamma',0,2))
    # ig from C and i from the triangle supply minus one.
    assert value(left,w)==F(1,2)
    q1=first_fan(mul(*map(graviton,range(3))))
    assert not plus((classical(left),4),(q1,1))
    assert not plus((classical(right),4),(q1,1))
    w,nested,fan=second_order_channel()
    print('Left and right first-order component monomials',len(left),len(right))
    print('Both first-order Ward equations checked on',len(q1),'nonzero source monomials.')
    print('Left and right two-ghost coefficients in g^2/pi^4 units',
          value(nested,w)/(8*32**2),-value(fan,w)/(8*32**2))

    # The pure-ghost cubic term may be dropped from either bracketing without
    # changing closure. At order two the differential has the common factor
    # i/(32 pi^2)^2 when the zero component of C denotes ig c.
    q2=second_fan_on_inputs()
    for o1,o2 in zip((left,right),second_order()):
        assert not plus((classical(o2),1),(first_fan_displaced(o1),-F(1,2)),(q2,F(1,48)))
        assert not first_fan_displaced(o2)
        assert max(sum(k[0]=='beta' for k in w) for w in o1)==2
        assert max(sum(k[0]=='beta' for k in w) for w in o2)==1
        print('Two-loop Ward identity on',len(o2),'product and',len(q2),'source monomials.')
