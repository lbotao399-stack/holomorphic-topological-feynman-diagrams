"""Pull back the two bulk edges in Oh–Zhou Fig. 12 and integrate the printed scalar contour."""
import sympy as sp


def wedge(a,b):
    out={}
    for I,c in a.items():
        for J,d in b.items():
            if set(I)&set(J):
                continue
            sign=(-1)**sum(i>j for i in I for j in J)
            key=tuple(sorted(I+J))
            out[key]=out.get(key,0)+sign*c*d
    return {I:sp.expand(c) for I,c in out.items() if sp.expand(c)!=0}


def basis():
    # Slots of the pulled-back one-forms: dt_2, dzbar_2, dwbar_2, dwbar_1, dwbar_3.
    return [{(i,):1} for i in range(5)]


def difference(a,b):
    return {i:a.get(i,0)-b.get(i,0) for i in a.keys()|b.keys() if a.get(i,0)!=b.get(i,0)}


def propagator(values,differentials):
    t,z,w=values
    T,Z,W=differentials
    pairs=((z,W,T),(-w,Z,T),(t,Z,W))
    out={}
    for scalar,u,v in pairs:
        for ij,c in wedge(u,v).items():
            out[ij]=out.get(ij,0)+scalar*c
    return out


def exterior():
    t,z,w,w1,w3=sp.symbols('t zbar wbar wbar1 wbar3')
    T,Z,W,D1,D3=basis()
    coords=((0,0,w1),(t,z,w),(0,0,w3))
    diffs=(({}, {},D1),(T,Z,W),({}, {},D3))
    edges=[]
    for i,j in ((0,1),(1,2)):
        edges.append(propagator(tuple(x-y for x,y in zip(coords[i],coords[j])),
                                tuple(difference(x,y) for x,y in zip(diffs[i],diffs[j]))))
    mixed=[{(i,):v for (i,j),v in edge.items() if j==last}
           for edge,last in zip(edges,(3,4))]
    assert mixed==[{(1,):t,(0,):-z},{(1,):-t,(0,):z}]
    assert not wedge(*mixed)
    assert not {k:v for k,v in wedge(*edges).items() if 3 in k and 4 in k}


def circle():
    a,s=sp.symbols('a s',positive=True)
    # Set R^2=1 by homogeneity and epsilon^2=a^2/(1-a^2).
    positive=a**2/(1-a**2)
    second=2*(s**2-a**2)**2*(5-3*s**2)/(s**4*(1-a**2)**2)
    primitive=sp.integrate(sp.expand(2*(s**2-a**2)**2*(5-3*s**2)/s**4),s)/(1-a**2)**2
    assert sp.cancel(primitive.diff(s)-second)==0
    ratio=sp.factor(1-3*positive*(primitive.subs(s,1)-primitive.subs(s,a))/2)
    assert sp.factor(ratio-(1-a)**2*(5*a+1)/(1+a)**3)==0
    assert sp.factor(ratio.diff(a)-24*a*(a-1)/(a+1)**4)==0
    return ratio


if __name__=='__main__':
    exterior()
    print('wedge at the common bulk vertex = 0; printed circle scalar / first term =',circle())
