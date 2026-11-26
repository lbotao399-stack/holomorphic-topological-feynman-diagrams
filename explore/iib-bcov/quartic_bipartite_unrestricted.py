"""Tropical support of the original five-component heat cofactor on subdivided K3,3."""
from itertools import product
from math import prod
from random import Random
from flint import nmod_mat


PRIME=1009
OPPOSITE=(4,5,7,8)


def interpolate(samples):
    # Newton divided differences, followed by exact conversion to the monomial basis.
    xs=[x for x,y in samples]; coefficients=[y for x,y in samples]
    for level in range(1,len(xs)):
        for i in range(len(xs)-1,level-1,-1):
            coefficients[i]=(coefficients[i]-coefficients[i-1])*pow(xs[i]-xs[i-level],-1,PRIME)%PRIME
    out=[coefficients[-1]]
    for x,c in zip(reversed(xs[:-1]),reversed(coefficients[:-1])):
        out=[(c-x*out[0])%PRIME]+[(a-x*b)%PRIME for a,b in zip(out,out[1:])]+[out[-1]]
    return out


class Subdivision:
    def __init__(self,counts,seed):
        if len(counts)!=9 or sum(counts)!=12:
            raise ValueError("A bare four-loop twelve-wave placement is required.")
        self.counts=counts; random=Random(seed)
        frame=((1,0),(0,1),(-1,-1))
        directions=[tuple(u[r]*v[s] for r,s in product(range(2),repeat=2)) for u in frame for v in frame]
        self.parent=[];self.fractions=[];self.edges=[];new=6
        for i,m in enumerate(counts):
            vertices=[i//3,*range(new,new+m),3+i%3];new+=m
            lengths=[random.randrange(1,16) for j in range(m+1)]
            for e,t in zip(zip(vertices,vertices[1:]),lengths):
                self.edges.append(e);self.parent.append(i)
                self.fractions.append(t*pow(sum(lengths),-1,PRIME)%PRIME)
        self.cycles=[[directions[i][j] for i in self.parent] for j in range(4)]
        self.b=nmod_mat([[int(i==v)-int(j==v) for v in range(1,new)] for i,j in self.edges],PRIME)
        waves=[[random.randrange(-7,8) for k in range(5)] for i in range(11)]
        waves[0][4]=waves[1][4]=0
        waves.append([-sum(col) for col in zip(*waves)])
        self.source=nmod_mat([[0]*5 for _ in range(5)]+waves,PRIME)
        self.base=[random.randrange(2,70) for e in range(9)]
        self.core_cycles=nmod_mat(list(zip(*directions)),PRIME)

    def kirchhoff(self,t):
        c=self.core_cycles
        return int((c*nmod_mat([[t[i] if i==j else 0 for j in range(9)] for i in range(9)],PRIME)*c.transpose()).det())

    def sample(self,t):
        moment=self.kirchhoff(t)
        if not moment:return None
        segments=[t[e]*s%PRIME for e,s in zip(self.parent,self.fractions)]
        conductance=nmod_mat([[pow(segments[i],-1,PRIME) if i==j else 0 for j in range(21)] for i in range(21)],PRIME)
        flows=conductance*self.b
        potential=(self.b.transpose()*flows).solve(self.source)
        curr=flows*potential
        # Omit segment zero and divide that 20-by-20 minor by its time.
        rows=[[int(curr[j,v])*cy[j] for j in range(1,21)] for v in range(5) for cy in self.cycles]
        cofactor=int(nmod_mat(rows,PRIME).det())*pow(segments[0],-1,PRIME)%PRIME
        voltage=int(potential[5,4]-potential[6,4])
        jacobian=prod(pow(x,n,PRIME) for x,n in zip(t,self.counts))%PRIME
        numerator=jacobian*cofactor**2*voltage*pow(moment,9,PRIME)%PRIME
        return numerator,cofactor*pow(moment,4,PRIME)%PRIME

    def weight_degree(self,lower=-47,upper=47):
        # At fixed path fractions U^4*c has degree 15 and U*voltage has degree 5.
        # Twelve path-volume factors give the safe two-sided Laurent bound 47.
        data=[];held=[]
        for s in range(1,PRIME):
            t=[base*pow(s,1 if i==0 else -1 if i in OPPOSITE else 0,PRIME)%PRIME
               for i,base in enumerate(self.base)]
            sample=self.sample(t)
            if sample is None:continue
            pair=(s,sample[0]*pow(s,-lower,PRIME)%PRIME)
            if len(data)<upper-lower+1:data.append(pair)
            else:held.append(pair)
            if len(held)==3:break
        poly=interpolate(data)
        for s,val in held:assert sum(c*pow(s,j,PRIME) for j,c in enumerate(poly))%PRIME==val
        return max((lower+j for j,c in enumerate(poly) if c),default=None)

    def kirchhoff_divisor(self):
        # Vary one complete core time at fixed fractions. The current determinant
        # has degree at most 15, enough to test divisibility without a degree refinement.
        data=[];full=[]
        for s in range(1,PRIME):
            t=[s,*self.base[1:]];sample=self.sample(t)
            if sample:
                if len(data)<16:data.append((s,sample[1]))
                full.append((s,sample[0]))
            if len(full)==48:break
        poly=interpolate(data)
        full_poly=interpolate(full)
        zero=self.kirchhoff([0,*self.base[1:]]);one=self.kirchhoff([1,*self.base[1:]])
        if one==zero:raise ValueError("Kirchhoff root at infinity on this fixture")
        root=-zero*pow(one-zero,-1,PRIME)%PRIME
        return (root,sum(c*pow(root,i,PRIME) for i,c in enumerate(poly))%PRIME,
                sum(c*pow(root,i,PRIME) for i,c in enumerate(full_poly))%PRIME)


if __name__=='__main__':
    fixtures=(((4,1,2,1,1,0,2,0,1),711,-1,(376,449)),
              ((3,2,2,1,1,0,2,0,1),723,0,(360,836)),
              ((2,1,2,1,1,2,1,1,1),709,-5,(925,827)))
    for counts,seed,expected,residue in fixtures:
        graph=Subdivision(counts,seed)
        degree=graph.weight_degree()
        test=graph.kirchhoff_divisor()
        assert degree==expected and test[:2]==residue and test[2]
        print(counts,seed,'tropical numerator degree',degree,'Kirchhoff root, current, full numerator',test)
