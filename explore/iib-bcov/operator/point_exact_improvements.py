"""Exact changes of the rank-five Euler representative on the two integer probes."""
from functools import cache
from itertools import combinations, permutations


@cache
def determinant(rows):
    if len(rows)==1:
        return rows[0][0]
    return sum((-1)**j*x*determinant(tuple(row[:j]+row[j+1:] for row in rows[1:]))
               for j,x in enumerate(rows[0]) if x)


def improvements(plane):
    frame=tuple(tuple(int(i==j) for i in range(5)) for j in range(5))
    wave=tuple(tuple(p)+(0,0,0) for p in plane)
    scalar=(*wave,tuple(-sum(p[i] for p in wave)-int(i>=2) for i in range(5)))
    vector_momentum=frame[2]
    weights=((1,2,4),(1,3,3))
    answers=[0,0]
    # Only fusion of the four-vector with a scalar leaves the required ranks
    # (one trivector and one vector). The other seven scalars stay at the mark.
    for i,p in enumerate(scalar):
        merged=tuple(a+b for a,b in zip(frame[4],p))
        remaining=scalar[:i]+scalar[i+1:]
        for q,s,t,u in permutations(remaining,4):
            # Transverse four-vector e1^e2^e3^e4, vector e5.
            first=merged[4]*determinant(tuple(x[:4] for x in (merged,vector_momentum,q,s)))
            if not first:
                continue
            ds=tuple(determinant((merged,vector_momentum,q,x,y))
                     for x,y in ((s,t),(s,u),(t,u)))
            if not all(ds):
                continue
            for j,(a,b,c) in enumerate(weights):
                answers[j]+=first*ds[0]**a*ds[1]**b*ds[2]**c
    # Independent expansion of every five-row determinant in the two-plane.
    def bracket(x,y):
        return x[0]*y[1]-x[1]*y[0]
    reduced=[sum(bracket(p,s)*bracket(s,t)**a*bracket(s,u)**b*bracket(t,u)**c
                 for p,s,t,u in permutations(plane,4)) for a,b,c in weights]
    assert answers==reduced
    return tuple(answers)


def scaling_check(plane):
    frame=tuple(tuple(int(i==j) for i in range(5)) for j in range(5))
    vector=frame[2]
    for fused,p in enumerate(plane):
        spectator=tuple(tuple(q)+(0,0,0) for j,q in enumerate(plane) if j!=fused)
        determinants=[]
        tensors=[]
        for strength in (0,1,2):
            merged=(strength*p[0],strength*p[1],0,0,1)
            balance=tuple(-sum(q[j] for q in spectator)-merged[j]-vector[j]-frame[3][j]
                          for j in range(5))
            rows=(merged,vector,balance,*spectator)
            determinants.append(tuple(determinant(t) for t in combinations(rows,5)))
            tensors.append(tuple(determinant(tuple(q[:4] for q in (merged,*t)))*v[4]
                                 for t in combinations(rows,3) for v in rows))
        assert determinants[0]==determinants[1]==determinants[2]
        assert not any(tensors[0]) and any(tensors[1])
        assert tuple(2*x for x in tensors[1])==tensors[2]


def check():
    rows=(improvements(((1,0),(1,0),(2,0),(0,2),(0,2),(1,2),(-2,-1))),
          improvements(((2,0),(2,0),(0,-2),(0,1),(1,1),(2,2),(1,-1))))
    assert rows==((191032,-236568),(-1207384,-30648))
    assert determinant(rows)==-291483166848
    scaling_check(((1,0),(1,0),(2,0),(0,2),(0,2),(1,2),(-2,-1)))
    print('Exact representative changes:',rows,'determinant',determinant(rows))


if __name__=='__main__':
    check()
