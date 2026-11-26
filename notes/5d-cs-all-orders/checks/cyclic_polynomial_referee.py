"""Independent polygon-polynomial differential, with exact finite-field ranks.

Coordinates are vertex differences relative to the first vertex, not trace words.
No repository implementation is imported. The calculations test homogeneous
momentum degrees 2, 3, 4, 5 and cyclic cohomological degrees 2 and 3 over
the prime fields of orders 1000003 and 1000033. They are a bounded independent
check, not a replacement for the characteristic-zero all-degree proof.
"""
from functools import lru_cache

@lru_cache(None)
def monomials(n, degree):
    if n == 0:
        return ((),) if degree == 0 else ()
    if n == 1:
        return ((degree,),)
    return tuple((a,) + b for a in range(degree + 1)
                 for b in monomials(n - 1, degree - a))

def multiply(left, right):
    out = {}
    for a, x in left.items():
        for b, y in right.items():
            c = tuple(u+v for u,v in zip(a,b))
            out[c] = out.get(c,0)+x*y
    return {k:v for k,v in out.items() if v}

def substitute(exponents, images, n):
    out = {(0,)*n:1}
    for power, entries in zip(exponents, images):
        linear = {}
        for index, coeff in entries:
            mon = [0]*n
            mon[index] = 1
            linear[tuple(mon)] = coeff
        for _ in range(power):
            out = multiply(out, linear)
    return out

def maps(j,m):
    source = monomials(2*j,m)
    dest = monomials(2*(j+1),m)
    sindex = {v:i for i,v in enumerate(source)}
    dindex = {v:i for i,v in enumerate(dest)}
    cycle = [[0]*len(source) for _ in source]
    delta = [[0]*len(source) for _ in dest]
    # Rotated polygon is (x_j,0,x_1,...,x_{j-1}); subtract x_j.
    rotate = []
    for vertex in range(j):
        for coordinate in range(2):
            terms = [(2*(j-1)+coordinate,-1)]
            if vertex:
                terms.append((2*(vertex-1)+coordinate,1))
            rotate.append(terms)
    for column, mon in enumerate(source):
        for term, value in substitute(mon,rotate,2*j).items():
            cycle[sindex[term]][column] += (-1)**j*value
        cycle[column][column] -= 1
        # Face deleting x_0 translates the remaining polygon by -x_1.
        first = [[(2*(v+1)+c,1),(c,-1)]
                 for v in range(j) for c in range(2)]
        for term,value in substitute(mon,first,2*(j+1)).items():
            delta[dindex[term]][column] += value
        for omitted in range(1,j+2):
            retain = [v for v in range(1,j+2) if v!=omitted]
            images = [[(2*(v-1)+c,1)] for v in retain for c in range(2)]
            for term,value in substitute(mon,images,2*(j+1)).items():
                delta[dindex[term]][column] += (-1)**omitted*value
    return cycle,delta

def rank(matrix,p):
    if not matrix: return 0
    rows = [[x%p for x in row] for row in matrix]
    pivot = 0
    for c in range(len(rows[0])):
        chosen = next((r for r in range(pivot,len(rows)) if rows[r][c]),None)
        if chosen is None: continue
        rows[pivot],rows[chosen] = rows[chosen],rows[pivot]
        inv = pow(rows[pivot][c],-1,p)
        rows[pivot] = [x*inv%p for x in rows[pivot]]
        for r in range(pivot+1,len(rows)):
            value=rows[r][c]
            if value:
                rows[r]=[(x-value*y)%p for x,y in zip(rows[r],rows[pivot])]
        pivot+=1
        if pivot==len(rows): break
    return pivot

def run(m,p):
    rows=[]
    for j in range(1,4):
        a,b=maps(j,m)
        ra=rank(a,p)
        rb=rank(a+b,p)-ra
        rows.append((j,len(a)-ra,rb))
    # CE arities 3 and 4 correspond to cyclic j=2 and j=3.
    h2=rows[1][1]-rows[1][2]-rows[0][2]
    h3=rows[2][1]-rows[2][2]-rows[1][2]
    print('degree',m,'prime',p,'(j,cyclic dimension,d rank)',rows,
          'HC^2,HC^3',h2,h3,flush=True)
    assert (h2,h3)==((1,0) if m==2 else (0,0))

if __name__=='__main__':
    print('Independent polygon-polynomial check; exact ranks over finite fields.')
    print('Momentum degrees: 2, 3, 4, 5. Cohomological degrees: 2 and 3.')
    print('Primes: 1000003 and 1000033. Degree 2 is a nonzero positive control.')
    print('This bounded check does not replace the all-degree complex-number proof.')
    for m in (2,3,4,5):
        for p in (1000003,1000033): run(m,p)
