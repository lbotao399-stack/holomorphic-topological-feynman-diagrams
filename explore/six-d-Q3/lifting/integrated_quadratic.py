"""Classify invariant quadratic Fourier symbols for the two odd source coordinates."""
from sympy import Matrix


def product(i,j):
    if i&j:
        return 0,0
    return i|j, -1 if (i&2 and j&1) else 1


def constraints(order,parity):
    entries=[(i,j) for i in range(4) for j in range(4)
             if (i.bit_count()+j.bit_count())%2==parity]
    columns=[]
    for target,source in entries:
        def operator(i):
            return (target,1) if i==source else (0,0)
        column=[]
        for i in range(4):
            for j in range(4):
                ij,ij_sign = product(i,j)
                o,mo = operator(ij)
                u,mu = operator(i)
                v,mv = operator(j)
                ouj,uj = product(u,j)
                oiv,iv = product(i,v)
                sign = (-1)**(parity*i.bit_count())
                # Integrated skew-adjointness at p+q=0, picking out theta_1 theta_2.
                column.append(mu*uj*(ouj==3)+sign*(-1)**order*mv*iv*(oiv==3))
                for a,b in ((1,0),(0,1),(1,1)):
                    result=[0]*4
                    result[o] += (a+b)**order*ij_sign*mo
                    result[ouj] -= a**order*mu*uj
                    result[oiv] -= sign*b**order*mv*iv
                    column.extend(result)
        columns.append(column)
    return Matrix(columns).T


if __name__=='__main__':
    for order in range(7):
        dimensions=tuple(8-constraints(order,parity).rank() for parity in (0,1))
        assert dimensions==((3,2) if order==0 else (2,2) if order==1 else (0,0))
        print(order,dimensions)
