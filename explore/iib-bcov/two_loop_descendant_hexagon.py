"""The one-loop insertion in the nine-wave genus-two BCOV Ward equation on C^5."""
from fractions import Fraction
from itertools import combinations, permutations
from math import factorial, lcm

from flint import fmpz_mat
from marked_theta import green_moments


def add_rows(rows):
    return tuple(map(sum,zip(*rows)))


def det(rows):
    return int(fmpz_mat(rows).det())


def moments(coeff):
    """Integral of G(s_i,s_j)*(sum_k coeff[k]*s_k)^2 on the five-cube."""
    if len(coeff)!=5 or sum(coeff):
        raise ValueError('The five determinant coefficients must sum to zero.')
    total = sum(x*x for x in coeff)
    return [[(Fraction(total,72)-Fraction(a*a,180) if i==j else
              Fraction(total,144)-Fraction(a*a+b*b-a*b,360))
             for j,b in enumerate(coeff)] for i,a in enumerate(coeff)]


def group_potentials(group, single):
    volume = det(single)
    if not volume:
        return [(Fraction(0),)*5]*5
    coeff = tuple(det((*group,p)) for p in single)
    green = moments(coeff)
    # Positive six-field Hessian coefficient 1/24; rooted circle has reflection factor 1/2.
    return [tuple(Fraction(volume**2,48)*sum(m*x for m,x in zip(row,coordinate))
                  for coordinate in zip(*single)) for row in green]


def potentials(momenta, descendant):
    """All slot potentials relative to the marked genus-one vertex, in Li-Wang flux sign."""
    if len(momenta)!=9 or any(len(q)!=5 for q in momenta) or any(add_rows(momenta)):
        raise ValueError('Nine conserved fivefold momenta are required.')
    if not 0<=descendant<9:
        raise ValueError('Expected a marked wave label 0 through 8.')
    den = lcm(*(Fraction(x).denominator for p in momenta for x in p))
    integer = tuple(tuple(int(den*Fraction(x)) for x in p) for p in momenta)
    labels = [i for i in range(9) if i!=descendant]
    out = [[Fraction(0)]*5 for _ in momenta]
    for others in combinations(labels,3):
        group = (descendant,*others)
        active = tuple(i for i in labels if i not in others)
        rows = group_potentials(tuple(integer[i] for i in group),tuple(integer[i] for i in active))
        for i,row in zip(active,rows):
            for j,x in enumerate(row):
                out[i][j] += x
    return [tuple(x/den**21 for x in row) for row in out]


def compositions(total,count):
    if count==1:
        yield (total,)
        return
    for i in range(total+1):
        for row in compositions(total-i,count-1):
            yield (i,*row)


def simplex_test(group,single):
    """Literal simplex cubature at five ordered waves, including both internal momenta."""
    ans = [[Fraction(0)]*5 for _ in range(5)]
    # Grundmann-Moller rule of degree five on the projective five-simplex.
    for j in range(3):
        den = 10-2*j
        weight = Fraction((-1)**j*den**5,16*factorial(j)*factorial(10-j))
        partial = [[0]*5 for _ in range(5)]
        for parts in compositions(2-j,6):
            gaps = tuple(2*i+1 for i in parts)
            positions = [sum(gaps[:i]) for i in range(1,6)]
            for ordered in permutations(positions):
                current = tuple(sum(x*t for x,t in zip(col,ordered)) for col in zip(*single))
                other = tuple(den*sum(col)-x for x,col in zip(current,zip(*single)))
                # group[0] is the descendant. Use its five primary partners, including both currents.
                value, remainder = divmod(det((*group[1:],current,other)),den)
                assert not remainder
                for a,x in enumerate(ordered):
                    for b,y in enumerate(ordered):
                        partial[a][b] += value**2*(den*min(x,y)-x*y)
        for a,b in ((a,b) for a in range(5) for b in range(5)):
            ans[a][b] += weight*Fraction(partial[a][b],den**4)
    assert ans==moments(tuple(det((*group,p)) for p in single))


def check():
    single = tuple(tuple(int(i==j) for j in range(5)) for i in range(5))
    three = ((2,-1,1,0,1),(0,1,-2,1,0),(-1,0,2,1,-2))
    group = (*three,tuple(-x for x in add_rows((*single,*three))))
    simplex_test(group,single)
    coeff = tuple(det((*group,p)) for p in single)
    green = moments(coeff)
    for i,j in combinations(range(5),2):
        reference = sum(a*b*green_moments(tuple(int(h==k)+int(h==l) for h in range(5)),i,j)[1]
                        for k,a in enumerate(coeff) for l,b in enumerate(coeff))
        assert reference==green[i][j]
    p = (*group,*single)
    out = potentials(p,0)
    assert out[0]==(0,)*5 and any(any(row) for row in out)
    perm = (8,5,6,2,0,1,7,3,4)
    assert potentials(tuple(p[i] for i in perm),perm.index(0))==[out[i] for i in perm]
    half = tuple(tuple(Fraction(x,2) for x in row) for row in p)
    assert potentials(half,0)==[tuple(x/2**21 for x in row) for row in out]
    zero_mark = ((0,)*5,*p[1:-1],add_rows((p[0],p[-1])))
    assert not any(any(row) for row in potentials(zero_mark,0))
    print('Raw five-simplex moment matrix passed; 56 groups, permutation, scaling and dilaton limits passed.')
    print('Full finite-momentum potential at wave 5:',out[5])


if __name__=='__main__':
    check()
