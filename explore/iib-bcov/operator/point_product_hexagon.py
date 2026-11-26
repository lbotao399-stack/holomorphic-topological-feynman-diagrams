"""One-loop marked products and the required correction of two protected BCOV point classes."""
from fractions import Fraction
from functools import cache
from itertools import permutations
from math import factorial, prod

from point_products import exterior


def determinant(rows):
    if len(rows)==1:
        return rows[0][0]
    return sum((-1)**j*x*determinant([row[:j]+row[j+1:] for row in rows[1:]])
               for j,x in enumerate(rows[0]) if x)


def currents(momenta,times):
    # Edge e runs from vertex e to e+1, with the auxiliary edge 5 from A to B.
    jumps = [[sum(p[i] for p in momenta[:e]) for i in range(5)] for e in range(6)]
    mean = [Fraction(sum(t*p[i] for t,p in zip(times,jumps)),sum(times)) for i in range(5)]
    return [[x-y for x,y in zip(p,mean)] for p in jumps]


def exterior_kernel(flow,fields=None):
    # With two scalar ends only the four bulk theta variables survive.
    # With traced marks the moving fifth vertex also has a theta variable.
    last = 4 if fields is None else len(fields)
    value = {0:1}
    for edge in range(5 if fields is None else last):
        delta = []
        for coordinate in range(5):
            linear = {}
            if edge:
                linear[1<<(5*(edge-1)+coordinate)] = 1
            if edge<last:
                linear[1<<(5*edge+coordinate)] = -1
            delta.append(linear)
        propagator = {}
        for missing in range(5):
            term = {0:(-1)**missing*flow[edge][missing]}
            for i in range(5):
                if i!=missing:
                    term = exterior(term,delta[i])
            for mask,coefficient in term.items():
                propagator[mask] = propagator.get(mask,0)+coefficient
        value = exterior(value,propagator)
        if edge:
            if fields is not None:
                value = exterior(value,{sum(1<<(5*(edge-1)+c) for c in fields[edge-1]):1})
            filled = (1<<(5*edge))-1
            value = {k:v for k,v in value.items() if k&filled==filled}
    if fields is not None:
        value = exterior(value,{sum(1<<(5*(last-1)+c) for c in fields[-1]):1})
    return value.get((1<<(5*last))-1,0)


def certificate():
    frame = [[int(i==j) for i in range(5)] for j in range(5)]
    mixed = [[2,-1,1,0,1],[1,2,0,0,-1],[0,-1,3,0,1],[-1,0,1,2,0],[0,2,1,1,1]]
    for momenta,times in [(frame,(1,2,3,4,5,6)),(mixed,(7,2,1,4,6,3))]:
        flow = currents(momenta,times)
        volume = determinant(momenta)
        physical = exterior_kernel(flow)
        assert physical==determinant(flow[:5])==-Fraction(times[5],sum(times))*volume
        heat = -volume
        assert heat*physical==Fraction(times[5],sum(times))*volume**2
        # The normalized first simplex moment is 1/6, its volume 1/5!.
        ordered = Fraction(volume**2,6*factorial(5))
        assert sum(ordered for _ in permutations(momenta))==Fraction(volume**2,6)

    # Reversing the inputs leaves the moving mark with no surviving letter.
    zero_mark = [*frame[:4],[0]*5]
    assert determinant(zero_mark)==0
    assert exterior_kernel(currents(zero_mark,(5,7,3,1,8,2)))==0
    print('Rooted hexagon has squared polarized reduced value det(p)^4/36; reversed scalar projection is zero.')


def traced_marks(vector_axis=0,marked_fields=None):
    # Put the five active momenta e^i successively on the path from the
    # rooted mark. Wave zero carries e2 e3 e4 e5 and wave four carries e1.
    # At the root the two constant spectators fill the polyvector trace,
    # so the contracted field there is scalar.
    full = Fraction()
    for wave_order in permutations(range(5)):
        positions = [wave_order.index(wave)+1 for wave in range(5)]
        sign = (-1)**sum(i>j for at,i in enumerate(wave_order) for j in wave_order[at+1:])
        # Column i fills the missing edge direction i with the external
        # polyvector of wave 4 (i=0) or 0 (i=1,...,4).
        if marked_fields is None:
            columns = [(4,vector_axis),*((0,i) for i in (1,2,3,4))]
        else:
            columns = [(label,i) for label,blade in marked_fields for i in blade]
        source_sites = [positions[label] for label,_ in columns]
        coordinates = [i for _,i in columns]
        # Each current component is an integer jump minus one resistance tail.
        # Expand the determinant in its five columns before simplex integration.
        average = Fraction()
        for subset in range(32):
            tails = sorted([positions[c] for i,c in enumerate(coordinates) if subset>>i&1],
                           reverse=True)
            matrix = [[int(edge<source_sites[i])*(-1 if subset>>i&1
                       else int(edge>=positions[c])) for i,c in enumerate(coordinates)]
                      for edge in range(5)]
            moment = prod(Fraction(6-site+j,6+j) for j,site in enumerate(tails))
            average += determinant(matrix)*moment
        full -= sign*average/factorial(5)
    return full


def traced_cubature(vector_axis=None):
    # The ordinary-product cross term has six physical edges and needs
    # degree-seven cubature; the auxiliary marked product needs degree five.
    degree = 7 if vector_axis is None else 5
    half = (degree-1)//2
    def compositions(total,length):
        if length==1:
            yield (total,)
            return
        for i in range(total+1):
            for rest in compositions(total-i,length-1):
                yield (i,*rest)

    nodes = []
    for j in range(half+1):
        den = 6+2*half-2*j
        weight = Fraction((-1)**j*den**degree,4**half*factorial(j)*factorial(6+2*half-j))
        for beta in compositions(half-j,6):
            nodes.append((tuple(2*x+1 for x in beta),den,weight))
    for powers in compositions(degree,6):
        expected = Fraction(prod(factorial(k) for k in powers),factorial(5+degree))
        assert sum(weight*Fraction(prod(x**k for x,k in zip(t,powers)),den**degree)
                   for t,den,weight in nodes)==expected

    frame = [tuple(int(i==j) for i in range(5)) for j in range(5)]
    total = Fraction()
    for sequence,order in enumerate(permutations(range(5))):
        momenta = [frame[i] for i in order]
        if vector_axis is None:
            fields = [(0,1) if i==4 else () for i in order]+[(1,2,3,4)]
        else:
            fields = [{0:(1,2,3,4),4:(vector_axis,)}.get(i,()) for i in order]
        heat = -determinant(momenta)
        size = len(fields)
        positions = [order.index(i)+1 for i in range(5)]
        columns = [(j,c) for j,blade in enumerate(fields,1) for c in blade]
        for node,(times,den,weight) in enumerate(nodes):
            flow = [[(den if positions[i]<=e else 0)-sum(times[positions[i]:])
                     for i in range(5)] for e in range(size)]
            physical = determinant([[flow[e][c] if e<j else 0 for j,c in columns]
                                    for e in range(size)])
            if node==0 and sequence in (6,17,97,119):
                assert physical==exterior_kernel(flow,fields)
            total += weight*heat*Fraction(physical,den**size)
    return total


def ordinary_leibniz_defect():
    # Split the two ends of the physical cycle between the two point traces.
    # The rooted rank-two and rank-three spectators fill E1, leaving a scalar
    # end. E0 receives the constant four-vector. The only active nonscalar
    # is the fifth wave with bivector e1 e2 and momentum e^5.
    integer_minor = cache(determinant)
    full = Fraction()
    coordinates = (0,1,1,2,3,4)
    for order in permutations(range(5)):
        positions = [order.index(i)+1 for i in range(5)]
        heat = -(-1)**sum(i>j for at,i in enumerate(order) for j in order[at+1:])
        sites = (positions[4],positions[4],6,6,6,6)
        for subset in range(64):
            matrix = tuple(tuple(int(edge<sites[i])*(-1 if subset>>i&1
                           else int(edge>=positions[c])) for i,c in enumerate(coordinates))
                           for edge in range(6))
            minor = integer_minor(matrix)
            if minor:
                tails = sorted([positions[c] for i,c in enumerate(coordinates) if subset>>i&1],
                               reverse=True)
                moment = prod(Fraction(6-site+j,6+j) for j,site in enumerate(tails))
                full += heat*minor*moment/factorial(5)
    return full


if __name__=='__main__':
    certificate()
    factors = [traced_marks(i) for i in range(5)]
    assert factors==[Fraction(1,60),*[Fraction(1,45)]*3,Fraction(1,30)]
    assert traced_cubature(0)==factors[0] and traced_cubature(1)==factors[1]
    print('E0 and E1 traced hexagon, signed heat times exterior factors:',factors)
    defect = ordinary_leibniz_defect()
    assert defect==traced_cubature()==-Fraction(1,24)
    print('Ordinary product one-loop Leibniz defect, oriented current moment:',defect)
    fusions = []
    for label in (1,2,3):
        remaining = tuple(i for i in (1,2,3,4) if i!=label)
        fusions.append((-1)**(label-1)*traced_marks(marked_fields=[(label,remaining),(4,(0,1))]))
    assert fusions==[Fraction(1,72)]*3 and defect+sum(fusions)==0
    print('Classical variation through the three nonzero fused waves:',fusions,sum(fusions))
