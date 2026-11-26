"""Lower bounds for closed even classes from polynomial products of quadratic traces."""

from collections import Counter, defaultdict
from itertools import product

from relative import Sector, Span, kernel, map_row, primitive


def quadratics(sector):
    per_weight = defaultdict(list)
    for i in range(0,len(sector.degrees),3):
        for j in range(i,len(sector.degrees),3):
            weight = tuple(a+b for a,b in zip(sector.degrees[i], sector.degrees[j]))
            if any(a>b for a,b in zip(weight, sector.weight)) or (i==j and sector.odd[i]):
                continue
            trace = Counter()
            for c,d,coeff in ((0,2,1),(1,1,2),(2,0,1)):
                sign, mon = sector.order((i+c,j+d))
                trace[mon] += coeff*sign
            per_weight[weight].append(primitive(trace))
    gen = {}
    for w, rows in per_weight.items():
        gens = []
        for relation in kernel(map_row(row, sector.differential) for row in rows):
            row = map_row(relation, rows.__getitem__)
            gens.append(primitive({m:v for m,v in row.items() if all(i%3==1 for i in m)}))
        if gens:
            gen[w] = gens
    return gen


def multiply(left, right, sector):
    result = Counter()
    for m, a in left.items():
        for n, b in right.items():
            sign, mon = sector.order(m+n)
            result[mon] += a*b*sign
    return primitive(result)


def count(bound, by_length=False):
    sector = Sector(bound)
    gens = quadratics(sector)
    # Process one generating quadratic at a time. The finite graded span
    # before multiplication plus its new multiples is the next subalgebra.
    zero = (0,)*5
    algebra = defaultdict(Span)
    algebra[zero].insert({():1})
    weights = sorted(product(*(range(x+1) for x in bound)), key=lambda w:(sum(w),w))
    for weight, polynomials in sorted(gens.items(), key=lambda item:(sum(item[0]),item[0])):
        relevant = [(w,tuple(i-j for i,j in zip(w,weight))) for w in weights
                    if all(i>=j for i,j in zip(w,weight))]
        for generator in polynomials:
            for target,lower in relevant:
                if lower in algebra:
                    space = algebra[target]
                    for row in algebra[lower].pivots.values():
                        space.insert(multiply(generator, row, sector))
    if by_length:
        return {w:Counter(len(mon) for mon in space.pivots) for w,space in algebra.items()}
    return {w:len(space.pivots) for w,space in algebra.items()}


if __name__ == "__main__":
    import sys
    from index_scan import index
    bound = tuple(map(int,sys.argv[1:])) if len(sys.argv)>1 else (2,1,1,2,2)
    counts = count(bound)
    euler = index(bound)
    gaps = sorted(((w,n,euler[w]) for w,n in counts.items() if euler[w]<n),
                  key=lambda item:(sum(item[0]),item[0]))
    print(len(counts), "product sectors,", sum(counts.values()), "classes;", len(gaps), "forced odd sectors", gaps[:12])
