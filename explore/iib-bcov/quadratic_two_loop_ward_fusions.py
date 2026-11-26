"""Fivefold genus-two Ward contributions on the correlated level-three background."""
from argparse import ArgumentParser
from collections import Counter
from fractions import Fraction
from itertools import product
from math import factorial, prod

from genus_two_graph_ward import projected_ward_difference, vertex
from two_loop_moving_descendant import BASIS


def fusion_orbits():
    # A pair counts donors from one transverse axis absorbed at each of the two A3 marks.
    options=((0,0),(1,0),(0,1),(1,1),(0,2),(2,0))
    multiplicities=Counter()
    for assignment in product(options,repeat=4):
        if any(sum(column)>3 for column in zip(*assignment)):
            continue
        swapped=((b,a) for a,b in assignment)
        orbit=min(tuple(sorted(assignment)),tuple(sorted(swapped)))
        multiplicities[orbit]+=prod(2//(factorial(2-a-b)*factorial(a)*factorial(b)) for a,b in assignment)
    return multiplicities


def ward_value(orbit):
    first,second=zip(*orbit)
    levels=(3-sum(first),3-sum(second))
    waves=((0,*first),(0,*second),BASIS[0],(-1,-2,-2,-2,-2),
           *(BASIS[i] for i,(a,b) in enumerate(orbit,1) for _ in range(2-a-b)))
    assert not any(map(sum,zip(*waves)))
    # The constant vector is a zero-charge primary in the scalar trace. The
    # string equation lowers each surviving descendant once. Its root pairing
    # is -1, which cancels the genus-two heat sign.
    result=Fraction()
    for i,k in enumerate(levels):
        if k:
            lowered=tuple(p-int(j==i) for j,p in enumerate(levels))
            result+=vertex(waves,(*lowered,*(0 for _ in waves[2:])))
    return result


SECTORS=(None,-Fraction(13483,1944),Fraction(14981,2160),-Fraction(42011,12960),
         Fraction(10757,15552),-Fraction(211,4320),Fraction())


def unfused_sample(step):
    """Move the two descendant momenta oppositely along the second minus third axis."""
    transfer=(0,step,-step,0,0)
    recoil=tuple(-x for x in transfer)
    waves=(transfer,recoil,BASIS[0],(-1,-2,-2,-2,-2),
           *(p for p in BASIS[1:] for _ in range(2)))
    # At precursor levels (2,2) the axis exchange swaps the two marked sites.
    # Their raised vertices agree. The difference of their Ward rows is
    # therefore twice their common vertex times the first marked momentum.
    return projected_ward_difference(waves,(2,2)+(0,)*10,0,1,BASIS[1])/(2*step)


def unfused_value(samples):
    # Each five-dimensional determinant is affine in this rank-one change.
    # Swapping the second and third axes acts trivially by the fourth determinant power.
    v1,v2,v3=samples
    return (15*v1-6*v2+v3)/10


def check(selected):
    orbits=fusion_orbits()
    for n in selected:
        active=[(q,w) for q,w in orbits.items() if sum(a+b for a,b in q)==n]
        value=sum((w*ward_value(q) for q,w in active),Fraction())
        # Identical donor waves contribute (2!)^4 and the identical hard marks another 2!.
        value*=Fraction((-1)**n,2**5)
        assert value==SECTORS[n]
        print('Fusions, orbits, labeled assignments and heat coefficient:',n,len(active),sum(w for _,w in active),value,flush=True)
    assert sum(SECTORS[1:])==-Fraction(67361,25920)


if __name__=='__main__':
    cli=ArgumentParser(description=__doc__)
    cli.add_argument('fusions',type=int,nargs='*',metavar='1..6')
    cli.add_argument('--unfused-sample',choices=(1,2,3),type=int)
    cli.add_argument('--unfused',action='store_true',help='evaluate the complete zero-fusion sector')
    args=cli.parse_args()
    if any(j not in range(1,7) for j in args.fusions):
        cli.error('Use fusion counts from one to six.')
    if args.unfused_sample:
        print(args.unfused_sample,unfused_sample(args.unfused_sample),flush=True)
    if args.unfused:
        value=unfused_value([unfused_sample(s) for s in (1,2,3)])
        assert value/16 + sum(SECTORS[1:]) == Fraction(1,36)
        print('Unfused vertex and complete fivefold response:',value,Fraction(1,36),flush=True)
    if args.fusions or not (args.unfused or args.unfused_sample):
        check(args.fusions or range(1,7))
