"""Sum the original K3,3 weights and assemble the fifty-shape elementary-field response."""
from argparse import ArgumentParser
from fractions import Fraction as F
from gzip import open as gzopen
from itertools import islice, permutations, product
from pathlib import Path

from quartic_prism_physical_ibp import FULL_TOTAL, verify_physical_ranges, zeta_interval

HERE = Path(__file__).resolve().parent
ONE_DOUBLE = ((F(11144483303,16895050672500000),F(197793221359,36248733480960000000)),
              (-F(186919,340540200000),-F(507937,111878807040000)), (F(0),F(0)))
TWO_DOUBLES = (-F(211,8214877440000),-F(1,6490767360000))
WEIGHT_FIVE = (F(109867352441461633,21778907950800000),F(13733372208504593,597364332364800000))
BIPARTITE = (
    (-F(4417821905516992411874825280217,58088006361319334400000000000),
     -F(4350142347269738714354694912071,11152897221373312204800000000000)),
    (-F(25047509866224749020693,123134594952600000000),
     -F(822126351867227342863,788061407696640000000)),
    (F(249372365272156031,806626220400000),F(2209968587039007449,1393850108851200000)),
)


def physical_channel(period_file,weights_file=HERE/"quartic_k33.weights.gz",row_count=465564):
    """Read one certified orbit column and match its three-by-three edge matrix to each weight."""
    import numpy as np
    from flint import fmpq

    with period_file.open() as f:
        values = {int(k):fmpq(v) for k,v in (line.split() for line in f)}
    if len(values) != 695485:
        raise ValueError("Expected the complete positive-orbit column")
    # Raw project edges use increasing vertex pairs. Interior keys use the minimum
    # row-major exponent matrix over graph isomorphisms, with a nibble per edge.
    edges = sorted(tuple(sorted((i,j))) for i in (0,2,4) for j in (1,3,5))
    maps = []
    for first,second in (((0,2,4),(1,3,5)),((1,3,5),(0,2,4))):
        for rows, cols in product(permutations(first),permutations(second)):
            order = [edges.index(tuple(sorted((i,j)))) for i in rows for j in cols]
            maps.append(order)
    radix = np.asarray([16**i for i in reversed(range(9))],dtype=np.int64)
    transforms = np.zeros((9,len(maps)),dtype=np.int64)
    for i,indices in enumerate(maps):
        transforms[indices,i] = radix
    totals = [fmpq(0),fmpq(0)]
    count = 0
    with gzopen(weights_file,"rt") as stream:
        while lines := list(islice(stream,2048)):
            data = [line.split() for line in lines]
            exponents = np.asarray([r[:9] for r in data],dtype=np.int64)
            keys = (exponents @ transforms).min(axis=1)
            for key,row in zip(keys,data):
                moment = values[int(key)]
                for i in range(2):
                    totals[i] += moment*fmpq(row[9+i])
            count += len(lines)
    assert count == row_count
    return tuple(F(str(value)) for value in totals)


def total_channels(bipartite):
    return tuple(tuple(x+y+z+w for x,y,z,w in zip(bipartite[i],FULL_TOTAL[i],ONE_DOUBLE[i],
                       TWO_DOUBLES if i == 0 else (0,0))) for i in range(3))


def residue(value, prime):
    return value.numerator * pow(value.denominator,-1,prime) % prime


def field_check():
    prime = 999999937
    observed = ((691097625,727592239),(224776989,616432627),(852726112,25526974))
    total = [tuple(residue(x,prime) for x in row) for row in total_channels(observed)]
    assert total == [(450469298,683141421),(396949588,272572399),(790650093,875668526)]
    assert total[2] == tuple(residue(x,prime) for x in WEIGHT_FIVE)
    assert (total[1][0]*total[2][1] - total[1][1]*total[2][0]) % prime == 323157369
    return total


def interval(coefficients, places=24):
    lo = hi = coefficients[0]
    for c,s in zip(coefficients[1:],(3,5)):
        a,b = zeta_interval(s,128)
        lo += min(c*a,c*b)
        hi += max(c*a,c*b)
    scale = 10**places
    return F((lo*scale).__floor__(),scale), F((hi*scale).__ceil__(),scale)


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("channels",nargs="*",type=Path,help="certified interior tables in 1,zeta(3),zeta(5) order")
    args = parser.parse_args()
    residues = field_check()
    print("full period residues in channel order",residues,flush=True)
    if args.channels:
        if len(args.channels) != 3:
            parser.error("Provide all three period tables")
        bipartite = []
        for path in args.channels:
            bipartite.append(physical_channel(path))
            print("K3,3", *bipartite[-1], flush=True)
        assert tuple(bipartite) == BIPARTITE
    values = total_channels(BIPARTITE)
    assert values[2] == WEIGHT_FIVE
    assert all(tuple(residue(x,999999937) for x in row)==r for row,r in zip(values,residues))
    verify_physical_ranges(values)
    for name,column in zip(("K4","J4"),zip(*values)):
        print(name,*column,"interval",*interval(column))
    elementary = tuple(F(9072,11)*j for k,j in values)
    print("reduced Q4 phi",*elementary,"interval",*interval(elementary))
    print("reduced one-flux vector per unit flux",*(40320*x for x in elementary))
