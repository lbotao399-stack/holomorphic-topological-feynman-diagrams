"""Closed homogeneous positive-orbit solution with zero four-wheel boundary."""
from argparse import ArgumentParser
from fractions import Fraction
from math import factorial, prod
from pathlib import Path
from struct import Struct


def value(key,q):
    powers=tuple((key>>(4*(8-i))) & 15 for i in range(9))
    rows=[sum(powers[i:i+3]) for i in (0,3,6)]
    columns=[sum(powers[i::3]) for i in range(3)]
    limit=2*q-1
    if limit not in rows or limit not in columns:
        return 0
    i,j=rows.index(limit),columns.index(limit)
    if powers[3*i+j]!=q-1:
        return 0
    weight=prod(factorial(limit-1-s) for s in rows+columns if s!=limit)
    denominator=prod(factorial(a-1) for a in powers)
    return Fraction(weight*factorial(q-3)*factorial(q-2)**3,
                    denominator*factorial(2*q-5)**2)


def check(q,base,reference=None):
    from math import lcm
    values={}
    known={}
    with Path(str(base)+'.nodes').open() as stream:
        for idx,line in enumerate(stream):
            key=int(line)
            result=value(key,q)
            if result:
                values[idx]=result
                if reference:known[key]=result
    if reference:
        with reference.open() as stream:
            assert {int(k):Fraction(x) for k,x in map(str.split,stream)}==known
    scale=lcm(*(x.denominator for x in values.values()))
    scaled={idx:int(scale*v) for idx,v in values.items()}
    record=Struct('<ib');checked=0
    with Path(str(base)+'.bin').open('rb') as stream:
        while header:=stream.read(1):
            n=header[0];data=stream.read(n*record.size)
            if len(data)!=n*record.size:raise ValueError('short identity')
            assert sum(c*scaled.get(idx,0) for idx,c in record.iter_unpack(data))==0
            checked+=1
    return len(scaled),checked


if __name__=='__main__':
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('q',type=int)
    parser.add_argument('base',type=Path)
    parser.add_argument('--reference',type=Path)
    args=parser.parse_args()
    print('nonzero orbits, exact identities',*check(args.q,args.base,args.reference))
