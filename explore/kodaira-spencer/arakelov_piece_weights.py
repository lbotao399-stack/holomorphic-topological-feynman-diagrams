"""Merge the sorted projective weights of disjoint Arakelov numerator pieces."""
from argparse import ArgumentParser
from contextlib import ExitStack
from gzip import open as gzip_open
from heapq import merge
from itertools import groupby
from pathlib import Path
from tempfile import NamedTemporaryFile

from flint import fmpq


def rows(stream, multiplicity=1):
    previous = None
    for line in stream:
        *raw, value, other = line.split()
        powers = tuple(map(int,raw))
        if len(powers)!=9 or sum(powers)!=48 or not 0<min(powers)<=max(powers)<12 or other!='0':
            raise ValueError('Invalid genus-four scalar moment row.')
        if previous is not None and previous>=powers:
            raise ValueError('Input moment rows must be in strictly increasing order.')
        previous = powers
        yield powers,multiplicity*fmpq(value)


def merge_weights(inputs, output, multiplicities=None):
    inputs = [Path(p).resolve() for p in inputs]
    if not inputs or len(set(inputs))!=len(inputs) or Path(output).resolve() in inputs:
        raise ValueError('Use distinct input and output files.')
    if multiplicities is None:
        multiplicities = [1]*len(inputs)
    if len(multiplicities)!=len(inputs) or any(type(m) is not int or m<1 for m in multiplicities):
        raise ValueError('Provide one positive orbit multiplicity for each input.')
    count = 0
    with NamedTemporaryFile(dir=Path(output).parent,delete=False) as staged:
        pending = Path(staged.name)
    try:
        with ExitStack() as files, gzip_open(pending,'wt') as result:
            sources = [rows(files.enter_context(gzip_open(p,'rt')),m)
                       for p,m in zip(inputs,multiplicities)]
            combined = merge(*sources,key=lambda x:x[0])
            for powers,terms in groupby(combined,key=lambda x:x[0]):
                value = sum((c for _,c in terms),fmpq(0))
                if value:
                    print(*powers,value,0,file=result)
                    count += 1
        pending.replace(output)
        return count
    finally:
        pending.unlink(missing_ok=True)


if __name__=='__main__':
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    parser.add_argument('inputs',type=Path,nargs='+')
    parser.add_argument('--multiplicities',type=int,nargs='+',
                        help='Size of the represented signed-derivative orbit per input file.')
    args=parser.parse_args()
    print(merge_weights(args.inputs,args.output,args.multiplicities),'nonzero folded moments')
