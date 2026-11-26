"""Exact three-line integration of unique K3,3 vertices onto prism position periods."""
from argparse import ArgumentParser
from fractions import Fraction
from math import factorial, prod
from pathlib import Path
from subprocess import run


def star_charts(key,q,unique_only=False):
    powers=tuple((key>>(4*(8-j)))&15 for j in range(9))
    matrix=tuple(powers[i:i+3] for i in (0,3,6))
    for rows in (matrix,tuple(zip(*matrix))):
        for i,triangle in enumerate(rows):
            deficit=q-sum(triangle)
            if deficit>=0 and (not unique_only or not deficit):
                yield deficit,triangle,rows[(i+1)%3],rows[(i+2)%3]


def prism_powers(q,first,second,edges):
    # Points 0 and 1 are the two remaining vertices on the original side.
    return (edges[2],q-first[0],q-first[1],q-second[1],q-second[0],
            q-first[2],edges[1],edges[0],q-second[2])


def chart_expansion(q,deficit,triangle,first,second):
    scale=Fraction(prod(factorial(q-a-1) for a in first+second),
                   prod(factorial(a-1) for a in first+second))
    terms=[]
    for k in range(deficit+1):
        for l in range(deficit+1-k):
            composition=(k,l,deficit-k-l)
            edges=tuple(a+deficit-j for a,j in zip(triangle,composition))
            moment=factorial(deficit)//prod(factorial(j) for j in composition)
            moment*=prod(factorial(b-1)//factorial(a-1) for b,a in zip(edges,triangle))
            powers=prism_powers(q,first,second,edges)
            assert min(powers)>0 and max(powers)<q and sum(powers)==4*q
            terms.append((powers,scale*moment))
    return terms


def star_expansions(key,q,unique_only=False):
    for chart in star_charts(key,q,unique_only):
        yield chart[0],chart_expansion(q,*chart)


def preferred_star(key,q):
    # At a fixed deficit the first prism term recovers all three original rows.
    # This chooses exactly the lexicographically least expansion without
    # constructing any of the alternatives.
    def order(chart):
        d,triangle,first,second=chart
        return d,prism_powers(q,first,second,(triangle[0]+d,triangle[1]+d,triangle[2]))
    chart=min(star_charts(key,q),key=order)
    return chart_expansion(q,*chart)


def stars(key,q):
    for _,terms in star_expansions(key,q,unique_only=True):
        yield terms[0]


def emit(keys, q, queries, instructions):
    lookup = {}
    count = 0
    with keys.open() as source, queries.open('w') as target, instructions.open('w') as recipe:
        for line in source:
            key = int(line)
            for position, scale in stars(key,q):
                if position not in lookup:
                    lookup[position] = len(lookup)
                    print(*position, file=target)
                print(key,lookup[position],scale,file=recipe)
                count += 1
    print('star tests',count,'different position queries',len(lookup),flush=True)


def collect(instructions, periods, output):
    with periods.open() as source:
        values = [tuple(map(Fraction,line.split())) for line in source]
    assert all(len(v)==3 for v in values)
    known = {}
    with instructions.open() as source:
        for line in source:
            key, index, scale = line.split()
            key = int(key)
            value = tuple(Fraction(scale)*v for v in values[int(index)])
            if key in known and known[key] != value:
                raise ArithmeticError(f'incompatible star reductions at {key}')
            known[key] = value
    with output.open('w') as target:
        for key, row in sorted(known.items()):
            print(key,*row,file=target)
    print('exact positive bipartite orbit anchors',len(known),flush=True)


def check_general(q, channels, native, prism_prefix):
    """Compare a full prism expansion at every available surplus degree."""
    from contextlib import ExitStack
    with ExitStack() as streams:
        entries = [streams.enter_context(Path(str(channels)+str(j)+'.tsv').open()) for j in range(3)]
        samples = {}
        for records in zip(*entries,strict=True):
            key = int(records[0].split()[0])
            matrix = [(key>>(4*(8-i))) & 15 for i in range(9)]
            deficits = {q-sum(matrix[i*3:i*3+3]) for i in range(3)}
            deficits.update(q-sum(matrix[i::3]) for i in range(3))
            if all(h<0 or h in samples for h in deficits):
                continue
            exact = tuple(Fraction(line.split()[1]) for line in records)
            for h,terms in star_expansions(key,q):
                samples.setdefault(h,(terms,exact))
    if set(samples)!=set(range(q-2)):
        raise ArithmeticError('incomplete surplus degree coverage')
    queries=[p for h in sorted(samples) for p,c in samples[h][0]]
    data=''.join(' '.join(map(str,row))+'\n' for row in queries)
    proc=run([str(native),'--all-periods',str(q),'--query',str(prism_prefix)],
             text=True,input=data,capture_output=True,check=True)
    responses=iter(proc.stdout.splitlines())
    for h,(terms,expected) in sorted(samples.items()):
        result=[Fraction(0) for _ in expected]
        for _,factor in terms:
            ans=tuple(map(Fraction,next(responses).split()))
            for j,v in enumerate(ans):
                result[j]+=factor*v
        assert tuple(result)==expected,(h,result,expected)
    assert next(responses,None) is None
    print('surplus degrees',sorted(samples),'prism terms',len(queries),flush=True)


if __name__ == '__main__':
    parser=ArgumentParser(description=__doc__)
    ops=parser.add_subparsers(dest='op',required=True)
    send=ops.add_parser('emit')
    send.add_argument('q',type=int)
    for k in ('nodes','queries','instructions'):
        send.add_argument(k,type=Path)
    receive=ops.add_parser('collect')
    for k in ('instructions','periods','output'):
        receive.add_argument(k,type=Path)
    check=ops.add_parser('general')
    check.add_argument('q',type=int)
    for k in ('channel_prefix','prism_binary','prism_prefix'):
        check.add_argument(k,type=Path)
    args=parser.parse_args()
    if args.op=='emit':
        emit(args.nodes,args.q,args.queries,args.instructions)
    elif args.op=='collect':
        collect(args.instructions,args.periods,args.output)
    else:
        check_general(args.q,args.channel_prefix,args.prism_binary,args.prism_prefix)
