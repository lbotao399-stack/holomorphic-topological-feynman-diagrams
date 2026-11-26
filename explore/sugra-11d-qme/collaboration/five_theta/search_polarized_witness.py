"""Search the five theta kernels for a nonzero component with transverse physical polarizations."""
from argparse import ArgumentParser
from fractions import Fraction
from itertools import permutations, product
from pathlib import Path
from random import Random
from tempfile import TemporaryDirectory
from replay import CASES, Commands, build, calculate_chunk

waves=(1,2,3,4,5,-1,-2,-4,-3,-5)

def sign(a):return (-1)**sum(x>y for i,x in enumerate(a) for y in a[i+1:])

def accepted(mask,polar):
    m=(1+mask.bit_count())//2
    # Two momentum columns per axis occur in D; every gamma also supplies its own momentum.
    for j in range(5):
        npol=sum(1 if mask>>l&1 else -1 for l,p in enumerate(polar) if p==j)
        ext=sum(bool(mask>>i&1) for i,w in enumerate(waves) if abs(w)==j+1)
        if m-npol-ext<0:return False
    # A signed permutation stabilizing all wave and polarization lines must act
    # consistently on the m volume tensors and on the external odd inputs.
    for perm in permutations(range(5)):
        for signs in product((-1,1),repeat=5):
            vp=[waves.index((1 if p>0 else -1)*signs[abs(p)-1]*(1+perm[abs(p)-1])) for p in waves]
            if any((mask>>n&1)!=(mask>>v&1) or perm[polar[n]]!=polar[v] for n,v in enumerate(vp)):continue
            ext=sign([v for n,v in enumerate(vp) if not mask>>n&1]);det=sign(perm)
            for x in signs:det*=x
            for p in polar:ext*=signs[p]
            if ext!=det**m:return False
    return True


def matches(shape):
    path=[j for j,n in enumerate(shape) for _ in range(n-1)];current=[-1]*10
    def gen(i):
        if i==5:yield tuple(current);return
        v=current.index(-1);current[v]=i
        for k in range(v+1,10):
            if path[v]==path[k] or current[k]!=-1:continue
            current[k]=i;yield from gen(i+1);current[k]=-1
        current[v]=-1
    words=sorted(gen(0));ranks={w:i for i,w in enumerate(words)};groups=[];seen=set()
    permuts=[p for p in permutations(range(3)) if all(shape[j]==shape[k] for j,k in enumerate(p))]
    for i,w in enumerate(words):
        if i in seen:continue
        paths=[[v for v,p in zip(w,path) if p==j] for j in range(3)];orbit=set()
        for p in permuts:
            for direction in (1,-1):
                a=[j for k in p for j in paths[k][::direction]];table={};canon=tuple(table.setdefault(l,len(table)) for l in a);orbit.add(ranks[canon])
        seen|=orbit;groups.append(sorted(orbit))
    return groups



def candidate(rng,order):
    # The two trials alternate. Start from the prior transverse unit fixture and
    # change four polarization slots, keeping every momentum-polarization pairing zero.
    mask=1023^(1<<rng.randrange(10)) if order==5 else sum(1<<i for i in rng.sample(range(10),3))
    while True:
        polar=[1,2,3,4,1,1,0,0,0,0]
        for slot in rng.sample(range(10),4):
            polar[slot]=rng.choice([i for i in range(5) if i!=abs(waves[slot])-1])
        if accepted(mask,polar):return mask,','.join(str(i+1) for i in polar)


def search(scratch,commands,attempts):
    rng=Random(68022)
    pending=[5,2]
    orbits={shape:matches(shape) for shape in CASES}
    for attempt in range(attempts):
        for order in list(pending):
            commands.check()
            mask,axes=candidate(rng,order)
            print(f'trial {attempt}, m={order}, --gamma-mask {mask:#x} --polarizations {axes}',flush=True)
            stored={}
            def at(shape,rank):
                if (shape,rank) not in stored:
                    stored[shape,rank]=calculate_chunk(scratch,shape,rank,rank+1,(4,),commands,mask,axes)[3]
                return stored[shape,rank]
            probes=[at(s,groups[len(groups)//2][0]) for s,groups in orbits.items()]
            # This screening heuristic may pass over nonzero evaluations. Any retained
            # result below uses all 768 matchings through complete automorphism orbits.
            if sum(x!=0 for x in probes)<3:continue
            weighted=Fraction()
            for shape,(count,automorphisms,_) in CASES.items():
                total=Fraction()
                assert sum(map(len,orbits[shape]))==count
                for index,group in enumerate(orbits[shape]):
                    value=at(shape,group[0])
                    # Reversing all paths and permuting equal-length paths relabels the
                    # integrated graph after all physical placements have been summed.
                    if index<2 and len(group)>1 and at(shape,group[-1])!=value:
                        raise ValueError(f'automorphism mismatch: {shape}, {group}')
                    total+=len(group)*value
                print(f'{shape} raw / pi = {total}',flush=True)
                weighted+=total/automorphisms
            print(f'five-graph sum / pi = {weighted}',flush=True)
            if weighted:pending.remove(order)
        if not pending:return True
    return False


if __name__=='__main__':
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('--trials',type=int,default=1,help='maximum trials per coupling order')
    args=parser.parse_args()
    if args.trials<1:parser.error('trials must be positive')
    commands=Commands()
    with commands.signals(), TemporaryDirectory(prefix='theta-witness-') as tmp:
        build(Path(tmp),commands)
        found=search(Path(tmp),commands,args.trials)
        commands.check()
    if not found:raise SystemExit('At least one coupling order has no nonzero evaluation within these trials.')
