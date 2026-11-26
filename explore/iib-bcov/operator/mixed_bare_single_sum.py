"""Sum original hard-label voltages and bubble moments before the tetrahedron cut."""

from collections import Counter,defaultdict
from fractions import Fraction as F
from functools import cache
from hashlib import sha256
from itertools import combinations,permutations
from math import gcd
from pathlib import Path
from tempfile import NamedTemporaryFile
import argparse
import fcntl
import json

from mixed_bare_ranks import WAVES
from mixed_bare_single import chamber_data,chamber_polynomials,shape_chambers
from mixed_bare_single_core import read_records
from mixed_bare_single_fold import project


@cache
def signed_voltages():
    # The Green potential reduces to -1/2 times the resistance contracted
    # against the source. For eight exchanges the mean times four therefore
    # equals one quarter of the integer signed resistance coefficients.
    rows = defaultdict(Counter)
    for order in permutations(range(6)):
        word = tuple((0,1,1,2,2,0)[i] for i in order)
        plus,minus = order.index(0),order.index(1)
        for site,label in enumerate(order):
            if site!=plus:
                rows[word][tuple(sorted((site,plus)))] -= WAVES[label][3]
            if site!=minus:
                rows[word][tuple(sorted((site,minus)))] += WAVES[label][3]
    assert all(c%4==0 for row in rows.values() for c in row.values())
    return {word:[(i,j,c//4) for (i,j),c in row.items() if c] for word,row in rows.items()}


def streamed_period(data,word,weight):
    """Keep six signed minors and a single soft cofactor in memory."""
    from mixed_bare_ranks import hard_waves
    from quartic_seven_plane_double import numerator,tree_flow
    from quartic_seven_plane_single import CYCLES

    graph,ctx,_,lengths,resistances,_,_ = data
    _,hard = hard_waves(word)
    charges = [(0,)*4]*6+[(int(a),int(b+c),int(d),int(c)) for a,b,c,d in hard]
    flow = tree_flow(graph,charges)
    pairs = list(combinations(range(4),2))
    minors = [numerator(graph,charges,flow,[(3,i),(3,j)]) for i,j in pairs]
    sectors = {}
    for i,j in combinations(range(9),2):
        ray = [CYCLES[a][i]*CYCLES[b][j]-CYCLES[a][j]*CYCLES[b][i] for a,b in pairs]
        content = gcd(*ray)
        if not content:
            continue
        if next(v for v in ray if v)<0:
            content = -content
        ray = tuple(v//content for v in ray)
        sectors[ray] = sectors.get(ray,0)+content**2*lengths[i]*lengths[j]
    raw = sum(c*resistances[6+i][6+j] for i,j,c in signed_voltages()[word])
    voltage = ctx.from_dict({e:int(c) for e,c in raw.terms()})
    pair = (F(),F())
    count = 0
    for ray,position in sectors.items():
        root = sum((v*f for v,f in zip(ray,minors) if v),ctx.constant(0))
        if not root:
            continue
        content = int(root.content())
        root = root//content
        # On one color the signed voltage is common to all nonnegative
        # squares. Each absolute summand is bounded by the full integrand.
        value = project(data,[(root,position)],voltage,sparse=True)
        pair = tuple(a+b*content**2*weight/2 for a,b in zip(pair,value))
        count += 1
    return pair,[word] if count else [],count,count


def grouped_period(data,rows,statistics=None,separate_roots=False):
    if separate_roots:
        if len(rows)!=1:
            raise ValueError('Independent soft-square projections require a single signed voltage.')
        return streamed_period(data,*rows[0])
    ctx,resistances = data[1],data[4]
    buckets = defaultdict(list)
    roots = 0
    nonzero = []
    for word,weight in rows:
        squares,formula = chamber_polynomials(data,word)
        if not squares:
            continue
        nonzero.append(word)
        raw = sum(c*resistances[6+i][6+j] for i,j,c in signed_voltages()[word])
        voltage = ctx.from_dict({k:int(c) for k,c in raw.terms()})
        assert voltage==formula and weight.denominator==1
        for root,soft in squares:
            roots += 1
            key = len(root),int(root(*range(2,17)))
            for row in buckets[key]:
                if row[0]==root:
                    row[1] += weight.numerator*soft*voltage
                    break
            else:
                buckets[key].append([root,weight.numerator*soft*voltage])
    # All colors now share one basis of bubble and exterior polynomials.
    # Equal triples across colors contract only once inside project().
    terms = [row for bucket in buckets.values() for row in bucket if row[1]]
    pair = project(data,terms,ctx.constant(1),statistics=statistics) if terms else (F(),F())
    return tuple(x/2 for x in pair),nonzero,roots,len(terms)


def save(record):
    bank = Path(__file__).with_suffix('.jsonl')
    with bank.open('a+') as stream:
        fcntl.flock(stream,fcntl.LOCK_EX)
        stream.seek(0)
        for line in stream:
            existing = json.loads(line)
            if existing['paths']==record['paths']:
                assert existing['period']==record['period'] and existing['nonzero']==record['nonzero']
                return
        stream.write(json.dumps(record,separators=(',',':'))+'\n')
        stream.flush()


def integrate(marks,rows,batch_size,write=False,checkpoint=None,reuse_colors=False,separate_roots=False):
    if write and (reuse_colors or checkpoint is not None):
        bank = Path(__file__).with_suffix('.jsonl')
        if bank.exists():
            with bank.open() as stream:
                fcntl.flock(stream,fcntl.LOCK_SH)
                for line in stream:
                    record = json.loads(line)
                    if tuple(record['paths'])==marks:
                        print('Already integrated',marks,flush=True)
                        return record
    data = chamber_data(marks)
    saved = read_records()
    total = [F(),F()]
    nonzero = []
    weights = dict(rows)
    if reuse_colors:
        remaining = []
        for word,weight in rows:
            if (marks,word) not in saved:
                remaining.append((word,weight))
                continue
            previous,a,b = saved[marks,word]
            assert previous==weight
            if a or b or chamber_polynomials(data,word)[0]:
                nonzero.append(word)
            total[0] += a
            total[1] += b
        rows = remaining
        # A dense tile can exceed the worker bound. Retain all disjoint
        # completed tiles even if a retry changes their size or row order.
        if checkpoint is not None and checkpoint.exists():
            available = dict(rows)
            path_word = ''.join(map(str,marks))
            files = sorted(checkpoint.glob(f'{path_word}-*.json'))
            for file in files:
                entry = json.loads(file.read_text())
                inputs = [(tuple(map(int,w)),F(v)) for w,v in entry['inputs']]
                if not all(w in available for w,v in inputs):
                    continue
                assert all(available[w]==v for w,v in inputs)
                words = [tuple(map(int,w)) for w in entry['nonzero']]
                assert len(words)==len(set(words)) and set(words)<=dict(inputs).keys()
                for w,v in inputs:
                    del available[w]
                nonzero.extend(words)
                total = [a+F(b) for a,b in zip(total,entry['period'],strict=True)]
            rows = [(w,v) for w,v in rows if w in available]
    nroots = ngroups = 0
    for start in range(0,len(rows),batch_size):
        block = rows[start:start+batch_size]
        inputs = [[''.join(map(str,w)),str(v)] for w,v in block]
        filename = None
        if checkpoint is not None:
            checkpoint.mkdir(exist_ok=True,parents=True)
            path_word = ''.join(map(str,marks))
            stem = f'{path_word}-{batch_size}-{start}'
            original = checkpoint/f'{stem}.json'
            selection = sha256(json.dumps(inputs).encode()).hexdigest()[:16]
            exact = checkpoint/f'{stem}-{selection}.json'
            filename = exact if reuse_colors else original
            alternate = original if reuse_colors else exact
            if not filename.exists() and alternate.exists():
                if json.loads(alternate.read_text())['inputs']==inputs:
                    filename = alternate
        if filename is not None and filename.exists():
            entry = json.loads(filename.read_text())
            assert entry['inputs']==inputs
            pair = tuple(map(F,entry['period']))
            words = [tuple(map(int,w)) for w in entry['nonzero']]
            roots,groups = entry['counts']
            assert len(pair)==2 and len(set(words))==len(words) and set(words)<=dict(block).keys()
            assert roots>=groups and roots>=len(words)
        else:
            pair,words,roots,groups = grouped_period(data,block,separate_roots=separate_roots)
            if filename is not None:
                entry = dict(inputs=inputs,period=list(map(str,pair)),
                             nonzero=[''.join(map(str,w)) for w in words],counts=[roots,groups])
                with NamedTemporaryFile(mode='w',dir=filename.parent,delete=False) as out:
                    json.dump(entry,out,separators=(',',':'))
                Path(out.name).replace(filename)
        if all((marks,w) in saved for w in words):
            direct = tuple(sum((saved[marks,w][j] for w in words),F()) for j in (1,2))
            assert pair==direct
        nonzero.extend(words)
        nroots += roots
        ngroups += groups
        total = [a+b for a,b in zip(total,pair)]
        print('paths',marks,'color representatives',start+len(block),'of',len(rows),flush=True)
    labeled = sum(8*weights[w] for w in nonzero)
    assert labeled.denominator==1
    record = dict(paths=list(marks),nonzero=len(nonzero),labeled=labeled.numerator,
                  period=list(map(str,total)))
    if write:
        save(record)
    print(record,'squares before and after aggregation',nroots,ngroups,flush=True)
    return record


def coverage():
    from mixed_bare_symbolic import classified_orbits

    physical,_ = classified_orbits(2)
    by_shape = defaultdict(dict)
    for shape,word,multiplicity,weight in physical:
        by_shape[shape][word] = (multiplicity,weight)
    assert len(physical)==3995 and len(by_shape)==102
    saved = read_records()
    for (shape,word),(weight,a,b) in list(saved.items()):
        if word in by_shape[shape]:
            assert by_shape[shape][word][1]==weight
        else:
            assert not a and not b
            del saved[shape,word]
    aggregate = {}
    bank = Path(__file__).with_suffix('.jsonl')
    if bank.exists():
        with bank.open() as stream:
            fcntl.flock(stream,fcntl.LOCK_SH)
            for line in stream:
                record = json.loads(line)
                shape = tuple(record['paths'])
                assert shape not in aggregate and record['nonzero']==len(by_shape[shape])
                assert record['labeled']==sum(n for n,w in by_shape[shape].values())
                aggregate[shape] = tuple(map(F,record['period']))
    total = [F(),F()]
    represented = complete = 0
    classes = defaultdict(list)
    for shape,words in by_shape.items():
        rows = [saved[shape,w] for w in words if (shape,w) in saved]
        pair = [sum((r[j] for r in rows),F()) for j in (1,2)]
        if shape in aggregate:
            if len(rows)==len(words):
                assert tuple(pair)==aggregate[shape]
            pair = aggregate[shape]
            represented += len(words)
            complete += 1
        else:
            represented += len(rows)
            complete += len(rows)==len(words)
        total = [a+b for a,b in zip(total,pair)]
        # The two parallel arms and the two adjacent connectors form an
        # intrinsic subset invariant under every core automorphism.
        arms = sorted((shape.count(0),shape.count(1)),reverse=True)
        signature = (*arms,shape.count(2)+shape.count(3))
        classes[signature].append((shape in aggregate or len(rows)==len(words),
                                   len(words),sum(n for n,w in words.values()),*pair))
    print('Individual signed periods:',len(saved),'of',len(physical))
    print('Covered by exact sums:',represented,'with',complete,'complete placements out of',len(by_shape))
    print('Subtotal without double counting:',*total)
    for signature,rows in sorted(classes.items()):
        if all(r[0] for r in rows):
            print('Complete parallel/connector population',signature,'placements',len(rows),
                  'chambers, labelings, rational and zeta(3) coefficients',
                  *(sum(r[j] for r in rows) for j in range(1,5)))
    assert {signature[2] for signature in classes}=={0,1}
    connector = [r for signature,rows in classes.items() if signature[2] for r in rows]
    if all(r[0] for r in connector):
        count,labeled,a,b = (sum(r[j] for r in connector) for j in range(1,5))
        assert (len(connector),count,labeled)==(40,1527,95472)
        from quartic_prism_physical_ibp import zeta_interval
        low,high = sorted(a+b*z for z in zeta_interval(3))
        assert 1572661674 < 10**20*low <= 10**20*high < 1572661675
        print('All nonempty connectors:',a,b)
    if complete==102:
        assert represented==3995 and sum(r[2] for rows in classes.values() for r in rows)==248928
        # Absorb the one-double term after the two-double, theta and pentagon
        # contributions have been inserted in the mixed D3 Ward equation.
        a,b = total
        fixed = F(3515607603589,597906382071398400000),-F(469153,98075577600000)
        known = fixed[0]+F(3,16)*a,fixed[1]+F(3,16)*b
        lower,upper = sorted(a+b*z for z in zeta_interval(3))
        assert -7513250331 < 10**20*lower <= 10**20*upper < -7513250330
        lower,upper = sorted(known[0]+known[1]*z for z in zeta_interval(3))
        assert 11563257517 < 10**20*lower <= 10**20*upper < 11563257518
        print('Mixed Ward coefficients before the two primitive cores:',*known)
    return len(saved),represented,complete,tuple(total)


def check_projection():
    from fourfold_quartic_single import project_squares
    from mixed_bare_single import check_polynomials

    word = (0,0,2,1,1,2)
    tests = (((0,1,2,8,8,8),1,494516),((0,1,2,4,5,8),1,56790),
             ((0,0,0,4,7,8),7,146622))
    for path,count,residue in tests:
        data = chamber_data(path)
        squares,voltage = chamber_polynomials(data,word)
        assert len(squares)==count and check_polynomials(data,word,squares,voltage)[2]==residue
        # This projection multiplies the exterior polynomials explicitly;
        # the grouped route contracts their five path moments at integer scale.
        pair = project_squares(data,squares,voltage)
        actual,words,_,_ = grouped_period(data,[(word,F(8))])
        assert words==[word] and actual==tuple(4*x for x in pair)
        assert streamed_period(data,word,F(8))[:2]==(actual,words)
        print('Direct electrical and exterior-polynomial control:',path,*pair)
    path = (0,0,0,2,4,7)
    rows = shape_chambers()[path]
    data = chamber_data(path)
    separate = Counter()
    direct = [F(),F()]
    for row in rows:
        sizes = {}
        pair,_,_,_ = grouped_period(data,[row],sizes)
        separate.update(sizes)
        direct = [a+b for a,b in zip(direct,pair)]
    combined = {}
    pair,words,roots,groups = grouped_period(data,rows,combined)
    assert pair==tuple(direct) and (len(words),roots,groups)==(51,51,45)
    assert separate['exterior_products']==735 and combined['exterior_products']==101
    assert separate['exterior_rows']==576857 and combined['exterior_rows']==85001
    print('All 51 signed colors:',pair,'exterior products and rows before and after grouping:',
          separate['exterior_products'],combined['exterior_products'],
          separate['exterior_rows'],combined['exterior_rows'])


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('indices',nargs='*',type=int)
    parser.add_argument('--block',type=int,default=4)
    parser.add_argument('--write',action='store_true')
    parser.add_argument('--checkpoint',type=Path,help='Resume complete color blocks from this directory')
    parser.add_argument('--reuse-colors',action='store_true',help='Consume saved individual chamber integrals before grouping the rest')
    parser.add_argument('--separate-roots',action='store_true',help='Evaluate one soft square at a time when --block 1 is used')
    parser.add_argument('--controls',action='store_true',help='Check direct currents and explicit exterior tensor products on three placements')
    args = parser.parse_args()
    if args.block<1:
        parser.error('The block size must be positive.')
    if args.separate_roots and args.block!=1:
        parser.error('--separate-roots requires --block 1.')
    if args.controls:
        check_projection()
    elif args.indices:
        table = list(shape_chambers().items())
        for index in args.indices:
            if args.checkpoint is None:
                integrate(*table[index],args.block,args.write,None,args.reuse_colors,args.separate_roots)
                continue
            args.checkpoint.mkdir(parents=True,exist_ok=True)
            path_word = ''.join(map(str,table[index][0]))
            with (args.checkpoint/f'{path_word}.lock').open('a') as lock:
                try:
                    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
                except BlockingIOError:
                    print('Active worker on',table[index][0],flush=True)
                    continue
                integrate(*table[index],args.block,args.write,args.checkpoint,args.reuse_colors,args.separate_roots)
    else:
        coverage()
