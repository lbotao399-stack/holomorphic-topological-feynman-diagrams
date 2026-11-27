"""Independent exact replay of the complete fixed-weight GL(2) calculation.

Basis completeness uses a multigraded product enumeration.  Ranks use
descending-pivot elimination, independently of the production kernel
construction.  Q0 columns are recomputed by the original Wick engine.
Q2 uses its order-two/order-three differential-operator expansion; a
complete local commutator check proves descent to the entire complex.
"""
from collections import defaultdict
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import sys
import time

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import gl2_cohomology as production
from flint import fmpq


def add(a,b,c=1):
    for k,v in b.items():
        q=a.get(k,0)+c*v
        if q:a[k]=q
        else:a.pop(k,None)


def rank(vectors):
    pivots={}
    for vector in vectors:
        row={k:production.exact(v) for k,v in vector.items() if v}
        while row:
            pivot=max(row)
            if pivot not in pivots:
                lead=row[pivot]
                pivots[pivot]={k:v/lead for k,v in row.items()}
                break
            add(row,pivots[pivot],-row[pivot])
    return len(pivots)


def polynomial(raw):
    return {tuple(map(tuple,m)):fmpq(c) for m,c in raw}


def dimensions(complex_):
    """Coefficient of each weight in the exact bosonic/fermionic product."""
    weight=complex_.weight
    counts={(0,0,0,0,0):1}
    for derivative in range(weight//4+1):
        for name,dimension in (("c",4),("B",4),("p0",2),("s0",2)):
            for colour in range(dimension):
                letter=(name,colour,derivative)
                w=complex_.jw(letter)
                if w>weight:continue
                r=complex_.rw(letter); a=complex_.af(letter); h=complex_.charges(letter)
                odd=complex_.algebra.parity(letter)
                new=defaultdict(int)
                for (j,R,af,h1,h2),multiplicity in counts.items():
                    maximum=1 if odd else (weight-j)//w
                    for exponent in range(maximum+1):
                        if j+exponent*w>weight:break
                        key=(j+exponent*w,R+exponent*r,af+exponent*a,
                             h1+exponent*h[0],h2+exponent*h[1])
                        new[key]+=multiplicity
                counts=dict(new)
    return {(R,af):n for (j,R,af,h1,h2),n in counts.items()
            if j==weight and h1==h2==0}


class FullOperator:
    def __init__(self,complex_):
        self.complex=complex_; self.algebra=complex_.algebra
        self.engine=complex_.engine
        self.cache={}

    def q2_monomial(self,m):
        if m in self.cache:return self.cache[m]
        out={}
        for size,kind in ((2,"tad"),(3,"fan")):
            for selected in combinations(range(len(m)),size):
                part=tuple(m[i] for i in selected)
                if sum(self.complex.af(x) for x in part)<2:continue
                remaining=tuple(i for i in range(len(m)) if i not in selected)
                sign=self.algebra.perm_sign(m,selected+remaining)
                rest=tuple(m[i] for i in remaining)
                for image,c in self.engine.Q2_mono(part,only=kind).items():
                    s,target=self.algebra.canon(image+rest)
                    if s:add(out,{target:production.exact(c)},sign*s)
        self.cache[m]=out
        return out

    def q2(self,p):
        out={}
        for m,c in p.items():add(out,self.q2_monomial(m),c)
        return out


def verify(path):
    started=time.monotonic()
    data=json.loads(path.read_text())
    assert data["complete"] and data["theory"]=="GL(2)+1F"
    assert data["r"]=="1/2" and data["CS"]=="0" and data["W"]=="0"
    for relative,digest in data["implementation_sha256"].items():
        assert sha256((production.ROOT/relative).read_bytes()).hexdigest()==digest,relative
    complex_=production.Complex(data["Jw"])
    operator=FullOperator(complex_)
    sectors={(s["Rw"],s["antifields"]):s for s in data["sectors"]}
    expected=dimensions(complex_)
    assert expected=={key:s["dimension"] for key,s in sectors.items()}
    bases={key:[tuple(map(tuple,m)) for m in s["basis"]] for key,s in sectors.items()}
    indices={key:{m:i for i,m in enumerate(b)} for key,b in bases.items()}
    local={2:set(),3:set()}; all_columns={}; reps={}; tested=0
    for (R,af),sector in sorted(sectors.items()):
        b=bases[(R,af)]
        assert len(set(b))==len(b)==expected[(R,af)]
        target=indices.get((R+4,af),{})
        cols=[]
        for j,m in enumerate(b):
            assert complex_.algebra.canon(m)==(1,m)
            assert sum(complex_.jw(x) for x in m)==data["Jw"]
            assert sum(complex_.rw(x) for x in m)==R
            assert sum(complex_.af(x) for x in m)==af
            assert tuple(sum(complex_.charges(x)[i] for x in m) for i in range(2))==(0,0)
            image=complex_.engine.Q0_mono(m)
            vector={target[n]:production.exact(c) for n,c in image.items()}
            recorded={i:fmpq(c) for i,c in sector["Q0_columns"][j]}
            assert vector==recorded,((R,af),j)
            cols.append(vector); tested+=1
            for size in (2,3):
                for selected in combinations(range(len(m)),size):
                    part=tuple(m[i] for i in selected)
                    if sum(complex_.af(x) for x in part)>=2:local[size].add(part)
        assert rank(cols)==sector["Q0_rank"]
        all_columns[(R,af)]=cols
        reps[(R,af)]=[polynomial(r["polynomial"]) for r in sector["representatives"]]
        print("verified Q0",R,af,len(b),"seconds",round(time.monotonic()-started,1),flush=True)
    square_checks=0
    for (R,af),cols in all_columns.items():
        next_columns=all_columns.get((R+4,af),[])
        for column in cols:
            residual={}
            for i,c in column.items():add(residual,next_columns[i],c)
            assert not residual,(R,af,column)
            square_checks+=1
    images={}; maps={}; nreps=0
    for (R,af),sector in sorted(sectors.items()):
        previous=all_columns.get((R-4,af),[])
        b_rank=rank(previous)
        assert b_rank==sector["boundary_rank"]
        vectors=[{indices[(R,af)][m]:c for m,c in p.items()} for p in reps[(R,af)]]
        assert sector["cohomology_dimension"]==sector["dimension"]-sector["Q0_rank"]-b_rank
        assert rank(previous+vectors)==b_rank+len(vectors)
        target_key=(R+4,af-2); map_columns=[]
        for j,record in enumerate(sector["representatives"]):
            p=reps[(R,af)][j]
            assert not complex_.q0(p)
            actual=operator.q2(p)
            assert actual==polynomial(record["Q2_image"]),((R,af),j)
            assert not complex_.q0(actual)
            h={i:fmpq(c) for i,c in record["cohomology_coordinates"]}
            boundary={i:fmpq(c) for i,c in record["boundary_preimage_coordinates"]}
            witness={}
            for i,c in h.items():add(witness,reps[target_key][i],c)
            previous_target=(R,af-2)
            for i,c in boundary.items():
                for k,v in all_columns[previous_target][i].items():
                    add(witness,{bases[target_key][k]:v},c)
            assert witness==actual
            map_columns.append(h); nreps+=1
        assert rank(map_columns)==sector["d2_rank"]
        maps[(R,af)]=map_columns
    for (R,af),columns in maps.items():
        next_=maps.get((R+4,af-2),[])
        for column in columns:
            square={}
            for i,c in column.items():add(square,next_[i],c)
            assert not square
    nonzero=0; checked=0
    for size in (2,3):
        for m in sorted(local[size]):
            left=complex_.q0(operator.q2_monomial(m))
            right=operator.q2(complex_.q0_monomial(m))
            residual=dict(left); add(residual,right)
            assert not residual,(m,residual)
            nonzero+=bool(left or right); checked+=1
            if checked%100==0:
                print("local commutator",checked,"of",sum(map(len,local.values())),
                      "seconds",round(time.monotonic()-started,1),flush=True)
    report={"schema":"gl2-cohomology-independent-replay-v1","result":"passed",
            "certificate_sha256":sha256(path.read_bytes()).hexdigest(),
            "complete_cochain_basis":tested,"complete_cohomology_basis":nreps,
            "Q0_columns_recomputed_by_original_Wick_engine":tested,
            "Q0_squared_columns_independently_checked":square_checks,
            "multigraded_product_counts_equal_complete_bases":True,
            "Q2_images_and_boundary_witnesses_checked":nreps,
            "local_commutator_inputs":{str(k):len(v) for k,v in local.items()},
            "nonzero_commutator_compositions_cancelled":nonzero,
            "full_commutator_justification":"Q0 has differential order one; Q2 has orders two and three. Their graded commutator has order at most three. It vanishes on constants and generators and on all occurring two/three-letter submonomials. Thus it vanishes on the complete declared complex.",
            "d2_square_zero":True,"seconds":round(time.monotonic()-started,2),
            "verifier_sha256":sha256(Path(__file__).read_bytes()).hexdigest()}
    output=path.with_name(path.stem+"_verification.json")
    output.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print(json.dumps(report,indent=2))


if __name__=="__main__":
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else HERE/"gl2_jw20.json")
