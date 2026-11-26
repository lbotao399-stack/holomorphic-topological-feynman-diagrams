"""Reconstruct the entire labelled core numerator and integrate its two period channels."""
from fractions import Fraction as F
from itertools import combinations_with_replacement
from math import ceil, comb, factorial, gcd, lcm, prod
from pathlib import Path
from tempfile import NamedTemporaryFile
import json

import numpy as np
from flint import fmpq, fmpz, fmpq_mpoly_ctx, nmod_mat, nmod_mpoly_ctx
from core_cubature import evaluations, kirchhoff, rational_table, supports
from fourfold_bulk import CORES, symmetries
from ordered_period_bounds import determinant_masses
from trace_reconstruction import load_rule


def position_denominator():
    # The squared cofactor has total degree ten and single-position degree
    # four. The interpolated voltage adds two degrees on distinct slots.
    denominator=1
    for degrees in combinations_with_replacement(range(6),10):
        if sum(degrees)<=12 and degrees.count(5)<=2:
            denominator=lcm(denominator,prod(a+1 for a in degrees))
    for a in range(5):
        for b in range(5):
            bridge=(F(1,a+2)+F(1,b+2))/(a+b+3)
            for rest in combinations_with_replacement(range(5),8):
                if a+b+sum(rest)<=10:
                    denominator=lcm(denominator,(bridge/prod(k+1 for k in rest)).denominator)
    return denominator


def coefficient_bounds(name,source,jumps,probes):
    orbits=supports(CORES[name],12)
    powers=np.array([k for orbit in orbits for k in orbit],dtype=np.int8)
    indices=np.cumsum([0]+[len(orbit) for orbit in orbits[:-1]])
    path_den=[2*prod(factorial(k) for k in orbit[0]) for orbit in orbits]
    cube_den=position_denominator()
    denominators=[gcd(cube_den,d) for d in path_den]
    masses=determinant_masses(name)
    common=lcm(*(mass.denominator for mass in masses.values()))
    total=np.zeros(len(orbits),dtype=object)
    local=np.array([factorial(5)**(a//5)*factorial(a%5) for a in range(11)],dtype=np.int64)
    for counts,mass in masses.items():
        difference=powers-np.asarray(counts,dtype=np.int8)
        valid=np.all((difference>=0)&(difference<=5*(np.asarray(counts)+1)),axis=1)
        factors=np.zeros(len(powers),dtype=np.int64)
        factors[valid]=np.prod(local[difference[valid]],axis=1)
        maximum=np.maximum.reduceat(factors,indices).astype(object)
        total+=(common*mass).numerator*maximum
    halves=[F(sum(max(v,0) for v in axis),2) for axis in zip(*source)]
    scale=prod(v**6 for v in halves)*prod(sum((x/h)**2 for x,h in zip(row,halves)) for row in jumps)
    content=prod(gcd(*row) for row in jumps)**2
    ceilings=[]
    for left,right,axis in probes:
        if source[left][axis] or source[right][axis]:raise ValueError('Expected a transverse probe.')
        weight=sum(abs(p[axis]) for p in source)
        factor=2*scale*comb(16,4)*weight/(common*content)
        ceilings.append([ceil(factor*m*d/q) for m,d,q in zip(total,denominators,path_den)])
    return denominators,content,ceilings


def numerator_rows(cache,name,fixture,p,probes,orbits):
    rule=load_rule(cache/f'core-rule-{name}-{p}.tsv',name,p)
    values={}
    for line in (cache/f'primary-{name}-{fixture}-{p}.tsv').read_text().splitlines():
        row=tuple(map(int,line.split()))
        if len(row)!=6+len(probes) or row[:6] in values or not all(0<=v<p for v in row[6:]):
            raise ValueError('Invalid numerator checkpoint.')
        values[row[:6]]=row[6:]
    nodes=[r[:6] for r in rule]
    if set(nodes)!=set(values):raise ValueError('Incomplete numerator checkpoint.')
    matrix=evaluations(nodes,orbits,12,p).transpose()
    rhs=nmod_mat([list(values[t]) for t in nodes],p)
    result=matrix.solve(rhs)
    assert matrix*result==rhs
    return result.tolist()


def reduced_periods(core,power,quotient):
    orbits=supports(core,power)
    coefficients=quotient.to_dict()
    assert not set(coefficients).difference(e for orbit in orbits for e in orbit)
    moments=rational_table(core,power)
    values=[F(0),F(0)]
    for orbit,row in zip(orbits,moments):
        c=F(str(coefficients.get(orbit[0],0)))
        for k in (0,1):values[k]+=c*row[k]
    return [x/len(symmetries(core)) for x in values]


def common_core_factor(cache,name):
    core=CORES[name]
    orbits=supports(core,12)
    ctx=fmpq_mpoly_ctx.get(tuple(f't{i}' for i in range(6)))
    modctx=nmod_mpoly_ctx.get(tuple(f't{i}' for i in range(6)),1009)
    t=ctx.gens()
    u=sum(prod(t[i] for i in cut) for cut in kirchhoff(core))
    extra=ctx.constant(1)
    if name=='double-square':extra=(t[0]*t[1]*t[4]*t[5])**2*(t[0]+t[1])*(t[4]+t[5])
    residual_gcd=None
    for fixture in ('four','three','soft-two'):
        record=json.loads((cache/f'primary-{name}-{fixture}-polynomial.json').read_text())
        for row,period in zip(record['orbit_coefficients'],record['values']):
            coefficients=list(map(F,row))
            poly=ctx.from_dict({e:fmpq(v.numerator,v.denominator)
                                for v,orbit in zip(coefficients,orbits) if v for e in orbit})
            reduced,remainder=divmod(poly,u**2)
            assert not remainder and reduced%u
            assert reduced_periods(core,10,reduced)==list(map(F,period))
            quotient,remainder=divmod(reduced,extra)
            assert not remainder
            if residual_gcd is None or residual_gcd.total_degree()>0:
                residual=modctx.from_dict({tuple(e):int(v.numerator)*pow(int(v.denominator),-1,1009)%1009
                                           for e,v in quotient.terms()})
                residual_gcd=residual if residual_gcd is None else residual_gcd.gcd(residual)
    assert residual_gcd.total_degree()==0
    return dict(kirchhoff_power=2,extra_factor='(t0*t1*t4*t5)^2*(t0+t1)*(t4+t5)' if name=='double-square' else '1',
                reduced_period_denominator_power=10,residual_gcd_prime=1009)


def reconstruct_polynomials(cache,name,fixture,source,jumps,probes):
    core=CORES[name]
    orbits=supports(core,12)
    den,content,bounds=coefficient_bounds(name,source,jumps,probes)
    nums=[[0]*len(probes) for _ in orbits]
    modulus=1
    candidate=999983
    while modulus<=2*max(map(max,bounds)):
        while not fmpz(candidate).is_prime():candidate-=2
        p=candidate
        candidate-=2
        rows=numerator_rows(cache,name,fixture,p,probes,orbits)
        modulus_inverse=pow(modulus,-1,p)
        content_inverse=pow(content,-1,p)
        for i,row in enumerate(rows):
            for j,v in enumerate(row):
                lift=(int(v)*den[i]*content_inverse-nums[i][j])*modulus_inverse%p
                nums[i][j]+=modulus*lift
        modulus*=p
        print(name,fixture,'polynomial CRT bits',modulus.bit_length(),flush=True)
    polys=[[] for _ in probes]
    for i,row in enumerate(nums):
        for j,x in enumerate(row):
            v=x-modulus if 2*x>modulus else x
            assert abs(v)<=bounds[j][i]
            polys[j].append(F(v*content,den[i]))
    # Check every reconstructed coefficient against an independently sampled prime.
    held=numerator_rows(cache,name,fixture,1009,probes,orbits)
    for i,row in enumerate(held):
        for j,x in enumerate(row):
            v=polys[j][i]
            assert v.numerator*pow(v.denominator,-1,1009)%1009==int(x)

    periods={tuple(map(int,r[:6])):tuple(map(F,r[6:])) for r in
             (line.split() for line in (cache/f'core-periods-{name}.tsv').read_text().splitlines())}
    totals=[]
    factors=[]
    ctx=fmpq_mpoly_ctx.get(tuple(f't{i}' for i in range(6)))
    u=sum(prod(ctx.gen(i) for i in cut) for cut in kirchhoff(core))
    for poly in polys:
        totals.append([sum(x*len(o)*periods[o[0]][k] for x,o in zip(poly,orbits))/len(symmetries(core))
                       for k in (0,1)])
        flat=ctx.from_dict({e:fmpq(value.numerator,value.denominator)
                           for orbit,value in zip(orbits,poly) if value for e in orbit})
        mult=0
        while flat:
            quotient,remainder=divmod(flat,u)
            if remainder:break
            mult+=1
            flat=quotient
        factors.append(mult)
        if mult:
            assert reduced_periods(core,12-mult,flat)==totals[-1]
    return dict(core=name,fixture=fixture,probes=probes,modulus=str(modulus),
                verification_prime=1009,content=content,position_denominator=position_denominator(),
                coefficient_bounds=[max(r) for r in bounds],kirchhoff_multiplicities=factors,
                orbit_coefficients=[[str(x) for x in poly] for poly in polys],
                values=[[str(v) for v in row] for row in totals])


def fixture_data(fixture):
    from trace_reconstruction import TESTS,waves
    if fixture=='soft-two':
        from soft_trace_measurement import momenta
        source=momenta(0)
        jumps=source[:]
        jumps[7]=(0,0,1,1,0)
        probes=((2,3,0),)
    else:
        source=waves()[fixture]
        jumps=source
        probes=TESTS[fixture]
    return source,jumps,probes


def verify_periods(cache,name,fixture,record):
    source,jumps,probes=fixture_data(fixture)
    den,content,bounds=coefficient_bounds(name,source,jumps,probes)
    assert record['position_denominator']==position_denominator() and record['content']==content
    assert record['coefficient_bounds']==[max(row) for row in bounds]
    assert int(record['modulus'])>2*max(record['coefficient_bounds'])
    core=CORES[name]
    orbits=supports(core,12)
    periods={tuple(map(int,r[:6])):tuple(map(F,r[6:])) for r in
             (line.split() for line in (cache/f'core-periods-{name}.tsv').read_text().splitlines())}
    assert len(record['orbit_coefficients'])==len(probes)
    for j,row in enumerate(record['orbit_coefficients']):
        coefficients=list(map(F,row))
        assert len(coefficients)==len(orbits)
        for x,d,b in zip(coefficients,den,bounds[j]):
            v=x*d/content
            assert v.denominator==1 and abs(v)<=b
        values=[sum(x*len(o)*periods[o[0]][k] for x,o in zip(coefficients,orbits))/len(symmetries(core))
                for k in (0,1)]
        assert [str(x) for x in values]==record['values'][j]


def run_fixture(cache,name,fixture):
    source,jumps,probes=fixture_data(fixture)
    record=reconstruct_polynomials(cache,name,fixture,source,jumps,probes)
    path=cache/f'primary-{name}-{fixture}-polynomial.json'
    with NamedTemporaryFile('w',dir=cache,delete=False) as out:
        json.dump(record,out,indent=2)
        print(file=out)
        tmp=Path(out.name)
    tmp.replace(path)
    return record


if __name__=='__main__':
    from argparse import ArgumentParser
    cli=ArgumentParser()
    cli.add_argument('cache',type=Path)
    cli.add_argument('core',choices=('tetrahedron','double-square'))
    cli.add_argument('fixture',choices=('soft-two','three','four'))
    args=cli.parse_args()
    record=run_fixture(args.cache,args.core,args.fixture)
    print(json.dumps({k:v for k,v in record.items() if k!='orbit_coefficients'},indent=2))
