"""Oriented GKW bitriangle Gaussian pushforward, exact at all positive times.

The order of all differentials is part of the input.  The independent
polynomial identities compare the original position-space propagators
with the repository cycle-Pfaffian density, retaining all three sources.
"""
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path[:0] = [str(ROOT / "code"), str(ROOT / "explore/circuit-general/cores")]
from flint import fmpq_mpoly_ctx
from htk.graph import Graph
from htk import brute
import numer

EDGES = ((0,1),(1,2),(1,3),(2,0),(3,2))
NAMES = tuple(["t%d" % i for i in range(5)] + ["barz%d" % i for i in range(3)]
              + ["u%d" % i for i in range(3)] + ["lambda%d" % i for i in range(3)])


def wedge(a,b):
    result={}
    for mask,x in a.items():
        for other,y in b.items():
            if mask&other:
                continue
            inversions=sum((other & ((1<<i)-1)).bit_count()
                           for i in range(mask.bit_length()) if mask & (1<<i))
            key=mask|other
            value=result.get(key,0)+(-1)**inversions*x*y
            if value:
                result[key]=value
            else:
                result.pop(key,None)
    return result


def rational(x):
    if hasattr(x,"to_dict"):
        d=x.to_dict()
        assert all(not any(e) for e in d)
        return sum((rational(c) for c in d.values()),Fraction(0))
    if hasattr(x,"p"):
        return Fraction(int(x.p),int(x.q))
    return Fraction(x)


def terms(p):
    return [[list(map(int,e)),str(rational(c))] for e,c in sorted(p.to_dict().items())]


def derive():
    ctx=fmpq_mpoly_ctx.get(NAMES)
    v=ctx.gens(); times=v[:5]; bar=v[5:8]; top=v[8:11]; source=v[11:14]
    graph=Graph(4,list(EDGES)); incidence=graph.incidence(); cycles=graph.cycle_matrix()
    full_sources=((1,0,0),(0,1,0),(0,0,1),(-1,-1,-1))
    tree=numer.tree_current_polys(graph,full_sources)
    expand=lambda p:ctx.from_dict({tuple(e)+(0,)*9:c for e,c in p.items()})
    current=[[expand(tree[e][j]) for j in range(3)] for e in range(5)]
    matrix=[[sum(cycles[i][e]*cycles[j][e]*times[e] for e in range(5))
             for j in range(2)] for i in range(2)]
    kirchhoff=matrix[0][0]*matrix[1][1]-matrix[0][1]*matrix[1][0]
    adj=((matrix[1][1],-matrix[0][1]),(-matrix[1][0],matrix[0][0]))
    # Exact grounded voltage response: U phi_1=t_2(U y_2),
    # U phi_2=-t_4(U y_4), U phi_0=t_0(U y_0)+U phi_1.
    response=[[times[0]*current[0][j]+times[2]*current[2][j] for j in range(3)],
              [times[2]*current[2][j] for j in range(3)],
              [-times[4]*current[4][j] for j in range(3)]]
    assert all(response[i][j]==response[j][i] for i in range(3) for j in range(3))
    assert all(sum(incidence[e][i]*current[e][j] for e in range(5))
               == (kirchhoff if i==j else 0) for i in range(3) for j in range(3))
    potential=[sum(response[i][j]*source[j] for j in range(3)) for i in range(3)]
    y=[sum(current[e][j]*source[j] for j in range(3)) for e in range(5)]
    product_times=times[0]*times[1]*times[2]*times[3]*times[4]
    # Ordered coordinate differentials:
    # dbarz0,du0,dbarz1,du1,dbarz2,du2,dt0,...,dt4.
    # Propagator forms are F_e=t_e dbarz_e-barz_e dt_e and
    # G_e=t_e du_e-u_e dt_e/2, so dY_e ds_e=t_e^(-7/2) F_e G_e.
    forms=[]
    for e in range(5):
        f={1<<(2*i):incidence[e][i]*times[e] for i in range(3) if incidence[e][i]}
        g={1<<(2*i+1):incidence[e][i]*times[e] for i in range(3) if incidence[e][i]}
        f[1<<(6+e)]=-sum(incidence[e][i]*bar[i] for i in range(3))
        g[1<<(6+e)]=-sum(incidence[e][i]*top[i] for i in range(3))/2
        forms.extend((f,g))
    raw={0:ctx.constant(1)}
    for form in forms:
        raw=wedge(raw,form)
    lambda_product=source[0]*sum(source)
    position_controls=[]
    for omitted in range(5):
        mask=((1<<11)-1)^(1<<(6+omitted))
        polynomial=raw.get(mask,ctx.constant(0))
        averaged=ctx.constant(0)
        for exponent,coefficient in polynomial.to_dict().items():
            assert sum(exponent[5:8])==2 and sum(exponent[8:11])==2
            positions=[i for i in range(3) for _ in range(int(exponent[8+i]))]
            monomial=ctx.from_dict({tuple(exponent[:5])+(0,)*9:coefficient})
            for i in range(3):
                monomial*=potential[i]**int(exponent[5+i])
            monomial*=response[positions[0]][positions[1]]
            averaged+=monomial
        # E(raw)=averaged/(2 U^3): two barz means and one real covariance.
        expected=-(-1)**omitted*times[omitted]*product_times**2*kirchhoff**2*lambda_product/4
        assert averaged==expected, (omitted,averaged-expected)
        position_controls.append({"omitted_edge":omitted,
            "raw_wedge_coefficient":terms(polynomial),
            "gaussian_numerator":terms(averaged),
            "exact_identity":"-(-1)^e*t_e*(product t)^2*U^2*lambda0*(lambda0+lambda1+lambda2)/4"})
    # Separate current-space Pfaffian construction, before any Gaussian
    # pushforward or simplification of the position-space forms.
    dmatrix=[[{1<<e:ctx.constant(cycles[i][e]*cycles[j][e])
               for e in range(5) if cycles[i][e]*cycles[j][e]}
              for j in range(2)] for i in range(2)]
    top_form={}
    for a in range(2):
        for b in range(2):
            for mask,p in wedge(dmatrix[0][a],dmatrix[1][b]).items():
                value=top_form.get(mask,0)-adj[a][b]*p/2
                if value:
                    top_form[mask]=value
                else:
                    top_form.pop(mask,None)
    psi=[{1<<e:cycles[i][e]*y[e] for e in range(5) if cycles[i][e]}
         for i in range(2)]
    cycle_form=wedge(wedge(top_form,psi[0]),psi[1])
    cycle_controls=[]
    for omitted in range(5):
        coefficient=cycle_form.get(((1<<5)-1)^(1<<omitted),ctx.constant(0))
        expected=(-1)**omitted*times[omitted]*kirchhoff**2*lambda_product/2
        assert coefficient==expected, (omitted,coefficient-expected)
        cycle_controls.append({"omitted_edge":omitted,"coefficient":terms(coefficient)})
    unequal=[]
    for values in ((1,1,1,1,1),(2,3,5,7,11),(3,2,7,5,13)):
        for lam in ((1,0,0),(1,2,3),(-2,3,4)):
            density=brute.density(graph,values,[[a] for a in lam],T=1,check_basic=True)
            u=brute.symanzik_U(graph,values)
            expected=Fraction(lam[0]*sum(lam),2*u*u)
            assert density.coeff==expected and density.radicand==u
            unequal.append({"times":values,"sources":lam,"U":str(u),
                            "htk_coefficient":str(expected),"epsilon2":-1})
    # Zero-loop calibration with a single oriented edge: the two scaled
    # forms t*dbarz and t*du give t^2; the normalized Gaussian and
    # propagator factors multiply to t^(-2).  Thus I_0=K_0=1.
    data={"schema":"gkw-epsilon2-oriented-pushforward-v1","edges":EDGES,
          "root":3,"variables":NAMES,"epsilon2":-1,
          "primary_source":"https://arxiv.org/html/2403.13049v1",
          "primary_equations":{"spacetime_orientation":49,"full_ordered_integral":59,
                               "super_propagator":60,"printed_bitriangle":109},
          "differential_order":["dbarz0","du0","dbarz1","du1","dbarz2","du2",
                                "dt0","dt1","dt2","dt3","dt4"],
          "holomorphic_volume_order":["dz0/(2*pi*i)","dz1/(2*pi*i)","dz2/(2*pi*i)"],
          "holomorphic_volume_reordering_sign":1,
          "projective_form":"i_E(dt0 wedge dt1 wedge dt2 wedge dt3 wedge dt4)",
          "ordered_scaled_edge_forms":[[[mask.bit_length()-1,str(p)] for mask,p in sorted(f.items())]
                                       for f in forms],
          "kirchhoff":terms(kirchhoff),"grounded_response_numerators":[[terms(p) for p in row] for row in response],
          "position_space_cofactors":position_controls,"cycle_space_cofactors":cycle_controls,
          "gaussian_prefactor_after_scaled_forms":"1/(pi*(product t)^2*U^(3/2))",
          "gkw_density":"-lambda0*(lambda0+lambda1+lambda2)/(8*pi*U^(5/2))",
          "htk_density":"lambda0*(lambda0+lambda1+lambda2)/(2*U^(5/2))",
          "conversion":"I_GKW=-(4*pi)^(-1)*K_htk",
          "zero_loop_calibration":"I0=K0=1",
          "unequal_time_original_graph_controls":unequal}
    paths=[Path(__file__),ROOT/"code/htk/graph.py",ROOT/"code/htk/brute.py",ROOT/"code/htk/linalg.py",
           ROOT/"explore/circuit-general/cores/numer.py"]
    data["implementation_sha256"]={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths}
    return data


if __name__=="__main__":
    data=derive()
    path=HERE/"epsilon2_sign.json"
    path.write_text(json.dumps(data,indent=2,sort_keys=True)+"\n")
    print("epsilon2=-1: all five position and cycle cofactors agree identically; nine original-graph controls; zero-loop calibration +1.")
