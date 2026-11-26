"""Bounded exact regression of the K=0 fixed-site kernels; C++17 and SymPy."""
import itertools
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile

from sympy import QQ
from sympy.polys.rings import ring

HERE = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(HERE.parents[1] / 'verify'))
from vcommon import PQR, theta_graph
from v_integrand import grouped_integrand

DRIVER = r'''
#include "multi_gamma.hpp"
#include <array>
#include <iostream>
#include <set>

void check_supports(const twoloop::detail::Graph& graph) {
    // Independently orient each edge and read the external fields from the indegrees.
    std::array<std::set<int>,1024> compatible;
    for (int mask=0;mask<8192;++mask) {
        std::array<int,12> incoming{};
        for (int e=0;e<13;++e) ++incoming[(mask>>e&1)?graph.head[e]:graph.tail[e]];
        if (incoming[0]%2!=1 || incoming[11]%2!=1) continue;
        int sites=0,cubics=(incoming[0]==3)+(incoming[11]==3);
        for(int v=1;v<=10;++v){
            if(incoming[v]!=1) sites|=1<<(v-1);
            cubics+=(incoming[v]==2);
        }
        if(2*cubics!=1+__builtin_popcount(sites)) throw std::logic_error("coupling order");
        compatible[sites].insert(mask);
    }
    for (unsigned sites=0;sites<1024;++sites) {
        if (__builtin_popcount(sites)%2==0) {
            if (!compatible[sites].empty()) throw std::logic_error("even support");
            continue;
        }
        std::set<int> generated;
        for(unsigned start:{0U,3U,5U,6U}) {
            auto flow=twoloop::detail::multi_directions(graph,twoloop::GammaSites(sites),start);
            int mask=0;for(int e=0;e<13;++e) if(flow.forward[e]) mask|=1<<e;
            generated.insert(mask);
        }
        if(generated.size()!=4 || generated!=compatible[sites]) throw std::logic_error("flow mismatch");
    }
}
int main(){
    twoloop::Lengths lengths; twoloop::Slots slots; unsigned mask;int check;
    if(!(std::cin>>lengths[0]>>lengths[1]>>lengths[2]>>mask>>check))return 2;
    for(auto& leg:slots)std::cin>>leg.momentum>>leg.vector_axis;
    if(check)check_supports(twoloop::detail::Graph(lengths,slots,0));
    auto p=twoloop::multi_gamma_numerator(lengths,slots,twoloop::GammaSites(mask));
    if(__builtin_popcount(mask)==1 && p!=twoloop::numerator(lengths,slots,__builtin_ctz(mask)))return 1;
    for(auto [key,value]:p)std::cout<<key<<' '<<value<<'\n';
}
'''


def momenta(lengths, rng):
    paths = [[], [], []]
    axes = list(range(1, 6))
    rng.shuffle(axes)
    for i, j in itertools.combinations(range(3), 2):
        for _ in range(6 - lengths[3-i-j]):
            signed_axis = axes.pop() * rng.choice((-1, 1))
            paths[i].append(signed_axis)
            paths[j].append(-signed_axis)
    for path in paths:
        rng.shuffle(path)
    return [p for path in paths for p in path]


def orientations(nv, edges):
    grouped = {}
    for bits in itertools.product((False, True), repeat=len(edges)):
        flow = [(v, u) if reverse else (u, v) for (u, v), reverse in zip(edges, bits)]
        inc = [0]*nv
        for _, head in flow:
            inc[head] += 1
        if inc[0] not in (1, 3) or inc[-1] not in (1, 3):
            continue
        sites = sum(1 << (v-1) for v in range(1, nv-1) if inc[v] != 1)
        types = {v: 'c' if inc[v] == (3 if v in (0, nv-1) else 2) else 'b' for v in range(nv)}
        grouped.setdefault(sites, []).append((flow, types))
    return grouped


def oracle(nv, edges, paths, flows, momenta, polarizations):
    _, *vars_ = ring('b1,b2,b3,b4,b5,c1,c2,c3,c4,c5', QQ)
    b, c = vars_[:5], vars_[5:]
    current = []
    zero = [0]*5
    ext = {v: dict(lam=zero, u=zero, vg=zero) for v in (0, nv-1)}
    for (vertices, _), y in zip(paths, ([-x-z for x, z in zip(b, c)], b, c)):
        current.append(y)
        for v in vertices:
            m, pol = momenta[v-1], polarizations[v-1]
            lam = [(1 if m>0 else -1)*int(i == abs(m)-1) for i in range(5)]
            axis = [int(i == pol) for i in range(5)]
            ext[v] = dict(lam=lam, u=axis, vg=axis)
            y = [x+k for x, k in zip(y, lam)]
            current.append(y)
    value = sum(grouped_integrand(nv, edges, paths, f, t, current, ext) for f, t in flows)
    return {sum(k << (4*i) for i, k in enumerate(mon)): int(coeff) for mon, coeff in value.items()}


def generated(exe, lengths, sites, momenta, polarizations, check):
    args = [*lengths, sites, int(check), *itertools.chain.from_iterable(zip(momenta, polarizations))]
    process = subprocess.run([exe], input=' '.join(map(str, args))+'\n', text=True, check=True, capture_output=True)
    return dict(tuple(map(int, line.split())) for line in process.stdout.splitlines())


def main():
    rng = random.Random(271608)
    checked = 0
    with tempfile.TemporaryDirectory() as directory:
        exe = Path(directory)/'multi_gamma'
        subprocess.run([os.environ.get('CXX','c++'), '-std=c++17', '-O2', '-fsanitize=signed-integer-overflow',
                        '-fno-sanitize-recover=all', '-I', HERE, '-x','c++','-','-o', exe],input=DRIVER,text=True,check=True)
        for lengths in PQR:
            nv, edges, paths = theta_graph(*lengths)
            grouped = orientations(nv, edges)
            for m in range(1, 6):
                for attempt in range(500):
                    sites = sum(1 << j for j in rng.sample(range(10), 2*m-1))
                    moment = momenta(lengths, rng)
                    polar = [rng.choice([i for i in range(5) if i!=abs(q)-1]) for q in moment]
                    poly = generated(exe, lengths, sites, moment, polar, m==1 and attempt==0)
                    if poly:
                        break
                else:
                    raise AssertionError(('missing nonzero control', lengths, m))
                assert len(grouped[sites])==4
                assert oracle(nv, edges, paths, grouped[sites], moment, polar)==poly, (lengths,m,sites)
                # The exact gauge polarization on any gamma slot annihilates its vertex.
                v = next(i for i in range(10) if sites & (1<<i))
                polar[v] = abs(moment[v])-1
                assert not generated(exe, lengths, sites, moment, polar, False)
                checked+=1
                print(lengths,'m',m,'support',hex(sites),len(poly),'terms',flush=True)
    print(checked, 'nonzero polynomial identities, 25 longitudinal zeros, 2560 support sets exhaustively checked.')

if __name__=='__main__':
    main()
