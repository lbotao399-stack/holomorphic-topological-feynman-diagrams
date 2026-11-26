"""Integrable degree-fifteen prism monomials, cut in twelve dimensions."""
from collections import defaultdict
from functools import cache
from itertools import product
from pathlib import Path
import re
import subprocess
import sys

import networkx as nx

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
from graphs import PRISM

GRAPH = nx.Graph(PRISM)
INDEX = {frozenset(e): i for i, e in enumerate(PRISM)}
MAPS = tuple(nx.isomorphism.GraphMatcher(GRAPH, GRAPH).isomorphisms_iter())
PERMS = tuple(tuple(INDEX[frozenset((g[u], g[v]))] for u, v in PRISM) for g in MAPS)
CUTS = (0, 1)
BATCH = 128


@cache
def margins():
    rows = []
    for mask in range(1, 511):
        subset = tuple(e for e in range(9) if mask & (1 << e))
        g = nx.Graph(PRISM[e] for e in subset)
        loops = len(subset) - len(g) + nx.number_connected_components(g)
        if 6*loops >= len(subset):
            rows.append((subset, 6*loops - len(subset)))
    return rows


@cache
def support():
    seen = set()
    keep = []
    for ex in product(range(4), repeat=9):
        if sum(ex) != 15 or ex in seen:
            continue
        seen.update(tuple(ex[e] for e in p) for p in PERMS)
        if all(sum(ex[e] for e in subset) > threshold for subset, threshold in margins()):
            keep.append(ex)
    return keep


def edge_order(cut):
    return (cut,) + tuple(i for i in range(9) if i != cut)


def vertices(cut, edges=PRISM):
    incident = [[] for _ in range(6)]
    for i, k in enumerate(edge_order(cut)):
        a, b = edges[k]
        label = f'p{i}' if i else 'Q'
        incident[a].append(label)
        incident[b].append('-' + label)
    return '*'.join('vx(' + ','.join(row) + ')' for row in incident)


def form(folder, executable, name, expression, cut=0, option='', edges=PRISM):
    job = folder/(name+'.frm')
    job.write_text('#-\n#include- forcer.h\nS tag;\nL F=' + vertices(cut,edges) + '*(' + expression
                   + ');\n#call Forcer(' + option + ')\n#write <' + name + '.out> "%E",F\n.end\n')
    with (folder/(name+'.log')).open('w') as log:
        subprocess.run([str(executable), '-q', job.name], cwd=folder, stdout=log,
                       stderr=subprocess.STDOUT, check=True)
    return (folder/(name+'.out')).read_text()


def tagged(terms):
    return '\n+'.join(f'tag^{i+2}*({term})' for i, term in terms)


def read_vectors(text):
    text = ''.join(text.split()).replace('\\', '')
    terms = list(re.finditer(r'([^+]+?)\*tag\^(\d+)\*rat\(([^,]+),([^()]+)\)', text))
    assert '+'.join(m.group(0) for m in terms) == text
    rows = defaultdict(dict)
    for match in terms:
        master, tag, top, bottom = match.groups()
        vector = rows[int(tag)-2]
        assert master not in vector
        vector[master] = (top, bottom)
    return rows


def monomial_jobs(folder, executable):
    exponents = support()
    for start in range(0, len(exponents), BATCH):
        for cut in CUTS:
            order = edge_order(cut)
            terms = []
            for i, ex in enumerate(exponents[start:start+BATCH], start):
                powers = max(tuple(ex[p[e]] for e in order) for p in PERMS)
                denominator = '*'.join(f'1/(p{j}.p{j})^{k+1}' for j,k in enumerate(powers) if j)
                terms.append((i, denominator))
            form(folder, executable, f'period{start}_{cut}', tagged(terms), cut)
        print('two cuts', min(start+BATCH, len(exponents)), '/', len(exponents), flush=True)
