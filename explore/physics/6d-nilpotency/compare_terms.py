"""Composition-by-composition comparison of two run_nilp logs (same seed): E_1 = s E_2 for s = +1 or -1.
Usage: python3 compare_terms.py LOG1 LOG2 [prime]"""
import re
import sys

P = int(sys.argv[3]) if len(sys.argv) > 3 else 2147483647


def read(fn):
    E = {}
    for line in open(fn):
        m = re.match(r'\s+<hatQ\[(\S+)\] hatQ\[(\S+)\] O>: \[([0-9, ]+)\]', line)
        if m:
            E[(m.group(1), m.group(2))] = [int(x) for x in m.group(3).split(',')]
    return E


E1, E2 = read(sys.argv[1]), read(sys.argv[2])
keys = sorted(set(E1) & set(E2))
plus = all(E1[k] == E2[k] for k in keys)
minus = all(all((a + b) % P == 0 for a, b in zip(E1[k], E2[k])) for k in keys)
nz = sum(1 for k in keys for x in E1[k] if x)
print('compositions compared: %d (nonzero values %d); E1 == +E2: %s; E1 == -E2: %s' % (len(keys), nz, plus, minus))
