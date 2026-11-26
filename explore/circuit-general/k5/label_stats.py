"""Label (constant) content of the K5 tables per kind: which subsets of {1, z3, z5} occur."""
import json, sys, os
from fractions import Fraction as F
HERE = os.path.dirname(os.path.abspath(__file__))
for q in map(int, sys.argv[1:]):
    d = json.load(open(os.path.join(HERE, 'tables', 'K5_q%d.json' % q)))
    st = {}
    for a, kd, v in d['values']:
        key = (kd, tuple(sorted(v, key=lambda s: (len(s), s))))
        st[key] = st.get(key, 0) + 1
    tot = {}
    for (kd, labs), c in st.items():
        tot[kd] = tot.get(kd, 0) + c
    print('q=%d:' % q, ', '.join('%s %s: %d' % (kd, '+'.join(labs) if labs else '0', c) for (kd, labs), c in sorted(st.items())))
