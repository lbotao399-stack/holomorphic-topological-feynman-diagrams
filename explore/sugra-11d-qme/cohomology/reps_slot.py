"""Representatives of H(C_quot)^{w=1} on one slot (n, t, N) and the w=0 class of their beta_0-coefficient.
usage: python3 reps_slot.py n N t1 t2 ... -> appends to reps_slot.log"""
import sys
from cohom import Space, Cohom
from jets import Echelon, vec
from reps import reps_C

n = int(sys.argv[1]); N = int(sys.argv[2]); t = tuple(int(x) for x in sys.argv[3:3 + n])
G = -1 - sum(t)
out = open('reps_slot.log', 'a')


def log(s):
    print(s, flush=True)
    out.write(s + '\n')
    out.flush()


S1 = Space(n, 1, G + n + 2)
S0 = Space(n, 0, G + n + 3)
C1 = Cohom(S1)
A, A0 = S1.A, S0.A
h = C1.HC(t)
log('n=%d w=1 G=%d t=%s : H(C_quot) = %s' % (n, G, t, {N_: v[1] for N_, v in h.items() if v[1]}))
zs = reps_C(C1, t, N)
b0 = A.L('b', 0, (0,) * n)
for z in zs:
    nb = sum(1 for m in z if b0 in m)
    types = {}
    for m in z:
        key = tuple(sorted(A.info[x][0] + ('' if A.info[x][0] == 'X' else str(sum(A.info[x][2]))) for x in m if A.info[x][0] != 'X'))
        types[key] = types.get(key, 0) + 1
    log(' N=%d representative: %d terms; terms by non-X content: %s' % (N, len(z), types))
    log('   z = ' + A.pstr(z, 12))
    c = {}
    for m, v in z.items():
        if b0 in m:
            j = m.index(b0)
            sg = -1 if sum(1 for x in m[:j] if A.par[x]) % 2 else 1
            r = m[:j] + m[j + 1:]
            c[r] = c.get(r, 0) + sg * v
    c = {m: v for m, v in c.items() if v}
    if not c:
        log('   beta_0-coefficient: 0')
        continue
    c0 = {}
    for m, v in c.items():
        ids = [A0.L(*A.info[x]) for x in m]
        odd = [x for x in ids if A0.par[x]]
        inv = sum(1 for i in range(len(odd)) for j in range(i + 1, len(odd)) if odd[i] > odd[j])
        mm = tuple(sorted(ids))
        c0[mm] = c0.get(mm, 0) + (-1 if inv % 2 else 1) * v
    closed = not A0.Q(c0)
    src = S0.slot(t).get(N - 2, [])
    mid = S0.slot(t).get(N - 1, [])
    idx = {m: j for j, m in enumerate(mid)}
    E = Echelon()
    for m in src:
        E.insert(vec(A0.Qmono(m), idx))
    r0 = E.rank
    E.insert(vec(c0, idx))
    log('   beta_0-coefficient: %d terms, Q c = 0: %s, c exact in C^{w=0}: %s' % (len(c0), closed, E.rank == r0))
