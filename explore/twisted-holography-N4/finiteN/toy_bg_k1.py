import sys, time
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
from holography import bg_words as bw, bg_complex as bc
which = sys.argv[1]; N = int(sys.argv[2]); nmax = int(sys.argv[3])
if which == 'theta':
    alg, _ = bw.ctheta_letters()
else:
    alg, _ = bw.cz_letters(nmax)
bg = bc.BG1(alg)
# checks
for n in range(1, min(nmax, 6)+1):
    for ber in (N, N+1):
        t, b = bc.check_d2(bg, (n,), 3, ber, absolute=False)
        if b: print('d2 fail rel', n, ber, b, t)
# absolute check with zero modes: build monomials including P0, MM0, R0 by hand
import itertools
P0=(bc.P,()); MM0=(bc.MM,()); R0=(bc.R,()); M0=(bc.M,())
def dd(m, Nn):
    v = bg.d(list(m), Nn); acc={}
    for k,c in v.items():
        for k2,c2 in bg.d(list(k), Nn).items():
            acc[k2]=acc.get(k2,0)+c*c2
    return {k:x for k,x in acc.items() if x}
ch = bg.chains((3,), 1)
cnt=0; bad=0
for g, monos in ch.items():
    for m in monos[:50]:
        for extra in ([P0],[MM0],[R0],[P0,R0],[MM0,M0],[P0,MM0]):
            sg, key = bg.canon(list(m)+extra)
            if not sg: continue
            cnt+=1
            if dd(key, 5): bad+=1
print('absolute d^2 tests', cnt, 'bad', bad)
# D_r anticommutes with d on C'_{N+1}
cnt=0; bad=0
for n in range(1, min(nmax,5)+1):
    ch = bg.chains((n,), N+1)
    for g, monos in ch.items():
        for m in monos[:100]:
            a = {}
            for k,c in bg.Dr(m).items():
                for k2,c2 in bg.d(k, N).items(): a[k2]=a.get(k2,0)+c*c2
            for k,c in bg.d(m, N).items():
                if any(bg.is_zero_mode_b(x) for x in k): continue
                for k2,c2 in bg.Dr(k).items(): a[k2]=a.get(k2,0)+c*c2
            a={k:x for k,x in a.items() if x}
            cnt+=1
            if a: bad+=1
print('dDr+Drd tests', cnt, 'bad', bad)
for n in range(0, nmax+1):
    t=time.time()
    res = bc.homology(bg, (n,), N)
    hs = {g: v[3] for g, v in res.items() if v[3]}
    chi = sum((-1)**(g % 2) * h for g, h in hs.items())
    print('x^%d'%n, 'H_g:', hs, 'chi', chi, '%.2fs'%(time.time()-t), flush=True)
