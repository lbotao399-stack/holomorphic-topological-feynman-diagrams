import sys, time
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
from holography import bg_words as bw, bg_complex_k as bk, bg_complex as bc, bg_sizes as bs
alg, keys = bw.n4_letters(12)
# cyclic counts vs enumeration
for lam in [(2,1,0,0,0),(1,1,1,0,0),(2,1,1,0,0),(1,1,0,1,0),(2,2,0,0,0),(1,0,0,2,1)]:
    cc = bs.cyclic_class_counts(alg, lam)
    enum = {}
    for w in alg.words_upto(lam):
        s, cw = alg.canon_cyclic(w)
        if s and cw == w:
            pass
    cyc = set()
    for w in alg.words_upto(lam):
        s, cw = alg.canon_cyclic(w)
        if s: cyc.add(cw)
    for w in cyc:
        key = (alg.wword(w), len(w))
        enum[key] = enum.get(key, 0) + 1
    print(lam, 'cyclic ok' if enum == cc else ('cyclic MISMATCH', enum, cc))
# chain dims vs enumeration, k = 2, N = 1, 2
for N in (1, 2):
    bg = bk.BGk(alg, 2)
    for lam in [(1,0,0,0,0),(2,0,0,0,0),(1,1,0,0,0),(1,0,0,1,0),(2,1,0,0,0),(1,1,1,0,0),(0,0,0,1,1),(3,0,0,0,0)]:
        t0 = time.time()
        ch = bg.chains(lam, (2*N, 2*N+1))
        Tch = (N, N, N, N)
        ghs = sorted(set(g for (c, g) in ch))
        extra = [((N+1, N-1), (N, N)), ((N, N), (N+1, N-1)), ((N+1, N-1), (N+1, N-1)), ((N+1, N), (N, N+1)), ((N, N+1), (N+1, N))]
        d = bs.chain_dims(alg, lam, 2, ghs, N, extra_charges=extra)
        bad = 0
        for (Lc, Rc, g), v in d.items():
            key = (tuple(Lc) + tuple(Rc), g)
            e = len(ch.get(key, []))
            if e != round(v):
                bad += 1
                print('  mismatch', lam, Lc, Rc, g, e, v)
        print('N', N, lam, 'k2 dims', 'ok' if not bad else 'BAD', {g: d[((N,N),(N,N),g)] for g in ghs}, '%.1fs' % (time.time()-t0), flush=True)
# k = 1
bg1 = bc.BG1(alg)
for N in (1, 2):
    for lam in [(2,1,0,0,0),(1,1,1,0,0),(2,2,0,0,0),(2,1,1,0,0),(1,1,0,1,0)]:
        ch = bg1.chains(lam, N)
        ghs = sorted(ch)
        d = bs.chain_dims(alg, lam, 1, ghs, N)
        ok = all(len(ch[g]) == round(d[((N,),(N,),g)]) for g in ghs)
        print('k1 N', N, lam, 'ok' if ok else ('BAD', {g: len(ch[g]) for g in ghs}, d))
