import sys, time
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
from holography import bg_words as bw, bg_complex_k as bk
which = sys.argv[1]; N = int(sys.argv[2]); k = int(sys.argv[3]); nmax = int(sys.argv[4])
alg = bw.ctheta_letters()[0] if which == 'theta' else bw.cz_letters(nmax)[0]
bg = bk.BGk(alg, k)
for n in range(0, nmax+1):
    t=time.time()
    res = bk.homology(bg, (n,), N)
    hs = {g: v[3] for g, v in res.items() if v[3]}
    dims = sum(v[0] for v in res.values())
    print('x^%d'%n, 'H_g:', hs, 'dimV_T', dims, '%.2fs'%(time.time()-t), flush=True)
