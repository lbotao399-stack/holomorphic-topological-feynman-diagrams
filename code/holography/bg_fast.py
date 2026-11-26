"""Fast homology of the Budzik-Gaiotto k-meson complex (conventions of bg_complex_k / bg_complex).

The C++ kernel native/bgk.cpp enumerates the chain blocks V_{charge,g}, applies Q + Delta_BV and the
b-operators (sign rules of bg_complex_k.BGk), and computes ranks mod p by sparse elimination.
For k = 1 the complex is that of bg_complex.BG1 (b-operator R_0).

homology(alg, lam, k, N, p, ghs=None, mode=1, tables=None) -> {g: (dimV_T, rankS_g, rankd_g, dim H_g)}
    mode 0: all b-operators of bg_complex_k.b_ops (P0_xy, MM0_xy for x != y, all R0_xy);
    mode 1: {P0_{x,x+1}, MM0_{x,x+1}, R0_{0,k-1}}, same b-image at charge T (raising operators of
            gl_k + gl_k; highest-weight component of r = <R0_xy>).
    S_g is eliminated block by block (blocks = T-monomial and multiset of words, preserved by every
    b-operator); rank dbar_{g-1} uses d of the non-pivot monomials of S_{g-1} reduced modulo S_g.
    tables = small_tables(lam): HA^! model (E_1 page of the meson-number filtration); output = E_2.
Kernel environment flags: BGK_NOBLOCK (global elimination), BGK_TRANSPOSE (column-wise elimination),
BGK_DIAG / BGK_DIAG_ONLY (block statistics), BGK_BLOCKLOG (per-block timing), BGK_FULLD.
"""
import os
import subprocess
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, "native", "bgk")
SRC = os.path.join(HERE, "native", "bgk.cpp")
P1 = 2147483647
P2 = 2147483629


def build(force=False):
    """Compile the kernel if missing or older than the source (atomic replace: running kernels keep their image)."""
    if force or not os.path.exists(BIN) or os.path.getmtime(BIN) < os.path.getmtime(SRC):
        tmp = f"{BIN}.tmp{os.getpid()}"
        subprocess.check_call(["clang++", "-O3", "-std=c++17", "-march=native", "-o", tmp, SRC])
        os.replace(tmp, BIN)


def sector_tables(alg, lam):
    """Words (index 0 = empty) in Python tuple order, Q on words, cyclic classes, Q on traces."""
    lam = tuple(lam)
    words = [()] + sorted(set(alg.words_upto(lam)))
    widx = {w: i for i, w in enumerate(words)}
    cyc = set()
    canon = []
    for w in words:
        if not w:
            canon.append((0, None))
            continue
        s, cw = alg.canon_cyclic(w)
        canon.append((s, cw))
        if s:
            cyc.add(cw)
    cyc = sorted(cyc)
    cidx = {c: i for i, c in enumerate(cyc)}
    WQ = []
    for w in words:
        if not w:
            WQ.append([])
            continue
        WQ.append([(widx[nw], c) for nw, c in sorted(alg.Qword(w).items())])
    CQ = []
    for cw in cyc:
        acc = {}
        for nw, c in alg.Qword(cw).items():
            s, c2 = alg.canon_cyclic(nw)
            if s:
                acc[cidx[c2]] = acc.get(cidx[c2], 0) + c * s
        CQ.append(sorted((a, b) for a, b in acc.items() if b))
    return dict(words=words, widx=widx, canon=canon, cyc=cyc, cidx=cidx, WQ=WQ, CQ=CQ)


def write_input(alg, lam, k, N, p, ghs, mode, path, tables=None):
    tb = tables or sector_tables(alg, lam)
    words, canon, cyc, cidx = tb["words"], tb["canon"], tb["cyc"], tb["cidx"]
    out = [f"{k} {N} {p}", " ".join(map(str, lam)), str(len(words))]
    if tb.get("prodmode"):
        for w in words:
            out.append("0")
        for i in range(len(words)):
            out.append(" ".join(map(str, list(tb["wt"][i]) + [tb["par"][i], tb["gh"][i]])))
    else:
        for w in words:
            out.append(" ".join(map(str, [len(w)] + list(w))))
        for w in words:
            out.append(" ".join(map(str, list(alg.wword(w)) + [alg.pword(w), alg.ghword(w)])))
    for q in tb["WQ"]:
        out.append(" ".join(map(str, [len(q)] + [x for pair in q for x in pair])))
    for s, cw in canon:
        out.append(f"{s} {cidx[cw] if s else 0}")
    for i in range(len(words)):
        out.append(str(tb["mes"][i] if "mes" in tb else 1))
    out.append(str(len(cyc)))
    for cw in cyc:
        out.append(str(tb["widx"][cw]))
    for q in tb["CQ"]:
        out.append(" ".join(map(str, [len(q)] + [x for pair in q for x in pair])))
    out.append(str(1 if tb.get("prodmode") else 0))
    out.append(str(mode))
    out.append(" ".join(map(str, [len(ghs)] + list(ghs))))
    with open(path, "w") as fh:
        fh.write("\n".join(out) + "\n")


def small_tables(lam):
    """HA^! model (E_1 page of the meson-number filtration) for a sector with n+ = n- = 0:
    words = monomials x^a of C[x1,x2,x3] (gh 0, even; product = weight sum, Q = 0);
    T generators = reduced HC(C[x]): HC_0 class of x^a (= target of the self-contraction of P_{x^a}, MM_{x^a})
    and spectator classes HC_n, n >= 1, binom(s(a)-1, n) at weight a, gh -n, parity n mod 2."""
    import itertools
    from math import comb
    lam = tuple(lam)
    assert lam[3] == 0 and lam[4] == 0
    mons = [a for a in itertools.product(*[range(x + 1) for x in lam[:3]]) if sum(a)]
    mons.sort(key=lambda a: (sum(a), a))
    words = [("e",)]
    wt, par, gh, mes = [(0, 0, 0, 0, 0)], [0], [0], [1]
    for a in mons:
        words.append(("m",) + a)
        wt.append(tuple(a) + (0, 0)); par.append(0); gh.append(0); mes.append(1)
    spect = []
    for a in mons:
        s = sum(1 for x in a if x)
        for n in range(1, s):
            for j in range(comb(s - 1, n)):
                words.append(("h", n, j) + a)
                wt.append(tuple(a) + (0, 0)); par.append(n & 1); gh.append(-n); mes.append(0)
                spect.append(len(words) - 1)
    widx = {w: i for i, w in enumerate(words)}
    # cyclic classes: HC_0 classes first (one per monomial), then spectators
    cyc = [words[i] for i in range(1, len(mons) + 1)] + [words[i] for i in spect]
    cidx = {c: i for i, c in enumerate(cyc)}
    canon = [(0, None)] + [(1, words[i]) for i in range(1, len(mons) + 1)] + [(0, None)] * len(spect)
    WQ = [[] for _ in words]
    CQ = [[] for _ in cyc]
    return dict(words=words, widx=widx, canon=canon, cyc=cyc, cidx=cidx, WQ=WQ, CQ=CQ,
                prodmode=True, wt=wt, par=par, gh=gh, mes=mes)


def run_cols(alg, lam, k, N, p, cols, mode=1, tables=None, log=None, workdir=None):
    """Run the kernel on column spaces cols; returns ({g: (dimV, rankS, rankd_{g-1})}, raw output).
    The kernel output is streamed to `log` (appended) while it runs."""
    build()
    fd, path = tempfile.mkstemp(suffix=".in", dir=workdir)
    os.close(fd)
    fo, opath = tempfile.mkstemp(suffix=".out", dir=workdir)
    os.close(fo)
    try:
        write_input(alg, lam, k, N, p, cols, mode, path, tables)
        if log:
            with open(log, "a") as fh:
                fh.write(f"# bgk k={k} N={N} p={p} lam={tuple(lam)} cols={list(cols)} mode={mode}\n")
        with open(log if log else opath, "a") as fh:
            res = subprocess.run([BIN, path], stdout=fh, stderr=subprocess.PIPE, text=True)
        if res.returncode != 0:
            raise RuntimeError(f"bgk failed ({res.returncode}): {res.stderr[-2000:]}")
        txt = open(log if log else opath).read()
        if log:
            txt = txt[txt.rindex("# bgk"):]
    finally:
        os.remove(path)
        os.remove(opath)
    out = {}
    for line in txt.splitlines():
        if line.startswith("COL"):
            parts = line.split()
            g = int(parts[2])
            out[g] = (int(parts[4]), int(parts[6]), int(parts[8]))
    return out, txt


def homology(alg, lam, k, N, p=P1, ghs=None, mode=1, tables=None, log=None):
    """{g: (dimV_T, rankS_g, rankd_g, dim H_g)} for g in ghs (default: all g with V_T != 0).
    tables = small_tables(lam) gives the HA^! model (E_2 page of the meson-number spectral sequence)."""
    if ghs is None:
        from . import bg_sizes as bs
        L = 2 * sum(lam[:3]) + 3 * sum(lam[3:])
        gmin = -2 * (2 * sum(lam[:3]) + sum(lam[3:])) - 2
        gmax = sum(lam[3:])
        cand = list(range(gmin, gmax + 1))
        d = bs.chain_dims(alg, lam, k, cand, N) if k in (1, 2) else None
        ghs = [g for g in cand if d is None or round(d[((N,) * k, (N,) * k, g)]) > 0]
    ghs = sorted(ghs)
    # consecutive column spaces g0-1 .. g1+1: the d-rows of column g are a basis of V_{g-1}/S_{g-1}
    cols = list(range(ghs[0] - 1, ghs[-1] + 2))
    r, _ = run_cols(alg, lam, k, N, p, cols, mode, tables, log)
    out = {}
    for g in ghs:
        dimV, rS, rdprev = r[g]
        rd = r[g + 1][2]
        out[g] = (dimV, rS, rd, dimV - rS - rd - rdprev)
    return out
