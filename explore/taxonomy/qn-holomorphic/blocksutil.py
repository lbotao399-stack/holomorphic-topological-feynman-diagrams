"""Partition of the basis of g_tot into constant-parity blocks (finest available from the constructors' block labels)."""


def blocks_for(sl):
    N = sl.N
    names = sorted(sl.blocks, key=lambda b: (len(sl.blocks[b]), b))
    cover = {}
    for b in names:
        ids = sl.blocks[b]
        if len({int(sl.pi[i]) for i in ids}) != 1:
            continue
        if any(i in cover for i in ids):
            continue
        for i in ids:
            cover[i] = b
    out = {}
    for i in range(N):
        b = cover.get(i, 'par%d' % int(sl.pi[i]))
        out.setdefault(b, []).append(i)
    return out
