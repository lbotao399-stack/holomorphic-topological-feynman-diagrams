def relabel_any(nv, edges, rng):
    """random vertex permutation, random edge reversals, random edge order, random external order"""
    pi = list(range(nv)); rng.shuffle(pi)
    E = []
    for a, b in edges:
        a2, b2 = pi[a], pi[b]
        if rng.random() < 0.5:
            a2, b2 = b2, a2
        E.append((a2, b2))
    rng.shuffle(E)
    deg = [0] * nv
    for a, b in E:
        deg[a] += 1; deg[b] += 1
    ext = [v for v in range(nv) if deg[v] == 2]
    rng.shuffle(ext)
    return E, ext
