/* Relation rows for the K6 family at level q (integer coefficients).
   (X)  for each class b and ordered distinct (x, x1, c):
        (b_xc - b_xx1) J(b) - sum_i b_{x xi} [ J(b - e_{x x1} + e_{x xi} + e_{x1 c} - e_{xi c})
                                             - J(b + e_{x xi} - e_{x1 xi} + e_{x1 c} - e_{xc}) ] = 0
        emitted only if every term with nonzero coefficient is convergent.
   (DS) for each level-(q-1) class b and vertex c:
        (q-1) J_{q-1}(b) - sum_T prod_{e in T} b_e J_q(b + 1_T, completed at c) = 0,
        T spanning trees of K6 - c; emitted only if every term with nonzero coefficient is convergent.
   Input: classes_q.txt (level q, one canonical vector per line, first 15 ints), optionally classes_{q-1}.txt.
   Output (stdout): one row per line: 'X' or 'D', then pairs 'coef idx' (idx = level-q class index;
   for D rows the level-(q-1) class is written as 'coef L<idx>').
   usage: xrows q classes_q.txt [classes_qm1.txt] */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int EI[6][6];
static int ED[15][2];
static int PM[720][15];
static int NP = 0;
static int q;
static int TRIE[20][3];

typedef struct { int8_t k[15]; int used; int idx; } Slot;
static Slot *H;
static size_t HCAP;

static uint64_t hsh(const int8_t *k) {
    uint64_t h = 1469598103934665603ULL;
    for (int i = 0; i < 15; i++) { h ^= (uint8_t)k[i]; h *= 1099511628211ULL; }
    return h;
}
static void hput(const int8_t *k, int idx) {
    uint64_t h = hsh(k) & (HCAP - 1);
    while (H[h].used) { if (!memcmp(H[h].k, k, 15)) return; h = (h + 1) & (HCAP - 1); }
    memcpy(H[h].k, k, 15); H[h].used = 1; H[h].idx = idx;
}
static int hget(const int8_t *k) {
    uint64_t h = hsh(k) & (HCAP - 1);
    while (H[h].used) { if (!memcmp(H[h].k, k, 15)) return H[h].idx; h = (h + 1) & (HCAP - 1); }
    return -1;
}

static void perms(void) {
    int s[6] = {0, 1, 2, 3, 4, 5};
    for (;;) {
        for (int k = 0; k < 15; k++) PM[NP][EI[s[ED[k][0]]][s[ED[k][1]]]] = k;
        NP++;
        int i = 4;
        while (i >= 0 && s[i] > s[i + 1]) i--;
        if (i < 0) break;
        int j = 5;
        while (s[j] < s[i]) j--;
        int t = s[i]; s[i] = s[j]; s[j] = t;
        for (int a = i + 1, c = 5; a < c; a++, c--) { t = s[a]; s[a] = s[c]; s[c] = t; }
    }
}

static void canon(const int *b, int8_t *best) {
    int m = b[0];
    for (int k = 1; k < 15; k++) if (b[k] < m) m = b[k];
    int8_t cur[15];
    int have = 0;
    for (int p = 0; p < NP; p++) {
        if (b[PM[p][0]] != m) continue;
        int better = !have, decided = !have;
        for (int k = 0; k < 15; k++) {
            int8_t x = (int8_t)b[PM[p][k]];
            cur[k] = x;
            if (!decided) {
                if (x < best[k]) { better = 1; decided = 1; }
                else if (x > best[k]) { better = 0; decided = 1; break; }
            }
        }
        if (better) { memcpy(best, cur, 15); have = 1; }
    }
}

static int convergent(const int *b, int qq) {
    for (int k = 0; k < 15; k++) if (b[k] >= qq) return 0;
    for (int t = 0; t < 20; t++) if (b[TRIE[t][0]] + b[TRIE[t][1]] + b[TRIE[t][2]] >= 2 * qq) return 0;
    return 1;
}

static int read_classes(const char *fn, int8_t **out) {
    FILE *f = fopen(fn, "r");
    if (!f) { perror(fn); exit(1); }
    int cap = 1024, n = 0;
    int8_t *a = malloc(cap * 15);
    char line[512];
    while (fgets(line, sizeof line, f)) {
        int v[16];
        int r = sscanf(line, "%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d", &v[0], &v[1], &v[2], &v[3], &v[4], &v[5],
                       &v[6], &v[7], &v[8], &v[9], &v[10], &v[11], &v[12], &v[13], &v[14]);
        if (r != 15) continue;
        if (n == cap) { cap *= 2; a = realloc(a, cap * 15); }
        for (int k = 0; k < 15; k++) a[n * 15 + k] = (int8_t)v[k];
        n++;
    }
    fclose(f);
    *out = a;
    return n;
}

/* spanning trees of K5 on vertex set (complement of c): list of 4-edge sets */
static int TREES[6][125][4];
static int NT[6];

static int find(int *par, int v) { while (par[v] != v) v = par[v] = par[par[v]]; return v; }

static void trees(void) {
    for (int c = 0; c < 6; c++) {
        int es[10], ne = 0;
        for (int k = 0; k < 15; k++) if (ED[k][0] != c && ED[k][1] != c) es[ne++] = k;
        NT[c] = 0;
        for (int m = 0; m < (1 << 10); m++) {
            if (__builtin_popcount(m) != 4) continue;
            int par[6] = {0, 1, 2, 3, 4, 5}, ok = 1, t[4], nt = 0;
            for (int i = 0; i < 10 && ok; i++) if (m >> i & 1) {
                int a = find(par, ED[es[i]][0]), bb = find(par, ED[es[i]][1]);
                if (a == bb) ok = 0; else par[a] = bb;
                t[nt++] = es[i];
            }
            if (ok) { memcpy(TREES[c][NT[c]], t, sizeof t); NT[c]++; }
        }
    }
}

int main(int argc, char **argv) {
    q = atoi(argv[1]);
    int k = 0;
    for (int i = 0; i < 6; i++) for (int j = i + 1; j < 6; j++) { EI[i][j] = EI[j][i] = k; ED[k][0] = i; ED[k][1] = j; k++; }
    perms();
    trees();
    int t = 0;
    for (int a = 0; a < 6; a++) for (int bb = a + 1; bb < 6; bb++) for (int cc = bb + 1; cc < 6; cc++) {
        TRIE[t][0] = EI[a][bb]; TRIE[t][1] = EI[a][cc]; TRIE[t][2] = EI[bb][cc]; t++;
    }
    int8_t *C;
    int n = read_classes(argv[2], &C);
    HCAP = 1; while (HCAP < (size_t)n * 4) HCAP <<= 1;
    H = calloc(HCAP, sizeof(Slot));
    for (int i = 0; i < n; i++) hput(C + 15 * i, i);
    int8_t cb[15];
    long nrows = 0, missing = 0;
    /* X rows */
    for (int i = 0; i < n; i++) {
        int b[15];
        for (int z = 0; z < 15; z++) b[z] = C[15 * i + z];
        for (int x = 0; x < 6; x++) for (int x1 = 0; x1 < 6; x1++) for (int c = 0; c < 6; c++) {
            if (x == x1 || x == c || x1 == c) continue;
            int coefs[8], idxs[8], nn = 0, ok = 1;
            coefs[nn] = b[EI[x][c]] - b[EI[x][x1]]; idxs[nn] = i; nn++;
            for (int xi = 0; xi < 6 && ok; xi++) {
                if (xi == x || xi == x1 || xi == c) continue;
                int w = b[EI[x][xi]];
                if (w == 0) continue;
                int t1[15], t2[15];
                memcpy(t1, b, sizeof b); memcpy(t2, b, sizeof b);
                t1[EI[x][x1]] -= 1; t1[EI[x][xi]] += 1; t1[EI[x1][c]] += 1; t1[EI[xi][c]] -= 1;
                t2[EI[x][xi]] += 1; t2[EI[x1][xi]] -= 1; t2[EI[x1][c]] += 1; t2[EI[x][c]] -= 1;
                if (!convergent(t1, q) || !convergent(t2, q)) { ok = 0; break; }
                canon(t1, cb); int j1 = hget(cb);
                canon(t2, cb); int j2 = hget(cb);
                if (j1 < 0 || j2 < 0) { missing++; ok = 0; break; }
                coefs[nn] = -w; idxs[nn] = j1; nn++;
                coefs[nn] = w; idxs[nn] = j2; nn++;
            }
            if (!ok) continue;
            printf("X");
            for (int z = 0; z < nn; z++) printf(" %d %d", coefs[z], idxs[z]);
            printf("\n");
            nrows++;
        }
    }
    long nds = 0;
    if (argc > 3) {
        int8_t *L;
        int nl = read_classes(argv[3], &L);
        for (int i = 0; i < nl; i++) {
            int b[15];
            for (int z = 0; z < 15; z++) b[z] = L[15 * i + z];
            for (int c = 0; c < 6; c++) {
                int ok = 1;
                long coefs[130]; int idxs[130]; int nn = 0;
                for (int ti = 0; ti < NT[c] && ok; ti++) {
                    long coef = 1;
                    for (int z = 0; z < 4; z++) coef *= b[TREES[c][ti][z]];
                    if (coef == 0) continue;
                    int nb[15];
                    memcpy(nb, b, sizeof b);
                    for (int z = 0; z < 4; z++) nb[TREES[c][ti][z]] += 1;
                    for (int v = 0; v < 6; v++) {
                        if (v == c) continue;
                        int s = 0;
                        for (int w = 0; w < 6; w++) if (w != v) s += nb[EI[v][w]];
                        nb[EI[v][c]] += 2 * q - s;
                    }
                    if (!convergent(nb, q)) { ok = 0; break; }
                    canon(nb, cb);
                    int j = hget(cb);
                    if (j < 0) { missing++; ok = 0; break; }
                    int found = -1;
                    for (int z = 0; z < nn; z++) if (idxs[z] == j) { found = z; break; }
                    if (found >= 0) coefs[found] += coef; else { coefs[nn] = coef; idxs[nn] = j; nn++; }
                }
                if (!ok) continue;
                printf("D %d L%d", q - 1, i);
                for (int z = 0; z < nn; z++) if (coefs[z]) printf(" %ld %d", -coefs[z], idxs[z]);
                printf("\n");
                nds++;
            }
        }
    }
    fprintf(stderr, "q=%d classes %d X rows %ld DS rows %ld missing-term rows %ld\n", q, n, nrows, nds, missing);
    return 0;
}
