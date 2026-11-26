/* Enumerate S6 classes of balanced convergent K6 weights at level q.
   beta in Z^15 (edges (i<j) lexicographic), sum over the 5 edges at each vertex = 2q,
   beta_e <= q-1, triangle sums <= 2q-1, beta_e >= 4-2q.
   Symmetry breaking: beta_01 = global minimum.  Canonical form = lexicographic minimum over S6
   of beta'[k] = beta[PM[p][k]] (same convention as k6.py).  Output: one canonical vector per line
   followed by the number of labelled vectors with beta_01 = min that map to it.
   usage: enum6 q > file */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int EI[6][6];
static int ED[15][2];
static int PM[720][15];
static int NP = 0;
static int q, lo, hi;
static int b[15];
static int TRI[20][3];

/* hash set of canonical forms */
typedef struct { int8_t k[15]; int8_t used; int64_t cnt; } Slot;
static Slot *H;
static size_t HN = 0, HCAP = 1 << 22;

static uint64_t hsh(const int8_t *k) {
    uint64_t h = 1469598103934665603ULL;
    for (int i = 0; i < 15; i++) { h ^= (uint8_t)k[i]; h *= 1099511628211ULL; }
    return h;
}

static void grow(void);
static void insert(const int8_t *k) {
    if (HN * 2 > HCAP) grow();
    uint64_t h = hsh(k) & (HCAP - 1);
    while (H[h].used) {
        if (memcmp(H[h].k, k, 15) == 0) { H[h].cnt++; return; }
        h = (h + 1) & (HCAP - 1);
    }
    memcpy(H[h].k, k, 15); H[h].used = 1; H[h].cnt = 1; HN++;
}
static void grow(void) {
    Slot *old = H; size_t oc = HCAP;
    HCAP *= 2; H = calloc(HCAP, sizeof(Slot)); HN = 0;
    for (size_t i = 0; i < oc; i++) if (old[i].used) {
        uint64_t h = hsh(old[i].k) & (HCAP - 1);
        while (H[h].used) h = (h + 1) & (HCAP - 1);
        H[h] = old[i]; HN++;
    }
    free(old);
}

static void perms(void) {
    int s[6] = {0, 1, 2, 3, 4, 5};
    /* Heap's algorithm would reorder; use lexicographic next_permutation */
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

static void canon_and_insert(void) {
    int m = b[0];
    int8_t best[15], cur[15];
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
    insert(best);
}

/* fill order: vertex v, edges (v,w) w>v */
static int ord[15];
static int vstart[7];
static int trilast[15][20][3];
static int ntrilast[15];
static int pos_of[15];

static int tri_ok(int k) {
    for (int t = 0; t < ntrilast[k]; t++)
        if (b[trilast[k][t][0]] + b[trilast[k][t][1]] + b[trilast[k][t][2]] >= 2 * q) return 0;
    return 1;
}

static int vmin; /* = b[0] once set */
static void rec(int v);

static void fill(int v, int i, int n, int rem) {
    int k = ord[vstart[v] + i];
    int l = (k == 0) ? lo : vmin;
    if (i == n - 1) {
        if (rem >= l && rem <= hi) {
            b[k] = rem;
            if (k == 0) vmin = rem;
            if (tri_ok(k)) rec(v + 1);
        }
        return;
    }
    int left = n - i - 1;
    int lmin = (k == 0) ? lo : vmin;
    int a0 = rem - hi * left; if (a0 < lmin) a0 = lmin;
    int a1 = rem - lmin * left; if (a1 > hi) a1 = hi;
    /* for k == 0 the remaining edges must be >= b[0] */
    for (int x = a0; x <= a1; x++) {
        b[k] = x;
        if (k == 0) {
            vmin = x;
            if (rem - x < x * left || rem - x > hi * left) continue;
        }
        if (tri_ok(k)) fill(v, i + 1, n, rem - x);
    }
}

static void rec(int v) {
    if (v == 5) {
        int s = 0;
        for (int w = 0; w < 5; w++) s += b[EI[w][5]];
        if (s == 2 * q) canon_and_insert();
        return;
    }
    int fixed = 0;
    for (int w = 0; w < v; w++) fixed += b[EI[w][v]];
    int n = 5 - v;
    fill(v, 0, n, 2 * q - fixed);
}

int main(int argc, char **argv) {
    q = atoi(argv[1]);
    lo = 4 - 2 * q; hi = q - 1;
    int k = 0;
    for (int i = 0; i < 6; i++) for (int j = i + 1; j < 6; j++) { EI[i][j] = EI[j][i] = k; ED[k][0] = i; ED[k][1] = j; k++; }
    perms();
    int c = 0;
    for (int v = 0; v < 6; v++) { vstart[v] = c; for (int w = v + 1; w < 6; w++) ord[c++] = EI[v][w]; }
    vstart[6] = c;
    for (int i = 0; i < 15; i++) pos_of[ord[i]] = i;
    int t = 0;
    for (int a = 0; a < 6; a++) for (int bb = a + 1; bb < 6; bb++) for (int cc = bb + 1; cc < 6; cc++) {
        int e[3] = {EI[a][bb], EI[a][cc], EI[bb][cc]};
        int last = e[0];
        for (int z = 1; z < 3; z++) if (pos_of[e[z]] > pos_of[last]) last = e[z];
        int n = ntrilast[last]++;
        trilast[last][n][0] = e[0]; trilast[last][n][1] = e[1]; trilast[last][n][2] = e[2];
        t++;
    }
    H = calloc(HCAP, sizeof(Slot));
    vmin = lo;
    rec(0);
    for (size_t i = 0; i < HCAP; i++) if (H[i].used) {
        for (int j = 0; j < 15; j++) printf("%d ", H[i].k[j]);
        printf("%lld\n", (long long)H[i].cnt);
    }
    fprintf(stderr, "q=%d classes %zu\n", q, HN);
    return 0;
}
