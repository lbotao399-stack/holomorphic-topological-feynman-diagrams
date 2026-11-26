/* Modular sparse solver for one level of the K6 family.
   Input files:
     rows file (xrows output): 'X c j c j ...' or 'D qm1 L<i> c j c j ...' (D: coefficient of J_{q-1}(L_i) is qm1)
     lower values file: first line 'nl nlab', then nl lines of nlab integers (values mod p, -1 = unknown class)
   Unknown columns 0..n-1 (n given).  Row: sum_j c_j J(j) + qm1 * sum_lab val(L_i, lab) * [lab] = 0.
   Rows referencing a lower class with unknown value are skipped.
   Elimination mod p with pivot = smallest column; then full back-substitution.
   Output: line 'rank R'; then for each column j: 'j det v_1 ... v_nlab' (det = 1 if determined).
   env MS_ROWS_OUT=file: write the indices of the rows that produced pivots;
   env MS_ROWS_IN=file: process only the rows with these indices (a pivot set found with another prime).
   usage: modsolve p n rowsfile lowerfile */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef int64_t i64;
static i64 P;
static int n, nlab;

static i64 md(i64 a) { a %= P; return a < 0 ? a + P : a; }
static i64 pw(i64 a, i64 e) { i64 r = 1; a = md(a); while (e) { if (e & 1) r = (__int128)r * a % P; a = (__int128)a * a % P; e >>= 1; } return r; }
static i64 inv(i64 a) { return pw(a, P - 2); }

typedef struct { int len, cap; int *col; i64 *val; i64 *rhs; } Row;

static long selcap = 0;
static Row *piv;      /* piv[j]: pivot row with leading column j (len 0 = none) */
static int *haspiv;

static void row_init(Row *r, int cap) { r->len = 0; r->cap = cap; r->col = malloc(cap * sizeof(int)); r->val = malloc(cap * sizeof(i64)); r->rhs = calloc(nlab, sizeof(i64)); }
static void row_free(Row *r) { free(r->col); free(r->val); free(r->rhs); }

/* r <- r - f * s  (both sorted by column) */
static int *tc; static i64 *tv; static int tcap = 0;
static void axpy(Row *r, i64 f, const Row *s) {
    int need = r->len + s->len;
    if (need > tcap) { tcap = need * 2; tc = realloc(tc, tcap * sizeof(int)); tv = realloc(tv, tcap * sizeof(i64)); }
    int i = 0, k = 0, m = 0;
    while (i < r->len || k < s->len) {
        if (k >= s->len || (i < r->len && r->col[i] < s->col[k])) { tc[m] = r->col[i]; tv[m] = r->val[i]; m++; i++; }
        else if (i >= r->len || s->col[k] < r->col[i]) { tc[m] = s->col[k]; tv[m] = md(-(__int128)f * s->val[k] % P); if (tv[m]) m++; k++; }
        else { i64 v = md(r->val[i] - (__int128)f * s->val[k] % P); if (v) { tc[m] = r->col[i]; tv[m] = v; m++; } i++; k++; }
    }
    if (m > r->cap) { r->cap = m * 2; r->col = realloc(r->col, r->cap * sizeof(int)); r->val = realloc(r->val, r->cap * sizeof(i64)); }
    memcpy(r->col, tc, m * sizeof(int)); memcpy(r->val, tv, m * sizeof(i64)); r->len = m;
    for (int l = 0; l < nlab; l++) r->rhs[l] = md(r->rhs[l] - (__int128)f * s->rhs[l] % P);
}

static int cmpcol(const void *a, const void *b) { return (*(const int *)a) - (*(const int *)b); }

int main(int argc, char **argv) {
    P = atoll(argv[1]);
    n = atoi(argv[2]);
    FILE *fl = fopen(argv[4], "r");
    int nl;
    if (fscanf(fl, "%d %d", &nl, &nlab) != 2) return 1;
    i64 *low = malloc((size_t)nl * nlab * sizeof(i64));
    char *lowok = malloc(nl);
    for (int i = 0; i < nl; i++) {
        lowok[i] = 1;
        for (int l = 0; l < nlab; l++) { long long x; if (fscanf(fl, "%lld", &x) != 1) return 1; if (x < 0) lowok[i] = 0; low[(size_t)i * nlab + l] = x < 0 ? 0 : x; }
    }
    fclose(fl);
    piv = calloc(n, sizeof(Row));
    haspiv = calloc(n, sizeof(int));
    char *sel = NULL;
    if (getenv("MS_ROWS_IN")) {
        FILE *fi = fopen(getenv("MS_ROWS_IN"), "r");
        long cap = 1 << 20, x;
        sel = calloc(cap, 1);
        while (fscanf(fi, "%ld", &x) == 1) { while (x >= cap) { sel = realloc(sel, cap * 2); memset(sel + cap, 0, cap); cap *= 2; } sel[x] = 1; }
        fclose(fi);
        /* rows beyond cap: not selected */
        sel = realloc(sel, cap + 1); 
        selcap = cap;
    }
    FILE *fo = getenv("MS_ROWS_OUT") ? fopen(getenv("MS_ROWS_OUT"), "w") : NULL;
    FILE *fr = fopen(argv[3], "r");
    char *line = NULL; size_t lcap = 0;
    long nrows = 0, nskip = 0, rank = 0, ninc = 0;
    int *cols = malloc(4096 * sizeof(int)); i64 *vals = malloc(4096 * sizeof(i64));
    while (getline(&line, &lcap, fr) > 0) {
        char *s = line;
        int isD = (s[0] == 'D');
        s += 1;
        Row r; row_init(&r, 16);
        int skip = 0;
        if (isD) {
            long qm1 = strtol(s, &s, 10);
            while (*s == ' ') s++;
            s++; /* 'L' */
            long li = strtol(s, &s, 10);
            if (!lowok[li]) skip = 1;
            else for (int l = 0; l < nlab; l++) r.rhs[l] = md((__int128)qm1 * low[(size_t)li * nlab + l] % P);
        }
        int m = 0;
        for (;;) {
            char *e;
            long c = strtol(s, &e, 10);
            if (e == s) break;
            s = e;
            long j = strtol(s, &e, 10);
            s = e;
            cols[m] = (int)j; vals[m] = md(c); m++;
        }
        nrows++;
        if (sel && (nrows - 1 >= selcap || !sel[nrows - 1])) skip = 2;
        if (skip) { if (skip == 1) nskip++; row_free(&r); continue; }
        /* sort and merge duplicates */
        int idx[4096];
        for (int i = 0; i < m; i++) idx[i] = i;
        /* simple insertion sort by column */
        for (int i = 1; i < m; i++) { int t = idx[i], k = i - 1; while (k >= 0 && cols[idx[k]] > cols[t]) { idx[k + 1] = idx[k]; k--; } idx[k + 1] = t; }
        for (int i = 0; i < m; i++) {
            int c = cols[idx[i]]; i64 v = vals[idx[i]];
            if (r.len && r.col[r.len - 1] == c) { r.val[r.len - 1] = md(r.val[r.len - 1] + v); if (!r.val[r.len - 1]) r.len--; }
            else if (v) { if (r.len == r.cap) { r.cap *= 2; r.col = realloc(r.col, r.cap * sizeof(int)); r.val = realloc(r.val, r.cap * sizeof(i64)); } r.col[r.len] = c; r.val[r.len] = v; r.len++; }
        }
        /* reduce */
        while (r.len) {
            int c = r.col[0];
            if (haspiv[c]) { i64 f = r.val[0]; axpy(&r, f, &piv[c]); }
            else {
                i64 f = inv(r.val[0]);
                for (int i = 0; i < r.len; i++) r.val[i] = (__int128)r.val[i] * f % P;
                for (int l = 0; l < nlab; l++) r.rhs[l] = (__int128)r.rhs[l] * f % P;
                piv[c] = r; haspiv[c] = 1; rank++;
                if (fo) fprintf(fo, "%ld\n", nrows - 1);
                break;
            }
        }
        if (!r.len) {
            int nz = 0;
            for (int l = 0; l < nlab; l++) if (r.rhs[l]) nz = 1;
            if (nz) ninc++;
            row_free(&r);
        }
    }
    fclose(fr);
    if (fo) fclose(fo);
    /* back substitution: descending pivot columns */
    for (int c = n - 1; c >= 0; c--) {
        if (!haspiv[c]) continue;
        Row *r = &piv[c];
        int changed = 1;
        while (changed) {
            changed = 0;
            for (int i = 1; i < r->len; i++) {
                int cc = r->col[i];
                if (haspiv[cc]) { i64 f = r->val[i]; axpy(r, f, &piv[cc]); changed = 1; break; }
            }
        }
    }
    printf("rank %ld rows %ld skipped %ld inconsistent %ld\n", rank, nrows, nskip, ninc);
    for (int c = 0; c < n; c++) {
        int det = haspiv[c] && piv[c].len == 1;
        printf("%d %d", c, det);
        for (int l = 0; l < nlab; l++) printf(" %lld", det ? (long long)md(-piv[c].rhs[l]) : 0LL);
        printf("\n");
    }
    return 0;
}
