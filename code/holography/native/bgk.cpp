// bgk.cpp -- Budzik-Gaiotto k-meson complex (conventions of holography/bg_complex_k.py), homology mod p.
//
// Input (text, written by holography/bg_fast.py):
//   k N p
//   lam[5]
//   W                      number of words (index 0 = empty word), Python tuple order
//   W lines: len l_1 .. l_len            (product mode: 0)
//   W lines: wt[5] par gh
//   W lines: nQ (widx coeff)*            Q_r on the word
//   W lines: csign cidx                  canonical graded-cyclic class of the word (0 0 if it vanishes)
//   W lines: mesflag                     1 if meson generators exist on the word
//   C                      number of trace generators T (cyclic classes), sorted
//   C lines: widx                        word carrying the weight/parity/gh of the trace
//   C lines: nQ (cidx coeff)*            Q on the trace (cyclic signs included)
//   prodmode               0: concatenation of letter strings; 1: commutative monomials (product = weight sum)
//   mode                   0: all b-operators (bg_complex_k.b_ops); 1: {P0_{x,x+1}, MM0_{x,x+1}, R0_{0,k-1}}
//   ng g_1 .. g_ng         column spaces V_{T,g}, processed in increasing order
// Output: "COL g dimV rankS rankd_prev ..." per column; rankd_prev = rank of dbar_{g-1} computed from the
//   non-pivot monomials of S_{g-1} (-1 if column g-1 was not processed).
// Generators: T: id = cidx (0..C-1); meson (t, w, a, b), t in {M,P,MM,R} = {1,2,3,4}:
//   id = C + ((t-1)*W + w)*k*k + a*k + b  (Python order (t, w, a, b) == id order).
// b-modes: t in {P,MM,R} with w = 0 (never in V; results containing them are dropped).
// S_g is eliminated per skeleton block (T-monomial, multiset of words), which every b-operator preserves.
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <chrono>
#include <array>
#include <map>

using namespace std;
typedef uint32_t gid;
typedef int64_t i64;
typedef uint64_t u64;

static int K, NN;
static u64 PR;
static int LAM[5];
static int W, C;
static vector<vector<uint8_t>> WL;      // word letters
static vector<array<int,5>> WWT;
static vector<int> WPAR, WGH;
static vector<vector<pair<int,i64>>> WQ;
static vector<int> WCS, WCI;            // canonical cyclic sign, index
static vector<int> CW;                  // cyclic -> word index
static vector<vector<pair<int,i64>>> CQ;
static unordered_map<string,int> WIDX;

enum { T_ = 0, M_ = 1, P_ = 2, MM_ = 3, R_ = 4 };

// generator decoding
static inline int gtype(gid g) { if ((int)g < C) return T_; return 1 + (int)((g - C) / ((u64)W * K * K)); }
static inline int gword(gid g) { if ((int)g < C) return CW[g]; return (int)(((g - C) / (K * K)) % W); }
static inline int gfa(gid g) { if ((int)g < C) return 0; return (int)(((g - C) % (K * K)) / K); }
static inline int gfb(gid g) { if ((int)g < C) return 0; return (int)((g - C) % K); }
static inline gid mkgen(int t, int w, int a, int b) { return (gid)(C + (((u64)(t - 1) * W + w) * K * K) + a * K + b); }
static inline int gpar(gid g) {
    int t = gtype(g); int pw = WPAR[gword(g)];
    return (t == T_ || t == M_ || t == R_) ? pw : 1 - pw;
}
static inline int ggh(gid g) { static const int sh[5] = {0, 0, 1, 1, 2}; return WGH[gword(g)] - sh[gtype(g)]; }
static inline bool isbmode(gid g) { int t = gtype(g); return (t == P_ || t == MM_ || t == R_) && gword(g) == 0; }
static inline void gcharge(gid g, int *ch) { // ch[0..k-1] = Lc, ch[k..2k-1] = Rc (added)
    int t = gtype(g); if (t == T_) return;
    int a = gfa(g), b = gfb(g);
    if (t == M_) { ch[a] += 1; ch[K + b] += 1; }
    else if (t == P_) { ch[a] += 1; ch[b] -= 1; }
    else if (t == MM_) { ch[K + a] -= 1; ch[K + b] += 1; }
    else { ch[K + a] -= 1; ch[b] -= 1; }
}
static inline bool START_PB(int t) { return t == M_ || t == P_; }
static inline bool START_UB(int t) { return t == MM_ || t == R_; }
static inline bool END_PSI(int t) { return t == M_ || t == MM_; }
static inline bool END_U(int t) { return t == P_ || t == R_; }
static inline int mtype(bool spb, bool epsi) { return spb ? (epsi ? M_ : P_) : (epsi ? MM_ : R_); }

static int PRODMODE = 0;
static size_t PROGRESS = 0;               // 1: words are monomials of a commutative even algebra, product = weight sum
static vector<int> WMES;                // 1 if mesons exist on the word
static map<array<int,5>, int> WBYWT;    // product mode: meson-word by weight
static int concat(int wa, int wb) { // word index of wa + wb, -1 if absent
    if (wa == 0) return wb;
    if (wb == 0) return wa;
    if (PRODMODE) {
        array<int,5> w; for (int j = 0; j < 5; j++) w[j] = WWT[wa][j] + WWT[wb][j];
        auto it = WBYWT.find(w);
        return it == WBYWT.end() ? -1 : it->second;
    }
    string s; s.reserve(WL[wa].size() + WL[wb].size());
    for (auto x : WL[wa]) s.push_back((char)x);
    for (auto x : WL[wb]) s.push_back((char)x);
    auto it = WIDX.find(s);
    return it == WIDX.end() ? -1 : it->second;
}

// sort with Koszul sign; returns 0 if a repeated odd generator
static int sort_sign(vector<gid> &v) {
    int sign = 1;
    for (size_t i = 1; i < v.size(); i++) {
        gid x = v[i]; int px = gpar(x);
        size_t m = i;
        while (m > 0 && v[m - 1] > x) {
            if (px && gpar(v[m - 1])) sign = -sign;
            v[m] = v[m - 1]; m--;
        }
        v[m] = x;
    }
    for (size_t i = 1; i < v.size(); i++) if (v[i] == v[i - 1] && gpar(v[i])) return 0;
    return sign;
}

// ------------------------------------------------------------------ monomial store with hash
struct Block {
    vector<gid> data; vector<u64> off;  // monomial i = data[off[i] .. off[i+1])
    vector<int64_t> table; u64 mask = 0; size_t n = 0;
    Block() { off.push_back(0); }
    static u64 hsh(const gid *p, size_t len) {
        u64 h = 1469598103934665603ULL ^ (len * 0x9E3779B97F4A7C15ULL);
        for (size_t i = 0; i < len; i++) { h ^= p[i] + 0x9E3779B97F4A7C15ULL + (h << 6) + (h >> 2); h *= 1099511628211ULL; }
        return h ^ (h >> 29);
    }
    void add(const vector<gid> &m) { data.insert(data.end(), m.begin(), m.end()); off.push_back(data.size()); n++; }
    void build() {
        u64 cap = 16; while (cap < 2 * n + 16) cap <<= 1;
        table.assign(cap, -1); mask = cap - 1;
        for (size_t i = 0; i < n; i++) {
            u64 h = hsh(&data[off[i]], off[i + 1] - off[i]) & mask;
            while (table[h] >= 0) h = (h + 1) & mask;
            table[h] = (int64_t)i;
        }
    }
    int64_t find(const vector<gid> &m) const {
        if (n == 0) return -1;
        u64 h = hsh(m.data(), m.size()) & mask;
        while (true) {
            int64_t i = table[h];
            if (i < 0) return -1;
            size_t len = off[i + 1] - off[i];
            if (len == m.size() && memcmp(&data[off[i]], m.data(), len * sizeof(gid)) == 0) return i;
            h = (h + 1) & mask;
        }
    }
    void get(size_t i, vector<gid> &m) const { m.assign(data.begin() + off[i], data.begin() + off[i + 1]); }
};

// ------------------------------------------------------------------ enumeration
struct Target { vector<int> ch; int g; Block *blk; };
static vector<Target> TARGETS;
static vector<gid> EG;          // enumeration generators (no M_0, no b-modes), sorted by (class, id)
static vector<int> ECLS, CLS_END;
static vector<array<int,5>> CLS_WT;
static vector<int> EGH, EPAR; static vector<vector<int>> ECH;
static int GHMIN_T, GHMAX_T, NPOS;

static void leaf(vector<gid> &cur, const int *ch, int gh) {
    for (auto &t : TARGETS) {
        if (t.g != gh) continue;
        // M_0 completion: rows r_a = Lc_t[a]-Lc[a], cols c_b = Rc_t[b]-Rc[b]
        int r[3], c[3]; bool ok = true; int sr = 0, sc = 0;
        for (int a = 0; a < K; a++) { r[a] = t.ch[a] - ch[a]; c[a] = t.ch[K + a] - ch[K + a]; if (r[a] < 0 || c[a] < 0) ok = false; sr += r[a]; sc += c[a]; }
        if (!ok || sr != sc) continue;
        // enumerate KxK nonneg matrices with row sums r, col sums c
        vector<int> mvals(K * K, 0);
        vector<gid> mono;
        auto emit = [&]() {
            mono = cur;
            for (int a = 0; a < K; a++) for (int b = 0; b < K; b++) for (int j = 0; j < mvals[a * K + b]; j++) mono.push_back(mkgen(M_, 0, a, b));
            sort(mono.begin(), mono.end());
            t.blk->add(mono);
        };
        if (K == 1) { mvals[0] = r[0]; emit(); }
        else if (K == 2) {
            for (int m11 = 0; m11 <= min(r[0], c[0]); m11++) {
                int m12 = r[0] - m11, m21 = c[0] - m11, m22 = r[1] - m21;
                if (m12 < 0 || m21 < 0 || m22 < 0 || m12 + m22 != c[1]) continue;
                mvals = {m11, m12, m21, m22}; emit();
            }
        } else if (K == 3) {
            for (int m00 = 0; m00 <= min(r[0], c[0]); m00++)
            for (int m01 = 0; m01 <= min(r[0] - m00, c[1]); m01++) {
                int m02 = r[0] - m00 - m01; if (m02 < 0 || m02 > c[2]) continue;
                for (int m10 = 0; m10 <= min(r[1], c[0] - m00); m10++)
                for (int m11 = 0; m11 <= min(r[1] - m10, c[1] - m01); m11++) {
                    int m12 = r[1] - m10 - m11; if (m12 < 0 || m12 > c[2] - m02) continue;
                    int m20 = c[0] - m00 - m10, m21 = c[1] - m01 - m11, m22 = c[2] - m02 - m12;
                    if (m20 < 0 || m21 < 0 || m22 < 0 || m20 + m21 + m22 != r[2]) continue;
                    mvals = {m00, m01, m02, m10, m11, m12, m20, m21, m22}; emit();
                }
            }
        }
    }
}

static long long NLEAF = 0;
static void dfs(size_t pos, int *rem, int *ch, int gh, vector<gid> &cur) {
    bool zero = true; int units = 0;
    for (int i = 0; i < 5; i++) { if (rem[i]) zero = false; units += rem[i]; }
    if (zero) { NLEAF++; leaf(cur, ch, gh); return; }
    int npos_rem = rem[3] + rem[4];
    if (gh + npos_rem < GHMIN_T) return;          // gh can rise by at most #z-units
    // charge pruning: Lc_x can drop by at most 1 per further generator (<= units)
    for (int a = 0; a < 2 * K; a++) if (ch[a] > NN + units + 1) return;
    size_t i = pos;
    while (i < EG.size()) {
        int c = ECLS[i];
        const auto &w = CLS_WT[c];
        bool fits = true;
        for (int j = 0; j < 5; j++) if (w[j] > rem[j]) { fits = false; break; }
        if (!fits) { i = CLS_END[c]; continue; }
        gid g = EG[i];
        for (int j = 0; j < 5; j++) rem[j] -= w[j];
        for (int a = 0; a < 2 * K; a++) ch[a] += ECH[i][a];
        cur.push_back(g);
        dfs(EPAR[i] ? i + 1 : i, rem, ch, gh + EGH[i], cur);
        cur.pop_back();
        for (int a = 0; a < 2 * K; a++) ch[a] -= ECH[i][a];
        for (int j = 0; j < 5; j++) rem[j] += w[j];
        i++;
    }
}

// ------------------------------------------------------------------ operators
typedef vector<pair<vector<gid>, i64>> Terms;

static void emit_term(Terms &out, vector<gid> v, i64 coef) {
    int s = sort_sign(v);
    if (!s) return;
    out.push_back({std::move(v), coef * s});
}

static void Qop(const vector<gid> &gens, Terms &out) {
    int pre = 0;
    for (size_t pos = 0; pos < gens.size(); pos++) {
        gid g = gens[pos];
        i64 s0 = pre ? -1 : 1;
        int t = gtype(g);
        if (t == T_) {
            for (auto &pc : CQ[g]) { vector<gid> nv = gens; nv[pos] = (gid)pc.first; emit_term(out, nv, s0 * pc.second); }
        } else {
            int w = gword(g), a = gfa(g), b = gfb(g);
            i64 sg = START_PB(t) ? -1 : 1;
            for (auto &pc : WQ[w]) { vector<gid> nv = gens; nv[pos] = mkgen(t, pc.first, a, b); emit_term(out, nv, s0 * sg * pc.second); }
        }
        pre ^= gpar(g);
    }
}

// Delta_BV with Tr 1 = N; only >= 0 keeps contractions involving that position
static void Dop(const vector<gid> &gens, int N, int only, Terms &out) {
    int n = (int)gens.size();
    vector<int> pars(n), prefix(n + 1, 0), tt(n), ww(n), fa(n), fb(n);
    for (int i = 0; i < n; i++) { pars[i] = gpar(gens[i]); prefix[i + 1] = prefix[i] + pars[i]; tt[i] = gtype(gens[i]); ww[i] = gword(gens[i]); fa[i] = gfa(gens[i]); fb[i] = gfb(gens[i]); }
    // type 1: d_psibar d_u : psibar at start of gens[a], u at end of gens[b]
    for (int b = 0; b < n; b++) {
        if (!END_U(tt[b])) continue;
        for (int a = 0; a < n; a++) {
            if (!START_PB(tt[a]) || fa[a] != fb[b]) continue;
            if (only >= 0 && only != a && only != b) continue;
            i64 s = (prefix[a] & 1) ? -1 : 1;
            if (a == b) {
                if (ww[a] == 0) {
                    if (N) { vector<gid> rest; for (int i = 0; i < n; i++) if (i != a) rest.push_back(gens[i]); emit_term(out, rest, s * N); }
                } else {
                    int sc = WCS[ww[a]];
                    if (sc) { vector<gid> nv = gens; nv[a] = (gid)WCI[ww[a]]; emit_term(out, nv, s * sc); }
                }
                continue;
            }
            int pA = 1 - pars[a], pB = pars[b];
            int nw = concat(ww[b], ww[a]);
            if (nw < 0) { fprintf(stderr, "concat missing\n"); exit(2); }
            gid ng = mkgen(mtype(START_PB(tt[b]), END_PSI(tt[a])), nw, fa[b], fb[a]);
            i64 sm;
            if (a < b) { int mid = prefix[b] - prefix[a + 1]; sm = ((pA * (mid + pB)) & 1) ? -1 : 1; }
            else { int mid = prefix[a] - prefix[b + 1]; sm = ((pA * mid) & 1) ? -1 : 1; }
            vector<gid> nv; nv.reserve(n - 1);
            for (int i = 0; i < n; i++) { if (i == a) continue; nv.push_back(i == b ? ng : gens[i]); }
            emit_term(out, nv, s * sm);
        }
    }
    // type 2: d_ubar d_psi : psi at end of gens[a], ubar at start of gens[b]
    for (int a = 0; a < n; a++) {
        if (!END_PSI(tt[a])) continue;
        i64 s = ((prefix[a] + pars[a] - 1) & 1) ? -1 : 1;
        for (int b = 0; b < n; b++) {
            if (!START_UB(tt[b]) || fb[a] != fa[b]) continue;
            if (only >= 0 && only != a && only != b) continue;
            if (a == b) {
                if (ww[a] == 0) {
                    if (N) { vector<gid> rest; for (int i = 0; i < n; i++) if (i != a) rest.push_back(gens[i]); emit_term(out, rest, s * N); }
                } else {
                    int sc = WCS[ww[a]];
                    if (sc) { vector<gid> nv = gens; nv[a] = (gid)WCI[ww[a]]; emit_term(out, nv, s * sc); }
                }
                continue;
            }
            int pA = 1 - pars[a], pB = pars[b];
            int nw = concat(ww[a], ww[b]);
            if (nw < 0) { fprintf(stderr, "concat missing\n"); exit(2); }
            gid ng = mkgen(mtype(START_PB(tt[a]), END_PSI(tt[b])), nw, fa[a], fb[b]);
            i64 sm;
            if (b > a) { int mid = prefix[b] - prefix[a + 1]; sm = ((pB * mid) & 1) ? -1 : 1; }
            else { int mid = prefix[a] - prefix[b + 1]; sm = ((pB * (mid + pA)) & 1) ? -1 : 1; }
            vector<gid> nv; nv.reserve(n - 1);
            for (int i = 0; i < n; i++) { if (i == b) continue; nv.push_back(i == a ? ng : gens[i]); }
            emit_term(out, nv, s * sm);
        }
    }
}

static void rho(gid beta, const vector<gid> &mono, Terms &out) {
    vector<gid> gens; gens.reserve(mono.size() + 1);
    gens.push_back(beta); gens.insert(gens.end(), mono.begin(), mono.end());
    Terms raw; Dop(gens, 0, 0, raw);
    for (auto &t : raw) {
        bool bad = false;
        for (auto g : t.first) if (isbmode(g)) { bad = true; break; }
        if (!bad) out.push_back(std::move(t));
    }
}

static void dop(const vector<gid> &mono, Terms &out) { Qop(mono, out); Dop(mono, NN, -1, out); }

// ------------------------------------------------------------------ elimination mod p
struct Span {
    u64 p; size_t ncols;
    vector<int64_t> pivstart; vector<uint32_t> plen;    // per column: row start in pool, length
    vector<uint32_t> pcol; vector<uint32_t> pval;       // pool (leading entry included, value 1)
    vector<u64> acc; vector<uint32_t> heap; vector<char> inheap;
    size_t rank = 0; u64 nops = 0;
    Span(u64 p_, size_t nc) : p(p_), ncols(nc), pivstart(nc, -1), plen(nc, 0), acc(nc, 0), inheap(nc, 0) {}
    void push(uint32_t c) { if (!inheap[c]) { inheap[c] = 1; heap.push_back(c); push_heap(heap.begin(), heap.end(), greater<uint32_t>()); } }
    uint32_t pop() { pop_heap(heap.begin(), heap.end(), greater<uint32_t>()); uint32_t c = heap.back(); heap.pop_back(); inheap[c] = 0; return c; }
    static u64 inv(u64 a, u64 p) { u64 r = 1, e = p - 2; a %= p; while (e) { if (e & 1) r = (__uint128_t)r * a % p; a = (__uint128_t)a * a % p; e >>= 1; } return r; }
    // normal form modulo the span: eliminate every pivot column, keep the others (returned sorted)
    void nf(const vector<pair<uint32_t, i64>> &row, vector<pair<uint32_t, u64>> &out) {
        out.clear();
        for (auto &e : row) {
            i64 v = e.second % (i64)p; if (v < 0) v += p;
            if (!v) continue;
            acc[e.first] = (acc[e.first] + (u64)v) % p;
            push(e.first);
        }
        while (!heap.empty()) {
            uint32_t c = pop();
            u64 a = acc[c];
            if (!a) continue;
            acc[c] = 0;
            if (pivstart[c] >= 0) {
                u64 f = p - a;
                int64_t st = pivstart[c]; uint32_t L = plen[c];
                nops += L;
                for (uint32_t j = 1; j < L; j++) {
                    uint32_t cc = pcol[st + j];
                    acc[cc] = (acc[cc] + f * (u64)pval[st + j]) % p;
                    push(cc);
                }
            } else out.push_back({c, a});
        }
    }
    bool add(const vector<pair<uint32_t, i64>> &row) {
        for (auto &e : row) {
            i64 v = e.second % (i64)p; if (v < 0) v += p;
            if (!v) continue;
            acc[e.first] = (acc[e.first] + (u64)v) % p;
            push(e.first);
        }
        while (!heap.empty()) {
            uint32_t c = pop();
            u64 a = acc[c];
            if (!a) continue;
            if (pivstart[c] >= 0) {
                u64 f = p - a; acc[c] = 0;
                int64_t st = pivstart[c]; uint32_t L = plen[c];
                nops += L;
                for (uint32_t j = 1; j < L; j++) {
                    uint32_t cc = pcol[st + j];
                    acc[cc] = (acc[cc] + f * (u64)pval[st + j]) % p;
                    push(cc);
                }
                continue;
            }
            // new pivot at c
            u64 ia = inv(a, p);
            int64_t st = (int64_t)pcol.size();
            pcol.push_back(c); pval.push_back(1); acc[c] = 0;
            vector<uint32_t> rest;
            while (!heap.empty()) { uint32_t cc = pop(); if (acc[cc]) rest.push_back(cc); }
            sort(rest.begin(), rest.end());
            for (auto cc : rest) { pcol.push_back(cc); pval.push_back((uint32_t)((acc[cc] * ia) % p)); acc[cc] = 0; }
            pivstart[c] = st; plen[c] = (uint32_t)(pcol.size() - st);
            rank++;
            return true;
        }
        return false;
    }
};

static void terms_to_row(const Terms &tm, const Block &blk, vector<pair<uint32_t, i64>> &row) {
    row.clear();
    for (auto &t : tm) {
        int64_t i = blk.find(t.first);
        if (i < 0) { fprintf(stderr, "term not in target block (len %zu)\n", t.first.size()); exit(3); }
        row.push_back({(uint32_t)i, t.second});
    }
    sort(row.begin(), row.end());
    size_t m = 0;
    for (size_t j = 0; j < row.size(); j++) {
        if (m && row[m - 1].first == row[j].first) row[m - 1].second += row[j].second;
        else row[m++] = row[j];
    }
    row.resize(m);
    size_t q = 0;
    for (size_t j = 0; j < row.size(); j++) if (row[j].second % (i64)PR) row[q++] = row[j];
    row.resize(q);
}

static double now() { return chrono::duration<double>(chrono::steady_clock::now().time_since_epoch()).count(); }

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: bgk input [progress_every]\n"); return 1; }
    if (argc > 2) PROGRESS = (size_t)atoll(argv[2]);
    FILE *f = fopen(argv[1], "r");
    if (!f) { perror("open"); return 1; }
    long long pp;
    if (fscanf(f, "%d %d %lld", &K, &NN, &pp) != 3) return 1;
    PR = (u64)pp;
    for (int i = 0; i < 5; i++) if (fscanf(f, "%d", &LAM[i]) != 1) return 1;
    if (fscanf(f, "%d", &W) != 1) return 1;
    WL.resize(W); WWT.resize(W); WPAR.resize(W); WGH.resize(W); WQ.resize(W); WCS.resize(W); WCI.resize(W);
    for (int i = 0; i < W; i++) {
        int len; if (fscanf(f, "%d", &len) != 1) return 1;
        WL[i].resize(len);
        for (int j = 0; j < len; j++) { int x; if (fscanf(f, "%d", &x) != 1) return 1; WL[i][j] = (uint8_t)x; }
        string s(WL[i].begin(), WL[i].end());
        WIDX[s] = i;
    }
    for (int i = 0; i < W; i++) { for (int j = 0; j < 5; j++) if (fscanf(f, "%d", &WWT[i][j]) != 1) return 1; if (fscanf(f, "%d %d", &WPAR[i], &WGH[i]) != 2) return 1; }
    for (int i = 0; i < W; i++) { int n; if (fscanf(f, "%d", &n) != 1) return 1; for (int j = 0; j < n; j++) { int w; long long c; if (fscanf(f, "%d %lld", &w, &c) != 2) return 1; WQ[i].push_back({w, c}); } }
    for (int i = 0; i < W; i++) if (fscanf(f, "%d %d", &WCS[i], &WCI[i]) != 2) return 1;
    WMES.assign(W, 1);
    for (int i = 0; i < W; i++) if (fscanf(f, "%d", &WMES[i]) != 1) return 1;
    if (fscanf(f, "%d", &C) != 1) return 1;
    CW.resize(C); CQ.resize(C);
    for (int i = 0; i < C; i++) if (fscanf(f, "%d", &CW[i]) != 1) return 1;
    for (int i = 0; i < C; i++) { int n; if (fscanf(f, "%d", &n) != 1) return 1; for (int j = 0; j < n; j++) { int c2; long long cc; if (fscanf(f, "%d %lld", &c2, &cc) != 2) return 1; CQ[i].push_back({c2, cc}); } }
    if (fscanf(f, "%d", &PRODMODE) != 1) return 1;
    if (PRODMODE) for (int i = 1; i < W; i++) if (WMES[i]) WBYWT[WWT[i]] = i;
    int mode; if (fscanf(f, "%d", &mode) != 1) return 1;
    int ng; if (fscanf(f, "%d", &ng) != 1) return 1;
    vector<int> GS(ng); for (int i = 0; i < ng; i++) if (fscanf(f, "%d", &GS[i]) != 1) return 1;
    fclose(f);
    // WCI holds cyclic index -> generator id is cidx itself
    double t0 = now();
    // enumeration generators: T's and mesons with nonempty words, grouped by weight class
    vector<pair<array<int,5>, gid>> lst;
    for (int c = 0; c < C; c++) { auto w = WWT[CW[c]]; bool ok = true; for (int j = 0; j < 5; j++) if (w[j] > LAM[j]) ok = false; if (ok) lst.push_back({w, (gid)c}); }
    for (int t = M_; t <= R_; t++) for (int w = 1; w < W; w++) {
        if (!WMES[w]) continue;
        bool ok = true; for (int j = 0; j < 5; j++) if (WWT[w][j] > LAM[j]) ok = false;
        if (!ok) continue;
        for (int a = 0; a < K; a++) for (int b = 0; b < K; b++) lst.push_back({WWT[w], mkgen(t, w, a, b)});
    }
    sort(lst.begin(), lst.end());
    {
        int cls = -1; array<int,5> last = {-1, -1, -1, -1, -1};
        for (size_t i = 0; i < lst.size(); i++) {
            if (lst[i].first != last) { cls++; CLS_WT.push_back(lst[i].first); last = lst[i].first; if (cls > 0) CLS_END.push_back((int)i); }
            EG.push_back(lst[i].second); ECLS.push_back(cls);
            EGH.push_back(ggh(lst[i].second)); EPAR.push_back(gpar(lst[i].second));
            vector<int> ch(2 * K, 0); gcharge(lst[i].second, ch.data()); ECH.push_back(ch);
        }
        CLS_END.push_back((int)lst.size());
    }
    // b-operators
    struct BOp { gid beta; vector<int> ch; int dg; };
    vector<BOp> ops;
    auto addop = [&](int t, int x, int y, int dg) { BOp o; o.beta = mkgen(t, 0, x, y); o.ch.assign(2 * K, 0); gcharge(o.beta, o.ch.data()); o.dg = dg; ops.push_back(o); };
    if (mode == 0) {
        for (int x = 0; x < K; x++) for (int y = 0; y < K; y++) {
            if (x != y) { addop(P_, x, y, 0); addop(MM_, x, y, 0); }
            addop(R_, x, y, -1);
        }
    } else {
        for (int x = 0; x + 1 < K; x++) { addop(P_, x, x + 1, 0); addop(MM_, x, x + 1, 0); }
        addop(R_, 0, K - 1, -1);
    }
    vector<int> Tch(2 * K, NN);
    // blocks: for each column space g: V_{T,g} (cols), V_{T,g-1} (d rows), sources for S_g
    // collect all (charge, gh) blocks
    map<pair<vector<int>, int>, Block*> blocks;
    auto need = [&](const vector<int> &ch, int g) { auto key = make_pair(ch, g); if (!blocks.count(key)) blocks[key] = new Block(); return blocks[key]; };
    for (int g : GS) {
        need(Tch, g); need(Tch, g - 1);
        for (auto &o : ops) { vector<int> sc(2 * K); for (int a = 0; a < 2 * K; a++) sc[a] = Tch[a] - o.ch[a]; need(sc, g - o.dg); }
    }
    GHMIN_T = 1 << 30; GHMAX_T = -(1 << 30);
    for (auto &kv : blocks) { TARGETS.push_back({kv.first.first, kv.first.second, kv.second}); GHMIN_T = min(GHMIN_T, kv.first.second); GHMAX_T = max(GHMAX_T, kv.first.second); }
    {
        int rem[5]; for (int j = 0; j < 5; j++) rem[j] = LAM[j];
        int ch[6] = {0, 0, 0, 0, 0, 0};
        vector<gid> cur;
        dfs(0, rem, ch, 0, cur);
    }
    for (auto &kv : blocks) kv.second->build();
    fprintf(stdout, "ENUM %.2fs leaves %lld gens %zu\n", now() - t0, NLEAF, EG.size());
    for (auto &kv : blocks) {
        fprintf(stdout, "BLOCK ch");
        for (int x : kv.first.first) fprintf(stdout, " %d", x);
        fprintf(stdout, " g %d dim %zu\n", kv.first.second, kv.second->n);
    }
    fflush(stdout);
    vector<gid> mono; Terms tm; vector<pair<uint32_t, i64>> row;
    sort(GS.begin(), GS.end());
    if (getenv("BGK_DIAG")) {
        for (int g : GS) {
            Block *col = blocks[make_pair(Tch, g)];
            map<vector<gid>, size_t> cnt;
            for (size_t i = 0; i < col->n; i++) {
                col->get(i, mono);
                vector<gid> key; vector<gid> wds;
                for (auto x : mono) { if (gtype(x) == T_) key.push_back(x); else wds.push_back((gid)gword(x)); }
                sort(wds.begin(), wds.end()); key.push_back(0xFFFFFFFFu); key.insert(key.end(), wds.begin(), wds.end());
                cnt[key]++;
            }
            vector<size_t> sz; for (auto &kv : cnt) sz.push_back(kv.second);
            sort(sz.rbegin(), sz.rend());
            double s2 = 0; for (auto x : sz) s2 += (double)x * x;
            fprintf(stdout, "DIAG g %d blocks %zu max %zu %zu %zu %zu sum_sq %.3g\n", g, sz.size(), sz.size() > 0 ? sz[0] : 0, sz.size() > 1 ? sz[1] : 0, sz.size() > 2 ? sz[2] : 0, sz.size() > 3 ? sz[3] : 0, s2);
        }
        fflush(stdout);
        if (getenv("BGK_DIAG_ONLY")) return 0;
    }
    map<int, vector<size_t>> NP;     // non-pivot columns of S_g in V_{T,g}: a basis of V_{T,g}/S_g
    bool TRANS = getenv("BGK_TRANSPOSE") != nullptr;
    bool BLK = getenv("BGK_NOBLOCK") == nullptr && !TRANS;
    auto skel_key = [&](const vector<gid> &m) {
        vector<gid> key, wds;
        for (auto x : m) { if (gtype(x) == T_) key.push_back(x); else wds.push_back((gid)gword(x)); }
        sort(wds.begin(), wds.end()); key.push_back(0xFFFFFFFFu); key.insert(key.end(), wds.begin(), wds.end());
        return key;
    };
    for (int g : GS) {
      if (BLK) {
        double t1 = now();
        Block *col = blocks[make_pair(Tch, g)];
        Block *rows = blocks[make_pair(Tch, g - 1)];
        // skeleton groups of the column space
        map<vector<gid>, uint32_t> sid; vector<vector<uint32_t>> grp;
        vector<uint32_t> csk(col->n), cloc(col->n);
        for (size_t i = 0; i < col->n; i++) {
            col->get(i, mono); auto key = skel_key(mono);
            auto it = sid.find(key);
            uint32_t k_;
            if (it == sid.end()) { k_ = (uint32_t)grp.size(); sid[key] = k_; grp.push_back({}); } else k_ = it->second;
            csk[i] = k_; cloc[i] = (uint32_t)grp[k_].size(); grp[k_].push_back((uint32_t)i);
        }
        // source groups per operator
        vector<vector<vector<uint32_t>>> sgrp(ops.size(), vector<vector<uint32_t>>(grp.size()));
        vector<Block*> srcb(ops.size());
        for (size_t oi = 0; oi < ops.size(); oi++) {
            auto &o = ops[oi];
            vector<int> sc(2 * K); for (int a = 0; a < 2 * K; a++) sc[a] = Tch[a] - o.ch[a];
            Block *src = blocks[make_pair(sc, g - o.dg)]; srcb[oi] = src;
            for (size_t i = 0; i < src->n; i++) {
                src->get(i, mono); auto key = skel_key(mono);
                auto it = sid.find(key);
                if (it == sid.end()) continue;      // its images vanish or leave V_T: checked below
                sgrp[oi][it->second].push_back((uint32_t)i);
            }
        }
        // d rows from a basis of V_{g-1}/S_{g-1}, bucketed by skeleton of their columns
        vector<size_t> drows;
        bool dskip = false;
        if (NP.count(g - 1)) drows = NP[g - 1]; else if (getenv("BGK_FULLD")) { for (size_t i = 0; i < rows->n; i++) drows.push_back(i); } else dskip = true;
        vector<vector<pair<uint32_t, i64>>> drow(drows.size());
        vector<vector<pair<uint32_t, uint32_t>>> dbucket(grp.size());   // (drow index, position start) - store per skeleton list of drow ids
        vector<vector<uint32_t>> dby(grp.size());
        for (size_t r = 0; r < drows.size(); r++) {
            rows->get(drows[r], mono); tm.clear(); dop(mono, tm);
            terms_to_row(tm, *col, drow[r]);
            uint32_t last = 0xFFFFFFFFu;
            vector<uint32_t> seen;
            for (auto &e : drow[r]) seen.push_back(csk[e.first]);
            sort(seen.begin(), seen.end()); seen.erase(unique(seen.begin(), seen.end()), seen.end());
            for (auto k_ : seen) dby[k_].push_back((uint32_t)r);
            (void)last;
        }
        // global NP numbering and reduced d rows
        vector<size_t> np; vector<int64_t> npidx(col->n, -1);
        vector<vector<pair<uint32_t, i64>>> dred(drows.size());
        size_t rS = 0, nrowS = 0; double tmax = 0; size_t bmax = 0; u64 opsum = 0;
        vector<pair<uint32_t, i64>> lrow; vector<pair<uint32_t, u64>> nfo;
        for (size_t k_ = 0; k_ < grp.size(); k_++) {
            double tb = now();
            size_t nl = grp[k_].size();
            Span sp(PR, nl);
            for (size_t oi = 0; oi < ops.size(); oi++) {
                for (auto i : sgrp[oi][k_]) {
                    srcb[oi]->get(i, mono); tm.clear(); rho(ops[oi].beta, mono, tm);
                    terms_to_row(tm, *col, row);
                    if (row.empty()) continue;
                    lrow.clear();
                    for (auto &e : row) {
                        if (csk[e.first] != k_) { fprintf(stderr, "S-row leaves its skeleton block\n"); exit(4); }
                        lrow.push_back({cloc[e.first], e.second});
                    }
                    sort(lrow.begin(), lrow.end());
                    sp.add(lrow); nrowS++;
                }
            }
            rS += sp.rank;
            for (size_t c = 0; c < nl; c++) if (sp.pivstart[c] < 0) { npidx[grp[k_][c]] = (int64_t)np.size(); np.push_back(grp[k_][c]); }
            for (auto r : dby[k_]) {
                lrow.clear();
                for (auto &e : drow[r]) if (csk[e.first] == k_) lrow.push_back({cloc[e.first], e.second});
                sort(lrow.begin(), lrow.end());
                sp.nf(lrow, nfo);
                for (auto &e : nfo) dred[r].push_back({(uint32_t)npidx[grp[k_][e.first]], (i64)e.second});
            }
            opsum += sp.nops;
            double dt = now() - tb; if (dt > tmax) { tmax = dt; bmax = nl; }
            if (getenv("BGK_BLOCKLOG") && nl >= 2000) { fprintf(stdout, "  block %zu size %zu rows %zu rank %zu pool %zu ops %llu %.2fs\n", k_, nl, (size_t)0, sp.rank, sp.pcol.size(), (unsigned long long)sp.nops, dt); fflush(stdout); }
        }
        // sources not matched to a column skeleton must have zero images in V_T (rho preserves skeletons)
        double t2 = now();
        Span spd(PR, np.size());
        for (auto &r : dred) { sort(r.begin(), r.end()); spd.add(r); }
        size_t rd = spd.rank;
        NP[g] = np;
        fprintf(stdout, "COL g %d dimV %zu rankS %zu rankd_prev %lld  (S rows %zu, d rows %zu, blocks %zu, max block %zu took %.2fs, ops %llu, %.2fs + %.2fs) [B]\n", g, col->n, rS, dskip ? -1LL : (long long)rd, nrowS, drows.size(), grp.size(), bmax, tmax, (unsigned long long)opsum, t2 - t1, now() - t2);
        fflush(stdout);
        continue;
      }
        double t1 = now();
        Block *col = blocks[make_pair(Tch, g)];
        Block *rows = blocks[make_pair(Tch, g - 1)];
        if (TRANS) {
            // rows of [S; D] in CSR, then column-wise incremental rank (rank of the transpose)
            vector<size_t> rst; vector<uint32_t> rc; vector<i64> rv; rst.push_back(0);
            size_t nrowS = 0;
            for (auto &o : ops) {
                vector<int> sc(2 * K); for (int a = 0; a < 2 * K; a++) sc[a] = Tch[a] - o.ch[a];
                Block *src = blocks[make_pair(sc, g - o.dg)];
                for (size_t i = 0; i < src->n; i++) {
                    src->get(i, mono); tm.clear(); rho(o.beta, mono, tm);
                    terms_to_row(tm, *col, row);
                    if (row.empty()) continue;
                    for (auto &e : row) { rc.push_back(e.first); rv.push_back(e.second); }
                    rst.push_back(rc.size()); nrowS++;
                }
            }
            size_t nrowd = 0;
            vector<size_t> drows;
            if (NP.count(g - 1)) drows = NP[g - 1]; else { for (size_t i = 0; i < rows->n; i++) drows.push_back(i); }
            for (size_t i : drows) {
                rows->get(i, mono); tm.clear(); dop(mono, tm);
                terms_to_row(tm, *col, row);
                for (auto &e : row) { rc.push_back(e.first); rv.push_back(e.second); }
                rst.push_back(rc.size()); nrowd++;
            }
            size_t nr = nrowS + nrowd;
            // CSC
            vector<size_t> cst(col->n + 1, 0);
            for (auto c : rc) cst[c + 1]++;
            for (size_t c = 0; c < col->n; c++) cst[c + 1] += cst[c];
            vector<uint32_t> cr(rc.size()); vector<i64> cv(rc.size());
            { vector<size_t> pos(cst.begin(), cst.end() - 1);
              for (size_t r = 0; r < nr; r++) for (size_t j = rst[r]; j < rst[r + 1]; j++) { size_t q = pos[rc[j]]++; cr[q] = (uint32_t)r; cv[q] = rv[j]; } }
            vector<size_t>().swap(rst); vector<uint32_t>().swap(rc); vector<i64>().swap(rv);
            double t2 = now();
            vector<size_t> np; size_t rS = 0;
            { Span sp(PR, nrowS);
              vector<pair<uint32_t, i64>> v;
              for (size_t c = 0; c < col->n; c++) {
                  v.clear();
                  for (size_t j = cst[c]; j < cst[c + 1]; j++) if (cr[j] < nrowS) v.push_back({cr[j], cv[j]});
                  if (!sp.add(v)) np.push_back(c);
              }
              rS = sp.rank; }
            double t3 = now();
            size_t rSD = 0;
            { Span sp(PR, nr);
              vector<pair<uint32_t, i64>> v;
              for (size_t c = 0; c < col->n; c++) {
                  v.clear();
                  for (size_t j = cst[c]; j < cst[c + 1]; j++) v.push_back({cr[j], cv[j]});
                  sp.add(v);
              }
              rSD = sp.rank; }
            NP[g] = np;
            fprintf(stdout, "COL g %d dimV %zu rankS %zu rankd_prev %zu  (S rows %zu, d rows %zu, build %.2fs, S %.2fs, SD %.2fs) [T]\n", g, col->n, rS, rSD - rS, nrowS, nrowd, t2 - t1, t3 - t2, now() - t3);
            fflush(stdout);
            continue;
        }
        Span sp(PR, col->n);
        size_t nrowS = 0;
        for (auto &o : ops) {
            vector<int> sc(2 * K); for (int a = 0; a < 2 * K; a++) sc[a] = Tch[a] - o.ch[a];
            Block *src = blocks[make_pair(sc, g - o.dg)];
            for (size_t i = 0; i < src->n; i++) {
                src->get(i, mono); tm.clear(); rho(o.beta, mono, tm);
                terms_to_row(tm, *col, row);
                if (!row.empty()) { sp.add(row); nrowS++; }
                if (PROGRESS && (nrowS % PROGRESS) == 0) { fprintf(stdout, "  prog g %d rows %zu rank %zu pool %zu ops %llu %.1fs\n", g, nrowS, sp.rank, sp.pcol.size(), (unsigned long long)sp.nops, now() - t1); fflush(stdout); }
            }
        }
        size_t rS = sp.rank;
        vector<size_t> np;
        for (size_t c = 0; c < col->n; c++) if (sp.pivstart[c] < 0) np.push_back(c);
        double t2 = now();
        size_t nrowd = 0;
        if (NP.count(g - 1)) {
            for (size_t i : NP[g - 1]) {
                rows->get(i, mono); tm.clear(); dop(mono, tm);
                terms_to_row(tm, *col, row);
                if (!row.empty()) sp.add(row);
                nrowd++;
            }
        } else {
            for (size_t i = 0; i < rows->n; i++) {
                rows->get(i, mono); tm.clear(); dop(mono, tm);
                terms_to_row(tm, *col, row);
                if (!row.empty()) sp.add(row);
                nrowd++;
            }
        }
        size_t rd = sp.rank - rS;
        NP[g] = np;
        fprintf(stdout, "COL g %d dimV %zu rankS %zu rankd_prev %zu  (S rows %zu, d rows %zu, %.2fs + %.2fs, pool %zu, ops %llu)\n", g, col->n, rS, rd, nrowS, nrowd, t2 - t1, now() - t2, sp.pcol.size(), (unsigned long long)sp.nops);
        fflush(stdout);
    }
    return 0;
}
