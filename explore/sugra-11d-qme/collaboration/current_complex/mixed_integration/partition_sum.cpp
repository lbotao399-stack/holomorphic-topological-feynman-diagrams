// The four conductance exponents also determine the two path factorials.
// This keeps the inner contraction integral. GMP is needed only for the final
// bins.
#include <algorithm>
#include <array>
#include <cstdint>
#include <gmpxx.h>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <unordered_map>
#include <vector>
using namespace std;
namespace {
using V = array<int, 5>;
constexpr int max_degree = 12, fields = 1 << 20;
int residue_prime = 0;
#ifdef ARBITRARY_PRECISION
using Wide = mpz_class;
using Coefficient = mpz_class;
Wide checked_mul(const Wide &a, const Wide &b) { return a * b; }
Wide checked_add(const Wide &a, const Wide &b) { return a + b; }
const Wide &exact(const Wide &a) { return a; }
void add_product(Wide &out, const Wide &a, const Wide &b) {
  mpz_addmul(out.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
}
void add_scaled(Wide &out, const Wide &a, long n) {
  if (n > 0)
    mpz_addmul_ui(out.get_mpz_t(), a.get_mpz_t(), n);
  if (n < 0)
    mpz_submul_ui(out.get_mpz_t(), a.get_mpz_t(), -n);
}
#else
using Wide = __int128_t;
using Coefficient = int;
Wide checked_add(Wide a, Wide b) {
  Wide out;
  if (__builtin_add_overflow(a, b, &out))
    throw overflow_error("128-bit sum");
  return out;
}
Wide checked_mul(Wide a, Wide b) {
  Wide out;
  if (__builtin_mul_overflow(a, b, &out))
    throw overflow_error("128-bit product");
  return out;
}
mpz_class exact(Wide a) {
  bool neg = a < 0;
  __uint128_t mag = neg ? -__uint128_t(a) : __uint128_t(a);
  mpz_class z;
  mpz_import(z.get_mpz_t(), 1, -1, sizeof(mag), 0, 0, &mag);
  return neg ? -z : z;
}
void add_product(Wide &out, Wide a, Wide b) {
  out = checked_add(out, checked_mul(a, b));
}
void add_scaled(Wide &out, Wide a, long n) { add_product(out, a, n); }
#endif
long fact(int n) {
  long a = 1;
  for (int i = 2; i <= n; ++i)
    a *= i;
  return a;
}
long choose(int n, int k) { return fact(n) / fact(k) / fact(n - k); }
int degree(int code) {
  int d = 0;
  while (code) {
    d += code & 15;
    code >>= 4;
  }
  return d;
}
struct Index {
  vector<int> code, id = vector<int>(fields, -1);
  Index() {
    for (int d = 0; d <= max_degree; ++d)
      enumerate(0, d, 0);
  }
  void enumerate(int i, int n, int key) {
    if (i == 4) {
      key += n << (4 * i);
      id[key] = int(code.size());
      code.push_back(key);
      return;
    }
    for (int j = 0; j <= n; ++j)
      enumerate(i + 1, n - j, key + (j << (4 * i)));
  }
} indices;
using Prefix = array<V, 5>;
vector<Wide> moments(const Prefix &waves) {
  vector<Wide> h(indices.code.size());
  h[0] = 1;
  V q{};
  for (auto p : waves) {
    for (int i = 0; i < 5; ++i)
      q[i] += p[i];
    for (size_t pos = 1; pos < h.size(); ++pos) {
      int key = indices.code[pos];
      for (int i = 0; i < 5; ++i)
        if ((key >> (4 * i)) & 15) {
          int lower = indices.id[key - (1 << (4 * i))];
          add_scaled(h[pos], h[lower], q[i]);
        }
    }
  }
  for (size_t pos = 0; pos < h.size(); ++pos) {
    int k = indices.code[pos];
    while (k) {
      h[pos] = checked_mul(h[pos], fact(k & 15));
      k >>= 4;
    }
  }
  return h;
}
struct Chain {
  vector<pair<int, Coefficient>> coeff;
  vector<Wide> values;
  vector<bool> ready;
  const Wide &moment(const vector<vector<Wide>> &mom, int a) {
    if (values.empty()) {
      values.resize(indices.code.size());
      ready.resize(indices.code.size());
    }
    auto &v = values[a];
    if (!ready[a]) {
      ready[a] = true;
      for (const auto &[n, c] : coeff)
        add_product(v, mom[n][a], c);
      if (residue_prime)
        v = (v % residue_prime + residue_prime) % residue_prime;
    }
    return v;
  }
};
struct Split {
  int key, deg, binomial;
};
const vector<Split> &splits(int code) {
  static unordered_map<int, vector<Split>> cache;
  auto [it, added] = cache.try_emplace(code);
  auto &result = it->second;
  if (added) {
    result.push_back({0, 0, 1});
    for (int i = 0; i < 5; ++i) {
      int n = (code >> (4 * i)) & 15, m = int(result.size());
      for (int j = 1; j <= n; ++j)
        for (int s = 0; s < m; ++s) {
          auto t = result[s];
          result.push_back({t.key + (j << (4 * i)), t.deg + j,
                            t.binomial * int(choose(n, j))});
        }
    }
  }
  return result;
}
#ifdef CARTESIAN_REFERENCE
struct Term {
  uint16_t b, c, bin, binomial;
};
const vector<Term> &expansion(int x, int y) {
  static unordered_map<uint64_t, vector<Term>> cache;
  auto [it, added] = cache.try_emplace(uint64_t(x) + (uint64_t(y) << 20));
  auto &v = it->second;
  if (added) {
    int dx = degree(x), dy = degree(y);
    if (dx + dy > max_degree)
      throw logic_error("degree");
    const auto &sx = splits(x), &sy = splits(y);
    v.reserve(sx.size() * sy.size());
    for (auto a : sx)
      for (auto b : sy) {
        int ib = indices.id[x + y - a.key - b.key],
            ic = indices.id[a.key + b.key];
        int bin =
            (dx - a.deg) + (a.deg << 4) + (b.deg << 8) + ((dy - b.deg) << 12);
        if (ib < 0 || ic < 0)
          throw logic_error("missing exponent");
        v.push_back({uint16_t(ib), uint16_t(ic), uint16_t(bin),
                     uint16_t(a.binomial * b.binomial)});
      }
  }
  return v;
}
#endif
long odd(int start, int n) {
  long a = 1;
  for (int i = 0; i < n; ++i)
    a *= start + 2 * i;
  return a;
}
mpq_class weight(int key) {
  int a = key & 15, b = (key >> 4) & 15, c = (key >> 8) & 15, d = key >> 12;
  mpq_class w;
  for (int i = 0; i <= a; ++i)
    for (int j = 0; j <= c; ++j) {
      mpz_class num = mpz_class(choose(a, i)) * choose(c, j) * odd(1, i + b) *
                      odd(1, j + d);
      if ((a - i + c - j) & 1)
        num = -num;
      mpq_class m(num, odd(13, i + b + j + d));
      m.canonicalize();
      w += m;
    }
  return w * fact(5) / fact(5 + a + d) * fact(5) / fact(5 + b + c);
}
using Groups = array<map<int, Chain>, 25>;
struct End {
  int js, kt, b, c, weight;
};
// Each derivative collects all monomials allocating the same powers to the
// other path.
struct LocalPower {
  int exponent, degree, binomial;
  Chain *chain;
};
struct ByDegree {
  array<Wide, max_degree + 1> values{};
  unsigned nonzero = 0;
};
struct Derivative {
  int exponent, degree;
  vector<LocalPower> powers;
  ByDegree evaluate(const vector<vector<Wide>> &mom, int incoming) const {
    ByDegree out;
    for (auto t : powers) {
      int i = indices.id[t.exponent + incoming];
      if (i < 0)
        throw logic_error("combined path degree");
      add_scaled(out.values[t.degree], t.chain->moment(mom, i), t.binomial);
    }
    for (int i = 0; i <= max_degree; ++i) {
      if (residue_prime)
        out.values[i] %= residue_prime;
      if (out.values[i])
        out.nonzero |= 1u << i;
    }
    return out;
  }
};
vector<Derivative> derivatives(map<int, Chain> &group, int offset) {
  map<int, vector<LocalPower>> table;
  for (auto &[power, chain] : group) {
    int total = power + offset;
    for (auto s : splits(total)) {
      int left = total - s.key;
      table[s.key].push_back({left, degree(left), s.binomial, &chain});
    }
  }
  vector<Derivative> out;
  for (auto &[key, row] : table)
    out.push_back({key, degree(key), std::move(row)});
  return out;
}
void contract(vector<Wide> &bins, const vector<End> &ends, Groups &b, Groups &c,
              const array<vector<vector<Wide>>, 2> &mom) {
  array<map<pair<int, int>, vector<Derivative>>, 2> cache;
  for (auto e : ends) {
    auto [lb, nb] = cache[0].try_emplace({e.js, e.b});
    auto [rc, nc] = cache[1].try_emplace({e.kt, e.c});
    if (nb)
      lb->second = derivatives(b[e.js], e.b);
    if (nc)
      rc->second = derivatives(c[e.kt], e.c);
    for (const auto &left : lb->second)
      for (const auto &right : rc->second) {
        auto l = left.evaluate(mom[0], right.exponent);
        if (!l.nonzero)
          continue;
        auto r = right.evaluate(mom[1], left.exponent);
        int base = (left.degree << 4) + (right.degree << 12);
        for (unsigned lm = l.nonzero; lm; lm &= lm - 1) {
          int a = __builtin_ctz(lm);
          Wide x = checked_mul(e.weight, l.values[a]);
          for (unsigned rm = r.nonzero; rm; rm &= rm - 1) {
            int c = __builtin_ctz(rm), bin = base + a + (c << 8);
            add_product(bins[bin], x, r.values[c]);
          }
        }
      }
  }
}
} // namespace
#ifndef TWIST_PARTITION_LIBRARY
int main(int argc, char **argv) {
  if (argc > 2)
    throw invalid_argument("usage: partition_sum [prime]");
  if (argc == 2) {
    residue_prime = stoi(argv[1]);
    if (residue_prime <= 37 ||
        !mpz_probab_prime_p(mpz_class(residue_prime).get_mpz_t(), 10))
      throw invalid_argument("prime greater than 37 required");
  }
  ios::sync_with_stdio(false);
  cin.tie(nullptr);
  array<array<Groups, 4>, 2> groups;
  array<vector<vector<Wide>>, 2> mom;
  array<vector<End>, 4> ends;
  int det = 0;
  char type;
  while (cin >> type) {
    if (type == 'D')
      cin >> det;
    else if (type == 'R') {
      int h, n;
      cin >> h >> n;
      Prefix p{};
      for (auto &row : p)
        for (auto &v : row)
          cin >> v;
      if (!cin || h < 1 || h > 2 || int(mom[h - 1].size()) != n)
        throw logic_error("input order");
      mom[h - 1].push_back(moments(p));
    } else if (type == 'P') {
      int h, f, i, j, key, n;
      Coefficient c;
      if (!(cin >> h >> f >> i >> j >> key >> n >> c) || h < 1 || h > 2 ||
          f < 0 || f >= 4 || i < 0 || i >= 5 || j < 0 || j >= 5 || key < 0 ||
          key >= fields || indices.id[key] < 0 || n < 0 ||
          size_t(n) >= mom[h - 1].size())
        throw invalid_argument("invalid chain coefficient");
      groups[h - 1][f][5 * i + j][key].coeff.push_back({n, c});
    } else if (type == 'E') {
      int f, j, s, k, t, b, c, w;
      cin >> f >> j >> s >> k >> t >> b >> c >> w;
      if (!cin || f < 0 || f >= 4 || min({j, s, k, t}) < 0 ||
          max({j, s, k, t}) >= 5 || b < 0 || c < 0 || b >= fields ||
          c >= fields || indices.id[b] < 0 || indices.id[c] < 0)
        throw invalid_argument("invalid endpoint");
      ends[f].push_back({5 * j + s, 5 * k + t, b, c, w});
    } else
      throw logic_error("unknown record");
  }
  if (!cin.eof() || mom[0].empty() || mom[1].empty())
    throw invalid_argument("incomplete input");
  vector<Wide> bins(1 << 16);
  for (int f = 0; f < 4; ++f) {
#ifdef CARTESIAN_REFERENCE
    for (auto e : ends[f])
      for (auto &[x, left] : groups[0][f][e.js])
        for (auto &[y, right] : groups[1][f][e.kt]) {
          for (auto t : expansion(x + e.b, y + e.c)) {
            auto a = left.moment(mom[0], t.b);
            if (!a)
              continue;
            auto b = right.moment(mom[1], t.c);
            auto v =
                checked_mul(checked_mul(a, b), long(e.weight) * t.binomial);
            bins[t.bin] = checked_add(bins[t.bin], v);
          }
        }
#else
    contract(bins, ends[f], groups[0][f], groups[1][f], mom);
#endif
    for (auto &h : groups)
      for (auto &g : h[f])
        g.clear();
  }
  mpq_class answer;
  for (int i = 0; i < int(bins.size()); ++i)
    if (bins[i])
      answer += exact(bins[i]) * weight(i);
  if (residue_prime) {
    mpz_class p(residue_prime), inv;
    if (!mpz_invert(inv.get_mpz_t(), answer.get_den_mpz_t(), p.get_mpz_t()))
      throw logic_error("nonunit");
    mpz_class reduced = (answer.get_num() * inv) % p;
    answer = (reduced + p) % p;
    cout << "orders " << mom[0].size() * mom[1].size() << ", signed Wick "
         << answer << ", determinant-weighted "
         << ((det * answer.get_num() % p + p) % p) << " mod " << p << "\n";
    return 0;
  }
  cout << "orders " << mom[0].size() * mom[1].size() << ", signed Wick "
       << answer << ", determinant-weighted " << det * answer << "\n";
}

#endif
