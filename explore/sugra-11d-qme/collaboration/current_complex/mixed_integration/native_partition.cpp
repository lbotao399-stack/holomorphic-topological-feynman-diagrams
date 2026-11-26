// Finite-field path sum for all polarizations on a fixed path partition.
// Dense real matrix multiplication is exact here: every dot product is a
// sum of at most 120 nonnegative integer products of residues below 1009.
#define TWIST_PARTITION_LIBRARY
#include "partition_sum.cpp"
#define ACCELERATE_NEW_LAPACK
#include <Accelerate/Accelerate.h>
#include <numeric>
namespace native {
#ifndef TWIST_RESIDUE_PRIME
#define TWIST_RESIDUE_PRIME 1009
#endif
constexpr int prime = TWIST_RESIDUE_PRIME;
static_assert(prime > 37 && 120ULL * (prime - 1) * (prime - 1) < (1ULL << 53));
using Entry = conditional_t<(prime < 65536), uint16_t, uint32_t>;
int inverse_prime(int v) {
  long long a = v, b = 1;
  for (int n = prime - 2; n; n >>= 1) {
    if (n & 1)
      b = b * a % prime;
    a = a * a % prime;
  }
  return int(b);
}
#ifdef UNCUT_CURRENT
constexpr int path_bound = 12, mean_bound = 12;
#else
constexpr int path_bound = 4, mean_bound = 8;
#endif
const size_t mean_count =
    upper_bound(indices.code.begin(), indices.code.end(), mean_bound,
                [](int d, int key) { return d < degree(key); }) -
    indices.code.begin();
int red(long long x) { return int((x % prime + prime) % prime); }
using Mat = array<long long, 25>;
using Poly = map<int, Mat>;
using Fields = array<V, 5>;
struct Group {
  int cell, power, shift;
  array<Entry, 120> coeff{};
};
struct Rows {
  array<map<pair<int, int>, int>, 25> ids;
  vector<Group> groups;
  void add(const Poly &matrix, int order, int w, int ell = 0, int shift = 0) {
    for (const auto &[power, m] : matrix)
      if (degree(power + ell) + degree(shift) <= path_bound)
        for (int i = 0; i < 25; ++i)
          if (m[i]) {
            auto [it, fresh] = ids[i].try_emplace(make_pair(power + ell, shift),
                                                  int(groups.size()));
            if (fresh)
              groups.push_back({i, power + ell, shift, {}});
            auto &c = groups[it->second].coeff[order];
            c = red(c + static_cast<long long>(w) * m[i]);
          }
  }
};
int dot(const V &a, const V &b) {
  return inner_product(a.begin(), a.end(), b.begin(), 0);
}
V plus(V a, const V &b) {
  for (int i = 0; i < 5; ++i)
    a[i] += b[i];
  return a;
}
void scaled(Mat &z, const Mat &x, int c) {
  if (c)
    for (int i = 0; i < 25; ++i)
      z[i] += c * x[i];
}
void right_product(Mat &z, const Mat &a, const Mat &b) {
  for (int i = 0; i < 5; ++i)
    for (int k = 0; k < 5; ++k)
      if (a[5 * i + k])
        for (int j = 0; j < 5; ++j)
          z[5 * i + j] += a[5 * i + k] * b[5 * k + j];
}
int parity(const V &a) {
  int s = 1;
  for (int i = 0; i < 5; ++i)
    for (int j = 0; j < i; ++j)
      if (a[j] > a[i])
        s = -s;
  return s;
}
struct Waves {
  array<int, 5> labels;
  Fields p;
  vector<Fields> polar;
  vector<double> moment;
  vector<array<int, 5>> orders;
  int active = 5;
  void make_moments(int limit) {
    array<int, 5> perm{0, 1, 2, 3, 4};
    if (limit == 120) {
      auto it = stable_partition(perm.begin(), perm.end(), [&](int j) {
        return labels[j] == 9 || p[j] != V{};
      });
      active = int(it - perm.begin());
    }
    do {
      orders.push_back(perm);
      if (int(orders.size()) == limit)
        break;
    } while (next_permutation(perm.begin(), perm.begin() + active));
    size_t m = mean_count;
    moment.resize(m * orders.size());
    array<int, 13> scale{};
    for (int d = 0; d <= 12; ++d) {
      int v = red(fact(active + d));
      scale[d] = red(fact(5) * inverse_prime(v));
    }
    map<Fields, int> used;
    for (size_t n = 0; n < orders.size(); ++n) {
      Fields gaps{};
      V sum{};
      for (int s = 0; s < active; ++s) {
        sum = plus(sum, p[orders[n][s]]);
        gaps[s] = sum;
      }
      double *target = moment.data() + n * m;
      auto [it, fresh] = used.try_emplace(gaps, int(n));
      if (!fresh) {
        copy_n(moment.data() + it->second * m, m, target);
        continue;
      }
      vector<int> h(m);
      h[0] = 1;
      for (auto q : gaps)
        for (size_t j = 1; j < m; ++j) {
          int key = indices.code[j];
          long long s = h[j];
          for (int i = 0; i < 5; ++i)
            if (q[i] && ((key >> (4 * i)) & 15))
              s += long(q[i]) * h[indices.id[key - (1 << (4 * i))]];
          h[j] = red(s);
        }
      for (size_t j = 0; j < m; ++j) {
        int k = indices.code[j], a = h[j];
        while (k) {
          a = red(long(a) * fact(k & 15));
          k >>= 4;
        }
        target[j] = red(long(a) * scale[degree(indices.code[j])]);
      }
    }
  }
  map<pair<int, int>, int> zero_factors(int polar_id, int kind,
                                        bool marked) const {
    map<pair<int, int>, int> result = {{{0, 0}, 1}};
    const int sign = (marked ? !kind : bool(kind)) ? -1 : 1;
    for (int slot = active; slot < 5; ++slot) {
      auto u = polar[polar_id][orders.front()[slot]];
      map<pair<int, int>, int> next;
      for (auto [powers, c] : result)
        for (int i = 0; i < 5; ++i)
          if (u[i]) {
            int k = 1 << (4 * i), v = c * sign * u[i];
            next[{powers.first + k, powers.second}] += v;
            next[{powers.first, powers.second + k}] += v;
          }
      result = std::move(next);
    }
    return result;
  }
  Poly chain(int order, int polar_id, int kind) const {
    bool marked = find(labels.begin(), labels.end(), 9) != labels.end();
    // On an ordinary path kind=0 is reverse, 1 is forward.
    // On a marked path kind=0 is source, 1 is sink; before gamma these are
    // reverse and forward respectively, and after gamma the direction flips.
    bool f = bool(kind);
    Poly a;
    for (int i = 0; i < 5; ++i)
      a[0][6 * i] = 1;
    V prefix{};
    for (int slot = 0; slot < active; ++slot) {
      int local = orders[order][slot];
      const auto &q = p[local], &u = polar[polar_id][local];
      bool gamma = labels[local] == 9;
      Mat constant{};
      array<Mat, 5> linear{};
      if (gamma && kind == 1) {
        V ix{0, 1, 2, 3, 4};
        do {
          int i = ix[0], j = ix[1], g = ix[2], k = ix[3], l = ix[4];
          int v = parity(ix) * u[g] * q[l];
          linear[k][5 * i + j] += v;
          constant[5 * i + j] += v * prefix[k];
        } while (next_permutation(ix.begin(), ix.end()));
      } else if (gamma) {
        for (int i = 0; i < 5; ++i)
          for (int j = 0; j < 5; ++j)
            constant[5 * i + j] = u[i] * q[j] - q[i] * u[j];
      } else {
        int trace = (f ? -1 : 1) * dot(u, prefix);
        for (int i = 0; i < 5; ++i)
          for (int j = 0; j < 5; ++j)
            constant[5 * i + j] =
                -(f ? u[i] * q[j] : q[i] * u[j]) + (i == j ? trace : 0);
      }
      Poly b;
      for (const auto &[key, x] : a) {
        right_product(b[key], x, constant);
        if (gamma && kind == 1) {
          for (int i = 0; i < 5; ++i)
            if (any_of(linear[i].begin(), linear[i].end(),
                       [](int v) { return v != 0; }))
              right_product(b[key + (1 << (4 * i))], x, linear[i]);
        } else if (!gamma)
          for (int i = 0; i < 5; ++i)
            if (u[i])
              scaled(b[key + (1 << (4 * i))], x, (f ? -1 : 1) * u[i]);
      }
      for (auto it = b.begin(); it != b.end();) {
        bool changed = false;
        for (auto &x : it->second) {
          x = red(x);
          changed |= x != 0;
        }
        if (!changed)
          it = b.erase(it);
        else
          ++it;
      }
      a = std::move(b);
      prefix = plus(prefix, q);
      if (marked && gamma)
        f = !f;
    }
    return a;
  }
};
struct Value {
  int power, direct_degree;
  const Entry *v;
};
struct Dense {
  array<vector<Value>, 25> by_cell;
  vector<Entry> data;
  Dense(const Rows &rows, const Waves &w) {
    int n = int(w.orders.size()), m = int(mean_count),
        size = int(rows.groups.size());
    map<array<int, 3>, int> unique;
    for (const auto &g : rows.groups)
      unique.try_emplace({g.cell, g.power, degree(g.shift)},
                         int(unique.size()));
    data.resize(unique.size() * m);
    for (auto [pair, slot] : unique)
      by_cell[pair[0]].push_back(
          {pair[1], pair[2], data.data() + size_t(slot) * m});
    map<int, vector<int>> shifts;
    constexpr int batch = 96;
    vector<double> input(batch * n), out(batch * m);
    for (int j = 0; j < size; j += batch) {
      int count = min(batch, size - j);
      for (int i = 0; i < count; ++i)
        for (int k = 0; k < n; ++k)
          input[i * n + k] = rows.groups[j + i].coeff[k];
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, count, m, n, 1.,
                  input.data(), n, w.moment.data(), m, 0., out.data(), m);
      for (int i = 0; i < count; ++i) {
        auto &r = rows.groups[j + i];
        auto *dst = data.data() +
                    size_t(unique.at({r.cell, r.power, degree(r.shift)})) * m;
        if (!r.shift) {
          for (int k = 0; k < m; ++k)
            dst[k] = (dst[k] + static_cast<uint64_t>(out[i * m + k])) % prime;
        } else {
          auto [it, fresh] = shifts.try_emplace(r.shift);
          if (fresh) {
            auto &ids = it->second;
            ids.assign(m, -1);
            int maxd = mean_bound - degree(r.shift);
            for (int k = 0; k < m && degree(indices.code[k]) <= maxd; ++k)
              ids[k] = indices.id[indices.code[k] + r.shift];
          }
          const auto &q = it->second;
          for (int k = 0; k < m && q[k] >= 0; ++k)
            dst[k] =
                (dst[k] + static_cast<uint64_t>(out[i * m + q[k]])) % prime;
        }
      }
    }
  }
};
struct Power {
  int exp, degree, bin;
  const Entry *mom;
};
struct Degree {
  array<long long, path_bound + 1> v{};
  unsigned mask = 0;
};
struct Diff {
  int key, degree;
  vector<Power> term;
  Degree eval(int incoming) const {
    Degree out;
    for (auto t : term)
      out.v[t.degree] += long(t.bin) * t.mom[indices.id[t.exp + incoming]];
    for (int i = 0; i <= path_bound; ++i)
      if (out.v[i]) {
        out.v[i] %= prime;
        if (out.v[i])
          out.mask |= 1u << i;
      }
    return out;
  }
};
vector<Diff> split_path(const vector<Value> &rows, int offset) {
  map<int, vector<Power>> p;
  for (auto r : rows) {
    int total = offset + r.power;
    if (degree(total) + r.direct_degree > path_bound)
      continue;
    for (auto s : splits(total))
      p[s.key].push_back(
          {total - s.key, degree(total) - s.deg, s.binomial, r.v});
  }
  vector<Diff> out;
  for (auto &[k, v] : p)
    out.push_back({k, degree(k), std::move(v)});
  return out;
}
int rational_weight(int bin) {
  static map<int, int> cache;
  auto [it, fresh] = cache.try_emplace(bin);
  if (fresh) {
    int a = bin & 15, b = (bin >> 4) & 15, c = (bin >> 8) & 15, d = bin >> 12;
    mpq_class q =
        weight(bin) * fact(5 + a + d) * fact(5 + b + c) / fact(5) / fact(5);
    mpz_class z, p(prime);
    if (!mpz_invert(z.get_mpz_t(), q.get_den_mpz_t(), p.get_mpz_t()))
      throw logic_error("denominator");
    mpz_class v = q.get_num() * z;
    it->second = mpz_fdiv_ui(v.get_mpz_t(), prime);
  }
  return it->second;
}
#if defined(DIRECT_GLUE) || defined(UNCUT_CURRENT)
int glue(const Dense &b, const Dense &c, const vector<End> &ends) {
  array<map<pair<int, int>, vector<Diff>>, 2> cache;
  vector<__int128_t> bins(1 << 16);
  for (auto e : ends) {
    auto [lb, nb] = cache[0].try_emplace({e.js, e.b});
    auto [rc, nc] = cache[1].try_emplace({e.kt, e.c});
    if (nb)
      lb->second = split_path(b.by_cell[e.js], e.b);
    if (nc)
      rc->second = split_path(c.by_cell[e.kt], e.c);
    for (const auto &l : lb->second)
      for (const auto &r : rc->second) {
        const auto &a = l.eval(r.key);
        if (!a.mask)
          continue;
        const auto &v = r.eval(l.key);
        int base = (l.degree << 4) + (r.degree << 12);
        for (unsigned ai = a.mask; ai; ai &= ai - 1) {
          int d = __builtin_ctz(ai);
          long x = long(e.weight) * a.v[d];
          for (unsigned ci = v.mask; ci; ci &= ci - 1) {
            int k = __builtin_ctz(ci);
            bins[base + d + (k << 8)] += x * v.v[k];
          }
        }
      }
  }
  int total = 0;
  for (int k = 0; k < int(bins.size()); ++k)
    if (bins[k])
      total = red(total + long(red(static_cast<long long>(bins[k] % prime))) *
                              rational_weight(k));
  return total;
}
#else
struct Cuts {
  map<pair<int, int>, int> index;
  vector<vector<Diff>> derivatives;
  vector<int> exponents;
  void add(const Dense &path, int cell, int off) {
    auto [it, fresh] = index.try_emplace({cell, off}, int(index.size()));
    if (fresh)
      derivatives.push_back(split_path(path.by_cell[cell], off));
  }
  void select() {
    vector<bool> seen(indices.code.size());
    for (const auto &row : derivatives)
      for (const auto &d : row) {
        int i = indices.id[d.key];
        if (!seen[i]) {
          exponents.push_back(d.key);
          seen[i] = true;
        }
      }
    sort(exponents.begin(), exponents.end());
  }
};
int glue(const Dense &b, const Dense &c, const vector<End> &ends) {
  Cuts l, r;
  for (auto e : ends) {
    l.add(b, e.js, e.b);
    r.add(c, e.kt, e.c);
  }
  l.select();
  r.select();
  if (l.exponents.empty() || r.exponents.empty())
    return 0;
  vector<int> start(indices.code.size()), right_id(indices.code.size());
  int count = 0;
  for (int x : l.exponents) {
    start[indices.id[x]] = count;
    count += (path_bound - degree(x) + 1) * int(r.exponents.size());
  }
  for (int i = 0; i < int(r.exponents.size()); ++i)
    right_id[indices.id[r.exponents[i]]] = i;
  vector<Entry> a(l.derivatives.size() * count),
      z(r.derivatives.size() * count);
  for (int h = 0; h < int(l.derivatives.size()); ++h)
    for (const auto &d : l.derivatives[h]) {
      int len = path_bound - d.degree + 1;
      auto *target = a.data() + size_t(h) * count + start[indices.id[d.key]];
      for (int j = 0; j < int(r.exponents.size()); ++j) {
        auto v = d.eval(r.exponents[j]);
        for (int k = 0; k < len; ++k)
          target[j * len + k] = v.v[k];
      }
    }
  int w[5][5][5][5] = {};
  for (int x = 0; x <= 4; ++x)
    for (int y = 0; y <= 4; ++y)
      for (int i = 0; i + x <= 4; ++i)
        for (int j = 0; j + y <= 4; ++j)
          w[i][x][j][y] = rational_weight(i + (x << 4) + (j << 8) + (y << 12));
  for (int h = 0; h < int(r.derivatives.size()); ++h)
    for (const auto &d : r.derivatives[h]) {
      int j = right_id[indices.id[d.key]];
      for (int x : l.exponents) {
        int dx = degree(x), len = path_bound - dx + 1;
        auto v = d.eval(x);
        if (!v.mask)
          continue;
        auto *target =
            z.data() + size_t(h) * count + start[indices.id[x]] + j * len;
        for (int i = 0; i < len; ++i) {
          int64_t sum = 0;
          for (int k = 0; k + d.degree <= path_bound; ++k)
            sum += long(w[i][dx][k][d.degree]) * v.v[k];
          target[i] = sum % prime;
        }
      }
    }
  // At degree four there are at most 31752 columns. At the certified maximum
  // prime the scalar product fits in a signed 64-bit integer.
  if ((uint64_t(count) * (prime - 1)) >
      uint64_t(numeric_limits<int64_t>::max()) / (prime - 1))
    throw overflow_error("dense contraction");
  __int128_t answer = 0;
  for (auto e : ends) {
    const auto *left = a.data() + size_t(l.index.at({e.js, e.b})) * count;
    const auto *right = z.data() + size_t(r.index.at({e.kt, e.c})) * count;
    int64_t sum = 0;
    for (int i = 0; i < count; ++i)
      sum += int64_t(left[i]) * right[i];
    answer += __int128_t(e.weight) * sum;
  }
  return red(static_cast<long long>(answer % prime));
}
#endif
struct Geometry {
  int id, det, limit;
  array<Waves, 2> w;
  vector<vector<int>> weights;
  array<vector<End>, 4> ends;
  bool read() {
    if (!(cin >> id >> det >> limit))
      return false;
    if (limit < 1 || limit > 120)
      throw invalid_argument("order count");
    for (auto &path : w) {
      for (auto &l : path.labels)
        cin >> l;
      for (auto &q : path.p)
        for (auto &v : q)
          cin >> v;
      int n;
      cin >> n;
      path.polar.resize(n);
      for (auto &u : path.polar)
        for (auto &a : u)
          for (auto &v : a)
            cin >> v;
    }
    if (find(w[0].labels.begin(), w[0].labels.end(), 9) != w[0].labels.end())
      throw invalid_argument("gamma belongs on C");
    weights.resize(w[0].polar.size());
    for (auto &row : weights) {
      row.resize(w[1].polar.size());
      for (auto &v : row)
        cin >> v;
    }
    for (auto &list : ends) {
      int n;
      cin >> n;
      list.resize(n);
      for (auto &e : list)
        cin >> e.js >> e.kt >> e.b >> e.c >> e.weight;
    }
    if (!cin)
      throw invalid_argument("partial case");
    return true;
  }
  int solve() {
    if (!det)
      return 0;
    for (auto &path : w)
      path.make_moments(limit);
    int classes = int(w[0].polar.size());
    array<vector<array<Rows, 2>>, 2> rows;
    for (auto &h : rows)
      h.resize(classes);
    for (int h = 0; h < 2; ++h) {
      auto &path = w[h];
      vector<array<map<pair<int, int>, int>, 2>> zero(path.polar.size());
      for (int u = 0; u < int(path.polar.size()); ++u)
        for (int t = 0; t < 2; ++t)
          zero[u][t] = path.zero_factors(u, t, bool(h));
      for (int n = 0; n < int(path.orders.size()); ++n) {
        int g = 0;
        if (h) {
          while (path.labels[path.orders[n][g]] != 9)
            ++g;
        }
        int external = h && ((4 - g) & 1) ? -1 : 1;
        for (int u = 0; u < int(path.polar.size()); ++u)
          for (int t = 0; t < 2; ++t) {
            auto m = path.chain(n, u, t);
            for (auto [powers, c] : zero[u][t])
              if (c) {
                if (!h)
                  rows[0][u][t].add(m, n, c, powers.first, powers.second);
                else
                  for (int l = 0; l < classes; ++l)
                    if (weights[l][u])
                      rows[1][l][t].add(m, n, external * weights[l][u] * c,
                                        powers.first, powers.second);
              }
          }
      }
    }
    int result = 0;
    for (int l = 0; l < classes; ++l) {
      Dense br(rows[0][l][0], w[0]), bf(rows[0][l][1], w[0]),
          cs(rows[1][l][0], w[1]), ce(rows[1][l][1], w[1]);
      result = red(result + glue(br, cs, ends[0]) + glue(bf, cs, ends[1]) +
                   glue(br, ce, ends[2]) + glue(bf, ce, ends[3]));
    }
    return red(long(det) * result);
  }
};
} // namespace native
#ifndef TWIST_NATIVE_PARTITION_LIBRARY
int main(int argc, char **argv) {
  if (argc == 2 && string(argv[1]) == "--prime") {
    cout << native::prime << "\n";
    return 0;
  }
  if (argc != 1)
    throw invalid_argument("usage: native_partition [--prime]");
  ios::sync_with_stdio(false);
  cin.tie(nullptr);
  while (true) {
    native::Geometry g;
    if (!g.read())
      break;
    cout << g.id << ' ' << g.solve() << '\n' << flush;
  }
}
#endif
