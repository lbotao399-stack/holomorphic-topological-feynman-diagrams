#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <stdexcept>
#include <unordered_map>
#include <vector>
using namespace std;
#ifndef EXTREME_PRIME
#define EXTREME_PRIME 251
#endif
constexpr int P = EXTREME_PRIME;
int mod(int x) { return (x % P + P) % P; }
int inv(int x) {
  int v = 1;
  for (int k = 0; k < P - 2; k++)
    v = v * x % P;
  return v;
}
using Weight = array<int, 5>;
// Divided-power monomial z^[a] times the vector basis element v.
struct Mono {
  Weight a;
  int v;
  bool operator<(const Mono &b) const { return tie(v, a) < tie(b.v, b.a); }
};
using Poly = map<Mono, int>;
struct Jet {
  vector<Poly> basis;
  vector<Weight> weights;
  map<Mono, int> pos;
  // Among equal-weight monomials with nonzero divergence, subtract the
  // component in direction 4. All resulting basis vectors are integral.
  Jet(int d) {
    vector<Mono> all;
    Weight a{};
    auto rec = [&](auto &&rec, int i, int left) -> void {
      if (i == 4) {
        a[i] = left;
        for (int k = 0; k < 5; k++)
          all.push_back({a, k});
        return;
      }
      for (int j = 0; j <= left; j++) {
        a[i] = j;
        rec(rec, i + 1, left - j);
      }
    };
    rec(rec, 0, d);
    for (const auto &m : all) {
      if (m.a[m.v] && m.v == 4)
        continue;
      Poly p{{m, 1}};
      if (m.a[m.v]) {
        Mono h = m;
        h.a[m.v]--;
        h.a[4]++;
        h.v = 4;
        p[h] = -1;
      }
      pos[m] = basis.size();
      basis.push_back(p);
      Weight w = m.a;
      for (int &x : w)
        x = -x;
      w[m.v]++;
      weights.push_back(w);
    }
  }
  vector<pair<int, int>> coords(const Poly &p) const {
    vector<pair<int, int>> c;
    for (auto [m, x] : p) {
      auto it = pos.find(m);
      if (it != pos.end() && x)
        c.push_back({it->second, x});
    }
    return c;
  }
};
using Word = uint64_t;
using Vec = unordered_map<Word, int>;
using Slots = array<int, 8>;
Word encode(const Slots &a) {
  Word m = 0;
  for (int x : a)
    m = (m << 7) | x;
  return m;
}
Slots decode(Word x) {
  Slots a{};
  for (int k = 7; k >= 0; k--) {
    a[k] = x & 127;
    x >>= 7;
  }
  return a;
}
// The bracket takes a pair of quadratic jets to a cubic jet. Matrix
// transposes below are algebraic dual maps in the chosen divided-power bases.
struct Family {
  Jet low{2}, cubic{3};
  vector<pair<int, int>> action[5][5][70], transposed[5][5][70],
      brackets[70][70];
  vector<array<int, 3>> preimages[160];
  Family() {
    for (int k = 0; k < 70; k++)
      for (int i = 0; i < 5; i++)
        for (int j = 0; j < 5; j++)
          if (i != j) {
            Poly p;
            for (auto [m, c] : low.basis[k]) {
              if (m.v == j) {
                Mono t = m;
                t.v = i;
                p[t] += c;
              }
              if (m.a[i]) {
                Mono t = m;
                t.a[i]--;
                t.a[j]++;
                p[t] -= c * t.a[j];
              }
            }
            action[i][j][k] = low.coords(p);
            for (auto [h, c] : action[i][j][k])
              transposed[i][j][h].push_back({k, c});
          }
    for (int k = 0; k < 70; k++)
      for (int l = k + 1; l < 70; l++) {
        Poly p;
        auto add = [&](int s, int t, int sg) {
          for (auto [m, c] : low.basis[s])
            for (auto [n, d] : low.basis[t]) {
              if (!n.a[m.v])
                continue;
              Mono b = n;
              b.a[m.v]--;
              int mult = c * d * sg;
              for (int h = 0; h < 5; h++) {
                int y = b.a[h], x = m.a[h];
                b.a[h] += x;
                for (int v = 1; v <= y; v++)
                  mult = mult * (x + v) / v;
              }
              p[b] += mult;
            }
        };
        add(k, l, 1);
        add(l, k, -1);
        brackets[k][l] = cubic.coords(p);
        for (auto [b, c] : brackets[k][l])
          preimages[b].push_back({k, l, c});
      }
  }
  template <typename T>
  unordered_map<Word, T> apply(const unordered_map<Word, T> &v, int i, int j,
                               bool transpose = false,
                               bool reduce = true) const {
    const auto &ops = transpose ? transposed[i][j] : action[i][j];
    unordered_map<Word, T> out;
    out.reserve(v.size() * 2);
    for (auto [w, c] : v) {
      auto a = decode(w);
      for (int k = 0; k < 8 && a[k] != 127; k++)
        for (auto [newslot, t] : ops[a[k]]) {
          auto b = a;
          int z = k, s = 1;
          for (; z > 0 && newslot < b[z - 1]; z--) {
            b[z] = b[z - 1];
            s = -s;
          }
          for (; z < 7 && newslot > b[z + 1]; z++) {
            b[z] = b[z + 1];
            s = -s;
          }
          if ((z > 0 && b[z - 1] == newslot) || (z < 7 && b[z + 1] == newslot))
            continue;
          b[z] = newslot;
          out[encode(b)] += s * t * c;
        }
    }
    for (auto it = out.begin(); it != out.end();) {
      if (reduce)
        it->second = mod(it->second % P);
      if (!it->second)
        it = out.erase(it);
      else
        ++it;
    }
    return out;
  }
};
// Enumerate the exterior weight space using coordinatewise suffix bounds.
struct Seeds {
  mt19937 rng{123456};
  vector<Word> seeds;
  Seeds(const Family &f, Weight w, int arity = 8) {
    vector<int> order(70);
    iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(),
         [&](int i, int j) { return f.low.weights[i] > f.low.weights[j]; });
    int lo[71][9][5]{}, hi[71][9][5]{};
    for (int k = 69; k >= 0; k--)
      for (int n = 1; n <= min(8, 70 - k); n++)
        for (int i = 0; i < 5; i++) {
          int v = f.low.weights[order[k]][i];
          lo[k][n][i] = lo[k + 1][n - 1][i] + v;
          hi[k][n][i] = hi[k + 1][n - 1][i] + v;
          if (n < 70 - k) {
            lo[k][n][i] = min(lo[k][n][i], lo[k + 1][n][i]);
            hi[k][n][i] = max(hi[k][n][i], hi[k + 1][n][i]);
          }
        }
    Slots sel;
    sel.fill(127);
    auto walk = [&](auto &&self, int j, int n, Weight rem) -> void {
      if (n == 0) {
        auto a = sel;
        sort(a.begin(), a.end());
        seeds.push_back(encode(a));
        return;
      }
      for (int k = j; k <= 70 - n; k++) {
        Weight next;
        bool ok = true;
        for (int i = 0; i < 5; i++) {
          next[i] = rem[i] - f.low.weights[order[k]][i];
          if (next[i] < lo[k + 1][n - 1][i] || next[i] > hi[k + 1][n - 1][i]) {
            ok = false;
            break;
          }
        }
        if (ok) {
          sel[arity - n] = order[k];
          self(self, k + 1, n - 1, next);
        }
      }
    };
    walk(walk, 0, arity, w);
    shuffle(seeds.begin(), seeds.end(), rng);
  }
};
// Ordered root factors of the extremal projector. Each sample below also
// passes an explicit annihilation check for all four simple raising maps.
Vec project(Family &f, Weight w, Vec v) {
  for (int j = 1; j < 5; j++)
    for (int i = 0; i < j; i++) {
      Vec up = v;
      int coeff = 1;
      for (int n = 1;; n++) {
        up = f.apply(up, i, j);
        if (up.empty())
          break;
        Vec down = up;
        for (int k = 0; k < n; k++)
          down = f.apply(down, j, i);
        coeff = mod(-coeff * inv(mod(n * (w[i] - w[j] + j - i + n))));
        for (auto [m, c] : down) {
          int t = mod(v[m] + coeff * c);
          if (t)
            v[m] = t;
          else
            v.erase(m);
        }
      }
    }
  return v;
}
struct Probe {
  vector<pair<Word, int>> entries;
  int eval(const Vec &v) const {
    int64_t x = 0;
    for (auto [m, c] : entries) {
      auto t = v.find(m);
      if (t != v.end())
        x += c * t->second;
    }
    return mod(x % P);
  }
};
struct Rank {
  vector<vector<int>> pivots;
  int rank = 0;
  unsigned digest = 0;
  Rank(int width) : pivots(width) {}
  bool insert(vector<int> x) {
    for (int k = 0; k < (int)x.size(); k++)
      if (x[k]) {
        if (pivots[k].empty()) {
          digest = digest * 1009u + unsigned(x[k]);
          int sc = inv(x[k]);
          for (size_t j = k; j < x.size(); j++)
            x[j] = x[j] * sc % P;
          pivots[k] = std::move(x);
          rank++;
          return true;
        }
        int sc = P - x[k];
        for (size_t j = k; j < x.size(); j++)
          x[j] = (x[j] + pivots[k][j] * sc) % P;
      }
    return false;
  }
};
// A single coordinate of the boundary in wedge^6 R1 tensor R2.
Probe select_probe(const Family &f, Seeds &s) {
  for (;;) {
    auto a = decode(s.seeds[s.rng() % s.seeds.size()]);
    int i = s.rng() % 8, j = s.rng() % 8;
    if (j == i)
      continue;
    if (j < i)
      swap(i, j);
    const auto &choices = f.brackets[a[i]][a[j]];
    if (choices.empty())
      continue;
    int v = choices[s.rng() % choices.size()].first;
    vector<int> rest;
    for (int t = 0; t < 8; t++)
      if (t != i && t != j)
        rest.push_back(a[t]);
    Probe probe;
    for (auto [k, l, c] : f.preimages[v]) {
      if (binary_search(rest.begin(), rest.end(), k) ||
          binary_search(rest.begin(), rest.end(), l))
        continue;
      Slots all;
      copy(rest.begin(), rest.end(), all.begin());
      all[6] = k;
      all[7] = l;
      sort(all.begin(), all.end());
      int parity = (lower_bound(all.begin(), all.end(), k) - all.begin()) +
                   (lower_bound(all.begin(), all.end(), l) - all.begin());
      probe.entries.push_back({encode(all), parity % 2 ? -c : c});
    }
    return probe;
  }
}
using Probes = array<vector<Probe>, 3>;
// Pull a coordinate of wedge^7 R1 back through contraction and through
// each lowering word in the dual quadratic/high bracket.
Probes high_fusion(const Family &f, const array<Weight, 3> &ws,
                   const string &filename) {
  struct FusionTerm {
    int input, sector, coeff;
    vector<int> lowering;
  };
  array<vector<FusionTerm>, 3> terms;
  ifstream in(filename);
  if (!in)
    throw runtime_error("missing fusion coefficients");
  int m, q, h, c, n;
  while (in >> m >> q >> h >> c >> n) {
    FusionTerm x{q, h, mod(c), {}};
    for (int k = 0; k < n; k++) {
      int l;
      in >> l;
      x.lowering.push_back(l);
    }
    terms.at(m).push_back(x);
  }
  Probes probes;
  for (int m = 0; m < 3; m++) {
    Weight w = ws[m];
    w[0]++;
    Seeds seeds(f, w, 7);
    for (int row = 0; row < min(64, int(seeds.seeds.size())); row++) {
      auto a = decode(seeds.seeds[row]);
      array<Vec, 3> coeffs;
      for (const auto &x : terms[m]) {
        int pos = lower_bound(a.begin(), a.begin() + 7, x.input) - a.begin();
        if (a[pos] == x.input)
          continue;
        auto b = a;
        b[7] = x.input;
        sort(b.begin(), b.end());
        Vec v{{encode(b), 1}};
        for (auto it = x.lowering.rbegin(); it != x.lowering.rend(); ++it)
          v = f.apply(v, *it + 1, *it, true);
        for (auto [k, c] : v)
          coeffs[x.sector][k] += ((pos % 2) ? -1 : 1) * c * x.coeff;
      }
      for (h = 0; h < 3; h++) {
        Probe probe;
        for (auto [k, c] : coeffs[h])
          if (mod(c))
            probe.entries.push_back({k, mod(c)});
        probes[h].push_back(std::move(probe));
      }
    }
  }
  return probes;
}
#ifndef EXTREME_LIBRARY
int main(int ac, char **av) {
  static_assert(sizeof(int) >= 4);
  Family f;
  array<Weight, 3> ws{
      {{5, -3, -3, -3, -4}, {4, -2, -3, -3, -4}, {4, -3, -3, -3, -3}}};
  array<int, 3> bounds{14, 145, 16};
  if (ac != 2)
    throw runtime_error("expected dual-bracket coefficient file");
  auto high = high_fusion(f, ws, av[1]);
  const int width_qq = 512;
  Rank combined(3 * width_qq + high[0].size());
  for (int idx = 0; idx < 3; idx++) {
    auto w = ws[idx];
    Seeds s(f, w);
    int target = bounds[idx], width = target + 100;
    vector<Probe> probes;
    for (int k = 0; k < width_qq; k++)
      probes.push_back(select_probe(f, s));
    Rank chains(width), boundaries(probes.size());
    for (int k = 0; k < target + 100; k++) {
      Vec test;
      for (int j = 0; j < 4; j++)
        test[s.seeds[s.rng() % s.seeds.size()]] = 1 + s.rng() % (P - 1);
      auto v = project(f, w, test);
      for (int i = 0; i < 4; i++)
        if (!f.apply(v, i, i + 1).empty())
          throw logic_error("nonhighest projection");
      vector<int> x(width);
      for (int i = 0; i < width; i++)
        x[i] = v[s.seeds[i]];
      chains.insert(x);
      x.clear();
      for (auto &probe : probes)
        x.push_back(probe.eval(v));
      vector<int> all(3 * width_qq);
      copy(x.begin(), x.end(), all.begin() + idx * width_qq);
      for (auto &probe : high[idx])
        all.push_back(probe.eval(v));
      boundaries.insert(x);
      combined.insert(all);
      if (k % 30 == 0 || chains.rank == target)
        cerr << idx << " samples " << k + 1 << " chains " << chains.rank
             << " brackets " << boundaries.rank << " both " << combined.rank
             << " hash " << combined.digest << endl;
      if (chains.rank == target)
        break;
    }
    if (chains.rank != target)
      throw runtime_error("failed to span");
  }
  if (combined.rank != 174)
    throw runtime_error("projected rank certificate failed");
  cout << "combined rank " << combined.rank << endl;
}

#endif
