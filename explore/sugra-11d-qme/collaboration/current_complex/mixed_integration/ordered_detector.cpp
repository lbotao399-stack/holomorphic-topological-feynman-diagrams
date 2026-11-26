// Exact finite-field integration of one ordered theta graph with three paths.
#define TWIST_NATIVE_PARTITION_LIBRARY
#include "native_partition.cpp"

namespace ordered {
using Key = uint64_t;
// The fibers consist of small sparse polynomials. Contiguous entries avoid
// allocating a node per monomial. Large fibers use a linear-probing index.
Key hash_key(Key value) {
  value ^= value >> 30;
  value *= 0xbf58476d1ce4e5b9ULL;
  value ^= value >> 27;
  value *= 0x94d049bb133111ebULL;
  return value ^ (value >> 31);
}

struct Poly {
  vector<pair<Key, int>> terms;
  vector<unsigned> buckets;
  bool empty() const { return terms.empty(); }
  auto begin() const { return terms.begin(); }
  auto end() const { return terms.end(); }

  unsigned &bucket(Key key) {
    size_t mask = buckets.size() - 1;
    size_t index = hash_key(key) & mask;
    while (buckets[index] && terms[buckets[index] - 1].first != key)
      index = (index + 1) & mask;
    return buckets[index];
  }
  void rebuild(size_t capacity) {
    buckets.assign(capacity, 0);
    for (unsigned i = 0; i < terms.size(); ++i)
      bucket(terms[i].first) = i + 1;
  }
  int &operator[](Key key) {
    if (buckets.empty() && terms.size() < 12) {
      for (auto &term : terms)
        if (term.first == key)
          return term.second;
      terms.emplace_back(key, 0);
      return terms.back().second;
    }
    if (terms.size() * 3 >= buckets.size() * 2)
      rebuild(max<size_t>(32, 2 * buckets.size()));
    unsigned &slot = bucket(key);
    if (!slot) {
      terms.emplace_back(key, 0);
      slot = terms.size();
    }
    return terms[slot - 1].second;
  }
};
using Fiber = array<Poly, 25>;
using Triple = array<Poly, 125>;
using native::prime;
using native::red;
constexpr Key mask = (Key(1) << 20) - 1;
Key pack(int b, int c) { return Key(b) + (Key(c) << 20); }
void add(Poly &dst, Key key, int64_t v) {
  int &c = dst[key];
  c = red(c + v);
}
void times(Poly &dst, const Poly &x, const Poly &y, int scale = 1) {
  if (!scale)
    return;
  for (auto [i, a] : x)
    if (a)
      for (auto [j, b] : y)
        if (b)
          add(dst, i + j, int64_t(scale) * a * b);
}
void trim(Poly &p) {
  auto end = remove_if(p.terms.begin(), p.terms.end(),
                       [](const auto &term) { return !term.second; });
  if (end == p.terms.end())
    return;
  p.terms.erase(end, p.terms.end());
  if (!p.buckets.empty()) {
    if (p.terms.size() < 12)
      p.buckets.clear();
    else
      p.rebuild(p.buckets.size());
  }
}
struct Split3 {
  int keys[3], bown, cown, coefficient;
};
const vector<Split3> &parts(int key) {
  static map<int, vector<Split3>> cache;
  auto [it, fresh] = cache.try_emplace(key);
  if (fresh)
    for (auto a : splits(key))
      for (auto b : splits(key - a.key))
        it->second.push_back({{a.key, b.key, key - a.key - b.key},
                              b.deg,
                              degree(key) - a.deg - b.deg,
                              a.binomial * b.binomial});
  return it->second;
}
struct Job {
  int id, det, marked = -1;
  array<int, 3> shape;
  array<native::Waves, 3> waves;
  vector<array<int, 4>> probes;
  array<array<Triple, 2>, 4> endpoints{};
  array<int, 4> signs{};
  vector<int> angular;
  bool read() {
    if (!(cin >> id >> det))
      return false;
    for (int &n : shape) {
      cin >> n;
      if (n < 1 || n > 6)
        throw invalid_argument("path length");
    }
    for (int h = 0; h < 3; ++h) {
      auto &w = waves[h];
      w.active = shape[h] - 1;
      w.labels.fill(-1);
      for (int i = 0; i < w.active; ++i) {
        cin >> w.labels[i];
        if (w.labels[i] == 9)
          marked = h;
      }
      for (int i = 0; i < w.active; ++i)
        for (int &v : w.p[i])
          cin >> v;
      int slots;
      cin >> slots;
      w.polar.resize(slots);
      for (auto &f : w.polar)
        for (int i = 0; i < w.active; ++i)
          for (int &x : f[i])
            cin >> x;
      w.make_moments(1);
      int scale = red(fact(w.active) * native::inverse_prime(fact(5)));
      for (double &x : w.moment)
        x = red(int64_t(x) * scale);
    }
    if (marked < 0)
      throw invalid_argument("missing one-form");
    int n;
    cin >> n;
    probes.resize(n);
    for (auto &probe : probes)
      for (int &v : probe)
        cin >> v;
    for (int t = 0; t < 4; ++t) {
      cin >> signs[t];
      for (Triple &end : endpoints[t]) {
        int count;
        cin >> count;
        for (int j = 0; j < count; ++j) {
          int cell, b, c, v;
          cin >> cell >> b >> c >> v;
          add(end.at(cell), pack(b, c), v);
        }
      }
    }
    if (!cin)
      throw invalid_argument("partial geometry");
    return true;
  }
  Fiber fiber(int h, int index, int flow, int order = 0) const {
    auto &w = waves[h];
    int kind =
        h == marked
            ? int(flow >= 2)
            : (flow < 2 ? flow : ((h == (marked + 1) % 3) == (flow == 2)));
    Fiber matrix;
    for (auto &[k, raw] : w.chain(order, index, kind)) {
      int sign = degree(k) & 1 ? -1 : 1;
      for (int t = 0; t < 25; ++t)
        if (raw[t]) {
          if (h == 1)
            add(matrix[t], pack(k, 0), raw[t]);
          else if (h == 2)
            add(matrix[t], pack(0, k), raw[t]);
          else
            for (auto part : splits(k))
              add(matrix[t], pack(part.key, k - part.key),
                  sign * part.binomial * raw[t]);
        }
    }
    return matrix;
  }
  static Triple apply(const Triple &tensor, const Fiber &matrix, int h) {
    Triple next;
    const int step = h == 0 ? 25 : 5;
    for (int i = 0; i < 125; ++i)
      if (!tensor[i].empty()) {
        int k = (i / step) % 5;
        for (int j = 0; j < 5; ++j)
          times(next[i + (j - k) * step], tensor[i], matrix[5 * k + j]);
      }
    for (auto &x : next)
      trim(x);
    return next;
  }
  Poly polynomial() const {
    Poly result;
    for (int t = 0; t < 4; ++t) {
      array<vector<Fiber>, 3> mats;
      for (int h = 0; h < 3; ++h)
        for (int i = 0; i < int(waves[h].polar.size()); ++i)
          mats[h].push_back(fiber(h, i, t));
      map<pair<int, int>, Fiber> combined;
      for (auto [a, b, c, weight] : probes) {
        auto &target = combined[{a, b}];
        for (int j = 0; j < 25; ++j)
          for (auto [k, v] : mats[2][c][j])
            add(target[j], k, int64_t(weight) * v);
      }
      vector<Triple> first;
      for (auto &matrix : mats[0])
        first.push_back(apply(endpoints[t][0], matrix, 0));
      for (auto &[ab, third] : combined) {
        auto left = apply(first[ab.first], mats[1][ab.second], 1);
        Fiber sides;
        for (int s = 0; s < 25; ++s)
          for (int i = 0; i < 5; ++i)
            for (int j = 0; j < 5; ++j)
              times(sides[5 * i + j], left[5 * s + i],
                    endpoints[t][1][5 * s + j]);
        for (int j = 0; j < 25; ++j) {
          trim(sides[j]);
          trim(third[j]);
          times(result, sides[j], third[j], signs[t]);
        }
      }
    }
    trim(result);
    return result;
  }
  static int slot(int b, int c, int x, int y) {
    return ((b * 9 + c) * 9 + x) * 9 + y;
  }
  void make_angular() {
    int moment[9][9] = {};
    for (int i = 0; i <= 8; ++i)
      for (int j = 0; i + j <= 8; ++j) {
        int m = 1, d = 1;
        for (int k = 0; k < i; ++k)
          m = red(m * (13 - 2 * shape[1] + 2 * k));
        for (int k = 0; k < j; ++k)
          m = red(m * (13 - 2 * shape[2] + 2 * k));
        for (int k = 0; k < i + j; ++k)
          d = red(d * (13 + 2 * k));
        moment[i][j] = red(int64_t(m) * native::inverse_prime(d));
      }
    angular.resize(9 * 9 * 9 * 9);
    for (int b = 0; b <= 8; ++b)
      for (int c = 0; b + c <= 8; ++c)
        for (int x = 0; x <= b; ++x)
          for (int y = 0; y <= c; ++y) {
            int v = 0;
            for (int i = 0; i <= x; ++i)
              for (int j = 0; j <= y; ++j)
                v += ((x + y - i - j) % 2 ? -1 : 1) * choose(x, i) *
                     choose(y, j) * moment[b - x + i][c - y + j];
            angular[slot(b, c, x, y)] = red(v);
          }
  }
  int solve() {
    make_angular();
    auto f = polynomial();
    int result = 0;
    for (auto [k, v] : f) {
      int b = k & mask, c = k >> 20, db = degree(b), dc = degree(c);
      if (db + dc > 8)
        throw logic_error("total degree bound");
      int64_t integral = 0;
      for (const auto &s : parts(b))
        for (const auto &t : parts(c)) {
          int m = angular[slot(db, dc, s.bown, t.cown)];
          for (int h = 0; h < 3 && m; ++h)
            m = red(int64_t(m) *
                    int(waves[h].moment[indices.id[s.keys[h] + t.keys[h]]]));
          integral += int64_t(m) * s.coefficient * t.coefficient;
        }
      result = red(result + int64_t(v) * red(integral));
    }
    return red(int64_t(result) * det);
  }
};
} // namespace ordered
#ifndef TWIST_ORDERED_DETECTOR_LIBRARY
int main() {
  ios::sync_with_stdio(false);
  cin.tie(nullptr);
  while (true) {
    ordered::Job j;
    if (!j.read())
      break;
    cout << j.id << ' ' << j.solve() << '\n' << flush;
  }
  return 0;
}

#endif
