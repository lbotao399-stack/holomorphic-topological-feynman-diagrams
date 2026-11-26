// Filtered vertical CE differential on four vector inputs and one marked
// function. Ordinary plane waves are projected to Taylor order >=2 at each
// vector input.
#include <algorithm>
#include <array>
#include <gmpxx.h>
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>
#include <vector>

using namespace std;
using Vec = array<int, 3>;
int prime = 251;
bool exact = false;
int mod(int64_t x) { return exact ? int(x) : (x % prime + prime) % prime; }
int inv(int x) {
  for (int t = 1; t < prime; ++t)
    if (x * t % prime == 1)
      return t;
  throw runtime_error("nonunit pivot");
}
Vec sum(Vec x, Vec y) {
  for (int k = 0; k < 3; ++k)
    x[k] = mod(x[k] + y[k]);
  return x;
}
Vec scale(Vec x, int c) {
  for (int &v : x)
    v = mod(v * c);
  return x;
}
int dot(Vec a, Vec b) {
  return mod(inner_product(a.begin(), a.end(), b.begin(), 0));
}
Vec cross(Vec a, Vec b) {
  Vec c{};
  for (int i = 0; i < 3; ++i) {
    int j = (i + 1) % 3, k = (i + 2) % 3;
    c[i] = mod(a[j] * b[k] - a[k] * b[j]);
  }
  return c;
}

struct Dual {
  int v = 0, d = 0;
};
mt19937 rng{37791};
Vec random_vec() {
  return {int(rng() % 5) - 2, int(rng() % 5) - 2, int(rng() % 5) - 2};
}

template <int Q> struct Waves {
  // Covectors are the vector momenta, the function momentum, and the test
  // one-form.
  array<Vec, Q + 2> cov{}, dcov{};
  array<Vec, Q> pol{}, dpol{};
};
template <int Q> struct Evaluations {
  static constexpr int n = Q + 2;
  array<Dual, n * n * n> det{};
  array<Dual, n * Q> pair{};
  explicit Evaluations(const Waves<Q> &w) {
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j)
        for (int k = 0; k < n; ++k) {
          if (i == j || i == k || j == k)
            continue;
          auto cp = cross(w.cov[j], w.cov[k]);
          det[n * n * i + n * j + k] = {
              dot(w.cov[i], cp),
              mod(dot(w.dcov[i], cp) +
                  dot(w.cov[i], cross(w.dcov[j], w.cov[k])) +
                  dot(w.cov[i], cross(w.cov[j], w.dcov[k])))};
        }
      for (int j = 0; j < Q; ++j)
        pair[n * j + i] = {
            dot(w.cov[i], w.pol[j]),
            mod(dot(w.dcov[i], w.pol[j]) + dot(w.cov[i], w.dpol[j]))};
    }
  }
};
template <int Q> struct Product {
  array<int, 3> det;
  array<int, Q> pair;
  int sign;
};
template <int Q> using Column = vector<Product<Q>>;

template <int Q>
Column<Q> alternate(array<array<int, 3>, 3> det, array<int, Q> pair) {
  constexpr int n = Q + 2;
  Column<Q> terms;
  array<int, n> ord;
  iota(ord.begin(), ord.end(), 0);
  do {
    int sign = 1;
    for (int i = 0; i < Q; ++i)
      for (int j = 0; j < i; ++j)
        if (ord[j] > ord[i])
          sign = -sign;
    Product<Q> p{};
    p.sign = sign;
    for (int i = 0; i < 3; ++i)
      p.det[i] = n * n * ord[det[i][0]] + n * ord[det[i][1]] + ord[det[i][2]];
    for (int i = 0; i < Q; ++i)
      p.pair[i] = n * ord[i] + ord[pair[i]];
    terms.push_back(p);
  } while (next_permutation(ord.begin(), ord.begin() + Q));
  return terms;
}
template <int Q> int value(const Evaluations<Q> &w, const Column<Q> &c) {
  int total = 0;
  for (const auto &p : c) {
    int v = p.sign;
    for (auto i : p.det)
      v = v * w.det[i].v % prime;
    for (auto i : p.pair)
      v = v * w.pair[i].v % prime;
    total += v;
  }
  return mod(total);
}
struct Echelon {
  vector<vector<int>> piv;
  int rank = 0;
  explicit Echelon(int n) : piv(n) {}
  bool add(vector<int> row) {
    for (size_t i = 0; i < row.size(); ++i)
      if (int v = row[i]) {
        if (piv[i].empty()) {
          int c = inv(v);
          for (size_t j = i; j < row.size(); ++j)
            row[j] = row[j] * c % prime;
          piv[i] = std::move(row);
          ++rank;
          return true;
        }
        for (size_t j = i; j < row.size(); ++j)
          row[j] = mod(row[j] - v * piv[i][j]);
      }
    return false;
  }
};
template <int Q> Waves<Q> random_waves() {
  Waves<Q> w;
  for (auto &p : w.cov)
    p = random_vec();
  for (int i = 0; i < Q; ++i)
    w.pol[i] = cross(w.cov[i], random_vec());
  return w;
}

template <int Q>
vector<Column<Q>> basis(const vector<Waves<Q>> &waves,
                        const vector<int>& dimensions, bool scalar = false,
                        int first_degree = 1) {
  constexpr int n = Q + 2;
  vector<Evaluations<Q>> probes;
  for (const auto &w : waves)
    probes.emplace_back(w);
  vector<Echelon> ranks(dimensions.size(), Echelon(int(waves.size())));
  const int available = scalar ? n - 1 : n;
  int dimension = accumulate(dimensions.begin(), dimensions.end(), 0);
  vector<array<int, 3>> triples;
  for (int i = 0; i < n; ++i)
    for (int j = i + 1; j < n; ++j)
      for (int k = j + 1; k < available; ++k)
        triples.push_back({i, j, k});
  vector<Column<Q>> cols;
  for (int attempt = 0; attempt < 5000000 && int(cols.size()) < dimension;
       ++attempt) {
    array<int, n> degrees{};
    array<array<int, 3>, 3> det;
    array<int, Q> inner;
    for (auto &t : det) {
      t = triples[rng() % triples.size()];
      for (int i : t)
        ++degrees[i];
    }
    if (degrees[Q + 1] > int(!scalar))
      continue;
    for (int i = 0; i < Q; ++i) {
      inner[i] = (i + 1 + rng() % (available - 1)) % available;
      ++degrees[inner[i]];
    }
    int s = degrees[Q];
    int grade = s - first_degree;
    if (degrees[Q + 1] != int(!scalar) || grade < 0 || grade >= int(dimensions.size()) ||
        ranks[grade].rank == dimensions[grade] ||
        *min_element(degrees.begin(), degrees.begin() + Q) < 2)
      continue;
    Column<Q> c = alternate<Q>(det, inner);
    vector<int> row;
    for (const auto &w : probes)
      row.push_back(value(w, c));
    if (ranks[grade].add(std::move(row)))
      cols.push_back(std::move(c));
  }
  if (int(cols.size()) != dimension)
    throw runtime_error("basis incomplete");
  return cols;
}

// Original momenta have components <=2 and transverse polarizations <=8.
// After fusion determinants are bounded by 384. At most one ordinary
// pairing can use the fused polarization (bound 9216); the other three
// pairings are bounded by 96. The product rule, 24 permutations and all
// bracket terms keep every signed integer row entry below 2^90.
// Each differentiated determinant or pairing is taken by the product rule.
using Wide = __int128_t;
template <int Q>
Wide integer_value(const Evaluations<Q> &w, const Column<Q> &c, bool low) {
  Wide total = 0;
  for (const auto &p : c) {
    Wide v = 1, d = 0;
    auto mult = [&](Dual a) {
      d = d * a.v + v * a.d;
      v *= a.v;
    };
    for (int i : p.det)
      mult(w.det[i]);
    for (int i : p.pair)
      mult(w.pair[i]);
    total += p.sign * (low ? v + d : v);
  }
  return total;
}
template <int Q>
vector<Wide> differential(const vector<Column<Q>> &cols,
                          const Waves<Q + 1> &input) {
  const auto &p = input.cov;
  const auto &u = input.pol;
  Vec f = p[Q + 1], eta = p[Q + 2];
  vector<Wide> row(cols.size());
  auto add = [&](const Waves<Q> &w, int sign, bool low) {
    Evaluations<Q> e(w);
    for (size_t i = 0; i < cols.size(); ++i) {
      row[i] += sign * integer_value(e, cols[i], low);
    }
  };
  for (int i = 0; i < Q + 1; ++i) {
    Waves<Q> w;
    w.cov[Q] = f;
    w.cov[Q + 1] = eta;
    int n = 0;
    for (int k = 0; k < Q + 1; ++k)
      if (k != i) {
        w.cov[n] = p[k];
        w.pol[n++] = u[k];
      }
    // Y_i(f) at momentum f + t p_i, with its constant and linear Y_i jet
    // removed.
    int sgn = (i % 2 ? 1 : -1) * dot(u[i], f);
    w.cov[Q] = sum(f, p[i]);
    add(w, sgn, false);
    w.cov[Q] = f;
    w.dcov[Q] = p[i];
    add(w, -sgn, true);

    for (int j = i + 1; j < Q + 1; ++j) {
      Waves<Q> b;
      b.cov[Q] = f;
      b.cov[Q + 1] = eta;
      n = 1;
      for (int k = 0; k < Q + 1; ++k)
        if (k != i && k != j) {
          b.cov[n] = p[k];
          b.pol[n++] = u[k];
        }
      Vec a = scale(u[j], dot(p[j], u[i])), z = scale(u[i], -dot(p[i], u[j]));
      sgn = (i + j) % 2 ? -1 : 1;
      b.cov[0] = sum(p[i], p[j]);
      b.pol[0] = sum(a, z);
      add(b, sgn, false);
      b.cov[0] = p[j];
      b.dcov[0] = p[i];
      b.pol[0] = a;
      b.dpol[0] = z;
      add(b, -sgn, true);
      b.cov[0] = p[i];
      b.dcov[0] = p[j];
      b.pol[0] = z;
      b.dpol[0] = a;
      add(b, -sgn, true);
    }
  }
  return row;
}
mpz_class big(Wide x) {
  __uint128_t u = x < 0 ? -x : x;
  uint64_t limbs[] = {uint64_t(u), uint64_t(u >> 64)};
  mpz_class z;
  mpz_import(z.get_mpz_t(), 2, -1, sizeof(uint64_t), 0, 0, limbs);
  return x < 0 ? -z : z;
}
int bareiss(vector<vector<mpz_class>> &m) {
  int r = 0;
  mpz_class prev = 1;
  for (size_t j = 0; j < m.front().size(); ++j) {
    size_t k = r;
    while (k < m.size() && m[k][j] == 0)
      ++k;
    if (k == m.size())
      continue;
    swap(m[k], m[r]);
    auto pivot = m[r][j];
    for (size_t i = r + 1; i < m.size(); ++i) {
      for (size_t c = j + 1; c < m[i].size(); ++c) {
        m[i][c] = pivot * m[i][c] - m[i][j] * m[r][c];
        if (!mpz_divisible_p(m[i][c].get_mpz_t(), prev.get_mpz_t()))
          throw runtime_error("nonintegral Bareiss division");
        mpz_divexact(m[i][c].get_mpz_t(), m[i][c].get_mpz_t(),
                     prev.get_mpz_t());
      }
      m[i][j] = 0;
    }
    prev = pivot;
    ++r;
  }
  return r;
}
#ifndef MIXED_FUNCTION_COMPONENTS
int main(int argc, char **argv) {
  if (argc > 1)
    prime = stoi(argv[1]);
  if (prime != 251 && prime != 257)
    throw runtime_error("choose 251 or 257");
  vector<Waves<4>> four;
  vector<Waves<5>> five;
  // Generate the integer polarizations before any finite-field reduction.
  exact = true;
  while (four.size() < 110)
    four.push_back(random_waves<4>());
  while (five.size() < 200)
    five.push_back(random_waves<5>());
  exact = false;
  auto cols = basis(four, {57, 49, 26});
  auto targets = basis(five, {84, 45, 9});
  Echelon samples(targets.size());
  for (auto w : five) {
    Evaluations<5> e(w);
    vector<int> row;
    for (const auto &c : targets)
      row.push_back(value(e, c));
    samples.add(std::move(row));
  }
  if (samples.rank != 138)
    throw runtime_error("evaluation on target noninjective");
  exact = true;
  vector<vector<mpz_class>> m;
  for (auto w : five) {
    vector<mpz_class> row;
    for (auto x : differential(cols, w))
      row.push_back(big(x));
    m.push_back(std::move(row));
  }
  int rank = bareiss(m);
  cout << "basis " << cols.size() << ", injective evaluation rank modulo "
       << prime << ": " << samples.rank << ", rational filtered CE rank "
       << rank << "\n";
  if (rank != 88)
    throw runtime_error("unexpected filtered rank");
}
#endif
