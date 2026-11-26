// Exact rank of the second relative CE differential on vector currents and
// scalar outputs. All evaluations and pivots are in F_251. See README.md for
// the rank over Q.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;
constexpr int P = 251, projection_width = 2600;
int mod(int x) { return (x % P + P) % P; }
int power(int a, int n) {
  int b = 1;
  for (; n; n /= 2, a = a * a % P)
    if (n % 2)
      b = b * a % P;
  return b;
}
using Coordinates = array<int, 5>;
using SixSlots = array<int, 6>;
using SevenSlots = array<int, 7>;
int dot(const Coordinates &a, const Coordinates &b) {
  return mod(inner_product(a.begin(), a.end(), b.begin(), 0));
}
int determinant(array<Coordinates, 5> m) {
  int ans = 1;
  for (int i = 0; i < 5; i++) {
    int j = i;
    while (j < 5 && m[j][i] == 0)
      j++;
    if (j == 5)
      return 0;
    if (i != j) {
      swap(m[i], m[j]);
      ans = P - ans;
    }
    ans = ans * m[i][i] % P;
    int inv = power(m[i][i], P - 2);
    for (j = i + 1; j < 5; j++) {
      int f = m[j][i] * inv % P;
      for (int h = i + 1; h < 5; h++)
        m[j][h] = mod(m[j][h] - f * m[i][h]);
    }
  }
  return ans;
}
int encode(const SixSlots &f) {
  int a = 0;
  for (int i = 0; i < 6; i++)
    a = a * 5 + f[i] - (f[i] > i);
  return a;
}
struct Pairing {
  SixSlots f, indegree{};
};
struct Determinants {
  int a, b, c, d;
  SixSlots omissions{};
};
vector<Pairing> contraction_maps;
vector<Determinants> determinant_products;
map<SixSlots, vector<int>> omissions;
// Each orbit sum has one signed projection column. projection_width+1 denotes
// an orbit with an odd stabilizer.
vector<int> projected_columns;
int determinant_index[6][6][6][7];
bool current_output;
// D_a D_b E_cd, with a<=b and c<d. Each D omits one momentum, E omits two and
// inserts gamma.
int signed_orbit_columns() {
  auto append_product = [](int a, int b, int c, int d) {
    Determinants h{a, b, c, d};
    for (int i : {a, b, c, d})
      if (i >= 0)
        h.omissions[i]++;
    determinant_index[a][b][c][d + 1] = determinant_products.size();
    omissions[h.omissions].push_back(determinant_products.size());
    determinant_products.push_back(h);
  };
  for (int a = 0; a < 6; a++)
    for (int b = a; b < 6; b++)
      for (int c = 0; c < 6; c++) {
        if (current_output) {
          for (int d = c + 1; d < 6; d++)
            append_product(a, b, c, d);
        } else if (c >= b)
          append_product(a, b, c, -1);
      }
  for (int index = 0; index < 15625; index++) {
    int a = index;
    Pairing p;
    for (int i = 5; i >= 0; i--) {
      int v = a % 5;
      a /= 5;
      p.f[i] = v + (v >= i);
      p.indegree[p.f[i]]++;
    }
    contraction_maps.push_back(p);
  }
  projected_columns.assign(
      contraction_maps.size() * determinant_products.size(), 0);
  SixSlots identity{0, 1, 2, 3, 4, 5};
  vector<SixSlots> perms;
  SixSlots perm = identity;
  do {
    perms.push_back(perm);
  } while (next_permutation(perm.begin(), perm.end()));
  mt19937 rng(10091);
  int count = 0;
  for (int fi = 0; fi < (int)contraction_maps.size(); fi++)
    for (int di = 0; di < (int)determinant_products.size(); di++) {
      int key = fi * determinant_products.size() + di;
      if (projected_columns[key])
        continue;
      const auto &h = determinant_products[di];
      const auto &m = contraction_maps[fi];
      bool positive = true;
      for (int i = 0; i < 6; i++)
        if (3 + m.indegree[i] - h.omissions[i] < 2)
          positive = false;
      if (!positive)
        continue;
      vector<pair<int, int>> keys;
      bool cancels = false;
      for (const auto &p : perms) {
        SixSlots f;
        for (int i = 0; i < 6; i++)
          f[p[i]] = p[m.f[i]];
        int a = p[h.a], b = p[h.b], c = p[h.c],
            d = current_output ? p[h.d] : -1;
        if (a > b)
          swap(a, b);
        int sign = 1;
        if (!current_output) {
          if (b > c)
            swap(b, c);
          if (a > b)
            swap(a, b);
        } else if (c > d) {
          swap(c, d);
          sign = -1;
        }
        int k = encode(f) * determinant_products.size() +
                determinant_index[a][b][c][d + 1];
        keys.push_back({k, sign});
        if (k == key && sign == -1)
          cancels = true;
      }
      int col = 1 + rng() % projection_width;
      if (rng() % 2)
        col = -col;
      for (auto [k, s] : keys)
        projected_columns[k] = cancels ? projection_width + 1 : s * col;
      if (!cancels)
        count++;
    }
  return count;
}
// x[k] is the coefficient of s^k t^(degree-k).
struct Poly {
  int degree = 0;
  array<int, 11> x{1};
};
Poly multiply_linear(Poly p, int v, int w, bool linear = true) {
  if (linear) {
    p.x[p.degree + 1] = p.x[p.degree] * w % P;
    for (int j = p.degree; j >= 1; j--)
      p.x[j] = (p.x[j] * v + p.x[j - 1] * w) % P;
    p.degree++;
  } else {
    for (int j = 1; j <= p.degree; j++)
      p.x[j] = p.x[j] * v % P;
  }
  p.x[0] = p.x[0] * v % P;
  return p;
}
struct Fixture {
  array<Coordinates, 7> p, u;
  Coordinates gamma;
  SevenSlots degree;
};
vector<SevenSlots> test_orders;
Fixture sample(mt19937 &rng) {
  Fixture v;
  v.degree = test_orders[rng() % test_orders.size()];
  for (int j = 6; j > 0; j--)
    swap(v.degree[j], v.degree[rng() % (j + 1)]);
  for (int t = 0; t < 7; t++) {
    auto &p = v.p[t];
    auto &u = v.u[t];
    p.fill(0);
    u.fill(0);
    for (int j = 0; j < 3; j++)
      p[rng() % 5] = 1 + rng() % (P - 1);
    for (int j = 0; j < 3; j++)
      u[rng() % 5] = 1 + rng() % (P - 1);
    int k = 0;
    while (!p[k])
      k++;
    u[k] = 0;
    u[k] = mod(-dot(p, u) * power(p[k], P - 2));
  }
  for (auto &x : v.gamma)
    x = rng() % P;
  return v;
}
vector<int> differential_row(const Fixture &v) {
  vector<int> result(projection_width);
  int dots[7][7];
  for (int a = 0; a < 7; a++)
    for (int b = 0; b < 7; b++)
      dots[a][b] = dot(v.u[a], v.p[b]);
  for (int i = 0; i < 7; i++)
    for (int j = i + 1; j < 7; j++) {
      if (!dots[i][j] && !dots[j][i])
        continue;
      SixSlots rest{};
      for (int k = 0, l = 1; k < 7; k++)
        if (k != i && k != j)
          rest[l++] = k;
      int weight = ((i + j) % 2) ? -1 : 1;
      // Signed minors of six momenta; E omits two and has gamma in the last
      // row.
      int single_omissions[2][6], double_omissions[2][6][6];
      for (int side = 0; side < 2; side++) {
        rest[0] = side ? i : j;
        for (int a = 0; a < 6; a++) {
          array<Coordinates, 5> m;
          int s = 0;
          for (int t = 0; t < 6; t++)
            if (t != a)
              m[s++] = v.p[rest[t]];
          single_omissions[side][a] = mod((a % 2 ? -1 : 1) * determinant(m));
          if (!current_output)
            continue;
          for (int b = a + 1; b < 6; b++) {
            s = 0;
            for (int t = 0; t < 6; t++)
              if (t != a && t != b)
                m[s++] = v.p[rest[t]];
            m[4] = v.gamma;
            double_omissions[side][a][b] =
                mod(((a + b - 1) % 2 ? -1 : 1) * determinant(m));
          }
        }
      }
      // Precompute the determinant products as binary forms in the two fused
      // momenta.
      vector<Poly> minor_polynomials(determinant_products.size());
      for (int di = 0; di < (int)determinant_products.size(); di++) {
        const auto &h = determinant_products[di];
        Poly p;
        for (int a : {h.a, h.b})
          p = multiply_linear(p, single_omissions[0][a], single_omissions[1][a],
                              a != 0);
        if (!current_output) {
          minor_polynomials[di] = multiply_linear(
              p, single_omissions[0][h.c], single_omissions[1][h.c], h.c != 0);
          continue;
        }
        minor_polynomials[di] =
            multiply_linear(p, double_omissions[0][h.c][h.d],
                            double_omissions[1][h.c][h.d], h.c != 0);
      }
      for (int fi = 0; fi < (int)contraction_maps.size(); fi++) {
        const auto &m = contraction_maps[fi];
        const auto &f = m.f;
        int sp = 1;
        for (int k = 1; k < 6 && sp; k++)
          if (f[k])
            sp = sp * dots[rest[k]][rest[f[k]]] % P;
        if (!sp)
          continue;
        // Spectator orders determine the required multiset of omitted
        // indices.
        SixSlots o{};
        o[0] = current_output ? 4 : 3;
        bool good = true;
        for (int k = 1; k < 6; k++) {
          o[k] = 3 + m.indegree[k] - v.degree[rest[k]];
          o[0] -= o[k];
          if (o[k] < 0 || o[k] > 3)
            good = false;
        }
        if (!good || o[0] < 0 || o[0] > 3)
          continue;
        auto found = omissions.find(o);
        if (found == omissions.end())
          continue;
        Poly p;
        p.x[0] = sp;
        int a = rest[f[0]];
        p = multiply_linear(p, dots[i][j] * dots[j][a] % P,
                            mod(-dots[j][i] * dots[i][a]));
        for (int k = 1; k < 6; k++)
          if (!f[k])
            p = multiply_linear(p, dots[rest[k]][j], dots[rest[k]][i]);
        for (int di : found->second) {
          int col = projected_columns[fi * determinant_products.size() + di];
          if (!col || col == projection_width + 1)
            continue;
          const auto &dp = minor_polynomials[di];
          int val = 0;
          for (int k = 0; k <= p.degree; k++) {
            int h = v.degree[i] - k;
            if (h >= 0 && h <= dp.degree)
              val += p.x[k] * dp.x[h];
          }
          result[abs(col) - 1] += weight * (col < 0 ? -1 : 1) * (val % P);
        }
      }
    }
  for (auto &v : result)
    v = mod(v);
  return result;
}
// Direct torus interpolation of the fused momentum power, using unexpanded
// determinants. The other six-argument momentum degrees are checked before
// evaluation.
int direct_entry(const Fixture &v, int column) {
  vector<pair<int, int>> terms;
  for (int f = 0; f < (int)contraction_maps.size(); f++)
    for (int h = 0; h < (int)determinant_products.size(); h++)
      if (abs(projected_columns[f * determinant_products.size() + h]) ==
          column + 1)
        terms.push_back({f, h});
  int root = 2;
  while (power(root, 25) != 1 || power(root, 5) == 1)
    ++root;
  int answer = 0;
  for (int i = 0; i < 7; i++)
    for (int j = i + 1; j < 7; j++) {
      SixSlots order{};
      for (int s = 0, slot = 1; s < 7; s++)
        if (s != i && s != j)
          order[slot++] = s;
      vector<pair<int, int>> allowed;
      for (auto [fi, di] : terms) {
        bool ok = true;
        for (int k = 1; k < 6; k++)
          if (3 + contraction_maps[fi].indegree[k] -
                  determinant_products[di].omissions[k] !=
              v.degree[order[k]])
            ok = false;
        if (ok)
          allowed.push_back({fi, di});
      }
      if (allowed.empty())
        continue;
      for (int l = 0, scale = 1; l < 25; l++, scale = scale * root % P) {
        array<Coordinates, 6> p, u;
        for (int a = 0; a < 5; a++) {
          p[0][a] = mod(scale * v.p[i][a] + v.p[j][a]);
          u[0][a] = mod(dot(v.u[i], v.p[j]) * v.u[j][a] -
                        scale * dot(v.u[j], v.p[i]) * v.u[i][a]);
        }
        for (int s = 1; s < 6; s++) {
          p[s] = v.p[order[s]];
          u[s] = v.u[order[s]];
        }
        auto minor = [&](int a, int b) {
          array<Coordinates, 5> m;
          int s = 0;
          for (int k = 0; k < 6; k++)
            if (k != a && k != b)
              m[s++] = p[k];
          int power_of_minus_one = a;
          if (b >= 0) {
            m[4] = v.gamma;
            power_of_minus_one += b - 1;
          }
          return mod((power_of_minus_one % 2 ? -1 : 1) * determinant(m));
        };
        int value = 0;
        for (auto [fi, di] : allowed) {
          auto h = determinant_products[di];
          int term = minor(h.a, -1) * minor(h.b, -1) % P * minor(h.c, h.d) % P;
          for (int k = 0; k < 6; k++)
            term = term * dot(u[k], p[contraction_maps[fi].f[k]]) % P;
          value =
              mod(value +
                  (projected_columns[fi * determinant_products.size() + di] < 0
                       ? -term
                       : term));
        }
        int projection =
            power(scale, mod(P - 1 - v.degree[i])) * power(25, P - 2) % P;
        answer = mod(answer + ((i + j) % 2 ? -1 : 1) * value * projection);
      }
    }
  return answer;
}
int main(int argc, char **argv) {
  if (argc != 2 ||
      (std::string(argv[1]) != "current" && std::string(argv[1]) != "scalar"))
    throw runtime_error("expected current or scalar");
  current_output = std::string(argv[1]) == "current";
  const int samples = 8000, target = current_output ? 1801 : 555;
  cerr << "columns " << signed_orbit_columns() << endl;
  SevenSlots d;
  auto gen = [&](auto &&g, int i, int left, int lo) -> void {
    if (i == 7) {
      if (!left)
        test_orders.push_back(d);
      return;
    }
    for (int k = lo; k * (7 - i) <= left; k++) {
      d[i] = k;
      g(g, i + 1, left - k, k);
    }
  };
  gen(gen, 0, current_output ? 21 : 22, 2);
  mt19937 rng(15123);
  vector<vector<int>> pivots(projection_width);
  int rank = 0, interpolations = 0;
  uint32_t hash = 0;
  for (int s = 1; s <= samples; s++) {
    auto fixture = sample(rng);
    auto v = differential_row(fixture);
    // Compare the coefficient convolution to a 25-point interpolation before
    // row reduction.
    if (s <= 32) {
      auto pos = find_if(v.begin(), v.end(), [](int x) { return x != 0; });
      if (pos != v.end()) {
        if (direct_entry(fixture, int(pos - v.begin())) != *pos)
          throw runtime_error("direct interpolation disagrees");
        interpolations++;
      }
    }
    for (int k = 0; k < projection_width; k++)
      if (v[k]) {
        if (pivots[k].empty()) {
          int inv = power(v[k], P - 2);
          for (int j = k; j < projection_width; j++)
            v[j] = v[j] * inv % P;
          pivots[k] = std::move(v);
          rank++;
          hash = hash * 1000003u + (s * projection_width + k);
          break;
        }
        int f = P - v[k];
        for (int j = k; j < projection_width; j++)
          v[j] = (v[j] + f * pivots[k][j]) % P;
      }
    if (s % 1000 == 0 || s == samples || rank == target)
      cerr << "samples " << s << " rank " << rank << " hash " << hash << endl;
    if (rank == target)
      break;
  }
  cerr << "independent interpolations " << interpolations << endl;
  if (rank != target || interpolations == 0)
    throw runtime_error("rank certificate incomplete");
}
