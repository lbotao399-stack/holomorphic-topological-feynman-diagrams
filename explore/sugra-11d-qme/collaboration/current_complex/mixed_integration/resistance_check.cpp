// Check the summed resistance pole directly in the original four cubic Wick
// words.
#define SUGRA_PHYSICAL_WORDS_ONLY
#include "../../five_theta/check_physical.cpp"
#include <algorithm>
#include <random>
using namespace std;
namespace {
constexpr int ring = 49, max_degree = 11;
int reduce(int x) { return (x % ring + ring) % ring; }
int inner(V a, V b) {
  int x = 0;
  for (int i = 0; i < 5; ++i)
    x += a[i] * b[i];
  return reduce(x);
}
V plus(V a, V b, int c = 1) {
  for (int i = 0; i < 5; ++i)
    a[i] = reduce(a[i] + c * b[i]);
  return a;
}
using K = array<V, 14>;
using Mat = array<int, 25>;
using Path = array<Mat, 6>;
using Endpoint = array<Tensor, 3>;
using Observed = array<Path, max_degree + 1>;

Tensor literal(const Vertex &vertex, const K &momenta, const V &polar) {
  Tensor a{};
  for (auto w : vertex.words) {
    int n = vertex.cubic ? w.scale / inv6 * 41 : w.scale;
    if (w.external_axis >= 0)
      n *= polar[w.external_axis];
    n = reduce(n * momenta[w.first_slot][w.first_axis]);
    if (w.second_slot >= 0)
      n *= momenta[w.second_slot][w.second_axis];
    a[w.element] += n;
  }
  for (auto &n : a)
    n = reduce(n);
  return a;
}
struct Fixture {
  array<V, 10> p{}, u{};
  array<V, 3> ell{}, total{};
  array<array<V, 2>, 3> obs{};
  array<array<int, 5>, 2> lists{{{0, 1, 2, 3, 9}, {4, 5, 6, 7, 8}}};
  explicit Fixture(mt19937 &rng) {
    auto vec = [&]() {
      V v{};
      for (auto &a : v)
        a = rng() % ring;
      return v;
    };
    for (int j = 0; j < 9; ++j) {
      p[j] = vec();
      V v = vec();
      // A random skew matrix applied to p gives an exactly transverse
      // polarization.
      for (int a = 0; a < 5; ++a) {
        int b = (a + 1) % 5;
        u[j][a] += v[a] * p[j][b];
        u[j][b] -= v[a] * p[j][a];
      }
      for (auto &a : u[j])
        a = reduce(a);
      if (inner(p[j], u[j]))
        throw logic_error("polarization");
      p[9] = plus(p[9], p[j], -1);
    }
    u[9] = vec();
    ell[1] = vec();
    ell[2] = vec();
    ell[0] = plus(plus(ell[0], ell[1], -1), ell[2], -1);
    for (int h = 1; h <= 2; ++h) {
      for (int l : lists[h - 1])
        total[h] = plus(total[h], p[l]);
      for (auto &a : obs[h])
        a = vec();
    }
  }
  Endpoint endpoint(const Flow &flow, bool right) const {
    int v = right ? 11 : 0;
    array<Tensor, 3> vals{};
    for (int t = 0; t < 3; ++t) {
      K k{};
      for (int h = 0; h < 3; ++h) {
        int e = right ? edgeids[h].back() : edgeids[h].front();
        k[e] = plus(right ? total[h] : V{}, ell[h], t);
        if (!right)
          k[e] = plus(V{}, k[e], -1);
      }
      vals[t] = literal(*flow.local[v], k, {});
    }
    Endpoint coeff{};
    for (int i = 0; i < 125; ++i) {
      coeff[0][i] = vals[0][i];
      coeff[2][i] = reduce(25 * (vals[2][i] - 2 * vals[1][i] + vals[0][i]));
      coeff[1][i] = reduce(vals[1][i] - vals[0][i] - coeff[2][i]);
    }
    return coeff;
  }
  Observed path(const FlowPlans &plans, int flow_id, int h) const {
    Observed result{};
    auto labels = lists[h - 1];
    do {
      int g =
          h == 1 ? find(labels.begin(), labels.end(), 9) - labels.begin() : 0;
      auto &flow = plans.by_sites[1 << g].at(flow_id);
      int order_sign = (h == 1 && (g & 1)) ? -1 : 1;
      Path matrix{};
      for (int j = 0; j < 5; ++j)
        matrix[0][6 * j] = 1;
      array<array<int, 7>, 7> moments{};
      moments[0][0] = 1;
      V prefix{};
      for (int pos = 0; pos < 5; ++pos) {
        int l = labels[pos], v = vertices[h][pos];
        array<Mat, 2> factor{};
        for (int t = 0; t < 2; ++t) {
          K k{};
          k[13] = p[l];
          auto in = plus(prefix, ell[h], t);
          k[edgeids[h][pos]] = in;
          k[edgeids[h][pos + 1]] = plus(V{}, plus(in, p[l]), -1);
          auto val = literal(*flow.local[v], k, u[l]);
          for (int k = 0; k < 25; ++k)
            factor[t][k] = reduce(val[k] - (t ? factor[0][k] : 0));
        }
        Path next{};
        for (int d = 0; d <= pos; ++d)
          for (int i = 0; i < 5; ++i)
            for (int j = 0; j < 5; ++j)
              for (int k = 0; k < 5; ++k)
                for (int q = 0; q < 2; ++q)
                  next[d + q][5 * i + j] +=
                      matrix[d][5 * i + k] * factor[q][5 * k + j];
        for (auto &a : next)
          for (auto &x : a)
            x = reduce(x);
        matrix = next;
        prefix = plus(prefix, p[l]);
        int sx = inner(prefix, obs[h][0]), sy = inner(prefix, obs[h][1]);
        for (int i = 0; i <= 6; ++i)
          for (int j = 0; j <= 6; ++j) {
            if (i)
              moments[i][j] += sx * moments[i - 1][j];
            if (j)
              moments[i][j] += sy * moments[i][j - 1];
            moments[i][j] = reduce(moments[i][j]);
          }
      }
      for (int d = 0; d <= max_degree; ++d) {
        int n = order_sign * moments[d / 2][d - d / 2];
        for (int a = 0; a <= 5; ++a)
          for (int j = 0; j < 25; ++j)
            result[d][a][j] += n * matrix[a][j];
      }
    } while (next_permutation(labels.begin(), labels.end()));
    for (auto &p : result)
      for (auto &m : p)
        for (int &v : m)
          v = reduce(v);
    return result;
  }
};
int coefficient(const Endpoint &left, const Endpoint &right, const Path &b,
                const Path &c, int d) {
  long long acc = 0;
  for (int a = 0; a < 3; ++a)
    for (int e = 0; e < 3; ++e)
      for (int j = 0; j < 6; ++j) {
        int k = d - a - e - j;
        if (k < 0 || k > 5)
          continue;
        for (int n = 0; n < 125; ++n)
          if (left[a][n]) {
            int i = n / 25, r = n / 5 % 5, s = n % 5;
            for (int u = 0; u < 5; ++u)
              for (int v = 0; v < 5; ++v)
                acc += 1LL * left[a][n] * right[e][25 * i + 5 * u + v] *
                       b[j][5 * r + u] * c[k][5 * s + v];
          }
      }
  return acc % ring;
}
} // namespace
int main() {
  topology({1, 6, 6});
  FlowPlans plans(1);
  mt19937 rng(711637);
  int single = 0, doubled = 0, nonzero = 0;
  for (int test = 0; test < 6; ++test) {
    Fixture f(rng);
    for (int n = 0; n < 4; ++n) {
      const auto &flow = plans.by_sites[1][n];
      for (int g = 1; g < 5; ++g) {
        const auto &other = plans.by_sites[1 << g][n];
        if (flow.wick != other.wick || flow.local[0] != other.local[0] ||
            flow.local[11] != other.local[11])
          throw logic_error("Wick orientation varies with position");
      }
      auto left = f.endpoint(flow, false), right = f.endpoint(flow, true);
      auto b = f.path(plans, n, 1), c = f.path(plans, n, 2);
      for (int i = 0; i <= 10; ++i)
        for (int j = 0; j <= 10 && i + j <= 11; ++j) {
          int blocks = (5 + i) / 7 + (5 + j) / 7;
          int v = coefficient(left, right, b[i], c[j], i + j);
          if (blocks == 2) {
            ++single;
            if (v % 7)
              throw logic_error("undivided path double pole");
          }
          if (blocks == 3) {
            ++doubled;
            if (v % 49)
              throw logic_error("undivided path triple pole");
          }
          if (blocks < 2 && v % 7)
            ++nonzero;
        }
    }
  }
  if (!nonzero)
    throw logic_error("zero controls");
  cout << single << " single-factor checks, " << doubled
       << " double-factor checks, " << nonzero << " nonzero controls\n";
}
