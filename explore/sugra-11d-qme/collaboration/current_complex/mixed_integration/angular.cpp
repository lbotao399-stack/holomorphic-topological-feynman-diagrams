// Degree (9,2) symbol.  The two-dimensional conductance functional is evaluated
// mod 49.
#define SUGRA_PHYSICAL_WORDS_ONLY
#include "../../five_theta/check_physical.cpp"
#include <numeric>

using namespace std;

int residue(int x) { return (x % 49 + 49) % 49; }
Scalar angular(const FlowPlans &plans, int g, int sing, int src, int ref,
               int wB, int wC) {
  array<V, 12> p{};
  p[g + 1][2] = 1;
  p[src + 1][3] = 1;
  p[ref + 1][2] = -1;
  p[ref + 1][3] = -1;
  array<V, 13> y{};
  int w[3] = {1 - wB - wC, wB, wC};
  for (int k = 0; k < 3; ++k) {
    V current{};
    current[0] = w[k] - (k == 1);
    current[1] = w[k] - (k == 2);
    for (size_t i = 0; i < edgeids[k].size(); ++i) {
      y[edgeids[k][i]] = current;
      if (i < vertices[k].size())
        for (int a = 0; a < 5; ++a)
          current[a] += p[vertices[k][i]][a];
    }
  }
  // The physical words use inv6=841 at 1009. Lift their signed cubic
  // coefficient to +/-1/6 and then use 41, the inverse of six modulo 49.
  array<array<V, 14>, 12> k{};
  for (int v = 0; v < 12; ++v) {
    k[v][13] = p[v];
    for (int e = 0; e < 13; ++e)
      if (::ends[e][0] == v || ::ends[e][1] == v)
        for (int a = 0; a < 5; ++a)
          k[v][e][a] = (::ends[e][1] == v ? 1 : -1) * y[e][a];
  }
  int total = 0;
  for (auto &flow : plans.by_sites[1 << g]) {
    array<Tensor, 12> tv{};
    for (int v = 0; v < 12; ++v) {
      int axis = v == g + 1 ? 4 : v == sing + 1 ? 1 : 0;
      for (auto &word : flow.local[v]->words) {
        if (word.external_axis >= 0 && word.external_axis != axis)
          continue;
        int val = k[v][word.first_slot][word.first_axis];
        if (word.second_slot >= 0)
          val *= k[v][word.second_slot][word.second_axis];
        int coefficient =
            flow.local[v]->cubic ? (word.scale / inv6) * 41 : word.scale;
        tv[v][word.element] += coefficient * val;
      }
      for (auto &v0 : tv[v])
        v0 = residue(v0);
    }
    Tensor acc = tv[0];
    for (int h = 0; h < 3; ++h) {
      int stride = h == 0 ? 25 : h == 1 ? 5 : 1;
      for (int v : vertices[h]) {
        Tensor next{};
        for (int n = 0; n < 125; ++n)
          if (acc[n]) {
            int i = n / stride % 5;
            for (int j = 0; j < 5; ++j)
              next[n + (j - i) * stride] += acc[n] * tv[v][5 * i + j];
          }
        for (auto &a : next)
          a %= 49;
        acc = next;
      }
    }
    int sum = 0;
    for (int i = 0; i < 125; ++i)
      sum += acc[i] * tv[11][i];
    total += flow.wick * sum;
  }
  return residue(total);
}
int main() {
  topology({1, 6, 6});
  FlowPlans plans(1);
  vector<array<int, 3>> samples;
  array<int, 3> p;
  while (cin >> p[0] >> p[1] >> p[2])
    samples.push_back(p);
  if (samples.empty())
    return 1;
  int nnz = 0;
  for (int g = 0; g < 10; ++g) {
    int ref = g == 0 ? 1 : 0;
    for (int j = 0; j < 10; ++j)
      if (j != g)
        for (int i = 0; i < 10; ++i)
          if (i != g && i != ref) {
            int sum = 0;
            for (auto [x, y, w] : samples)
              sum += w * angular(plans, g, j, i, ref, x, y);
            int v = residue(sum);
            if (v % 7)
              throw runtime_error("unnormalized 7-adic pole");
            v /= 7;
            if (v)
              ++nnz;
            if (v)
              cout << g << " " << j << " " << i << " " << ref << " " << v
                   << "\n";
          }
  }
  cerr << samples.size()
       << " interpolation nodes; nonzero angular coefficients " << nnz
       << " / 720\n";
  return nnz != 0;
}
