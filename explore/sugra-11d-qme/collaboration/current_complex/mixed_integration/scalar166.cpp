// Integrate every Schwinger parameter at a single scalar test in the detector.
#define main mixed_graph_point_program
#include "../mixed_graph_point.cpp"
#undef main
#include <map>
using namespace std;
int coefficient_modulus = 1009;
int reduce_mod(long long x) {
  return int((x % coefficient_modulus + coefficient_modulus) %
             coefficient_modulus);
}
int exact_det(array<array<int, 10>, 10> a) {
  array<array<long long, 10>, 10> t{};
  for (int i = 0; i < 10; ++i)
    for (int j = 0; j < 10; ++j)
      t[i][j] = a[i][j];
  long long prev = 1, sign = 1;
  for (int k = 0; k < 9; ++k) {
    int j = k;
    while (j < 10 && t[j][k] == 0)
      ++j;
    if (j == 10)
      return 0;
    if (j != k) {
      swap(t[j], t[k]);
      sign = -sign;
    }
    long long v = t[k][k];
    for (int i = k + 1; i < 10; ++i)
      for (int h = k + 1; h < 10; ++h)
        t[i][h] = (t[i][h] * v - t[i][k] * t[k][h]) / prev;
    prev = v;
  }
  return int(sign * t[9][9]);
}
struct Term {
  int pos, a, b, coeff;
};
struct Local {
  vector<Term> terms;
};
struct Prepared {
  array<array<Local, 12>, 4> local;
  array<int, 4> sign{};
  Prepared(const Ordered &graph, const Fields &f, const array<V, 10> &p) {
    for (int n = 0; n < 4; ++n) {
      const auto &flow = graph.flows[n];
      sign[n] = flow.wick;
      for (int v = 0; v < 12; ++v) {
        map<array<int, 3>, int> terms;
        for (auto w : flow.local[v]->words) {
          int c = flow.local[v]->cubic ? w.scale / inv6 *
                                             (coefficient_modulus == 49  ? 41
                                              : coefficient_modulus == 7 ? 6
                                                                         : 841)
                                       : w.scale;
          int l = v >= 1 && v <= 10 ? graph.order[v - 1] : -1;
          if (w.external_axis >= 0)
            c *= l == 9 ? int(w.external_axis == 2 + f.eta)
                        : f.polar[l][w.external_axis];
          if (!c)
            continue;
          array<int, 2> keys{65, 65};
          for (int i = 0; i < 2; ++i) {
            int slot = i ? w.second_slot : w.first_slot;
            int axis = i ? w.second_axis : w.first_axis;
            if (slot < 0)
              continue;
            if (slot == 13)
              c *= p[l][axis];
            else {
              keys[i] = 5 * slot + axis;
              if (::ends[slot][0] == v)
                c = -c;
            }
          }
          sort(keys.begin(), keys.end());
          terms[{w.element, keys[0], keys[1]}] += c;
        }
        for (auto [key, c] : terms)
          if (reduce_mod(c)) {
            auto [pos, a, b] = key;
            local[n][v].terms.push_back({pos, a, b, reduce_mod(c)});
          }
      }
    }
  }
  int evaluate(const array<int, 66> &y) const {
    int ans = 0;
    for (int n = 0; n < 4; ++n) {
      array<Tensor, 12> tv{};
      for (int v = 0; v < 12; ++v) {
        for (auto &w : local[n][v].terms)
          tv[v][w.pos] += reduce_mod(1LL * w.coeff * y[w.a] * y[w.b]);
        for (int &a : tv[v])
          a %= coefficient_modulus;
      }
      array<Matrix, 2> mats{};
      for (int h = 1; h <= 2; ++h) {
        Matrix m{};
        for (int i = 0; i < 5; ++i)
          m[i * 6] = 1;
        for (int v : vertices[h]) {
          Matrix next{};
          for (int i = 0; i < 5; ++i)
            for (int k = 0; k < 5; ++k)
              if (m[5 * i + k])
                for (int j = 0; j < 5; ++j)
                  next[5 * i + j] += m[5 * i + k] * tv[v][5 * k + j];
          for (int &a : next)
            a %= coefficient_modulus;
          m = next;
        }
        mats[h - 1] = m;
      }
      for (int i = 0; i < 5; ++i) {
        Matrix mid{};
        for (int b = 0; b < 5; ++b)
          for (int c = 0; c < 5; ++c)
            if (int x = tv[0][25 * i + 5 * b + c])
              for (int h = 0; h < 5; ++h)
                mid[5 * h + c] += x * mats[0][5 * b + h];
        for (int &a : mid)
          a %= coefficient_modulus;
        for (int j = 0; j < 5; ++j)
          for (int k = 0; k < 5; ++k)
            if (int end = tv[11][25 * i + 5 * j + k]) {
              int t = 0;
              for (int c = 0; c < 5; ++c)
                t += mid[5 * j + c] * mats[1][5 * c + k];
              ans = reduce_mod(ans +
                               1LL * sign[n] * end * (t % coefficient_modulus));
            }
      }
    }
    return ans;
  }
};
struct Rule {
  array<array<int, 6>, 2> x;
  int w;
};
int main(int argc, char **argv) {
  if (argc != 8 && argc != 9) {
    cerr << "usage: scalar166 resist.txt conduct.txt order branch mask "
            "selector eta\n";
    return 2;
  }
  if (argc == 9)
    coefficient_modulus = stoi(argv[8]);
  if (coefficient_modulus != 49 && coefficient_modulus != 7 &&
      coefficient_modulus != 1009)
    throw runtime_error("unsupported residue ring");
  vector<Rule> rule;
  ifstream r(argv[1]);
  Rule rn{};
  while (r >> rn.x[0][0]) {
    for (int i = 1; i < 12; ++i)
      r >> rn.x[i / 6][i % 6];
    r >> rn.w;
    rule.push_back(rn);
  }
  vector<array<int, 3>> conduct;
  ifstream a(argv[2]);
  array<int, 3> q{};
  while (a >> q[0] >> q[1] >> q[2])
    conduct.push_back(q);
  if (rule.empty() || conduct.empty())
    throw runtime_error("empty cubature rule");
  array<int, 3> shape{1, 6, 6};
  topology(shape);
  FlowPlans flows(1);
  array<int, 10> labels{};
  string s = argv[3];
  if (s.size() != 10)
    throw runtime_error("bad order");
  for (int i = 0; i < 10; ++i)
    labels[i] = s[i] - '0';
  auto check = labels;
  sort(check.begin(), check.end());
  for (int i = 0; i < 10; ++i)
    if (check[i] != i)
      throw runtime_error("invalid label permutation");
  Ordered graph(shape, labels, flows);
  int branch = stoi(argv[4]);
  Fields field(branch, stoi(argv[5]), stoi(argv[6]), stoi(argv[7]), 1);
  array<V, 10> p{};
  array<int, 9> scales{1, 2, 3, 1, 1, 1, 1, 1, 1};
  for (int i = 0; i < 9; ++i)
    for (int j = 0; j < 5; ++j)
      p[i][j] =
          scales[i] * field.direction[i][j] + (i == 3 ? field.function[j] : 0);
  for (int i = 0; i < 9; ++i)
    for (int j = 0; j < 5; ++j)
      p[9][j] -= p[i][j];
  array<array<int, 10>, 10> rows{};
  for (int i = 0; i < 10; ++i)
    for (int j = 0; j < 5; ++j) {
      if (graph.paths[i] != 2)
        rows[j][i] = p[i][j];
      if (graph.paths[i] != 0)
        rows[j + 5][i] = p[i][j];
    }
  int det = reduce_mod(graph.sign * exact_det(rows));
  Prepared prepared(graph, field, p);
  array<V, 13> prefix{};
  for (int h = 0; h < 3; ++h) {
    V sum{};
    for (size_t e = 0; e < edgeids[h].size(); ++e) {
      prefix[edgeids[h][e]] = sum;
      if (e < vertices[h].size())
        for (int j = 0; j < 5; ++j)
          sum[j] += p[labels[vertices[h][e] - 1]][j];
    }
  }
  int total = 0;
  for (auto &r : rule) {
    array<V, 3> m{};
    for (int h = 1; h < 3; ++h)
      for (int e = 0; e < 6; ++e)
        for (int j = 0; j < 5; ++j)
          m[h][j] += r.x[h - 1][e] * prefix[edgeids[h][e]][j];
    V all{};
    for (int j = 0; j < 5; ++j)
      all[j] = reduce_mod(m[1][j] + m[2][j]);
    int angular = 0;
    for (auto [wB, wC, mass] : conduct) {
      if ((mass * r.w) % coefficient_modulus == 0)
        continue;
      array<int, 3> w{reduce_mod(1 - wB - wC), wB, wC};
      array<int, 66> y{};
      y[65] = 1;
      for (int h = 0; h < 3; ++h)
        for (int e : edgeids[h])
          for (int j = 0; j < 5; ++j)
            y[5 * e + j] = reduce_mod(w[h] * all[j] - m[h][j] + prefix[e][j]);
      angular = reduce_mod(angular + mass * prepared.evaluate(y));
    }
    total = reduce_mod(total + angular * r.w);
  }
  cout << "D with external sign " << det << ", supplied-rule Wick " << total
       << ", product " << reduce_mod(det * total) << " mod "
       << coefficient_modulus << "; evaluations "
       << rule.size() * conduct.size() << "\n";
}
