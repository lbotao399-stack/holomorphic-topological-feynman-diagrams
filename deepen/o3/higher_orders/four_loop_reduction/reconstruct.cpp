// Exact evaluation of the chain-integrated numerator on the quadratic
// interpolation grid.
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>
using namespace std;
using Int = long long;
using Wide = __int128;
using Nine = array<int, 9>;
struct Term {
  array<int, 4> j;
  array<int, 4> h;
  int sign;
};
struct Pair {
  int a, b;
  bool right = false;
};
struct Product {
  int left, right;
  Int coefficient;
};
struct Evaluation {
  array<int, 4> j;
  int sign;
  vector<Product> products;
};
const int qsize = 15;
struct Graph {
  array<array<int, 4>, 9> c;
  array<array<Int, 5>, 18> offset;
  array<array<Int, qsize>, qsize> moment;
  array<array<int, 5>, 5> symmetric{};
  vector<array<int, 4>> bases;
  struct Cofactor {
    int i, j, k;
    array<int, 9> minor;
  };
  vector<Cofactor> cofactors;
  vector<Pair> pairs;
  vector<Evaluation> terms;
  int orientation;
  Graph() {
    for (auto &e : c)
      for (auto &v : e)
        cin >> v;
    for (auto &row : offset)
      for (auto &v : row)
        cin >> v;
    cin >> orientation;
    int t = 0;
    for (int i = 0; i < 5; i++)
      for (int j = i; j < 5; j++)
        symmetric[i][j] = symmetric[j][i] = t++;
    for (auto &row : moment)
      for (auto &v : row)
        cin >> v;
    auto minor = [&](array<int, 4> cols) {
      array<int, 4> p{0, 1, 2, 3};
      int total = 0;
      do {
        int value = 1;
        for (int i = 0; i < 4; i++) {
          value *= c[cols[i]][p[i]];
          for (int j = 0; j < i; j++)
            if (p[j] > p[i])
              value = -value;
        }
        total += value;
      } while (next_permutation(p.begin(), p.end()));
      if (total < -1 || total > 1)
        throw runtime_error("non-unimodular cycle minor");
      return total;
    };
    for (int i = 0; i < 9; i++)
      for (int j = i + 1; j < 9; j++)
        for (int k = j + 1; k < 9; k++) {
          Cofactor cof{i, j, k};
          for (int e = 0; e < 9; e++)
            cof.minor[e] = minor({e, i, j, k});
          cofactors.push_back(cof);
          for (int l = k + 1; l < 9; l++)
            if (minor({i, j, k, l}))
              bases.push_back({i, j, k, l});
        }
    map<pair<int, int>, int> cache;
    auto pair_id = [&](int a, int b, bool right) {
      if (a > b)
        swap(a, b);
      auto [it, fresh] = cache.emplace(make_pair(a, b), pairs.size());
      if (fresh)
        pairs.push_back({a, b});
      pairs[it->second].right |= right;
      return it->second;
    };
    int n;
    cin >> n;
    for (int e = 0; e < n; e++) {
      Term raw;
      for (auto &i : raw.j)
        cin >> i;
      for (auto &i : raw.h)
        cin >> i;
      cin >> raw.sign;
      Evaluation term{raw.j, raw.sign};
      for (int m = 0; m < 16; m++) {
        Int c0;
        cin >> c0;
        if (!c0)
          continue;
        int indices[4];
        for (int k = 0; k < 4; k++)
          indices[k] = 2 * raw.h[k] + ((m >> k) & 1);
        term.products.push_back({pair_id(indices[0], indices[1], false),
                                 pair_id(indices[2], indices[3], true), c0});
      }
      terms.push_back(term);
    }
    if (!cin)
      throw runtime_error("incomplete graph inputs");
  }
  Wide value(Nine r) const {
    Int U = 0;
    for (auto t : bases)
      U += r[t[0]] * r[t[1]] * r[t[2]] * r[t[3]];
    Int D[9][9]{};
    for (auto &cof : cofactors) {
      int m = r[cof.i] * r[cof.j] * r[cof.k];
      for (int e = 0; e < 9; e++)
        if (cof.minor[e])
          for (int f = 0; f < 9; f++)
            D[e][f] += m * cof.minor[e] * cof.minor[f];
    }
    array<array<Int, 5>, 18> h{};
    Int largest = 0;
    for (int i = 0; i < 9; i++)
      for (int k = 0; k < 2; k++)
        for (int t = 0; t < 5; t++) {
          auto &z = h[2 * i + k][t];
          z = U * offset[2 * i + k][t];
          for (int j = 0; j < 9; j++)
            z -= D[i][j] * r[j] * offset[2 * j + k][t];
          largest = max(largest, abs(z));
        }
    // Bounds for the quadratic products (2 h^2), the moment transform
    // (15*8!*2 h^2) and their wide inner products are enforced at each node.
    if (largest > 1000000)
      throw runtime_error("horizontal coefficient out of 64-bit range");
    vector<array<Int, qsize>> q(pairs.size()), weighted(pairs.size());
    for (size_t n = 0; n < pairs.size(); n++) {
      auto &pair = pairs[n];
      auto &left = h[pair.a];
      auto &right = h[pair.b];
      for (int i = 0; i < 5; i++)
        for (int j = 0; j < 5; j++)
          q[n][symmetric[i][j]] += left[i] * right[j];
      if (pair.right)
        for (int i = 0; i < qsize; i++)
          for (int j = 0; j < qsize; j++)
            weighted[n][i] += moment[i][j] * q[n][j];
    }
    Wide value = 0;
    for (auto &term : terms) {
      const auto [a, b, c, d] = term.j;
      Int haf = D[a][b] * D[c][d] + D[a][c] * D[b][d] + D[a][d] * D[b][c];
      if (!haf)
        continue;
      Wide pol = 0;
      for (auto &p : term.products) {
        Wide m = 0;
        for (int i = 0; i < qsize; i++)
          m += Wide(q[p.left][i]) * weighted[p.right][i];
        Wide v;
        if (__builtin_mul_overflow(m, Wide(p.coefficient), &v) ||
            __builtin_add_overflow(pol, v, &pol))
          throw overflow_error("horizontal moment");
      }
      Wide v;
      if (__builtin_mul_overflow(pol, Wide(term.sign * haf), &v) ||
          __builtin_add_overflow(value, v, &value))
        throw overflow_error("Pfaffian sum");
    }
    Wide divisor = Wide(r[8]) * U * U * U;
    if (value % divisor)
      throw runtime_error("nonintegral primitive times 4*8!");
    return orientation * (value / divisor);
  }
};
void print(Wide v) {
  if (v < 0) {
    cout << '-';
    v = -v;
  }
  char out[50];
  int n = 0;
  do {
    out[n++] = char('0' + v % 10);
    v /= 10;
  } while (v);
  while (n)
    cout << out[--n];
  cout << '\n';
}
int main() {
  Graph g;
  Nine r;
  r.fill(1);
  for (int k = 0; k < 6561; k++) {
    int j = k;
    for (int i = 7; i >= 0; i--) {
      r[i] = 1 + j % 3;
      j /= 3;
    }
    print(g.value(r));
  }
}
