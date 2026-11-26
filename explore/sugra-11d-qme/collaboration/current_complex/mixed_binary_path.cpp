// Summed squarefree path moments on the actual forward, reverse and marked
// fibers.
#include <algorithm>
#include <array>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>
using V = std::array<int, 5>;
using Mat = std::array<int, 25>;
int mul[16][16];
std::mt19937 gen(92731);
void make_field() {
  for (int a = 0; a < 16; ++a)
    for (int b = 0; b < 16; ++b) {
      int v = 0, k = a;
      for (int j = 0; j < 4; ++j) {
        if (b & (1 << j))
          v ^= k;
        k <<= 1;
        if (k & 16)
          k ^= 19;
      }
      mul[a][b] = v;
    }
}
int inverse(int v) {
  for (int i = 1; i < 16; ++i)
    if (mul[v][i] == 1)
      return i;
  throw std::runtime_error("nonunit");
}
V plus(V a, V b) {
  for (int k = 0; k < 5; ++k)
    a[k] ^= b[k];
  return a;
}
int dot(V p, V u) {
  int s = 0;
  for (int k = 0; k < 5; ++k)
    s ^= mul[p[k]][u[k]];
  return s;
}
Mat product(const Mat &a, const Mat &b) {
  Mat c{};
  for (int i = 0; i < 5; ++i)
    for (int j = 0; j < 5; ++j)
      for (int k = 0; k < 5; ++k)
        c[5 * i + j] ^= mul[a[5 * i + k]][b[5 * k + j]];
  return c;
}
struct Tests {
  std::array<V, 5> momenta{}, u{};
  V start{}, slope{};
  Tests() {
    for (int i = 0; i < 5; ++i) {
      for (int k = 0; k < 5; ++k)
        momenta[i][k] = gen() % 16;
      int first = i, second = (i + 2) % 5;
      u[i][first] = momenta[i][second];
      u[i][second] = momenta[i][first];
    }
    for (int k = 0; k < 5; ++k) {
      start[k] = gen() % 16;
      slope[k] = gen() % 16;
    }
  }
  Mat transfer(int which, int i, int prefix, V k) const {
    Mat m{};
    auto p = momenta[i];
    auto pol = u[i];
    if (i == 0 && which >= 2) {
      if (which == 2) {
        for (int a = 0; a < 5; ++a)
          for (int b = 0; b < 5; ++b)
            m[5 * a + b] = mul[pol[a]][p[b]] ^ mul[pol[b]][p[a]];
      } else {
        std::array<int, 5> perm{0, 1, 2, 3, 4};
        do {
          m[5 * perm[0] + perm[1]] ^=
              mul[mul[pol[perm[2]]][k[perm[3]]]][p[perm[4]]];
        } while (std::next_permutation(perm.begin(), perm.end()));
      }
      return m;
    }
    bool fwd = which == 0 || (which == 2 && (prefix & 1)) ||
               (which == 3 && !(prefix & 1));
    for (int a = 0; a < 5; ++a)
      for (int b = 0; b < 5; ++b) {
        m[5 * a + b] = fwd ? mul[pol[a]][p[b]] : mul[p[a]][pol[b]];
        if (a == b)
          m[5 * a + b] ^= dot(k, pol);
      }
    return m;
  }
  Mat path(int which, int d, int t) const {
    V initial = start;
    for (int k = 0; k < 5; ++k)
      initial[k] ^= mul[t][slope[k]];
    std::array<V, 32> Q{};
    for (int s = 1; s < 32; ++s)
      Q[s] = plus(Q[s & (s - 1)], momenta[__builtin_ctz(unsigned(s))]);
    std::array<std::array<Mat, 5>, 32> trans{};
    for (int s = 0; s < 32; ++s)
      for (int j = 0; j < 5; ++j)
        if (!(s & (1 << j)))
          trans[s][j] = transfer(which, j, s, plus(initial, Q[s]));
    std::vector<Mat> state(1 << (5 + d));
    for (int i = 0; i < 5; ++i)
      state[0][5 * i + i] = 1;
    for (int s = 1; s < (1 << (5 + d)); ++s) {
      int slots = s & 31;
      for (int j = 0; j < 5; ++j)
        if (slots & (1 << j)) {
          int prev = s ^ (1 << j);
          auto m = product(state[prev], trans[prev & 31][j]);
          for (int k = 0; k < 25; ++k)
            state[s][k] ^= m[k];
        }
      for (int i = 0; i < d; ++i)
        if (s & (1 << (5 + i))) {
          auto &v = state[s ^ (1 << (5 + i))];
          int q = Q[slots][i];
          for (int k = 0; k < 25; ++k)
            state[s][k] ^= mul[q][v[k]];
        }
    }
    return state.back();
  }
};
int degree(const std::array<Mat, 8> &samples) {
  auto v = samples;
  int d = -1;
  for (int k = 0; k < 8; ++k) {
    if (std::any_of(v[0].begin(), v[0].end(), [](int n) { return n != 0; }))
      d = k;
    for (int j = 0; j < 7 - k; ++j) {
      int inv = inverse(j ^ (j + k + 1));
      for (int i = 0; i < 25; ++i)
        v[j][i] = mul[v[j][i] ^ v[j + 1][i]][inv];
    }
  }
  return d;
}
int main() {
  make_field();
  for (int type = 0; type < 4; ++type) {
    std::array<int, 3> maxdeg{-1, -1, -1};
    for (int trial = 0; trial < 4; ++trial) {
      Tests fixture;
      for (int d = 3; d <= 5; ++d) {
        std::array<Mat, 8> samples;
        for (int t = 0; t < 8; ++t)
          samples[t] = fixture.path(type, d, t);
        int actual = degree(samples);
        maxdeg[d - 3] = std::max(maxdeg[d - 3], actual);
        if (actual > __builtin_popcount(5 + d))
          throw std::runtime_error("opposite-current degree");
      }
    }
    std::cout << "fiber " << type << ", observed degrees " << maxdeg[0] << " "
              << maxdeg[1] << " " << maxdeg[2] << "\n";
  }
}
