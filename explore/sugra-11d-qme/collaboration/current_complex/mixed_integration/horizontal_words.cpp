// Literal cubic Wick words at fixed active Schwinger data and one linear field.
#define SUGRA_PHYSICAL_WORDS_ONLY
#include "../../five_theta/check_physical.cpp"
#include <memory>
#include <numeric>
#include <random>

using Shape = std::array<int, 3>;

int inverse_field(int x) {
  int n = P - 2, result = 1;
  for (x = mod(x); n; n >>= 1, x = mod(x * x))
    if (n & 1)
      result = mod(result * x);
  return result;
}

int loop_determinant(const std::array<V, 10> &columns,
                     const std::array<int, 10> &path) {
  std::array<std::array<int, 10>, 10> a{};
  for (int j = 0; j < 10; ++j)
    for (int i = 0; i < 5; ++i) {
      if (path[j] != 2)
        a[i][j] = mod(columns[j][i]);
      if (path[j] != 0)
        a[i + 5][j] = mod(columns[j][i]);
    }
  int answer = 1;
  for (int i = 0; i < 10; ++i) {
    int k = i;
    while (k < 10 && !a[k][i])
      ++k;
    if (k == 10)
      return 0;
    if (k != i) {
      std::swap(a[k], a[i]);
      answer = -answer;
    }
    int pivot = a[i][i];
    answer = mod(answer * pivot);
    for (int j = i + 1; j < 10; ++j) {
      int s = mod(a[j][i] * inverse_field(pivot));
      for (int m = i; m < 10; ++m)
        a[j][m] = mod(a[j][m] - s * a[i][m]);
    }
  }
  return answer;
}

struct Case {
  static constexpr int length = 12;
  Shape size, w = {2, -2, 1};
  unsigned gamma_labels;
  std::array<std::vector<int>, 3> ordered, time;
  std::array<V, 10> p{}, col{};
  std::array<int, 10> polar{};
  std::array<V, 3> initial{};
  Case(Shape n, unsigned mask, std::mt19937 &rng)
      : size(n), gamma_labels(mask) {
    std::array<int, 9> labels = {0, 1, 2, 3, 4, 5, 6, 7, 9};
    std::shuffle(labels.begin(), labels.end(), rng);
    for (int l = 0; l < 9; ++l) {
      polar[l] = rng() % 5;
      for (int i = 0; i < 5; ++i)
        if (i != polar[l]) {
          col[l][i] = int(rng() % 3) - 1;
          if (l != 8) {
            p[l][i] = length * col[l][i];
            p[9][i] -= p[l][i];
          }
        }
    }
    polar[9] = rng() % 5;
    col[9] = p[9];
    int j = 0;
    for (int h = 0; h < 3; ++h)
      for (int k = 0; k < n[h]; ++k) {
        ordered[h].push_back(labels[j++]);
        time[h].push_back(length - 1 - 2 * k);
      }
    std::array<V, 3> means{};
    V total{};
    for (int h = 0; h < 3; ++h)
      for (int j = 0; j < n[h]; ++j)
        for (int i = 0; i < 5; ++i) {
          int value = p[ordered[h][j]][i] * time[h][j] / length;
          means[h][i] += value;
          total[i] += value;
        }
    for (int h = 0; h < 3; ++h)
      for (int i = 0; i < 5; ++i)
        initial[h][i] = w[h] * total[i] - means[h][i];
  }
  int inserted(int h) {
    if (size[h] == 5)
      return 0;
    Shape shape = size;
    for (auto &n : shape)
      ++n;
    ++shape[h];
    vertices = {};
    edgeids = {};
    topology(shape);
    FlowPlans plans(__builtin_popcount(gamma_labels));
    int contribution = 0;
    for (int gap = 0; gap <= size[h]; ++gap) {
      auto order = ordered;
      order[h].insert(order[h].begin() + gap, 8);
      std::array<V, 12> waves{};
      std::array<int, 12> axes{};
      std::array<V, 13> current{};
      std::array<int, 10> paths{};
      std::vector<int> full, odd;
      int gamma_sites = 0;
      for (int r = 0; r < 3; ++r) {
        V k = initial[r];
        for (size_t j = 0; j < edgeids[r].size(); ++j) {
          current[edgeids[r][j]] = k;
          if (j < order[r].size()) {
            int v = vertices[r][j], l = order[r][j];
            paths[l] = r;
            waves[v] = p[l];
            axes[v] = polar[l];
            full.push_back(l);
            if ((gamma_labels >> l) & 1)
              gamma_sites |= 1 << (v - 1);
            else
              odd.push_back(l);
            for (int i = 0; i < 5; ++i)
              k[i] += p[l][i];
          }
        }
      }
      int det = loop_determinant(col, paths);
      int primitive =
          raw_numerator(plans.by_sites[gamma_sites], current, axes, waves);
      int width = (gap ? time[h][gap - 1] : length) -
                  (gap == size[h] ? 0 : time[h][gap]);
      int ext = signature(full) * signature(odd);
      contribution = mod(contribution + long(ext * width) * det * primitive);
    }
    return mod(contribution * inverse_field(w[h]));
  }
};
int main() {
  std::mt19937 rng(528264);
  for (int m = 1; m <= 5; ++m) {
    int checks = 0, nonzero = 0;
    for (int i = 0; i <= 5; ++i)
      for (int j = 0; j <= 5; ++j) {
        int k = 9 - i - j;
        if (k < 0 || k > 5)
          continue;
        for (int r = 0; r < 12; ++r) {
          unsigned mask = 1 << 9;
          // Keep label eight as the linear vector and label nine as a one-form.
          // The other ordinary legs already obey transverse polarization.
          std::array<int, 8> extra{0, 1, 2, 3, 4, 5, 6, 7};
          if (m > 1)
            std::shuffle(extra.begin(), extra.end(), rng);
          for (int t = 0; t < 2 * m - 2; ++t)
            mask |= 1 << extra[t];
          Case test({i, j, k}, mask, rng);
          int sum = 0;
          for (int h = 0; h < 3; ++h) {
            int v = test.inserted(h);
            sum = mod(sum + v);
            nonzero += v != 0;
          }
          if (sum)
            throw std::logic_error("one-linear-field correlation");
          ++checks;
        }
      }
    if (!nonzero)
      throw std::logic_error("vacuous sample");
    std::cout << "g^" << m << " " << checks << " active geometries; " << nonzero
              << " nonzero path choices; all summed zero\n";
  }
}
