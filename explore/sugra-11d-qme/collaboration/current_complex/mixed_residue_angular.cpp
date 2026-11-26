// Exact cancellation of the top conductance symbol of theta166.
#define SUGRA_PHYSICAL_WORDS_ONLY
#define EXACT_INTEGER
#include "../five_theta/check_physical.cpp"
#include <numeric>

V addv(V a, V b, int factor = 1) {
  for (int i = 0; i < 5; ++i)
    a[i] = mod(a[i] + factor * b[i]);
  return a;
}
std::array<long long, 4> raw(const FlowPlans &plans, int j, int wb, int wc,
                             V prefix, V path_sum, V mean) {
  const V m{1, 0, 0, 0, 0}, pg{0, 0, 1, 0, 0};
  const int gamma = vertices[1][j];
  std::array<V, 12> p{};
  p[gamma] = pg;
  if (j > 0)
    p[vertices[1][0]] = prefix;
  else if (prefix != V{})
    throw std::runtime_error("nonempty gamma prefix");
  V tail = addv(addv(path_sum, prefix, -1), pg, -1);
  if (j < 4)
    p[vertices[1][j + 1]] = tail;
  else if (tail != V{})
    throw std::runtime_error("nonempty gamma suffix");
  p[vertices[2][0]] = addv({}, path_sum, -1);
  // Every external vector is e_0; the external one-form is e^1.
  std::array<int, 12> pol{};
  pol[gamma] = 1;
  std::array<V, 13> current{};
  for (int k = 0; k < 3; ++k) {
    int w = k == 0 ? -wb - wc : k == 1 ? wb : wc;
    V v = addv(addv({}, m, w), mean, k == 0 ? 0 : k == 1 ? -1 : 1);
    for (size_t l = 0; l < edgeids[k].size(); ++l) {
      current[edgeids[k][l]] = v;
      if (l < vertices[k].size())
        v = addv(v, p[vertices[k][l]]);
    }
  }
  std::array<long long, 4> values{};
  int i = 0;
  for (const auto &flow : plans.by_sites[1 << (gamma - 1)])
    values[i++] = raw_numerator({flow}, current, pol, p);
  return values;
}
int main() {
  topology({1, 6, 6});
  FlowPlans plans(1);
  V pg{0, 0, 1, 0, 0}, a{0, 0, 0, 1, 0}, b{0, 0, 0, 0, 1};
  int total = 0, nonzero = 0;
  for (int j = 0; j < 5; ++j)
    for (int pair = 0; pair < 3; ++pair) {
      if ((j == 0 || j == 4) && pair != 2)
        continue;
      V prefix = pair < 2 ? a : V{},
        sum = pair == 0   ? b
              : pair == 2 ? a
                          : V{},
        mean = pair ? b : V{};
      if (j == 4)
        prefix = addv(sum, pg, -1);
      for (int t = -5; t <= 5; ++t) {
        auto values = raw(plans, j, t, 1, prefix, sum, mean);
        if (std::accumulate(values.begin(), values.end(), 0LL) != 0)
          throw std::runtime_error("degree-ten symbol survives");
        for (auto v : values)
          if (v)
            ++nonzero;
        ++total;
      }
    }
  if (nonzero == 0)
    throw std::runtime_error("all individual orientation tests vanished");
  std::cout
      << total << " integer interpolation tests, " << nonzero
      << " nonzero individual orientations, zero summed degree-ten symbol\n";
}
