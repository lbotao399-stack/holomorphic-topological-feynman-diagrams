// Evidence matrices for homogeneous jets of the marked transverse function.
#define MIXED_FUNCTION_COMPONENTS
#include "mixed_function.cpp"
#include <fstream>
#include <string>

namespace {
void line(int label, const vector<Wide> &row) {
  std::cout << label;
  for (auto x : row)
    std::cout << ' ' << big(x);
  std::cout << '\n';
}
void append(vector<Wide> &row, const Evaluations<5> &e,
            const vector<Column<5>> &cols) {
  for (const auto &col : cols)
    row.push_back(integer_value(e, col, false));
}
Waves<5> fields(bool single, int mask, int selector, int eta) {
  Waves<5> w;
  w.cov = {{{1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 1, 1}, {0, 0, 1}, {}, {}}};
  if (single) {
    w.cov[1] = {1, 1, 1};
    w.cov[3] = {1, 1, -1};
    w.cov[5] = {1, 0, -1};
  } else
    w.cov[5][selector % 3] = selector == 3 ? 2 : 1;
  for (int j = 0; j < 5; j++) {
    int axis = j / 2, other = (axis + 1 + ((mask >> j) & 1)) % 3;
    w.pol[j][other] = 1;
    w.pol[j][axis] = -w.cov[j][other];
  }
  w.cov[6][eta] = 1;
  return w;
}
} // namespace
#ifndef MIXED_PROBE_JETS_LIBRARY
int main(int argc, char **argv) {
  if (argc != 3)
    throw std::runtime_error("usage: mixed_probe_jets old.tsv single.tsv");
  exact = true;
  vector<Waves<4>> four;
  vector<Waves<5>> five;
  vector<Waves<6>> six;
  while (four.size() < 110)
    four.push_back(random_waves<4>());
  while (five.size() < 110)
    five.push_back(random_waves<5>());
  while (six.size() < 65)
    six.push_back(random_waves<6>());
  exact = false;
  auto inc = basis(
      four, {57, 49,
             26}); // keep the original RNG order for both differential matrices
  auto middle = basis(five, {84, 0, 0});
  auto target = basis(five, {84, 45, 9});
  exact = true;
  for (const auto &w : six)
    line(0, differential(middle, w));
  for (auto w : five) {
    array<vector<Wide>, 3> vals;
    auto f = w.cov[5];
    for (int z = 1; z <= 3; z++) {
      w.cov[5] = scale(f, z);
      auto &v = vals[z - 1];
      v = differential(inc, w);
      Evaluations<5> e(w);
      append(v, e, middle);
      append(v, e, target);
      for (auto &a : v) {
        if (a % z)
          throw std::logic_error("function origin");
        a /= z;
      }
    }
    array<vector<Wide>, 3> coefficients;
    for (size_t j = 0; j < vals[0].size(); j++) {
      Wide a = vals[0][j], b = vals[1][j], c = vals[2][j];
      Wide high = c - 2 * b + a;
      if (high % 2)
        throw std::logic_error("second difference");
      high /= 2;
      Wide mid = b - a - 3 * high;
      coefficients[0].push_back(a - mid - high);
      coefficients[1].push_back(mid);
      coefficients[2].push_back(high);
    }
    for (int s = 0; s < 3; s++)
      line(s + 1, coefficients[s]);
  }
  for (int dataset = 0; dataset < 2; dataset++) {
    std::ifstream tests(argv[dataset + 1]);
    mpz_class weight;
    int mask, selector = 0, eta, count = 0;
    if (!tests)
      throw std::runtime_error("missing detector weights");
    while (tests >> weight >> mask) {
      if (!dataset)
        tests >> selector;
      tests >> eta;
      if (!tests || mask < 0 || mask >= 32 || eta < 0 || eta >= 3)
        throw std::runtime_error("invalid detector row");
      auto w = fields(dataset, mask, selector, eta);
      auto row = differential(inc, w);
      append(row, Evaluations<5>(w), middle);
      std::cout << (dataset + 4) << ' ' << weight;
      for (auto x : row)
        std::cout << ' ' << big(x);
      std::cout << '\n';
      count++;
    }
    if (count != (dataset ? 94 : 89) || !tests.eof())
      throw std::runtime_error("incomplete detector weights");
  }
}
#endif
