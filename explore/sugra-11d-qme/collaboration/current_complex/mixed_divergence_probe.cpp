// Test strict transverse conservation on the one-mark current line.
#define MIXED_PROBE_JETS_LIBRARY
#include "mixed_probe_jets.cpp"

int main(int argc, char **argv) {
  bool graded = argc == 3 && string(argv[2]) == "--graded";
  if (argc != 2 && !graded)
    throw runtime_error("expected single detector table");
  exact = true;
  vector<Waves<5>> five;
  vector<Waves<6>> six;
  for (int k = 0; k < 172; ++k) {
    five.push_back(random_waves<5>());
    six.push_back(random_waves<6>());
  }
  exact = false;
  auto cols = basis(five, {84, 45, 9});
  vector<Waves<4>> four;
  for (int i = 0; i < 100; ++i)
    four.push_back(random_waves<4>());
  auto scalars = basis(four, {31, 36, 22, 10, 1}, true);
  auto scalar_five = basis(five, {75, 51, 21, 3}, true);
  exact = true;
  cout << 0;
  for (auto &col : cols) {
    auto &r = col.front();
    int degree = 0;
    for (int j : r.det) {
      for (int k = 0; k < 3; ++k) {
        degree += (j % 7 == 5);
        j /= 7;
      }
    }
    for (int j : r.pair)
      degree += (j % 7 == 5);
    cout << " " << degree;
  }
  cout << "\n";
  for (auto &w : six) {
    line(1, differential(cols, w));
    line(5, differential(scalar_five, w));
  }
  for (auto w : five) {
    Vec f = w.cov[5], p{};
    for (int i = 0; i < 5; ++i)
      p = sum(p, w.cov[i]);
    array<vector<Wide>, 4> vals;
    array<int, 4> scales{1, -1, 2, -2};
    for (int k = 0; k < (graded ? 4 : 1); ++k) {
      w.cov[5] = scale(f, scales[k]);
      w.cov[6] = sum(p, w.cov[5]);
      append(vals[k], Evaluations<5>(w), cols);
      auto sc = differential(scalars, w);
      if (k == 0) {
        line(2, vals[0]);
        line(4, sc);
      }
      vals[k].insert(vals[k].end(), sc.begin(), sc.end());
    }
    if (graded) {
      array<vector<Wide>, 4> coeffs;
      for (size_t j = 0; j < vals[0].size(); ++j) {
        auto a = vals[0][j], b = vals[1][j], c = vals[2][j], d = vals[3][j];
        Wide c3 = (c - d - 2 * (a - b)) / 12, c4 = (c + d - 4 * (a + b)) / 24;
        Wide c1 = (a - b) / 2 - c3, c2 = (a + b) / 2 - c4;
        for (int i = 0; i < 4; ++i) {
          int t = scales[i];
          if (t * (c1 + t * (c2 + t * (c3 + t * c4))) != vals[i][j])
            throw logic_error("marked degree");
        }
        coeffs[0].push_back(c1);
        coeffs[1].push_back(c2);
        coeffs[2].push_back(c3);
        coeffs[3].push_back(c4);
      }
      for (int s = 0; s < 4; ++s)
        line(11 + s, coeffs[s]);
    }
  }
  vector<mpz_class> wgt(cols.size());
  ifstream fd(argv[1]);
  int mask, eta, count = 0;
  mpz_class weight;
  while (fd >> weight >> mask >> eta) {
    auto w = fields(true, mask, 0, eta);
    Evaluations<5> e(w);
    for (size_t j = 0; j < cols.size(); ++j)
      wgt[j] += weight * big(integer_value(e, cols[j], false));
    ++count;
  }
  if (count != 94)
    throw runtime_error("incomplete detector");
  cout << 3;
  for (auto &x : wgt)
    cout << " " << x;
  cout << "\n";
  return 0;
}
