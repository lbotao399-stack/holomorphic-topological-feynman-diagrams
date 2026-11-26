// Positive-jet tests on the transverse five-vector factor of the graph.
#define MIXED_FUNCTION_COMPONENTS
#include "mixed_function.cpp"
#include <fstream>

struct AxisTest {
  mpz_class coefficient;
  int mask, f, eta;
  Waves<5> waves() const {
    Waves<5> w;
    w.cov = {{{1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 1, 1}, {0, 0, 1}, {}, {}}};
    for (int j = 0; j < 5; ++j) {
      int axis = j / 2, k = (axis + 1 + ((mask >> j) & 1)) % 3;
      w.pol[j][k] = 1;
      w.pol[j][axis] = -w.cov[j][k];
    }
    w.cov[5][f % 3] = (f == 3 ? 2 : 1);
    w.cov[6][eta] = 1;
    return w;
  }
};

vector<vector<int>> kernel(const Echelon &rows) {
  vector<vector<int>> result;
  for (size_t n = 0; n < rows.piv.size(); ++n)
    if (rows.piv[n].empty()) {
      vector<int> v(rows.piv.size());
      v[n] = 1;
      for (int i = int(v.size()) - 1; i >= 0; --i)
        if (!rows.piv[i].empty()) {
          for (size_t j = i + 1; j < v.size(); ++j)
            v[i] = mod(v[i] - rows.piv[i][j] * v[j]);
        }
      result.push_back(std::move(v));
    }
  return result;
}
int main(int argc, char **argv) {
  if (argc > 1)
    prime = stoi(argv[1]);
  if (prime != 251 && prime != 257)
    throw runtime_error("choose 251 or 257");
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
  auto incoming = basis(four, {57, 49, 26});
  auto degree_five = basis(five, {84, 0, 0});
  exact = true;
  Echelon outgoing(84);
  for (const auto &w : six) {
    vector<int> row;
    for (Wide x : differential(degree_five, w))
      row.push_back(int((x % prime + prime) % prime));
    exact = false;
    outgoing.add(std::move(row));
    exact = true;
  }
  if (outgoing.rank != 37)
    throw runtime_error("outgoing rank");
  exact = false;
  auto closed = kernel(outgoing);
  ifstream input(argc > 2 ? argv[2] : "mixed_graph_probe.tsv");
  if (!input)
    throw runtime_error("cannot open mixed graph probe");
  vector<mpz_class> boundary(incoming.size());
  vector<int> values(degree_five.size());
  AxisTest t;
  int count = 0;
  while (input >> t.coefficient >> t.mask >> t.f >> t.eta) {
    if (t.mask < 0 || t.mask >= 32 || (t.f != 0 && t.f != 1 && t.f != 3) ||
        t.eta < 0 || t.eta >= 3)
      throw runtime_error("invalid probe row");
    exact = true;
    auto w = t.waves();
    auto row = differential(incoming, w);
    for (size_t j = 0; j < row.size(); ++j)
      boundary[j] += t.coefficient * big(row[j]);
    exact = false;
    Evaluations<5> eval(w);
    int c = mpz_fdiv_ui(t.coefficient.get_mpz_t(), prime);
    for (size_t j = 0; j < values.size(); ++j)
      values[j] = mod(values[j] + c * value(eval, degree_five[j]));
    ++count;
  }
  if (!input.eof() || count != 89)
    throw runtime_error("incomplete probe");
  if (any_of(boundary.begin(), boundary.end(), [](auto &v) { return v != 0; }))
    throw runtime_error("nonzero integer boundary");
  for (auto &v : closed) {
    int pairing = 0;
    for (size_t j = 0; j < v.size(); ++j)
      pairing = mod(pairing + v[j] * values[j]);
    if (pairing) {
      cout << "89 terms on one five-vector momentum list, 132 rational "
              "boundary evaluations zero, closed pairing modulo "
           << prime << ": " << pairing << "\n";
      return 0;
    }
  }
  throw runtime_error("probe killed all closed cochains");
}
