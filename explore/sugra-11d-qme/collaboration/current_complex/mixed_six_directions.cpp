// Exact scale and marked-momentum coefficients on three and four ordinary lines.
#define MIXED_SIX_COMPONENTS
#include "mixed_six_probe.cpp"
#include <map>

using Powers = uint64_t;
Powers variable(int i) { return Powers{1} << (4 * i); }

struct SymbolTerm {
  int sign;
  Powers scales = 0;
  array<int, 9> ordinary{};
  array<int, 3> marked{};
  int ordinary_size = 0, marked_size = 0;
};

vector<vector<SymbolTerm>> symbols(const vector<Column<6>> &columns) {
  vector<vector<SymbolTerm>> result;
  for (const auto &column : columns) {
    vector<SymbolTerm> terms;
    for (const auto &product : column) {
      SymbolTerm term{product.sign};
      auto add = [&](int index, initializer_list<int> covectors) {
        bool marked = false;
        for (int i : covectors) {
          if (i == 6) marked = true;
          else term.scales += variable(i);
        }
        if (marked) term.marked.at(term.marked_size++) = index;
        else term.ordinary.at(term.ordinary_size++) = index;
      };
      for (int i : product.det) add(i, {i / 64, i / 8 % 8, i % 8});
      for (int i : product.pair) add(512 + i, {i % 8});
      terms.push_back(term);
    }
    result.push_back(std::move(terms));
  }
  return result;
}

void chart(const array<Vec,6> &directions, int mask, int tag,
           const vector<vector<SymbolTerm>> &terms,
           const vector<Column<6>> &columns) {
  Waves<6> w;
  for (int i = 0; i < 6; ++i) {
    w.cov[i] = directions[i];
    int pivot = 0;
    while (!directions[i][pivot]) ++pivot;
    int axis = (pivot + 1 + ((mask >> i) & 1)) % 3;
    w.pol[i][axis] = directions[i][pivot];
    w.pol[i][pivot] = -directions[i][axis];
  }
  // Each elementary determinant and pairing is linear in the marked momentum
  // if it contains that momentum. Three coordinate evaluations recover it.
  array<array<int,560>,3> values{};
  for (int k = 0; k < 3; ++k) {
    w.cov[6] = {};
    w.cov[6][k] = 1;
    Evaluations<6> e(w);
    for (size_t i = 0; i < e.det.size(); ++i) values[k][i] = e.det[i].v;
    for (size_t i = 0; i < e.pair.size(); ++i) values[k][512+i] = e.pair[i].v;
  }

  map<Powers, vector<int>> coefficients;
  for (size_t j = 0; j < terms.size(); ++j) {
    for (const auto &t : terms[j]) {
      int value = t.sign;
      for (int i = 0; i < t.ordinary_size && value; ++i)
        value *= values[0][t.ordinary[i]];
      if (!value) continue;
      auto expand = [&](auto &self, int i, int coeff, Powers powers) -> void {
        if (i == t.marked_size) {
          auto &row = coefficients[powers];
          if (row.empty()) row.resize(terms.size());
          row[j] += coeff;
          return;
        }
        for (int k = 0; k < 3; ++k)
          if (int z = values[k][t.marked[i]])
            self(self, i + 1, coeff * z, powers + variable(6 + k));
      };
      expand(expand, 0, value, t.scales);
    }
  }

  const array<int,9> test{1,-2,3,2,-1,3,-2,1,1};
  vector<Wide> check(terms.size());
  for (const auto &[powers,row] : coefficients) {
    if (none_of(row.begin(), row.end(), [](int x) { return x; })) continue;
    int monomial = 1;
    for (int i = 0; i < 9; ++i)
      for (int e = (powers >> (4 * i)) & 15; e; --e) monomial *= test[i];
    cout << tag;
    for (size_t j = 0; j < terms.size(); ++j) {
      check[j] += monomial * row[j];
      cout << ' ' << row[j];
    }
    cout << '\n';
  }
  for (int i = 0; i < 6; ++i) w.cov[i] = scale(w.cov[i], test[i]);
  copy(test.begin() + 6, test.end(), w.cov[6].begin());
  Evaluations<6> evaluation(w);
  for (size_t j = 0; j < terms.size(); ++j)
    if (integer_value(evaluation, columns[j], false) != check[j])
      throw logic_error("scale polynomial failed direct substitution");
}

int main() {
  exact = true;
  vector<Waves<5>> five;
  vector<Waves<6>> six;
  for (int i = 0; i < 110; ++i) {
    five.push_back(nonzero_waves<5>());
    six.push_back(nonzero_waves<6>());
  }
  for (int i = 0; i < 45; ++i) nonzero_waves<7>(); // match mixed_six_probe's basis
  exact = false;
  auto source = basis(five, {75,51,21,3}, true);
  auto target = basis(six, {70,28,7}, true);
  exact = true;
  for (int i = 0; i < 85; ++i) six.push_back(nonzero_waves<6>());
  for (const auto &w : six) {
    Evaluations<6> e(w);
    vector<Wide> row;
    for (const auto &c : target) row.push_back(integer_value(e, c, false));
    auto boundary = differential(source, w);
    row.insert(row.end(), boundary.begin(), boundary.end());
    line(0, row);
  }
  const auto terms = symbols(target);
  Vec x{1,0,0}, y{0,1,0}, z{0,0,1}, u{1,1,1};
  for (int mask = 0; mask < 64; ++mask) {
    chart({x,x,y,y,z,z}, mask, 1, terms, target);
    chart({x,x,y,y,z,u}, mask, 2, terms, target);
  }
}
