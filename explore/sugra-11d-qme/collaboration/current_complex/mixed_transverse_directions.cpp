// Exact restriction of the five-vector transverse-current basis to three lines.
#define MIXED_FUNCTION_COMPONENTS
#include "mixed_function.cpp"
#include <map>

struct Monomial {
  int coefficient = 0;
  uint32_t powers = 0;
};

void direction_chart(const array<int, 5> &axes, int type, int mask, int eta,
                     const vector<Column<5>> &cols) {
  // Four bits per exponent for the five wave scales followed by the three
  // components of the marked function momentum. Every exponent is at most 5.
  array<int, 7> axis{};
  copy(axes.begin(), axes.end(), axis.begin());
  axis[5] = -1;
  axis[6] = eta;
  array<int, 5> pol{};
  for (int i = 0; i < 5; ++i)
    pol[i] = (axis[i] + 1 + ((mask >> i) & 1)) % 3;
  array<Monomial, 343> determinants{};
  array<Monomial, 35> pairings{};
  auto factor = [](int input, int coordinate) -> uint32_t {
    if (input == 6)
      return 0;
    int variable = input == 5 ? 5 + coordinate : input;
    return 1u << (4 * variable);
  };
  for (int i = 0; i < 7; ++i)
    for (int j = 0; j < 7; ++j)
      for (int k = 0; k < 7; ++k) {
        if (i == j || i == k || j == k)
          continue;
        bool has_function = i == 5 || j == 5 || k == 5;
        for (int coordinate = 0; coordinate < (has_function ? 3 : 1); ++coordinate) {
          int x = i == 5 ? coordinate : axis[i];
          int y = j == 5 ? coordinate : axis[j];
          int z = k == 5 ? coordinate : axis[k];
          if (x == y || y == z || x == z)
            continue;
          auto &value = determinants[49 * i + 7 * j + k];
          if (value.coefficient)
            throw logic_error("more than one monomial in an axial determinant");
          value = {(y - x + 3) % 3 == 1 ? 1 : -1,
                   factor(i, x) + factor(j, y) + factor(k, z)};
        }
      }
  for (int j = 0; j < 5; ++j)
    for (int i = 0; i < 7; ++i)
      if (i == 5 || axis[i] == pol[j])
        pairings[7 * j + i] = {1, factor(i, pol[j])};

  map<uint32_t, vector<int>> coefficients;
  for (size_t j = 0; j < cols.size(); ++j)
    for (const auto &term : cols[j]) {
      int sign = term.sign;
      uint32_t powers = 0;
      for (auto i : term.det) {
        sign *= determinants[i].coefficient;
        powers += determinants[i].powers;
      }
      for (auto i : term.pair) {
        sign *= pairings[i].coefficient;
        powers += pairings[i].powers;
      }
      if (sign) {
        auto &row = coefficients[powers];
        if (row.empty())
          row.resize(cols.size());
        row[j] += sign;
      }
    }
  const array<int, 8> fixture{1, -2, 3, 2, -1, 3, -2, 1};
  vector<Wide> check(cols.size());
  for (const auto &[powers, row] : coefficients) {
    if (none_of(row.begin(), row.end(), [](int x) { return x != 0; }))
      continue;
    int monomial = 1;
    for (int j = 0; j < 8; ++j)
      for (int exponent = (powers >> (4 * j)) & 15; exponent; --exponent)
        monomial *= fixture[j];
    cout << 1 + type;
    for (size_t j = 0; j < cols.size(); ++j) {
      check[j] += monomial * row[j];
      cout << ' ' << row[j];
    }
    cout << '\n';
  }
  Waves<5> w;
  for (int i = 0; i < 5; ++i) {
    w.cov[i][axis[i]] = fixture[i];
    w.pol[i][pol[i]] = 1;
  }
  copy(fixture.begin() + 5, fixture.end(), w.cov[5].begin());
  w.cov[6][eta] = 1;
  Evaluations<5> eval(w);
  for (size_t j = 0; j < cols.size(); ++j)
    if (check[j] != integer_value(eval, cols[j], false))
      throw logic_error("monomial substitution");
}

void four_lines(const vector<Column<4>> &primitives,
                const vector<Column<5>> &currents) {
  const array<array<int, 8>, 5> fixtures{{
      {-3, 2, 3, 3,-3, 2,-1, 1}, { 1, 2,-3, 2, 2, 1,-1,-2},
      {-1,-3, 3, 2,-2, 0, 0, 2}, { 1, 1, 1, 1, 1, 0, 1, 0},
      { 1,-1,-2, 2, 1, 0, 1,-1}
  }};
  for (const auto &data : fixtures)
    for (int mask = 0; mask < 32; ++mask) {
      if (data[0] == data[1] && (mask & 3) != 2)
        continue;
      for (int eta = 0; eta < 3; ++eta) {
        Waves<5> w;
        w.cov[0] = w.cov[1] = {1, 0, 0};
        w.cov[2] = {0, 1, 0};
        w.cov[3] = {0, 0, 1};
        w.cov[4] = {1, 1, 1};
        for (int i = 0; i < 5; ++i)
          w.cov[i] = scale(w.cov[i], data[i]);
        w.pol[0][1 + (mask & 1)] = 1;
        w.pol[1][1 + ((mask >> 1) & 1)] = 1;
        w.pol[2][mask & 4 ? 2 : 0] = 1;
        w.pol[3][mask & 8 ? 1 : 0] = 1;
        w.pol[4] = mask & 16 ? Vec{1, 0, -1} : Vec{1, -1, 0};
        copy(data.begin() + 5, data.end(), w.cov[5].begin());
        w.cov[6][eta] = 1;
        Evaluations<5> eval(w);
        cout << 3;
        for (Wide x : differential(primitives, w))
          cout << ' ' << big(x);
        for (const auto &c : currents)
          cout << ' ' << big(integer_value(eval, c, false));
        cout << '\n';
      }
    }
  for (int i = 0; i < 130; ++i) {
    cout << 4;
    for (Wide x : differential(currents, random_waves<6>()))
      cout << ' ' << big(x);
    cout << '\n';
  }
}

int main() {
  exact = true;
  vector<Waves<4>> four;
  vector<Waves<5>> five;
  while (four.size() < 110)
    four.push_back(random_waves<4>());
  while (five.size() < 238)
    five.push_back(random_waves<5>());
  exact = false;
  auto primitives = basis(four, {57, 49, 26});
  auto currents = basis(five, {84, 45, 9});
  exact = true;
  for (const auto &w : five) {
    cout << 0;
    Evaluations<5> eval(w);
    for (const auto &c : currents)
      cout << ' ' << big(integer_value(eval, c, false));
    for (Wide x : differential(primitives, w))
      cout << ' ' << big(x);
    cout << '\n';
  }
  for (int mask = 0; mask < 32; ++mask)
    for (int eta = 0; eta < 3; ++eta) {
      direction_chart({0, 0, 1, 1, 2}, 0, mask, eta, currents);
      direction_chart({0, 0, 0, 1, 2}, 1, mask, eta, currents);
    }
  four_lines(primitives, currents);
}
