// Integer wave evaluations for the six-transverse-vector scalar class.
#define MIXED_FUNCTION_COMPONENTS
#include "mixed_function.cpp"

template <int Q> Waves<Q> nonzero_waves() {
  for (;;) {
    auto w = random_waves<Q>();
    if (w.cov[Q] != Vec{} &&
        all_of(w.pol.begin(), w.pol.end(), [](Vec u) { return u != Vec{}; }))
      return w;
  }
}

Waves<6> fields(int mask, int selector) {
  Waves<6> w;
  w.cov = {{{1, 0, 0}, {1, 1, 1}, {0, 1, 0}, {1, 1, -1},
            {0, 0, 1}, {1, 0, 1}, {1, 0, -1}, {}}};
  if (selector == 1)
    w.cov[6] = {0, 1, 1};
  for (int i = 0; i < 6; ++i) {
    int axis = i / 2, other = (axis + 1 + ((mask >> i) & 1)) % 3;
    w.pol[i][other] = 1;
    w.pol[i][axis] = -w.cov[i][other];
  }
  return w;
}

void line(int tag, const vector<Wide> &row) {
  cout << tag;
  for (auto x : row)
    cout << ' ' << big(x);
  cout << '\n';
}

Waves<6> conserving_fields(int mask, bool five_lines) {
  const array<Vec, 6> directions = five_lines
      ? array<Vec,6>{{{1,0,0},{1,0,0},{0,1,0},{0,0,1},{1,1,1},{1,-1,2}}}
      : array<Vec,6>{{{1,0,0},{0,1,0},{0,0,1},{0,1,1},{1,1,0},{1,0,1}}};
  const array<int,6> signs = five_lines ? array<int,6>{-1,1,-1,-1,1,-1}
                                        : array<int,6>{1,1,-1,1,-1,-1};
  Waves<6> w;
  for (int i = 0; i < 6; ++i) {
    auto p = directions[i];
    int pivot = 0;
    while (!p[pivot]) ++pivot;
    int axis = (pivot + 1 + ((mask >> i) & 1)) % 3;
    w.cov[i] = scale(p, signs[i]);
    w.cov[6] = sum(w.cov[6], scale(w.cov[i], -1));
    w.pol[i][axis] = p[pivot];
    w.pol[i][pivot] = -p[axis];
  }
  return w;
}

#ifndef MIXED_SIX_COMPONENTS
int main(int argc, char **argv) {
  bool balanced = argc == 2 && string(argv[1]) == "--balanced";
  bool linear = argc == 2 && string(argv[1]) == "--linear";
  bool geometries = argc == 2 && string(argv[1]) == "--geometries";
  if (argc > 1 && !balanced && !linear && !geometries)
    return 2;
  vector<Waves<5>> five;
  vector<Waves<6>> six;
  vector<Waves<7>> seven;
  exact = true;
  for (int i = 0; i < 110; ++i) {
    five.push_back(nonzero_waves<5>());
    six.push_back(nonzero_waves<6>());
  }
  for (int i = 0; i < 45; ++i)
    seven.push_back(nonzero_waves<7>());

  exact = false;
  auto source = basis(five, {75, 51, 21, 3}, true);
  auto target = basis(six, {70, 28, 7}, true);
  if (linear) {
    vector<Wide> degrees;
    for (const auto &column : target) {
      int degree = 0;
      for (int k : column.front().det) {
        for (int j=0;j<3;++j) {
          degree += (k % 8 == 6);
          k /= 8;
        }
      }
      for (int k : column.front().pair)
        degree += (k % 8 == 6);
      degrees.push_back(degree);
    }
    line(3, degrees);
  }
  exact = true;
  for (const auto &w : seven)
    line(0, differential(target, w));
  if (balanced) {
    for (int i = 0; i < 34; ++i) {
      auto w = nonzero_waves<6>();
      w.cov[6] = {};
      for (int j = 0; j < 6; ++j)
        w.cov[6] = sum(w.cov[6], scale(w.cov[j], -1));
      auto values = differential(source, w);
      Evaluations<6> e(w);
      for (const auto &col : target)
        values.push_back(integer_value(e, col, false));
      line(1, values);
    }
    return 0;
  }
  for (int selector = 0; selector < (geometries ? 4 : 2); ++selector) {
    for (int mask = 0; mask < 64; ++mask) {
      auto w = selector < 2 ? fields(mask, selector) : conserving_fields(mask, selector == 3);
      auto values = differential(source, w);
      Evaluations<6> e(w);
      for (const auto &col : target)
        values.push_back(integer_value(e, col, false));
      line(selector + 1, values);
    }
  }
}
#endif
