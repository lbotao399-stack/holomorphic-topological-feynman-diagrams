// Finite-field lower bounds on the two final marked-scalar differentials.
#define MIXED_FUNCTION_COMPONENTS
#include "mixed_function.cpp"

template <int Q> Waves<Q> nonzero_waves() {
  for (;;) {
    auto w = random_waves<Q>();
    if (w.cov[Q] != Vec{} && all_of(w.pol.begin(), w.pol.end(),
                                    [](const Vec &u) { return u != Vec{}; }))
      return w;
  }
}

template <int Q, int MarkedDegree = -1> int lower_bound(int target, int probes, int dimension) {
  constexpr int n = Q + int(MarkedDegree != 0);
  exact = true;
  vector<Waves<Q + 1>> samples;
  for (int j = 0; j < probes; ++j)
    samples.push_back(nonzero_waves<Q + 1>());
  vector<array<int, 3>> triples;
  for (int a = 0; a < n; ++a)
    for (int b = a + 1; b < n; ++b)
      for (int c = b + 1; c < n; ++c)
        triples.push_back({a, b, c});
  exact = false;
  vector<Evaluations<Q>> source_probes;
  for (int j = 0; j < dimension + 16; ++j)
    source_probes.emplace_back(nonzero_waves<Q>());
  Echelon source(dimension + 16), rank(probes);
  int tested = 0;
  for (int attempt = 0; attempt < 1500000 && rank.rank < target; ++attempt) {
    array<int, n> degree{};
    array<array<int, 3>, 3> det;
    array<int, Q> pair;
    for (auto &t : det) {
      t = triples[rng() % triples.size()];
      for (auto i : t)
        ++degree[i];
    }
    for (int i = 0; i < Q; ++i) {
      pair[i] = (i + 1 + rng() % (n - 1)) % n;
      ++degree[pair[i]];
    }
    if constexpr (MarkedDegree < 0) {
      if (!degree[Q]) continue;
    } else if constexpr (MarkedDegree > 0) {
      if (degree[Q] != MarkedDegree) continue;
    }
    if (*min_element(degree.begin(), degree.begin() + Q) < 2)
      continue;
    vector<Column<Q>> one{alternate<Q>(det, pair)};
    exact = false;
    vector<int> evals;
    for (auto &e : source_probes)
      evals.push_back(value(e, one[0]));
    if (!source.add(std::move(evals)))
      continue;
    vector<int> values;
    for (auto &s : samples) {
      exact = true;
      auto row = differential(one, s);
      values.push_back(int((row[0] % prime + prime) % prime));
    }
    exact = false;
    rank.add(std::move(values));
    ++tested;
    if (tested == dimension)
      break;
  }
  if (rank.rank != target)
    throw std::runtime_error("rank certificate incomplete");
  cout << "arity " << Q << " columns " << tested << " certified rank at least "
       << rank.rank << " mod " << prime << '\n'
       << flush;
  return rank.rank;
}
int main(int argc, char **argv) {
  if (argc == 2 && string(argv[1]) == "--unmarked") {
    lower_bound<3, 0>(2, 8, 2);
    lower_bound<4, 0>(14, 26, 16);
    lower_bound<7, 0>(10, 24, 44);
    lower_bound<8, 0>(1, 3, 11);
    return 0;
  }
  if (argc == 2 && string(argv[1]) == "--current") {
    lower_bound<3, 1>(4, 12, 4);
    lower_bound<4, 1>(27, 48, 31);
    lower_bound<5, 1>(47, 85, 75);
    lower_bound<6, 1>(22, 45, 70);
    lower_bound<7, 1>(3, 12, 25);
    return 0;
  }
  if (argc != 1) return 2;
  lower_bound<2>(3, 10, 3);
  lower_bound<3>(27, 40, 30);
  lower_bound<7>(3, 10, 31);
  lower_bound<6>(28, 40, 105);
  return 0;
}
