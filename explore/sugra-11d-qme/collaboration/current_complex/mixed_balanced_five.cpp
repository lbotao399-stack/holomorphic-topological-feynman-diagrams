// Integer values of the five-vector current, its balanced restriction, and its first normal jet.
#define MIXED_PROBE_JETS_LIBRARY
#include "mixed_probe_jets.cpp"

int main(int argc, char **argv) {
  if (argc != 2)
    throw std::runtime_error("expected the 94-row detector weights");
  exact = true;
  std::vector<Waves<5>> five;
  std::vector<Waves<6>> six;
  for (int i = 0; i < 169; ++i) {
    five.push_back(random_waves<5>());
    six.push_back(random_waves<6>());
  }
  exact = false;
  auto marked = basis(five, {84, 45, 9});
  auto unmarked = basis(five, {75}, false, 0);
  exact = true;

  // Restrict the unmarked target a second time, balancing its five vector momenta.
  for (auto waves : five) {
    Vec total{};
    for (int j = 0; j < 4; ++j)
      total = sum(total, waves.cov[j]);
    waves.cov[4] = scale(total, -1);
    waves.pol[4] = cross(waves.cov[4], random_vec());
    std::vector<Wide> row;
    append(row, Evaluations<5>(waves), unmarked);
    line(3, row);
  }
  for (const auto &waves : six)
    line(0, differential(marked, waves));
  for (auto waves : five) {
    Vec total{};
    for (int j = 0; j < 5; ++j)
      total = sum(total, waves.cov[j]);
    waves.cov[5] = scale(total, -1);
    std::vector<Wide> values;
    append(values, Evaluations<5>(waves), marked);
    line(1, values);
    waves.dcov[5] = random_vec();
    Evaluations<5> first(waves);
    for (size_t j = 0; j < marked.size(); ++j)
      values[j] = integer_value(first, marked[j], true) - values[j];
    line(4, values);
  }

  std::ifstream src(argv[1]);
  mpz_class coefficient;
  int mask, eta, count = 0;
  std::vector<mpz_class> target(marked.size());
  while (src >> coefficient >> mask >> eta) {
    Evaluations<5> eval(fields(true, mask, 0, eta));
    for (size_t j = 0; j < marked.size(); ++j)
      target[j] += coefficient * big(integer_value(eval, marked[j], false));
    ++count;
  }
  if (count != 94 || !src.eof())
    throw std::runtime_error("incomplete detector weights");
  std::cout << 2;
  for (const auto &x : target)
    std::cout << " " << x;
  std::cout << '\n';
}
