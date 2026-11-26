// Integral five-vector transverse cycles for the mixed four-plus-five current.
#define MIXED_PROBE_JETS_LIBRARY
#include "../mixed_probe_jets.cpp"

namespace {
Waves<5> small_wave(int choice, int mask, int eta) {
  Waves<5> w;
  w.cov = {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {-1, -1, -1}, {-1, -1, 1}, {}, {}}};
  if (choice == 2)
    w.cov[4] = {1, 1, -1};
  const std::array<Vec, 3> functions{{{1, 0, 0}, {1, 0, -1}, {0, 1, -1}}};
  w.cov[5] = functions.at(choice);
  w.cov[6][eta] = 1;
  for (int i = 0; i < 5; ++i) {
    int pivot = 0;
    while (!w.cov[i][pivot])
      ++pivot;
    int axis = (pivot + 1 + ((mask >> i) & 1)) % 3;
    w.pol[i][axis] = w.cov[i][pivot];
    w.pol[i][pivot] = -w.cov[i][axis];
  }
  return w;
}
void accumulate(vector<mpz_class> &result, const mpz_class &weight,
                const Waves<5> &input, const vector<Column<4>> &primitives,
                const vector<Column<5>> &target) {
  auto row = differential(primitives, input);
  append(row, Evaluations<5>(input), target);
  for (size_t i = 0; i < result.size(); ++i)
    result[i] += weight * big(row[i]);
}
} // namespace

int main(int argc, char **argv) {
  if (argc != 3)
    throw std::invalid_argument("usage: four_small_cycle old-weights small-weights");
  exact = true;
  vector<Waves<4>> four;
  vector<Waves<5>> five;
  while (four.size() < 110)
    four.push_back(random_waves<4>());
  while (five.size() < 110)
    five.push_back(random_waves<5>());
  exact = false;
  auto primitives = basis(four, {57, 49, 26});
  auto target = basis(five, {84, 45, 9});
  exact = true;
  for (int j = 0; j < 100; ++j)
    line(0, differential(target, random_waves<6>()));
  std::ifstream old(argv[1]), small(argv[2]);
  if (!old || !small)
    throw std::runtime_error("missing polarization weights");
  vector<mpz_class> original(270), alternative(270);
  mpz_class weight;
  int mask, eta, choice;
  while (old >> weight >> mask >> eta)
    accumulate(original, weight, fields(true, mask, 0, eta), primitives, target);
  while (small >> weight >> choice >> mask >> eta)
    accumulate(alternative, weight, small_wave(choice, mask, eta), primitives, target);
  for (const auto *row : {&original, &alternative}) {
    std::cout << (row == &original ? 1 : 2);
    for (const auto &x : *row)
      std::cout << ' ' << x;
    std::cout << '\n';
  }
}
