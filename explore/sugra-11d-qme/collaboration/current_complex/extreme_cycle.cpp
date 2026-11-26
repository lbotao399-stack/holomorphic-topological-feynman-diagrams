// Recover the unique surviving highest vector and verify its integral boundary.
#define EXTREME_LIBRARY
#include "extreme_profile.cpp"
using Exact = unordered_map<Word, int64_t>;
const Weight highest{5, -3, -3, -3, -4};
void residue() {
  Family f;
  Seeds s(f, highest);
  constexpr int width = 320, limit = 34;
  vector<Probe> probes;
  for (int i = 0; i < width; i++)
    probes.push_back(select_probe(f, s));
  vector<vector<int>> piv(width), sources(width);
  vector<Vec> seeds;
  for (int j = 0; j < limit; j++) {
    Vec seed;
    for (int k = 0; k < 4; k++)
      seed[s.seeds[s.rng() % s.seeds.size()]] = 1 + s.rng() % (P - 1);
    seeds.push_back(seed);
    auto v = project(f, highest, seed);
    vector<int> row, word(limit);
    word[j] = 1;
    for (auto &p : probes)
      row.push_back(p.eval(v));
    bool independent = false;
    for (int i = 0; i < width; i++)
      if (row[i]) {
        if (piv[i].empty()) {
          int c = inv(row[i]);
          for (int &x : row)
            x = x * c % P;
          for (int &x : word)
            x = x * c % P;
          piv[i] = row;
          sources[i] = word;
          independent = true;
          break;
        }
        int c = P - row[i];
        for (int k = i; k < width; k++)
          row[k] = (row[k] + c * piv[i][k]) % P;
        for (int k = 0; k <= j; k++)
          word[k] = (word[k] + c * sources[i][k]) % P;
      }
    if (independent)
      continue;
    Vec combined;
    for (int k = 0; k <= j; k++)
      if (word[k])
        for (auto [m, c] : seeds[k])
          combined[m] = mod(combined[m] + word[k] * c);
    Vec cycle = project(f, highest, combined);
    if (cycle.empty())
      continue;
    map<Word, int> ordered(cycle.begin(), cycle.end());
    int scalar = inv(ordered.begin()->second);
    for (auto [m, c] : ordered)
      cout << m << " " << c * scalar % P << "\n";
    return;
  }
  throw runtime_error("no kernel residue");
}
void verify(const string &file) {
  Family f;
  Exact cycle;
  Word key;
  int64_t coefficient;
  while (cin >> key >> coefficient) {
    auto s = decode(key);
    Weight w{};
    for (int k = 0; k < 8; k++) {
      if (s[k] >= 70 || (k && s[k] <= s[k - 1]))
        throw runtime_error("invalid exterior monomial");
      for (int t = 0; t < 5; t++)
        w[t] += f.low.weights[s[k]][t];
    }
    if (w != highest || !coefficient || !cycle.emplace(key, coefficient).second)
      throw runtime_error("invalid highest-weight vector");
  }
  if (cycle.size() != 32676)
    throw runtime_error("unexpected cycle support");
  for (int i = 0; i < 4; i++)
    if (!f.apply(cycle, i, i + 1, false, false).empty())
      throw runtime_error("raising map failed");
  Exact qq;
  for (auto [m, c] : cycle) {
    auto a = decode(m);
    for (int i = 0; i < 8; i++)
      for (int j = i + 1; j < 8; j++) {
        Word rest = 0;
        for (int k = 0; k < 8; k++)
          if (k != i && k != j)
            rest = rest * 128 + a[k];
        for (auto [b, d] : f.brackets[a[i]][a[j]])
          qq[rest * 256 + b] += ((i + j) % 2 ? -1 : 1) * c * d;
      }
  }
  for (auto [m, c] : qq)
    if (c)
      throw runtime_error("quadratic bracket failed");
  cout << "quadratic/cubic output coordinates " << qq.size() << " all zero\n";
  // Evaluate the high-bracket intertwiner on every coordinate of each of its
  // three highest vectors. Equivariance then kills every input to that map.
  struct Term {
    int input, coeff;
    vector<int> word;
  };
  array<vector<Term>, 3> terms;
  ifstream in(file);
  if (!in)
    throw runtime_error("missing dual-bracket coefficients");
  int h, i, t, c, n;
  while (in >> h >> i >> t >> c >> n) {
    Term term{i, c, {}};
    for (int j = 0; j < n; j++) {
      int root;
      in >> root;
      term.word.push_back(root);
    }
    if (!t)
      terms.at(h).push_back(std::move(term));
  }
  array<Weight, 3> weights{
      {{6, -3, -3, -3, -4}, {5, -2, -3, -3, -4}, {5, -3, -3, -3, -3}}};
  for (int h = 0; h < 3; h++) {
    Seeds outputs(f, weights[h], 7);
    if (terms[h].empty())
      throw runtime_error("empty bracket channel");
    for (auto m : outputs.seeds) {
      auto a = decode(m);
      int64_t value = 0;
      for (auto &term : terms[h]) {
        int pos = lower_bound(a.begin(), a.begin() + 7, term.input) - a.begin();
        if (a[pos] == term.input)
          continue;
        auto b = a;
        b[7] = term.input;
        sort(b.begin(), b.end());
        Exact dual{{encode(b), 1}};
        for (auto it = term.word.rbegin(); it != term.word.rend(); ++it)
          dual = f.apply(dual, *it + 1, *it, true, false);
        for (auto [k, c] : dual) {
          auto hit = cycle.find(k);
          if (hit != cycle.end())
            value += (pos % 2 ? -1 : 1) * hit->second * c * term.coeff;
        }
      }
      if (value)
        throw runtime_error("quadratic/high bracket failed");
    }
    cout << "quadratic/high output " << h << " " << outputs.seeds.size()
         << " coordinates all zero\n";
  }
}
int main(int argc, char **argv) {
  if (argc == 2 && string(argv[1]) == "residue")
    residue();
  else if (argc == 3 && string(argv[1]) == "verify")
    verify(argv[2]);
  else
    throw runtime_error("use residue or verify DUAL_BRACKET_FILE");
}
