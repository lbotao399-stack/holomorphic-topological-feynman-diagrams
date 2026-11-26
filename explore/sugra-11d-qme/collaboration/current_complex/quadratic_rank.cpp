// Rank lower bounds for the quadratic bracket on one isotypic component.
#define EXTREME_LIBRARY
#include "extreme_profile.cpp"

// Each exterior output coordinate, including its cubic label, has a fixed
// signed image in the narrow test space. This map can only lower rank.
uint64_t mix_coordinate(uint64_t x) {
  x ^= x >> 30;
  x *= 0xbf58476d1ce4e5b9ULL;
  x ^= x >> 27;
  x *= 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}
int main() {
  Family f;
  int q, dim, amount, trace;
  Weight w;
  while (cin >> q >> dim >> amount >> trace >> w[0] >> w[1] >> w[2] >> w[3] >>
         w[4]) {
    Seeds seeds(f, w, q);
    if (int(seeds.seeds.size()) != trace || q < 2 || q > 8)
      throw runtime_error("inconsistent weight dimension");
    const int width = dim + 48;
    Rank inputs(width), outputs(width);
    vector<vector<pair<int, int>>> columns;
    columns.reserve(seeds.seeds.size());
    for (auto s : seeds.seeds) {
      auto a = decode(s);
      vector<pair<int, int>> col;
      for (int i = 0; i < q; i++)
        for (int j = i + 1; j < q; j++) {
          Word rest = 0;
          for (int k = 0; k < q; k++)
            if (k != i && k != j)
              rest = rest * 128 + a[k];
          for (auto [b, c] : f.brackets[a[i]][a[j]]) {
            auto hash = mix_coordinate(rest * 256 + b);
            col.push_back(
                {(hash / 2) % width, ((i + j + (hash & 1)) % 2 ? -1 : 1) * c});
          }
        }
      columns.push_back(std::move(col));
    }
    for (int n = 0; n < dim + 40; n++) {
      Vec seed;
      for (int r = 0; r < min(16, trace); r++)
        seed[seeds.seeds[seeds.rng() % trace]] = 1 + seeds.rng() % (P - 1);
      auto v = project(f, w, seed);
      for (int i = 0; i < 4; i++)
        if (!f.apply(v, i, i + 1).empty())
          throw runtime_error("failed highest vector");
      vector<int> x(width);
      for (int i = 0; i < min(width, trace); i++)
        x[i] = v[seeds.seeds[i]];
      inputs.insert(x);
      vector<int64_t> out(width);
      for (int i = 0; i < trace; i++) {
        auto it = v.find(seeds.seeds[i]);
        if (it == v.end())
          continue;
        for (auto [t, c] : columns[i])
          out[t] += it->second * c;
      }
      for (int i = 0; i < width; i++)
        x[i] = mod(out[i] % P);
      outputs.insert(x);
      if (inputs.rank == dim)
        break;
    }
    if (inputs.rank != dim)
      throw runtime_error("highest weights did not span");
    cout << q;
    for (int x : w)
      cout << " " << x;
    cout << " " << dim << " " << amount << " " << outputs.rank << endl;
  }
}
