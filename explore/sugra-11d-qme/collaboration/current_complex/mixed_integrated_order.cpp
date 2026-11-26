// Shape 166. Sum the complete positive-jet detector at one labeled vertex
// order.
#define main mixed_graph_point_program
#include "mixed_graph_point.cpp"
#undef main
#include <map>
#include <unordered_map>

namespace {
using Key = uint64_t;
using Poly = std::unordered_map<Key, mpz_class>;
using PMat = array<Poly, 25>;
using PTensor = array<Poly, 125>;
constexpr Key mask5 = (1ULL << 20) - 1;
const array<int, 10> labels{0, 1, 4, 6, 7, 2, 3, 5, 8, 9};
struct Probe {
  int mask, f, eta;
  mpz_class weight;
};
Key unit(int axis) { return Key{1} << (4 * axis); }
void plus(Poly &a, const Poly &b, const mpz_class &w = 1) {
  if (!w)
    return;
  for (const auto &[k, v] : b)
    mpz_addmul(a[k].get_mpz_t(), w.get_mpz_t(), v.get_mpz_t());
}
void times(Poly &a, const Poly &b, const Poly &c, int w = 1) {
  if (!w)
    return;
  for (const auto &[k, v] : b)
    if (v) {
      mpz_class t = v * w;
      for (const auto &[l, x] : c)
        if (x)
          mpz_addmul(a[k + l].get_mpz_t(), t.get_mpz_t(), x.get_mpz_t());
    }
}
void clean(Poly &a) {
  for (auto i = a.begin(); i != a.end();)
    if (!i->second)
      i = a.erase(i);
    else
      ++i;
}
int degree(Key k) {
  int d = 0;
  while (k) {
    d += k & 15;
    k >>= 4;
  }
  return d;
}
int bin[16][16];
mpz_class fac[20], odd[14], angular_den[14], measure_den, path_den[14][14];
void factorials() {
  fac[0] = odd[0] = angular_den[0] = 1;
  for (int i = 1; i < 20; i++)
    fac[i] = fac[i - 1] * i;
  for (int i = 1; i < 14; i++) {
    odd[i] = odd[i - 1] * (2 * i - 1);
    angular_den[i] = angular_den[i - 1] * (11 + 2 * i);
  }
  for (int i = 0; i < 16; i++) {
    bin[i][0] = bin[i][i] = 1;
    for (int j = 1; j < i; j++)
      bin[i][j] = bin[i - 1][j - 1] + bin[i - 1][j];
  }
  mpz_class common = 1;
  for (int d = 0; d <= 13; d++)
    for (int e = 0; e + d <= 13; e++) {
      mpz_class den = fac[5 + d] * fac[5 + e] / (fac[5] * fac[5]);
      mpz_lcm(common.get_mpz_t(), common.get_mpz_t(), den.get_mpz_t());
    }
  for (int d = 0; d <= 13; d++)
    for (int e = 0; e + d <= 13; e++)
      path_den[d][e] = common * fac[5] * fac[5] / (fac[5 + d] * fac[5 + e]);
  measure_den = common * angular_den[13];
}
struct Split {
  Key key;
  int d, c;
};
std::unordered_map<Key, vector<Split>> split_cache;
const vector<Split> &splits(Key k) {
  auto f = split_cache.find(k);
  if (f != split_cache.end())
    return f->second;
  vector<Split> result{{0, 0, 1}};
  for (int i = 0; i < 5; i++) {
    int a = (k >> (4 * i)) & 15;
    if (!a)
      continue;
    int len = result.size();
    for (int j = 1; j <= a; j++)
      for (int t = 0; t < len; t++) {
        auto part = result[t];
        result.push_back(
            {part.key + j * unit(i), part.d + j, bin[a][j] * part.c});
      }
  }
  return split_cache.emplace(k, std::move(result)).first->second;
}
std::map<array<int, 4>, mpz_class> conductances;
const mpz_class &conductance(int a, int b, int c, int d) {
  const array<int, 4> key{a, b, c, d};
  auto found = conductances.find(key);
  if (found != conductances.end())
    return found->second;
  mpz_class v = 0;
  for (int i = 0; i <= a; i++)
    for (int j = 0; j <= c; j++) {
      mpz_class term = odd[i + b] * odd[j + d] *
                       (angular_den[13] / angular_den[i + b + j + d]) *
                       bin[a][i] * bin[c][j];
      if ((a - i + c - j) & 1)
        v -= term;
      else
        v += term;
    }
  return conductances.emplace(key, std::move(v)).first->second;
}
struct MomentPair {
  Key left, right;
  mpz_class weight;
};
const vector<MomentPair> &mean_expansion(Key key) {
  static std::unordered_map<Key, vector<MomentPair>> cache;
  auto [entry, added] = cache.try_emplace(key);
  if (!added)
    return entry->second;
  Key b = key & mask5, c = key >> 20;
  int db = degree(b), dc = degree(c);
  std::unordered_map<Key, mpz_class> collected;
  for (auto u : splits(b))
    for (auto z : splits(c)) {
      Key left = b - u.key + c - z.key, right = u.key + z.key;
      int dr = u.d + z.d, dl = db + dc - dr;
      collected[left + (right << 20)] +=
          conductance(db - u.d, u.d, z.d, dc - z.d) * u.c * z.c *
          path_den[dl][dr];
    }
  auto &terms = entry->second;
  for (auto &[pair, weight] : collected)
    if (weight != 0)
      terms.push_back({pair & mask5, pair >> 20, std::move(weight)});
  return terms;
}
struct PathMoment {
  array<V, 6> q{};
  array<std::unordered_map<Key, mpz_class>, 6> homogeneous;
  std::unordered_map<Key, mpz_class> moments;
  mpz_class hom(int j, Key k) {
    if (!k)
      return 1;
    if (!j)
      return 0;
    auto &row = homogeneous[j];
    auto found = row.find(k);
    if (found != row.end())
      return found->second;
    mpz_class t = hom(j - 1, k);
    for (int i = 0; i < 5; i++)
      if (q[j][i] && ((k >> (4 * i)) & 15))
        t += q[j][i] * hom(j, k - unit(i));
    row[k] = t;
    return t;
  }
  const mpz_class &operator()(Key k) {
    auto f = moments.find(k);
    if (f != moments.end())
      return f->second;
    mpz_class v = hom(5, k);
    for (int i = 0; i < 5; i++)
      v *= fac[(k >> (4 * i)) & 15];
    return moments.emplace(k, std::move(v)).first->second;
  }
};
struct Geometry {
  array<V, 10> p{};
  array<array<Poly, 5>, 13> currents{};
  array<PathMoment, 3> moments;

  Geometry(const Fields &fields, const V &sc, int switches) {
    for (int i = 0; i < 9; i++) {
      int t = i < 4 ? sc[i] : (switches >> (i - 4)) & 1;
      for (int a = 0; a < 5; a++) {
        p[i][a] = t * fields.direction[i][a] +
                  (i == 3 && (switches & 32) ? fields.function[a] : 0);
        p[9][a] -= p[i][a];
      }
    }
    for (int h = 0; h < 3; h++)
      for (int j = 0; j < edgeids[h].size(); j++) {
        int e = edgeids[h][j];
        for (int i = 0; i < 5; i++) {
          auto &t = currents[e][i];
          if (h != 2)
            t[unit(i)] = h ? 1 : -1;
          if (h != 1)
            t[unit(i + 5)] = h ? 1 : -1;
          if (int z = moments[h].q[j][i])
            t[0] = z;
          if (j < vertices[h].size())
            moments[h].q[j + 1][i] =
                moments[h].q[j][i] + p[labels[vertices[h][j] - 1]][i];
        }
      }
  }
  int sample(const Poly &numerator) const {
    V mb{}, mc{};
    const int den = inverse(21), six = inverse(6);
    for (int i = 0; i < 5; i++)
      for (int j = 1; j <= 5; j++) {
        mb[i] = mod(mb[i] + (j + 1) * den % P * moments[1].q[j][i]);
        mc[i] = mod(mc[i] + (j + 1) * den % P * moments[2].q[j][i]);
      }
    array<array<int, 14>, 10> powers{};
    for (int i = 0; i < 10; i++) {
      int j = i % 5, t = mod((i < 5 ? 2 : 3) * six % P * (mb[j] + mc[j]) -
                             (i < 5 ? mb[j] : mc[j]));
      powers[i][0] = 1;
      for (int k = 1; k < 14; k++)
        powers[i][k] = mod(powers[i][k - 1] * t);
    }
    int value = 0;
    for (const auto &[k, exact] : numerator) {
      int v = mpz_fdiv_ui(exact.get_mpz_t(), P);
      for (int i = 0; i < 10; i++)
        v = mod(v * powers[i][(k >> (4 * i)) & 15]);
      value = mod(value + v);
    }
    return value;
  }
  mpz_class integrate(const Poly &numerator) {
    mpz_class result = 0;
    for (const auto &[k, v] : numerator)
      if (v != 0) {
        mpz_class sum = 0;
        for (const auto &term : mean_expansion(k)) {
          const auto &l = moments[1](term.left);
          if (l == 0)
            continue;
          const auto &r = moments[2](term.right);
          if (r == 0)
            continue;
          sum += l * r * term.weight;
        }
        mpz_addmul(result.get_mpz_t(), sum.get_mpz_t(), v.get_mpz_t());
      }
    return result;
  }
  PTensor literal(const Vertex &vertex, int v, const V &polar) const {
    PTensor out;
    for (const auto &w : vertex.words) {
      int n = (vertex.cubic ? w.scale / inv6 : w.scale) *
              (w.external_axis < 0 ? 1 : polar[w.external_axis]);
      if (!n)
        continue;
      const Poly *args[2]{};
      int s = 0;
      for (auto [id, a] : {std::pair{w.first_slot, w.first_axis},
                           std::pair{w.second_slot, w.second_axis}}) {
        if (id < 0)
          continue;
        if (id == 13)
          n *= p[labels[v - 1]][a];
        else {
          if (::ends[id][0] == v)
            n = -n;
          args[s++] = &currents[id][a];
        }
      }
      auto &slot = out[w.element];
      if (s == 0)
        slot[0] += n;
      if (s == 1)
        plus(slot, *args[0], n);
      if (s == 2)
        times(slot, *args[0], *args[1], n);
    }
    for (auto &x : out) {
      clean(x);
      if (vertex.cubic)
        for (auto &[key, coeff] : x) {
          if (!mpz_divisible_ui_p(coeff.get_mpz_t(), 6))
            throw std::logic_error("cubic normalization");
          coeff /= 6;
        }
    }
    return out;
  }
  PMat path(const Flow &flow, const Fields &field, int h) const {
    PMat m;
    for (int i = 0; i < 5; i++)
      m[i * 6][0] = 1;
    for (int j = 0; j < 5; j++) {
      int v = vertices[h][j], label = labels[v - 1];
      PMat out;
      if (label == 9) {
        V eta{};
        eta[2 + field.eta] = 1;
        auto g = literal(*flow.local[v], v, eta);
        for (int i = 0; i < 5; i++)
          for (int r = 0; r < 5; r++)
            if (!m[5 * i + r].empty())
              for (int c = 0; c < 5; c++)
                times(out[5 * i + c], m[5 * i + r], g[5 * r + c]);
      } else {
        const auto &u = field.polar[label], &q = p[label];
        // The primitive sign of the scalar diagonal is fixed by the gamma edge.
        // In a forward transfer the differentiated gamma is on the entering
        // edge.
        bool forward = flow.local[v]->words[0].first_slot == edgeids[h][j];
        int sign = forward ? -1 : 1;
        Poly trace;
        for (int a = 0; a < 5; a++)
          plus(trace, currents[edgeids[h][j]][a], sign * u[a]);
        for (int i = 0; i < 5; i++)
          for (int r = 0; r < 5; r++)
            if (!m[5 * i + r].empty()) {
              times(out[5 * i + r], m[5 * i + r], trace);
              for (int c = 0; c < 5; c++)
                plus(out[5 * i + c], m[5 * i + r],
                     forward ? -u[r] * q[c] : -q[r] * u[c]);
            }
      }
      for (auto &x : out)
        clean(x);
      m = std::move(out);
    }
    return m;
  }
  Poly polynomial(const FlowPlans &plans, int branch, int f,
                  const vector<Probe> &probes) const {
    Poly total;
    for (const auto &flow : plans.by_sites[1 << 9]) {
      auto left = literal(*flow.local[0], 0, {}),
           right = literal(*flow.local[11], 11, {});
      array<PMat, 32> sides, other;
      array<bool, 32> done{}, other_done{};
      // Masks on the two paths are disjoint. Combine the C-path detectors
      // before multiplying B.
      array<PMat, 32> weighted;
      for (const auto &probe : probes)
        if (probe.f == f && probe.weight) {
          int b = probe.mask & 13, c = (probe.mask & 18),
              key = (c & 2) + ((c >> 2) & 4) + 8 * probe.eta;
          if (!done[b]) {
            auto mat = path(flow, Fields(branch, b, f, 0, 1), 1);
            for (int k = 0; k < 5; k++)
              for (int t = 0; t < 5; t++) {
                auto &out = sides[b][5 * k + t];
                for (int i = 0; i < 5; i++)
                  for (int j = 0; j < 5; j++)
                    if (!left[25 * i + 5 * j + k].empty())
                      for (int s = 0; s < 5; s++)
                        if (!right[25 * i + 5 * s + t].empty() &&
                            !mat[5 * j + s].empty()) {
                          Poly pair;
                          times(pair, left[25 * i + 5 * j + k],
                                right[25 * i + 5 * s + t]);
                          clean(pair);
                          times(out, pair, mat[5 * j + s]);
                        }
                clean(out);
              }
            done[b] = true;
          }
          if (!other_done[key]) {
            other[key] = path(flow, Fields(branch, c, f, probe.eta, 1), 2);
            other_done[key] = true;
          }
          for (int n = 0; n < 25; n++)
            plus(weighted[b][n], other[key][n], probe.weight);
        }
      for (int b = 0; b < 32; b++)
        if (done[b])
          for (int n = 0; n < 25; n++) {
            clean(weighted[b][n]);
            times(total, sides[b][n], weighted[b][n], flow.wick);
          }
    }
    clean(total);
    return total;
  }
};
int residual(int branch, int f, const V &sc, int s) {
  const int k = (s & 32) ? (f == 3 ? 2 : f == 1 ? -1 : 1) : 0;
  int v = -(s & 1) + ((s >> 2) & 1) + ((s >> 3) & 1);
  return (branch ? k * (sc[1] - sc[0]) : k * sc[1]) +
         (branch ? -2 : 1) * sc[3] * v;
}
} // namespace
int main(int argc, char **argv) {
  if (argc < 2 || argc > 4)
    throw std::runtime_error(
        "usage: mixed_integrated_order probe.tsv [first-node [last-node]]");
  topology({1, 6, 6});
  FlowPlans plans(1);
  factorials();
  std::ifstream in(argv[1]);
  vector<Probe> rows;
  mpz_class w;
  int m, f, e;
  while (in >> w >> m >> f >> e)
    rows.push_back({m, f, e, w});
  if (rows.size() != 89)
    throw std::runtime_error("89 detector rows required");
  Geometry control(Fields(0, 0, 0, 0, 1), {1, 2, 3, 1, 0}, 63);
  auto raw = control.polynomial(plans, 0, 0, {{0, 0, 0, 1}});
  auto check = control.integrate(raw);
  if (check * mpz_class("25758708393600") !=
      measure_den * mpz_class("4170228116417"))
    throw std::runtime_error("exact scalar integration control");
  mpz_class total = 0;
  int point = 0, id = 0;
  int lo = argc > 2 ? std::stoi(argv[2]) : 0,
      hi = argc > 3 ? std::stoi(argv[3]) : 24;
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 4; k++, id++) {
        if (id < lo || id >= hi)
          continue;
        V sc{1, 2 * i - 1, 2 * j - 2, 2 * k - 3, 0};
        mpz_class value = 0;
        int sample = 0;
        for (int branch = 0; branch < 2; branch++)
          for (int selector : {0, 1, 3})
            for (int switches = 0; switches < 64; switches++) {
              int det = residual(branch, selector, sc, switches);
              if (!det)
                continue;
              Geometry g(Fields(branch, 0, selector, 0, 1), sc, switches);
              auto numerator = g.polynomial(plans, branch, selector, rows);
              int weight = (branch ? 17280 : -1313280) * det;
              if (__builtin_popcount((unsigned)switches) & 1)
                weight = -weight;
              value += weight * g.integrate(numerator);
              sample = mod(sample + 1L * weight * g.sample(numerator));
            }
        int factor = ((6 - i - j - k) & 1 ? -1 : 1) * 2 * bin[2][j] * bin[3][k];
        total += value * factor;
        point = mod(point + sample * factor);
        std::cout << "node " << id << " completed\n" << std::flush;
      }
  mpq_class exact(total, measure_den * 1536);
  exact.canonicalize();
  int sample = mod(point * inverse(mod(1536)));
  if (lo == 0 && hi == 24 &&
      (sample != 554 ||
       mod(mpz_fdiv_ui(exact.get_num_mpz_t(), P) *
           inverse(mpz_fdiv_ui(exact.get_den_mpz_t(), P))) != 735))
    throw std::runtime_error("independent finite-field controls");
  std::cout << "nodes " << lo << " " << hi << " normalized integral " << exact
            << " point residue " << sample << "\n";
}
