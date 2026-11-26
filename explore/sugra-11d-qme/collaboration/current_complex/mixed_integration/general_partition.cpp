// Three-path kernel. Contract every order of the longest path into its
// normalized resistance moments before joining the other two path tensors.
#define TWIST_ORDERED_DETECTOR_LIBRARY
#include "ordered_detector.cpp"

namespace complete {
using namespace ordered;
// Both dense dot products have at most 120 terms. For a degree-eight
// angular monomial the number of triple assignments is at most 3^8.
static_assert(32768ULL * (prime - 1) * (prime - 1) + prime < (1ULL << 63));
using Residue = native::Entry;
using Product = std::conditional_t<(prime <= 32768), int, int64_t>;
struct SplitBatch {};
using Rows = std::array<std::map<int, std::vector<int>>, 25>;
using Joint = std::map<std::pair<int, int>, Rows>;
struct Driver : Job {
  std::array<int, 3> requested;
  static int gamma_sign(const native::Waves &w, int order) {
    for (int i = 0; i < w.active; ++i) {
      int j = w.orders[order][i];
      if (w.labels[j] == 9)
        return ((i + j) & 1) ? -1 : 1;
    }
    return 1;
  }
  int reversal = 1;
  void all_orders() {
    for (int h = 0; h < 3; ++h) {
      auto &w = waves[h];
      int count = w.active;
      // The native path enumerator assumes at most five active waves and
      // uses the six-gap weight. Restore the actual number of path vertices.
      for (int i = 0; i < count; ++i)
        if (w.p[i] == V{})
          throw std::runtime_error("all-active kernel");
      for (int i = count; i < 5; ++i)
        w.p[i] = V{};
      w.orders.clear();
      w.moment.clear();
      w.make_moments(requested[h] >= fact(count) ? 120 : requested[h]);
      if (w.active != count)
        throw std::logic_error("path size");
      int scale = red(fact(count) * native::inverse_prime(fact(5)));
      for (double &v : w.moment)
        v = red(int64_t(v) * scale);
    }
    bool full = true;
    for (int h = 0; h < 3; ++h)
      full &= int(waves[h].orders.size()) == fact(waves[h].active);
    if (full) {
      int h = waves[0].active >= 2 ? 0 : 1;
      auto &path = waves[h];
      if (path.active < 2)
        throw std::logic_error("missing reversible path");
      int kept = 0;
      for (int q = 0; q < int(path.orders.size()); ++q) {
        if (path.orders[q][0] > path.orders[q][path.active - 1])
          continue;
        if (kept != q) {
          path.orders[kept] = path.orders[q];
          std::copy_n(path.moment.data() + q * native::mean_count,
                      native::mean_count,
                      path.moment.data() + kept * native::mean_count);
        }
        ++kept;
      }
      if (2 * kept != int(path.orders.size()))
        throw std::logic_error("reversal orbit");
      path.orders.resize(kept);
      path.moment.resize(kept * native::mean_count);
      reversal = 2;
    }
  }
  int third_kind(int flow) const {
    if (marked == 2)
      return int(flow >= 2);
    if (flow < 2)
      return flow;
    return (2 == (marked + 1) % 3) == (flow == 2);
  }
  Joint third(int kind) {
    auto &w = waves[2];
    Joint combined;
    std::vector<Rows> basis(w.polar.size());
    for (int u = 0; u < int(w.polar.size()); ++u) {
      std::map<std::pair<int, int>, std::vector<double>> groups;
      int m = native::mean_count, n = w.orders.size();
      for (int i = 0; i < n; ++i) {
        int sign = gamma_sign(w, i);
        for (auto &[exp, matrix] : w.chain(i, u, kind))
          for (int j = 0; j < 25; ++j)
            if (matrix[j]) {
              auto &terms = groups[{j, exp}];
              if (terms.empty())
                terms.resize(n);
              terms[i] = red(sign * matrix[j]);
            }
      }
      if (groups.empty())
        continue;
      std::vector<double> vals, averages(groups.size() * m);
      for (auto &[key, r] : groups)
        vals.insert(vals.end(), r.begin(), r.end());
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, groups.size(), m,
                  n, 1, vals.data(), n, w.moment.data(), m, 0, averages.data(),
                  m);
      size_t pos = 0;
      for (auto &[key, r] : groups) {
        auto &out = basis[u][key.first][key.second];
        out.resize(m);
        for (int j = 0; j < m; ++j)
          out[j] = uint64_t(averages[pos++]) % prime;
      }
    }
    for (auto [ia, ib, ic, weight] : probes) {
      auto &dst = combined[{ia, ib}];
      for (int j = 0; j < 25; ++j)
        for (auto &[key, values] : basis[ic][j]) {
          auto &row = dst[j][key];
          if (row.empty())
            row.resize(values.size());
          for (size_t k = 0; k < values.size(); ++k)
            row[k] = red(row[k] + int64_t(weight) * values[k]);
        }
    }
    return combined;
  }
  struct CTerm {
    int exponent, degree, factor, vector;
  };
  struct Lane {
    int ordinal = 0;
    std::vector<CTerm> terms;
  };
  struct Contribution {
    int vector, row, factor;
  };
  struct RowIndex {
    // The exponent tables are sorted by degree. After fixing the first
    // five-variable index, all allowed second indices form an initial segment.
    static int left_offset(int exponent) {
      static const auto offsets = [] {
        std::array<int, 1287> value{};
        int position = 0;
        for (size_t i = 0; i < value.size(); ++i) {
          value[i] = position;
          int remaining = 8 - degree(indices.code[i]);
          for (size_t j = 0; j < value.size() && degree(indices.code[j]) <= remaining; ++j)
            ++position;
        }
        if (position != 43758)
          throw std::logic_error("joint degree-eight dimension");
        return value;
      }();
      return offsets[indices.id[exponent]];
    }

    std::vector<Key> keys;
    std::vector<Residue> source;
    std::vector<int> position = std::vector<int>(43758);
    std::vector<std::vector<Contribution>> groups;
    size_t records = 0;
    void add(Key key, int monomial, Contribution c) {
      if (++records > 6000000)
        throw SplitBatch{};
      int &slot = position[monomial];
      if (!slot) {
        keys.push_back(key);
        groups.emplace_back();
        slot = keys.size();
      }
      groups[slot - 1].push_back(c);
    }
  };
  int batch() {
    all_orders();
    make_angular();
    using Paths = std::map<std::pair<int, int>, std::array<Lane, 25>>;
    std::array<Paths, 4> lanes;
    std::vector<std::vector<Residue>> dictionary;
    std::unordered_map<uint64_t, std::vector<int>> buckets;
    int lane_count = 0;
    std::array<int, 2> representative{-1, -1};
    for (int t = 0; t < 4; ++t) {
      int kind = third_kind(t);
      int previous = representative[kind];
      if (previous != -1) {
        lanes[t] = lanes[previous];
        for (auto &[ab, cells] : lanes[t])
          for (auto &lane : cells) {
            lane.ordinal = lane_count++;
            if (signs[t] != signs[previous])
              for (auto &term : lane.terms)
                term.factor = red(-term.factor);
          }
        continue;
      }
      representative[kind] = t;
      auto original = third(kind);
      for (auto &[ab, cells] : original) {
        auto &l = lanes[t][ab];
        for (auto &lane : l)
          lane.ordinal = lane_count++;
        for (int cell = 0; cell < 25; ++cell)
          for (auto &[k, v] : cells[cell]) {
            auto first =
                std::find_if(v.begin(), v.end(), [](int x) { return bool(x); });
            if (first == v.end())
              continue;
            int scale = *first, inv = native::inverse_prime(scale);
            uint64_t hash = 0;
            std::vector<Residue> normalized(v.size());
            for (int i = 0; i < int(v.size()); ++i) {
              normalized[i] = red(int64_t(v[i]) * inv);
              hash = hash * 1099511628211ULL + normalized[i];
            }
            auto &candidates = buckets[hash];
            int id = -1;
            for (int c : candidates)
              if (dictionary[c] == normalized) {
                id = c;
                break;
              }
            if (id < 0) {
              id = dictionary.size();
              dictionary.push_back(std::move(normalized));
              candidates.push_back(id);
            }
            l[cell].terms.push_back({k, degree(k), red(signs[t] * scale), id});
          }
      }
    }
    int result = 0;
    const int n = waves[1].orders.size(), m = native::mean_count;
    std::vector<int> bmoment(m * n);
    for (int i = 0; i < m; ++i)
      for (int o = 0; o < n; ++o)
        bmoment[i * n + o] = int(waves[1].moment[o * m + i]);
    auto integrate_block = [&](int ao, int begin, int end_lane) {
      RowIndex index;
      for (int t = 0; t < 4; ++t) {
        const auto &end = endpoints[t][1];
        using Link = std::pair<int, Poly>;
        // Contract both endpoint tensors with path A once per A order.
        // A link runs from a path-B matrix entry to a path-C matrix entry.
        std::vector<std::array<std::vector<Link>, 25>> links(waves[0].polar.size());
        for (int u = 0; u < int(waves[0].polar.size()); ++u) {
          auto left = apply(endpoints[t][0], fiber(0, u, t, ao), 0);
          for (int k = 0; k < 5; ++k)
            for (int l = 0; l < 5; ++l)
              for (int b = 0; b < 5; ++b)
                for (int j = 0; j < 5; ++j) {
                  Poly link;
                  for (int a = 0; a < 5; ++a)
                    times(link, left[25 * a + 5 * b + k], end[25 * a + 5 * j + l]);
                  trim(link);
                  if (!link.empty())
                    links[u][5 * k + l].emplace_back(5 * b + j, std::move(link));
                }
        }
        std::vector<std::vector<Fiber>> second(n);
        for (int bo = 0; bo < n; ++bo)
          for (int u = 0; u < int(waves[1].polar.size()); ++u)
            second[bo].push_back(fiber(1, u, t, bo));
        std::vector<int> source_row(43758, -1), base_ids;
        std::vector<Key> base_keys;
        std::vector<int64_t> side(43758, 0);
        std::vector<std::pair<int, Key>> touched;
        // For fixed C exponent, all allowed B exponents form a prefix:
        // the number of five-variable monomials of degree <= d is C(d+5,5).
        constexpr std::array<int, 9> limits{1, 6, 21, 56, 126, 252, 462, 792, 1287};
        for (const auto &[ab, group] : lanes[t]) {
          for (int cell = 0; cell < 25; ++cell) {
            const auto &lane = group[cell];
            if (lane.ordinal < begin || lane.ordinal >= end_lane || lane.terms.empty())
              continue;
            base_ids.clear();
            base_keys.clear();
            for (int bo = 0; bo < n; ++bo) {
              touched.clear();
              for (const auto &[slot, link] : links[ab.first][cell]) {
                const auto &matrix = second[bo][ab.second][slot];
                if (matrix.empty())
                  continue;
                for (auto [key, factor] : link) {
                  int left = key & mask, right = key >> 20;
                  int dr = degree(right);
                  if (dr > 8)
                    continue;
                  int offset = RowIndex::left_offset(right);
                  for (auto [power, coeff] : matrix) {
                    if (!coeff)
                      continue;
                    int i = indices.id[left + power];
                    if (i >= limits[8 - dr])
                      continue;
                    int pos = offset + i;
                    if (!side[pos])
                      touched.emplace_back(pos, key + power);
                    // At fixed C exponent a specified B monomial has at most
                    // 108 factorizations of degree <= 8, in each of 25 slots.
                    side[pos] += int64_t(factor) * coeff;
                  }
                }
              }
              int sigma = gamma_sign(waves[0], ao) * gamma_sign(waves[1], bo);
              for (auto [pos, key] : touched) {
                int val = red(sigma * side[pos]);
                side[pos] = 0;
                if (!val)
                  continue;
                int &row = source_row[pos];
                if (row < 0) {
                  if (index.source.size() + n > 16000000)
                    throw SplitBatch{};
                  row = index.source.size();
                  index.source.resize(row + n, 0);
                  base_ids.push_back(pos);
                  base_keys.push_back(key);
                }
                index.source[row + bo] = val;
              }
            }
            for (size_t i = 0; i < base_ids.size(); ++i) {
              Key base = base_keys[i];
              int pos = base_ids[i], row = source_row[pos];
              source_row[pos] = -1;
              int room = 8 - degree(int(base & mask)) - degree(int(base >> 20));
              int y = base >> 20, origin = RowIndex::left_offset(base & mask);
              for (const auto &term : lane.terms) {
                if (term.degree > room)
                  continue;
                Key key = base + (Key(term.exponent) << 20);
                int monomial = origin + indices.id[y + term.exponent];
                index.add(key, monomial, {term.vector, row, term.factor});
              }
            }
          }
        }
      }
      // Each common kernel acts on all B orders in one block. Only subpowers
      // of the combined moment occur on path C.
      std::array<int, 1287> local{};
      std::vector<int64_t> accumulated;
      std::vector<double> right, left, mixed;
      const double *ma = waves[0].moment.data() + ao * m;
      std::vector<int64_t> coefficients(dictionary.size() * n);
      std::vector<int> last_monomial(dictionary.size(), -1), patterns;
      for (int pos = 0; pos < int(index.keys.size()); ++pos) {
        Key k = index.keys[pos];
        int x = k & mask, y = k >> 20, dx = degree(x), dy = degree(y);
        const auto &sub = splits(x + y);
        int count = sub.size();
        accumulated.assign(count * n, 0);
        for (int i = 0; i < count; ++i)
          local[indices.id[sub[i].key]] = i;
        for (const auto &s : parts(x))
          for (const auto &t : parts(y)) {
            int av = int(ma[indices.id[s.keys[0] + t.keys[0]]]);
            if (!av)
              continue;
            int a = angular[slot(dx, dy, s.bown, t.cown)];
            int factor = red(int64_t(av) * a * s.coefficient * t.coefficient);
            if (!factor)
              continue;
            const int *bm =
                bmoment.data() + n * indices.id[s.keys[1] + t.keys[1]];
            auto *target = accumulated.data() +
                           n * local[indices.id[s.keys[2] + t.keys[2]]];
            for (int b = 0; b < n; ++b)
              target[b] += Product(factor) * bm[b];
          }
        right.resize(count * n);
        for (int i = 0; i < count * n; ++i)
          right[i] = accumulated[i] % prime;
        patterns.clear();
        const auto &records = index.groups[pos];
        // Fold repeated third-path vectors one monomial at a time. At the
        // largest supported prime each chunk's 32768 products fit in int64.
        for (size_t start = 0; start < records.size(); start += 32768) {
          size_t stop = std::min(start + 32768, records.size());
          for (size_t j = start; j < stop; ++j) {
            const auto &entry = records[j];
            auto *dest = coefficients.data() + size_t(entry.vector) * n;
            if (last_monomial[entry.vector] != pos) {
              last_monomial[entry.vector] = pos;
              patterns.push_back(entry.vector);
              std::fill_n(dest, n, 0);
            }
            const auto *weights = index.source.data() + entry.row;
            for (int b = 0; b < n; ++b)
              dest[b] += int64_t(entry.factor) * weights[b];
          }
          for (int v : patterns) {
            auto *dest = coefficients.data() + size_t(v) * n;
            for (int b = 0; b < n; ++b)
              dest[b] %= prime;
          }
        }
        int active = 0;
        for (int v : patterns) {
          const auto *c = coefficients.data() + size_t(v) * n;
          if (std::any_of(c, c + n, [](int64_t x) { return x != 0; }))
            patterns[active++] = v;
        }
        if (!active)
          continue;
        left.resize(active * count);
        mixed.resize(active * n);
        for (int r = 0; r < active; ++r) {
          const auto &v = dictionary[patterns[r]];
          for (int i = 0; i < count; ++i)
            left[r * count + i] = v[indices.id[sub[i].key]];
        }
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, active, n, count,
                    1., left.data(), count, right.data(), n, 0., mixed.data(),
                    n);
        int64_t total = 0;
        for (int r = 0; r < active; ++r) {
          // Reduce each row for large primes. Its 120 products already fit
          // below 2^53 by the dense-multiplication bound.
          if constexpr (prime > 32768)
            total %= prime;
          const auto *c = coefficients.data() + size_t(patterns[r]) * n;
          for (int b = 0; b < n; ++b)
            total += int64_t(uint64_t(mixed[r * n + b]) % prime) * c[b];
        }
        result = red(result + total);
      }
    };
    for (int ao = 0; ao < int(waves[0].orders.size()); ++ao) {
      std::vector<std::pair<int, int>> pending{{0, lane_count}};
      while (!pending.empty()) {
        auto [begin, end] = pending.back();
        pending.pop_back();
        try {
          integrate_block(ao, begin, end);
        } catch (const SplitBatch &) {
          // A block splits before it contributes to the answer.
          int mid = (begin + end) / 2;
          if (mid == begin)
            throw std::length_error("single fiber block");
          pending.push_back({mid, end});
          pending.push_back({begin, mid});
        }
      }
    }
    return red(int64_t(det) * reversal * result);
  }
};
} // namespace complete
int main(int argc, char **argv) {
  if (argc == 2 && std::string(argv[1]) == "--prime") {
    std::cout << native::prime << '\n';
    return 0;
  }
  std::array<int, 3> limits{120, 120, 120};
  if (argc == 4) {
    for (int i = 0; i < 3; ++i)
      limits[i] = std::stoi(argv[i + 1]);
  } else if (argc != 1)
    return 2;
  for (int v : limits)
    if (v < 1 || v > 120)
      return 2;
  std::ios::sync_with_stdio(false);
  std::cin.tie(nullptr);
  for (;;) {
    complete::Driver j;
    if (!j.read())
      break;
    j.requested = limits;
    std::cout << j.id << ' ' << j.batch() << '\n' << std::flush;
  }
}
