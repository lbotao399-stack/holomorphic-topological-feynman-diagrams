// Compare every primitive Wick word under the order-four graph symmetry.
#define SUGRA_PHYSICAL_WORDS_ONLY
#include "../five_theta/check_physical.cpp"
#include <numeric>

using Map = std::map<std::array<int, 6>, int>;
using Edges = std::array<int, 13>;
std::vector<int> incident(int v) {
  std::vector<int> result;
  for (int e = 0; e < 13; e++)
    if (ends[e][0] == v || ends[e][1] == v)
      result.push_back(e);
  return result;
}
Map relabel(const Vertex &vertex, int old_v, int new_v, const Edges &edges) {
  Map words;
  auto old = incident(old_v), next = incident(new_v);
  for (auto w : vertex.words) {
    int pos = w.element, idx = 0;
    std::array<int, 3> axis{};
    for (int k = int(old.size()) - 1; k >= 0; k--) {
      int n = std::find(next.begin(), next.end(), edges[old[k]]) - next.begin();
      if (n >= next.size())
        throw std::logic_error("incidence");
      axis[n] = pos % 5;
      pos /= 5;
    }
    for (size_t k = 0; k < old.size(); k++)
      idx = 5 * idx + axis[k];
    auto slot = [&](int e) { return e < 0 || e == 13 ? e : edges[e]; };
    int a = slot(w.first_slot), i = w.first_axis, b = slot(w.second_slot),
        j = w.second_axis;
    if (b >= 0 && std::pair{b, j} < std::pair{a, i}) {
      std::swap(a, b);
      std::swap(i, j);
    }
    words[{idx, w.external_axis, a, i, b, j}] += w.scale;
  }
  for (auto it = words.begin(); it != words.end();)
    if (!it->second)
      it = words.erase(it);
    else
      ++it;
  return words;
}
int compare(Map a, const Map &b) {
  if (a == b)
    return 1;
  for (auto &[k, v] : a)
    v = -v;
  return a == b ? -1 : 0;
}
int small_factorial(int n) {
  int v = 1;
  for (int i = 1; i <= n; i++)
    v *= i;
  return v;
}
// Exact conditional integral for each set of zero external momenta. The
// active insertion times are held fixed, so the other path stays arbitrary.
int zero_waves(const FlowPlans &plans) {
  int checks = 0;
  for (int marked : {1, 2})
    for (int h : {1, 2})
      for (int mask = 1; mask < 32; mask++) {
        if (h == marked && (mask & 1))
          continue;
        int k = __builtin_popcount(unsigned(mask)), d = 5 - k;
        std::array<int, 6> gap{};
        for (int j = 0; j <= d; j++)
          gap[j] = 1 + j % 3;
        auto value = [](int label, int gap) {
          return 1 + (2 * label + 3 * gap + label * gap) % 5;
        };
        std::array<int, 5> base{0, 1, 2, 3, 4};
        std::stable_partition(base.begin(), base.end(),
                              [&](int j) { return !(mask & (1 << j)); });
        auto signed_flow = [&](std::array<int, 5> order, int flow_id) {
          int pos = h == marked ? int(std::find(order.begin(), order.end(), 0) -
                                      order.begin())
                                : 4;
          int gamma = vertices[marked][pos] - 1;
          auto &flow = plans.by_sites[1 << gamma][flow_id];
          int sign = flow.wick * ((9 - gamma) & 1 ? -1 : 1);
          for (int j = 0; j < 5; j++)
            if (mask & (1 << order[j])) {
              auto &words = flow.local[vertices[h][j]]->words;
              sign *= words[0].first_slot == edgeids[h][j] ? -1 : 1;
            }
          return sign;
        };
        for (int flow_id = 0; flow_id < 4; flow_id++) {
          int reference = signed_flow(base, flow_id);
          long long right = reference * small_factorial(k);
          for (int label = 0; label < 5; label++)
            if (mask & (1 << label)) {
              int sum = 0;
              for (int j = 0; j <= d; j++)
                sum += gap[j] * value(label, j);
              right *= sum;
            }
          long long left = 0;
          std::array<int, 5> order{0, 1, 2, 3, 4};
          do {
            int last = -1, index = 0;
            std::array<int, 6> count{};
            bool sorted = true;
            long long integrand = 1;
            for (int label : order) {
              if (mask & (1 << label)) {
                integrand *= gap[index] * value(label, index);
                count[index]++;
              } else {
                if (label < last)
                  sorted = false;
                last = label;
                index++;
              }
            }
            if (!sorted)
              continue;
            int sign = signed_flow(order, flow_id);
            if (sign != reference)
              throw std::logic_error(
                  "constant wave crossed gamma with wrong sign");
            int den = 1;
            for (int n : count)
              den *= small_factorial(n);
            left += sign * integrand * (small_factorial(k) / den);
          } while (std::next_permutation(order.begin(), order.end()));
          if (left != right)
            throw std::logic_error("conditional zero-wave moment");
          checks++;
        }
      }
  return checks;
}
int main() {
  topology({1, 6, 6});
  FlowPlans plans(1);
  Edges identity;
  std::iota(identity.begin(), identity.end(), 0);
  int matched = 0;
  std::cout << zero_waves(plans) << " conditional constant-wave integrations\n";
  for (bool swap : {false, true})
    for (bool reverse : {false, true}) {
      Edges es{};
      std::array<int, 12> vert{};
      vert[0] = reverse ? 11 : 0;
      vert[11] = reverse ? 0 : 11;
      for (int h = 0; h < 3; h++) {
        auto dest = h && swap ? 3 - h : h;
        int n = edgeids[h].size();
        for (int j = 0; j < n; j++)
          es[edgeids[h][j]] = edgeids[dest][reverse ? n - 1 - j : j];
        for (int j = 0; j + 1 < n; j++)
          vert[vertices[h][j]] = vertices[dest][reverse ? n - 2 - j : j];
      }
      const auto &old = plans.by_sites[1 << 9],
                 &next = plans.by_sites[1 << (vert[10] - 1)];
      std::set<int> seen;
      for (auto &f : old) {
        int ways = 0;
        for (int j = 0; j < 4; j++) {
          const auto &candidate = next[j];
          int s = f.wick * candidate.wick;
          for (int v = 0; v < 12 && s; v++)
            s *= compare(
                relabel(*f.local[v], v, vert[v], es),
                relabel(*candidate.local[vert[v]], vert[v], vert[v], identity));
          if (s) {
            if (s != 1)
              throw std::logic_error("orientation sign");
            seen.insert(j);
            ways++;
          }
        }
        if (ways != 1)
          throw std::logic_error("Wick matching");
        matched += ways;
      }
      if (seen.size() != 4)
        throw std::logic_error("flow permutation");
      int external = (9 - (vert[10] - 1)) % 2 ? -1 : 1;
      int determinant = swap ? -1 : 1;
      if (external * determinant != 1)
        throw std::logic_error("external and determinant signs");
      std::cout << "swap " << swap << " reverse " << reverse << " gamma slot "
                << vert[10] - 1 << " determinant " << determinant
                << " external " << external << "\n";
    }
  std::cout << matched << " global primitive Wick-word identities\n";
}
