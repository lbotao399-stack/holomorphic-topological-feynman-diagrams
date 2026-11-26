// Evaluate one ordered, four-orientation theta integrand on the complete cycle.
#define SUGRA_PHYSICAL_WORDS_ONLY
#include "../five_theta/check_physical.cpp"
#include <gmpxx.h>
#include <numeric>

using std::array;
using std::vector;
struct Jet {
  array<int, 12> c{};
  Jet(int n = 0) { c[0] = mod(n); }
  bool zero() const {
    return std::all_of(c.begin(), c.end(), [](int n) { return n == 0; });
  }
};
Jet scale(const Jet &x, int s) {
  Jet z;
  for (int i = 0; i < 12; ++i)
    z.c[i] = mod(1LL * x.c[i] * s);
  return z;
}
void add(Jet &x, const Jet &y, int s = 1) {
  for (int i = 0; i < 12; ++i)
    x.c[i] = mod(x.c[i] + 1LL * s * y.c[i]);
}
Jet product(const Jet &x, const Jet &y) {
  Jet z;
  for (int i = 0; i < 12; ++i)
    if (x.c[i])
      for (int j = 0; j < 12; ++j)
        if (y.c[j] && (i & j & 3) == 0 && i / 4 + j / 4 < 3)
          z.c[i + j] += x.c[i] * y.c[j];
  for (int &v : z.c)
    v %= P;
  return z;
}
int inverse(int v) {
  for (int i = 1; i < P; ++i)
    if (i * v % P == 1)
      return i;
  throw std::runtime_error("nonunit denominator");
}
using Cov = array<Jet, 5>;
using JTensor = array<Jet, 125>;
using JMatrix = array<Jet, 25>;
using Momentum = array<Cov, 10>;
int determinant(array<array<int, 10>, 10> rows) {
  int sign = 1;
  for (int i = 0; i < 10; ++i) {
    int j = i;
    while (j < 10 && rows[j][i] == 0)
      ++j;
    if (j == 10)
      return 0;
    if (j != i) {
      std::swap(rows[i], rows[j]);
      sign = -sign;
    }
    int pivot = rows[i][i];
    sign = mod(sign * pivot);
    int u = inverse(pivot);
    for (int k = i + 1; k < 10; ++k)
      if (rows[k][i]) {
        int v = mod(u * rows[k][i]);
        for (int h = i; h < 10; ++h)
          rows[k][h] = mod(rows[k][h] - v * rows[i][h]);
      }
  }
  return mod(sign);
}

struct Fields {
  array<V, 9> direction{}, polar{};
  V function{};
  int f, eta, mask;
  mpz_class coefficient;
  Fields(int branch, int bits, int selector, int covector,
         const mpz_class &weight)
      : f(selector), eta(covector), mask(bits), coefficient(weight) {
    direction[0][0] = direction[1][1] = 1;
    direction[2][0] = 1;
    direction[2][1] = branch;
    direction[3][0] = branch;
    direction[3][1] = branch ? -1 : 1;
    for (int j = 0; j < 4; ++j) {
      polar[j][0] = -direction[j][1];
      polar[j][1] = direction[j][0];
    }
    for (int j = 0; j < 5; ++j) {
      int axis = 2 + j / 2, other = 2 + (j / 2 + 1 + ((bits >> j) & 1)) % 3;
      direction[j + 4][axis] = 1;
      if (j == 1 || j == 3)
        direction[j + 4][axis + 1] = 1;
      polar[j + 4][other] = 1;
      polar[j + 4][axis] = -direction[j + 4][other];
    }
    function[2 + f % 3] = (f == 3 ? 2 : 1);
  }
  Momentum momenta(int switches) const {
    Momentum p{};
    for (int i = 0; i < 9; ++i) {
      Jet s(1);
      if (i < 3) {
        s = Jet();
        s.c[1 << i] = 1;
      }
      if (i >= 4)
        s = Jet((switches >> (i - 4)) & 1);
      for (int j = 0; j < 5; ++j)
        p[i][j] = scale(s, direction[i][j]);
    }
    if (switches & 32)
      for (int j = 0; j < 5; ++j)
        add(p[3][j], Jet(function[j]));
    for (int i = 0; i < 9; ++i)
      for (int j = 0; j < 5; ++j)
        add(p[9][j], p[i][j], -1);
    return p;
  }
};

struct Ordered {
  array<int, 3> shape;
  array<int, 10> order;
  array<int, 10> paths{};
  const vector<Flow> &flows;
  int gamma_slot, sign;
  Ordered(array<int, 3> lens, array<int, 10> ids, const FlowPlans &source)
      : shape(lens), order(ids),
        flows(source.by_sites[1 << (std::find(ids.begin(), ids.end(), 9) -
                                    ids.begin())]) {
    gamma_slot = std::find(ids.begin(), ids.end(), 9) - ids.begin();
    sign = (9 - gamma_slot) % 2 ? -1 : 1;
    for (int k = 0; k < 3; ++k)
      for (int v : vertices[k])
        paths[order[v - 1]] = k;
  }
  Jet det_residual(const Fields &field, const Momentum &p, int switches) const {
    // The eight pure label columns have had their scalar scales removed.
    Jet result;
    for (int d : {0, 1, 2, 4}) {
      array<array<int, 10>, 10> a{};
      for (int i = 0; i < 10; ++i)
        for (int k = 0; k < 5; ++k) {
          int v = i == 9 ? p[9][k].c[d] : field.direction[i][k];
          if (i == 3 && (switches & 32))
            v += field.function[k];
          if (paths[i] != 2)
            a[k][i] = mod(v);
          if (paths[i] != 0)
            a[k + 5][i] = mod(v);
        }
      result.c[d] = determinant(a);
    }
    return result;
  }
  array<Cov, 13> currents(const Momentum &p) const {
    // Ordered edge resistance fractions proportional to 1,2,... within
    // each path. The conductance fractions are (1,2,3)/6.
    array<Cov, 3> tail{};
    Cov all{};
    for (int k = 0; k < 3; ++k) {
      int denom = shape[k] * (shape[k] + 1) / 2;
      int t = denom;
      for (int j = 0; j < shape[k] - 1; ++j) {
        t -= j + 1;
        const auto &cov = p[order[vertices[k][j] - 1]];
        for (int l = 0; l < 5; ++l)
          add(tail[k][l], cov[l], mod(t * inverse(denom)));
      }
      for (int l = 0; l < 5; ++l)
        add(all[l], tail[k][l]);
    }
    array<Cov, 13> ys{};
    for (int k = 0; k < 3; ++k) {
      Cov y;
      for (int l = 0; l < 5; ++l) {
        y[l] = scale(all[l], mod((k + 1) * inverse(6)));
        add(y[l], tail[k][l], -1);
      }
      for (int j = 0; j < shape[k]; ++j) {
        ys[edgeids[k][j]] = y;
        if (j < shape[k] - 1)
          for (int l = 0; l < 5; ++l)
            add(y[l], p[order[vertices[k][j] - 1]][l]);
      }
    }
    return ys;
  }
};

JMatrix chain_jet(const array<JTensor, 12> &tv, int path) {
  JMatrix mat{};
  for (int i = 0; i < 5; ++i)
    mat[5 * i + i] = 1;
  for (int v : vertices[path]) {
    JMatrix next{};
    for (int i = 0; i < 5; ++i)
      for (int k = 0; k < 5; ++k)
        if (!mat[5 * i + k].zero())
          for (int j = 0; j < 5; ++j)
            if (!tv[v][5 * k + j].zero())
              add(next[5 * i + j], product(mat[5 * i + k], tv[v][5 * k + j]));
    mat = next;
  }
  return mat;
}
Jet wick(const Ordered &order, const Fields &field, const Momentum &p,
         const array<Cov, 13> &ys) {
  array<array<Cov, 14>, 12> signed_y{};
  array<V, 12> pol{};
  for (int v = 0; v < 12; ++v) {
    if (v >= 1 && v <= 10) {
      int i = order.order[v - 1];
      signed_y[v][13] = p[i];
      if (i == 9)
        pol[v][2 + field.eta] = 1;
      else
        pol[v] = field.polar[i];
    }
    for (int e = 0; e < 13; ++e)
      if (ends[e][0] == v || ends[e][1] == v) {
        int s = ends[e][1] == v ? 1 : -1;
        for (int k = 0; k < 5; ++k)
          signed_y[v][e][k] = scale(ys[e][k], s);
      }
  }
  Jet total;
  for (const auto &flow : order.flows) {
    array<JTensor, 12> tv{};
    for (int v = 0; v < 12; ++v)
      for (const auto &w : flow.local[v]->words) {
        int c = w.scale;
        if (w.external_axis >= 0)
          c *= pol[v][w.external_axis];
        if (!c)
          continue;
        auto t = signed_y[v][w.first_slot][w.first_axis];
        if (w.second_slot >= 0)
          t = product(t, signed_y[v][w.second_slot][w.second_axis]);
        add(tv[v][w.element], t, c);
      }
    JTensor acc = tv[0];
    for (int k = 0; k < 3; ++k) {
      auto mat = chain_jet(tv, k);
      JTensor next{};
      int stride = k == 0 ? 25 : k == 1 ? 5 : 1;
      for (int i = 0; i < 125; ++i)
        if (!acc[i].zero()) {
          int index = i / stride % 5;
          for (int j = 0; j < 5; ++j)
            if (!mat[5 * index + j].zero())
              add(next[i + (j - index) * stride],
                  product(acc[i], mat[5 * index + j]));
        }
      acc = next;
    }
    for (int i = 0; i < 125; ++i)
      if (!acc[i].zero() && !tv[11][i].zero())
        add(total, product(acc[i], tv[11][i]), flow.wick);
  }
  return total;
}
// Independent extraction of x1 x2 x3^2 from a polynomial of planar degree <=7
// after the three pure planar momenta have been factored out. Interpolate a
// scalar T, then polarize its quartic term at (iT,jT,kT).
void check_extraction(const Ordered &graph, const Fields &fields) {
  int switches = 63;
  auto jets = fields.momenta(switches);
  auto determinant = graph.det_residual(fields, jets, switches);
  int direct =
      product(determinant, wick(graph, fields, jets, graph.currents(jets)))
          .c[11];
  int total = 0;
  for (int t = 0; t <= 7; ++t) {
    array<int, 8> lagrange{};
    lagrange[0] = 1;
    int degree = 0, den = 1;
    for (int z = 0; z <= 7; ++z)
      if (t != z) {
        ++degree;
        for (int h = degree; h >= 0; --h)
          lagrange[h] = mod((h ? lagrange[h - 1] : 0) - z * lagrange[h]);
        den = mod(den * (t - z));
      }
    int coefficient = mod(lagrange[4] * inverse(den));
    for (int i = 0; i <= 1; ++i)
      for (int j = 0; j <= 1; ++j)
        for (int k = 0; k <= 2; ++k) {
          auto numeric = jets;
          for (auto &p : numeric)
            for (auto &e : p)
              e = Jet(e.c[0] + t * (i * e.c[1] + j * e.c[2] + k * e.c[4]));
          auto d = graph.det_residual(fields, numeric, switches).c[0];
          auto n = wick(graph, fields, numeric, graph.currents(numeric)).c[0];
          int sign = (i + j + k) % 2 ? -1 : 1;
          total = mod(total +
                      1LL * sign * (k == 1 ? 2 : 1) * coefficient * d % P * n);
        }
  }
  if (mod(total * inverse(2)) != direct)
    throw std::runtime_error("independent jet interpolation failed");
}
int main(int argc, char **argv) {
  if (argc != 6) {
    std::cerr << "usage: mixed_graph_point p q r ten-label-order probe.tsv\n";
    return 2;
  }
  array<int, 3> shape = {std::stoi(argv[1]), std::stoi(argv[2]),
                         std::stoi(argv[3])};
  array<int, 10> order{};
  std::string word = argv[4];
  if (word.size() != 10)
    throw std::runtime_error("ten vertex labels required");
  for (int i = 0; i < 10; ++i)
    order[i] = word[i] - '0';
  auto sorted = order;
  std::sort(sorted.begin(), sorted.end());
  for (int i = 0; i < 10; ++i)
    if (sorted[i] != i)
      throw std::runtime_error("invalid vertex assignment");
  topology(shape);
  FlowPlans physical(1);
  Ordered graph(shape, order, physical);
  std::ifstream tsv(argv[5]);
  int mask, f, eta, count = 0, result = 0;
  mpz_class weight;
  while (tsv >> weight >> mask >> f >> eta) {
    for (int branch = 0; branch < 2; ++branch) {
      Fields fields(branch, mask, f, eta, weight);
      if (count == 0)
        check_extraction(graph, fields);
      int projected = 0;
      for (int switches = 0; switches < 64; ++switches) {
        auto p = fields.momenta(switches);
        auto d = graph.det_residual(fields, p, switches);
        if (d.zero())
          continue;
        auto raw = wick(graph, fields, p, graph.currents(p));
        int sign = __builtin_popcount(unsigned(switches)) % 2 ? -1 : 1;
        projected += sign * product(d, raw).c[11];
      }
      int scale = branch ? mod(5 * inverse(6)) : mod(-190 * inverse(3));
      result = mod(result + 1LL * mod(projected) * scale % P *
                                mpz_fdiv_ui(weight.get_mpz_t(), P));
    }
    ++count;
  }
  if (!count || !tsv.eof())
    throw std::runtime_error("invalid weights");
  result = mod(1LL * graph.sign * 20736 * result);
  std::cout << "projected ordered four-flow integrand modulo " << P << ": "
            << result << " (terms " << count << ")\n";
  return 0;
}
