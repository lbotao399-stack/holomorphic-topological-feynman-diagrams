// Recover each of the four primitive Wick contractions using two SL5
// transports.
#define SUGRA_PHYSICAL_WORDS_ONLY
#include "../five_theta/check_physical.cpp"
#include <numeric>
#include <random>

namespace {
Matrix identity() {
  Matrix v{};
  for (int i = 0; i < 5; i++)
    v[i * 6] = 1;
  return v;
}
Matrix tr(Matrix a) {
  Matrix t{};
  for (int i = 0; i < 5; i++)
    for (int j = 0; j < 5; j++)
      t[i * 5 + j] = a[j * 5 + i];
  return t;
}
int trace(Matrix a) {
  int t = 0;
  for (int i = 0; i < 5; i++)
    t += a[6 * i];
  return mod(t);
}
int inner(V a, V b) {
  long x = 0;
  for (int i = 0; i < 5; i++)
    x += long(a[i]) * b[i];
  return mod(x);
}
V addv(V a, V b, int c = 1) {
  for (int i = 0; i < 5; i++)
    a[i] = mod(a[i] + c * b[i]);
  return a;
}
V mv(Matrix a, V x) {
  V out{};
  for (int i = 0; i < 5; i++)
    for (int j = 0; j < 5; j++)
      out[i] = mod(out[i] + a[5 * i + j] * x[j]);
  return out;
}
Matrix mm(Matrix a, Matrix b) {
  Matrix out{};
  for (int i = 0; i < 5; i++)
    for (int j = 0; j < 5; j++)
      for (int k = 0; k < 5; k++)
        out[i * 5 + j] = mod(out[i * 5 + j] + a[i * 5 + k] * b[k * 5 + j]);
  return out;
}
Matrix wedge(V x, V y) {
  Matrix a;
  for (int i = 0; i < 5; i++)
    for (int j = 0; j < 5; j++)
      a[5 * i + j] = mod(x[i] * y[j] - y[i] * x[j]);
  return a;
}
Matrix bivector(V x, V y, V z) {
  Matrix a{};
  V perm{0, 1, 2, 3, 4};
  do {
    auto [i, j, k, l, m] = perm;
    a[5 * i + j] =
        mod(a[5 * i + j] + long(signature(perm)) * x[k] * y[l] * z[m]);
  } while (std::next_permutation(perm.begin(), perm.end()));
  return a;
}
int inverse0(int x) {
  if (!x)
    throw std::logic_error("zero transfer scalar");
  int n = P - 2, y = 1;
  while (n) {
    if (n & 1)
      y = mod(y * x);
    x = mod(x * x);
    n >>= 1;
  }
  return y;
}
struct Transport {
  Matrix g = identity(), dual = identity();
  int scale = 1;
  void append(V p, V u, V k) {
    int alpha = inner(k, u);
    scale = mod(scale * alpha);
    int ia = inverse0(alpha);
    Matrix a = identity(), b = a;
    for (int i = 0; i < 5; i++)
      for (int j = 0; j < 5; j++) {
        a[5 * i + j] = mod(a[5 * i + j] + long(ia) * u[i] * p[j]);
        b[5 * i + j] = mod(b[5 * i + j] - long(ia) * p[i] * u[j]);
      }
    g = mm(g, a);
    dual = mm(dual, b);
  }
};
struct Fixture {
  int gamma;
  std::array<int, 12> axis{};
  std::array<V, 12> p{};
  std::array<V, 3> initial{};
  Fixture(int j, std::mt19937 &rng) : gamma(j + 1) {
    for (int i = 1; i <= 10; i++)
      if (i != gamma) {
        axis[i] = rng() % 5;
        for (int k = 0; k < 5; k++)
          if (k != axis[i]) {
            p[i][k] = int(rng() % 7) - 3;
            p[gamma][k] -= p[i][k];
          }
      }
    axis[gamma] = rng() % 5;
    for (int i = 0; i < 5; i++) {
      initial[1][i] = 12 + rng() % 5;
      initial[2][i] = 13 + rng() % 5;
      initial[0][i] = -initial[1][i] - initial[2][i];
    }
  }
  std::array<V, 13> currents(int vary = 0, int scale = 1) const {
    auto init = initial;
    if (vary) {
      for (int i = 0; i < 5; i++) {
        init[vary][i] *= scale;
        init[0][i] = -init[1][i] - init[2][i];
      }
    }
    std::array<V, 13> y{};
    for (int h = 0; h < 3; h++) {
      V k = init[h];
      for (int j = 0; j < int(edgeids[h].size()); j++) {
        y[edgeids[h][j]] = k;
        if (j < int(vertices[h].size()))
          for (int i = 0; i < 5; i++)
            k[i] += p[vertices[h][j]][i];
      }
    }
    return y;
  }
  int marked() const {
    for (int h = 0; h < 3; h++)
      if (std::find(vertices[h].begin(), vertices[h].end(), gamma) !=
          vertices[h].end())
        return h;
    throw std::logic_error("one-form site");
  }
  std::array<int, 4> reduced(const std::array<V, 13> &y) const {
    int b = marked(), a = (b + 2) % 3, c = (b + 1) % 3;
    Transport reference, left, right, other;
    for (int h = 0; h < 3; h++)
      for (int j = 0; j < int(vertices[h].size()); j++) {
        int v = vertices[h][j];
        if (v == gamma)
          continue;
        V u{};
        u[axis[v]] = 1;
        auto &segment = h == a      ? reference
                        : h == c    ? other
                        : v < gamma ? left
                                    : right;
        segment.append(p[v], u, y[edgeids[h][j]]);
      }
    V A = addv({}, y[edgeids[a].front()]), B = addv({}, y[edgeids[b].front()]),
      C = addv({}, y[edgeids[c].front()]),
      E = mv(reference.dual, y[edgeids[b].back()]),
      F = mv(reference.dual, y[edgeids[c].back()]);
    V eta{};
    eta[axis[gamma]] = 1;
    eta = mv(left.dual, eta);
    V q = mv(left.dual, p[gamma]);
    auto refinv = tr(reference.dual), reftrans = tr(reference.g);
    Matrix G = mm(mm(left.g, right.g), refinv),
           D = mm(mm(left.dual, right.dual), reftrans), H = mm(other.g, refinv),
           J = mm(other.dual, reftrans);
    if (mv(D, E) != addv(B, q) || mv(tr(H), C) != F)
      throw std::logic_error("transported momenta");
    Matrix S = mm(wedge(eta, q), G), T = mm(bivector(eta, B, q), D);
    std::array<int, 4> r{};
    auto JE = mv(J, E), htB = mv(tr(H), B);
    auto SJ = mm(S, tr(J)), sH = mm(tr(S), H);
    V per{0, 1, 2, 3, 4};
    do {
      auto [i, j, k, a, b] = per;
      long s = signature(per);
      r[0] = mod(r[0] - s * B[a] * C[b] % P *
                            mod(S[j * 5 + i] * JE[k] - E[i] * SJ[5 * j + k]));
      r[1] = mod(r[1] + s * E[a] * F[b] % P *
                            mod(S[5 * i + j] * htB[k] - B[i] * sH[5 * j + k]));
    } while (std::next_permutation(per.begin(), per.end()));
    auto t = [&](V first, V second) { return inner(first, mv(T, second)); };
    auto htA = mv(tr(H), A), htF = mv(tr(H), F), ja = mv(J, A), jc = mv(J, C);
    r[2] = mod(-(-t(htA, F) + t(F, htA) + trace(H) * t(A, F) - t(A, htF)));
    r[3] = mod(-t(C, ja) + trace(J) * t(C, A) + t(ja, C) - t(jc, A));
    for (auto &x : r)
      x = mod(long(x) * left.scale % P * right.scale % P * other.scale % P *
              reference.scale);
    return r;
  }
};
} // namespace
int main() {
  const std::array<std::array<int, 3>, 5> shapes{
      {{1, 6, 6}, {2, 5, 6}, {3, 4, 6}, {3, 5, 5}, {4, 4, 5}}};
  std::mt19937 rng(83591);
  for (int n = 0; n < 5; n++) {
    vertices = {};
    edgeids = {};
    topology(shapes[n]);
    FlowPlans plans(1);
    const int shape_sign = (n == 1 || n == 3) ? -1 : 1;
    int count = 0, nonzero = 0;
    std::array<int, 2> highest{-1, -1}, highest_joint{-1, -1};
    for (int trial = 0; trial < 50; trial++) {
      Fixture z(trial % 10, rng);
      if (n == 0)
        for (int h = 1; h <= 2; h++) {
          std::array<int, 9> samples{};
          for (int j = 0; j < 9; j++)
            samples[j] = raw_numerator(plans.by_sites[1 << (z.gamma - 1)],
                                       z.currents(h, j - 4), z.axis, z.p);
          for (int d = 0; d < 9; d++) {
            if (samples[0])
              highest[h - 1] = std::max(highest[h - 1], d);
            for (int j = 0; j < 8 - d; j++)
              samples[j] = mod(samples[j + 1] - samples[j]);
          }
        }
      for (int mode = 0; mode < 2; mode++) {
        std::array<int, 14> values{};
        for (int j = 0; j < 14; j++) {
          auto input = z;
          for (int h = 1; h <= 2; h++)
            for (int k = 0; k < 5; k++) {
              int direction =
                  mode ? h * (1 + k % 3) : ((2 * k + h * k + h) % 5 - 2);
              input.initial[h][k] = z.initial[h][k] + (j - 7) * direction;
              input.initial[0][k] = -input.initial[1][k] - input.initial[2][k];
            }
          values[j] = raw_numerator(plans.by_sites[1 << (z.gamma - 1)],
                                    input.currents(), input.axis, input.p);
        }
        for (int d = 0; d < 14; d++) {
          if (values[0])
            highest_joint[mode] = std::max(highest_joint[mode], d);
          for (int j = 0; j < 13 - d; j++)
            values[j] = mod(values[j + 1] - values[j]);
        }
      }
      auto y = z.currents();
      auto reduced = z.reduced(y);
      std::array<bool, 4> found{};
      int b = z.marked(), a = (b + 2) % 3, c = (b + 1) % 3;
      std::array<int, 3> order{a, b, c};
      for (const auto &flow : plans.by_sites[1 << (z.gamma - 1)]) {
        int id = flow.local[0]->cubic                                       ? 0
                 : flow.local[11]->cubic                                    ? 1
                 : flow.local[0]->words[0].first_slot == edgeids[a].front() ? 2
                                                                            : 3;
        int external = (10 - z.gamma) & 1 ? -1 : 1;
        int primitive = external * flow.wick;
        for (int h = 0; h < 3; h++)
          for (int j = 0; j < int(vertices[h].size()); j++) {
            int v = vertices[h][j];
            if (v != z.gamma &&
                flow.local[v]->words[0].first_slot == edgeids[h][j])
              primitive = -primitive;
          }
        for (int v : {0, 11})
          if (!flow.local[v]->cubic) {
            int differentiated = flow.local[v]->words[0].first_slot, last = -1;
            for (int h : order)
              if (differentiated !=
                  (v ? edgeids[h].back() : edgeids[h].front())) {
                if (h < last)
                  primitive = -primitive;
                last = h;
              }
          }
        if (primitive != shape_sign * (id == 2 ? -1 : 1))
          throw std::logic_error("external and local signs");
        int raw = mod(external * raw_numerator({flow}, y, z.axis, z.p));
        if (found[id] || raw != mod(shape_sign * reduced[id]))
          throw std::logic_error("literal orientation mismatch");
        found[id] = true;
        nonzero += raw != 0;
        count++;
      }
    }
    if (n == 0) {
      if (highest != std::array<int, 2>{4, 4})
        throw std::logic_error("sharp separate degrees");
      std::cout << "166 separate degrees " << highest[0] << ' ' << highest[1]
                << "\n";
    }
    if (highest_joint[0] > 8 || highest_joint[1] > 7)
      throw std::logic_error("universal degrees");
    std::cout << "joint and collinear degrees " << highest_joint[0] << " "
              << highest_joint[1] << " ";
    std::cout << shapes[n][0] << shapes[n][1] << shapes[n][2] << " shape sign "
              << shape_sign << ", matched primitive orientations " << count
              << ", nonzero " << nonzero << "\n";
  }
}
