// Complete the planar dual to six weighted points and group equivalent cuts.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <gmpxx.h>
#include <limits>
#include <string>
#include <vector>

#include "quartic_prism_wheel_charts.hpp"

using Integer = __int128_t;
struct Row {
    std::uint64_t key;
    Integer k, j;
};

Integer scaled(const std::string &entry, long long scale) {
    auto slash = entry.find('/');
    auto numerator = std::stoll(entry.substr(0, slash));
    auto denominator = slash == std::string::npos ? 1 : std::stoll(entry.substr(slash+1));
    if (denominator <= 0 || scale % denominator) throw std::invalid_argument("bank denominator");
    return Integer(numerator) * (scale / denominator);
}

std::string decimal(Integer x) {
    if (!x) return "0";
    bool negative = x < 0;
    if (negative) x = -x;
    std::string result;
    while (x) { result.push_back('0' + int(x % 10)); x /= 10; }
    if (negative) result.push_back('-');
    std::reverse(result.begin(), result.end());
    return result;
}

mpq_class rational(Integer x) { return mpq_class(decimal(x)); }
const mpq_class &rational(const mpq_class &x) { return x; }

struct RationalRow {
    std::uint64_t key;
    mpq_class k, j;
};

// A scalar three-line uniqueness integral leaves a five-point completion.
// Opposite tetrahedron edges exchange the dual position and momentum graphs.
struct ScalarCut {
    std::array<int,6> powers{};
    mpq_class ratio{};
};

ScalarCut scalar_cut(const Edges &a, int q, const std::vector<mpz_class> &factorial) {
    for (int v=0; v<6; ++v) {
        std::vector<int> neighbors;
        for (int i=0; i<6; ++i)
            if (a[v][i]) neighbors.push_back(i);
        if (neighbors.size()!=3 ||
            std::any_of(neighbors.begin(),neighbors.end(),[&](int u){return a[v][u]<=0;})) continue;
        Edges b=a;
        mpq_class ratio=1;
        for (int h=0;h<3;++h) {
            const int i=neighbors[(h+1)%3], j=neighbors[(h+2)%3], power=a[v][neighbors[h]];
            b[i][j]=b[j][i]=a[i][j]+q-power;
            ratio*=mpq_class(factorial[q-power-1],factorial[power-1]);
        }
        ratio.canonicalize();
        for (int inf=0; inf<6;++inf) if (inf!=v) {
            ScalarCut cut;
            int k=6;
            for (int i=0;i<6;++i) if (i!=inf && i!=v)
                for (int j=i+1;j<6;++j) if (j!=inf && j!=v) cut.powers[--k]=b[i][j];
            if (*std::min_element(cut.powers.begin(),cut.powers.end())<=0 ||
                *std::max_element(cut.powers.begin(),cut.powers.end())>=q) continue;
            cut.ratio=ratio;
            return cut;
        }
    }
    return {};
}

// Sending one point to infinity and another to the origin leaves four angular
// vertices. Two absent opposite pairs give a single Gegenbauer summation.
std::uint64_t wheel_key(const Edges &g) {
    for (int i=0; i<6; ++i) for (int j=i+1; j<6; ++j) {
        if (g[i][j]) continue;
        for (int k=i+1; k<6; ++k) for (int l=k+1; l<6; ++l) {
            if (k==j || l==j || g[k][l]) continue;
            for (int c=0; c<6; ++c) {
                if (c==i || c==j || c==k || c==l) continue;
                // Rim order i,k,j,l. Zero at (0,0) uses the existing wheel-bank packing.
                Powers a{0,g[i][k],g[j][l],g[l][i],g[c][i],g[c][l],
                         g[k][j],g[c][k],g[c][j]};
                if (*std::min_element(a.begin()+1,a.end())<=0) continue;
                std::uint64_t packed=0;
                for (int power : a) packed=(packed<<4)|power;
                return packed;
            }
        }
    }
    return 0;
}

template<class Bank>
void reflect_bank(const Bank &bank, int q, long scale, bool emit, bool keep_bubbles=false, bool tensor_charts=false, int wheel_charts=0) {
    std::vector<mpz_class> factorial(q,1);
    for (int i = 1; i < q; ++i) factorial[i] = factorial[i-1]*i;
    std::vector<bool> seen(bank.size(),false);
    long bubbles = 0, groups = 0, keep_k = 0, keep_j = 0, cancelled = 0, largest = 0;
    for (std::size_t start = 0; start < bank.size(); ++start) {
        if (seen[start]) continue;
        auto initial = graph_key(complete(decode(bank[start].key,q),q));
        // Scalar factors relate the seed cut to each visited cut.
        std::map<Key,mpq_class> factors;
        std::vector<Key> queue{initial};
        factors[initial] = 1;
        bool bubble = false;
        ScalarCut terminal;
        Key bubble_key{};
        mpq_class bubble_ratio;
        int bubble_score=1000000;
        std::uint64_t wheel=0;
        mpq_class wheel_ratio;
        std::vector<std::size_t> members;
        mpq_class k = 0, j = 0;
        for (std::size_t i = 0; i < queue.size(); ++i) {
            const Key current = queue[i];
            const Edges a = unpack(current);
            const mpq_class &factor = factors.at(current);
            auto code = canonical(a,q);
            auto it = std::lower_bound(bank.begin(),bank.end(),code,
                [](const auto &r, std::uint64_t c) { return r.key < c; });
            if (it != bank.end() && it->key == code) {
                std::size_t index = it-bank.begin();
                if (seen[index]) throw std::logic_error("overlapping reflected components");
                seen[index] = true;
                members.push_back(index);
                k += rational(it->k)/factor;
                j += rational(it->j)/factor;
            }
            if (wheel_charts && !wheel) {
                if(wheel_charts<0) wheel=numerator_wheel_key(a,q,-wheel_charts-1,false);
                else wheel=wheel_charts>=2 ? numerator_wheel_key(a,q,wheel_charts-1) : wheel_key(a);
                if (wheel) wheel_ratio=factor;
            }
            const bool boundary=has_bubble(a);
            bubble |= boundary;
            if (tensor_charts && boundary) {
                for (const auto &row : a) {
                    if (std::count_if(row.begin(),row.end(),[](int x){return x>0;})!=3) continue;
                    int score=0;
                    for (int w : row) score-=100*std::min(w,0);
                    for (const auto &other : a)
                        for (int w : other) score-=std::min(w,0);
                    if (score<bubble_score) {
                        bubble_score=score;
                        bubble_key=current;
                        bubble_ratio=factor;
                    }
                }
            }
            if (keep_bubbles && boundary && terminal.ratio==0) {
                terminal=scalar_cut(a,q,factorial);
                terminal.ratio*=factor;
            }
            for (int v = 0; v < 6; ++v) {
                std::vector<int> neighbors;
                for (int u = 0; u < 6; ++u) if (a[v][u]) neighbors.push_back(u);
                if (neighbors.size() != 4 ||
                    std::any_of(neighbors.begin(),neighbors.end(),
                                [&](int u) { return a[v][u] <= 0 || a[v][u] >= q; })) continue;
                Edges b = a;
                mpq_class ratio = 1;
                for (int u : neighbors) {
                    int power = a[v][u];
                    b[v][u] = b[u][v] = q-power;
                    ratio *= mpq_class(factorial[q-power-1],factorial[power-1]);
                    for (int w : neighbors) if (w > u)
                        b[u][w] = b[w][u] = a[u][w] + power + a[v][w] - q;
                }
                ratio.canonicalize();
                Key child = graph_key(b);
                mpq_class new_factor = factor*ratio;
                auto [pos, fresh] = factors.emplace(child,new_factor);
                if (fresh) queue.push_back(child);
                else if (new_factor != pos->second) throw std::logic_error("star factor mismatch");
            }
        }
        largest = std::max(largest,long(queue.size()));
        if (bubble) {
            bubbles += members.size();
            if (!keep_bubbles) continue;
        }
        ++groups;
        keep_k += k != 0;
        keep_j += j != 0;
        cancelled += k == 0 && j == 0;
        if (emit && (k != 0 || j != 0)) {
            if (wheel_charts && !wheel) continue;
            for (int a : decode(bank[start].key,q)) std::cout << a << ' ';
            std::cout << k/scale << ' ' << j/scale;
            if (wheel_charts) {
                std::cout << ' ' << wheel << ' ' << wheel_ratio;
            } else if (tensor_charts) {
                if (!bubble) throw std::invalid_argument("requested chart has no bubble");
                for (int w : bubble_key) std::cout << ' ' << w;
                std::cout << ' ' << bubble_ratio;
            } else if (keep_bubbles) {
                std::cout << ' ' << bubble;
                if (terminal.ratio!=0) {
                    for (int a : terminal.powers) std::cout << ' ' << a;
                    std::cout << ' ' << terminal.ratio;
                }
            }
            std::cout << '\n';
        }
    }
    std::cerr << bubbles << ' ' << groups << ' ' << keep_k << ' '
              << keep_j << ' ' << cancelled << ' ' << largest << '\n';
}

// Group a sparse bank with arbitrary rational coefficients, retaining every
// rational, weight-three and weight-five contribution of tensor bubbles.
int full_components(int q, bool tensor_charts=false, int wheel_charts=0) {
    std::map<std::uint64_t, RationalRow> rows;
    Powers a;
    long count = 0;
    int margin = q;
    while (std::cin >> a[0]) {
        for (int i = 1; i < 9; ++i) if (!(std::cin >> a[i])) return 1;
        mpq_class k, j;
        if (!(std::cin >> k >> j)) return 1;
        const Edges g = complete(a,q);
        for (int i=0;i<6;++i) {
            int degree=0;
            for(int w : g[i]) degree+=w;
            if (degree!=2*q) throw std::invalid_argument("completed degree");
            for (int j=i+1;j<6;++j) {
                margin=std::min(margin,q-g[i][j]);
                for (int h=j+1;h<6;++h)
                    margin=std::min(margin,2*q-g[i][j]-g[i][h]-g[j][h]);
            }
        }
        if (margin<=0) throw std::invalid_argument("divergent prism cut");
        const auto key=canonical(g,q);
        if (key==std::numeric_limits<std::uint64_t>::max()) return 1;
        auto &row=rows[key];
        row.key=key;
        row.k+=k;
        row.j+=j;
        ++count;
    }
    if (!std::cin.eof()) return 1;
    std::vector<RationalRow> bank;
    for (const auto &pair : rows)
        if (pair.second.k!=0 || pair.second.j!=0) bank.push_back(pair.second);
    std::cerr << "rows " << count << " completed " << rows.size() << " nonzero "
              << bank.size() << " minimum margin " << margin << '\n';
    reflect_bank(bank,q,1,true,true,tensor_charts,wheel_charts);
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 3 || argc > 4) return 2;
    if (argc==4 && (std::string(argv[2])=="--numerator-wheel-charts" ||
                    std::string(argv[2])=="--formal-numerator-wheel-charts")) {
        int q=std::stoi(argv[1]),degree=std::stoi(argv[3]);
        int mode=std::string(argv[2])=="--numerator-wheel-charts" ? degree+1 : -(degree+1);
        return q>=3 && q<=16 && degree>0 && degree<=2*q-4 ? full_components(q,false,mode) : 2;
    }
    if (argc==3 && (std::string(argv[2])=="--full-components" || std::string(argv[2])=="--bubble-charts" || std::string(argv[2])=="--wheel-charts" || std::string(argv[2])=="--linear-wheel-charts" || std::string(argv[2])=="--quadratic-wheel-charts")) {
        int q=std::stoi(argv[1]);
        int wheels=std::string(argv[2])=="--quadratic-wheel-charts" ? 3 :
                   std::string(argv[2])=="--linear-wheel-charts" ? 2 : int(std::string(argv[2])=="--wheel-charts");
        return q>=2 && q<=(wheels?16:24) ? full_components(q,std::string(argv[2])=="--bubble-charts",wheels) : 2;
    }
    auto scale = std::stoll(argv[2]);
    std::string mode = argc == 4 ? argv[3] : "";
    if (scale <= 0 || (mode != "" && mode != "--weights" && mode != "--reflect"
                                      && mode != "--reflected-weights")) return 2;
    int q = std::stoi(argv[1]);
    if (q < 2 || q > 24) return 2;
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    Powers a;
    std::vector<Row> active;
    std::array<long, 3> erased{};
    long count = 0;
    while (std::cin >> a[0]) {
        for (int i = 1; i < 9; ++i) if (!(std::cin >> a[i])) return 1;
        std::string k, j;
        if (!(std::cin >> k >> j)) return 1;
        const Edges graph = complete(a, q);
        if (has_bubble(graph)) {
            ++erased[0];
            erased[1] += k != "0";
            erased[2] += j != "0";
        } else {
            auto code = canonical(graph, q);
            if (code == std::numeric_limits<std::uint64_t>::max()) return 1;
            active.push_back({code, scaled(k, scale), scaled(j, scale)});
        }
        ++count;
    }
    std::sort(active.begin(), active.end(), [](const Row &a, const Row &b) {
        return a.key < b.key;
    });
    std::vector<Row> combined;
    std::size_t used = 0, distinct = 0, nonzero_k = 0, nonzero_j = 0, cancelled = 0;
    while (used < active.size()) {
        Row row = active[used++];
        while (used < active.size() && active[used].key == row.key) {
            row.k += active[used].k;
            row.j += active[used++].j;
        }
        ++distinct;
        nonzero_k += row.k != 0;
        nonzero_j += row.j != 0;
        cancelled += !row.k && !row.j;
        if (mode == "--reflect" || mode == "--reflected-weights") combined.push_back(row);
        if (mode == "--weights" && (row.k || row.j)) {
            Powers powers{};
            for (int i = 8; i >= 0; --i) { powers[i] = row.key % q; row.key /= q; }
            for (int power : powers) std::cout << power << ' ';
            std::cout << decimal(row.k) << ' ' << decimal(row.j) << '\n';
        }
    }
    std::cerr << count << ' ' << erased[0] << ' ' << erased[1] << ' '
              << erased[2] << ' ' << active.size() << ' ' << distinct << ' '
              << nonzero_k << ' ' << nonzero_j << ' ' << cancelled << '\n';
    if (!combined.empty()) reflect_bank(combined,q,scale,mode == "--reflected-weights");
}
