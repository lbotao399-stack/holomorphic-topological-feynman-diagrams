// Exact S_9 orbit evaluations modulo 1009 for the scalar quartic's transverse tests.
#include <algorithm>
#include <array>
#include <iostream>
#include <numeric>
#include <vector>

using Five = std::array<int, 5>;
constexpr int prime = 1009;

int key(const Five& ids) {
    int base = 0;
    for (int i : ids) base = 9 * base + i;
    return base;
}

int main(int argc, char** argv) {
    const bool scalar_slots = argc == 2;
    std::array<Five, 9> p;
    for (auto& v : p) for (int& x : v) if (!(std::cin >> x)) return 1;
    int count;
    std::cin >> count;
    std::vector<std::array<int, 4>> orbits(count);
    for (auto& o : orbits) for (int& j : o) std::cin >> j;
    std::vector<int> table(59049);
    std::vector<Five> subsets;
    for (int i = 0; i < 5; ++i) for (int j = i+1; j < 6; ++j)
    for (int k = j+1; k < 7; ++k) for (int l = k+1; l < 8; ++l)
    for (int m = l+1; m < 9; ++m) {
        Five ids{i,j,k,l,m}, perm{0,1,2,3,4};
        int determinant = 0;
        do {
            int term = 1;
            for (int a = 0; a < 5; ++a) {
                for (int b = 0; b < a; ++b) if (perm[b] > perm[a]) term = -term;
                term *= p[ids[a]][perm[a]];
            }
            determinant += term;
        } while (std::next_permutation(perm.begin(), perm.end()));
        subsets.push_back(ids);
        do {
            int sign = 1;
            for (int a = 0; a < 5; ++a) for (int b = 0; b < a; ++b)
                if (ids[b] > ids[a]) sign = -sign;
            table[key(ids)] = ((determinant * sign) % prime + prime) % prime;
        } while (std::next_permutation(ids.begin(), ids.end()));
    }
    std::vector<std::array<long long, 81>> slots(count);
    std::array<int,9> g;
    std::iota(g.begin(), g.end(), 0);
    do {
        std::array<int,126> values;
        for (int i = 0; i < 126; ++i) {
            int index = 0;
            for (int j : subsets[i]) index = 9*index + g[j];
            values[i] = table[index];
        }
        int cell = 9*g[0]+g[1];
        for (int c = 0; c < count; ++c) {
            auto o = orbits[c];
            slots[c][cell] += 1LL * values[o[0]] * values[o[1]] * values[o[2]] * values[o[3]];
        }
    } while (std::next_permutation(g.begin(), g.end()));
    for (const auto& col : slots) {
        if (scalar_slots) {
            for (int k = 0; k < 5; ++k) {
                long long v = 0;
                for (int i = 0; i < 9; ++i) for (int j = 0; j < 9; ++j)
                    v += (col[9*i+j] % prime)*p[i][k];
                std::cout << (v % prime + prime) % prime << ' ';
            }
            std::cout << '\n';
            continue;
        }
        for (int i = 0; i < 9; ++i) {
            int pivot = 0;
            while (!p[i][pivot]) ++pivot;
            for (int k = 0; k < 5; ++k) if (k != pivot) {
                long long v = 0;
                for (int j = 0; j < 9; ++j)
                    v += (col[9*i+j] % prime)*(p[i][pivot]*p[j][k]-p[i][k]*p[j][pivot]);
                std::cout << (v % prime + prime) % prime << ' ';
            }
        }
        std::cout << '\n';
    }
}
