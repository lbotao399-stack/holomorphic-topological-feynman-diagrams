#pragma once
#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <utility>

using Edges = std::array<std::array<int, 6>, 6>;
using Powers = std::array<int, 9>;

// Finite points 0,1 are the triangle faces, 2,3,4 are the rectangular faces.
constexpr std::array<std::array<int, 2>, 9> faces = {{
    {2,3}, {0,2}, {0,3}, {1,3}, {1,2}, {0,4}, {2,4}, {3,4}, {1,4}
}};

Edges complete(const Powers &powers, int q) {
    Edges a{};
    for (int i = 0; i < 9; ++i) {
        const auto [u, v] = faces[i];
        a[u][v] = a[v][u] = powers[i];
    }
    for (int i = 0; i < 5; ++i) {
        int left = 2*q;
        for (int j = 0; j < 5; ++j) left -= a[i][j];
        a[i][5] = a[5][i] = left;
    }
    return a;
}

bool has_bubble(const Edges &a) {
    for (const auto &row : a)
        if (std::count_if(row.begin(), row.end(), [](int x) { return x > 0; }) <= 3)
            return true;
    return false;
}

// A positive convergent prism chart has one missing pair and eight cut lines.
// Every original row is an eligible chart. Encode its nine powers in base q.
std::uint64_t canonical(const Edges &a, int q) {
    std::uint64_t best = std::numeric_limits<std::uint64_t>::max();
    for (int infinity = 0; infinity < 6; ++infinity) {
        int zeros = 0;
        std::array<int, 2> apices{};
        bool eligible = true;
        for (int i = 0; i < 6; ++i) for (int j = i+1; j < 6; ++j) {
            if (i == infinity || j == infinity) continue;
            int power = a[i][j];
            if (power < 0 || power >= q) eligible = false;
            if (power == 0) { apices = {i,j}; ++zeros; }
        }
        if (!eligible || zeros != 1) continue;
        std::array<int, 5> p = {apices[0], apices[1], 0, 0, 0};
        int next = 2;
        for (int i = 0; i < 6; ++i)
            if (i != infinity && i != p[0] && i != p[1]) p[next++] = i;
        do {
            for (int orientation = 0; orientation < 2; ++orientation) {
                std::uint64_t code = 0;
                for (const auto [i,j] : faces) code = code*q + a[p[i]][p[j]];
                best = std::min(best, code);
                std::swap(p[0], p[1]);
            }
        } while (std::next_permutation(p.begin()+2, p.end()));
    }
    return best;
}

using Key = std::array<signed char,14>;
constexpr std::array<std::array<int,2>,14> key_pairs = {{
    {0,2},{1,2},{0,3},{1,3},{0,4},{1,4},{0,5},{1,5},
    {2,3},{2,4},{2,5},{3,4},{3,5},{4,5}
}};

Key graph_key(const Edges &a) {
    Key best;
    best.fill(127);
    for (int x = 0; x < 6; ++x) for (int y = 0; y < 6; ++y) {
        if (x == y || a[x][y]) continue;
        std::array<int,6> p{x,y,0,0,0,0};
        int next = 2;
        for (int z = 0; z < 6; ++z) if (z != x && z != y) p[next++] = z;
        auto pair_less = [&](int i, int j) {
            return std::make_pair(a[x][i],a[y][i]) < std::make_pair(a[x][j],a[y][j]);
        };
        std::sort(p.begin()+2,p.end(),[&](int i,int j) {
            return pair_less(i,j) || (!pair_less(j,i) && i < j);
        });
        do {
            Key k{};
            for (int i = 0; i < 14; ++i) k[i] = a[p[key_pairs[i][0]]][p[key_pairs[i][1]]];
            if (k < best) best = k;
            // Enumerate only permutations within tied blocks of neighbor weights.
            bool more = false;
            int end = 6;
            while (end > 2 && !more) {
                int first = end-1;
                while (first > 2 && !pair_less(p[first-1],p[first])) --first;
                more = std::next_permutation(p.begin()+first,p.begin()+end);
                end = first;
            }
            if (!more) break;
        } while (true);
    }
    return best;
}

Edges unpack(const Key &key) {
    Edges a{};
    for (int i = 0; i < 14; ++i) {
        const auto [u,v] = key_pairs[i];
        a[u][v] = a[v][u] = key[i];
    }
    return a;
}

Powers decode(std::uint64_t code, int q) {
    Powers p{};
    for (int i = 8; i >= 0; --i) { p[i] = code % q; code /= q; }
    return p;
}

