// Evaluate P = U^6 D on the interpolation grid using the tree-current
// formula and the two 4-by-4 minors of the cycle matrix. All divisions
// and fixed-width integer operations affecting the determinant are checked.
#include <array>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

using Wide = __int128_t;
struct Tree {
    std::array<int, 4> missing;
    int current[9][2];
};
struct Minor {
    int sign;
    std::array<int, 8> columns;
};

int main(int argc, char** argv) {
    assert(argc == 3 || argc == 4);
    std::ifstream input(argv[1]);
    std::vector<Tree> trees(argc == 4 ? std::stoi(argv[3]) : 81);
    for (auto& tree : trees) {
        for (auto& edge : tree.missing) input >> edge;
        for (auto& row : tree.current) input >> row[0] >> row[1];
    }
    assert(input);
    std::vector<Minor> minors;
    Minor minor;
    while (input >> minor.sign) {
        for (auto& edge : minor.columns) input >> edge;
        assert(input);
        minors.push_back(minor);
    }

    std::ofstream output(argv[2], std::ios::binary);
    for (int n = 0; n < (1 << 18); n++) {
        int t[9];
        int key = n;
        for (int edge = 8; edge >= 0; edge--) {
            t[edge] = 1 + key % 4;
            key /= 4;
        }
        int u = 0;
        int y[9][2]{};
        for (const auto& tree : trees) {
            int weight = 1;
            for (int edge : tree.missing) weight *= t[edge];
            u += weight;
            for (int edge = 1; edge < 9; edge++) {
                for (int i = 0; i < 2; i++) y[edge][i] += weight * tree.current[edge][i];
            }
        }

        Wide determinant = 0;
        for (const auto& term : minors) {
            Wide product = term.sign;
            for (int j = 0; j < 8; j++) {
                const bool overflow = __builtin_mul_overflow(product, Wide(y[term.columns[j]][j/4]), &product);
                assert(!overflow);
            }
            const bool overflow = __builtin_add_overflow(determinant, product, &determinant);
            assert(!overflow);
        }
        const Wide denominator = Wide(u) * u * u * u * t[0];
        assert(determinant % denominator == 0);
        const auto p = static_cast<std::int64_t>(determinant / denominator);
        assert(Wide(p) * denominator == determinant);
        output.write(reinterpret_cast<const char*>(&p), sizeof(p));
    }
    assert(output);
}
