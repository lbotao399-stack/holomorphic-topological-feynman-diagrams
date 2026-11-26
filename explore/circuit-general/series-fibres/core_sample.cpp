// Exact finite-field evaluation of a series-integrated reciprocal numerator.
// The Python wrapper validates the graph/source conventions and supplies the
// exact contracted determinant data. No floating-point operations occur.
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

using word = std::uint64_t;

static inline word mul(word a, word b, word prime) { return (a * b) % prime; }
static inline word add(word a, word b, word prime) {
    word value = a + b;
    return value >= prime ? value - prime : value;
}
static inline word sub(word a, word b, word prime) {
    return a >= b ? a - b : a + prime - b;
}
static word power(word value, word exponent, word prime) {
    word result = 1;
    while (exponent) {
        if (exponent & 1) result = mul(result, value, prime);
        value = mul(value, value, prime);
        exponent >>= 1;
    }
    return result;
}

// All matrices are row-major. The sources/constant elimination multiplier
// may be rational before reduction modulo the prime. The wrapper rejects
// primes dividing any of their denominators.
extern "C" int htk_core_samples(
    std::size_t count, int n, int H, int m, int blocks, int gauge,
    const std::int64_t *conductances, word prime,
    const std::int64_t *cycles, const std::int64_t *last_flow,
    const std::int64_t *interior_sources, const int *block_chains,
    const std::int64_t *lift, word determinant_factor,
    int node_count, const std::int64_t *nodes, const std::int64_t *weights,
    int threads, std::int64_t *output) {
    if (n < 1 || n > 6 || H < 1 || m < 1 || blocks < 0 || gauge < 0 || gauge >= m
        || prime <= 5 || prime >= (word(1) << 31) || node_count < 1 || threads < 1)
        return 1;
    const int dimension = m - 1;
    const int components = (blocks + 1) * H;
    const int width = n + components;
    const std::size_t square = std::size_t(dimension) * dimension;
    std::vector<int> chain_block(m, -1);
    for (int b = 0; b < blocks; ++b) chain_block[block_chains[b]] = b;

    #pragma omp parallel num_threads(threads)
    {
        std::vector<word> resistance(m), input_prefix(m + 1);
        std::vector<word> augmented(std::size_t(n) * width);
        std::vector<word> means(std::size_t(m) * components);
        std::vector<word> coefficients(std::size_t(blocks + 1) * square);
        std::vector<word> matrices(std::size_t(node_count) * square);
        std::vector<word> determinants(node_count), pivots(node_count), inverse_pivots(node_count);
        std::vector<word> pivot_prefix(node_count + 1);

        #pragma omp for schedule(static)
        for (std::size_t point = 0; point < count; ++point) {
            const std::int64_t *g = conductances + point * m;
            bool singular = false;
            input_prefix[0] = 1;
            word reciprocal_factor = determinant_factor;
            for (int c = 0; c < m; ++c) {
                if (!g[c]) singular = true;
                input_prefix[c + 1] = mul(input_prefix[c], g[c], prime);
                reciprocal_factor = mul(reciprocal_factor, power(g[c], n - 1, prime), prime);
            }
            if (singular) { output[point] = -1; continue; }
            word inverse_product = power(input_prefix[m], prime - 2, prime);
            for (int c = m - 1; c >= 0; --c) {
                resistance[c] = mul(inverse_product, input_prefix[c], prime);
                inverse_product = mul(inverse_product, g[c], prime);
            }

            // One core inverse, simultaneously applied to the constant source
            // and every coefficient of a chain-fraction source variation.
            std::fill(augmented.begin(), augmented.end(), 0);
            for (int c = 0; c < m; ++c) {
                for (int r = 0; r < n; ++r) {
                    word factor = mul(resistance[c], cycles[c * n + r], prime);
                    for (int s = 0; s < n; ++s)
                        augmented[r * width + s] = add(augmented[r * width + s],
                            mul(factor, cycles[c * n + s], prime), prime);
                    for (int i = 0; i < H; ++i)
                        augmented[r * width + n + i] = sub(augmented[r * width + n + i],
                            mul(factor, last_flow[c * H + i], prime), prime);
                    int b = chain_block[c];
                    if (b >= 0)
                        for (int i = 0; i < H; ++i)
                            augmented[r * width + n + (b + 1) * H + i] =
                                mul(factor, interior_sources[b * H + i], prime);
                }
            }
            word U = 1;
            for (int column = 0; column < n; ++column) {
                int pivot = column;
                while (pivot < n && !augmented[pivot * width + column]) ++pivot;
                if (pivot == n) { singular = true; break; }
                if (pivot != column) {
                    for (int j = column; j < width; ++j)
                        std::swap(augmented[pivot * width + j], augmented[column * width + j]);
                    U = U ? prime - U : 0;
                }
                word value = augmented[column * width + column];
                U = mul(U, value, prime);
                word inverse_value = power(value, prime - 2, prime);
                for (int j = column; j < width; ++j)
                    augmented[column * width + j] = mul(augmented[column * width + j], inverse_value, prime);
                for (int r = 0; r < n; ++r) {
                    if (r == column) continue;
                    word factor = augmented[r * width + column];
                    if (!factor) continue;
                    for (int j = column; j < width; ++j)
                        augmented[r * width + j] = sub(augmented[r * width + j],
                            mul(factor, augmented[column * width + j], prime), prime);
                }
            }
            if (singular || !U) { output[point] = -1; continue; }

            for (int c = 0; c < m; ++c) {
                for (int k = 0; k <= blocks; ++k) {
                    for (int i = 0; i < H; ++i) {
                        word value = k == 0 ? last_flow[c * H + i] : 0;
                        for (int r = 0; r < n; ++r)
                            value = add(value, mul(cycles[c * n + r],
                                augmented[r * width + n + k * H + i], prime), prime);
                        if (k && chain_block[c] == k - 1)
                            value = sub(value, interior_sources[(k - 1) * H + i], prime);
                        means[c * components + k * H + i] = value;
                    }
                }
            }
            std::fill(coefficients.begin(), coefficients.end(), 0);
            int kept = 0;
            for (int c = 0; c < m; ++c) {
                if (c == gauge) continue;
                for (int row = 0; row < dimension; ++row)
                    for (int k = 0; k <= blocks; ++k) {
                        word value = 0;
                        for (int i = 0; i < H; ++i)
                            value = add(value, mul(lift[(c * dimension + row) * H + i],
                                means[c * components + k * H + i], prime), prime);
                        coefficients[std::size_t(k) * square + row * dimension + kept] = value;
                    }
                ++kept;
            }
            for (int node = 0; node < node_count; ++node) {
                for (std::size_t entry = 0; entry < square; ++entry) {
                    word value = coefficients[entry];
                    for (int k = 0; k < blocks; ++k)
                        value = add(value, mul(nodes[node * blocks + k],
                            coefficients[std::size_t(k + 1) * square + entry], prime), prime);
                    matrices[std::size_t(node) * square + entry] = value;
                }
                determinants[node] = 1;
            }

            // Batch inversion of the cubature determinants' pivots: one
            // exponentiation per column, rather than one per node/column.
            for (int column = 0; column < dimension; ++column) {
                pivot_prefix[0] = 1;
                for (int node = 0; node < node_count; ++node) {
                    word *matrix = matrices.data() + std::size_t(node) * square;
                    word value = 1;
                    if (determinants[node]) {
                        int pivot = column;
                        while (pivot < dimension && !matrix[pivot * dimension + column]) ++pivot;
                        if (pivot == dimension) determinants[node] = 0;
                        else {
                            if (pivot != column) {
                                for (int j = column; j < dimension; ++j)
                                    std::swap(matrix[pivot * dimension + j], matrix[column * dimension + j]);
                                determinants[node] = prime - determinants[node];
                            }
                            value = matrix[column * dimension + column];
                        }
                    }
                    pivots[node] = value;
                    pivot_prefix[node + 1] = mul(pivot_prefix[node], value, prime);
                }
                word inverse = power(pivot_prefix[node_count], prime - 2, prime);
                for (int node = node_count - 1; node >= 0; --node) {
                    inverse_pivots[node] = mul(inverse, pivot_prefix[node], prime);
                    inverse = mul(inverse, pivots[node], prime);
                }
                for (int node = 0; node < node_count; ++node) {
                    if (!determinants[node]) continue;
                    determinants[node] = mul(determinants[node], pivots[node], prime);
                    word *matrix = matrices.data() + std::size_t(node) * square;
                    for (int row = column + 1; row < dimension; ++row) {
                        word factor = mul(matrix[row * dimension + column], inverse_pivots[node], prime);
                        if (!factor) continue;
                        for (int j = column + 1; j < dimension; ++j)
                            matrix[row * dimension + j] = sub(matrix[row * dimension + j],
                                mul(factor, matrix[column * dimension + j], prime), prime);
                        matrix[row * dimension + column] = 0;
                    }
                }
            }
            word integral = 0;
            for (int node = 0; node < node_count; ++node)
                integral = add(integral, mul(weights[node], determinants[node], prime), prime);
            reciprocal_factor = mul(reciprocal_factor, power(U, n, prime), prime);
            reciprocal_factor = mul(reciprocal_factor, g[gauge], prime);
            output[point] = mul(reciprocal_factor, integral, prime);
        }
    }
    return 0;
}
