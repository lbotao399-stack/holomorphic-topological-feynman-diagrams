// Sparse rank lower bound over F_65521. Each vector passed in uses CSR row order.
#include <algorithm>
#include <cstdint>
#include <random>
#include <vector>

namespace {
using U = std::uint32_t;
using W = std::uint64_t;
constexpr U prime = 65521;

U inverse(U x) {
    U answer = 1;
    for (U power = prime - 2; power; power >>= 1) {
        if (power & 1) answer = W(answer) * x % prime;
        x = W(x) * x % prime;
    }
    return answer;
}
}

extern "C" int krylov_rank(int nrows, int ncols, const int* offsets, const int* cols,
                          const U* entries, int goal, U seed) {
    if (goal <= 0 || ncols == 0 || nrows == 0) return 0;
    std::vector<int> transpose_offsets, transpose_cols;
    std::vector<U> transpose_entries;
    if (ncols > nrows) {
        transpose_offsets.assign(ncols + 1, 0);
        for (int j = 0; j < offsets[nrows]; ++j) ++transpose_offsets[cols[j] + 1];
        for (int j = 0; j < ncols; ++j) transpose_offsets[j+1] += transpose_offsets[j];
        auto next = transpose_offsets;
        transpose_cols.resize(offsets[nrows]);
        transpose_entries.resize(offsets[nrows]);
        for (int i = 0; i < nrows; ++i) {
            for (int j = offsets[i]; j < offsets[i+1]; ++j) {
                int k = next[cols[j]]++;
                transpose_cols[k] = i;
                transpose_entries[k] = entries[j];
            }
        }
        std::swap(nrows, ncols);
        offsets = transpose_offsets.data();
        cols = transpose_cols.data();
        entries = transpose_entries.data();
    }

    std::mt19937 rng(seed);
    auto random_nonzero = [&]() -> U { return 1 + rng() % (prime - 1); };
    std::vector<U> left(ncols), v(ncols), row_diagonal(nrows), col_diagonal(ncols);
    for (auto* array : {&left, &v, &row_diagonal, &col_diagonal}) {
        std::generate(array->begin(), array->end(), random_nonzero);
    }
    std::vector<W> sums(ncols);
    std::vector<U> sequence, polynomial{1}, previous{1};
    int length = 0, shift = 1;
    U last_inverse = 1;
    for (int step = 0; step < 2*goal; ++step) {
        std::fill(sums.begin(), sums.end(), W(0));
        for (int i = 0; i < nrows; ++i) {
            W dot = 0;
            for (int j = offsets[i]; j < offsets[i+1]; ++j) dot += W(entries[j]) * v[cols[j]];
            U weight = (dot % prime) * row_diagonal[i] % prime;
            for (int j = offsets[i]; j < offsets[i+1]; ++j) sums[cols[j]] += W(entries[j]) * weight;
        }
        W sum = 0;
        for (int j = 0; j < ncols; ++j) {
            v[j] = (sums[j] % prime) * col_diagonal[j] % prime;
            sum += W(left[j]) * v[j];
        }
        sequence.push_back(sum % prime);
        W difference = sequence.back();
        for (int i = 1; i <= length; ++i) difference += W(polynomial[i]) * sequence[step-i];
        U delta = difference % prime;
        if (delta == 0) {
            ++shift;
            continue;
        }
        const bool grow = 2*length <= step;
        std::vector<U> saved;
        if (grow) saved = polynomial;
        U factor = W(delta) * last_inverse % prime;
        if (polynomial.size() < previous.size() + shift) polynomial.resize(previous.size() + shift);
        for (std::size_t i = 0; i < previous.size(); ++i) {
            auto& out = polynomial[i+shift];
            out = (out + prime - W(factor) * previous[i] % prime) % prime;
        }
        if (grow) {
            length = step + 1 - length;
            if (length >= goal) return length;
            previous.swap(saved);
            last_inverse = inverse(delta);
            shift = 1;
        } else {
            ++shift;
        }
    }
    return length;
}
