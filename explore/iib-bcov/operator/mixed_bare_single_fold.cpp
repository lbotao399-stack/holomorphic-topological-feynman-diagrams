// Exterior ordered-simplex moments of two integer polynomials on five tetrahedron paths.
#include <array>
#include <cstdint>
#include <limits>
#include <vector>

namespace {
using Integer = __int128;
struct Term { uint64_t monomial; unsigned path_degree; int64_t coefficient; };
std::array<std::vector<uint64_t>, 6> moment_tables;
std::vector<Integer> sums(1u << 20);
std::vector<unsigned char> active(1u << 20);

const std::vector<uint64_t>& path_moments(unsigned size) {
    auto& table = moment_tables[size];
    if (!table.empty()) return table;
    const unsigned end = 1u << (4*size);
    table.resize(end);
    std::array<uint64_t,16> factorial{1};
    for (int j=1; j<16; ++j) factorial[j] = j*factorial[j-1];
    for (unsigned i=0; i<end; ++i) {
        uint64_t product=1;
        unsigned degree=0;
        for (unsigned j=0; j<size; ++j) {
            const unsigned exponent=(i>>(4*j))&15;
            degree+=exponent;
            if (degree>15) { product=0; break; }
            product*=factorial[exponent];
        }
        table[i]=product;
    }
    return table;
}

std::vector<Term> terms(unsigned count, const uint64_t* monomials, const int64_t* coefficients,
                        const unsigned* sizes, std::array<unsigned,5>& maxima) {
    std::vector<Term> result;
    result.reserve(count);
    for (unsigned i=0;i<count;++i) {
        unsigned key=0,offset=0;
        for (unsigned k=0;k<5;++k) {
            unsigned degree=0;
            for (unsigned j=0;j<sizes[k];++j) degree+=(monomials[i]>>(4*offset++))&15;
            if (degree>maxima[k]) maxima[k]=degree;
            key|=degree<<(4*k);
        }
        result.push_back({monomials[i],key,coefficients[i]});
    }
    return result;
}
}

extern "C" int fold(unsigned nl,const uint64_t* lm,const int64_t* lc,
                    unsigned nr,const uint64_t* rm,const int64_t* rc,const unsigned* sizes,
                    unsigned capacity,uint32_t* keys,uint64_t* low,int64_t* high) {
    std::array<unsigned,5> left_max{},right_max{};
    auto left=terms(nl,lm,lc,sizes,left_max);
    auto right=terms(nr,rm,rc,sizes,right_max);
    struct Split { unsigned shift,mask; const uint64_t* moments; };
    std::array<Split,5> splits{};
    unsigned split_count=0,offset=0;
    for (unsigned i=0;i<5;++i) {
        if (left_max[i]+right_max[i]>15 || sizes[i]>5 || !sizes[i]) return -1;
        if (sizes[i]>1)
            splits[split_count++]={4*offset,(1u<<(4*sizes[i]))-1,path_moments(sizes[i]).data()};
        offset+=sizes[i];
    }
    unsigned used=0;
    int error=0;
    for (const auto& a:left) {
        for (const auto& b:right) {
            const unsigned degree=a.path_degree+b.path_degree;
            if (!active[degree]) {
                if (used==capacity) { error=-2;break; }
                active[degree]=1;
                keys[used++]=degree;
            }
            Integer value=static_cast<Integer>(a.coefficient)*b.coefficient;
            const auto combined=a.monomial+b.monomial;
            for (unsigned j=0;j<split_count;++j) {
                const auto& s=splits[j];
                if (__builtin_mul_overflow(value,s.moments[(combined>>s.shift)&s.mask],&value)) {
                    error=-3;break;
                }
            }
            if (error || __builtin_add_overflow(sums[degree],value,&sums[degree])) {error=-3;break;}
        }
        if (error) break;
    }
    for (unsigned i=0;i<used;++i) {
        const auto value=sums[keys[i]];
        low[i]=static_cast<uint64_t>(value);
        high[i]=static_cast<int64_t>(value>>64);
        sums[keys[i]]=0;
        active[keys[i]]=0;
    }
    return error ? error : static_cast<int>(used);
}
