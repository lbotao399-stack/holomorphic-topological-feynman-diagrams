// Integer angular moments of sparse polynomial products. A checked failure asks
// the caller to use its arbitrary-precision implementation for this block.
#include <cstdint>
#include <limits>
#include <algorithm>
#include <unordered_map>
#include <vector>

using Wide = __int128_t;
using Unsigned = __uint128_t;

// The right polynomials stay fixed while every pair of cofactor rows is tested.
// Reuse their factorial-moment vector whenever the same left monomial occurs.
struct MomentTable {
    int n,count,limit;
    int64_t factorial[21]{1},rising[21][21]{};
    std::vector<int> starts,slots,ends,highest;
    std::vector<uint8_t> shifts;
    std::vector<int64_t> coeffs;
    std::unordered_map<uint64_t,std::vector<int64_t>> cache;
    MomentTable(int nv,int columns,const int32_t *off,const uint8_t *powers,const int64_t *weights)
        :n(nv),count(columns),limit(std::min(131072,24000000/(64+8*columns))),
         starts(off[columns]+1),ends(off,off+columns+1),highest(nv),coeffs(weights,weights+off[columns]) {
        for(int i=1;i<=20;++i)factorial[i]=factorial[i-1]*i;
        for(int i=0;i<=20;++i)for(int j=0;i+j<=20;++j)rising[i][j]=factorial[i+j]/factorial[i];
        for(int t=0;t<ends.back();++t){
            starts[t]=slots.size();
            for(int k=0;k<n;++k)if(int p=powers[t*n+k]){
                highest[k]=std::max(highest[k],p);slots.push_back(k);shifts.push_back(p);
            }
        }
        starts.back()=slots.size();
    }
    bool evaluate(int na,const uint8_t* powers,const int64_t* coeff,uint64_t *out){
        std::vector<Wide> totals(count);
        std::vector<int64_t> work(count);
        for(int i=0;i<na;++i){
            const auto ex=powers+i*n;
            uint64_t key=0;Wide scale=coeff[i];
            for(int k=0;k<n;++k){
                if(ex[k]+highest[k]>20||
                   __builtin_mul_overflow(scale,static_cast<Wide>(factorial[ex[k]]),&scale))return false;
                key=(key<<5)|ex[k];
            }
            auto hit=cache.find(key);
            const int64_t* row=work.data();
            if(hit!=cache.end())row=hit->second.data();
            else {
                for(int j=0;j<count;++j){
                    int64_t sum=0;
                    for(int t=ends[j];t<ends[j+1];++t){
                        int64_t weight=coeffs[t];
                        for(int s=starts[t];s<starts[t+1];++s)
                            if(__builtin_mul_overflow(weight,rising[ex[slots[s]]][shifts[s]],&weight))return false;
                        if(__builtin_add_overflow(sum,weight,&sum))return false;
                    }
                    work[j]=sum;
                }
                if(int(cache.size())<limit)cache.emplace(key,work);
            }
            for(int j=0;j<count;++j)if(row[j]){
                Wide term;
                if(__builtin_mul_overflow(scale,static_cast<Wide>(row[j]),&term)||
                   __builtin_add_overflow(totals[j],term,&totals[j]))return false;
            }
        }
        for(int j=0;j<count;++j){
            auto v=static_cast<Unsigned>(totals[j]);
            out[2*j]=static_cast<uint64_t>(v);out[2*j+1]=static_cast<uint64_t>(v>>64);
        }
        return true;
    }
};

extern "C" void* ordered_moment_table(int n,int count,const int32_t *ends,const uint8_t *powers,const int64_t *coeffs){
    if(n<=0||n>12||count<=0)return nullptr;
    return new MomentTable(n,count,ends,powers,coeffs);
}
extern "C" void free_ordered_moment_table(void* table){delete static_cast<MomentTable*>(table);}
extern "C" int table_moments(void* table,int na,const uint8_t* powers,const int64_t* coeffs,uint64_t* output){
    return static_cast<MomentTable*>(table)->evaluate(na,powers,coeffs,output);
}

// Factors have at most 12 variables in five-bit lanes. The Python adapter
// checks the largest componentwise exponent before adding the packed keys.
extern "C" int product_moments(void* table,int nx,const uint64_t* x,const int64_t* a,
        int ny,const uint64_t* y,const int64_t* b,uint64_t* out){
    auto* angular=static_cast<MomentTable*>(table);
    if(static_cast<int64_t>(nx)*ny>4000000)return 0;
    std::unordered_map<uint64_t,int64_t> terms;
    for(int i=0;i<nx;++i)for(int j=0;j<ny;++j){
        int64_t value;
        if(__builtin_mul_overflow(a[i],b[j],&value))return 0;
        auto& coeff=terms[x[i]+y[j]];
        if(__builtin_add_overflow(coeff,value,&coeff))return 0;
    }
    std::vector<uint8_t> powers;
    std::vector<int64_t> coeffs;
    powers.reserve(angular->n*terms.size());
    coeffs.reserve(terms.size());
    for(auto [key,value]:terms)if(value){
        for(int k=angular->n-1;k>=0;--k)
            powers.push_back((key>>(5*k))&31);
        coeffs.push_back(value);
    }
    return angular->evaluate(coeffs.size(),powers.data(),coeffs.data(),out);
}

extern "C" int ordered_moments(
    int n, int na, const uint8_t* ae, const int64_t* ac, int count,
    const int32_t* ends, const uint8_t* be, const int64_t* bc,
    uint64_t* output) {
    if (n <= 0 || na < 0 || count < 0) return 0;
    constexpr int max_index = 20;
    int64_t factorial[max_index+1]{1};
    for (int i=1; i<=max_index; ++i) factorial[i] = factorial[i-1]*i;
    int64_t rising[max_index+1][max_index+1]{};
    for (int i=0; i<=max_index; ++i)
        for (int j=0; i+j<=max_index; ++j)
            rising[i][j] = factorial[i+j]/factorial[i];

    const int nb = ends[count];
    std::vector<int> starts(nb+1), slots;
    std::vector<uint8_t> shifts;
    std::vector<int> highest(n);
    for (int t=0; t<nb; ++t) {
        starts[t] = static_cast<int>(slots.size());
        for (int k=0; k<n; ++k) {
            const int p = be[t*n+k];
            if (p > highest[k]) highest[k] = p;
            if (p) {
                slots.push_back(k);
                shifts.push_back(p);
            }
        }
    }
    starts[nb] = static_cast<int>(slots.size());
    std::vector<Wide> sums(count);
    for (int i=0; i<na; ++i) {
        const uint8_t* ex = ae+i*n;
        Wide scale = ac[i];
        for (int k=0; k<n; ++k) {
            if (ex[k]+highest[k] > max_index ||
                __builtin_mul_overflow(scale,static_cast<Wide>(factorial[ex[k]]),&scale)) return 0;
        }
        for (int j=0; j<count; ++j) {
            int64_t sum = 0;
            for (int q=ends[j]; q<ends[j+1]; ++q) {
                int64_t weight = bc[q];
                for (int h=starts[q]; h<starts[q+1]; ++h) {
                    const int64_t f = rising[ex[slots[h]]][shifts[h]];
                    if (__builtin_mul_overflow(weight,f,&weight)) return 0;
                }
                if (__builtin_add_overflow(sum,weight,&sum)) return 0;
            }
            Wide term;
            if (__builtin_mul_overflow(scale,static_cast<Wide>(sum),&term) ||
                __builtin_add_overflow(sums[j],term,&sums[j])) return 0;
        }
    }
    for (int j=0; j<count; ++j) {
        const Unsigned v = static_cast<Unsigned>(sums[j]);
        output[2*j] = static_cast<uint64_t>(v);
        output[2*j+1] = static_cast<uint64_t>(v >> 64);
    }
    return 1;
}
