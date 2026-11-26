// Exact Dirichlet integration of sparse polynomials in the first path currents.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <gmpxx.h>
#include <iostream>
#include <stdexcept>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;
using i128 = __int128;
using Order = std::array<int,5>;
constexpr u64 factorial17=355687428096000ULL;
constexpr u64 path_scale=factorial17/120;
constexpr u64 weight_scale=21322087106945625ULL;
constexpr unsigned mask=(1u<<20)-1;

struct Split { unsigned powers,degree,multiple; };
std::array<unsigned,13> binomial_row(unsigned n) {
    std::array<unsigned,13> row={1};
    for(unsigned k=1;k<=n;k++) row[k]=row[k-1]*(n+1-k)/k;
    return row;
}
std::vector<Split> split(unsigned powers) {
    std::vector<Split> result={{0,0,1}};
    for(unsigned i=0;i<5;++i) {
        unsigned n=(powers>>(4*i))&15;
        auto choose=binomial_row(n);
        std::size_t count=result.size();
        for(unsigned k=1;k<=n;++k)
            for(std::size_t j=0;j<count;++j) {
                auto item=result[j];
                result.push_back({item.powers+(k<<(4*i)),item.degree+k,item.multiple*choose[k]});
            }
    }
    return result;
}

u64 odd_product(unsigned first, unsigned length) {
    u64 value=1;
    for(unsigned n=0;n<length;++n) value*=first+2*n;
    return value;
}
long scaled_weight(unsigned a,unsigned u,unsigned v,unsigned d) {
    // Mean of (wB-1)^a (-wB)^u wC^v (1-wC)^d, scaled by (13*15*...*35).
    auto ca=binomial_row(a), cd=binomial_row(d);
    i128 sum=0;
    for(unsigned i=0;i<=a;i++) for(unsigned j=0;j<=d;j++) {
        unsigned n=i+u, m=j+v;
        i128 value=(i128)(weight_scale/odd_product(13,n+m))*odd_product(1,n)*odd_product(1,m);
        value*=ca[i]*cd[j];
        if((a-i+u+j)%2) sum-=value; else sum+=value;
    }
    if(sum>weight_scale || sum<-(i128)weight_scale) throw std::logic_error("weight bound");
    return (long)sum;
}
struct PathMoments {
    Order order;
    std::vector<u64> values = std::vector<u64>(1u<<20);
    u64 operator()(unsigned powers) {
        if(values[powers]) return values[powers];
        unsigned tail=0; u64 den=1;
        for(int j=4;j>=0;--j) {
            tail += (powers>>(4*order[j]))&15;
            den *= 5-j+tail;
        }
        if(factorial17%den) throw std::logic_error("unsupported moment degree");
        return values[powers]=factorial17/den;
    }
};

mpz_class big(u128 n) {
    const u64 words[]={(u64)n,(u64)(n>>64)};
    mpz_class result;
    mpz_import(result.get_mpz_t(),2,-1,sizeof(u64),0,0,words);
    return result;
}
unsigned degree(unsigned key) {
    unsigned n=0;for(int i=0;i<5;i++) n+=(key>>(4*i))&15;
    return n;
}

mpz_class integrated_term(u64 key, PathMoments& path_b, PathMoments& path_c) {
    unsigned b=key&mask, c=key>>20;
    unsigned kb=degree(b),kc=degree(c);
    if(kb+kc>12) throw std::logic_error("maximum insertion degree exceeded");
    auto bs=split(b),cs=split(c);
    // Binomial multiplicities sum to 2^(kb+kc). Each path moment is <= path_scale.
    u128 groups[13][13]={};
    for(const auto& us:bs) for(const auto& vs:cs) {
        unsigned s=b-us.powers+vs.powers, t=c-vs.powers+us.powers;
        groups[us.degree][vs.degree]+=(u128)path_b(s)*path_c(t)*us.multiple*vs.multiple;
    }
    mpz_class result;
    static long weights[13][13][13][13]={};
    for(unsigned u=0;u<=kb;++u) for(unsigned v=0;v<=kc;++v) {
        if(!groups[u][v]) continue;
        long& w=weights[kb-u][u][v][kc-v];
        if(!w) w=scaled_weight(kb-u,u,v,kc-v);
        result+=big(groups[u][v])*w;
    }
    return result;
}

int main(int argc,char**argv) {
    if(argc<2){std::cerr<<"Usage: moments polynomial-chunk-file ...\n";return 2;}
    mpz_class sum;
    mpz_class scale=mpz_class((unsigned long)path_scale)*((unsigned long)path_scale)
                    *((unsigned long)weight_scale);
    PathMoments path_b{{0,1,2,3,4}};
    std::array<bool,120> seen{};
    int processed=0;
    for(int arg=1;arg<argc;++arg){
        std::ifstream input(argv[arg]);if(!input) return 2;
        int rank; size_t count;
        while(input>>rank>>count){
            if(rank<0 || rank>=120 || seen[rank]) throw std::logic_error("invalid relative order");
            seen[rank]=true;
            Order order={0,1,2,3,4};
            for(int t=0;t<rank;++t) std::next_permutation(order.begin(),order.end());
            PathMoments path_c{order};
            mpz_class partial;
            for(size_t j=0;j<count;j++){
                u64 key;long coef;
                if(!(input>>key>>coef)) throw std::logic_error("truncated input");
                partial+=integrated_term(key,path_b,path_c)*coef;
            }
            sum+=partial; ++processed;
            mpq_class value(partial,scale*158400);value.canonicalize();
            std::cout<<"order "<<rank<<" terms "<<count<<" kernel/pi "<<value<<"\n";
        }
        if(!input.eof()) throw std::logic_error("invalid polynomial data");
    }
    mpq_class total(sum,scale*158400);total.canonicalize();
    std::cout<<"orders "<<processed<<" signed kernel/pi "<<total<<"\n";
}
