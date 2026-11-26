// Exact weighted nonplanar master projection of the retained K3,3 cut orbits.
#include <gmpxx.h>
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

using Exponents = std::array<int,8>;
using Integer = mpz_class;
using Rational = mpq_class;

struct Entry {
    Exponents exponent{};
    int coefficient{},offset{};
    uint64_t packed{};
};

struct Target {
    Exponents exponent{};
    std::array<Rational,2> weight{};
    int order{},flat{};
};

struct State {
    uint64_t packed{};
    int flat{};
    uint8_t through{};
};

struct Coefficient {
    // Arithmetic modulo 2^192. The proven one-norm bounds for every S_r
    // below are below 2^191, so centered readback is exact.
    uint64_t words[3]{};
    bool zero() const {return !(words[0]|words[1]|words[2]);}
    void add(const Coefficient& input,int multiplier){
        __int128 carry=0;
        for(int j=0;j<3;++j){
            const __int128 t=static_cast<__int128>(input.words[j])*multiplier+words[j]+carry;
            words[j]=static_cast<uint64_t>(t);
            carry=t>>64;
        }
    }
    Integer value() const {
        Integer v;
        mpz_import(v.get_mpz_t(),3,-1,sizeof(uint64_t),0,0,words);
        if(words[2]>>63) v-=Integer(1)<<192;
        return v;
    }
};

uint64_t pack(const Exponents& ex){
    uint64_t v=0;
    for(int j=0;j<8;++j) v|=uint64_t(ex[j])<<(5*j);
    return v;
}

int main(int argc,char** argv){
    const bool inspect=argc>1 && std::string(argv[1])=="--inspect";
    std::array<std::array<int,2>,2> pairs{};
    std::array<std::array<std::vector<Entry>,2>,2> operations;
    for(int endpoint=0;endpoint<2;++endpoint){
        std::cin>>pairs[endpoint][0]>>pairs[endpoint][1];
        for(auto& entries:operations[endpoint]){
            int count; std::cin>>count; entries.resize(count);
            for(auto& entry:entries){
                for(auto& x:entry.exponent) std::cin>>x;
                std::cin>>entry.coefficient;
            }
        }
    }
    std::vector<Target> targets;
    Exponents largest{};
    int last=0;
    for(Target row;std::cin>>row.exponent[0];){
        for(int j=1;j<8;++j) std::cin>>row.exponent[j];
        std::cin>>row.weight[0]>>row.weight[1];
        row.order=std::accumulate(row.exponent.begin(),row.exponent.end(),-20);
        if(row.order<0) throw std::logic_error("low-degree test target");
        last=std::max(last,row.order);
        for(int j=0;j<8;++j) largest[j]=std::max(largest[j],row.exponent[j]);
        targets.push_back(row);
    }
    if(targets.empty()) throw std::runtime_error("empty projection input");
    Exponents stride{};
    int volume=1;
    for(int j=7;j>=0;--j){
        stride[j]=volume;
        if(largest[j]>15) throw std::logic_error("packed exponent overflow");
        volume*=largest[j]+1;
    }
    auto encode=[&](const Exponents& a){
        return std::inner_product(a.begin(),a.end(),stride.begin(),0);
    };
    std::vector<int8_t> needed(volume,-1);
    for(auto& row:targets){
        row.flat=encode(row.exponent);
        needed[row.flat]=std::max<int>(needed[row.flat],row.order);
    }
    // Downward-closed support of each required coefficient array.
    for(int j=7;j>=0;--j){
        int width=stride[j],block=width*(largest[j]+1);
        for(int start=0;start<volume;start+=block)
            for(int at=start+block-1;at>=start+width;--at)
                needed[at-width]=std::max(needed[at-width],needed[at]);
    }
    std::vector<int> index(volume,-1);
    std::vector<State> states;
    states.reserve(std::count_if(needed.begin(),needed.end(),[](int x){return x>=0;}));
    Exponents a{};
    for(int i=0;i<volume;++i){
        if(needed[i]>=0){
            index[i]=int(states.size());
            states.push_back({pack(a),i,uint8_t(needed[i])});
        }
        for(int j=7;j>=0;--j){if(++a[j]<=largest[j])break;a[j]=0;}
    }
    Exponents bias{},mask{};
    for(int j=0;j<8;++j){bias[j]=15-largest[j];mask[j]=16;}
    const uint64_t cap=pack(bias),overflow=pack(mask);
    std::cerr<<"targets "<<targets.size()<<", rectangular cells "<<volume
             <<", downward closure "<<states.size()<<"\n";
    if(inspect)return 0;
    std::vector<Coefficient> first(states.size()),second(states.size()),result(states.size());
    std::array<Rational,2> total{};
    for(int endpoint=0;endpoint<2;++endpoint){
        auto& ops=operations[endpoint];
        std::array<int,2> norm{};
        for(int j=0;j<2;++j){
            auto& terms=ops[j];
            for(auto& t:terms){
                norm[j]+=abs(t.coefficient);
                t.packed=pack(t.exponent);t.offset=encode(t.exponent);
            }
        }
        std::fill(first.begin(),first.end(),Coefficient{});
        std::fill(second.begin(),second.end(),Coefficient{});
        second[index[0]].words[0]=1;
        Integer n0=0,n1=1;
        Integer factorial=1;
        for(int r=1;r<=last;++r){
            std::fill(result.begin(),result.end(),Coefficient{});
            const std::array<int,2> scales{21-2*r,4*(r-1)*(21-r)};
            Integer bound=abs(scales[0])*norm[0]*n1+abs(scales[1])*norm[1]*n0;
            if(bound>=(Integer(1)<<191)) throw std::overflow_error("signed coefficient bound");
            n0=n1;n1=bound;factorial*=r;
            for(int j=0;j<2;++j){
                const int scale=scales[j];
                if(!scale)continue;
                const auto& data=j?first:second;
                for(size_t i=0;i<states.size();++i){
                    const auto& s=states[i];
                    if(s.through<r || data[i].zero())continue;
                    for(const auto& t:ops[j]){
                        if((s.packed+cap+t.packed)&overflow)continue;
                        const int next=s.flat+t.offset;
                        if(needed[next]<r)continue;
                        result[index[next]].add(data[i],scale*t.coefficient);
                    }
                }
            }
            for(const auto& row:targets)if(row.order==r){
                const auto& pair=pairs[endpoint];
                int offset=stride[pair[0]]+stride[pair[1]];
                Integer binomial=1,value=0;
                for(int l=0;l<=std::min(row.exponent[pair[0]],row.exponent[pair[1]]);++l){
                    value+=binomial*result[index[row.flat-l*offset]].value();
                    binomial*=-2*(19-2*r-2*l);
                    if(!mpz_divisible_ui_p(binomial.get_mpz_t(),l+1)) throw std::logic_error("binomial");
                    binomial/=l+1;
                }
                Rational scaled(value,(Integer(1)<<r)*factorial);
                scaled.canonicalize();
                for(int channel=0;channel<2;++channel)total[channel]+=scaled*row.weight[channel];
            }
            first.swap(second);second.swap(result);
        }
        std::cerr<<"endpoint "<<endpoint<<" through homogeneous degree "<<last<<"\n";
    }
    Rational normalization(1,Integer(8)*(Integer(1)<<38));
    for(int j=0;j<10;++j) normalization*=Rational(8,(2*j+1)*(2*j+2)*(2*j+3));
    for(auto& v:total){v*=normalization;v.canonicalize();std::cout<<v<<"\n";}
}
