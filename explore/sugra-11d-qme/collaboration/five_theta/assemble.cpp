#include "matching.hpp"
#include "multi_gamma.hpp"
#include "gamma_options.hpp"
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace {
using namespace twoloop;
constexpr std::array<int,10> waves={1,2,3,4,5,-1,-2,-4,-3,-5};
int wave_index(int momentum) {
    for(int label=0;label<10;++label)if(waves[label]==momentum)return label;
    throw std::logic_error("unknown wave");
}

struct Transported {
    Slots slots{};
    GammaSites sites;
    int anchor_slot=-1, factor=1;
};
Transported transport(const Word& word,const Word& sides,const GammaOptions& options,
                      const std::array<int,5>& permutation,const std::array<int,5>& inverse,
                      int sign_mask) {
    const int count_gamma=__builtin_popcount(options.mask), m=(count_gamma+1)/2;
    std::array<int,5> signs{};
    int det=permutation_sign(permutation);
    for(int i=0;i<5;++i){signs[i]=(sign_mask&(1<<i))?-1:1;det*=signs[i];}
    Transported transformed;
    transformed.factor=determinant(word)*((m&1)?det:1);
    std::vector<int> odd;odd.reserve(10-count_gamma);
    const int anchor=__builtin_ctz(options.mask);
    for(int v=0;v<10;++v) {
        int i=word[v];
        int label=wave_index(sides[v]*signs[i]*(1+permutation[i]));
        int axis=inverse[options.polarizations[label]];
        transformed.slots[v]={sides[v]*(i+1),axis};
        transformed.factor*=signs[axis];
        if((options.mask>>label)&1U)transformed.sites.set(v);else odd.push_back(label);
        if(label==anchor)transformed.anchor_slot=v;
    }
    if(int(odd.size())!=10-count_gamma || transformed.anchor_slot<0)throw std::logic_error("external field mismatch");
    transformed.factor*=permutation_sign(odd);
    return transformed;
}

Polynomial aggregate(const Lengths& lengths,const Word& word,int wanted_slot,const GammaOptions& options) {
    std::array<int,5> permutation={0,1,2,3,4};
    const auto sides=momentum_signs(word);
    detail::Collected factored;
    int assignments=0;
    Transported transformed;
    do {
        std::array<int,5> inverse{};
        for(int i=0;i<5;++i)inverse[permutation[i]]=i;
        for(int mask=0;mask<32;++mask) {
            transformed=transport(word,sides,options,permutation,inverse,mask);
            if(wanted_slot>=0 && transformed.anchor_slot!=wanted_slot)continue;
            detail::Graph graph(lengths,transformed.slots,transformed.anchor_slot);
            detail::collect_gamma_sites(graph,transformed.sites,transformed.factor,factored);
            ++assignments;
        }
    }while(std::next_permutation(permutation.begin(),permutation.end()));
    if(assignments!=(wanted_slot<0?3840:384))throw std::logic_error("assignment count mismatch");
    detail::Graph graph(lengths,transformed.slots,0);
    return detail::expand(graph,factored);
}
} // namespace
int main(int argc,char** argv) {
    GammaOptions options;
    if(!options.parse(argc,argv) || options.positional.size()<6 || options.positional.size()>7){
        std::cerr<<"Usage: assemble p q r first end output [gamma_slot] [--gamma-mask bits] [--polarizations i1,...,i10]\n"
                 <<"Physical wave order: +1,+2,+3,+4,+5,-1,-2,-4,-3,-5; gamma_slot pins the first selected wave.\n";return 2;
    }
    const auto& args=options.positional;
    const twoloop::Lengths lengths={std::stoi(args[0]),std::stoi(args[1]),std::stoi(args[2])};
    const int first=std::stoi(args[3]),end=std::stoi(args[4]),slot=args.size()==7?std::stoi(args[6]):-1;
    auto patterns=twoloop::matchings(lengths);
    if(first<0||first>=end||end>int(patterns.size())||slot< -1||slot>=10)return 2;
    std::ofstream out(args[5]); if(!out)return 2;
    for(int rank=first;rank<end;++rank){
        auto start=std::chrono::steady_clock::now();
        auto p=aggregate(lengths,patterns[rank],slot,options);
        std::vector<std::pair<std::uint64_t,std::int64_t>> terms;
        for(auto term:p)if(term.second)terms.push_back(term);
        std::sort(terms.begin(),terms.end());
        out<<rank<<' '<<terms.size()<<'\n';for(auto [key,coefficient]:terms)out<<key<<' '<<coefficient<<'\n';
        out.flush();if(!out)return 1;
        if(rank==first || (rank+1)%10==0)std::cerr<<"rank "<<rank<<'/'<<patterns.size()<<" terms "<<terms.size()<<" seconds "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<'\n';
    }
    return 0;
}
