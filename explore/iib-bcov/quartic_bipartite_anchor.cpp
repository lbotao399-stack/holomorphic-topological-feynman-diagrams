// Lift all three q=3 convergent projective K3,3 periods by multiplication by U^(q-3).
// Rows in the output are packed row-major powers followed by three exact orbit sums.
#define BCOV_K33_ORBIT_LIBRARY
#include "quartic_bipartite_orbit_ibp.cpp"

std::vector<Matrix> cotree_shifts(){
    std::vector<Matrix> shifts;
    for(int mask=0;mask<(1<<9);++mask){
        if(__builtin_popcount(unsigned(mask))!=4)continue;
        int seen=1;
        for(int step=0;step<5;++step){
            for(int j=0;j<9;++j)if(!(mask&(1<<j))){
                int ends=(1<<(j/3))|(1<<(3+j%3));
                if(seen&ends)seen|=ends;
            }
        }
        if(seen!=63)continue;
        Matrix shift{};
        for(int j=0;j<9;++j)shift[j]=(mask>>j)&1;
        shifts.push_back(shift);
    }
    if(shifts.size()!=81)throw std::logic_error("spanning tree count");
    return shifts;
}

int main(int argc,char**argv){
    if(argc!=3)throw std::invalid_argument("q and tsv output required");
    int dimension=std::stoi(argv[1]);
    if(dimension<3||dimension>12)throw std::invalid_argument("expected 3 <= q <= 12");
    auto shifts=cotree_shifts();
    auto seeds=Domain(3).keys;
    if(seeds.size()!=3)throw std::logic_error("unexpected terminal rank");
    using Multiplicities=std::array<uint64_t,3>;
    std::unordered_map<Key,Multiplicities> current;
    for(int j=0;j<3;++j)current[seeds[j]][j]=1;
    std::cerr<<"terminal keys "<<seeds[0]<<" "<<seeds[1]<<" "<<seeds[2]<<"\n";
    uint64_t norm=1;
    for(int q=3;q<dimension;++q){
        std::unordered_map<Key,Multiplicities> next;
        next.reserve(current.size()*4);
        for(const auto& [key,counts]:current){
            auto powers=unpack(key);
            for(const auto& shift:shifts){
                Matrix p;
                for(int i=0;i<9;++i)p[i]=powers[i]+shift[i];
                auto& target=next[canonical(p)];
                for(int j=0;j<3;++j)target[j]+=counts[j];
            }
        }
        current.swap(next);
        norm*=81; // 81^9 < 2^58, also bounding every intermediate orbit sum.
        Multiplicities mass{};
        for(const auto& [key,row]:current){
            if(!finite(unpack(key),q+1))throw std::logic_error("nonconvergent lift");
            for(int j=0;j<3;++j)mass[j]+=row[j];
        }
        if(mass!=Multiplicities{norm,norm,norm})throw std::logic_error("lost orbit mass");
        std::cerr<<"q "<<q+1<<" supported orbits "<<current.size()<<std::endl;
    }
    std::vector<Key> ordered;
    for(const auto& kv:current)ordered.push_back(kv.first);
    std::sort(ordered.begin(),ordered.end());
    std::ofstream output(argv[2]);
    if(!output)throw std::runtime_error("cannot open anchor output");
    for(auto key:ordered){
        auto row=current.at(key);
        output<<key<<' '<<row[0]<<' '<<row[1]<<' '<<row[2]<<'\n';
    }
}
