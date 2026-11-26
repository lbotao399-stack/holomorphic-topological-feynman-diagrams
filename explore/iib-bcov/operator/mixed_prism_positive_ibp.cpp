// Positive five-term face identities in a single triangular-prism chart.
#include "../quartic_prism_graph.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

struct PositiveIndex {
    std::vector<uint64_t> keys;
    std::vector<int> ids;
    size_t mask;
    explicit PositiveIndex(size_t n){
        size_t cap=1;
        while(cap*3/4<n)cap*=2;
        keys.assign(cap,0);ids.resize(cap);mask=cap-1;
    }
    size_t position(uint64_t k) const {
        uint64_t h=k;
        h=(h^(h>>30))*0xbf58476d1ce4e5b9ULL;
        h=(h^(h>>27))*0x94d049bb133111ebULL;
        size_t i=(h^(h>>31))&mask;
        while(keys[i] && keys[i]!=k)i=(i+1)&mask;
        return i;
    }
    void insert(uint64_t k,int id){
        auto i=position(k);
        if(keys[i])throw std::logic_error("Repeated finite prism");
        keys[i]=k;ids[i]=id;
    }
    int get(uint64_t k) const {auto i=position(k);return keys[i]?ids[i]:-1;}
};

std::array<std::array<int,9>,12> maps;
int labels[5][5];
uint64_t finite_key(const Powers& p){
    const int first=std::min({p[0],p[6],p[7]});
    uint64_t best=UINT64_MAX;
    for(const auto &g:maps)if(p[g[0]]==first){
        uint64_t k=0;
        for(int j:g)k=(k<<4)|p[j];
        best=std::min(best,k);
    }
    return best;
}

struct Row{
    std::array<int,7> id{};
    std::array<int8_t,7> c{};
    void add(int key,int n){
        if(!n)return;
        for(size_t i=0;i<id.size();++i)
            if(c[i] && id[i]==key){c[i]+=n;return;}
        for(size_t i=0;i<id.size();++i)
            if(!c[i]){id[i]=key;c[i]=n;return;}
        throw std::logic_error("Too many terms in a prism face equation");
    }
};

int main(int argc,char** argv){
    if(argc<7 || argc>9)throw std::invalid_argument("q angular.keys finite.aliases anchors.keys scratch plan [--all-finite [scaled.bin]]");
    const int q=std::stoi(argv[1]);
    const bool all_finite=argc>=8 && std::string(argv[7])=="--all-finite";
    if(argc>=8 && !all_finite)throw std::invalid_argument("Unknown prism face mode");
    std::ifstream seeds(argv[2]),aliases(argv[3],std::ios::binary),anchors(argv[4]);
    std::ofstream plan(argv[6]);
    if(!seeds||!aliases||!anchors||!plan)throw std::runtime_error("Cannot open prism inputs");
    std::vector<uint64_t> keys;
    uint64_t k;
    while(seeds>>k)keys.push_back(k);
    if(!std::is_sorted(keys.begin(),keys.end()))throw std::logic_error("Expected sorted angular keys");
    PositiveIndex index(std::filesystem::file_size(argv[3])/16);
    uint64_t angular;
    while(aliases.read(reinterpret_cast<char*>(&k),8)){
        if(!aliases.read(reinterpret_cast<char*>(&angular),8))throw std::runtime_error("Truncated aliases");
        auto it=std::lower_bound(keys.begin(),keys.end(),angular);
        if(it==keys.end()||*it!=angular)throw std::logic_error("Missing angular class");
        index.insert(k,int(it-keys.begin()));
    }
    std::vector<bool> known(keys.size(),false);
    while(anchors>>k){
        auto it=std::lower_bound(keys.begin(),keys.end(),k);
        if(it!=keys.end()&&*it==k)known[it-keys.begin()]=true;
    }
    const auto initial=std::count(known.begin(),known.end(),true);
    std::vector<uint64_t> scaled;
    uint64_t limbs=0;
    if(argc==9){
        std::ifstream input(argv[8],std::ios::binary);
        if(!input.read(reinterpret_cast<char*>(&limbs),8) || !limbs || limbs>32)
            throw std::runtime_error("Invalid fixed-width exact period table");
        scaled.resize(keys.size()*3*limbs);
        if(!input.read(reinterpret_cast<char*>(scaled.data()),scaled.size()*8) || input.peek()!=EOF)
            throw std::runtime_error("Inconsistent fixed-width exact period table");
    }
    for(int i=0;i<9;++i){auto [u,v]=faces[i];labels[u][v]=labels[v][u]=i;}
    std::array<int,5> perm{0,1,2,3,4};
    int count=0;
    do{
        for(int orientation=0;orientation<2;++orientation){
            for(int i=0;i<9;++i){auto [u,v]=faces[i];maps[count][i]=labels[perm[u]][perm[v]];}
            ++count;std::swap(perm[0],perm[1]);
        }
    }while(std::next_permutation(perm.begin()+2,perm.end()));
    uint64_t extra=0,rows=0;
    auto propagate=[&](const Row &r){
        if(!scaled.empty()){
            // The caller bounds the absolute rational numerator of every row
            // by less than 2^(64*limbs), so a zero residue proves equality.
            for(int col=0;col<3;++col){
                __int128 carry=0;
                for(size_t j=0;j<limbs;++j){
                    for(int i=0;i<7;++i)if(r.c[i]){
                        uint64_t w=scaled[(size_t(r.id[i])*3+col)*limbs+j];
                        carry+=static_cast<__int128>(w)*r.c[i];
                    }
                    if(uint64_t(carry)!=0)throw std::logic_error("Violated exact prism face equation");
                    carry>>=64;
                }
            }
            return true;
        }
        int missing=-1;
        for(int j=0;j<7;++j)if(r.c[j]&&!known[r.id[j]]){
            if(missing>=0)return false;
            missing=j;
        }
        if(missing>=0){
            ++extra;known[r.id[missing]]=true;
            plan<<keys[r.id[missing]]<<" "<<int(r.c[missing]);
            for(int i=0;i<7;++i)if(i!=missing && r.c[i])
                plan<<" "<<keys[r.id[i]]<<" "<<int(r.c[i]);
            plan<<"\n";
        }
        return true;
    };
    std::filesystem::path old=std::filesystem::path(argv[5])/"positive-face-rows-a";
    auto next=std::filesystem::path(argv[5])/"positive-face-rows-b";
    {
        std::ofstream disk(old,std::ios::binary);
        auto generate=[&](size_t id,uint64_t key){
            Powers p;
            for(int i=8;i>=0;--i){p[i]=key&15;key>>=4;}
            for(int v=0;v<(all_finite?5:2);++v){
                int inf=2*q;
                const int first=v<2?2:0;
                for(int j=first;j<=4;++j)if(j!=v)inf-=p[labels[v][j]];
                for(int r=2;r<=4;++r){
                    if(r==v)continue;
                    const int vr=labels[v][r];
                    if(p[vr]==1)continue;
                    Row row;row.add(id,inf-p[vr]);bool ok=true;
                    for(int j=first;j<=4 && ok;++j)if(j!=r && j!=v){
                        int vj=labels[v][j],jr=labels[r][j];
                        if(p[vj]>=q-1||p[jr]==1){ok=false;break;}
                        for(int side=0;side<2;++side){
                        auto h=p;h[vj]++;h[side?jr:vr]--;
                            int neighbor=index.get(finite_key(h));
                            if(neighbor<0){ok=false;break;}
                            row.add(neighbor,side?p[vj]:-p[vj]);
                        }
                    }
                    if(!ok)continue;
                    ++rows;
                    if(!propagate(row))disk.write(reinterpret_cast<const char*>(&row),sizeof(row));
                }
            }
        };
        if(all_finite){
            for(size_t i=0;i<index.keys.size();++i)if(index.keys[i])generate(index.ids[i],index.keys[i]);
        }else for(size_t id=0;id<keys.size();++id)generate(id,keys[id]);
    }
    std::cerr<<(all_finite?"Convergent five- and seven-term rows ":"Convergent five-term rows ")<<rows<<", initial angular classes "<<initial
             <<", substitutions on generation "<<extra<<"\n";
    if(!scaled.empty()){
        std::filesystem::remove(old);
        std::cerr<<"Exact rational row checks in all three columns "<<rows<<"\n";
        return 0;
    }
    for(int pass=1;;++pass){
        uint64_t before=extra,pending=0;
        {
            std::ifstream input(old,std::ios::binary);std::ofstream output(next,std::ios::binary);
            Row row;
            while(input.read(reinterpret_cast<char*>(&row),sizeof(row)))
                if(!propagate(row)){output.write(reinterpret_cast<const char*>(&row),sizeof(row));++pending;}
        }
        std::cerr<<"Closure pass "<<pass<<", new "<<extra-before<<", pending "<<pending<<"\n";
        std::filesystem::remove(old);std::swap(old,next);
        if(extra==before)break;
    }
    std::filesystem::remove(old);
    std::cerr<<"New angular periods "<<extra<<", remaining "<<keys.size()-initial-extra<<"\n";
}
