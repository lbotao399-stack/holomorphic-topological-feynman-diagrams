// Sparse rank certificate for the finite positive-orbit equations with prescribed wheel faces.
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

constexpr int prime=1000000007;
using Term=std::pair<int,int>;
using Row=std::vector<Term>;

int inverse(int n){
    int val=1;
    for(int k=prime-2;k;k>>=1){if(k&1)val=int(int64_t(val)*n%prime);n=int(int64_t(n)*n%prime);}
    return val;
}
void normalize(Row& r){
    int factor=inverse(r.back().second);
    for(auto& t:r)t.second=int(int64_t(t.second)*factor%prime);
}
bool eliminate(Row& row,const Row& pivot,int cap){
    Row v;v.reserve(std::min(cap+1,int(row.size()+pivot.size())));
    unsigned i=0,j=0;
    int scale=prime-row.back().second;
    while(i<row.size()||j<pivot.size()){
        if(j==pivot.size()||(i<row.size()&&row[i].first<pivot[j].first))v.push_back(row[i++]);
        else if(i==row.size()||row[i].first>pivot[j].first){
            v.emplace_back(pivot[j].first,int(int64_t(scale)*pivot[j].second%prime));++j;
        }else{
            int c=int((row[i].second+int64_t(scale)*pivot[j].second)%prime);
            if(c)v.emplace_back(row[i].first,c);
            ++i;++j;
        }
        if(int(v.size())>cap)return false;
    }
    row.swap(v);return true;
}
template<class F> void scan(const std::string& base,F visit){
    std::ifstream f(base+".bin",std::ios::binary);
    if(!f)throw std::runtime_error("equation file missing");
    Row row;uint8_t n;
    while(f.read(reinterpret_cast<char*>(&n),1)){
        row.clear();
        for(int j=0;j<n;++j){
            int id;int8_t c;
            f.read(reinterpret_cast<char*>(&id),4);f.read(reinterpret_cast<char*>(&c),1);
            if(id>=0)row.emplace_back(id,c>0?c:prime+c);
        }
        std::sort(row.begin(),row.end());
        if(!row.empty())visit(row);
    }
}
int main(int argc,char**argv){
    if(argc<3)throw std::invalid_argument("equation basename and support limits required");
    std::string base=argv[1];
    std::ifstream nodes(base+".nodes");
    uint64_t key;int n=0;std::vector<uint64_t> keys;
    while(nodes>>key){++n;keys.push_back(key);}
    std::vector<Row> pivots(n);int rank=0;
    scan(base,[&](Row& r){
        auto& v=pivots[r.back().first];
        if(!v.empty() && v.size()<=r.size())return;
        if(v.empty())++rank;
        v=r;normalize(v);
    });
    int start=2,target=n-1;
    if(std::string(argv[2])=="--equal" || std::string(argv[2])=="--low"){
        std::string mode=argv[2];int q=std::stoi(argv[3]);
        for(int i=0;i<n;++i){
            uint64_t a=keys[i];int sums[6]{};
            for(int j=8;j>=0;--j){sums[j/3]+=a&15;sums[3+j%3]+=a&15;a>>=4;}
            if(std::any_of(sums,sums+6,[&](int s){return mode=="--low"?s<=q:s==q;})){
                if(pivots[i].empty())++rank;
                pivots[i]={{i,1}};
            }
        }
        start=4;target=n;
    }
    std::cerr<<"initial rank "<<rank<<" missing "<<n-rank<<std::endl;
    for(int pass=start;pass<argc && rank<target;++pass){
        int cap=std::stoi(argv[pass]);
        scan(base,[&](Row& r){
            while(!r.empty()){
                auto& v=pivots[r.back().first];
                if(v.empty()){v=r;normalize(v);++rank;break;}
                if(r.size()<v.size()){normalize(r);r.swap(v);}
                if(!eliminate(r,v,cap))break;
            }
        });
        std::cerr<<"support cap "<<cap<<" rank "<<rank<<" missing "<<n-rank<<std::endl;
    }
    int free_count=n-rank;
    if(free_count>0 && free_count<=128 && rank<target){
        std::vector<int> evaluations(size_t(n)*free_count);
        int at=0;
        for(int i=0;i<n;++i){
            int* dest=evaluations.data()+size_t(i)*free_count;
            const auto& current=pivots[i];
            if(current.empty()){dest[at++]=1;continue;}
            for(int k=0;k<free_count;++k){
                __uint128_t sum=0;
                for(unsigned j=0;j+1<current.size();++j){
                    auto [idx,c]=current[j];
                    sum+=uint64_t(c)*evaluations[size_t(idx)*free_count+k];
                }
                dest[k]=(prime-uint64_t(sum%prime))%prime;
            }
        }
        std::vector<std::vector<int>> last(free_count);
        int extra=0;
        scan(base,[&](Row& row){
            if(rank+extra>=target)return;
            std::vector<int> v(free_count);
            for(int k=0;k<free_count;++k){
                uint64_t sum=0;
                for(auto [i,c]:row)sum+=uint64_t(c)*evaluations[size_t(i)*free_count+k];
                v[k]=sum%prime;
            }
            for(int k=0;k<free_count;++k)if(v[k]){
                if(last[k].empty()){
                    int inv=inverse(v[k]);for(int j=k;j<free_count;++j)v[j]=int(int64_t(v[j])*inv%prime);
                    last[k]=std::move(v);++extra;break;
                }
                int c=prime-v[k];for(int j=k;j<free_count;++j)
                    v[j]=int((int64_t(c)*last[k][j]+v[j])%prime);
            }
        });
        rank+=extra;
        std::cerr<<"free-coordinate residual rank "<<extra<<" remaining "<<n-rank<<std::endl;
    }
    size_t width=0,total=0;for(const auto& p:pivots){width=std::max(width,p.size());total+=p.size();}
    std::cout<<"positive variables "<<n<<" prime "<<prime<<" rank lower bound "<<rank
             <<" nonzeros "<<total<<" max row "<<width<<std::endl;
}
