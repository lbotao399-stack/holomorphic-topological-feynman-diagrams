// Direct evaluation of the two original cubic Wick words on physical external data.
// No generated numerator, graph tensor, or assembler module is imported.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

using V=std::array<int,5>; using Order=V;
using Matrix=std::array<int,25>; using Tensor=std::array<int,125>;
constexpr int P=1009,inv6=841; // 6*841 = 1 modulo 1009.
int mod(long long x){return (x%P+P)%P;}
template<class T> int signature(const T& a){
    int s=1;for(unsigned i=0;i<a.size();++i)for(unsigned j=i+1;j<a.size();++j)if(a[i]>a[j])s=-s;return s;
}
struct Dword{int element,external_axis,scale,first_slot,first_axis,second_slot,second_axis;};
struct Vertex{std::vector<Dword> words;};
using Endpoints=std::array<std::array<int,2>,13>;
const Endpoints ends={{{0,11},{0,1},{1,2},{2,3},{3,4},{4,5},{5,11},{0,6},{6,7},{7,8},{8,9},{9,10},{10,11}}};
struct Flow{std::array<Vertex,12> local;int wick;};

// Local slot ids are sorted edge numbers followed by external-slot sentinel 13.
Vertex make_vertex(int v,const std::vector<int>& incid,const Endpoints& field,bool cubic){
    auto slots=incid;if(slots.size()==2)slots.push_back(13);
    Vertex result;
    std::array<int,3> perm={0,1,2};
    do {
        // perm maps letter position in the source word to local fiber position.
        if(!cubic){
            bool admissible=true;
            for(int l=0;l<3;++l){
                int e=slots[perm[l]];
                if(e!=13 && field[e][l<2?0:1]!=v)admissible=false;
            }
            if(!admissible)continue;
            int mu0=slots[perm[0]],mu1=slots[perm[1]];
            int wick=(mu1!=13 && (mu0==13 || mu0<mu1))?-1:1;
            for(int a=0;a<5;++a)for(int b=0;b<5;++b){
                std::array<int,3> axes{};
                axes[perm[0]]=a;axes[perm[1]]=b;axes[perm[2]]=b;
                int ix=incid.size()==3?25*axes[0]+5*axes[1]+axes[2]:5*axes[0]+axes[1];
                result.words.push_back({ix,incid.size()==3?-1:axes[2],wick,slots[perm[2]],a,-1,0});
            }
        } else {
            Order ie={0,1,2,3,4};
            do{
                std::array<int,3> axes{};
                axes[perm[0]]=ie[0];axes[perm[1]]=ie[2];axes[perm[2]]=ie[4];
                int ix=incid.size()==3?25*axes[0]+5*axes[1]+axes[2]:5*axes[0]+axes[1];
                result.words.push_back({ix,incid.size()==3?-1:axes[2],signature(ie)*inv6,
                                        slots[perm[1]],ie[1],slots[perm[2]],ie[3]});
            }while(std::next_permutation(ie.begin(),ie.end()));
        }
    }while(std::next_permutation(perm.begin(),perm.end()));
    return result;
}
std::vector<Flow> make_flows(int gamma){
    std::vector<Flow> out;
    for(int mask=0;mask<8192;++mask){
        Endpoints f=ends;
        std::array<int,12> incoming{},degree{};
        std::array<std::vector<int>,12> incid;
        std::array<int,13> tailorder;
        for(int e=0;e<13;++e){
            if(mask>>e&1)std::swap(f[e][0],f[e][1]);
            tailorder[e]=f[e][0]; ++incoming[f[e][1]];
            for(int v:f[e]){++degree[v];incid[v].push_back(e);}
        }
        bool valid=true; int cubics=0, gamma_count=0;
        for(int v=0;v<12;++v){
            if(degree[v]==3 && incoming[v]!=1 && incoming[v]!=3)valid=false;
            if(incoming[v]==degree[v])++cubics;
            if(degree[v]==2 && incoming[v]!=1){++gamma_count;if(v!=gamma)valid=false;}
        }
        if(!valid || cubics!=1 || gamma_count!=1)continue;
        Flow plan;plan.wick=signature(tailorder);
        for(int v=0;v<12;++v)plan.local[v]=make_vertex(v,incid[v],f,incoming[v]==degree[v]);
        out.push_back(std::move(plan));
    }
    if(out.size()!=4)throw std::runtime_error("orientation count");
    return out;
}
Matrix chain(const std::array<Tensor,12>& tv,int offset){
    Matrix mat{},next{};
    for(int i=0;i<5;++i)mat[5*i+i]=1;
    for(int v=offset;v<offset+5;++v){
        for(int i=0;i<5;++i)for(int j=0;j<5;++j){
            int sum=0;for(int k=0;k<5;++k)sum+=mat[5*i+k]*tv[v][5*k+j];
            next[5*i+j]=sum%P;
        }
        mat=next;
    }
    return mat;
}
int raw_numerator(const std::vector<Flow>& flows,const std::array<V,13>& y,
                  const std::array<int,12>& u,const std::array<V,12>& p){
    std::array<std::array<V,14>,12> signed_y{};
    for(int v=0;v<12;++v){
        signed_y[v][13]=p[v];
        for(int e=0;e<13;++e)for(int i=0;i<5;++i)
            signed_y[v][e][i]=(ends[e][1]==v?1:-1)*y[e][i];
    }
    long long total=0;
    for(auto& flow:flows){
        std::array<Tensor,12> tv{};
        for(int v=0;v<12;++v){
            for(const auto& w:flow.local[v].words){
                if(w.external_axis!=-1 && w.external_axis!=u[v])continue;
                int t=signed_y[v][w.first_slot][w.first_axis];
                if(w.second_slot>=0)t*=signed_y[v][w.second_slot][w.second_axis];
                // Input currents are small signed integers. Only accumulated vertex entries use residues.
                tv[v][w.element]+=w.scale*t;
            }
            for(auto& a:tv[v])a=mod(a);
        }
        Matrix b=chain(tv,1),c=chain(tv,6);
        long long contraction=0;
        for(int a=0;a<5;++a)for(int i=0;i<5;++i)for(int j=0;j<5;++j){
            int t=tv[0][25*a+5*i+j];if(t==0)continue;
            for(int k=0;k<5;++k)for(int l=0;l<5;++l)
                contraction+=(long long)t*tv[11][25*a+5*k+l]*b[5*i+k]*c[5*j+l];
        }
        total+=flow.wick*mod(contraction);
    }
    return mod(total);
}
int main(int argc,char** argv){
    if(argc!=5){std::cerr<<"need gamma index (1..3), polynomial path, seed, max relative order (-1 for all)\n";return 2;}
    int gamma=std::stoi(argv[1]), seed=std::stoi(argv[3]),maxrank=std::stoi(argv[4]);
    auto flows=make_flows(gamma);
    V b{},c{};std::array<int,10> xy{};
    for(int i=0;i<5;++i){xy[i]=b[i]=1+((i+seed)*(i+2)+seed)%17;xy[5+i]=c[i]=-1-((3*i+seed)*(i+5)+2*seed)%19;}
    Order sorted={0,1,2,3,4};std::map<Order,int> rank;int r=0;
    do{rank[sorted]=r++;}while(std::next_permutation(sorted.begin(),sorted.end()));
    std::array<int,120> reduced{},wrong{},valid{},expect{};std::array<bool,120> seen{};
    std::array<int,9> odd={0,1,2,3,4,5,6,7,8};
    // Labels and polarizations in fixed exterior input order.
    const std::array<int,9> moment={2,3,4,5,-1,-2,-4,-3,-5};
    const std::array<int,9> polar={2,3,4,1,1,0,0,0,0};
    do{
        int cursor=0;
        std::array<int,10>waves{},axis{};
        for(int v=0;v<10;++v){
            if(v==gamma-1){waves[v]=1;axis[v]=1;}
            else{waves[v]=moment[odd[cursor]];axis[v]=polar[odd[cursor++]];}
        }
        V bp{},cp{},where{};int md=1;
        int occupied=0;
        for(int i=0;i<5;++i){bp[i]=std::abs(waves[i])-1;cp[i]=std::abs(waves[5+i])-1;where[bp[i]]=i;
            occupied|=1<<bp[i];md*=waves[i]>0?1:-1;md*=waves[5+i]>0?1:-1;}
        if(occupied!=31)continue;
        md*=signature(bp)*signature(cp);
        Order rel{};for(int i=0;i<5;++i)rel[i]=where[cp[i]];
        int rr=rank.at(rel);if(maxrank!=-1 && rr>maxrank)continue;
        std::array<V,12> p{};std::array<int,12> u{};
        for(int i=0;i<10;++i){p[i+1][std::abs(waves[i])-1]=waves[i]>0?1:-1;u[i+1]=axis[i];}
        std::array<V,13> y{};
        for(int i=0;i<5;++i){
            int ps=waves[i]>0?1:-1;
            y[1][bp[i]]=ps*b[i];y[7][bp[i]]=ps*c[i];
            y[0][bp[i]]=-ps*(b[i]+c[i]);
        }
        for(int i=0;i<5;++i)for(int j=0;j<5;++j){
            y[2+i][j]=y[1+i][j]+p[1+i][j];
            y[8+i][j]=y[7+i][j]+p[6+i][j];
        }
        int value=raw_numerator(flows,y,u,p);
        reduced[rr]=mod(reduced[rr]+signature(odd)*md*value);
        wrong[rr]=mod(wrong[rr]+md*value);
        ++valid[rr];
    }while(std::next_permutation(odd.begin(),odd.end()));
    std::ifstream in(argv[2]);if(!in)return 2;
    int block,nterms;
    while(in>>block>>nterms){
        if(seen[block])throw std::runtime_error("duplicate polynomial block");
        seen[block]=true;
        for(int k=0;k<nterms;++k){
            uint64_t key;long long coeff;if(!(in>>key>>coeff))throw std::runtime_error("truncated block");
            int val=mod(coeff);
            for(int i=0;i<10;++i){int power=(key>>(4*i))&15;while(power--)val=mod(val*xy[i]);}
            expect[block]=mod(expect[block]+val);
        }
    }
    int blocks=0,wrong_fail=0,nonzero=0;
    for(int rr=0;rr<120;++rr)if(valid[rr]){
        if(!seen[rr] || valid[rr]!=384 || reduced[rr]!=expect[rr]){
            std::cerr<<"mismatch slot "<<gamma<<" rank "<<rr<<" count "<<valid[rr]<<" direct "<<reduced[rr]<<" polynomial "<<expect[rr]<<"\n";return 1;}
        ++blocks;wrong_fail+=wrong[rr]!=expect[rr];nonzero+=expect[rr]!=0;
    }
    std::cout<<"B"<<gamma<<" point "<<seed<<": "<<blocks<<" blocks x 384 assignments x 4 Wick orientations; "<<nonzero<<" nonzero blocks; dropping exterior sign fails "<<wrong_fail<<" blocks, F_1009\n";
}
