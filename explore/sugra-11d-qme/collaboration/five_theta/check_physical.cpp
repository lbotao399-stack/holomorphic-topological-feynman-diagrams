// Direct evaluation of the two original cubic Wick words on physical external data.
// No generated numerator, graph tensor, or assembler module is imported.
#include "gamma_options.hpp"
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <sstream>
#include <chrono>
#include <filesystem>
#include <set>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

using V=std::array<int,5>; using Order=V;
#ifdef EXACT_INTEGER
using Scalar=long long;
constexpr int inv6=1;
constexpr char mode[]="integer";
inline Scalar mod(long long x){return x;}
#else
using Scalar=int;
#ifdef SUGRA_PHYSICAL_PRIME
constexpr int P=SUGRA_PHYSICAL_PRIME;
#else
constexpr int P=1009;
#endif
constexpr int inv6=(P%6==1?5*P+1:P+1)/6;
constexpr char mode[]="F1009";
inline Scalar mod(long long x){return (x%P+P)%P;}
#endif
using Matrix=std::array<Scalar,25>; using Tensor=std::array<Scalar,125>;
template<class T> int signature(const T& a){
    int s=1;for(unsigned i=0;i<a.size();++i)for(unsigned j=i+1;j<a.size();++j)if(a[i]>a[j])s=-s;return s;
}
struct Dword{int element,external_axis,scale,first_slot,first_axis,second_slot,second_axis;};
struct Vertex{std::vector<Dword> words; bool cubic;};
using Endpoints=std::array<std::array<int,2>,13>;
Endpoints ends; std::array<std::vector<int>,3> vertices, edgeids;
void topology(std::array<int,3> lengths){int e=0,next=1;for(int p=0;p<3;++p){int left=0;for(int j=0;j<lengths[p];++j){int right=j+1==lengths[p]?11:next++;ends[e]={left,right};edgeids[p].push_back(e++);if(right!=11)vertices[p].push_back(right);left=right;}}}
struct Flow{std::array<const Vertex*,12> local;int wick;};

// Local slot ids are sorted edge numbers followed by external-slot sentinel 13.
Vertex make_vertex(int v,const std::vector<int>& incid,const Endpoints& field,bool cubic){
    auto slots=incid;if(slots.size()==2)slots.push_back(13);
    Vertex result;result.cubic=cubic;
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
// Cache the literal words for each local edge orientation. Different gamma supports
// share these vertices; only four small orientation records are needed per support.
struct FlowPlans {
    std::map<std::pair<int,int>,Vertex> vertices;
    std::array<std::vector<Flow>,1024> by_sites;
    FlowPlans(const FlowPlans&)=delete;
    FlowPlans& operator=(const FlowPlans&)=delete;
    explicit FlowPlans(int gamma_count) {
        std::array<std::vector<int>,12> incident;
        for(int e=0;e<13;++e)for(int v:ends[e])incident[v].push_back(e);
        for(int mask=0;mask<8192;++mask) {
            Endpoints field=ends;
            std::array<int,12> incoming{};
            std::array<int,13> tails;
            for(int e=0;e<13;++e) {
                if(mask>>e&1)std::swap(field[e][0],field[e][1]);
                tails[e]=field[e][0];++incoming[field[e][1]];
            }
            if(!(incoming[0]&1) || !(incoming[11]&1))continue;
            unsigned sites=0;int cubics=0;
            for(int v=0;v<12;++v) {
                if(incoming[v]==int(incident[v].size()))++cubics;
                if(incident[v].size()==2 && incoming[v]!=1)sites|=1U<<(v-1);
            }
            if(__builtin_popcount(sites)!=gamma_count)continue;
            if(2*cubics!=1+gamma_count)throw std::logic_error("wrong coupling order");
            Flow flow;flow.wick=signature(tails);
            for(int v=0;v<12;++v) {
                int local_mask=0;
                for(int e:incident[v]) local_mask=2*local_mask+int(field[e][1]==v);
                const auto key=std::make_pair(v,local_mask);
                auto found=vertices.find(key);
                if(found==vertices.end())found=vertices.emplace(key,make_vertex(v,incident[v],field,incoming[v]==int(incident[v].size()))).first;
                flow.local[v]=&found->second;
            }
            by_sites[sites].push_back(flow);
        }
        for(unsigned sites=0;sites<1024;++sites)
            if(__builtin_popcount(sites)==gamma_count && by_sites[sites].size()!=4)throw std::runtime_error("orientation count");
    }
};
Matrix chain(const std::array<Tensor,12>& tv,int path){
    Matrix mat{},next{};
    for(int i=0;i<5;++i)mat[5*i+i]=1;
    for(int v:vertices[path]){
        for(int i=0;i<5;++i)for(int j=0;j<5;++j){
            Scalar sum=0;for(int k=0;k<5;++k)sum+=mat[5*i+k]*tv[v][5*k+j];
            next[5*i+j]=mod(sum);
        }
        mat=next;
    }
    return mat;
}
Scalar raw_numerator(const std::vector<Flow>& flows,const std::array<V,13>& y,
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
            for(const auto& w:flow.local[v]->words){
                if(w.external_axis!=-1 && w.external_axis!=u[v])continue;
                int t=signed_y[v][w.first_slot][w.first_axis];
                if(w.second_slot>=0)t*=signed_y[v][w.second_slot][w.second_axis];
                // Input currents are small signed integers. Only accumulated vertex entries use residues.
                tv[v][w.element]+=w.scale*t;
            }
            for(auto& a:tv[v]) {
#ifdef EXACT_INTEGER
                if(flow.local[v]->cubic) {
                    if(a%6)throw std::logic_error("cubic word normalization");
                    a/=6;
                }
#endif
                a=mod(a);
            }
        }
        Matrix A=chain(tv,0),B=chain(tv,1),C=chain(tv,2);
        // Apply one path matrix per tensor mode. Reduce before multiplying the next matrix.
        Tensor T{},U{};
        for(int b=0;b<5;++b)for(int c=0;c<5;++c)for(int i=0;i<5;++i){
            Scalar s=0;for(int a=0;a<5;++a)s+=tv[0][25*a+5*b+c]*A[5*a+i];T[25*i+5*b+c]=mod(s);}
        for(int i=0;i<5;++i)for(int c=0;c<5;++c)for(int j=0;j<5;++j){
            Scalar s=0;for(int b=0;b<5;++b)s+=T[25*i+5*b+c]*B[5*b+j];U[25*i+5*j+c]=mod(s);}
        long long contraction=0;
        for(int i=0;i<5;++i)for(int j=0;j<5;++j)for(int k=0;k<5;++k){
            Scalar s=0;for(int c=0;c<5;++c)s+=U[25*i+5*j+c]*C[5*c+k];
            contraction+=static_cast<long long>(mod(s))*tv[11][25*i+5*j+k];}
        total+=flow.wick*mod(contraction);
    }
    return mod(total);
}
#ifndef SUGRA_PHYSICAL_WORDS_ONLY
using Word=std::array<int,10>;
std::vector<Word> pairings(){
    // Enumerate perfect matchings by pairing the first unpaired slot with each later slot.
    Word assigned{};assigned.fill(-1);std::vector<Word> words;
    std::array<int,12> path{};for(int p=0;p<3;++p)for(int v:vertices[p])path[v]=p;
    auto rec=[&](auto&& self,int label)->void{
        if(label==5){words.push_back(assigned);return;}
        int a=0;while(assigned[a]>=0)++a;assigned[a]=label;
        for(int b=a+1;b<10;++b)if(assigned[b]<0 && path[a+1]!=path[b+1]){
            assigned[b]=label;self(self,label+1);assigned[b]=-1;
        }
        assigned[a]=-1;
    };rec(rec,0);std::sort(words.begin(),words.end());return words;
}
// Integer determinant by grouping the columns by physical spatial coordinate.
// Each 2x2 block consists of the corresponding two holomorphic momenta.
int physical_determinant(const Word& moment,const Word& path){
    std::array<int,10> togroup{};int cursor=0,det=1;
    for(int axis=1;axis<=5;++axis){
        int lo=-1,hi=-1;for(int v=0;v<10;++v)if(std::abs(moment[v])==axis){if(lo<0)lo=v;else hi=v;}
        togroup[cursor++]=lo;togroup[cursor++]=hi;
        int p=path[lo],q=path[hi];int a=(moment[lo]>0?1:-1),b=(moment[hi]>0?1:-1);
        det*=a*b*((p!=2)*(q!=0)-(p!=0)*(q!=2));
    }
    return det*signature(togroup); // Interleaving the two blocks of 5 rows has 10 inversions.
}
std::array<Scalar,2> evaluate(const Word& pattern,const V& b,const V& c,const FlowPlans& plans,const GammaOptions& options){
    // Fill the pair of physical momentum slots axis-by-axis. The positive eigenwave and
    // its negative partner share an axis, have independent physical polarizations and
    // carry fixed physical Grassmann labels. This loop does not act on any tensor index.
    const int labels[11] = {9,7,8,6,5,-9,0,1,2,3,4}; // momentum m indexes [m+5]
    std::array<std::array<int,2>,5> pair{};std::array<int,5> cnt{};
    for(int v=0;v<10;++v)pair[pattern[v]][cnt[pattern[v]]++]=v;
    Word path{};for(int p=0;p<3;++p)for(int v:vertices[p])path[v-1]=p;
    std::array<int,5> choice={0,1,2,3,4};Scalar sum=0,wrong=0; int count=0;
    do{
        for(int signs=0;signs<32;++signs){
            Word waves{};
            std::array<int,12> u{};std::array<V,12> mom{};
            std::array<V,13> currents{};
            for(int axis=0;axis<5;++axis){
                int j=choice[axis],s=(signs>>axis&1)?-1:1;
                int plus=pair[j][s==1?0:1],minus=pair[j][s==1?1:0];
                waves[plus]=axis+1;waves[minus]=-axis-1;
                mom[plus+1][axis]=1;mom[minus+1][axis]=-1;
                u[plus+1]=options.polarizations[labels[axis+6]];u[minus+1]=options.polarizations[labels[4-axis]];
                currents[edgeids[1][0]][axis]=s*b[j];currents[edgeids[2][0]][axis]=s*c[j];
                currents[edgeids[0][0]][axis]=-s*(b[j]+c[j]);
            }
            std::vector<int> odd;unsigned sites=0;
            for(int v=0;v<10;++v) {
                int label=labels[waves[v]+5];
                if((options.mask>>label)&1U)sites|=1U<<v;else odd.push_back(label);
            }
            for(int p=0;p<3;++p)for(unsigned j=0;j<vertices[p].size();++j)
                for(int i=0;i<5;++i)currents[edgeids[p][j+1]][i]=currents[edgeids[p][j]][i]+mom[vertices[p][j]][i];
            int d=physical_determinant(waves,path);if(d!=1 && d!=-1)throw std::logic_error("invalid Dtheta");
            Scalar term=d*raw_numerator(plans.by_sites[sites],currents,u,mom);
            sum+=signature(odd)*term;wrong+=term;++count;
        }
    }while(std::next_permutation(choice.begin(),choice.end()));
    if(count!=3840)throw std::logic_error("orbits");return {mod(sum),mod(wrong)};
}
int main(int argc,char** argv){
    GammaOptions options;
    if(!options.parse(argc,argv) || options.positional.size()<5) {
        std::cerr<<"Usage: check_physical p q r seed polynomial... [--gamma-mask bits] [--polarizations i1,...,i10]\n";return 2;
    }
    const auto& args=options.positional;
    topology({std::stoi(args[0]),std::stoi(args[1]),std::stoi(args[2])});
    int seed=std::stoi(args[3]);V b{},c{};std::array<int,10> xy{};
    for(int i=0;i<5;++i){xy[i]=b[i]=1+((i+seed)*(i+2)+seed)%17;xy[5+i]=c[i]=-1-((3*i+seed)*(i+5)+2*seed)%19;}
#ifdef EXACT_INTEGER
    b={1,-2,3,1,2};c={2,1,-1,-2,3};
    std::copy(b.begin(),b.end(),xy.begin());std::copy(c.begin(),c.end(),xy.begin()+5);
#endif
    auto patterns=pairings();FlowPlans plans(__builtin_popcount(options.mask));
    for(std::size_t arg=4;arg<args.size();++arg){std::ifstream in(args[arg]);if(!in)throw std::runtime_error("missing polynomial");
        int rank,count;int blocks=0,dropped=0,nonzero=0;
        while(in>>rank>>count){Scalar expected=0;for(int n=0;n<count;++n){uint64_t key;long long coefficient;if(!(in>>key>>coefficient))throw std::logic_error("truncated polynomial");
            Scalar v=mod(coefficient);for(int j=0;j<10;++j)for(int k=0;k<int((key>>(4*j))&15);++k)v=mod(v*xy[j]);expected=mod(expected+v);}
            auto actual=evaluate(patterns.at(rank),b,c,plans,options);
            if(expected!=actual[0]){std::cerr<<args[arg]<<" rank "<<rank<<" actual "<<actual[0]<<" expected "<<expected<<"\n";return 1;}
            ++blocks;dropped+=actual[1]!=expected;nonzero+=expected!=0;
        }
        if(!blocks)throw std::logic_error("empty polynomial");
        std::cout<<std::filesystem::path(args[arg]).filename().string()<<' '<<blocks<<" blocks x 3840 physical assignments x 4 flows; "<<mode<<" seed "<<seed<<"; "<<nonzero<<" nonzero, no exterior sign fails "<<dropped<<'\n'<<std::flush;
    }
}
#endif
