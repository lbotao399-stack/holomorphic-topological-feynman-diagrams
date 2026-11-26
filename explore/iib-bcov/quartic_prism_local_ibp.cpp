// Reduce a finite set of completed cuts to scalar and numerator wheels with
// strictly convergent position-space integration-by-parts identities.
#include "quartic_prism_wheel_charts.hpp"
#include <gmpxx.h>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <unordered_map>
#include <set>
#include <fstream>

struct Hash {
    std::size_t operator()(const Key &key) const {
        std::size_t hash=1;
        for(int power : key) hash=31*hash+unsigned(power+32);
        return hash;
    }
};
constexpr int unresolved_degree=100;
using Ratio=mpq_class;
// factor is the period of the seed divided by the period of this chart.
struct Entry{int idx;Ratio factor;};
struct Component{std::vector<Key> charts;int degree;uint64_t wheel;Ratio factor;};
std::unordered_map<Key,Entry,Hash> labels;
std::vector<Component> nodes;
std::vector<mpz_class> fact;
std::vector<uint64_t> known_prisms;
std::array<std::vector<uint64_t>,7> known_wheels;
int q;

int register_orbit(const Key &seed) {
    if(auto it=labels.find(seed);it!=labels.end())return it->second.idx;
    int idx=nodes.size();
    labels[seed]={idx,1};
    std::vector<Key> queue{seed};
    int degree=unresolved_degree;
    uint64_t wheel=0;
    Ratio terminal=0;
    bool wheel_cached=false;
    for(size_t i=0;i<queue.size();++i){
        auto a=unpack(queue[i]);
        Ratio f=labels.at(queue[i]).factor;
        if(degree>=0 && !known_prisms.empty()) {
            auto code=canonical(a,q);
            if(code!=std::numeric_limits<uint64_t>::max()) {
                uint64_t candidate=0;
                for(auto exponent:decode(code,q))candidate=(candidate<<4)|exponent;
                if(std::binary_search(known_prisms.begin(),known_prisms.end(),candidate)){
                    degree=-1;wheel=candidate;terminal=f;
                }
            }
        }
        if(degree>=0 && !wheel_cached){
            for(int d=0;d<std::min(7,q);++d){
                if(known_wheels[d].empty())continue;
                auto w=numerator_wheel_key(a,q,d,false,&known_wheels[d]);
                if(w){degree=d;wheel=w;terminal=f;wheel_cached=true;break;}
            }
        }
        for(int d=0;!wheel_cached && d<std::min({degree,7,q});++d) {
            auto w=numerator_wheel_key(a,q,d,false);
            if(w){degree=d;wheel=w;terminal=f;break;}
        }
        // Four-point reflection preserves the period up to its gamma factor.
        for(int v=0;v<6;++v){
            std::vector<int> ns;
            for(int j=0;j<6;++j)if(a[v][j])ns.push_back(j);
            if(ns.size()!=4 || std::any_of(ns.begin(),ns.end(),[&](int u){return a[v][u]<1||a[v][u]>=q;}))continue;
            auto b=a; Ratio factor=f;
            for(int j:ns){
                int t=a[v][j];
                b[v][j]=b[j][v]=q-t;
                factor*=Ratio(fact[q-t-1],fact[t-1]);
                for(int l:ns)if(l>j)b[j][l]=b[l][j]=a[j][l]+t+a[v][l]-q;
            }
            factor.canonicalize();
            auto k=graph_key(b);auto [it,newkey]=labels.emplace(k,Entry{idx,factor});
            if(newkey)queue.push_back(k);
            else if(it->second.idx!=idx||it->second.factor!=factor)throw std::logic_error("bad orbit");
        }
    }
    nodes.push_back({std::move(queue),degree,wheel,terminal});
    return idx;
}

void shift(Edges &graph,int i,int j,int increment) {
    graph[i][j]+=increment;
    graph[j][i]+=increment;
}
using Equation=std::map<int,Ratio>;

// Integrate the divergence at v with reference point ref after sending inf to infinity.
bool row(const Key &key,int v,int ref,int inf,Equation &out) {
    auto a=unpack(key);
    int c=a[v][inf]-a[v][ref];
    std::vector<std::pair<Key,int>> terms{{key,c}};
    for(int j=0;j<6;++j)if(j!=v&&j!=ref&&j!=inf&&a[v][j]){
        for(int side=0;side<=1;++side){
            auto h=a;shift(h,v,j,1);shift(h,ref,inf,1);
            shift(h,side?j:v,ref,-1);shift(h,side?v:j,inf,-1);
            if(!convergent(h,q))return false;
            terms.push_back({graph_key(h),side?a[v][j]:-a[v][j]});
        }
    }
    for(auto [k,c]:terms){
        int id=register_orbit(k);out[id]+=c/labels.at(k).factor;
        if(!out[id])out.erase(id);
    }
    return true;
}

int main(int argc,char **argv){
    if(argc!=3 && argc!=4)throw std::invalid_argument("q, output file and optional sorted prism keys required");
    q=std::stoi(argv[1]);
    if(q<3 || q>13)throw std::invalid_argument("unsupported wheel dimension");
    const bool strict=std::string(argv[2])=="--canonical-convergent";
    if(strict || std::string(argv[2])=="--canonical-powers"){
        // Regrounding a completed conformal period only permutes its vertices.
        std::ios::sync_with_stdio(false);
        std::cin.tie(nullptr);
        uint64_t input;
        while(std::cin>>input){
            Powers powers;
            for(int k=8;k>=0;--k){powers[k]=input&15;input>>=4;}
            auto graph=complete(powers,q);
            if(strict){
                int total=0;
                for(int p:powers)total+=p;
                if(input || total!=4*q || *std::min_element(powers.begin(),powers.end())<1 ||
                    !convergent(graph,q)){
                    std::cerr<<"Outside the positive convergent prism domain:";
                    for(int p:powers)std::cerr<<" "<<p;
                    std::cerr<<"\n";
                    return 1;
                }
            }
            auto canonical_key=canonical(graph,q);
            if(canonical_key==std::numeric_limits<uint64_t>::max())
                throw std::logic_error("no positive prism chart");
            uint64_t answer=0;
            for(auto exponent:decode(canonical_key,q))answer=(answer<<4)|exponent;
            std::cout<<answer<<"\n";
        }
        return 0;
    }
    fact.assign(q,1);
    for(int i=1;i<q;++i)fact[i]=fact[i-1]*i;
    if(argc==4){
        std::ifstream known(argv[3],std::ios::binary|std::ios::ate);
        if(!known)throw std::runtime_error("cannot open scalar prism key set");
        const auto size=known.tellg();
        if(size<0 || size%sizeof(uint64_t))throw std::runtime_error("invalid scalar prism key set");
        known_prisms.resize(size/sizeof(uint64_t));known.seekg(0);
        known.read(reinterpret_cast<char*>(known_prisms.data()),size);
        if(!known || !std::is_sorted(known_prisms.begin(),known_prisms.end()))
            throw std::runtime_error("invalid order or truncation in scalar prism keys");
        std::ifstream wheels(std::string(argv[3])+".wheels",std::ios::binary);
        uint64_t item;
        while(wheels.read(reinterpret_cast<char*>(&item),sizeof(item))){
            const auto d=item>>32;
            if(d>=known_wheels.size())throw std::runtime_error("invalid numerator-wheel degree");
            known_wheels[d].push_back(item&0xffffffff);
        }
        for(const auto &keys:known_wheels)
            if(!std::is_sorted(keys.begin(),keys.end()))throw std::runtime_error("unordered numerator-wheel keys");
    }
    std::vector<int> target;
    std::set<int> distinct;
    std::vector<Entry> original;
    Powers p;
    while(std::cin>>p[0]){
        for(int j=1;j<9;++j)if(!(std::cin>>p[j]))throw std::runtime_error("truncated seed");
        auto k=graph_key(complete(p,q));
        int id=register_orbit(k);original.push_back(labels.at(k));
        if(nodes[id].degree==unresolved_degree && distinct.insert(id).second)target.push_back(id);
    }
    struct Rule {Key source; unsigned char v,a,b; std::vector<int> unknown;};
    std::vector<Rule> rows;std::map<int,std::vector<int>> reverse;
    std::vector<int> todo;
    std::vector<int> unresolved,order;
    std::set<int> examined;
    std::map<int,int> pivots;
    auto generate=[&](int id) {
        if(!examined.insert(id).second)return;
        // Registering an IBP neighbor can reallocate nodes.
        const auto charts=nodes[id].charts;
        for(auto &k:charts)for(int v=0;v<6;++v)for(int a=2;a<6;++a)for(int b=a+1;b<6;++b){
            if(a==v||b==v)continue;
            Equation r;if(!row(k,v,a,b,r) || r.empty())continue;
            std::vector<int> nonterm;
            for(const auto &[j,z]:r)
                if(nodes[j].degree==unresolved_degree && !pivots.count(j))nonterm.push_back(j);
            if(nonterm.empty())continue;
            int index=rows.size();
            for(int j:nonterm)reverse[j].push_back(index);
            if(nonterm.size()==1)todo.push_back(index);
            unresolved.push_back(nonterm.size());
            rows.push_back({k,(unsigned char)v,(unsigned char)a,(unsigned char)b,std::move(nonterm)});
        }
    };
    for(int id:target)generate(id);
    size_t cursor=0;
    while(true) {
        // Unit pivots can unlock older equations as well as newly generated ones.
        while(cursor<todo.size()){
            int ix=todo[cursor++];if(unresolved[ix]!=1)continue;
            int id=-1;
            for(int j:rows[ix].unknown)if(!pivots.count(j)){id=j;break;}
            if(id<0)throw std::logic_error("lost pivot");
            pivots[id]=ix;order.push_back(id);
            for(int next:reverse[id])if(--unresolved[next]==1)todo.push_back(next);
        }
        // Expand the nearest unevaluated neighborhoods of surviving seed periods.
        std::set<int> visited,frontier,extension;
        for(int id:target)if(!pivots.count(id))frontier.insert(id);
        while(!frontier.empty() && extension.empty()){
            std::set<int> next;
            for(int i:frontier)if(visited.insert(i).second){
                if(!examined.count(i)){extension.insert(i);continue;}
                for(int ix:reverse[i])for(int j:rows[ix].unknown)
                    if(!pivots.count(j) && !visited.count(j))next.insert(j);
            }
            frontier.swap(next);
        }
        if(extension.empty())break;
        std::cout<<"additional integration-by-parts neighborhoods "<<extension.size()<<std::endl;
        for(int i:extension)generate(i);
    }
    std::set<int> needed;
    for(int id:target)if(pivots.count(id))needed.insert(id);
    int count=needed.size();
    std::ofstream f(argv[2]);
    if(!f)throw std::runtime_error("cannot write reduction");
    std::set<int> terminal;
    for(size_t i=0;i<original.size();++i){
        auto entry=original[i];
        if(nodes[entry.idx].degree!=unresolved_degree)terminal.insert(entry.idx);
        else if(!pivots.count(entry.idx))continue;
        f<<"I "<<i<<" "<<entry.idx<<" "<<1/entry.factor<<"\n";
    }
    // Follow dependencies backwards so the evaluator needs only reachable leaves.
    for(auto it=order.rbegin();it!=order.rend();++it)if(needed.count(*it)){
        int id=*it;
        auto &equation=rows[pivots[id]];
        Equation r;
        if(!row(equation.source,equation.v,equation.a,equation.b,r))
            throw std::logic_error("failed reconstruction of pivot");
        auto c=r.at(id);r.erase(id);
        for(auto [j,z]:r) {
            Ratio coef=-z/c;
            f<<"P "<<id<<" "<<j<<" "<<coef<<"\n";
            if(nodes[j].degree==unresolved_degree)needed.insert(j);
            else terminal.insert(j);
        }
    }
    for(int j:terminal)
        f<<"T "<<j<<" "<<nodes[j].degree<<" "<<nodes[j].wheel<<" "<<nodes[j].factor<<"\n";
    std::cout<<"closed unresolved seed orbits "<<count<<" of "<<target.size()
             <<", necessary pivots "<<needed.size()<<", candidate equations "<<rows.size()
             <<", terminal orbits "<<terminal.size()<<std::endl;
}
