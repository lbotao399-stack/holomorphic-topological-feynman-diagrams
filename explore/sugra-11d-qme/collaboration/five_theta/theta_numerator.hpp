#pragma once

#include <array>
#include <cstdint>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace twoloop {
struct Leg { int momentum; int vector_axis; };
using Lengths = std::array<int,3>;
using Slots = std::array<Leg,10>;
using Polynomial = std::unordered_map<std::uint64_t,std::int64_t>;

namespace detail {
struct Powers {
    std::uint64_t plain = 0, shifted = 0;
    bool operator==(const Powers& p) const { return plain == p.plain && shifted == p.shifted; }
};
struct Hash {
    std::size_t operator()(Powers p) const {
        std::uint64_t v = p.plain + 0x9e3779b97f4a7c15ULL * p.shifted;
        v ^= v >> 30; v *= 0xbf58476d1ce4e5b9ULL; v ^= v >> 27;
        return v;
    }
};
struct Term {
    Powers p;
    int sign = 1;
};
inline Term multiply(Term a, Term b) {
    return {{a.p.plain + b.p.plain, a.p.shifted + b.p.shifted}, a.sign*b.sign};
}
using Entries = std::array<std::array<std::vector<Term>,5>,5>;
using Collected = std::unordered_map<Powers,std::int64_t,Hash>;

inline int epsilon(const std::array<int,5>& a) {
    int sign=1;
    for (int i=0;i<5;++i) for (int j=i+1;j<5;++j) {
        if (a[i]==a[j]) return 0;
        if (a[i]>a[j]) sign=-sign;
    }
    return sign;
}
// Paths have geometric orientation N=0 -> N'=11 and are ordered A,B,C.
// For intermediate products keep b_i,c_i,a_i in that order. Four bits per
// exponent suffice at total degree at most 13. Expand a=-b-c only after Wick summation.
struct Graph {
    const Slots& slots;
    std::array<int,3> first{}, last{};
    std::array<int,13> path{}, tail{}, head{};
    std::array<std::array<int,5>,13> offsets{};
    std::array<int,15> shifts{};
    int gamma, gamma_path=0;
    static int variables(int path) { return ((path+2)%3)*5; }
    Graph(const Lengths& lengths,const Slots& inputs,int gamma_slot): slots(inputs),gamma(gamma_slot+1) {
        if (gamma_slot<0 || gamma_slot>=10 || lengths[0]+lengths[1]+lengths[2]!=13)
            throw std::invalid_argument("invalid theta data");
        int edge=0,next=1;
        for (int p=0;p<3;++p) {
            if (lengths[p]<1 || lengths[p]>6) throw std::invalid_argument("invalid path length");
            first[p]=edge;
            int previous=0;
            std::array<int,5> prefix{};
            for (int j=0;j<lengths[p];++j,++edge) {
                int v=(j+1==lengths[p])?11:next++;
                tail[edge]=previous;head[edge]=v;path[edge]=p;offsets[edge]=prefix;previous=v;
                if (v==11) continue;
                if (v==gamma) gamma_path=p;
                auto leg=slots[v-1];
                if (!leg.momentum || leg.momentum < -5 || leg.momentum>5 || leg.vector_axis<0 || leg.vector_axis>4)
                    throw std::invalid_argument("invalid leg");
                int axis=std::abs(leg.momentum)-1;
                if(prefix[axis]) throw std::invalid_argument("duplicate path momentum axis");
                prefix[axis]=leg.momentum>0?1:-1;
            }
            last[p]=edge-1;
            for (int i=0;i<5;++i) shifts[variables(p)+i]=prefix[i];
        }
        for (int i=0;i<5;++i) if (shifts[i]+shifts[i+5]+shifts[i+10])
            throw std::invalid_argument("nonzero total momentum");
    }
    Term current(int e,int i,int sign) const {
        const auto bit=std::uint64_t{1}<<(4*(variables(path[e])+i));
        return offsets[e][i] ? Term{{0,bit},sign} : Term{{bit,0},sign};
    }
};
struct Directions {
    std::array<bool,13> forward{};
    int sign=1;
    void set_wick_sign(const Graph& g) {
        sign=1;
        for (int i=0;i<13;++i) for (int j=i+1;j<13;++j) {
            int vi=forward[i]?g.tail[i]:g.head[i], vj=forward[j]?g.tail[j]:g.head[j];
            if(vi>vj) sign=-sign;
        }
    }
};
struct Flow : Directions {
    int c_vertex=0;
    Flow(const Graph& g,int which) {
        c_vertex=(which<2)?(which==0?0:11):g.gamma;
        for (int p=0;p<3;++p) {
            bool before_gamma=true;
            for(int e=g.first[p];e<=g.last[p];++e) {
                if(p==g.gamma_path) forward[e]=(which<2) ? !before_gamma : before_gamma;
                else if(which<2) forward[e]=(which==1);
                else forward[e]=((p==(g.gamma_path+1)%3)==(which==2));
                if(g.head[e]==g.gamma) before_gamma=false;
            }
        }
        set_wick_sign(g);
    }
};
inline Entries transfer(const Graph& g,const Directions& f,int path) {
    Entries matrix;
    for(int i=0;i<5;++i) matrix[i][i].push_back({{},1});
    for(int e=g.first[path];e<g.last[path];++e) {
        const int v=g.head[e],u=g.slots[v-1].vector_axis;
        const int momentum=g.slots[v-1].momentum, axis=std::abs(momentum)-1;
        Entries next;
        for(int i=0;i<5;++i) for (int j=0;j<5;++j) {
            Term factor;
            if(f.forward[e] && !f.forward[e+1]) {
                // Signed internal momenta at c are y and -(y+lambda).
                int missing=-1,sign=0;
                for(int k=0;k<5;++k) {
                    sign=epsilon({i,j,u,k,axis});
                    if(sign){missing=k;break;}
                }
                if(!sign)continue;
                factor=g.current(e,missing,sign*(momentum>0?1:-1));
            } else if(!f.forward[e] && f.forward[e+1]) {
                int sign=int(i==u && j==axis)-int(i==axis && j==u);
                if (!sign) continue;
                factor={{},sign*(momentum>0?1:-1)};
            } else {
                // Transfer u J^T - (J.u)1; transpose for reverse flow.
                bool forward=f.forward[e];
                int row=forward?i:j,col=forward?j:i;
                if((row==col && row==u) || (row!=col && row!=u))continue;
                int sign=forward?1:-1;
                factor=g.current(e+int(!forward),row==col?u:col,row==col?-sign:sign);
            }
            for (int left=0;left<5;++left) for (auto entry:matrix[left][i])
                next[left][j].push_back(multiply(entry,factor));
        }
        matrix=std::move(next);
    }
    return matrix;
}
struct EndpointEntry { int a,b,c; Term term; };
inline std::vector<EndpointEntry> endpoint(const Graph& g,const Directions& f,bool root) {
    std::vector<EndpointEntry> entries;
    const auto& edges=root?g.last:g.first;
    int incident_sign=root?1:-1;
    bool cubic=f.forward[edges[0]]==root && f.forward[edges[1]]==root && f.forward[edges[2]]==root;
    int h=-1;
    if(!cubic)for(int k=0;k<3;++k)if(f.forward[edges[k]]==root)h=k;
    for(int a=0;a<5;++a)for(int b=0;b<5;++b)for(int c=0;c<5;++c) {
        if(cubic) {
            // k_A+k_B+k_C=0, so the endpoint tensor is -epsilon * k_B wedge k_C.
            if(a==b || a==c || b==c)continue;
            int i=-1,j=-1;
            for(int k=0;k<5;++k)if(k!=a&&k!=b&&k!=c){if(i<0)i=k;else j=k;}
            int s=-epsilon({a,b,c,i,j});
            entries.push_back({a,b,c,multiply(g.current(edges[1],i,s),g.current(edges[2],j,1))});
            entries.push_back({a,b,c,multiply(g.current(edges[1],j,-s),g.current(edges[2],i,1))});
        } else {
            int idx[3]={a,b,c};int t1=h==0?1:0,t2=h==2?1:2;
            if(idx[t1]==idx[t2])continue;
            if(idx[h]==idx[t1]) entries.push_back({a,b,c,g.current(edges[h],idx[t2],incident_sign)});
            else if(idx[h]==idx[t2]) entries.push_back({a,b,c,g.current(edges[h],idx[t1],-incident_sign)});
        }
    }
    return entries;
}
inline void collect(const Graph& g,const Directions& f,Collected& coefficients) {
    const auto A=transfer(g,f,0),B=transfer(g,f,1),C=transfer(g,f,2);
    const auto ln=endpoint(g,f,false),rn=endpoint(g,f,true);
    for(auto l:ln)for(auto r:rn) {
        const auto& as=A[l.a][r.a];const auto& bs=B[l.b][r.b];const auto& cs=C[l.c][r.c];
        if(as.empty() || bs.empty() || cs.empty())continue;
        auto end=multiply(l.term,r.term);end.sign*=f.sign;
        for(auto a:as) {
            auto first=multiply(end,a);
            for(auto b:bs) {
                auto second=multiply(first,b);
                for(auto c:cs){auto total=multiply(second,c);coefficients[total.p]+=total.sign;}
            }
        }
    }
}
// Compute the polynomial in 15 plain formal first currents before eliminating a.
// Keeping each substitution linear permits expansion once per surviving monomial.
inline void expand_shifted(Powers key,std::int64_t value,const Graph& g,Polynomial& out) {
    if(!key.shifted){out[key.plain]+=value;return;}
    const int axis=__builtin_ctzll(key.shifted)/4,shift=4*axis;
    const auto bit=std::uint64_t{1}<<shift;
    int degree=(key.shifted>>shift)&15;key.shifted-=degree*bit;
    int choose=1;
    for(int k=0;k<=degree;++k) {
        int sign=g.shifts[axis]<0 && ((degree-k)&1)?-1:1;
        expand_shifted({key.plain+k*bit,key.shifted},value*choose*sign,g,out);
        choose=choose*(degree-k)/(k+1);
    }
}
inline void eliminate_a(std::uint64_t key,std::int64_t value,Polynomial& out) {
    auto a_key=key>>40;
    if(!a_key){out[key]+=value;return;}
    int axis=__builtin_ctzll(a_key)/4,shift=axis*4;
    const auto bbit=std::uint64_t{1}<<shift,cbit=bbit<<20,abit=bbit<<40;
    int degree=(a_key>>shift)&15;key-=degree*abit;
    if(degree&1)value=-value;
    int choose=1;
    for(int k=0;k<=degree;++k) {
        eliminate_a(key+k*bbit+(degree-k)*cbit,value*choose,out);
        choose=choose*(degree-k)/(k+1);
    }
}
inline Polynomial expand(const Graph& graph,const Collected& factored) {
    Polynomial intermediate, result;
    for(auto [key,value]:factored)if(value)expand_shifted(key,value,graph,intermediate);
    for(auto [key,value]:intermediate)if(value)eliminate_a(key,value,result);
    for(auto it=result.begin();it!=result.end();) {
        if(!it->second)it=result.erase(it);else ++it;
    }
    return result;
}
} // namespace detail

inline Polynomial numerator(const Lengths& lengths,const Slots& slots,int gamma_index) {
    detail::Graph graph(lengths,slots,gamma_index);
    detail::Collected factored;
    for(int orientation=0;orientation<4;++orientation)detail::collect(graph,detail::Flow(graph,orientation),factored);
    return detail::expand(graph,factored);
}
} // namespace twoloop
