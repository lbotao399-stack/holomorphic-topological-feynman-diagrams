#pragma once

#include <array>
#include <cstdint>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace theta166 {
struct Leg { int momentum; int vector_axis; };
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

// The graph's geometric orientation points N -> N' along all three paths.
// The four rows are the field directions on A, on first and remaining B edges,
// and on C. true denotes geometric orientation. The c vertex is the only sink
// of three gamma slots; in rows 0 and 2 two of these slots are internal.
struct Flow { bool a,b_first,b_rest,c; int c_vertex; int gamma_vertex=1; };
inline constexpr std::array<Flow,4> flows = {{{true,true,false,false,1},
    {true,false,true,true,11},{false,true,false,true,1},{false,false,true,false,0}}};
inline int edge_end(int e, bool head) {
    if (e == 0) return head ? 11 : 0;
    if (e <= 6) return head ? (e == 6 ? 11 : e) : e - 1;
    return head ? (e == 12 ? 11 : e - 1) : (e == 7 ? 0 : e - 2);
}
inline bool aligned(const Flow& f, int e) {
    if (e == 0) return f.a;
    if (e <= f.gamma_vertex) return f.b_first;
    return e <= 6 ? f.b_rest : f.c;
}
inline int tail_parity(const Flow& f) {
    int s = 1;
    for (int i=0;i<13;++i) for (int j=i+1;j<13;++j)
        if (edge_end(i,!aligned(f,i)) > edge_end(j,!aligned(f,j))) s = -s;
    return s;
}
inline int epsilon(const std::array<int,5>& a) {
    int sign=1;
    for (int i=0;i<5;++i) for (int j=i+1;j<5;++j) {
        if (a[i]==a[j]) return 0;
        if (a[i]>a[j]) sign=-sign;
    }
    return sign;
}
struct Fixture {
    const Slots& slots;
    std::array<std::array<int,5>,13> offsets{};
    std::array<int,10> shifts{};
    explicit Fixture(const Slots& inputs): slots(inputs) {
        for (int path=0;path<2;++path) {
            std::array<int,5> q{};
            const int first=1+6*path;
            for (int j=0;j<5;++j) {
                auto leg = slots[5*path+j];
                if (leg.momentum==0 || leg.momentum < -5 || leg.momentum>5 ||
                    leg.vector_axis<0 || leg.vector_axis>4) throw std::invalid_argument("invalid leg");
                int axis = (leg.momentum>0 ? leg.momentum : -leg.momentum)-1;
                if (q[axis]!=0) throw std::invalid_argument("a theta166 momentum determinant vanishes");
                q[axis] = leg.momentum>0 ? 1 : -1;
                offsets[first+j+1]=q;
            }
            for (int i=0;i<5;++i) shifts[path*5+i]=q[i];
        }
        for (int i=0;i<5;++i) if (shifts[i] != -shifts[5+i])
            throw std::invalid_argument("the two empty endpoints require zero total momentum");
    }
    Term variable(int axis,int offset,int sign) const {
        const auto flag=std::uint64_t{1}<<(4*axis);
        return {{offset ? 0 : flag,offset ? flag : 0},sign};
    }
    Term path_current(int e,int i,int multiplier) const {
        const int axis = e<=6 ? i : i+5;
        return variable(axis,offsets[e][i],multiplier);
    }
    std::vector<Term> current(int e,int i,int sign) const {
        // y_A=-b-c. All other edge currents consist of one (possibly shifted) factor.
        if (e==0) return {variable(i,0,-sign),variable(i+5,0,-sign)};
        return {path_current(e,i,sign)};
    }
};

inline Entries pass_chain(const Fixture& g,int first,int last,int offset,bool forward) {
    Entries prev;
    for (int i=0;i<5;++i) prev[i][i].push_back({{},1});
    for (int vertex=first;vertex<=last;++vertex) {
        Entries next;
        int u=g.slots[vertex-1].vector_axis;
        int e=vertex+offset+int(!forward);
        int current_sign=forward ? 1 : -1;
        for (int incoming=0;incoming<5;++incoming) for (int outgoing=0;outgoing<5;++outgoing) {
            // Forward transfer is u J^T - (Ju)1, reverse its transpose.
            int row=forward?incoming:outgoing, col=forward?outgoing:incoming;
            if (row==col && row==u) continue;
            if (row!=col && row!=u) continue;
            auto factor = g.path_current(e,row==col?u:col,row==col?-current_sign:current_sign);
            for (int left=0;left<5;++left) for (Term t:prev[left][incoming])
                next[left][outgoing].push_back(multiply(t,factor));
        }
        prev=std::move(next);
    }
    return prev;
}
inline Entries b_chain(const Fixture& g,const Flow& f) {
    const int vertex=f.gamma_vertex;
    const auto left=pass_chain(g,1,vertex-1,0,f.b_first);
    const auto right=pass_chain(g,vertex+1,5,0,f.b_rest);
    Entries result;
    const int moment=g.slots[vertex-1].momentum, gamma=g.slots[vertex-1].vector_axis;
    const int lambda_axis=(moment>0?moment:-moment)-1, lambda_sign=moment>0?1:-1;
    for (int a=0;a<5;++a) for (int k=0;k<5;++k) {
        Term factor;
        if (f.c_vertex==vertex) {
            // c(a, k, gamma): the two signed edge momenta are y and -(y+lambda).
            bool nonzero=false;
            for (int i=0;i<5;++i) {
                int s=epsilon({a,k,gamma,i,lambda_axis});
                if (s) {factor=g.path_current(vertex,i,s*lambda_sign);nonzero=true;break;}
            }
            if (!nonzero) continue;
        } else {
            int s=lambda_sign*(int(a==gamma && k==lambda_axis)-int(a==lambda_axis && k==gamma));
            if (!s) continue;
            factor={{},s};
        }
        for (int i=0;i<5;++i) for (auto l:left[i][a]) {
            auto first=multiply(l,factor);
            for (int j=0;j<5;++j) for (auto r:right[k][j])
                result[i][j].push_back(multiply(first,r));
        }
    }
    return result;
}
struct EndpointEntry { int a,b,c; Term term; };
inline std::vector<EndpointEntry> endpoint(const Fixture& g,const Flow& f,bool root) {
    std::vector<EndpointEntry> entries;
    const std::array<int,3> edges = root ? std::array<int,3>{0,6,12} : std::array<int,3>{0,1,7};
    const int current_sign=root ? 1 : -1;
    const bool is_c=f.c_vertex==(root?11:0);
    int h=-1;
    if (!is_c) for (int i=0;i<3;++i) if (aligned(f,edges[i])==root) h=i;
    for (int a=0;a<5;++a) for (int b=0;b<5;++b) for (int c=0;c<5;++c) {
        std::array<int,3> idx={a,b,c};
        if (is_c) {
            // -epsilon * k_B wedge k_C, since k_A + k_B + k_C = 0.
            for (int i=0;i<5;++i) for (int j=0;j<5;++j) {
                int s=-epsilon({a,b,c,i,j});
                if (s) entries.push_back({a,b,c,multiply(g.path_current(edges[1],i,s),g.path_current(edges[2],j,1))});
            }
        } else {
            int t1= h==0 ? 1 : 0, t2=h==2?1:2;
            if (idx[t1]==idx[t2]) continue;
            int axis, s;
            if (idx[h]==idx[t1]) {axis=idx[t2];s=current_sign;}
            else if (idx[h]==idx[t2]) {axis=idx[t1];s=-current_sign;}
            else continue;
            for (auto term:g.current(edges[h],axis,s)) entries.push_back({a,b,c,term});
        }
    }
    return entries;
}
inline void collect(const Fixture& g,const Flow& f,Collected& coefficients) {
    const auto B=b_chain(g,f), C=pass_chain(g,6,10,1,f.c);
    const auto at_n=endpoint(g,f,false), at_root=endpoint(g,f,true);
    const int wick=tail_parity(f);
    for (auto l:at_n) for (auto r:at_root) {
        if (l.a!=r.a) continue;
        const auto endpoints=multiply(l.term,r.term);
        for (auto tb:B[l.b][r.b]) {
            auto left=multiply(endpoints,tb);
            for (auto tc:C[l.c][r.c]) {
                auto t=multiply(left,tc);
                coefficients[t.p]+=wick*t.sign;
            }
        }
    }
}
inline void expand(Powers key,std::int64_t value,const Fixture& g,Polynomial& out) {
    if (!key.shifted) {out[key.plain]+=value;return;}
    const int axis=__builtin_ctzll(key.shifted)/4, shift=4*axis;
    const auto bit=std::uint64_t{1}<<shift;
    int degree=(key.shifted>>shift)&15;
    key.shifted-=degree*bit;
    int binomial=1;
    for(int k=0;k<=degree;++k) {
        int sign=(g.shifts[axis]==-1 && ((degree-k)&1)) ? -1 : 1;
        expand({key.plain+k*bit,key.shifted},value*binomial*sign,g,out);
        binomial=binomial*(degree-k)/(k+1);
    }
}
} // namespace detail

inline Polynomial numerator(const Slots& slots,int gamma_slot=0) {
    if (gamma_slot<0 || gamma_slot>4) throw std::invalid_argument("gamma_slot must lie on B");
    detail::Fixture fixture(slots);
    detail::Collected coefficients;
    for (auto flow:detail::flows) {
        flow.gamma_vertex=gamma_slot+1;
        if (flow.c_vertex==1) flow.c_vertex=flow.gamma_vertex;
        detail::collect(fixture,flow,coefficients);
    }
    Polynomial result;
    for (auto [key,value]:coefficients) if (value) detail::expand(key,value,fixture,result);
    for (auto it=result.begin();it!=result.end();) {
        if (!it->second) it=result.erase(it); else ++it;
    }
    return result;
}
} // namespace theta166
