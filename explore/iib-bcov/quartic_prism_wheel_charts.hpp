#pragma once
#include "quartic_prism_graph.hpp"
#include <vector>

inline bool convergent(const Edges &g,int q) {
    for (int i=0;i<6;++i) for(int j=i+1;j<6;++j) {
        if(g[i][j]>=q) return false;
        for(int k=j+1;k<6;++k) if(g[i][j]+g[i][k]+g[j][k]>=2*q) return false;
    }
    return true;
}

// Put |x_k-x_l|^(2*degree) on diagonal (0,2). Both extreme radial
// monomials must converge before separating the angular expansion.
inline std::uint64_t numerator_wheel_key(const Edges &g,int q,int degree,bool require_convergence=true,
                                        const std::vector<std::uint64_t>* available=nullptr) {
    for(int i=0;i<6;++i) for(int j=i+1;j<6;++j) if(!g[i][j]) {
        for(int k=0;k<6;++k) for(int l=k+1;l<6;++l) {
            if(k==i || k==j || l==i || l==j || g[k][l]!=-degree) continue;
            for(int c=0;c<6;++c) {
                if(c==i || c==j || c==k || c==l) continue;
                std::array<int,8> powers{g[c][k],g[c][i],g[c][l],g[c][j],
                                         g[j][k],g[k][i],g[i][l],g[l][j]};
                if(*std::min_element(powers.begin(),powers.end())<=0) continue;
                int inf=15-i-j-k-l-c;
                bool valid=true;
                for(int v:{k,l}) {
                    Edges h=g;
                    h[k][l]=h[l][k]=0;
                    h[v][c]-=degree;h[c][v]-=degree;
                    int other=k+l-v;
                    h[inf][other]-=degree;h[other][inf]-=degree;
                    h[inf][c]+=degree;h[c][inf]+=degree;
                    valid &= convergent(h,q);
                }
                if(require_convergence && !valid) continue;
                std::uint64_t packed=0;
                for(int w:powers) packed=(packed<<4)|w;
                if(available && !std::binary_search(available->begin(),available->end(),packed))continue;
                return packed;
            }
        }
    }
    return 0;
}

