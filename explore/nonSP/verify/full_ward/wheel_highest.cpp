// Highest odd-zeta coefficient of an integer-power projective wheel, modulo a prime.
// Input: n q prime, then rows of n cyclic spoke powers and n rim powers (i,i+1).
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
#include <vector>

using Z = int64_t;
Z mod;
Z residue(Z v) { return (v%mod+mod)%mod; }
Z times(Z x,Z y) { return x*y%mod; }
Z power(Z x,Z n) {
    Z y=1;
    for(;n;n/=2,x=times(x,x)) if(n%2) y=times(y,x);
    return y;
}
struct Local { Z lead=1, slope=0; int poles=0; };
struct Ordering { std::vector<int> sign; Z all=0; std::vector<Z> gap; };

struct Wheel {
    int n,q;
    std::vector<Z> fac,inv;
    std::vector<Local> multiplicity;
    std::vector<std::vector<std::vector<Local>>> angular;
    std::vector<Ordering> orderings;
    Z reciprocal(int v) const {
        assert(v && std::abs(v)<int(inv.size()));
        return v<0 ? mod-inv[-v] : inv[v];
    }
    Wheel(int count,int dimension) : n(count),q(dimension),fac(2*q+1,1),inv(12*n*q+10) {
        assert(n>=3 && q>=2 && mod>int(inv.size()) && mod<(int64_t(1)<<31));
        for(int k=1;k<int(inv.size());++k) inv[k]=power(k,mod-2);
        for(int k=1;k<int(fac.size());++k) fac[k]=times(k,fac[k-1]);
        for(int h=0;h<q-1;++h) {
            Local d; d.lead=times(2,power(fac[2*q-2],mod-2));
            // The multiplicity has roots -(q-2),...,(q-2) and a second zero at the origin.
            for(int j=2-q;j<=q-2;++j) if(h!=j) {
                d.lead=times(d.lead,residue(h-j));
                d.slope=residue(d.slope+reciprocal(h-j));
            }
            if(h) {
                d.lead=times(d.lead,h);
                d.slope=residue(d.slope+reciprocal(h));
            }
            multiplicity.push_back(d);
            angular.emplace_back();
            for(int j=0;j<=q-2;++j) {
                angular[h].emplace_back();
                for(int k=0;k<=j;++k) {
                    Local e;
                    e.lead=times(fac[j],power(times(fac[k],fac[j-k]),mod-2));
                    if(k%2) e.lead=residue(-e.lead);
                    for(int z=0;z<=j;++z) {
                        int a=h+k-j+z;
                        if(!a) ++e.poles;
                        else {
                            e.lead=times(e.lead,reciprocal(a));
                            e.slope=residue(e.slope-reciprocal(a));
                        }
                    }
                    angular[h][j].push_back(e);
                }
            }
        }
        std::vector<int> p(n),pos(n),cross(n-1),sign(n);
        std::iota(p.begin(),p.end(),0);
        std::map<std::vector<int>,Ordering> groups;
        do {
            for(int i=0;i<n;++i) pos[p[i]]=i;
            for(int i=0;i<n;++i) sign[i]=pos[i]<pos[(i+1)%n] ? 1 : -1;
            Z weight=1;
            for(int j=0;j<n-1;++j) {
                cross[j]=0;
                for(int i=0;i<n;++i) cross[j]+=(pos[i]<=j)!=(pos[(i+1)%n]<=j);
                weight=times(weight,reciprocal(cross[j]));
            }
            auto& g=groups[sign];
            if(g.gap.empty()) {g.sign=sign;g.gap.resize(n*n);}
            g.all=residue(g.all+weight);
            for(int j=0;j<n-1;++j) {
                Z& v=g.gap[p[j]*n+p[j+1]];
                v=residue(v+weight*cross[j]);
            }
        } while(std::next_permutation(p.begin(),p.end()));
        for(auto& kv:groups) orderings.push_back(std::move(kv.second));
    }
    Z coefficient(const std::vector<int>& sp,const std::vector<int>& rim) const {
        assert(std::accumulate(sp.begin(),sp.end(),0)+std::accumulate(rim.begin(),rim.end(),0)==n*q);
        std::vector<int> d(n),c(n);
        Z constant=power(2,n-1),answer=0;
        for(int i=0;i<n;++i) {
            assert(sp[i]>0 && sp[i]<q && rim[i]>0 && rim[i]<q);
            d[i]=rim[i]-1;
            constant=times(constant,fac[q-sp[i]-1]);
            c[i]=(i?c[i-1]:0)+2*sp[i]+rim[(i+n-1)%n]+rim[i]-2*q;
        }
        assert(c.back()==0);
        const int first=std::min_element(d.begin(),d.end())-d.begin();
        for(int h=0;h<=*std::max_element(d.begin(),d.end());++h) {
            const auto& mult=multiplicity[h];
            Z contribution=0;
            for(const auto& g:orderings) {
                // A zero gap at every cut forces sign[i]*p[i]+c[i] constant.
                for(int k=0;k<=d[first];++k) {
                    int w=g.sign[first]*(h+2*k-d[first])+c[first];
                    Local row;
                    for(int i=0;i<n;++i) {
                        int twice=g.sign[i]*(w-c[i])-h+d[i];
                        if(twice<0 || twice>2*d[i] || twice%2) {row.lead=0;break;}
                        const auto& e=angular[h][d[i]][twice/2];
                        row.lead=times(row.lead,e.lead);
                        row.poles+=e.poles;
                        row.slope=residue(row.slope+e.slope);
                    }
                    int excess=row.poles+(n-1)-(h?1:2)-(2*n-3);
                    if(excess<0 || !row.lead) continue;
                    assert(excess<=1);
                    Z v=times(row.lead,g.all);
                    if(excess) v=times(v,residue(mult.slope+row.slope));
                    contribution=residue(contribution+v);
                }
                // Exactly one nonzero gap can reach weight 2n-3 only at h>0 with n angular poles.
                if(!h || d[first]<h) continue;
                for(int u=0;u<n;++u) for(int v=0;v<n;++v) {
                    Z weight=g.gap[u*n+v];
                    if(!weight) continue;
                    auto arc=[&](int a,int end,int k) {
                        int w=g.sign[a]*(h+2*k-d[a])+c[a];
                        Z z=1;
                        for(int i=a;i!=end;i=(i+1)%n) {
                            int twice=g.sign[i]*(w-c[i])-h+d[i];
                            if(twice<0 || twice>2*(d[i]-h) || twice%2) return std::pair<int,Z>{w,0};
                            z=times(z,angular[h][d[i]][twice/2].lead);
                        }
                        return std::pair<int,Z>{w,z};
                    };
                    for(int ku=0;ku<=d[u]-h;++ku) {
                        auto [a,x]=arc(u,v,ku);
                        if(!x) continue;
                        for(int kv=0;kv<=d[v]-h;++kv) {
                            auto [b,y]=arc(v,u,kv);
                            if(y && a!=b) contribution=residue(contribution+
                                times(times(weight,x),times(y,reciprocal(a-b))));
                        }
                    }
                }
            }
            answer=residue(answer+times((h?2:1)*mult.lead%mod,contribution));
        }
        return times(constant,answer);
    }
};

int main() {
    int n,q;
    if(!(std::cin>>n>>q>>mod)) return 1;
    Wheel wheel(n,q);
    std::vector<int> spoke(n),rim(n);
    while(std::cin>>spoke[0]) {
        for(int i=1;i<n;++i) std::cin>>spoke[i];
        for(auto& r:rim) std::cin>>r;
        if(!std::cin) return 1;
        std::cout<<wheel.coefficient(spoke,rim)<<'\n';
    }
}
