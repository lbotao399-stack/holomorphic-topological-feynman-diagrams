// Transpose the time interpolation: apply each exact projective integration
// weight to all source probes at the same time point. No floating point.
#include <array>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include <gmpxx.h>

using Wide=__int128_t;
using Source=std::array<std::array<int,2>,5>;
constexpr int points=65536;
constexpr int digits=4;
using Digits=std::array<int32_t,digits>;
using Accumulator=std::array<std::array<Wide,digits>,3>;
struct Tree {std::array<int,4> missing; int flow[8][5];};
struct Minor {int sign; int edges[8];};

mpz_class integer(Wide v) {
    bool neg=v<0;
    __uint128_t x=neg ? -static_cast<__uint128_t>(v) : static_cast<__uint128_t>(v);
    mpz_class result;
    mpz_import(result.get_mpz_t(),1,1,sizeof(x),0,0,&x);
    return neg ? -result : result;
}

int main(int argc,char**argv) {
    const bool grid=argc==5 && std::string(argv[4])=="--grid";
    assert(argc==4 || grid);
    std::ifstream graph(argv[1]),probes(argv[2]);
    std::ofstream out(argv[3], grid ? std::ios::binary : std::ios::out);
    int nt,nm,np;
    graph>>nt>>nm;probes>>np;
    assert(graph&&probes&&out&&nt>0&&nm>0&&np>0);
    std::vector<Tree> trees(nt);
    for(auto& tree:trees) {
        for(int& e:tree.missing)graph>>e;
        for(auto& row:tree.flow)for(int& v:row)graph>>v;
    }
    std::vector<Minor> minors(nm);
    for(auto& m:minors) {
        graph>>m.sign;
        for(int& e:m.edges){graph>>e;--e;assert(e>=0&&e<8);}
    }
    std::vector<std::array<Digits,3>> weights(points);
    std::array<mpz_class,3> denom;
    std::array<int,3> active{};
    if(!grid)for(int j=0;j<3;j++){
        graph>>denom[j];
        for(auto& entry:weights) {
            mpz_class w;graph>>w;
            int s=sgn(w);w=abs(w);
            for(int k=0;k<digits;k++) {
                entry[j][k]=s*int(mpz_get_ui(w.get_mpz_t())&((1u<<30)-1));
                if(entry[j][k])active[j]=std::max(active[j],k+1);
                w>>=30;
            }
            assert(w==0);
        }
    }
    std::vector<Source> source(np);
    for(auto& p:source)for(auto& row:p)for(int& v:row){probes>>v;assert(abs(v)<=4);}
    assert(graph&&probes);
    std::vector<Accumulator> a(np);
    for(int code=0;code<points;code++) {
        int t[9]={1};int n=code;
        for(int e=8;e>0;e--){t[e]=1+n%4;n/=4;}
        int u=0,flow[8][5]{};
        for(const auto& tree:trees) {
            int w=1;for(int e:tree.missing)w*=t[e];
            u+=w;
            for(int e=0;e<8;e++)for(int j=0;j<5;j++)flow[e][j]+=w*tree.flow[e][j];
        }
        const Wide divisor=Wide(u)*u*u*u;
        for(int i=0;i<np;i++) {
            int y[8][2]{};
            for(int e=0;e<8;e++)for(int j=0;j<5;j++)
                for(int d=0;d<2;d++)y[e][d]+=flow[e][j]*source[i][j][d];
            Wide det=0;
            for(const auto& m:minors) {
                // Products of two int32 flows fit int64, and four fit Wide
                // under the six-source bound. Guard the eight-flow product.
                Wide x=Wide(int64_t(y[m.edges[0]][0])*y[m.edges[1]][0])*
                       (int64_t(y[m.edges[2]][0])*y[m.edges[3]][0]);
                Wide z=Wide(int64_t(y[m.edges[4]][1])*y[m.edges[5]][1])*
                       (int64_t(y[m.edges[6]][1])*y[m.edges[7]][1]);
                Wide term,next;
                bool failed=__builtin_mul_overflow(x,z,&term);
                assert(!failed);
                failed=m.sign>0 ? __builtin_add_overflow(det,term,&next)
                                : __builtin_sub_overflow(det,term,&next);
                assert(!failed);det=next;
            }
            Wide exact=det/divisor;
            int64_t p=static_cast<int64_t>(exact);
            assert(Wide(p)==exact&&exact*divisor==det);
            if(grid) {
                out.write(reinterpret_cast<const char*>(&p),sizeof(p));
                continue;
            }
            // Each of at most 2^16 terms is bounded by 2^63 times 2^30.
            for(int j=0;j<3;j++)for(int k=0;k<active[j];k++)
                a[i][j][k]+=Wide(p)*weights[code][j][k];
        }
    }
    if(!grid)for(const auto& row:a) {
        for(int j=0;j<3;j++) {
            mpz_class v=0;
            for(int k=digits-1;k>=0;k--){v<<=30;v+=integer(row[j][k]);}
            mpq_class coefficient(v,denom[j]);coefficient.canonicalize();out<<coefficient<<' ';
        }
        out<<'\n';
    }
    assert(out);
}
