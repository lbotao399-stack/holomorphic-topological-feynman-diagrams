// Time-grid values of the five-loop source numerator in the prime field.
#include <array>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <vector>
#ifndef GRID_PRIME
#define GRID_PRIME 1000003
#endif
constexpr int prime=GRID_PRIME;
using Values=std::array<int,21>; // U followed by two source numerators per nonpivot edge
struct Tree {int mask; Values v;};
struct Minor {int sign,left,right;};

static int reduce(int64_t x) {int r=int(x%prime); return r<0?r+prime:r;}

static std::array<int,32> products(const Values& v,int start) {
    std::array<int,32> a{};
    a[0]=1;
    for(int i=0;i<5;++i) {
        int x=v[1+2*(start+i)],y=v[2+2*(start+i)];
        for(int j=0;j<(1<<i);++j) {
            a[j+(1<<i)]=reduce(int64_t(a[j])*x);
            a[j]=reduce(int64_t(a[j])*y);
        }
    }
    return a;
}

static void suffix(int n,const Values* v,const std::vector<int>& inv,
                   const std::vector<Minor>& minors,std::vector<int32_t>& output) {
    if(n==0) {
        auto l=products(v[0],0),r=products(v[0],5);
        int64_t det=0;
        for(const auto& m:minors)det+=int64_t(m.sign)*l[m.left]*r[m.right];
        assert(v[0][0]>0&&size_t(v[0][0])<inv.size());
        output.push_back(reduce(int64_t(reduce(det))*inv[v[0][0]]));
        return;
    }
    Values next[16];
    int length=1<<(n-1);
    for(int t=1;t<=5;++t) {
        for(int b=0;b<length;++b)
            for(int f=0;f<21;++f)next[b][f]=v[2*b][f]+t*v[2*b+1][f];
        suffix(n-1,next,inv,minors,output);
    }
}

int main(int argc,char** argv) {
    assert(argc==3);
    std::ifstream input(argv[1]);
    int nt,nm;input>>nt>>nm;
    assert(input&&nt>0&&nm>0);
    std::vector<Tree> trees(nt);
    for(auto& t:trees) {
        input>>t.mask;
        for(int& x:t.v){input>>x;assert(x>=-8&&x<=8);}
    }
    std::vector<Minor> minors(nm);
    for(auto& m:minors)input>>m.sign>>m.left>>m.right;
    assert(input);
    int maximum=nt*3125;
    assert(maximum<prime);
    std::vector<int> inv(maximum+1);
    inv[1]=1;
    for(int i=2;i<=maximum;++i)inv[i]=prime-int(int64_t(prime/i)*inv[prime%i]%prime);
    for(int& v:inv){int sq=reduce(int64_t(v)*v);v=reduce(int64_t(reduce(int64_t(sq)*sq))*v);}
    std::ofstream out(argv[2],std::ios::binary);
    for(int prefix=0;prefix<3125;++prefix) {
        int code=prefix,t[5];
        for(int i=4;i>=0;--i){t[i]=1+code%5;code/=5;}
        int scale[32]={1};
        for(int i=0;i<5;++i)for(int j=0;j<(1<<i);++j)scale[j+(1<<i)]=scale[j]*t[i];
        Values coeff[32]{};
        for(const auto& tree:trees) {
            auto& target=coeff[tree.mask>>5];
            int weight=scale[tree.mask&31];
            for(int i=0;i<21;++i)target[i]+=weight*tree.v[i];
        }
        std::vector<int32_t> result;result.reserve(3125);
        suffix(5,coeff,inv,minors,result);
        out.write(reinterpret_cast<const char*>(result.data()),result.size()*sizeof(int32_t));
    }
    assert(out);
}
