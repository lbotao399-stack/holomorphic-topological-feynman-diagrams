// Eight finite momenta and a soft covector, or nine finite slots, in a symmetric trace.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
#include <vector>
using Five=std::array<int,5>;
int key(const Five& a) {int v=0;for(int x:a)v=9*v+x;return v;}
struct Order { std::array<int8_t,126> image; int soft, penultimate, antepenultimate; };
int main() {
    int count,tests,linear;
    if(!(std::cin>>count>>tests>>linear)) return 1;
    int modulus=1009;
    if(linear==3 || linear==6 || linear==7) {std::cin>>modulus;linear=linear==3?0:linear-2;}
    if(linear==9)std::cin>>modulus; // Leading degree two in each of three soft slots.
    if(modulus>1000000 || modulus<2) return 2;
    std::vector<std::array<int,6>> words(count);
    for(auto& w:words) for(int& j:w) std::cin>>j;
    std::vector<Five> subsets;
    std::array<int,59049> signs{};
    for(int a=0;a<5;++a) for(int b=a+1;b<6;++b) for(int c=b+1;c<7;++c)
    for(int d=c+1;d<8;++d) for(int e=d+1;e<9;++e) {
        Five ids{a,b,c,d,e};subsets.push_back(ids);
        do {
            int s=1;
            for(int i=0;i<5;++i)for(int j=0;j<i;++j)if(ids[j]>ids[i])s=-s;
            signs[key(ids)]=s*int(subsets.size());
        }while(std::next_permutation(ids.begin(),ids.end()));
    }
    bool derivativeMode=linear==2 || linear==8;
    std::array<std::vector<int>,9> active;
    std::vector<std::array<int,9>> degrees(count);
    for(int t=0;t<9;++t)for(int c=0;c<count;++c){
        int degree=0;for(int b:words[c])for(int j:subsets[b])degree+=(j==t);
        degrees[c][t]=degree;
        if((linear!=1 || degree==1) && (linear<4 || degree==2))active[t].push_back(c);
    }
    std::vector<Order> orbit;orbit.reserve(362880);
    std::array<int,9> g;std::iota(g.begin(),g.end(),0);
    do{
        Order ord;
        for(int i=0;i<9;++i){
            if(g[i]==8)ord.soft=i;
            if(g[i]==7)ord.penultimate=i;
            if(g[i]==6)ord.antepenultimate=i;
        }
        if(active[ord.soft].empty())continue;
        for(int h=0;h<126;++h){
            int z=0;for(int i:subsets[h])z=9*z+g[i];
            ord.image[h]=signs[z];
        }
        orbit.push_back(ord);
    }while(std::next_permutation(g.begin(),g.end()));
    // In symbolic mode each minor is supplied as zero or a signed factor ID.
    // Collect equal commutative products before expanding any polynomial.
    if(linear==-1){
        std::array<int,253> factors{};
        for(int h=1;h<=126;++h){
            int id; if(!(std::cin>>id))return 1;
            factors[126+h]=id;
            factors[126-h]=-id;
        }
        for(auto const& word:words){
            std::map<std::array<int,6>,int> multiplicities;
            for(auto const& g:orbit){
                std::array<int,6> product;
                int sign=1;
                for(int j=0;j<6;++j){
                    int id=factors[126+g.image[word[j]]];
                    if(id==0){ sign=0; break; }
                    if(id<0)sign=-sign;
                    product[j]=std::abs(id);
                }
                if(sign){
                    std::sort(product.begin(),product.end());
                    multiplicities[product]+=sign;
                }
            }
            for(auto const& [product,multiple]:multiplicities)if(multiple){
                std::cout<<multiple;
                for(int id:product)std::cout<<' '<<id;
                std::cout<<'\n';
            }
            std::cout<<"end\n";
        }
        return 0;
    }
    for(int t=0;t<tests;++t){
        std::array<int,253> d{};
        for(int h=1;h<=126;++h){int v;std::cin>>v;v%=modulus;d[126+h]=v;d[126-h]=-v;}
        std::array<int,253> dd{};
        if(derivativeMode)for(int h=1;h<=126;++h){int v;std::cin>>v;v%=modulus;dd[126+h]=v;dd[126-h]=-v;}
        std::vector<long long> sums(count);
        for(const auto& g:orbit){
            std::array<int,126> v;
            for(int h=0;h<126;++h)v[h]=d[126+g.image[h]];
            std::array<int,126> dv;
            if(derivativeMode) for(int h=0;h<126;++h)dv[h]=dd[126+g.image[h]];
            for(int c:active[g.soft]){
                if((linear==5 || linear==9) && degrees[c][g.penultimate]!=2)continue;
                if(linear==9 && degrees[c][g.antepenultimate]!=2)continue;
                auto w=words[c];
                long long x=1LL*v[w[0]]*v[w[1]]*v[w[2]]%modulus;
                long long y=1LL*v[w[3]]*v[w[4]]*v[w[5]]%modulus;
                if(derivativeMode){
                    long long dx=1LL*dv[w[0]]*v[w[1]]*v[w[2]]+1LL*v[w[0]]*dv[w[1]]*v[w[2]]
                                  +1LL*v[w[0]]*v[w[1]]*dv[w[2]];
                    long long dy=1LL*dv[w[3]]*v[w[4]]*v[w[5]]+1LL*v[w[3]]*dv[w[4]]*v[w[5]]
                                  +1LL*v[w[3]]*v[w[4]]*dv[w[5]];
                    sums[c]+=(dx%modulus)*y+x*(dy%modulus);
                }else sums[c]+=x*y;
            }
        }
        for(auto x:sums)std::cout<<(x%modulus+modulus)%modulus<<' ';
        std::cout<<'\n';
    }
}
