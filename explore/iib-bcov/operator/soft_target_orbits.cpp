// Pair triple-minor columns by exact, bounded integer Gram products on the S_8 stabilizer.
#include <Accelerate/Accelerate.h>
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>
using Five=std::array<int,5>;
constexpr int prime=1009;
int bit(int x){return 1<<x;}
int main(){
    int test_count,k,count;if(!(std::cin>>test_count>>k>>count))return 1;
    std::vector<std::array<int,3>> triples(k);
    std::vector<std::array<int,2>> words(count);
    for(auto& w:triples)for(int& i:w)std::cin>>i;
    for(auto& w:words)for(int& i:w)std::cin>>i;
    std::vector<Five> subsets;
    std::array<int,1024> subset_id{};
    for(int a=0;a<6;++a)for(int b=a+1;b<7;++b)for(int c=b+1;c<8;++c)
    for(int d=c+1;d<9;++d)for(int e=d+1;e<10;++e){
        subsets.push_back({a,b,c,d,e});subset_id[bit(a)+bit(b)+bit(c)+bit(d)+bit(e)]=subsets.size();
    }
    auto map_subset=[&](const Five& s,const std::array<int,10>& perm) {
        int mask=0,sign=1;
        for(int h:s){int v=perm[h];if(__builtin_popcount(mask>>(v+1))%2)sign=-sign;mask+=bit(v);}
        return sign*subset_id[mask];
    };
    std::array<int,10> perm;std::iota(perm.begin(),perm.end(),0);
    std::vector<std::array<int16_t,252>> stabilizer;
    stabilizer.reserve(40320);
    do{
        std::array<int16_t,252> m;
        for(int s=0;s<252;++s)m[s]=map_subset(subsets[s],perm);
        stabilizer.push_back(m);
    }while(std::next_permutation(perm.begin()+2,perm.end()));
    std::vector<double> features(40320*k),gram(k*k);
    for(int t=0;t<test_count;++t){
        std::array<Five,10> p;for(auto& r:p)for(int& x:r)std::cin>>x;
        std::array<int,505> minors{};
        for(int i=1;i<=252;++i){int x;std::cin>>x;minors[252+i]=x;minors[252-i]=-x;}
        std::vector<std::array<long long,50>> results(count);
        for(int i=0;i<10;++i)for(int j=0;j<10;++j)if(i!=j){
            perm[0]=i;perm[1]=j;int pos=2;
            for(int h=0;h<10;++h)if(h!=i&&h!=j)perm[pos++]=h;
            std::array<int,505> d{};
            for(int s=0;s<252;++s){
                int v=minors[252+map_subset(subsets[s],perm)];
                d[253+s]=v;d[251-s]=-v;
            }
            double* out=features.data();
            for(auto& m:stabilizer){
                std::array<int,252> v;
                for(int h=0;h<252;++h)v[h]=d[252+m[h]];
                for(auto& w:triples)*out++=(v[w[0]]*v[w[1]]*v[w[2]])%prime;
            }
            cblas_dsyrk(CblasRowMajor,CblasUpper,CblasTrans,k,40320,1.0,features.data(),k,0.0,gram.data(),k);
            for(int w=0;w<count;++w){
                auto ab=words[w];
                long long x=static_cast<long long>(gram[k*ab[0]+ab[1]])%prime;
                for(int h=0;h<5;++h)results[w][5*i+h]+=x*p[j][h];
            }
        }
        for(const auto& col:results){
            for(auto x:col)std::cout<<(x%prime+prime)%prime<<' ';
            std::cout<<'\n';
        }
        std::cout.flush();
    }
}
