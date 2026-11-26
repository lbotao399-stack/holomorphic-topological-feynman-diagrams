// Signed orbit sums of six five-row determinants on eight differentiated momenta.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>
using Five=std::array<int,5>;
int key(const Five& ids) {
    int n=0;for(int i:ids)n=8*n+i;return n;
}
int main() {
    int count, tests;if(!(std::cin>>count>>tests)) return 1;
    std::vector<std::array<int,6>> words(count);
    for(auto& w:words) for(int& x:w) std::cin>>x;
    std::vector<Five> subsets;
    std::array<int,32768> lookup{};
    for(int a=0;a<4;++a) for(int b=a+1;b<5;++b) for(int c=b+1;c<6;++c)
    for(int d=c+1;d<7;++d) for(int e=d+1;e<8;++e) {
        Five q{a,b,c,d,e};subsets.push_back(q);
        do {
            int sign=1;
            for(int i=0;i<5;++i) for(int j=0;j<i;++j) if(q[j]>q[i]) sign=-sign;
            lookup[key(q)]=sign*int(subsets.size());
        } while(std::next_permutation(q.begin(),q.end()));
    }
    std::vector<std::array<int8_t,56>> permuted;
    std::array<int,8> perm;std::iota(perm.begin(),perm.end(),0);
    do {
        std::array<int8_t,56> m;
        for(int i=0;i<56;++i) {
            int k=0;for(int j:subsets[i])k=8*k+perm[j];
            m[i]=lookup[k];
        }
        permuted.push_back(m);
    } while(std::next_permutation(perm.begin(),perm.end()));
    for(int t=0;t<tests;++t) {
        std::array<int,113> d{};
        for(int i=1;i<=56;++i) {
            int x;std::cin>>x;x%=1009;d[56+i]=x;d[56-i]=-x;
        }
        std::vector<long long> sums(count);
        for(const auto& perm:permuted) {
            std::array<int,56> v;
            for(int i=0;i<56;++i)v[i]=d[56+perm[i]];
            for(int i=0;i<count;++i) {
                const auto& word=words[i];
                long long a=1LL*v[word[0]]*v[word[1]]*v[word[2]]%1009;
                long long b=1LL*v[word[3]]*v[word[4]]*v[word[5]]%1009;
                sums[i]+=a*b;
            }
        }
        for(long long x:sums)std::cout<<(x%1009+1009)%1009<<' ';
        std::cout<<'\n';
    }
}
