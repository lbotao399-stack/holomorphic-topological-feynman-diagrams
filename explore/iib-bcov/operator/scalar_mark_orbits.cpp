// Polarized scalar invariant of four signed five-row minors on nine conserved slots.
#include <algorithm>
#include <array>
#include <iostream>
#include <numeric>
#include <set>
#include <vector>
using Five=std::array<int,5>;
using Four=std::array<int,4>;
int key(const Five& a) {
    int n=0;
    for(int v:a) n=9*n+v;
    return n;
}
void print(__int128 n) {
    if(n<0) {std::cout<<'-'; n=-n;}
    if(n>9) print(n/10);
    std::cout<<char('0'+n%10);
}
int main() {
    int count, tests;
    if(!(std::cin>>count>>tests)) return 1;
    std::vector<Four> words(count);
    for(auto& w:words) for(int& j:w) std::cin>>j;
    std::vector<Five> subsets;
    std::array<int,59049> signs{};
    for(int a=0;a<5;++a) for(int b=a+1;b<6;++b) for(int c=b+1;c<7;++c)
    for(int d=c+1;d<8;++d) for(int e=d+1;e<9;++e) {
        Five q{a,b,c,d,e}; subsets.push_back(q);
        do {
            int s=1;
            for(int i=0;i<5;++i) for(int j=0;j<i;++j) if(q[j]>q[i]) s=-s;
            signs[key(q)]=s*int(subsets.size());
        } while(std::next_permutation(q.begin(),q.end()));
    }
    std::vector<int> used;
    for(const auto& w:words) for(int j:w) used.push_back(j);
    std::sort(used.begin(),used.end());
    used.erase(std::unique(used.begin(),used.end()),used.end());
    // A stabilizer with odd sign annihilates its whole symmetrization.
    std::vector<std::set<Four>> orbits(count);
    std::array<int,9> g;
    std::iota(g.begin(),g.end(),0);
    do {
        std::array<int,126> d{};
        for(int j:used) {
            Five v;
            for(int k=0;k<5;++k) v[k]=g[subsets[j][k]];
            d[j]=signs[key(v)];
        }
        for(int n=0;n<count;++n) {
            int sign=1;Four w;
            for(int k=0;k<4;++k) {
                int z=d[words[n][k]]; w[k]=std::abs(z);
                if(z<0) sign=-sign;
            }
            std::sort(w.begin(),w.end());
            w[0]*=sign;
            orbits[n].insert(w);
        }
    } while(std::next_permutation(g.begin(),g.end()));
    for(int t=0;t<tests;++t) {
        std::vector<long long> d(127);
        for(int j=1;j<=126;++j) std::cin>>d[j];
        for(const auto& orbit:orbits) {
            __int128 sum=0;
            for(auto w:orbit) {
                int x=w[0];sum+=(x<0?-1:1)*(__int128)d[std::abs(x)]*d[w[1]]*d[w[2]]*d[w[3]];
            }
            print(sum);std::cout<<' ';
        }
        std::cout<<'\n';
    }
}
