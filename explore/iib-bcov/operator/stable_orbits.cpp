// Polarize four alternating determinants at up to thirteen scalar slots.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>
using std::array;
using std::vector;
constexpr int modulus=1009;
using Integer=__int128;
using Vec=array<long long,5>;

void write_integer(Integer x) {
    if(x<0) {std::cout<<'-';x=-x;}
    if(x>=10) write_integer(x/10);
    std::cout<<char('0'+x%10);
}

struct Wave {
    vector<Vec> p;
    vector<long long> minors;
};

struct Orbit {
    int n,all,mode,free_leg=-1,base_sign=1;
    vector<int> flags,order;
    array<int,4> masks{};
    const vector<Wave>& waves;
    vector<array<Integer,5>> total;
    Orbit(int labels,int modular,const vector<int>& f,const vector<Wave>& w)
          : n(labels),all((1<<n)-1),mode(modular),flags(f),waves(w),total(w.size()) {
        order.resize(f.size());std::iota(order.begin(),order.end(),0);
        // Place equal incidence types consecutively, first placing the types that fill many volumes.
        std::sort(order.begin(),order.end(),[&](int a,int b) {
            auto score=[&](int i) {return std::pair<int,int>{__builtin_popcount(f[i]&15),f[i]};};
            return score(a)>score(b);
        });
        for(int i=0;i<(int)order.size();++i) {
            if(f[order[i]]&16) free_leg=i;
            for(int j=0;j<i;++j) if(order[j]>order[i])
                if(__builtin_popcount(f[order[j]]&f[order[i]]&15)%2) base_sign=-base_sign;
        }
    }
    void sum(int slot,int taken,int prev,int free_index,int sign) {
        if(slot==(int)order.size()) {
            for(int t=0;t<(int)waves.size();++t) {
                const auto& d=waves[t].minors;
                long long a=d[masks[0]],b=d[masks[1]],c=d[masks[2]],e=d[masks[3]];
                if(!a||!b||!c||!e) continue;
                Integer val=Integer(sign)*a*b*c*e;
                if(mode) val%=modulus;
                if(free_index<0) total[t][0]+=val;
                else for(int h=0;h<5;++h) total[t][h]+=val*waves[t].p[free_index][h];
            }
            return;
        }
        int f=flags[order[slot]],avail=all^taken;
        if(slot && flags[order[slot-1]]==f) avail&=~((1<<(prev+1))-1);
        while(avail) {
            int bit=avail&(-avail),j=__builtin_ctz(bit);avail-=bit;
            int parity=0;
            for(int h=0;h<4;++h) if((f>>h)&1) {
                parity+=__builtin_popcount(unsigned(masks[h])>>(j+1));
                masks[h]+=bit;
            }
            sum(slot+1,taken+bit,j,(f&16)?j:free_index,(parity%2)?-sign:sign);
            for(int h=0;h<4;++h) if((f>>h)&1) masks[h]-=bit;
        }
    }
};
int main() {
    int n,testcount,count,modular;
    if(!(std::cin>>n>>testcount>>count>>modular)||n>13||n<5) return 1;
    vector<vector<int>> shapes(count);
    for(auto& f:shapes) {
        int size;std::cin>>size;f.resize(size);
        for(int& x:f) std::cin>>x;
    }
    vector<Wave> waves(testcount);
    for(auto& w:waves) {
        w.p.resize(n);w.minors.assign(1<<n,0);
        for(auto& v:w.p) for(auto& x:v) std::cin>>x;
        for(int i=0;i<n;++i) for(int j=i+1;j<n;++j) for(int k=j+1;k<n;++k)
        for(int l=k+1;l<n;++l) for(int m=l+1;m<n;++m) {
            int mask=(1<<i)+(1<<j)+(1<<k)+(1<<l)+(1<<m);std::cin>>w.minors[mask];
            if(modular) w.minors[mask]%=modulus;
        }
    }
    for(const auto& f:shapes) {
        Orbit x(n,modular,f,waves);x.sum(0,0,0,-1,x.base_sign);
        for(const auto& row:x.total) {
            for(int h=0;h<(x.free_leg<0?1:5);++h) {
                Integer val=row[h];if(modular) val=(val%modulus+modulus)%modulus;
                write_integer(val);std::cout<<' ';
            }
        }
        std::cout<<'\n';
    }
}
