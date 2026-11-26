// Integrate the first-edge currents for all five theta graphs by exact integer moments.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <gmpxx.h>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <unordered_map>

using u64=std::uint64_t;
using u128=unsigned __int128;
using i128=__int128;
using Word=std::array<int,10>;
using Shape=std::array<int,3>;
constexpr u64 wscale=21322087106945625ULL;  // 13*15*...*35
constexpr unsigned mask=(1u<<20)-1;

std::vector<Word> matchings(const Shape& n) {
    std::vector<Word> list;
    Word word{}, seen{}; std::array<unsigned,3> occupied{};
    auto extend=[&](const auto& self, int j, int labels)->void {
        if(j==10) { if(labels==5) list.push_back(word); return; }
        const int path=j<n[0]-1 ? 0 : j<n[0]+n[1]-2 ? 1 : 2;
        for(int a=0;a<std::min(5,labels+1);++a) {
            if(seen[a]==2 || (occupied[path]&(1u<<a))) continue;
            word[j]=a; ++seen[a]; occupied[path]^=1u<<a;
            self(self,j+1,std::max(labels,a+1));
            --seen[a]; occupied[path]^=1u<<a;
        }
    };
    extend(extend,0,0); return list;
}

u64 factorial(unsigned n){u64 f=1;for(unsigned i=2;i<=n;++i) f*=i;return f;}
u64 odd_product(unsigned first,unsigned n){u64 f=1;for(unsigned i=0;i<n;++i)f*=first+2*i;return f;}
unsigned deg(unsigned x) {
    return (x&15)+((x>>4)&15)+((x>>8)&15)+((x>>12)&15)+(x>>16);
}
unsigned choose(unsigned n,unsigned k) {
    unsigned a=1;for(unsigned i=1;i<=k;++i) a=a*(n+1-i)/i;return a;
}

struct Split{unsigned code,degree,coef;};
std::vector<Split> split(unsigned x) {
    std::vector<Split> result={{0,0,1}};
    for(unsigned i=0;i<5;i++){
        int n=(x>>(4*i))&15;size_t old=result.size();
        for(int k=1;k<=n;k++) for(size_t j=0;j<old;j++){
            auto s=result[j];result.push_back({s.code+((unsigned)k<<(4*i)),s.degree+k,s.coef*choose(n,k)});
        }
    }
    return result;
}

struct Path {
    std::vector<int> axes;
    std::unordered_map<unsigned,u64> cache;
    unsigned plus=0, minus=0;
    u64 fact=1;
    u64 scale=1;
    void init() {
        if(!axes.empty()) {
            fact=factorial(axes.size()+12);
            scale=fact/factorial(axes.size());
            // Nonnegative path exponents have total at most twelve.
            cache.reserve(choose(static_cast<unsigned>(axes.size())+12, static_cast<unsigned>(axes.size())));
        }
    }
    u64 moment(unsigned s,unsigned t) {
        if(axes.empty()) return 1;
        unsigned x=(s&plus)|(t&minus);
        auto found=cache.find(x);
        if(found!=cache.end()) return found->second;
        unsigned tail=0;u64 denom=1;
        for(size_t j=0;j<axes.size();j++) {
            tail+=(x>>(4*axes[axes.size()-1-j]))&15;
            denom*=1+j+tail;
        }
        if(fact%denom) throw std::logic_error("nonintegral path scale");
        return cache.emplace(x,fact/denom).first->second;
    }
};
struct Rule {
    std::array<Path,3> path;
    long weights[13][13][13][13]={};
    unsigned alphaB,alphaC; // twice the two Dirichlet parameters
    Rule(const Shape& n, const Word& word):alphaB(13-2*n[1]),alphaC(13-2*n[2]) {
        unsigned seen=0;int j=0;
        for(int p=0;p<3;++p) {
            for(int i=0;i<n[p]-1;++i,++j) {
                unsigned a=word[j], bit=1u<<a;
                path[p].axes.push_back(a);
                (seen&bit?path[p].minus:path[p].plus)|=15u<<(4*a);
                seen|=bit;
            }
            path[p].init();
        }
    }
    long weight(unsigned a,unsigned b,unsigned c,unsigned d) {
        // E[(wB-1)^a wB^b (wC-1)^c wC^d] times common scale.
        long& cached=weights[a][b][c][d];if(cached)return cached;
        i128 value=0;
        for(unsigned i=0;i<=a;++i)for(unsigned j=0;j<=c;++j){
            unsigned k=i+b, l=j+d;
            i128 term=(i128)(wscale/odd_product(13,k+l))*odd_product(alphaB,k)*odd_product(alphaC,l);
            term*=choose(a,i)*choose(c,j);
            value+=(a-i+c-j)%2?-term:term;
        }
        if(value>(i128)wscale || value<-(i128)wscale) throw std::logic_error("conductance moment out of bounds");
        return cached=(long)value;
    }
    mpz_class num(u64 key) {
        const unsigned b=key&mask,c=key>>20,k=deg(b),l=deg(c);
        if(k+l>12)throw std::logic_error("excessive current degree");
        i128 groups[13][13]={};
        auto us=split(b),vs=split(c);
        for(auto& u:us) for(auto& v:vs) {
            // u and v count the powers assigned to the negative momentum in each axis.
            unsigned s=b-u.code+c-v.code, t=u.code+v.code;
            unsigned a=deg(((b-u.code)&path[1].plus)|(u.code&path[1].minus));
            unsigned h=deg(((c-v.code)&path[2].plus)|(v.code&path[2].minus));
            // Bound: 2^12 times the three path scales is < 2^127 for all five shapes.
            i128 term=(i128)path[0].moment(s,t)*path[1].moment(s,t)*path[2].moment(s,t)*u.coef*v.coef;
            groups[a][h]+=((u.degree+v.degree)%2)?-term:term;
        }
        mpz_class answer;
        for(unsigned a=0;a<=k;a++)for(unsigned h=0;h<=l;h++)if(groups[a][h]){
            u128 mag=groups[a][h]<0?-groups[a][h]:groups[a][h];
            u64 words[]={(u64)mag,(u64)(mag>>64)};mpz_class n;
            mpz_import(n.get_mpz_t(),2,-1,sizeof(u64),0,0,words);
            long w=weight(a,k-a,h,l-h);
            answer+=n*(groups[a][h]<0?-w:w);
        }
        return answer;
    }
};
int scalar_denominator(Shape n) {
    const Shape shapes[]={{1,6,6},{2,5,6},{3,4,6},{3,5,5},{4,4,5}};
    int den[]={158400,285120,332640,798336,997920};
    for(int i=0;i<5;i++)if(n==shapes[i])return den[i];
    throw std::runtime_error("shape not in the two-loop catalogue");
}
int main(int argc,char**argv){
    if(argc<5){std::cerr<<"Usage: moments p q r chunks...\n";return 2;}
    Shape shape={std::stoi(argv[1]),std::stoi(argv[2]),std::stoi(argv[3])};
    int scalar=scalar_denominator(shape);
    auto list=matchings(shape);
    mpz_class sum,scale;
    std::vector<bool> seen(list.size());
    int nblocks=0;
    for(int f=4;f<argc;f++){
        std::ifstream in(argv[f]);if(!in) return 2;
        int rank;size_t nterms;
        while(in>>rank>>nterms){
            if(rank<0 || (size_t)rank>=list.size() || seen[rank])throw std::runtime_error("invalid rank");
            seen[rank]=true; Rule rule(shape,list[rank]);
            scale=(unsigned long)wscale;
            for(auto& p:rule.path) scale*=(unsigned long)p.scale;
            scale*=scalar;
            mpz_class subtotal;
            for(size_t i=0;i<nterms;i++){
                u64 key; long coefficient;
                if(!(in>>key>>coefficient)) throw std::runtime_error("truncated polynomial");
                subtotal += rule.num(key)*coefficient;
            }
            sum+=subtotal;++nblocks;
            mpq_class result(subtotal,scale);result.canonicalize();
            std::cout<<"rank "<<rank<<" terms "<<nterms<<" kernel/pi "<<result<<"\n";
        }
        if(!in.eof())throw std::runtime_error("parse error");
    }
    if(!nblocks) return 2;
    mpq_class result(sum,scale);result.canonicalize();
    std::cout<<"ranks "<<nblocks<<" signed kernel/pi "<<result<<"\n";
}
