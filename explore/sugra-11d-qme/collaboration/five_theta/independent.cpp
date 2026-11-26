// Exact positive-expansion Dirichlet integration of five 11d two-loop theta graphs.
// The adjacent replay.py builds this executable with the installed C++ compiler and GMP.
// independent p q r block_file ...  (only complete "rank nterms" blocks)
#include <gmpxx.h>
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>

using Word=unsigned long;
using Wide=unsigned __int128;
using Trio=std::array<int,3>;
using Labels=std::array<int,10>;
using Axes=std::array<int,5>;
constexpr int max_degree=13;
const std::array<Trio,5> topologies={{{1,6,6},{2,5,6},{3,4,6},{3,5,5},{4,4,5}}};

Word fact(int n) {
    Word ret=1;
    for(int i=2;i<=n;++i) ret*=i;
    return ret;
}
Word step2(int x,int n) {
    Word ret=1;
    for(int i=0;i<n;++i) ret*=x+2*i;
    return ret;
}

struct Topology {
    Trio edges,interiors;
    Labels path{};
    std::vector<Labels> matchings;
    Word wden=step2(13,max_degree);
    Word bins[14][14]{}, weight[14][14][14][14]{};
    mpz_class scale;
    mpq_class scalar;
    Topology(Trio lengths):edges(lengths) {
        if(std::find(topologies.begin(),topologies.end(),edges)==topologies.end())
            throw std::invalid_argument("unknown theta graph");
        int pos=0;
        scale=wden;
        Word kernel_num=1,kernel_den=step2(1,6);
        for(int p=0;p<3;++p) {
            int n=interiors[p]=edges[p]-1;
            for(int i=0;i<n;++i) path[pos++]=p;
            scale*=fact(n+max_degree)/fact(n);
            kernel_num*=step2(1,6-edges[p]);
            kernel_den*=fact(n);
        }
        scalar=mpq_class(kernel_num,kernel_den);scalar.canonicalize();
        Labels lab{};
        Axes counts{},first{};
        enumerate(lab,counts,first,0,0);
        const Word expected=32*fact(5)*fact(interiors[0])*fact(interiors[1])*fact(interiors[2])
            /(fact(6-edges[0])*fact(6-edges[1])*fact(6-edges[2])*32*fact(5));
        if(matchings.size()!=expected) throw std::logic_error("wrong match count");
        for(int n=0;n<=max_degree;++n) {
            bins[n][0]=bins[n][n]=1;
            for(int k=1;k<n;++k) bins[n][k]=bins[n-1][k-1]+bins[n-1][k];
        }
        for(int a=0;a<=max_degree;++a) for(int b=0;a+b<=max_degree;++b)
        for(int c=0;a+b+c<=max_degree;++c) for(int d=0;a+b+c+d<=max_degree;++d) {
            int n=a+b+c+d;
            Wide val=0;
            // E[(wA+wC)^a wB^b wC^c (wA+wB)^d]. All terms nonnegative.
            for(int i=0;i<=a;++i) for(int j=0;j<=d;++j)
                val+=Wide(bins[a][i])*bins[d][j]*step2(13-2*edges[0],a+d-i-j)
                    *step2(13-2*edges[1],b+j)*step2(13-2*edges[2],c+i);
            val*=step2(13+2*n,max_degree-n);
            if(val>wden) throw std::logic_error("Dirichlet moment exceeds one");
            weight[a][b][c][d]=static_cast<Word>(val);
        }
    }
    void enumerate(Labels& lab,Axes& count,Axes& first,int pos,int next) {
        if(pos==10){matchings.push_back(lab);return;}
        for(int i=0;i<5 && i<=next;++i) {
            if(count[i]==2 || (count[i] && first[i]==path[pos])) continue;
            if(count[i]==0) first[i]=path[pos];
            ++count[i];lab[pos]=i;
            enumerate(lab,count,first,pos+1,std::max(next,i+1));
            --count[i];
        }
    }
};

struct Block {
    const Topology& graph;
    Axes first{},second{},b{},c{};
    Labels powers{};
    std::array<std::vector<Word>,3> cached;
    mpz_class acc,term;
    Block(const Topology& g,int rank):graph(g) {
        first.fill(-1);
        for(int j=0;j<10;++j) {
            const int axis=g.matchings.at(rank)[j];
            (first[axis]<0 ? first[axis] : second[axis])=j;
        }
        for(int p=0;p<3;++p) cached[p].resize(1UL<<(4*g.interiors[p]));
    }
    Word pmoment(int p,int start) {
        unsigned key=0;int n=graph.interiors[p];
        for(int j=0;j<n;++j) key|=unsigned(powers[start+j])<<(4*j);
        Word& result=cached[p][key];
        if(result)return result;
        Word den=1;int tail=0;
        for(int j=n-1;j>=0;--j) {tail+=powers[start+j];den*=n-j+tail;}
        // The denominators form a subset of 1,...,n+max_degree.
        const Word scale_n=fact(n+max_degree);
        if(scale_n%den)throw std::logic_error("ordered tails moment denominator");
        result=scale_n/den;
        return result;
    }
    void choose(int i,int comp_b,int plain_b,int plain_c,int comp_c,int neg,Word mult) {
        if(i==5) {
            term=graph.weight[comp_b][plain_b][plain_c][comp_c];
            int pos=0;
            for(int p=0;p<3;++p){term*=pmoment(p,pos);pos+=graph.interiors[p];}
            if(neg%2) mpz_submul_ui(acc.get_mpz_t(),term.get_mpz_t(),mult);
            else mpz_addmul_ui(acc.get_mpz_t(),term.get_mpz_t(),mult);
            return;
        }
        const int plus=first[i],minus=second[i];
        const int u=graph.path[plus],v=graph.path[minus];
        for(int j=0;j<=b[i];++j) for(int k=0;k<=c[i];++k) {
            // First wave is +e_i, second wave -e_i. For Q=B,C,
            // ell_Q component = (w_Q-delta_Q,U)*S + (delta_Q,V-w_Q)*T.
            // S,T are the two tail fractions. j and k count S choices.
            powers[plus]=j+k;
            powers[minus]=b[i]+c[i]-j-k;
            const int cb=(u==1 ? j:0)+(v==1 ? b[i]-j:0);
            const int cc=(u==2 ? k:0)+(v==2 ? c[i]-k:0);
            const int local_neg=(u==1 ? j:0)+(v!=1 ? b[i]-j:0)
                               +(u==2 ? k:0)+(v!=2 ? c[i]-k:0);
            choose(i+1,comp_b+cb,plain_b+b[i]-cb,plain_c+c[i]-cc,comp_c+cc,
                   neg+local_neg,mult*graph.bins[b[i]][j]*graph.bins[c[i]][k]);
        }
    }
    mpz_class monomial(Word key) {
        for(int i=0;i<5;++i) {b[i]=(key>>(4*i))&15;c[i]=(key>>(4*(i+5)))&15;}
        if(std::accumulate(b.begin(),b.end(),0)+std::accumulate(c.begin(),c.end(),0)>max_degree)
            throw std::logic_error("unsupported moment degree");
        acc=0;choose(0,0,0,0,0,0,1);return acc;
    }
};

int main(int argc,char**argv) {
    if(argc<5){std::cerr<<"independent p q r blocks...\n";return 1;}
    Topology g({std::stoi(argv[1]),std::stoi(argv[2]),std::stoi(argv[3])});
    if(Block(g,0).monomial(0)!=g.scale) return 2;
    std::cerr<<"matchings "<<g.matchings.size()<<" unit-kernel/pi "<<g.scalar<<"\n";
    mpz_class total=0;
    for(int i=4;i<argc;++i) {
        std::ifstream in(argv[i]);
        if(!in) throw std::runtime_error("cannot open block file");
        int rank,n;
        while(in>>rank>>n) {
            if(n<0 || rank<0 || size_t(rank)>=g.matchings.size())return 3;
            Block block(g,rank);mpz_class subtotal=0;
            for(int t=0;t<n;++t) {
                Word key;mpz_class coef;
                if(!(in>>key>>coef)) throw std::runtime_error("truncated block");
                subtotal+=coef*block.monomial(key);
            }
            total+=subtotal;
            mpq_class mean(subtotal,g.scale);mean.canonicalize();
            std::cout<<rank<<" "<<n<<" "<<mean<<" "<<mean*g.scalar<<"\n";
        }
        if(!in.eof())throw std::runtime_error("invalid input");
    }
    mpq_class mean(total,g.scale);mean.canonicalize();
    std::cout<<"SUM "<<mean<<"\nKERNEL_OVER_PI "<<mean*g.scalar<<"\n";
}
