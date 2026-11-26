// Exact integer integration of canonical-order current polynomials.
// Build and verify with replay.py in this directory.
// Input: blocks "rank number_of_terms", then "key signed_integer_coefficient".
#include <gmpxx.h>
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>

using Order = std::array<int,5>;
using Word = unsigned long;
using Wide = unsigned __int128;
constexpr Word factorial18=6402373705728000UL;
constexpr Word path_den=factorial18/120;
constexpr Word weight_den=788917222956988125UL;
Word binom[14][14],weights[14][14][14][14];
std::array<Order,120> orders;

Word step2(int start,int n) {
    Word ans=1;
    for(int i=0;i<n;++i) ans*=start+2*i;
    return ans;
}

void prepare() {
    Order order={0,1,2,3,4};
    for(auto& x:orders){x=order; std::next_permutation(order.begin(),order.end());}
    for(int n=0;n<=13;++n) {
        binom[n][0]=binom[n][n]=1;
        for(int j=1;j<n;++j) binom[n][j]=binom[n-1][j-1]+binom[n-1][j];
    }
    for(int a=0;a<=13;++a) for(int b=0;a+b<=13;++b)
    for(int c=0;a+b+c<=13;++c) for(int d=0;a+b+c+d<=13;++d) {
        // |moment| = E[(wA+wC)^a wB^b wC^c (wA+wB)^d].
        // i and j count factors wC and wB selected from the binomials.
        const int degree=a+b+c+d;
        Wide sum=0;
        for(int i=0;i<=a;++i) for(int j=0;j<=d;++j)
            sum+=Wide(binom[a][i])*binom[d][j]*step2(11,a+d-i-j)
                 *step2(1,b+j)*step2(1,c+i);
        sum*=step2(13+2*degree,13-degree);
        if(sum>weight_den) throw std::logic_error("conductance moment exceeds one");
        weights[a][b][c][d]=static_cast<Word>(sum);
    }
}

Word path_scaled(const Order& powers,const Order& order) {
    int sum=0; Word den=1;
    for(int j=4;j>=0;--j) {
        sum+=powers[order[j]];
        den*=5-j+sum;
    }
    if(factorial18%den) throw std::logic_error("ordered ratio denominator");
    return factorial18/den;
}

struct Integral {
    Order b={},c={},sp={},tp={};
    const Order& order;
    mpz_class subtotal;
    mpz_class term;
    Integral(const Order& o):order(o){}

    void choose(int i,int bs,int bt,int cs,int ct,Word multiplicity) {
        if(i==5) {
            term=weights[bs][bt][cs][ct];
            term*=path_scaled(sp,orders[0]);
            term*=path_scaled(tp,order);
            mpz_addmul_ui(subtotal.get_mpz_t(),term.get_mpz_t(),multiplicity);
            return;
        }
        // b_i=-(1-wB) S_i - wB T_i ; c_i=wC S_i+(1-wC) T_i.
        for(int j=0;j<=b[i];++j) for(int k=0;k<=c[i];++k) {
            sp[i]=j+k;tp[i]=b[i]+c[i]-j-k;
            choose(i+1,bs+j,bt+b[i]-j,cs+k,ct+c[i]-k,
                   multiplicity*binom[b[i]][j]*binom[c[i]][k]);
        }
    }

    mpz_class monomial(Word key) {
        for(int i=0;i<5;++i) {b[i]=(key>>(4*i))&15; c[i]=(key>>(4*(i+5)))&15;}
        int nb=std::accumulate(b.begin(),b.end(),0);
        int nc=std::accumulate(c.begin(),c.end(),0);
        if(nb+nc>13) throw std::logic_error("degree greater than 13");
        subtotal=0;
        choose(0,0,0,0,0,1);
        return nb%2 ? -subtotal:subtotal;
    }
};

int main(int argc,char**argv) {
    prepare();
    mpz_class den=mpz_class(path_den)*path_den*weight_den;
    if(Integral(orders[0]).monomial(0)!=den) return 3;
    if(argc<2) {std::cerr<<"independent_moments block_file(s)\n"; return 1;}
    mpz_class total=0;
    for(int f=1;f<argc;++f) {
        std::ifstream in(argv[f]);
        if(!in) throw std::runtime_error("cannot open aggregate");
        int rank,n;
        while(in>>rank>>n) {
            if(rank<0||rank>=120||n<0) return 2;
            Integral integrator(orders[rank]);
            mpz_class sum=0;
            for(int i=0;i<n;++i) {
                Word key; mpz_class coefficient;
                if(!(in>>key>>coefficient)) throw std::runtime_error("truncated block");
                sum+=coefficient*integrator.monomial(key);
            }
            total+=sum;
            mpq_class value(sum,den); value.canonicalize();
            std::cout<<rank<<" "<<n<<" "<<value<<"\n";
        }
        if(!in.eof()) throw std::logic_error("invalid input");
    }
    mpq_class scalar(total,den); scalar.canonicalize();
    std::cout<<"SUM "<<scalar<<"\nKERNEL_OVER_PI "<<scalar/158400<<"\n";
}
