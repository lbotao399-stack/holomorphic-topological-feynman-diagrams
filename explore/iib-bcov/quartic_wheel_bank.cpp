// Modular finite angular sampling of K3,3 wheel faces, with cached tree messages.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <sstream>
#include <vector>
using U=uint32_t;
using L=uint64_t;
constexpr int width=12;
using List=std::array<U,width>;
int q,prime;
U mul(U a,U b){return L(a)*b%prime;}
U powmod(U a,int n){U r=1;for(;n;n>>=1){if(n&1)r=mul(r,a);a=mul(a,a);}return r;}
U inv(int a){return powmod((a%prime+prime)%prime,prime-2);}
U binom(int a,int b){U c=1;for(int j=1;j<=b;++j)c=mul(mul(c,a-j+1),inv(j));return c;}
struct Layout{
    std::array<int,3> masks;
    std::array<std::array<int,2>,3> crossings;
    bool star;
    int center,outer1,mid,outer2; // path outer1-center-mid-outer2, or a star centered at center.
    int h1,h2,h3;
};
std::vector<Layout> layouts(){
    std::array<int,4> order{0,1,2,3};
    std::vector<Layout> out;
    do{
        auto [i,j,k,l]=order;
        if((i-j)%2==0 && i>j)continue;
        int a=1<<i,b=1<<j,c=a|b|(1<<k);
        Layout row{};row.masks={a,(i-j)%2?a|b:b,c};
        int edge_id[4][4]{},degree[4]{};
        for(int h=0;h<3;++h){
            int mask=row.masks[h];int n=0;
            for(int v=0;v<4;++v)if(((mask>>v)^(mask>>((v+3)%4)))&1)row.crossings[h][n++]=v;
            if(n!=2)throw std::logic_error("noninteger radial root");
            auto [v,w]=row.crossings[h];edge_id[v][w]=edge_id[w][v]=h+1;++degree[v];++degree[w];
        }
        row.center=int(std::max_element(degree,degree+4)-degree);row.star=degree[row.center]==3;
        int count=0;
        for(int v=0;v<4;++v)if(degree[v]==1 && edge_id[row.center][v]){
            if(count==0){row.outer1=v;row.h1=edge_id[row.center][v]-1;}
            if(count==1){row.mid=v;row.h2=edge_id[row.center][v]-1;}
            if(count==2){row.outer2=v;row.h3=edge_id[row.center][v]-1;}
            ++count;
        }
        if(!row.star){
            row.mid=0;while(row.mid==row.center||degree[row.mid]!=2)++row.mid;
            row.outer2=6-row.center-row.outer1-row.mid;
            row.h2=edge_id[row.center][row.mid]-1;row.h3=edge_id[row.mid][row.outer2]-1;
        }
        out.push_back(row);
    }while(std::next_permutation(order.begin(),order.end()));
    if(out.size()!=20)throw std::logic_error("invalid ordering partition");
    return out;
}
struct Chamber{int a,b;};
struct Face{std::array<Chamber,20> terms;};
struct Tables{
    int slots, upper, span,stub_count,path_count;
    std::vector<U> leaf,step;
    std::vector<List> stub,other;
    std::vector<int> active_other,active_stubs,active_steps,lookup;
    std::vector<uint8_t> need_stub,need_step;
    explicit Tables(int degree):slots((degree+1)*q),upper(3*q-3),span(upper+q),
             stub_count(slots*slots*upper),path_count(stub_count*upper),leaf(slots*span),step(stub_count*span),
             stub(stub_count),lookup(2*path_count,-1),need_stub(stub_count),need_step(stub_count*span){}
    int sid(int center,int outside,int h)const{return (center*slots+outside)*upper+h;}
    int oid(bool star,int b,int c,int u,int v)const{return (star?path_count:0)+((b*slots+c)*upper+u)*upper+v;}
    void require_stub(int id){
        if(!need_stub[id]){need_stub[id]=1;active_stubs.push_back(id);}
    }
    Face face(L packed,const std::vector<Layout>& all){
        int a[3][3],zi=-1,zj=-1;
        for(int i=2;i>=0;--i)for(int j=2;j>=0;--j){a[i][j]=packed&15;packed>>=4;if(!a[i][j]){zi=i;zj=j;}}
        if(zi<0)throw std::logic_error("face contains no contracted line");
        int i=(zi+1)%3,k=(zi+2)%3,j=(zj+1)%3,l=(zj+2)%3;
        std::array<int,4> u{a[i][j],a[k][j],a[k][l],a[i][l]};
        std::array<int,4> b{a[i][zj],a[zi][j],a[k][zj],a[zi][l]};
        return sector(u,b,all,0,0,0);
    }
    Face sector(const std::array<int,4> &u,const std::array<int,4> &b,
                const std::vector<Layout>& all,int higher,int radial0,int radial2,int gap=1) {
        std::array<int,4> modes=b;
        for(int v=0;v<4;++v) if((higher>>v)&1) modes[v]+=gap*q;
        Face output{};
        for(unsigned t=0;t<all.size();++t){
            auto row=all[t];int h[3]{};
            for(int e=0;e<3;++e){
                int mask=row.masks[e];
                for(int v=0;v<4;++v)if(mask>>v&1){
                    h[e]+=q-u[v];if(mask>>((v+3)%4)&1)h[e]-=b[v];
                }
                int extra=(mask&1 ? radial0 : 0)+(mask&4 ? radial2 : 0);
                for(int v:row.crossings[e]) extra+=gap*((higher>>v)&1);
                if(extra%2) throw std::logic_error("half-integer radial exponent");
                h[e]+=extra/2;
                if (h[e]<=0 || h[e]>=upper) {
                    throw std::logic_error("nonconvergent radial chamber");
                }
            }
            int c1=sid(modes[row.center],modes[row.outer1],h[row.h1]);
            require_stub(c1);
            int c2=oid(row.star,modes[row.mid],modes[row.outer2],h[row.h2],h[row.h3]);
            if(lookup[c2]<0){
                lookup[c2]=int(active_other.size());
                active_other.push_back(c2);
                other.emplace_back();
                if(!row.star){
                    int id=sid(modes[row.mid],modes[row.outer2],h[row.h3]);
                    require_stub(id);
                    for(int k=0;k<q-1;++k){
                        int address=id*span+h[row.h2]+k;
                        if(!need_step[address]){need_step[address]=1;active_steps.push_back(address);}
                    }
                }
            }
            output.terms[t]={c1,lookup[c2]};
        }
        return output;
    }
    U prepare(int n){
        std::vector<U> reciprocal(n+6*q+3),fact(n+3*q+1,1);
        for(unsigned j=1;j<reciprocal.size();++j)reciprocal[j]=inv(int(j));
        for(unsigned j=1;j<fact.size();++j)fact[j]=mul(j,fact[j-1]);
        std::vector<List> weights(slots);
        const int lam=q-1;
        for(int id=1;id<slots;++id){
            const int b=id%q, degree=n+id/q;
            if(!b) continue;
            U v=mul(binom(b+degree-1,degree),inv(binom(lam+degree-1,degree)));
            for(int k=0;k<q-b;++k){
                weights[id][k]=v;
                if(k<lam-b)v=mul(mul(mul(mul(v,b+degree+k),prime-(lam-b-k)),
                                       reciprocal[k+1]),reciprocal[degree+lam+1+k]);
            }
        }
        for(int b=1;b<slots;++b)if(b%q)for(int s=1;s<span;++s){
            L sum=0;for(int k=0;k<q-(b%q);++k)sum+=L(weights[b][k])*reciprocal[n+s+k];
            leaf[b*span+s]=sum%prime;
        }
        for(int id:active_stubs){
            int h=id%upper,c=(id/upper)%slots,b=id/(upper*slots);
            for(int k=0;k<q-(b%q);++k)stub[id][k]=mul(weights[b][k],leaf[c*span+h+k]);
        }
        for(int address:active_steps){
            int id=address/span,s=address%span,b=id/(upper*slots);
            L sum=0;for(int k=0;k<q-(b%q);++k)sum+=L(stub[id][k])*reciprocal[n+s+k];
            step[address]=sum%prime;
        }
        for(int id:active_other){
            int index=id%path_count,v=index%upper;index/=upper;
            int u=index%upper;index/=upper;
            int c=index%slots,b=index/slots;
            for(int k=0;k<q-1;++k)
                other[lookup[id]][k]=id>=path_count ? mul(leaf[b*span+u+k],leaf[c*span+v+k])
                                             : step[sid(b,c,v)*span+u+k];
        }
        U g=mul(powmod(lam,3),inv(powmod(fact[lam],3)));
        g=mul(g,binom(2*lam-1+n,n));return mul(g,inv(powmod(n+lam,3)));
    }
};
// The four shift indices form a path or a star after the radial order is fixed.
U shell(const Tables &data,const Face &face) {
    L value=0;
    for(auto t:face.terms) {
        const List& a=data.stub[t.a];const List& b=data.other[t.b];
        L dot=0;for(int k=0;k<q-1;++k) dot+=L(a[k])*b[k];
        value+=dot%prime;
    }
    return value%prime;
}

struct Transition {
    int radial0,jump,mask;
};
// Each pair (radial degree at endpoint zero, spin jump) has a common radial
// integral. Sum all multinomial terms feeding it before multiplying the shell.
std::vector<Transition> transitions(int degree) {
    std::vector<Transition> out;
    for(int r=0;r<=2*degree;++r) {
        int top=std::min(r,2*degree-r);
        for(int j=r%2;j<=top;j+=2) {
            if(!j) out.push_back({r,0,0});
            else for(int mask : {6,9}) out.push_back({r,j,mask});
        }
    }
    return out;
}

std::vector<U> transition_weights(unsigned n,int degree,const std::vector<Transition>& sectors) {
    const int lam=q-1,center=degree;
    std::vector<std::vector<U>> powers(degree+1,std::vector<U>(2*degree+1));
    powers[0][center]=1;
    for(int p=0;p<degree;++p) for(int h=-p;h<=p;h+=2) {
        const int spin=int(n)+h;
        if (spin<0) continue;
        const U weight=powers[p][center+h];
        const U factor=mul(weight,inv(spin+lam));
        U &up=powers[p+1][center+h+1];
        up=(up+prime-mul(factor,spin+1))%prime;
        if(spin) {
            U &down=powers[p+1][center+h-1];
            down=(down+prime-mul(factor,spin+2*lam-1))%prime;
        }
    }
    std::vector<U> ratio(degree+1,1);
    U norm=1;
    for(int j=1;j<=degree;++j){
        norm=mul(norm,mul(n+2*lam+j-1,inv(n+j)));
        ratio[j]=mul(norm,powmod(mul(n+lam,inv(n+lam+j)),2));
    }
    std::vector<U> result;
    for(auto sector:sectors){
        int r=sector.radial0,j=sector.jump,top=std::min(r,2*degree-r);
        L coefficient=0;
        for(int p=j;p<=top;p+=2) {
            const int i=(r-p)/2;
            U choose=mul(binom(degree,p),binom(degree-p,i));
            coefficient+=mul(choose,powers[p][center+j]);
        }
        result.push_back(mul(coefficient%prime,ratio[j]));
    }
    return result;
}

int main(int argc,char**argv){
    const std::string option=argc==7 ? argv[6] : "";
    const bool fixed_spin=option.rfind("--fixed-spin=",0)==0;
    const bool general=fixed_spin || option.rfind("--degree=",0)==0;
    const int chord=general ? std::stoi(option.substr(fixed_spin ? 13 : 9)) : (option=="--quadratic" ? 2 : int(option=="--signed"));
    if(argc!=6 && !chord && !general)throw std::invalid_argument("q prime kernel-file wheel-keys output-table [--signed|--quadratic|--degree=N]");
    q=std::stoi(argv[1]);prime=std::stoi(argv[2]);
    if(q>width+1 || q<2 || prime>1000000007 || prime<2)throw std::invalid_argument("wheel field or dimension");
    if(general && (chord<0 || chord>=q)) throw std::invalid_argument("invalid numerator degree");
    Tables data(chord);
    auto sectors=general ? transitions(chord) : std::vector<Transition>{};
    auto all=layouts();
    std::ifstream file(argv[4]);
    std::vector<Face> faces;
    for(L key;file>>key;) {
        if (!chord && !general) {faces.push_back(data.face(key,all));continue;}
        std::array<int,4> u{},b{};
        for(int i=3;i>=0;--i){b[i]=key&15;key>>=4;}
        for(int i=3;i>=0;--i){u[i]=key&15;key>>=4;}
        if(general){
            for(auto v:sectors) faces.push_back(data.sector(u,b,all,v.mask,v.radial0,2*chord-v.radial0,v.jump));
        } else if (chord==1) {
            faces.push_back(data.sector(u,b,all,0,2,0));
            faces.push_back(data.sector(u,b,all,0,0,2));
            faces.push_back(data.sector(u,b,all,6,1,1));
            faces.push_back(data.sector(u,b,all,9,1,1));
        } else {
            faces.push_back(data.sector(u,b,all,0,4,0));
            faces.push_back(data.sector(u,b,all,0,0,4));
            faces.push_back(data.sector(u,b,all,0,2,2));
            for (int mask : {6,9}) {
                faces.push_back(data.sector(u,b,all,mask,3,1));
                faces.push_back(data.sector(u,b,all,mask,1,3));
                faces.push_back(data.sector(u,b,all,mask,2,2,2));
            }
        }
    }
    if(faces.empty())throw std::runtime_error("empty wheel list");
    int channels=general ? 0 : (chord==2 ? 8 : 7);
    std::vector<std::vector<U>> kernel;
    std::ifstream input(argv[3]);
    std::string line;
    while(std::getline(input,line)){
        std::istringstream row(line);
        std::vector<U> projection;
        for(U value;row>>value;) projection.push_back(value);
        if(!channels) channels=int(projection.size());
        if(!channels || int(projection.size())!=channels) throw std::runtime_error("inconsistent kernel columns");
        kernel.push_back(projection);
    }
    if(kernel.empty())throw std::runtime_error("empty interpolation kernel");
    std::cerr<<"faces "<<faces.size()<<" nodes "<<kernel.size()<<" active messages "<<data.active_other.size()<<std::endl;
    const int stride=general ? int(sectors.size()) : (chord==2 ? 9 : (chord ? 4 : 1));
    std::vector<std::vector<U>> values(faces.size()/stride,std::vector<U>(channels));
    for(unsigned n=0;n<kernel.size()+(fixed_spin ? chord : 0);++n) {
        U factor=data.prepare(n);
        const int lam=q-1;
        U ratio=chord ? prime-mul(mul(n+lam,n+2*lam),powmod(inv(n+lam+1),2)) : 0;
        U diagonal=0,second=0;
        if (chord==2) {
            U up=mul(mul(n+1,n+2*lam),inv(n+lam+1));
            U down=n ? mul(mul(n,n+2*lam-1),inv(n+lam-1)) : 0;
            diagonal=(2+mul((up+down)%prime,inv(n+lam)))%prime;
            second=mul(mul(mul(n+lam,n+2*lam),n+2*lam+1),
                       mul(inv(n+lam+1),powmod(inv(n+lam+2),2)));
        }
        auto factors=general ? transition_weights(n,chord,sectors) : std::vector<U>{};
        if(!fixed_spin)for(auto& w:kernel[n])w=mul(w,factor);
        for(unsigned i=0;i<values.size();++i) {
            const unsigned j=stride*i;
            if(fixed_spin) {
                // Regulate the angular momentum m on rim edge zero. The tail
                // interpolation variable is m-chord, and the first chord
                // shells contribute finite rational terms.
                std::vector<U> grouped(chord+1);
                for(int k=0;k<stride;++k) {
                    int shift=(sectors[k].mask&1) ? sectors[k].jump : 0;
                    grouped[shift]=(grouped[shift]+mul(factors[k],shell(data,faces[j+k])))%prime;
                }
                for(int shift=0;shift<=chord;++shift) {
                    const int index=int(n)+shift-chord;
                    U term=mul(factor,grouped[shift]);
                    if(index<0) values[i][0]=(values[i][0]+term)%prime;
                    else if(unsigned(index)<kernel.size()) for(int h=0;h<channels;++h)
                        values[i][h]=(values[i][h]+L(term)*kernel[index][h])%prime;
                }
                continue;
            }
            U value=0;
            if (general) {
                L total=0;
                for (int k=0;k<stride;++k) total+=mul(factors[k],shell(data,faces[j+k]));
                value=total%prime;
            } else value=shell(data,faces[j]);
            if (!general && chord==1) {
                U cross=(shell(data,faces[j+2])+shell(data,faces[j+3]))%prime;
                value=(L(value)+shell(data,faces[j+1])+mul(ratio,cross))%prime;
            }
            if (!general && chord==2) {
                L linear=0,quadratic=0;
                for (unsigned k : {j+3,j+6}) {
                    linear+=shell(data,faces[k])+L(shell(data,faces[k+1]));
                    quadratic+=shell(data,faces[k+2]);
                }
                value=(L(value)+shell(data,faces[j+1])+mul(diagonal,shell(data,faces[j+2]))+
                       mul((2*L(ratio))%prime,linear%prime)+mul(second,quadratic%prime))%prime;
            }
            for(int h=0;h<channels;++h)values[i][h]=(values[i][h]+L(value)*kernel[n][h])%prime;
        }
        if((n+1)%25==0)std::cerr<<"sample "<<n+1<<std::endl;
    }
    std::ofstream out(argv[5]);
    for(const auto &v:values) {
        if(chord || general) for(int j=0;j<channels;++j) out<<v[j]<<(j+1==channels?'\n':' ');
        else out<<v[0]<<' '<<v[3]<<' '<<v[5]<<'\n';
    }
}
