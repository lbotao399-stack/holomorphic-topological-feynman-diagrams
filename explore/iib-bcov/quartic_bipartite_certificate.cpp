// Rational reconstruction at free direct rows, exact triangular evaluation, all-row check.
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#include <gmpxx.h>

using Key=uint64_t;
using Triple=std::array<mpq_class,3>;
struct Equation {
    int ids[7]{};
    int8_t coefficients[7]{};
    int8_t count=0;
};
template<class Function>
void equations(const std::string& base,Function visit){
    std::ifstream input(base+".bin",std::ios::binary);
    if(!input)throw std::runtime_error("missing orbit identities");
    int n;
    while((n=input.get())!=EOF){
        if(n<1||n>7)throw std::runtime_error("invalid row size");
        char bytes[35]; input.read(bytes,5*n);
        if(!input)throw std::runtime_error("truncated orbit identity");
        Equation row;
        for(int i=0;i<n;++i){
            int32_t id;
            std::memcpy(&id,bytes+5*i,4);
            int8_t c=bytes[5*i+4];
            if(c){row.ids[row.count]=id;row.coefficients[row.count++]=c;}
        }
        if(row.count)visit(row);
    }
}
std::vector<Key> keys(const std::string& path){
    std::ifstream input(path); if(!input)throw std::runtime_error("missing orbit keys");
    std::vector<Key> out;
    for(Key key; input>>key;)out.push_back(key);
    return out;
}
void read_rationals(std::istream& s,Triple& row){
    for(auto& x:row)if(!(s>>x))throw std::runtime_error("truncated rational row");
}
unsigned residue(const mpq_class& x,unsigned p){
    mpz_class inverse;
    mpz_class modulus=p;
    if(!mpz_invert(inverse.get_mpz_t(),x.get_den_mpz_t(),modulus.get_mpz_t()))
        throw std::runtime_error("field meets a denominator");
    return uint64_t(mpz_fdiv_ui(x.get_num_mpz_t(),p))*inverse.get_ui()%p;
}
bool rational(mpq_class& out,mpz_class x,const mpz_class& modulus){
    mpz_class bound=modulus/2;
    mpz_sqrt(bound.get_mpz_t(),bound.get_mpz_t());
    mpz_class a=modulus,t=0,u=1,quot,rem;
    while(x>bound){
        mpz_tdiv_qr(quot.get_mpz_t(),rem.get_mpz_t(),a.get_mpz_t(),x.get_mpz_t());
        a.swap(x); x.swap(rem);
        rem=t-quot*u; t.swap(u); u.swap(rem);
    }
    if(!u||abs(u)>bound)return false;
    mpz_gcd(rem.get_mpz_t(),x.get_mpz_t(),u.get_mpz_t());
    if(rem!=1)return false;
    out=mpq_class(x,u);out.canonicalize();return true;
}

int main(int argc,char** argv){
    if(argc!=7)throw std::invalid_argument("equation basename, exact wheels, star anchors or -, manifest, channel 0/1/2, output");
    const std::string base=argv[1]; const int ch=std::stoi(argv[5]);
    if(ch<0||ch>2)throw std::invalid_argument("invalid channel");
    auto nodes=keys(base+".nodes"), wheels=keys(base+".wheel");
    std::vector<Equation> pivots(nodes.size());
    equations(base,[&](const Equation& row){
        int lead=*std::max_element(row.ids,row.ids+row.count);
        if(lead<0)return;
        if(!pivots[lead].count || row.count<pivots[lead].count)pivots[lead]=row;
    });
    std::vector<mpq_class> value(nodes.size()+wheels.size());
    std::vector<bool> fixed(nodes.size());
    std::ifstream boundary(argv[2]); Triple row;
    if(!boundary)throw std::runtime_error("missing exact contact table");
    mpz_class denominator=1;
    auto set=[&](size_t i,const mpq_class& x){
        value[i]=x;
        mpz_lcm(denominator.get_mpz_t(),denominator.get_mpz_t(),x.get_den_mpz_t());
    };
    for(size_t i=0;i<wheels.size();++i){
        Key key; if(!(boundary>>key)||key!=wheels[i])throw std::runtime_error("contact key mismatch");
        read_rationals(boundary,row);set(nodes.size()+i,row[ch]);
    }
    if(std::string(argv[3])!="-"){
        std::ifstream stars(argv[3]);if(!stars)throw std::runtime_error("missing star anchors");
        for(Key key;stars>>key;){
            read_rationals(stars,row);auto it=std::lower_bound(nodes.begin(),nodes.end(),key);
            if(it==nodes.end()||*it!=key)throw std::logic_error("unknown star");
            size_t i=it-nodes.begin();fixed[i]=true;set(i,row[ch]);
        }
    }
    std::vector<int> free;
    for(int i=0;i<int(nodes.size());++i)if(!fixed[i]&&!pivots[i].count)free.push_back(i);
    std::vector<mpz_class> lift(free.size());
    std::vector<std::pair<unsigned,std::string>> fields;
    std::ifstream manifest(argv[4]); unsigned p;std::string path;
    while(manifest>>p>>path)fields.emplace_back(p,path);
    if(fields.size()<2)throw std::runtime_error("need CRT fields and an independent field");
    mpz_class modulus=1;
    for(size_t j=0;j<fields.size();++j){
        auto [prime,file]=fields[j];
        std::ifstream input(file);if(!input)throw std::runtime_error("missing modular interior table");
        mpz_class inv;
        mpz_class field=prime;
        if(!mpz_invert(inv.get_mpz_t(),modulus.get_mpz_t(),field.get_mpz_t()))
            throw std::runtime_error("repeated field");
        unsigned multiplier=inv.get_ui();
        size_t k=0; std::array<unsigned,3> residues{};
        for(int i=0;i<int(nodes.size());++i){
            if(!(input>>residues[0]>>residues[1]>>residues[2]))throw std::logic_error("short prime table");
            if(k==free.size()||free[k]!=i)continue;
            unsigned r=residues[ch];
            if(j+1<fields.size()){
                unsigned diff=(r+prime-mpz_fdiv_ui(lift[k].get_mpz_t(),prime))%prime;
                mpz_addmul_ui(lift[k].get_mpz_t(),modulus.get_mpz_t(),uint64_t(diff)*multiplier%prime);
            }else{
                mpz_class x=(lift[k]*denominator)%modulus;
                mpz_class centered=x;if(2*x>modulus)centered-=modulus;
                mpq_class candidate(centered,denominator);candidate.canonicalize();
                if(residue(candidate,prime)!=r){
                    if(rational(candidate,x,modulus)){candidate/=denominator;}
                    if(residue(candidate,prime)!=r &&
                       !(rational(candidate,lift[k],modulus)&&residue(candidate,prime)==r)){
                        std::cerr<<"insufficient CRT at orbit "<<nodes[i]<<std::endl;return 2;
                    }
                }
                value[i]=candidate;
            }
            ++k;
        }
        if(j+1<fields.size())modulus*=prime;
    }
    std::cerr<<"free direct values "<<free.size()<<" CRT bits "<<mpz_sizeinbase(modulus.get_mpz_t(),2)<<std::endl;
    auto position=[&](int id)->size_t{return id<0?nodes.size()+size_t(-1-id):size_t(id);};
    for(size_t i=0;i<nodes.size();++i)if(!fixed[i]&&pivots[i].count){
        auto& eq=pivots[i];int lead=0;mpq_class v=0;
        for(int j=0;j<eq.count;++j){
            if(eq.ids[j]==int(i))lead=eq.coefficients[j];
            else v+=value[position(eq.ids[j])]*int(eq.coefficients[j]);
        }
        if(!lead)throw std::logic_error("zero diagonal");value[i]=-v/lead;
    }
    for(const auto& v:value)
        mpz_lcm(denominator.get_mpz_t(),denominator.get_mpz_t(),v.get_den_mpz_t());
    std::vector<mpz_class> integer;
    integer.reserve(value.size());
    for(const auto& v:value)integer.emplace_back((denominator/v.get_den())*v.get_num());
    unsigned long checked=0;mpz_class residual;
    equations(base,[&](const Equation& eq){
        residual=0;
        for(int j=0;j<eq.count;++j){
            int c=eq.coefficients[j];auto& v=integer[position(eq.ids[j])];
            if(c>0)mpz_addmul_ui(residual.get_mpz_t(),v.get_mpz_t(),c);
            else mpz_submul_ui(residual.get_mpz_t(),v.get_mpz_t(),-c);
        }
        if(residual!=0)throw std::logic_error("exact positive or face identity failed at row "+std::to_string(checked));
        ++checked;
    });
    std::ofstream output(argv[6]);if(!output)throw std::runtime_error("cannot save checked channel");
    for(size_t i=0;i<nodes.size();++i)output<<nodes[i]<<" "<<value[i]<<"\n";
    std::cout<<"channel "<<ch<<" exact rows "<<checked<<" common denominator "<<denominator<<std::endl;
}
