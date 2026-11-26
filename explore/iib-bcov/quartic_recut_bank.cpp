// Out-of-core single-endpoint convolution for the K3,3 planar and nine-cut projections.
#include <gmpxx.h>
#include <algorithm>
#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <stdexcept>
#include <vector>

using Exponents=std::array<int,8>;
using Integer=mpz_class;
using Rational=mpq_class;
namespace fs=std::filesystem;

struct Coefficient {
    std::array<uint64_t,4> words{};
    bool zero()const{return !(words[0]|words[1]|words[2]|words[3]);}
    void add(const Coefficient& source,int factor){
        __int128 carry=0;
        for(int j=0;j<4;++j){
            auto t=static_cast<__int128>(source.words[j])*factor+words[j]+carry;
            words[j]=static_cast<uint64_t>(t);carry=t>>64;
        }
    }
    Integer value()const{
        Integer v;
        mpz_import(v.get_mpz_t(),4,-1,8,0,0,words.data());
        if(words[3]>>63)v-=Integer(1)<<256;
        return v;
    }
};
struct State {uint32_t packed;int flat;uint32_t through;};
struct Active {State state;Coefficient value;};
struct Entry {uint32_t packed;int offset,coefficient;};
struct Target {int flat,order,pair_depth;int64_t k,j;int count;};
struct Tile {int first,second;std::vector<State> states;std::vector<Target> targets;};
using Ops=std::map<std::pair<int,int>,std::vector<Entry>>;
std::pair<int,int> first_two(const Exponents& e){return {e[0],e[1]};}
uint32_t pack(const Exponents& e){
    uint32_t v=0;
    for(int k=2;k<8;++k)v|=uint32_t(e[k])<<(5*(k-2));
    return v;
}
template<class T> void write_vector(const fs::path& path,const std::vector<T>& data){
    std::ofstream out(path,std::ios::binary);
    out.write(reinterpret_cast<const char*>(data.data()),data.size()*sizeof(T));
    if(!out)throw std::runtime_error("could not write coefficient tile");
}
template<class T> std::vector<T> read_vector(const fs::path& path,size_t count){
    std::vector<T> data(count);
    std::ifstream in(path,std::ios::binary);
    in.read(reinterpret_cast<char*>(data.data()),count*sizeof(T));
    if(!in)throw std::runtime_error("incomplete coefficient tile");
    return data;
}

int main(int argc,char**argv){
    if(argc!=4)throw std::invalid_argument("operators  packed-weights  temporary-directory");
    std::ifstream weights(argv[2],std::ios::binary),operators(argv[1]);
    const fs::path work=argv[3];
    if(!fs::is_directory(work) || !fs::is_empty(work))throw std::invalid_argument("empty temporary directory required");
    std::array<std::vector<std::pair<Exponents,int>>,2> raw_ops;
    for(auto& terms:raw_ops){
        int count;operators>>count;
        if(count<=0)throw std::runtime_error("missing Gram polynomial");
        terms.resize(count);
        for(auto& [e,c]:terms){for(auto& k:e)operators>>k;operators>>c;}
    }
    Exponents largest{},stride{};
    auto decode=[](uint32_t key){
        Exponents e;
        for(int k=7;k>=0;--k){e[k]=key&15;key>>=4;}return e;
    };
    auto for_targets=[&](auto emit){
        weights.clear();weights.seekg(0);
        for(uint32_t key;weights.read(reinterpret_cast<char*>(&key),4);){
            int64_t k,j;uint8_t count;
            weights.read(reinterpret_cast<char*>(&k),8);weights.read(reinterpret_cast<char*>(&j),8);
            weights.read(reinterpret_cast<char*>(&count),1);
            if(!weights)throw std::runtime_error("incomplete packed weight");
            emit(decode(key),k,j,count);
        }
    };
    for_targets([&](const Exponents& e,int64_t,int64_t,int){for(int k=0;k<8;++k)largest[k]=std::max(largest[k],e[k]);});
    int volume=1,last=0;
    for(int k=7;k>=0;--k){stride[k]=volume;volume*=largest[k]+1;}
    const int plane=stride[1];
    auto encode=[&](const Exponents& e){return std::inner_product(e.begin(),e.end(),stride.begin(),0);};
    auto tile_id=[&](int i,int j){return i*(largest[1]+1)+j;};
    std::vector<Tile> tiles(volume/plane);
    std::vector<int8_t> needed(volume,-1);
    for_targets([&](const Exponents& e,int64_t k,int64_t j,int count){
        int r=std::accumulate(e.begin(),e.end(),-20),at=encode(e);
        if(r<0||r>19)throw std::logic_error("physical target outside Gram degree");
        needed[at]=std::max<int>(needed[at],r); last=std::max(last,r);
        tiles[at/plane].targets.push_back({at%plane,r,std::min(e[5],e[6]),k,j,count});
    });
    for(int k=7;k>=0;--k){
        int width=stride[k],block=width*(largest[k]+1);
        for(int b=0;b<volume;b+=block)for(int j=b+block-1;j>=b+width;--j)
            needed[j-width]=std::max(needed[j-width],needed[j]);
    }
    auto meta=[&](int tile){return work/("states-"+std::to_string(tile));};
    auto page=[&](int r,int tile){return work/(std::to_string(r)+"-"+std::to_string(tile));};
    std::vector<size_t> sizes(tiles.size());
    Exponents lower{};
    for(size_t i=0;i<tiles.size();++i){
        Tile& tile=tiles[i];tile.first=i/(largest[1]+1);tile.second=i%(largest[1]+1);
        lower={};
        for(int j=0;j<plane;++j){
            int degree=needed[i*plane+j];
            if(degree>=0)tile.states.push_back({pack(lower),j,uint32_t(degree)});
            for(int k=7;k>=2;--k){if(++lower[k]<=largest[k])break;lower[k]=0;}
        }
        sizes[i]=tile.states.size();write_vector(meta(i),tile.states);
        std::vector<State>().swap(tile.states);
    }
    std::vector<int8_t>().swap(needed);
    std::cerr<<"tiles "<<tiles.size()<<" compressed cells "<<std::accumulate(sizes.begin(),sizes.end(),size_t(0))
             <<" largest tile "<<*std::max_element(sizes.begin(),sizes.end())<<std::endl;
    std::array<Ops,2> operations;
    std::array<int,2> one_norm{};
    for(int h=0;h<2;++h)for(const auto&[e,c]:raw_ops[h]){
        operations[h][first_two(e)].push_back({pack(e),encode(e)%plane,c});
        one_norm[h]+=abs(c);
    }
    Exponents cap{},mask{};
    for(int k=2;k<8;++k){cap[k]=15-largest[k];mask[k]=16;}
    const uint32_t bias=pack(cap),overflow=pack(mask);
    std::vector<std::vector<size_t>> active(last+1,std::vector<size_t>(tiles.size()));
    std::array<Rational,4> total{}; Integer factorial=1,n0=0,n1=1;
    for(int r=0;r<=last;++r){
        std::array<Integer,4> accum{};
        if(r>0)factorial*=r;
        const std::array<int,2> scale{21-2*r,4*(r-1)*(21-r)};
        if(r>0){
            Integer bound=abs(scale[0])*one_norm[0]*n1+abs(scale[1])*one_norm[1]*n0;
            if(bound>=(Integer(1)<<255))throw std::overflow_error("centered integer convolution");
            n0=n1;n1=bound;
        }
        for(size_t t=0;t<tiles.size();++t){
            if(!sizes[t])continue;
            auto dest=read_vector<State>(meta(t),sizes[t]);
            std::vector<Coefficient> result(sizes[t]);
            std::vector<int> index(plane,-1);
            for(size_t j=0;j<dest.size();++j)if(int(dest[j].through)>=r)index[dest[j].flat]=int(j);
            if(r==0 && t==0)result[index[0]].words[0]=1;
            for(int h=0;h<2;++h){
                if(r<h+1 || !scale[h])continue;
                for(auto& [pair,terms]:operations[h]){
                    int x=tiles[t].first-pair.first,y=tiles[t].second-pair.second;
                    if(x<0||y<0)continue;
                    int source=tile_id(x,y);if(!sizes[source])continue;
                    if(!active[r-h-1][source])continue;
                    auto data=read_vector<Active>(page(r-h-1,source),active[r-h-1][source]);
                    for(size_t j=0;j<data.size();++j){
                        if(int(data[j].state.through)<r)continue;
                        auto packed=data[j].state.packed+bias;
                        auto flat=data[j].state.flat;
                        for(const Entry& term:terms){
                            if((packed+term.packed)&overflow)continue;
                            int k=index[flat+term.offset];
                            if(k>=0)result[k].add(data[j].value,scale[h]*term.coefficient);
                        }
                    }
                }
            }
            for(const Target& row:tiles[t].targets)if(row.order==r){
                Integer binomial=1,value=0;
                for(int j=0;j<=row.pair_depth;++j){
                    value+=binomial*result[index[row.flat-j*(stride[5]+stride[6])]].value();
                    binomial*=-2*(19-2*r-2*j);
                    if(!mpz_divisible_ui_p(binomial.get_mpz_t(),j+1))throw std::logic_error("nonintegral conic coefficient");
                    binomial/=j+1;
                }
                if(!value)continue;
                for(int f=0;f<2;++f){
                    const int mult=f ? row.count>>4 : row.count&15;
                    if(!mult)continue;
                    Integer weighted=mult*value;
                    accum[2*f]+=weighted*static_cast<long>(row.k);
                    accum[2*f+1]+=weighted*static_cast<long>(row.j);
                }
            }
            if(r<last){
                std::vector<Active> output;
                for(size_t j=0;j<result.size();++j)if(!result[j].zero() && int(dest[j].through)>r)
                    output.push_back({dest[j],result[j]});
                active[r][t]=output.size();
                if(!output.empty())write_vector(page(r,t),output);
            }
        }
        for(int j=0;j<4;++j){Rational v(accum[j],(Integer(1)<<r)*factorial);v.canonicalize();total[j]+=v;}
        if(r>=2)for(size_t i=0;i<tiles.size();++i)fs::remove(page(r-2,i));
        std::cerr<<"degree "<<r<<" of "<<last<<" live coefficients "<<std::accumulate(active[r].begin(),active[r].end(),size_t(0))<<std::endl;
    }
    Rational normalization(1,Integer(8*1925)*(Integer(1)<<38));
    for(int j=0;j<10;++j)normalization*=Rational(8,(2*j+1)*(2*j+2)*(2*j+3));
    for(auto v:total){v*=normalization;v.canonicalize();std::cout<<v<<"\n";}
}
