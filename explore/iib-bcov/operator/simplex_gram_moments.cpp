// First apply every soft moment form to one hard row, then pair the result
// with all remaining hard rows. The Python caller proves a common bound
// on the absolute sum of the expanded terms before entering this kernel.
#include <algorithm>
#include <cstdint>
#include <vector>
#include <unordered_map>
using Wide=__int128_t;
using UWide=__uint128_t;

struct Gram {
    int rows, monomials, columns;
    std::vector<int32_t> ends, indices;
    std::vector<int64_t> coefficients;
    std::vector<Wide> gram, partial;
    int bits=0,digits=0;
    std::vector<double> exact_digits,partial_digits;
};
extern "C" void* simplex_gram(int n,int m,int w,int num, const uint64_t* hx,const uint64_t* wx,
        const int32_t* he,const int32_t* hi,const int64_t* hc,
        const int32_t* we,const int32_t* wi,const int64_t* wc,
        const uint64_t* sums,const uint64_t* lifted){
    auto* t=new Gram{n,m,w,{he,he+n+1},{hi,hi+he[n]},{hc,hc+he[n]},
                     std::vector<Wide>(int64_t(m)*m*w),std::vector<Wide>(m*w)};
    std::unordered_map<uint64_t,Wide> moments;
    moments.reserve(num);
    for(int i=0;i<num;i++)moments[sums[i]]=(UWide(lifted[2*i+1])<<64)|lifted[2*i];
    for(int a=0;a<m;a++)for(int b=a;b<m;b++){
        auto* x=&t->gram[(int64_t(a)*m+b)*w];
        auto* y=&t->gram[(int64_t(b)*m+a)*w];
        uint64_t key=hx[a]+hx[b];
        for(int v=0;v<w;v++){
            Wide total=0;
            for(int j=we[v];j<we[v+1];j++)total+=wc[j]*moments.at(key+wx[wi[j]]);
            x[v]=y[v]=total;
        }
    }
    return t;
}
extern "C" void free_simplex_gram(void* p){delete static_cast<Gram*>(p);}
// The caller bounds the absolute sum in each digit by 2^52, including
// every partial sparse dot. Binary64 therefore carries exact integers.
extern "C" void digit_simplex_gram(void* p,int bits,int digits){
    auto& t=*static_cast<Gram*>(p);
    t.bits=bits;t.digits=digits;
    t.exact_digits.resize(t.gram.size()*digits);
    t.partial_digits.resize(t.monomials*t.columns*digits);
    UWide mask=(UWide(1)<<bits)-1;
    for(size_t i=0;i<t.gram.size();i++){
        Wide v=t.gram[i];
        UWide whole=v<0?-v:v;
        for(int d=0;d<digits;d++){
            double x=static_cast<uint64_t>(whole&mask);
            t.exact_digits[i*digits+d]=v<0?-x:x;
            whole>>=bits;
        }
    }
    std::vector<Wide>().swap(t.gram);
    std::vector<Wide>().swap(t.partial);
}

static void digit_row(Gram& t,int i,uint64_t* result){
    int m=t.monomials,w=t.columns,d=t.digits;
    int stride=w*d;
    auto& partial=t.partial_digits;
    std::fill(partial.begin(),partial.end(),0);
    for(int k=t.ends[i];k<t.ends[i+1];k++){
        const double* g=t.exact_digits.data()+int64_t(t.indices[k])*m*stride;
        double weight=t.coefficients[k];
        for(int j=0;j<m*stride;j++)partial[j]+=weight*g[j];
    }
    std::vector<double> total(stride);
    Wide base=Wide(1)<<t.bits;
    for(int j=i;j<t.rows;j++){
        std::fill(total.begin(),total.end(),0);
        for(int k=t.ends[j];k<t.ends[j+1];k++){
            const double* row=partial.data()+t.indices[k]*stride;
            double weight=t.coefficients[k];
            for(int v=0;v<stride;v++)total[v]+=weight*row[v];
        }
        for(int v=0;v<w;v++){
            Wide n=0;
            for(int a=d-1;a>=0;a--)n=n*base+static_cast<int64_t>(total[v*d+a]);
            auto u=static_cast<UWide>(n);
            result[2*(j*w+v)]=uint64_t(u);result[2*(j*w+v)+1]=uint64_t(u>>64);
        }
    }
}

extern "C" void simplex_gram_row(void* p,int i,uint64_t* result){
    auto& t=*static_cast<Gram*>(p);
    if(t.digits){digit_row(t,i,result);return;}
    int m=t.monomials,w=t.columns;
    std::fill(t.partial.begin(),t.partial.end(),0);
    for(int k=t.ends[i];k<t.ends[i+1];k++){
        const auto* g=t.gram.data()+int64_t(t.indices[k])*m*w;
        int64_t weight=t.coefficients[k];
        for(int j=0;j<m*w;j++)t.partial[j]+=weight*g[j];
    }
    std::vector<Wide> total(w);
    for(int j=i;j<t.rows;j++){
        std::fill(total.begin(),total.end(),0);
        for(int k=t.ends[j];k<t.ends[j+1];k++){
            const auto* row=t.partial.data()+t.indices[k]*w;
            int64_t weight=t.coefficients[k];
            for(int v=0;v<w;v++)total[v]+=weight*row[v];
        }
        for(int v=0;v<w;v++){
            auto n=UWide(total[v]);
            result[2*(j*w+v)]=uint64_t(n);
            result[2*(j*w+v)+1]=uint64_t(n>>64);
        }
    }
}

// Core exponents have nine five-bit lanes. Using a contiguous sparse table
// avoids node allocations in the two metric convolutions. No key can equal
// the all-ones sentinel.
template<class Value> class MetricMap {
    static constexpr uint64_t empty=UINT64_MAX;
    std::vector<uint64_t> keys;
    std::vector<Value> values;
    size_t used=0;
    size_t locate(uint64_t key) const {
        uint64_t x=key^(key>>21);
        x*=UINT64_C(0x9e3779b97f4a7c15);
        x^=x>>29;
        size_t i=x&(keys.size()-1);
        while(keys[i]!=empty && keys[i]!=key)i=(i+1)&(keys.size()-1);
        return i;
    }
    void grow(){
        std::vector<uint64_t> oldkeys(keys.size()*2,empty);
        std::vector<Value> oldvalues(oldkeys.size());
        oldkeys.swap(keys);oldvalues.swap(values);
        for(size_t k=0;k<oldkeys.size();k++)if(oldkeys[k]!=empty){
            size_t i=locate(oldkeys[k]);
            keys[i]=oldkeys[k];values[i]=oldvalues[k];
        }
    }
public:
    MetricMap():keys(1024,empty),values(1024){}
    size_t size() const {return used;}
    bool has(uint64_t key) const{return keys[locate(key)]!=empty;}
    Value& add(uint64_t key,Value initial,bool& created){
        size_t i=locate(key);
        created=(keys[i]==empty);
        if(created){
            if((used+1)*5>=keys.size()*4){grow();i=locate(key);}
            keys[i]=key;values[i]=initial;used++;
        }
        return values[i];
    }
    Value& operator[](uint64_t key){bool unused;return add(key,0,unused);}
    const Value& at(uint64_t key) const{return values[locate(key)];}
    struct Iterator{
        const MetricMap* t;
        size_t i;
        void skip(){while(i<t->keys.size()&&t->keys[i]==empty)i++;}
        auto operator*()const{return std::make_pair(t->keys[i],t->values[i]);}
        Iterator& operator++(){i++;skip();return *this;}
        bool operator!=(const Iterator& b)const{return i!=b.i;}
    };
    Iterator begin()const{Iterator it{this,0};it.skip();return it;}
    Iterator end()const{return {this,keys.size()};}
};
using Polynomial=MetricMap<Wide>;

struct MomentBundle{
    int columns;
    MetricMap<size_t> ids;
    std::vector<uint64_t> keys;
    std::vector<Wide> values;
};
extern "C" void free_metric_moments(void* p){delete static_cast<MomentBundle*>(p);}
extern "C" void* empty_metric_moments(int columns){return new MomentBundle{columns,{},{},{}};}

extern "C" void* metric_moment_pairs(void* table,const uint64_t* hc,int first,int last,uint64_t* count){
    auto& t=*static_cast<Gram*>(table);
    auto* out=new MomentBundle{t.columns,{},{},{}};
    auto& ids=out->ids;
    auto& sums=out->keys;
    for(int i=first;i<last;i++)for(int j=i;j<t.rows;j++){
        uint64_t key=hc[i]+hc[j];
        bool created;ids.add(key,sums.size(),created);
        if(created){
            sums.push_back(key);
            if(sums.size()*t.columns>22000000){delete out;return nullptr;}
        }
    }
    int w=t.columns;
    auto& paired=out->values;
    paired.resize(sums.size()*w);
    std::vector<uint64_t> row(2*t.rows*w);
    for(int i=first;i<last;i++){
        simplex_gram_row(table,i,row.data());
        for(int j=i;j<t.rows;j++){
            auto* ptr=paired.data()+ids.at(hc[i]+hc[j])*w;
            const auto* src=row.data()+2*j*w;
            for(int v=0;v<w;v++){
                Wide n=static_cast<Wide>((UWide(src[2*v+1])<<64)|src[2*v]);
                ptr[v]+=(i==j ? n : 2*n);
            }
        }
    }
    *count=sums.size();
    return out;
}

extern "C" void rescale_metric_moments(void* p,int64_t scale){
    for(auto& n:static_cast<MomentBundle*>(p)->values)n*=scale;
}
extern "C" int fit_metric_shifts(void* target,void* source,int ns,const uint64_t* shifts,size_t cap){
    auto& out=*static_cast<MomentBundle*>(target);
    auto& src=*static_cast<MomentBundle*>(source);
    if(out.keys.size()+src.keys.size()*ns<=cap)return 1;
    MetricMap<uint8_t> extra;
    for(auto k:src.keys)for(int j=0;j<ns;j++)if(!out.ids.has(k+shifts[j])){
        extra[k+shifts[j]]=1;
        if(out.keys.size()+extra.size()>cap)return 0;
    }
    return 1;
}
extern "C" size_t add_metric_shifts(void* output,void* source,int ns,
        const uint64_t* shifts,const int64_t* weights){
    auto& out=*static_cast<MomentBundle*>(output);
    auto& src=*static_cast<MomentBundle*>(source);
    int w=src.columns;
    for(size_t i=0;i<src.keys.size();i++)for(int j=0;j<ns;j++){
        uint64_t key=src.keys[i]+shifts[j];
        bool created;
        size_t q=out.ids.add(key,out.keys.size(),created);
        if(created){out.keys.push_back(key);out.values.resize(out.keys.size()*w);}
        auto* dest=out.values.data()+q*w;
        auto* row=src.values.data()+i*w;
        int64_t c=weights[j];
        for(int k=0;k<w;k++)dest[k]+=c*row[k];
    }
    return out.keys.size();
}

extern "C" void* apply_voltage(void* p,int ns,const uint64_t* sc,const int32_t* ends,
        const int32_t* indices,const int64_t* coeff,size_t first,size_t last,uint64_t* count){
    auto& bundle=*static_cast<MomentBundle*>(p);
    const auto& sums=bundle.keys;
    const auto& paired=bundle.values;
    int w=bundle.columns;
    auto *out=new Polynomial;
    for(size_t i=first;i<last;i++){
        const auto* ptr=paired.data()+i*w;
        for(int j=0;j<ns;j++){
            Wide n=0;
            for(int k=ends[j];k<ends[j+1];k++)n+=coeff[k]*ptr[indices[k]];
            if(n){
                (*out)[sums[i]+sc[j]]+=n;
                if(out->size()>6500000){delete out;return nullptr;}
            }
        }
    }
    *count=out->size();return out;
}

extern "C" void* simplex_square(void* table,const uint64_t* hc,int ns,const uint64_t* sc,
        const int32_t* ends,const int32_t* indices,const int64_t* coeff,
        uint64_t* sizes,int first,int last){
    void* pairs=metric_moment_pairs(table,hc,first,last,sizes);
    if(!pairs)return nullptr;
    auto* out=apply_voltage(pairs,ns,sc,ends,indices,coeff,0,sizes[0],sizes+1);
    free_metric_moments(pairs);
    return out;
}

extern "C" size_t simplex_terms(void* p,uint64_t* keys,uint64_t* coeffs){
    auto& poly=*static_cast<Polynomial*>(p);
    size_t n=0;
    for(auto [key,c]:poly)if(c){
        keys[n]=key;
        UWide u=static_cast<UWide>(c);
        coeffs[2*n]=uint64_t(u);coeffs[2*n+1]=uint64_t(u>>64);
        n++;
    }
    return n;
}
extern "C" void free_simplex_terms(void* p){delete static_cast<Polynomial*>(p);}

static void* fold_core(void* in,uint64_t jac,const int32_t* order,
        const int32_t* metrics,int32_t* checks,uint64_t* count,bool prism){
    auto* out=new Polynomial;
    int pows[2][9][32];
    for(int s=0;s<2;s++)for(int e=0;e<9;e++){
        pows[s][e][0]=1;
        for(int k=1;k<32;k++)pows[s][e][k]=pows[s][e][k-1]*metrics[s*9+e]%1009;
    }
    for(auto [monomial,val]:*static_cast<Polynomial*>(in))if(val){
        uint64_t key=monomial+jac;
        int v[9];
        for(int i=0;i<9;i++)v[i]=(key>>(5*i))&31;
        for(int s=0;s<2;s++){
            int n=int(val%1009);
            for(int j=0;j<9;j++)n=n*pows[s][j][v[j]]%1009;
            checks[s]=(checks[s]+n+1009)%1009;
        }
        uint64_t best=UINT64_MAX;
        if(prism){
            for(int g=0;g<12;g++){
                uint64_t candidate=0;
                for(int j=0;j<9;j++)candidate=(candidate<<5)|v[order[9*g+j]];
                best=std::min(best,candidate);
            }
        }else for(int trans=0;trans<2;trans++)for(int a=0;a<3;a++)for(int b=0;b<3;b++)if(a!=b){
            int c=3-a-b,cols[3];
            for(int j=0;j<3;j++){
                int i0=trans?3*j+a:3*a+j,i1=trans?3*j+b:3*b+j,i2=trans?3*j+c:3*c+j;
                cols[j]=(v[order[i0]]<<10)|(v[order[i1]]<<5)|v[order[i2]];
            }
            if(cols[0]>cols[1])std::swap(cols[0],cols[1]);
            if(cols[1]>cols[2])std::swap(cols[1],cols[2]);
            if(cols[0]>cols[1])std::swap(cols[0],cols[1]);
            uint64_t cand=0;
            for(int shift:{10,5,0})for(int n:cols)cand=(cand<<5)|((n>>shift)&31);
            best=std::min(best,cand);
        }
        (*out)[best]+=val;
    }
    *count=out->size();
    return out;
}

extern "C" void* simplex_bipartite_fold(void* in,uint64_t jac,const int32_t* order,
        const int32_t* metrics,int32_t* checks,uint64_t* count){
    return fold_core(in,jac,order,metrics,checks,count,false);
}
extern "C" void* simplex_prism_fold(void* in,uint64_t jac,const int32_t* order,
        const int32_t* metrics,int32_t* checks,uint64_t* count){
    return fold_core(in,jac,order,metrics,checks,count,true);
}
