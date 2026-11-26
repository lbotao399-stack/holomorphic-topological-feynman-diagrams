// Transpose of a selected triangular positive-orbit system; original wheel rows retain their order.
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

using Term=std::pair<int,int>;
using Row=std::vector<Term>;
static int prime;
int mod(int64_t a){a%=prime; return int(a<0?a+prime:a);}
int inv(int v){
    if(!v)throw std::runtime_error("zero divisor");
    int result=1;
    for(int e=prime-2;e;e/=2){if(e%2)result=mod(int64_t(result)*v);v=mod(int64_t(v)*v);}
    return result;
}
int coefficient(const Row& row,int index){
    auto it=std::lower_bound(row.begin(),row.end(),Term{index,0});
    return it!=row.end()&&it->first==index?it->second:0;
}
void scale(Row& row,int s){for(auto& t:row)t.second=mod(int64_t(s)*t.second);}
bool add(Row& a,const Row& b,int factor,int cap){
    Row r;r.reserve(std::min<int>(a.size()+b.size(),cap+2));
    size_t i=0,j=0;
    while(i<a.size()||j<b.size()){
        if(j==b.size()||(i<a.size()&&a[i].first<b[j].first))r.push_back(a[i++]);
        else if(i==a.size()||a[i].first>b[j].first){r.emplace_back(b[j].first,mod(int64_t(factor)*b[j].second));++j;}
        else{int v=mod(a[i].second+int64_t(factor)*b[j].second);if(v)r.emplace_back(a[i].first,v);++i;++j;}
        if(int(r.size())>cap)return false;
    }
    a.swap(r);return true;
}
struct Bank{
    const uint8_t* data;size_t size;std::vector<uint32_t> offset;int n=0,wheels=0;
    explicit Bank(const std::string& base){
        for(auto [ext,count]:{std::pair{".nodes",&n},std::pair{".wheel",&wheels}}){
            std::ifstream s(base+ext);uint64_t key;while(s>>key)++*count;
        }
        int fd=open((base+".bin").c_str(),O_RDONLY);if(fd<0||!n)throw std::runtime_error("missing input bank");
        struct stat st{};fstat(fd,&st);size=st.st_size;
        data=static_cast<const uint8_t*>(mmap(nullptr,size,PROT_READ,MAP_PRIVATE,fd,0));close(fd);
        if(data==MAP_FAILED||size>UINT32_MAX)throw std::runtime_error("bank mapping failed");
        for(size_t at=0;at<size;){offset.push_back(at);at+=1+5*data[at];}
    }
    ~Bank(){munmap((void*)data,size);}
    template<class F>void scanrow(int id,F f)const{
        auto p=data+offset.at(id);int n=*p++;
        for(int k=0;k<n;++k){int32_t i;memcpy(&i,p,4);int c=int8_t(p[4]);p+=5;f(i,c);}
    }
    Row interior(int id)const{
        Row r;scanrow(id,[&](int i,int c){if(i>=0&&c)r.emplace_back(i,mod(c));});
        std::sort(r.begin(),r.end());return r;
    }
};
// Negative indices encode original rows by bitwise complement. Other indices name earlier normalized pivots.
struct Expression{int first;std::vector<int> targets;};
// A stored row is normalization * (first + sum_j edges[j].coefficient * edges[j].row).
struct DAG{
    struct Node{int first,offset,count,normalization;};
    std::vector<Node> nodes;
    std::vector<Term> edges;
    int save(int first,const Row& uses,int norm){
        int id=nodes.size();nodes.push_back({first,int(edges.size()),int(uses.size()),norm});
        edges.insert(edges.end(),uses.begin(),uses.end());return id;
    }
};
struct Pivot{Row row;int token=-1;};
struct Factor{
    const Bank& b;std::vector<Pivot> piv;DAG dag;std::vector<int> terminal,free;
    explicit Factor(const Bank& bank):b(bank),piv(b.n){}
    int add_pivot(Row& row,int first,const Row& uses){
        int factor=inv(row.back().second);scale(row,factor);
        return dag.save(first,uses,factor);
    }
    void discover(){
        std::vector<int> shortest(b.n,-1),width(b.n,100);
        for(int r=0;r<int(b.offset.size());++r){
            int lead=-1,sz=0;b.scanrow(r,[&](int i,int c){if(i>=0&&c){lead=std::max(lead,i);++sz;}});
            if(lead>=0&&sz<width[lead]){shortest[lead]=r;width[lead]=sz;}
        }
        int rank=0;
        for(int i=0;i<b.n;++i)if(shortest[i]>=0){
            auto& p=piv[i];p.row=b.interior(shortest[i]);p.token=add_pivot(p.row,~shortest[i],{});++rank;
        }
        for(int cap:{8,16,32,64,128,256,512}){
            for(int r=0;r<int(b.offset.size());++r){
                Row row=b.interior(r);int first=~r;Row uses;
                while(!row.empty()){
                    auto& p=piv[row.back().first];
                    if(p.row.empty()){
                        p.token=add_pivot(row,first,uses);p.row=std::move(row);++rank;break;
                    }
                    if(row.size()<p.row.size()){
                        int old=p.token;p.token=add_pivot(row,first,uses);p.row.swap(row);
                        uses.clear();first=old;
                    }
                    int mult=prime-row.back().second;uses.emplace_back(p.token,mult);
                    if(!add(row,p.row,mult,cap))break;
                }
            }
            std::cerr<<"cap "<<cap<<" free "<<b.n-rank<<" history "<<dag.nodes.size()<<" edges "<<dag.edges.size()<<std::endl;
            if(b.n-rank<=128)break;
        }
        finish({});
    }
    // The triangular pivots leave a small space of provisional free coordinates.
    // Select original equations on that space, then complete with coordinate rows.
    void finish(const std::vector<int>& prescribed){
        std::vector<int> holes;for(int i=0;i<b.n;++i)if(piv[i].row.empty())holes.push_back(i);
        int h=holes.size();if(h>128)throw std::runtime_error("too many unresolved variables");
        std::vector<int> shadow(size_t(b.n)*h);int next=0;
        for(int i=0;i<b.n;++i){
            auto dest=shadow.data()+size_t(i)*h;const auto& row=piv[i].row;
            if(row.empty()){dest[next++]=1;continue;}
            for(int k=0;k<h;++k){
                __uint128_t sum=0;
                for(size_t j=0;j+1<row.size();++j)sum+=uint64_t(row[j].second)*shadow[size_t(row[j].first)*h+k];
                dest[k]=mod(-int64_t(sum%prime));
            }
        }
        std::vector<std::vector<int>> small(h);terminal.clear();free.clear();
        auto accept=[&](int id,const Row& raw){
            std::vector<int> v(h);
            for(int k=0;k<h;++k){__uint128_t sum=0;for(auto [i,c]:raw)sum+=uint64_t(c)*shadow[size_t(i)*h+k];v[k]=int(sum%prime);}
            for(int k=0;k<h;++k)if(v[k]){
                if(small[k].empty()){int scale=inv(v[k]);for(int j=k;j<h;++j)v[j]=mod(int64_t(v[j])*scale);
                    small[k]=std::move(v);terminal.push_back(id);return true;}
                int scale=prime-v[k];for(int j=k;j<h;++j)v[j]=mod(v[j]+int64_t(scale)*small[k][j]);
            }
            return false;
        };
        if(prescribed.empty()){
            for(int r=0;r<int(b.offset.size())&&int(terminal.size())<h;++r)accept(r,b.interior(r));
        }else for(int r:prescribed)if(r>=0&&!accept(r,b.interior(r)))throw std::runtime_error("exceptional terminal row");
        for(int k=0;k<h;++k)if(accept(~holes[k],{{holes[k],1}}))free.push_back(holes[k]);
        std::cerr<<"residual interior dimension "<<free.size()<<" of "<<h<<std::endl;
        terminal_shadow.assign(h,std::vector<int>(h));
        for(int j=0;j<h;++j){int id=terminal[j];auto row=id>=0?b.interior(id):Row{{~id,1}};
            for(int k=0;k<h;++k){__uint128_t v=0;for(auto [i,c]:row)v+=uint64_t(c)*shadow[size_t(i)*h+k];terminal_shadow[k][j]=v%prime;}
        }
    }
    std::vector<std::vector<int>> terminal_shadow;
    void transfer(std::vector<int>& w,std::vector<int>* multipliers=nullptr)const{
        for(int i=b.n-1;i>=0;--i)if(w[i]&&!piv[i].row.empty()){
            int c=w[i];if(multipliers)(*multipliers)[piv[i].token]=mod((*multipliers)[piv[i].token]+int64_t(c));
            for(auto [j,k]:piv[i].row)w[j]=mod(w[j]-int64_t(c)*k);
        }
    }
    std::vector<int> solve_terminal(std::vector<int> rhs,bool transpose=false)const{
        auto m=terminal_shadow;int h=m.size();
        if(transpose)for(int i=0;i<h;++i)for(int j=0;j<i;++j)std::swap(m[i][j],m[j][i]);
        for(int i=0;i<h;++i){
            int j=i;while(j<h&&!m[j][i])++j;if(j==h)throw std::runtime_error("singular terminal basis");
            std::swap(m[i],m[j]);std::swap(rhs[i],rhs[j]);int s=inv(m[i][i]);
            for(int k=i;k<h;++k)m[i][k]=mod(int64_t(m[i][k])*s);rhs[i]=mod(int64_t(rhs[i])*s);
            for(int j=0;j<h;++j)if(j!=i&&m[j][i]){int c=m[j][i];
                for(int k=i;k<h;++k)m[j][k]=mod(m[j][k]-int64_t(c)*m[i][k]);rhs[j]=mod(rhs[j]-int64_t(c)*rhs[i]);}
        }
        return rhs;
    }
    // Split the input covector into terminal rows and triangular pivot rows.
    // Expanding the recorded operations backwards gives original row multipliers.
    // Their positive-column contraction must reproduce the input with no remainder.
    std::vector<int> adjoint(const std::vector<int>& input,std::vector<int>& residual)const{
        auto w=input;transfer(w);std::vector<int> small;
        for(int i=0;i<b.n;++i)if(piv[i].row.empty())small.push_back(w[i]);
        auto y=solve_terminal(small);std::vector<int> raw(b.offset.size()),history(dag.nodes.size());
        w=input;residual.assign(b.n,0);
        for(size_t r=0;r<terminal.size();++r){
            int i=terminal[r];if(i>=0){raw[i]=y[r];b.scanrow(i,[&](int j,int c){if(j>=0)w[j]=mod(w[j]-int64_t(y[r])*c);});}
            else{w[~i]=mod(w[~i]-int64_t(y[r]));residual[~i]=y[r];}
        }
        transfer(w,&history);
        if(std::any_of(w.begin(),w.end(),[](int v){return v!=0;}))throw std::runtime_error("transposed free reduction failed");
        auto accumulate=[&](int i,int v){auto& cell=i<0?raw[~i]:history[i];cell=mod(int64_t(cell)+v);};
        for(int i=int(history.size())-1;i>=0;--i)if(history[i]){
            auto n=dag.nodes[i];int base=mod(int64_t(history[i])*n.normalization);
            accumulate(n.first,base);for(int j=0;j<n.count;++j){auto [token,c]=dag.edges[n.offset+j];accumulate(token,mod(int64_t(c)*base));}
        }
        w=residual;std::vector<int> wheel(b.wheels);
        for(int i=0;i<int(raw.size());++i)if(raw[i])b.scanrow(i,[&](int j,int c){
            if(j>=0)w[j]=mod(w[j]+int64_t(raw[i])*c);
            else wheel[~j]=mod(wheel[~j]-int64_t(raw[i])*c);
        });
        if(w!=input)throw std::runtime_error("original-row adjoint certificate failed");
        return wheel;
    }
};
template<class S>void io(S& stream,int& value){
    if constexpr(std::is_same_v<S,std::ifstream>){stream.read((char*)&value,4);if(!stream)throw std::runtime_error("short recipe");}
    else stream.write((char*)&value,4);
}
void save(const Factor& f,const std::string& path){
    std::ofstream s(path,std::ios::binary);int n=f.dag.nodes.size();io(s,n);
    for(auto node:f.dag.nodes){
        io(s,node.first);io(s,node.count);
        for(int j=0;j<node.count;++j){int id=f.dag.edges[node.offset+j].first;io(s,id);}
    }
    for(auto& p:f.piv){int i=p.token;io(s,i);}
    n=f.terminal.size();io(s,n);for(auto i:f.terminal)io(s,i);
}
// Only the sequence of eliminated row IDs crosses prime fields. Recalculate
// every modular multiplier from the current field, including the leading inverse.
void load(Factor& f,const std::string& path){
    std::ifstream s(path,std::ios::binary);int n;io(s,n);
    std::vector<Expression> recipe(n);std::vector<int> references(n);
    for(auto& node:recipe){
        io(s,node.first);int k;io(s,k);node.targets.resize(k);if(node.first>=0)++references[node.first];
        for(int& j:node.targets){io(s,j);++references[j];}
    }
    std::vector<int> final(f.b.n);for(int& i:final){io(s,i);if(i>=0)++references[i];}
    int h;io(s,h);std::vector<int> terminal(h);for(int& i:terminal)io(s,i);
    std::vector<Row> live(n);
    auto drop=[&](int id){if(id>=0&&!--references[id])Row().swap(live[id]);};
    for(int i=0;i<n;++i){
        auto& e=recipe[i];Row row=e.first<0?f.b.interior(~e.first):live[e.first];drop(e.first);
        Row uses;uses.reserve(e.targets.size());
        for(int j:e.targets){
            const auto& p=live[j];int idx=p.back().first;int c=coefficient(row,idx);
            if(c){add(row,p,prime-c,1000000);uses.emplace_back(j,prime-c);}
            drop(j);
        }
        if(row.empty())throw std::runtime_error("exceptional recipe pivot");int id=f.add_pivot(row,e.first,uses);
        if(id!=i)throw std::runtime_error("recipe index mismatch");if(references[i])live[i].swap(row);
        std::vector<int>().swap(e.targets);
    }
    for(int i=0;i<f.b.n;++i)if(final[i]>=0){
        f.piv[i].token=final[i];f.piv[i].row=std::move(live[final[i]]);
        if(f.piv[i].row.back().first!=i)throw std::runtime_error("exceptional leading position");
    }
    f.finish(terminal);
}
int main(int argc,char** argv){
    if(argc<5)throw std::runtime_error("discover|replay bank prime structural-recipe [physical/anchor input]");
    prime=std::stoi(argv[3]);Bank b(argv[2]);Factor f(b);
    if(std::string(argv[1])=="discover"){f.discover();save(f,argv[4]);}
    else if(std::string(argv[1])=="replay"||std::string(argv[1])=="solve"||std::string(argv[1])=="boundary")load(f,argv[4]);
    else throw std::runtime_error("unknown operation");
    if(argc<6)return 0;
    if(std::string(argv[1])=="boundary"){
        if(argc!=7)throw std::runtime_error("need an exact boundary table and an output table");
        extern void solve_boundary(const Factor&,const char*,const char*,const char*);
        solve_boundary(f,argv[2],argv[5],argv[6]);
        return 0;
    }
    if(argc==6){extern void check_terminal(const Bank&,const Factor&,const char*);check_terminal(b,f,argv[5]);return 0;}
    extern void measure(const Bank&,const Factor&,char**,int);
    measure(b,f,argv,argc);
    return 0;
}
#include <unordered_map>
#include <zlib.h>
int residue(std::string s){
    int sign=1,r=0;size_t p=0;if(s[0]=='-'){sign=-1;++p;}
    for(;p<s.size()&&s[p]!='/';++p)r=mod(10LL*r+(s[p]-'0'));
    r=mod(sign*r);return p==s.size()?r:mod(int64_t(r)*inv(residue(s.substr(p+1))));
}
using Key=uint64_t;
Key canonical_graph(const std::array<int,9>& ex){
    static const int order[9]={0,1,2,3,5,6,4,7,8};std::array<int,9> p{};
    for(int i=0;i<9;++i)p[i]=ex[order[i]];
    Key best=UINT64_MAX;
    for(int tr=0;tr<2;++tr){
        int a[3]={0,1,2};
        do{
            int col[3]={};
            for(int j=0;j<3;++j)for(int i=0;i<3;++i)col[j]=(col[j]<<4)|p[tr?3*j+a[i]:3*a[i]+j];
            std::sort(col,col+3);Key k=0;for(int shift:{8,4,0})for(int x:col)k=(k<<4)|((x>>shift)&15);
            best=std::min(best,k);
        }while(std::next_permutation(a,a+3));
    }
    return best;
}
#include <sstream>
using Values=std::array<int,4>;
Values plus(Values a,const Values& b,int factor){
    for(int j=0;j<4;++j)a[j]=mod(a[j]+int64_t(factor)*b[j]);return a;
}
std::vector<Values> forward(const Factor& f,const std::vector<std::array<int,3>>& wheel){
    const auto& b=f.b;
    if(f.free.size()>1)throw std::runtime_error("forward solver requires at most one free mode");
    auto rhs=[&](int id){
        Values v{};b.scanrow(id,[&](int j,int c){if(j<0)for(int k=0;k<3;++k)v[k]=mod(v[k]-int64_t(c)*wheel[~j][k]);});return v;
    };
    std::vector<Values> history(f.dag.nodes.size()),out(b.n);
    for(size_t j=0;j<history.size();++j){
        auto e=f.dag.nodes[j];auto v=e.first<0?rhs(~e.first):history[e.first];
        for(int k=0;k<e.count;++k){auto [i,c]=f.dag.edges[e.offset+k];v=plus(v,history[i],c);}
        for(int k=0;k<4;++k)history[j][k]=mod(int64_t(e.normalization)*v[k]);
    }
    std::vector<int> holes;
    auto substitute=[&](){
        for(int i=0;i<b.n;++i){
            const auto& p=f.piv[i];if(p.row.empty())continue;
            auto v=history[p.token];for(size_t k=0;k+1<p.row.size();++k){auto [j,c]=p.row[k];v=plus(v,out[j],prime-c);}
            out[i]=v;
        }
    };
    substitute();for(int i=0;i<b.n;++i)if(f.piv[i].row.empty())holes.push_back(i);
    std::vector<Values> terminal;
    for(int r:f.terminal){
        Values v{};
        if(r<0)v[3]=1; // one free coordinate, normalized to one in the homogeneous column
        else{v=rhs(r);b.scanrow(r,[&](int i,int c){if(i>=0)v=plus(v,out[i],-c);});}
        terminal.push_back(v);
    }
    for(int k=0;k<4;++k){
        std::vector<int> a;for(const auto& v:terminal)a.push_back(v[k]);
        auto h=f.solve_terminal(a,true);for(size_t j=0;j<holes.size();++j)out[holes[j]][k]=h[j];
    }
    substitute();
    for(int r=0;r<int(b.offset.size());++r){
        auto v=rhs(r);b.scanrow(r,[&](int i,int c){if(i>=0)v=plus(v,out[i],-c);});
        if(v!=Values{})throw std::runtime_error("original-row solution residual at "+std::to_string(r));
    }
    return out;
}
void check_terminal(const Bank& b,const Factor& f,const char* values){
    if(b.n!=3||f.free.size()!=1)throw std::runtime_error("expected the q=3 control");
    std::ifstream input(values);std::vector<std::array<int,3>> wheel(b.wheels);
    for(auto& v:wheel)if(!(input>>v[0]>>v[1]>>v[2]))throw std::runtime_error("short control table");
    auto s=forward(f,wheel);
    const std::string ref[3][3]={{"0","-7/3","10/3"},{"0","10/3","-10/3"},{"0","-4","5"}};
    std::array<int,3> free{};
    for(int k=0;k<3;++k)free[k]=residue(ref[f.free.front()][k]);
    for(int i=0;i<3;++i){
        std::vector<int> target(3),rem;target[i]=1;auto c=f.adjoint(target,rem);
        for(int k=0;k<3;++k){
            int d=mod(int64_t(rem[f.free.front()])*free[k]);
            for(int j=0;j<b.wheels;++j)d=mod(d+int64_t(c[j])*wheel[j][k]);
            if(d!=mod(s[i][k]+int64_t(s[i][3])*free[k])||d!=residue(ref[i][k]))throw std::runtime_error("terminal direct/dual discrepancy");
        }
    }
    std::cout<<"All three q=3 terminal periods pass in each of three channels, in both directions."<<std::endl;
}
void solve_boundary(const Factor& f,const char* base,const char* path,const char* output){
    if(!f.free.empty())throw std::runtime_error("boundary data do not determine every interior value");
    std::ifstream keys(std::string(base)+".wheel"),known(path);
    std::vector<std::array<int,3>> values(f.b.wheels);
    for(auto& v:values){
        Key expected,key;
        if(!(keys>>expected)||!(known>>key)||key!=expected)throw std::runtime_error("boundary key mismatch");
        std::string coeff;
        for(int& x:v){if(!(known>>coeff))throw std::runtime_error("short boundary table");x=residue(coeff);}
    }
    std::ofstream dest(output);
    for(const auto& a:forward(f,values))dest<<a[0]<<" "<<a[1]<<" "<<a[2]<<"\n";
    if(!dest)throw std::runtime_error("cannot write interior field table");
}
void measure(const Bank& b,const Factor& f,char**argv,int argc){
    if(argc!=9)throw std::runtime_error("additional arguments: physical.weights.gz anchor.tsv wheel_values.tsv output-prefix");
    std::unordered_map<Key,int> lookup;std::vector<Key> keys;
    std::ifstream nodes(std::string(argv[2])+".nodes");Key key;while(nodes>>key){lookup[key]=keys.size();keys.push_back(key);}
    std::vector<std::vector<int>> weights(6,std::vector<int>(b.n));
    gzFile gz=gzopen(argv[5],"rb");if(!gz)throw std::runtime_error("missing weights");
    char line[4096];while(gzgets(gz,line,sizeof(line))){
        std::istringstream s(line);std::array<int,9> ex;for(int& i:ex)s>>i;
        std::string k,j;s>>k>>j;if(!s)throw std::runtime_error("bad physical row");int id=lookup.at(canonical_graph(ex));
        weights[0][id]=mod(weights[0][id]+int64_t(residue(k)));weights[1][id]=mod(weights[1][id]+int64_t(residue(j)));
    }
    gzclose(gz);
    int fac[12]={1};for(int i=1;i<12;++i)fac[i]=mod(int64_t(i)*fac[i-1]);
    // The determinant identities use C(a) = Gamma(12) / prod Gamma(a_e) times
    // the projective period. The physical .weights.gz already multiplies C(a).
    auto normalized=[&](Key key){int r=inv(fac[11]);for(int i=0;i<9;++i){r=mod(int64_t(r)*fac[(key&15)-1]);key>>=4;}return r;};
    std::ifstream a(argv[6]);std::string nums[3];int count=0;
    while(a>>key>>nums[0]>>nums[1]>>nums[2]){
        int i=lookup.at(key),factor=normalized(key);for(int c=0;c<3;++c)weights[c+2][i]=mod(int64_t(factor)*residue(nums[c]));++count;
    }
    if(count!=b.n)throw std::runtime_error("incomplete anchors");
    auto calibration=canonical_graph({9,6,5,6,4,5,4,5,4});weights[5][lookup.at(calibration)]=normalized(calibration);
    std::vector<std::array<int,3>> wheel(b.wheels);std::ifstream vals(argv[7]);
    for(auto& x:wheel)if(!(vals>>x[0]>>x[1]>>x[2]))throw std::runtime_error("short wheel table");
    std::array<std::array<int,3>,6> measured{};std::vector<std::vector<int>> free_weights(6);
    auto solution=forward(f,wheel);
    for(int j=0;j<6;++j){
        Values direct{};
        for(int i=0;i<b.n;++i)if(weights[j][i])direct=plus(direct,solution[i],weights[j][i]);
        if(std::string(argv[1])=="solve"){
            for(int k=0;k<3;++k)measured[j][k]=direct[k];
            if(!f.free.empty())free_weights[j].push_back(direct[3]);
            continue;
        }
        std::vector<int> residual;auto c=f.adjoint(weights[j],residual);
        for(int i:f.free)free_weights[j].push_back(residual[i]);
        std::array<__uint128_t,3> sums{};
        int support=0;for(int i=0;i<b.wheels;++i)if(c[i]){++support;for(int k=0;k<3;++k)sums[k]+=uint64_t(c[i])*wheel[i][k];}
        std::cout<<"functional "<<j<<" wheel_support "<<support<<" wheel";for(int k=0;k<3;++k){
            measured[j][k]=sums[k]%prime;std::cout<<" "<<measured[j][k];}
        std::cout<<" residual";for(int v:free_weights[j])std::cout<<" "<<v;std::cout<<std::endl;
        for(int k=0;k<3;++k)if(direct[k]!=measured[j][k])throw std::runtime_error("dual and primal projected values disagree");
        if(!f.free.empty()&&direct[3]!=free_weights[j][0])throw std::runtime_error("dual and primal free coefficients disagree");
        if(j<2){std::ofstream out(std::string(argv[8])+"."+(j?"J":"K")+".bin",std::ios::binary);
            out.write((const char*)c.data(),sizeof(int)*c.size());}
    }
    if(f.free.size()>1)throw std::runtime_error("need multiple independent calibration rows");
    const std::string references[4][3]={{"0","-7/6","5/3"},{"0","5/3","-5/3"},{"0","-2","5/2"},
        {"-186964375676833/9579551657195520000","-53055501437/1016349037704000","1693/21339318"}};
    std::array<int,3> unknown{};bool anchored=f.free.empty();
    for(int j=2;j<6;++j)if(!anchored&&free_weights[j][0]){
        for(int k=0;k<3;++k)unknown[k]=mod(int64_t(mod(residue(references[j-2][k])-measured[j][k]))*inv(free_weights[j][0]));
        anchored=true;
    }
    if(!anchored)throw std::runtime_error("all supplied anchors miss the free period");
    for(int j=0;j<6;++j)for(int k=0;k<3;++k){
        if(!f.free.empty())measured[j][k]=mod(measured[j][k]+int64_t(free_weights[j][0])*unknown[k]);
        if(j>=2&&measured[j][k]!=residue(references[j-2][k]))throw std::runtime_error("calibration discrepancy at "+std::to_string(j)+","+std::to_string(k));
    }
    if(measured[0][2]!=residue("249372365272156031/806626220400000")||
       measured[1][2]!=residue("2209968587039007449/1393850108851200000"))
        throw std::runtime_error("known physical zeta(5) coefficient disagrees");
    for(int j=0;j<2;++j){std::cout<<(j?"J":"K");for(auto v:measured[j])std::cout<<" "<<v;std::cout<<std::endl;}
    std::ofstream null(std::string(argv[8])+".homogeneous");
    for(size_t i=0;i<solution.size();++i)if(solution[i][3])null<<i<<" "<<solution[i][3]<<"\n";
    std::ofstream interior(std::string(argv[8])+".nodes.tsv");
    for(auto v:solution){for(int k=0;k<3;++k)interior<<mod(v[k]+int64_t(v[3])*unknown[k])<<(k==2?'\n':' ');}
    std::ofstream result(std::string(argv[8])+".periods");for(int j=0;j<2;++j)result<<measured[j][0]<<" "<<measured[j][1]<<" "<<measured[j][2]<<"\n";
}
