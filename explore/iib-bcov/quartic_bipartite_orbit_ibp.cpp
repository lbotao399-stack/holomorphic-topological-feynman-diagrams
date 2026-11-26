// Convergence-preserving determinant symmetries on positive K3,3 orbit periods.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <vector>

using Matrix=std::array<int,9>;
using Key=uint64_t;
constexpr int field=1000000007;
const int permutations[6][3]={{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};

Key pack(const Matrix& p){Key out=0;for(int v:p)out=(out<<4)|v;return out;}
Matrix unpack(Key key){Matrix p;for(int i=8;i>=0;--i){p[i]=key&15;key>>=4;}return p;}
bool finite(const Matrix& p,int q){
    for(int k=0;k<9;++k)if(p[k]<0||p[k]>=q)return false;
    std::array<int,3> rows{},cols{};
    for(int i=0;i<3;++i)for(int j=0;j<3;++j){rows[i]+=p[3*i+j];cols[j]+=p[3*i+j];}
    for(int i=0;i<3;++i)if(rows[i]>=2*q || cols[i]>=2*q)return false;
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)
        if(4*q-rows[i]-cols[j]+p[3*i+j]<=q)return false;
    return true;
}

Key canonical(const Matrix& p){
    Key best=uint64_t(1)<<40;
    for(int transpose=0;transpose<2;++transpose){
        int c[3][3];
        for(int i=0;i<3;++i)for(int j=0;j<3;++j)c[i][j]=transpose?p[3*j+i]:p[3*i+j];
        for(const auto& ord:permutations){
            std::array<int,3> cols{};
            for(int j=0;j<3;++j)cols[j]=(c[ord[0]][j]<<8)|(c[ord[1]][j]<<4)|c[ord[2]][j];
            if(cols[0]>cols[1])std::swap(cols[0],cols[1]);
            if(cols[1]>cols[2])std::swap(cols[1],cols[2]);
            if(cols[0]>cols[1])std::swap(cols[0],cols[1]);
            Key key=0;
            for(int k:{8,4,0})for(int v:cols)key=(key<<4)|((v>>k)&15);
            best=std::min(best,key);
        }
    }
    return best;
}

struct Domain{
    int q;
    std::vector<Key> keys;
    std::unordered_map<Key,int> ids;
    explicit Domain(int level):q(level){
        Matrix p{};
        auto enumerate=[&](auto&& self,int k,int remaining)->void{
            if(k==8){
                if(remaining<1 || remaining>=q)return;
                p[k]=remaining;
                if(finite(p,q) && canonical(p)==pack(p))keys.push_back(pack(p));
                return;
            }
            int lo=std::max(1,remaining-(8-k)*(q-1));
            int hi=std::min(q-1,remaining-(8-k));
            if(k==1 || k==2)lo=std::max(lo,p[k-1]);
            // Sorted rows and columns are necessary for canonical row-major order.
            if(k==3 || k==6)lo=std::max(lo,p[k-3]);
            if(k==4 || k==7){
                if(p[k-1]==p[k-4])lo=std::max(lo,p[k-3]);
                if(p[0]==p[1] && (k==4 || p[3]==p[4]))lo=std::max(lo,p[k-1]);
            }
            if(k==5 && p[1]==p[2])lo=std::max(lo,p[4]);
            for(int v=lo;v<=hi;++v){p[k]=v;self(self,k+1,remaining-v);}
        };
        enumerate(enumerate,0,4*q);
        for(unsigned i=0;i<keys.size();++i)ids[keys[i]]=int(i);
    }
};

struct Row{std::array<int,7> i{},c{};int size=0;};
void add(Row& row,int index,int coeff){
    for(int j=0;j<row.size;++j)if(row.i[j]==index){row.c[j]+=coeff;return;}
    row.i[row.size]=index;row.c[row.size++]=coeff;
}
bool generate(const Domain& d,int source,int rule,Row& row,std::unordered_map<Key,int>& contacts,
              std::vector<Key>* contact_keys=nullptr){
    Matrix p=unpack(d.keys[source]);
    int last=rule%3,left=(last+1)%3,right=(last+2)%3;
    auto edge=[&](int x,int y){return rule<3?3*x+y:3*y+x;};
    int diag=0;
    std::array<std::pair<Key,int>,7> all{};
    int n=0;
    auto term=[&](Key key,int v){
        for(int i=0;i<n;++i)if(all[i].first==key){all[i].second+=v;return;}
        all[n++]={key,v};
    };
    for(int col=0;col<3;++col){
        int x=edge(left,col),y=edge(right,col),z=edge(last,col);
        diag+=p[x]-p[y];
        for(auto [at,sign]:{std::pair<int,int>{x,1},{y,-1}}){
            --p[at];++p[z];term(canonical(p),sign*(p[z]-1));++p[at];--p[z];
        }
    }
    term(d.keys[source],diag);
    for(int k=0;k<n;++k)if(all[k].second){
        auto [key,coeff]=all[k];
        auto found=d.ids.find(key);
        if(found!=d.ids.end())add(row,found->second,coeff);
        else{
            auto m=unpack(key);
            if(!finite(m,d.q))return false;
            if(*std::min_element(m.begin(),m.end())!=0)throw std::logic_error("incomplete positive orbit enumeration");
            auto [it,fresh]=contacts.emplace(key,int(contacts.size()));
            if(fresh&&contact_keys)contact_keys->push_back(key);
            add(row,-1-it->second,coeff);
        }
    }
    return row.size>0;
}

#ifndef BCOV_K33_ORBIT_LIBRARY
int main(int argc,char**argv){
    if(argc<2 || argc>3)throw std::invalid_argument("q  optional-output-basename");
    Domain d(std::stoi(argv[1]));
    std::cerr<<"positive orbits "<<d.keys.size()<<std::endl;
    std::ofstream out;
    if(argc==3){out.open(std::string(argv[2])+".bin",std::ios::binary);}
    std::unordered_map<Key,int> contacts;
    std::vector<uint8_t> mask(d.keys.size());
    unsigned long rows=0;
    for(int i=0;i<int(d.keys.size());++i)for(int rule=0;rule<6;++rule){
        Row row;
        if(!generate(d,i,rule,row,contacts))continue;
        int leading=-1;
        for(int j=0;j<row.size;++j)leading=std::max(leading,row.i[j]);
        if(leading<0)continue;
        mask[leading]=1;++rows;
        if(out.is_open()){
            uint8_t n=row.size;out.write(reinterpret_cast<char*>(&n),1);
            for(int j=0;j<row.size;++j){
                out.write(reinterpret_cast<char*>(&row.i[j]),4);
                int8_t c=row.c[j];out.write(reinterpret_cast<char*>(&c),1);
            }
        }
    }
    std::cerr<<"wheel contacts "<<contacts.size()<<" valid rows "<<rows<<" directly missing leads "
             <<std::count(mask.begin(),mask.end(),0)<<std::endl;
    if(out.is_open()){
        std::ofstream nodes(std::string(argv[2])+".nodes");
        for(Key key:d.keys)nodes<<key<<"\n";
        std::vector<Key> wheel(contacts.size());
        for(auto [k,i]:contacts)wheel[i]=k;
        std::ofstream faces(std::string(argv[2])+".wheel");
        for(Key key:wheel)faces<<key<<"\n";
    }
}
#endif
