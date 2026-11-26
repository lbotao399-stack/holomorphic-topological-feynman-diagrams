// The scalar CE map on seven balanced positive vector jets, modulo total divergence.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <stdexcept>
#include <vector>
using namespace std;
#ifndef BALANCED_PRIME
#define BALANCED_PRIME 251
#endif
constexpr int p=BALANCED_PRIME;
template<class T> void shuffle_reproducibly(vector<T>& values, mt19937& rng) {
    for(size_t n=values.size();n>1;n--)std::swap(values[n-1],values[rng()%n]);
}
using Vec=array<int,5>;
using Perm=array<int,7>;
using Support=array<uint8_t,7>;
using Marks=array<uint8_t,3>;
int reduce(int x){return (x%p+p)%p;}
int inverse(int a){int b=1;for(int i=p-2;i;i/=2,a=a*a%p)if(i&1)b=b*a%p;return b;}
int dot(const Vec&a,const Vec&b){return inner_product(a.begin(),a.end(),b.begin(),0)%p;}
int determinant(vector<Vec> m,int skip=-1){
    int d=1,n=m.size(),j=0;
    for(int i=0;i<n;i++,j++){
        if(j==skip)j++;
        int k=i;while(k<n&&!m[k][j])k++;
        if(k==n)return 0;
        if(k!=i){swap(m[i],m[k]);d=-d;}
        d=d*m[i][j]%p;
        int inv=inverse(m[i][j]);
        for(k=i+1;k<n;k++){
            int c=m[k][j]*inv%p;
            for(int t=j+1;t<5;t++)m[k][t]=reduce(m[k][t]-c*m[i][t]);
        }
    }
    return reduce(d);
}
struct Waves{vector<Vec> momenta,polarizations;};
Waves sample(int n,mt19937&rng){
    Waves test{vector<Vec>(n),vector<Vec>(n)};
    for(int j=0;j<5;j++){
        int total=0;
        for(int i=0;i<n-1;i++){test.momenta[i][j]=rng()%p;total+=test.momenta[i][j];}
        test.momenta.back()[j]=reduce(-total);
    }
    for(int i=0;i<n;i++){
        vector<int> others;
        for(int j=0;j<n;j++)if(i!=j)others.push_back(j);
        shuffle_reproducibly(others,rng);
        vector<Vec> normals={test.momenta[i],test.momenta[others[0]],test.momenta[others[1]],test.momenta[others[2]]};
        for(int j=0;j<5;j++)test.polarizations[i][j]=reduce((j%2?-1:1)*determinant(normals,j));
        if(dot(test.polarizations[i],test.momenta[i]))throw runtime_error("transversality");
    }
    return test;
}
Waves fuse(const Waves& w,int i,int j){
    Waves out{{Vec{}},{Vec{}}};
    int a=dot(w.polarizations[i],w.momenta[j]),b=dot(w.polarizations[j],w.momenta[i]);
    for(int k=0;k<5;k++){
        out.momenta[0][k]=(w.momenta[i][k]+w.momenta[j][k])%p;
        out.polarizations[0][k]=reduce(a*w.polarizations[j][k]-b*w.polarizations[i][k]);
    }
    for(int k=0;k<int(w.momenta.size());k++)if(k!=i&&k!=j){
        out.momenta.push_back(w.momenta[k]);out.polarizations.push_back(w.polarizations[k]);
    }
    return out;
}
struct Group{Perm f;vector<Support> contractions;vector<Marks> marks;};
vector<Perm> permutations;
vector<pair<int,int>> endpoints;
vector<array<uint8_t,21>> edge_images;
vector<Group> groups;
void prepare(){
    Perm perm={0,1,2,3,4,5,6};
    for(int i=0;i<7;i++)for(int j=i+1;j<7;j++)endpoints.push_back({i,j});
    do{
        permutations.push_back(perm);
        array<uint8_t,21> edges;
        for(int e=0;e<21;e++)edges[e]=perm[endpoints[e].first]*7+perm[endpoints[e].second];
        edge_images.push_back(edges);
    }while(next_permutation(perm.begin(),perm.end()));
    // Enumerate one representative of each conjugacy class of maps without fixed points.
    set<int> seen;mt19937 rng(1009);
    for(int index=0;index<279936;index++){
        if(seen.count(index))continue;
        Perm f{},indegree{};int code=index;
        for(int i=6;i>=0;i--){int v=code%6;code/=6;f[i]=v+(v>=i);indegree[f[i]]++;}
        Group g{f,{},{}};
        for(auto perm:permutations){
            Perm image{};Support es;
            for(int i=0;i<7;i++){image[perm[i]]=perm[f[i]];es[i]=perm[i]*7+perm[f[i]];}
            int h=0;for(int i=0;i<7;i++)h=6*h+image[i]-(image[i]>i);
            seen.insert(h);g.contractions.push_back(es);
        }
        vector<Marks> valid;
        for(int a=0;a<21;a++)for(int b=a;b<21;b++)for(int c=b;c<21;c++){
            Marks marks={(uint8_t)a,(uint8_t)b,(uint8_t)c};Perm omitted{};
            for(int e:marks){omitted[endpoints[e].first]++;omitted[endpoints[e].second]++;}
            bool ok=true;for(int i=0;i<7;i++)if(3+indegree[i]-omitted[i]<2)ok=false;
            if(ok)valid.push_back(marks);
        }
        shuffle_reproducibly(valid,rng);valid.resize(min(size_t(48),valid.size()));
        g.marks=valid;groups.push_back(std::move(g));
    }
    shuffle_reproducibly(groups,rng);
}
struct Column{int group,mark;};
struct Rank{
    map<int,vector<int>> pivots;int determinant=1;
    void insert(vector<int> row){
        for(const auto&[k,v]:pivots){
            int scale=row[k];if(!scale)continue;
            for(int j=k;j<int(v.size());j++)row[j]=reduce(row[j]-scale*v[j]);
        }
        auto it=find_if(row.begin(),row.end(),[](int a){return a;});if(it==row.end())return;
        int k=it-row.begin();determinant=determinant*row[k]%p;int u=inverse(row[k]);
        for(int j=k;j<int(row.size());j++)row[j]=row[j]*u%p;
        pivots.emplace(k,std::move(row));
    }
};
void evaluate(const Waves& w,const vector<vector<pair<int,int>>>&selected,vector<int>&result,int sign){
    array<int,49> pairing{},minor{};
    for(int i=0;i<7;i++)for(int j=0;j<7;j++)pairing[7*i+j]=dot(w.polarizations[i],w.momenta[j]);
    for(auto [i,j]:endpoints){
        vector<Vec> m;
        for(int k=0;k<7;k++)if(k!=i&&k!=j)m.push_back(w.momenta[k]);
        // Sign of [sorted complement,i,j], for i<j.
        int delta=reduce((2*((i+j)%2)-1)*determinant(m));
        minor[i*7+j]=delta;minor[j*7+i]=reduce(-delta);
    }
    for(int x=0;x<int(groups.size());x++){
        if(selected[x].empty())continue;
        const auto&g=groups[x];vector<int64_t> totals(selected[x].size());
        for(int n=0;n<5040;n++){
            auto es=g.contractions[n];
            int v=1;
            for(int e:es){if(!pairing[e]){v=0;break;}v=v*pairing[e]%p;}
            if(!v)continue;
            array<int,21> delta;
            for(int e=0;e<21;e++)delta[e]=minor[edge_images[n][e]];
            for(int j=0;j<int(selected[x].size());j++){
                auto marks=g.marks[selected[x][j].first];
                totals[j]+=int64_t(v)*delta[marks[0]]*delta[marks[1]]*delta[marks[2]];
            }
        }
        for(int j=0;j<int(selected[x].size());j++){
            int& a=result[selected[x][j].second];a=reduce(a+sign*(totals[j]%p));
        }
    }
}
int main(int argc,char**argv){
    prepare();mt19937 rng(12348);
    int limit=argc>1?stoi(argv[1]):365;
    vector<vector<pair<int,int>>> selected(groups.size());
    vector<Column> available;
    for(int g=0;g<int(groups.size());g++)for(int m=0;m<int(groups[g].marks.size());m++){
        selected[g].push_back({m,available.size()});available.push_back({g,m});
    }
    Rank initial;
    for(int i=0;i<limit;i++){
        vector<int> row(available.size());evaluate(sample(7,rng),selected,row,1);initial.insert(std::move(row));
    }
    cerr<<"raw probes "<<limit<<" rank "<<initial.pivots.size()<<" groups "<<groups.size()<<endl;
    selected.assign(groups.size(),{});int count=0;
    for(const auto&[i,row]:initial.pivots){
        auto col=available[i];selected[col.group].push_back({col.mark,count++});
    }
    Rank differential;
    for(int test=0;test<limit;test++){
        auto waves=sample(8,rng);vector<int> row(count);
        for(int i=0;i<8;i++)for(int j=i+1;j<8;j++)
            if(dot(waves.polarizations[i],waves.momenta[j])||dot(waves.polarizations[j],waves.momenta[i]))
                evaluate(fuse(waves,i,j),selected,row,(i+j)%2?-1:1);
        differential.insert(std::move(row));
        if((test+1)%50==0)cerr<<"differential probes "<<test+1<<" rank "<<differential.pivots.size()<<endl;
        if(differential.pivots.size()==340)break;
    }
    if(limit>=365 && (count!=349 || differential.pivots.size()!=340))
        throw runtime_error("missing balanced rank certificate");
    cout<<"raw "<<count<<" differential "<<differential.pivots.size()<<" pivot product "<<differential.determinant<<endl;
}
