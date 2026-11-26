// Propagate convergent bipartite identities to vertices admitting a finite star expansion.
#define BCOV_K33_ORBIT_LIBRARY
#include "../quartic_bipartite_orbit_ibp.cpp"
#include <gmpxx.h>
#include <limits>
#include <cstring>
#include <queue>

int smallest_vertex(Key key) {
    const auto p=unpack(key);
    int v=100;
    for(int i=0;i<3;++i)
        v=std::min({v,p[3*i]+p[3*i+1]+p[3*i+2],p[i]+p[3+i]+p[6+i]});
    return v;
}

int verify_positive(const Domain& domain){
    using Columns=std::array<mpq_class,3>;
    std::vector<Columns> values(domain.keys.size()),wheels;
    std::vector<bool> read(domain.keys.size());
    std::unordered_map<Key,int> contacts;
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    Key key;
    Columns columns;
    while(std::cin>>key){
        for(auto& v:columns)if(!(std::cin>>v))throw std::runtime_error("truncated period row");
        auto entry=domain.ids.find(key);
        if(entry==domain.ids.end()){
            auto powers=unpack(key);
            if(std::count(powers.begin(),powers.end(),0)!=1||!finite(powers,domain.q))
                throw std::logic_error("unexpected contracted graph");
            if(!contacts.emplace(key,wheels.size()).second)throw std::logic_error("duplicate contracted graph");
            wheels.push_back(columns);
        }else{
            if(read[entry->second])throw std::logic_error("duplicate positive graph");
            read[entry->second]=true;
            values[entry->second]=columns;
        }
    }
    if(std::find(read.begin(),read.end(),false)!=read.end())throw std::logic_error("missing positive graph");
    long verified=0,positive_only=0,missing_cut=0;
    for(unsigned i=0;i<read.size();++i)for(int rule=0;rule<6;++rule){
        Row row;
        if(!generate(domain,i,rule,row,contacts))continue;
        bool missing=false,internal=true;
        for(int j=0;j<row.size;++j)if(row.i[j]<0){
            internal=false;
            if(static_cast<size_t>(-1-row.i[j])>=wheels.size())missing=true;
        }
        if(missing){++missing_cut;continue;}
        for(unsigned c=0;c<3;++c){
            mpq_class residual;
            for(int j=0;j<row.size;++j){
                int id=row.i[j];
                residual+=row.c[j]*(id>=0?values[id][c]:wheels[-1-id][c]);
            }
            if(residual!=0)throw std::logic_error("nonzero period residual");
        }
        ++verified;
        positive_only+=internal;
    }
    std::cout<<"positive "<<values.size()<<" contracted "<<wheels.size()<<" verified "<<verified
             <<" positive-only "<<positive_only<<" missing contracted input "<<missing_cut<<"\n";
    return 0;
}

Key canonical_scalar_wheel(Key packed,int q){
    int u[4],b[4],a[6][6]{};
    for(int i=3;i>=0;--i){b[i]=packed&15;packed>>=4;}
    for(int i=3;i>=0;--i){u[i]=packed&15;packed>>=4;}
    for(int i=0;i<4;++i){
        int prev=(i+3)%4;
        a[4][i]=a[i][4]=u[i];
        a[prev][i]=a[i][prev]=b[i];
    }
    for(int i=0;i<5;++i){
        int remaining=2*q;
        for(int j=0;j<5;++j)remaining-=a[i][j];
        a[5][i]=a[i][5]=remaining;
    }
    Key best=std::numeric_limits<Key>::max();
    for(int center=0;center<6;++center)for(int inf=0;inf<6;++inf){
        if(center==inf)continue;
        std::array<int,4> rim{};int h=0;
        for(int v=0;v<6;++v)if(v!=center&&v!=inf)rim[h++]=v;
        // Fix one rim vertex, then place its opposite on a missing diagonal.
        for(int other=1;other<4;++other){
            if(a[rim[0]][rim[other]])continue;
            int k=1;while(k==other)++k;
            int l=6-other-k;
            if(a[rim[k]][rim[l]])continue;
            int order[4]={rim[0],rim[k],rim[other],rim[l]};
            for(int shift=0;shift<4;++shift)for(int direction:{1,3}){
                Key code=0;bool valid=true;
                for(int type=0;type<2;++type)for(int j=0;j<4;++j){
                    int x=order[(shift+direction*j)%4];
                    int prev=order[(shift+direction*(j+3))%4];
                    int e=type?a[x][prev]:a[center][x];
                    if(e<1||e>=q)valid=false;
                    code=(code<<4)|unsigned(e&15);
                }
                if(valid)best=std::min(best,code);
            }
        }
    }
    if(best==std::numeric_limits<Key>::max())throw std::logic_error("no convergent scalar wheel frame");
    return best;
}

struct CompletedSystem {
    std::vector<Key> keys,raw_wheels;
    std::vector<int> wheel_ids;
    std::unordered_map<Key,int> ids,scalar_ids;
    size_t positive;
    std::vector<unsigned char> matrix;
    CompletedSystem(const std::string& base,const std::string& aliases){
        std::ifstream core(base+".nodes"),cuts(base+".wheel"),labels(aliases);
        if(!core||!cuts||!labels)throw std::runtime_error("missing complete orbit input");
        for(Key key;core>>key;){
            ids.emplace(key,keys.size());keys.push_back(key);
        }
        positive=keys.size();
        for(Key key;cuts>>key;){
            Key raw,canonical;
            if(!(labels>>raw>>canonical)||raw!=key)throw std::runtime_error("unordered scalar-wheel registration");
            auto [it,fresh]=scalar_ids.emplace(canonical,keys.size());
            if(fresh)keys.push_back(canonical);
            int id=it->second;
            ids.emplace(raw,id);
            raw_wheels.push_back(raw);wheel_ids.push_back(id);
        }
        std::ifstream stream(base+".bin",std::ios::binary|std::ios::ate);
        if(!stream)throw std::runtime_error("missing strict identities");
        matrix.resize(stream.tellg());stream.seekg(0);
        stream.read(reinterpret_cast<char*>(matrix.data()),matrix.size());
    }
    Row equation(size_t offset)const{
        Row row;
        int n=matrix[offset++];
        for(int j=0;j<n;++j){
            int32_t index;
            std::memcpy(&index,matrix.data()+offset,4);
            int id=index<0?wheel_ids[-1-index]:index;
            add(row,id,int8_t(matrix[offset+4]));
            offset+=5;
        }
        return row;
    }
};

void plan_closure(const CompletedSystem& system,const std::string& prefix,bool find_frontier=false){
    struct ShortRow {uint32_t offset,first;unsigned char length;};
    std::vector<ShortRow> rows;
    std::vector<int> columns,incidence_count(system.keys.size());
    std::vector<bool> known(system.keys.size());
    for(Key key;std::cin>>key;)known[system.ids.at(key)]=true;
    if(find_frontier)std::fill(known.begin()+system.positive,known.end(),true);
    for(size_t offset=0;offset<system.matrix.size();){
        const auto row=system.equation(offset);
        uint32_t first=columns.size();
        for(int j=0;j<row.size;++j)if(row.c[j]&&!known[row.i[j]]){
            columns.push_back(row.i[j]);++incidence_count[row.i[j]];
        }
        if(first!=columns.size())rows.push_back({uint32_t(offset),first,uint8_t(columns.size()-first)});
        offset+=1+5*system.matrix[offset];
    }
    std::vector<int> reverse(columns.size()),begin(known.size()+1),cursor;
    for(size_t i=0;i<known.size();++i)begin[i+1]=begin[i]+incidence_count[i];
    cursor=begin;
    std::vector<int> single,degree(known.size());
    std::vector<unsigned char> missing(rows.size());
    for(size_t i=0;i<rows.size();++i){
        const auto row=rows[i];
        missing[i]=row.length;
        if(row.length==1)single.push_back(i);
        for(int k=0;k<row.length;++k){
            int v=columns[row.first+k];
            reverse[cursor[v]++]=i;
            degree[v]+=row.length==2;
        }
    }
    using Choice=std::pair<int,int>;
    std::priority_queue<Choice> choices;
    for(size_t i=system.positive;i<known.size();++i)if(!known[i]&&incidence_count[i])
        choices.emplace(degree[i],i);
    size_t next=0,forced=0,derived=0;
    std::ofstream anchors,pivots;
    if(!find_frontier){
        anchors.open(prefix+".angular");pivots.open(prefix+".pivots",std::ios::binary);
        if(!anchors||!pivots)throw std::runtime_error("cannot open closure outputs");
    }
    auto learn=[&](int target){
        known[target]=true;
        for(int p=begin[target];p<begin[target+1];++p){
            int j=reverse[p];auto remain=--missing[j];
            if(remain==1)single.push_back(j);
            if(remain==2){
                const auto row=rows[j];
                for(int h=0;h<row.length;++h){
                    int v=columns[row.first+h];
                    if(!known[v]&&v>=int(system.positive))choices.emplace(++degree[v],v);
                }
            }
        }
    };
    while(true){
        while(next<single.size()){
            int r=single[next++];
            if(missing[r]!=1)continue;
            const auto row=rows[r];
            int target=-1;
            for(int j=0;j<row.length;++j){
                int v=columns[row.first+j];
                if(!known[v]){target=v;break;}
            }
            if(target<0)throw std::logic_error("missing pivot target");
            if(!find_frontier){
                pivots.write(reinterpret_cast<const char*>(&row.offset),4);
                pivots.write(reinterpret_cast<const char*>(&target),4);
            }
            learn(target);++derived;
        }
        while(!choices.empty()&&known[choices.top().second])choices.pop();
        if(choices.empty())break;
        int target=choices.top().second;choices.pop();
        anchors<<system.keys[target]<<"\n";
        learn(target);++forced;
    }
    if(find_frontier){
        std::ofstream frontier(prefix+".wheel");size_t count=0;
        for(size_t i=0;i<system.positive;++i)if(!known[i]){frontier<<system.keys[i]<<"\n";++count;}
        std::cerr<<"extra positive angular stars "<<count<<std::endl;
        return;
    }
    size_t blind=0;
    std::ofstream unsupported(prefix+".unused");
    for(size_t i=0;i<known.size();++i)if(!known[i]){
        if(i<system.positive||incidence_count[i])throw std::logic_error("incomplete identity closure");
        unsupported<<system.keys[i]<<"\n";++blind;
    }
    std::cerr<<"positive "<<system.positive<<" canonical scalar wheels "<<known.size()-system.positive
             <<" single pivots "<<derived<<" direct angular anchors "<<forced<<" unused scalar wheels "<<blind<<std::endl;
}

void replay_closure(const CompletedSystem& system,const std::string& prefix,int channel,
                    const std::string& angular,const std::string& output){
    if(channel<0||channel>=3)throw std::invalid_argument("period channel");
    std::vector<mpq_class> values(system.keys.size());
    std::vector<bool> known(values.size());
    auto save=[&](int id,const mpq_class& value){
        if(known[id]&&values[id]!=value)throw std::logic_error("unequal exact anchors");
        known[id]=true;values[id]=value;
    };
    auto read=[&](std::istream& input,bool canonical){
        Key key;std::array<mpq_class,3> cols;
        while(input>>key){
            for(auto& val:cols)if(!(input>>val))throw std::runtime_error("truncated known columns");
            save(canonical?system.scalar_ids.at(key):system.ids.at(key),cols[channel]);
        }
    };
    read(std::cin,false);
    std::ifstream direct(angular);if(!direct)throw std::runtime_error("missing scalar anchors");
    read(direct,true);
    std::ifstream pivots(prefix+".pivots",std::ios::binary);
    if(!pivots)throw std::runtime_error("missing triangular recipe");
    uint32_t offset;int target;long count=0;
    while(pivots.read(reinterpret_cast<char*>(&offset),4)){
        if(!pivots.read(reinterpret_cast<char*>(&target),4))throw std::runtime_error("truncated pivot");
        auto row=system.equation(offset);int diagonal=0;mpq_class residual;
        for(int k=0;k<row.size;++k)if(row.c[k]){
            if(row.i[k]==target){diagonal=row.c[k];continue;}
            if(!known[row.i[k]])throw std::logic_error("unordered pivot");
            residual-=row.c[k]*values[row.i[k]];
        }
        if(!diagonal)throw std::logic_error("zero diagonal");
        save(target,residual/diagonal);++count;
    }
    long verified=0,positive_only=0;
    for(size_t pos=0;pos<system.matrix.size();){
        const auto row=system.equation(pos);
        mpq_class residual;
        bool internal=true;
        for(int k=0;k<row.size;++k)if(row.c[k]){
            int id=row.i[k];
            if(!known[id])throw std::logic_error("unknown in strict original row");
            internal &= size_t(id)<system.positive;
            residual+=row.c[k]*values[id];
        }
        if(residual!=0)throw std::logic_error("nonzero residual");
        ++verified;positive_only+=internal;
        pos+=1+5*system.matrix[pos];
    }
    std::ofstream out(output);
    if(!out)throw std::runtime_error("cannot write completed table");
    for(size_t i=0;i<system.positive;++i){
        if(!known[i])throw std::logic_error("missing positive period");
        out<<system.keys[i]<<" "<<values[i]<<"\n";
    }
    for(size_t i=0;i<system.raw_wheels.size();++i)if(known[system.wheel_ids[i]])
        out<<system.raw_wheels[i]<<" "<<values[system.wheel_ids[i]]<<"\n";
    std::cerr<<"channel "<<channel<<" solved "<<count<<" verified "<<verified<<" positive-only "<<positive_only<<std::endl;
}

int main(int argc,char **argv) {
    if(argc>=5&&(std::string(argv[1])=="plan"||std::string(argv[1])=="frontier"||std::string(argv[1])=="resolve")){
        std::ios::sync_with_stdio(false);std::cin.tie(nullptr);
        CompletedSystem system(argv[2],argv[3]);
        if(std::string(argv[1])=="frontier"&&argc==5)plan_closure(system,argv[4],true);
        else if(std::string(argv[1])=="plan"&&argc==5)plan_closure(system,argv[4]);
        else if(argc==8)replay_closure(system,argv[4],std::stoi(argv[5]),argv[6],argv[7]);
        else throw std::invalid_argument("plan base aliases recipe, or resolve base aliases recipe channel angular output");
        return 0;
    }
    if(argc!=2&&argc!=3)throw std::invalid_argument("half dimension and optional positive-system file prefix");
    const int q=std::stoi(argv[1]);
    if(argc==3&&std::string(argv[2])=="--canonical-wheel"){
        std::ios::sync_with_stdio(false);std::cin.tie(nullptr);
        for(Key key;std::cin>>key;)std::cout<<canonical_scalar_wheel(key,q)<<"\n";
        return 0;
    }
    if(argc==3&&std::string(argv[2])=="--contact-map"){
        std::ios::sync_with_stdio(false);std::cin.tie(nullptr);
        for(Key key;std::cin>>key;){
            auto powers=unpack(key);int zero=std::find(powers.begin(),powers.end(),0)-powers.begin();
            if(zero>=9||std::count(powers.begin(),powers.end(),0)!=1)throw std::logic_error("expected one zero line");
            int i=zero/3,j=zero%3,k=(i+1)%3,l=(i+2)%3,m=(j+1)%3,n=(j+2)%3;
            Key wheel=0;
            for(int e:{3*k+m,3*l+m,3*l+n,3*k+n,3*k+j,3*i+m,3*l+j,3*i+n})wheel=(wheel<<4)|powers[e];
            std::cout<<key<<" "<<canonical_scalar_wheel(wheel,q)<<"\n";
        }
        return 0;
    }
    Domain domain(q);
    if(argc==3&&std::string(argv[2])=="--verify-positive")return verify_positive(domain);
    if(argc==3&&std::string(argv[2])=="--keys"){
        for(const auto key:domain.keys)std::cout<<key<<"\n";
        return 0;
    }
    const int n=domain.keys.size();
    std::vector<int> smallest,pivot(n,-1),order,queue;
    for(auto k:domain.keys)smallest.push_back(smallest_vertex(k));
    std::vector<std::vector<int>> incident(n);
    struct Rule {int source,rule,count;};
    std::vector<Rule> rules;
    std::unordered_map<Key,int> contacts;
    std::vector<int> unknown_index(n,-1);
    std::unordered_map<Key,int> leaves;
    std::ofstream system;
    if(argc==3) {
        std::ofstream nodes(std::string(argv[2])+".nodes");
        system.open(std::string(argv[2])+".bin",std::ios::binary);
        int index=0;
        for(int i=0;i<n;++i)if(smallest[i]>q) {
            unknown_index[i]=index++;
            nodes<<domain.keys[i]<<"\n";
        }
    }
    std::vector<Key> contact_keys;
    for(int i=0;i<n;++i)if(smallest[i]>=q)
        for(int rule=0;rule<6;++rule) {
            Row row;
            if(!generate(domain,i,rule,row,contacts,&contact_keys))continue;
            const int id=rules.size();
            int count=0;
            for(int j=0;j<row.size;++j)if(row.i[j]>=0&&smallest[row.i[j]]>q) {
                incident[row.i[j]].push_back(id);
                ++count;
            }
            if(count)rules.push_back({i,rule,count});
            if(count==1)queue.push_back(id);
            if(count&&system.is_open()) {
                uint8_t sz=row.size;
                system.write(reinterpret_cast<char*>(&sz),1);
                for(int j=0;j<row.size;++j) {
                    const int col=row.i[j];int idx=col<0?-1:unknown_index[col];
                    if(idx<0) {
                        Key key=col<0?contact_keys[-1-col]:domain.keys[col];
                        idx=~leaves.emplace(key,int(leaves.size())).first->second;
                    }
                    int8_t coeff=row.c[j];
                    system.write(reinterpret_cast<char*>(&idx),4);
                    system.write(reinterpret_cast<char*>(&coeff),1);
                }
            }
        }
    if(system.is_open()) {
        system.close();
        std::vector<Key> codes(leaves.size());
        for(auto [key,i]:leaves)codes[i]=key;
        std::ofstream table(std::string(argv[2])+".wheel");
        for(Key key:codes)table<<key<<"\n";
        std::cerr<<"finite positive-star and contracted boundary orbits "<<codes.size()<<std::endl;
    }
    for(size_t i=0;i<queue.size();++i) {
        const int id=queue[i];
        const auto &r=rules[id];
        if(r.count!=1)continue;
        Row row;
        if(!generate(domain,r.source,r.rule,row,contacts))throw std::logic_error("lost integral relation");
        int target=-1;
        for(int j=0;j<row.size;++j)
            if(row.i[j]>=0&&smallest[row.i[j]]>q&&pivot[row.i[j]]<0)target=row.i[j];
        if(target<0)throw std::logic_error("lost unit pivot");
        pivot[target]=id;
        order.push_back(target);
        for(int ix:incident[target])if(--rules[ix].count==1)queue.push_back(ix);
    }
    const int unknown=std::count_if(smallest.begin(),smallest.end(),[q](int v){return v>q;});
    std::cerr<<"dimension "<<2*q<<" nodes "<<n<<" unmatched "<<unknown
             <<" single-pivot reductions "<<order.size()<<" rows "<<rules.size()<<std::endl;
    std::unordered_map<int,mpq_class> weights;
    Key key;mpq_class w;
    while(std::cin>>key>>w)weights[domain.ids.at(key)]+=w;
    for(auto it=order.rbegin();it!=order.rend();++it) {
        auto pos=weights.find(*it);
        if(pos==weights.end()||!sgn(pos->second))continue;
        const mpq_class weight=pos->second;
        weights.erase(pos);
        const auto &r=rules[pivot[*it]];
        Row row;
        if(!generate(domain,r.source,r.rule,row,contacts))throw std::logic_error("lost pivot equation");
        int diagonal=0;
        for(int j=0;j<row.size;++j)if(row.i[j]==*it)diagonal=row.c[j];
        if(!diagonal)throw std::logic_error("zero pivot coefficient");
        for(int j=0;j<row.size;++j)if(row.i[j]!=*it)
            weights[row.i[j]]-=weight*row.c[j]/diagonal;
    }
    std::vector<Key> wheels(contacts.size());
    for(auto [k,i]:contacts)wheels[i]=k;
    for(const auto &[i,v]:weights)if(sgn(v)) {
        if(i<0)std::cout<<"W "<<wheels[-1-i];
        else std::cout<<(smallest[i]<=q?"S ":"U ")<<domain.keys[i];
        std::cout<<" "<<v<<"\n";
    }
}
