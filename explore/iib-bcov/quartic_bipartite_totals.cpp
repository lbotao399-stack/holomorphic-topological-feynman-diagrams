// Integrate the two physical numerator banks and the three independently known U^9 lifts.
#define BCOV_K33_ORBIT_LIBRARY
#include "quartic_bipartite_orbit_ibp.cpp"
#include <gmpxx.h>
#include <sstream>
#include <zlib.h>

std::pair<int,std::array<mpq_class,2>> weighted_bank(
        const char* path,const std::unordered_map<Key,mpq_class>& values){
    gzFile gz=gzopen(path,"rb");
    if(!gz)throw std::runtime_error("missing physical graph weights");
    std::array<mpq_class,2> sums{};
    int count=0;
    char buffer[1024];
    while(gzgets(gz,buffer,sizeof(buffer))){
        std::istringstream line(buffer);
        Matrix ordinary,p;
        for(auto& x:ordinary)if(!(line>>x))throw std::runtime_error("missing edge power");
        const int order[9]={0,1,2,3,5,6,4,7,8};
        for(int j=0;j<9;++j)p[j]=ordinary[order[j]];
        const auto& exact=values.at(canonical(p));
        for(auto& s:sums){
            mpq_class weight;
            if(!(line>>weight))throw std::runtime_error("short physical weights");
            s+=weight*exact;
        }
        std::string extra;
        if(line>>extra)throw std::runtime_error("unexpected weight column");
        ++count;
    }
    if(gzclose(gz)!=Z_OK||!count)throw std::logic_error("invalid physical bank");
    return {count,sums};
}

int main(int argc,char** argv){
    if(argc<5)throw std::invalid_argument("checked channel table, U^9 orbit counts, physical weights.gz, channel id, optional additional banks");
    int ch=std::stoi(argv[4]);
    if(ch<0||ch>2)throw std::invalid_argument("invalid period column");
    std::unordered_map<Key,mpq_class> values;
    std::ifstream certified(argv[1]);
    Key key;
    mpq_class val;
    while(certified>>key>>val){
        if(!values.emplace(key,val).second)throw std::runtime_error("duplicate certified value");
    }
    if(values.size()!=695485)throw std::logic_error("incomplete q=12 period table");
    unsigned long fact[12]{1};
    for(int j=1;j<12;++j)fact[j]=fact[j-1]*j;
    std::array<mpq_class,3> anchors{};
    std::ifstream calibration(argv[2]);size_t n=0;
    while(calibration>>key){
        auto p=unpack(key);
        auto ref=values.at(key);
        for(int a:p)ref*=fact[a-1];
        for(int i=0;i<3;++i){
            unsigned long count;
            if(!(calibration>>count))throw std::runtime_error("short U^9 table");
            anchors[i]+=count*ref;
        }
        ++n;
    }
    const mpq_class reference[3][3]={{0,mpq_class(-7,6),mpq_class(5,3)},
                                     {0,mpq_class(5,3),mpq_class(-5,3)},
                                     {0,-2,mpq_class(5,2)}};
    if(n!=values.size())throw std::logic_error("incomplete absolute anchors");
    for(int i=0;i<3;++i)if(anchors[i]/fact[11]!=reference[i][ch])
        throw std::logic_error("exact dimension-lift anchor failure");
    auto [count,sums]=weighted_bank(argv[3],values);
    if(count!=465564)throw std::logic_error("incomplete physical bank");
    if(ch==2 && (sums[0]!=mpq_class("249372365272156031/806626220400000")
             || sums[1]!=mpq_class("2209968587039007449/1393850108851200000")))
        throw std::logic_error("independent nine-cut weight-five test failed");
    std::cout<<"q12 exact three integral anchors and "<<count<<" physical orbits: "
             <<ch<<" "<<sums[0]<<" "<<sums[1]<<std::endl;
    for(int i=5;i<argc;++i){
        auto [rows,extra]=weighted_bank(argv[i],values);
        std::cout<<"bank\t"<<argv[i]<<"\t"<<rows<<"\t"<<extra[0]<<"\t"<<extra[1]<<std::endl;
    }
}
