// Enumerate positive convergent scalar prism powers and completed angular classes.
#include "../quartic_prism_wheel_charts.hpp"
#include <iostream>
#include <stdexcept>

int main(int argc,char** argv){
    if(argc!=2 && argc!=3)throw std::invalid_argument("Expected the positive half-dimension q");
    const bool aliases=argc==3 && std::string(argv[2])=="--finite-aliases";
    if(argc==3 && !aliases)throw std::invalid_argument("Unknown domain output");
    const int q=std::stoi(argv[1]);
    if(q<3||q>13)throw std::invalid_argument("Expected 3 <= q <= 13");
    std::ios::sync_with_stdio(false);
    std::vector<uint64_t> classes;
    uint64_t labeled=0,finite=0;
    // The rectangular permutation group acts as S3 on the three rung powers.
    for(int x=1;x<q;++x)for(int y=x;y<q;++y)for(int z=y;z<q;++z){
        if(x+y+z>=2*q)continue;
        const int rim[3][3]={{0,x,y},{x,0,z},{y,z,0}};
        std::vector<std::array<int,3>> permutations;
        std::array<int,3> perm{0,1,2};
        do{permutations.push_back(perm);}while(std::next_permutation(perm.begin(),perm.end()));
        for(int a=1;a<q;++a)for(int b=1;b<q;++b)for(int c=1;c<q;++c){
            const int bottom=4*q-x-y-z-a-b-c;
            if(a+b+c<=q||bottom<=q||a+b+x>=2*q||a+c+y>=2*q||b+c+z>=2*q)continue;
            for(int d=a;d<q;++d)for(int e=1;e<q;++e){
                const int f=bottom-d-e;
                if(f<1||f>=q||d+e+x>=2*q||d+f+y>=2*q||e+f+z>=2*q)continue;
                if(a+d+x+y<=q||b+e+x+z<=q||c+f+y+z<=q)continue;
                const std::array<int,6> spokes{a,b,c,d,e,f};
                int stabilizer=0;
                uint64_t finite_key=std::numeric_limits<uint64_t>::max();
                bool least=true;
                for(const auto& p:permutations){
                    int i=p[0],j=p[1],k=p[2];
                    int X=rim[i][j],Y=rim[i][k],Z=rim[j][k];
                    for(int t=0;t<=3;t+=3){
                        std::array<int,6> alternative{spokes[t+i],spokes[t+j],spokes[t+k],
                            spokes[3-t+i],spokes[3-t+j],spokes[3-t+k]};
                        if(X==x && Y==y && Z==z){
                            if(alternative<spokes)least=false;
                            if(alternative==spokes)++stabilizer;
                        }
                        const Powers ordered{X,alternative[0],alternative[1],alternative[4],
                            alternative[3],alternative[2],Y,Z,alternative[5]};
                        uint64_t code=0;
                        for(int power:ordered)code=code*q+power;
                        finite_key=std::min(finite_key,code);
                    }
                }
                if(!least)continue;
                ++finite;
                labeled+=12/stabilizer;
                const auto graph=complete(Powers{x,a,b,e,d,c,y,z,f},q);
                if(!convergent(graph,q))throw std::logic_error("Missed a completed collision");
                const auto angular_key=canonical(graph,q);
                if(aliases){
                    uint64_t raw=0,angular=0;
                    for(auto w:decode(finite_key,q))raw=(raw<<4)|w;
                    for(auto w:decode(angular_key,q))angular=(angular<<4)|w;
                    std::cout.write(reinterpret_cast<const char*>(&raw),sizeof(raw));
                    std::cout.write(reinterpret_cast<const char*>(&angular),sizeof(angular));
                }
                if(angular_key!=finite_key)continue;
                uint64_t key=0;
                for(auto w:decode(finite_key,q))key=(key<<4)|w;
                classes.push_back(key);
            }
        }
    }
    std::sort(classes.begin(),classes.end());
    if(std::adjacent_find(classes.begin(),classes.end())!=classes.end())
        throw std::logic_error("Duplicate completed scalar orbit");
    if(!aliases)for(uint64_t key:classes)std::cout<<key<<"\n";
    std::cerr<<"Positive labeled prisms, finite graph orbits, and angular classes: "
        <<labeled<<" "<<finite<<" "<<classes.size()<<"\n";
}
