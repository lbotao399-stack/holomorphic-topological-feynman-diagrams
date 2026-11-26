// Exact Fermat odd projection before integrating the two remaining hard points.
#define PARTIAL_COTREES_LIBRARY
#include "partial_cotrees_mod.cpp"

namespace {
int tangent_voltage(const Bundle& b){
    const Field& f=b.f;
    const auto& g=b.g;
    const auto& q=b.q;
    const Group& p=b.groups[0];
    const Group& r=b.groups[1];
    const Env env=b.environment(0,1);
    const int n=p.masks.size(), ne=g.edges.size(), nv=g.green.size();
    // Both hard singleton matrices are symmetric. Fold their two off-diagonal
    // copies before applying the common three-normal exterior contraction.
    Vec ids,swapped;
    for(int i=0;i<n;++i)for(int j=i;j<n;++j){
        ids.push_back(i*n+j);swapped.push_back(j*n+i);
    }
    const int count=ids.size();
    std::vector<std::vector<std::pair<int,int>>> coefficient(count);
    for(int i=0;i<count;++i)for(int j=0;j<count;++j){
        int a=ids[i],c=swapped[i],u=ids[j],v=swapped[j];
        int h=env.coefficients[a][u];
        if(a!=c)h=f.add(h,env.coefficients[c][u]);
        if(u!=v)h=f.add(h,env.coefficients[a][v]);
        if(a!=c&&u!=v)h=f.add(h,env.coefficients[c][v]);
        if(h)coefficient[i].push_back({j,h});
    }
    Mat mass=zero(nv,count),greenMass=zero(nv,count);
    Vec total(count);
    for(int j=0;j<count;++j)total[j]=r.mean[ids[j]];
    for(int e=0;e<ne;++e)for(int h=0;h<5;++h){
        auto [left,right]=g.edges[e];
        int a=f.mul(q.weights[h],f.sub(1,q.x[h])),d=f.mul(q.weights[h],q.x[h]);
        const Vec& matrix=r.conditional[0][e][h];
        for(int j=0;j<count;++j){
            int value=matrix[ids[j]];
            mass[left][j]=f.add(mass[left][j],f.mul(a,value));
            mass[right][j]=f.add(mass[right][j],f.mul(d,value));
        }
    }
    for(int i=0;i<nv;++i)for(int j=0;j<count;++j){
        uint64_t sum=0;for(int k=0;k<nv;++k)sum+=(uint64_t)g.green[i][k]*mass[k][j];
        greenMass[i][j]=sum%f.p;
    }
    // With the recoil grounded, reciprocity gives -G(p,q)+4G(p,r)-3G(q,r)
    // for hard scalar r=(1,2). The regular Green part interpolates linearly;
    // the same-edge bridge has the preintegrated two-position moments q.bridge.
    Vec integrated(count);int result=0;
    for(int e=0;e<ne;++e)for(int h=0;h<5;++h){
        auto [left,right]=g.edges[e];
        int s=q.x[h],a=f.sub(1,s);
        int drive=f.mul(4,f.add(f.mul(a,g.green[left][4]),f.mul(s,g.green[right][4])));
        for(int j=0;j<count;++j){
            uint64_t bridge=0;
            for(int k=0;k<5;++k)bridge+=(uint64_t)q.bridge[h][k]*r.conditional[0][e][k][ids[j]];
            int v=f.sub(f.mul(drive,total[j]),f.mul(3,greenMass[4][j]));
            v=f.sub(v,f.add(f.mul(a,greenMass[left][j]),f.mul(s,greenMass[right][j])));
            integrated[j]=f.sub(f.mul(q.weights[h],v),f.mul(g.t[e],bridge%f.p));
        }
        const Vec& fixed=p.conditional[0][e][h];
        for(int i=0;i<count;++i)if(fixed[ids[i]]){
            uint64_t value=0;
            for(auto [j,c]:coefficient[i])value+=(uint64_t)c*integrated[j];
            result=f.add(result,f.mul(fixed[ids[i]],value%f.p));
        }
    }
    return f.mul(b.factor,result);
}

int layout(const Core& core,const Field& f,const Quadrature& q,const Vec& times,
           const std::vector<Vec>& byPath,const Vec& pos,bool literal) {
    std::vector<Edge> edges;
    std::vector<Three> cycles;
    Vec lengths;
    for(int e=0;e<6;++e){
        int at=core.edges[e][0],previous=0;
        auto emit=[&](int to,int s){
            edges.push_back({at,to});
            cycles.push_back(core.cycle[e]);
            lengths.push_back(f.mul(times[e],f.sub(s,previous)));
            at=to;previous=s;
        };
        for(int h:byPath[e])emit(4+h,pos[h]);
        emit(core.edges[e][1],1);
    }
    // Hard waves (1,0),(0,1),(1,2), and their recoil. The last two have
    // positions 4 and 5; the grounded potential at the recoil is zero.
    Metric metric(f,edges,cycles,lengths,6,5);
    Mat source=zero(6,5);
    source[4][0]=1;source[4][1]=2;
    source[5][0]=f.val(-2);source[5][1]=f.val(-3);
    std::vector<std::pair<int,int>> owner={{0,0},{1,0},{-1,4},{-1,5}};
    std::vector<Three> probes={{0,3,1},{1,3,0}};
    std::vector<Group> groups(5);
    groups[0].strengths=groups[1].strengths={1};
    for(int k=2;k<5;++k)groups[k].strengths={1,1};
    auto jacobi=[&](Bundle& b){
        return f.sub(f.mul(2,b.potential(0,1)),f.mul(3,b.potential(1,0)));
    };
    if(literal){
        int sum=0;
        for(int mask=0;mask<8;++mask){
            for(int k=2;k<5;++k){
                int v=(mask&(1<<(k-2)))?-1:1;
                groups[k].strengths={1,f.val(v)};
                source[5][k]=f.val(-1-v);
            }
            Bundle b(metric,q,source,groups,owner,probes);
            int v=jacobi(b);
            sum=(pc(mask)%2)?f.sub(sum,v):f.add(sum,v);
        }
        return f.mul(sum,f.inv(64));
    }
    Bundle b(metric,q,source,groups,owner,probes,nullptr,nullptr,7);
    Vec plus=b.groups[2].mean;
    b.groups[2].strengths[1]=f.val(-1);
    b.affineMoments(2);
    for(int i=0;i<(int)plus.size();++i)b.groups[2].mean[i]=f.sub(plus[i],b.groups[2].mean[i]);
    b.groups[3]=b.groups[4]=b.groups[2];
    return f.mul(tangent_voltage(b),f.inv(64));
}

int average(const Core& core,const Field& f,const Quadrature& q,const Vec& times){
    // Each evaluator process has one fixed prime.
    static const auto separate=cubeRules(f,{6,6});
    const auto& innerRule=separate[0];
    static const auto outerRule=rule(f,7);
    std::vector<Vec> byPath(6);
    Vec pos(2);
    int sum=0;
    for(int i=0;i<6;++i)for(int j=0;j<6;++j){
        byPath[i].push_back(0);byPath[j].push_back(1);
        int weight=f.mul(times[i],times[j]);
        if(i!=j){
            for(int a=0;a<4;++a)for(int b=0;b<4;++b){
                pos[0]=separate[0].first[a];pos[1]=separate[1].first[b];
                int w=f.mul(separate[0].second[a],separate[1].second[b]);
                sum=f.add(sum,f.mul(f.mul(weight,w),layout(core,f,q,times,byPath,pos,false)));
            }
        }else{
            // Total position degree on one core path is at most six. The
            // ordering Jacobian supplies one additional outer power.
            for(int outer:{1,0}){
                byPath[i]={1-outer,outer};
                for(int a=0;a<(int)outerRule.first.size();++a){
                    pos[outer]=outerRule.first[a];
                    for(int b=0;b<(int)innerRule.first.size();++b){
                        pos[1-outer]=f.mul(pos[outer],innerRule.first[b]);
                        int w=f.mul(pos[outer],f.mul(outerRule.second[a],innerRule.second[b]));
                        sum=f.add(sum,f.mul(f.mul(weight,w),layout(core,f,q,times,byPath,pos,false)));
                    }
                }
            }
        }
        byPath[i].clear();byPath[j].clear();
    }
    int u=core.kirchhoff(f,times);
    for(int k=0;k<7;++k)sum=f.mul(sum,u);
    return sum;
}
}

int main(int argc,char** argv){
    try{
        if(argc!=3)throw std::runtime_error("usage: fermat_jacobi_core PRIME CORE");
        Field f(std::stoul(argv[1]));Core core(argv[2]);Quadrature q(f);Vec t(6);
        while(std::cin>>t[0]){
            for(int i=1;i<6;++i)if(!(std::cin>>t[i]))return 2;
            for(int& x:t)x=f.val(x);
            std::vector<Vec> paths(6);paths[1]={0};paths[5]={1};
            Vec pos={f.inv(3),f.inv(4)};
            if(layout(core,f,q,t,paths,pos,true)!=layout(core,f,q,t,paths,pos,false))
                throw std::runtime_error("direct eight-sign polarization disagrees");
            std::cout<<average(core,f,q,t)<<std::endl;
        }
    }catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}
}
