// Exact three-loop squared-current averages on a fixed metric cubic core.
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

using Vec=std::vector<int>;
using Mat=std::vector<Vec>;
using Edge=std::array<int,2>;
using Three=std::array<int,3>;
using Wave=std::array<int,5>;
struct Field {
    unsigned p; uint64_t reciprocal;
    explicit Field(unsigned prime):p(prime),reciprocal(prime?UINT64_MAX/prime:0){
        if(prime<11 || prime>10000000)throw std::runtime_error("prime exceeds supported 64-bit moment bounds");
        for(unsigned i=2;i*i<=prime;++i)if(prime%i==0)throw std::runtime_error("modulus is composite");
    }
    int add(int a,int b)const {unsigned z=a+b;return z>=p?z-p:z;}
    int sub(int a,int b)const {return a>=b?a-b:p+a-b;}
    int mul(int a,int b)const {
        uint64_t n=(uint64_t)a*b, q=((__uint128_t)n*reciprocal)>>64;
        unsigned r=n-q*p;return r>=p?r-p:r;
    }
    int val(long long a)const {a%=p;return a<0?a+p:a;}
    int inv(int a)const {
        if(!a)throw std::runtime_error("noninvertible field denominator");
        unsigned e=p-2;int b=1;for(;e;e>>=1,a=mul(a,a))if(e&1)b=mul(a,b);return b;
    }
    int sum(const Vec& a)const{int z=0;for(int v:a)z=add(z,v);return z;}
};
Mat zero(int n,int m){return Mat(n,Vec(m));}
int pc(int x){return __builtin_popcount((unsigned)x);}
int wedgeSign(int a,int b){int z=0;for(;a;a&=a-1)z+=pc((a&-a)-1 & b);return z&1?-1:1;}
Three cross(Three a,Three b){return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};}
int dot(Three a,Three b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
struct Metric {
    const Field& f;std::vector<Edge> edges;std::vector<Three> cycles;Vec t;Mat green,flow;
    Metric(const Field& ff,std::vector<Edge> ed,std::vector<Three> cc,Vec times,int nv,int root):
        f(ff),edges(ed),cycles(cc),t(times),green(zero(nv,nv)),flow(zero(ed.size(),nv)){
        Mat a=zero(nv,2*nv);Vec conductance;
        for(int i=0;i<(int)t.size();++i){
            auto [v,w]=edges[i];int y=f.inv(t[i]);conductance.push_back(y);
            a[v][v]=f.add(a[v][v],y);a[w][w]=f.add(a[w][w],y);
            a[v][w]=f.sub(a[v][w],y);a[w][v]=f.sub(a[w][v],y);
        }
        for(int v=0;v<nv;++v){a[v][root]=a[root][v]=0;a[v][nv+v]=1;}
        a[root][root]=1;
        for(int k=0;k<nv;++k){
            int i=k;while(i<nv&&!a[i][k])++i;
            if(i==nv)throw std::runtime_error("singular Kirchhoff matrix");
            std::swap(a[i],a[k]);int h=f.inv(a[k][k]);
            for(int j=k;j<2*nv;++j)a[k][j]=f.mul(a[k][j],h);
            for(i=0;i<nv;++i)if(i!=k&&a[i][k]){
                h=a[i][k];for(int j=k;j<2*nv;++j)a[i][j]=f.sub(a[i][j],f.mul(h,a[k][j]));
            }
        }
        for(int v=0;v<nv;++v)for(int w=0;w<nv;++w)
            if(v!=root&&w!=root)green[v][w]=a[v][nv+w];
        for(int i=0;i<(int)t.size();++i)for(int w=0;w<nv;++w)
            flow[i][w]=f.mul(conductance[i],f.sub(green[edges[i][0]][w],green[edges[i][1]][w]));
    }
};
struct Group {
    Vec strengths,labels,masks; int degree;
    // Conditional outer products at five interpolation nodes of a distinguished insertion.
    std::vector<std::vector<Mat>> conditional;
    Vec mean;
};
struct Quadrature {
    Vec x,weights;Mat bridge;
    explicit Quadrature(const Field& f):x(5),weights(5),bridge(zero(5,5)){
        Mat lagrange(5);
        for(int k=0;k<5;++k)x[k]=f.mul(k,f.inv(4));
        for(int k=0;k<5;++k){
            Vec q{1};int den=1;
            for(int j=0;j<5;++j)if(j!=k){
                Vec v(q.size()+1);
                for(int h=0;h<(int)q.size();++h){v[h]=f.sub(v[h],f.mul(x[j],q[h]));v[h+1]=f.add(v[h+1],q[h]);}
                q=v;den=f.mul(den,f.sub(x[k],x[j]));
            }
            for(int h=0;h<5;++h){q[h]=f.mul(q[h],f.inv(den));weights[k]=f.add(weights[k],f.mul(q[h],f.inv(h+1)));}
            lagrange[k]=q;
        }
        for(int i=0;i<5;++i)for(int j=0;j<5;++j)for(int a=0;a<5;++a)for(int b=0;b<5;++b){
            int m=f.sub(f.mul(f.add(f.inv(a+2),f.inv(b+2)),f.inv(a+b+3)),f.mul(f.inv(a+2),f.inv(b+2)));
            bridge[i][j]=f.add(bridge[i][j],f.mul(m,f.mul(lagrange[i][a],lagrange[j][b])));
        }
    }
};
struct Env {
    Vec ma,mb;Mat coefficients;
    Vec against(const Field& f,const Vec& b)const{
        Vec out(coefficients.size());
        for(int r=0;r<(int)out.size();++r){
            uint64_t z=0;
            for(int j=0;j<(int)b.size();++j){
                z+=(uint64_t)coefficients[r][j]*b[j];
            }
            out[r]=z%f.p;
        }
        return out;
    }
};
struct Bundle {
    const Field& f;const Metric& g;const Quadrature& q;
    int d=0,factor=0;Vec pivot,free; const Three* softCycle;const Wave* softWave;
    Mat fixed,charges;std::array<bool,5> marked{};mutable std::map<std::pair<int,int>,Env> environments;std::vector<Group> groups;std::vector<std::pair<int,int>> owner;
    Bundle(const Metric& metric,const Quadrature& quad,const Mat& source,
           std::vector<Group> removed,const std::vector<std::pair<int,int>>& slots,const std::vector<Three>& probes,
           const Three* color=nullptr,const Wave* soft=nullptr,int calculate=31):
        f(metric.f),g(metric),q(quad),softCycle(color),softWave(soft),charges(source),groups(removed),owner(slots){
        int e=g.t.size(), nv=g.green.size(), cols=e-1+int(soft!=nullptr);
        for(auto [i,j,k]:probes){marked[k]=true;if(owner[i].first>=0)marked[owner[i].first]=true;
            if(owner[j].first>=0)marked[owner[j].first]=true;}
        for(int k=0;k<5;++k)if(groups[k].strengths.empty()){
            for(int a=0;a<3;++a){
                Vec row(cols);
                for(int h=0;h<e-1;++h){
                    for(int v=0;v<nv;++v)row[h]=f.add(row[h],f.mul(source[v][k],g.flow[h][v]));
                    row[h]=f.mul(row[h],f.val(g.cycles[h][a]));
                }
                if(soft)row.back()=f.mul(f.val((*color)[a]),f.val((*soft)[k]));
                fixed.push_back(row);
            }
        }
        int r=0;factor=1;
        for(int h=0;h<cols;++h){
            int j=r;while(j<(int)fixed.size()&&!fixed[j][h])++j;
            if(j==(int)fixed.size()){free.push_back(h);continue;}
            std::swap(fixed[j],fixed[r]);int val=fixed[r][h];factor=f.mul(factor,f.mul(val,val));
            val=f.inv(val);for(int a=h;a<cols;++a)fixed[r][a]=f.mul(fixed[r][a],val);
            for(j=0;j<(int)fixed.size();++j)if(j!=r&&fixed[j][h]){
                val=fixed[j][h];for(int a=h;a<cols;++a)fixed[j][a]=f.sub(fixed[j][a],f.mul(val,fixed[r][a]));
            }
            pivot.push_back(h);++r;
        }
        if(r!=(int)fixed.size()){factor=0;return;}
        d=free.size();int den=f.inv(g.t.back());factor=f.mul(factor,f.mul(den,den));
        for(int k=0;k<5;++k)if((calculate&(1<<k))&&!groups[k].strengths.empty())groupMoments(k);
    }
    // Integrate the affine edge response and its two-row minors in a shifted Legendre basis.
    void affineMoments(int k){
        Group& gk=groups[k];
        const int m=gk.strengths.size(),ne=g.t.size(),nv=g.green.size();
        const int n=gk.masks.size(),size=n*n;
        const int half=f.inv(2),twelfth=f.inv(12),quadraticVariance=f.inv(180);
        // Each nonempty color supplies at most two residual rows, so d <= 10.
        if(d>10)throw std::logic_error("residual exterior dimension exceeds twice the number of colors");
        using Row=std::array<int,45>;
        using Pair=std::array<Row,2>;
        Vec base(ne-1);
        for(int e=0;e<ne-1;++e)for(int v=0;v<nv;++v)
            base[e]=f.add(base[e],f.mul(charges[v][k],g.flow[e][v]));
        auto linear=[&](Three normal,int i,int j,Row& mid,Pair& slope){
            mid.fill(0);
            for(auto& s:slope)s.fill(0);
            for(int l=0;l<ne-1;++l){
                int sign=dot(normal,g.cycles[l]);if(!sign)continue;
                int value=base[l],tangent=0;
                for(int h=0;h<m;++h){
                    auto [a,b]=g.edges[h?j:i];
                    int charge=gk.strengths[h];
                    value=f.add(value,f.mul(charge,g.flow[l][a]));
                    int delta=f.mul(charge,f.sub(g.flow[l][b],g.flow[l][a]));
                    tangent=f.add(tangent,delta);slope[h][l]=f.mul(f.val(sign),delta);
                }
                mid[l]=f.mul(f.val(sign),f.add(value,f.mul(half,tangent)));
            }
            if(softWave)mid[ne-1]=f.mul(f.val(dot(normal,*softCycle)),f.val((*softWave)[k]));
            if(!pivot.empty()){
                auto project=[&](Row& row){
                    Row raw=row;
                    for(int a=0;a<d;++a){
                        int value=raw[free[a]];
                        for(int j=0;j<(int)pivot.size();++j)value=f.sub(value,f.mul(fixed[j][free[a]],raw[pivot[j]]));
                        row[a]=value;
                    }
                };
                project(mid);
                for(int i=0;i<m;++i)project(slope[i]);
            }
        };
        using Wide=int64_t;
        std::vector<Wide> means(size);
        std::vector<Wide> conditional(marked[k]&&m==2?m*ne*3*size:0);
        if(marked[k]&&m==1)gk.conditional.assign(1,std::vector<Mat>(ne,zero(5,size)));
        auto outer=[&](Wide* out,const Row& a,const Row& b,int weight){
            for(int i=0;i<n;++i){
                Wide x=f.mul(weight,a[i]);
                for(int j=i;j<n;++j)out[i*n+j]+=x*b[j];
            }
        };
        auto symmetric=[&](Wide* out,const Row& a,const Row& b,int weight){
            for(int i=0;i<n;++i){
                Wide x=f.mul(weight,a[i]),y=f.mul(weight,b[i]);
                for(int j=i;j<n;++j)out[i*n+j]+=x*b[j]+y*a[j];
            }
        };
        int jump=1;for(int z:gk.strengths)jump=f.mul(jump,z);
        jump=f.mul(jump,jump);
        for(int i=0;i<ne;++i)for(int j=0;j<(m==2?ne:1);++j){
            int vol=f.mul(jump,g.t[i]);
            if(m==2){
                Three normal=cross(g.cycles[i],g.cycles[j]);
                if(normal==Three{})continue;
                vol=f.mul(vol,g.t[j]);int variance=f.mul(vol,twelfth);
                Row mid;Pair slope;linear(normal,i,j,mid,slope);
                outer(means.data(),mid,mid,vol);
                for(auto& s:slope)outer(means.data(),s,s,variance);
                if(marked[k])for(int h=0;h<2;++h){
                    Wide* c=&conditional[(h*ne+(h?j:i))*3*size];
                    outer(c,mid,mid,vol);outer(c,slope[1-h],slope[1-h],variance);
                    symmetric(c+size,mid,slope[h],vol);outer(c+2*size,slope[h],slope[h],vol);
                }
            }else{
                Three cyc=g.cycles[i];int pivot=0;while(pivot<3&&!cyc[pivot])++pivot;
                if(pivot==3)continue;
                Row midpoint[2];Pair slopes[2];int r=0;
                for(int h=0;h<3;++h)if(h!=pivot){
                    Three normal{};normal[h]=cyc[pivot];normal[pivot]=-cyc[h];
                    linear(normal,i,-1,midpoint[r],slopes[r]);++r;
                }
                Row mean{},derivative{},quadratic{};
                const Row &a=midpoint[0],&b=midpoint[1],&u=slopes[0][0],&v=slopes[1][0];
                for(int y=1,l=0;y<d;++y)for(int x=0;x<y;++x,++l){
                    auto wedge=[&](const Row& e,const Row& h){return f.sub(f.mul(e[x],h[y]),f.mul(e[y],h[x]));};
                    quadratic[l]=wedge(u,v);
                    derivative[l]=f.add(wedge(a,v),wedge(u,b));
                    mean[l]=f.add(wedge(a,b),f.mul(twelfth,quadratic[l]));
                }
                outer(means.data(),mean,mean,vol);
                outer(means.data(),derivative,derivative,f.mul(vol,twelfth));
                outer(means.data(),quadratic,quadratic,f.mul(vol,quadraticVariance));
                if(marked[k]){
                    for(int h=0;h<5;++h){
                        int z=f.sub(q.x[h],half),second=f.sub(f.mul(z,z),twelfth);
                        Row value{};for(int l=0;l<n;++l)value[l]=f.add(mean[l],f.add(f.mul(z,derivative[l]),f.mul(second,quadratic[l])));
                        Vec& out=gk.conditional[0][i][h];
                        for(int a=0;a<n;++a)for(int b=a;b<n;++b)out[a*n+b]=out[b*n+a]=f.mul(vol,f.mul(value[a],value[b]));
                    }
                }
            }
        }
        for(int a=0;a<n;++a)for(int b=a;b<n;++b)gk.mean[a*n+b]=gk.mean[b*n+a]=means[a*n+b]%f.p;
        if(marked[k]&&m==2){
            gk.conditional.assign(m,std::vector<Mat>(ne,zero(5,size)));
            for(int side=0;side<2;++side)for(int i=0;i<ne;++i){
                Wide* c=&conditional[(side*ne+i)*3*size];
                for(int a=0;a<n;++a)for(int b=a;b<n;++b){
                    int l=a*n+b,constant=c[l]%f.p,first=c[size+l]%f.p,second=c[2*size+l]%f.p;
                    for(int h=0;h<5;++h){
                        int z=f.sub(q.x[h],half),value=f.add(constant,f.mul(z,f.add(first,f.mul(z,second))));
                        auto& out=gk.conditional[side][i][h];out[l]=out[b*n+a]=value;
                    }
                }
            }
        }
    }
    void groupMoments(int k){
        Group& b=groups[k];int m=b.strengths.size(),e=g.t.size(),nv=g.green.size();b.degree=3-m;
        for(int mask=0;mask<(1<<d);++mask)if(pc(mask)==b.degree)b.masks.push_back(mask);
        int l=b.masks.size(),size=l*l;b.mean.assign(size,0);
        if(m==3){
            Vec mass(e);int total=0;
            int scale=1;for(int a:b.strengths)scale=f.mul(scale,f.mul(a,a));
            for(int i=0;i<e;++i)for(int j=i+1;j<e;++j)for(int h=j+1;h<e;++h){
                int det=dot(g.cycles[i],cross(g.cycles[j],g.cycles[h]));if(!det)continue;
                int w=f.mul(scale,f.mul(f.val(det*det),f.mul(g.t[i],f.mul(g.t[j],g.t[h]))));
                total=f.add(total,f.mul(6,w));
                for(int edge:{i,j,h})mass[edge]=f.add(mass[edge],f.mul(2,w));
            }
            b.mean[0]=total;
            if(marked[k]){
                b.conditional.assign(3,std::vector<Mat>(e,zero(5,1)));
                for(auto& label:b.conditional)for(int i=0;i<e;++i)for(Vec& row:label[i])row[0]=mass[i];
            }
            return;
        }
        affineMoments(k);
    }
    Env environment(int a,int b)const;
    int potential(int target,int load,int difference=-1)const;
};
Env Bundle::environment(int a,int b)const{
    int n=1<<d;
    struct Split {int left,right,sign;};
    struct Schedule {Vec masks;std::vector<std::vector<Split>> options;};
    static std::map<Three,Schedule> plans;
    Mat dp{{1}};Vec masks{0};int degree=0;
    for(int k=0;k<5;++k)if(k!=a&&k!=b&&!groups[k].strengths.empty()){
        const Group& x=groups[k];int z=x.masks.size();Three key{d,degree,x.degree};
        if(!plans.count(key)){
            Schedule schedule;
            for(int m=0;m<n;++m)if(pc(m)==degree+x.degree){
                std::vector<Split> splits;
                for(int j=0;j<z;++j)if((x.masks[j]&m)==x.masks[j]){
                    int left=m^x.masks[j];
                    int id=std::lower_bound(masks.begin(),masks.end(),left)-masks.begin();
                    splits.push_back({id,j,wedgeSign(left,x.masks[j])});
                }
                schedule.masks.push_back(m);schedule.options.push_back(splits);
            }
            plans[key]=std::move(schedule);
        }
        const auto& schedule=plans.at(key);int size=schedule.masks.size();Mat next=zero(size,size);
        for(int i=0;i<size;++i)for(int j=i;j<size;++j){
            int64_t sum=0;
            for(auto [l,r,sign]:schedule.options[i]){
                const int* lhs=dp[l].data();const int* rhs=&x.mean[r*z];
                for(auto w:schedule.options[j])sum+=(int64_t)(sign*w.sign)*lhs[w.left]*rhs[w.right];
            }
            next[i][j]=next[j][i]=f.val(sum);
        }
        dp=std::move(next);masks=schedule.masks;degree+=x.degree;
    }
    Vec indices(n);for(int i=0;i<(int)masks.size();++i)indices[masks[i]]=i;
    Env out;
    out.ma=a<0?Vec{0}:groups[a].masks;out.mb=b<0?Vec{0}:groups[b].masks;
    int na=out.ma.size(),nb=out.mb.size();out.coefficients=zero(na*na,nb*nb);
    for(int i=0;i<na;++i)for(int k=0;k<nb;++k){
        int u=out.ma[i],v=out.mb[k];if(u&v)continue;
        int rem=(n-1)^(u|v),left=wedgeSign(rem,u)*wedgeSign(rem|u,v);
        for(int j=0;j<na;++j)for(int l=0;l<nb;++l){
            int x=out.ma[j],y=out.mb[l];if(x&y)continue;
            int end=(n-1)^(x|y);int right=wedgeSign(end,x)*wedgeSign(end|x,y);
            int val=dp[indices[rem]][indices[end]];out.coefficients[i*na+j][k*nb+l]=left==right?val:f.sub(0,val);
        }
    }
    return out;
}
int inner(const Field& f,const Vec& a,const Vec& b){
    uint64_t z=0;for(int i=0;i<(int)a.size();++i){z+=(uint64_t)a[i]*b[i];}return z%f.p;
}
int Bundle::potential(int target,int load,int difference)const{
    if(!factor)return 0;
    int a=owner.at(target).first,pos=owner.at(target).second;int b=groups[load].strengths.empty()?-1:load;
    if(a>=0&&a==b)throw std::runtime_error("the requested target is longitudinal to this component");
    auto key=std::make_pair(a,b);
    if(!environments.count(key))environments[key]=environment(a,b);
    const Env& env=environments.at(key);int nv=g.green.size(),ne=g.t.size();
    std::vector<Mat> delta;
    if(difference>=0){
        delta=groups[a].conditional[pos];
        for(int e=0;e<ne;++e)for(int h=0;h<5;++h)for(int j=0;j<(int)delta[e][h].size();++j)
            delta[e][h][j]=f.sub(delta[e][h][j],groups[a].conditional[difference][e][h][j]);
    }
    auto observed=[&](int e,int h)->const Vec&{return difference<0?groups[a].conditional[pos][e][h]:delta[e][h];};
    int la=env.ma.size()*env.ma.size(),lb=env.mb.size()*env.mb.size();
    Vec dual=env.against(f,b<0?Vec{1}:groups[b].mean);
    Mat measureA=zero(nv,la),measureB=zero(nv,lb);
    if(a<0)measureA[pos][0]=1;
    int bridge=0;
    for(int e=0;e<ne;++e){
        auto [left,right]=g.edges[e];Mat tuned;
        for(int h=0;h<5;++h){
            int w[2]={f.mul(q.weights[h],f.sub(1,q.x[h])),f.mul(q.weights[h],q.x[h])};
            if(a>=0){
                const Vec& v=observed(e,h);
                for(int x=0;x<la;++x){
                    measureA[left][x]=f.add(measureA[left][x],f.mul(w[0],v[x]));
                    measureA[right][x]=f.add(measureA[right][x],f.mul(w[1],v[x]));
                }
            }
            if(b>=0){
                Vec v(lb);
                for(int j=0;j<(int)groups[b].strengths.size();++j){
                    const Vec& c=groups[b].conditional[j][e][h];int strength=groups[b].strengths[j];
                    for(int x=0;x<lb;++x)v[x]=f.add(v[x],f.mul(strength,c[x]));
                }
                for(int x=0;x<lb;++x){
                    measureB[left][x]=f.add(measureB[left][x],f.mul(w[0],v[x]));
                    measureB[right][x]=f.add(measureB[right][x],f.mul(w[1],v[x]));
                }
                if(a>=0)tuned.push_back(env.against(f,v));
            }
        }
        if(a>=0&&b>=0)for(int h=0;h<5;++h)for(int j=0;j<5;++j){
            int term=inner(f,observed(e,h),tuned[j]);
            bridge=f.add(bridge,f.mul(g.t[e],f.mul(q.bridge[h][j],term)));
        }
    }
    int answer=bridge;
    for(int w=0;w<nv;++w){
        Vec bvertex=b>=0?env.against(f,measureB[w]):Vec(la);
        for(int j=0;j<la;++j)bvertex[j]=f.add(bvertex[j],f.mul(charges[w][load],dual[j]));
        for(int v=0;v<nv;++v)
            answer=f.add(answer,f.mul(g.green[v][w],inner(f,measureA[v],bvertex)));
    }
    return f.mul(factor,answer);
}
// In the final quadratic-soft probe the soft jump is carried only by colors
// 2 and 3. Contract its two occurrences while wedging those two spectator blocks.
struct PairingPlan {
    struct Term {int left,right,sign;};
    std::vector<std::vector<Term>> terms;Vec masks;
    PairingPlan(const Vec& left,const Vec& right,int dim,int degree,bool ghost){
        for(int m=0;m<(1<<dim);++m)if(pc(m)==degree && (!ghost||(m&(1<<(dim-1))))){
            std::vector<Term> opts;
            for(int i=0;i<(int)left.size();++i)if((m&left[i])==left[i]){
                int other=m^left[i];
                int j=std::lower_bound(right.begin(),right.end(),other)-right.begin();
                if(j<(int)right.size()&&right[j]==other)opts.push_back({i,j,wedgeSign(left[i],other)});
            }
            masks.push_back(ghost?m^(1<<(dim-1)):m);terms.push_back(opts);
        }
    }
    void add(const Field& f,const Vec& left,int nl,const Vec& right,int nr,int weight,Vec& out)const{
        int n=masks.size();
        for(int i=0;i<n;++i)for(int j=i;j<n;++j){
            int64_t sum=0;
            for(auto a:terms[i]){
                const int* x=&left[a.left*nl];const int* y=&right[a.right*nr];
                for(auto b:terms[j])sum+=(int64_t)(a.sign*b.sign)*x[b.left]*y[b.right];
            }
            int z=f.add(out[i*n+j],f.mul(weight,f.val(sum)));out[i*n+j]=out[j*n+i]=z;
        }
    }
};
struct SoftPairing {
    PairingPlan withGhost,remaining;
    static Vec masks(int n,int degree){
        Vec result;for(int i=0;i<(1<<n);++i)if(pc(i)==degree)result.push_back(i);return result;
    }
    SoftPairing():withGhost(masks(8,2),masks(8,2),8,4,true),
                  remaining(withGhost.masks,masks(8,2),7,5,false){}
    Env finish(const Bundle& b,const Vec& first)const{
        const Field& f=b.f;const Group& extra=b.groups[4];
        Vec top(21*21);remaining.add(f,first,35,extra.mean,28,1,top);
        Vec index(128);for(int i=0;i<21;++i)index[remaining.masks[i]]=i;
        Env env;env.ma=b.groups[1].masks;env.mb=b.groups[0].masks;env.coefficients=zero(64,64);
        for(int i=0;i<7;++i)for(int k=0;k<7;++k)if(i!=k){
            int comp=127^(1<<i)^(1<<k),s=wedgeSign(comp,1<<i)*wedgeSign(comp|(1<<i),1<<k);
            for(int j=0;j<7;++j)for(int l=0;l<7;++l)if(j!=l){
                int other=127^(1<<j)^(1<<l);
                int t=wedgeSign(other,1<<j)*wedgeSign(other|(1<<j),1<<l);
                int value=top[21*index[comp]+index[other]];
                env.coefficients[8*i+j][8*k+l]=s==t?value:f.sub(0,value);
            }
        }
        return env;
    }
};
struct Core {
    std::vector<Edge> edges;std::vector<Three> cycle;
    explicit Core(const std::string& name){
        if(name=="tetrahedron"){
            edges={{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
            cycle={{1,1,0},{-1,0,1},{0,-1,-1},{1,0,0},{0,1,0},{0,0,1}};
        }else if(name=="double-square"){
            edges={{0,1},{0,1},{0,3},{1,2},{2,3},{2,3}};
            cycle={{-1,1,1},{1,0,0},{0,-1,-1},{0,1,1},{0,1,0},{0,0,1}};
        }else throw std::runtime_error("unknown core");
    }
    int kirchhoff(const Field& f,const Vec& times)const{
        int z=0;
        for(int a=0;a<6;++a)for(int b=a+1;b<6;++b)for(int c=b+1;c<6;++c){
            int v=dot(cycle[a],cross(cycle[b],cycle[c]));
            z=f.add(z,f.mul(f.val(v*v),f.mul(times[a],f.mul(times[b],times[c]))));
        }
        return z;
    }
};
std::vector<Wave> makeFixture(const std::string& name){
    std::vector<Wave> waves;
    auto axis=[&](int k,int a){Wave x{};x[k]=a;waves.push_back(x);};
    if(name=="control"){
        axis(0,1);axis(0,1);axis(0,1);axis(1,1);axis(1,1);axis(2,1);axis(2,1);axis(3,1);axis(4,1);
    }else if(name=="four"||name=="equal"){
        for(int k=0;k<4;++k){axis(k,1);axis(k,name=="equal"?1:k+2);}axis(4,1);
    }else if(name=="three"){
        axis(0,1);axis(0,2);axis(1,1);axis(1,3);axis(2,1);axis(2,-1);axis(3,1);axis(4,1);
        waves.push_back({0,-2,-1,-1,1});
    }else if(name=="soft-two"||name.rfind("soft-full=",0)==0){
        axis(0,1);axis(0,2);axis(1,1);axis(1,3);axis(2,1);axis(3,1);axis(4,1);
        int t=name=="soft-two"?1:std::stoi(name.substr(10));
        waves.push_back({0,0,t,t,0});waves.push_back({0,0,1,2,1});
    }else if(name=="two"){
        axis(0,1);axis(0,2);axis(1,1);axis(1,3);axis(2,1);axis(3,1);axis(4,1);
        waves.push_back({1,0,1,-2,-2});waves.push_back({-2,0,1,1,-2});
    }else throw std::runtime_error("unknown source fixture");
    Wave balanced{};for(auto x:waves)for(int k=0;k<5;++k)balanced[k]-=x[k];waves.push_back(balanced);
    return waves;
}
struct Oracle {
    const Field& f;const Core& core;const Quadrature& q;const Vec times;std::vector<Wave> waves;
    Vec live;std::vector<Group> groups;std::vector<std::pair<int,int>> owner;
    std::vector<Three> probes;int softLabel;std::vector<std::pair<Three,int>> softDirections;
    Oracle(const Field& field,const Core& graph,const Quadrature& quad,Vec t,std::vector<Wave> p,
           std::vector<Three> tests,int keep=0,int soft=-1):f(field),core(graph),q(quad),times(t),waves(p),groups(5),owner(p.size()),probes(tests),softLabel(soft){
        if(soft>=0){
            for(int k=0;k<5;++k)waves.back()[k]+=p[soft][k];
            for(auto [i,j,k]:tests)if(i==soft||j==soft||p[soft][k])throw std::runtime_error("soft probe carries a voltage insertion");
            Mat gram=zero(3,3);
            for(int e=0;e<6;++e)for(int i=0;i<3;++i)for(int j=0;j<3;++j)
                gram[i][j]=f.add(gram[i][j],f.mul(t[e],f.val(core.cycle[e][i]*core.cycle[e][j])));
            for(int step=0;step<3;++step){
                int i=0;while(i<3&&!gram[i][i])++i;
                if(i==3){softDirections.clear();for(int e=0;e<6;++e)softDirections.push_back({core.cycle[e],t[e]});break;}
                Three col{gram[0][i],gram[1][i],gram[2][i]};int weight=f.inv(col[i]);
                softDirections.push_back({col,weight});
                for(int a=0;a<3;++a)for(int b=0;b<3;++b)
                    gram[a][b]=f.sub(gram[a][b],f.mul(weight,f.mul(col[a],col[b])));
            }
        }
        for(auto [i,j,k]:tests)if(p.at(i).at(k)||p.at(j).at(k))throw std::runtime_error("probe must be transverse on both observed slots");
        for(int v=0;v<(int)p.size();++v){
            if(v==soft){owner[v]={-2,0};continue;}
            int n=0,g=-1;for(int h=0;h<5;++h)if(p[v][h]){++n;g=h;}
            if(n==1&&v+1!=(int)p.size()&&!(keep&(1<<v))){
                owner[v]={g,(int)groups[g].strengths.size()};groups[g].strengths.push_back(f.val(p[v][g]));
            }else{owner[v]={-1,(int)live.size()+4};live.push_back(v);}
        }
        for(const auto& g:groups)if(g.strengths.size()>3)throw std::runtime_error("expected at most three pure insertions per coordinate");
    }
    Vec layout(const std::vector<Vec>& byPath,const Vec& fraction)const{
        std::vector<Edge> ed;std::vector<Three> cycles;Vec t;Mat src=zero(4+live.size(),5);
        for(int h=0;h<(int)live.size();++h)for(int k=0;k<5;++k)src[4+h][k]=f.val(waves[live[h]][k]);
        for(int e=0;e<6;++e){
            int at=core.edges[e][0],previous=0;
            auto emit=[&](int vertex,int pos){
                ed.push_back({at,vertex});cycles.push_back(core.cycle[e]);t.push_back(f.mul(times[e],f.sub(pos,previous)));
                at=vertex;previous=pos;
            };
            for(int h:byPath[e])emit(4+h,fraction[h]);
            emit(core.edges[e][1],1);
        }
        Metric metric(f,ed,cycles,t,src.size(),4+live.size()-1);
        auto evaluate=[&](const Three* color,const Wave* soft){
            Bundle bundle(metric,q,src,groups,owner,probes,color,soft);
            int cache[10][5];for(auto& v:cache)std::fill(v,v+5,-1);
            auto at=[&](int v,int k){if(cache[v][k]<0)cache[v][k]=bundle.potential(v,k);return cache[v][k];};
            Vec result;
            for(auto [i,j,k]:probes){
                if(probes.size()==1 && owner[i].first>=0 && owner[i].first==owner[j].first)
                    result.push_back(bundle.potential(i,k,owner[j].second));
                else result.push_back(f.sub(at(i,k),at(j,k)));
            }
            return result;
        };
        if(softLabel==7 && probes==std::vector<Three>{{2,3,0}} && waves[softLabel]==Wave{0,0,1,1,0}){
            static const SoftPairing plan;
            Vec contracted(35*35);
            std::unique_ptr<Bundle> common;
            for(const auto& [color,weight]:softDirections){
                auto b=std::make_unique<Bundle>(metric,q,src,groups,owner,probes,&color,&waves[softLabel],
                                                common?((1<<2)|(1<<3)):31);
                if(!b->factor)return Vec{0};
                plan.withGhost.add(f,b->groups[2].mean,28,b->groups[3].mean,28,weight,contracted);
                if(!common)common=std::move(b);
            }
            common->environments[{1,0}]=plan.finish(*common,contracted);
            return Vec{common->potential(2,0,owner[3].second)};
        }
        if(softLabel<0)return evaluate(nullptr,nullptr);
        Vec result(probes.size());
        for(const auto& [color,weight]:softDirections){
            Vec v=evaluate(&color,&waves[softLabel]);
            for(int i=0;i<(int)v.size();++i)result[i]=f.add(result[i],f.mul(weight,v[i]));
        }
        return result;
    }
    Vec average()const;
    Vec orderedAverage()const;
};
// Interior rational interpolation is used because every subdivided metric is nonsingular.
std::pair<Vec,Vec> rule(const Field& f,int degree){
    int n=degree+1;Vec nodes(n),weights(n);
    for(int i=0;i<n;++i)nodes[i]=f.mul(i+1,f.inv(n+1));
    for(int i=0;i<n;++i){
        Vec poly{1};int denom=1;
        for(int j=0;j<n;++j)if(i!=j){
            Vec next(poly.size()+1);
            for(int h=0;h<(int)poly.size();++h){
                next[h]=f.sub(next[h],f.mul(nodes[j],poly[h]));next[h+1]=f.add(next[h+1],poly[h]);
            }
            poly=next;denom=f.mul(denom,f.sub(nodes[i],nodes[j]));
        }
        int w=0;for(int h=0;h<n;++h)w=f.add(w,f.mul(poly[h],f.inv(h+1)));
        weights[i]=f.mul(w,f.inv(denom));
    }
    return {nodes,weights};
}
std::vector<std::pair<Vec,Vec>> cubeRules(const Field& f,const Vec& degrees){
    // A polynomial P_n+a P_{n-1} with n distinct roots gives a degree-(2n-2) rule.
    std::vector<std::pair<Vec,Vec>> rules;Vec used(f.p);
    for(int degree:degrees){
        int points=(degree+2)/2;
        std::vector<int> heads(f.p,-1),size(f.p),next(f.p,-1);
        bool found=false;
        for(int x=2;x<(int)f.p&&!found;++x){
            if(used[x])continue;
            int y=f.sub(f.mul(2,x),1),y2=f.mul(y,y),y3=f.mul(y2,y);
            int p2=f.sub(f.mul(3,y2),1),p3=f.sub(f.mul(5,y3),f.mul(3,y));
            int p4=f.add(3,f.sub(f.mul(35,f.mul(y2,y2)),f.mul(30,y2)));
            int lower=points==4?p3:p2,upper=points==4?p4:p3;
            if(!lower)continue;
            int a=f.mul(f.sub(0,upper),f.inv(lower));
            next[x]=heads[a];heads[a]=x;if(++size[a]!=points)continue;
            Vec nodes;for(int i=x;i>=0;i=next[i]){nodes.push_back(i);used[i]=1;}
            Vec weights(points);
            for(int i=0;i<points;++i){
                Vec polynomial{1};int denominator=1;
                for(int j=0;j<points;++j)if(i!=j){
                    Vec out(polynomial.size()+1);
                    for(int h=0;h<(int)polynomial.size();++h){
                        out[h]=f.sub(out[h],f.mul(nodes[j],polynomial[h]));out[h+1]=f.add(out[h+1],polynomial[h]);
                    }
                    polynomial=out;denominator=f.mul(denominator,f.sub(nodes[i],nodes[j]));
                }
                int weight=0;for(int h=0;h<points;++h)weight=f.add(weight,f.mul(polynomial[h],f.inv(h+1)));
                weights[i]=f.mul(weight,f.inv(denominator));
            }
            Vec powers(points,1);
            for(int k=0;k<=degree;++k){
                if(inner(f,weights,powers)!=f.inv(k+1))throw std::runtime_error("inexact core-position cubature");
                for(int j=0;j<points;++j)powers[j]=f.mul(powers[j],nodes[j]);
            }
            rules.push_back({nodes,weights});found=true;
        }
        if(!found)throw std::runtime_error("insufficient finite-field position rules");
    }
    return rules;
}
Vec Oracle::average()const{
    int m=live.size();std::vector<Vec> byPath(6);Vec pos(m),sum(probes.size());
    Vec degrees(m,4);
    for(auto [i,j,k]:probes)for(int h=0;h<m;++h)
        if(owner[i].first<0||owner[j].first<0||waves[live[h]][k])degrees[h]=6;
    static std::map<std::pair<int,Vec>,std::vector<std::pair<Vec,Vec>>> saved;
    auto key=std::make_pair((int)f.p,degrees);
    if(!saved.count(key))saved[key]=cubeRules(f,degrees);
    auto& grid=saved[key];
    for(auto [i,j,k]:probes)if((owner[i].first<0||owner[j].first<0)&&m>1){
        // A retained target can have a piecewise-polynomial bridge to a retained source.
        // Sorting the simplex positions then requires the ordered rule.
        bool active=false;
        for(int h:live)if(h!=i&&h!=j&&waves[h][k])active=true;
        if(active)return orderedAverage();
    }
    std::function<void(int,int)> assign=[&](int j,int weight){
        if(j==m){Vec v=layout(byPath,pos);for(int i=0;i<(int)v.size();++i)sum[i]=f.add(sum[i],f.mul(weight,v[i]));return;}
        for(int e=0;e<6;++e){
            byPath[e].push_back(j);
            for(int h=0;h<(int)grid[j].first.size();++h){
                pos[j]=grid[j].first[h];
                assign(j+1,f.mul(weight,f.mul(times[e],grid[j].second[h])));
            }
            byPath[e].pop_back();
        }
    };
    assign(0,1);
    return sum;
}
Vec Oracle::orderedAverage()const{
    int m=live.size();std::vector<Vec> byPath(6);Vec pos(m),sum(probes.size());
    std::vector<std::pair<Vec,Vec>> grid(m);
    for(int j=0;j<m;++j)grid[j]=rule(f,6+j);
    std::function<void(int,int)> integrate=[&](int e,int weight){
        if(e==6){Vec v=layout(byPath,pos);for(int i=0;i<(int)v.size();++i)sum[i]=f.add(sum[i],f.mul(weight,v[i]));return;}
        Vec& sites=byPath[e];if(sites.empty()){integrate(e+1,weight);return;}
        do{
            // Integrate 0 < s_1 < ... < s_l < 1 using nested open unit variables.
            std::function<void(int,int,int)> descend=[&](int j,int upper,int wt){
                if(j<0){integrate(e+1,wt);return;}
                const auto& [nodes,weights]=grid[j];
                for(int h=0;h<(int)nodes.size();++h){
                    pos[sites[j]]=f.mul(upper,nodes[h]);
                    int local=f.mul(upper,f.mul(times[e],weights[h]));
                    descend(j-1,pos[sites[j]],f.mul(wt,local));
                }
            };
            descend(sites.size()-1,1,weight);
        }while(std::next_permutation(sites.begin(),sites.end()));
    };
    std::function<void(int)> assign=[&](int j){
        if(j==m){integrate(0,1);return;}
        for(int e=0;e<6;++e){byPath[e].push_back(j);assign(j+1);byPath[e].pop_back();}
    };
    assign(0);
    return sum;
}
#ifndef PARTIAL_COTREES_LIBRARY
int main(int argc,char** argv){
    try{
        if(argc==3&&std::string(argv[1])=="--check"){
            Field f(std::stoul(argv[2]));Core core("tetrahedron");Quadrature q(f);
            Oracle first(f,core,q,{5,9,13,17,10,11},makeFixture("three"),{{0,6,4}},(1<<6)|(1<<7));
            Vec pos;for(auto [a,b]:std::vector<Edge>{{2,5},{4,9},{6,13},{8,17}})pos.push_back(f.mul(a,f.inv(b)));
            int value=first.layout({{0},{1},{2},{3},{},{}},pos).at(0);
            std::cout<<"rational-layout "<<value<<std::endl;
            Oracle base(f,core,q,{1,2,3,4,5,6},makeFixture("four"),{{0,8,1}});
            Oracle extra(f,core,q,{1,2,3,4,5,6},makeFixture("four"),{{0,8,1}},1<<8);
            int v=base.average().at(0),w=extra.average().at(0);
            std::cout<<"one or two retained labels "<<v<<" "<<w<<std::endl;if(v!=w)return 3;
            return 0;
        }
        if(argc<7||(argc-4)%3){std::cerr<<"usage: partial-cotrees PRIME CORE {four|equal|three|two} LEFT RIGHT LOAD\n";return 2;}
        Field f(std::stoul(argv[1]));Core core(argv[2]);Quadrature q(f);
        auto waves=makeFixture(argv[3]);std::vector<Three> tests;
        for(int i=4;i<argc;i+=3)tests.push_back({std::stoi(argv[i]),std::stoi(argv[i+1]),std::stoi(argv[i+2])});
        Vec t(6);
        while(std::cin>>t[0]){
            for(int i=1;i<6;++i)if(!(std::cin>>t[i]))return 2;
            for(int& z:t)z=f.val(z);
            int u=core.kirchhoff(f,t);if(!u)throw std::runtime_error("zero Kirchhoff polynomial at core node");
            Oracle oracle(f,core,q,t,waves,tests,0,std::string(argv[3])=="soft-two"?7:-1);
            Vec values=oracle.average();int power=1;for(int i=0;i<7;++i)power=f.mul(power,u);
            for(int i=0;i<(int)values.size();++i)std::cout<<(i?" ":"")<<f.mul(values[i],power);
            std::cout<<std::endl;
        }
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
#endif
