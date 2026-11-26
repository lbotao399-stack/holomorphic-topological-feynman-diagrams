// Completed planar cut charts and convergence-safe position-space IBP in even dimension.
#include "quartic_prism_graph.hpp"
#include <chrono>
#include <fstream>
#include <set>
#include <random>
#include <cstring>
#include <map>
#include <gmpxx.h>
#include <iostream>
#include <stdexcept>
#include <unordered_map>
#include <vector>

using Factors = std::array<short,4>;
constexpr std::array<int,4> primes{2,3,5,7};

Factors change(Factors a, const Factors &b, int sign=1) {
    for (int i=0; i<4; ++i) a[i] += sign*b[i];
    return a;
}

struct Hash {
    std::size_t operator()(const Key &a) const {
        std::uint64_t h = 1469598103934665603ULL;
        for (unsigned char c : a) h = (h^c)*1099511628211ULL;
        return h;
    }
};
struct Label {
    int index;
    Factors ratio; // The integral of the canonical graph divided by its representative's integral.
    int record = -1;
};

bool convergent(const Edges &a, int q) {
    for (int i=0;i<6;++i) for (int j=i+1;j<6;++j) {
        if (a[i][j]>=q) return false;
        for (int k=j+1;k<6;++k) if (a[i][j]+a[i][k]+a[j][k]>=2*q) return false;
    }
    return true;
}

struct Domain {
    int q;
    bool keep_boundaries;
    std::vector<Factors> star_factor;
    std::unordered_map<Key,Label,Hash> labels;
    std::vector<Key> representatives;

    explicit Domain(int dim,bool full=false):q(dim),keep_boundaries(full),star_factor(q) {
        if (q!=6 && q!=12) throw std::invalid_argument("q must be 6 or 12");
        for (int i=1;i<q;++i) {
            auto factorial_factor=[](int n) {
                Factors r{};
                for (int k=0;k<4;++k) for (int p=primes[k];p<=n;p*=primes[k]) r[k] += n/p;
                return r;
            };
            star_factor[i]=change(factorial_factor(q-i-1),factorial_factor(i-1),-1);
        }
        labels.reserve(q==6 ? 6000 : 4000000);
    }

    void register_chart(const Key &seed) {
        if (labels.find(seed)!=labels.end()) return;
        std::vector<Key> queue{seed};
        labels.emplace(seed,Label{-2,{}});
        bool vanishes=false;
        for (std::size_t t=0;t<queue.size();++t) {
            const Key g=queue[t];
            const Edges a=unpack(g);
            Factors factor=labels.at(g).ratio;
            vanishes |= !keep_boundaries && has_bubble(a);
            for (int v=0;v<6;++v) {
                std::array<int,5> neighbors{};
                int n=0;
                bool positive=true;
                for (int u=0;u<6;++u) if (a[v][u]) {
                    neighbors[n++]=u;
                    positive &= a[v][u]>0 && a[v][u]<q;
                }
                if (n!=4 || !positive) continue;
                Edges b=a;
                Factors f=factor;
                for (int j=0;j<n;++j) {
                    int u=neighbors[j], w=a[v][u];
                    b[u][v]=b[v][u]=q-w;
                    f=change(f,star_factor[w],-1);
                    for (int k=0;k<j;++k) {
                        int z=neighbors[k];
                        b[u][z]=b[z][u]=a[u][z]+w+a[v][z]-q;
                    }
                }
                Key child=graph_key(b);
                auto [it,fresh]=labels.emplace(child,Label{-2,f});
                if (fresh) queue.push_back(child);
                else if (it->second.index!=-2 || it->second.ratio!=f)
                    throw std::logic_error("inconsistent reflection");
            }
        }
        int idx=vanishes ? -1 : int(representatives.size());
        if (!vanishes) representatives.push_back(seed);
        for (const Key &k:queue) labels.at(k).index=idx;
    }

    void enumerate() {
        std::vector<std::array<int,4>> stars,sorted;
        for (int a=3-q;a<q;++a) for (int b=3-q;b<q;++b) for (int c=3-q;c<q;++c) {
            int d=2*q-a-b-c;
            if (d<3-q || d>=q) continue;
            std::array<int,4> row{a,b,c,d};
            stars.push_back(row);
            if (std::is_sorted(row.begin(),row.end())) sorted.push_back(row);
        }
        for (const auto &a:sorted) for (const auto &b:stars) {
            bool ordered=true;
            for (int i=1;i<4;++i) if (a[i-1]==a[i] && b[i-1]>b[i]) ordered=false;
            if (!ordered) continue;
            Edges g{};
            for (int i=0;i<4;++i) {
                g[0][i+2]=g[i+2][0]=a[i];g[1][i+2]=g[i+2][1]=b[i];
            }
            for (int x=4-2*q;x<q;++x) for (int y=4-2*q;y<q;++y) {
                int z=a[3]+b[3]-x-y;
                std::array<int,6> c{x,y,2*q-a[0]-b[0]-x-y,z,
                                    2*q-a[1]-b[1]-x-z,2*q-a[2]-b[2]-y-z};
                if (*std::max_element(c.begin(),c.end())>=q ||
                    *std::min_element(c.begin(),c.end())<4-2*q) continue;
                int e=0;
                for (int i=2;i<6;++i) for (int j=i+1;j<6;++j)
                    g[i][j]=g[j][i]=c[e++];
                if (convergent(g,q)) register_chart(graph_key(g));
            }
        }
    }
};

constexpr int prime = 1000000007;
int modpow(long long a, int power,int modulus=prime) {
    long long answer=1;
    while (power) {
        if (power&1) answer=answer*a%modulus;
        a=a*a%modulus;power>>=1;
    }
    return int(answer);
}

struct Term {
    int index, scalar;
    Factors ratio;
};
struct Rule {
    Key source;
    unsigned char v,ref,inf;
};
struct ModularTerm { int idx,coef; };
struct Best {
    Rule rule{};
    int support = 100;
    long long score = 0;
};

void shift(Edges &a,int i,int j,int delta) { a[i][j]+=delta;a[j][i]+=delta; }

template<class Add> bool visit_row(const Domain &dom,const Edges &g,const Label &parent,
                                   int v,int ref,int inf,Add add) {
    add(parent,g[v][inf]-g[v][ref]);
    for (int j=0;j<6;++j) {
        if (j==v || j==inf || j==ref || !g[v][j]) continue;
        for (int side=0;side<2;++side) {
            Edges h=g;
            shift(h,v,j,1);shift(h,ref,inf,1);
            shift(h,side ? j : v,ref,-1);shift(h,side ? v : j,inf,-1);
            if (!convergent(h,dom.q)) return false;
            if (!dom.keep_boundaries && has_bubble(h)) continue;
            auto found=dom.labels.find(graph_key(h));
            if (found==dom.labels.end()) throw std::logic_error("missing shifted chart");
            add(found->second,side ? g[v][j] : -g[v][j]);
        }
    }
    return true;
}

bool make_row(const Domain &dom,const Edges &g,const Label &parent,
              int v,int ref,int inf,std::array<Term,9> &terms,int &n) {
    n=0;
    auto add=[&](const Label &l,int w) {
        if(l.index>=0 && w)terms[n++]={l.index,w,l.ratio};
    };
    return visit_row(dom,g,parent,v,ref,inf,add) && n>0;
}

class ModFactors {
    static constexpr int limit=256;
    std::array<std::array<int,2*limit+1>,4> values{};
public:
    int modulus;
    explicit ModFactors(int p=prime):modulus(p) {
        for (int k=0;k<4;++k) {
            int inv=modpow(primes[k],modulus-2,modulus);
            values[k][limit]=1;
            for (int j=1;j<=limit;++j) {
                values[k][limit+j]=long(values[k][limit+j-1])*primes[k]%modulus;
                values[k][limit-j]=long(values[k][limit-j+1])*inv%modulus;
            }
        }
    }
    int operator()(const Factors &f) const {
        long value=1;
        for (int i=0;i<4;++i) {
            if (f[i]<-limit || f[i]>limit) throw std::logic_error("prime exponent bound");
            value=value*values[i][limit+f[i]]%modulus;
        }
        return value;
    }
};

int collect(const std::array<Term,9> &terms,int n,const ModFactors &factor,
            std::array<ModularTerm,9> &values) {
    int size=0;
    int modulus=factor.modulus;
    for (int j=0;j<n;++j) {
        int i=0;
        while (i<size && values[i].idx!=terms[j].index) ++i;
        if (i==size) values[size++]={terms[j].index,0};
        int increment=long(terms[j].scalar)*factor(terms[j].ratio)%modulus;
        int result=values[i].coef+increment;
        if (result<0) result+=modulus;
        else if (result>=modulus) result-=modulus;
        values[i].coef=result;
    }
    return std::remove_if(values.begin(),values.begin()+size,
                          [](auto v){return v.coef==0;})-values.begin();
}

void pivots(const Domain &dom, const char *output) {
    std::vector<Best> best(dom.representatives.size());
    std::array<Term,9> terms;
    std::array<ModularTerm,9> row;
    ModFactors factor;
    long rows=0, seen=0;
    for (const auto &[key,parent]:dom.labels) {
        Edges g=unpack(key);
        for (int v=0;v<6;++v) for (int ref=2;ref<6;++ref) for (int inf=2;inf<6;++inf) {
            if (v==ref || v==inf || ref==inf) continue;
            int n=0;
            if (!make_row(dom,g,parent,v,ref,inf,terms,n)) continue;
            int size=collect(terms,n,factor,row);
            if (!size) continue;
            ++rows;
            int pivot=0;long long score=0;
            for (int i=0;i<size;++i) { pivot=std::max(pivot,row[i].idx);score+=row[i].idx; }
            Best &r=best[pivot];
            if (size<r.support || (size==r.support && score<r.score)) {
                r.rule={key,(unsigned char)v,(unsigned char)ref,(unsigned char)inf};
                r.support=size; r.score=score;
            }
        }
        ++seen;
        if (seen%500000==0) std::cerr << "charts checked " << seen << std::endl;
    }
    long rank=0;
    std::ofstream file(output);
    if (!file) throw std::runtime_error("cannot write pivot cache");
    for (int k=0;k<int(best.size());++k) {
        const auto &p=best[k];
        if (p.support==100) continue;
        ++rank;
        file << k << ' ';
        for (int value:p.rule.source) file << value << ' ';
        file << int(p.rule.v) << ' ' << int(p.rule.ref) << ' ' << int(p.rule.inf) << '\n';
    }
    std::cout << "nonzero rows " << rows << " direct triangular pivots " << rank
              << " remaining " << best.size()-rank << std::endl;
}

struct ExactRow {
    std::array<Term,9> terms;
    int n=0;
};
struct FastRow {
    std::array<ModularTerm,9> terms;
    int n=0;
};
FastRow convert(const ExactRow &r, const ModFactors &f) {
    FastRow m;
    m.n=collect(r.terms,r.n,f,m.terms);
    return m;
}
Rule read_rule(std::istream &in) {
    Rule r;
    int x=0,v,ref,inf;
    for(auto &c:r.source) {in>>x;c=x;}
    in>>v>>ref>>inf; r.v=v;r.ref=ref;r.inf=inf;
    if(!in) throw std::runtime_error("bad rule file");
    return r;
}
void write_rule(std::ostream &out,const Rule &r) {
    for(int a:r.source) out << a << ' ';
    out << int(r.v) << ' ' << int(r.ref) << ' ' << int(r.inf) << '\n';
}
ExactRow generate(const Domain &d,const Rule &r) {
    ExactRow result;
    if (!make_row(d,unpack(r.source),d.labels.at(r.source),r.v,r.ref,r.inf,result.terms,result.n))
        throw std::runtime_error("invalid rule");
    return result;
}

struct SeedSystem {
    std::vector<ExactRow> rows;
    std::vector<int> free;
    std::vector<Rule> candidates;
    ExactRow anchor;
    mpq_class anchor_value;
    SeedSystem(const Domain &d,const char *filename):rows(d.representatives.size()) {
        std::ifstream input(filename);
        if (!input) throw std::runtime_error("missing pivot file");
        int k;
        while(input>>k) rows.at(k)=generate(d,read_rule(input));
        for(int i=0;i<int(rows.size());++i) if(!rows[i].n) free.push_back(i);
        Powers a = d.q==6? Powers{1,1,2,1,3,4,4,4,4}:Powers{4,4,5,9,6,4,5,5,6};
        auto l=d.labels.at(graph_key(complete(a,d.q)));
        if(l.index<0)throw std::logic_error("zero anchor");
        anchor.terms[0]={l.index,1,l.ratio};anchor.n=1;
        anchor_value = mpq_class(d.q==6?"-1087/2268":"1/41191991672832");
    }

    void sample(const Domain &d,int additional,const std::vector<int> *remaining=nullptr) {
        // Include equations at every remaining variable, then a reproducible
        // spread of other completed charts. Save only rows selected for rank.
        std::set<Key> sources;
        for(int k:(remaining?*remaining:free)) sources.insert(d.representatives[k]);
        std::mt19937 rng(1909);
        std::vector<Key> extra;
        int seen=0;
        for(const auto &[g,l]:d.labels) {
            ++seen;
            if(int(extra.size())<additional) extra.push_back(g);
            else {int slot=rng()%seen;if(slot<additional)extra[slot]=g;}
        }
        sources.insert(extra.begin(),extra.end());
        for(const Key &g:sources) {
            Edges m=unpack(g);
            for(int v=0;v<6;++v) for(int b=2;b<6;++b) for(int i=2;i<6;++i) {
                if(v==b||v==i||b==i)continue;
                Rule r{g,(unsigned char)v,(unsigned char)b,(unsigned char)i};
                ExactRow entry;
                if(make_row(d,m,d.labels.at(g),v,b,i,entry.terms,entry.n))
                    candidates.push_back(r);
            }
        }
    }
};

std::vector<FastRow> normalized(const SeedSystem &s,const ModFactors &f) {
    std::vector<FastRow> out(s.rows.size());
    for(int k=0;k<int(out.size());++k) {
        if(!s.rows[k].n)continue;
        FastRow m=convert(s.rows[k],f);
        int top=0;
        for(int i=0;i<m.n;++i) {
            if(m.terms[i].idx>k)throw std::logic_error("nontriangular pivot");
            if(m.terms[i].idx==k)top=m.terms[i].coef;
            else out[k].terms[out[k].n++]=m.terms[i];
        }
        if(!top)throw std::logic_error("singular pivot prime");
        int inv=f.modulus-modpow(top,f.modulus-2,f.modulus);
        for(int j=0;j<out[k].n;++j)out[k].terms[j].coef=long(out[k].terms[j].coef)*inv%f.modulus;
    }
    return out;
}


using Sparse = std::vector<ModularTerm>;
void add_sparse(Sparse &r,int index,int value,int p) {
    if(!value)return;
    auto it=std::lower_bound(r.begin(),r.end(),index,
                            [](const ModularTerm &a,int i){return a.idx<i;});
    if(it==r.end() || it->idx!=index) {r.insert(it,{index,value});return;}
    it->coef+=value;
    if(it->coef>=p)it->coef-=p;
    if(!it->coef)r.erase(it);
}

struct Added {int pivot; Rule rule;};

struct SparseSystem {
    std::vector<FastRow> primary;
    std::vector<Sparse> secondary;
    std::vector<int> free;
    const SeedSystem &seed;
    ModFactors field;
    SparseSystem(const SeedSystem &s,int modulus):primary(normalized(s,ModFactors(modulus))),
        secondary(s.rows.size()),free(s.free),seed(s),field(modulus) {}
    bool has_rule(int k) const {return seed.rows[k].n || !secondary[k].empty();}

    // Each substitution replaces the largest active variable by smaller ones.
    Sparse reduce(const ExactRow &eq,int support_limit,int step_limit,int target=-1) const {
        auto m=convert(eq,field);
        Sparse r(m.terms.begin(),m.terms.begin()+m.n);
        std::sort(r.begin(),r.end(),[](auto a,auto b){return a.idx<b.idx;});
        int steps=0;
        const int p=field.modulus;
        while(!r.empty()) {
            const int k=r.back().idx,c=r.back().coef;
            if(!has_rule(k) || k==target) return r;
            if(++steps>step_limit || int(r.size())>support_limit)return {};
            r.pop_back();
            auto add=[&](const ModularTerm &t){add_sparse(r,t.idx,long(t.coef)*c%p,p);};
            if(seed.rows[k].n)for(int i=0;i<primary[k].n;++i)add(primary[k].terms[i]);
            else for(const auto &t:secondary[k])add(t);
        }
        return {};
    }
    void install(Sparse r) {
        int k=r.back().idx,p=field.modulus;
        int inv=p-modpow(r.back().coef,p-2,p);
        r.pop_back();
        if(r.empty())throw std::logic_error("unexpected zero period in full sector");
        for(auto &t:r)t.coef=long(t.coef)*inv%p;
        secondary[k]=std::move(r);
    }
    void refresh() {
        free.clear();for(int k:seed.free)if(!has_rule(k))free.push_back(k);
    }
};

void sparse_pass(const Domain &d,const SeedSystem &s,const char *save,int extra) {
    SparseSystem sys(s,prime);
    std::ofstream file(save);
    if(!file)throw std::runtime_error("cannot store sparse rows");
    std::vector<Key> pool;
    pool.reserve(d.labels.size());
    for(const auto &[k,l]:d.labels)pool.push_back(k);
    std::mt19937 rng(1909);
    std::shuffle(pool.begin(),pool.end(),rng);
    int rank=0;
    for(int pass=0;pass<3;++pass) {
        auto todo=sys.free;
        std::shuffle(todo.begin(),todo.end(),rng);
        std::vector<Key> candidates;
        for(int k:todo)candidates.push_back(d.representatives[k]);
        int cap=std::min<int>((pass+1)*extra,pool.size());
        candidates.insert(candidates.end(),pool.begin()+std::min<int>(pass*extra,pool.size()),pool.begin()+cap);
        int max_terms=std::array<int,3>{16,32,72}[pass];
        int max_steps=std::array<int,3>{14,40,120}[pass];
        for(const Key &g:candidates) {
            auto parent=d.labels.at(g);
            Edges m=unpack(g);
            for(int v=0;v<6;++v)for(int ref=2;ref<6;++ref)for(int inf=2;inf<6;++inf) {
                if(v==ref || v==inf || ref==inf)continue;
                ExactRow eq;
                if(!make_row(d,m,parent,v,ref,inf,eq.terms,eq.n))continue;
                Sparse r=sys.reduce(eq,max_terms,max_steps);
                if(r.empty())continue;
                int k=r.back().idx;
                if(sys.has_rule(k))throw std::logic_error("duplicate pivot");
                file << k << ' ';
                write_rule(file,{g,(unsigned char)v,(unsigned char)ref,(unsigned char)inf});
                sys.install(r);++rank;
            }
        }
        sys.refresh();
        file.flush();
        std::cout << "sparse pass " << pass+1 << " added " << rank
                  << " remaining " << sys.free.size() << std::endl;
        if(sys.free.size()<=3)break;
    }
}
// Evaluate each retained equation on all free basis vectors in small blocks.
std::vector<int> residuals(const SeedSystem &s, const std::vector<FastRow> &tri,
                           const std::vector<FastRow> &checks,int modulus,
                           const SparseSystem *extended=nullptr) {
    constexpr int width=24;
    const auto &free=extended?extended->free:s.free;
    int nf=free.size();
    std::vector<int> output(std::size_t(checks.size())*nf);
    std::vector<int> columns(std::size_t(tri.size())*width);
    auto evaluate=[&](const ModularTerm *terms,int count,int *dest,int k) {
        std::array<std::uint64_t,width> sums{};
        for(int j=0;j<count;++j) {
            if(j && j%8==0) for(int a=0;a<k;++a)sums[a]%=modulus;
            const auto [idx,c]=terms[j];
            const int *p=columns.data()+std::size_t(idx)*width;
            for(int a=0;a<k;++a)sums[a]+=std::uint64_t(c)*p[a];
        }
        for(int a=0;a<k;++a)dest[a]=sums[a]%modulus;
    };
    for(int start=0;start<nf;start+=width) {
        int k=std::min(nf-start,width), which=0;
        for(int i=0;i<int(tri.size());++i) {
            int *dest=columns.data()+std::size_t(i)*width;
            if(which<nf && free[which]==i) {
                std::fill(dest,dest+k,0);
                if(start<=which && which<start+k)dest[which-start]=1;
                ++which;
            } else if(extended && !extended->secondary[i].empty()){
                const auto &r=extended->secondary[i];evaluate(r.data(),r.size(),dest,k);
            } else evaluate(tri[i].terms.data(),tri[i].n,dest,k);
        }
        for(int i=0;i<int(checks.size());++i)
            evaluate(checks[i].terms.data(),checks[i].n,output.data()+std::size_t(i)*nf+start,k);
    }
    return output;
}

struct Echelon {
    std::vector<std::vector<int>> basis;
    int modulus;
    Echelon(int n,int p):basis(n),modulus(p) {}
    bool add(std::vector<int> row) {
        for(int i=0;i<int(basis.size());++i) {
            int c=row[i];if(!c)continue;
            if(basis[i].empty()) {
                int inv=modpow(c,modulus-2,modulus);
                for(int j=i;j<int(row.size());++j)row[j]=long(row[j])*inv%modulus;
                basis[i]=std::move(row); return true;
            }
            for(int j=i;j<int(row.size());++j)
                row[j]=(row[j]+long(modulus-c)*basis[i][j])%modulus;
        }
        if(row.size()>basis.size() && row.back())throw std::logic_error("inconsistent equations");
        return false;
    }
};

void choose(const Domain &d,SeedSystem &s,const char *filename,int extra) {
    s.sample(d,extra);
    const ModFactors f;
    auto tri=normalized(s,f);
    std::vector<FastRow> checks;
    for(const Rule &r:s.candidates)checks.push_back(convert(generate(d,r),f));
    std::vector<int> mat=residuals(s,tri,checks,prime);
    int n=s.free.size();
    Echelon solve(n,prime);
    std::ofstream out(filename);
    if(!out)throw std::runtime_error("cannot write independent rows");
    int rank=0;
    for(int i=0;i<int(checks.size());++i) {
        const int *row=mat.data()+std::size_t(i)*n;
        if(solve.add(std::vector<int>(row,row+n))) {write_rule(out,s.candidates[i]);++rank;}
        if(rank==n-1)break;
    }
    std::cout << "candidate rows " << checks.size() << " residual rank " << rank
              << " of " << n << std::endl;
}


struct BoundarySystem {
    std::vector<std::pair<int,ExactRow>> added;
    std::array<ExactRow,3> anchors;
    std::array<std::array<mpq_class,3>,3> anchor_values{};
    BoundarySystem(const Domain &d,const SeedSystem &s,const char *file) {
        std::ifstream f(file);if(!f)throw std::runtime_error("cannot open sparse pivots");
        int k;while(f>>k)added.emplace_back(k,generate(d,read_rule(f)));
        const int q=d.q,b=2*q/3;
        Edges triple{},kite{};
        for(int i=0;i<3;++i)for(int j=3;j<6;++j)shift(triple,i,j,b);
        for(int i=0;i<4;++i)for(int j=i+1;j<4;++j)shift(kite,i,j,q/2);
        shift(kite,0,1,q-2*b);shift(kite,0,4,b);shift(kite,1,4,b);
        for(int i=0;i<5;++i) {
            int left=2*q;for(int j=0;j<5;++j)left-=kite[i][j];
            shift(kite,i,5,left);
        }
        for(int i=0;i<2;++i) {
            const auto &l=d.labels.at(graph_key(i?kite:triple));
            anchors[i].n=1;anchors[i].terms[0]={l.index,1,l.ratio};
        }
        anchors[2]=s.anchor;
        if(q==6) {
            anchor_values[0][0]=mpq_class("1/10077696");
            anchor_values[1][0]=mpq_class("-173/2488320");
            anchor_values[1][1]=mpq_class("1/17280");
            anchor_values[2][0]=mpq_class("83175731/705438720");
            anchor_values[2][1]=mpq_class("14713/46656");
        } else {
            anchor_values[0][0]=mpq_class("1/208215748530929664000000000");
            anchor_values[1][0]=mpq_class("-3541028833/9363576310742384640000000000");
            anchor_values[1][1]=mpq_class("2273/7224981721251840000000");
            anchor_values[2][0]=mpq_class("-531372573405073/88991429257295623618560000000");
            anchor_values[2][1]=mpq_class("-80588077/5044865604155080704000");
        }
        anchor_values[2][2]=s.anchor_value;
    }
    void restore(SparseSystem &sys) const {
        for(const auto &[k,raw]:added) {
            auto r=sys.reduce(raw,10000,10000,k);
            if(r.empty()||r.back().idx!=k)throw std::logic_error("incompatible sparse pivot prime");
            sys.install(r);
        }
        sys.refresh();
    }
};

// Null coordinates on the last triangular variables, extended in chart order.
std::vector<int> null_vector(const Echelon &solver,const std::vector<int> &parameters) {
    const int n=solver.basis.size(),p=solver.modulus;
    std::vector<int> result(n);
    int slot=parameters.size();
    for(int i=n-1;i>=0;--i) {
        const auto &r=solver.basis[i];
        if(r.empty()) {result[i]=parameters.at(--slot);continue;}
        long sum=0;
        for(int j=i+1;j<n;++j)sum=(sum+long(r[j])*result[j])%p;
        result[i]=sum?p-sum:0;
    }
    if(slot)throw std::logic_error("null parameter count");
    return result;
}
template<class T> void for_rhs(const SparseSystem &sys,int i,T add) {
    if(sys.seed.rows[i].n) {
        const auto &r=sys.primary[i];
        for(int j=0;j<r.n;++j)add(r.terms[j]);
    } else for(const auto &t:sys.secondary[i])add(t);
}
std::vector<int> extend(const SparseSystem &sys,const std::vector<int> &base) {
    int slot=0;
    std::vector<int> out(sys.primary.size());
    for(int i=0;i<int(out.size());++i) {
        if(slot<int(sys.free.size()) && sys.free[slot]==i)out[i]=base[slot++];
        else {
            long sum=0;
            for_rhs(sys,i,[&](const ModularTerm &t){sum=(sum+long(t.coef)*out[t.idx])%sys.field.modulus;});
            out[i]=sum;
        }
    }
    return out;
}

int complete_rank(const Domain &dom,const SparseSystem &sys,Echelon &first,
                  std::ofstream &file,int rank) {
    int n=sys.free.size(),k=n-rank;
    if(k==3)return rank;
    std::vector<std::vector<int>> columns;
    for(int j=0;j<k;++j) {
        std::vector<int> param(k);param[j]=1;
        columns.push_back(extend(sys,null_vector(first,param)));
    }
    Echelon extra(k,prime);
    std::mt19937 rng(220901);
    auto probe=[&]() {
        std::vector<int> param(n-rank);
        for(int &v:param)v=1+rng()%(prime-1);
        return extend(sys,null_vector(first,null_vector(extra,param)));
    };
    auto trial=probe();
    long tested=0;
    std::vector<Key> pool;
    pool.reserve(dom.labels.size());
    for(const auto &[key,label]:dom.labels)pool.push_back(key);
    std::shuffle(pool.begin(),pool.end(),rng);
    for(const auto &key:pool) {
        const auto &parent=dom.labels.at(key);const auto g=unpack(key);
        for(int v=0;v<6;++v)for(int ref=2;ref<6;++ref)for(int inf=2;inf<6;++inf) {
            if(v==ref || v==inf || ref==inf)continue;
            ExactRow raw;
            if(!make_row(dom,g,parent,v,ref,inf,raw.terms,raw.n))continue;
            ++tested;
            auto row=convert(raw,sys.field);
            auto evaluate=[&](const std::vector<int> &vec){
                std::uint64_t sum=0;
                for(int j=0;j<row.n;++j)sum+=std::uint64_t(row.terms[j].coef)*vec[row.terms[j].idx];
                return int(sum%prime);
            };
            if(!evaluate(trial))continue;
            std::vector<int> r(k);
            for(int j=0;j<k;++j)r[j]=evaluate(columns[j]);
            if(!extra.add(r))throw std::logic_error("null witness failed");
            write_rule(file,{key,(unsigned char)v,(unsigned char)ref,(unsigned char)inf});
            file.flush();
            ++rank;
            if(rank%10==0 || rank==n-3)
                std::cout<<"streamed equations "<<tested<<" residual rank "<<rank<<" of "<<n<<std::endl;
            if(rank==n-3)return rank;
            trial=probe();
        }
    }
    return rank;
}

void choose_full(const Domain &d,SeedSystem &s,const BoundarySystem &b,
                 const char *target,int extra) {
    SparseSystem sys(s,prime);b.restore(sys);
    s.sample(d,extra,&sys.free);
    int n=sys.free.size();std::ofstream file(target);
    if(!file)throw std::runtime_error("cannot write full rank rows");
    if(n==3) {std::cout<<"three sparse coordinates"<<std::endl;return;}
    std::vector<FastRow> checks;
    for(const auto &r:s.candidates)checks.push_back(convert(generate(d,r),sys.field));
    auto mat=residuals(s,sys.primary,checks,prime,&sys);
    Echelon solver(n,prime);int rank=0;
    for(int i=0;i<int(checks.size());++i) {
        const int *p=mat.data()+std::size_t(i)*n;
        if(solver.add(std::vector<int>(p,p+n))) {write_rule(file,s.candidates[i]);++rank;}
        if(rank==n-3)break;
    }
    std::cout<<"full candidate equations "<<checks.size()<<" rank "<<rank<<" of "<<n<<std::endl;
    file.flush();
    std::vector<int>().swap(mat);std::vector<FastRow>().swap(checks);
    std::vector<Rule>().swap(s.candidates);
    if(complete_rank(d,sys,solver,file,rank)!=n-3)throw std::logic_error("insufficient full rank");
}
std::vector<int> modular_values(const SeedSystem &s,const std::vector<ExactRow> &conditions,
                                int modulus) {
    ModFactors f(modulus);
    auto tri=normalized(s,f);
    std::vector<FastRow> checks;
    for(const auto &r:conditions)checks.push_back(convert(r,f));
    checks.push_back(convert(s.anchor,f));
    int nf=s.free.size(),n=s.rows.size();
    auto matrix=residuals(s,tri,checks,modulus);
    Echelon pivot(nf,modulus);
    if(int(checks.size())!=nf)throw std::logic_error("incomplete final system");
    for(int i=0;i<nf;++i) {
        const int *a=matrix.data()+std::size_t(i)*nf;
        std::vector<int> r(a,a+nf);
        int rhs=0;
        if(i+1==nf) {
            int num=mpz_fdiv_ui(s.anchor_value.get_num_mpz_t(),modulus);
            int den=mpz_fdiv_ui(s.anchor_value.get_den_mpz_t(),modulus);
            if(!den)throw std::logic_error("singular anchor denominator");
            rhs=long(num)*modpow(den,modulus-2,modulus)%modulus;
        }
        r.push_back(rhs);
        if(!pivot.add(r))throw std::logic_error("dependent final equation at this prime");
    }
    std::vector<int> base(nf),values(n);
    for(int i=nf-1;i>=0;--i) {
        const auto &r=pivot.basis[i];
        long long d=r.back();
        for(int j=i+1;j<nf;++j)d=(d-long(r[j])*base[j])%modulus;
        base[i]=int((d+modulus)%modulus);
    }
    int a=0;
    for(int i=0;i<n;++i) {
        if(a<nf && s.free[a]==i) {values[i]=base[a++];continue;}
        std::uint64_t d=0;
        for(int j=0;j<tri[i].n;++j) {
            const auto [idx,c]=tri[i].terms[j];
            d+=std::uint64_t(c)*values[idx];
        }
        values[i]=d%modulus;
    }
    return values;
}

bool reconstruct(mpq_class &out,const mpz_class &residue,const mpz_class &mod,
                 const mpz_class &bound) {
    if(residue==0) {out=0;return true;}
    mpz_class a=mod,b=residue, t=0,u=1, quotient,rem,temp;
    while(b>bound) {
        mpz_tdiv_qr(quotient.get_mpz_t(),rem.get_mpz_t(),a.get_mpz_t(),b.get_mpz_t());
        a.swap(b);b.swap(rem);
        temp=t-quotient*u; t.swap(u);u.swap(temp);
    }
    if(!u || abs(u)>bound)return false;
    mpz_gcd(temp.get_mpz_t(),b.get_mpz_t(),u.get_mpz_t());
    if(temp!=1)return false;
    out=mpq_class(b,u);out.canonicalize();return true;
}

class ExactFactors {
    std::map<Factors,mpq_class> cache;
public:
    const mpq_class &get(const Factors &f) {
        auto [it,fresh]=cache.emplace(f,0);
        if(fresh) {
            mpz_class num=1,den=1,tmp;
            for(int j=0;j<4;++j) {
                mpz_ui_pow_ui(tmp.get_mpz_t(),primes[j],std::abs(f[j]));
                (f[j]<0?den:num)*=tmp;
            }
            it->second=mpq_class(num,den);
        }
        return it->second;
    }
    mpq_class evaluate(const ExactRow &row,const std::vector<mpq_class> &values) {
        mpq_class sum;
        for(int j=0;j<row.n;++j) {
            const auto &t=row.terms[j];
            sum+=get(t.ratio)*t.scalar*values[t.index];
        }
        return sum;
    }
};

std::vector<mpq_class> solve_exact(const Domain &d,SeedSystem &s,const char *last) {
    std::ifstream src(last);
    if(!src)throw std::runtime_error("cannot open final equations");
    std::vector<ExactRow> conditions;
    while(src>>std::ws && src.peek()!=EOF)conditions.push_back(generate(d,read_rule(src)));
    mpz_class modulus=1;
    std::vector<mpz_class> residues(s.rows.size());
    int p=prime;
    for(int pass=1;pass<=60;++pass) {
        auto values=modular_values(s,conditions,p);
        int inv=modpow(mpz_fdiv_ui(modulus.get_mpz_t(),p),p-2,p);
        for(std::size_t i=0;i<residues.size();++i) {
            int r=mpz_fdiv_ui(residues[i].get_mpz_t(),p);
            auto t=long((values[i]-r+p)%p)*inv%p;
            mpz_addmul_ui(residues[i].get_mpz_t(),modulus.get_mpz_t(),t);
        }
        modulus*=p;
        std::vector<mpq_class> exact(residues.size());
        mpz_class bound=modulus/2;
        mpz_sqrt(bound.get_mpz_t(),bound.get_mpz_t());
        bool complete=true;
        for(std::size_t i=0;i<residues.size();++i)
            if(!reconstruct(exact[i],residues[i],modulus,bound)) {complete=false;break;}
        if(complete) {
            ExactFactors factors;
            complete=factors.evaluate(s.anchor,exact)==s.anchor_value;
            for(const auto &r:s.rows) if(complete && r.n && factors.evaluate(r,exact)!=0)complete=false;
            for(const auto &r:conditions) if(complete && factors.evaluate(r,exact)!=0)complete=false;
            if(complete) {
                int nonzero=0;unsigned long numbits=0,denbits=0;
                for(const auto &v:exact) {
                    nonzero+=v!=0;
                    numbits=std::max(numbits,mpz_sizeinbase(v.get_num_mpz_t(),2));
                    denbits=std::max(denbits,mpz_sizeinbase(v.get_den_mpz_t(),2));
                }
                std::cout << "exact pass " << pass << " nonzero " << nonzero
                          << " max numerator and denominator bits " << numbits << ' ' << denbits
                          << std::endl;
                return exact;
            }
        }
        std::cout << "CRT pass " << pass << " (" << mpz_sizeinbase(modulus.get_mpz_t(),2)
                  << " bits)" << std::endl;
        do {
            p-=2;
            bool yes=true;
            for(int i=3;long(i)*i<=p;i+=2) if(p%i==0) {yes=false;break;}
            if(yes)break;
        } while(true);
    }
    throw std::runtime_error("CRT budget exceeded");
}

void physical_weights(const Domain &d,const std::vector<mpq_class> &values) {
    ExactFactors factors;
    std::array<mpq_class,2> sums{};
    long rows=0,nonzero=0;
    Powers a;
    std::string k,j;
    while(std::cin>>a[0]) {
        for(int i=1;i<9;++i)if(!(std::cin>>a[i]))throw std::runtime_error("truncated powers");
        if(!(std::cin>>k>>j))throw std::runtime_error("truncated weights");
        ++rows;
        Edges g=complete(a,d.q);
        auto l=d.labels.at(graph_key(g));
        if(l.index<0 || values[l.index]==0)continue;
        ++nonzero;
        mpq_class coeff=values[l.index]*factors.get(l.ratio);
        sums[0]+=mpq_class(k)*coeff;sums[1]+=mpq_class(j)*coeff;
    }
    std::cout << "weight sums " << rows << ' ' << nonzero << ' ' << sums[0] << ' '
              << sums[1] << std::endl;
}

std::vector<mpq_class> read_table(const Domain &domain,const char *path) {
    std::ifstream f(path);
    if(!f)throw std::runtime_error("cannot open coefficients");
    std::vector<mpq_class> values;
    for(const Key &rep:domain.representatives) {
        for(int a:rep) {
            int b;if(!(f>>b) || a!=b)throw std::logic_error("mismatched coefficient chart");
        }
        mpq_class value;
        if(!(f>>value))throw std::runtime_error("missing coefficient");
        values.push_back(value);
    }
    std::string extra;if(f>>extra)throw std::logic_error("unexpected coefficient");
    return values;
}

void check_all(const Domain &d,const std::vector<mpq_class> &exact) {
    const int p=1000000009;
    ModFactors factors(p);
    std::vector<int> values;
    for(const auto &v:exact) {
        auto den=mpz_fdiv_ui(v.get_den_mpz_t(),p);
        if(!den)throw std::logic_error("validation prime divides denominator");
        values.push_back(mpz_fdiv_ui(v.get_num_mpz_t(),p)*modpow(den,p-2,p)%p);
    }
    std::array<Term,9> terms;long checks=0;
    for(const auto &[key,parent]:d.labels) {
        Edges g=unpack(key);
        for(int v=0;v<6;++v) for(int ref=2;ref<6;++ref) for(int inf=2;inf<6;++inf) {
            if(v==ref||v==inf||ref==inf)continue;
            int n=0;
            if(!make_row(d,g,parent,v,ref,inf,terms,n))continue;
            long long value=0;
            for(int j=0;j<n;++j) {
                auto &t=terms[j];
                value+=long(factors(t.ratio))*values[t.index]%p*t.scalar;
            }
            if(value%p)throw std::logic_error("independent IBP validation failed");
            ++checks;
        }
    }
    std::cout << "independent modular checks " << checks << " prime " << p << std::endl;
}

int reduce_rational(const mpq_class &v,int p) {
    int den=mpz_fdiv_ui(v.get_den_mpz_t(),p);
    if(!den)throw std::logic_error("prime divides anchor denominator");
    return mpz_fdiv_ui(v.get_num_mpz_t(),p)*modpow(den,p-2,p)%p;
}

std::vector<std::array<int,3>> seed_values(const SeedSystem &s,const BoundarySystem &b,
                                            const std::vector<ExactRow> &conditions,int p) {
    SparseSystem sys(s,p);b.restore(sys);
    const int n=sys.free.size();
    if(int(conditions.size())+3!=n)throw std::logic_error("wrong final rank");
    std::vector<FastRow> checks;
    for(const auto &r:conditions)checks.push_back(convert(r,sys.field));
    for(const auto &r:b.anchors)checks.push_back(convert(r,sys.field));
    auto mat=residuals(s,sys.primary,checks,p,&sys);
    Echelon solver(n,p);
    for(int i=0;i<n;++i) {
        auto ptr=mat.data()+std::size_t(i)*n;
        std::vector<int> row(ptr,ptr+n);
        for(int j=0;j<3;++j)
            row.push_back(i<n-3?0:reduce_rational(b.anchor_values[i-(n-3)][j],p));
        if(!solver.add(row))throw std::logic_error("dependent anchor or singular modular minor");
    }
    std::vector<std::array<int,3>> out(s.free.size());
    for(int ch=0;ch<3;++ch) {
        std::vector<int> base(n);
        for(int i=n-1;i>=0;--i) {
            const auto &r=solver.basis[i];
            long sum=r[n+ch];
            for(int j=i+1;j<n;++j)sum=(sum-long(r[j])*base[j])%p;
            base[i]=(sum+p)%p;
        }
        auto values=extend(sys,base);
        for(int k=0;k<int(s.free.size());++k)out[k][ch]=values[s.free[k]];
    }
    return out;
}

void exact_channel(const Domain &dom,const SeedSystem &s,const BoundarySystem &b,
                   const std::vector<ExactRow> &last,
                   const std::vector<std::array<mpq_class,3>> &free,int ch,const std::string &path) {
    ExactFactors factors;
    std::vector<mpq_class> values(s.rows.size());
    int slot=0;
    for(int i=0;i<int(values.size());++i) {
        const auto &r=s.rows[i];
        if(!r.n) {values[i]=free.at(slot++)[ch];continue;}
        mpq_class top,sum;
        std::map<int,mpq_class> upper;
        for(int j=0;j<r.n;++j) {
            const auto &t=r.terms[j];
            const mpq_class c=t.scalar*factors.get(t.ratio);
            if(t.index==i)top+=c;
            else if(t.index<i)sum+=c*values[t.index];
            else upper[t.index]+=c;
        }
        for(const auto &[idx,value]:upper)if(value!=0)throw std::logic_error("exact row out of order");
        if(top==0)throw std::logic_error("zero exact leading coefficient");
        values[i]=-sum/top;
    }
    for(const auto &r:s.rows)if(r.n && factors.evaluate(r,values)!=0)
        throw std::logic_error("exact direct pivot check");
    for(const auto &[k,r]:b.added)if(factors.evaluate(r,values)!=0)
        throw std::logic_error("exact sparse pivot check");
    for(const auto &r:last)if(factors.evaluate(r,values)!=0)
        throw std::logic_error("exact final pivot check");
    for(int a=0;a<3;++a)if(factors.evaluate(b.anchors[a],values)!=b.anchor_values[a][ch])
        throw std::logic_error("exact normalization");
    std::ofstream out(path);
    if(!out)throw std::runtime_error("cannot create full period table");
    for(int i=0;i<int(values.size());++i) {
        for(int j:dom.representatives[i])out<<j<<' ';
        out<<values[i]<<'\n';
    }
    if(!out)throw std::runtime_error("cannot store full period table");
    std::cout<<"verified full channel "<<ch<<" classes "<<values.size()<<std::endl;
}

// Reconstruct the free coordinates of the direct triangular system only. All
// other coordinates are computed over Q by direct substitution.
void full_solve(const Domain &d,SeedSystem &s,const BoundarySystem &b,
                const char *final,const char *prefix) {
    std::ifstream src(final);
    if(!src)throw std::runtime_error("missing final pivot file");
    std::vector<ExactRow> last;
    while(src>>std::ws && src.peek()!=EOF)last.push_back(generate(d,read_rule(src)));
    std::vector<std::array<mpz_class,3>> residues(s.free.size());
    std::vector<std::array<mpq_class,3>> exact(s.free.size());
    mpz_class mod=1;
    bool candidate=false;
    for(int pass=1,p=prime;pass<=60;++pass) {
        const auto vals=seed_values(s,b,last,p);
        if(candidate) {
            for(int i=0;i<int(vals.size()) && candidate;++i)for(int j=0;j<3;++j)
                if(reduce_rational(exact[i][j],p)!=vals[i][j]) {candidate=false;break;}
            if(candidate) {
                for(int ch=0;ch<3;++ch)exact_channel(d,s,b,last,exact,ch,std::string(prefix)+std::to_string(ch));
                return;
            }
        }
        int inv=modpow(mpz_fdiv_ui(mod.get_mpz_t(),p),p-2,p);
        for(int i=0;i<int(vals.size());++i)for(int j=0;j<3;++j) {
            auto &z=residues[i][j];
            int old=mpz_fdiv_ui(z.get_mpz_t(),p);
            auto t=long((vals[i][j]-old+p)%p)*inv%p;
            mpz_addmul_ui(z.get_mpz_t(),mod.get_mpz_t(),t);
        }
        mod*=p;
        mpz_class bound=mod/2;mpz_sqrt(bound.get_mpz_t(),bound.get_mpz_t());
        candidate=true;
        for(int i=0;i<int(vals.size()) && candidate;++i)for(int j=0;j<3;++j)
            if(!reconstruct(exact[i][j],residues[i][j],mod,bound)){candidate=false;break;}
        std::cout<<"full CRT pass "<<pass<<" bits "<<mpz_sizeinbase(mod.get_mpz_t(),2)
                 <<" candidate "<<candidate<<std::endl;
        for(bool next=false;!next;) {
            p-=2;next=true;
            for(int i=3;long(i)*i<=p;i+=2)if(p%i==0){next=false;break;}
        }
    }
    throw std::logic_error("full CRT exhaustion");
}

void compare_cuts(const Domain &d,const std::vector<mpq_class> &values) {
    ExactFactors factors;
    Powers powers;
    int n=0;
    while(std::cin>>powers[0]) {
        for(int i=1;i<9;++i)if(!(std::cin>>powers[i]))throw std::runtime_error("short cut");
        mpq_class expected;
        if(!(std::cin>>expected))throw std::runtime_error("missing expected cut");
        const auto &l=d.labels.at(graph_key(complete(powers,d.q)));
        if(l.index<0 || values[l.index]*factors.get(l.ratio)!=expected)
            throw std::logic_error("independent cut disagrees on row "+std::to_string(n));
        ++n;
    }
    std::cout<<"independent cuts "<<n<<std::endl;
}

void group_weights(const Domain &d,const char *output) {
    ExactFactors factors;
    std::vector<std::array<mpq_class,2>> bank(d.representatives.size());
    long rows=0,groups=0;
    Powers a;
    while(std::cin>>a[0]) {
        for(int i=1;i<9;++i)if(!(std::cin>>a[i]))throw std::runtime_error("short prism powers");
        const auto &l=d.labels.at(graph_key(complete(a,d.q)));
        if(l.index<0)throw std::logic_error("missing full chart");
        auto &w=bank[l.index];
        mpq_class k,j;
        if(!(std::cin>>k>>j))throw std::runtime_error("missing prism weight");
        const auto &scale=factors.get(l.ratio);
        w[0]+=k*scale;w[1]+=j*scale;++rows;
    }
    std::ofstream f(output);if(!f)throw std::runtime_error("cannot store grouped bank");
    for(int i=0;i<int(bank.size());++i)if(bank[i][0]!=0 || bank[i][1]!=0) {
        f<<i<<' '<<bank[i][0]<<' '<<bank[i][1]<<'\n';++groups;
    }
    std::cout<<"physical rows "<<rows<<" reflected groups "<<groups<<std::endl;
}

void group_sum(const std::vector<mpq_class> &values) {
    std::array<mpq_class,2>sums{};long n=0;int index;
    while(std::cin>>index) {
        mpq_class k,j;if(!(std::cin>>k>>j))throw std::runtime_error("missing grouped weights");
        sums[0]+=values.at(index)*k;sums[1]+=values.at(index)*j;++n;
    }
    std::cout<<"grouped sum "<<n<<' '<<sums[0]<<' '<<sums[1]<<std::endl;
}

// Encode the three rational coefficients as one integer per completed chart.
// Each radix exceeds the largest absolute residual of its coefficient, so
// vanishing of the packed integer forces all three rational residuals to vanish.
void check_integer_rows(Domain &d,const std::string &prefix) {
    ExactFactors factors;
    std::map<Factors,mpz_class> multipliers;
    for(const auto &[key,l]:d.labels)multipliers.try_emplace(l.ratio,0);
    mpz_class common=1,largest=0;
    for(const auto &[ratio,m]:multipliers)
        mpz_lcm(common.get_mpz_t(),common.get_mpz_t(),factors.get(ratio).get_den_mpz_t());
    for(auto &[ratio,m]:multipliers) {
        const auto &r=factors.get(ratio);
        m=(common/r.get_den())*r.get_num();
        if(m>largest)largest=m;
    }
    std::vector<mpz_class> packed(d.representatives.size());
    for(int ch=0;ch<3;++ch) {
        auto values=read_table(d,(prefix+std::to_string(ch)).c_str());
        mpz_class den=1,maxabs=0;
        for(const auto &v:values)mpz_lcm(den.get_mpz_t(),den.get_mpz_t(),v.get_den_mpz_t());
        std::vector<mpz_class> integer(values.size());
        for(int i=0;i<int(values.size());++i) {
            integer[i]=(den/values[i].get_den())*values[i].get_num();
            const auto magnitude=abs(integer[i]);
            if(magnitude>maxabs)maxabs=magnitude;
        }
        mpz_class bound=27*d.q*largest*maxabs;
        const int shift=mpz_sizeinbase(bound.get_mpz_t(),2);
        for(int i=0;i<int(values.size());++i)packed[i]=(packed[i]<<shift)+integer[i];
        std::cout<<"integer channel "<<ch<<" residual bound bits "<<shift<<std::endl;
    }
    std::vector<mpz_class> charts;
    charts.reserve(d.labels.size());
    for(auto &[key,l]:d.labels) {
        l.record=charts.size();
        charts.emplace_back(packed[l.index]*multipliers.at(l.ratio));
    }
    std::vector<mpz_class>().swap(packed);
    long checked=0;
    mpz_class residual;
    for(const auto &[key,parent]:d.labels) {
        const auto g=unpack(key);
        for(int v=0;v<6;++v)for(int ref=2;ref<6;++ref)for(int inf=2;inf<6;++inf) {
            if(v==ref || v==inf || ref==inf)continue;
            int n=0;
            std::array<std::pair<int,int>,9> raw;
            auto collect=[&](const Label &l,int c){if(c)raw[n++]={l.record,c};};
            if(!visit_row(d,g,parent,v,ref,inf,collect)||!n)continue;
            residual=0;
            for(int j=0;j<n;++j) {
                auto [idx,c]=raw[j];
                if(c>0)mpz_addmul_ui(residual.get_mpz_t(),charts[idx].get_mpz_t(),c);
                else mpz_submul_ui(residual.get_mpz_t(),charts[idx].get_mpz_t(),-c);
            }
            if(residual!=0)throw std::logic_error("exact integer IBP check failed");
            ++checked;
        }
    }
    std::cout<<"all exact IBP checks "<<checked<<std::endl;
}

// Each input row uses the nine position powers in the order of faces[].
void query_periods(const Domain &d,const std::string &prefix) {
    std::vector<Label> labels;
    Powers a;
    while(std::cin>>a[0]) {
        for(int i=1;i<9;++i)if(!(std::cin>>a[i]))throw std::runtime_error("short prism query");
        int degree=0;for(int k:a)degree+=k;
        const auto chart=complete(a,d.q);
        if(degree!=4*d.q || !convergent(chart,d.q))throw std::invalid_argument("invalid prism query");
        labels.push_back(d.labels.at(graph_key(chart)));
    }
    ExactFactors factors;
    std::vector<std::array<mpq_class,3>> answers(labels.size());
    for(int ch=0;ch<3;++ch) {
        const auto values=read_table(d,(prefix+std::to_string(ch)).c_str());
        for(int i=0;i<int(labels.size());++i) {
            const auto &l=labels[i];
            answers[i][ch]=values.at(l.index)*factors.get(l.ratio);
        }
    }
    for(const auto &r:answers)std::cout<<r[0]<<' '<<r[1]<<' '<<r[2]<<'\n';
}

int run_program(int argc,char **argv) {
    bool full=argc>2 && std::string(argv[1])=="--all-periods";
    if(full) {--argc;++argv;}
    if (argc<2 || argc>8) return 2;
    Domain domain(std::stoi(argv[1]),full);
    domain.enumerate();
    if(argc==4 && full && std::string(argv[2])=="--query") {query_periods(domain,argv[3]);return 0;}
    std::cout << "charts " << domain.labels.size() << " variables "
              << domain.representatives.size() << std::endl;
    if(argc==4 && full && std::string(argv[2])=="--check-integers") {check_integer_rows(domain,argv[3]);return 0;}
    if (argc==3) pivots(domain,argv[2]);
    if (argc==4 && std::string(argv[2])=="--group") {group_weights(domain,argv[3]);return 0;}
    if (argc==4) {
        auto exact=read_table(domain,argv[3]);
        if(std::string(argv[2])=="--check-all") check_all(domain,exact);
        else if(std::string(argv[2])=="--sum") physical_weights(domain,exact);
        else if(std::string(argv[2])=="--compare") compare_cuts(domain,exact);
        else if(std::string(argv[2])=="--group-sum") group_sum(exact);
        else return 2;
    }
    if (argc==7 && std::string(argv[3])=="--full-solve") {
        SeedSystem seed(domain,argv[2]);BoundarySystem b(domain,seed,argv[4]);
        full_solve(domain,seed,b,argv[5],argv[6]);return 0;
    }
    if (argc==7 && std::string(argv[3])=="--final-rank") {
        SeedSystem seed(domain,argv[2]);BoundarySystem b(domain,seed,argv[4]);
        choose_full(domain,seed,b,argv[5],std::stoi(argv[6]));return 0;
    }
    if (argc==6 && std::string(argv[3])=="--sparse") {
        SeedSystem s(domain,argv[2]);
        sparse_pass(domain,s,argv[4],std::stoi(argv[5]));
        return 0;
    }
    if (argc==6) {
        if (std::string(argv[4])!="--solve") return 2;
        SeedSystem s(domain,argv[2]);
        auto values=solve_exact(domain,s,argv[3]);
        if(std::string(argv[5])!="-") {
            std::ofstream f(argv[5]);
            if(!f)throw std::runtime_error("cannot write coefficients");
            for(std::size_t i=0;i<values.size();++i) {
                for(int v:domain.representatives[i]) f << v << ' ';
                f << values[i] << '\n';
            }
        }
        physical_weights(domain,values);
    }
    if (argc==5) {SeedSystem s(domain,argv[2]);choose(domain,s,argv[3],std::stoi(argv[4]));}
    return 0;
}

int main(int argc,char **argv) {
    try {
        return run_program(argc,argv);
    } catch(const std::exception &error) {
        std::cerr<<"prism reduction failed: "<<error.what()<<std::endl;
        return 1;
    }
}
