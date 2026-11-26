// Exact zeroth Hepp moment of both orbit-summed four-loop BCOV kernels.
// Input: nine core edges, followed by nine powers and the two cut-normalized weights per row.
#include <gmpxx.h>
#include <array>
#include <algorithm>
#include <climits>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

constexpr unsigned edges = 9, full = (1u << edges) - 1;
using Integer = mpz_class;
using Rational = mpq_class;

struct Moments {
    std::array<int,512> cycles{}, split{}, power_sum{};
    std::array<uint64_t,512> digits{};
    std::array<Integer,512> scratch{};
    std::array<bool,512> visited{};
    std::array<std::unordered_map<uint64_t,Integer>,512> cache;
    std::array<Integer,61> division{};
    Integer scale = 1, denominator;
    uint64_t code = 0;
    int q;

    explicit Moments(const std::array<std::pair<int,int>,edges>& graph,int power=12):q(power) {
        for (unsigned mask=0; mask<=full; ++mask) {
            std::array<int,6> parents{0,1,2,3,4,5};
            auto root = [&](int i) { while (parents[i]!=i) i=parents[i]; return i; };
            for (unsigned i=0; i<edges; ++i) {
                if (!(mask & (1u<<i))) continue;
                auto [v,w] = graph[i];
                int a=root(v), b=root(w);
                if (a==b) ++cycles[mask]; else parents[a]=b;
                digits[mask] |= uint64_t(15)<<(4*i);
            }
            for (unsigned a=(mask-1)&mask; a; a=(a-1)&mask)
                if (cycles[a]+cycles[a^mask]==cycles[mask]) { split[mask]=a; break; }
        }
        if (cycles[full]!=4) throw std::invalid_argument("wrong graph genus");
        for (int i=1; i<=4*q; ++i) mpz_lcm_ui(scale.get_mpz_t(),scale.get_mpz_t(),i);
        for (int i=1; i<=4*q; ++i) division[i]=scale/i;
        mpz_pow_ui(denominator.get_mpz_t(),scale.get_mpz_t(),8);
        scratch[0]=1;
    }

    const Integer& volume(unsigned mask) {
        if (!mask || visited[mask]) return scratch[mask];
        visited[mask]=true;
        auto& out=scratch[mask];
        uint64_t key=code&digits[mask];
        bool save=__builtin_popcount(mask)<=6;
        if (save) {
            auto found=cache[mask].find(key);
            if (found!=cache[mask].end()) return out=found->second;
        }
        unsigned part=split[mask];
        if (part && mask!=full) out=volume(part)*volume(mask^part);
        else {
            out=0;
            for (unsigned s=mask; s; s&=s-1) out+=volume(mask^(s&-s));
            if (mask!=full) {
                int margin=power_sum[mask]-q*cycles[mask];
                if (margin<=0) throw std::invalid_argument("nonpositive convergence margin");
                out*=division[margin];
            }
        }
        if (save && cache[mask].size()<1536) cache[mask].emplace(key,out);
        return out;
    }

    const Integer& operator()(const std::array<int,edges>& powers) {
        code=0;
        visited.fill(false);
        for (unsigned i=0; i<edges; ++i) {
            if (powers[i]<1 || powers[i]>=q) throw std::invalid_argument("invalid propagator power");
            code |= uint64_t(powers[i])<<(4*i);
        }
        for (unsigned mask=1; mask<=full; ++mask) {
            unsigned bit=mask&-mask;
            power_sum[mask]=power_sum[mask^bit]+powers[__builtin_ctz(bit)];
        }
        if (power_sum[full]!=4*q) throw std::invalid_argument("nonprojective exponents");
        return volume(full);
    }
};

using SmallRow = std::array<__int128_t,2>;
using Polynomial = std::unordered_map<uint64_t,SmallRow>;
using Graph = std::array<std::pair<int,int>,edges>;
struct Permutation {
    std::array<unsigned,edges> target;
    std::vector<uint64_t> shifts;
    uint64_t operator()(uint64_t code) const {
        uint64_t result=0;
        for (unsigned i=0; i<edges; ++i) result |= ((code>>(4*i))&15)<<(4*target[i]);
        return result;
    }
};

std::vector<Permutation> symmetries(const Graph& graph,const Moments& moment) {
    std::vector<uint64_t> trees;
    for (unsigned mask=0; mask<=full; ++mask)
        if (__builtin_popcount(mask)==5 && !moment.cycles[mask]) {
            uint64_t code=0;
            for (unsigned i=0; i<edges; ++i) if (!(mask&(1u<<i))) code|=uint64_t(1)<<(4*i);
            trees.push_back(code);
        }
    std::vector<Permutation> result;
    std::array<int,6> order{0,1,2,3,4,5};
    do {
        Permutation p;
        bool valid=true;
        for (unsigned i=0; i<edges; ++i) {
            auto [a,b]=graph[i];
            auto pair=std::minmax(order[a],order[b]);
            auto it=std::find(graph.begin(),graph.end(),std::pair<int,int>(pair));
            if (it==graph.end()) { valid=false; break; }
            p.target[i]=it-graph.begin();
        }
        if (!valid) continue;
        for (auto t:trees) p.shifts.push_back(p(t));
        result.push_back(std::move(p));
    } while (std::next_permutation(order.begin(),order.end()));
    if (!((trees.size()==75 && result.size()==12) || (trees.size()==81 && result.size()==72)))
        throw std::invalid_argument("unexpected primitive graph symmetry");
    return result;
}

void add_product(Polynomial& next,const std::vector<Permutation>& group,uint64_t code,const SmallRow& row) {
    size_t count=group.front().shifts.size();
    std::array<uint64_t,81> keys;
    keys.fill(UINT64_MAX);
    for (const auto& p:group) {
        uint64_t base=p(code);
        for (size_t i=0; i<count; ++i) keys[i]=std::min(keys[i],base+p.shifts[i]);
    }
    for (size_t i=0; i<count; ++i) {
        auto& target=next[keys[i]];
        for (unsigned j=0; j<2; ++j)
            if (__builtin_add_overflow(target[j],row[j],&target[j]))
                throw std::overflow_error("orbit coefficient exceeds signed 128 bits");
    }
}

int main(int argc,char** argv) {
    int order=argc==3 ? std::stoi(argv[1]) : 0;
    if ((argc!=1 && argc!=3) || order<0 || order>3)
        throw std::invalid_argument("expected optional raw moment order 1..3 and coefficient denominator");
    Integer coefficient_den=argc==3 ? Integer(argv[2]) : Integer(1);
    std::array<std::pair<int,int>,edges> graph;
    for (auto& [a,b]:graph)
        if (!(std::cin>>a>>b) || a<0 || b<0 || a>=6 || b>=6)
            throw std::invalid_argument("expected nine edges on six vertices");
    for (auto& [a,b]:graph) if (a>b) std::swap(a,b);
    Moments moment(graph);
    auto group=symmetries(graph,moment);
    Polynomial polynomial;
    std::array<int,edges> powers;
    std::array<Integer,12> factorial;
    factorial[0]=1;
    for (int i=1; i<=11; ++i) factorial[i]=factorial[i-1]*i;
    std::array<Rational,2> sum{}, absolute{};
    size_t count=0;
    while (std::cin>>powers[0]) {
        for (unsigned i=1; i<edges; ++i)
            if (!(std::cin>>powers[i])) throw std::invalid_argument("partial power row");
        uint64_t code=0;
        for (unsigned j=0; j<edges; ++j) {
            if (powers[j]<=0 || powers[j]>=12) throw std::invalid_argument("invalid input powers");
            code |= uint64_t(powers[j])<<(4*j);
        }
        Integer normalization_den=1;
        for (int p:powers) normalization_den*=factorial[p-1];
        Rational scaled(factorial[11],normalization_den);
        scaled.canonicalize();
        SmallRow small{};
        if (!order) scaled*=moment(powers);
        for (unsigned j=0; j<2; ++j) {
            std::string weight;
            if (!(std::cin>>weight)) throw std::invalid_argument("incomplete BCOV weights");
            Rational term=Rational(weight)*scaled;
            if (order) {
                term*=coefficient_den;
                if (term.get_den()!=1 || !mpz_fits_slong_p(term.get_num_mpz_t()))
                    throw std::overflow_error("invalid common denominator or initial bound");
                small[j]=term.get_num().get_si();
            } else {
                sum[j]+=term;
                absolute[j]+=abs(term);
            }
        }
        if (order) add_product(polynomial,group,code,small);
        ++count;
    }
    if (order) {
        for (int k=1; k<=order; ++k) {
            Moments shifted(graph,12+k);
            std::array<Integer,2> result{};
            for (const auto& [key,row]:polynomial) {
                for (unsigned i=0; i<edges; ++i) powers[i]=(key>>(4*i))&15;
                const auto& value=shifted(powers);
                for (unsigned j=0; j<2; ++j) {
                    if (row[j]>=LONG_MIN && row[j]<=LONG_MAX) result[j]+=static_cast<long>(row[j])*value;
                    else {
                        Integer multiple;
                        __uint128_t magnitude=row[j]<0 ? -static_cast<__uint128_t>(row[j]) : row[j];
                        mpz_import(multiple.get_mpz_t(),1,1,sizeof(magnitude),0,0,&magnitude);
                        if (row[j]<0) multiple=-multiple;
                        result[j]+=multiple*value;
                    }
                }
            }
            std::cout<<count<<' '<<polynomial.size();
            for (const auto& v:result) {
                Rational out(v,coefficient_den*shifted.denominator);
                out.canonicalize();
                std::cout<<' '<<out;
            }
            std::cout<<std::endl;
            if (k<order) {
                Polynomial next;
                next.reserve(2*polynomial.size());
                for (auto it=polynomial.begin(); it!=polynomial.end();) {
                    if (it->second[0] || it->second[1]) add_product(next,group,it->first,it->second);
                    it=polynomial.erase(it);
                }
                polynomial.swap(next);
            }
        }
        return 0;
    }
    std::cout<<count<<'\n';
    for (unsigned i=0; i<2; ++i)
        std::cout<<sum[i]/moment.denominator<<' '<<absolute[i]/moment.denominator<<'\n';
}
