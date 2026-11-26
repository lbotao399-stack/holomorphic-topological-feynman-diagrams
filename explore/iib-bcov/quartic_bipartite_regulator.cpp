// Check canceled divergent terms along the common analytic regulator a_e -> a_e+z.
#define BCOV_K33_ORBIT_LIBRARY
#include "quartic_bipartite_orbit_ibp.cpp"

int main(int argc, char** argv) {
    if (argc != 2) throw std::invalid_argument("pole order q required");
    Domain domain(std::stoi(argv[1]));
    std::unordered_map<Key,int> contacts;
    uint64_t accepted=0, canceled=0, exceptional=0;
    for (int i=0; i<int(domain.keys.size()); ++i) {
        for (int rule=0; rule<6; ++rule) {
            Row row;
            if (!generate(domain,i,rule,row,contacts)) continue;
            bool positive=false;
            for (int j=0; j<row.size; ++j) positive |= row.i[j]>=0;
            if (!positive) continue;
            ++accepted;
            const int last=rule%3;
            auto slot=[&](int r,int c){return rule<3 ? 3*r+c : 3*c+r;};
            Matrix a=unpack(domain.keys[i]);
            std::map<Key,std::pair<int,int>> divergent;
            for (int col=0; col<3; ++col) {
                const int dest=slot(last,col);
                for (int step : {1,2}) {
                    const int src=slot((last+step)%3,col);
                    const int sign=step==1 ? 1 : -1;
                    --a[src]; ++a[dest];
                    if (!finite(a,domain.q)) {
                        auto& coeff=divergent[canonical(a)];
                        coeff.first+=sign*(a[dest]-1);
                        coeff.second+=sign;
                    }
                    ++a[src]; --a[dest];
                }
            }
            for (const auto& [key,coeff] : divergent) {
                ++canceled;
                if (coeff.first != 0) throw std::logic_error("nonzero divergent coefficient");
                if (coeff.second != 0) {
                    if (++exceptional <= 4)
                        std::cerr << "regulator residual " << domain.keys[i] << " rule " << rule
                                  << " at " << key << " derivative " << coeff.second << '\n';
                }
            }
        }
    }
    std::cout << "accepted " << accepted << " canceled divergent orbits " << canceled
              << " nonzero regulator residuals " << exceptional << '\n';
    return exceptional != 0;
}
