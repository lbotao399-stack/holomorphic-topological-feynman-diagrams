#pragma once

#include "theta_numerator.hpp"
#include <bitset>

namespace twoloop {
using GammaSites = std::bitset<10>;

namespace detail {
// Each external gamma reverses the direction along its path (a b source or c sink).
// An external mu retains the direction. At N either zero or two edges leave.
inline Directions multi_directions(const Graph& g,const GammaSites& sites,unsigned outgoing_at_N) {
    Directions flow;
    for (int path=0;path<3;++path) {
        bool forward=(outgoing_at_N>>path)&1U;
        for (int edge=g.first[path];edge<=g.last[path];++edge) {
            flow.forward[edge]=forward;
            if (edge!=g.last[path] && sites.test(g.head[edge]-1)) forward=!forward;
        }
    }
    flow.set_wick_sign(g);
    return flow;
}
// Delay expansion until after summing all transported polynomials in one pairing.
inline void collect_gamma_sites(const Graph& graph,const GammaSites& sites,int sign,Collected& factored) {
    for (unsigned outgoing : {0U,3U,5U,6U}) {
        auto flow=multi_directions(graph,sites,outgoing);
        flow.sign*=sign;
        collect(graph,flow,factored);
    }
}
} // namespace detail

// K=0. Bit j selects an external gamma at the (j+1)-st interior vertex of A,B,C.
// Return the sum of the four oriented numerators stripped of g^m, D_theta,
// the Schwinger measure, graph multiplicity and external polarization permutations.
inline Polynomial multi_gamma_numerator(const Lengths& lengths,const Slots& slots,const GammaSites& sites) {
    if (sites.count()%2!=1) throw std::invalid_argument("expected an odd number of external gamma sites");
    int first=0;
    while (!sites.test(first)) ++first;
    detail::Graph graph(lengths,slots,first);
    detail::Collected factored;
    detail::collect_gamma_sites(graph,sites,1,factored);
    // Of the 12+m holomorphic derivatives in the vertices, one external momentum
    // can be factored at each of the 2m-1 external gamma vertices by conservation.
    // The remaining current degree is 13-m <= 12, within the four-bit exponent format.
    return detail::expand(graph,factored);
}
} // namespace twoloop
