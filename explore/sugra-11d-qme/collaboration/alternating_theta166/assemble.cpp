// Exterior coefficient for the transverse input specified in transverse_theta166.md.
#include "numerator.hpp"
#include <algorithm>
#include <array>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using Order = std::array<int, 5>;

template<class Indices>
int sign_of(const Indices& a) {
    int sign = 1;
    for (size_t i=0; i<a.size(); ++i)
        for (size_t j=i+1; j<a.size(); ++j)
            if (a[i]>a[j]) sign = -sign;
    return sign;
}

struct Wave {
    int momentum;       // +(i+1) or -(i+1)
    int vector_axis;    // zero-based polarization axis (dual polarization for gamma)
    int odd_position;   // in exterior product of the nine vector inputs, -1 for gamma
};

const std::array<Wave,10> waves = {{
    {1,1,-1}, {2,2,0}, {3,3,1}, {4,4,2}, {5,1,3},
    {-1,1,4}, {-2,0,5}, {-4,0,6}, {-3,0,7}, {-5,0,8}
}};

const Wave& wave(int momentum) {
    for (const auto& w : waves) if (w.momentum == momentum) return w;
    throw std::logic_error("unknown wave");
}

int main(int argc, char** argv) {
    if (argc != 5) {
        std::cerr << "Usage: assemble gamma-slot[0..2] first-order[0..119] end-exclusive output\n";
        return 2;
    }
    const int gamma = std::stoi(argv[1]), first = std::stoi(argv[2]), end = std::stoi(argv[3]);
    if(gamma < 0 || gamma > 2 || first < 0 || end > 120 || first >= end) return 2;
    std::ofstream output(argv[4]); if(!output) return 2;
    // Relative order is the sequence of B-axis positions seen along C.
    Order order = {0,1,2,3,4};
    for (int rank=0; rank<end; ++rank,std::next_permutation(order.begin(),order.end())) {
        if(rank<first) continue;
        theta166::Polynomial aggregate;
        size_t count=0;
        Order b_axes={0,1,2,3,4};
        do {
            if (b_axes[gamma]!=0) continue;
            for (int mask=0; mask<16; ++mask) {
                Order sign={},position={};
                int det=sign_of(b_axes);
                for (int i=0; i<5; ++i) {
                    const int phys=b_axes[i];
                    sign[i]=phys==0 || !(mask & (1 << (phys-1))) ? 1 : -1;
                    position[phys]=i;
                    det*=sign[i];
                }
                theta166::Slots slots{};
                std::vector<int> odd;
                // Each cubic gamma vertex has one volume tensor; D_theta has two.
                int prefactor=det * (-sign_of(order));
                for(int v=0;v<10;++v) {
                    int i=v<5 ? v : order[v-5];
                    int side=v<5 ? 1 : -1;
                    const auto& physical=wave(side*sign[i]*(b_axes[i]+1));
                    int axis=position[physical.vector_axis];
                    slots[v]={side*(i+1), axis};
                    prefactor*=sign[axis];
                    if (physical.odd_position>=0) odd.push_back(physical.odd_position);
                }
                prefactor*=sign_of(odd);
                auto polynomial=theta166::numerator(slots,gamma);
                for(const auto& [key,coeff]:polynomial) aggregate[key]+=prefactor*coeff;
                ++count;
            }
        } while(std::next_permutation(b_axes.begin(),b_axes.end()));
        if(count!=384) throw std::logic_error("wrong number of labellings");
        std::vector<std::pair<uint64_t,int64_t>> terms;
        for(const auto& item:aggregate) if(item.second) terms.push_back(item);
        std::sort(terms.begin(),terms.end());
        output<<rank<<" "<<terms.size()<<"\n";
        for(const auto& [key,coef]:terms) output<<key<<" "<<coef<<"\n";
        if(!output) return 1;
        if(rank==first || (rank+1)%10==0)
            std::cerr<<"slot "<<gamma<<" order "<<rank<<" terms "<<terms.size()<<"\n";
    }
}
