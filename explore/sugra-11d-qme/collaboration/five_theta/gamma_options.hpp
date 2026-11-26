#pragma once
#include <array>
#include <exception>
#include <string>
#include <vector>

// Bit i selects the ith wave in the fixed physical fixture, independently of
// the vertex to which that wave is assigned. Decimal and 0x-prefixed masks work.
struct GammaOptions {
    unsigned mask=1;
    // Zero-based internally, in fixture wave order +1,+2,+3,+4,+5,-1,-2,-4,-3,-5.
    std::array<int,10> polarizations={1,2,3,4,1,1,0,0,0,0};
    std::vector<std::string> positional;
    bool parse(int argc,char** argv) {
        bool supplied=false, axes_supplied=false;
        for(int i=1;i<argc;++i) {
            std::string flag=argv[i];
            if(flag=="--polarizations") {
                if(axes_supplied || ++i==argc)return false;
                axes_supplied=true;
                std::string axes=argv[i];
                if(axes.size()!=19)return false;
                for(int k=0;k<10;++k) {
                    char item=axes[2*k];
                    if(item<'1' || item>'5' || (k<9 && axes[2*k+1]!=','))return false;
                    polarizations[k]=item-'1';
                }
                continue;
            }
            if(flag!="--gamma-mask") {positional.emplace_back(flag);continue;}
            if(supplied || ++i==argc)return false;
            supplied=true;std::string value=argv[i];
            try {
                if(value.empty() || value[0]=='-')return false;
                bool hex=value.size()>2 && value[0]=='0' && (value[1]=='x' || value[1]=='X');
                std::size_t end;auto bits=std::stoul(value,&end,hex?16:10);
                if(end!=value.size() || bits>=1024 || !(__builtin_popcount(static_cast<unsigned>(bits))&1))return false;
                mask=static_cast<unsigned>(bits);
            } catch(const std::exception&) {return false;}
        }
        return true;
    }
};
