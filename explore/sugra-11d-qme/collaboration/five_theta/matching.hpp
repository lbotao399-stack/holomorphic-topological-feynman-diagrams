#pragma once
#include "theta_numerator.hpp"
#include <algorithm>
#include <array>
#include <vector>

namespace twoloop {
using Word=std::array<int,10>;
template<class T> inline int permutation_sign(const T& sequence) {
    int sign=1;
    for(std::size_t i=0;i<sequence.size();++i)for(std::size_t j=i+1;j<sequence.size();++j)
        if (sequence[i]>sequence[j]) sign=-sign;
    return sign;
}
inline std::array<int,10> momentum_signs(const Word& w) {
    std::array<int,10> signs{};std::array<bool,5> seen{};
    for (int j=0;j<10;++j) {signs[j]=seen[w[j]]?-1:1;seen[w[j]]=true;}
    return signs;
}
inline int determinant(const Word& w) {return -permutation_sign(w);}
inline int direct_determinant(const Word& w,const Lengths& lengths) {
    int matrix[10][10]{};
    auto coeff=momentum_signs(w);
    for(int v=0;v<10;++v) {
        int path=v<lengths[0]-1?0:(v<lengths[0]+lengths[1]-2?1:2);
        if(path<2)matrix[w[v]][v]=coeff[v];
        if(path>0)matrix[w[v]+5][v]=coeff[v];
    }
    int sign=1, prev=1;
    for(int j=0;j<9;++j) {
        int k=j;while(k<10&&!matrix[k][j])++k;
        if(k==10)return 0;
        if(k!=j){std::swap(matrix[k],matrix[j]);sign=-sign;}
        int pivot=matrix[j][j];
        for(int a=j+1;a<10;++a)for(int b=j+1;b<10;++b)
            matrix[a][b]=(matrix[a][b]*pivot-matrix[a][j]*matrix[j][b])/prev;
        prev=pivot;
    }
    return sign*matrix[9][9];
}
inline std::vector<Word> matchings(const Lengths& lengths) {
    if(lengths[0]+lengths[1]+lengths[2]!=13 || *std::min_element(lengths.begin(),lengths.end())<1 ||
       *std::max_element(lengths.begin(),lengths.end())>6)throw std::invalid_argument("invalid lengths");
    std::vector<Word> result;
    Word word{};std::array<int,5> counts{},pathmask{};
    auto visit=[&](auto&& self,int pos,int introduced)->void {
        if(pos==10){
            if(introduced!=5)return;
            if (direct_determinant(word,lengths)!=determinant(word))throw std::logic_error("determinant parity mismatch");
            result.push_back(word);return;
        }
        int path=pos<lengths[0]-1?0:(pos<lengths[0]+lengths[1]-2?1:2);
        for(int axis=0;axis<std::min(introduced+1,5);++axis){
            if(counts[axis]==2 || (pathmask[axis]&(1<<path)))continue;
            word[pos]=axis;++counts[axis];pathmask[axis]|=1<<path;
            self(self,pos+1,std::max(introduced,axis+1));
            --counts[axis];pathmask[axis]^=1<<path;
        }
    };
    visit(visit,0,0);
    return result;
}
} // namespace twoloop
