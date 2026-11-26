// Bounded GL_5 Weyl-character extraction for soft momentum ideals on a conserved hyperplane.
#include <algorithm>
#include <array>
#include <iostream>
#include <numeric>
#include <vector>
using A=std::array<int,5>;
using I=__int128;
void output(I n){if(n<0){std::cout<<'-';n=-n;}if(n>=10)output(n/10);std::cout<<char('0'+n%10);}
struct Character {
    A lambda,cap,stride;
    int count=1;
    std::vector<I> p;
    Character(A shape):lambda(shape) {
        for(int i=0;i<5;++i){cap[i]=lambda[i]+4-i;stride[i]=count;count*=1+cap[i];}
        p.assign(count,0);p[0]=1;
    }
    void cycle(int len,int minDegree) {
        auto old=p;
        for(int axis=0;axis<5;++axis){
            int step=stride[axis],block=(cap[axis]+1)*step,delta=len*step;
            for(int b=0;b<count;b+=block)for(int j=b+delta;j<b+block;++j)p[j]+=p[j-delta];
        }
        if(minDegree==0)return;
        for(int i=0;i<count;++i)p[i]-=old[i];
        if(minDegree==1)return;
        for(int axis=0;axis<5;++axis){
            int step=stride[axis],block=(cap[axis]+1)*step,delta=len*step;
            for(int b=0;b<count;b+=block)for(int j=b+delta;j<b+block;++j)p[j]-=old[j-delta];
        }
    }
    I irreducible(bool conserve) {
        if(conserve)for(int axis=0;axis<5;++axis){
            int dx=stride[axis],block=(cap[axis]+1)*dx;
            for(int b=0;b<count;b+=block)for(int j=b+block-1;j>=b+dx;--j)p[j]-=p[j-dx];
        }
        A perm{0,1,2,3,4};I sum=0;
        do {
            int index=0,sign=1;bool ok=true;
            for(int i=0;i<5;++i){
                int e=cap[i]-(4-perm[i]);if(e<0){ok=false;break;}index+=e*stride[i];
                for(int j=0;j<i;++j)if(perm[j]>perm[i])sign=-sign;
            }
            if(ok)sum+=sign*p[index];
        }while(std::next_permutation(perm.begin(),perm.end()));
        return sum;
    }
};
int main(){
    int weight,kind,conserve,n;
    while(std::cin>>weight>>kind>>conserve>>n){
        A shape{weight,weight,weight,weight,weight};
        for(int i=0;i<kind;++i)++shape[i];
        Character c(shape);
        for(int j=0;j<n;++j){int len,m;std::cin>>len>>m;c.cycle(len,m);}
        output(c.irreducible(conserve));std::cout<<'\n';
    }
}
