"""Bounded symbolic transport, literal Wick and exact moment checks for custom unit polarizations."""
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
from random import Random
from tempfile import TemporaryDirectory
import os
import shlex
import subprocess
import sys

sys.dont_write_bytecode=True
from replay import Commands, gmp_flags
from verify_multi_gamma import oracle, orientations, theta_graph

HERE=Path(__file__).resolve().parent
WAVES=(1,2,3,4,5,-1,-2,-4,-3,-5)
DEFAULT='2,3,4,5,2,2,1,1,1,1'
SAMPLE=r'''
#define main assemble_entry_point
#include "assemble.cpp"
#undef main
int main(int argc,char** argv){
    GammaOptions options;if(!options.parse(argc,argv))return 2;
    twoloop::Lengths shape;std::array<int,5> permutation,inverse;int rank,signs;
    std::cin>>shape[0]>>shape[1]>>shape[2]>>rank>>signs;
    for(int i=0;i<5;++i){std::cin>>permutation[i];inverse[permutation[i]]=i;}
    auto word=twoloop::matchings(shape).at(rank);auto sides=twoloop::momentum_signs(word);
    auto data=transport(word,sides,options,permutation,inverse,signs);
    for(int x:word)std::cout<<x;std::cout<<'\n';
    auto p=twoloop::multi_gamma_numerator(shape,data.slots,data.sites);
    for(auto [key,coefficient]:p)std::cout<<key<<' '<<data.factor*coefficient<<'\n';
}
'''

def run(args,stdin=None):
    return subprocess.run(args,input=stdin,text=True,capture_output=True,check=True).stdout


def parity(seq):
    return (-1)**sum(a>b for i,a in enumerate(seq) for b in seq[i+1:])


def check_transported(exe):
    rng=Random(1103241)
    for shape,rank,mask in [((2,5,6),13,0x130),((4,4,5),54,0x3e),((3,4,6),2,0x82)]:
        if mask.bit_count()%2==0:mask|=1
        nv,edges,paths=theta_graph(*shape)
        choices=orientations(nv,edges,)
        for _ in range(100):
            polars=[rng.choice([i for i in range(5) if not (mask>>l&1) or i!=abs(q)-1]) for l,q in enumerate(WAVES)]
            perm=rng.sample(range(5),5);reflection=rng.randrange(32)
            signs=[-1 if reflection>>j&1 else 1 for j in range(5)]
            cli=['--gamma-mask',hex(mask),'--polarizations',','.join(str(i+1) for i in polars)]
            args=' '.join(map(str,(*shape,rank,reflection,*perm)))+'\n'
            word,*rows=run([exe,*cli],args).splitlines()
            actual=dict(tuple(map(int,row.split())) for row in rows)
            if actual:break
        else:raise AssertionError('no nonzero symbolic control')
        seen=set();physical=[]
        for x in map(int,word):
            physical.append((-1 if x in seen else 1)*signs[x]*(perm[x]+1));seen.add(x)
        labels=[WAVES.index(m) for m in physical];slots=sum(1<<v for v,l in enumerate(labels) if mask>>l&1)
        own={v+1:[int(abs(m)-1==i)*(1 if m>0 else -1) for i in range(5)] for v,m in enumerate(physical)}
        from vcommon import D_theta
        factor=int(D_theta(own,paths))*parity([l for l in labels if not(mask>>l&1)])
        full=oracle(nv,edges,paths,choices[slots],physical,[polars[l] for l in labels])
        # The reference did all contractions in the physical basis. Pull its two
        # first-edge current vectors back only now with the actual signed permutation.
        expected={}
        for key,coeff in full.items():
            newkey=0;s=factor
            for offset in (0,5):
                for j,k in enumerate(perm):
                    power=(key>>(4*(offset+k)))&15
                    newkey|=power<<(4*(offset+j));s*=signs[j]**power
            expected[newkey]=s*coeff
        assert actual==expected,(shape,rank,mask)
        print(shape,'mask',hex(mask),len(actual),'custom transported coefficients matched',flush=True)


def main():
    with TemporaryDirectory(prefix='theta-polarizations-') as tmp:
        base=Path(tmp);flags=gmp_flags(Commands());cxx=shlex.split(os.environ.get('CXX','c++'))
        for name in ('assemble','check_physical','moments','independent'):
            extra=flags if name in ('moments','independent') else []
            run([*cxx,'-std=c++17','-O3','-DEXACT_INTEGER','-fsanitize=signed-integer-overflow','-fno-sanitize-recover=all',HERE/(name+'.cpp'),*extra,'-o',base/name])
        run([*cxx,'-std=c++17','-O2','-I',HERE,'-x','c++','-','-o',base/'sample'],SAMPLE)
        check_transported(base/'sample')
        poly=base/'sample.poly';assembly=[base/'assemble',*'346','8','9',poly];direct=[base/'check_physical',*'346','4',poly]
        for bad in ('','1,2','0,2,3,4,5,1,2,3,4,5','6,2,3,4,5,1,2,3,4,5',DEFAULT+',',DEFAULT.replace(',',';',1)):
            for cmd in (assembly,direct):assert subprocess.run([*cmd,'--polarizations',bad],capture_output=True).returncode==2
        for cmd in (assembly,direct):assert subprocess.run([*cmd,'--polarizations',DEFAULT,'--polarizations',DEFAULT],capture_output=True).returncode==2
        run([base/'assemble',*'166','0','1',poly,'--polarizations',DEFAULT]);assert sha256(poly.read_bytes()).hexdigest()=='de375d908e51af4036095e36005b1efa9ecf071012f9315b00543201bd8613d6'
        axes='2,5,2,2,2,2,5,1,1,5';settings=['--gamma-mask','0x130','--polarizations',axes]
        run([*assembly,*settings]);assert '; 1 nonzero,' in run([*direct,*settings]);assert len(poly.read_text().splitlines())==25164
        expected=Fraction(-21889056259327,306013455715968000000)
        for tool in ('moments','independent'):
            assert Fraction(run([base/tool,*'346',poly]).split()[-1])==expected
        axes=axes.split(',');axes[4]='5'  # gamma at physical +e^5 has covector e^5
        long=['--polarizations',','.join(axes),'--gamma-mask','0x130'];run([*assembly,*long]);assert poly.read_text()=='8 0\n'
        assert '; 0 nonzero,' in run([*direct,*long]);print('custom integrated kernel/pi',expected,'; longitudinal gamma zero')

if __name__=='__main__':main()
