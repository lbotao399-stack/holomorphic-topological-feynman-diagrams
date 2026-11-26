"""Exact ordered-fixture regression for the physical --gamma-mask option (C++17)."""
from collections import Counter
from hashlib import sha256
from pathlib import Path
from tempfile import TemporaryDirectory
import os
import shlex
import subprocess

HERE=Path(__file__).resolve().parent
# One nonzero matching block for each of four higher coupling orders and five graphs.
CASES={
    '166': ((29,0x130),(66,0x174),(16,0xdf),(106,0x3df)),
    '256': ((69,0x89),(46,0x25a),(35,0x1af),(63,0x3fe)),
    '346': ((107,0x320),(87,0x2cc),(38,0xfe),(46,0x3fe)),
    '355': ((32,0x224),(143,0x11d),(131,0x1bb),(68,0x3df)),
    '445': ((121,0x228),(22,0xe6),(86,0x3b9),(21,0x3df)),
}

def run(args):
    return subprocess.run(args,check=True,text=True,capture_output=True).stdout


def coefficients(file):
    lines=file.read_text().splitlines();rank,n=map(int,lines[0].split())
    polynomial={int(key):int(c) for key,c in (row.split() for row in lines[1:])}
    assert len(polynomial)==n
    return rank,polynomial


def main():
    with TemporaryDirectory(prefix='gamma-mask-cli-') as tmp:
        base=Path(tmp);assembler=base/'assemble';physical=base/'check_physical';integer=base/'integer'
        for source,target,extra in [('assemble',assembler,[]),('check_physical',physical,[]),('check_physical',integer,['-DEXACT_INTEGER'])]:
            run([*shlex.split(os.environ.get('CXX','c++')),'-O3','-std=c++17',*extra,'-fsanitize=signed-integer-overflow',
                 '-fno-sanitize-recover=all',HERE/(source+'.cpp'),'-o',target])
        output=base/'case.poly'
        for tail,digest in [((), 'de375d908e51af4036095e36005b1efa9ecf071012f9315b00543201bd8613d6'),
                            (('0',),'71b0c284cd9ab4113e2d42b29ba5ffe448343715146281706489de81c8974bca')]:
            for flags in ((),('--gamma-mask','0x1')):
                run([assembler,*'166','0','1',output,*tail,*flags])
                assert sha256(output.read_bytes()).hexdigest()==digest
        for exe,args in [(assembler,[*'166','0','1',output]),(physical,[*'166','4',output])]:
            for invalid in ('0','3','-1','1024','0x1000','not-a-number'):
                assert subprocess.run([exe,*args,'--gamma-mask',invalid],capture_output=True).returncode==2
        for shape,cases in CASES.items():
            for m,(rank,mask) in enumerate(cases,2):
                assert mask.bit_count()==2*m-1
                bounds=[str(rank),str(rank+1)]
                run([assembler,*shape,*bounds,output,'--gamma-mask',hex(mask)])
                block,poly=coefficients(output);assert block==rank and poly
                for checker in (physical,integer):
                    result=run([checker,'--gamma-mask',str(mask),*shape,'4',output])
                    assert '1 blocks x 3840 physical assignments x 4 flows;' in result and '; 1 nonzero,' in result,result
                if shape=='166' and m==2:
                    combined=Counter()
                    for slot in range(10):
                        run([assembler,*shape,*bounds,output,str(slot),'--gamma-mask',hex(mask)])
                        combined.update(coefficients(output)[1])
                    assert {k:v for k,v in combined.items() if v}==poly
                print(shape,'m',m,'mask',hex(mask),'rank',rank,'pointwise exact and F1009 checks passed',flush=True)
    print('Default and explicit g^1 masks reproduced prior bytes. Ten anchor-slot polynomials sum to the unrestricted block.')

if __name__=='__main__':main()
