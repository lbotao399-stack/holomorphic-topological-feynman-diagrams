"""Project a complete mixed prism placement in the original nine core lengths."""
from argparse import ArgumentParser
from contextlib import ExitStack
from hashlib import sha256
from heapq import merge
from itertools import groupby
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory
import gzip
import json
import os
import shlex
import sys

from flint import fmpq
from mixed_bare_symbolic import classified_orbits


def weighted_rows(source):
    previous=-1
    for line in source:
        encoded,value=line.split()
        key=int(encoded)
        assert key>previous
        previous=key
        yield key,value


def merge_weights(filenames,destination):
    signature=sha256()
    count=0
    pending=destination.with_suffix('.partial')
    with ExitStack() as stack,gzip.open(pending,'wb') as fout:
        inputs=[weighted_rows(stack.enter_context(gzip.open(path,'rt'))) for path in filenames]
        for key,rows in groupby(merge(*inputs,key=lambda row:row[0]),key=lambda row:row[0]):
            value=sum((fmpq(v) for _,v in rows),fmpq(0))
            if value:
                line=f'{key} {value}\n'.encode()
                signature.update(line);fout.write(line);count+=1
    pending.replace(destination)
    return count,signature.hexdigest()


def bank_info(path):
    signature=sha256()
    count=0
    with gzip.open(path,'rb') as stream:
        for block in iter(lambda:stream.read(1048576),b''):
            signature.update(block)
            count+=block.count(b'\n')
    return count,signature.hexdigest()


def project_placement(word,folder):
    from quartic_prism_physical_ibp import compile_program

    marks=tuple(map(int,word))
    selected=[(color,n) for m,color,n,w in classified_orbits(3)[0] if m==marks]
    if not selected:raise ValueError('the prism placement must be a nonzero canonical representative')
    folder=Path(folder)
    folder.mkdir(parents=True,exist_ok=True)
    source=Path(__file__).with_name('mixed_bipartite_simplex_square.py')
    rows={}
    files=[]
    for color,n in selected:
        tag=''.join(map(str,color))
        file=folder/f'{word}-{tag}.gz'
        if not file.exists():
            run([sys.executable,'-Bu',str(source),'--core','prism',word,tag,str(file)],check=True)
        count,signature=bank_info(file)
        rows[tag]=dict(nonzero_hard_orders=n,core_orbits=count,weights_sha256=signature)
        files.append(file)
        print('Prism color',word,tag,count,flush=True)
    bank=folder/f'{word}-direct.gz'
    with TemporaryDirectory(dir=folder) as temp:
        program=Path(temp)/'merge'
        compile_program('operator/mixed_prism_weight_merge.cpp',program,('-lz',))
        pending=Path(temp)/'merged.gz'
        merged=run([program,pending,*files],capture_output=True,text=True,check=True)
        count,signature=bank_info(pending)
        if count!=int(merged.stdout):
            raise ArithmeticError('Merged prism bank row count differs from written count')
        pending.replace(bank)
    result=dict(colors=rows,core_orbits=count,nonzero_hard_orders=sum(n for c,n in selected),
                weights_sha256=signature)
    pending=folder/(word+'-weights.partial')
    pending.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    pending.replace(folder/(word+'-weights.json'))
    print('Complete prism numerator',word,len(rows),count,signature,flush=True)


if __name__=='__main__':
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('word')
    parser.add_argument('scratch',type=Path)
    args=parser.parse_args()
    with TemporaryDirectory() as tmp:
        if not os.environ.get('MIXED_SIMPLEX_GRAM_LIB'):
            lib=Path(tmp)/'prism-moments.so'
            run([*shlex.split(os.environ.get('CXX','c++')),'-O3','-std=c++17','-shared','-fPIC',
                 str(Path(__file__).with_name('simplex_gram_moments.cpp')),'-o',str(lib)],check=True)
            os.environ['MIXED_SIMPLEX_GRAM_LIB']=str(lib)
        project_placement(args.word,args.scratch)
