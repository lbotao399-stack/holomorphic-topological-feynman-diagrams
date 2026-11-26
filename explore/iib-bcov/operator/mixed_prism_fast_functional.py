"""Integrate weighted prism powers with the complete positive q=13 scalar table."""
from argparse import ArgumentParser
from array import array
from functools import cache
from hashlib import sha256
from itertools import permutations
from math import factorial,lcm,prod
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory
import fcntl
import gzip
import json
import os

import numpy as np
from flint import fmpq
from mixed_bipartite_scalar_functional import file_hash,ratio
from mixed_prism_scalar_functional import ARCHIVE,METADATA,FACE_ORDER,batches,compile_program,pack,unpack


FINITE=np.dtype([('key','<u8'),('angular','<u4')])


def angular_dtype(widths):
    return np.dtype([('key','<u8')]+[(f'c{j}',f'V{n}') for j,n in enumerate(widths)])


def build(folder,source):
    shape=source['angular_classes']
    common=[1]*3
    bits=[0]*3
    checksum=sha256()
    with gzip.open(ARCHIVE,'rb') as stream:
        for line in stream:
            checksum.update(line)
            for j,c in enumerate(line.split()[1:]):
                n,d=ratio(c);common[j]=lcm(common[j],d);bits[j]=max(bits[j],abs(n).bit_length())
    if checksum.hexdigest()!=source['uncompressed_sha256']:raise ArithmeticError('Angular archive checksum mismatch')
    widths=[(b+d.bit_length()+8)//8 for b,d in zip(bits,common,strict=True)]
    dtype=angular_dtype(widths)
    with TemporaryDirectory(dir=folder) as tmp:
        scratch=Path(tmp);angular=scratch/'angular.bin';finite=scratch/'finite.bin'
        count=0
        with angular.open('wb') as out,gzip.open(ARCHIVE,'rb') as stream:
            for lines in batches(stream,16384):
                block=np.empty(len(lines),dtype=dtype)
                for i,line in enumerate(lines):
                    key,*terms=line.split();cols=[]
                    for c,d,width in zip(terms,common,widths,strict=True):
                        n,v=ratio(c)
                        cols.append((n*(d//v)).to_bytes(width,'little',signed=True))
                    block[i]=(int(key),*cols)
                block.tofile(out);count+=len(block)
        if count!=shape:raise ArithmeticError('Incomplete angular archive')
        angular_rows=np.memmap(angular,mode='r',dtype=dtype)
        keys=np.ascontiguousarray(angular_rows['key'])
        if not np.all(keys[:-1]<keys[1:]):raise ArithmeticError('Angular powers out of order')
        compile_program('operator/mixed_prism_positive_domain.cpp',scratch/'domain')
        pairs=scratch/'pairs.bin'
        with pairs.open('wb') as out:
            run([scratch/'domain','13','--finite-aliases'],stdout=out,check=True)
        raw=np.memmap(pairs,mode='r',dtype='<u8').reshape((-1,2))
        if len(raw)!=source['finite_prism_orbits']:raise ArithmeticError('Incomplete finite prism census')
        order=np.argsort(raw[:,0],kind='stable')
        previous=0
        with finite.open('wb') as out:
            for first in range(0,len(order),65536):
                rows=raw[order[first:first+65536]]
                ids=np.searchsorted(keys,rows[:,1])
                if np.any(ids>=shape) or not np.array_equal(keys[ids],rows[:,1]):raise ArithmeticError('Missing angular orbit')
                if rows[0,0]<=previous or not np.all(rows[:-1,0]<rows[1:,0]):raise ArithmeticError('Repeated finite prism power')
                previous=rows[-1,0]
                block=np.empty(len(rows),dtype=FINITE)
                block['key']=rows[:,0];block['angular']=ids
                block.tofile(out)
        info={'source':source['uncompressed_sha256'],'denominators':common,'widths':widths,
              'angular_sha256':file_hash(angular),'finite_sha256':file_hash(finite)}
        manifest=scratch/'index.json';manifest.write_text(json.dumps(info,sort_keys=True)+'\n')
        angular.replace(folder/angular.name);finite.replace(folder/finite.name);manifest.replace(folder/manifest.name)
    return info


@cache
def open_index(folder):
    folder=Path(folder).resolve();folder.mkdir(exist_ok=True,parents=True)
    source=json.loads(METADATA.read_text())
    lock=os.open(folder,os.O_RDONLY)
    try:
        fcntl.flock(lock,fcntl.LOCK_EX)
        path=folder/'index.json'
        info=json.loads(path.read_text()) if path.is_file() else {}
        ang,fin=folder/'angular.bin',folder/'finite.bin'
        if (info.get('source')!=source['uncompressed_sha256'] or not ang.is_file() or not fin.is_file()
                or file_hash(ang)!=info.get('angular_sha256') or file_hash(fin)!=info.get('finite_sha256')):
            info=build(folder,source)
        angular=np.memmap(ang,mode='r',dtype=angular_dtype(info['widths']),shape=source['angular_classes'])
        finite=np.memmap(fin,mode='r',dtype=FINITE,shape=source['finite_prism_orbits'])
    finally:os.close(lock)
    return angular,finite,info['denominators']


@cache
def finite_permutations():
    groups=[]
    for triangle in ((0,1),(1,0)):
        for rectangular in permutations((2,3,4)):
            move=triangle+rectangular
            groups.append(tuple(FACE_ORDER.index(tuple(sorted((move[i],move[j])))) for i,j in FACE_ORDER))
    return tuple(groups)


def select(packed,finite):
    field=finite['key']
    found=np.searchsorted(field,packed)
    def missing(ids,inputs):
        return (ids==len(field)) | (field[np.minimum(ids,len(field)-1)]!=inputs)
    mask=missing(found,packed)
    if np.any(mask):
        shifts=np.arange(32,-1,-4,dtype=np.uint64)
        powers=(packed[mask,None]>>shifts)&15
        canonical=packed[mask].copy()
        for g in finite_permutations():
            canonical=np.minimum(canonical,np.sum(powers[:,g]<<shifts,axis=1))
        indices=np.searchsorted(field,canonical)
        if np.any(missing(indices,canonical)):raise ValueError('Power outside the positive convergent prism domain')
        found[mask]=indices
    return finite['angular'][found]


def integer_columns(folder,keys):
    angular,finite,denominators=open_index(Path(folder).resolve())
    ids=select(np.asarray(keys,dtype=np.uint64),finite)
    distinct,positions=np.unique(ids,return_inverse=True)
    rows=angular[distinct]
    aliases=rows['key'].tolist()
    periods={key:tuple(int.from_bytes(row[f'c{j}'],'little',signed=True) for j in range(3))
             for key,row in zip(aliases,rows,strict=True)}
    return [aliases[i] for i in positions],periods,denominators


def measure(folder,path,raw=False):
    angular,finite,denominators=open_index(Path(folder).resolve())
    totals=[0]*3;common=1
    opener=gzip.open if path.suffix=='.gz' else open
    with opener(path,'rb') as stream:
        for lines in batches(stream,32768):
            keys=array('Q');weights=[]
            for line in lines:
                terms=line.split()
                if len(terms)==2:
                    key=int(terms[0])
                    if not 0<=key<1<<36:raise ValueError('Expected nine packed prism powers')
                elif len(terms)==10:
                    powers=tuple(map(int,terms[:9]))
                    if min(powers)<1 or max(powers)>12:raise ValueError('Invalid prism power')
                    key=pack(powers)
                else:raise ValueError('Expected a packed angular weight or nine powers and a weight')
                keys.append(key)
                n,d=ratio(terms[-1])
                if raw:
                    powers=unpack(key,9)
                    if min(powers)<1:raise ValueError('Zero prism power')
                    n*=prod(factorial(p-1) for p in powers);d*=factorial(12)
                weights.append((n,d))
            ids=select(np.frombuffer(keys,dtype=np.uint64),finite)
            old=common
            for _,d in weights:common=lcm(common,d)
            factors=[n*(common//d) for n,d in weights]
            selected=angular[ids]
            for j in range(3):
                totals[j]=totals[j]*(common//old)+sum(c*int.from_bytes(x,'little',signed=True)
                    for c,x in zip(factors,selected[f'c{j}'],strict=True))
    return tuple(str(fmpq(n,common*d)) for n,d in zip(totals,denominators,strict=True))


if __name__=='__main__':
    p=ArgumentParser(description=__doc__)
    p.add_argument('cache',type=Path)
    p.add_argument('banks',nargs='+',type=Path)
    p.add_argument('--raw',action='store_true')
    args=p.parse_args()
    for path in args.banks:print('Integrated prism',path,*measure(args.cache,path,args.raw),flush=True)
