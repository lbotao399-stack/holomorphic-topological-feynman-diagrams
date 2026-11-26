"""Integer dot products of physical weights with the complete positive K3,3 scalar table."""
from array import array
from functools import cache
from hashlib import sha256
from math import lcm
from pathlib import Path
from tempfile import TemporaryDirectory
import fcntl
import gzip
import json
import os

import numpy as np
from flint import fmpq
from mixed_bipartite_star_bank import ARCHIVE,COMPLETION,verify_archive,verify_member


def ratio(s):
    n,sep,d=s.partition(b'/')
    return int(n),int(d) if sep else 1


def positive_rows():
    # The canonical 3×3 array is row and column sorted up to transpose, hence
    # its first entry is its minimum. A leading zero denotes a contraction.
    for source in (ARCHIVE,COMPLETION):
        with gzip.open(source,'rb') as f:
            for line in f:
                key,*values=line.split()
                k=int(key)
                if k>>32:yield k,[ratio(s) for s in values]


def file_hash(path):
    with path.open('rb') as f:
        h=sha256()
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
        return h.hexdigest()


def build(folder,source):
    meta=json.loads(Path(__file__).with_name('mixed_bipartite_star_bank.json').read_text())
    verify_archive()
    more=meta['completion']
    verify_member(COMPLETION,more['extra_positive']+more['extra_contracted'],more['sha256'])
    size=meta['positive_orbits']
    denom=[1,1,1];numer=[0,0,0];n=0
    for _,cols in positive_rows():
        n+=1
        for j,(v,d) in enumerate(cols):
            denom[j]=lcm(denom[j],d);numer[j]=max(numer[j],abs(v).bit_length())
    assert n==size
    widths=[(d.bit_length()+b+8)//8 for d,b in zip(denom,numer)]
    dtype=np.dtype([('key','<u8')]+[(f'c{j}',f'V{k}') for j,k in enumerate(widths)])
    with TemporaryDirectory(dir=folder) as tmp:
        scratch=Path(tmp)
        original=np.memmap(scratch/'unsorted',dtype=dtype,mode='w+',shape=size)
        for i,(k,cols) in enumerate(positive_rows()):
            original[i]=tuple([k]+[(v*(a//d)).to_bytes(w,'little',signed=True)
                                   for (v,d),a,w in zip(cols,denom,widths,strict=True)])
        permutation=np.argsort(original['key'],kind='stable')
        target=scratch/'positive.bin'
        keys=original['key'][permutation]
        assert np.all(keys[:-1]<keys[1:])
        with target.open('wb') as out:
            for start in range(0,size,32768):original[permutation[start:start+32768]].tofile(out)
        del original
        info=dict(source=source,rows=size,widths=widths,denominators=denom,sha256=file_hash(target))
        (scratch/'index.json').write_text(json.dumps(info,sort_keys=True)+'\n')
        target.replace(folder/'positive.bin')
        (scratch/'index.json').replace(folder/'index.json')
        return info


@cache
def open_table(folder):
    folder=Path(folder)
    folder.mkdir(exist_ok=True,parents=True)
    source_meta=json.loads(Path(__file__).with_name('mixed_bipartite_star_bank.json').read_text())
    signature=[source_meta['archive_sha256'],source_meta['completion']['sha256']]
    fd=os.open(folder,os.O_RDONLY)
    try:
        fcntl.flock(fd,fcntl.LOCK_EX)
        index=folder/'index.json'
        info=json.loads(index.read_text()) if index.exists() else {}
        path=folder/'positive.bin'
        if info.get('source')!=signature or not path.is_file() or file_hash(path)!=info.get('sha256'):
            info=build(folder,signature)
        dtype=np.dtype([('key','<u8')]+[(f'c{j}',f'V{k}') for j,k in enumerate(info['widths'])])
        rows=np.memmap(path,dtype=dtype,mode='r',shape=info['rows'])
    finally:os.close(fd)
    return rows,info['denominators']


def measure(folder,source):
    rows,denom=open_table(Path(folder).resolve())
    totals=[0,0,0];common=1
    with gzip.open(source,'rb') as f:
        while True:
            keys=array('Q');weights=[]
            for line in f:
                k,c=line.split();keys.append(int(k));weights.append(ratio(c))
                if len(keys)==16384:break
            if not keys:break
            lookup=np.frombuffer(keys,dtype=np.uint64)
            ids=np.searchsorted(rows['key'],lookup)
            if np.any(ids==len(rows)) or not np.array_equal(rows['key'][ids],lookup):
                raise ValueError('physical numerator includes a power outside the positive scalar domain')
            old=common
            for _,d in weights:common=lcm(common,d)
            for j in range(3):totals[j]*=common//old
            factors=[n*(common//d) for n,d in weights]
            selected=rows[ids]
            for j in range(3):
                totals[j]+=sum(c*int.from_bytes(v,'little',signed=True)
                               for c,v in zip(factors,selected[f'c{j}'],strict=True))
    return tuple(str(fmpq(c,common*d)) for c,d in zip(totals,denom,strict=True))


if __name__=='__main__':
    import sys
    for file in sys.argv[2:]:print(file,*measure(Path(sys.argv[1]),Path(file)),flush=True)
