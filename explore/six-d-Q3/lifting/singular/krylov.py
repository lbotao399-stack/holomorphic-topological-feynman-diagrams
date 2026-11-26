"""Rational cohomology certificates from sparse Krylov lower bounds."""

import ctypes as ct
from collections import defaultdict
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory
from time import monotonic
import numpy as np

from balanced import balanced
from diagonal_products import count
from relative import Sector, map_row

PRIME = 65521


def build(folder):
    path = Path(folder) / "rank.so"
    subprocess.run(["c++","-std=c++17","-O3","-shared","-fPIC",str(Path(__file__).with_suffix(".cpp")),
                    "-o", str(path)], check=True)
    fn = ct.CDLL(str(path)).krylov_rank
    fn.argtypes = [ct.c_int,ct.c_int,np.ctypeslib.ndpointer(np.int32),
                   np.ctypeslib.ndpointer(np.int32),np.ctypeslib.ndpointer(np.uint32),ct.c_int,ct.c_uint]
    fn.restype = ct.c_int
    return fn


def check(weight, rank_fn):
    products = count(weight, by_length=True).get(weight,{})
    sec = Sector(weight)
    groups = defaultdict(list)
    anchors = defaultdict(dict)
    for letters,tensors in sec.singlet_blocks():
        length = len(letters)
        for free,tensor in tensors:
            mon = tuple(f+i for f,i in zip(letters,free))
            anchors[length][mon] = (len(anchors[length]), tensor[free])
        groups[length].append((letters,tensors))

    previous = 0
    for length in range(sum(weight)+1):
        required = len(anchors[length])-previous-products.get(length,0)
        if required==0:
            previous=0
            continue
        target = anchors[length+1]
        projection = {m:(i,pow(c,-1,PRIME)) for m,(i,c) in target.items()}
        offsets, columns, values = [0], [], []
        for letters, tensors in groups.pop(length,[]):
            cache = {}
            def field_image(colors):
                if colors not in cache:
                    mon = tuple(i+c for i,c in zip(letters,colors))
                    cache[colors] = {projection[m][0]:v*projection[m][1]
                                     for m,v in sec.differential(mon).items() if m in projection}
                return cache[colors]
            for _,tensor in tensors:
                for key,val in map_row(tensor,field_image).items():
                    residue = val % PRIME
                    if residue:
                        columns.append(key)
                        values.append(residue)
                offsets.append(len(columns))
        args = (len(offsets)-1,len(target),np.asarray(offsets,dtype=np.int32),
                np.asarray(columns,dtype=np.int32),np.asarray(values,dtype=np.uint32),required)
        certified=0
        for seed in (1237,7519,16993):
            certified=max(certified,rank_fn(*args,seed))
            if certified==required:
                break
        if certified != required:
            return False, length, required, certified
        previous=certified
    return True, dict(sorted(products.items()))


if __name__ == "__main__":
    import sys
    total = int(sys.argv[1])
    begin = int(sys.argv[2]) if len(sys.argv)>2 else 2
    end = int(sys.argv[3]) if len(sys.argv)>3 else total
    with TemporaryDirectory(prefix="twisted-ce-krylov-") as directory:
        rank_fn=build(directory)
        for h in range(begin,end+1):
            weight=balanced(h,3)+balanced(total-h,2)
            start=monotonic()
            print(weight,check(weight,rank_fn),round(monotonic()-start,2),flush=True)
