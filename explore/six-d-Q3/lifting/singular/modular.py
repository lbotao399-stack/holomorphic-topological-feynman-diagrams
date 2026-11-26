"""Certified rational bounds via bit-parallel arithmetic mod the prime 127."""

from collections import Counter, defaultdict
from diagonal_products import count
from relative import Sector, map_row


class PackedSpan:
    """One base-2^15 lane per coordinate; 2^7=1 modulo 127."""
    def __init__(self, dim):
        self.prime = 127
        self.width = 15
        self.ones = ((1 << (self.width*dim))-1) // ((1 << self.width)-1)
        self.low = self.ones * self.prime
        self.high = self.ones << 7
        self.pivots = {}

    def add_multiple(self, a, multiple, b):
        # A lane reaches at most 126**2+126 < 2**14 before reduction.
        s = a + multiple*b
        s = (s & self.low) + ((s >> 7) & self.low)
        carry = ((s+self.ones) & self.high) >> 7
        return s - carry*self.prime

    def insert(self, row):
        bits = sum((v%self.prime) << (self.width*i) for i,v in row.items())
        while bits:
            last = bits.bit_length()-1
            shift = (last//self.width)*self.width
            val = (bits >> shift) & self.prime
            basis = self.pivots.get(shift)
            if basis is None:
                self.pivots[shift] = self.add_multiple(0,pow(val,-1,self.prime),bits)
                return True
            bits = self.add_multiple(bits,self.prime-val,basis)
        return False


def cohomology_upper_bound(weight):
    sec = Sector(weight)
    blocks = defaultdict(list)
    anchors = defaultdict(dict)
    for letters, tensors in sec.singlet_blocks():
        length = len(letters)
        for free, tensor in tensors:
            anchor = tuple(f+i for f,i in zip(letters, free))
            anchors[length][anchor] = (len(anchors[length]),tensor[free])
        blocks[length].append((letters,tensors))
    counts = {}
    boundary_pivots = set()
    for length in range(sum(weight)+1):
        target = anchors[length+1]
        mod = PackedSpan(len(target))
        proj = {m:(i,pow(n,-1,mod.prime)) for m,(i,n) in target.items()}
        sparse_rows = []
        source_number = 0
        for letters,tensors in blocks.pop(length,[]):
            cache = {}
            def image(colors):
                if colors not in cache:
                    mon = tuple(f+i for f,i in zip(letters,colors))
                    cache[colors] = {proj[m][0]: n*proj[m][1]
                                     for m,n in sec.differential(mon).items() if m in proj}
                return cache[colors]
            for _,tensor in tensors:
                if source_number not in boundary_pivots:
                    row = {k:v%mod.prime for k,v in map_row(tensor,image).items() if v%mod.prime}
                    if row:
                        sparse_rows.append(row)
                source_number += 1
        # Q²=0 permits deletion of a basis for the previous boundary
        # space from the current domain. Put the sparsest target coordinates
        # at the pivot end of the packed vector to limit fill-in.
        frequencies = Counter(k for row in sparse_rows for k in row)
        target_order = sorted(frequencies, key=frequencies.get, reverse=True)
        position = {old:new for new,old in enumerate(target_order)}
        for row in sorted(sparse_rows,key=len):
            mod.insert({position[k]:v for k,v in row.items()})
        rank = len(mod.pivots)
        dim = len(anchors[length])-rank-len(boundary_pivots)
        if dim:
            counts[length] = dim
        boundary_pivots = {target_order[shift//mod.width] for shift in mod.pivots}
    return counts


if __name__ == '__main__':
    from balanced import balanced
    from time import monotonic
    import sys
    degree = int(sys.argv[1])
    begin = int(sys.argv[2]) if len(sys.argv)>2 else 2
    end = int(sys.argv[3]) if len(sys.argv)>3 else degree
    for h in range(begin,end+1):
        weight = balanced(h,3)+balanced(degree-h,2)
        start = monotonic()
        upper = cohomology_upper_bound(weight)
        gen = count(weight).get(weight,0)
        assert gen <= sum(upper.values()), (weight, upper, gen)
        print(weight, upper, "products",gen,"gap",sum(upper.values())-gen,
              round(monotonic()-start,2), flush=True)
