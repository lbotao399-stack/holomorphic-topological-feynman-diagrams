"""Certify all normalized four-wheel periods in the K3,3 contact order."""
from argparse import ArgumentParser
from contextlib import ExitStack
from fractions import Fraction
from hashlib import file_digest, sha256
from itertools import zip_longest
from math import gcd, isqrt, lcm, prod
from pathlib import Path
from subprocess import run
from tempfile import NamedTemporaryFile, TemporaryDirectory

from quartic_wheel_odd import denominator_bound, kernel, scaled_numerator_bound


def fields(q):
    primes = [1000000007, 999999937]
    modulus = prod(primes)
    candidate = primes[-1] - 2
    while modulus <= 2*scaled_numerator_bound(q):
        if all(candidate % i for i in range(3,isqrt(candidate)+1,2)):
            primes.append(candidate)
            modulus *= candidate
        candidate -= 2
    return primes


def field_table(bank, keys, cache, previous, q, prime):
    with keys.open("rb") as stream:
        key_hash = file_digest(stream,"sha256").hexdigest()
    name = f"quartic_wheel_values_{q}_{prime}_{key_hash[:16]}.tsv"
    for directory in (cache, previous):
        if directory is not None and (directory/name).is_file():
            return directory/name
    weights = cache/f"kernel_{q}_{prime}.tsv"
    with NamedTemporaryFile(dir=cache, mode="w",delete=False) as f:
        f.write("".join(" ".join(map(str,row))+"\n" for row in kernel(q,prime)))
        staging = Path(f.name)
    staging.replace(weights)
    with NamedTemporaryFile(dir=cache,delete=False) as f:
        pending = Path(f.name)
    try:
        run([str(bank),str(q),str(prime),str(weights),str(keys),str(pending)],check=True)
        pending.replace(cache/name)
    finally:
        pending.unlink(missing_ok=True)
    return cache/name


def reconstruct(q, primes, files, keys, output):
    denom = denominator_bound(q)
    limit = scaled_numerator_bound(q)
    modulus = prod(primes)
    assert modulus > 2*limit and gcd(denom,modulus) == 1
    factors = [(modulus//p) * ((pow(modulus//p,-1,p) * denom) % p) for p in primes]
    rows = max_bits = 0
    common = [1,1,1]
    digest = sha256()
    with ExitStack() as stack:
        temporary = Path(stack.enter_context(TemporaryDirectory(dir=output.parent))) / output.name
        out = stack.enter_context(temporary.open("wb"))
        inputs = [stack.enter_context(path.open("rb")) for path in (keys,*files)]
        for key, *lines in zip_longest(*inputs):
            if key is None or None in lines:
                raise ValueError("Wheel keys and field tables have different lengths")
            columns = [tuple(map(int,line.split())) for line in lines]
            if any(len(row) != 3 for row in columns):
                raise ValueError("Expected three channel entries")
            nums = []
            for i, residues in enumerate(zip(*columns)):
                v = sum(c*r for c,r in zip(factors,residues)) % modulus
                if 2*v > modulus:
                    v -= modulus
                if abs(v) > limit:
                    raise ValueError(f"Height bound exceeded at {int(key)} in column {i}")
                value = Fraction(v,denom)
                common[i] = lcm(common[i],value.denominator)
                max_bits = max(max_bits, abs(value.numerator).bit_length())
                nums.append(str(value))
            record = (str(int(key))+" "+" ".join(nums)+"\n").encode()
            out.write(record)
            digest.update(record)
            rows += 1
        if rows == 0:
            raise ValueError("Empty wheel bank")
        out.close()
        temporary.replace(output)
    print("certified faces",rows,"common denominators",*common,"maximum reduced numerator bits",max_bits,
          "sha256",digest.hexdigest(),flush=True)


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("q", type=int)
    parser.add_argument("bank", type=Path)
    parser.add_argument("keys", type=Path)
    parser.add_argument("cache", type=Path)
    parser.add_argument("--previous", type=Path)
    args = parser.parse_args()
    args.cache.mkdir(parents=True,exist_ok=True)
    primes = fields(args.q)
    print("fields",*primes,flush=True)
    files = []
    for p in primes:
        files.append(field_table(args.bank.resolve(),args.keys,args.cache,args.previous,args.q,p))
        print("field complete",p,flush=True)
    reconstruct(args.q,primes,files,args.keys,args.cache/"wheel_exact.tsv")
