"""Rebuild all five alternating theta-functionals in scratch, with exact and modular checks."""
from argparse import ArgumentParser
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from contextlib import contextmanager
from threading import Lock
from fractions import Fraction
from pathlib import Path
from tempfile import TemporaryDirectory
import os
import re
import shlex
import shutil
import signal
import subprocess
import time

SOURCES = Path(__file__).resolve().parent
# Shape, number of signed coordinate matchings, automorphism order, raw kernel / pi.
CASES = {
    (1, 6, 6): (120, 4, Fraction(-12644, 892812375)),
    (2, 5, 6): (120, 2, Fraction(-278323, 44895708000)),
    (3, 4, 6): (120, 2, Fraction(-4141, 1662804000)),
    (3, 5, 5): (192, 4, Fraction(-519697, 314269956000)),
    (4, 4, 5): (216, 4, Fraction(-35669, 62853991200)),
}
DEFAULT_AXES = "2,3,4,5,2,2,1,1,1,1"
# Raw coefficients divided by pi, in CASES order, for this ordered basis fixture.
SECTOR_VALUES = {
    1: tuple(case[2] for case in CASES.values()),
    0x130: tuple(map(Fraction, (
        '0', '0', '0',
        '0', '0',
    ))),
    0x174: tuple(map(Fraction, (
        '-7593791/114279984000', '-899281/27935107200', '-142790357/7542478944000',
        '-718847/46558512000', '-44809553/3771239472000',
    ))),
    0xdf: tuple(map(Fraction, (
        '-644333/6722352000', '-10918933/209513304000', '-79947643/1885619736000',
        '-60008659/1885619736000', '-1402801/44895708000',
    ))),
    0x3df: tuple(map(Fraction, (
        '-491/3007368000', '-2861/62853991200', '-10627/36972936000',
        '2113/34918884000', '9473/12324312000',
    ))),
}
CUSTOM_VALUES = {
    (0x250,"4,3,4,5,4,2,1,1,1,2"): tuple(map(Fraction, (
        "-68432557/2628439632000", "-437149/34358688000", "-507388727/57825671904000",
        "-9979231/1807052247000", "-379363129/86738507856000",
    ))),
    (0x3fb,"5,3,4,5,3,2,1,3,4,1"): tuple(map(Fraction, (
        "4020047/28569996000", "96368917/1257079824000", "115393697/1885619736000",
        "86067409/1885619736000", "9109651/209513304000",
    ))),
}
CHUNK = 4


class Commands:
    """One replay's cancellation state. Each caller owns and reaps its child group."""
    def __init__(self):
        self.lock = Lock()
        self.failure = None
        self.signal_number = None

    def check(self):
        if self.signal_number is not None:
            raise SystemExit(128+self.signal_number)
        if self.failure is not None:
            raise self.failure

    def fail(self, exc):
        with self.lock:
            if self.failure is None:
                self.failure = exc

    @contextmanager
    def signals(self):
        def request_stop(number, frame):
            self.signal_number = number
        old = {s: signal.signal(s, request_stop) for s in (signal.SIGINT, signal.SIGTERM)}
        try:
            yield
        finally:
            for number, handler in old.items():
                signal.signal(number, handler)

    @staticmethod
    def clean_group(child):
        # start_new_session makes the child's PID its process-group ID. A compiler's
        # helpers may outlive its leader or hold its output pipes open after an error.
        def send(number):
            try:
                os.killpg(child.pid, number)
                return True
            except ProcessLookupError:
                return False
        if send(signal.SIGTERM):
            deadline = time.monotonic()+1.0
            while time.monotonic()<deadline:
                child.poll()
                if not send(0):
                    break
                time.sleep(0.03)
            else:
                send(signal.SIGKILL)
        child.wait(timeout=3)
        child.stdout.close()
        child.stderr.close()

    def run(self, args, *, optional=False):
        with self.lock:
            self.check()
            child = subprocess.Popen(args, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                     start_new_session=True)
        try:
            while True:
                self.check()
                try:
                    out, err = child.communicate(timeout=0.15)
                    break
                except subprocess.TimeoutExpired:
                    # A failed leader may leave an inherited pipe open in a helper.
                    if child.poll() not in (None, 0):
                        if optional:
                            return None
                        raise RuntimeError(f"{args}: exit {child.returncode}")
            self.check()
            if child.returncode:
                if optional:
                    return None
                raise RuntimeError(f"{args}: exit {child.returncode}\n{err}\n{out}")
            return out
        except Exception as exc:
            self.fail(exc)
            raise
        finally:
            self.clean_group(child)


def gmp_flags(commands):
    if shutil.which("pkg-config"):
        out = commands.run(["pkg-config", "--cflags", "--libs", "gmpxx"], optional=True)
        if out is not None:
            return shlex.split(out)
    if shutil.which("brew"):
        out = commands.run(["brew", "--prefix", "gmp"], optional=True)
        if out is not None:
            p = Path(out.strip())
            return [f"-I{p / 'include'}", f"-L{p / 'lib'}", "-lgmpxx", "-lgmp"]
    return ["-lgmpxx", "-lgmp"]


def build(scratch, commands):
    flags = gmp_flags(commands)
    compiler = shlex.split(os.environ.get("CXX", "c++"))
    for name in ("assemble", "moments", "independent", "check_physical"):
        args = [*compiler, "-O3", "-std=c++17", SOURCES / (name+".cpp")]
        if name in ("moments", "independent"):
            args.extend(flags)
        commands.run([*args, "-o", scratch / name])


def pairings(available):
    if not available:
        yield ()
        return
    first = available[0]
    for i in range(1, len(available)):
        second = available[i]
        for tail in pairings(available[1:i]+available[i+1:]):
            yield ((first, second),)+tail


def sign(entries):
    inversions = sum(a > b for i, a in enumerate(entries) for b in entries[i+1:])
    return (-1)**inversions


def check_collisions(shape):
    # Group determinant rows by coordinate and columns by paired slots.
    # Each pair's second momentum is negative in source order.
    row_order = [r for axis in range(5) for r in (axis,axis+5)]
    path = [p for p,edges in enumerate(shape) for _ in range(edges-1)]
    column = ((1,0), (1,1), (0,1))
    negative = positive = collision = 0
    for pairs in pairings(tuple(range(10))):
        col_order = [v for pair in pairs for v in pair]
        determinant = sign(row_order)*sign(col_order)
        axes = [0]*10
        for axis, (i,j) in enumerate(pairs):
            axes[i] = axes[j] = axis
            a,b = column[path[i]]
            c,d = column[path[j]]
            determinant *= b*c-a*d
        hit = any(path[i] == path[j] for i,j in pairs)
        if hit:
            if determinant != 0:
                raise ValueError("same-path momentum collision has nonzero determinant")
            collision += 1
        else:
            if determinant != -sign(axes):
                raise ValueError("oriented determinant identity failed")
            negative += determinant < 0
            positive += determinant > 0
    if negative+positive != CASES[shape][0] or negative+positive+collision != 945:
        raise ValueError("coordinate matching counts failed")
    return negative,positive,collision


def calculate_chunk(scratch,shape,first,last,seeds,commands,gamma_mask=1,polarizations=None):
    run=commands.run
    prefix = ''.join(map(str,shape))
    args = list(map(str,shape))
    path = scratch / f"{prefix}_{first}-{last}.poly"
    mask_args = [] if gamma_mask == 1 else ["--gamma-mask", str(gamma_mask)]
    if polarizations is not None:
        mask_args += ["--polarizations", polarizations]
    run([scratch/"assemble",*args,str(first),str(last),path,*mask_args])
    a = run([scratch/"moments",*args,path]).splitlines()
    b = run([scratch/"independent",*args,path]).splitlines()
    n = last-first
    if len(a)!=n+1 or len(b)!=n+2:
        raise ValueError(f"{path.name}: incomplete integral output")
    total = Fraction(0)
    for index in range(n):
        left,right = a[index].split(),b[index].split()
        if (int(left[1]),int(right[0]))!=(first+index,first+index) or int(left[3])!=int(right[1]):
            raise ValueError(f"{path.name}: duplicate, missing, or inconsistent rank")
        primary,independent = Fraction(left[-1]),Fraction(right[-1])
        if primary!=independent:
            raise ValueError(f"{path.name}: exact disagreement at rank {first+index}")
        total += primary
    if total!=Fraction(a[-1].split()[-1]) or total!=Fraction(b[-1].split()[-1]):
        raise ValueError(f"{path.name}: inconsistent integral total")
    controls=[]
    for seed in seeds:
        out=run([scratch/"check_physical",*args,str(seed),path,*mask_args]).strip()
        match=re.fullmatch(r"\S+ (\d+) blocks x 3840 physical assignments x 4 flows; F1009 seed (\d+); (\d+) nonzero, no exterior sign fails (\d+)",out)
        if match is None or tuple(map(int,match.groups()[:2]))!=(n,seed):
            raise ValueError(f"{path.name}: incomplete physical-coordinate check: {out}")
        controls.append(tuple(map(int,match.groups()[2:])))
    path.unlink()
    return shape,first,last,total,controls


def replay(scratch,jobs,seeds,commands,gamma_mask=1,polarizations=None):
    build(scratch,commands)
    for shape in CASES:
        commands.check()
        neg,pos,zero=check_collisions(shape)
        print(f"{shape}: determinant −1,+1,0 counts {neg},{pos},{zero}",flush=True)
    bounds=[(s,i,min(i+CHUNK,n)) for i in range(0,max(c[0] for c in CASES.values()),CHUNK)
            for s,(n,_,_) in CASES.items() if i<n]
    totals={s:Fraction(0) for s in CASES}
    ranks={s:set() for s in CASES}
    controls={seed:[0,0] for seed in seeds}
    def task(bound):
        try:
            commands.check()
            return calculate_chunk(scratch,*bound,seeds,commands,gamma_mask,polarizations)
        except Exception as exc:
            commands.fail(exc)
            raise

    pending=set()
    todo=iter(bounds)
    completed=0
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        try:
            while True:
                commands.check()
                while len(pending)<jobs:
                    bound=next(todo,None)
                    if bound is None:
                        break
                    commands.check()
                    pending.add(pool.submit(task,bound))
                if not pending:
                    break
                ready,pending=wait(pending,timeout=0.15,return_when=FIRST_COMPLETED)
                for future in ready:
                    s,i,j,total,local=future.result()
                    chunk_ranks=set(range(i,j))
                    if ranks[s]&chunk_ranks:
                        raise ValueError("duplicate chunk")
                    ranks[s].update(chunk_ranks)
                    totals[s]+=total
                    for seed,counts in zip(seeds,local):
                        for k in range(2):
                            controls[seed][k]+=counts[k]
                    completed+=1
                    if completed%24==0 or completed==len(bounds):
                        print(f"{completed}/{len(bounds)} chunks checked",flush=True)
        except BaseException as exc:
            commands.fail(exc)
            for future in pending:
                future.cancel()
            raise
    weighted=Fraction(0)
    reference = (SECTOR_VALUES.get(gamma_mask) if polarizations in (None,DEFAULT_AXES)
                 else CUSTOM_VALUES.get((gamma_mask,polarizations)))
    for i,(shape,(count,automorphisms,_)) in enumerate(CASES.items()):
        if ranks[shape]!=set(range(count)):
            raise ValueError(f"{shape}: missing matching ranks")
        if reference is not None and totals[shape]!=reference[i]:
            raise ValueError(f"{shape}: got {totals[shape]}, expected {reference[i]}")
        weighted+=totals[shape]/automorphisms
        print(f"{shape} raw / pi = {totals[shape]}; symmetry-weighted / pi = {totals[shape]/automorphisms}")
    for seed,(nonzero,rejected) in controls.items():
        print(f"F1009 seed {seed}: all 768 signed blocks agree; {nonzero} nonzero, {rejected} change without the exterior sign.")
    if gamma_mask == 1 and polarizations in (None,DEFAULT_AXES) and weighted!=Fraction(-73,8648640):
        raise ValueError(f"five-graph sum mismatch: {weighted}")
    print(f"Five-graph symmetry-weighted sum / pi = {weighted}")


if __name__=="__main__":
    parser=ArgumentParser(description=__doc__)
    parser.add_argument("--jobs",type=int,choices=range(1,7),default=min(4,os.cpu_count() or 1),
                        help="maximum simultaneous chunks, 1–6 (default at most 4)")
    parser.add_argument("--seeds",type=int,nargs="+",default=[4],
                        help="deterministic physical-coordinate finite-field points (default 4)")
    parser.add_argument("--gamma-mask", type=lambda value: int(value, 0), default=1,
                        help="Gamma inputs as a 10-bit mask in source wave order (default 1).")
    parser.add_argument("--polarizations", help="Ten coordinate axes, each 1–5, comma-separated in source wave order.")
    options=parser.parse_args()
    if options.polarizations is not None and not re.fullmatch(r"[1-5](,[1-5]){9}", options.polarizations):
        parser.error("polarizations must be ten comma-separated digits from 1 through 5")
    if not 0 < options.gamma_mask < 1024 or options.gamma_mask.bit_count() % 2 != 1:
        parser.error("gamma mask must select an odd number of the ten external waves")
    if len(set(options.seeds))!=len(options.seeds) or any(s<0 or s>1000 for s in options.seeds):
        parser.error("seeds must be distinct integers from 0 through 1000")
    commands=Commands()
    with commands.signals(), TemporaryDirectory(prefix="five-theta-") as tmp:
        replay(Path(tmp),options.jobs,tuple(options.seeds),commands,options.gamma_mask,options.polarizations)
        commands.check()
