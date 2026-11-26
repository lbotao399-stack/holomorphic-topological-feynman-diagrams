"""Sum all one-form-on-C partitions of the single-momentum mixed detector."""

import argparse
import fcntl
import os
import signal
import subprocess
import sys
from itertools import combinations
from pathlib import Path

from single_projection import cases, determinant_table
from detector_node import single_fields

ROOT = Path(__file__).resolve().parent
STORED = ROOT / "single_partition_02345.tsv"


def partitions():
    # The determinant has one mixed column and one momentum-conserving gamma.
    # All partitions share the external fields at each of the centered nodes.
    nodes = [
        (br, columns[3], columns[9])
        for _, br, scales, switches, _ in cases()
        for _, columns in [single_fields(br, 0, 0, 0, scales, switches)]
    ]
    active = []
    for left in combinations(range(9), 5):
        order = left + tuple(i for i in range(10) if i not in left)
        tables = [determinant_table(order, branch) for branch in (0, 1)]
        count = sum(
            any(tables[br]) and sum(c * p[i] * q[j] for i, j, c in tables[br]) != 0
            for br, p, q in nodes
        )
        if count:
            active.append((count, "".join(map(str, order))))
    if len(active) != 75 or sum(n for n, _ in active) != 130328:
        raise ValueError("unexpected single-momentum determinant support")
    return sorted(active)


def read_ledger(lines, counts, prime):
    saved = {}
    for row in lines:
        if row.startswith("#"):
            if row.strip() != f"# mixed_single_weights.tsv 166 {prime}":
                raise ValueError("incompatible shape ledger")
            continue
        order, n, v = row.split()
        if order in saved or counts.get(order) != int(n) or not 0 <= int(v) < prime:
            raise ValueError("invalid completed partition")
        saved[order] = int(v)
    return saved


def run(args):
    tasks = partitions()
    counts = {order: n for n, order in tasks}
    ledger = args.directory / "partitions.tsv"
    if args.status:
        if not ledger.exists():
            print(0, "of", len(tasks), "partitions")
            return
        lines = ledger.read_text().splitlines()
        prime = int(lines[0].split()[-1])
        done = read_ledger(lines, counts, prime)
        pending = list(args.directory.glob("*.partial.tsv"))
        current = sum(max(0, len(p.read_text().splitlines()) - 1) for p in pending)
        print(
            len(done),
            "of",
            len(tasks),
            "partitions;",
            sum(counts[w] for w in done) + current,
            "of",
            sum(counts.values()),
            "nonzero nodes",
        )
        return
    solver = str(args.solver.resolve())
    prime = int(subprocess.check_output([solver, "--prime"], text=True))
    args.directory.mkdir(parents=True, exist_ok=True)
    with ledger.open("a+") as f:
        fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
        f.seek(0)
        existing = f.readlines()
        saved = read_ledger(existing, counts, prime)
        if not existing:
            print(f"# mixed_single_weights.tsv 166 {prime}", file=f, flush=True)
        for count, order in tasks:
            if order in saved:
                continue
            available = order == "0234516789" and prime == 1009 and STORED.is_file()
            partial = STORED if available else args.directory / (order + ".partial.tsv")
            cmd = [
                sys.executable,
                str(ROOT / "single_projection.py"),
                solver,
                "--order",
                order,
                "--checkpoint",
                str(partial),
            ]
            job = subprocess.Popen(
                cmd, stdout=subprocess.PIPE, text=True, start_new_session=True
            )
            try:
                stdout, _ = job.communicate()
                if job.returncode:
                    raise RuntimeError(f"partition {order} exited {job.returncode}")
            finally:
                if job.poll() is None:
                    os.killpg(job.pid, signal.SIGTERM)
                    try:
                        job.wait(timeout=8)
                    except subprocess.TimeoutExpired:
                        os.killpg(job.pid, signal.SIGKILL)
                        job.wait()
            words = stdout.split()
            if (
                words[:5] != ["partition", order, "surviving", "cases", str(count)]
                or words[5:7] != ["orders", "14400"]
                or words[-2:] != ["mod", str(prime)]
            ):
                raise ValueError(stdout)
            value = int(words[-3])
            if available and value != 302:
                raise ValueError("independent stored partition")
            saved[order] = value
            print(order, count, value, file=f, flush=True)
            os.fsync(f.fileno())
            if not available:
                partial.unlink()
            print(
                "partition",
                order,
                "residue",
                value,
                "completed",
                len(saved),
                "of",
                len(tasks),
                flush=True,
            )
        normalized = sum(saved.values()) * pow(2, -1, prime) % prime
        print(
            "Full normalized 166 with graph symmetry",
            normalized,
            "coefficient of pi",
            normalized * pow(158400, -1, prime) % prime,
            "mod",
            prime,
            flush=True,
        )


def interrupted(signum, frame):
    raise KeyboardInterrupt


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("solver", type=Path)
    p.add_argument("--directory", type=Path, default=ROOT / "single_shape_166")
    p.add_argument("--status", action="store_true")
    signal.signal(signal.SIGTERM, interrupted)
    run(p.parse_args())
