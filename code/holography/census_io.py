"""Read-only access to the LLC census data and to the Chang-Lin appendix tables transcribed in the
census code.  Paths are fixed; files are opened read-only."""
import ast
import csv
from fractions import Fraction

REPO = "/Users/libotao/Desktop/研究项目_非Git/All_loop_arxiv_2026-11-10/work/repo"
CODE = "/Users/libotao/Desktop/研究项目_非Git/BPS枚举/python版"


def su2_levels():
    """data/N4_SU2_levels.csv -> {j: dict}."""
    out = {}
    with open(f"{REPO}/data/N4_SU2_levels.csv", newline="") as fh:
        for row in csv.DictReader(fh):
            j = int(row["j"])
            out[j] = {k: (Fraction(v) if "/" in v else int(v)) for k, v in row.items()}
    return out


def su2_nongraviton_blocks():
    rows = []
    with open(f"{REPO}/data/N4_SU2_nongraviton_blocks.csv", newline="") as fh:
        for row in csv.DictReader(fh):
            rows.append(row)
    return rows


def chang_lin_tables():
    """reference_tables.py (transcription of the appendix Tables 3-5 of Chang-Lin): parsed with ast,
    not executed.  Returns {rank: {field: list}}."""
    src = open(f"{CODE}/reference_tables.py", encoding="utf-8").read()
    tree = ast.parse(src)
    out = {}
    for node in tree.body:
        if isinstance(node, ast.Assign) and isinstance(node.value, ast.Call):
            name = node.targets[0].id
            if not name.startswith("TABLE"):
                continue
            kw = {k.arg: ast.literal_eval(k.value) for k in node.value.keywords}
            out[kw["gauge_rank"]] = kw
    return out
