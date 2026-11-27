"""Exploratory modular rank test with exact color and spectator retention.

This is not the rational calculation or a nilpotency proof over Q.
Arguments are those of the inherited poolnilp.py driver.
"""
from pathlib import Path
import runpy
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path[:0] = [str(ROOT / "explore/taxonomy/qn-holomorphic"), str(HERE)]
import fmode
import exact_colour


def bounded_integer_colour(*args, **kwargs):
    result = exact_colour.colour_values(*args, **kwargs)
    assert all(-(1 << 62) < value < (1 << 62) for value in result.flat)
    return np.array(result, dtype=np.int64)


fmode.colour_values = bounded_integer_colour
runpy.run_path(str(ROOT / "explore/taxonomy/b4-round2/poolnilp.py"), run_name="__main__")
print("exact-color statistics", exact_colour.STATISTICS, flush=True)
