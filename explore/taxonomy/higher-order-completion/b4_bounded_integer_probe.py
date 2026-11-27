"""Corrected pooled driver with proved intermediate integer bounds.

Kernel integration is modulo the recorded prime; color contractions are
over the integers, retaining every original spectator color.
"""
from pathlib import Path
import runpy
import sys
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "qn-holomorphic"))
import safe_integer_contraction
import exact_colour

safe_integer_contraction.install(exact_colour)
original_colour = exact_colour.colour_values
controls = set()


def controlled_colour(*arguments, **keywords):
    actual = original_colour(*arguments, **keywords)
    bindings, original_inputs = arguments[2], arguments[4]
    spectator = bool(set(original_inputs) - set(bindings[0]) - set(bindings[1]))
    if spectator not in controls and actual.any():
        accelerated = exact_colour.np
        exact_colour.np = np
        try:
            expected = original_colour(*arguments, **keywords)
        finally:
            exact_colour.np = accelerated
        assert np.array_equal(actual, expected), "complete integer color tensor differs"
        controls.add(spectator)
        print("complete Python-integer tensor control", actual.shape,
              "spectator", spectator, "entries", actual.size, flush=True)
    return actual


exact_colour.colour_values = controlled_colour
runpy.run_path(str(Path(__file__).with_name("b4_modular_probe.py")), run_name="__main__")
print("certified integer bounds", safe_integer_contraction.STATISTICS, flush=True)
print("complete object-color controls", sorted(controls), flush=True)
