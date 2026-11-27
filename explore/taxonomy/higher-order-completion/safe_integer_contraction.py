"""Bounded machine-integer acceleration of exact tensor contractions.

Each intermediate is checked *before* execution against an absolute sum
bound. If that bound reaches 2**62, Python integers perform that step.
There is no floating-point contraction or modular reconstruction.
"""
from math import prod
from numbers import Integral
import numpy as np

STATISTICS = {"bounded_integer_steps": 0, "python_integer_steps": 0,
              "largest_machine_step_bound": 0}


def maximum(array):
    if not array.size:
        return 0
    if array.dtype == np.dtype("int64"):
        return max(abs(int(array.min())), abs(int(array.max())))
    return max(abs(int(value)) for value in array.flat)


def einsum(specification, *operands, optimize=False, **kwargs):
    if kwargs or "..." in specification or "->" not in specification:
        raise ValueError("the certified contraction requires explicit integer indices")
    left, output = specification.split("->")
    labels = left.split(",")
    arrays = [np.asarray(array) for array in operands]
    dimensions = {}
    for indices, array in zip(labels, arrays):
        if array.dtype.kind not in "iu" and not (
                array.dtype.kind == "O" and all(isinstance(value, Integral) for value in array.flat)):
            raise TypeError("color contractions require exact integer operands")
        if len(indices) != array.ndim:
            raise ValueError("contraction rank mismatch")
        for index, size in zip(indices, array.shape):
            if index in dimensions and dimensions[index] != size:
                raise ValueError("broadcast dimensions are not part of this certificate")
            dimensions[index] = size
    if len(arrays) != len(labels) or len(set(output)) != len(output):
        raise ValueError("invalid explicit contraction")
    if isinstance(optimize, (list, tuple)) and optimize[0] == "einsum_path":
        path = optimize[1:]
    else:
        path = np.einsum_path(specification, *arrays, optimize="greedy")[0][1:]
    for positions in path:
        positions = tuple(positions)
        chosen_labels = [labels[i] for i in positions]
        chosen_arrays = [arrays[i] for i in positions]
        remaining_labels = [label for i, label in enumerate(labels) if i not in positions]
        remaining_arrays = [array for i, array in enumerate(arrays) if i not in positions]
        present = set("".join(chosen_labels))
        retained = present & set("".join(remaining_labels) + output)
        result_labels = output if not remaining_labels else "".join(sorted(retained))
        summed = present - retained
        maxima = [maximum(array) for array in chosen_arrays]
        bound = prod(maxima) * prod(dimensions[i] for i in summed)
        operation = ",".join(chosen_labels) + "->" + result_labels
        if bound < 1 << 62 and all(value < 1 << 62 for value in maxima):
            cast = [np.asarray(array, dtype=np.int64) for array in chosen_arrays]
            result = np.einsum(operation, *cast, optimize=False)
            STATISTICS["bounded_integer_steps"] += 1
            STATISTICS["largest_machine_step_bound"] = max(STATISTICS["largest_machine_step_bound"], bound)
        else:
            cast = [np.asarray(array, dtype=object) for array in chosen_arrays]
            result = np.einsum(operation, *cast, optimize=False)
            STATISTICS["python_integer_steps"] += 1
        arrays = remaining_arrays + [np.asarray(result)]
        labels = remaining_labels + [result_labels]
    if len(arrays) != 1 or labels != [output]:
        raise ArithmeticError("the complete contraction path was not executed")
    return arrays[0]


class NumpyWithCertifiedContractions:
    def __getattr__(self, name):
        return einsum if name == "einsum" else getattr(np, name)


def install(exact_colour):
    exact_colour.np = NumpyWithCertifiedContractions()
