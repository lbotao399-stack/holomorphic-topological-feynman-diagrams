"""Bounded, parameter-aware memoization for the unchanged scalar recurrence.

The original cut evaluator clears its reduction cache on *every* request.
Here identical requests at fixed q and Laurent precision reuse those exact
values.  A change of either parameter clears both recurrence caches.  The
context restores every patched binding and never modifies a repository file.
"""
from contextlib import contextmanager
from functools import lru_cache
from fractions import Fraction


@contextmanager
def fixed_parameter_cache(maxsize=20000):
    import laurent
    import prop

    original_set_q = prop.set_q
    original_reduce = prop.reduce_graph
    original_gamma_bubble = prop.G_value
    cached_reduce = lru_cache(maxsize=maxsize)(original_reduce.__wrapped__)
    cached_bubble = lru_cache(maxsize=maxsize)(original_gamma_bubble)
    parameter = None
    stats = {"parameter_changes": [], "maxsize": maxsize}

    def set_q(q):
        nonlocal parameter
        key = (Fraction(q), laurent.ORDER)
        if key != parameter:
            if parameter is not None and key[1] != parameter[1]:
                laurent.gamma_ratio.cache_clear()
            cached_reduce.cache_clear()
            cached_bubble.cache_clear()
            prop.Q = Fraction(q)
            parameter = key
            stats["parameter_changes"].append([str(key[0]), key[1]])

    prop.reduce_graph = cached_reduce
    prop.G_value = cached_bubble
    prop.set_q = set_q
    try:
        yield stats
    finally:
        stats["last_reduction_cache"] = cached_reduce.cache_info()._asdict()
        stats["last_bubble_cache"] = cached_bubble.cache_info()._asdict()
        cached_reduce.cache_clear()
        cached_bubble.cache_clear()
        prop.reduce_graph = original_reduce
        prop.G_value = original_gamma_bubble
        prop.set_q = original_set_q
        # The original cache may belong to a previous q.  Do not expose it
        # after the context with a different global dimension parameter.
        original_reduce.cache_clear()
