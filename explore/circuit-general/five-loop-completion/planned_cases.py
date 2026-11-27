"""Finish budget-limited bare cases with exact alternative trees or signed zeros."""
import argparse
from collections import deque
import hashlib
import json
from pathlib import Path
import signal
import time
from types import SimpleNamespace

from bulk3_periods import Bulk3Periods, HERE, ROOT, current_minor, encoded, polynomial_record
from certify import ROWS, atomic_json, hashes, selected, zero_witness
from htk.graph import Graph
from htk.treeflow import tree_flow


def candidate_trees(graph, sources):
    lookup = {frozenset(edge): i for i, edge in enumerate(graph.edges)}
    adjacency = {v: [] for v in range(graph.nv)}
    for edge, (left, right) in enumerate(graph.edges):
        adjacency[left].append(right)
        adjacency[right].append(left)
    candidates = {tuple(graph.spanning_tree())}
    for root in range(graph.nv):
        seen, queue, tree = {root}, deque([root]), []
        while queue:
            vertex = queue.popleft()
            for neighbor in sorted(adjacency[vertex]):
                if neighbor not in seen:
                    seen.add(neighbor)
                    queue.append(neighbor)
                    tree.append(lookup[frozenset((vertex, neighbor))])
        candidates.add(tuple(sorted(tree)))
    ranked = []
    for tree in sorted(candidates):
        flow = [[int(value) for value in row] for row in tree_flow(graph, sources, tree)]
        proxy = SimpleNamespace(n=graph.loops, N=graph.ne, C=graph.cycle_matrix(tree))
        score = []
        def capture(matrix, context):
            score.extend((sum(len(value) for row in matrix for value in row),
                          sum(bool(value) for row in matrix for value in row)))
            return context.constant(0)
        current_minor.numerator(proxy, sources, flow, determinant=capture)
        ranked.append((tuple(score), tree))
    return sorted(ranked)


def complete(index, seconds, currents=False):
    path = HERE / "cases" / ("bulk3_%03d.json" % index)
    result = json.loads(path.read_text())
    if result["dependencies_sha256"] != hashes():
        raise ArithmeticError("The original bounded case has changed dependencies.")
    if result["bare_status"] == "complete" and (not currents or result["first_current_status"] == "complete"):
        print("Already complete", index, flush=True)
        return
    row = ROWS[index]
    witness = zero_witness(row)
    if witness and not currents:
        certificate_path = HERE / "symmetry_certificate.json"
        certificate = json.loads(certificate_path.read_text())
        if index not in [entry["index"] for entry in certificate["complete_bare_zeros"]]:
            raise ArithmeticError("The whole source-polynomial zero is not certified independently.")
        result.update(method="complete_signed_automorphism", bare_status="complete", bare_coefficient={},
                      independent_zero_automorphism=witness, complete_finite_targets=[],
                      independent_signed_zero_certificate_sha256=hashlib.sha256(certificate_path.read_bytes()).hexdigest())
        result.pop("error", None)
        result["planning_dependencies_sha256"] = {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (Path(__file__), certificate_path, HERE / "symmetry_controls.py")}
        atomic_json(path, result)
        print("Completed signed zero", index, flush=True)
        return
    graph = Graph(row["nv"], row["edges"])
    sources = [[0, 0, 0] for _ in range(graph.nv)]
    for i, vertex in enumerate(row["S0"][:3]):
        sources[vertex][i] = 1
    sources[row["S0"][-1]] = [-1, -1, -1]
    plans = candidate_trees(graph, sources)
    attempts = []
    for score, tree in plans[:4]:
        engine = None
        started = time.monotonic()
        def expired(*unused):
            raise TimeoutError("Alternative complete plan exceeded the declared time budget.")
        signal.signal(signal.SIGALRM, expired)
        signal.alarm(seconds)
        try:
            engine = Bulk3Periods(row, tree=tree)
            bare = engine.integrate()
            if witness and bare:
                raise ArithmeticError("The complete bare value violates the independent signed-zero theorem.")
            values = None
            if currents:
                values = engine.complete_first_currents(
                    bare, progress=lambda edge: print("Planned case", index, "whole chord", edge, "done", flush=True))
            result.update(bare_status="complete", bare_coefficient=encoded(bare),
                          source_cycle_tree=list(tree), method="complete_planned_linear_minor",
                          complete_bare_numerator=polynomial_record(engine._bare),
                          original_grounded_laplacian_density_checks=engine.density_checks,
                          complete_series_chains=[list(chain.edges) for chain in engine._chains],
                          scalar_core_edges=[list(edge) for edge in engine._core],
                          complete_finite_targets=engine.integrations,
                          complete_tree_planning={"considered_trees": len(plans), "selected_score": score,
                                                  "failed_complete_plans": attempts},
                          planning_dependencies_sha256={
                              str(Path(__file__).relative_to(ROOT)): hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
            if currents:
                result.update(first_current_status="complete",
                              all_original_first_current_coefficients=[
                                  [encoded(value) for value in item] for item in values],
                              integrated_original_source_identities=3 * row["nv"],
                              complete_affine_current_polynomial_identities=3 * len(row["edges"]))
            result.pop("error", None)
            atomic_json(path, result)
            print("Completed", index, "terms", len(engine._bare), "seconds", time.monotonic() - started, flush=True)
            return
        except Exception as error:
            attempts.append({"tree": list(tree), "score": score,
                             "failure": type(error).__name__ + ": " + str(error)})
            print("Plan failed", index, attempts[-1], flush=True)
        finally:
            signal.alarm(0)
            if engine is not None:
                engine.close()
    raise ArithmeticError("All bounded complete plans failed for graph %d." % index)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--indices", default="19,25,39")
    parser.add_argument("--seconds", type=int, default=180)
    parser.add_argument("--currents", action="store_true", help="finish all original first-current vectors too")
    args = parser.parse_args()
    for index in selected(args.indices):
        complete(index, args.seconds, args.currents)
