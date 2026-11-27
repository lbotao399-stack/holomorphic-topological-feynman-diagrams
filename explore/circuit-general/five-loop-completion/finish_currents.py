"""Finish whole original-current tensors by exact fibres and complete relations."""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import signal
import time

from bulk3_periods import Bulk3Periods, HERE, ROOT, encoded, polynomial_record
from certify import ROWS, atomic_json, hashes, selected
from chain_currents import ChainCurrentPeriods
from current_relations import CurrentRelations
from planned_cases import candidate_trees
from htk.graph import Graph


def complete(index, seconds):
    path = HERE / "cases" / ("bulk3_%03d.json" % index)
    result = json.loads(path.read_text())
    if result["dependencies_sha256"] != hashes():
        raise ArithmeticError("The existing case has changed exact dependencies.")
    if result["first_current_status"] == "complete":
        print("Already complete", index, flush=True)
        return
    row = ROWS[index]
    graph = Graph(row["nv"], row["edges"])
    sources = [[0, 0, 0] for _ in range(graph.nv)]
    for i, vertex in enumerate(row["S0"][:3]):
        sources[vertex][i] = 1
    sources[row["S0"][-1]] = [-1, -1, -1]
    plans = candidate_trees(graph, sources)
    relations = CurrentRelations(row)
    print("Whole tensor", index, "requires", len(relations.free), "independent entries", flush=True)
    attempts = []
    for score, tree in plans[:4]:
        engine = None
        started = time.monotonic()
        def expired(*unused):
            raise TimeoutError("The whole tensor plan exceeded its declared time budget.")
        signal.signal(signal.SIGALRM, expired)
        signal.alarm(seconds)
        try:
            engine = Bulk3Periods(row, tree=tree)
            bare = engine.integrate()
            if result["bare_status"] == "complete" and encoded(bare) != result["bare_coefficient"]:
                raise ArithmeticError("The alternative tree changed a complete bare value.")
            if "complete_bare_numerator" in result and polynomial_record(engine._bare) != result["complete_bare_numerator"]:
                raise ArithmeticError("The alternative tree changed the complete original numerator.")
            old_targets = {(entry["edge"], entry["source_direction"]): entry
                           for entry in result.get("complete_finite_targets", [])
                           if entry["edge"] is not None}
            pushed = ChainCurrentPeriods(engine)
            independent, reused = {}, []
            for position in relations.free:
                edge, component = divmod(position, 3)
                if (edge, component) in old_targets:
                    entry = old_targets[edge, component]
                    independent[position] = {key: Fraction(value) for key, value in entry["value"].items()}
                    reused.append(entry)
                else:
                    independent[position] = pushed.integrate(edge, component)
                print("Whole tensor", index, "entry", edge, component, "done;",
                      round(time.monotonic() - started, 3), "seconds", flush=True)
            values = relations.reconstruct(bare, independent)
            paths = [Path(__file__), HERE / "chain_currents.py", HERE / "current_relations.py",
                     HERE / "planned_cases.py"]
            result.update(bare_status="complete", bare_coefficient=encoded(bare),
                          first_current_status="complete", source_cycle_tree=list(tree),
                          method="complete_current_relations_and_exact_series_fibres",
                          complete_bare_numerator=polynomial_record(engine._bare),
                          original_grounded_laplacian_density_checks=engine.density_checks,
                          complete_series_chains=[list(chain.edges) for chain in engine._chains],
                          scalar_core_edges=[list(edge) for edge in engine._core],
                          complete_finite_targets=engine.integrations + reused + pushed.integrations,
                          complete_current_relations=relations.record(),
                          all_original_first_current_coefficients=[
                              [encoded(value) for value in item] for item in values],
                          integrated_original_source_identities=3 * row["nv"],
                          complete_loop_adjugate_polynomial_identities=25,
                          complete_tree_planning={"considered_trees": len(plans), "selected_score": score,
                                                  "failed_complete_plans": attempts},
                          planning_dependencies_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                                                        for p in paths},
                          completion_seconds=round(time.monotonic() - started, 3))
            result.pop("error", None)
            atomic_json(path, result)
            print("Completed every original current vector for", index, flush=True)
            return
        except Exception as error:
            attempts.append({"tree": list(tree), "score": score,
                             "failure": type(error).__name__ + ": " + str(error)})
            print("Whole tensor plan failed", index, attempts[-1], flush=True)
        finally:
            signal.alarm(0)
            if engine is not None:
                engine.close()
    raise ArithmeticError("All bounded whole-tensor plans failed for graph %d." % index)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--indices", default="47,48,51,53-57")
    parser.add_argument("--seconds", type=int, default=600)
    args = parser.parse_args()
    for index in selected(args.indices):
        complete(index, args.seconds)
