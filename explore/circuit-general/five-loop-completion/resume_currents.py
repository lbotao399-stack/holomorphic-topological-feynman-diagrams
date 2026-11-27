"""Checkpoint every complete current integral using an already verified source tree."""
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


def complete(index, seconds):
    path = HERE / "cases" / ("bulk3_%03d.json" % index)
    result = json.loads(path.read_text())
    if result["dependencies_sha256"] != hashes():
        raise ArithmeticError("The complete bare input has changed exact dependencies.")
    if result["first_current_status"] == "complete":
        print("Already complete", index, flush=True)
        return
    if result["bare_status"] != "complete" or "source_cycle_tree" not in result:
        raise ValueError("This replay needs a complete bare polynomial and its successful source tree.")
    paths = [Path(__file__), HERE / "chain_currents.py", HERE / "current_relations.py"]
    dependencies = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    checkpoint_path = HERE / "partial_currents" / ("bulk3_%03d.json" % index)
    checkpoint = {"index": index, "complete_original_numerator": result["complete_bare_numerator"],
                  "dependencies_sha256": {**hashes(), **dependencies}, "complete_individual_targets": []}
    if checkpoint_path.exists():
        saved = json.loads(checkpoint_path.read_text())
        if any(saved[key] != checkpoint[key] for key in ("index", "complete_original_numerator", "dependencies_sha256")):
            raise ArithmeticError("An individual-target checkpoint has stale exact dependencies.")
        checkpoint = saved
    old = {(entry["edge"], entry["source_direction"]): entry
           for entry in result.get("complete_finite_targets", []) if entry["edge"] is not None}
    old.update({(entry["edge"], entry["source_direction"]): entry
                for entry in checkpoint["complete_individual_targets"]})
    started, engine = time.monotonic(), None
    def expired(*unused):
        raise TimeoutError("The complete graph budget expired; every finished whole target is checkpointed.")
    signal.signal(signal.SIGALRM, expired)
    signal.alarm(seconds)
    try:
        engine = Bulk3Periods(ROWS[index], tree=result["source_cycle_tree"])
        if polynomial_record(engine._bare) != result["complete_bare_numerator"]:
            raise ArithmeticError("The complete original numerator changed.")
        bare = {key: Fraction(value) for key, value in result["bare_coefficient"].items()}
        relations, pushed = CurrentRelations(ROWS[index]), ChainCurrentPeriods(engine)
        values, used = {}, []
        print("Checkpointed whole tensor", index, "requires", len(relations.free), "entries", flush=True)
        for position in relations.free:
            edge, component = divmod(position, 3)
            if (edge, component) not in old:
                pushed.integrate(edge, component)
                entry = pushed.integrations[-1]
                checkpoint["complete_individual_targets"].append(entry)
                atomic_json(checkpoint_path, checkpoint)
                old[edge, component] = entry
            entry = old[edge, component]
            values[position] = {key: Fraction(value) for key, value in entry["value"].items()}
            used.append(entry)
            print("Checkpointed", index, "entry", edge, component, "done;",
                  round(time.monotonic() - started, 3), "seconds", flush=True)
        currents = relations.reconstruct(bare, values)
        bare_targets = [entry for entry in result.get("complete_finite_targets", []) if entry["edge"] is None]
        if len(bare_targets) != 1:
            raise ArithmeticError("The existing complete bare target is missing or duplicated.")
        result.update(first_current_status="complete",
                      method="complete_current_relations_and_checkpointed_exact_series_fibres",
                      complete_finite_targets=bare_targets + used,
                      complete_current_relations=relations.record(),
                      all_original_first_current_coefficients=[
                          [encoded(value) for value in spinor] for spinor in currents],
                      integrated_original_source_identities=3 * ROWS[index]["nv"],
                      complete_loop_adjugate_polynomial_identities=25 if pushed._ready else 0,
                      planning_dependencies_sha256={**result.get("planning_dependencies_sha256", {}), **dependencies},
                      completion_seconds=round(time.monotonic() - started, 3))
        result.pop("error", None)
        atomic_json(path, result)
        print("Completed every original current vector for", index, flush=True)
        checkpoint_path.unlink(missing_ok=True)
    finally:
        signal.alarm(0)
        if engine is not None:
            engine.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--indices", default="54-57")
    parser.add_argument("--seconds", type=int, default=900)
    args = parser.parse_args()
    for index in selected(args.indices):
        complete(index, args.seconds)
