"""Independent exact tensor with one derivative on every contracted input.

The original inputs are beta1, beta2, beta3, b, in that index order.
All sixteen entries of the rank-four tensor are retained. This direct
period route supplies an independent oracle for the rational boundary
evaluator; it does not use any boundary identity to evaluate a graph.
"""
from fractions import Fraction
from itertools import product
from pathlib import Path
import argparse
import json
import sys
from time import perf_counter

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'universal-product-jets'))
import routing

jets = routing.jets
OUTPUT = HERE / 'mixed_input_control.json'
VECTORS = ((2, 1), (-1, 4), (3, -2))
GRAPHS = (((0, 1), (0, 3), (0, 4), (1, 2), (1, 4), (2, 3), (3, 4)),
          ((0, 1), (0, 3), (0, 4), (1, 2), (2, 3), (2, 4), (3, 4)),
          ((0, 1), (0, 3), (0, 4), (1, 2), (1, 3), (2, 4), (3, 4)))


def evaluate():
    ledger = routing.certify()['all_original_matching_maps']
    tensors = []
    for edges in GRAPHS:
        record = next(row for row in ledger if row['commutator_slot'] == 1 and row['ordered_edges'] == edges)
        locations = tuple(record['full_source_vertices'][source] for source in
                          ('lambda_gamma4', 'lambda_c', 'lambda_gamma5'))
        rows = jets.base._sources(locations, VECTORS)
        original_edges = tuple(row['edge_position'] for row in record['all_original_input_current_maps'])
        entries = []
        engine = jets.universal.wick.CoreCurrentPeriods(jets.universal.wick.Graph(5, edges), rows,
                                                       marked=(3, 4), check=2)
        try:
            for indices in product((0, 1), repeat=4):
                powers = [0] * 14
                for edge, index in zip(original_edges, indices):
                    powers[2 * edge + index] += 1
                value = jets.base._periods(engine.integrate({tuple(powers): Fraction(1)}, allow_series=False))
                entries.append({'ordered_free_spinor_entries': indices, 'complete_period_value': value})
        finally:
            engine.close()
        tensors.append({'ordered_edges': edges, 'complete_balanced_sources': rows,
                        'original_edges_beta1_beta2_beta3_b': original_edges, 'complete_rank_four_tensor': entries})
        print('Completed direct mixed-input graph', len(tensors), '/', len(GRAPHS), flush=True)
    combinations = []
    for partner in (1, 2):
        entries = []
        for left, right in zip(tensors[0]['complete_rank_four_tensor'], tensors[partner]['complete_rank_four_tensor']):
            value = dict(left['complete_period_value'])
            for label, coefficient in right['complete_period_value'].items():
                value[label] = value.get(label, Fraction(0)) + coefficient
            value = {label: coefficient for label, coefficient in value.items() if coefficient}
            if set(value) - {'1'}:
                raise ArithmeticError('a period survived in a complete mixed-input graph-pair entry')
            entries.append({'ordered_free_spinor_entries': left['ordered_free_spinor_entries'],
                            'complete_rational_value': value.get('1', Fraction(0))})
        combinations.append({'graph_positions': (0, partner), 'complete_rank_four_tensor': entries})
    return jets.base._jsonable({'schema': 1, 'source_spinors_gamma4_c_gamma5': VECTORS,
            'complete_tensor_index_order': ['beta1', 'beta2', 'beta3', 'b'],
            'scope': 'all 16 entries of the entire rank-four tensor, every original source and routed edge retained',
            'common_chord_input_currents': ['y_AB-y_WA-lambda_gamma5', 'y_WA', 'y_WB',
                'lambda_gamma4+lambda_c+lambda_gamma5-y_WB-y_AB'],
            'physical_head_factors_not_applied': '(-2)^4; this is the complete raw geometric current tensor',
            'complete_original_graph_tensors': tensors, 'complete_rational_pair_tensors': combinations})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    start = perf_counter()
    data = evaluate()
    if args.write:
        OUTPUT.write_text(json.dumps(data, indent=2) + '\n', encoding='utf-8')
        print('Wrote the complete independent mixed-input tensor oracle.')
    elif json.loads(OUTPUT.read_text()) != data:
        raise ArithmeticError('a complete mixed-input tensor entry differs')
    else:
        print('Every original graph and rational pair in the mixed-input tensor agrees.')
    print('Elapsed %.3f seconds.' % (perf_counter() - start), flush=True)


if __name__ == '__main__':
    main()
