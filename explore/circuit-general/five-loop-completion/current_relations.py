"""Complete original-current reconstruction from graph symmetries and sources.

These are identities on the entire 16 by 3 coefficient array, including the
affine bare term in every original vertex conservation equation.  Gaussian
elimination only chooses which integrals to evaluate; all output entries
are reconstructed and checked against every original equation.
"""
from fractions import Fraction

from flint import fmpq_mat
import networkx as nx

from bulk3_periods import add


def parity(values):
    return (-1) ** sum(x > y for i, x in enumerate(values) for y in values[i + 1:])


class CurrentRelations:
    def __init__(self, row):
        sources = row["S0"]
        frame = [[int(i == j) for j in range(3)] for i in range(3)] + [[-1, -1, -1]]
        edges = list(map(tuple, row["edges"]))
        lookup = {tuple(sorted(edge)): i for i, edge in enumerate(edges)}
        relations = set()
        self.automorphisms = []
        for vertex in range(row["nv"]):
            for a in range(3):
                line = [0] * 49
                for edge, (left, right) in enumerate(edges):
                    line[3 * edge + a] = int(left == vertex) - int(right == vertex)
                line[48] = -(frame[sources.index(vertex)][a] if vertex in sources else 0)
                relations.add(tuple(line))
        graph = nx.Graph(edges)
        for mapping in nx.algorithms.isomorphism.GraphMatcher(graph, graph).isomorphisms_iter():
            inverse = {mapping[v]: v for v in mapping}
            matrix = [frame[sources.index(inverse[v])] for v in sources[:3]]
            reversals = sum(mapping[left] > mapping[right] for left, right in edges)
            factor = (parity([mapping[v] for v in range(row["nv"])])
                      * (-1) ** reversals
                      * parity([sources.index(mapping[v]) for v in sources]) ** 5)
            self.automorphisms.append({"vertex_permutation": [mapping[v] for v in range(row["nv"])],
                                       "source_matrix": matrix, "bare_multiplier": factor})
            for edge, (left, right) in enumerate(edges):
                image = lookup[tuple(sorted((mapping[left], mapping[right])))]
                edge_sign = 1 if (mapping[left], mapping[right]) == edges[image] else -1
                for a in range(3):
                    line = [0] * 49
                    for b in range(3):
                        line[3 * image + b] += matrix[b][a]
                    line[3 * edge + a] -= factor * edge_sign
                    if any(line):
                        relations.add(tuple(line))
        self.rows = sorted(relations)
        reduced, self.rank = fmpq_mat(self.rows).rref()
        self.reduced = [[Fraction(int(reduced[i, j].p), int(reduced[i, j].q)) for j in range(49)]
                        for i in range(self.rank)]
        self.pivots = [next(j for j, value in enumerate(line) if value) for line in self.reduced]
        self.free = [j for j in range(48) if j not in self.pivots]
        self.bare_forced_zero = 48 in self.pivots

    def check(self, bare, values):
        entries = [value for row in values for value in row] + [bare]
        if len(entries) != 49:
            raise ValueError("Every original edge and all three source directions are required.")
        for row in self.rows:
            residual = {}
            for value, coefficient in zip(entries, row):
                if coefficient:
                    residual = add(residual, value, coefficient)
            if residual:
                raise ArithmeticError("A complete graph-symmetry or affine source relation failed.")

    def reconstruct(self, bare, independent):
        if set(independent) != set(self.free):
            raise ValueError("Every independent integral, and no truncated subset, is required.")
        if self.bare_forced_zero and bare:
            raise ArithmeticError("The complete source relations force a zero bare kernel.")
        values = [None] * 48 + [bare]
        for index, value in independent.items():
            values[index] = value
        for pivot, row in zip(self.pivots, self.reduced):
            if pivot == 48:
                continue
            result = {}
            for index in self.free + [48]:
                if row[index]:
                    result = add(result, values[index], -row[index])
            values[pivot] = result
        if any(value is None for value in values):
            raise ArithmeticError("The full reconstruction left an undetermined original coefficient.")
        result = [values[3 * e:3 * e + 3] for e in range(16)]
        self.check(bare, result)
        return result

    def record(self):
        return {"unknown_order": "original edge major; source directions0,1,2; final column is bare coefficient",
                "automorphisms": self.automorphisms,
                "complete_original_relation_rows": [
                    {str(j): str(value) for j, value in enumerate(row) if value} for row in self.rows],
                "exact_reduced_rows": [
                    {str(j): str(value) for j, value in enumerate(row) if value} for row in self.reduced],
                "rank": self.rank, "independent_current_entries": [divmod(j, 3) for j in self.free],
                "bare_forced_zero": self.bare_forced_zero}
