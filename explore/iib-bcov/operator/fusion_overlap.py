"""Gluing all thirteen anchor fusions of a generic fourteen-wave configuration."""
from itertools import combinations
from math import comb

import numpy as np
from flint import nmod_mat

from fusion_fiber import Fiber, PRIME, fixture, kernel_basis, modular_matrix, modular_rank


def check(seed=70992):
    original = np.asarray(fixture(14,seed),dtype=np.int64)
    anchor = 13
    faces = list(combinations(range(13),2))
    common_choices = {}
    triple_rows = {}
    difference = np.zeros((len(faces)*33,13*513),dtype=np.int16)
    for i in range(13):
        singles = [s for s in range(13) if s!=i]
        groups = [1<<s for s in singles]+[(1<<i) | (1<<anchor)]
        momenta = [tuple(original[s]) for s in singles]+[tuple(original[i]+original[anchor])]
        fiber = Fiber(momenta)
        kernel,dim = modular_matrix(fiber.trace_map()).nullspace()
        assert dim==513
        keys = {number:tuple(sorted(sum(groups[s] for s in part) for part in partition))
                for partition,number in fiber.partition_ids.items()}
        for slot,j in enumerate(singles):
            edge = tuple(sorted((i,j)))
            block = fiber.pairs.index((slot,12))
            columns,values = fiber.blocks[block]
            order = sorted(range(len(columns)),key=lambda k:keys[columns[k]])
            columns = columns[order]
            canonical = [keys[k] for k in columns]
            if edge not in common_choices:
                face_map = np.asarray(values[order].T,dtype=np.int64)
                echelon,rank = nmod_mat(face_map.tolist(),PRIME).rref()
                assert rank==33
                pivots = np.argmax(np.array(echelon.tolist(),dtype=np.int64)[:rank]!=0,axis=1)
                free = np.setdiff1d(np.arange(66),pivots)
                # Each of the eleven quadruple labels can be prescribed
                # independently in the face kernel; this fixes the rank of
                # the triple-overlap conditions on all discrepancies.
                face_kernel = kernel_basis(face_map)
                assert np.array_equal(face_kernel[free],np.eye(33,dtype=np.int64))
                four = [k for k,parts in enumerate(canonical) if max(p.bit_count() for p in parts)==4]
                assert len(four)==11 and modular_rank(face_kernel[four])==11
                known = (1<<anchor)|(1<<i)|(1<<j)
                for index in four:
                    quadruple = next(p for p in canonical[index] if p.bit_count()==4)
                    k = (quadruple ^ known).bit_length()-1
                    triple_rows[edge,k] = face_kernel[index]
                common_choices[edge] = canonical,free
            else:
                previous,free = common_choices[edge]
                assert previous==canonical
            start = 33*faces.index(edge)
            sign = 1 if i<j else -1
            for k,col in enumerate(columns[free]):
                difference[start+k,i*513:(i+1)*513] = [sign*int(kernel[int(col),h]) for h in range(513)]
        del kernel
        print('anchor',i+1,'of 13',flush=True)
    for size in (10,11,12,13):
        selected = [pos for pos,(_,last) in enumerate(faces) if last<size]
        part = difference.reshape(len(faces),33,-1)[selected,:,:size*513].reshape(-1,size*513)
        rank = modular_rank(part)
        compatible = 33*len(selected)-comb(size,3)
        assert rank >= {10:1365,11:1650,12:1957,13:2284}[size]
        assert rank <= compatible
        print(size,'parents, difference rank',rank,'compatible discrepancies',compatible,
              'residual',compatible-rank,flush=True)
    # Independent equations on the discrepancies come first from the
    # quadruple block common to each triple of parents. Extend their row
    # space to the complete left kernel of the compatibility matrix.
    triples = list(combinations(range(13),3))
    compatibility = np.zeros((len(triples),len(faces),33),dtype=np.int64)
    for r,(i,j,k) in enumerate(triples):
        for edge,m,sign in (((i,j),k,1),((i,k),j,-1),((j,k),i,1)):
            compatibility[r,faces.index(edge)] = sign*triple_rows[edge,m]
    compatibility = compatibility.reshape(len(triples),-1)
    assert modular_rank(compatibility)==len(triples)
    dual,dim = modular_matrix(difference).transpose().nullspace()
    certificates = np.array([[int(dual[i,j]) for j in range(dim)]
                             for i in range(difference.shape[0])],dtype=np.int64)
    form,_ = modular_matrix(compatibility).rref()
    form = np.array(form.tolist(),dtype=np.int64)
    pivots = np.argmax(form!=0,axis=1)
    certificates = (certificates-form.T@certificates[pivots])%PRIME
    basis,independent = modular_matrix(certificates.T).rref()
    basis = np.array(basis.tolist()[:independent],dtype=np.int64)
    assert independent == compatible-rank
    assert not np.any((basis@difference.astype(np.int64))%PRIME)
    print('independent consistency tests modulo triple overlaps',independent,
          'pivot coordinates',np.argmax(basis!=0,axis=1).tolist(),flush=True)
    return basis


if __name__=='__main__':
    check()
