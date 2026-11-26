"""Unimodular cycle and lattice bases for exact mixed four-loop current minors."""
from copy import copy
from functools import lru_cache
from itertools import combinations
from math import lcm
import os

import networkx as nx
import numpy as np
from flint import fmpq_mat,fmpz_mat
from mixed_prism_empty_stars import normalized_roots,soft_polynomials
from quartic_seven_plane_double import numerator,polynomial_det


def reduce_rows(array):
    n,m,variables=array.shape
    matrix=fmpz_mat(array.reshape(n,m*variables).tolist())
    reduced,transform=matrix.lll(transform=True)
    orientation=int(transform.det())
    assert abs(orientation)==1 and transform*matrix==reduced
    return np.array(reduced.tolist(),dtype=object).reshape(n,m,variables),orientation


def sparse_matrix(matrix,ctx):
    n=len(matrix)
    variables=ctx.gens()
    array=np.zeros((n,n,len(variables)),dtype=object)
    for i,row in enumerate(matrix):
        for j,entry in enumerate(row):
            if entry:
                for ex,c in entry.terms():
                    powers=list(map(int,ex))
                    assert sum(powers)==1
                    array[i,j,powers.index(1)]=int(c)
    orientation=1
    previous=None
    for _ in range(5):
        array,sign=reduce_rows(array)
        orientation*=sign
        array,sign=reduce_rows(array.transpose(1,0,2))
        array=array.transpose(1,0,2)
        orientation*=sign
        score=np.count_nonzero(np.any(array,axis=2)),np.count_nonzero(array)
        if score==previous:break
        previous=score
    rows=[[sum(int(c)*variables[k] for k,c in enumerate(e) if c) for e in row] for row in array]
    return rows,orientation


def laplace_determinant(matrix,ctx):
    n=len(matrix)
    full=(1<<n)-1
    row=[sum(1<<j for j,p in enumerate(entries) if p) for entries in matrix]
    col=[sum(1<<i for i in range(n) if row[i]>>j&1) for j in range(n)]
    def bits(mask):
        while mask:
            v=mask&-mask
            yield v.bit_length()-1
            mask-=v
    @lru_cache(maxsize=1024)
    def minor(left,right):
        if not left:return ctx.constant(1)
        # Expanding on the least populated surviving row or column keeps every
        # product linear in one factor, including the final determinant step.
        count,indices,slot,transpose=min(
            [(int.bit_count(row[i]&right),row[i]&right,i,False) for i in bits(left)]+
            [(int.bit_count(col[j]&left),col[j]&left,j,True) for j in bits(right)])
        total=ctx.constant(0)
        for item in bits(indices):
            i,j=(item,slot) if transpose else (slot,item)
            value=minor(left^(1<<i),right^(1<<j))
            if value:
                parity=((left&((1<<i)-1)).bit_count()+(right&((1<<j)-1)).bit_count())%2
                term=matrix[i][j]*value
                total=total-term if parity else total+term
        return total
    try:
        return minor(full,full)
    finally:
        minor.cache_clear()


def sparse_determinant(matrix,ctx):
    rows,sign=sparse_matrix(matrix,ctx)
    return sign*(laplace_determinant if os.environ.get('MIXED_LAPLACE_DET')=='1' else polynomial_det)(rows,ctx)


def block_work(matrix):
    n=len(matrix)
    graph=nx.Graph()
    graph.add_nodes_from(range(2*n))
    graph.add_edges_from((i,n+j) for i,row in enumerate(matrix) for j,p in enumerate(row) if p)
    matching=nx.bipartite.maximum_matching(graph,range(n))
    if len(matching)!=2*n:return 0
    directed=nx.DiGraph()
    directed.add_nodes_from(range(n))
    directed.add_edges_from((i,matching[j]) for i in range(n) for j in graph[i])
    return sum(2**len(ids)*sum(len(matrix[i][matching[j]-n]) for i in ids for j in ids
                    if matrix[i][matching[j]-n]) for ids in nx.strongly_connected_components(directed))


def adapted_minors(probe,data):
    original=probe.cycles
    possibilities=[]
    for choice in combinations(range(9),4):
        block=fmpq_mat([[row[j] for j in choice] for row in original])
        if abs(block.det())==1:possibilities.append((len(set(choice).intersection(probe.marks)),choice,block))
    coverage=max(n for n,_,_ in possibilities)
    pairs=tuple(combinations(range(4),2))
    best={}
    for n,choice,block in possibilities:
        if n!=coverage:continue
        cycles=[list(map(int,row)) for row in (block.inv()*fmpq_mat(original)).tolist()]
        graph=copy(probe.graph)
        graph.C=[[r[j] for j,path in enumerate(data[5]) for _ in path] for r in cycles]
        for i,j in pairs:
            # In the new cycle coordinates the two unit jumps are columns i,j
            # of the original chord block. Its second minors give their coordinates.
            direction=tuple(int(block[a,i]*block[b,j]-block[a,j]*block[b,i]) for a,b in pairs)
            cost=[]
            def inspect(matrix,ctx):
                cost.append(block_work(sparse_matrix(matrix,ctx)[0]))
                return ctx.constant(0)
            numerator(graph,probe.sheared,probe.flow,[(3,i),(3,j)],determinant=inspect)
            work=sum(cost)
            if direction not in best or work<best[direction][0]:
                best[direction]=work,graph,i,j
    independent=[]
    jobs=[]
    for direction,(cost,graph,i,j) in sorted(best.items(),key=lambda entry:entry[1][0]):
        if fmpz_mat([*independent,direction]).rank()==len(independent):continue
        independent.append(direction)
        jobs.append((cost,graph,i,j))
        if len(independent)==6:break
    assert len(independent)==6
    values=[numerator(graph,probe.sheared,probe.flow,[(3,i),(3,j)],determinant=sparse_determinant)
            if cost else data[1].constant(0) for cost,graph,i,j in jobs]
    mins={}
    for pair,row in zip(pairs,fmpq_mat(independent).inv().tolist(),strict=True):
        denominator=lcm(*(int(c.denominator) for c in row))
        integer=sum((int(c*denominator)*v for c,v in zip(row,values) if c),data[1].constant(0))
        mins[pair]=integer//denominator
        assert denominator*mins[pair]==integer
    return mins


def adapted_roots(probe,data):
    mins=adapted_minors(probe,data)
    soft=soft_polynomials(probe.cycles,data[3],data[1])
    for root,rays,scale in normalized_roots(mins,soft):
        yield root,scale*soft[rays]


def check():
    from mixed_bare_primitive_pole import PrimitiveChamber
    from mixed_prism_empty_stars import prism_data
    from mixed_prism_rung_sum import sparse_flow
    fixtures=(((0,0,0,5,5,5),(0,0,2,1,1,2)),((0,0,5,5,8,8),(0,1,0,2,1,2)),
              ((0,0,1,1,3,3),(1,2,1,0,0,2)))
    for marks,color in fixtures:
        probe=PrimitiveChamber(4,marks,color)
        _,probe.sheared,probe.flow=sparse_flow(probe.graph,color,probe.sheared)
        data,_=prism_data(marks,4)
        mins=adapted_minors(probe,data)
        for shift in (1,4):
            times=[1+(13*k*k+shift*k+7)%37 for k in range(15)]
            weighted=[data[1].constant(t) for t in times]
            for (i,j),poly in mins.items():
                old=numerator(probe.graph,probe.sheared,probe.flow,[(3,i),(3,j)],edge_times=weighted)
                assert old==data[1].constant(poly(*times))
        print('Six signed bivector coordinates at two metrics:',marks,color,flush=True)


if __name__=='__main__':
    check()
