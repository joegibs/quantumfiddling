# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 12:45:57 2022

@author: jogib
"""

from quimb import *

data = [1, 2j, -3]
#%%
print(qu(data, qtype='dop'))

#%%
import itertools
from operator import add
import numpy as np
from quimb import *
#%%
def ham_heis_2D(n, m, j=1.0, bz=0.0, cyclic=False,
                sparse=True):

    dims = [[2] * m] * n  # shape (n, m)

    # generate tuple of all site coordinates
    sites = tuple(itertools.product(range(n), range(m)))

    # generate neighbouring pairs of coordinates
    def gen_pairs():
        for i, j, in sites:
            above, right = (i + 1) % n, (j + 1) % m
            # ignore wraparound coordinates if not cyclic
            if cyclic or above != 0:
                yield ((i, j), (above, j))
            if cyclic or right != 0:
                yield ((i, j), (i, right))

    # generate all pairs of coordinates and directions
    pairs_ss = tuple(itertools.product(gen_pairs(), 'xyz'))

    # make XX, YY and ZZ interaction from pair_s
    #     e.g. arg ([(3, 4), (3, 5)], 'z')
    def interactions(pair_s):
        pair, s = pair_s
        Sxyz = spin_operator(s, sparse=True)
        return ikron([j * Sxyz, Sxyz], dims, inds=pair)

    # function to make Z field at ``site``
    def fields(site):
        Sz = spin_operator('z', sparse=True)
        return ikron(bz * Sz, dims, inds=[site])

    # combine all terms
    all_terms = itertools.chain(map(interactions, pairs_ss),
                                map(fields, sites) if bz != 0.0 else ())
    H = sum(all_terms)

    # can improve speed of e.g. eigensolving if known to be real
    if isreal(H):
        H = H.real

    if not sparse:
        H = qarray(H.A)

    return H
#%%
n, m = 4, 5
dims = [[2] * m] * n

for row in dims:
    print(row)
#%%
H = ham_heis_2D(4, 5, cyclic=False)
#%%
H = H + 0.2 * ikron(spin_operator('Z', sparse=True), dims, [(1, 2)])

#%%
Sz = spin_operator('Z', stype='coo')
Sz_ij = [[ikron(Sz, dims, [(i, j)])
          for j in range(m)]
         for i in range(n)]

#%%
m_ij = [[expec(Sz_ij[i][j], gs)
         for j in range(m)]
        for i in range(n)]

%matplotlib inline
import matplotlib.pyplot as plt

plt.imshow(m_ij)
plt.colorbar()

