# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 13:01:27 2022

@author: jogib
"""

import quimb as qu
import quimb.tensor as qtn
import matplotlib.pyplot as plt
import numpy as np

#%%
b = qu.bell_state("psi-")
k = qu.bell_state("psi+")
mat = qu.pauli("X") & qu.pauli("Y")

#%%
"""
from example
"""
data = qu.bell_state("psi-").reshape(2, 2)
inds = ("k0", "k1")
tags = ("KET",)

ket = qtn.Tensor(data=data, inds=inds, tags=tags)
ket
# ket.draw()
X = qtn.Tensor(qu.pauli("X"), inds=("k0", "b0"), tags=["PAULI", "X", "0"])
Y = qtn.Tensor(qu.pauli("Y"), inds=("k1", "b1"), tags=["PAULI", "Y", "1"])

bra = qtn.Tensor(qu.bell_state("psi+").reshape(2, 2), inds=("b0", "b1"), tags=["BRA"])
#%% evaluate things
TN = ket.H & X & Y & bra
print(TN)
TN.draw(color=["KET", "PAULI", "BRA"], figsize=(4, 4), show_inds="all")
TN.trace("b1", "b0")
#%%
"""
simple bell state prep
"""
init = qtn.Tensor(
    qu.computational_state("01").reshape(2, 2), inds=("k0", "k1"), tags=("KET", "INIT",)
)
had_on_1 = qtn.Tensor(qu.hadamard(), inds=("k0", "b0"), tags=("HADAMARD", "0"))
ident_on_2 = qtn.Tensor(qu.identity(2), inds=("k1", "b1"), tags=("IDENTITY", "1"))
cn = qtn.Tensor(qu.CNOT(), inds=("b0", "b1"), tags=("CNOT", "0", "1"))

TN = cn & had_on_1 & ident_on_2 & init
print(TN)
#%%
N = 2
circ = qtn.Circuit(N)
circ.apply_gate("H", 0)
circ.apply_gate("CNOT", 0, 1)
#%%

# init two quibit
circ = qtn.Circuit(N=2, psi0=MPS_computational_state("01"))

#%% section 2.3
L = 5

# create the nodes, by default just the scalar 1.0
tensors = [qtn.Tensor() for _ in range(L)]

for i in range(L):
    # add the physical indices, each of size 2
    tensors[i].new_ind(f'k{i}', size=2)


tensors[0].new_bond(tensors[1], size=7)
tensors[0].new_bond(tensors[2], size=7)
tensors[1].new_bond(tensors[2], size=7)
tensors[1].new_bond(tensors[3], size=7)
tensors[2].new_bond(tensors[3], size=7)
tensors[3].new_bond(tensors[4], size=7)

mps = qtn.TensorNetwork(tensors)
mps.draw()
#%%
print(qu.bell_state(3))
#%% 4 split with svd to get to a mps
# create a tensor with 5 legs
t = qtn.rand_tensor([2, 3, 4, 5, 6], inds=['a', 'b', 'c', 'd', 'e'])
t.draw(initial_layout='kamada_kawai', figsize=(3, 3))
# split the tensor, by grouping some indices as 'left'
tn = t.split(['a', 'c', 'd'])
tn.draw(figsize=(3, 3))


#%%
a=qu.qu([1],qtype='dop')
b=qu.qu([0],qtype='dop')
print(a@b)
a=qu.qu([[1,0],[0,0]],qtype='dop')
b=qu.qu([[0,0],[0,0]],qtype='dop')
print(a@b)
a=qu.qu([[1,0],[0,1]],qtype='dop')
b=qu.qu([[0,1],[0,0]],qtype='dop')
print(a.H@b)