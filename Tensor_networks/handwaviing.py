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
circ = qtn.Circuit(N=2, psi0=qtn.MPS_computational_state("01"))

#%% section 2.3
L = 5

# create the nodes, by default just the scalar 1.0
tensors = [qtn.Tensor() for _ in range(L)]

for i in range(L):
    # add the physical indices, each of size 2
    tensors[i].new_ind(f"k{i}", size=2)


tensors[0].new_bond(tensors[1], size=7)
tensors[0].new_bond(tensors[2], size=7)
tensors[1].new_bond(tensors[2], size=7)
tensors[1].new_bond(tensors[3], size=7)
tensors[2].new_bond(tensors[3], size=7)
tensors[3].new_bond(tensors[4], size=7)

mps = qtn.TensorNetwork(tensors)
mps.draw()

#%% contractions
L = 6

# create the nodes, by default just the scalar 1.0
tensors = [qtn.Tensor(tags=[str(i)]) for i in range(L)]

for i in range(L):
    # add the physical indices, each of size 2
    tensors[i].new_ind(f"k{i}", size=2)

for i in range(0, L - 2, 2):
    tensors[i].new_bond(tensors[i + 1], size=7)
    tensors[i].new_bond(tensors[i + 2], size=7)
tensors[L - 2].new_bond(tensors[L - 1], size=7)
for i in range(1, L - 2, 2):
    print(i, i + 2)
    tensors[i].new_bond(tensors[i + 2], size=7)

mps = qtn.TensorNetwork(tensors)
mps.draw()
mp = mps.contract(tags=["0", "2"])
mp.draw()
#%%
"""
Problems
"""
# Section 1
#%% 1
a = np.array([[b ** 2 - 2 * a for b in range(3)] for a in range(3)])
b = np.array(
    [[[(-((3) ** a)) * g + d for a in range(3)] for g in range(3)] for d in range(3)]
)
c = np.array([[e for i in range(3)] for e in range(3)])
d = np.array([[[b * g * e for b in range(3)] for g in range(3)] for e in range(3)])

tensors = [
    qtn.Tensor(data=a, inds=["b", "a"], tags=["A"]),
    qtn.Tensor(data=b, inds=["a", "g", "d"], tags=["B"]),
    qtn.Tensor(data=c, inds=["d", "e"], tags=["C"]),
    qtn.Tensor(data=d, inds=["e", "b", "g"], tags=["D"]),
]
mps = qtn.TensorNetwork(tensors)
mps.draw()
print(mps ^ ...)
"""
need to check indicies simpler examples work as expected idk
"""
#%%
s = 0
for a in range(3):
    for b in range(3):
        for g in range(3):
            for d in range(3):
                for e in range(3):
                    s = s + (b ** 2 - 2 * a) * ((-((3) ** a)) * g + d) * g * b * e ** 2
print(s)
#%% simpler example to indicie check

a = np.array([[a + b for a in range(3)] for b in range(3)])
b = np.array([[a + 3 ** b for a in range(3)] for b in range(3)])
tensors = [
    qtn.Tensor(data=a, inds=["a", "b"], tags=["A"]),
    qtn.Tensor(data=b, inds=["a", "b"], tags=["B"]),
]
mps = qtn.TensorNetwork(tensors)
mps.draw()
print(mps ^ ...)

s = 0
for a in range(3):
    for b in range(3):
        s = s + (b + a) * (3 ** b + a)
print(s)
#%%
a = np.array([[a + b for a in range(3)] for b in range(3)])
b = np.array([[c + 3 ** b for c in range(3)] for b in range(3)])
b = np.array([[a + 3 ** c for a in range(3)] for c in range(3)])
tensors = [
    qtn.Tensor(data=a, inds=["a", "b"], tags=["A"]),
    qtn.Tensor(data=b, inds=["c", "b"], tags=["B"]),
    qtn.Tensor(data=b, inds=["a", "c"], tags=["C"]),
]
mps = qtn.TensorNetwork(tensors)
mps.draw()
print(mps ^ ...)

s = 0
for a in range(3):
    for b in range(3):
        for c in range(3):
            s = s + (b + a) * (3 ** b + c) * (a + 3 ** c)
print(s)
#%% 2 table

a = np.array([[b ** 2 - 2 * a for b in range(3)] for a in range(3)])
b = np.array(
    [[[(-((3) ** a)) * g + d for a in range(3)] for g in range(3)] for d in range(3)]
)

L = 5
tensors = [qtn.Tensor(data=a, inds=["b", "a"], tags=["A"])]


mps = qtn.TensorNetwork(tensors)
print(mps.all_inds())
mps.add_tensor(qtn.Tensor(data=a, inds=["g", "a"], tags=["B"]))
mp = mps ^ ...
print(mp.inds)
mps.add_tensor(qtn.Tensor(data=b, inds=["g", "d", "e"], tags=["C"]))
mp = mps ^ ...
print(mp.inds)
mps.add_tensor(qtn.Tensor(data=a, inds=["b", "d"], tags=["D"]))
mp = mps ^ ...
print(mp.inds)
mps.add_tensor(qtn.Tensor(data=[1], inds=["e"], tags=["e"]))
mp = mps ^ ...
print(mp)

"""
data doesn't really matter just filling so that i have appropriatly sized arrays'
"""

#%%
"""
qi examples
"""
#%%
print(qu.bell_state(3))
#%% 4 split with svd to get to a mps
# create a tensor with 5 legs
inds = ["a", "b", "c", "d", "e"]
t = qtn.rand_tensor([2, 3, 4, 5, 6], inds=inds)
t.draw(initial_layout="kamada_kawai", figsize=(3, 3))
# split the tensor, by grouping some indices as 'left'

tn = t.split(["a"])
tn.draw(figsize=(3, 3))
tb = tn.split(["b"])
tb.draw(figsize=(3, 3))

#%% 4.2
a = qu.qu([1], qtype="dop")
b = qu.qu([0], qtype="dop")
print(a @ b)
a = qu.qu([[1, 0], [0, 0]], qtype="dop")
b = qu.qu([[0, 0], [0, 0]], qtype="dop")
print(a @ b)

print(a.H @ b)
#%%
a = np.array([[1, 0], [0, 0]])
b = np.array([[0, 0], [0, 1]])
mps= qtn.tensor_1d.MatrixProductState([a,b])

#%%
L = 2

# create the nodes, by default just the scalar 1.0
tensors = [qtn.Tensor() for _ in range(L)]

for i in range(L):
    # add the physical indices, each of size 2
    tensors[i].new_ind(f"k{i}", size=2)

    # add bonds between neighbouring tensors, of size 7
    if i != L - 1:
        tensors[i].new_bond(tensors[(i + 1) % L], size=1)
a = qu.qu([1], qtype="dop")
b = qu.qu([0], qtype="dop")
tensors[0].modify(data=a)
tensors[1].modify(data=b)
mps = qtn.TensorNetwork(tensors)
mps.draw()


#%%
p = qtn.MPS_rand_state(L=20, bond_dim=1)
print(f"Site tags: '{p.site_tag_id}', site inds: '{p.site_ind_id}'")
