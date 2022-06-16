# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 13:01:27 2022

@author: jogib
"""

import quimb as qu
import quimb.tensor as qtn
import matplotlib.pyplot as plt

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
