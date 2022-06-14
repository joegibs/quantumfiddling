# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 19:16:27 2022

@author: jogib
"""

import itertools
from operator import add
import numpy as np
from quimb import *
import matplotlib.pyplot as plt
#%%
import random
import quimb as qu
import quimb.tensor as qtn
#%%

# 10 qubits and tag the initial wavefunction tensors
circ = qtn.Circuit(N=10)

# initial layer of hadamards
for i in range(10):
    circ.apply_gate('H', i, gate_round=0)

# 8 rounds of entangling gates
for r in range(1, 9):

    # even pairs
    for i in range(0, 10, 2):
        circ.apply_gate('CNOT', i, i + 1, gate_round=r)

    # odd pairs
    for i in range(1, 9, 2):
        circ.apply_gate('CZ', i, i + 1, gate_round=r)


# final layer of hadamards
for i in range(10):
    circ.apply_gate('H', i, gate_round=r + 1)

circ

#%%
num_elems = 10

circ = qtn.Circuit(N=num_elems)

for i in range(num_elems):
    circ.apply_gate('H', i, gate_round=0)

for r in range(1, 10):

    # even pairs
    for i in range(0, num_elems, 2):
        circ.apply_gate('CNOT', i, i + 1, gate_round=r)
    for i in range(num_elems):
        circ.apply_gate('RZ', 0, i, gate_round=r)
    # odd pairs
    for i in range(1, num_elems-1, 2):
        circ.apply_gate('CZ', i, i + 1, gate_round=r)
    for i in range(num_elems):
        circ.apply_gate('RZ', 0, i, gate_round=r)
        
for i in range(num_elems):
    circ.apply_gate('H', i, gate_round=r + 1)