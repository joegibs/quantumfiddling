# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 07:51:00 2022

@author: jogib
"""

import quimb as qu
import quimb.tensor as qtn
import matplotlib.pyplot as plt
import numpy as np
import itertools

#%%

# init all the tensors
n = 4
lis = ["X", "Y", "Z", "I"]
combos = itertools.combinations_with_replacement(lis, n)
tensors = [qtn.Tensor(tags=["".join(x for x in i)]) for i in combos]

# make connections
def anticommute_check(A, B):
    ch_arr = []
    for i in range(len(A)):
        a_t = qu.pauli(A[i])
        b_t = qu.pauli(B[i])
        check = a_t @ b_t + b_t @ a_t
        ch_arr.append(not np.any(check))
    return np.any(ch_arr)


for i in itertools.combinations([i for i in range(len(tensors))], 2):
    A = list(tensors[i[0]].tags)[0]
    B = list(tensors[i[1]].tags)[0]

    if anticommute_check(A, B):
        # print(A,B)
        tensors[i[0]].new_bond(tensors[i[1]], size=1)

mps = qtn.TensorNetwork(tensors)
mps.delete("".join("I" for i in range(n)))
#%% #plot
mps.draw(
    iterations=0,
    initial_layout="shell",
    edge_color="black",
    edge_alpha=1.0,
    edge_scale=0.5,
)
