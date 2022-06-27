# -*- coding: utf-8 -*-
"""
Created on Sun Jun 26 16:15:11 2022

@author: jogib
"""

from quimb import *
import matplotlib.pyplot as plt

#%%
class rabi_ham:
    def __init__(self, L):
        # self.h_interaction = qarray(np.array([[1,0],[0,2]]))
        self.h_field = ham_ising(L, sparse=True, jz=0.0, bx=1.0, cyclic=False)

    def __call__(self, t):
        return np.array([[3, 0], [0, 3]]) + np.array([[0, 3], [3, 0]])


#%%
L = 1
# our initial state
psi0 = rand_product_state(L)
# instantiate the ham object, it's __call__ method will be used by Evolution
fn_ham_t = rabi_ham(L)
#%%
compute = {
    "time": lambda t, p: t,
    "state1": lambda t, p, ham: p[0].H @ p[0],
    "state2": lambda t, p, ham: p[1].H @ p[1],
}

evo = Evolution(psi0, fn_ham_t, progbar=True, compute=compute, int_small_step=True)

evo.update_to(5)
#%%
plt.plot(
    evo.results["time"],
    evo.results["state1"],
    evo.results["time"],
    evo.results["state2"],
)

#%%
#%%
b0 = 1
b1 = 1
w = 0.1


class rabi_ham:
    def __init__(self, L):
        # self.h_interaction = qarray(np.array([[1,0],[0,2]]))
        self.h_field = ham_ising(L, sparse=True, jz=0.0, bx=1.0, cyclic=False)

    def __call__(self, t):
        return (
            b0 * pauli("Z")
            + b1 * cos(w * t) * pauli("X")
            + b1 * sin(w * t) * pauli("Y")
        )


#%%
L = 1
# our initial state
psi0 = rand_product_state(L)
# instantiate the ham object, it's __call__ method will be used by Evolution
fn_ham_t = rabi_ham(L)
#%%
compute = {
    "time": lambda t, p: t,
    "state1": lambda t, p, ham: p[0].H @ p[0],
    "state2": lambda t, p, ham: p[1].H @ p[1],
}

evo = Evolution(psi0, fn_ham_t, progbar=True, compute=compute, int_small_step=True)

evo.update_to(20)
#%%
plt.plot(
    evo.results["time"],
    evo.results["state1"],
    evo.results["time"],
    evo.results["state2"],
)
