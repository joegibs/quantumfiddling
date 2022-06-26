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
        return np.array([[3,0],[0,3]])+np.array([[0,3],[3,0]])

#%%
L = 1

# our initial state
psi0 = np.array([1,0])#rand_product_state(L)

# instantiate the ham object, it's __call__ method will be used by Evolution
fn_ham_t = rabi_ham(L)
#%%
compute = {
    'time': lambda t, p: t,
    'state1': lambda t, p, ham: p[0].H@p[0],
    'state2': lambda t, p, ham: p[1].H@p[1],
    # 'entropy': lambda t, p: entropy_subsys(p, dims=[2] * L, sysa=range(L // 2))
}

evo = Evolution(psi0, fn_ham_t, progbar=True, compute=compute,int_small_step=True)

evo.update_to(5)
#%%
plt.plot(evo.results['time'], evo.results['state1'],evo.results['time'], evo.results['state2'])