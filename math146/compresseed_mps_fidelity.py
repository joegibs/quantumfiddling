#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 09:50:48 2024

@author: joeg
"""

%config InlineBackend.figure_formats = ['svg']

from quimb.tensor import *
from quimb import *
import numpy as np
#%%

def generate_initial_vector(n=10,n_gates=50):
    # the initial state
    cyclic = False
    psi = MPS_computational_state("1"+"0"*(n-1),cyclic=cyclic, tags='KET', dtype='complex128')    
    # the gates
    gates = [rand_uni(4) for _ in range(n_gates)]
    u_tags = [f'U{i}' for i in range(n_gates)]
    
    for U, t in zip(gates, u_tags):
        # generate a random coordinate
        i = np.random.randint(0, n - int(not cyclic))
        
        # apply the next gate to the coordinate
        #     propagate_tags='sites' (the default in fact) specifies that the
        #     new gate tensor should inherit the site tags from tensors it acts on
        psi.gate_(U, where=[i, i + 1], tags=t, propagate_tags='sites')
        
    return psi.to_dense()

from scipy.linalg import logm
def log2m(A):
    """Compute matrix logarithm base 2"""
    return np.log2(np.exp(1)) * logm(A)

def KL(a, b):
    """relative entropy of a wrt b"""
 
    a = np.matrix(a)
    b = np.matrix(b)

    if min(a.shape)==1 and min(b.shape)==1:
        A = np.diagflat(a)
        B = np.diagflat(b)
    else:
        A = a
        B = b
            
    return np.trace( A @ log2m(A) ) - np.trace( A @ log2m (B) )

def KL(a,b):
    tot=0
    for i in range(len(a)):
        A = a[i]*np.conjugate(a[i])
        B = b[i]*np.conjugate(b[i])
        tot+= A*np.log2(A/B)
    return tot

#%%
n=10
init= generate_initial_vector(n,100)

full =MatrixProductState.from_dense(init,[2]*n)
compressed = MatrixProductState.from_dense(init,[2]*n,max_bond=2)
fidelity(compressed.to_dense(),full.to_dense(),squared=True)



#%%
n=10
fids=[]
for trial in range(0,10):
    fids_trial=[]
    init= generate_initial_vector(n,1000)
    full =MatrixProductState.from_dense(init,[2]*n)
    
    bonds = range(32,0,-1)
    for max_bond in bonds:
        compressed = MatrixProductState.from_dense(init,[2]*n,max_bond=max_bond)
        fids_trial.append(fidelity(compressed.to_dense(),full.to_dense(),squared=True))
    fids.append(fids_trial)
    
plt.plot(bonds,np.mean(fids,axis=0))
