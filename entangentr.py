# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 15:01:29 2023

@author: jogib
"""
import numpy as np
from quimb import *
#%%
R=np.array([[1,1],[1,-1]])
state=qu(bell_state(0),qtype='dop')

print(entropy(ptr(state,[2]*2,1)))

xstate = 1/4*kron(R,R)@state@kron(R,R).H

print(entropy(ptr(xstate,[2]*2,1)))
#%%
Rbell = 1/np.sqrt(2)*qu(np.array([[1,0,0,1],[0,1,1,0],[0,1,-1,0],[1,0,0,-1]]),qtype='dop')
bstate=Rbell@state@Rbell.H
print(entropy(ptr(bstate,[2]*2,1)))

#%%
def ent_vec(x):
    ent = sum([i * np.log(i) for i in x])
    return ent


