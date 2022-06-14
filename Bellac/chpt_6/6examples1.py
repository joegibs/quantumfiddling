# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 07:33:23 2022

@author: jogib
"""

import numpy as np
import random
#%% util funcs
def normalized_2():
    a=random.random()
    b=random.random()
    nor=np.sqrt(np.power(a,2)+np.power(b,2))
    return(a/nor,b/nor)

#%% util consts
had = 1/np.sqrt(2)*np.array([[1,1],[1,-1]])
cnot = [[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]]


#%%
spin_1=np.array(normalized_2())
spin_2=np.array(normalized_2())

#%%
# uu
a=np.array([1,0,0,0])
aa=np.outer(a,a)
had4=np.kron(had,np.identity(2))

ha=np.dot(had4,a)
be = np.dot(cnot,ha)
bep = np.outer(be,be)
haha = np.outer(ha,ha)

np.trace(haha.reshape(2,2,2,2),axis1=0,axis2=2)
np.trace(bep.reshape(2,2,2,2),axis1=0,axis2=2)
# np.trace(ha)

#%%
haha2=np.dot(np.dot(had4,aa),had4)
bep2 =np.dot(np.dot(cnot,haha2),cnot)
np.trace(bep.reshape(2,2,2,2),axis1=0,axis2=2)

#%% 6.5.8
'''init'''
thet = np.array(normalized_2())
ang = np.arccos(thet[0])
thet_perp = np.array([-np.sin(ang),np.cos(ang)])

Psi = 1/np.sqrt(2)*(np.kron(thet,thet_perp)-np.kron(thet_perp,thet))

#rotation about theta in z axis
#exp(i sigmaz theta/2)



#%% 

init = np.kron(np.dot(had,[0,1]),np.dot(had,[1,0]))
np.dot(np.kron(had,had),np.kron([0,1],[1,0]))

def fx(balanced=True):
    if balanced:
        return 1,1
    elif:
        return 1,0
