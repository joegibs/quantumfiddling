# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 12:53:07 2023

@author: jogib
"""

import numpy as np
from quimb import *
import matplotlib.pyplot as plt

#%%

def funct(*kwargs):
    for a in kwargs:
        print(a)

#%%
vec = np.random.rand(5)
# vec = vec/sum(vec)
vec = np.random.dirichlet(np.ones(5))
# print(ent_vec(vec))

class sites():
    def __init__(self,length=2,init=0):
        self.length=2
        if init:
            self.probs=np.array(qu([[np.random.dirichlet(np.ones(2))] for i in range(self.length)]))
        else:
            self.probs=np.array(qu([[[0,1]] for i in range(self.length)]))
        
    #measures    
    def ent_vec(self,x):
        #S = p log(p)
        ent = -sum([i * np.log2(i) if i != 0 else 0 for i in np.real(x)])
        return ent
    def ent_site(self,site):
        #entropy for one site
        return self.ent_vec(self.probs[site][0])
    def ent_joint(self,sites):
        joint_vec = kron(*[self.probs[i] for i in sites])[0]
        return self.ent_vec(joint_vec)
    def mut_inf(self,sites):
        tot = sum([self.ent_site(i) for i in sites])
        tot-= self.ent_joint(sites)
        return tot
    #operations
    def op_markov(self, sites, rand=1):
        joint_vec = qu(kron(*[self.probs[i] for i in sites])[0],qtype='ket')
        if rand:
            mark = self.gen_markov()
        else:
            mark = self.set_markov()
        fin = mark@joint_vec
        sep = self.sep_joint_vec(fin)
        self.probs[sites]=np.reshape(sep,np.shape(self.probs[sites]))
    def op_meas(self,site):
        choice = np.random.choice([0,1],p=np.real(a.probs[site][0]))
        se =[0,0]
        se[choice] = 1
        self.probs[site][0]=se
        
    #utilities
    # def set_markov(self):
    #     return np.array([[1,0,0,0],[0,1,0,0],[0,0,0.5,0.5],[0,0,0.5,0.5]])
    def set_markov(self):
        return np.array([[0.5,0,0,0.5],[0,0,0,0],[0,0,0,0],[0.5,0,0,0.5]])
    def gen_markov(self):
        return np.transpose(np.random.dirichlet([1]*4, 4))
    def sep_joint_vec(self,vec):
        return np.array([[vec[0]+vec[1],vec[2]+vec[3]],[vec[0]+vec[2],vec[1]+vec[3]]])
    
    #%%
    
a=sites(init=1)
print("before ent ",a.ent_joint([0,1]))
print("ent site 1 ",a.ent_site(1))
print("before mut ",a.mut_inf([0,1]))

for i in range(5):
    a.op_markov([0,1],rand=0)
print("probs sum ",np.sum(a.probs,axis=2))
print("after ent ",a.ent_joint([0,1]))
print("ent site 1 ",a.ent_site(1))
print("after mut ",a.mut_inf([0,1]))
#%%
data = []
data_corr = []
for i in range(1000):
    a=sites()
    for i in range(5):
        a.op_markov([0,1],rand=1)
    data.append(a.ent_joint([0,1]))
    data_corr.append(a.mut_inf([0,1]))
plt.hist([data,data_corr])
plt.legend(["ent","mut"])
#%%
#%%
data = []
data_corr=[]
for i in range(1000):
    a=sites(init=1)
    tst = a.ent_joint([0,1])
    tst_corr = a.mut_inf([0,1])
    for i in range(1):
        a.op_meas(1)
    data.append(tst-a.ent_joint([0,1]))
    data_corr.append(tst_corr-a.mut_inf([0,1]))
plt.hist([data,data_corr])
plt.legend(["diff ent","diff mut"])
#%%
a=sites()

ss=[0,1]
joint_vec = qu(kron(*[a.probs[i] for i in ss])[0],qtype='ket')
mark = np.transpose(a.gen_markov())
fin = mark@joint_vec
sep = a.sep_joint_vec(fin)
a.probs[ss]=np.reshape(sep,np.shape(a.probs[ss]))