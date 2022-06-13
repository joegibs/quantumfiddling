# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 13:19:49 2022

@author: jogib
"""
import itertools
from operator import add
import numpy as np
from quimb import *
import matplotlib.pyplot as plt
#%%
class circuit():
    '''
    contains the bits and mixing functions

    Returns
    -------
    None.

    '''
    def __init__(self, num_elems):
        self.num_elems = num_elems 
        
        ''' this need to be updated for different inits'''
        self.dop = computational_state("".join(['0' for x in range(self.num_elems)]),
                                       sparse=False)
        self.dims = [2] * self.num_elems
        
        self.step_num = 0 
        self.pairs = []
        self.target = 0
        
        
        
    def gen_step(self,eoo):
        self.gen_pairs(eoo)
        for i in self.pairs:
            self.do_operation(i)
        self.step_num = self.step_num + 1
        
    def gen_pairs(self,eoo):
        self.pairs = []
        #some control over where indexing starts
        if eoo: 
            i =0
        else:
            i=1
        # get pairs
        while i+2<= self.num_elems:
            self.pairs.append([i,i+1])
            i=i+2
            
    def do_operation(self,pair):
          # Needs some work
          had = ikron(hadamard(),[2]*self.num_elems,0)
          step1 = had@ self.dop
          self.dop = self.dop#pkron(CNOT(),[2]*self.num_elems,pair)@step1
        
    def mutinfo(self,target = 0):
        #this is mem bad
        arr=[ptr(self.dop, dims=self.dims, keep=[target,x]) for x in range(np.size(self.dims))]
        mi = [mutinf(arr[x] if x != target else purify(arr[x])) for x in range(np.size(arr))]
        return mi
        
    def mut_info_array_gen(self, num_steps, target=0):
        arr=[]
        for i in range(num_steps):
            circ.gen_step(1-i%2)
            arr.append(circ.mutinfo(target))
        return arr
            
            
    def print_arr(self):
        print()
        print(self.dop)
#%%
circ=circuit(11)
arr = circ.mut_info_array_gen(5,0)



#%%
plt.imshow(arr)
plt.title("Mutual Information with site 1")
plt.ylabel("step number")
plt.xlabel("site number")
plt.colorbar()
