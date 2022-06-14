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
    def __init__(self, num_elems, gate = 'bell', init = 'up'):
        '''
        num_elems: number of elements in the chain
        gate: type of gate used
            -"bell" hadamard then cnot to make a bell state
            -"haar" haar random unitary operators
        init: initialization of qqubits
        
        architecture: arrangement of gates
            -"brick": alternating pairs
            -staircase: each step only has one gate
        '''
        
        self.num_elems = num_elems 
        
        ''' this need to be updated for different inits'''
        if init == 'up':
            self.dop = computational_state("".join(['0' for x in range(self.num_elems)]),
                                           sparse=False)
        elif init == 'rand':
            self.dop = rand_product_state(self.num_elems)
            
        self.dims = [2] * self.num_elems
        self.gate=gate
        
        self.step_num = 0 
        self.pairs = []
        self.target = 0
        
        
        
    def gen_step(self,eoo):
        self.gen_pairs(eoo)
        print(self.pairs)
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
        if self.gate == 'bell':
            print(pair[0])
            had = ikron(hadamard(),[2]*self.num_elems,pair[0])
            step1 = had@ self.dop
            self.dop = pkron(CNOT(),[2]*self.num_elems,pair)@step1#qu(pkron(CNOT(),[2]*self.num_elems,pair)@step1,qtype='dop')
            
        elif self.gate == 'haar':
            haar = ikron(rand_uni(2),[2]*self.num_elems,pair)
            self.dop = haar@ self.dop
        
    def mutinfo(self,target = 0):
        #this is mem bad
        arr =[ptr(self.dop, dims=self.dims, keep=[target,x]).round(4) for x in range(np.size(self.dims))]
        mi = [mutinf(arr[x] if x != target else purify(arr[x])) for x in range(np.size(arr))]
        return mi
        
    def mut_info_array_gen(self, num_steps, target=0):
        arr=[circ.mutinfo(target)]
        for i in range(num_steps):
            circ.gen_step(1-i%2)
            arr.append(circ.mutinfo(target))
        return arr
            
#%%
circ=circuit(6,gate='bell',init='up')
arr = circ.mut_info_array_gen(6,1)



#%%
plt.imshow(arr)
plt.title("Mutual Information with site 1")
plt.ylabel("step number")
plt.xlabel("site number")
plt.colorbar()
#%%
dop = computational_state("".join(['0' for x in range(2)]),
                               sparse=False)
had = ikron(hadamard(),[2]*2,0)
step1 = had@ dop
end = pkron(CNOT(),[2]*2,[0,1])@step1
# had2 = ikron(hadamard(),[2]*3,1)
# step2 = had@ end
# en2 = pkron(CNOT(),[2]*3,[1,2])@step2
#%%
step1 = had@ end
end = pkron(CNOT(),[2]*2,[0,1])@step1
#%%
step1 = had@ end
end = pkron(CNOT(),[2]*2,[0,1])@step1