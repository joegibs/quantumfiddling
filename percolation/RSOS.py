# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt

#%%
class surface:
    def __init__(
            self,
            d=1,
            length=5,
            ):
        self.grid=np.zeros([length]*d)
        self.length=length
        self.d=d
    
    def show_grid(self):
        print(self.grid)
    
    def add_height(self,*args):
        self.grid[args]+=1
    
    def rand_growth(self,number):
        for i in range(number):
            self.add_height(*np.random.choice(self.length,self.d))
    def sos_growth(self,gamma,steps,N):
        for count in range(steps):
            i = np.random.choice(np.arange(self.length))
            if i == 0:
                if abs(self.grid[i] - self.grid[i+1])<N and abs(self.grid[i] - self.grid[-1])<N:
                    self.grid[i]+=np.random.choice([1,0],p=[gamma,1-gamma])
            elif i == self.length-1:
                if abs(self.grid[i] - self.grid[0])<N and abs(self.grid[i] - self.grid[i-1])<N:
                    self.grid[i]+=np.random.choice([1,0],p=[gamma,1-gamma])
            else:
                if abs(self.grid[i] - self.grid[i+1])<N and abs(self.grid[i] - self.grid[i-1])<N:
                    self.grid[i]+=np.random.choice([1,0],p=[gamma,1-gamma])
                    
class sos_model:
    def __init__(
            self,
            length=5,
            ):
        self.length = length
        self.sites = np.array([1,0]*length)
        self.active_sites = np.array
    def print_sites(self):
        print(self.sites)
    def sos_growth(self,steps,N=1):
        for count in range(steps):
            i = np.random.choice(np.arange(len(self.sites)))
            if len(self.sites) - 1 ==i:
                if self.sites[i-1]-self.sites[i]==1 and self.sites[0] -self.sites[i]==1:
                    self.sites+=2
            elif self.sites[i-1]-self.sites[i]==1 and self.sites[i+1] -self.sites[i]==1:
                self.sites[i]+=2