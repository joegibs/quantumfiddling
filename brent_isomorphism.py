# -*- coding: utf-8 -*-
"""
Created on Fri Dec 23 14:24:45 2022

@author: jogib
"""

import networkx as nx
from networkx.algorithms import isomorphism
#%%
a=np.array([[0., 1., 1.],
       [1., 0., 0.],
       [1., 0., 0.]])
b=np.array([[0., 1., 0.],
       [1., 0., 1.],
       [0., 1., 0.]])
test = nx.is_isomorphic(nx.from_numpy_array(a), nx.from_numpy_array(b))
#%%
import numpy as np
from itertools import product
#%%
N=3
rows_g = product([0,1],repeat=N)
rows= [i for i in rows_g]
#%%
matricies_g = product(rows,repeat=N)
matricies = [i for i in matricies_g if np.linalg.det(i)%2==1]

#%%
def remove_iso(a):
    arr = [a[0]]
    for mat in a[1:]:
        print(arr)
        chk = nx.from_numpy_array(np.array(mat),create_using=nx.DiGraph)
        app = 0
        for i in arr:
            if nx.is_isomorphic(chk, nx.from_numpy_array(np.array(i),create_using=nx.DiGraph)):
                app = 1
                pass
        if not app:
            arr.append(mat)
    return arr
# def remove_iso(a,j):
#     if j >= len(a):
#         return a
#     chk = nx.from_numpy_array(np.array(a[j]),create_using=nx.DiGraph)
#     i = j+1
#     while len(a) != i:
#         if  nx.is_isomorphic(chk, nx.from_numpy_array(np.array(a[i]),create_using=nx.DiGraph)):
#             a.pop(i)
#         else: 
#             i+=1
#     return a
#%%
val = len(matricies)
stop_val=1
chk = remove_iso(matricies)
j=1
#%%
while val != 8:
    print(len(chk))
    val = len(chk)
    j+=1
    if j>len(chk):
        break
    chk = remove_iso(chk,j)
    stop_val = len(chk)
print(chk)
#%%
test = nx.is_isomorphic(nx.from_numpy_array(np.array(chk[0]),create_using=nx.Graph), nx.from_numpy_array(np.array(chk[1]),create_using=nx.Graph))
#%%
#%%
def check_rows(arr_1,arr_2):
    return frozenset(chk[0]) == frozenset(chk[7])
def check_cols(arr_1,arr_2):
    return frozenset(chk[0]) == frozenset(chk[7])

#%%
G = nx.DiGraph()
G.add_node(1)
G.add_node(2)
G.add_edge(1, 1)
G.add_edge(2, 2)
G.add_edge(1,2)
G2 = nx.DiGraph()
G2.add_node(1)
G2.add_node(2)
G2.add_edge(1, 1)
G2.add_edge(2, 2)
G2.add_edge(2, 1)
ahh = nx.to_numpy_array(G2)