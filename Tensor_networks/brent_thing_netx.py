# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 10:44:42 2022

@author: jogib
"""

import networkx as nx
import itertools
import numpy as np
import matplotlib.pyplot as plt
import copy
import functools

#%%
pauli = {
    "I": np.array(np.mat("1+0j,0;0,1")),
    "Z": np.array(np.mat("1+0j,0;0,-1")),
    "X": np.array(np.mat("0+0j,1;1,0")),
    "Y": np.array(np.mat("0,-1j;1j,0")),
}
#%%
# init all the tensors
n = 2
lis = ["X", "Y", "Z", "I"]
combos = itertools.product(lis, repeat=n)
names = ["".join(x for x in i) for i in combos]
#%%


def anticommute_check(A, B):
    ch_arr = []
    for i in range(len(A)):
        a_t = pauli[A[i]]
        b_t = pauli[B[i]]
        check = np.dot(a_t, b_t) + np.dot(b_t, a_t)
        ch_arr.append(not np.any(check))
    return np.count_nonzero(ch_arr) % 2


def gen_graph(n=2):
    lis = ["X", "Y", "Z", "I"]
    combos = itertools.product(lis, repeat=n)
    d = {i: "".join(j for j in x) for i, x in enumerate(combos)}
    G = nx.empty_graph(len(d.values()))
    G = nx.relabel_nodes(G, d)
    # add edges
    for i in itertools.combinations(list(d.values()), 2):

        A = i[0]
        B = i[1]
        if anticommute_check(A, B):
            G.add_edge(i[0], i[1], **{"color": "tab:blue", "width": 0.2})
    return G


#%%
n = 4
G = gen_graph(n)
G.remove_node("".join("I" for i in range(n)))
print(nx.algorithms.approximation.large_clique_size(G))

#%%
pos = nx.circular_layout(G)
fig, ax = plt.subplots(figsize=(16, 16))
node_opts = {"node_size": 500, "node_color": "w", "edgecolors": "k", "linewidths": 2.0}
nx.draw_networkx_nodes(G, pos, **node_opts)
nx.draw_networkx_labels(G, pos, font_size=10)
nx.draw_networkx_edges(G, pos, width=0.2)
plt.show()

#%% cliques need to fix
cliques=[]
num_cliq = 0
check=0

max_xs=nx.algorithms.approximation.large_clique_size(G)
for i in nx.find_cliques(G):
    if len(i) == 9:#max_xs:
        cliques.append(set(i))

# for i in nx.find_cliques(G):
#     if not any([set(i).issubset(b) for b in cliques]):
#         if len(i) == 3:
#             num_cliq = num_cliq +1

print(len(cliques))
#%%
def gen_paulis(lis):
    """
    Parameters
    ----------
    lis : list
        list of the form ['xyz','xxi'...] so that each the first char in
        each element of the list is pauli corresponding to the char.

    Returns
    -------
    paul_lis : lis
        [[pauli(x),pauli(y) pauli(z)],].

    """
    pl=[[pauli[j] for j in i] for i in list(lis)]
    pauli_lis = [[row[i] for row in pl] for i in range(len(list(lis)[0]))]
    return pauli_lis

def multiply_lis(lis):
    asd = [functools.reduce(np.dot, i) for i in lis]
    return asd

def find_pauli(lis):
    
    st=''
    for i in lis:
        if i[0][0]!=0:
            if i[0][0]+i[1][1]==0:
                st=st+'Z'
            else:
                st=st+'I'
        elif i[0][1]+i[1][0]==0:
            st=st+'Y'
        else:
            st=st+'X'
    return st
    
for i in itertools.combinations(cliques[0],2):
    chk = copy.copy(cliques[0])
    for j in i: chk.remove(j)
    print(i,find_pauli(multiply_lis(gen_paulis(chk))))
    
#%%sub cliques
basis_sets=[]
count=0
for i in cliques:
#pick 2
    for j in itertools.combinations(i,2):
        count+=1
        chk = copy.copy(i)
        for k in j: chk.remove(k)
        se=set(j)
        se.add(find_pauli(multiply_lis(gen_paulis(chk))))
        
        # if not any([se.issubset(b) for b in basis_sets]):
        basis_sets.append(se)
#multiply the rest
        #gen_paulis
        #multiply
        #find pauli return to 'x''y'z' etc



#%% count max identity
max_xs=nx.algorithms.approximation.large_clique_size(G)
for b in nx.find_cliques(G):
    if len(b) == 3:
        count = 0;
        for i in b:
            if ({'I'}.issubset(set(i)) and len(set(i))==2):
                count+=1
        if count == 3:
            print(b)
#%%
clique = next(nx.find_cliques(G))

for i in itertools.combinations(clique, 2):
    G[i[0]][i[1]]["color"] = "tab:red"
    G[i[0]][i[1]]["width"] = 2.5

pos = nx.spiral_layout(G)
fig, ax = plt.subplots(figsize=(16, 16))
node_opts = {"node_size": 500, "node_color": "w", "edgecolors": "k", "linewidths": 2.0}
nx.draw_networkx_nodes(G, pos, **node_opts)
nx.draw_networkx_labels(G, pos, font_size=10)
# nx.draw_networkx_edges(G, pos, width=.01)


edge_colors = [edgedata["color"] for _, _, edgedata in G.edges(data=True)]
edge_width = [edgedata["width"] for _, _, edgedata in G.edges(data=True)]
nx.draw_networkx_edges(G, pos, width=edge_width, edge_color=edge_colors)
plt.show()
