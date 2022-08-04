# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 10:44:42 2022

@author: jogib
"""

import networkx as nx
import itertools
import numpy as np
import matplotlib.pyplot as plt

#%%
pauli = {
    "I": np.array(np.mat("1,0;0,1")),
    "Z": np.array(np.mat("1,0;0,-1")),
    "X": np.array(np.mat("0,1;1,0")),
    "Y": np.array(1j * np.mat("0,-1;1,0")),
}
#%%
# init all the tensors
n = 4
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
            G.add_edge(i[0], i[1], **{"color": "tab:blue", "width": 0.05})
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

#%% cliques
for i in nx.find_cliques(G):
    if len(i) == 9:
        print(i)
#%%
clique = next(nx.find_cliques(G))

for i in itertools.combinations(clique, 2):
    G[i[0]][i[1]]["color"] = "tab:red"
    G[i[0]][i[1]]["width"] = 7

pos = nx.circular_layout(G)
fig, ax = plt.subplots(figsize=(32, 32))
node_opts = {"node_size": 500, "node_color": "w", "edgecolors": "k", "linewidths": 2.0}
nx.draw_networkx_nodes(G, pos, **node_opts)
nx.draw_networkx_labels(G, pos, font_size=10)
# nx.draw_networkx_edges(G, pos, width=.01)


edge_colors = [edgedata["color"] for _, _, edgedata in G.edges(data=True)]
edge_width = [edgedata["width"] for _, _, edgedata in G.edges(data=True)]
nx.draw_networkx_edges(G, pos, width=edge_width, edge_color=edge_colors)
plt.show()
