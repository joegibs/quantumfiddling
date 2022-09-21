"""
Functions used to generate and count all possible cliques of anticommuting pauli strings

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
    "X": np.array(np.mat("0+0j,1;1,0")),
    "Y": np.array(np.mat("0,-1j;1j,0")),
}
#%%
# init all the tensors
n = 2
lis = ["X", "Y"]
combos = itertools.product(lis, repeat=n)
names = ["".join(x for x in i) for i in combos]
#%%


def anticommute_check(A: str, B: str):
    """
    Checks if two strings anti commute
    
    Parameters
    ----------
    A : str
    B : str
    Returns
    -------
    bool

    """
    ch_arr = []
    for i in range(len(A)):
        a_t = pauli[A[i]]
        b_t = pauli[B[i]]
        check = np.dot(a_t, b_t) + np.dot(b_t, a_t)
        ch_arr.append(not np.any(check))
    return np.count_nonzero(ch_arr) % 2


def gen_graph(n=2):
    lis = ["X", "Y"]
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


def clique_counter(G, clique_size: int):
    """
    

    Parameters
    ----------
    G : networkx graph
    clique_size : int
        size of clique we are counting

    Returns
    -------
    size : TYPE
        DESCRIPTION.

    """
    size = 0
    check = 0

    if nx.algorithms.approximation.large_clique_size(G) < clique_size:
        return 0
    for i in nx.find_cliques(G):
        if len(i) == clique_size:
            size = size + 1
    return size


#%%
n = 2
G = gen_graph(n)
# G.remove_node("".join("I" for i in range(n)))
print(nx.algorithms.approximation.large_clique_size(G))
print(clique_counter(G, 3))

#%%
# pos = nx.circular_layout(G)
# fig, ax = plt.subplots(figsize=(16, 16))
# node_opts = {"node_size": 500, "node_color": "w", "edgecolors": "k", "linewidths": 2.0}
# nx.draw_networkx_nodes(G, pos, **node_opts)
# nx.draw_networkx_labels(G, pos, font_size=10)
# nx.draw_networkx_edges(G, pos, width=0.2)
# plt.show()

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
    pl = [[pauli[j] for j in i] for i in list(lis)]
    pauli_lis = [[row[i] for row in pl] for i in range(len(list(lis)[0]))]
    return pauli_lis


def multiply_lis(lis):
    list_prod = [functools.reduce(np.dot, i) for i in lis]
    return list_prod


def find_pauli(lis):

    st = ""
    for i in lis:
        if i[0][0] != 0:
            if i[0][0] + i[1][1] == 0:
                st = st + "Z"
            else:
                st = st + "I"
        elif i[0][1] + i[1][0] == 0:
            st = st + "Y"
        else:
            st = st + "X"
    return st

def order_of_repeats(cliques):
    total = {}
    for element in cliques:
        each = {}
        for st in element:
            count = {}
            for s in st:
              if s in count:
                count[s] += 1
                each[s] += 1
              else:
                count[s] = 1
                each[s] = 1
        order = max(list(each.values()))
        if order in total:
            total[order]+=1
        else:
            total[order]=1
    return total,each,count
# #%%
# clique = next(nx.find_cliques(G))

# for i in itertools.combinations(clique, 2):
#     G[i[0]][i[1]]["color"] = "tab:red"
#     G[i[0]][i[1]]["width"] = 2.5

# pos = nx.circular_layout(G)
# fig, ax = plt.subplots(figsize=(16, 16))
# node_opts = {"node_size": 500, "node_color": "w", "edgecolors": "k", "linewidths": 2.0}
# nx.draw_networkx_nodes(G, pos, **node_opts)
# nx.draw_networkx_labels(G, pos, font_size=10)
# # nx.draw_networkx_edges(G, pos, width=.01)


# edge_colors = [edgedata["color"] for _, _, edgedata in G.edges(data=True)]
# edge_width = [edgedata["width"] for _, _, edgedata in G.edges(data=True)]
# nx.draw_networkx_edges(G, pos, width=edge_width, edge_color=edge_colors)
# plt.show()

#%%sub cliques
#%% cliques need to fix

# multiply the rest
# gen_paulis
# multiply
# find pauli return to 'x''y'z' etc

#%%
"""
7->5
get all combos of size 3
generate 5 cliques
"""

num_cliq = 0
check = 0

n = 50
G = gen_graph(n)
# G.remove_node("".join("I" for i in range(n)))
# print(nx.algorithms.approximation.large_clique_size(G))
# print(clique_counter(G, 3))
#%%
cliques = []
for i in nx.find_cliques(G):
    # if len(i) == 5:  # max_xs:
    cliques.append(set(i))
print(len(cliques))
#%%
cliques = []
max_xs = nx.algorithms.approximation.large_clique_size(G)
for i in nx.find_cliques(G):
    if len(i) == 3:  # max_xs:
        cliques.append(set(i))

basis_5 = set()
count = 0
for i in cliques:
    # pick 4
    for j in itertools.combinations(i, 4):
        count += 1
        chk = copy.copy(i)
        for k in j:
            chk.remove(k)
        se = set(j)
        se.add(find_pauli(multiply_lis(gen_paulis(chk))))

        # if not any([se.issubset(b) for b in basis_sets]):
        basis_5.add(frozenset(se))

#%%
"""
5->3
"""
basis_3 = set()
count = 0
for i in basis_5:
    # pick 4
    for j in itertools.combinations(set(i), 2):
        count += 1
        chk = copy.copy(set(i))
        for k in j:
            chk.remove(k)
        se = set(j)
        se.add(find_pauli(multiply_lis(gen_paulis(chk))))

        # if not any([se.issubset(b) for b in basis_sets]):
        basis_3.add(frozenset(se))

#%% dummy position
count = 0
sets = []
for i in basis_5:
    Icount = 0
    if "I" not in list(i)[0]:
        pass

    else:
        pos = list(i)[0].index("I")
        arr = [1 if x[pos] == "I" else 0 for x in list(i)]
        if all(arr):
            sets.append(i)
            count += 1
#%% 
def add_triang(cliques):
    new_triang=set()
    n=len(list(cliques)[0])
    for tria in cliques:
        trian= list(tria)
        for i in range(n+1):
            st=frozenset([x[:i]+'I'+x[i:] for x in trian])
            new_triang.add(st)
            for p in {'X','Y','Z'}:
                one = trian[0][:i]+'I'+trian[0][i:]
                two = trian[1][:i]+p+trian[1][i:]
                three=trian[2][:i]+p+trian[2][i:]
                st=frozenset([one,two,three])
                new_triang.add(st)
                one = trian[0][:i]+p+trian[0][i:]
                two = trian[1][:i]+'I'+trian[1][i:]
                three=trian[2][:i]+p+trian[2][i:]
                st=frozenset([one,two,three])
                new_triang.add(st)
                one = trian[0][:i]+p+trian[0][i:]
                two = trian[1][:i]+p+trian[1][i:]
                three=trian[2][:i]+'I'+trian[2][i:]
                st=frozenset([one,two,three])
                new_triang.add(st)
    return new_triang
        
