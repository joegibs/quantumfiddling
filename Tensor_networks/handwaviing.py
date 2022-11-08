# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 13:01:27 2022

@author: jogib
"""

import quimb as qu
import quimb.tensor as qtn
import matplotlib.pyplot as plt
import numpy as np
from quimb import *
import itertools

from opt_einsum import contract
#%%
b = qu.bell_state("psi-")
k = qu.bell_state("psi+")
mat = qu.pauli("X") & qu.pauli("Y")

#%%
"""
from example
"""
data = qu.bell_state("psi-").reshape(2, 2)
inds = ("k0", "k1")
tags = ("KET",)

ket = qtn.Tensor(data=data, inds=inds, tags=tags)
ket
# ket.draw()
X = qtn.Tensor(qu.pauli("X"), inds=("k0", "b0"), tags=["PAULI", "X", "0"])
Y = qtn.Tensor(qu.pauli("Y"), inds=("k1", "b1"), tags=["PAULI", "Y", "1"])

bra = qtn.Tensor(qu.bell_state("psi+").reshape(2, 2), inds=("b0", "b1"), tags=["BRA"])
#%% evaluate things
TN = ket.H & X & Y & bra
print(TN)
TN.draw(color=["KET", "PAULI", "BRA"], figsize=(4, 4), show_inds="all")
TN.trace("b1", "b0")
#%%
"""
simple bell state prep
"""
init = qtn.Tensor(
    qu.computational_state("01").reshape(2, 2), inds=("k0", "k1"), tags=("KET", "INIT",)
)
had_on_1 = qtn.Tensor(qu.hadamard(), inds=("k0", "b0"), tags=("HADAMARD", "0"))
ident_on_2 = qtn.Tensor(qu.identity(2), inds=("k1", "b1"), tags=("IDENTITY", "1"))
cn = qtn.Tensor(qu.CNOT(), inds=("b0", "b1"), tags=("CNOT", "0", "1"))

TN = cn & had_on_1 & ident_on_2 & init
print(TN)
#%%
N = 2
circ = qtn.Circuit(N)
circ.apply_gate("H", 0)
circ.apply_gate("CNOT", 0, 1)
#%%

# init two quibit
circ = qtn.Circuit(N=2, psi0=qtn.MPS_computational_state("01"))

#%% section 2.3
L = 5

# create the nodes, by default just the scalar 1.0
tensors = [qtn.Tensor() for _ in range(L)]

for i in range(L):
    # add the physical indices, each of size 2
    tensors[i].new_ind(f"k{i}", size=2)


tensors[0].new_bond(tensors[1], size=7)
tensors[0].new_bond(tensors[2], size=7)
tensors[1].new_bond(tensors[2], size=7)
tensors[1].new_bond(tensors[3], size=7)
tensors[2].new_bond(tensors[3], size=7)
tensors[3].new_bond(tensors[4], size=7)

mps = qtn.TensorNetwork(tensors)
mps.draw()

#%% contractions
L = 6

# create the nodes, by default just the scalar 1.0
tensors = [qtn.Tensor(tags=[str(i)]) for i in range(L)]

for i in range(L):
    # add the physical indices, each of size 2
    tensors[i].new_ind(f"k{i}", size=2)

for i in range(0, L - 2, 2):
    tensors[i].new_bond(tensors[i + 1], size=7)
    tensors[i].new_bond(tensors[i + 2], size=7)
tensors[L - 2].new_bond(tensors[L - 1], size=7)
for i in range(1, L - 2, 2):
    print(i, i + 2)
    tensors[i].new_bond(tensors[i + 2], size=7)

mps = qtn.TensorNetwork(tensors)
mps.draw()
mp = mps.contract(tags=["0", "2"])
mp.draw()
#%%
"""
Problems
"""
# Section 1
#%% 1
a = np.array([[b ** 2 - 2 * a for b in range(3)] for a in range(3)])
b = np.array(
    [[[(-((3) ** a)) * g + d for a in range(3)] for g in range(3)] for d in range(3)]
)
c = np.array([[e for i in range(3)] for e in range(3)])
d = np.array([[[b * g * e for b in range(3)] for g in range(3)] for e in range(3)])

tensors = [
    qtn.Tensor(data=a, inds=["a", "b"], tags=["A"]),
    qtn.Tensor(data=b, inds=["d", "g", "a"], tags=["B"]),
    qtn.Tensor(data=c, inds=["e", "d"], tags=["C"]),
    qtn.Tensor(data=d, inds=["e", "g", "b"], tags=["D"]),
]

mps = qtn.TensorNetwork(tensors)
mps.draw()

print(mps ^ ...)
# print(a1@b1@c1@d1)
"""
Indicies are troublesome be very careful
"""
#%%
print(contract("ab,gda,ed,bge", a,b,c,d))
print(contract("ab,dga,ed,bge", a,b,c,d))

#%%
s = 0
for a in range(3):
    for b in range(3):
        for d in range(3):
            for e in range(3):
                for g in range(3):
                    s = s + (b ** 2 - 2 * a) * ((-((3) ** a)) * g + d) * g * b * e ** 2
print(s)
#%% simpler example to indicie check

a = np.array([[a + b for a in range(3)] for b in range(3)])
b = np.array([[a + 3 ** b for a in range(3)] for b in range(3)])
tensors = [
    qtn.Tensor(data=a, inds=["a", "b"], tags=["A"]),
    qtn.Tensor(data=b, inds=["a", "b"], tags=["B"]),
]
mps = qtn.TensorNetwork(tensors)
mps.draw()
print(mps ^ ...)

s = 0
for a in range(3):
    for b in range(3):
        s = s + (b + a) * (3 ** b + a)
print(s)
#%%
a = np.array([[a + b for a in range(3)] for b in range(3)])
b = np.array([[c + 3 ** b for c in range(3)] for b in range(3)])
b = np.array([[a + 3 ** c for a in range(3)] for c in range(3)])
tensors = [
    qtn.Tensor(data=a, inds=["a", "b"], tags=["A"]),
    qtn.Tensor(data=b, inds=["c", "b"], tags=["B"]),
    qtn.Tensor(data=b, inds=["a", "c"], tags=["C"]),
]
mps = qtn.TensorNetwork(tensors)
mps.draw()
print(mps ^ ...)

s = 0
for a in range(3):
    for b in range(3):
        for c in range(3):
            s = s + (b + a) * (3 ** b + c) * (a + 3 ** c)
print(s)
#%% 2 table

a = np.array([[b ** 2 - 2 * a for b in range(3)] for a in range(3)])
b = np.array(
    [[[(-((3) ** a)) * g + d for a in range(3)] for g in range(3)] for d in range(3)]
)

L = 5
tensors = [qtn.Tensor(data=a, inds=["b", "a"], tags=["A"])]


mps = qtn.TensorNetwork(tensors)
print(mps.all_inds())
mps.add_tensor(qtn.Tensor(data=a, inds=["g", "a"], tags=["B"]))
mp = mps ^ ...
print(mp.inds)
mps.add_tensor(qtn.Tensor(data=b, inds=["g", "d", "e"], tags=["C"]))
mp = mps ^ ...
print(mp.inds)
mps.add_tensor(qtn.Tensor(data=a, inds=["b", "d"], tags=["D"]))
mp = mps ^ ...
print(mp.inds)
mps.add_tensor(qtn.Tensor(data=[1], inds=["e"], tags=["e"]))
mp = mps ^ ...
print(mp)

"""
data doesn't really matter just filling so that i have appropriatly sized arrays'
"""
#%%
"""
3 colorings
not sure what the data looks like exactly 
The value of e is 1 if and only if all indices are identical, and zero otherwise,
whilst n has value 1 if and only if all legs differ and 0 otherwise
"""
L = 13

# create the nodes, by default just the scalar 1.0



def flat_num(arr):
    return 9*arr[0]+arr[1]+3*(arr[2]-1)

for i in range(6):
    dat=np.ones_like(tensors[i].data)
    arr = [0 if x not in (0,13,26) else 1 for x in range(np.product(np.shape(tensors[i].data))) ]
    tensors[i].modify(data = np.reshape(arr,np.shape(tensors[i].data)))
    # sizes = [tensors[i].ind_size(x) for x in tensors[i].inds]
    
#     tensors[i].modify(data=dat)
for i in range(6,13):
    dat=np.ones_like(tensors[i].data)-(np.ones_like(tensors[i].data)-np.identity(3))

    arr = [0 if x in (0,13,26) else 1 for x in range(np.product(np.shape(tensors[i].data))) ]
    tensors[i].modify(data = dat)#np.reshape(arr,np.shape(tensors[i].data)))


mps = qtn.TensorNetwork(tensors)
mps.draw(color=['e','n'])
print(mps^...)
#%%
"""
different coloring alg from https://arxiv.org/pdf/1708.00006.pdf
"""
L = 2

def levi_cevita_tensor(dim):   
    arr=np.zeros(tuple([dim for _ in range(dim)]))
    for x in itertools.permutations(tuple(range(dim))):
        mat = np.zeros((dim, dim), dtype=np.int32)
        for i, j in zip(range(dim), x):
            print(i,j, x)
            mat[i, j] = 1
        arr[x]=int(np.linalg.det(mat))
    return arr
#%%
# create the nodes, by default just the scalar 1.0

tensors = [qtn.Tensor(tags='e') for i in range(L)]

tensors[0].new_bond(tensors[1], size=6)
tensors[0].new_bond(tensors[1], size=6)
tensors[0].new_bond(tensors[1], size=6)
tensors[0].new_bond(tensors[1], size=6)
tensors[0].new_bond(tensors[1], size=6)
tensors[0].new_bond(tensors[1], size=6)

for i in range(2):
    
    tensors[i].modify(data = levi_cevita_tensor(6))

mps = qtn.TensorNetwork(tensors)
# mps.draw(color=['e','n'])
print(mps^...)
#%%
a=np.nditer(tensors[3].data)
np.res
for elem in a:
    print(a.iterindex,elem)
#%%
'''levi_cevita'''
def levi_cevita_tensor(dim):
    
    arr=np.zeros(tuple([dim for _ in range(dim)]))
    for x in itertools.permutations(tuple(range(dim))):
        mat = np.zeros((dim, dim), dtype=np.int32)
        for i, j in zip(range(n), x):
            mat[i, j] = 1
        arr[x]=int(np.linalg.det(mat))
    return arr
    
#%%
        
perm=[0,1,3,2]
n = len(perm)
if n != len(set(perm)):  # infer there are repeated elements
    print(0)
else:
    mat = np.zeros((n, n), dtype=np.int32)
    for i, j in zip(range(n), perm):
        print(i,j)
        mat[i, j] = 1
    print(int(np.linalg.det(mat)))
#%%
"""
qi examples
"""
#%%
print(qu.bell_state(3))
#%% 4 split with svd to get to a mps
# create a tensor with 5 legs
inds = ["a", "b", "c", "d", "e"]
t = qtn.rand_tensor([2, 3, 4, 5, 6], inds=inds)
t.draw(initial_layout="kamada_kawai", figsize=(3, 3))
# split the tensor, by grouping some indices as 'left'

tn = t.split(["a"])
tn.draw(figsize=(3, 3))
tb = tn.split(["b"])
tb.draw(figsize=(3, 3))

#%% 4.2
a = qu.qu([1], qtype="dop")
b = qu.qu([0], qtype="dop")
print(a @ b)
a = qu.qu([[1, 0], [0, 0]], qtype="dop")
b = qu.qu([[0, 0], [0, 0]], qtype="dop")
print(a @ b)

print(a.H @ b)

#%%
a = np.array([[1, 0], [0, 0]])
b = np.array([[0, 0], [0, 1]])
mps= qtn.tensor_1d.MatrixProductState([a,b])

#%%
L = 2

# create the nodes, by default just the scalar 1.0
tensors = [qtn.Tensor() for _ in range(L)]

for i in range(L):
    # add the physical indices, each of size 2
    tensors[i].new_ind(f"k{i}", size=2)

    # add bonds between neighbouring tensors, of size 7
    if i != L - 1:
        tensors[i].new_bond(tensors[(i + 1) % L], size=1)
a = qu.qu([1], qtype="dop")
b = qu.qu([0], qtype="dop")
tensors[0].modify(data=a)
tensors[1].modify(data=b)
mps = qtn.TensorNetwork(tensors)
mps.draw()


#%%
p = qtn.MPS_rand_state(L=20, bond_dim=1)
print(f"Site tags: '{p.site_tag_id}', site inds: '{p.site_ind_id}'")
#%%
"""
MPS
"""
#%%
L = 1

# create the nodes, by default just the scalar 1.0
tensors = [qtn.Tensor() for _ in range(L)]

tensors[0].new_ind(f"1", size=1)
tensors[0].new_ind(f"2", size=1)
tensors[0].new_ind(f"3", size=1)
a= np.array([[[1,0],[0.5,-0.5]],[[0,1],[0.5,-0.5]]])
b = np.array([[1,0],[0,-1]])
#contract("ijk,km->jm",a,b)
tensors[0].modify(data=a)

# tensors[1].modify(data=b)
mps = qtn.TensorNetwork(tensors)
mps.draw()
#%%
a0=np.array([[1,0],[0,1]])
a1=
#%%
a=[[1,0],[0,1]]
b=[[[1,0],[0,1]],[[0,1],[0,0]]]
d=[[1,0],[0,1]]
print(contract("ai,ibk,dk",a,b,d))
#%%3.2 the solution is dumb
L = 6

# create the nodes, by default just the scalar 1.0
tensors = [qtn.Tensor() for _ in range(L)]


tensors[0].new_ind(f"k1", size=2)
tensors[0].new_ind(f"k2", size=2)
tensors[0].new_ind(f"k3", size=2)
tensors[1].new_ind(f"k3", size=2)
tensors[1].new_ind(f"k4", size=2)
tensors[1].new_ind(f"k5", size=2)

tensors[2].new_ind(f"k1", size=2)
tensors[3].new_ind(f"k5", size=2)
tensors[4].new_ind(f"k2", size=2)
tensors[5].new_ind(f"k4", size=2)
    
a0=np.array([[1,0],[0,1]])
a1=np.array([[0,1],[1,0]])
# a=(a1@a0@a1@a1@a1).reshape(2,2,2)
# a=a.reshape(4,4,4)
b0=np.array([0,1])
b1=np.array([1,0])
a= np.array([[[1., 0.],
         [0., 1.]],
 
        [[0., 1.],
         [1., 0.]]])
tensors[0].modify(data=a)
tensors[1].modify(data=a)
tensors[2].modify(data=b0)
tensors[3].modify(data=b0)
tensors[4].modify(data=b1)
tensors[5].modify(data=b0)
tensors[0].modify(tags=['KET'])
tensors[1].modify(tags=['KET'])
tensors[2].modify(tags=['KET'])
# tensors[3].modify(data=b0)
mps = qtn.TensorNetwork(tensors)
mps.draw()
print((mps^...))
#%%
a=[[[1,0],[0,1]],[[1,0],[0,1]]]
b=[[[1,0],[0,1]],[[0,1],[0,0]]]
d=[[1,0],[0,-1]]
contract("jai,ibk,jk",a,b,d)
#%%
L=5
a= np.array([[[1., 0.],[0., 1.]],[[0., 1.],[1., 0.]]])
arr = [a for _ in range(L)]
arr.insert(0,np.array([[1,0],[0,1]]))
arr.append(np.array([[0., 1.],[1., 0.]]))
p=qtn.tensor_1d.MatrixProductState(arr)
tensors=[i for i in p.tensors]

b0=np.array([0,1])
b1=np.array([1,0])
b=[b0,b1]

for i in range(len(tensors)):
    tensors.append(qtn.Tensor())
    tensors[-1].new_ind(f"k{i}", size=2)
for i,j in enumerate("0000000"):
    tensors[i+L+2].modify(data=b[int(j)],tags=(j))
mps = qtn.TensorNetwork(tensors)
mps.draw()

print(mps^...)
#%%
A = qtn.MPO_rand_herm(7, bond_dim=3, tags=['HAM'])
pH = p.H

# This inplace modifies the indices of each to form overlap
p.align_(A, pH)

(pH & A & p).draw(color='HAM', iterations=20, initial_layout='kamada_kawai')
#%%
peps = qtn.PEPS.rand(Lx=3, Ly=3, bond_dim=2, seed=666)
peps.show()
norm = peps.H & peps
norm.draw(color=norm.site_tags, legend=False, figsize=(4, 4))
#%%
%%time
a=norm.contract_boundary(max_bond=32)

#%%
L = 44
zeros = '0' * ((L - 2) // 3)
binary = zeros + '1' + zeros + '1' + zeros
print('psi0:', f"|{binary}>")

psi0 = qtn.MPS_computational_state(binary)
psi0.show()

H = qtn.ham_1d_heis(L)

tebd = qtn.TEBD(psi0, H)

# Since entanglement will not grow too much, we can set quite
#     a small cutoff for splitting after each gate application
tebd.split_opts['cutoff'] = 1e-12

# times we are interested in
ts = np.linspace(0, 80, 101)

mz_t_j = []  # z-magnetization
be_t_b = []  # block entropy
sg_t_b = []  # schmidt gap

# range of bonds, and sites
js = np.arange(0, L)
bs = np.arange(1, L)

for psit in tebd.at_times(ts, tol=1e-3):
    mz_j = []
    be_b = []
    sg_b = []

    # there is one more site than bond, so start with mag
    #     this also sets the orthog center to 0
    mz_j += [psit.magnetization(0)]

    for j in range(1, L):
        # after which we only need to move it from previous site
        mz_j += [psit.magnetization(j, cur_orthog=j - 1)]
        be_b += [psit.entropy(j, cur_orthog=j)]
        sg_b += [psit.schmidt_gap(j, cur_orthog=j)]

    mz_t_j += [mz_j]
    be_t_b += [be_b]
    sg_t_b += [sg_b]
    #%%
plt.figure(figsize=(12, 7))

# plot the magnetization
ax1 = plt.subplot(131)
plt.pcolormesh(js, ts, np.real(mz_t_j), vmin=-0.5, vmax=0.5)
plt.set_cmap('RdYlBu')
plt.colorbar()
plt.title('Z-Magnetization')
plt.xlabel('Site')
plt.ylabel('time [ $Jt$ ]')

# plot the entropy
ax2 = plt.subplot(132, sharey=ax1)
plt.pcolormesh(bs, ts, be_t_b)
plt.setp(ax2.get_yticklabels(), visible=False)
plt.set_cmap('viridis'), plt.colorbar()
plt.title('Block Entropy')
plt.xlabel('Bond')

# plot the schmidt gap
ax3 = plt.subplot(133, sharey=ax1)
plt.pcolormesh(bs, ts, sg_t_b, vmin=0, vmax=1)
plt.setp(ax3.get_yticklabels(), visible=False)
plt.set_cmap('magma_r')
plt.colorbar()
plt.title('Schmidt Gap')
plt.xlabel('Bond')

plt.show()
# #%%
# Lx = 3
# Ly = 3
# zeros = '0' * ((L - 2) // 3)
# binary = zeros + '1' + zeros + '1' + zeros
# print('psi0:', f"|{binary}>")

# psi0 = qtn.PEPS.rand(Lx,Ly,bond_dim=2)
# psi0.show()

# H = qtn.ham_2d_heis(Lx,Ly)

# tebd = qtn.TEBD2D(psi0, H)

# # Since entanglement will not grow too much, we can set quite
# #     a small cutoff for splitting after each gate application
# # tebd.split_opts['cutoff'] = 1e-12

# # times we are interested in
# ts = np.linspace(0, 80, 11)

# mz_t_j = []  # z-magnetization
# be_t_b = []  # block entropy
# sg_t_b = []  # schmidt gap

# # range of bonds, and sites
# js = np.arange(0, L)
# bs = np.arange(1, L)

# for time in ts:
#     tebd.evolve(1,7.2)
#     psit = tebd.get_state()
#     mz_j = []
#     be_b = []
#     sg_b = []

#     # there is one more site than bond, so start with mag
#     #     this also sets the orthog center to 0
#     mz_j += [psit.magnetization(0)]

#     for j in range(1, L):
#         # after which we only need to move it from previous site
#         mz_j += [psit.magnetization(j, cur_orthog=j - 1)]
#         be_b += [psit.entropy(j, cur_orthog=j)]
#         sg_b += [psit.schmidt_gap(j, cur_orthog=j)]

#     mz_t_j += [mz_j]
#     be_t_b += [be_b]
#     sg_t_b += [sg_b]
#%%
L=2
mps = qtn.MPS_computational_state('00',)
XX = pauli('X') & pauli('X')

YY = pauli('Y') & pauli('Y')
Z = pauli('Z')
H = qtn.LocalHam1D(2,H2=XX+YY,H1=Z)

tebd = qtn.TEBD(mps, H)

# Since entanglement will not grow too much, we can set quite
#     a small cutoff for splitting after each gate application
tebd.split_opts['cutoff'] = 1e-12

# times we are interested in
ts = np.linspace(0, 10, 10)

mz_t_j = []  # z-magnetization
be_t_b = []  # block entropy
sg_t_b = []  # schmidt gap

# range of bonds, and sites
js = np.arange(0, L)
bs = np.arange(1, L)

for psit in tebd.at_times(ts, tol=1e-3):
    print(psit.arrays)
#%%
L=3
XX = pauli('X') & pauli('X')

YY = pauli('Y') & pauli('Y')

fT = np.array([[0,0],[1,0]])
f = np.array([[0,1],[0,0]])
cT=1/(np.sqrt(2))*(f+fT)
c= (-fT+f)/(1j*np.sqrt(2))

XX = pauli('X') & pauli('X')
YY = pauli('Y') & pauli('Y')

H = qtn.LocalHam1D(3,H2=c@cT)

#%%
'''more or less correct structure
need to check ham
correlation func'''
# need this one thing
cT = np.array([[0,0],[1,0]])
c = np.array([[0,1],[0,0]])

L=4
J=1
mps = qtn.MPS_rand_state(L,bond_dim=2,cyclic=True)

H2 = {(i, i + 1): kron(cT,c) + kron(cT,c) for i in range(1,L)}
H2[(0,1)]=kron(cT,c)+kron(cT,c)
H2[(L-1,0)]=kron(cT,c)+kron(cT,c)
H1 = {i: cT@c for i in range(L)}

H = qtn.LocalHam1D(L=L, H2=H2,H1 = H1,cyclic=True)

tebd = qtn.TEBD(mps, H)

# Since entanglement will not grow too much, we can set quite
#     a small cutoff for splitting after each gate application
tebd.split_opts['cutoff'] = 1e-12

# times we are interested in
ts = np.linspace(0, 4, 5)

corr_b = []
# range of bonds, and sites
js = np.arange(0, L)
bs = np.arange(1, L)


for psit in tebd.at_times(ts, tol=1e-3):

    corr = []
    # there is one more site than bond, so start with mag
    #     this also sets the orthog center to 0

    for j in range(1, L):
        # after which we only need to move it from previous site
        # be_b += [psit.entropy(j, cur_orthog=j)]
        corr += [psit.correlation(identity(2),0,j).real]
    corr_b += [corr]

    #%%
plt.plot(corr)
#%%
plt.figure(figsize=(12, 7))

# plot the entropy
ax1 = plt.subplot(131)
plt.pcolormesh(bs, ts, corr_b)
plt.setp(ax1.get_yticklabels(), visible=False)
plt.set_cmap('viridis'), plt.colorbar()
plt.title('correlation')
plt.xlabel('Bond')


plt.show()