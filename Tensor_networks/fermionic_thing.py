import quimb as qu
import quimb.tensor as qtn
import matplotlib.pyplot as plt
import numpy as np
from quimb import *
import itertools

from opt_einsum import contract
#%%
def corrlate(psit, A, i, j, B=None, **expec_opts):
    if B is None:
        B = A
    bra = psit.H

    pA = psit.gate(A, i, contract=True)
    cA = qtn.expec_TN_1D(bra, pA, **expec_opts)

    pB = psit.gate(B, j, contract=True)
    cB = qtn.expec_TN_1D(bra, pB, **expec_opts)

    pAB = pA.gate_(B, j, contract=True)
    cAB = qtn.expec_TN_1D(bra, pAB, **expec_opts)

    return cAB 
#%%
'''more or less correct structure
need to check ham
correlation func'''
# need this one thing
cT = np.array([[0,0],[1,0]])
c = np.array([[0,1],[0,0]])
fT = np.array([[0,0],[1,0]])
f = np.array([[0,1],[0,0]])
cT=1/(np.sqrt(2))*(f+fT)
c= (-fT+f)/(1j*np.sqrt(2))

L=10
J=1
Jp=0.05*J

st=''
for i in range(L):
    if i%2:
        st+='0'
    else:
        st+='1'
        
mps = qtn.MPS_computational_state(st)
# mps = qtn.MPS_rand_state(L,bond_dim=2)
H2 = {(i, i + 1): J*(kron(cT,c) + kron(c,cT)) for i in range(1,L-1)}
H2[(0,1)]=Jp*(kron(cT,c)+kron(c,cT))
# H2[(L-1,1)]=kron(cT,c)+kron(c,cT)
H1 = {i: cT@c for i in range(L)}

H = qtn.LocalHam1D(L=L, H2=H2,H1 = H1,cyclic=0)
# H = qtn.ham_1d_heis(L)

tebd = qtn.TEBD(mps, H)

# Since entanglement will not grow too much, we can set quite
#     a small cutoff for splitting after each gate application
tebd.split_opts['cutoff'] = 1e-12

# times we are interested in
ts = np.linspace(0, 5, 101)

corr_b = []
# range of bonds, and sites
js = np.arange(0, L)
bs = np.arange(1, L)


for psit in tebd.at_times(ts, tol=1e-3):

    corr = []
    # there is one more site than bond, so start with mag
    #     this also sets the orthog center to 0

    for j in range(1, L,10):
    #     # after which we only need to move it from previous site
    #     # be_b += [psit.entropy(j, cur_orthog=j)]
    #     corr += [1/Jp*np.abs(corrlate(psit,cT,0,j,B=c))]
    # corr_b += [corr]
#%%
for j in range(1, L):
    corr += [1/Jp*np.abs(corrlate(psit,cT,0,j,B=c))]
#%%
x=np.linspace(1,L,L-1)*Jp
plt.plot(x,corr)
#%%
plt.figure(figsize=(12, 7))

# plot the entropy
ax1 = plt.subplot(131)
plt.pcolormesh(bs, ts, (corr_b))
plt.setp(ax1.get_yticklabels(), visible=False)
plt.set_cmap('viridis'), plt.colorbar()
plt.title('correlation')
plt.xlabel('Bond')


plt.show()

#%%
'''

'''
#%%
L=4
XX = pauli('X') & pauli('X')
YY = pauli('Y') & pauli('Y')


mps = qtn.MPS_rand_state(L,bond_dim=2,cyclic=0)
H = qtn.LocalHam1D(L=L, H2=XX + YY,cyclic = 0)
builder = qtn.SpinHam1D(1/2,cyclic=0)
builder += 0.5, '+', '-'
builder += 0.5, '-', '+'
builder += 1.0, 'Z', 'Z'

H=builder.build_mpo(L)

tebd = qtn.TEBD(mps, H)

# Since entanglement will not grow too much, we can set quite
#     a small cutoff for splitting after each gate application
tebd.split_opts['cutoff'] = 1e-12

# times we are interested in
ts = np.linspace(0, 10, 101)

corr_b = []
# range of bonds, and sites
js = np.arange(0, L)
bs = np.arange(1, L)


for psit in tebd.at_times(ts, tol=1e-3):

    # corr = []
    # there is one more site than bond, so start with mag
    #     this also sets the orthog center to 0

    # for j in range(1, L):
        # after which we only need to move it from previous site
        # be_b += [psit.entropy(j, cur_orthog=j)]
    corr += [psit.correlation(identity(2),0,L-1).real]
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