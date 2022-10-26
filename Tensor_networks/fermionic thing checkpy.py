

from quimb import *
from quimb.tensor import *

from tqdm import tqdm
#%%
L=1000
H = MPO_ham_heis(L)
# H.draw()
#%%
H2 = {(i, i + 1): J*(kron(cT,c) + kron(c,cT)) for i in range(1,L-1)}
H2[(0,1)]=Jp*(kron(cT,c)+kron(c,cT))
# H2[(L-1,1)]=kron(cT,c)+kron(c,cT)
H1 = {i: cT@c for i in range(L)}

H = qtn.LocalHam1D(L=L, H2=H2,H1 = H1,cyclic=0)
#%%
E_exact = heisenberg_energy(L)
E_exact
#%%
dmrg = DMRG2(H)
#%%
dmrg.TN_energy.draw(color=['_KET', '_HAM', '_BRA'])  # might be slow as uses force repulsion
#%%
dmrg.solve(max_sweeps=40, verbosity=1, cutoffs=1e-6)
#%%
(dmrg.energy - E_exact) / abs(E_exact)

#%%
dmrg.opts

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
cT = np.array([[0,0],[1,0]])
c = np.array([[0,1],[0,0]])
fT = np.array([[0,0],[1,0]])
f = np.array([[0,1],[0,0]])
# cT=1/(np.sqrt(2))*(f+fT)
# c= (-fT+f)/(1j*np.sqrt(2))

corr = np.empty(L)
for j in tqdm(range(1, L)):
    #     # after which we only need to move it from previous site
    corr[j-1] = np.abs(corrlate(dmrg.state,cT,0,j,B=c))

#%%
plt.loglog([i for i in range(1,L+1)],corr.round(7),'o')