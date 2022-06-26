# -*- coding: utf-8 -*-
"""
Created on Sat Jun 25 17:55:49 2022

@author: jogib
"""

import quimb as qu

n = 18
H = qu.ham_heis(n, sparse=True).real
psi0 = qu.rand_product_state(n)

#%%

# check normalization
print(qu.expectation(psi0, psi0))
#%%
print(qu.expec(H, psi0))
#%%
print(# find the initial variance in energy
psi0.H @ H @ H @ psi0

)
#%%%%time

en_low, en_high = qu.bound_spectrum(H)

print("Highest energy:", en_high)
print("Lowest energy:", en_low)

#%%
def compute(t, pt):
    """Perform computation at time ``t`` with state ``pt``.
    """
    dims = [2] * n
    lns = [qu.logneg_subsys(pt, dims, i, i + 1) for i in range(n - 1)]
    mis = [qu.mutinf_subsys(pt, dims, i, i + 1) for i in range(n - 1)]
    return t, lns, mis
#%%
evo = qu.Evolution(psi0, H, compute=compute, progbar=True)
#%%
evo.update_to(5)
#%%
%matplotlib inline
import matplotlib.pyplot as plt

ts, lns, mis = zip(*evo.results)

fig, axs = plt.subplots(2, 1, sharex=True)
axs[0].plot(ts, lns);
axs[0].set_title("Logarithmic Negativity")
axs[1].plot(ts, mis);
axs[1].set_title("Mutual Information")

plt.show()
#%%
print(qu.expec(H, evo.pt))