# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 13:19:02 2022

@author: jogib
"""

import quimb as qu
import quimb.tensor as qtn
import numpy as np
import matplotlib.pyplot as plt
#%%
L = 44
zeros = '0' * ((L - 2) // 3)
binary = zeros + '1' + zeros + '1' + zeros
print('psi0:', f"|{binary}>")

H = qtn.ham_1d_heis(L)
psi0 = qtn.MPS_computational_state(binary)
#%%
tebd = qtn.TEBD(psi0, H)

# Since entanglement will not grow too much, we can set quite
#     a small cutoff for splitting after each gate application
tebd.split_opts['cutoff'] = 1e-12
#%%
# times we are interested in
ts = np.linspace(0, 80, 101)

mz_t_j = []  # z-magnetization
be_t_b = []  # block entropy
sg_t_b = []  # schmidt gap

# range of bonds, and sites
js = np.arange(0, L)
bs = np.arange(1, L)

#%%
# generate the state at each time in ts
#     and target error 1e-3 for whole evolution
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

#%%
import numpy as np, matplotlib.pyplot as plt
from matplotlib import animation
from numpy.random import random as rand

#%matplotlib qt5
    
#%%
fig, ax = plt.subplots()
ims = []
for i in range(60):
    x += np.pi / 15.
    y += np.pi / 20.
    im = ax.pcolor(js, ts[0:2], np.real(mz_t_j[i:i+2]), vmin=-0.5, vmax=0.5,)
    
    if i == 0:
        ax.pcolor(js, ts[0:2], np.real(mz_t_j[i:i+2]), vmin=-0.5, vmax=0.5)
        ax.set_ylim(0,0.4)
        ax.get_yaxis().set_visible(False)
        ax.set_title('Z-Magnetization')
        ax.set_xlabel('Site')
        ax.set_aspect(30)
    ims.append([im])

ani = animation.ArtistAnimation(fig, ims, interval=50, blit=True,
                                repeat_delay=1000)

# To save the animation, use e.g.
#
# ani.save("movie.mp4")
#
# or
#
# writer = animation.FFMpegWriter(
#     fps=15, metadata=dict(artist='Me'), bitrate=1800)
# ani.save("movie.mp4", writer=writer)

plt.show()