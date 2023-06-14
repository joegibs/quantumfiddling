# -*- coding: utf-8 -*-
"""
Created on Fri May 26 10:19:05 2023

@author: jogib
"""

#%%
from pylab import *

from pylab import *
L = 20
p = 0.4
z = rand(L,L)
m = z<p
imshow(m, origin='lower')
show()
#%%
from scipy import ndimage
L = 100
p = 0.5
z = rand(L,L)
m = z<p
lw, num = ndimage.label(m)

b = arange(lw.max() + 1)
shuffle(b[1:])
shuffledLw = b[lw]
area = ndimage.sum(m, lw, index=arange(lw.max() + 1))
areaImg = area[lw]
imshow(shuffledLw, origin='lower',cmap='gist_heat_r')
colorbar()
# imshow(shuffledLw, origin='lower',cmap='gnuplot')#'gist_heat_r'

#%% image for different probabilities
L = 50
pv = [0.2,0.3,0.4,0.5,0.6,0.7]
z = rand(L,L*2)
for i in range(len(pv)):
    p = pv[i]
    m = z<p
    lw, num = ndimage.label(m)
    b = arange(lw.max() + 1)
    shuffle(b[1:])
    shuffledLw = b[lw]
    area = ndimage.sum(m, lw, index=arange(lw.max() + 1))
    areaImg = area[lw]
    subplot(2,3,i+1)
    tit = 'p='+str(p)
    imshow(shuffledLw, origin='lower')
    title(tit)
    axis()
    
#%%
lw,num = ndimage.label(z)
perc_x = intersect1d(lw[0,:],lw[-1,:])
perc = perc_x[where(perc_x>0)]


L = 100
colors = ['tab:orange','r','k']
ci=0
for L in [10,50,100]:
    p = linspace(0.0,1.0,200)
    nx = len(p)
    Ni = zeros(nx)
    P = zeros(nx)
    N = 1000
    for i in range(N):
        z = rand(L,L)
        for ip in range(nx):
            m = z<p[ip]
            lw, num = ndimage.label(m)
            perc_x = intersect1d(lw[0,:],lw[-1,:])
            perc = perc_x[where(perc_x>0)]
            if (len(perc)>0):
                Ni[ip] = Ni[ip] + 1
    Pi = Ni/N
    plot(p,Pi,colors[ci])
    ci+=1
    xlabel('Probability of open site $p$')
    ylabel('Probability of an spanning cluster $\Pi(p)$')
    legend([10,50,100])
