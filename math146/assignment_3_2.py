# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 10:04:48 2022

@author: jogib
"""

#%%
import numpy as np
import matplotlib.pyplot as plt

from scipy.linalg import eigh
from scipy.stats import vonmises as vonmises
import scipy.integrate as integrate
from scipy.special import iv as bessel

#%%
def phi_i(l,theta):
    return np.exp(1j*l*theta)
#%%
L=8
kappa=8
def gen_A_L(L,kappa):
    #Toeplitz matrix of the expansion coefficients
    A_L=np.zeros((2*L+1,2*L+1))
    for ij in range(2*L,-2*L-1,-1):
        ans = bessel(np.abs(ij), kappa) / bessel(0, kappa)
        A_L+=np.diag([ans]*(2*L-abs(ij)+1),ij)
    chk = -(2*(np.indices((2*L+1,2*L+1)).sum(axis=0) % 2))+np.ones((2*L+1,2*L+1))
    return A_L#np.multiply(A_L,chk)
A_L=gen_A_L(L,kappa)

eigval,eigvecs=np.linalg.eigh(A_L)

#%%
length=200
theta=np.linspace(0,1,length)*2*np.pi
tot = np.zeros_like(theta,dtype = 'complex')
ets=[]
for i,eigvec in enumerate(eigvecs.T):
    eta_j = np.zeros_like(theta,dtype='complex')
    for x,j in enumerate(eigvec):
        eta_j += -j.real *  phi_i(x-L,theta)
    ets.append(eta_j.real)

#%% 
"""
Question 1, your code had a fancy slide i just have a variable to select which eigen val
"""
x=16

plt.plot(theta,ets[x])
plt.show()

#%%
"""
Question 2: point wise evaluation plotted with a shift
"""
ep = 10/180 
thp = np.linspace(0, 2*np.pi, length)
thp2 = thp + np.pi/4
point_w = np.zeros(length, dtype = 'complex')

def p_theta_thetap(th, thp):
    return np.sqrt(np.exp(np.cos(th-thp)/ep)/bessel(0, ep**(-1)))

for i in range(length):
    xi_0 = np.zeros(2*L+1, dtype = 'complex')    
    for j in range(-L, L):
        for an in thp:
            xi_0[j+L] += (p_theta_thetap(theta[i], an) * phi_i(j, an))
    
    xi_norm = np.linalg.norm(xi_0)
    xi = xi_0/xi_norm
    point_w[i] = np.matmul(np.conjugate(xi), np.matmul(A_L, xi))/(2*np.pi)
    
plt.plot(np.mod(thp2,2*np.pi), point_w.real)
plt.plot(np.mod(thp2,2*np.pi),vonmises(kappa).pdf(thp))
