#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 19:03:07 2024

@author: joeg
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
from scipy.stats import vonmises as vonmises
from scipy.special import iv as bessel

from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation

from quimb import *
import quimb.tensor as qtn

#%%
def composition(theta, T, f, **kwargs):
    thetan = T(theta, **kwargs)
    zn = f(thetan, **kwargs)
    return thetan, zn


def T_blank(theta, **kwargs):
    return np.mod(theta, 2 * np.pi)


def T_rotation(theta, **kwargs):
    if "alpha" in kwargs.keys():
        return np.mod((theta + kwargs["alpha"]), 2 * np.pi)
    else:
        return np.mod(theta, 2 * np.pi)


def T_double_map(theta, **kwargs):
    return np.mod(theta * 2, 2 * np.pi)


def f_ang(theta, **kwargs):
    return theta


def f_vonmis(theta, **kwargs):
    if "kappa" in kwargs.keys():
        return vonmises(kwargs["kappa"]).pdf(np.mod(theta, np.pi * 2))
    else:
        return vonmises(1).pdf(theta)
#%%
def phi_i(l, theta):
    return np.exp(1j * l * theta)
    # return np.cos(l*theta)+1j*np.sin(l*theta)


def U_ij(i, j, T, X=[0, 2 * np.pi], **kwargs):
    # print(i,j)
    f = lambda x: (np.conjugate(phi_i(i, x)) * phi_i(j, T(x, **kwargs)))
    a = f(np.linspace(X[0], X[1], 300))
    return np.trapz(a, np.linspace(X[0], X[1], 300))

def basis_expansion(g, basis_func, theta):
    expan = np.zeros_like(theta, dtype="complex")
    for ij in np.ndindex(np.shape(g)):
        # print(ij)
        expan += g[ij] * basis_func(lplus1_space(ij[0]), theta)
    return expan
def lplus1_space(i):
    return int(i - L)
def calculate_g(U, f, n=1):
    g = U @ f
    for i in range(1, n):
        g = U @ g
    return g
def gen_u_and_f(alpha,kappa,L,T):
    U = np.zeros((2 * L + 1, 2 * L + 1),dtype='complex')

    """
    calculate U_ij,
    """
    for ij in np.ndindex(np.shape(U)):
        # print(lplus1_space(ij[0]),lplus1_space(ij[1]))
        ans = U_ij(lplus1_space(ij[0]), lplus1_space(ij[1]), T, **{"alpha": alpha})
        U[ij] = (ans) / (np.pi * 2)
    """
    Calculate fl
    """
    f = np.zeros((2 * L + 1, 1))
    for ij in np.ndindex(np.shape(f)):
        # print(abs(lplus1_space(ij[0])))
        ans = bessel(abs(lplus1_space(ij[0])), kappa) / bessel(0, kappa)
        f[ij] = ans
    f= f/np.linalg.norm(f)
    return U,f

#%%

alpha = -2*np.pi/2
kappa = 100
L=8
angle = np.linspace(0, 2*np.pi, 300)

_, colors = composition(angle, T_rotation, f_vonmis, **{"alpha": alpha, "kappa": kappa})

U,f = gen_u_and_f(alpha,kappa,L,T_rotation)
g = calculate_g(U,f)

colorsd = basis_expansion(g, phi_i, angle).real/(np.pi*2)
f=f/np.linalg.norm(f)
colorsf = basis_expansion(f, phi_i, angle).real/(np.pi*2)
fig, axs = plt.subplots(2)
fig.suptitle(f'Discritized vs Continuous transformation, alpha = {alpha}, kappa = {kappa}')

axs[0].plot(angle, colors)
axs[0].set_title('Continuous')
axs[0].set_xlabel('')
axs[1].plot(angle, colorsd)
axs[1].set_title('Discrete')
axs[1].plot(angle, colorsf)
axs[1].set_title('Discrete')

plt.show()

#%% ok continuous to discrete works can i f into an mps
mps_f = qtn.MatrixProductState.from_dense(f,[2*L+1])
mps_U = qtn.MatrixProductOperator.from_dense(U,(2*L+1))