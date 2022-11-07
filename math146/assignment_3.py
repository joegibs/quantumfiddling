# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
from scipy.stats import vonmises as vonmises
import scipy.integrate as integrate
from scipy.special import jv as bessel
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
T=T_rotation
f = f_vonmis
kwargs={"alpha": np.pi / 2, "kappa": 1}
angle = np.linspace(0, 2 * np.pi, 300)

radius = 1

yarg, colors = composition(angle, T, f, **kwargs)
#%%
"""
define some functions for discritization
"""
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
        
    return U,f
#%%
"""
Generate the A_L matrix
"""
L=2
kappa=3
def gen_A_L(L,kappa):
    #Toeplitz matrix of the expansion coefficients
    A_L=np.zeros((2*L+1,2*L+1))
    for ij in range(-2*L,2*L+1):
        ans = bessel(np.abs(ij), kappa) / bessel(0, kappa)
        A_L+=np.diag([ans]*(2*L-abs(ij)+1),ij)
    return A_L
A_L=gen_A_L(L,kappa)

eigval,eigvecs=np.linalg.eigh(A_L)
#%%
"""
Next, make histograms of the set of eigenvalues aj of AL (i.e., the spectrum of the operator
AL ∈B(HL)) for various values of L
"""
# def plot_loghist(x, bins):
#   hist, bins = np.histogram(x, bins=bins)
#   # print(bins)
#   logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
#   plt.hist(x, bins=logbins)
#   plt.xscale('log')

# plot_loghist(eigval,100)

plt.hist(eigval,bins=100)
#%%
"""
Compare your histograms with histograms for the values f(θ)
of f (i.e., the spectrum of f as an element of L∞(μ)) obtained by sampling θ from the Lebesgue
measure on S1.
"""
samp = np.random.rand(500)*2*np.pi
def f_theta(theta):
    thing = np.exp(kappa*np.cos(theta)) / bessel(0, kappa)
    return thing
thing = f_theta(samp)
plt.hist(thing, bins= 100)
#%%
"""
In addition, for representative L and j, reconstruct the functions ηj = PL
l=−L uljφl
on S1 associated with eigenvector uj. What do you observe as L increases
"""
length = 200
theta = np.linspace(0,2*np.pi,length)
tot = np.zeros(length,dtype='complex')

for i,eigvec in enumerate(eigvecs[:,:]):
    expan = np.zeros_like(theta, dtype="complex")
    # print(i,eigvec)
    for j,val in enumerate(eigvec):
        # print(j, val)
        expan += val * phi_i(lplus1_space(j), theta)
    tot=tot+expan
    # if i==0:
    #     plt.plot(theta,np.real(expan))
    #     plt.plot(theta,f_vonmis(theta,**{'kappa':kappa}))

#%%
plt.plot(theta,np.real(tot))
plt.plot(theta,f_vonmis(theta,**{'kappa':kappa}))

#%%
"""
Write code that implements RL for the von Mises kernel from (3) or another kernel of your
choice. 
"""
