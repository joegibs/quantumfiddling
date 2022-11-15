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
ep = 1/180 
thp = np.linspace(0, 2*np.pi, length)
thp2 = thp + np.pi
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



#%%%
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
alpha = 1

U,f = gen_u_and_f(alpha,kappa,L,T_rotation)
New_G = np.conjugate(U) @ A_L @ U
ep = 1/180 
thp = np.linspace(0, 2*np.pi, length)
thp2 = thp + np.pi
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
    point_w[i] = np.matmul(np.conjugate(xi), np.matmul(New_G, xi))/(2*np.pi)
    
plt.plot(np.mod(thp2,2*np.pi), point_w.real,'o')
plt.plot(np.mod(thp2,2*np.pi),vonmises(kappa).pdf(thp))
