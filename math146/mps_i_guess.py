#%%
import numpy as np
import matplotlib.pyplot as plt

from scipy.linalg import eigh
from scipy.stats import vonmises as vonmises
import scipy.integrate as integrate
from scipy.special import iv as bessel

from quimb import *
import quimb.tensor as qtn

#%%
def phi_i(l,theta):
    return np.exp(1j*l*theta)
#%%
L=2
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
alpha = 1.5

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
#%%
"""
Need to split the A_L and the U matricies into MPS and MPOs
"""

"""
get next power of two
"""
p2 = 8
"""
pad a_L with zeros
"""
A_L_pad = np.lib.pad(A_L,(p2-(2*L+1),0))
U_pad = np.lib.pad(U,(p2-(2*L+1),0))
"""
create a padded identity matrix
"""
I_=np.identity(p2-(2*L+1))
I_pad = np.lib.pad(I_,(0,(2*L+1)))
"""
merge!
"""
A_L_TOT = A_L_pad
U_TOT = U_pad

"""
SVD till i cry
"""
# def return_decomp(A_L_TOT):
#     size=np.prod(np.shape(A_L_TOT))
#     length = np.log(size)/np.log(2)
#     arrs= []
    
#     #reshape to 2^2 x 2^(p2-2)
#     mat = A_L
#     for i in range(length-1,0,-1):
#         np.reshape(mat,(2,2**i))
#         u,s,v=svd(mat)
#         arrs.append(u)
#         mat=s@v
#     return arrs
size=np.prod(np.shape(A_L_TOT))
length_1 = np.log(size)/np.log(2)
dims = [2]*int(length_1)
psi_A = A_L_TOT.reshape(1,p2**2)
mps_A = qtn.MatrixProductState.from_dense(psi_A, dims)

size=np.prod(np.shape(U_TOT))
length_1 = np.log(size)/np.log(2)
dims = [2]*int(length_1)
psi_U = U_TOT.reshape(1,p2**2)
mps_U = qtn.MatrixProductState.from_dense(psi_U, dims)

H = qtn.MPO_identity(6,tags=["IDENT"])
H = qtn.MPO_ham_heis(6,tags=["IDENT"])

mps_A.align_(H)
tot = (H&mps_A)
tot.draw(color='IDENT')
fin = tot^...
new_g = fin.fuse({'b0':['b0','b1','b2'],'b1':['b3','b4','b5']}).data

new_g = new_g[3:,3:]
for i in range(length):
    xi_0 = np.zeros(2*L+1, dtype = 'complex')    
    for j in range(-L, L):
        for an in thp:
            xi_0[j+L] += (p_theta_thetap(theta[i], an) * phi_i(j, an))
    
    xi_norm = np.linalg.norm(xi_0)
    xi = xi_0/xi_norm
    point_w[i] = np.matmul(np.conjugate(xi), np.matmul(new_g, xi))/(2*np.pi)
    
plt.plot(np.mod(thp2,2*np.pi), point_w.real,'o')
plt.plot(np.mod(thp2,2*np.pi),vonmises(kappa).pdf(thp))
#%%
"""
get coeff
"""
g = [bessel(np.abs(ij), kappa) / bessel(0, kappa) for ij in range(5)]
"""
pad
"""
psi = np.lib.pad(g,(3,0))
mps_A = qtn.MatrixProductState.from_dense(psi, dims)