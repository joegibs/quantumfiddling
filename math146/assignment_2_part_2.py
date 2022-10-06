'''
2-d maps
'''

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
from scipy.stats import vonmises as vonmises
from scipy.special import iv as bessel

from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation

from opt_einsum import contract

#%%

def lplus1_space(i):
    return int(i - L)
def iter_2l(X):
    return  [i - L for i in X]

def calculate_g(U, f, n=1):
    g = U @ f
    for i in range(1, n):
        g = U @ g
    return g
#%%
def init_space(n):
    x=np.linspace(0,np.pi*2,n,endpoint=False)
    y=np.linspace(0,np.pi*2,n,endpoint=False)
    
    theta = np.meshgrid(x,y)
    return theta

def composition(theta,T,f,**kwargs):
    
    thetan=T(theta,**kwargs)
    zn=f(thetan,**kwargs)
    
    return thetan,zn

def T_rotation(theta,**kwargs):
    thetan=np.mod(np.array([theta[0]+kwargs['a1'],theta[1]+kwargs['a2']],dtype='object'),2*np.pi)
    return thetan

def T_skew_rotation(theta,**kwargs):
    thetan=np.array([np.mod(theta[0]+kwargs['a1'],np.pi*2),np.mod(theta[0]+theta[1],np.pi*2)])
    return thetan

def T_cat_map(theta,**kwargs):
    n=np.shape(theta[0])[0]
    thetan=np.empty_like(np.array(theta))
    
    for i in range(np.shape(theta[0])[0]):
        for j in range(np.shape(theta[0])[1]):
            v = np.array([[theta[0][i][j]],[theta[1][i][j]]])
            t = np.mod(np.matmul(kwargs['A'],v),1)
            thetan[0][i][j]=t[0][0]
            thetan[1][i][j]=t[1][0]
    return thetan

def f_add(theta):
    return np.add(theta[0],theta[1])

def f_bump(theta):
    xx=np.array(theta[0],dtype='float')
    yy=np.array(theta[1],dtype='float')
    return np.exp(1*(np.cos(xx)+np.cos(yy)))
def f_vonmis(theta, **kwargs):
    xx= np.array(theta[0],dtype='float')
    yy= np.array(theta[1],dtype='float')
    if "kappa" in kwargs.keys():
        return np.array(vonmises(kwargs["kappa"]).pdf(np.mod(xx, np.pi * 2))*vonmises(kwargs["kappa"]).pdf(np.mod(yy, np.pi * 2)))
def f_blank(theta):
    return np.zeros_like(theta[0])

#%%
def phi_i(l, theta):
    theta1=np.array(theta[0],dtype='float')
    theta2=np.array(theta[1],dtype='float')
    return np.exp(1j * (l[0] * theta1+l[1] * theta2))

def U_ij(index, T, X=[0, 2 * np.pi], **kwargs):
    f = lambda x: np.conjugate(phi_i([index[0],index[1]], theta)) * phi_i([index[2],index[3]], T(theta, **kwargs))
    n=100
    x=np.linspace(X[0], X[1],n,endpoint=False)
    y=np.linspace(X[0], X[1],n,endpoint=False)
    theta=np.meshgrid(x,y)
    a = f(theta)
    return np.trapz(np.trapz(a,x,axis=1),y)


#%%
def gen_u_and_f(L,T,**kwargs):
    U = np.zeros((2 * L + 1,2 * L + 1,2 * L + 1,2 * L + 1),dtype='complex')
    for ij in np.ndindex(np.shape(U)):
        # print(ij)
        ans = U_ij(iter_2l(ij), T, **kwargs)
        U[ij] = (ans) / ((np.pi * 2)**2)

    f = np.zeros((2 * L + 1,2 * L + 1))
    for ij in np.ndindex(np.shape(f)):
        ans1 = bessel(abs(lplus1_space(ij[0])), kwargs['kappa']) / bessel(0, kwargs['kappa'])
        ans2 =  bessel(abs(lplus1_space(ij[1])), kwargs['kappa']) / bessel(0, kwargs['kappa'])
        f[ij] = ans1*ans2
    
        
    return U,f

def basis_expansion(g, basis_func, theta):
    expan = np.zeros_like(theta[0], dtype="complex")
    for ij in np.ndindex(np.shape(g)):
        expan += g[ij] * basis_func([lplus1_space(ij[0]),lplus1_space(ij[1])], theta)
    return expan

#%%
n=100
a1 = np.pi*0
a2 = np.pi*0
L  = 1
kappa=1
T = T_rotation
theta=init_space(n)

kwargs = {'a1':a1,'a2':a2,'A':np.array([[2,1],[1,1]]),'kappa':kappa}

U,f = gen_u_and_f(L,T_rotation,**kwargs)
g=contract("ijkm,im",U,f)
zn=basis_expansion(g,phi_i,theta).real

h = plt.contourf(theta[0],theta[1], zn)
plt.axis('scaled')
plt.colorbar()
plt.show()
#%%
fig, axs = plt.subplots(1,2)

n=100
a1 =np.pi
a2 = 0*np.pi/2
L  = 3
kappa=5
T = T_skew_rotation
theta=init_space(n)

kwargs = {'a1':a1,'a2':a2,'A':np.array([[2,1],[1,1]]),'kappa':kappa}

_,zn=composition(theta,T,f_vonmis,**kwargs)
                                              
axs[0].contourf(theta[0],theta[1], zn)
axs[0].set_aspect("equal")
axs[0].set_title('Continuous')
# plt.colorbar()

U,f = gen_u_and_f(L,T,**kwargs)
g=contract("ijkm,jk->im",U,f)
zn2=basis_expansion(g,phi_i,theta).real

axs[1].contourf(theta[0],theta[1], zn2)
axs[1].set_aspect("equal")
axs[1].set_title('Discrete')

fig.suptitle(f'Kappa: {kappa}, L:{L}')
fig.show()

#%%
"""
animate for a couple different values
"""
plt.close("all")  # Clear anything left over from prior runs.

n=100
a1 = 1/10
a2 = 1/10
L  = 3
kappa=1
T = T_rotation
theta=init_space(n)

kwargs = {'a1':a1,'a2':a2,'A':np.array([[2,1],[1,1]]),'kappa':kappa}
Nt = 100

class data_z:
    def __init__(self, angle,Ud,fd, T,**kwargs):
        self.angle = angle
        self.anglen = angle
        self.T=T
        self.Ud = Ud
        self.fd = fd
        self.g =contract("ijkm,jk->im",self.Ud,self.fd)

        self.kwargs = kwargs
        self.znc = composition(self.angle,T,f_vonmis,**kwargs)
        self.znd = basis_expansion(self.g, phi_i, self.angle).real

    def step(self):
        self.g=contract("ijkm,jk->im",self.Ud,np.reshape(self.g,np.shape(self.fd)))
        self.znd = basis_expansion(self.g, phi_i, self.angle).real
        _, self.znc = composition(self.anglen, self.T, f_vonmis, **self.kwargs)
        self.anglen = self.T(self.anglen, **self.kwargs)

def animate(ii, lis, data):
    global cont
    data.step()
    zn1 = data.znc
    thetan= data.anglen
    zn2 = data.znd
    # print(np.shape(zn))
    for i in cont:
        for c in i.collections:
            c.remove()  # removes only the contours, leaves the rest intact
    cont[1].contourf(theta[0], theta[1], zn2)
    cont[0].contourf(theta[0], theta[1], zn1)
    fig.suptitle(f't = {ii}')

    return cont

fig, cont = plt.subplots(1,2, sharey=True)


_,zn=composition(theta,T,f_vonmis,**kwargs)
                                              
cont[0].contourf(theta[0],theta[1], zn)

U,f = gen_u_and_f(L,T,**kwargs)
g=contract("ijkm,jk->im",U,f)
zn2=basis_expansion(g,phi_i,theta).real

cont[1].contourf(theta[0], theta[1], zn2)    # first image on screen
cont[0].set_aspect("equal")
cont[1].set_aspect("equal")

mydata = data_z(theta,U,f, T,**kwargs)

mylis = [1, 2, 3]
anim = animation.FuncAnimation(
    fig, animate, fargs=(mylis, mydata), frames=Nt, interval=100
)

#%%
n=100
theta=init_space(n)
a1 = np.pi
a2 = np.pi
L  = 2
kappa=2
T = T_rotation
kwargs = {'a1':a1,'a2':a2,'A':np.array([[2,1],[1,1]]),'kappa':kappa}

U,f = gen_u_and_f(L,T,**kwargs)
g=contract("ijkm,jk->im",U,f)
zn=basis_expansion(g,phi_i,theta).real
# g=contract("ijkm,jk->im",U,np.reshape(g,np.shape(f)))
# g=contract("ijkm,jk->im",U,np.reshape(g,np.shape(f)))
# g=contract("ijkm,jk->im",U,np.reshape(g,np.shape(f)))
# g=contract("ijkm,jk->im",U,np.reshape(g,np.shape(f)))
# g=contract("ijkm,jk->im",U,np.reshape(g,np.shape(f)))

# zn=basis_expansion(g,phi_i,theta).real

# axs[1].contourf(theta[0],theta[1], zn)
h = plt.contourf(theta[0],theta[1], zn)
plt.axis('scaled')
plt.colorbar()
plt.show()