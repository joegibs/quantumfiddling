
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
from scipy.stats import vonmises as vonmises
import scipy.integrate as integrate
from scipy.special import jv as bessel

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
#%%

angle = np.linspace( 0 , 2 * np.pi , 150 ) 
 
radius = 1
x = radius * np.cos( angle ) 
y = radius * np.sin( angle ) 
_,colors = composition(angle,T_rotation,f_vonmis,**{'alpha':0.5,'kappa':5})
points = np.array([x, y]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)

fig, axs = plt.subplots(1, sharex=True, sharey=True)

norm = plt.Normalize(colors.min(), colors.max())
lc = LineCollection(segments, cmap='viridis', norm=norm)
lc.set_array(colors)
lc.set_linewidth(8)
line = axs.add_collection(lc)
fig.colorbar(line, ax=axs)

axs.set_xlim(-1.3, 1.3)
axs.set_ylim(-1.3, 1.3)

axs.set_aspect('equal')
# plt.show()
#%%
def plot(T,f,angle=None,**kwargs):
    angle = np.linspace( 0 , 2 * np.pi , 150 ) 
     
    radius = 1
    x = radius * np.cos( angle ) 
    y = radius * np.sin( angle ) 
    _,colors = composition(angle,T,f,**kwargs)
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    fig, axs = plt.subplots(1, sharex=True, sharey=True)

    norm = plt.Normalize(colors.min(), colors.max())
    lc = LineCollection(segments, cmap='viridis', norm=norm)
    lc.set_array(colors)
    lc.set_linewidth(8)
    line = axs.add_collection(lc)
    fig.colorbar(line, ax=axs)

    axs.set_xlim(-1.3, 1.3)
    axs.set_ylim(-1.3, 1.3)

    axs.set_aspect('equal')
    plt.show()
#%%
def composition(theta,T,f,**kwargs):
    thetan=T(theta,**kwargs)
    zn=f(thetan,**kwargs)
    return thetan,zn
def T_blank(theta,**kwargs):
    return theta
def T_rotation(theta,**kwargs):
    if 'alpha' in kwargs.keys():
        return np.mod(theta+kwargs['alpha'],np.pi*2)
    else:
        return theta+np.pi
def T_double_map(theta,**kwargs):
    return np.mod(theta*2,2*np.pi)
def f_ang(theta,**kwargs):
    return theta
def f_vonmis(theta,**kwargs):
    if 'kappa' in kwargs.keys():
        return vonmises(kwargs['kappa']).pdf(theta)
    else:
        return vonmises(1).pdf(theta)
#%%
plot(T_rotation,f_vonmis,**{'alpha':-np.pi/2,'kappa':5})
#%%
plt.close('all') # Clear anything left over from prior runs.

Nt=180
bound=1
alpha=1/20
angle = np.linspace( 0 , 2 * np.pi , 150 ) 


class data_z:
    def __init__(self,angle,T,f,**kwargs):
        self.angle=angle
        self.anglen=angle
        self.T = T
        self.f = f
        self.kwargs=kwargs
        self.zn=f(angle)
    def step(self):
        self.anglen,self.zn=composition(self.anglen,self.T,self.f,**self.kwargs)

# animation function
def animate(i,lis,data):
    global cont
    for c in cont.collections:
        c.remove()  # removes only the contours, leaves the rest intact
    data.step()
    radius = 1
    angle = data.angle
    x = radius * np.cos( angle ) 
    y = radius * np.sin( angle ) 
    
    _,colors = composition(data.anglen,data.T,data.f)
    points = np.transpose(np.array([x, y])).reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    norm = plt.Normalize(colors.min(), colors.max())
    lc = LineCollection(segments, cmap='viridis', norm=norm)
    lc.set_array(colors)
    lc.set_linewidth(8)
    line = cont.add_collection(lc)
    

    plt.title('t = %i' % (i))
    return cont

mydata=data_z(angle,T_rotation,f_vonmis,**{'alpha':alpha})


angle = np.linspace( 0 , 2 * np.pi , 150 ) 
 
radius = 1
x = radius * np.cos( angle ) 
y = radius * np.sin( angle ) 
_,colors = composition(mydata.angle,T_rotation,f_vonmis,**{'alpha':0,'kappa':5})
points = np.array([x, y]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)

fig, cont = plt.subplots(1, sharex=True, sharey=True)

norm = plt.Normalize(colors.min(), colors.max())
lc = LineCollection(segments, cmap='viridis', norm=norm)
lc.set_array(colors)
lc.set_linewidth(8)
line = cont.add_collection(lc)
fig.colorbar(line, ax=cont)

cont.set_xlim(-1.3, 1.3)
cont.set_ylim(-1.3, 1.3)

cont.set_aspect('equal')


mylis=[1,2,3]
anim = animation.FuncAnimation(fig, animate, fargs=(mylis,mydata), frames=Nt,interval = 1)


#%%
def phi_i(l,theta):
    return np.exp(1j*l*theta)
    # return np.cos(l*theta)+1j*np.sin(l*theta)
def U_ij(i,j,T,X=[0,2*np.pi],**kwargs):
    return integrate.quad(lambda x: (np.conjugate(phi_i(i,x))*phi_i(j,T(x,**kwargs))),X[0],X[1])

#%%
'''
not sure what n is, 
not sure n the iteration,

'''
L=4
kappa=5
U=np.zeros((2*L+1,2*L+1))

def lplus1_space(i):
    return int(i-L)
'''
calculate U_ij
'''
for ij in np.ndindex(np.shape(U)):
    ans = U_ij(lplus1_space(ij[0]),lplus1_space(ij[1]),T_double_map,**{'alpha':1})
    U[ij] = np.round((ans[0]+1j*ans[1])/(2*np.pi),6)
    
f=np.zeros((L,1))
'''
Calculate fl
'''
for ij in np.ndindex(np.shape(f)):
    ans = bessel(lplus1_space(ij[0]),kappa)/bessel(0,kappa)
    f[ij] = ans

def calculate_g(U,f,n=1):
    g=U@f
    for i in range(1,n):
        g=U@g
    return g
g=calculate_g(U,f)
def basis_expansion(g,basis_func,theta):
    expan=np.zeros_like(theta,dtype='complex')
    for ij in np.ndindex(np.shape(g)):
        # print(ij)
        expan+=g[ij]*basis_func(lplus1_space(ij[0]),theta)
    return expan
    
theta=np.linspace(0 , 2 * np.pi , 150 ) 
radius = 1
x = radius * np.cos( angle ) 
y = radius * np.sin( angle ) 
colors = basis_expansion(g,phi_i,theta).real
points = np.array([x, y]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)

fig, axs = plt.subplots(1, sharex=True, sharey=True)

norm = plt.Normalize(colors.min(), colors.max())
lc = LineCollection(segments, cmap='viridis', norm=norm)
lc.set_array(colors)
lc.set_linewidth(8)
line = axs.add_collection(lc)
fig.colorbar(line, ax=axs)

axs.set_xlim(-1.3, 1.3)
axs.set_ylim(-1.3, 1.3)

axs.set_aspect('equal')
#%%

'''
should just need to calculate g over and over?
'''
plt.close('all') # Clear anything left over from prior runs.

Nt=180
bound=1
alpha=1/20
angle = np.linspace( 0 , 2 * np.pi , 150 ) 


class data_z:
    def __init__(self,angle,T,f,**kwargs):
        self.angle=angle
        self.anglen=angle
        self.T = T
        self.f = f
        self.kwargs=kwargs
        self.zn=f(angle)
    def step(self):
        self.anglen,self.zn=composition(self.anglen,self.T,self.f,**self.kwargs)

# animation function
def animate(i,lis,data):
    global cont
    for c in cont.collections:
        c.remove()  # removes only the contours, leaves the rest intact
    data.step()
    radius = 1
    angle = data.angle
    x = radius * np.cos( angle ) 
    y = radius * np.sin( angle ) 
    
    _,colors = composition(data.anglen,data.T,data.f)
    points = np.transpose(np.array([x, y])).reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    norm = plt.Normalize(colors.min(), colors.max())
    lc = LineCollection(segments, cmap='viridis', norm=norm)
    lc.set_array(colors)
    lc.set_linewidth(8)
    line = cont.add_collection(lc)
    

    plt.title('t = %i' % (i))
    return cont

mydata=data_z(angle,T_rotation,f_vonmis,**{'alpha':alpha})


angle = np.linspace( 0 , 2 * np.pi , 150 ) 
 
radius = 1
x = radius * np.cos( angle ) 
y = radius * np.sin( angle ) 
_,colors = composition(mydata.angle,T_rotation,f_vonmis,**{'alpha':0,'kappa':5})
points = np.array([x, y]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)

fig, cont = plt.subplots(1, sharex=True, sharey=True)

norm = plt.Normalize(colors.min(), colors.max())
lc = LineCollection(segments, cmap='viridis', norm=norm)
lc.set_array(colors)
lc.set_linewidth(8)
line = cont.add_collection(lc)
fig.colorbar(line, ax=cont)

cont.set_xlim(-1.3, 1.3)
cont.set_ylim(-1.3, 1.3)

cont.set_aspect('equal')


mylis=[1,2,3]
anim = animation.FuncAnimation(fig, animate, fargs=(mylis,mydata), frames=Nt,interval = 1)
