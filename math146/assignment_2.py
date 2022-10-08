import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
from scipy.stats import vonmises as vonmises
import scipy.integrate as integrate
from scipy.special import iv as bessel

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation

#%%

# angle = np.linspace( 0 , 2 * np.pi , 150 )

# radius = 1
# x = radius * np.cos( angle )
# y = radius * np.sin( angle )
# _,colors = composition(angle,T_rotation,f_vonmis,**{'alpha':0.5*np.pi,'kappa':5})
# points = np.array([x, y]).T.reshape(-1, 1, 2)
# segments = np.concatenate([points[:-1], points[1:]], axis=1)

# fig, axs = plt.subplots(1, sharex=True, sharey=True)

# norm = plt.Normalize(colors.min(), colors.max())
# lc = LineCollection(segments, cmap='viridis', norm=norm)
# lc.set_array(colors)
# lc.set_linewidth(8)
# line = axs.add_collection(lc)
# fig.colorbar(line, ax=axs)

# axs.set_xlim(-1.3, 1.3)
# axs.set_ylim(-1.3, 1.3)

# axs.set_aspect('equal')
# # plt.show()
#%%
def plot(T, f, angle=None, **kwargs):
    angle = np.linspace(0, 2 * np.pi, 300)

    radius = 1

    yarg, colors = composition(angle, T, f, **kwargs)
    x = radius * np.cos(angle)
    y = radius * np.sin(angle)
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    fig, axs = plt.subplots(1, sharex=True, sharey=True)

    norm = plt.Normalize(colors.min(), colors.max())
    lc = LineCollection(segments, cmap="viridis", norm=norm)
    lc.set_array(colors)
    lc.set_linewidth(8)
    line = axs.add_collection(lc)
    fig.colorbar(line, ax=axs)

    axs.set_xlim(-1.3, 1.3)
    axs.set_ylim(-1.3, 1.3)

    axs.set_aspect("equal")
    plt.show()


#%%
def composition(theta, T, f, **kwargs):
    thetan = T(theta, **kwargs)
    zn = f(thetan, **kwargs)
    return thetan, zn


def T_blank(theta, **kwargs):
    return np.mod(theta, 2 * np.pi)


def T_rotation(theta, **kwargs):
    if "alpha" in kwargs.keys():
        return np.mod((theta - kwargs["alpha"]), 2 * np.pi)
    else:
        return np.mod(theta, 2 * np.pi)


def T_double_map(theta, **kwargs):
    return np.mod(theta * 2, 2 * np.pi)


def f_ang(theta, **kwargs):
    return theta


def f_vonmis(theta, **kwargs):
    if "kappa" in kwargs.keys():
        return vonmises(kwargs["kappa"]).pdf(np.mod(theta + np.pi, np.pi * 2))
    else:
        return vonmises(1).pdf(theta)


#%%
a = plot(T_rotation, f_vonmis, **{"alpha": np.pi / 2, "kappa": 1})
#%%
plt.close("all")  # Clear anything left over from prior runs.

Nt = 180
bound = 1
alpha = np.pi / 20
angle = np.linspace(0, 2 * np.pi, 150)


class data_z:
    def __init__(self, angle, T, f, **kwargs):
        self.angle = angle
        self.anglen = angle
        self.T = T
        self.f = f
        self.kwargs = kwargs
        self.z = f(angle)
        self.zn = f(angle)

    def step(self):
        _, self.zn = composition(self.anglen, self.T, self.f, **self.kwargs)
        self.anglen = self.T(self.anglen, **self.kwargs)


# animation function
def animate(i, lis, data):
    global cont
    for c in cont.collections:
        c.remove()  # removes only the contours, leaves the rest intact
    colors = data.z
    data.step()
    radius = 1
    angle = data.anglen
    x = radius * np.cos(angle)
    y = radius * np.sin(angle)

    points = np.transpose(np.array([x, y])).reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    norm = plt.Normalize(colors.min(), colors.max())
    lc = LineCollection(segments, cmap="viridis", norm=norm)
    lc.set_array(colors)
    lc.set_linewidth(8)
    line = cont.add_collection(lc)

    plt.title("t = %i" % (i))
    return cont


mydata = data_z(angle, T_rotation, f_vonmis, **{"alpha": alpha})


angle = np.linspace(0, 2 * np.pi, 150)

radius = 1
x = radius * np.cos(angle)
y = radius * np.sin(angle)
_, colors = composition(mydata.angle, T_rotation, f_vonmis, **{"alpha": 0, "kappa": 5})
points = np.array([x, y]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)

fig, cont = plt.subplots(1, sharex=True, sharey=True)

norm = plt.Normalize(colors.min(), colors.max())
lc = LineCollection(segments, cmap="viridis", norm=norm)
lc.set_array(colors)
lc.set_linewidth(8)
line = cont.add_collection(lc)
fig.colorbar(line, ax=cont)

cont.set_xlim(-1.3, 1.3)
cont.set_ylim(-1.3, 1.3)

cont.set_aspect("equal")


mylis = [1, 2, 3]
anim = animation.FuncAnimation(
    fig, animate, fargs=(mylis, mydata), frames=Nt, interval=1
)


#%%
def phi_i(l, theta):
    return np.exp(1j * l * theta)
    # return np.cos(l*theta)+1j*np.sin(l*theta)


def U_ij(i, j, T, X=[0, 2 * np.pi], **kwargs):
    # print(i,j)
    f = lambda x: (np.conjugate(phi_i(i, x)) * phi_i(j, T(x, **kwargs)))
    a = f(np.linspace(X[0], X[1], 300))
    return np.trapz(a, np.linspace(X[0], X[1], 300))


#%%
"""

"""
L = 8
kappa = 50
U = np.zeros((2 * L + 1, 2 * L + 1),dtype='complex')
alpha = np.pi/2*0



def lplus1_space(i):
    return int(i - L)


"""
calculate U_ij
"""
for ij in np.ndindex(np.shape(U)):
    # print(lplus1_space(ij[0]),lplus1_space(ij[1]))
    ans = U_ij(lplus1_space(ij[0]), lplus1_space(ij[1]), T_rotation, **{"alpha": alpha})
    U[ij] = (ans) / (np.pi * 2)

f = np.zeros((2 * L + 1, 1))
"""
Calculate fl
"""
for ij in np.ndindex(np.shape(f)):
    # print(abs(lplus1_space(ij[0])))
    ans = bessel(abs(lplus1_space(ij[0])), kappa) / bessel(0, kappa)
    f[ij] = ans


def calculate_g(U, f, n=1):
    g = U @ f
    for i in range(1, n):
        g = U @ g
    return g


g = calculate_g(U, f,5)


def basis_expansion(g, basis_func, theta):
    expan = np.zeros_like(theta, dtype="complex")
    for ij in np.ndindex(np.shape(g)):
        # print(ij)
        expan += g[ij] * basis_func(lplus1_space(ij[0]), theta)
    return expan


theta = np.linspace(0, 2 * np.pi, 300)
radius = 1
x = radius * np.cos(theta)
y = radius * np.sin(theta)
colors = basis_expansion(g, phi_i, theta).real
points = np.array([x, y]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)

fig, axs = plt.subplots(1, sharex=True, sharey=True)

norm = plt.Normalize(colors.min(), colors.max())
lc = LineCollection(segments, cmap="viridis", norm=norm)
lc.set_array(colors)
lc.set_linewidth(8)
line = axs.add_collection(lc)
fig.colorbar(line, ax=axs)

axs.set_xlim(-1.3, 1.3)
axs.set_ylim(-1.3, 1.3)

axs.set_aspect("equal")
# plt.plot(theta/np.pi,colors)
# plt.title(f'{alpha/np.pi}')
#%%
def basis_expansion(g, basis_func, theta):
    expan = np.zeros_like(theta, dtype="complex")
    for ij in np.ndindex(np.shape(g)):
        # print(ij)
        expan += g[ij] * basis_func(lplus1_space(ij[0]), theta)
    return expan
#%%
# '''
# should just need to calculate g over and over?
# '''
plt.close("all")  # Clear anything left over from prior runs.

Nt = 180
bound = 1
alpha = np.pi / 50
# angle = np.linspace(0, 2 * np.pi, 150)


class data_z:
    def __init__(self, angle, U, f, **kwargs):
        self.angle = angle
        self.anglen = angle
        self.U = U
        self.f = f
        self.g = self.U@self.f
        self.kwargs = kwargs
        # self.z = f(angle)
        self.zn = basis_expansion(self.g, phi_i, self.angle).real

    def step(self):
        self.g=self.U@self.g
        self.zn = basis_expansion(self.g, phi_i, self.angle).real
        # self.anglen = self.T(self.anglen, **self.kwargs)


# animation function
def animate(i, lis, data):
    global cont
    for c in cont.collections:
        c.remove()  # removes only the contours, leaves the rest intact
    
    data.step()

    colors = data.zn
    radius = 1
    angle = np.linspace(0, 2 * np.pi, 300)

    x = radius * np.cos(angle)
    y = radius * np.sin(angle)

    points = np.transpose(np.array([x, y])).reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    norm = plt.Normalize(colors.min(), colors.max())
    lc = LineCollection(segments, cmap="viridis", norm=norm)
    lc.set_array(colors)
    lc.set_linewidth(8)
    line = cont.add_collection(lc)

    plt.title("t = %i" % (i))
    return cont

L = 8
kappa = 2
U = np.zeros((2 * L + 1, 2 * L + 1),dtype='complex')


def lplus1_space(i):
    return int(i - L)
for ij in np.ndindex(np.shape(U)):
    # print(lplus1_space(ij[0]),lplus1_space(ij[1]))
    ans = U_ij(lplus1_space(ij[0]), lplus1_space(ij[1]), T_rotation, **{"alpha": alpha})
    U[ij] = (ans) / (np.pi * 2)

f = np.zeros((2 * L + 1, 1))
for ij in np.ndindex(np.shape(f)):
    # print(abs(lplus1_space(ij[0])))
    ans = bessel(abs(lplus1_space(ij[0])), kappa) / bessel(0, kappa)
    f[ij] = ans
def calculate_g(U, f, n=1):
    g = U @ f
    for i in range(1, n):
        g = U @ g
    return g



theta = np.linspace(0, 2 * np.pi, 300)

mydata = data_z(theta, U, f)

radius = 1
x = radius * np.cos(theta)
y = radius * np.sin(theta)
colors = basis_expansion(f, phi_i, theta).real
points = np.array([x, y]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)

fig, cont = plt.subplots(1, sharex=True, sharey=True)

norm = plt.Normalize(colors.min(), colors.max())
lc = LineCollection(segments, cmap="viridis", norm=norm)
lc.set_array(colors)
lc.set_linewidth(8)
line = cont.add_collection(lc)
fig.colorbar(line, ax=cont)

cont.set_xlim(-1.3, 1.3)
cont.set_ylim(-1.3, 1.3)

cont.set_aspect("equal")


mylis = [1, 2, 3]
anim = animation.FuncAnimation(
    fig, animate, fargs=(mylis, mydata), frames=Nt, interval=1
)

#%%
def lplus1_space(i):
    return int(i - L)
def gen_u_and_f(alpha,kappa,L):
    U = np.zeros((2 * L + 1, 2 * L + 1),dtype='complex')

    """
    calculate U_ij
    """
    for ij in np.ndindex(np.shape(U)):
        # print(lplus1_space(ij[0]),lplus1_space(ij[1]))
        ans = U_ij(lplus1_space(ij[0]), lplus1_space(ij[1]), T_rotation, **{"alpha": alpha})
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
# '''
# should just need to calculate g over and over?
# '''
plt.close("all")  # Clear anything left over from prior runs.

Nt = 180
bound = 1
alpha = np.pi / 50
# angle = np.linspace(0, 2 * np.pi, 150)


class data_z:
    def __init__(self, angle, U, f, **kwargs):
        self.angle = angle
        self.anglen = angle
        self.U = U
        self.f = f
        self.g = self.U@self.f
        self.kwargs = kwargs
        # self.z = f(angle)
        self.zn = basis_expansion(self.g, phi_i, self.angle).real

    def step(self):
        self.g=self.U@self.g
        self.zn = basis_expansion(self.g, phi_i, self.angle).real
        # self.anglen = self.T(self.anglen, **self.kwargs)


# animation function
def animate(ii, lis, data):
    global cont
    for i in range(len(data)):
        for c in cont[i].collections:
            c.remove()  # removes only the contours, leaves the rest intact
        
        data[i].step()
    
        colors = data[i].zn
        radius = 1
        angle = np.linspace(0, 2 * np.pi, 300)
    
        x = radius * np.cos(angle)
        y = radius * np.sin(angle)
    
        points = np.transpose(np.array([x, y])).reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
    
        norm = plt.Normalize(colors.min(), colors.max())
        lc = LineCollection(segments, cmap="viridis", norm=norm)
        lc.set_array(colors)
        lc.set_linewidth(12)
        line = cont[i].add_collection(lc)
    
    cont[0].set_title(f"\n kappa={kappav[0]}\n alpha = {np.round(alphav[0]/np.pi,3)}")
    cont[1].set_title(f"t = {ii}\n kappa={kappav[1]}\n alpha = {np.round(alphav[1]/np.pi,3)}")
    cont[2].set_title(f"\n kappa={kappav[2]}\n alpha = {np.round(alphav[2]/np.pi,3)}")
    return cont

L = 8
kappa = 2
U = np.zeros((2 * L + 1, 2 * L + 1),dtype='complex')


def lplus1_space(i):
    return int(i - L)
for ij in np.ndindex(np.shape(U)):
    # print(lplus1_space(ij[0]),lplus1_space(ij[1]))
    ans = U_ij(lplus1_space(ij[0]), lplus1_space(ij[1]), T_rotation, **{"alpha": alpha})
    U[ij] = (ans) / (np.pi * 2)

f = np.zeros((2 * L + 1, 1))
for ij in np.ndindex(np.shape(f)):
    # print(abs(lplus1_space(ij[0])))
    ans = bessel(abs(lplus1_space(ij[0])), kappa) / bessel(0, kappa)
    f[ij] = ans
def calculate_g(U, f, n=1):
    g = U @ f
    for i in range(1, n):
        g = U @ g
    return g



theta = np.linspace(0, 2 * np.pi, 300)
alphav=[np.pi/20,np.pi/30,np.pi/100]
kappav=[1,3,9]
uf = [gen_u_and_f(alphav[i],kappav[i],8) for i in range(3)]
mydata = [data_z(theta,uf[i][0] , uf[i][1]) for i in range(3)]

radius = 1
x = radius * np.cos(theta)
y = radius * np.sin(theta)
colors = basis_expansion(f, phi_i, theta).real
points = np.array([x, y]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)

fig, cont = plt.subplots(1,3, sharey=True)

norm = plt.Normalize(colors.min(), colors.max())
lc = LineCollection(segments, cmap="viridis", norm=norm)
lc.set_array(colors)
lc.set_linewidth(12)
line = cont[0].add_collection(lc)
# fig.colorbar(line, ax=cont)

for i in range(np.size(cont)):
    cont[i].set_xlim(-1.3, 1.3)
    cont[i].set_ylim(-1.3, 1.3)
    
    cont[i].set_aspect("equal")


mylis = [1, 2, 3]
anim = animation.FuncAnimation(
    fig, animate, fargs=(mylis, mydata), frames=Nt, interval=1
)

