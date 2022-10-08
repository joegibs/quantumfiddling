import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
from scipy.stats import vonmises as vonmises
from scipy.special import iv as bessel

from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation

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
'''
Start with continuous form
'''
a = plot(T_rotation, f_vonmis, **{"alpha": np.pi / 2, "kappa": 1})

#%%
''' animation '''
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
'''
plot
'''
alpha = 1.234
kappa = 3
L=4
T=T_rotation

U,f = gen_u_and_f(alpha,kappa,L,T)
g = calculate_g(U,f)
theta = np.linspace(0, 2 * np.pi, 300)
radius = 1
x = radius * np.cos(theta)
y = radius * np.sin(theta)
colors = basis_expansion(g, phi_i, theta).real
points = np.array([x, y]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)

fig, axs = plt.subplots(1)

norm = plt.Normalize(colors.min(), colors.max())
lc = LineCollection(segments, cmap="gist_ncar", norm=norm)
lc.set_array(colors)
lc.set_linewidth(13)
line = axs.add_collection(lc)
fig.colorbar(line, ax=axs)
axs.set_xlim(-1.3, 1.3)
axs.set_ylim(-1.3, 1.3)
axs.set_aspect("equal")
#%%
"""
animate for a couple different values
"""
plt.close("all")  # Clear anything left over from prior runs.

Nt = 180
bound = 1
alpha = np.pi / 50
T=T_rotation

# angle = np.linspace(0, 2 * np.pi, 150)

class data_z:
    def __init__(self, angle, U, f, **kwargs):
        self.angle = angle
        self.anglen = angle
        self.U = U
        self.f = f
        self.g = self.U@self.f
        self.kwargs = kwargs
        self.zn = basis_expansion(self.g, phi_i, self.angle).real

    def step(self):
        self.g=self.U@self.g
        self.zn = basis_expansion(self.g, phi_i, self.angle).real

def animate(ii, lis, data):
    global cont
    for i in range(len(data)):
        for c in cont[i].collections:
            c.remove()
        
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

theta = np.linspace(0, 2 * np.pi, 300)
alphav=[np.pi/20,np.pi/30,np.pi/100]
kappav=[1,3,9]
uf = [gen_u_and_f(alphav[i],kappav[i],8,T) for i in range(3)]
mydata = [data_z(theta,uf[i][0] , uf[i][1]) for i in range(3)]

radius = 1
x = radius * np.cos(theta)
y = radius * np.sin(theta)
colors = basis_expansion(f, phi_i, theta).real
points = np.array([x, y]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)

fig, cont = plt.subplots(1,3, sharey=True)

norm = plt.Normalize(colors.min(), colors.max())
lc = LineCollection(segments, cmap="flag", norm=norm)
lc.set_array(colors)
lc.set_linewidth(12)
line = cont[0].add_collection(lc)

for i in range(np.size(cont)):
    cont[i].set_xlim(-1.3, 1.3)
    cont[i].set_ylim(-1.3, 1.3)
    
    cont[i].set_aspect("equal")


mylis = [1, 2, 3]
anim = animation.FuncAnimation(
    fig, animate, fargs=(mylis, mydata), frames=Nt, interval=1
)
#%%
'''
less fun plots but easier to see the difference
'''
alpha = 1.2345
kappa = 100
L=8
angle = np.linspace(0, 2 * np.pi, 300)

_, colors = composition(angle, T_rotation, f_vonmis, **{"alpha": alpha, "kappa": kappa})

U,f = gen_u_and_f(alpha,kappa,L,T_rotation)
g = calculate_g(U,f)
colorsd = basis_expansion(g, phi_i, angle).real/(np.pi*2)

fig, axs = plt.subplots(2)
fig.suptitle(f'Discritized vs Continuous transformation, alpha = {alpha}, kappa = {kappa}')

axs[0].plot(angle, colors)
axs[0].set_title('Continuous')
axs[0].set_xlabel('')
axs[1].plot(angle, colorsd)
axs[1].set_title('Discrete')

plt.show()
#%%

'''
Same with doubling map
'''
alpha = 1.2345
kappa = 100
L=8
angle = np.linspace(0, 2 * np.pi, 300)

_, colors = composition(angle, T_double_map, f_vonmis, **{"alpha": alpha, "kappa": kappa})

U,f = gen_u_and_f(alpha,kappa,L,T_double_map)
g = calculate_g(U,f)
colorsd = basis_expansion(g, phi_i, angle).real/(np.pi*2)

fig, axs = plt.subplots(2)
fig.suptitle(f'Discritized vs Continuous transformation, alpha = {alpha}, kappa = {kappa}')

axs[0].plot(angle, colors)
axs[0].set_title('Continuous')
axs[0].set_xlabel('')
axs[1].plot(angle, colorsd)
axs[1].set_title('Discrete')

plt.show()