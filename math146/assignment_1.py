import numpy as np
import matplotlib.pyplot as plt
#%%
def init_space(n):
    x=np.linspace(0,1,n,endpoint=False)
    y=np.linspace(0,1,n,endpoint=False)
    
    xx,yy = np.meshgrid(x,y)
    return xx,yy

def composition(xx,yy,T,f,**kwargs):
    
    xn,yn=T(xx,yy,**kwargs)
    zn=f(xn,yn)
    
    return xn,yn,zn


def T_rotation(xx,yy,**kwargs):
    xn=(xx+kwargs['a1'])%1
    yn=(yy+kwargs['a2'])%1
    return xn,yn

def T_skew_rotation(xx,yy,**kwargs):
    xn=(xx+kwargs['a1'])%1
    yn=(yy+xx)%1
    return xn,yn

def T_cat_map(xx,yy,**kwargs):
    n=np.shape(xx)[0]
    xn=np.empty((n,n))
    yn=np.empty((n,n))
    
    for i in range(np.shape(xx)[0]):
        for j in range(np.shape(xx)[1]):
            v = np.array([[xx[i][j]],[yy[i][j]]])
            t = np.mod(np.matmul(kwargs['A'],v),1)
            xn[i][j]=t[0][0]
            yn[i][j]=t[1][0]
    return xn,yn

def f_add(xn,yn):
    return xn+yn

def f_bump(xx,yy):
    return np.exp(1*(np.cos(2*np.pi*xx)+np.cos(2*np.pi*yy)))

#%%
#%%
n=1000
xx,yy=init_space(n)
xn,yn,zn=composition(xx,yy,T_rotation,f_bump,**{'a1':1/2,'a2':1/2})
                                              
h = plt.contourf(xx,yy, zn)
plt.axis('scaled')
plt.colorbar()
plt.show()
#%%
n=1000
xx,yy=init_space(n)
xn,yn,zn=composition(xx,yy,T_cat_map,f_bump,**{'A':np.array([[2,1],[1,1]])})
          
                                    
h = plt.contourf(xx,yy, zn)
plt.axis('scaled')
plt.colorbar()
plt.show()

#%%
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
#%%
plt.close('all') # Clear anything left over from prior runs.

Nt=40
bound=1

xx,yy=init_space(100)

fig = plt.figure() # Must have figure object for movie.
ax = plt.axes(xlim=(0, bound), ylim=(0, bound))

class data_z:
    def __init__(self,xx,yy):
        self.xx=xx
        self.yy=yy
        self.xn=xx
        self.yn=yy
        self.zn=f_bump(xx,yy)
    def step(self):
        self.xn,self.yn,self.zn=composition(self.xn,self.yn,T_cat_map,f_bump,**{'A':np.array([[2,1],[1,1]])})

# animation function
def animate(i,lis,data):
    global cont
    for c in cont.collections:
        c.remove()  # removes only the contours, leaves the rest intact
        
    cont = plt.contourf(data.xx, data.yy, data.zn)
    data.step()
    plt.title('t = %i' % (i))
    return cont

mydata=data_z(xx,yy)
cont = plt.contourf(mydata.xx, mydata.yy, mydata.zn)
mylis=[1,2,3]
anim = animation.FuncAnimation(fig, animate, fargs=(mylis,mydata), frames=Nt,interval = 100)