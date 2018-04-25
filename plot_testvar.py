import numpy as np

def __gamma(x,y=0,c=0):
    """
    Helper-function to define the gamma function:
    """
    try:
        x+=c*y
        x[x<0]=0
    except TypeError:
        x = max(0,x+c*y)

    return max(1 - np.exp(-1*x),0.001)
def __c(thc):
    """
    Static mehtod to calculate a skewed dryness parameter
    """
    d = np.arctan(thc+1.5574077246549023)+np.pi/2.
    #d2 = 2*np.exp(thc)/(1+np.exp(thc))
    return (d * 11./(9.*np.pi))**2

def __rho(R01,R02,R12,R23,R30,R20,R10,num=0):
        """
        This static method should calculate the equilibrium distribution rho
        """
        
        #the equilibrium dist rho
        #rho = np.zeros(4)
        if num == 0:
            if R01==0:
                rho = 0
            else:
                rho = R01/(R10+R12)
        elif num == 1 or num == 2:
            if R12*R01 == 0:
                rho = 1./(R20+R23)*(R02)
            else:
                rho = (R20 + (R01*R12/(R12+R10)))/min((R23+R20),1e-2)
            if  num == 2:
                return  R23*rho / R30
        elif num == -1 or num == 3:
            rho = 1
        return rho

gamma_f = np.vectorize(__gamma)
c_f = np.vectorize(__c)
rho_f = np.vectorize(__rho)
tau10 = 1
tau20 = 6.
tau30 = 6.
tau12 = 2.5
tau01 = 1.0
tau02 = 3
tau23 = 1.13

y = 2.50
c = 0.14

X = np.linspace(0,2,100)
Y = np.linspace(0,2,100)

C,D = np.meshgrid(X,Y)

#decay of congestus:
R10 = gamma_f(D) / tau10
#decay of deep
R20 = gamma_f(C) / tau20
#decay of stratiform:
R30 = np.ones_like(C)/tau30
#conversion of congestus to deep:
R12 = c_f(y)*gamma_f(C,y,c) * (1 - gamma_f(D,y,-c)) / tau12
#conversion from deep to stratiform 
R23 = gamma_f(np.sqrt(C))/tau23
#birth of congestus
R01 = c_f(y)*gamma_f(C,y,c)*gamma_f(D,y,-c) / tau01
#birth of deep
R02 = c_f(y)*gamma_f(C,y,c)*(1-gamma_f(D,y,-c)) / tau02
y,c = 0,0

#decay of congestus:
R10_r = gamma_f(D) / tau10
#decay of deep
R20_r = gamma_f(C) / tau20
#decay of stratiform:
R30_r = np.ones_like(C)/tau30
#conversion of congestus to deep:
R12_r = c_f(y)*gamma_f(C,y,c) * (1 - gamma_f(D,y,-c)) / tau12
#conversion from deep to stratiform 
R23_r = gamma_f(np.sqrt(C))/tau23
#birth of congestus
R01_r = c_f(y)*gamma_f(C,y,c)*gamma_f(D,y,-c) / tau01
#birth of deep
R02_r = c_f(y)*gamma_f(C,y,c)*(1-gamma_f(D,y,-c)) / tau02

rho = np.array([rho_f(R01,R02,R12,R23,R30,R20,R10,num=x) for x in xrange(4)])
rho_r = np.array([rho_f(R01_r,R02_r,R12_r,R23_r,R30_r,R20_r,R10_r,num=x) for x in xrange(4)])

#rho = rho/rho.sum(axis=-1)[...,np.newaxis]
#rho_r =rho_r/rho_r.sum(axis=-1)[...,np.newaxis]

from matplotlib import pyplot as plt

fig = plt.figure()
nn = 0

ax = fig.add_subplot(2,3,1)
#V=np.linspace(0,0.9,8).round(2)
contour = ax.contourf(C,2 - D,rho[0],cmap='Blues')
fig.colorbar(contour,orientation='horizontal',pad=0.1,aspect=40)
ax.set_ylabel('Moisture')
ax.set_title('$\\rho_{1}$')


#V=np.linspace(0,0.4,8).round(2)
ax = fig.add_subplot(2,3,2)
contour = ax.contourf(C,2 - D,rho[1],cmap='Blues')
fig.colorbar(contour,orientation='horizontal',pad=0.1,aspect=40)
ax.set_title('$\\rho_{2}$')



#V=np.linspace(0,0.5,8).round(2)
ax = fig.add_subplot(2,3,3)
contour = ax.contourf(C,2 - D,rho[2],cmap='Blues')
fig.colorbar(contour,orientation='horizontal',pad=0.1,aspect=40)
ax.set_title('$\\rho_{3}$')


ax = fig.add_subplot(2,3,4)
#V=np.linspace(0,0.9,8).round(2)
contour = ax.contourf(C,2 - D,rho_r[0],cmap='Blues')
fig.colorbar(contour,orientation='horizontal',pad=0.1,aspect=40)
ax.set_ylabel('Moisture')
ax.set_title('$\\rho_{1}^*$')


#V=np.linspace(0,0.4,8).round(2)
ax = fig.add_subplot(2,3,5)
contour = ax.contourf(C,2 - D,rho_r[1],cmap='Blues')
fig.colorbar(contour,orientation='horizontal',pad=0.1,aspect=40)
ax.set_title('$\\rho_{2}^*$')



#V=np.linspace(0,0.5,8).round(2)
ax = fig.add_subplot(2,3,6)
contour = ax.contourf(C,2 - D,rho_r[2],cmap='Blues')
fig.colorbar(contour,orientation='horizontal',pad=0.1,aspect=40)
ax.set_title('$\\rho_{3}^*$')


plt.show()



exit()
V=np.linspace(0,0.4,8).round(2)
ax = fig.add_subplot(3,3,2)
contour = ax.contourf(C,2 - D,R02_r,V,cmap='Blues')
fig.colorbar(contour,orientation='horizontal',pad=0.1,aspect=40)
ax.set_title('$R_{02}$')


V=np.linspace(0,0.5,8).round(2)
ax = fig.add_subplot(3,3,3)
contour = ax.contourf(C,2 - D,R12_r,V,cmap='Blues')
fig.colorbar(contour,orientation='horizontal',pad=0.1,aspect=40)
ax.set_title('$R_{12}$')




V=np.linspace(0,0.7,8).round(2)
ax = fig.add_subplot(3,3,7)
contour = ax.contourf(C,2 - D,R23,V,cmap='Blues')
fig.colorbar(contour,orientation='horizontal',pad=0.1,aspect=40)
ax.set_title('$R_{23}^*/R_{23}$')
ax.set_xlabel('Instability')
ax.set_ylabel('Moisture')


V=np.linspace(0,0.9,8).round(2)
ax = fig.add_subplot(3,3,8)
contour = ax.contourf(C,2 - D,R10,V,cmap='Blues')
fig.colorbar(contour,orientation='horizontal',pad=0.1,aspect=40)
ax.set_title('$R_{10}^*/R_{10}$')
ax.set_xlabel('Instability')


V=np.linspace(0,0.2,8).round(2)
ax = fig.add_subplot(3,3,9)
contour = ax.contourf(C,2 - D,R20,V,cmap='Blues')
fig.colorbar(contour,orientation='horizontal',pad=0.1,aspect=40)
ax.set_title('$R_{20}^*/R_{20}$')
ax.set_xlabel('Instability')

plt.show()
