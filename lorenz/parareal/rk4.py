import numpy as np


# function that represents the Lorenz system
def lorenz(t, X, gamma): #X=(x,y,z)
    (sigma,b,r)=gamma
    (x,y,z)=X
    
    f_1 = sigma*(y-x)
    f_2 = x*(r-z)-y
    f_3 = x*y-b*z
    return np.array([f_1,f_2,f_3])

# Runge Kutta order 4
def RK4(X0,dt,t0,T,fct,gamma=None):
    X=np.array([X0])  
    
    t=t0 #=t_0
    while((t+dt)<=T or (np.isclose(t+dt,T))):
        K1=fct(t, X[-1],gamma)
        K2=fct(t+dt/2., X[-1] + 1./2. * K1 * dt,gamma)
        K3=fct(t+dt/2., X[-1] + 1./2. * K2 * dt,gamma)
        K4=fct(t+dt, X[-1]+ K3 * dt,gamma)
        
        X=np.append(X,[X[-1]+ dt/6.* (K1+2.*K2+2.*K3+K4)],axis=0) #X au t suivant
        t+=dt
    return X

gamma=(10.,8./3,28.) #(σ,b,r)
X0=[5.,5.,5.] #(x0,y0,z0)
t0=0.
T=200.
dt=0.001

sol = RK4(X0,dt,t0,T,lorenz,gamma)