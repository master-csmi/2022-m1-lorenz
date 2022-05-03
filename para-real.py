from cProfile import label
import numpy as np
import mpi4py

import matplotlib.pyplot as plt

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
    while(np.abs(T-t)>1e-9):
        K1=fct(t, X[-1],gamma)
        K2=fct(t+dt/2., X[-1] + 1./2. * K1 * dt,gamma)
        K3=fct(t+dt/2., X[-1] + 1./2. * K2 * dt,gamma)
        K4=fct(t+dt, X[-1]+ K3 * dt,gamma)
        
        X=np.append(X,[X[-1]+ dt/6.* (K1+2.*K2+2.*K3+K4)],axis=0)
        t+=dt
    return X


def diff_solutions(y_k,y_knext,eps=1e-9): #y_knext = y_{k=1}
    return np.max(np.abs(y_knext-y_k))

# P = number of processing units ; j in {0,...,P-1} 
def parareal_method(X0_t0,t0,T,P,fct,dt_G,dt_F,gamma=None): 
    # time between t_j and t_{j+1}
    dt_P = (T-t0)/P 
    # # dt for the fine integrator
    # NF = 100
    # dt_F = dt_P/NF
    # # dt for the coarse integrator
    # NG = 15
    # dt_G = dt_P/NG

    # we initialize the time intervals for each processing units
    times = [t0]
    for _ in range(1,P+1):
        times.append(times[-1] + dt_P)

    # fine integrator
    def F(initial_point,t_j1,t_j2):
        return RK4(initial_point,dt_F,t_j1,t_j2,fct,gamma) 
    # coarse intergator
    def G(initial_point,t_j1,t_j2):
        return RK4(initial_point,dt_G,t_j1,t_j2,fct,gamma)

    # Itération k=0

    # 1ère étape : calculer les X0 pour chaque tj (en séquentiel)
    # On utilise pour cela l'intégrateur grossier
    X0 = np.array([X0_t0]) 
    for j in range(1,P+1):
        X0_j = G(X0[-1],times[j-1],times[j])[-1]
        X0 = np.append(X0,[X0_j],axis=0)

    # 2ème étape : calculer la solution sur chaque sous-intervalle
    # On utilise l'intégrateur fin
    solk = np.array([X0_t0]) # solk = sol0
    for j in range(1,P+1):
        sol0_j = F(X0[j-1],times[j-1],times[j])

        timesk_j = np.arange(times[j-1],times[j]+1e-6,dt_F)
        plt.plot(timesk_j,sol0_j[:,0],label=f"fin pour j={j}")
        
        solk = np.append(solk,[sol0_j[-1]],axis=0)

    # Itérations suivantes (jusqu'à ce que la solution converge)

    #while
    #X0 = np.array([])
    #for j in range(1,P+1):
    #    X0_j = 
    #    X0 = np.append(X0,[X0_j],axis=0)

    x,y,z=X0[:,0],X0[:,1],X0[:,2]
    plt.plot(times,x,'.-k',label="grossier")
    x,y,z=solk[:,0],solk[:,1],solk[:,2]
    plt.plot(times,x,'.g',label="fin")
    plt.legend()
    plt.show()

    
gamma=(10.,8./3,9./10) #(σ,b,r)
X0=(2.,5.,5.) #(x0,y0,z0)
t0=0.
T=1
P=5
dt_G=0.1
dt_F=0.01

parareal_method(X0,t0,T,P,lorenz,dt_G,dt_F,gamma)
