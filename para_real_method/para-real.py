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


def sol_converge(y_k,y_knext,eps=1e-20): #y_knext = y_{k=1}
    return (np.max(np.abs(y_knext-y_k)) < eps)

def ecriture_sol_k(t,solk,k):
    f_sol_k = open("donnees_para_real/solution_"+str(k)+".csv", "w+")
    for i in range(len(solk)):
        f_sol_k.write(str(t[i])+", "+str(solk[i])+"\n")
    f_sol_k.close()    


def ecriture_csv(nom_fichier,x,y):
    fc = open(nom_fichier, "w+")
    for i in range(len(x)):
        fc.write(str(x[i])+", "+str(y[i])+"\n")
    fc.close()   

# P = number of processing units ; j in {0,...,P} 
def parareal_method(X0_t0,t0,T,P,fct,dt_G,dt_F,gamma=None): 
    # time between t_j and t_{j+1}
    dt_P = (T-t0)/P 

    # time between 0 and T for the fine integrator
    t = np.arange(t0,T+1e-6,(dt_F))

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

    #### Itération k=0 ####

    X0_k = np.empty((P+1,len(X0_t0))) 
    X0_k[0] = X0_t0

    grossier = np.empty((P,len(X0_t0))) #pas de valeur pour j=0
    fin = np.empty((P,len(X0_t0))) #pas de valeur pour j=0

    # 1ère étape : calculer les X0 pour chaque tj (en séquentiel)
    # On utilise pour cela l'intégrateur grossier
    for j in range(1,P+1):
        X0_j = G(X0_k[j-1],times[j-1],times[j])[-1]
        grossier[j-1] = X0_j
        X0_k[j] = X0_j

    ecriture_csv("donnees_para_real/X0_0.csv",times,X0_k[:,0])

    # 2ème étape : calculer la solution sur chaque sous-intervalle
    # On utilise l'intégrateur fin
    solk = np.array([X0_t0]) # solk = sol0
    for j in range(1,P+1):
        sol0_j = F(X0_k[j-1],times[j-1],times[j])
        
        fin[j-1] = sol0_j[-1]
        solk = np.concatenate((solk,sol0_j[1:]))
    
    ecriture_csv("donnees_para_real/solution_0.csv",t,solk[:,0])

    #### Itérations suivantes (jusqu'à ce que la solution converge) ####
    k=1
    # pour rentrer dans la boucle
    solkp = np.copy(solk)+np.ones((len(solk),len(X0))) 
    while(not sol_converge(solk,solkp)):
        solkp = np.array([X0_k[0]])
        for j in range(1,P+1):
            g_k=G(X0_k[j-1],times[j-1],times[j])[-1]
            X0_k[j] = g_k + (fin[j-1]-grossier[j-1])

            grossier[j-1] = g_k

            solk_j = F(X0_k[j-1],times[j-1],times[j])
            
            fin[j-1] = solk_j[-1]

            solkp = np.concatenate((solkp,solk_j[1:]))
        
        solk = np.copy(solkp)

        ecriture_csv("donnees_para_real/solution_"+str(k)+".csv",t,solk[:,0])
        ecriture_csv("donnees_para_real/X0_"+str(k)+".csv",times,X0_k[:,0])

        k+=1

# avec la méthode para-real
gamma=(10.,5./3,9./10) #(σ,b,r)
X0=[-1.,4.,5.] #(x0,y0,z0)
t0=0.
T=20
P=5
dt_G=0.1
dt_F=0.01

parareal_method(X0,t0,T,P,lorenz,dt_G,dt_F,gamma)

# avec RK4
dt = 0.001
t = np.arange(t0,T+1e-9,dt)
x=RK4(X0,dt,t0,T,lorenz,gamma)[:,0]

ecriture_csv("donnees_para_real/solution_rk4.csv",t,x)