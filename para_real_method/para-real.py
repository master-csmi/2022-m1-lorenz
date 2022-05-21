import numpy as np
import os
import pandas as panda
import mpi4py

#### TO CLEAN ####

# to delete old files
def supprimer_ancien_fichiers():
    if(os.path.isfile("donnees_para_real/solx.csv")):
        os.remove("donnees_para_real/solx.csv")
    if(os.path.isfile("donnees_para_real/soly.csv")):
        os.remove("donnees_para_real/soly.csv")
    if(os.path.isfile("donnees_para_real/solz.csv")):
        os.remove("donnees_para_real/solz.csv")

    if(os.path.isfile("donnees_para_real/sol_rk4.csv")):
        os.remove("donnees_para_real/sol_rk4.csv")
    
    if(os.path.isfile("donnees_para_real/init_pt_x.csv")):
        os.remove("donnees_para_real/init_pt_x.csv")
    if(os.path.isfile("donnees_para_real/init_pt_y.csv")):
        os.remove("donnees_para_real/init_pt_y.csv")
    if(os.path.isfile("donnees_para_real/init_pt_z.csv")):
        os.remove("donnees_para_real/init_pt_z.csv")

    print("Anciens fichiers supprimés")

supprimer_ancien_fichiers()

#### TO SOLVE ####

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

# to check if initial points converge
def sol_converge(x0_k,x0_knext,eps=1e-6):
    return (np.max(np.abs(x0_knext-x0_k)) < eps)

# P = number of processing units ; j in {0,...,P} 
# on part du principe que dt_G est un multiple de dt_F 
def parareal_method(X0_t0,t0,T,P,fct,dt_G,dt_F,gamma=None): 
    # time between t_j and t_{j+1}
    dt_P = (T-t0)/P 
    # nb points
    nb_pts = dt_P//dt_G

    # time between 0 and T for the fine integrator
    t = np.arange(t0,T+1e-6,dt_F)

    # we initialize the time intervals for each processing units
    times_exact = [t0]
    for j in range(1,P+1):
        times_exact.append(times_exact[-1] + dt_P)
    times_exact[-1]=T

    times = [t0]
    for i in range(1,P+1):
        exact=times_exact[i-1]+dt_P
        approche=times[i-1]+dt_G*nb_pts
        if(exact-approche<=abs(exact-(approche+dt_G))): #arrondi inférieur
            times.append(times[i-1] + dt_G*nb_pts)
        else: #arrondi supérieur
            times.append(times[i-1] + dt_G*(nb_pts+1))
    times[-1]=T

    solutions_x = panda.DataFrame(t,columns=['t'],dtype=np.float64)
    solutions_y = panda.DataFrame(t,columns=['t'],dtype=np.float64)
    solutions_z = panda.DataFrame(t,columns=['t'],dtype=np.float64)
    init_pts_x = panda.DataFrame(times,columns=['t'],dtype=np.float64)
    init_pts_y = panda.DataFrame(times,columns=['t'],dtype=np.float64)
    init_pts_z = panda.DataFrame(times,columns=['t'],dtype=np.float64)

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

    init_pts_x.insert(len(init_pts_x.columns),'k=0',X0_k[:,0])
    init_pts_y.insert(len(init_pts_y.columns),'k=0',X0_k[:,1])
    init_pts_z.insert(len(init_pts_z.columns),'k=0',X0_k[:,2])

    # 2ème étape : calculer la solution sur chaque sous-intervalle
    # On utilise l'intégrateur fin
    solk = np.array([X0_t0]) # solk = sol0
    for j in range(1,P+1):
        sol0_j = F(X0_k[j-1],times[j-1],times[j])
        
        fin[j-1] = sol0_j[-1]
        solk = np.concatenate((solk,sol0_j[1:]))

    solutions_x.insert(len(solutions_x.columns),'k=0',solk[:,0])
    solutions_y.insert(len(solutions_y.columns),'k=0',solk[:,1])
    solutions_z.insert(len(solutions_z.columns),'k=0',solk[:,2])

    #### Itérations suivantes (jusqu'à ce que la solution converge) ####
    k=1
    converge=False

    print(X0_k)

    # pour rentrer dans la boucle
    X0_kp = np.copy(X0_k)+np.ones((len(X0_k),len(X0))) 
    while(not converge):
        solkp = np.array([X0_k[0]])
        X0_kp = np.copy(X0_k)
        X0_kp[0] = X0_k[0]
        for j in range(1,P+1):
            print("j :",j)
            print("g_k :",G(X0_kp[j-1],times[j-1],times[j])[-1])
            print("fin[j-1] :",fin[j-1])
            print("grossier[j-1]",grossier[j-1])
            g_k=G(X0_kp[j-1],times[j-1],times[j])[-1]
            X0_kp[j] = g_k + (fin[j-1]-grossier[j-1])

            grossier[j-1] = g_k

            solk_j = F(X0_kp[j-1],times[j-1],times[j])
            
            fin[j-1] = solk_j[-1]

            solkp = np.concatenate((solkp,solk_j[1:]))
        
        converge = sol_converge(X0_k,X0_kp)
        print(X0_kp)
        X0_k = np.copy(X0_kp)

        init_pts_x.insert(len(init_pts_x.columns),'k='+str(k),X0_kp[:,0])
        init_pts_y.insert(len(init_pts_y.columns),'k='+str(k),X0_kp[:,1])
        init_pts_z.insert(len(init_pts_z.columns),'k='+str(k),X0_kp[:,2])

        solutions_x.insert(len(solutions_x.columns),'k='+str(k),solkp[:,0])
        solutions_y.insert(len(solutions_y.columns),'k='+str(k),solkp[:,1])
        solutions_z.insert(len(solutions_z.columns),'k='+str(k),solkp[:,2])

        k+=1

    init_pts_x.to_csv('donnees_para_real/init_pt_x.csv')
    init_pts_y.to_csv('donnees_para_real/init_pt_y.csv')
    init_pts_z.to_csv('donnees_para_real/init_pt_z.csv')
    solutions_x.to_csv('donnees_para_real/solx.csv')
    solutions_y.to_csv('donnees_para_real/soly.csv')
    solutions_z.to_csv('donnees_para_real/solz.csv')

    print(k," itérations")


gamma=(10.,8./3,28.) #(σ,b,r)
X0=[5.,5.,5.] #(x0,y0,z0)
t0=0.
T=20.

# avec la méthode para-real
P=3
dt_G=0.1
dt_F=0.01

parareal_method(X0,t0,T,P,lorenz,dt_G,dt_F,gamma)

# avec RK4
dt = 0.01
t = np.arange(t0,T+1e-9,dt)
sol = RK4(X0,dt,t0,T,lorenz,gamma)

sol_rk4 = panda.DataFrame(t,columns=['t'],dtype=np.float64)
sol_rk4.insert(len(sol_rk4.columns),'x',sol[:,0])
sol_rk4.insert(len(sol_rk4.columns),'y',sol[:,1])
sol_rk4.insert(len(sol_rk4.columns),'z',sol[:,2])
sol_rk4.to_csv('donnees_para_real/sol_rk4.csv')