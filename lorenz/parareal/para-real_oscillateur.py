import numpy as np
import os
import pandas as panda
from mpi4py import MPI

# python3 para-real.py | grep truc
# mpiexec -n 3 python3 para-real_parallele.py

comm = MPI.COMM_WORLD
P = comm.Get_size()
rank = comm.Get_rank()

#### TO CLEAN ####

# to delete old files
def supprimer_ancien_fichiers():
    if(os.path.isfile("donnees_para_real/solx_oscillateur.csv")):
        os.remove("donnees_para_real/solx_oscillateur.csv")

    if(os.path.isfile("donnees_para_real/sol_exacte_oscillateur.csv")):
        os.remove("donnees_para_real/sol_exacte_oscillateur.csv")
    
    if(os.path.isfile("donnees_para_real/init_pt_oscillateur.csv")):
        os.remove("donnees_para_real/init_pt_oscillateur.csv")

    print("Anciens fichiers supprimés")

if(rank == 0): 
    supprimer_ancien_fichiers()

#### TO SOLVE ####

# function that represents the Lorenz system
def oscillateur_harmonique(t, X, gamma): #X=(x,y,z)
    (x,v)=X
    w0=gamma[0]

    f_1 = v
    f_2 = -w0**2*x
    return np.array([f_1,f_2])

# solution exacte
def sol_exacte(t,gamma):
    (w0,x0,phi0)=gamma
    sol_x=x0*np.cos(w0*t+phi0)
    return sol_x

def d_sol_exacte(t,gamma):
    (w0,x0,phi0)=gamma
    sol_v=-x0*w0*np.sin(w0*t+phi0) # dérivée de sol_x
    return sol_v

# Runge Kutta order 4
def RK4(X0,dt,t0,T,fct,gamma=None):
    X=np.array([X0])  
    
    t=t0 #=t_0
    while((t+dt)<=T or (np.isclose(t+dt,T))):
    # while(not np.isclose(t,T)):
        K1=fct(t, X[-1],gamma)
        K2=fct(t+dt/2., X[-1] + 1./2. * K1 * dt,gamma)
        K3=fct(t+dt/2., X[-1] + 1./2. * K2 * dt,gamma)
        K4=fct(t+dt, X[-1]+ K3 * dt,gamma)
        
        X=np.append(X,[X[-1]+ dt/6.* (K1+2.*K2+2.*K3+K4)],axis=0)
        t+=dt
    return X

def compute_time(t0,T,dt_G):
    # time between t_j and t_{j+1}
    dt_P = (T-t0)/P 
    # nb points
    nb_pts = dt_P//dt_G  
    # t_j exact
    times_exact = [t0]
    for j in range(1,P+1):
        times_exact.append(times_exact[-1] + dt_P)
    times_exact[-1]=T
    # t_j approach (to be a multiple of dt_G)
    times = [t0]
    for i in range(1,P+1):
        exact=times_exact[i-1]+dt_P
        approche=times[i-1]+dt_G*nb_pts
        if(exact-approche<=abs(exact-(approche+dt_G))): #arrondi inférieur
            times.append(times[i-1] + dt_G*nb_pts)
        else: #arrondi supérieur
            times.append(times[i-1] + dt_G*(nb_pts+1))
    times[-1]=T
    return times

# to check if initial points converge
def sol_converge(x0_k,x0_knext,eps=1e-15):
    print(np.max(np.abs(x0_knext-x0_k))/np.max(np.abs(x0_knext)))
    return (np.max(np.abs(x0_knext-x0_k))/np.max(np.abs(x0_knext))) < eps

#to create and write in the csv files
def csv_files(t0,T,dt_F,times,nb_tj,k,solution,init_pts,reshape_size,fct,gamma=None):
    # time between t0 and T for the fine integrator ##
    t = np.arange(t0,T+1e-6,dt_F)

    # init
    init_pts_x = panda.DataFrame(times[:-1],columns=['t'],dtype=np.float64)

    init_pts_x.insert(len(init_pts_x.columns),'nb_pts',nb_tj)

    solutions_x = panda.DataFrame(t,columns=['t'],dtype=np.float64)

    # write
    for k in range(k):
        X0_k = init_pts[k,:]
        init_pts_x.insert(len(init_pts_x.columns),'k='+str(k),X0_k[:,0])

        sol_k = solution[k,:]
        sol_k = sol_k.reshape((-1,reshape_size))
        solutions_x.insert(len(solutions_x.columns),'k='+str(k),sol_k[:,0])

    # to convert dataframe to csv files
    init_pts_x.to_csv('donnees_para_real/init_pt_oscillateur.csv')
    solutions_x.to_csv('donnees_para_real/solx_oscillateur.csv')

    # avec RK4
    sol = sol_exacte(t,gamma)

    sol_ex = panda.DataFrame(t,columns=['t'],dtype=np.float64)
    sol_ex.insert(len(sol_ex.columns),'x',sol)
    sol_ex.to_csv('donnees_para_real/sol_exacte_oscillateur.csv')

    print("Fichiers csv créés")

def compute_sol_k(fin,X0_k_j):
    sol_k_j = fin[1:].flatten()
    if(rank==0):
        temp = np.array([X0_k_j])
        temp = np.append(temp,sol_k_j)
        sol_k_j=temp
    nb_tj = np.array(comm.gather(len(sol_k_j), 0))
    sol_k=None
    if(rank==0):
        sol_k = np.empty(sum(nb_tj))
    comm.Gatherv(sendbuf=sol_k_j, recvbuf=(sol_k, nb_tj), root=0)

    return nb_tj,sol_k


# P = number of processing units ; j in {0,...,P} 
def parareal_method(X0_t0,t0,T,fct,dt_G,dt_F,gamma=None,write_csv=True): 
    # fine integrator
    def F(initial_point,t_j1,t_j2):
        return RK4(initial_point,dt_F,t_j1,t_j2,fct,gamma) 
    # coarse intergator
    def G(initial_point,t_j1,t_j2):
        return RK4(initial_point,dt_G,t_j1,t_j2,fct,gamma)

    ## use by rank 0 ##

    times = []
    X0_k = np.empty((P,len(X0_t0))) 
    
    grossier_k = np.empty((P-1,len(X0_t0))) #pas de valeur pour j=0

    if(rank == 0):   
        ## we initialize the time intervals for each processing units ##
        times=compute_time(t0,T,dt_G)

        #### Itération k=0 ####

        # 1ère étape : calculer les X0 pour chaque tj (en séquentiel)
        # On utilise pour cela l'intégrateur grossier
        X0_k[0] = X0_t0
        X0_k_j = X0_t0
        for j in range(1,P):
            X0_j = G(X0_k[j-1],times[j-1],times[j])[-1]
            X0_k[j] = X0_j
            comm.send(X0_k[j], dest=j, tag=j)   
            grossier_k[j-1] = X0_j         

    if(rank>0):
        X0_k_j = comm.recv(source=0, tag=rank)
    t_j = comm.scatter(times[:-1], root=0)
    t_jp = comm.scatter(times[1:], root=0)

    # 2ème étape : calculer la solution sur chaque sous-intervalle
    # On utilise l'intégrateur fin

    fin = F(X0_k_j,t_j,t_jp)
    
    fin_k_j = fin[-1]
    fin_k = np.array(comm.gather(fin_k_j, 0))

    if(write_csv):
        nb_tj,sol_k = compute_sol_k(fin,X0_k_j)
        if(rank==0):
            solution = [sol_k]
            init_pts = [X0_k]

    #### Itérations suivantes (jusqu'à ce que la solution converge) ####
    k=1
    converge=False

    # pour rentrer dans la boucle
    X0_kp = np.copy(X0_k)+np.ones((len(X0_k),len(X0_t0))) 
    while(not converge):      
        X0_kp = np.zeros((P,len(X0_t0)))
        
        if(rank==0):
            X0_kp[0] = X0_k[0]
            X0_k_j = X0_k[0]
            for j in range(1,P):
                grossier_k_j = G(X0_kp[j-1],times[j-1],times[j])[-1]
                X0_kp[j] = grossier_k_j + fin_k[j-1] - grossier_k[j-1]
                comm.send(X0_k[j], dest=j, tag=j)
                grossier_k[j-1] = grossier_k_j

        if(rank>0):
            X0_k_j = comm.recv(source=0, tag=rank)

        fin = F(X0_k_j,t_j,t_jp)
    
        fin_k_j = fin[-1]
        fin_k = np.array(comm.gather(fin_k_j, 0))

        if(write_csv):
            nb_tj,sol_k = compute_sol_k(fin,X0_k_j)
            if(rank==0):
                solution = np.concatenate((solution,[sol_k]),axis=0)
                init_pts = np.concatenate((init_pts,[X0_kp]),axis=0)

        if(rank==0):
            converge = sol_converge(X0_k,X0_kp)
        comm.Barrier()
        converge = comm.bcast(converge, root=0)

        if(rank==0):
            X0_k = np.copy(X0_kp) 

        k+=1

    if(rank==0):
        print(k," itérations")
        if(csv_files):
            csv_files(t0,T,dt_F,times,nb_tj/2,k,solution,init_pts,len(X0_t0),fct,gamma)
            return solution[k-1,:].reshape((-1,len(X0_t0)))
        else:
            nb_tj,sol_k = compute_sol_k(fin,X0_k_j)
            return sol_k
        
    MPI.Finalize()


gamma=(5.,-1./5.,np.pi/2.) #=(w0,x0,phi0)
phi_0=[0.,1.] #(x0,y0,z0)
t0=0.
T=50.

dt_G=0.01
dt_F=0.001

sol=parareal_method(phi_0,t0,T,oscillateur_harmonique,dt_G,dt_F,gamma,True)