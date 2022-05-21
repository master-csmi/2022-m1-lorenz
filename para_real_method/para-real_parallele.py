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

if(rank == 0): 
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
    while(not np.isclose(t,T)):
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
def sol_converge(x0_k,x0_knext,eps=1e-6):
    return (np.max(np.abs(x0_knext-x0_k)) < eps)

#to create and write in the csv files
def csv_files(t0,T,dt_F,times,k,solution,init_pts,reshape_size,fct,gamma=None):
    # time between t0 and T for the fine integrator ##
    t = np.arange(t0,T+1e-6,dt_F)

    # init
    init_pts_x = panda.DataFrame(times[:-1],columns=['t'],dtype=np.float64)
    init_pts_y = panda.DataFrame(times[:-1],columns=['t'],dtype=np.float64)
    init_pts_z = panda.DataFrame(times[:-1],columns=['t'],dtype=np.float64)

    solutions_x = panda.DataFrame(t,columns=['t'],dtype=np.float64)
    solutions_y = panda.DataFrame(t,columns=['t'],dtype=np.float64)
    solutions_z = panda.DataFrame(t,columns=['t'],dtype=np.float64)

    # write
    for k in range(k):
        X0_k = init_pts[k,:]
        init_pts_x.insert(len(init_pts_x.columns),'k='+str(k),X0_k[:,0])
        init_pts_y.insert(len(init_pts_y.columns),'k='+str(k),X0_k[:,1])
        init_pts_z.insert(len(init_pts_z.columns),'k='+str(k),X0_k[:,2])

        sol_k = solution[k,:]
        sol_k = sol_k.reshape((-1,reshape_size))
        solutions_x.insert(len(solutions_x.columns),'k='+str(k),sol_k[:,0])
        solutions_y.insert(len(solutions_y.columns),'k='+str(k),sol_k[:,1])
        solutions_z.insert(len(solutions_z.columns),'k='+str(k),sol_k[:,2])

    # to convert dataframe to csv files
    init_pts_x.to_csv('donnees_para_real/init_pt_x.csv')
    init_pts_y.to_csv('donnees_para_real/init_pt_y.csv')
    init_pts_z.to_csv('donnees_para_real/init_pt_z.csv')
    solutions_x.to_csv('donnees_para_real/solx.csv')
    solutions_y.to_csv('donnees_para_real/soly.csv')
    solutions_z.to_csv('donnees_para_real/solz.csv')

    # avec RK4
    sol = RK4(init_pts[0,:][0],dt_F,t0,T,fct,gamma)

    sol_rk4 = panda.DataFrame(t,columns=['t'],dtype=np.float64)
    sol_rk4.insert(len(sol_rk4.columns),'x',sol[:,0])
    sol_rk4.insert(len(sol_rk4.columns),'y',sol[:,1])
    sol_rk4.insert(len(sol_rk4.columns),'z',sol[:,2])
    sol_rk4.to_csv('donnees_para_real/sol_rk4.csv')

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

    return sol_k


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
        sol_k = compute_sol_k(fin,X0_k_j)
        if(rank==0):
            solution = [sol_k]
            init_pts = [X0_k]

    #### Itérations suivantes (jusqu'à ce que la solution converge) ####
    k=1
    converge=False

    # pour rentrer dans la boucle
    X0_kp = np.copy(X0_k)+np.ones((len(X0_k),len(X0))) 
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
            sol_k = compute_sol_k(fin,X0_k_j)
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
        csv_files(t0,T,dt_F,times,k,solution,init_pts,len(X0_t0),fct,gamma)
        print(k," itérations")

    MPI.Finalize()


gamma=(10.,8./3,28.) #(σ,b,r)
X0=[5.,5.,5.] #(x0,y0,z0)
t0=0.
T=20.

dt_G=0.01
dt_F=0.001

parareal_method(X0,t0,T,lorenz,dt_G,dt_F,gamma,True)