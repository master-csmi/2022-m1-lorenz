import numpy as np
from mpi4py import MPI
from pyrsistent import v
from lorenz.parareal import utils

# python3 para-real.py | grep truc
# mpiexec -n 3 python3 para-real_parallele.py

comm = MPI.COMM_WORLD
P = comm.Get_size()
rank = comm.Get_rank()

def compute_sol_k(fin,X0_k_j):
    """Compute the solution at the current iteration and do a gather.

    Args:
        fin (numpy.ndarray): Fine solution on the current process at the 
            current iteration.
        X0_k_j (numpy.ndarray): Initial point on the current process at the 
            current iteration.

    Returns:
        tuple: Ndarray with the number of values of solution for each process
                and the entire solution between t0 and T.
    """    
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
def parareal_method(X0_t0,t0,T,prob,fct_res,fct_write,dt_G,dt_F,gamma=None,
        write_csv=True): 
    """Parareal method.

    Args:
        X0_t0 (list): Initial point.
        t0 (float): Starting time.
        T (float): Finish time.
        prob (function): Function which represent the ODE.
        fct_res (function): Function which represent the integrator.
        fct_write (function): Function which write the solution in csv files.
        dt_G (float): Coarse time step.
        dt_F (float): Time step.
        gamma (tuple, optional): System parameters. Defaults to None.
        write_csv (bool, optional): Boolean to true if we want to write the 
            solutions in a csv file. Defaults to True.

    Returns:
        numpy.ndarray: Solution at the last iteration.
    """    
    # fine integrator
    def F(initial_point,t_j1,t_j2):
        return utils.RK4(initial_point,dt_F,t_j1,t_j2,prob,gamma) 
    # coarse intergator
    def G(initial_point,t_j1,t_j2):
        return utils.RK4(initial_point,dt_G,t_j1,t_j2,prob,gamma)

    ## use by rank 0 ##

    times = []
    X0_k = np.empty((P,len(X0_t0))) 
    
    grossier_k = np.empty((P-1,len(X0_t0))) #pas de valeur pour j=0

    if(rank == 0):   
        ## we initialize the time intervals for each processing units ##
        times=utils.compute_time(t0,T,dt_G,P)

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
        print(type(nb_tj))
        print(nb_tj)
        print(type(sol_k))
        print(sol_k)
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
            converge = utils.sol_converge(X0_k,X0_kp)
        comm.Barrier()
        converge = comm.bcast(converge, root=0)

        if(rank==0):
            X0_k = np.copy(X0_kp) 

        k+=1

    if(rank==0):
        print(k," itérations")
        if(write_csv):
            fct_write(t0,T,dt_F,times,nb_tj/len(X0_t0),k,solution,
                init_pts,len(X0_t0),prob,fct_res,gamma)
            return solution[k-1,:].reshape((-1,len(X0_t0)))
        
    MPI.Finalize()
