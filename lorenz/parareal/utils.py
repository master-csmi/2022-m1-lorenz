import os
import numpy as np
import pandas as panda
import matplotlib.pyplot as plt

# to delete old files
def delete_old_files_lorenz():
    """Delete old csv files created for the lorenz system.
    """    
    if(os.path.isfile("data_parareal/solx.csv")):
        os.remove("data_parareal/solx.csv")
    if(os.path.isfile("data_parareal/soly.csv")):
        os.remove("data_parareal/soly.csv")
    if(os.path.isfile("data_parareal/solz.csv")):
        os.remove("data_parareal/solz.csv")

    if(os.path.isfile("data_parareal/sol_rk4.csv")):
        os.remove("data_parareal/sol_rk4.csv")
    
    if(os.path.isfile("data_parareal/init_pt_x.csv")):
        os.remove("data_parareal/init_pt_x.csv")
    if(os.path.isfile("data_parareal/init_pt_y.csv")):
        os.remove("data_parareal/init_pt_y.csv")
    if(os.path.isfile("data_parareal/init_pt_z.csv")):
        os.remove("data_parareal/init_pt_z.csv")

    print("Old files deleted")

def delete_olf_files_oscillator():
    """Delete old csv files created for the oscillator.
    """    
    if(os.path.isfile("data_parareal/solx_oscillateur.csv")):
        os.remove("data_parareal/solx_oscillateur.csv")

    if(os.path.isfile("data_parareal/sol_exacte_oscillateur.csv")):
        os.remove("data_parareal/sol_exacte_oscillateur.csv")
    
    if(os.path.isfile("data_parareal/init_pt_oscillateur.csv")):
        os.remove("data_parareal/init_pt_oscillateur.csv")

    print("Old files deleted")

def compute_time(t0,T,dt_G,P):
    """Compute time used to the system resolution.

    Args:
        t0 (float): Starting time.
        T (float): Finish time.
        dt_G (float): Coarse time step.
        P (int): Number of interval (= number of processes).

    Returns:
        list: Time used to the system resolution.
    """    
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
    """_summary_

    Args:
        x0_k (_type_): _description_
        x0_knext (_type_): _description_
        eps (_type_, optional): _description_. Defaults to 1e-15.

    Returns:
        _type_: _description_
    """    
    return (np.max(np.abs(x0_knext-x0_k))/np.max(np.abs(x0_knext))) < eps

#to create and write in the csv files
def csv_files_lorenz(t0,T,dt_F,times,nb_tj,k,solution,init_pts,reshape_size,fct,fct_res,gamma=None):
    # time between t0 and T for the fine integrator ##
    t = np.arange(t0,T+1e-6,dt_F)

    # init
    init_pts_x = panda.DataFrame(times[:-1],columns=['t'],dtype=np.float64)
    init_pts_x.insert(len(init_pts_x.columns),'nb_pts',nb_tj)
    init_pts_y = panda.DataFrame(times[:-1],columns=['t'],dtype=np.float64)
    init_pts_y.insert(len(init_pts_y.columns),'nb_pts',nb_tj)
    init_pts_z = panda.DataFrame(times[:-1],columns=['t'],dtype=np.float64)
    init_pts_z.insert(len(init_pts_z.columns),'nb_pts',nb_tj)

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
    init_pts_x.to_csv('data_parareal/init_pt_x.csv')
    init_pts_y.to_csv('data_parareal/init_pt_y.csv')
    init_pts_z.to_csv('data_parareal/init_pt_z.csv')
    solutions_x.to_csv('data_parareal/solx.csv')
    solutions_y.to_csv('data_parareal/soly.csv')
    solutions_z.to_csv('data_parareal/solz.csv')

    # avec RK4
    sol = fct_res(init_pts[0,:][0],dt_F,t0,T,fct,gamma)

    sol_rk4 = panda.DataFrame(t,columns=['t'],dtype=np.float64)
    sol_rk4.insert(len(sol_rk4.columns),'x',sol[:,0])
    sol_rk4.insert(len(sol_rk4.columns),'y',sol[:,1])
    sol_rk4.insert(len(sol_rk4.columns),'z',sol[:,2])
    sol_rk4.to_csv('data_parareal/sol_rk4.csv')

    print("Fichiers csv créés")

#to create and write in the csv files
def csv_files_oscillator(t0,T,dt_F,times,nb_tj,k,solution,init_pts,reshape_size,fct,fct_res,gamma=None):
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
    init_pts_x.to_csv('data_parareal/init_pt_oscillator.csv')
    solutions_x.to_csv('data_parareal/solx_oscillator.csv')

    # with exact solution
    sol = fct_res(t,gamma)

    sol_exacte = panda.DataFrame(t,columns=['t'],dtype=np.float64)
    sol_exacte.insert(len(sol_exacte.columns),'x',sol)
    sol_exacte.to_csv('data_parareal/sol_exacte_oscillator.csv')

    print("Fichiers csv créés")


# Runge Kutta order 4
def RK4(X0,dt,t0,T,fct,gamma=None):
    X=np.array([X0])  
    
    t=t0 #=t_0
    while((t+dt)<=T or (np.isclose(t+dt,T))):
        K1=fct(t, X[-1],gamma)
        K2=fct(t+dt/2., X[-1] + 1./2. * K1 * dt,gamma)
        K3=fct(t+dt/2., X[-1] + 1./2. * K2 * dt,gamma)
        K4=fct(t+dt, X[-1]+ K3 * dt,gamma)
        
        X=np.append(X,[X[-1]+ dt/6.* (K1+2.*K2+2.*K3+K4)],axis=0)
        t+=dt
    return X

def erreur(solx,solx_exacte):
    return np.max(np.abs(solx-solx_exacte))

def E_j_k(j,k,nb_pts,solx,solx_exacte,x0):
    suite_nb = [-1]
    for i in range(len(nb_pts)):
        suite_nb.append(suite_nb[-1]+nb_pts[i])
    suite_nb[0] = 0
    suite_nb[-1] += 1

    nb1=suite_nb[j]
    nb2=suite_nb[j+1]

    sol_k = solx[:,k]
    sol_k_j = sol_k[nb1:nb2]
    sol_ex_k = solx_exacte[nb1:nb2]
    err = erreur(sol_k_j,sol_ex_k)
    diff = np.abs(x0[j,k]-sol_ex_k[0])
    return (diff,err)

def cvg(t_ex,solx_exacte,t,solx,nb_iter,times,nb_pts,x0):
    nb_proc = len(nb_pts)
    dt_G=0.01
    diff = []
    err = []
    for k in range(nb_iter):
        diff_k = []
        err_k = []
        for j in range(nb_proc):
            diff_,err_=E_j_k(j,k,nb_pts,solx,solx_exacte,x0)
            diff_k.append(diff_)
            err_k.append(err_)
        diff.append(diff_k)
        err.append(err_k)

    diff = np.array(diff)
    err = np.array(err)

    delta_t = []
    for k in range(nb_iter):
        delta_t.append(dt_G**k)

    tab_k=np.arange(0,np.shape(diff)[0],1)

    plt.semilogy(tab_k,np.max(diff+err,axis=1),label="E_j_k")
    plt.semilogy(tab_k,delta_t,label="delta_t^k")
    plt.xticks(tab_k,tab_k)
    plt.title("Convergence based on iterations")
    plt.legend()
    plt.show()
