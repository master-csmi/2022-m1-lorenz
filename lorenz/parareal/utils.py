import os
import numpy as np
import pandas as panda

# to delete old files
def delete_old_files_lorenz():
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

    print("Old files deleted")

def delete_olf_files_oscillator():
    if(os.path.isfile("donnees_para_real/solx_oscillateur.csv")):
        os.remove("donnees_para_real/solx_oscillateur.csv")

    if(os.path.isfile("donnees_para_real/sol_exacte_oscillateur.csv")):
        os.remove("donnees_para_real/sol_exacte_oscillateur.csv")
    
    if(os.path.isfile("donnees_para_real/init_pt_oscillateur.csv")):
        os.remove("donnees_para_real/init_pt_oscillateur.csv")

    print("Anciens fichiers supprimés")

def compute_time(t0,T,dt_G,P):
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
    init_pts_x.to_csv('donnees_para_real/init_pt_x.csv')
    init_pts_y.to_csv('donnees_para_real/init_pt_y.csv')
    init_pts_z.to_csv('donnees_para_real/init_pt_z.csv')
    solutions_x.to_csv('donnees_para_real/solx.csv')
    solutions_y.to_csv('donnees_para_real/soly.csv')
    solutions_z.to_csv('donnees_para_real/solz.csv')

    # avec RK4
    sol = fct_res(init_pts[0,:][0],dt_F,t0,T,fct,gamma)

    sol_rk4 = panda.DataFrame(t,columns=['t'],dtype=np.float64)
    sol_rk4.insert(len(sol_rk4.columns),'x',sol[:,0])
    sol_rk4.insert(len(sol_rk4.columns),'y',sol[:,1])
    sol_rk4.insert(len(sol_rk4.columns),'z',sol[:,2])
    sol_rk4.to_csv('donnees_para_real/sol_rk4.csv')

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
    init_pts_x.to_csv('donnees_para_real/init_pt_oscillateur.csv')
    solutions_x.to_csv('donnees_para_real/solx_oscillateur.csv')

    # with exact solution
    sol = fct_res(t,gamma)

    sol_exacte = panda.DataFrame(t,columns=['t'],dtype=np.float64)
    sol_exacte.insert(len(sol_exacte.columns),'x',sol)
    sol_exacte.to_csv('donnees_para_real/sol_exacte_oscillateur.csv')

    print("Fichiers csv créés")


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

# function that represents the Lorenz system
def lorenz(t, X, gamma): #X=(x,y,z)
    (sigma,b,r)=gamma
    (x,y,z)=X
    
    f_1 = sigma*(y-x)
    f_2 = x*(r-z)-y
    f_3 = x*y-b*z
    return np.array([f_1,f_2,f_3])

# function that the oscillator
def oscillator(t, X, gamma): #X=(x,y,z)
    (x,v)=X
    w0=gamma[0]

    f_1 = v
    f_2 = -w0**2*x
    return np.array([f_1,f_2])

# exact solution
def sol_ex(t,gamma):
    (w0,x0,phi0)=gamma
    sol_x=x0*np.cos(w0*t+phi0)
    return sol_x