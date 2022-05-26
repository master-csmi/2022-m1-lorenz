import lorenz.parareal.parareal as lpp
import lorenz.parareal.utils as lpu
import lorenz.parareal.plot_sol as lp_plot
import numpy as np
import matplotlib.pyplot as plt
from mpi4py import MPI
import sys


def plot_all(var,n,nom_sol_ex,nom_sol,nom_pts):
    t_rk4,sol_rk4=lp_plot.read_sol_ex('data_parareal/'+nom_sol_ex)

    # Pour lire le fichier "sol[var].csv"
    t,sol,nb_iter=lp_plot.read_sol('data_parareal/'+nom_sol)

    # Pour lire le fichier "init_pt_x.csv"
    times,nb_pts,x0=lp_plot.read_init_pt('data_parareal/'+nom_pts,nb_iter)

    # Pour plot pour x
    lp_plot.plot_sol(var,t_rk4,sol_rk4[:,n],t,sol,nb_iter,times,x0)

    return t_rk4,sol_rk4[:,n],t,sol,nb_iter,times,nb_pts,x0


def erreur(solx,solx_exacte):
    return np.max(np.abs(solx-solx_exacte))

def val_k_n(k,n,nb_pts,solx,solx_exacte,x0):
    suite_nb = [-1]
    for i in range(len(nb_pts)):
        suite_nb.append(suite_nb[-1]+nb_pts[i])
    suite_nb[0] = 0
    suite_nb[-1] += 1

    nb1=suite_nb[n]
    nb2=suite_nb[n+1]

    sol_k = solx[:,k]
    sol_k_j = sol_k[nb1:nb2]
    sol_ex_k = solx_exacte[nb1:nb2]
    err = erreur(sol_k_j,sol_ex_k)
    diff = np.abs(x0[n,k]-sol_ex_k[0])
    return (diff,err)

def cvg(nom_sol_ex,nom_sol,nom_pts):
    t_ex,solx_exacte,t,solx,nb_iter,times,nb_pts,x0=plot_all('x',0,
            nom_sol_ex,nom_sol,nom_pts)

    nb_proc = len(nb_pts)
    dt_G=0.01
    diff = []
    err = []
    for k in range(nb_iter):
        diff_k = []
        err_k = []
        for n in range(nb_proc):
            diff_,err_=val_k_n(k,n,nb_pts,solx,solx_exacte,x0)
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

    plt.semilogy(tab_k,np.max(diff+err,axis=1),label="diff+err")
    plt.semilogy(tab_k,delta_t,label="delta_t^k")
    plt.legend()
    plt.show()


# mpiexec -n 3 python3 examples_parareal.py n_ex
 
nb_arg = len(sys.argv)

# display examples if no arguments
if(nb_arg==1):
    print("usage : [mpiexec -n 3] python3 examples_parareal.py n_ex")
    print("\t n_ex=0 : use rk4 on lorenz")
    print("\t n_ex=1 : use parareal method on lorenz system")
    print("\t n_ex=2 : use parareal method on oscillator")

comm = MPI.COMM_WORLD
P = comm.Get_size()
rank = comm.Get_rank()

# rk4 Lorenz

if(nb_arg>1 and sys.argv[1]=='0'):
    print("rk4 on lorenz")
    gamma=(10.,8./3,28.) #(σ,b,r)
    X0=[5.,5.,5.] #(x0,y0,z0)
    t0=0.
    T=200.
    dt=0.01

    sol = lpu.RK4(X0,dt,t0,T,lpu.lorenz,gamma)


## parareal

# Lorenz

if(nb_arg>1 and sys.argv[1]=='1'):
    if(rank==0):
        print("parareal method on lorenz system")
    gamma=(10.,8./3,28.) #(σ,b,r)
    X0=[5.,5.,5.] #(x0,y0,z0)
    t0=0.
    T=20.

    dt_G=0.1
    dt_F=0.01

    if(rank==0):
        lpu.delete_old_files_lorenz()
    sol=lpp.parareal_method(X0,t0,T,lpu.lorenz,lpu.RK4,lpu.csv_files_lorenz,
            dt_G,dt_F,gamma,True)

    if(rank==0):
        plot_all('x',0,'sol_rk4.csv','solx.csv','init_pt_x.csv')
        plot_all('y',1,'sol_rk4.csv','soly.csv','init_pt_y.csv')
        plot_all('z',2,'sol_rk4.csv','solz.csv','init_pt_z.csv')

# Oscillator

if(nb_arg>1 and sys.argv[1]=='2'):
    if(rank==0):
        print("parareal method on oscillator")
    gamma=(5.,-1./5.,np.pi/2.) #=(w0,x0,phi0)
    phi_0=[0.,1.] #(x0,y0,z0)
    t0=0.
    T=20.

    dt_G=0.01
    dt_F=0.001

    if(rank==0):
        lpu.delete_olf_files_oscillator()
    sol=lpp.parareal_method(phi_0,t0,T,lpu.oscillator,lpu.sol_ex,lpu.csv_files_oscillator,
            dt_G,dt_F,gamma,True)

    if(rank==0):
        cvg('sol_exacte_oscillator.csv','solx_oscillator.csv',
            'init_pt_oscillator.csv')