import lorenz.parareal.parareal as lpp
import lorenz.parareal.utils as lpu
import lorenz.parareal.plot_sol as lp_plot
import numpy as np
import matplotlib.pyplot as plt
from mpi4py import MPI
import sys

# mpiexec -n 3 python3 examples_parareal.py n_ex
 
nb_arg = len(sys.argv)

# display examples if no arguments
if(nb_arg==1):
    print("usage : [mpiexec -n 3] python3 examples_parareal.py n_ex")
    print("\t n_ex=0 : use rk4 on lorenz")
    print("\t n_ex=1 : use parareal method on lorenz system")
    print("\t n_ex=2 : use parareal method on oscillator")



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

    sol = lpu.RK4(X0,dt,t0,T,lorenz,gamma)

    # to show the 3D_sol
    # lp_plot.plot_3D(sol,X0)

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

    write=True

    if(rank==0):
        lpu.delete_old_files_lorenz()
    lpp.parareal_method(X0,t0,T,lorenz,lpu.RK4,lpu.csv_files_lorenz,
            dt_G,dt_F,gamma,write)

    if(rank==0 and write):
        lp_plot.plot_all('x',0,'sol_rk4.csv','solx.csv','init_pt_x.csv')
        lp_plot.plot_all('y',1,'sol_rk4.csv','soly.csv','init_pt_y.csv')
        lp_plot.plot_all('z',2,'sol_rk4.csv','solz.csv','init_pt_z.csv')

# Oscillator

if(nb_arg>1 and sys.argv[1]=='2'):
    if(rank==0):
        print("parareal method on oscillator")
    gamma=(5.,-1./5.,np.pi/2.) #=(w0,x0,phi0)
    pt_init=[0.,1.] #(x(0),v(0))
    t0=0.
    T=20.

    dt_G=0.01
    dt_F=0.001

    write=True

    if(rank==0):
        lpu.delete_olf_files_oscillator()
    lpp.parareal_method(pt_init,t0,T,oscillator,sol_ex,lpu.csv_files_oscillator,
            dt_G,dt_F,gamma,write)

    if(rank==0 and write):
        t_ex,solx_exacte,t,solx,nb_iter,times,nb_pts,x0=lp_plot.plot_all('x',0,
            'sol_exacte_oscillator.csv','solx_oscillator.csv',
            'init_pt_oscillator.csv')
        lpu.cvg(t_ex,solx_exacte,t,solx,nb_iter,times,nb_pts,x0)