import numpy as np
import matplotlib.pyplot as plt
import csv
import math

def read_sol_ex(file_name):
    """Read the file "file_name" where is the exact solution.

    Args:
        file_name (str): Name of the csv file where the exaction solution is.

    Returns:
        tuple: 
            Discretisation time between t0 and T for the exact solution
            and exact solution.
    """    
    t_rk4=[]
    sol_rk4=np.array([])
    with open(file_name, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        size=0
        for i,row in enumerate(reader):
            if(i==0):
                size=len(row)-2
            else:
                t_rk4.append(float(row[1]))
                sol_rk4 = np.append(sol_rk4,[float(item) for item in row[2:]])
        sol_rk4 = np.reshape(sol_rk4,[-1,size])
    
    return t_rk4,sol_rk4

def read_sol(file_name):
    """Read the file "file_name" where is the calculated solution.

    Args:
        file_name (str): Name of the csv file where the calculated solution is.

    Returns:
        tuple: 
            Discretisation time between t0 and T for the calculated solution
            and calculated solution.
    """ 
    t=[]
    sol_k=np.array([])
    with open(file_name, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for i,row in enumerate(reader):
            if(i==0):
                nb_iter=len(row)-2
            else:
                t.append(float(row[1]))
                sol_k = np.append(sol_k,[float(item) for item in row[2:]])
        sol_k = np.reshape(sol_k,[-1,nb_iter])

    return t,sol_k,nb_iter

def read_init_pt(file_name,dim):
    """Read the file "file_name" where are the initial points for each iteration.

    Args:
        file_name (str): Name of the csv file where the initial points are.
        dim (int): System dimension.

    Returns:
        tuple: 
            Times for the initial poins, number of value of the solution 
            for each process and initial points for each process and each iteration.
    """    
    times=[]
    nb_pts=[]
    x0=[]
    with open(file_name, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        next(reader, None)
        for row in reader:
            times.append(float(row[1]))
            nb_pts.append(math.floor(float(row[2])))
            x0 = np.append(x0,[float(item) for item in row[3:]])
        x0 = np.reshape(x0,[-1,dim])
    return times,nb_pts,x0

def plot_sol(entree,t_rk4,sol_rk4,t,sol,nb_iter,times,pt0):
    """Plot the solution of the parareal method for each iteration.

    Args:
        entree (str): Variable.
        t_rk4 (list): Discretisation time between t0 and T for the exact solution.
        sol_rk4 (numpy.ndarray): Exact solution.
        t (list): Discretisation time between t0 and T for the calculated solution.
        sol (numpy.ndarray): Calculated solution.
        nb_iter (int): Number of iteration (until the solution converge).
        times (list): Times used to the system resolution.
        pt0 (numpy.ndarray): Initial points for each process and each iteration.
    """   
    fig,axs=plt.subplots(int(np.ceil(nb_iter/3)),3,sharex=True,sharey=True,figsize=(10,6))

    fig.suptitle(entree+'(t)')
    fig.subplots_adjust(wspace=0.1,hspace=0.5)

    for k in range(nb_iter):
        ligne=k//3
        colonne=k%3

        if(nb_iter<=3):
            ax=axs[colonne]
        else:
            ax=axs[ligne,colonne]
        
        ax.set_title("k="+str(k))

        ax.plot(t_rk4,sol_rk4,linewidth=2,label="exact")
        ax.plot(t,sol[:,k],linewidth=0.8,label=entree,alpha=0.9)
        ax.plot(times,pt0[:,k],".k",label=entree+"0")

    if(nb_iter<=3):
        ax0=axs[0]
    else:
        ax0=axs[0,0]

    lines_labels = [ax0.get_legend_handles_labels()]
    lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
    fig.legend(lines, labels)

    plt.show()

def plot_all(var,n,nom_sol_ex,nom_sol,nom_pts):
    """Read the exact solution, calculated solution and initial points.
    Plot them all.

    Args:
        var (str): Variable.
        n (int): Position of the variable (between 0 and dimension).
        nom_sol_ex (str): Name of the csv file where the exaction solution is.
        nom_sol (str): Name of the csv file where the calculated solution is.
        nom_pts (str): Name of the csv file where the initial points are.

    Returns:
        tuple: 
            Discretisation time between t0 and T for the exact solution,
            exact solution, discretisation time between t0 and T for
            the calculated solution, calculated solution (only the variable
            concerned), number of iteration (until the solution converge),
            times for the initial poins, number of value of the solution
            for each process and initial points for each process
            and each iteration.
    """    
    t_rk4,sol_rk4=read_sol_ex('data_parareal/'+nom_sol_ex)
    t,sol,nb_iter=read_sol('data_parareal/'+nom_sol)
    times,nb_pts,x0=read_init_pt('data_parareal/'+nom_pts,nb_iter)

    plot_sol(var,t_rk4,sol_rk4[:,n],t,sol,nb_iter,times,x0)

    return t_rk4,sol_rk4[:,n],t,sol,nb_iter,times,nb_pts,x0

def plot_3D(solution,X0):
    """Plot in 3 dimension the solution and the initial point

    Args:
        solution (numpy.ndarray): Exact solution.
        X0 (numpy.ndarray): Initial point of the system.
    """    
    x=solution[:,0]
    y=solution[:,1]
    z=solution[:,2]
    fig=plt.figure(figsize=(18,5))

    ax3 = fig.add_subplot(1,3,3, projection='3d')
    ax3.plot(x,y,z)
    ax3.plot(X0[0],X0[1],X0[2],"r.")
    
    ax3.set_xlabel("x")
    ax3.set_ylabel("y")
    ax3.set_zlabel("z")
    
    plt.title("3D representation of the solution")
    plt.show()