import numpy as np
import matplotlib.pyplot as plt
import csv
import math

# Pour lire le fichier "sol_rk4.csv"
def read_sol_ex(nom_fichier):
    t_rk4=[]
    sol_rk4=np.array([])
    with open(nom_fichier, newline='') as csvfile:
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

def read_sol(nom_fichier):
    t=[]
    sol_k=np.array([])
    with open(nom_fichier, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for i,row in enumerate(reader):
            if(i==0):
                nb_iter=len(row)-2
            else:
                t.append(float(row[1]))
                sol_k = np.append(sol_k,[float(item) for item in row[2:]])
        sol_k = np.reshape(sol_k,[-1,nb_iter])

    return t,sol_k,nb_iter

def read_init_pt(nom_fichier,nb_iter):
    times=[]
    nb_pts=[]
    x0=[]
    with open(nom_fichier, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        next(reader, None)
        for row in reader:
            times.append(float(row[1]))
            nb_pts.append(math.floor(float(row[2])))
            x0 = np.append(x0,[float(item) for item in row[3:]])
        x0 = np.reshape(x0,[-1,nb_iter])
    return times,nb_pts,x0

def plot_sol(entree,t_rk4,sol_rk4,t,sol,nb_iter,times,pt0):
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
    t_rk4,sol_rk4=read_sol_ex('data_parareal/'+nom_sol_ex)

    # Pour lire le fichier "sol[var].csv"
    t,sol,nb_iter=read_sol('data_parareal/'+nom_sol)

    # Pour lire le fichier "init_pt_x.csv"
    times,nb_pts,x0=read_init_pt('data_parareal/'+nom_pts,nb_iter)

    # Pour plot pour x
    plot_sol(var,t_rk4,sol_rk4[:,n],t,sol,nb_iter,times,x0)

    return t_rk4,sol_rk4[:,n],t,sol,nb_iter,times,nb_pts,x0

def plot_3D(solution,X0):
    x=solution[:,0]
    y=solution[:,1]
    z=solution[:,2]
    fig=plt.figure(figsize=(18,5))
    
    #plot 3D
    ax3 = fig.add_subplot(1,3,3, projection='3d')
    ax3.plot(x,y,z)
    ax3.plot(X0[0],X0[1],X0[2],"r.")
    
    ax3.set_xlabel("x")
    ax3.set_ylabel("y")
    ax3.set_zlabel("z")
    # ax3.view_init(28,-95)
    
    plt.title("3D representation of the solution")
    plt.show()