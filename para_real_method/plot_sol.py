import numpy as np
import matplotlib.pyplot as plt
import csv

# Pour lire le fichier "sol_rk4.csv"
t_rk4=[]
solx_rk4=[]
soly_rk4=[]
solz_rk4=[]
with open('donnees_para_real/sol_rk4.csv', newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    next(reader, None)
    for row in reader:
        t_rk4.append(float(row[1]))
        solx_rk4.append(float(row[2]))
        soly_rk4.append(float(row[3]))
        solz_rk4.append(float(row[4]))

def read_sol(nom_fichier):
    t=[]
    sol_k=np.array([])
    with open('donnees_para_real/solx.csv', newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        legend=True
        nb_val=0
        for row in reader:
            if(legend):
                nb_iter=len(row)-2
                legend=False
            else:
                nb_val+=1
                t.append(float(row[1]))
                sol_k = np.append(sol_k,[float(item) for item in row[2:]])
        sol_k = np.reshape(sol_k,[-1,nb_iter])
    return t,sol_k,nb_iter

def read_init_pt(nom_fichier,nb_iter):
    times=[]
    x0=[]
    with open(nom_fichier, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        next(reader, None)
        for row in reader:
            times.append(float(row[1]))
            x0 = np.append(x0,[float(item) for item in row[2:]])
        x0 = np.reshape(x0,[-1,nb_iter])
    return times,x0

def plot(entree,t_rk4,sol_rk4,t,sol,times,pt0):
    fig,axs=plt.subplots(int(np.ceil(nb_iter/3)),3,sharex=True,sharey=True,figsize=(10,6))

    fig.suptitle(entree+'(t)')
    fig.subplots_adjust(wspace=0.1,hspace=0.5)

    for k in range(nb_iter):
        ligne=k//3
        colonne=k%3
        axs[ligne,colonne].set_title("k="+str(k))

        axs[ligne,colonne].plot(t_rk4,sol_rk4,linewidth=2,label="rk4")
        axs[ligne,colonne].plot(t,sol[:,k],linewidth=0.8,label=entree,alpha=0.9)
        axs[ligne,colonne].plot(times,pt0[:,k],".k",label=entree+"0")

    lines_labels = [axs[0,0].get_legend_handles_labels()]
    lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
    fig.legend(lines, labels)

    plt.show()



### Pour x ####

# Pour lire le fichier "solx.csv"
t,solx,nb_iter=read_sol('donnees_para_real/solx.csv')

# Pour lire le fichier "init_pt_x.csv"
times,x0=read_init_pt('donnees_para_real/init_pt_x.csv',nb_iter)

# Pour plot pour x
plot('x',t_rk4,solx_rk4,t,solx,times,x0)



### Pour y ####

# Pour lire le fichier "soly.csv"
t,soly,nb_iter=read_sol('donnees_para_real/soly.csv')

# Pour lire le fichier "init_pt_y.csv"
times,y0=read_init_pt('donnees_para_real/init_pt_y.csv',nb_iter)

# Pour plot pour y
plot('y',t_rk4,soly_rk4,t,soly,times,y0)



### Pour z ####

# Pour lire le fichier "solz.csv"
t,solz,nb_iter=read_sol('donnees_para_real/solz.csv')

# Pour lire le fichier "init_pt_z.csv"
times,z0=read_init_pt('donnees_para_real/init_pt_z.csv',nb_iter)

# Pour plot pour z
plot('z',t_rk4,solz_rk4,t,solz,times,z0)