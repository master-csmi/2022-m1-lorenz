import numpy as np
import matplotlib.pyplot as plt
import csv
import math

# Pour lire le fichier "sol_exacte.csv"
t_exact=[]
solx_exacte=[]
with open('donnees_para_real/sol_exacte_oscillateur.csv', newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    next(reader, None)
    for row in reader:
        t_exact.append(float(row[1]))
        solx_exacte.append(float(row[2]))

def read_sol(nom_fichier):
    t=[]
    sol_k=np.array([])
    with open(nom_fichier, newline='') as csvfile:
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



def plot(entree,t_exact,sol_exacte,t,sol,times,pt0):
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

        ax.plot(t_exact,sol_exacte,linewidth=2,label="exacte")
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



### Pour x ####

# Pour lire le fichier "solx.csv"
t,solx,nb_iter=read_sol('donnees_para_real/solx_oscillateur.csv')

# Pour lire le fichier "init_pt_x.csv"
times,nb_pts,x0=read_init_pt('donnees_para_real/init_pt_oscillateur.csv',nb_iter)

# Pour plot pour x
# plot('x',t_exact,solx_exacte,t,solx,times,x0)



#### CONVERGENCE

def erreur(solx,solx_exacte):
    return np.max(np.abs(solx-solx_exacte))

def cvg(k,n,nb_pts,solx,solx_exacte,x0):
    suite_nb = [0]
    for i in range(len(nb_pts)):
        suite_nb.append(suite_nb[-1]+nb_pts[i])
    nb1=suite_nb[n]
    nb2=suite_nb[n+1]

    sol_k = solx[:,k]
    sol_k_j = sol_k[nb1:nb2]
    sol_ex_k = solx_exacte[nb1:nb2]
    err = erreur(sol_k_j,sol_ex_k)
    diff = np.abs(x0[n,k]-sol_ex_k[0])
    return (diff,err)

nb_proc = len(nb_pts)
dt_G=0.01
diff = []
err = []
for k in range(nb_iter):
    diff_k = []
    err_k = []
    for n in range(nb_proc):
        diff_,err_=cvg(k,n,nb_pts,solx,solx_exacte,x0)
        diff_k.append(diff_)
        err_k.append(err_)
    diff.append(diff_k)
    err.append(err_k)

print("diff :",np.array(diff))
print("err :",np.array(err))

max_diff=np.max(diff,axis=1)
max_err=np.max(err,axis=1)
print("max_diff : ",max_diff)
print("max_err : ",max_err)

delta_t = []
for k in range(nb_iter):
    delta_t.append(dt_G**k)

tab_k=np.arange(0,np.shape(diff)[0],1)
print(tab_k)

plt.semilogy(tab_k,max_diff,label="|Y_k^n-y(T^n)|")
plt.semilogy(tab_k,max_err,label="max sur T^n T^{n+1}")
plt.semilogy(tab_k,delta_t,label="delta_t^k")
plt.legend()
plt.show()
