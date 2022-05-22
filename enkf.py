import numpy as np
from filterpy.kalman import EnsembleKalmanFilter
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import solve_ivp
import time


def plot_1_fig(lorenz1,lorenz2,ind,tab_temps,tab_cov,diff,tab_x_y_z,labelx,labely):
    fig=plt.figure(figsize=(18,5))
    ax1 = fig.add_subplot(1,3,1)
    ax1.plot(lorenz1[3],lorenz1[ind],label="observation (0.01)")
    ax1.plot(lorenz2[3],lorenz2[ind],label="Model (0.1)")
    ax1.plot(tab_temps,tab_x_y_z[:,ind],color='blue',label="apres l'assimilation de donnée")
    plt.fill_between(tab_temps,tab_x_y_z[:,ind]+(tab_cov[:,ind]+diff),tab_x_y_z[:,ind]-(tab_cov[:,ind]+diff), color='blue',alpha=0.2)
    plt.xlabel(labelx)
    plt.ylabel(labely)
    ax1.legend()
    
def plot(lorenz1,lorenz2,tab_temps,tab_x_y_z):
    fig=plt.figure(figsize=(18,6))
    ax1 = fig.add_subplot(1,3,1)
    ax1.plot(lorenz1[3],lorenz1[0],label="observation (0.01)")
    ax1.plot(lorenz2[3],lorenz2[0],label="Model (0.1)")
    ax1.plot(tab_temps,tab_x_y_z[:,0],label="apres l'assimilation de donnée")
    plt.xlabel("t")
    plt.ylabel("x")
    ax1.legend()
    
    ax2 = fig.add_subplot(1,3,2)
    ax2.plot(lorenz1[3],lorenz1[1],label="observation (0.01)")
    ax2.plot(lorenz2[3],lorenz2[1],label="Model (0.1)")
    ax2.plot(tab_temps,tab_x_y_z[:,1],label="apres l'assimilation de donnée")
    plt.xlabel("t")
    plt.ylabel("y")
    ax1.legend()
    
    ax3 = fig.add_subplot(1,3,3)
    ax3.plot(lorenz1[3],lorenz1[2],label="observation (0.01)")
    ax3.plot(lorenz2[3],lorenz2[2],label="Model (0.1)")
    ax3.plot(tab_temps,tab_x_y_z[:,2],label="apres l'assimilation de donnée")
    plt.xlabel("t")
    plt.ylabel("z")
    ax1.legend()
def f_orcillateur(t,X_n,w): #X_n=(x_n,y_n,z_n)
    (u,v)=X_n
    
    f_1 = v
    f_2 =-w**2*u
    return np.array([f_1,f_2])
def RK4_Lorenz_harmonique(w,X0,N,T): #we have N+1 discretization points
    dt=T/N
    X = np.zeros( (N+1, len(X0)) )
    T_tab=np.zeros(N+1)
    X[0] = X0
    t=0. #=t_0
    for n in range(1,N+1):
        
        K1=f_orcillateur(t, X[n-1],w)
        K2=f_orcillateur(t+dt/2., X[n-1] + 1./2. * K1 * dt,w)
        K3=f_orcillateur(t+dt/2., X[n-1] + 1./2. * K2 * dt,w)
        K4=f_orcillateur(t+dt, X[n-1]+ K3 * dt,w)
        T_tab[n]=t+dt
        X[n]=X[n-1]+ dt/6.* (K1+2.*K2+2.*K3+K4)
        t+=dt
        
    return X[:,0],X[:,1],T_tab



def assimilation_oscillateur_harmonique():
    def hx(x):
       return  np.array([x[0]])


    def fx(x, dt,w):
        def f(t,X_n,w): #X_n=(x_n,y_n,z_n)
            (u,v)=X_n
            f_1 = v
            f_2 =-w**2*u
            return np.array([f_1,f_2])
        K1=f_orcillateur(t, x,w)
        K2=f_orcillateur(t+dt/2., x + 1./2. * K1 * dt,w)
        K3=f_orcillateur(t+dt/2., x + 1./2. * K2 * dt,w)
        K4=f_orcillateur(t+dt, x+ K3 * dt,w)
        X_next=x+ dt/6.* (K1+2.*K2+2.*K3+K4)
        return X_next

    w=2
    P = 2*np.pi/w
    dt= P/20
    T = 3*P
    N = int(round(T/dt))

    x = np.array([2,0])#(x0,y0,z0)
    P = np.eye(2) * 2.

    f = EnsembleKalmanFilter (x=x, P=P, dim_z=1, dt=dt, N=40,
             hx=hx, fx=lambda x,dt:fx(x,dt,w))

    std_noise = np.eye(1)*0.001
    f.R *= std_noise # matrice de cov associer a la mesure
    f.Q=np.eye(2)*0.01
    #f.Q = Q_discrete_white_noise(2, dt, .01) #bruit blanc centree en 0

    def read_sensor(t,w):
         2*np.cos(w*t)
    t=0
    tab_etat=[]
    tab_temps=[]
    while (t<T):
        z = read_sensor(t,w)
        f.predict()
        f.update(z)
        tab_etat.append(f.x[0])
        tab_temps.append(t)
        t=t+dt
    return np.array(tab_etat),np.array(tab_temps)
        
def f(t_n,X_n,σ, b, r):
    (x,y,z)=X_n
    
    f_1 = σ*(y-x)
    f_2 = x*(r-z)-y
    f_3 = x*y-b*z
    return np.array([f_1,f_2,f_3])

def RK4_Lorenz(γ,X0,N,T): #we have N+1 discretization points
    (σ,b,r)=γ
    dt=T/N
    X = np.zeros( (N+1, len(X0)) )
    T_tab=np.zeros(N+1)
    X[0] = X0
    T_tab[0]=0
    
    t=0. #=t_0
    for n in range(1,N+1):
        
        K1=f(t, X[n-1],σ,b,r)
        K2=f(t+dt/2., X[n-1] + 1./2. * K1 * dt,σ,b,r)
        K3=f(t+dt/2., X[n-1] + 1./2. * K2 * dt,σ,b,r)
        K4=f(t+dt, X[n-1]+ K3 * dt,σ,b,r)
        T_tab[n]=t+dt
        X[n]=X[n-1]+ dt/6.* (K1+2.*K2+2.*K3+K4)
        t+=dt
        
    return X[:,0],X[:,1],X[:,2],T_tab

def assimilation_donnée(x,read_sensor,P,Q,R,T,dimz,dt,N,nb_echantillon,hx,fx,γ):
    f = EnsembleKalmanFilter (x=x, P=P, dim_z=dimz, dt=dt, N=nb_echantillon,
             hx=hx, fx=lambda x,dt:fx(x,t,dt,γ[0],γ[1],γ[2]))

    #std_noise = np.eye(3)*0.001
    f.R = R # matrice de cov associer a la mesure
    f.Q =Q   #bruit blanc  centree en 0
    t=0
    index=N
    tab_etat=[]
    tab_temps=[]
    tab_cov=[]
    tab_etat.append(f.x)
    tab_temps.append(t)
    tab_cov.append(f.P_post.diagonal())
    while (t<T-dt):
        z = read_sensor(index)
        f.predict()
        f.update(z)
        diag_cov=f.P_post.diagonal()
        tab_cov.append(diag_cov)
        index+=N
        t=t+dt
        tab_etat.append(f.x)
        tab_temps.append(t)
    return(np.array(tab_etat),np.array(tab_temps),np.array(tab_cov))

def fx(x,t, dt,σ,b,r):
    def f(t_n,X_n,σ, b, r): #X_n=(x_n,y_n,z_n)
        (x,y,z)=X_n
    
        f_1 = σ*(y-x)
        f_2 = x*(r-z)-y
        f_3 = x*y-b*z
        return np.array([f_1,f_2,f_3])
    K1=f(t, x,σ,b,r)
    K2=f(t+dt/2., x + 1./2. * K1 * dt,σ,b,r)
    K3=f(t+dt/2., x + 1./2. * K2 * dt,σ,b,r)
    K4=f(t+dt, x+ K3 * dt,σ,b,r)
    X_next=x+ dt/6.* (K1+2.*K2+2.*K3+K4)
    return X_next
    

