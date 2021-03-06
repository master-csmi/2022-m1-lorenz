import numpy as np
from filterpy.kalman import EnsembleKalmanFilter
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import solve_ivp
import time


def plot_1_fig(lorenz1,lorenz2,ind,tab_temps,tab_cov,diff,tab_x_y_z,labelx,labely):
    """
    Parameters:
        lorenz1: observations,an array of vectors where each vector is of size 3, representing [x,y,z] 
        lorenz2: model,an array of vectors where each vector is of size 3, representing [x,y,z] 
        ind:a real representing the index we want to plot (x, y or z)
        tab_temps:real type array representing time
        tab_cov:an array of vectors where each vector is of size 3, representing variances 
         (diag of the cov matrix associated to the state obtained through data assimilation)
        diff:represent the distance to create our standart deviation variance
        tab_x_y_z:represente an array of vectors where each vector is of size 3,array after data assimilation,
         representing [x,y,z] 
        labelx:legend of the x-curve
        labely:legend of the y-curve
    Returns:
        return nothing but plot 3 curves as a function of ind(either x,y or z) as a function of time

    
    """
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
    """
    Parameters:
        lorenz1: observations,an array of vectors where each vector is of size 3, representing [x,y,z] 
        lorenz2: model,an array of vectors where each vector is of size 3, representing [x,y,z] 
        tab_temps:real type array representing time
        tab_x_y_z:represente an array of vectors where each vector is of size 3,array after data assimilation,
         representing [x,y,z] 
    Returns:
        return nothing but makes three different plot, in each one we will have 
         3 curves (observations,model and stade after data assimilation) as a function of time

    
    """
    fig=plt.figure(figsize=(18,6))
    ax1 = fig.add_subplot(1,3,1)
    ax1.plot(lorenz1[3],lorenz1[0],label="observation (0.01)")
    ax1.plot(lorenz2[3],lorenz2[0],label="Model (0.1)")
    ax1.plot(tab_temps,tab_x_y_z[:,0],label="apres l'assimilation de donnée")
    plt.xlabel("t")
    plt.ylabel("x")
    ax1.legend(loc='upper right')
    
    ax2 = fig.add_subplot(1,3,2)
    ax2.plot(lorenz1[3],lorenz1[1],label="observation (0.01)")
    ax2.plot(lorenz2[3],lorenz2[1],label="Model (0.1)")
    ax2.plot(tab_temps,tab_x_y_z[:,1],label="apres l'assimilation de donnée")
    plt.xlabel("t")
    plt.ylabel("y")
    ax1.legend(loc='upper right')
    
    ax3 = fig.add_subplot(1,3,3)
    ax3.plot(lorenz1[3],lorenz1[2],label="observation (0.01)")
    ax3.plot(lorenz2[3],lorenz2[2],label="Model (0.1)")
    ax3.plot(tab_temps,tab_x_y_z[:,2],label="apres l'assimilation de donnée")
    plt.xlabel("t")
    plt.ylabel("z")
    ax1.legend(loc='upper right')
def f_orcillateur(t,X_n,w):
    """
    Parameters:
        X_n:un tableau de reel de taille 2 
        w: a real that represente the first parameter of the hamonic oscillator
    
    """
    (u,v)=X_n
    
    f_1 = v
    f_2 =-w**2*u
    return np.array([f_1,f_2])
def RK4_harmonique(w,X0,N,T): 
    """
    Parameters: 
        w: a real that represente the first parameter of the hamonic oscillator
        X0:a size 2 array with the initial constition, initial point.
        N: a real that represents the number of discritisation
        T: a real that represents the time interval
    Returns:
        this fonction return resolution with RK4 
        X[:,0]: array of x
        X[:,1]: array of y
        T_tab: array of the time
    
    """
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
def euler_explicit_harmonique(w,X0,N,T):
    """
    Parameters: 
        w: a real that represente the first parameter of the hamonic oscillator
        X0:a size 2 array with the initial constition, initial point.
        N: a real that represents the number of discritisation
        T: a real that represents the time interval
    Returns:
        this fonction return resolution with euler
        X[:,0]: array of x
        X[:,1]: array of y
        T_tab: array of the time
    
    
    """
    dt=T/N
    X = np.zeros( (N+1, len(X0)) )
    T_tab=np.zeros(N+1)
    X[0] = X0
    t=0. 
    for n in range(1,N+1):
        X[n]=[X[n-1,0]+(dt)*X[n-1,1],X[n-1,1]-w**2*X[n-1,0]*(dt)]
        T_tab[n]=t+dt
        t+=dt
    return X[:,0],X[:,1],T_tab

 

def f(t_n,X_n,σ, b, r):
    (x,y,z)=X_n
    
    f_1 = σ*(y-x)
    f_2 = x*(r-z)-y
    f_3 = x*y-b*z
    return np.array([f_1,f_2,f_3])

def RK4_Lorenz(γ,X0,N,T): 
    """
    Parameters: 
        γ: an array of real that represente the three parameter of the lorenz system (σ, b, r)
        X0:a size 3 array with the initial constition, initial point (x,y,z)
        N: a real that represents the number of discritisation
        T: a real that represents the time interval
    Returns:
        this fonction return resolution with RK4
        X[:,0]: array of x
        X[:,1]: array of y
        X[:,2]: array of z
        T_tab: array of the time
    
    
    """
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

def fx(x,t, dt,γ):
    def f(t_n,X_n,γ): 
        (x,y,z)=X_n
        f_1 = γ[0]*(y-x)
        f_2 = x*(γ[2]-z)-y
        f_3 = x*y-γ[1]*z
        return np.array([f_1,f_2,f_3])
    K1=f(t, x,γ)
    K2=f(t+dt/2., x + 1./2. * K1 * dt,γ)
    K3=f(t+dt/2., x + 1./2. * K2 * dt,γ)
    K4=f(t+dt, x+ K3 * dt,γ)
    X_next=x+ dt/6.* (K1+2.*K2+2.*K3+K4)
    return X_next

def assimilation_donnée(x,read_sensor,P,Q,R,T,dimz,dt,N,nb_echantillon,hx,fx,γ):
    """
    Parameters: 
        x:a size 3 array with the initial constition, initial point (x,y,z)
        read_sensor: a function that will read the observation every time
        P:covariance matrix  associated with the forecast state error
        Q:covariance matrix associated with the model error
        R:covariance matrix associated with the observations error
        T:a real that represents the time interval
        dimz: a real that represent the dimention of the observations
        dt:a real that represent the time for each discretisation
        N: a real that represents the index difference between the observations and the model
        nb_echantillon: a real that represents the index difference between the observations and the model
        hx:Measurement function. Convert state x into a measurement
        fx:State transition function
        γ: an array of real that represente the three parameter of the lorenz system (σ, b, r)
        N: a real that represents the number of discritisation
    Returns:
        this function do the data assimilation using the ensemble kalman filter
        -tab_etat: an array with the state mean after data assimilation
        -tab_temps: an array with the time
        -tab_cov: an array with the diagonal of the state covarience matrix pour each time
    """
    f = EnsembleKalmanFilter (x=x, P=P, dim_z=dimz, dt=dt, N=nb_echantillon,hx=hx, fx=lambda x,dt:fx(x,t,dt,γ))
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
    while (t<T-dt):#
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



