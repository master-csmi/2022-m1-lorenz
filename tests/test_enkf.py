import pytest
from lorenz.enkf import enkf
import numpy as np
from filterpy.kalman import EnsembleKalmanFilter
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import solve_ivp
import time



def test_enkf_1():
    w=2
    x = np.array([2,0])
    Pe = 2*np.pi/w
    dt= Pe/20
    T = 3*Pe
    N = int(round(T/dt))
    def fx_rk4(x,t, dt,w):
        def f(t,X_n,w): #X_n=(x_n,y_n,z_n)
            (u,v)=X_n
            f_1 = v
            f_2 =-w**2*u
            return np.array([f_1,f_2])
        K1=f(t, x,w)
        K2=f(t+dt/2., x + 1./2. * K1 * dt,w)
        K3=f(t+dt/2., x + 1./2. * K2 * dt,w)
        K4=f(t, x+ K3 * dt,w)
        X_next=x+ dt/6.* (K1+2.*K2+2.*K3+K4)
        return X_next

    P = np.eye(2) * 0.0
    R = np.eye(1)*0.001 # matrice de cov associer a la mesure
    Q=np.eye(2)*0.0
    def hx(x):
        return  np.array([x[0]])


    def read_sensor(t):
            2*np.cos(2*t)
    oscillateur=enkf.RK4_harmonique(w,x,N,T)

    tab_etat_oscillateur,tab_temps_oscillateur,tab_cov=enkf.assimilation_donnée(x,read_sensor,P,Q,R,T,1,dt,1,40,hx,fx_rk4,w)
    taille_tab_etat=len(tab_etat_oscillateur)
    compteur=0
    for i in range (taille_tab_etat):
        if ((tab_etat_oscillateur[i][0]-oscillateur[0][i])<1e-5):
            compteur+=1;
    if(compteur==taille_tab_etat):
        print("test oke pour enkf")
    else:
        print("test ne passe pas pour enkf")
