import numpy as np
import matplotlib.pyplot as plt
import csv
import math

k = 0

def read_sol_ex(filename):  
    t_rk4=[]
    sol_rk4=np.array([])
    with open(filename, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for i,row in enumerate(reader):
            if(i==0):
                dim = len(row)-1
            t_rk4.append(float(row[0]))
            sol_rk4 = np.append(sol_rk4,[float(item) for item in row[1:]])
        sol_rk4 = np.reshape(sol_rk4,[-1,dim])
    
    return t_rk4,sol_rk4

def read_sol(filename):
    t=[]
    sol_k=np.array([])
    with open(filename, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for i,row in enumerate(reader):
            if(i==0):
                dim = len(row)-1
            t.append(float(row[0]))
            sol_k = np.append(sol_k,[float(item) for item in row[1:]])
        sol_k = np.reshape(sol_k,[-1,dim])

    return t,sol_k,dim

def read_init_pt(filename):
    times=[]
    x0=np.array([])
    with open(filename, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for i,row in enumerate(reader):
            if(i==0):
                dim = len(row)-1
            times.append(float(row[0]))
            x0 = np.append(x0,[float(item) for item in row[1:]])
        x0 = np.reshape(x0,[-1,dim])

    return times,x0,dim

t,sol_k,dim = read_sol("data/solution_"+str(k)+".csv")
times,x0,dim = read_init_pt("data/init_pts_"+str(k)+".csv")
t_rk4,sol_rk4 = read_sol_ex("data/solution_ex.csv")


fig=plt.figure()

plt.plot(t,sol_k[:,0],label="x")
plt.plot(t_rk4,sol_rk4[:,0],label="x exact")
plt.plot(times,x0[:,0],".k",label="x0")

plt.title("k="+str(k))
plt.legend(loc='upper right')
plt.savefig("data/sol_"+str(k)+".png") 
plt.show()