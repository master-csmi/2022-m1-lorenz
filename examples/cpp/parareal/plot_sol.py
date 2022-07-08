import numpy as np
import matplotlib.pyplot as plt
import csv
import math

def read_sol(file_name):
    t=[]
    sol_k=np.array([])
    with open(file_name, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for i,row in enumerate(reader):
            if(i==0):
                dim = len(row)-1
            t.append(float(row[0]))
            sol_k = np.append(sol_k,[float(item) for item in row[1:]])
        sol_k = np.reshape(sol_k,[-1,dim])

    return t,sol_k,dim

def read_init_pt(file_name):
    times=[]
    x0=np.array([])
    with open(file_name, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for i,row in enumerate(reader):
            if(i==0):
                dim = len(row)-1
            times.append(float(row[0]))
            x0 = np.append(x0,[float(item) for item in row[1:]])
        x0 = np.reshape(x0,[-1,dim])

    return times,x0,dim

t,sol_k,dim = read_sol("data/solution_0.csv")
times,x0,dim = read_init_pt("data/solution_0.csv")

fig=plt.figure()

plt.plot(t,sol_k[:,0],label="x")
plt.plot(times,x0[:,0],".k",label="x0")

plt.show()



# print(t)
# print()
# print(sol_k)
# print()
# print(dim)
