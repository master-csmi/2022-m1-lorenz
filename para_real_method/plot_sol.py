import numpy as np
import matplotlib.pyplot as plt
import csv

k = 6

plt.title('Solution itération k='+str(k))
plt.xlabel("t")
plt.ylabel("x")

# Pour lire le fichier "solution_rk4.csv"
t=[]
x=[]
with open('donnees_para_real/solution_rk4.csv', newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    for row in reader:
        t.append(float(row[0]))
        x.append(float(row[1]))

plt.plot(t,x,linewidth=3,label="sol_rk4")

# Pour lire le fichier "solution_k.csv"
t=[]
sol_k=[]
with open('donnees_para_real/solution_'+str(k)+'.csv', newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    for row in reader:
        t.append(float(row[0]))
        sol_k.append(float(row[1]))

plt.plot(t,sol_k,label="sol_k")


# Pour lire le fichier "X0_k.csv"
t=[]
X0_k=[]
with open('donnees_para_real/X0_'+str(k)+'.csv', newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    for row in reader:
        t.append(float(row[0]))
        X0_k.append(float(row[1]))

plt.plot(t,X0_k,".",label="X0")

plt.legend()
plt.show()