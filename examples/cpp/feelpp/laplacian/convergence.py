import csv
import numpy as np
import matplotlib.pyplot as plt

filename = "circle/values.csv"

tab_h = []
normes_L2 = []
normes_H1 = []

# To read the csv file
file = open(filename)
csvreader = csv.reader(file)
for i,row in enumerate(csvreader):
    if(i!=0):
        tab_h += [(float)(row[0])]
        normes_L2 += [(float)(row[1])]
        normes_H1 += [(float)(row[2])]
file.close()
    
pente_L2 = np.polyfit(np.log(tab_h), np.log(normes_L2), 1)[0]
plt.loglog(tab_h,normes_L2,label="Norm L2 : slope="+str(np.round(pente_L2,2)))
print("Slope pour la norme L2 :",pente_L2)

pente_H1 = np.polyfit(np.log(tab_h), np.log(normes_H1), 1)[0]
plt.loglog(tab_h,normes_H1,label="Norm H1 : slope="+str(np.round(pente_H1,2)))
print("Slope pour la norme H1 :",pente_H1)

plt.legend()
plt.savefig("circle/cvg_laplacian.png") 
plt.show()