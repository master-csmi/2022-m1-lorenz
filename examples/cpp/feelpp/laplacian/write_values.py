import sys
import csv
import os

# python3 write_values.py n_ex

nb_arg = len(sys.argv)

write_file_name = "circle/values.csv"
values_file_name = "/home/flecourtier/feelppdb/laplacian.e/np_1/logs/values.csv"

# To delete old csv file
if(nb_arg>1 and sys.argv[1]=='1'):
    try:
        os.remove(write_file_name)
    except OSError as e:
        print(e)

h = sys.argv[2]

# to read the norm values
file = open(values_file_name)
csvreader = csv.reader(file)
for i,row in enumerate(csvreader):
    if(i==1):
        list_of_elem = [h] + row
        list_of_elem[1] = ' '+list_of_elem[1]
file.close()

header = ['h', ' L2', ' H1']

# Open file in append mode
with open(write_file_name, 'a+', newline='') as write_obj:
    csv_writer = csv.writer(write_obj)

    if(nb_arg>1 and sys.argv[1]=='1'):
        csv_writer.writerow(header)
    csv_writer.writerow(list_of_elem)