def create_circle(h,num):
    f = open("circle/circle_dirichlet_"+str(num)+".geo", "w")

    f.write('h = '+str(h)+';\n')
    f.write('\n')

    f.write('Point(1) = {0, 0, 0, h};\n')
    f.write('Point(2) = {1, 0, 0, h};\n')
    f.write('Point(3) = {-1, 0, 0, h};\n')
    f.write('Point(4) = {0, 1, 0, h};\n')
    f.write('Point(5) = {0, -1, 0, h};\n')
    f.write('\n')

    f.write('Circle(1) = {2, 1, 4};\n')
    f.write('Circle(2) = {4, 1, 3};\n')
    f.write('Circle(3) = {3, 1, 5};\n')
    f.write('Circle(5) = {5, 1, 2};\n')
    f.write('\n')

    f.write('Line Loop(6) = {1, 2, 3, 5};\n')
    f.write('Plane Surface(7) = {6};\n')
    f.write('\n')

    f.write('Physical Line("Dirichlet") = {1, 2, 3, 5};\n')
    f.write('Physical Surface(8) = {7};\n')
    f.write('\n')
    
    f.close()

tab_h = [0.4]
for i in range(4):
    tab_h.append(tab_h[-1]/2)

cercles = []
for i,h in enumerate(tab_h):
    print("Création du disque "+str(i+1)+" avec h = "+str(h))
    create_circle(h,i+1)