import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

rhit = []
zhit = []
phit = []
f = open('int_pts.out','r')
for line in f:
#    line = line.strip()
    data = line.strip().split()
    rhit.append(float(data[0]))
    zhit.append(float(data[1]))
    phit.append(float(data[2])*np.pi*2.)       
f.close()

xhit = rhit*np.cos(phit)
yhit = rhit*np.sin(phit)


rgeo = []
zgeo = []
pgeo = []
f = open('allparts.out','r')
line = f.readline()
line = f.readline()
#data = line.strip().split()
#print(data)
for line in f:
    line = line.strip()
    data = line.split()
    
f.close()

#fig = plt.figure()
#ax = Axes3D(fig)
#ax.scatter(xhit,yhit,zhit,marker=".",color="r")
#plt.show()


    
