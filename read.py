import numpy as np
import MDAnalysis
from numpy.linalg import norm
print "now"
cell = [21.6209, 21.6209, 21.6209]
delr = 360.0/1000.0
that=np.zeros((100,100))
this=np.zeros(100)
theta=np.zeros((100,100))
file1=np.zeros(100)
#hist=np.zeros(4)
#bin_edges=np.zeros(4)
count = 0
sumtheta=np.zeros((1000,1000))
### 

file = open("angletype", "r")
l= file.readlines()
num = len(l)

l = map(str.strip, l)
l = map(str.split, l)
l1 = sum(l, [])
print l1
numbin=1000

file.close()

a = np.genfromtxt("angletype", dtype=str)
 
print a
a.flatten(order="F")



#### loop to assign the atoms within the angles

#for i in range(0,num):
#    print lines[i]
#    this = lines[i].split()
#    print i, str(this[j])
#    print this[:]
#    for j in range(0,3):
#        print i, this[i]
#        that[i,j] = str(this[j])
#    print that[i, 0], str(that[i, 1]), str(that[i, 2]), str(that[i, 3])#


#u = MDAnalysis.Universe('chargeFFat.psf','0.xyz')

#for i in range(0,num):
#    for j in range(0,3):
#    print this[i] 
## loop over trajectory
#for ts in u.trajectory:
#    count = count+1
    
    ### loop over the different angles that will be broken

#    for p in range(0,6):    
#        A = u.select_atoms("bynum "+str(int(that[p, 0]))).center_of_geometry()    
#        B = u.select_atoms("bynum "+str(int(that[p, 1]))).center_of_geometry()
##        C = u.select_atoms("bynum "+str(int(that[p, 2]))).center_of_geometry()

#        BA = A - B
#        BC = C - B
#        for h in  range(0,3):
#            if BA[h] > cell[h]/2.0:
#                BA[h] = BA[h]-cell[h]
#        for h in  range(0,3):
#            if BA[h] < -cell[h]/2.0:
#                BA[h] = BA[h]+cell[h]       
#        for h in  range(0,3):
#            if BC[h] > cell[h]/2.0:
#                BC[h] = BC[h]-cell[h]
#        for h in  range(0,3):
#            if BC[h] < -cell[h]/2.0:
#                BC[h] = BC[h]+cell[h]
#        theta[count,p] = np.rad2deg(np.arccos(np.dot(BA, BC)/(norm(BA)*norm(BC))))
#        print theta[count,p]
