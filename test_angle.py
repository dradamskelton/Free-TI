import subprocess
import fileinput
import shutil
import numpy as np
import MDAnalysis
from MDAnalysis.lib import distances
import numpy.linalg
import scipy.integrate as integrate
from scipy.integrate import quad
from numpy.linalg import norm

numbin = 1000
file = open("test-dist.dat", "w")
count =0 

u = MDAnalysis.Universe('chargeFFat.psf','9.dcd')
num =  len(u.trajectory)
dist=np.zeros(num+1)


#A = u.select_atoms("bynum 732").center_of_geometry()
#B = u.select_atoms("bynum 751 ").center_of_geometry()
#C = u.select_atoms("bynum 749 ").center_of_geometry()


#BA = A - B
#BC = C - B
#theta = np.arccos(np.dot(BA, BC)/(norm(BA)*norm(BC)))
#print np.rad2deg(theta)

cell = [21.6209, 21.6209, 21.6209]
#two = [2.0, 2.0, 2.0]
#cell2 = cell*two
for ts in u.trajectory:
    A = u.select_atoms("bynum 732").center_of_geometry()
    B = u.select_atoms("bynum 751 ").center_of_geometry()
    C = u.select_atoms("bynum 749 ").center_of_geometry()
    BA = A - B
    BC = C - B
#    print BA
    for i in  range(0,3):
        if BA[i] > cell[i]/2.0:
            BA[i] = BA[i]-cell[i]
    for i in  range(0,3):
        if BA[i] < -cell[i]/2.0:
            BA[i] = BA[i]+cell[i]       
    for i in  range(0,3):
        if BC[i] > cell[i]/2.0:
            BC[i] = BC[i]-cell[i]
    for i in  range(0,3):
        if BC[i] < -cell[i]/2.0:
            BC[i] = BC[i]+cell[i]
              
    theta = np.arccos(np.dot(BA, BC)/(norm(BA)*norm(BC)))
    file.write (str(np.rad2deg(theta))+"\n")



    count = count+1
    dist[count]=distances.distance_array(sel.positions, sel1.positions, box = [21.6209, 21.6209, 21.6209, 90, 90, 90])
    hist, bin_edges = np.histogram(dist, bins=numbin, range=(1,5), density=True)

#for j in range(1, numbin, 1):
#    file.write(str(bin_edges[j])+" "+str(hist[j])+"\n")
