import subprocess
import fileinput
import shutil
import numpy as np
import MDAnalysis
from MDAnalysis.lib import distances
import numpy.linalg
import scipy.integrate as integrate
from scipy.integrate import quad

numbin = 1000
file = open("test-dist.dat", "w")
count =0 

u = MDAnalysis.Universe('chargeFFat.psf','5.dcd')
num =  len(u.trajectory)
dist=np.zeros(num+1)
sel = u.select_atoms("bynum 751 ")
sel1 = u.select_atoms("bynum 732")

for ts in u.trajectory:
    count = count+1
    dist[count]=distances.distance_array(sel.positions, sel1.positions, box = [21.6209, 21.6209, 21.6209, 90, 90, 90])
    hist, bin_edges = np.histogram(dist, bins=numbin, range=(1,5), density=True)

for j in range(1, numbin, 1):
    file.write(str(bin_edges[j])+" "+str(hist[j])+"\n")
