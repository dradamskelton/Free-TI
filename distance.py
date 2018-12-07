import numpy as np
import MDAnalysis
from MDAnalysis.lib import distances
import numpy.linalg
import pandas as pd
import matplotlib.pyplot as plt

dist=np.zeros(1000)
count =0 
u = MDAnalysis.Universe('chargeFFat.psf','4.dcd')
num =  len(u.trajectory)
kb = 306.432

dist=np.zeros(num+1)
du=np.zeros(num+1)


sel = u.select_atoms("bynum 751 ")
sel1 = u.select_atoms("bynum 732")




for ts in u.trajectory:
    count = count+1
    dist[count]=distances.distance_array(sel.positions, sel1.positions, box = [21.6209, 21.6209, 21.6209, 90, 90, 90])
    du[count] = 0.5*kb*((dist[count])**2)
    print du[count], dist[count]

hist, bin_edges = np.histogram(dist, bins=1000, range=(1.2,3.0))
hi_du, bin_du = np.histogram(du, bins=1000, range=(0,500))
#hist, bin_edges = np.histogram(dist, bins=1000, range=(1.2,3.0))

#for i in range(1, 1000, 1):
#    print bin_du[i], hi_du[i]

#print hist

#print bin_edges

#plt.plot(bin_edges,hist,"r--")
#plt.show()
