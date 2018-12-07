import numpy as np
import MDAnalysis
#from MDAnalysis.lib.distances import distance_array, self_distance_array
#from MDAnalysis.lib.c_distances import contact_matrix_no_pbc, contact_matrix_pbc
#from MDAnalysis.lib.NeighborSearch import AtomNeighborSearch
#from MDAnalysis.lib.distances import calc_bonds
from MDAnalysis.lib import distances
import numpy.linalg
import pandas as pd


dist=np.zeros(1000)
count =0 
u = MDAnalysis.Universe('chargeFFat.psf','4.dcd')

sel = u.select_atoms("bynum 751 ")
sel1 = u.select_atoms("bynum 732")
#protein1 = u.select_atoms("type CR")


#sel.write("test22.xyz")
#print (sel.center(sel.masses, compound='residues'))
#sel.distances
#print (distances.distance_array(sel.positions, sel1.positions, box = [21.6209, 21.6209, 21.6209, 90, 90, 90]))

#print (distances.distance_array(sel.positions, sel1.positions))
for ts in u.trajectory:
    count = count+1
    dist[count]=distances.distance_array(sel.positions, sel1.positions, box = [21.6209, 21.6209, 21.6209, 90, 90, 90])
    
#    print count, dist[0,0]

#print dist

hist, bin_edges = np.histogram(dist, bins=1000, range=(1.2,3.0))
print hist

print bin_edges
#print hist.size, bin_edges.size


#    print dist[0,0]*0.5*300
#    print (sel.positions[0,0])
#    print (sel.types[0], sel.positions[0,0], sel.positions[0,1], sel.positions[0,2])
#    print (sel.types[0], sel.positions[0])
#    print (sel.types[1], sel.positions[1,0], sel.positions[1,1], sel.positions[1,2])
