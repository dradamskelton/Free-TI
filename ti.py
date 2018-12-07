import subprocess
import fileinput
import shutil
import numpy as np
import MDAnalysis
from MDAnalysis.lib import distances
import numpy.linalg

#dist=np.zeros(1000)
count =0 
kb = 306.432
int_free =0
shutil.copyfile('again.data','type2.data')

free = open("free.dat", "w")
freeint = open("final_free.dat", "w")

for i in range(10, 11, 1):
    int_pr = 0
    count =0
    dusum = 0
    file = open(str(i)+"-dist.dat", "w")
    file1 = open(str(i)+"-du.dat", "w")
#    file.write(
    shutil.copyfile('again.data',str(i)+'.data')
    shutil.copyfile('sim.in',str(i)+".in")
    for line in fileinput.input(str(i)+".in", inplace=True):

        print line.replace("again.data",str(i)+'.data' ),
    
    new= (float(i)/10.0)*kb
    
    for line in fileinput.input(str(i)+'.data', inplace=True):

        print line.replace("ch", str(new)),
        
    subprocess.call(["mpirun", "-np", "4", "lammps-daily", "-in", str(i)+".in"])
    shutil.copyfile('log.lammps',str(i)+".lammps")
    shutil.copyfile('dump.dcd',str(i)+".dcd")
    
    u = MDAnalysis.Universe('chargeFFat.psf','dump.dcd')
    num =  len(u.trajectory)
    dist=np.zeros(num+1)
    du=np.zeros(num+1)
    sel = u.select_atoms("bynum 751 ")
    sel1 = u.select_atoms("bynum 732")
    for ts in u.trajectory:
        count = count+1
        dist[count]=distances.distance_array(sel.positions, sel1.positions, box = [21.6209, 21.6209, 21.6209, 90, 90, 90])
       # dist_diff = 
        du[count] = 0.5*kb*((dist[count]-1.5080)**2)
#        print du[count], dist[count]-1.5080
        dusum = du[count] + dusum
        
    hist, bin_edges = np.histogram(dist, bins=1000, range=(0,10.0))
    hi_du, bin_du = np.histogram(du, bins=1000, range=(0,100))
    
    for j in range(1, 1000, 1):
#        print bin_edges[i], hist[i]
        #print dist[0,0]
#        hist[j] = float(hist[j])/float(num)
        bin = 10.0/1000.0
        hi =  float(hist[j])/(float(num)*bin)
        hidu = float(hi_du[j])/(float(num)*bin)
        int_pr = (hidu *(bin*10)) + int_pr
#        bin = 10.0/1000.0
#        file.write(str(bin_edges[j])+" "+str(hist[j])+"\n")
        file.write(str(bin_edges[j])+" "+str(hi)+"\n")
        file1.write(str(bin_du[j])+" "+str(hidu)+"\n")
    avdusum = dusum/float(num)
    lam = float(i)/10.0
    free.write(str(lam)+" "+str(avdusum)+" "+str(int_pr)+"\n")
    int_free = (avdusum*0.1)+int_free

freeint.write(str(int_free))  
