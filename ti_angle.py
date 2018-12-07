import subprocess
import fileinput
import shutil
import numpy as np
import MDAnalysis
from MDAnalysis.lib import distances
import numpy.linalg
import scipy.integrate as integrate
from scipy.integrate import quad

### initiating constants and opening files 

#dist=np.zeros(1000)
count =0 
kb = 306.432
ka1 =  57.788
ka2 = 61.243 
int_free =0
free = open("free.dat", "w")
freeint = open("final_free.dat", "w")
numbin = 1000

shutil.copyfile('again.data','type2.data')

#### loop of lambda values ###

for i in range(1, 11, 1):
    int_pr = 0
    count =0
    dusum = 0
    file = open(str(i)+"-dist.dat", "w")
    file1 = open(str(i)+"-du.dat", "w")

    shutil.copyfile('again.data',str(i)+'.data')
    shutil.copyfile('sim.in',str(i)+".in")
    for line in fileinput.input(str(i)+".in", inplace=True):

        print line.replace("again.data",str(i)+'.data' ),
    
    new= (float(i)/10.0)*kb
    ang1 = (float(i)/10.0)*ka1
    ang2 = (float(i)/10.0)*ka2
    
####  changing the bond and angle parameters
    
    for line in fileinput.input(str(i)+'.data', inplace=True):

        print line.replace("ch", str(new)),

    for line in fileinput.input(str(i)+'.data', inplace=True):
        print line.replace("ang1", str(ang1)),
        
    for line in fileinput.input(str(i)+'.data', inplace=True):
        print line.replace("ang2", str(ang2)),    

    ######  running lammps  ###
 
    subprocess.call(["mpirun", "-np", "4", "lammps-daily", "-in", str(i)+".in"])
    shutil.copyfile('log.lammps',str(i)+".lammps")
    shutil.copyfile('dump.dcd',str(i)+".dcd")

    ## reading in  trajectory and setting up selections  
    
    u = MDAnalysis.Universe('chargeFFat.psf','dump.dcd')
    num =  len(u.trajectory)
#    sel = u.select_atoms("bynum 751 ")
#    sel1 = u.select_atoms("bynum 732")
    A = u.select_atoms("bynum 732").center_of_geometry()
    B = u.select_atoms("bynum 751 ").center_of_geometry()
    C = u.select_atoms("bynum 749 ").center_of_geometry()
    BA = A - B
    BC = C - B
    theta = np.arccos(np.dot(BA, BC)/(norm(BA)*norm(BC)))
#    print np.rad2deg(theta)
    
    dist=np.zeros(num+1)
    du=np.zeros(num+1)
    
    ### loop  for calculating  distance information from trajectory 
    
    for ts in u.trajectory:
        count = count+1
        A = u.select_atoms("bynum 732").center_of_geometry()
        B = u.select_atoms("bynum 751 ").center_of_geometry()
        C = u.select_atoms("bynum 749 ").center_of_geometry()

        
        dist[count]=distances.distance_array(sel.positions, sel1.positions, box = [21.6209, 21.6209, 21.6209, 90, 90, 90])

        du[count] = 0.5*kb*((dist[count]-1.5080)**2)

        dusum = du[count] + dusum
        
    hist, bin_edges = np.histogram(dist, bins=numbin, range=(1,5), density=True)
#    hi_du, bin_du = np.histogram(du, bins=1000, range=(0,100))

#### loop for manipulating histogram  ### 
    
    for j in range(1, numbin, 1):

        bin = 0.1/float(numbin)
        mid = (bin_edges[j]+ bin_edges[j+1])/2.0
        hi =  float(hist[j])*(0.5*kb*((mid-1.5080)**2))
#        hidu = float(hi_du[j])/(float(num)*bin)
        int_pr = (hi *(bin)) + int_pr

        file.write(str(bin_edges[j])+" "+str(hist[j])+"\n")
        file1.write(str(bin_edges[j])+" "+str(hi)+"\n")
    avdusum = dusum/float(num)
    lam = float(i)/10.0
    free.write(str(lam)+" "+str(avdusum)+"\n")
#    int_free = (avdusum*0.1)+int_free

#freeint.write(str(int_pr))  
