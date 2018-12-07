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

# reading in the angles that need to be broken

that=np.zeros((100,100))

file = open("angles", "r")

lines= file.readlines()
num = len(lines)

file.close()
for i in range(0,num):
#    print lines[i]
    this = lines[i].split()
    for j in range(0,4):
        that[i,j] = this[j]


### initiating constants and opening files 

#dist=np.zeros(1000)
count =0 
kb = 306.432
ka1 =  57.788
ka2 = 61.243
ka3 = 40.517
int_free =0
#free = open("free.dat", "w")
#freeint = open("final_free.dat", "w")
numbin = 1000
delr = 360.0/1000.0
sumtheta=np.zeros((1000,1000))
dusum=np.zeros(100)
file2=np.zeros(100)

shutil.copyfile('again.data','type2.data')
cell = [21.6209, 21.6209, 21.6209]

for p in range(0,6):
    file2 = open(str(p)+"-lam.dat", "w")
    file2.close()

#### loop of lambda values ###

for i in range(1, 11, 1):
    dusum=np.zeros(100)
    sumtheta=np.zeros((1000,1000))
    int_pr = 0
    count =0
#    dusum = 0
#    file = open(str(i)+"-dist.dat", "w")
#    file1 = open(str(i)+"-du.dat", "w")

    shutil.copyfile('again.data',str(i)+'.data')
    shutil.copyfile('sim.in',str(i)+".in")
    for line in fileinput.input(str(i)+".in", inplace=True):

        print line.replace("again.data",str(i)+'.data' ),
#    new = kb
    new= (float(i)/10.0)*kb
    ang1 = (float(i)/10.0)*ka1
#    ang1 = ka1
#    ang3 = ka3
    ang2 = (float(i)/10.0)*ka2
    ang3 = (float(i)/10.0)*ka3

    
####  changing the bond and angle parameters
    
    for line in fileinput.input(str(i)+'.data', inplace=True):

        print line.replace("ch", str(new)),

    for line in fileinput.input(str(i)+'.data', inplace=True):
        print line.replace("ang1", str(ang1)),
        
    for line in fileinput.input(str(i)+'.data', inplace=True):
        print line.replace("ang2", str(ang2)),

    for line in fileinput.input(str(i)+'.data', inplace=True):
        print line.replace("ang3", str(ang3)),    

    ######  running lammps  ###
 
    subprocess.call(["mpirun", "-np", "4", "lammps-daily", "-in", str(i)+".in"])
    shutil.copyfile('log.lammps',str(i)+".lammps")
    shutil.copyfile('dump.dcd',str(i)+".dcd")

    ## reading in  trajectory and setting up selections  
    
    u = MDAnalysis.Universe('chargeFFat.psf','dump.dcd')
    num =  len(u.trajectory)
#    sel = u.select_atoms("bynum 751 ")
#    sel1 = u.select_atoms("bynum 732")
 
    
    dist=np.zeros(num+1)
#    du=np.zeros(num+1)
    deg=np.zeros(num+1)
    theta=np.zeros((num+1,6))
    ### loop  for calculating  distance information from trajectory 
    
    for ts in u.trajectory:
        count = count+1

        ### loop for doing the different angles 
    
        for p in range(0,6):
            A = u.select_atoms("bynum "+str(int(that[p, 0]))).center_of_geometry()
            B = u.select_atoms("bynum "+str(int(that[p, 1]))).center_of_geometry()
            C = u.select_atoms("bynum "+str(int(that[p, 2]))).center_of_geometry()

            BA = A - B
            BC = C - B
            for h in  range(0,3):
                if BA[h] > cell[h]/2.0:
                    BA[h] = BA[h]-cell[h]
            for h in  range(0,3):
                if BA[h] < -cell[h]/2.0:
                    BA[h] = BA[h]+cell[h]       
            for h in  range(0,3):
                if BC[h] > cell[h]/2.0:
                    BC[h] = BC[h]-cell[h]
            for h in  range(0,3):
                if BC[h] < -cell[h]/2.0:
                    BC[h] = BC[h]+cell[h]
            theta[count,p] = np.rad2deg(np.arccos(np.dot(BA, BC)/(norm(BA)*norm(BC))))
            thetarad = np.arccos(np.dot(BA, BC)/(norm(BA)*norm(BC)))
  #          theta[count] = np.arccos(np.dot(BA, BC)/(norm(BA)*norm(BC)))
            bin_no = int(((theta[count,p])/delr) )
            sumtheta[bin_no,p] = sumtheta[bin_no,p] + 1
            if that[p, 3] == 1:
                du = 0.5*ka1*((thetarad - np.deg2rad(120.4190))**2)
            if that[p, 3] == 2:
                du = 0.5*ka2*((thetarad - np.deg2rad(109.6080))**2)
            if that[p, 3] == 3:
                du = 0.5*ka3*((thetarad - np.deg2rad(120.5710))**2)
                
 #           deg[count] = np.rad2deg(theta[count])
#    print np.rad2deg(theta)

#        dist[count]=distances.distance_array(sel.positions, sel1.positions, box = [21.6209, 21.6209, 21.6209, 90, 90, 90])

#        du[count] = 0.5*ka2*((theta[count]- np.deg2rad(109.6080))**2)

            dusum[p] = du + dusum[p]
        
    hist, bin_edges = np.histogram(deg, bins=numbin, range=(0,180), density=True)
#    hi_du, bin_du = np.histogram(du, bins=1000, range=(0,100))

#### loop for managing the output ##
    lam = float(i)/10.0
    for p in range(0,6):
        file1 = open(str(p)+"-"+str(i)+"test.dat", "w")
        file2 = open(str(p)+"-lam.dat", "a")

    #### loop for manipulating histogram   ### 
        for j in range(1, numbin, 1):
            file1.write(str(delr*float(j))+" "+str(sumtheta[j,p])+"\n")
        file2.write( str(lam)+" "+str(dusum[p]/float(num))+"\n")
        file2.close()
 #       bin = 0.1/float(numbin)
 #       mid = (bin_edges[j]+ bin_edges[j+1])/2.0
 #       hi =  float(hist[j])*(0.5*kb*((mid-1.5080)**2))
#        hidu = float(hi_du[j])/(float(num)*bin)
 ##       int_pr = (hi *(bin)) + int_pr
 
#        file.write(str(bin_edges[j])+" "+str(hist[j])+"\n")
 #       file1.write(str(bin_edges[j])+" "+str(hi)+"\n")
#    avdusum = dusum/float(num)
#    lam = float(i)/10.0
#    free.write(str(lam)+" "+str(avdusum)+"\n")
#    int_free = (avdusum*0.1)+int_free

#freeint.write(str(int_pr))  
