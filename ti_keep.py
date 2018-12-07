#from subprocess import Popen, PIPE, STDOUT
import subprocess
#import massedit
import fileinput
import shutil
import numpy as np
import MDAnalysis
from MDAnalysis.lib import distances
import numpy.linalg

shutil.copyfile('again.data','type2.data')


for i in range(0, 10, 1):
    
    shutil.copyfile('again.data',str(i)+'.data')
    shutil.copyfile('sim.in',str(i)+".in")
    for line in fileinput.input(str(i)+".in", inplace=True):

        print line.replace("again.data",str(i)+'.data' ),
    
    new= (float(i)/10.0)*306.432
    
    for line in fileinput.input(str(i)+'.data', inplace=True):

        print line.replace("ch", str(new)),
        
    subprocess.call(["mpirun", "-np", "4", "lammps-daily", "-in", str(i)+".in"])
    shutil.copyfile('log.lammps',str(i)+".lammps")
    shutil.copyfile('dump.dcd',str(i)+".dcd")
    
    u = MDAnalysis.Universe('chargeFFat.psf','dump.dcd')
    sel = u.select_atoms("bynum 751 ")
    sel1 = u.select_atoms("bynum 732")
    for ts in u.trajectory:
        dist=distances.distance_array(sel.positions, sel1.positions, box = [21.6209, 21.6209, 21.6209, 90, 90, 90])
        print dist[0,0]
        
    
#cmd = 'ls'
#p = Popen(cmd, shell=True, stdin=PIPE, stdout="file", stderr=STDOUT, close_fds=True)
#output = p.stdout.read()

#print output
#subprocess.call( ['ls'] stdout = open( 'logfile.log', 'w') )

#for i in range(0, 10, 1):
#outfile = open('test','w') #same with "w" or "a" as opening mode
#    outfile.write('Hello')
#    subprocess.Popen("sed", "{s/ch/"+str(i)+"/g}", "type.data",stdout=outfile)
#command = "sed", "{s/ch/10/g}", "type.data"
#subprocess.Popen(command.communicate)





#for i in range(0, 10, 1):


#subprocess.call(["lammps-daily", "-in", "sim.in"])
#    filenames = ['type.data']
#    cha = str(i)
#    subprocess.call(["sed", "{s/ch/"+str(i)+"/g}", "type.data"] stdout = open('logfile.log', 'w'))
#    massedit.edit_files(filenames, ["re.sub('ch', str(i), line)"])

#subprocess.call(["lammps-daily", "-in", "sim.in"])
#    subprocess.call(["lammps-daily", "-in", "sim.in"])
#    subprocess.call(["lammps-daily -in sim.in"])
#    subprocess.call(["cp", "log.lammps", str(i)+".log"])
 #   subprocess.call(["cp", "dump.xyz", str(i)+".xyz"])
    
    

