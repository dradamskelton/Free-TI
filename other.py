import MDAnalysis
#from MDAnalysis.tests.datafiles import PSF, DCD   # test trajectory
import numpy.linalg

#u = MDAnalysis.Universe(PSF,DCD)  # always start with a Universe

u = MDAnalysis.Universe('chargeFFat.psf','dump.xyz')

# can access via segid (4AKE) and atom name
# we take the first atom named N and the last atom named C
nterm = u.select_atoms("bynum 751")
cterm = u.select_atoms("bynum 732")

bb = u.select_atoms('protein and backbone')  # a selection (AtomGroup)

for ts in u.trajectory:     # iterate through all frames
    r = cterm.position - nterm.position # end-to-end vector from atom positions
    d = numpy.linalg.norm(r)  # end-to-end distance
    rgyr = bb.radius_of_gyration()  # method of AtomGroup
    print("frame = {0}: d = {1} A, Rgyr = {2} A".format(
          ts.frame, d, rgyr))
    
