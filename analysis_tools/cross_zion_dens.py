import sys
import os
import numpy as np

from grid import *
from cross_corr import *
import read_parameterfile as rp
import read_fields as rf

km_cm = 1.e5
Mpc_cm = 3.086e24
G = 6.67408e-8
mp_g = 1.673e-24


def round_down(num):
    if num < 0:
        return -np.ceil(abs(num))
    else:
        return np.int32(num)

inifile = sys.argv[1]
inputIsDouble = np.int32(sys.argv[2])
specie = np.int32(sys.argv[3])
outputfile = sys.argv[4]

lines = rp.read_inifile(inifile)

redshiftfile = rp.identify_string(lines, rp.redshiftfile_str, rp.splitting_str)

solve_he = rp.identify_int(lines, rp.solve_he_str, rp.splitting_str)

if(specie > 0 and solve_he == 0):
    sys.exit()

if(specie == 0):
    ionfile = rp.identify_string(lines, rp.ionfile_str, rp.splitting_str)
elif(specie == 1):
    ionfile = rp.identify_string(lines, rp.HeIIionfile_str, rp.splitting_str)
elif(specie == 2):
    ionfile = rp.identify_string(lines, rp.HeIIIionfile_str, rp.splitting_str)
else:
    sys.exit()
    
densfile = rp.identify_string(lines, rp.densfile_str, rp.splitting_str)
double_precision = rp.identify_int(lines, rp.double_precision_str, rp.splitting_str)
isPadded = rp.identify_int(lines, rp.padded_str, rp.splitting_str)
isPadded_factor = isPadded**(1./3.)
if(isPadded != 0):
    gridsize = np.int32(rp.identify_int(lines, rp.gridsize_str, rp.splitting_str)/isPadded_factor)
    boxsize = rp.identify_float(lines, rp.boxsize_str, rp.splitting_str)/isPadded_factor
else:
    gridsize = rp.identify_int(lines, rp.gridsize_str, rp.splitting_str)
    boxsize = rp.identify_float(lines, rp.boxsize_str, rp.splitting_str)

Omegab = rp.identify_float(lines, rp.omega_b_str, rp.splitting_str)
Omegam = rp.identify_float(lines, rp.omega_m_str, rp.splitting_str)
h = rp.identify_float(lines, rp.h_str, rp.splitting_str)

redshift, snap = np.loadtxt(redshiftfile, unpack='True', skiprows=0, usecols=(0,1))
snap = np.int32(snap)

print "\n-----------------------------------------------------------"
print "Computing cross correlation between zion and density fields"
print "-----------------------------------------------------------"

#----------------------------------------------
#----------------------------------------------
Hubble0 = h * 100. * km_cm / Mpc_cm
numdensity = 3.* Hubble0**2 * Omegab/(8. * np.pi * G * mp_g)

if(specie == 0):
    fp = open(outputfile + 'zion.dat', "rb")
    zion = np.fromfile(fp, count=gridsize*gridsize*gridsize, dtype=np.float64)
    fp.close()

elif(specie == 1):
    fp = open(outputfile + 'zionHeI.dat', "rb")
    zion = np.fromfile(fp, count=gridsize*gridsize*gridsize, dtype=np.float64)
    fp.close()

elif(specie == 2):
    fp = open(outputfile + 'zionHeII.dat', "rb")
    zion = np.fromfile(fp, count=gridsize*gridsize*gridsize, dtype=np.float64)
    fp.close()

else:
    sys.exit()

zion = np.reshape(zion, (gridsize, gridsize, gridsize))

counter = 0
for i in range(len(redshift)-1):    
    #-----------------------
    # density field
    #-----------------------
    if(snap[i] != 0):
        if(counter <10):
            dinfile = densfile + '_00' + str(counter)
        else:
            dinfile = densfile + '_0' + str(counter)
        counter = counter + 1
        if(os.path.isfile(dinfile) == False):
            print "!!! density field was not updated"
            continue
        
        dens = rf.read_dens(dinfile, isPadded, double_precision, gridsize)
    
        #----------------------------------------------
        # COMPUTING CROSS CORRELATION
        #----------------------------------------------

        compute_cross_corr(1.+zion, dens, boxsize, outputfile, 'zion', 'dens%02d'%(counter-1))

