import sys
import os
import numpy as np
from numpy.random import random
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as m
import glob
import re

from grid import *  
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

simulationtype = rp.identify_string(lines, rp.simulationtype_str, rp.splitting_str)
if(simulationtype == "EVOLVE_BY_SNAPSHOT"):
  snapshotstart = rp.identify_int(lines, rp.snapshotstart_str, rp.splitting_str)
else:
  snapshotstart = 0
  
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

print "\n---------------------------------"
print "Computing zion and dension fields"
print "---------------------------------"

#----------------------------------------------
#----------------------------------------------
Hubble0 = h * 100. * km_cm / Mpc_cm
numdensity = 3.* Hubble0**2 * Omegab/(8. * np.pi * G * mp_g)

threshold = 0.5
zion = np.zeros(gridsize**3)
dension = np.zeros(gridsize**3)

meanIon_list = np.zeros(len(redshift))
meanIon_mass_list = np.zeros(len(redshift))

ionfile_table = sorted(glob.glob(ionfile + '_*'))
densfile_table = sorted(glob.glob(densfile + '*'))

regex = re.compile(r'\d+')

dens_table_min = np.int32(regex.findall(densfile_table[0])[-1])
dens_table_max = np.int32(regex.findall(densfile_table[-1])[-1])

ion_table_min = np.int32(regex.findall(ionfile_table[0])[-1])
ion_table_max = np.int32(regex.findall(ionfile_table[-1])[-1])

print ion_table_min, ion_table_max
print dens_table_min, dens_table_max

if (len(densfile_table) > len(ionfile_table)):
    densfile_table = densfile_table[(ion_table_min-dens_table_min):ion_table_max+1]

assert(len(ionfile_table)==len(redshift)-1)
assert(len(densfile_table)<=len(redshift)-1)

counter = snapshotstart
for i in range(len(redshift)-1):
    z = redshift[i+1]

    #-----------------------
    # ionization field
    #-----------------------
    if(i + snapshotstart < 10):
        infile = ionfile + '_0' + str(i + snapshotstart)
    else:
        infile = ionfile + '_' + str(i + snapshotstart)
        
    infile = ionfile_table[i]
    if(os.path.isfile(infile) == False):
        continue
    
    ion = rf.read_ion(infile, isPadded, inputIsDouble, gridsize)
    print "z =", z, "\tXHII =", np.mean(ion, dtype=np.float64), "\t", infile
    ion = np.ravel(ion)
    
    #-----------------------
    # density field
    #-----------------------
    if(snap[i] != 0):
        if(counter < 10):
            dinfile = densfile + '_00' + str(counter)
        else:
            dinfile = densfile + '_0' + str(counter)
        
        print counter
        dinfile = densfile_table[counter]
        print dinfile
        counter = counter + 1
        if(os.path.isfile(dinfile) == False):
            print "!!! density field was not updated"
            continue
        
    dens = rf.read_dens(dinfile, isPadded, double_precision, gridsize)
    dens = np.ravel(dens)
    
    #-----------------------
    # zion and dension
    #-----------------------
    indices = np.intersect1d(np.where(ion>threshold)[0], np.where(zion==0.)[0])
    zion[indices] = z
    dension[indices] = dens[indices]*(1.+z)**3*numdensity

    meanIon_list[i] = np.mean(ion, dtype=np.float64)
    meanIon_mass_list[i] = np.mean(ion*dens, dtype=np.float64) / np.mean(dens, dtype=np.float64)
     
if(specie == 0):
    indices = np.where(zion==0.)[0]
    zion[indices] = z
    dension[indices] = dens[indices]*(1.+z)**3*numdensity

print np.amin(zion), np.amax(zion)
print z

if(specie == 0):
    fp = open(outputfile + 'zion.dat', "wb")
    fp.write(bytearray(zion))
    fp.close()

    fp = open(outputfile + 'dension.dat', "wb")
    fp.write(bytearray(dension))
    fp.close()

    np.savetxt(outputfile + 'meanIon_z.dat', np.c_[redshift, meanIon_list, meanIon_mass_list])
elif(specie == 1):
    fp = open(outputfile + 'zionHeI.dat', "wb")
    fp.write(bytearray(zion))
    fp.close()

    fp = open(outputfile + 'densionHeI.dat', "wb")
    fp.write(bytearray(dension))
    fp.close()

    np.savetxt(outputfile + 'meanIonHeI_z.dat', np.c_[redshift, meanIon_list, meanIon_mass_list])
elif(specie == 2):
    fp = open(outputfile + 'zionHeII.dat', "wb")
    fp.write(bytearray(zion))
    fp.close()

    fp = open(outputfile + 'densionHeII.dat', "wb")
    fp.write(bytearray(dension))
    fp.close()

    np.savetxt(outputfile + 'meanIonHeII_z.dat', np.c_[redshift, meanIon_list, meanIon_mass_list])
else:
    sys.exit()

zion = np.reshape(zion, (gridsize, gridsize, gridsize))
dension = np.reshape(dension, (gridsize, gridsize, gridsize))
dens = np.reshape(dens, (gridsize, gridsize, gridsize))


#----------------------------------------------
# Plotting fields
#----------------------------------------------
cut_slice = 64

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.rcParams.update({'figure.subplot.hspace':0.0, 'figure.subplot.bottom':0.1, 'figure.subplot.top':0.98, 'figure.subplot.left':0.08, 'figure.subplot.right':0.88})
plt.rcParams["figure.figsize"] = (7,6)

#----------------------------------------------
# zion field
#----------------------------------------------
fig = plt.figure()
ax = plt.gca()
image = ax.imshow(zion[:,:,cut_slice:cut_slice+1].mean(axis=-1),origin='lower',interpolation='nearest',  extent=[0,boxsize,0,boxsize], cmap="afmhot_r")

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = plt.colorbar(image, cax=cax)
cb.set_label('redshift at ionization   z$_{ion}$')

ax.set_xlabel('x  [h$^{-1}$ Mpc]')
ax.set_ylabel('y  [h$^{-1}$ Mpc]')
plt.savefig(outputfile + 'field_zion.png',format='png', dpi=512, transparent=True)

#----------------------------------------------
# dension field
#----------------------------------------------
fig = plt.figure()
ax = plt.gca()
image = ax.imshow(np.log10(dension)[:,:,cut_slice:cut_slice+1].mean(axis=-1),origin='lower',interpolation='nearest',  extent=[0,boxsize,0,boxsize], cmap="afmhot_r")

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = plt.colorbar(image, cax=cax)
cb.set_label('number density at ionization   n$_{ion}$ [cm$^{-3}$]')

ax.set_xlabel('x  [h$^{-1}$ Mpc]')
ax.set_ylabel('y  [h$^{-1}$ Mpc]')
plt.savefig(outputfile + 'field_dension.png',format='png', dpi=512, transparent=True)

