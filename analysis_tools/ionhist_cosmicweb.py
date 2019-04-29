import sys
import os
import numpy as np
from numpy.random import random

from grid import *
import read_parameterfile as rp
import read_fields as rf
    
inifile = sys.argv[1]
inputIsDouble = np.int32(sys.argv[2])
threshold = np.float32(sys.argv[3])
n = np.int32(sys.argv[4])
outputfile = sys.argv[5]

lines = rp.read_inifile(inifile)

redshiftfile = rp.identify_string(lines, rp.redshiftfile_str, rp.splitting_str) #sys.argv[4]

simulationtype = rp.identify_string(lines, rp.simulationtype_str, rp.splitting_str)
if(simulationtype == "EVOLVE_BY_SNAPSHOT"):
  snapshotstart = rp.identify_int(lines, rp.snapshotstart_str, rp.splitting_str)
else:
  snapshotstart = 0
  
ionfile = rp.identify_string(lines, rp.ionfile_str, rp.splitting_str) #sys.argv[1]
densfile = rp.identify_string(lines, rp.densfile_str, rp.splitting_str)
double_precision = rp.identify_int(lines, rp.double_precision_str, rp.splitting_str)
isPadded = rp.identify_int(lines, rp.padded_str, rp.splitting_str) #np.int32(sys.argv[3])
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
print "Computing vsfk ionization history"
print "---------------------------------"


outputfile = outputfile + '_th' + str(threshold) + '_res' + str(n) + '.dat'

f = open(outputfile,'w')

cosmic = None
counter = snapshotstart
for i in range(len(redshift)-1):
    z = redshift[i+1]
    
    if(i + snapshotstart < 10):
        infile = ionfile + '_0' + str(i + snapshotstart)
    else:
        infile = ionfile + '_' + str(i + snapshotstart)
    if(os.path.isfile(infile) == False):
        continue
    
    #ionization field
    ion = rf.read_ion(infile, isPadded, inputIsDouble, gridsize)
    ion.shape = (gridsize, gridsize, gridsize)
    mean_ionfrac = np.mean(ion, dtype=np.float64)
    print "z =", z, "\tXHII =", mean_ionfrac

    if(snap[i] != 0):
        if(counter < 10):
            dinfile = densfile + '_00' + str(counter)
        else:
            dinfile = densfile + '_0' + str(counter)
        counter = counter + 1
        if(os.path.isfile(dinfile) == False):
            print "!!! density field was not updated"
            continue
        
        dens = rf.read_dens(dinfile, isPadded, double_precision, gridsize)
        dens.shape = (gridsize, gridsize, gridsize)
        gridsize_n = gridsize/n
        newdens = dens.reshape((gridsize/n,n,gridsize/n,n,gridsize/n,n)).mean(axis=1).mean(axis=2).mean(axis=3).T
        
        #fourier transform of density field
        modesDens = ifftn(newdens)

        k, kmag, inv_k2 = make_k_values(boxsize,gridsize_n)

        modesPhi = -modesDens*inv_k2

        kz = np.tile(k,(len(k),len(k),1))
        ky = np.swapaxes(kz,1,2)
        kx = np.swapaxes(kz,0,2)

        Dxy = fftn(-kx*ky*modesPhi).real
        Dxz = fftn(-kx*kz*modesPhi).real
        Dyz = fftn(-ky*ky*modesPhi).real
        Dxx = fftn(-kx*kx*modesPhi).real
        Dyy = fftn(-ky*ky*modesPhi).real
        Dzz = fftn(-kz*kz*modesPhi).real

        ddPhi = np.array([[Dxx,Dxy,Dxz],[Dxy,Dyy,Dyz],[Dxz,Dyz,Dzz]])
        eigenval = np.linalg.eigh(ddPhi.T)[0]

        zeros_array = np.zeros((len(k),len(k),len(k),3))
        ones_array = np.ones((len(k),len(k),len(k),3))

        zeros_array_r = np.zeros((len(k),len(k),len(k)))
        ones_array_r = np.ones((len(k),len(k),len(k)))

        poseigenval = np.where(eigenval>threshold,ones_array,zeros_array)
        cosmic_n = np.add.reduce(poseigenval,3)

        cosmic = cosmic_n
        
    if(isinstance(cosmic, np.ndarray)):
        newion = (dens*ion).reshape((gridsize/n,n,gridsize/n,n,gridsize/n,n)).mean(axis=1).mean(axis=2).mean(axis=3).T / newdens

        void_ionfrac = np.sum(newion.T*np.where(cosmic==0.,ones_array_r,zeros_array_r),dtype=np.float64)/np.sum(np.where(cosmic==0.,ones_array_r,zeros_array_r),dtype=np.float64)
        sheet_ionfrac = np.sum(newion.T*np.where(cosmic==1.,ones_array_r,zeros_array_r),dtype=np.float64)/np.sum(np.where(cosmic==1.,ones_array_r,zeros_array_r),dtype=np.float64)
        filament_ionfrac = np.sum(newion.T*np.where(cosmic==2.,ones_array_r,zeros_array_r),dtype=np.float64)/np.sum(np.where(cosmic==2.,ones_array_r,zeros_array_r),dtype=np.float64)
        knot_ionfrac = np.sum(newion.T*np.where(cosmic==3.,ones_array_r,zeros_array_r),dtype=np.float64)/np.sum(np.where(cosmic==3.,ones_array_r,zeros_array_r),dtype=np.float64)

        void_ionfrac_m = np.sum(newion.T*newdens.T*np.where(cosmic==0.,ones_array_r,zeros_array_r),dtype=np.float64)/np.sum(newdens.T*np.where(cosmic==0.,ones_array_r,zeros_array_r),dtype=np.float64)
        sheet_ionfrac_m = np.sum(newion.T*newdens.T*np.where(cosmic==1.,ones_array_r,zeros_array_r),dtype=np.float64)/np.sum(newdens.T*np.where(cosmic==1.,ones_array_r,zeros_array_r),dtype=np.float64)
        filament_ionfrac_m = np.sum(newion.T*newdens.T*np.where(cosmic==2.,ones_array_r,zeros_array_r),dtype=np.float64)/np.sum(newdens.T*np.where(cosmic==2.,ones_array_r,zeros_array_r),dtype=np.float64)
        knot_ionfrac_m = np.sum(newion.T*newdens.T*np.where(cosmic==3.,ones_array_r,zeros_array_r),dtype=np.float64)/np.sum(newdens.T*np.where(cosmic==3.,ones_array_r,zeros_array_r),dtype=np.float64)

        line = str(z)+'\t'+str(mean_ionfrac)+'\t'+str(void_ionfrac)+'\t'+str(sheet_ionfrac)+'\t'+str(filament_ionfrac)+'\t'+str(knot_ionfrac)+'\t'+str(void_ionfrac_m)+'\t'+str(sheet_ionfrac_m)+'\t'+str(filament_ionfrac_m)+'\t'+str(knot_ionfrac_m)+'\n'
        f.write(line)
    else:
        print "There was no density field to compute the cosmic web!"
f.close()
