import sys
import os
import numpy as np
from numpy.random import random
import math as m
import matplotlib as m
m.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mt
from matplotlib import gridspec
from matplotlib import rc

from grid import *
import read_parameterfile as rp
import read_fields as rf
import observations as ob

threshold = 0.5

def compute_meanDens(infile, isPadded, inputIsDouble, gridsize):
    dens = rf.read_dens(infile, isPadded, inputIsDouble, gridsize)
    
    meanDens = np.mean(dens, dtype=np.float64)
    return meanDens
  
def compute_meanIon(infile, isPadded, inputIsDouble, gridsize):
    ion = rf.read_ion(infile, isPadded, inputIsDouble, gridsize)
    
    meanIon = np.mean(ion, dtype=np.float64)
    return meanIon

def compute_meanDensIon(infile, densfile, double_precision, isPadded, inputIsDouble, gridsize):
    ion = rf.read_ion(infile, isPadded, inputIsDouble, gridsize)
    dens = rf.read_dens(densfile, isPadded, double_precision, gridsize)

    #indexCells = np.where(ion > threshold)
    #meanDensInIon = np.mean(dens[indexCells], dtype=np.float64)/np.mean(dens, dtype=np.float64)
    meanDensInIon = np.mean(ion*dens, dtype=np.float64)/(np.mean(ion, dtype=np.float64) * np.mean(dens, dtype=np.float64))
    return meanDensInIon
  
def compute_meanDensNeutral(infile, densfile, double_precision, isPadded, inputIsDouble, gridsize):
    ion = rf.read_ion(infile, isPadded, inputIsDouble, gridsize)
    dens = rf.read_dens(densfile, isPadded, double_precision, gridsize)
    
    #indexCells = np.where(1. - ion > threshold)
    #meanDensInNeutral = np.mean(dens[indexCells], dtype=np.float64)/np.mean(dens, dtype=np.float64)
    meanDensInNeutral = np.mean((1.-ion)*dens, dtype=np.float64)/(np.mean(1.-ion, dtype=np.float64) * np.mean(dens, dtype=np.float64))
    return meanDensInNeutral

inifile = sys.argv[1]
inputIsDouble = np.int32(sys.argv[2])
outputfile = sys.argv[3]

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
    
redshift, snap = np.loadtxt(redshiftfile, unpack='True', skiprows=0, usecols=(0,1))

print "\n----------------------------"
print "Computing ionization history"
print "----------------------------"

hist_ion = np.zeros(len(redshift)-1)
hist_mass_ion = np.ones(len(redshift)-1)
hist_mass_neutral = np.ones(len(redshift)-1)

counter = snapshotstart
for i in range(len(redshift)-1):
    if(i + snapshotstart < 10):
        infile = ionfile + '_0' + str(i + snapshotstart)
    else:
        infile = ionfile + '_' + str(i + snapshotstart)
    
    print infile
    
    if(snap[i] != 0):
        if(counter <10):
            dinfile = densfile + '_00' + str(counter)
        else:
            dinfile = densfile + '_0' + str(counter)
        counter = counter + 1
            
    if(os.path.isfile(infile) == True):
        hist_ion[i] = compute_meanIon(infile, isPadded, inputIsDouble, gridsize)
    elif(i > 0):
        print "!!! XHII: replacing value at z=", redshift[i+1], "with previous value at", redshift[i]
        hist_ion[i] = hist_ion[i-1]
    else:
        print "!!! XHII: replacing value at z=", redshift[i+1], "with previous value 0"
        hist_ion[i] = 0.
        
    if(os.path.isfile(infile) == True and os.path.isfile(dinfile) == True):
        hist_mass_ion[i] = compute_meanDensIon(infile, dinfile, double_precision, isPadded, inputIsDouble, gridsize)
    elif(i > 0):
        print "!!! densHII: replacing value at z=", redshift[i+1], "with previous value at", redshift[i]
        hist_mass_ion[i] = hist_mass_ion[i-1] 
    else:
        print "!!! densHII: replacing value at z=", redshift[i+1], "with previous value 0"
        hist_mass_ion[i] = 1.
        
    if(os.path.isfile(infile) == True and os.path.isfile(dinfile) == True):
        hist_mass_neutral[i] = compute_meanDensNeutral(infile, dinfile, double_precision, isPadded, inputIsDouble, gridsize)
    elif(i > 0):
        print "!!! densHI: replacing value at z=", redshift[i+1], "with previous value at", redshift[i]
        hist_mass_neutral[i] = hist_mass_neutral[i-1] 
    else:
        print "!!! densHI: replacing value at z=", redshift[i+1], "with previous value 1"
        hist_mass_neutral[i] = 1.
        
 
print "z =",redshift[1:]
print "XHII =", hist_ion

np.savetxt(outputfile + '.dat',np.c_[redshift[1:], hist_ion, hist_mass_ion, hist_mass_neutral])

#----------------------------------------------
#----------------------------------------------
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.rcParams.update({'figure.subplot.hspace':0.0, 'figure.subplot.bottom':0.1, 'figure.subplot.top':0.95, 'figure.subplot.left':0.1, 'figure.subplot.right':0.97})
plt.rcParams["figure.figsize"] = (7,5)

fig, axes = plt.subplots(1,1, sharex='col', sharey='row')
plt.subplots_adjust(wspace=0., hspace=0.1)

majorLocator   = mt.MultipleLocator(0.1)
minorLocator   = mt.MultipleLocator(0.02)
ymajorLocator   = mt.MultipleLocator(0.2)
yminorLocator   = mt.MultipleLocator(0.05)

#----------------------------------------------

axes.plot(1.-hist_ion, hist_mass_neutral, linestyle='-', color='black', label='$\langle 1+\delta\\rangle_\mathrm{HI}$')
axes.plot(1.-hist_ion, hist_mass_ion, linestyle='--', color='black', label='$\langle 1+\delta\\rangle_\mathrm{HII}$')

axes.yaxis.set_major_locator(ymajorLocator)
#axes.yaxis.set_major_formatter(ymajorFormatter)
axes.yaxis.set_minor_locator(yminorLocator)
axes.xaxis.set_major_locator(majorLocator)
#axes.xaxis.set_major_formatter(majorFormatter)
axes.xaxis.set_minor_locator(minorLocator)

axes.set_xlim((0.,1.))
axes.set_ylim((0.,3.))
  
axes.set_xlabel('$\langle \chi_\mathrm{HI}\\rangle$')
axes.set_ylabel('Overdensity')

plt.legend(bbox_to_anchor=(0.75, 0.88), loc=2, borderaxespad=0., frameon=True, labelspacing=0.4, handlelength=3, prop={'size':11})

#----------------------------------------------
fig.savefig(outputfile + '.png', format='png', dpi=512)#, transparent=True)

