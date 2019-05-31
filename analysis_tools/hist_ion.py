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

def compute_meanIon(infile, isPadded, inputIsDouble, gridsize):
    ion = rf.read_ion(infile, isPadded, inputIsDouble, gridsize)
    
    meanIon = np.mean(ion, dtype=np.float64)
    return meanIon

def compute_meanMassIon(infile, densfile, double_precision, isPadded, inputIsDouble, gridsize):
    ion = rf.read_ion(infile, isPadded, inputIsDouble, gridsize)
    dens = rf.read_dens(densfile, isPadded, double_precision, gridsize)

    meanIon = np.mean(ion*dens, dtype=np.float64)/np.mean(dens, dtype=np.float64)
    return meanIon

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
    
solve_he = rp.identify_int(lines, rp.solve_he_str, rp.splitting_str)
HeIIionfile = rp.identify_string(lines, rp.HeIIionfile_str, rp.splitting_str)
HeIIIionfile = rp.identify_string(lines, rp.HeIIIionfile_str, rp.splitting_str)

if(solve_he == 1):
    HeIIionfile = rp.identify_string(lines, rp.HeIIionfile_str, rp.splitting_str)
    HeIIIionfile = rp.identify_string(lines, rp.HeIIIionfile_str, rp.splitting_str)

redshift, snap = np.loadtxt(redshiftfile, unpack='True', skiprows=0, usecols=(0,1))

print "\n----------------------------"
print "Computing ionization history"
print "----------------------------"

hist_ion = np.zeros(len(redshift)-1)
hist_mass_ion = np.zeros(len(redshift)-1)

if(solve_he == 1):
    hist_HeIIion = np.zeros(len(redshift)-1)
    hist_HeIIIion = np.zeros(len(redshift)-1)
    hist_mass_HeIIion = np.zeros(len(redshift)-1)
    hist_mass_HeIIIion = np.zeros(len(redshift)-1)

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
        hist_mass_ion[i] = compute_meanMassIon(infile, dinfile, double_precision, isPadded, inputIsDouble, gridsize)
    elif(i > 0):
        print "!!! XmHII: replacing value at z=", redshift[i+1], "with previous value at", redshift[i]
        hist_mass_ion[i] = hist_mass_ion[i-1] 
    else:
        print "!!! XmHII: replacing value at z=", redshift[i+1], "with previous value 0"
        hist_mass_ion[i] = 0.
        
    if(solve_he == 1):
        if(i + snapshotstart < 10):
            HeIIinfile = HeIIionfile + '_0' + str(i + snapshotstart)
            HeIIIinfile = HeIIIionfile + '_0' + str(i + snapshotstart)
        else:
            HeIIinfile = HeIIionfile + '_' + str(i + snapshotstart)
            HeIIIinfile = HeIIIionfile + '_' + str(i + snapshotstart)
            
        if(os.path.isfile(HeIIinfile) == True and os.path.isfile(dinfile) == True):
            hist_HeIIion[i] = compute_meanIon(HeIIinfile, isPadded, inputIsDouble, gridsize)
            hist_mass_HeIIion[i] = compute_meanMassIon(HeIIinfile, dinfile, double_precision, isPadded, inputIsDouble, gridsize)
        elif(i > 0):
            print "!!! replacing value at z=", redshift[i+1], "with previous value at", redshift[i]
            hist_HeIIion[i] = hist_HeIIion[i-1]
            hist_mass_HeIIion[i] = hist_mass_HeIIion[i-1] 
        else:
            print "!!! replacing value at z=", redshift[i+1], "with previous value 0"
            hist_HeIIion[i] = 0.
            hist_mass_HeIIion[i] = 0.
        
        if(os.path.isfile(HeIIIinfile) == True and os.path.isfile(dinfile) == True):
            hist_HeIIIion[i] = compute_meanIon(HeIIIinfile, isPadded, inputIsDouble, gridsize)
            hist_mass_HeIIIion[i] = compute_meanMassIon(HeIIIinfile, dinfile, double_precision, isPadded, inputIsDouble, gridsize)
        elif(i > 0):
            print "!!! replacing value at z=", redshift[i+1], "with previous value at", redshift[i]
            hist_HeIIIion[i] = hist_HeIIIion[i-1]
            hist_mass_HeIIIion[i] = hist_mass_HeIIIion[i-1] 
        else:
            print "!!! replacing value at z=", redshift[i+1], "with previous value 0"
            hist_HeIIIion[i] = 0.
            hist_mass_HeIIIion[i] = 0.
            
print "z =",redshift[1:]
print "XHII =", hist_ion
if(solve_he == 1):
    print "XHeII =", hist_HeIIion
    print "XHeIII =", hist_HeIIIion

if(solve_he == 1):
    np.savetxt(outputfile + '.dat',np.c_[redshift[1:], hist_ion, hist_HeIIion, hist_HeIIIion, hist_mass_ion, hist_mass_HeIIion, hist_mass_HeIIIion])
else:
    np.savetxt(outputfile + '.dat',np.c_[redshift[1:], hist_ion, hist_mass_ion])

#----------------------------------------------
#----------------------------------------------
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.rcParams.update({'figure.subplot.hspace':0.0, 'figure.subplot.bottom':0.1, 'figure.subplot.top':0.95, 'figure.subplot.left':0.1, 'figure.subplot.right':0.97})
plt.rcParams["figure.figsize"] = (7,5)
fig = plt.figure()

gs = gridspec.GridSpec(2, 1, height_ratios=[1,2])
plt.subplots_adjust(wspace=0.0, hspace=0.1)

#----------------------------------------------
ax1 = plt.subplot(gs[1])

redshift_range = redshift[0]-redshift[len(redshift)-1]
if(redshift_range <= 1.0):
    xmaxL = 0.1
    xminL = 0.02
elif(redshift_range <= 5.0):
    xmaxL = 0.5
    xminL = 0.1
elif(redshift_range <= 10.0):
    xmaxL = 1.
    xminL = 0.5
else:
    xmaxL = 2.
    xminL = 0.5

logchi_range = np.log10(1.-hist_ion[0]) - np.log10(1.-hist_ion[len(hist_ion)-1]) 

if(solve_he == 1):
    logchiHeII_range = np.log10(1.-hist_HeIIion[0]) - np.log10(1.-hist_ion[len(hist_HeIIion)-1])
    logchiHeIII_range = np.log10(hist_HeIIIion[0]) - np.log10(hist_ion[len(hist_HeIIIion)-1])
    
    if(logchiHeII_range > logchi_range):
        logchi_range = logchiHeII_range
    if(logchiHeIII_range > logchi_range):
        logchi_range = logchiHeIII_range

if(logchi_range <= 1.0):
    ymaxL = 0.2
    yminL = 0.05
    ystring = '%0.1f'
elif(logchi_range <= 3.0):
    ymaxL = 0.5
    yminL = 0.1
    ystring = '%0.1f'
else:
    ymaxL = 1.
    yminL = 0.2
    ystring = '%d'
    
majorLocator   = mt.MultipleLocator(xmaxL)
majorFormatter = mt.FormatStrFormatter('%0.1f')
minorLocator   = mt.MultipleLocator(xminL)

ymajorLocator   = mt.MultipleLocator(ymaxL)
ymajorFormatter = mt.FormatStrFormatter(ystring)
yminorLocator   = mt.MultipleLocator(yminL)

ax1.plot(redshift[1:], np.log10(1.-hist_ion), linestyle='-', color='black', label='$\langle\chi_\mathrm{HI}\\rangle$')
ax1.plot(redshift[1:], np.log10(1.-hist_mass_ion), linestyle='--', color='black', label='$\langle\chi_\mathrm{HI}\\rangle^{(m)}$')

if(solve_he == 1):
    plt.plot(redshift[1:], np.log10(1.-hist_HeIIion), linestyle='-', color='blue', label='$\langle\chi_\mathrm{HeI}\\rangle$')
    plt.plot(redshift[1:], np.log10(1.-hist_mass_HeIIion), linestyle='--', color='blue', label='$\langle\chi_\mathrm{HeI}\\rangle^{(m)}$')

    plt.plot(redshift[1:], np.log10(hist_HeIIIion), linestyle='-', color='green', label='$\langle\chi_\mathrm{HeIII}\\rangle$')
    plt.plot(redshift[1:], np.log10(hist_mass_HeIIIion), linestyle='--', color='green', label='$\langle\chi_\mathrm{HeIII}\\rangle^{(m)}$')

ax1.yaxis.set_major_locator(ymajorLocator)
ax1.yaxis.set_major_formatter(ymajorFormatter)
ax1.yaxis.set_minor_locator(yminorLocator)
ax1.xaxis.set_major_locator(majorLocator)
ax1.xaxis.set_major_formatter(majorFormatter)
ax1.xaxis.set_minor_locator(minorLocator)
		
ax1.set_xlabel('z')

ax1.set_ylabel('Log $\langle\chi\\rangle$')

plt.legend(bbox_to_anchor=(0.75, 0.88), loc=2, borderaxespad=0., frameon=True, labelspacing=0.4, handlelength=3, prop={'size':11})

#----------------------------------------------
ax0 = plt.subplot(gs[0], sharex=ax1)
plt.setp(ax0.get_xticklabels(), visible=False)

factor_range = hist_mass_ion[0]/hist_ion[0] - hist_mass_ion[len(hist_mass_ion)-1]/hist_ion[len(hist_ion)-1]

if(solve_he == 1):
    factor_min = np.min([hist_mass_ion[0]/hist_ion[0], hist_mass_ion[len(hist_mass_ion)-1]/hist_ion[len(hist_ion)-1], hist_mass_HeIIion[0]/hist_HeIIion[0], hist_mass_HeIIion[len(hist_mass_HeIIion)-1]/hist_HeIIion[len(hist_HeIIion)-1], hist_mass_HeIIIion[0]/hist_HeIIIion[0], hist_mass_HeIIIion[len(hist_mass_HeIIIion)-1]/hist_HeIIIion[len(hist_HeIIIion)-1]])
    factor_max = np.max([hist_mass_ion[0]/hist_ion[0], hist_mass_ion[len(hist_mass_ion)-1]/hist_ion[len(hist_ion)-1], hist_mass_HeIIion[0]/hist_HeIIion[0], hist_mass_HeIIion[len(hist_mass_HeIIion)-1]/hist_HeIIion[len(hist_HeIIion)-1], hist_mass_HeIIIion[0]/hist_HeIIIion[0], hist_mass_HeIIIion[len(hist_mass_HeIIIion)-1]/hist_HeIIIion[len(hist_HeIIIion)-1]])
    factor_range = factor_max - factor_min
        
if(factor_range <= 1.0):
    ymaxL = 0.2
    yminL = 0.05
    ystring = '%0.1f'
elif(factor_range <= 3.0):
    ymaxL = 0.5
    yminL = 0.1
    ystring = '%0.1f'
else:
    ymaxL = 1.
    yminL = 0.2
    ystring = '%d'
    
ymajorLocator   = mt.MultipleLocator(ymaxL)
ymajorFormatter = mt.FormatStrFormatter(ystring)
yminorLocator   = mt.MultipleLocator(yminL)

ax0.plot(redshift[1:], hist_mass_ion/hist_ion, linestyle='-', color='black')

if(solve_he == 1):
    ax0.plot(redshift[1:], hist_mass_HeIIion/hist_HeIIion, linestyle='-', color='blue')
    ax0.plot(redshift[1:], hist_mass_HeIIIion/hist_HeIIIion, linestyle='-', color='green')

ax0.yaxis.set_major_locator(ymajorLocator)
ax0.yaxis.set_major_formatter(ymajorFormatter)
ax0.yaxis.set_minor_locator(yminorLocator)

ax0.set_ylabel('$\langle\chi\\rangle^{(m)} / \langle\chi\\rangle$')

#----------------------------------------------
fig.savefig(outputfile + '.png', format='png', dpi=512)#, transparent=True)

