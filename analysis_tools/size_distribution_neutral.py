import sys
import os
import numpy as np
import matplotlib as m
m.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from statistics import *
import read_parameterfile as rp
import read_fields as rf

def round_down(num):
    if num < 0:
        return -np.ceil(abs(num))
    else:
        return np.int32(num)

inifile = sys.argv[1]
inputIsDouble = np.int32(sys.argv[2])
outputfile = sys.argv[3]

lines = rp.read_inifile(inifile)

redshiftfile = rp.identify_string(lines, rp.redshiftfile_str, rp.splitting_str)

simulationtype = rp.identify_string(lines, rp.simulationtype_str, rp.splitting_str)
if(simulationtype == "EVOLVE_BY_SNAPSHOT"):
  snapshotstart = rp.identify_int(lines, rp.snapshotstart_str, rp.splitting_str)
else:
  snapshotstart = 0
  
solve_he = rp.identify_int(lines, rp.solve_he_str, rp.splitting_str)

ionfile = rp.identify_string(lines, rp.ionfile_str, rp.splitting_str)
if(solve_he > 0):
    HeIIionfile = rp.identify_string(lines, rp.HeIIionfile_str, rp.splitting_str)
    HeIIIionfile = rp.identify_string(lines, rp.HeIIIionfile_str, rp.splitting_str)
    
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

print "\n-------------------------------------"
print "Computing neutral bubble distribution"
print "-------------------------------------"

#----------------------------------------------
#----------------------------------------------

if(solve_he > 0):
    numSpecie = 3
else:
    numSpecie = 1
    
for specie in range(numSpecie):
    if(specie == 0):
        inputfile = ionfile
        outputfile_tmp = outputfile + "_XHII"
    elif(specie == 1):
        inputfile = HeIIionfile
        outputfile_tmp = outputfile + "_XHeII"
    elif(specie == 2):
        inputfile = HeIIIionfile
        outputfile_tmp = outputfile + "_XHeIII"
    else:
        print "This specie does not exist."
        sys.exit()

            
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    plt.rcParams.update({'figure.subplot.hspace':0.0, 'figure.subplot.bottom':0.1, 'figure.subplot.top':0.95, 'figure.subplot.left':0.1, 'figure.subplot.right':0.87})
    plt.rcParams["figure.figsize"] = (7,5)
    fig = plt.figure()

    meanIon = np.zeros(len(redshift))

    maximum = 1.e-10

    threshold = 0.9
    Nrays = 100000

    for i in range(len(redshift)-1):
        z = redshift[i+1]
            
        if(i + snapshotstart < 10):
            infile = inputfile + '_0' + str(i + snapshotstart)
        else:
            infile = inputfile + '_' + str(i + snapshotstart)
        if(os.path.isfile(infile) == False):
            continue
    
        #ionization field
        ion = rf.read_ion(infile, isPadded, inputIsDouble, gridsize)
        meanIon[i] = np.mean(ion, dtype=np.float64)
        print "z =", z, "\tXHII =", meanIon[i]

        distance, histogram = calc_bubble_size_distribution(1.-ion, gridsize, boxsize, threshold, Nrays)
        
        if(distance[0] == None):
            continue
        if(i + snapshotstart < 10):
            outputfile_dat = outputfile_tmp + "_0" + str(i + snapshotstart) + ".dat"
        else:
            outputfile_dat = outputfile_tmp + "_" + str(i + snapshotstart) + ".dat"
        np.savetxt(outputfile_dat, np.c_[distance, histogram])
        
        plt.xscale('log')
        cmap = m.cm.get_cmap('jet_r')
        rgba = cmap(1.-meanIon[i])
        plt.plot(distance, histogram*distance, color=rgba)
        
        maximum_tmp = np.amax(histogram*distance)
        if(maximum_tmp > maximum):
            maximum = maximum_tmp
        
    #----------------------------------------------
    plt.xlabel('r  [ h$^{-1}$ Mpc ]')
    plt.ylabel('r dP / dr')

    plt.minorticks_on()

    if(maximum>0.5):
        maximum = 0.5
        
    plt.xlim([np.amin(distance), np.amax(distance)])
    plt.ylim([0, maximum])

    ax = fig.add_axes([0.89, 0.1, 0.03, 0.85])
    norm = m.colors.Normalize(vmin=0., vmax=1.)

    cb1 = m.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, orientation='vertical')
    cb1.set_label('$\langle\chi_\mathrm{HI}\\rangle$')

    outputfile_png = outputfile_tmp + '.png'
    fig.savefig(outputfile_png, format='png', dpi=512)#, transparent=True)
