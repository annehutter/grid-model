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


def get_X(X, redshift, z):
    if(z >= redshift[0]):
        result = 0.
    else:
        index = 0
        for i in range(len(X)):
            if(z > redshift[i]):
                index = i
                break
        if(index == 0):
            result = 1.
        else:
            result = X[index-1] + (X[index] - X[index-1])/(redshift[index] - redshift[index-1]) * (z - redshift[index-1])
    return result
  
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

hubble_h = rp.identify_float(lines, rp.h_str, rp.splitting_str)
omega_b = rp.identify_float(lines, rp.omega_b_str, rp.splitting_str)
omega_l = rp.identify_float(lines, rp.omega_l_str, rp.splitting_str)
omega_m = rp.identify_float(lines, rp.omega_m_str, rp.splitting_str)
Y = rp.identify_float(lines, rp.Y_str, rp.splitting_str)

clight = 3.0e10
G = 6.67408e-8
mp = 1.673e-24
sigmaT = 6.65e-25
Mpc_cm = 3.0857e24
H0 = hubble_h * 1.e7 / Mpc_cm

histionfile = os.path.dirname(os.path.abspath(outputfile)) + '/hist_ion.dat'

print "\n--------------------------------"
print "Computing electron optical depth"
print "--------------------------------"

if(solve_he == 1):
    redshift, XHII, XHeII, XHeIII, XmHII, XmHeII, XmHeIII = np.loadtxt(histionfile, unpack=True, usecols=(0,1,2,3,4,5,6))
else:
    redshift, XHII, XmHII = np.loadtxt(histionfile, unpack=True, usecols=(0,1,2))
    XHeII = XHII
    XHeIII = np.zeros(len(XHII))

prefactor = clight * sigmaT * H0 * omega_b / (4. * np.pi * G * mp)

dz = 0.2
zhigh = np.max(redshift) + 1
nbins = np.int32(zhigh / dz + 1)

z = dz * np.arange(nbins)
tau = np.zeros(nbins)

xHII = np.zeros(nbins)
xHeII = np.zeros(nbins)
xHeIII = np.zeros(nbins)

for i in range(nbins):
    xHII[i] = get_X(XHII, redshift, z[i])
    xHeII[i] = get_X(XHeII, redshift, z[i])
    if(z[i] <=3.):
        xHeIII[i] = 1.
    else:
        xHeIII[i] = 0.

for i in range(nbins-1):
    termi = (omega_m * (1. + z[i])**3 + omega_l)**0.5
    termf = (omega_m * (1. + z[i+1])**3 + omega_l)**0.5

    if(i > 0):
        tau[i] = tau[i-1] + prefactor * (xHII[i] * (1.-Y) + (xHeII[i] + 2. * xHeIII[i]) * 0.25 * Y) * (termf - termi)
    else:
        tau[i] = prefactor * (xHII[i] * (1.-Y) + (xHeII[i] + 2. * xHeIII[i]) * 0.25 * Y) * (termf - termi)

np.savetxt(outputfile + '.dat',np.c_[z, tau])

#----------------------------------------------
#----------------------------------------------
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.rcParams.update({'figure.subplot.hspace':0.0, 'figure.subplot.bottom':0.1, 'figure.subplot.top':0.95, 'figure.subplot.left':0.1, 'figure.subplot.right':0.97})
plt.rcParams["figure.figsize"] = (6,4)
fig = plt.figure()

plt.subplots_adjust(wspace=0.0, hspace=0.1)

#----------------------------------------------
ax1 = plt.subplot()

redshift_range = z[len(z)-1] - z[0]
print redshift_range
if(redshift_range <= 1.0):
    xmaxL = 0.1
    xminL = 0.02
elif(redshift_range <= 5.0):
    xmaxL = 0.5
    xminL = 0.1
else:
    xmaxL = 1.
    xminL = 0.5

print tau
print tau[0]
print tau[len(tau)-1]
tau_range = tau[len(tau)-2] - tau[0] 
print tau_range
if(tau_range <= 0.02):
    ymaxL = 0.005
    yminL = 0.001
    ystring = '%0.3f'
elif(tau_range <= 0.05):
    ymaxL = 0.01
    yminL = 0.002
    ystring = '%0.3f'
elif(tau_range <= 0.1):
    ymaxL = 0.01
    yminL = 0.005
    ystring = '%0.3f'
else:
    ymaxL = 0.1
    yminL = 0.05
    ystring = '%0.2f'
    
print ymaxL
print yminL
    
majorLocator   = mt.MultipleLocator(xmaxL)
majorFormatter = mt.FormatStrFormatter('%0.1f')
minorLocator   = mt.MultipleLocator(xminL)

ymajorLocator   = mt.MultipleLocator(ymaxL)
ymajorFormatter = mt.FormatStrFormatter(ystring)
yminorLocator   = mt.MultipleLocator(yminL)

ax1.plot(z, tau, linestyle='-', color='black', label='$\langle\chi_\mathrm{HI}\\rangle$')

ax1.yaxis.set_major_locator(ymajorLocator)
ax1.yaxis.set_major_formatter(ymajorFormatter)
ax1.yaxis.set_minor_locator(yminorLocator)
ax1.xaxis.set_major_locator(majorLocator)
ax1.xaxis.set_major_formatter(majorFormatter)
ax1.xaxis.set_minor_locator(minorLocator)
		
ax1.set_xlabel('Redshift z')

ax1.set_ylabel('Optical depth $\\tau$')

ax1.set_xlim((0., z[len(z)-2]))
#plt.legend(bbox_to_anchor=(0.75, 0.88), loc=2, borderaxespad=0., frameon=True, labelspacing=0.4, handlelength=3, prop={'size':11})

#----------------------------------------------
fig.savefig(outputfile + '.png', format='png', dpi=512)#, transparent=True)

