import sys
import os
import numpy as np
import matplotlib as m
m.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from grid import *
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

redshiftfile = rp.identify_string(lines, rp.redshiftfile_str, rp.splitting_str) #sys.argv[4]

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

print "\n-------------------------------"
print "Computing density power spectra"
print "-------------------------------"

#----------------------------------------------
#----------------------------------------------
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.rcParams.update({'figure.subplot.hspace':0.0, 'figure.subplot.bottom':0.1, 'figure.subplot.top':0.95, 'figure.subplot.left':0.1, 'figure.subplot.right':0.87})
plt.rcParams["figure.figsize"] = (7,5)
fig = plt.figure()

k = (2*np.pi)**3/(2*np.pi**2)

meanIon = np.zeros(len(redshift))

minimum = 1.e10
maximum = 1.e-10

counter = 0
for i in range(len(redshift)-1):
    z = redshift[i+1]
    
    if(i<10):
        infile = ionfile + '_0' + str(i)
    else:
        infile = ionfile + '_' + str(i)
    if(os.path.isfile(infile) == False):
        continue
    
    #ionization field
    ion = rf.read_ion(infile, isPadded, inputIsDouble, gridsize)

    meanIon[i] = np.mean(ion, dtype=np.float64)
    
    if(snap[i] != 0):
        if(counter <10):
            dinfile = densfile + '_00' + str(counter)
        else:
            dinfile = densfile + '_0' + str(counter)
        counter = counter + 1
        if(os.path.isfile(dinfile) == False):
            continue

    print "z =", z, "\tXHII =", meanIon[i], "\tfile =", dinfile, "\tnew/old file =", snap[i]
    dens = rf.read_dens(dinfile, isPadded, double_precision, gridsize)
    
    modesDensity = ifftn(dens)
    kmid_bins, powerspec, p_err = modes_to_pspec(modesDensity, boxsize=boxsize)

    if(minimum > np.min(powerspec[1:-1]*kmid_bins[1:-1]**3*k)):
        minimum = np.min(powerspec[1:-1]*kmid_bins[1:-1]**3*k)
    if(maximum < np.max(powerspec[1:-1]*kmid_bins[1:-1]**3*k)):
        maximum = np.max(powerspec[1:-1]*kmid_bins[1:-1]**3*k)

    plt.xscale('log')
    #plt.yscale('log')

    cmap = m.cm.get_cmap('jet_r')

    rgba = cmap(1.-meanIon[i])
    plt.plot(kmid_bins[1:-1], np.log10(powerspec[1:-1]*kmid_bins[1:-1]**3*k), color=rgba)
    
    if(i<10):
        outputfile_dat = outputfile+"_0"+str(i)+".dat"
    else:
        outputfile_dat = outputfile+"_"+str(i)+".dat"
    np.savetxt(outputfile_dat,np.c_[kmid_bins, powerspec, p_err])
    

#----------------------------------------------
plt.xlabel('k  [ h Mpc$^{-1}$]')
plt.ylabel('Log ( $\Delta^2_{\\rho_\mathrm{HI}}$ )')

plt.minorticks_on()

plt.xlim([10.**round_down(np.log10(kmid_bins[1])), kmid_bins[len(kmid_bins)-2]])
plt.ylim([np.log10(minimum), np.log10(maximum)])

ax = fig.add_axes([0.89, 0.1, 0.03, 0.85])
norm = m.colors.Normalize(vmin=0., vmax=1.)

cb1 = m.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, orientation='vertical')
cb1.set_label('$\langle\chi_\mathrm{HI}\\rangle$')

outputfile_png = outputfile + '.png'
fig.savefig(outputfile_png, format='png', dpi=512)#, transparent=True)
