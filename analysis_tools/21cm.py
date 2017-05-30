from numpy.random import random
import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib as m
from scipy import ndimage
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
print densfile
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
    z = redshift[i]
    print z, Omegab, Omegam, h

    T0 = 28.5*((1+z)/10.)**0.5*(Omegab/0.042*h/0.73)*(0.24/Omegam)**0.5

    if(i<10):
        infile = ionfile + '_0' + str(i)
    else:
        infile = ionfile + '_' + str(i)
        
    #ionization field
    ion = rf.read_ion(infile, isPadded, inputIsDouble, gridsize)
    meanIon[i] = np.mean(ion, dtype=np.float64)
    print meanIon[i]

    if(snap[i] != 0):
        if(counter <10):
            dinfile = densfile + '_00' + str(counter)
        else:
            dinfile = densfile + '_0' + str(counter)
        counter = counter + 1
        
    dens = rf.read_dens(dinfile, isPadded, double_precision, gridsize)

    Tb = (1.-ion)*dens#/(1.-meanIon)

    modes21 = ifftn(Tb)
    kmid_bins_21, powerspec_21, p_err_21 = modes_to_pspec(modes21, boxsize=boxsize)

    if(minimum > np.min(powerspec_21[1:-1]*kmid_bins_21[1:-1]**3*k*T0*T0)):
        minimum = np.min(powerspec_21[1:-1]*kmid_bins_21[1:-1]**3*k*T0*T0)
    if(maximum < np.max(powerspec_21[1:-1]*kmid_bins_21[1:-1]**3*k*T0*T0)):
        maximum = np.max(powerspec_21[1:-1]*kmid_bins_21[1:-1]**3*k*T0*T0)

    plt.xscale('log')
    #plt.yscale('log')

    cmap = m.cm.get_cmap('jet_r')

    rgba = cmap(1.-meanIon[i])
    plt.plot(kmid_bins_21[1:-1], np.log10(powerspec_21[1:-1]*kmid_bins_21[1:-1]**3*k*T0*T0), color=rgba)
    
    if(i<10):
        outputfile_dat = outputfile+"_0"+str(i)+".dat"
    else:
        outputfile_dat = outputfile+"_"+str(i)+".dat"
    np.savetxt(outputfile_dat, np.c_[kmid_bins_21, powerspec_21, p_err_21])
    

#----------------------------------------------
plt.xlabel('k  [ h Mpc$^{-1}$]')
plt.ylabel('Log ( $\Delta^2_{21\mathrm{cm}}$ )  [ mK$^2$ ]')

plt.minorticks_on()

plt.xlim([10.**round_down(np.log10(kmid_bins_21[1])), kmid_bins_21[len(kmid_bins_21)-2]])
plt.ylim([np.log10(minimum), np.log10(maximum)])

print minimum, maximum

ax = fig.add_axes([0.89, 0.1, 0.03, 0.85])
norm = m.colors.Normalize(vmin=0., vmax=1.)

cb1 = m.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, orientation='vertical')
cb1.set_label('$\langle\chi_\mathrm{HI}\\rangle$')

outputfile_png = outputfile + '.png'
fig.savefig(outputfile_png, format='png', dpi=512)#, transparent=True)
