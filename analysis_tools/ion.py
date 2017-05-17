from numpy.random import random
import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib as m
from scipy import ndimage
from grid import *
  
import read_parameterfile as rp

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
gridsize = rp.identify_int(lines, rp.gridsize_str, rp.splitting_str) #np.int32(sys.argv[6])
boxlength = rp.identify_float(lines, rp.boxsize_str, rp.splitting_str)

Omegab = rp.identify_float(lines, rp.omega_b_str, rp.splitting_str)
Omegam = rp.identify_float(lines, rp.omega_m_str, rp.splitting_str)
h = rp.identify_float(lines, rp.h_str, rp.splitting_str)

redshift, snap = np.loadtxt(redshiftfile, unpack='True', skiprows=0, usecols=(0,1))


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
    print z

    if(i<10):
        infile = ionfile + '_0' + str(i)
    else:
        infile = ionfile + '_' + str(i)
        
    #ionization field
    fi = open(infile, 'rb')
    
    if(isPadded == 0):
        if(inputIsDouble == 1):
            ion = np.fromfile(fi, count=gridsize*gridsize*gridsize, dtype=np.float64)
        else:
            ion = np.fromfile(fi, count=gridsize*gridsize*gridsize, dtype=np.float32)
        ion.shape = (gridsize,gridsize,gridsize)

    else:
        if(inputIsDouble == 1):
            ion = np.fromfile(fi, count=isPadded*isPadded*isPadded*gridsize*gridsize*gridsize, dtype=np.float64)
        else:
            ion = np.fromfile(fi, count=isPadded*isPadded*isPadded*gridsize*gridsize*gridsize, dtype=np.float32)
        ion.shape = (isPadded*gridsize,isPadded*gridsize,isPadded*gridsize)
        
        new_ion = np.array_split(ion, [gridsize], axis=2)
        new_ion = np.array_split(new_ion[0], [gridsize], axis=1)
        new_ion = np.array_split(new_ion[0], [gridsize], axis=0)
        print np.shape(new_ion[0])
        print len(new_ion)
        print  np.mean(new_ion[0], dtype=np.float64)
        ion = new_ion[0]


    meanIon[i] = np.mean(ion, dtype=np.float64)
    print meanIon[i]

    if(snap[i] != 0):
        if(counter <10):
            dinfile = densfile + '_00' + str(counter)
        else:
            dinfile = densfile + '_0' + str(counter)
        
    fd = open(dinfile, 'rb')
    
    if(isPadded == 0):
        if(double_precision == 1):
            dens = np.fromfile(fd, count=gridsize*gridsize*gridsize, dtype=np.float64)
        else:
            dens = np.fromfile(fd, count=gridsize*gridsize*gridsize, dtype=np.float32)
        dens.shape = (gridsize,gridsize,gridsize)
    else:
        if(double_precision == 1):
            dens = np.fromfile(fd, count=isPadded*isPadded*isPadded*gridsize*gridsize*gridsize, dtype=np.float64)
        else:
            dens = np.fromfile(fd, count=isPadded*isPadded*isPadded*gridsize*gridsize*gridsize, dtype=np.float32)
        dens.shape = (isPadded*gridsize,isPadded*gridsize,isPadded*gridsize)
        
        new_dens = np.array_split(dens, [gridsize], axis=2)
        new_dens = np.array_split(new_dens[0], [gridsize], axis=1)
        new_dens = np.array_split(new_dens[0], [gridsize], axis=0)
        dens = new_dens[0]
        
    meanDens = np.mean(dens, dtype=np.float64)
    
    modesIon = ifftn(ion*dens)
    kmid_bins, powerspec, p_err = modes_to_pspec(modesIon, boxsize=boxlength)

    if(minimum > np.min(powerspec[1:-1]*kmid_bins[1:-1]**3*k)):
        minimum = np.min(powerspec[1:-1]*kmid_bins[1:-1]**3*k)
    if(maximum < np.max(powerspec[1:-1]*kmid_bins[1:-1]**3*k)):
        maximum = np.max(powerspec[1:-1]*kmid_bins[1:-1]**3*k)

    plt.xscale('log')
    #plt.yscale('log')

    cmap = m.cm.get_cmap('jet_r')

    rgba = cmap(1.-meanIon[i])
    plt.plot(kmid_bins[1:-1], np.log10(powerspec[1:-1]*kmid_bins[1:-1]**3*k), color=rgba)
    
    
    np.savetxt("test_ion.dat",np.c_[kmid_bins, powerspec, p_err])
    

#----------------------------------------------
plt.xlabel('k  [ h Mpc$^{-1}$]')
plt.ylabel('Log ( $\Delta^2_{\\rho_\mathrm{HII}}$ )')

plt.minorticks_on()

plt.xlim([10.**round_down(np.log10(kmid_bins[1])), kmid_bins[len(kmid_bins)-2]])
plt.ylim([np.log10(minimum), np.log10(maximum)])

print minimum, maximum

ax = fig.add_axes([0.89, 0.1, 0.03, 0.85])
norm = m.colors.Normalize(vmin=0., vmax=1.)

cb1 = m.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, orientation='vertical')
cb1.set_label('$\langle\chi_\mathrm{HI}\\rangle$')

fig.savefig(outputfile, format='png', dpi=512)#, transparent=True)
