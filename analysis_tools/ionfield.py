import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib.colors as cl
from matplotlib.font_manager import FontProperties
from scipy import ndimage
from scipy import signal

import ntpath
import os

import read_parameterfile as rp
import read_fields as rf

Mpc = 3.086e24
km = 1.e5
yr = 365.*24.*60.*60.

def path_leaf(path):
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)

def plot_field(infile, i, toPlot, cut_slice, str_time, str_redshift, str_meanIon):
    fig = plt.figure()

    image = plt.imshow(toPlot[:,:,cut_slice:cut_slice+1].mean(axis=-1),origin='lower',interpolation='nearest',  extent=[0,boxsize,0,boxsize], cmap="afmhot_r")

    plt.clim([-8.,0.])
    
    ax = plt.gca()
    ax.set_xlabel('x  [h$^{-1}$ Mpc]')
    ax.set_ylabel('y  [h$^{-1}$ Mpc]')
    clb = plt.colorbar(image)#, ticks=[0,0.5,1])
    clb.set_label('Log ( $\chi_\mathrm{HI}$ )')
    
    font = FontProperties()
    font.set_weight('bold')
    t1 = plt.text(3., 75., str_time, fontsize=12, color="black", fontproperties=font, bbox={'facecolor':'white', 'alpha':0.7, 'edgecolor':'none', 'pad':2})
    t2 = plt.text(68., 75., str_redshift, fontsize=12, color="black", fontproperties=font, bbox={'facecolor':'white', 'alpha':0.7, 'edgecolor':'none', 'pad':2})
    t3 = plt.text(3., 3., str_meanIon, fontsize=12, color="black", fontproperties=font, bbox={'facecolor':'white', 'alpha':0.7, 'edgecolor':'none', 'pad':2})

    directory = rp.get_directory(inifile)
    ininame = os.path.splitext(path_leaf(inifile))[0]
    if(directory==''):
        newdir = 'results_' + ininame
    else:
        newdir = directory + '/results_' + ininame
    newdir = newdir + '/HI_fields/'
    basename = newdir + os.path.splitext(path_leaf(infile))[0]
        
    if(i<10):
        outfile = basename + '_0' + str(i) + '_x' + str(cut_slice) + '.png'
    else:
        outfile = basename + '_' + str(i) + '_x' + str(cut_slice) + '.png'
        
    plt.savefig(outfile,format='png', dpi=512, transparent=True)
    plt.close()

def plot_field_HeII(infile, i, toPlot, cut_slice, str_time, str_redshift, str_meanIon):
    fig = plt.figure()

    image = plt.imshow(toPlot[:,:,cut_slice:cut_slice+1].mean(axis=-1),origin='lower',interpolation='nearest',  extent=[0,boxsize,0,boxsize], cmap="afmhot_r")

    plt.clim([-8.,0.])
    
    ax = plt.gca()
    ax.set_xlabel('x  [h$^{-1}$ Mpc]')
    ax.set_ylabel('y  [h$^{-1}$ Mpc]')
    clb = plt.colorbar(image)#, ticks=[0,0.5,1])
    clb.set_label('Log ( $\chi_\mathrm{HeI}$ )')

    font = FontProperties()
    font.set_weight('bold')
    t1 = plt.text(3., 75., str_time, fontsize=12, color="black", fontproperties=font, bbox={'facecolor':'white', 'alpha':0.7, 'edgecolor':'none', 'pad':2})
    t2 = plt.text(68., 75., str_redshift, fontsize=12, color="black", fontproperties=font, bbox={'facecolor':'white', 'alpha':0.7, 'edgecolor':'none', 'pad':2})
    t3 = plt.text(3., 3., str_meanIon, fontsize=12, color="black", fontproperties=font, bbox={'facecolor':'white', 'alpha':0.7, 'edgecolor':'none', 'pad':2})
        
    directory = rp.get_directory(inifile)
    ininame = os.path.splitext(path_leaf(inifile))[0]
    if(directory==''):
        newdir = 'results_' + ininame
    else:
        newdir = directory + '/results_' + ininame
    newdir = newdir + '/HeI_fields/'
    basename = newdir + os.path.splitext(path_leaf(infile))[0]    
    
    if(i<10):
        outfile = basename + '_0' + str(i) + '_x' + str(cut_slice) + '.png'
    else:
        outfile = basename + '_' + str(i) + '_x' + str(cut_slice) + '.png'
        
    plt.savefig(outfile,format='png', dpi=512, transparent=True)
    plt.close()
    
def plot_field_HeIII(infile, i, toPlot, cut_slice, str_time, str_redshift, str_meanIon):
    fig = plt.figure()

    image = plt.imshow(toPlot[:,:,cut_slice:cut_slice+1].mean(axis=-1),origin='lower',interpolation='nearest',  extent=[0,boxsize,0,boxsize], cmap="afmhot_r")

    plt.clim([-8.,0.])
    
    ax = plt.gca()
    ax.set_xlabel('x  [h$^{-1}$ Mpc]')
    ax.set_ylabel('y  [h$^{-1}$ Mpc]')
    clb = plt.colorbar(image)#, ticks=[0,0.5,1])
    clb.set_label('Log ( $\chi_\mathrm{HeIII}$ )')

    font = FontProperties()
    font.set_weight('bold')
    t1 = plt.text(3., 75., str_time, fontsize=12, color="black", fontproperties=font, bbox={'facecolor':'white', 'alpha':0.7, 'edgecolor':'none', 'pad':2})
    t2 = plt.text(68., 75., str_redshift, fontsize=12, color="black", fontproperties=font, bbox={'facecolor':'white', 'alpha':0.7, 'edgecolor':'none', 'pad':2})
    t3 = plt.text(3., 3., str_meanIon, fontsize=12, color="black", fontproperties=font, bbox={'facecolor':'white', 'alpha':0.7, 'edgecolor':'none', 'pad':2})
        
    directory = rp.get_directory(inifile)
    ininame = os.path.splitext(path_leaf(inifile))[0]
    if(directory==''):
        newdir = 'results_' + ininame
    else:
        newdir = directory + '/results_' + ininame
    newdir = newdir + '/HeIII_fields/'
    basename = newdir + os.path.splitext(path_leaf(infile))[0]
    
    if(i<10):
        outfile = basename + '_0' + str(i) + '_x' + str(cut_slice) + '.png'
    else:
        outfile = basename + '_' + str(i) + '_x' + str(cut_slice) + '.png'
        
    plt.savefig(outfile,format='png', dpi=512, transparent=True)
    plt.close()

def time_from_redshift_flatuniverse(zmin, zmax):
    prefactor = 2./(3*H0*omega_l**0.5)
    tmp = (omega_l/omega_m)**0.5
    return prefactor*(np.arcsinh(tmp*(1.+zmin)**-1.5) - np.arcsinh(tmp*(1.+zmax)**-1.5))

def redshift_from_time_flatuniverse(zmax, time):
    tmp = (omega_l/omega_m)**0.5
    tmp2 = np.sinh(1.5*H0*omega_l**0.5*time + np.arcsinh(tmp*(1.+zmax)**-1.5))
    return (tmp/tmp2)**(2./3.) - 1.

inifile = sys.argv[1]
inputIsDouble = np.int32(sys.argv[2])
cut_slice = int(sys.argv[3])

lines = rp.read_inifile(inifile)

redshiftfile = rp.identify_string(lines, rp.redshiftfile_str, rp.splitting_str)

ionfile = rp.identify_string(lines, rp.ionfile_str, rp.splitting_str)
double_precision = rp.identify_int(lines, rp.double_precision_str, rp.splitting_str)
isPadded = rp.identify_int(lines, rp.padded_str, rp.splitting_str)
isPadded_factor = isPadded**(1./3.)
if(isPadded != 0):
    gridsize = np.int32(rp.identify_int(lines, rp.gridsize_str, rp.splitting_str)/isPadded_factor)
    boxsize = rp.identify_float(lines, rp.boxsize_str, rp.splitting_str)/isPadded_factor
else:
    gridsize = rp.identify_int(lines, rp.gridsize_str, rp.splitting_str)
    boxsize = rp.identify_float(lines, rp.boxsize_str, rp.splitting_str)

omega_b = rp.identify_float(lines, rp.omega_b_str, rp.splitting_str)
omega_m = rp.identify_float(lines, rp.omega_m_str, rp.splitting_str)
omega_l = rp.identify_float(lines, rp.omega_l_str, rp.splitting_str)
h = rp.identify_float(lines, rp.h_str, rp.splitting_str)

solve_he = rp.identify_int(lines, rp.solve_he_str, rp.splitting_str)
HeIIionfile = rp.identify_string(lines, rp.HeIIionfile_str, rp.splitting_str)
HeIIIionfile = rp.identify_string(lines, rp.HeIIIionfile_str, rp.splitting_str)

if(solve_he == 1):
    HeIIionfile = rp.identify_string(lines, rp.HeIIionfile_str, rp.splitting_str)
    HeIIIionfile = rp.identify_string(lines, rp.HeIIIionfile_str, rp.splitting_str)

redshift, snap = np.loadtxt(redshiftfile, unpack='True', skiprows=0, usecols=(0,1))

H0 = h*1.e2*km/Mpc
time = time_from_redshift_flatuniverse(redshift, 1.e10)/yr

#----------------------------------------------
#----------------------------------------------
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.rcParams.update({'figure.subplot.hspace':0.0, 'figure.subplot.bottom':0.1, 'figure.subplot.top':0.95, 'figure.subplot.left':0.1, 'figure.subplot.right':0.97})
plt.rcParams["figure.figsize"] = (6,5)

for i in range(len(redshift)-1):        
    
    str_redshift = "z = " + str("%3.1f"%redshift[i])
    str_time = "t = " + str("%4.0f"%(time[i]/1.e6)) + " Myr"
    
    if(i<10):
        infile = ionfile + '_0' + str(i)
    else:
        infile = ionfile + '_' + str(i)
        
    #ionization field
    ion = rf.read_ion(infile, isPadded, inputIsDouble, gridsize)
    toPlot = np.log10(1.-ion)
    str_mean = "$\langle\chi_\mathrm{HI}\\rangle$ =" + str("%4.2f"%(np.mean(1.-ion, dtype=np.float64)))
    plot_field(infile, i, toPlot, cut_slice, str_time, str_redshift, str_mean)
    
    if(solve_he == 1):
        if(i<10):
            HeIIinfile = HeIIionfile + '_0' + str(i)
            HeIIIinfile = HeIIIionfile + '_0' + str(i)
        else:
            HeIIinfile = HeIIionfile + '_' + str(i)
            HeIIIinfile = HeIIIionfile + '_' + str(i)
        
        ion_HeII = rf.read_ion(HeIIinfile, isPadded, inputIsDouble, gridsize)
        ion_HeIII = rf.read_ion(HeIIIinfile, isPadded, inputIsDouble, gridsize)
        
        toPlot_HeII = np.log10(1.-ion_HeII)
        str_mean = "$\langle\chi_\mathrm{HeI}\\rangle$ =" + str("%4.2f"%(np.mean(1.-ion_HeII, dtype=np.float64)))
        plot_field_HeII(HeIIinfile, i, toPlot_HeII, cut_slice, str_time, str_redshift, str_mean)
        
        toPlot_HeIII = np.log10(ion_HeIII)
        str_mean = "$\langle\chi_\mathrm{HeIII}\\rangle$ =" + str("%4.2f"%(np.mean(ion_HeIII, dtype=np.float64)))
        plot_field_HeIII(HeIIIinfile, i, toPlot_HeIII, cut_slice, str_time, str_redshift, str_mean)
    

