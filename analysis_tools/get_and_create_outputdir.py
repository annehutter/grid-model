import sys
import ntpath
import os

import read_parameterfile as rp

def path_leaf(path):
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)

#def create_result_dir(infile):
inifile = sys.argv[1]

lines = rp.read_inifile(inifile)
solve_he = rp.identify_int(lines, rp.solve_he_str, rp.splitting_str)

directory = rp.get_directory(inifile)
basename = os.path.splitext(path_leaf(inifile))[0]

if(directory==''):
    newdir = 'results_' + basename
else:
    newdir = directory + '/results_' + basename

if not os.path.exists(newdir):
    os.makedirs(newdir)
    
dir_21cm = newdir + '/21cm'
if not os.path.exists(dir_21cm):
    os.makedirs(dir_21cm)
    
dir_ion = newdir + '/ion'
if not os.path.exists(dir_ion):
    os.makedirs(dir_ion)

dir_HeIIion = newdir + '/HeIIion'
if not os.path.exists(dir_HeIIion):
    os.makedirs(dir_HeIIion)
    
dir_HeIIIion = newdir + '/HeIIIion'
if not os.path.exists(dir_HeIIIion):
    os.makedirs(dir_HeIIIion)
    
dir_iondens = newdir + '/iondens'
if not os.path.exists(dir_iondens):
    os.makedirs(dir_iondens)

dir_cross_iondens = newdir + '/cross_ion_dens'
if not os.path.exists(dir_cross_iondens):
    os.makedirs(dir_cross_iondens)
    
if(solve_he != 0):
    dir_HeIIiondens = newdir + '/HeIIiondens'
    if not os.path.exists(dir_HeIIiondens):
        os.makedirs(dir_HeIIiondens)
        
    dir_HeIIIiondens = newdir + '/HeIIIiondens'
    if not os.path.exists(dir_HeIIIiondens):
        os.makedirs(dir_HeIIIiondens)
    
dir_neutral = newdir + '/neutral'
if not os.path.exists(dir_neutral):
    os.makedirs(dir_neutral)
  
dir_neutraldens = newdir + '/neutraldens'
if not os.path.exists(dir_neutraldens):
    os.makedirs(dir_neutraldens)
    
dir_density = newdir + '/density'
if not os.path.exists(dir_density):
    os.makedirs(dir_density)
  
dir_HI = newdir + '/HI_fields'
if not os.path.exists(dir_HI):
    os.makedirs(dir_HI)
    
if(solve_he == 1):
    dir_HeI = newdir + '/HeI_fields'
    if not os.path.exists(dir_HeI):
        os.makedirs(dir_HeI)
        
    dir_HeIII = newdir + '/HeIII_fields'
    if not os.path.exists(dir_HeIII):
        os.makedirs(dir_HeIII)
    
dir_cross_ziondens = newdir + '/cross_zion_dens'
if not os.path.exists(dir_cross_ziondens):
    os.makedirs(dir_cross_ziondens)
    
dir_ion_bubbledistr = newdir + '/ion_bubbledistr'
if not os.path.exists(dir_ion_bubbledistr):
    os.makedirs(dir_ion_bubbledistr)
 
dir_ionhist_cosmicweb = newdir + '/ionhist_cosmicweb'
if not os.path.exists(dir_ionhist_cosmicweb):
    os.makedirs(dir_ionhist_cosmicweb)
    
print newdir
