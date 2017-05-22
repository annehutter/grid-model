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
    
dir_neutral = newdir + '/neutral'
if not os.path.exists(dir_neutral):
    os.makedirs(dir_neutral)
    
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
    
print newdir
