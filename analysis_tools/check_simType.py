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
simulationtype = rp.identify_string(lines, rp.simulationtype_str, rp.splitting_str)
solve_he = rp.identify_int(lines, rp.solve_he_str, rp.splitting_str)

if(simulationtype == 'EVOLVE_ALL' or simulationtype == 'EVOLVE_BY_SNAPSHOT'):
    if(simulationtype == "EVOLVE_BY_SNAPSHOT"):
        snapshotstart = rp.identify_int(lines, rp.snapshotstart_str, rp.splitting_str)
    else:
        snapshotstart = 0
        
    densfile = rp.identify_string(lines, rp.densfile_str, rp.splitting_str)
    if(snapshotstart < 10):
        dinfile = densfile + '_00' + str(snapshotstart)
    else:
        dinfile = densfile + '_0' + str(snapshotstart)
    
    if(os.path.isfile(dinfile) == True):
        evolvingDensity = 1
    else:
        evolvingDensity = 0
else:
    evolvingDensity = 0

print solve_he
print evolvingDensity
