import numpy as np
import sys

#inifile = '../iniFile128.ini'

def identify_string(lines, word, splitting_str):
    for i, line in enumerate(lines):
        if word in line:
            return line.split(splitting_str)[1]

def identify_float(lines, word, splitting_str):
    for i, line in enumerate(lines):
        if word in line:
            return np.float(line.split(splitting_str)[1])

def identify_int(lines, word, splitting_str):
    for i, line in enumerate(lines):
        if word in line:
            return np.int(line.split(splitting_str)[1])

def read_inifile(inifile):
    with open(inifile, 'r') as f:
        lines = f.read().split("\n")
    return lines
    
redshiftfile_str = 'redshiftFile'

solve_he_str = 'solveForHelium'
padded_str = 'paddedBox'

gridsize_str = 'gridsize'
boxsize_str = 'boxsize'

double_precision_str = 'inputFilesAreInDoublePrecision'

densfile_str = 'inputIgmDensityFile'
meansdens_str = 'meanDensity'
clumpfile_str = 'inputIgmClumpFile'

sourcefile_str = 'inputSourceFile'
nionfile_str = 'inputNionFile'

ionfile_str = 'output_XHII_file'

h_str = 'hubble_h'
omega_b_str = 'omega_b'
omega_l_str = 'omega_l'
omega_m_str = 'omega_m'
Y_str = 'Y'

HeIIionfile_str = 'output_XHeII_file'
HeIIIionfile_str = 'output_XHeIII_file'

splitting_str = ' = '

#print 'start'

#lines = read_inifile(inifile)

#ionfile = identify_string(lines, ionfile_str, splitting_str)
#densfile =identify_string(lines, densfile_str, splitting_str)

#redshiftfile = identify_string(lines, redshiftfile_str, splitting_str)

#h = identify_float(lines, h_str, splitting_str)
#omega_b = identify_float(lines, omega_b_str, splitting_str)
#omega_l = identify_float(lines, omega_l_str, splitting_str)


#print ionfile
#print h
        
