import numpy as np
import os
 
def read_ion(infile, isPadded, inputIsDouble, gridsize):
    fi = open(infile, 'rb')
    if(isPadded == 0):
        if(inputIsDouble == 1):
            ion = np.fromfile(fi, count=gridsize*gridsize*gridsize, dtype=np.float64)
        else:
            ion = np.fromfile(fi, count=gridsize*gridsize*gridsize, dtype=np.float32)
        ion.shape = (gridsize,gridsize,gridsize)
    else:
        if(inputIsDouble == 1):
            ion = np.fromfile(fi, count=isPadded*gridsize*gridsize*gridsize, dtype=np.float64)
        else:
            ion = np.fromfile(fi, count=isPadded*gridsize*gridsize*gridsize, dtype=np.float32)
        isPadded_factor = isPadded**(1./3.)
        ion.shape = (np.int32(isPadded_factor*gridsize),np.int32(isPadded_factor*gridsize),np.int32(isPadded_factor*gridsize))
        
        new_ion = np.array_split(ion, [gridsize], axis=2)
        new_ion = np.array_split(new_ion[0], [gridsize], axis=1)
        new_ion = np.array_split(new_ion[0], [gridsize], axis=0)
        ion = new_ion[0]

    fi.close()
    return ion


def read_dens(infile, isPadded, double_precision, gridsize):
    fd = open(infile, 'rb')
    
    if(isPadded == 0):
        if(double_precision == 1):
            dens = np.fromfile(fd, count=gridsize*gridsize*gridsize, dtype=np.float64)
        else:
            dens = np.fromfile(fd, count=gridsize*gridsize*gridsize, dtype=np.float32)
        dens.shape = (gridsize,gridsize,gridsize)
    else:
        if(double_precision == 1):
            dens = np.fromfile(fd, count=isPadded*gridsize*gridsize*gridsize, dtype=np.float64)
        else:
            dens = np.fromfile(fd, count=isPadded*gridsize*gridsize*gridsize, dtype=np.float32)
        isPadded_factor = isPadded**(1./3.)
        dens.shape = (np.int32(isPadded_factor*gridsize),np.int32(isPadded_factor*gridsize),np.int32(isPadded_factor*gridsize))
        
        
        new_dens = np.array_split(dens, [gridsize], axis=2)
        new_dens = np.array_split(new_dens[0], [gridsize], axis=1)
        new_dens = np.array_split(new_dens[0], [gridsize], axis=0)
        dens = new_dens[0]
        
    meanDens = np.mean(dens, dtype=np.float64)
    
    fd.close()
    return dens
