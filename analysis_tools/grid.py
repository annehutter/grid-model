from numpy.random import RandomState
from numpy import arange, pi, square, add, sqrt, empty, float64, mgrid, floor, reshape, int64, zeros, float32, bincount, digitize, array, remainder
from numpy.fft import fftn, ifftn
#from scipy.ndimage import map_coordinates

def make_k_values(boxsize, ngrid):
    """ 
    build the grid of |k| for the given box size, in a vaguely memory efficient way
    returns k (n,) array, inv_k2 (n,n,n) array, with 0 at i=j=k=0
    """
    # 1d component of k
    k1 = arange(ngrid)
    k1[1+ngrid/2:] -= ngrid
    k1 = k1 * (2 * pi / float(boxsize))
    k2 = square(k1)
    kmag = add.outer(add.outer(k2, k2), k2)
    inv_k2 = kmag.copy()
    kmag = sqrt(kmag)
    inv_k2[0,0,0] = 1.0
    inv_k2 = 1.0 / inv_k2
    inv_k2[0,0,0] = 0.0
    return k1, kmag, inv_k2
    


def powerspec_bins(ngrid, boxsize):
    """
    find power spectrum bins to for a cubic grid of size ngrid^3 of fourier modes.
    Assumes the FFT convention of 0, ..., n/2, -n/2+1, ..., -1
    ngrid   - num cells on side of cube
    boxsize - size of the box in real space

    returns kmin, kmax, kbins, kvol
    kmin  - the lower bound of the bin
    kmax  - the upper bound of the bin
    kbins - index (0, ..., m) of the bin of each cell
    kvol  - the volume in k space (number of cells of that type divided by k vol of each cell)
    """

    mid = ngrid/2
    # find the magnitude of the indices (i.e. ix**2+iy**2+iz**2 in the FFT convention)
    n1 = arange(ngrid)
    n1[1+mid:] -= ngrid
    n2 = square(n1)
    nmag = sqrt(add.outer(add.outer(n2, n2), n2)).ravel()
    
    nbins = (-1,) + tuple(arange(mid-1)+1.5) + (ngrid*2,)
    #print 'nbins', nbins
    kbins = digitize(nmag, nbins) - 1
    assert(kbins.min()==0)
    assert(kbins.max()==len(nbins)-2)
    
    # multiplier to go to k-space
    kmult = 2.0 * pi / boxsize

    kmin = (array(nbins) * kmult)[:-1]
    kmin[0] = 0

    kmax = (array(nbins) * kmult)[1:]
    kmax[-1] = mid * kmult * sqrt(3.0)
    
    kvol = bincount(kbins) * (kmult * kmult * kmult)
    return kmin, kmax, kbins, kvol
    
def build_displacement(boxsize, ngrid, power_func, seed=12345):
    """ 
    Build a grid of the displacement field using the given power spectrum

    boxsize    - Size of the box, units (e.g. Mpc/h) must be consistent with
                 the power spectrum function
    ngrid      - Integer size of the grid, e.g. 32 for a 32x32x32 grid 
    power_func - Power spectrum function to be used on each k-mode, units
                 should be inverse of the box size (e.g. h/Mpc)
    [seed]     - Optional seed for the random number generator.

    returns 
    disp_grid  - (3,ngrid,ngrid,ngrid) array of displacements, with
                 disp_grid[i] giving the displacement along the i-th axis etc.
    """
    mersenne = RandomState(seed)

    # make the k values 
    k, kmag, inv_k2 = make_k_values(boxsize, ngrid)
    dk = 2 * pi / boxsize # Distance between the modes (k)
    # amplitudes of the double integral of density (proportional to potential)
    mode_ampl = sqrt(power_func(kmag.ravel()) * dk * dk * dk) 
    mode_ampl.shape = (ngrid,ngrid,ngrid)
    white_real = mersenne.standard_normal(size=(ngrid,ngrid,ngrid))
    white_imag = mersenne.standard_normal(size=(ngrid,ngrid,ngrid))
    
    grf_k = (white_real + 1.0j*white_imag) * mode_ampl

    disp_fld = empty((3, ngrid, ngrid, ngrid), dtype=float64)

    k_shapes = ((ngrid,1,1), (ngrid,1), (ngrid,))

    for axis, k_shape in enumerate(k_shapes):
        ik_axis = 1.0j * reshape(k, k_shape) * inv_k2
        disp_fld[axis] = fftn(grf_k * ik_axis).real

    rms_disp = sqrt(square(disp_fld).sum(1).mean(dtype=float64))
    print 'RMS displacement', rms_disp
    return disp_fld

def interpolate_displacement_grid(uni_pts, disp_grid, boxsize, order=1):
    """ 
    Displace points via interpolation of displacement field 
    [order = 1] use linear interpolation (3=cubic spline)
    """
    assert(uni_pts.shape[0]==3)

    disps = empty((uni_pts.shape[1], uni_pts.shape[0]), dtype=float64)

    dx = boxsize / float(disp_grid.shape[1])

    coords = uni_pts * (1.0 / dx) # go from positions to indices
    print 'Indices in', coords.min(), coords.max(), 'should be in 0-%d'%(disp_grid.shape[1])
    # for each coordinate use the right displacement field
    for axis in range(3):
        disps[:,axis] = map_coordinates(disp_grid[axis], coords, None, order=order, mode='wrap') 
    
    return disps

def corner_idx_wt(corner, idx, wts, ngrid):
    if corner==0:
        return idx, 1.0 - wts
    else:
        return (idx + 1) % ngrid, wts

def cic_modes(pts, boxsize, ngrid, masses=None):
    """ find the power spec using cic """
    if pts.shape[0]!=3:
        raise Exception('cic_modes needs (3,n) array for point positions')
#    wrapped_pts = pts - boxsize * floor(pts * (1.0 / boxsize))
    wrapped_pts = remainder(pts, boxsize)
    dx = boxsize / float(ngrid)
    
    scaled_pts = wrapped_pts * (1.0/dx)
    f = floor(scaled_pts)
#    wts = scaled_pts - f
    wts = remainder(scaled_pts, 1.0)
    # idx of which grid cell each point lies in along each axis
    idx = f.astype(int64) 

    density = zeros(ngrid*ngrid*ngrid, dtype=float32)
    wt_tot = empty((pts.shape[1],), dtype=float32)
    ind_tot = empty(wt_tot.shape, dtype=int64)
    print 'Binning 8 corners',
    for corner in range(8):
        print corner,
        for axis in range(3):
            # is the cell corner to the left(=0) or right(=1) of the particle?
            leftright = (corner >> axis) & 1 
            # return the idx and 1d cloud-in-cell weight of this
            ind, wt = corner_idx_wt(leftright, idx[axis], wts[axis], ngrid)

            # 3d cloud-in-cell weights are the products of the 1d ones
            # indices into 3d array built similarly
            if axis==0:
                ind_tot[:] = ind
                wt_tot[:] = wt
            else:
                ind_tot *= ngrid
                ind_tot += ind
                wt_tot *= wt
        if masses is None:
            dens_part = bincount(ind_tot, weights=wt_tot)
        else:
            dens_part = bincount(ind_tot, weights=wt_tot * masses)
        density[:dens_part.size] += dens_part

    density *= 1.0 / density.mean(dtype=float64)  # Bug when the mean accumulator is less than float64
    density -= 1.0
    density.shape = (ngrid, ngrid, ngrid)
    
    modes = ifftn(density)

    return modes
    
def modes_to_pspec(modes, boxsize):
    """ 
    From a given set of fourier modes, compute the (binned) power spectrum 
    with errors.

    modes   - (n,n,n) array (from ifftn of overdensity values)
    boxsize - size of box (e.g. in Mpc/h)
    
    returns

    kmid    - k values of power spec
    pspec   - estimate of power spec at k
    perr    - error on estimate
    """

    ngrid = modes.shape[0]
    assert(modes.shape==(ngrid, ngrid,ngrid))
    kmin, kmax, kbins, kvol = powerspec_bins(ngrid, boxsize)

    wts = square(modes.ravel().real) + square(modes.ravel().imag)

    v1 = bincount(kbins, weights=wts) 
    powerspec = v1 * (1.0 / kvol)

    # work out error on power spectrum
    v2 = bincount(kbins, weights=square(wts))
    v0 = bincount(kbins)
    p_err = sqrt((v2*v0 - v1*v1)/(v0-1)) / kvol 

    kmid_bins = 0.5 * (kmin+kmax)

    return kmid_bins, powerspec, p_err
  
def two_modes_to_pspec(modes, modes2, boxsize):
    """ 
    From a given set of fourier modes, compute the (binned) power spectrum 
    with errors.

    modes   - (n,n,n) array (from ifftn of overdensity values)
    boxsize - size of box (e.g. in Mpc/h)
    
    returns

    kmid    - k values of power spec
    pspec   - estimate of power spec at k
    perr    - error on estimate
    """

    ngrid = modes.shape[0]
    assert(modes.shape==(ngrid, ngrid,ngrid))
    assert(modes2.shape==(ngrid, ngrid, ngrid))
    kmin, kmax, kbins, kvol = powerspec_bins(ngrid, boxsize)

    wts = modes.ravel().real*modes2.ravel().real + modes.ravel().imag*modes2.ravel().imag

    checksum = (modes.ravel().real*modes2.ravel().imag-modes2.ravel().real*modes.ravel().imag)
        
    v1 = bincount(kbins, weights=wts) 
    powerspec = v1 * (1.0 / kvol)
    
    v1check = bincount(kbins, weights=checksum)

    # work out error on power spectrum
    v2 = bincount(kbins, weights=square(wts))
    v0 = bincount(kbins)
    p_err = sqrt((v2*v0 - v1*v1)/(v0-1)) / kvol 

    kmid_bins = 0.5 * (kmin+kmax)

    return kmid_bins, powerspec, p_err
             
       
                    

