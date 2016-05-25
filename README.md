# grid-model
Code to compute ionization field from density fields and source catalogue

Requirements:
- Serial run: fftw3 library, gsl library (math, integrate)
- Parallel run with MPI: fftw3, fftw3-mpi and mpi libraries, gsl library (math, integrate)

Compilation:
make

Usage:
./cifog iniFile.ini

iniFile.ini contains the necessary input files and values for a run. For more details check example.


Parameters in iniFile.ini

- inputFilesAreInDoublePrecision: 0 for single, 1 for double precision of data files to be read in
- numSnapshots: number of outputs (output is automatically created for all redshifts where input files change)
- redshiftFile: redshifts of outputs, and 1 for new input files & 0 for no new input file
- inputIgmDensityFile: name of density file containing 3D density grid(if multiple, then of the basename and extensions _00i)

- gridsize: size of the grid (should be a power of 2)
- boxsize: comoving boxsize in Mpc/h
- size_linear_scale: comoving size in Mpc/h until which tophat kernel should be increased linearly
- first_increment_in_logscale: increment of tophat kernel beyond linear increase

- inputSourcesFile: (if existing) file containing the sources (x, y, z, Nion [s^-1], ID, fesc)
- inputNionFile: (if existing) name of file containing 3D grid of Nion [s^-1]
- evolutionTime: (if output for a single snapshot is chosen (calcIonHistory = 0)) t [Myrs]
- finalRedshift: final redshift (for calcIonHistory = 1), redshift of output (for calcIonHistory = 0)

- densityInOverdensity: set to 1 if density is in terms of overdensity, i.e. rho/mean(rho), otherwise 0
- meanDensity: assumed mean density, density is evolved as dens(z) = meanDensity*(1+z)^3
- useDefaultMeanDensity: set to 1 if default cosmological density value should be used (recommended), otherwise set to 0 "meanDensity" is used

- inputFilesAreComoving: set to 1 if input files are comoving, otherwise 0

- h: H = 100*h km/s/Mpc
- omega_b: baryon density parameter
- omega_m: matter density parameter
- omega_l: lambda density parameter
- sigma8: sigma8
- Y: mass fraction of Helium in the primordial gas (assumed to consist of H and He)

- output_XHII_file: output name for XHII fields

- useWebModel: set to 1 if web model as outlined in Sobacchi 2014 should be used, otherwise 0
- constantPhotHI: set to 1 if XHI should be calculated from a constant photoionization field, otherwise 0.
- photHI_bg: photoionization background value (used when constantPhotHI = 1)
- calcMeanFreePath: set to 1 if mfp is calculated as in Mirlada 2000, otherwise 0 (only applicable for constantPhotHI = 0)
- meanFreePathInIonizedMedium: mfp in physical Mpc (only applicable for calcMeanFreePath = 0)

- write_photHI_file: set to 1 if photoionization file should be written
- output_photHI_file: output name for photoionization field (comment: currently the one file is overwritten when new output is written, is going to be adapted in future versions)

- calcRecombinations: set to 1 if number of recombinations from redshift og first ionization should be calculated, otherwise 0
- recombinationTable: (table of recombination values, only change if you know exactly what you are doing! Below are the parameters of the table) nrec_values_batch_z6_20_0.01_f-9_2_0.1_d-4_4_0.1.dat
- zmin = 6.
- zmax = 20.
- dz = 0.01
- fmin = -9.
- fmax = 2.
- df = 0.1
- dcellmin = -4.
- dcellmax = 4.
- ddcell = 0.1

- calcIonHistory: set to 1 if ionization history should be calculated (either from a single snapshot from
- redshift_prevSnapshot to redshift, or from the given redshift_file), otherwise 0.

- readNrecFile: set to 1 if you want to provide a 3D grid of number of recombinations, otherwise 0.
- redshift_prevSnapshot = redshift to start the calculation, if no redshift_file is provided
- inputRecombFile = input name for the 3D grid of number of recombinations (comment: currently only one file can be read, is going to be adapted in future versions)
- outputRecombFile = output name for the number of recombination field (comment: currently the one file is overwritten when new output is written, is going to be adapted in future versions)
