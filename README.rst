Description
===========

Code to compute ionization field from density fields and source catalogues (or number of ionizing photon grids)


Why should you use it
=====================

1. **Modular** The code is written modular fashion, i.e. different modules such as ``solveForHelium`` and multiple ``photHImodel`` can be switched on or off or chosen.
2. **MPI Parallel** The code can be run on multiple cores and distributed memory.
3. **Residual HI fractions, recombinations & Helium** The code can compute residual HI fractions in ionized regions according to the chosen photoionization model; it accounts for HII, HeII and HeIII recombinations; it has the option to compute the HeII and HeIII ionization fields.

Installation
============

Pre-requisities
---------------

Serial run
``````````

1. fftw3 library: ``fftw3 >= 3.3.3``
2. gsl library (math, integrate): ``gsl >= 1.16``

Parallel run
````````````

1. MPI library
2. fftw3 & fftw3-mpi library: ``fftw3 >= 3.3.3``
3. gsl library (math, integrate): ``gsl >= 1.16``


Download & Build
----------------

::

    $ git clone https://github.com/annehutter/grid_model
    $ make

This will download the code and first test case from the github directory and compile the source code.

Execution
---------

The first test case can then be run by
::

    $ ./cifog iniFile128.ini

``iniFile128.ini`` contains all input parameters that are needed for any runs. For a different simulation the code does not need to be recompiled but just this parameter file iniFile.ini to be adapted.

Generation of the parameter file
````````````````````````````````
The parameter file can be adapted manually. However, since options can become complex, it is recommended to use program loacted under ``/create_parameter_file`` for the first runs, until options become more familiar.
::

    $ cd /create_inputfile
    $ make
    $ ./create_parameter_file
    
The program will guide you through the different options and generate the parameter file for you. It will also provide you with the information which other files are required.

Parameter file
''''''''''''''

**General**
...........

- ``calcIonHistory``: set to 1 if ionization history should be calculated (either from a single snapshot from redshift_prevSnapshot to redshift, or from the given redshift_file), otherwise 0.
- ``numSnapshots``: number of outputs (**Note**: output is automatically created for all redshifts where input files change)
- ``redshiftFile``: redshifts of in- and outputs: 1 for new input files & 0 for no new input file but output file
- ``redshift_prevSnapshot``: redshift to start the calculation, if no redshift_file is provided
- ``finalRedshift``: final redshift (for calcIonHistory = 1), or redshift of output (for calcIonHistory = 0)
- ``evolutionTime``: evolution time in Myrs if output for a single snapshot is chosen (calcIonHistory = 0) 

- ``size_linear_scale``: comoving size in h^{-1} Mpc until which tophat kernel should be increased linearly
- ``first_increment_in_logscale``: increment of the tophat kernel beyond linear increase
- ``max_scale``: maximum tophat kernel size in h^{-1} Mpc

- ``useDefaultMeanDensity``: set to 1 if default cosmological density value should be used (recommended), otherwise set to 0 if "meanDensity" is used

- ``useWebModel``: set to 1 if the residual HI fraction in ionized regions should be computed (this mode will require to choose a photHI model), otherwise 0
- ``constantPhotHI``: set to 1 if HI fraction should be calculated from a constant photoionization field, otherwise 0.
- ``calcMeanFreePath``: set to 1 if mfp is calculated from the size of the ionized regions and/or as in Miralda 2000, otherwise 0 (only applicable for constantPhotHI = 0)
- ``constantRecombinations``: set to 1 if rembination rate should be constant spatially, otherwise 0
- ``calcRecombinations``: set to 1 if number of recombinations should be calculated, otherwise 0
- ``solveForHelium``: set to 1 if HeII and HeIII fields should be computed, otherwise 0

- ``paddedBox``: set to 1 if a padded box is used, otherwise 0

**Input**
.........

- ``gridsize``: size of the grid (should be a power of 2)
- ``boxsize``: comoving boxsize in Mpc/h

- ``inputFilesAreInDoublePrecision``: 0 for single, 1 for double precision of data files to be read in
- ``inputFilesAreComoving``: set to 1 if input files are comoving, otherwise 0

- ``inputIgmDensityFile``: name of density file containing 3D density grid (if multiple then just the basename and neglecting extensions _00i)
- ``densityInOverdensity``: set to 1 if density is in terms of overdensity i.e. rho/mean(rho), otherwise 0
- ``meanDensity``: assumed mean density, density is evolved as dens(z) = meanDensity*(1+z)^3 (only effective when ``useDefaultMeanDensity=0``)

- ``inputIgmClumpFile``: name of clumping factor file, which is used to calculate the HI fraction at the listed outputs

- ``inputSourcesFile``: (if existing) file containing the sources (first line: #sources; every other line: x, y, z, Nion [s^-1], ID, fesc)
- ``inputNionFile``: (if existing) name of file containing 3D grid of Nion [s^-1]

**Output**
..........

- ``output_XHII_file``: basename for output of XHII fields
- ``write_photHI_file``: set to 1 if photoionization file should be written
- ``output_photHI_file``: basename for output of HI photoionization fields

**Cosmology**
.............

- ``h``: H = 100*h km/s/Mpc
- ``omega_b``: baryon density parameter
- ``omega_m``: matter density parameter
- ``omega_l``: lambda density parameter
- ``sigma8``: sigma8
- ``Y``: mass fraction of Helium in the primordial gas (assumed to consist of H and He)

**Photoionization**
...................

- ``photHI_bg_file``: name of file with a list of redshift, HI photoionization rates, HI photoheating rates, Q
- ``photHI_bg``: photoionization background value
- ``meanFreePathInIonizedMedium``: mfp in physical Mpc (only applicable for calcMeanFreePath = 0)
- ``sourceSlopeIndex``: spectral index of the spectrum of the ionizing sources, i.e. alpha for L_nu ~ nu^-alpha

**Recombinations**
..................

- ``dnrec_dt``: recombination rate when option ``constantRecombinations = 0`` is chosen.
- ``recombinationTable``: (table of recombination values, only change if you know exactly what you are doing! Below are the parameters of the table)
- ``zmin``: minimum redshift of recombination table
- ``zmax``: maximum redshift of recombination table
- ``dz``: increment in redshift in the recombination table
- ``fmin``: minimum factor (``= recombination rate/photionization rate in 10^{-12}s``) of recombination table
- ``fmax`` maximum factor (``= recombination rate/photionization rate in 10^{-12}s``) of recombination table
- ``df``: increment in factor in the recombination table
- ``dcellmin``: minimum dcell^{-1/3} of recombination table
- ``dcellmax``: minimum dcell^{-1/3} of recombination table
- ``ddcell``: increment in dcell^{-1/3} in the recombination table

**Helium**
..........

- ``inputSourcesHeIFile``: (if existing) file containing the sources (x, y, z, Nion_HeI [s^-1], ID, fesc)
- ``inputNionHeIFile``: (if existing) name of file containing 3D grid of Nion_HeI [s^-1]
- ``inputSourcesHeIFile``: (if existing) file containing the sources (x, y, z, Nion_HeII [s^-1], ID, fesc)
- ``inputNionHeIFile``: (if existing) name of file containing 3D grid of Nion_HeII [s^-1]

- ``output_XHeII_file``: output name for XHeII fields
- ``output_XHeIII_file``: output name for XHeIII fields


Options
=======

Helium
------

You can generate the corresponding input files of the ionizing photons of helium in **sourceFile format** by
::

    $ cd create_helium_nion_inputfiles/
    $ make
    $ ./create_helium_inputfiles

Before executing you may want to adjust the (in the directory) included iniFile, which lets you choose the in-and output names, the cosmology and the spectral shape of the sources.

HI photoionization models
-------------------------

1. ``photHI_model = 1``: This model assumes the photoionization rate to drop of as exp(-r/mfp)/r^2, whereas mfp is the mean mean free path of or in the ionized regions

2. ``photHI_model = 2``: This model computes the photoionization rate according to the mean free path of each cell. The mean free path corresponds to the filtering scale at which the cell became ionized.


