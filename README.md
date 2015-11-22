# grid-model
Code to compute ionization field from density fields and source catalogue

Requirements:
- Serial run: fftw3 library
- Parallel run with MPI: fftw3, fftw3-mpi and mpi libraries

Compilation:
make

Usage:
./cifog iniFile.ini

iniFile.ini contains the necessary input files and values for a run. For more details check example.
