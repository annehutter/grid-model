 #include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <complex.h>

#ifdef __MPI
#include <fftw3-mpi.h>
#include <mpi.h>
#else
#include <fftw3.h>
#endif

int file_exist(char *filename)
{
	FILE *file;
	if((file = fopen (filename, "rt"))){
		fclose(file);
		return 1;
	}else{
		return 0;
	}
}