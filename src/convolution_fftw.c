#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <complex.h>
#include <gsl/gsl_integration.h>

#ifdef __MPI
#include <fftw3-mpi.h>
#include <mpi.h>
#else
#include <fftw3.h>
#endif

#include "confObj.h"
#include "grid.h"

void convolve_fft(grid_t *thisGrid, fftw_complex *filter, fftw_complex *output, fftw_complex *input)
{
	int nbins;
	double factor;
#ifdef __MPI
	ptrdiff_t alloc_local, local_n0, local_0_start;
#else
	ptrdiff_t local_n0;
#endif
	
	fftw_complex *input_ft, *filter_ft;
	fftw_plan plan_input, plan_filter, plan_back;
	
	nbins = thisGrid->nbins;
	local_n0 = thisGrid->local_n0;
	
#ifdef __MPI
	local_0_start = thisGrid->local_0_start;
	
	alloc_local = fftw_mpi_local_size_3d(nbins, nbins, nbins, MPI_COMM_WORLD, &local_n0, &local_0_start);
	assert(local_0_start == thisGrid->local_0_start);
	assert(local_n0 == thisGrid->local_n0);
	
	input_ft = fftw_alloc_complex(alloc_local);
	filter_ft = fftw_alloc_complex(alloc_local);
	
	plan_input = fftw_mpi_plan_dft_3d(nbins, nbins, nbins, input, input_ft, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MPI_TRANSPOSED_OUT); //FFTW_MPI_TRANSPOSED_OUT
	plan_filter = fftw_mpi_plan_dft_3d(nbins, nbins, nbins, filter, filter_ft, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MPI_TRANSPOSED_OUT);
#else 
	input_ft = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
	if(input_ft == NULL)
	{
		fprintf(stderr, "input_ft in convolve_fft (filtering.c) could not be allocated\n");
		exit(EXIT_FAILURE);
	}
	filter_ft = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
	if(filter_ft == NULL)
	{
		fprintf(stderr, "filter_ft in convolve_fft (filtering.c) could not be allocated\n");
		exit(EXIT_FAILURE);
	}
	
	plan_input = fftw_plan_dft_3d(nbins, nbins, nbins, input, input_ft, FFTW_FORWARD, FFTW_ESTIMATE);
	plan_filter = fftw_plan_dft_3d(nbins, nbins, nbins, filter, filter_ft, FFTW_FORWARD, FFTW_ESTIMATE);
#endif
	
	fftw_execute(plan_input);
	fftw_execute(plan_filter);
	
	for(int i=0; i<local_n0; i++)
	{
		for(int j=0; j<nbins; j++)
		{
			for(int k=0; k<nbins; k++)
			{
				input_ft[i*nbins*nbins+j*nbins+k] = input_ft[i*nbins*nbins+j*nbins+k]*filter_ft[i*nbins*nbins+j*nbins+k];
			}
		}
	}
	
#ifdef __MPI
	plan_back = fftw_mpi_plan_dft_3d(nbins, nbins, nbins, input_ft, output, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MPI_TRANSPOSED_IN);
#else
	plan_back = fftw_plan_dft_3d(nbins, nbins, nbins, input_ft, output, FFTW_BACKWARD, FFTW_ESTIMATE);
#endif
	fftw_execute(plan_back);
	
	factor = 1./(nbins*nbins*nbins);
	
	for(int i=0; i<local_n0; i++)
	{
		for(int j=0; j<nbins; j++)
		{
			for(int k=0; k<nbins; k++)
			{
				output[i*nbins*nbins+j*nbins+k] = factor*output[i*nbins*nbins+j*nbins+k];
			}
		}
	}	
	
	fftw_destroy_plan(plan_input);
	fftw_destroy_plan(plan_filter);
	fftw_destroy_plan(plan_back);
	
	fftw_free(input_ft);
	fftw_free(filter_ft);
}
