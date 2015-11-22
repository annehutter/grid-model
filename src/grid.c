#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#ifdef __MPI
#include <fftw3-mpi.h>
#include <mpi.h>
#else
#include <fftw3.h>
#endif

#include "confObj.h"
#include "grid.h"

/* functions for grid -----------------------------------------------------------------------*/

grid_t *initGrid()
{
	grid_t *newGrid;
	newGrid = malloc(sizeof(grid_t));
	if(newGrid == NULL)
	{
		fprintf(stderr, "ERROR: initGrid Not enought memory to allocate Grid.\n");
                exit(EXIT_FAILURE);
	}
	
	newGrid->nbins = 0;
	newGrid->box_size =0.;
	
	newGrid->xmin = 0.;
	newGrid->ymin = 0.;
	newGrid->zmin = 0.;
	
	newGrid->igm_density = NULL;
	newGrid->halo_density = NULL;
	newGrid->igm_clump = NULL;
	newGrid->nion = NULL;
	
	newGrid->XHII = NULL;
	
	newGrid->local_n0 = 0;
	newGrid->local_0_start = 0;
	
	return newGrid;
}

void read_files_to_grid(grid_t *thisGrid, confObj_t thisInput)
{
	ptrdiff_t alloc_local, local_n0, local_0_start;
	int nbins;
	
	nbins = thisInput->grid_size;
	local_n0 = nbins;
	local_0_start = 0;
	
	thisGrid->nbins = nbins;
	thisGrid->box_size = thisInput->box_size;
	
	thisGrid->local_n0 = nbins;
	thisGrid->local_0_start = 0;
#ifdef __MPI
	fftw_mpi_init();
	
	alloc_local = fftw_mpi_local_size_3d(nbins, nbins, nbins, MPI_COMM_WORLD, &local_n0, &local_0_start);
	
	thisGrid->local_n0 = local_n0;
	thisGrid->local_0_start = local_0_start;
	
	thisGrid->igm_density = fftw_alloc_complex(alloc_local);
	thisGrid->halo_density = fftw_alloc_complex(alloc_local);
	thisGrid->igm_clump = fftw_alloc_complex(alloc_local);
	thisGrid->nion = fftw_alloc_complex(alloc_local);
	thisGrid->XHII = fftw_alloc_complex(alloc_local);
#else
	thisGrid->igm_density = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
	thisGrid->halo_density = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
	thisGrid->igm_clump = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
	thisGrid->nion = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
	thisGrid->XHII = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);

#endif
	read_grid(thisGrid->igm_density, nbins, local_n0, local_0_start, thisInput->igm_density_file);
	initialize_grid(thisGrid->halo_density, nbins, local_n0, local_0_start, 0.);
	read_grid(thisGrid->igm_clump, nbins, local_n0, local_0_start, thisInput->igm_clump_file);
// 	initialize_grid(thisGrid->igm_clump, nbins, local_n0, local_0_start, 1.);
	initialize_grid(thisGrid->nion, nbins, local_n0, local_0_start, 0.);
	initialize_grid(thisGrid->XHII, nbins, local_n0, local_0_start, 0.);
	
#ifdef __MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif
}

void read_grid(fftw_complex *toThisArray, int nbins, int local_n0, int local_0_start, char *filename)
{
	float *tmparray;
	FILE * fp;
		
	tmparray = (float*)malloc(sizeof(float)*local_n0*nbins*nbins);
#ifdef __MPI
	MPI_File mpifile;
	MPI_Offset offset;
	MPI_Status status;
	
	offset = (local_0_start*nbins*nbins*sizeof(float));
	
	MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY,MPI_INFO_NULL, &mpifile);
	MPI_File_read_at_all(mpifile,offset,tmparray, local_n0*nbins*nbins,MPI_FLOAT,&status);
	MPI_File_close(&mpifile);
#else
	fp = fopen(filename, "rb");
	fread(tmparray, sizeof(float), nbins*nbins*nbins, fp);
	fclose(fp);
#endif
	
	for(int i=0; i<local_n0; i++)
	{
		for(int j=0; j<nbins; j++)
		{
			for(int k=0; k<nbins; k++)
			{
				toThisArray[i*nbins*nbins+j*nbins+k] = (double)tmparray[i*nbins*nbins+j*nbins+k]+0.*I;
			}
		}
	}
	free(tmparray);
}

void initialize_grid(fftw_complex *thisArray, int nbins, int local_n0, int local_0_start, double value)
{
	for(int i=0; i<local_n0; i++)
	{
		for(int j=0; j<nbins; j++)
		{
			for(int k=0; k<nbins; k++)
			{
				thisArray[i*nbins*nbins+j*nbins+k] = value +0.*I;
			}
		}
	}
}

void deallocate_grid(grid_t *thisGrid)
{
	fftw_free(thisGrid->igm_density);
	fftw_free(thisGrid->halo_density);
	fftw_free(thisGrid->igm_clump);
	fftw_free(thisGrid->nion);
	fftw_free(thisGrid->XHII);
	free(thisGrid);
}

void write_grid_to_file_float(fftw_complex *thisArray, int nbins, int local_n0, int local_0_start, char *filename)
{
	float *tmparray;
	
	tmparray = (float*)malloc(sizeof(float)*local_n0*nbins*nbins);
	
	for(int i=0; i<local_n0; i++)
	{
		for(int j=0; j<nbins; j++)
		{
			for(int k=0; k<nbins; k++)
			{
				tmparray[i*nbins*nbins+j*nbins+k] = (float)creal(thisArray[i*nbins*nbins+j*nbins+k]);
// 				if(creal(thisArray[i*nbins*nbins+j*nbins+k])>1.e-13) printf("tmparray = %e\n",creal(thisArray[i*nbins*nbins+j*nbins+k]));
			}
		}
	}
		
#ifdef __MPI
	MPI_File mpifile;
	MPI_Offset offset;
	MPI_Status status;
	
	offset = (local_0_start*nbins*nbins*sizeof(float));

	MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &mpifile);
	MPI_File_write_at_all(mpifile, offset, tmparray, local_n0*nbins*nbins, MPI_FLOAT, &status);
	MPI_File_close(&mpifile);
#else
	FILE * fp;
	
	fp = fopen(filename, "wb");
	fwrite(tmparray, sizeof(float), nbins*nbins*nbins, fp);
	fclose(fp);
#endif
	
	free(tmparray);
}

void write_grid_to_file_double(fftw_complex *thisArray, int nbins, int local_n0, int local_0_start, char *filename)
{
	double *tmparray;
	
	tmparray = (double*)malloc(sizeof(double)*local_n0*nbins*nbins);
	
	for(int i=0; i<local_n0; i++)
	{
		for(int j=0; j<nbins; j++)
		{
			for(int k=0; k<nbins; k++)
			{
				tmparray[i*nbins*nbins+j*nbins+k] = (double)creal(thisArray[i*nbins*nbins+j*nbins+k]);
			}
		}
	}
	
#ifdef __MPI
	MPI_File mpifile;
	MPI_Offset offset;
	MPI_Status status;
	
	offset = (local_0_start*nbins*nbins*sizeof(double));
	
	MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &mpifile);
	MPI_File_write_at_all(mpifile, offset, tmparray, local_n0*nbins*nbins, MPI_DOUBLE, &status);
	MPI_File_close(&mpifile);
#else
	FILE * fp;
	
	fp = fopen(filename, "wb");
	fwrite(tmparray, sizeof(double), nbins*nbins*nbins, fp);
	fclose(fp);
#endif
	
	free(tmparray);
}

void save_to_file_XHII(grid_t *thisGrid, char *filename)
{
	write_grid_to_file_float(thisGrid->XHII, thisGrid->nbins, thisGrid->local_n0, thisGrid->local_0_start, filename);
}