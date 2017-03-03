#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <complex.h>
#include <assert.h>

#ifdef __MPI
#include <fftw3-mpi.h>
#include <mpi.h>
#else
#include <fftw3.h>
#endif

#include "utils.h"

#include "confObj.h"
#include "grid.h"
#include "input_grid.h"

/* read in / update sources or nion -------------------------------------------------------------*/
void read_update_igm_density(confObj_t simParam, grid_t *thisGrid, int snap)
{
	char igm_density_file[MAXLENGTH];
	char snap_string[8];
	
	for(int i=0; i<MAXLENGTH; i++) igm_density_file[i]='\0';
	if(snap >= 0)
	{
		sprintf(snap_string,"%03d",snap); 
		strcat(igm_density_file, simParam->igm_density_file);
		strcat(igm_density_file, "_");
		strcat(igm_density_file, snap_string);
		printf("\n reading %s\n", igm_density_file);
	}else{
		strcat(igm_density_file, simParam->igm_density_file);
		printf("\n reading %s\n", igm_density_file);
	}
  
	if(file_exist(igm_density_file) == 1)
	{
		read_array(thisGrid->igm_density, thisGrid, igm_density_file, simParam->input_doubleprecision);

		for(int i=0; i<thisGrid->nbins*thisGrid->nbins*thisGrid->local_n0; i++){
			if(creal(thisGrid->igm_density[i])<=0.){
				printf("density[%d] = %e\n",i,creal(thisGrid->igm_density[i]));
				thisGrid->igm_density[i] = 1.e-2 + 0.*I;
			}
		}
	}
	else if(file_exist(igm_density_file) == 1)
	{
		fprintf(stderr, "No density file available, or names are incorrect!\n");
		exit(EXIT_FAILURE);
	}
}

void read_update_igm_clump(confObj_t simParam, grid_t *thisGrid, int snap)
{
	char igm_clump_file[MAXLENGTH];
	char snap_string[8];
	
	for(int i=0; i<MAXLENGTH; i++) igm_clump_file[i]='\0';
	if(snap >= 0)
	{
		sprintf(snap_string,"%03d",snap); 
		strcat(igm_clump_file, simParam->igm_clump_file);
		strcat(igm_clump_file, "_");
		strcat(igm_clump_file, snap_string);
		printf("\n reading %s\n", igm_clump_file);
	}else{
		strcat(igm_clump_file, simParam->igm_clump_file);
		printf("\n reading %s\n", igm_clump_file);
	}
  
	if(file_exist(igm_clump_file) == 1)
	{
		read_array(thisGrid->igm_clump, thisGrid, igm_clump_file, simParam->input_doubleprecision);

		for(int i=0; i<thisGrid->nbins*thisGrid->nbins*thisGrid->local_n0; i++){
			if(creal(thisGrid->igm_clump[i])<=0.){
				printf("clump[%d] = %e\n",i,creal(thisGrid->igm_clump[i]));
				thisGrid->igm_clump[i] = 1. + 0.*I;
			}
		}
	}
	else if(file_exist(igm_clump_file) == 1)
	{
		printf("no clumping factor file exist; assume a clumping factor = 1\n");
	}
	
	double sum = 0.;
    double sum_total = 0.;
    int nbins = thisGrid->nbins;
    
	for(int i=0; i<thisGrid->local_n0; i++)
    {
        for(int j=0; j<nbins; j++)
		{
			for(int k=0; k<nbins; k++)
			{
                sum += creal(thisGrid->igm_clump[i*nbins*nbins+j*nbins+k]);
            }
        }
    }
    
    sum_total = sum;

#ifdef __MPI
	MPI_Allreduce(&sum, &sum_total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
    sum_total = sum_total/(nbins*nbins*nbins);
    printf(" mean clumping factor = %e\n", sum_total);
}
