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
	char igm_density_file[128];
	char snap_string[8];
	
	for(int i=0; i<128; i++) igm_density_file[i]='\0';
	if(snap >= 0)
	{
		sprintf(snap_string,"%03d",snap); 
		strcat(igm_density_file, simParam->igm_density_file);
		strcat(igm_density_file, "_");
		strcat(igm_density_file, snap_string);
	}else{
		strcat(igm_density_file, simParam->igm_density_file);
	}
  
	if(file_exist(simParam->igm_density_file) == 1)
	{
		read_igm_density(thisGrid, igm_density_file);
	}
	else if(file_exist(igm_density_file) == 1)
	{
		fprintf(stderr, "No source or nion file available, or names are incorrect!\n");
		exit(EXIT_FAILURE);
	}
}