#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <complex.h>

#ifdef __MPI
#include <fftw3-mpi.h>
#include <mpi.h>
#else
#include <fftw3.h>
#endif

#include "input_redshifts.h"
#include "phys_const.h"
#include "confObj.h"
#include "utils.h"

double *read_redshift_list(char *redshift_file, int num_snapshots)
{
	double *redshift_list;
	FILE * file;
	
	int counter;
	double redshift, existFile;
	char line[128];
	
	if(file_exist(redshift_file) == 1)
	{
        redshift_list = initRedshift_list(num_snapshots);

		file = fopen(redshift_file, "rt");
		counter = 0;
		while(fgets(line, 128, file) != NULL)
		{
		      /* get a line, up to 80 chars from fr.  done if NULL */
		      sscanf (line, "%le\t%le", &redshift, &existFile);
		      /* convert the string to a long int */
		      redshift_list[counter*2] = redshift;
		      redshift_list[counter*2+1] = existFile;
		      counter ++;
		      if(counter > num_snapshots)
		      {
				fprintf(stderr, "number of snapshots in input file is larger than in redshift list (input_redshifts.c).\n");
				exit(EXIT_FAILURE);
		      }
		}
		assert(counter == num_snapshots);
		fclose(file);  /* close the file prior to exiting the routine */  
		return redshift_list;
	}else{
		return  NULL;
	}
}

double *initRedshift_list(int num_snapshots)
{
	double *redshift_list;
	
	redshift_list = malloc(sizeof(double)*num_snapshots*2);
	if(redshift_list == NULL)
	{
		fprintf(stderr, "redshift_list in initRedshift_list (input_redshifts.c) could not be allocated\n");
		exit(EXIT_FAILURE);
	}
	
	return redshift_list;
}

void deallocateRedshift_list(double *redshift_list)
{
	if(redshift_list != NULL) 
    {
        free(redshift_list);
    }
}
