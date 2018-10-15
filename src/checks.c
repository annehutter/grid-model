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

#include "confObj.h"
#include "utils.h"
#include "checks.h"

int check_directory_exist(char *filename)
{
    char *directory = get_directory(filename);
    int result = -1;
    
    if(directory == NULL) result = 1;
    else
    {
        if(directory_exist(directory) == 0) result = 0;
        else result = 1;
        free(directory);
    }
    
    return result;
}

void check_output_directories_exist(confObj_t simParam)
{
    if(check_directory_exist(simParam->out_XHII_file) == 0)
    {
        printf("Directory for XHII output files does not exist\n");
        exit(EXIT_FAILURE);
    }
    if(simParam->write_photHI_file == 1)
    {
        if(check_directory_exist(simParam->out_photHI_file) == 0)
        {
            printf("Directory for photHI output files does not exist\n");
            exit(EXIT_FAILURE);
        }
    }
    if(simParam->solve_He == 1)
    {
        if(check_directory_exist(simParam->out_XHeII_file) == 0)
        {
            printf("Directory for XHeII output files does not exist\n");
            exit(EXIT_FAILURE);
        }
        if(check_directory_exist(simParam->out_XHeIII_file) == 0)
        {
            printf("Directory for XHeIII output files does not exist\n");
            exit(EXIT_FAILURE);
        }
    }
}