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
    newGrid->lin_scales = 0.;
    newGrid->inc_log_scales = 0.;
    newGrid->max_scale = 0.;
    
    newGrid->xmin = 0.;
    newGrid->ymin = 0.;
    newGrid->zmin = 0.;
    
    newGrid->igm_density = NULL;
    newGrid->igm_clump = NULL;
    
    //hydrogen
    newGrid->nion = NULL;
    newGrid->cum_nion = NULL;
    newGrid->cum_nrec = NULL;
    newGrid->cum_nabs = NULL;
    newGrid->frac_Q = NULL;
    
    newGrid->XHII = NULL;
    newGrid->nrec = NULL;
    
    newGrid->photHI = NULL;
    newGrid->mean_photHI = 0.;
    
    newGrid->mean_mfp = 0.;
    
    //helium
    newGrid->nion_HeI = NULL;
    newGrid->nion_HeII = NULL;
    newGrid->cum_nion_HeI = NULL;
    newGrid->cum_nion_HeII = NULL;
    newGrid->cum_nrec_HeI = NULL;
    newGrid->cum_nrec_HeII = NULL;
    newGrid->cum_nabs_HeI = NULL;
    newGrid->cum_nabs_HeII = NULL;
    newGrid->frac_Q_HeI = NULL;
    newGrid->frac_Q_HeII = NULL;
    
    newGrid->XHeII = NULL;
    newGrid->XHeIII = NULL;
    newGrid->nrec_HeI = NULL;
    newGrid->nrec_HeII = NULL;
    
    //domain decomposition
    newGrid->local_n0 = 0;
    newGrid->local_0_start = 0;
    
    return newGrid;
}

void read_files_to_grid(grid_t *thisGrid, confObj_t thisInput)
{
    int solve_He = thisInput->solve_He;
#ifdef __MPI
    ptrdiff_t alloc_local, local_n0, local_0_start;
#else
    ptrdiff_t local_n0;
#endif
    int nbins;
    
    nbins = thisInput->grid_size;
    local_n0 = nbins;
    
    thisGrid->nbins = nbins;
    thisGrid->box_size = thisInput->box_size;
    thisGrid->lin_scales = thisInput->lin_scales;
    thisGrid->inc_log_scales = thisInput->inc_log_scales;
    thisGrid->max_scale = thisInput->max_scale;
    
    thisGrid->local_n0 = nbins;
    thisGrid->local_0_start = 0;
#ifdef __MPI
    fftw_mpi_init();
    
    alloc_local = fftw_mpi_local_size_3d(nbins, nbins, nbins, MPI_COMM_WORLD, &local_n0, &local_0_start);
    
    thisGrid->local_n0 = local_n0;
    thisGrid->local_0_start = local_0_start;
    
    thisGrid->igm_density = fftw_alloc_complex(alloc_local);
    thisGrid->igm_clump = fftw_alloc_complex(alloc_local);
    
    thisGrid->nion = fftw_alloc_complex(alloc_local);
    thisGrid->cum_nion = fftw_alloc_complex(alloc_local);
    thisGrid->cum_nrec = fftw_alloc_complex(alloc_local);
    thisGrid->cum_nabs = fftw_alloc_complex(alloc_local);
    thisGrid->frac_Q = fftw_alloc_complex(alloc_local);
    
    thisGrid->XHII = fftw_alloc_complex(alloc_local);
    thisGrid->nrec = fftw_alloc_complex(alloc_local);
    thisGrid->photHI = fftw_alloc_complex(alloc_local);
    
    if(solve_He == 1){
        thisGrid->nion_HeI = fftw_alloc_complex(alloc_local);
        thisGrid->nion_HeII = fftw_alloc_complex(alloc_local);
        thisGrid->cum_nion_HeI = fftw_alloc_complex(alloc_local);
        thisGrid->cum_nion_HeII = fftw_alloc_complex(alloc_local);
        thisGrid->cum_nrec_HeI = fftw_alloc_complex(alloc_local);
        thisGrid->cum_nrec_HeII = fftw_alloc_complex(alloc_local);
        thisGrid->cum_nabs_HeI = fftw_alloc_complex(alloc_local);
        thisGrid->cum_nabs_HeII = fftw_alloc_complex(alloc_local);
        thisGrid->frac_Q_HeI = fftw_alloc_complex(alloc_local);
        thisGrid->frac_Q_HeII = fftw_alloc_complex(alloc_local);
        thisGrid->XHeII = fftw_alloc_complex(alloc_local);
        thisGrid->XHeIII = fftw_alloc_complex(alloc_local);
        thisGrid->nrec_HeI = fftw_alloc_complex(alloc_local);
        thisGrid->nrec_HeII = fftw_alloc_complex(alloc_local);
    }
#else
    thisGrid->igm_density = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
    thisGrid->igm_clump = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
    
    thisGrid->nion = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
    thisGrid->cum_nion = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
    thisGrid->cum_nrec = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
    thisGrid->cum_nabs = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
    thisGrid->frac_Q = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
    
    thisGrid->XHII = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
    thisGrid->nrec = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
    thisGrid->photHI = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
    
    if(solve_He == 1){
        thisGrid->nion_HeI = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
        thisGrid->nion_HeII = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
        thisGrid->cum_nion_HeI = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
        thisGrid->cum_nion_HeII = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
        thisGrid->cum_nrec_HeI = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
        thisGrid->cum_nrec_HeII = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
        thisGrid->cum_nabs_HeI = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
        thisGrid->cum_nabs_HeII = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
        thisGrid->frac_Q_HeI = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
        thisGrid->frac_Q_HeII = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
        thisGrid->XHeII = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
        thisGrid->XHeIII = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
        thisGrid->nrec_HeI = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
        thisGrid->nrec_HeII = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
    }

#endif

    initialize_grid(thisGrid->igm_density, nbins, local_n0, 0.);
    initialize_grid(thisGrid->igm_clump, nbins, local_n0, 1.);
    
    initialize_grid(thisGrid->nion, nbins, local_n0, 0.);
    initialize_grid(thisGrid->cum_nion, nbins, local_n0, 0.);
    initialize_grid(thisGrid->cum_nrec, nbins, local_n0, 0.);
    initialize_grid(thisGrid->cum_nabs, nbins, local_n0, 0.);
    initialize_grid(thisGrid->frac_Q, nbins, local_n0, 0.);
    
    initialize_grid(thisGrid->XHII, nbins, local_n0, 0.);
    initialize_grid(thisGrid->nrec, nbins, local_n0, 0.);
    initialize_grid(thisGrid->photHI, nbins, local_n0, 0.);
    
    if(solve_He == 1){
        initialize_grid(thisGrid->nion_HeI, nbins, local_n0, 0.);
        initialize_grid(thisGrid->nion_HeII, nbins, local_n0, 0.);
        initialize_grid(thisGrid->cum_nion_HeI, nbins, local_n0, 0.);
        initialize_grid(thisGrid->cum_nion_HeII, nbins, local_n0, 0.);
        initialize_grid(thisGrid->cum_nrec_HeI, nbins, local_n0, 0.);
        initialize_grid(thisGrid->cum_nrec_HeII, nbins, local_n0, 0.);
        initialize_grid(thisGrid->cum_nabs_HeI, nbins, local_n0, 0.);
        initialize_grid(thisGrid->cum_nabs_HeII, nbins, local_n0, 0.);
        initialize_grid(thisGrid->frac_Q_HeI, nbins, local_n0, 0.);
        initialize_grid(thisGrid->frac_Q_HeII, nbins, local_n0, 0.);
        
        initialize_grid(thisGrid->XHeII, nbins, local_n0, 0.);
        initialize_grid(thisGrid->XHeIII, nbins, local_n0, 0.);
        initialize_grid(thisGrid->nrec_HeI, nbins, local_n0, 0.);
        initialize_grid(thisGrid->nrec_HeII, nbins, local_n0, 0.);
    }
    
#ifdef __MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}

void read_array(fftw_complex *toThisArray, grid_t *thisGrid, char *filename, int double_precision)
{
#ifdef __MPI
    ptrdiff_t local_n0, local_0_start;
#else
    ptrdiff_t local_n0;
#endif
    int nbins;
    
    nbins = thisGrid->nbins;
    local_n0 = thisGrid->local_n0;
    
    if(double_precision == 1)
    {
#ifdef __MPI
    local_0_start = thisGrid->local_0_start;
    read_grid_doubleprecision(toThisArray, nbins, local_n0, local_0_start, filename);
#else
    read_grid_doubleprecision(toThisArray, nbins, local_n0, filename);
#endif
    }else{
#ifdef __MPI
    local_0_start = thisGrid->local_0_start;
    read_grid(toThisArray, nbins, local_n0, local_0_start, filename);
#else
    read_grid(toThisArray, nbins, local_n0, filename);
#endif
    }
}

#ifdef __MPI
void read_grid(fftw_complex *toThisArray, int nbins, int local_n0, int local_0_start, char *filename)
#else
void read_grid(fftw_complex *toThisArray, int nbins, int local_n0, char *filename)
#endif
{
    float *tmparray;
        
    tmparray = (float*)malloc(sizeof(float)*local_n0*nbins*nbins);
    if(tmparray == NULL)
    {
        fprintf(stderr, "tmparray in read_grid (grid.c) could not be allocated\n");
        exit(EXIT_FAILURE);
    }
#ifdef __MPI
    int success;
    int resultlen;
    char msg[MPI_MAX_ERROR_STRING];

    MPI_File mpifile;
    MPI_Offset offset;
    MPI_Status status;
    
    offset = (local_0_start*nbins*nbins*sizeof(float));
    
    success = MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY,MPI_INFO_NULL, &mpifile);
    if(success != MPI_SUCCESS)
    {
        MPI_Error_string(success, msg, &resultlen);
        fprintf(stderr, "MPI_File_open(): %s\n", msg);
        exit(-1);
    }
    MPI_File_read_at_all(mpifile,offset,tmparray, local_n0*nbins*nbins,MPI_FLOAT,&status);
    MPI_File_close(&mpifile);
#else
    FILE *fp;
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


#ifdef __MPI
void read_grid_doubleprecision(fftw_complex *toThisArray, int nbins, int local_n0, int local_0_start, char *filename)
#else
void read_grid_doubleprecision(fftw_complex *toThisArray, int nbins, int local_n0, char *filename)
#endif
{
    double *tmparray;
    
    tmparray = (double*)malloc(sizeof(double)*local_n0*nbins*nbins);
    if(tmparray == NULL)
    {
        fprintf(stderr, "tmparray in read_grid (grid.c) could not be allocated\n");
        exit(EXIT_FAILURE);
    }
#ifdef __MPI
    int success;
    int resultlen;
    char msg[MPI_MAX_ERROR_STRING];

    MPI_File mpifile;
    MPI_Offset offset;
    MPI_Status status;
    
    offset = (local_0_start*nbins*nbins*sizeof(double));
    
    success = MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY,MPI_INFO_NULL, &mpifile);
    if(success != MPI_SUCCESS)
    {
        MPI_Error_string(success, msg, &resultlen);
        fprintf(stderr, "MPI_File_open(): %s\n", msg);
        exit(-1);
    }
    MPI_File_read_at_all(mpifile,offset,tmparray, local_n0*nbins*nbins,MPI_DOUBLE,&status);
    MPI_File_close(&mpifile);
#else
    FILE *fp;
    fp = fopen(filename, "rb");
    fread(tmparray, sizeof(double), nbins*nbins*nbins, fp);
    fclose(fp);
#endif
    
    for(int i=0; i<local_n0; i++)
    {
        for(int j=0; j<nbins; j++)
        {
            for(int k=0; k<nbins; k++)
            {
                toThisArray[i*nbins*nbins+j*nbins+k] = tmparray[i*nbins*nbins+j*nbins+k]+0.*I;
            }
        }
    }
    free(tmparray);
}

void initialize_grid(fftw_complex *thisArray, int nbins, int local_n0, double value)
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

void deallocate_grid(grid_t *thisGrid, confObj_t thisInput)
{
    int solve_He = thisInput->solve_He;
    
    if(thisGrid->igm_density != NULL) fftw_free(thisGrid->igm_density);
    if(thisGrid->igm_clump != NULL) fftw_free(thisGrid->igm_clump);
    
    if(thisGrid->nion != NULL) fftw_free(thisGrid->nion);
    if(thisGrid->cum_nion != NULL) fftw_free(thisGrid->cum_nion);
    if(thisGrid->cum_nrec != NULL) fftw_free(thisGrid->cum_nrec);
    if(thisGrid->cum_nabs != NULL) fftw_free(thisGrid->cum_nabs);
    if(thisGrid->frac_Q != NULL) fftw_free(thisGrid->frac_Q);
    
    if(thisGrid->XHII != NULL) fftw_free(thisGrid->XHII);
    if(thisGrid->nrec != NULL) fftw_free(thisGrid->nrec);
    if(thisGrid->photHI != NULL) fftw_free(thisGrid->photHI);
    
    if(solve_He == 1)
    {
        if(thisGrid->nion_HeI != NULL) fftw_free(thisGrid->nion_HeI);
        if(thisGrid->nion_HeII != NULL) fftw_free(thisGrid->nion_HeII);
        if(thisGrid->cum_nion_HeI != NULL) fftw_free(thisGrid->cum_nion_HeI);
        if(thisGrid->cum_nion_HeII != NULL) fftw_free(thisGrid->cum_nion_HeII);
        if(thisGrid->cum_nrec_HeI != NULL) fftw_free(thisGrid->cum_nrec_HeI);
        if(thisGrid->cum_nrec_HeII != NULL) fftw_free(thisGrid->cum_nrec_HeII);
        if(thisGrid->cum_nabs_HeI != NULL) fftw_free(thisGrid->cum_nabs_HeI);
        if(thisGrid->cum_nabs_HeII != NULL) fftw_free(thisGrid->cum_nabs_HeII);
        if(thisGrid->frac_Q_HeI != NULL) fftw_free(thisGrid->frac_Q_HeI);
        if(thisGrid->frac_Q_HeII != NULL) fftw_free(thisGrid->frac_Q_HeII);
        if(thisGrid->XHeII != NULL) fftw_free(thisGrid->XHeII);
        if(thisGrid->XHeIII != NULL) fftw_free(thisGrid->XHeIII);
        if(thisGrid->nrec_HeI != NULL) fftw_free(thisGrid->nrec_HeI);
        if(thisGrid->nrec_HeII != NULL) fftw_free(thisGrid->nrec_HeII);
    }
    
    free(thisGrid);
}

#ifdef __MPI
void write_grid_to_file_float(fftw_complex *thisArray, int nbins, int local_n0, int local_0_start, char *filename)
#else
void write_grid_to_file_float(fftw_complex *thisArray, int nbins, int local_n0, char *filename)
#endif
{
    float *tmparray;
    
    tmparray = (float*)malloc(sizeof(float)*local_n0*nbins*nbins);
    if(tmparray == NULL)
    {
        fprintf(stderr, "tmparray in write_grid_to_file_float (grid.c) could not be allocated\n");
        exit(EXIT_FAILURE);
    }
    
    for(int i=0; i<local_n0; i++)
    {
        for(int j=0; j<nbins; j++)
        {
            for(int k=0; k<nbins; k++)
            {
                tmparray[i*nbins*nbins+j*nbins+k] = (float)creal(thisArray[i*nbins*nbins+j*nbins+k]);
            }
        }
    }
        
#ifdef __MPI
    MPI_File mpifile;
    MPI_Offset offset;
    MPI_Status status;
    
    offset = (local_0_start*nbins*nbins*sizeof(float));

    MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &mpifile);
    if (mpifile == MPI_FILE_NULL)
    {
        fprintf(stderr, "Could not open file %s\n", filename);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);     
    }
    MPI_File_write_at_all(mpifile, offset, tmparray, local_n0*nbins*nbins, MPI_FLOAT, &status);
    MPI_File_close(&mpifile);
#else
    FILE * fp;
    
    fp = fopen(filename, "wb");
    if (fp == NULL)
    {
        fprintf(stderr, "Could not open file %s\n", filename);
        exit(EXIT_FAILURE);
    } 
    fwrite(tmparray, sizeof(float), nbins*nbins*nbins, fp);
    fclose(fp);
#endif
    
    free(tmparray);
}

#ifdef __MPI
void write_grid_to_file_double(fftw_complex *thisArray, int nbins, int local_n0, int local_0_start, char *filename)
#else
void write_grid_to_file_double(fftw_complex *thisArray, int nbins, int local_n0, char *filename)
#endif
{
    double *tmparray;
    
    tmparray = (double*)malloc(sizeof(double)*local_n0*nbins*nbins);
    if(tmparray == NULL)
    {
        fprintf(stderr, "tmparray in write_grid_to_file_double (grid.c) could not be allocated\n");
        exit(EXIT_FAILURE);
    }
    
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
    if (mpifile == MPI_FILE_NULL)
    {
        fprintf(stderr, "Could not open file %s\n", filename);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE); 
    }
    MPI_File_write_at_all(mpifile, offset, tmparray, local_n0*nbins*nbins, MPI_DOUBLE, &status);
    MPI_File_close(&mpifile);
#else
    FILE * fp;
    
    fp = fopen(filename, "wb");
    if (fp == NULL)
    {
        fprintf(stderr, "Could not open file %s\n", filename);
        exit(EXIT_FAILURE);
    } 
    fwrite(tmparray, sizeof(double), nbins*nbins*nbins, fp);
    fclose(fp);
#endif
    
    free(tmparray);
}

void save_to_file(fftw_complex *thisArray, grid_t *thisGrid, char *filename)
{
#ifdef __MPI
    write_grid_to_file_double(thisArray, thisGrid->nbins, thisGrid->local_n0, thisGrid->local_0_start, filename);
#else
    write_grid_to_file_double(thisArray, thisGrid->nbins, thisGrid->local_n0, filename);
#endif
}

double get_mean_grid(fftw_complex *thisArray, int nbins, int local_n0)
{
    double sum = 0.;
#ifdef __MPI
    double sum_all = 0.;
#endif
    
    for(int i=0; i<local_n0; i++)
    {
        for(int j=0; j<nbins; j++)
        {
            for(int k=0; k<nbins; k++)
            {
                sum += creal(thisArray[i*nbins*nbins+j*nbins+k]);
            }
        }
    }
    
#ifdef __MPI
    MPI_Allreduce(&sum, &sum_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    sum = sum_all / (double)(nbins*nbins*nbins);
#else
    sum = sum / (double)(nbins*nbins*nbins);
#endif
    
    return sum;
}
