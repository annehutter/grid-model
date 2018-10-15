#include <stdio.h> 
#include <stdlib.h>
#include <complex.h>

#ifdef __MPI
#include <fftw3-mpi.h>
#include <mpi.h>
#else
#include <fftw3.h>
#endif

#include "utils.h"
#include "confObj.h"
#include "grid.h"
#include "restart.h"

#define MAXLEN 1024

int32_t save_restart_file(confObj_t simParam, grid_t *grid, int32_t cycle, int32_t snap, char *restartFile, int32_t myRank)
{
    FILE *outfile;
    char fname[MAXLEN];

    if (myRank == 0)
    { 
        snprintf(fname, MAXLEN, "%s_info", restartFile);
        outfile = fopen(fname, "wb");
        if (outfile == NULL)
        {
            fprintf(stderr, "Cannot open file %s\n", fname);
            return EXIT_FAILURE;
        }
        
        fwrite(&cycle, sizeof(int32_t), 1, outfile);
        fwrite(&snap, sizeof(int32_t), 1, outfile); 
        fwrite(&(simParam->grid_size), sizeof(int), 1, outfile);
        fwrite(&(simParam->num_snapshots), sizeof(int), 1, outfile);
        fwrite(&(simParam->input_doubleprecision), sizeof(int), 1, outfile);
        fwrite(&(simParam->photHI_model), sizeof(int), 1, outfile); 
        fwrite(&(simParam->solve_He), sizeof(int), 1, outfile);
        fwrite(&(simParam->mfp), sizeof(double), 1, outfile);
        fwrite(&(simParam->photHI_bg), sizeof(double), 1, outfile); 
        fwrite(&(grid->mean_mfp), sizeof(double), 1, outfile);
        fwrite(&(grid->mean_photHI), sizeof(double), 1, outfile);

        fclose(outfile);
    }

    snprintf(fname, MAXLEN, "%s_photHI", restartFile);
    save_to_file(grid->photHI, grid, fname);

    snprintf(fname, MAXLEN, "%s_XHII", restartFile);
    save_to_file(grid->XHII, grid, fname);

    snprintf(fname, MAXLEN, "%s_cum_nion", restartFile);
    save_to_file(grid->cum_nion, grid, fname);

    snprintf(fname, MAXLEN, "%s_cum_nrec", restartFile);
    save_to_file(grid->cum_nrec, grid, fname);

    snprintf(fname, MAXLEN, "%s_cum_nabs", restartFile);
    save_to_file(grid->cum_nabs, grid, fname);

    if (simParam->solve_He == 1)
    {
        snprintf(fname, MAXLEN, "%s_cum_nion_HeI", restartFile);
        save_to_file(grid->cum_nion_HeI, grid, fname);

        snprintf(fname, MAXLEN, "%s_cum_nion_HeII", restartFile);
        save_to_file(grid->cum_nion_HeII, grid, fname);

        snprintf(fname, MAXLEN, "%s_cum_nrec_HeI", restartFile);
        save_to_file(grid->cum_nrec_HeI, grid, fname);

        snprintf(fname, MAXLEN, "%s_cum_nrec_HeI", restartFile);
        save_to_file(grid->cum_nrec_HeI, grid, fname);

        snprintf(fname, MAXLEN, "%s_cum_nrec_HeII", restartFile);
        save_to_file(grid->cum_nrec_HeII, grid, fname);

        snprintf(fname, MAXLEN, "%s_cum_nabs_HeI", restartFile);
        save_to_file(grid->cum_nabs_HeI, grid, fname);

        snprintf(fname, MAXLEN, "%s_cum_nabs_HeII", restartFile);
        save_to_file(grid->cum_nabs_HeII, grid, fname);    
    }
    return EXIT_SUCCESS;
}

int32_t read_restart_file(confObj_t simParam, grid_t *grid, int32_t *start_cycle, int32_t *snap, char *restartFile) 
{
    FILE *infile;
    char fname[MAXLEN];
    int restart_gridsize, restart_numsnap, restart_doubleprecision, restart_photHImodel, restart_solveHe;

    snprintf(fname, MAXLEN, "%s_info", restartFile);
    infile = fopen(fname, "rb");
    if (infile == NULL)
    {
      fprintf(stderr, "Cannot open file %s\n", fname);
      return EXIT_FAILURE;
    }
    
    fread(start_cycle, sizeof(int32_t), 1, infile); // The cycle number in the restart file is the number it ended on.
    ++(*start_cycle); // So we want to start on the next one.

    fread(snap, sizeof(int32_t), 1, infile); // Since the snapshot doesn't necessarily have to change on each iteration, don't increment it here.
    fread(&restart_gridsize, sizeof(int), 1, infile);
    fread(&restart_numsnap, sizeof(int), 1, infile);
    fread(&restart_doubleprecision, sizeof(int), 1, infile);
    fread(&(restart_photHImodel), sizeof(int), 1, infile); 
    fread(&(restart_solveHe), sizeof(int), 1, infile);
    fread(&(simParam->mfp), sizeof(double), 1, infile);
    fread(&(simParam->photHI_bg), sizeof(double), 1, infile);
    fread(&(grid->mean_mfp), sizeof(double), 1, infile);
    fread(&(grid->mean_photHI), sizeof(double), 1, infile);

    fclose(infile);

    if (simParam->grid_size != restart_gridsize || simParam->num_snapshots != restart_numsnap || simParam->input_doubleprecision != restart_doubleprecision || simParam->photHI_model != restart_photHImodel || simParam->solve_He != restart_solveHe)
    {
        fprintf(stderr, "Some of the parameters read in from the restart file do not match those from the iniFile.\n\n");
        fprintf(stderr, "===============\n");
        fprintf(stderr, "Restart File Parameters\nGrid Size : %d\nNumber Snapshots : %d\nDouble Precision : %d\nPhotoionization Model : %d\nSolve He : %d\n", restart_gridsize, restart_numsnap, restart_doubleprecision, restart_photHImodel, restart_solveHe); 
        fprintf(stderr, "===============\n");
        fprintf(stderr, "iniFile Parameters\nGrid Size : %d\nNumber Snapshots : %d\nDouble Precision : %d\nPhotoionization Model : %d\nSolve He : %d\n", simParam->grid_size, simParam->num_snapshots, simParam->input_doubleprecision, simParam->photHI_model, simParam->solve_He); 
        fprintf(stderr, "===============\n");
        fprintf(stderr, "The iniFile must be unchanged from the one used to create the restart file.\n");
        return EXIT_FAILURE; 
    }

    snprintf(fname, MAXLEN, "%s_photHI", restartFile);
    read_array(grid->photHI, grid, fname, simParam->input_doubleprecision);

    snprintf(fname, MAXLEN, "%s_XHII", restartFile);
    read_array(grid->XHII, grid, fname, simParam->input_doubleprecision);

    snprintf(fname, MAXLEN, "%s_cum_nion", restartFile);
    read_array(grid->cum_nion, grid, fname, simParam->input_doubleprecision);

    snprintf(fname, MAXLEN, "%s_cum_nrec", restartFile);
    read_array(grid->cum_nrec, grid, fname, simParam->input_doubleprecision);

    snprintf(fname, MAXLEN, "%s_cum_nabs", restartFile);
    read_array(grid->cum_nabs, grid, fname, simParam->input_doubleprecision);

    if (simParam->solve_He == 1)
    {
        snprintf(fname, MAXLEN, "%s_cum_nion_HeI", restartFile);
        read_array(grid->cum_nion_HeI, grid, fname, simParam->input_doubleprecision);

        snprintf(fname, MAXLEN, "%s_cum_nion_HeII", restartFile);
        read_array(grid->cum_nion_HeII, grid, fname, simParam->input_doubleprecision);

        snprintf(fname, MAXLEN, "%s_cum_nrec_HeI", restartFile);
        read_array(grid->cum_nrec_HeI, grid, fname, simParam->input_doubleprecision);

        snprintf(fname, MAXLEN, "%s_cum_nrec_HeI", restartFile);
        read_array(grid->cum_nrec_HeI, grid, fname, simParam->input_doubleprecision);

        snprintf(fname, MAXLEN, "%s_cum_nrec_HeII", restartFile);
        read_array(grid->cum_nrec_HeII, grid, fname, simParam->input_doubleprecision);

        snprintf(fname, MAXLEN, "%s_cum_nabs_HeI", restartFile);
        read_array(grid->cum_nabs_HeI, grid, fname, simParam->input_doubleprecision);

        snprintf(fname, MAXLEN, "%s_cum_nabs_HeII", restartFile);
        read_array(grid->cum_nabs_HeII, grid, fname, simParam->input_doubleprecision);    
    }

    return EXIT_SUCCESS;
}