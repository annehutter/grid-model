/*
 *  confObj.c
 *  uvff
 *
 *  Created by 
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

/*--- Includes ----------------------------------------------------------*/
#include "confObj.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <inttypes.h>
#include <stdbool.h>
#include "xmem.h"

/*--- Defines for the Ini structure -------------------------------------*/


/*--- Prototypes of local functions -------------------------------------*/


/*--- Implementations of exported functios ------------------------------*/
extern confObj_t
confObj_new(parse_ini_t ini)
{
    confObj_t config;
    assert(ini != NULL);
    
    config = xmalloc(sizeof(struct confObj_struct));
    
    char *photHImodel = NULL;
    char *recombModel = NULL;
    
    //reading mandatory stuff
    
    //Type
    getFromIni(&(config->sim_type), parse_ini_get_string,
               ini, "simulationType", "Type");
    
    if(strcmp(config->sim_type, "FIXED_REDSHIFT") == 0)
    {
        printf("FIXED_REDSHIFT\n");
        config->calc_ion_history = 0;
        config->num_snapshots = 1;
        config->redshift_file = NULL;
        config->redshift_prev_snap = config->redshift;
        getFromIni(&(config->redshift), parse_ini_get_double,
                  ini, "redshift", "FixedRedshift");
        getFromIni(&(config->evol_time), parse_ini_get_double,
                  ini, "evolutionTime", "FixedRedshift");
    }
    else if(strcmp(config->sim_type, "EVOLVE_REDSHIFT") == 0)
    {
        printf("EVOLVE_REDSHIFT\n");
        config->calc_ion_history = 1;
        getFromIni(&(config->num_snapshots), parse_ini_get_int32,
                  ini, "numSnapshots", "EvolveRedshift");  
        config->redshift_file = NULL;
        getFromIni(&(config->redshift_prev_snap), parse_ini_get_double,
                  ini, "redshift_start", "EvolveRedshift");
        getFromIni(&(config->redshift), parse_ini_get_double,
                  ini, "redshift_end", "EvolveRedshift");
        config->evol_time = 0.;
    }
    else if(strcmp(config->sim_type, "EVOLVE_ALL") == 0)
    {
        printf("EVOLVE_ALL\n");
        config->calc_ion_history = 1;
        getFromIni(&(config->num_snapshots), parse_ini_get_int32,
                  ini, "numSnapshots", "EvolveAll");    
        getFromIni(&(config->redshift_file), parse_ini_get_string,
                  ini, "redshiftFile", "EvolveAll");
        config->redshift = 0.;
        config->redshift_prev_snap = 0.;
        config->evol_time = 0.;
    }
    else
    {
        printf("This simulation type is not supported!\nAborting.\n");
        exit(EXIT_FAILURE);
    }
    

    //Cosmology
    getFromIni(&(config->h), parse_ini_get_double,
               ini, "hubble_h", "Cosmology");
    getFromIni(&(config->omega_b), parse_ini_get_double,
               ini, "omega_b", "Cosmology");
    getFromIni(&(config->omega_m), parse_ini_get_double,
               ini, "omega_m", "Cosmology");
    getFromIni(&(config->omega_l), parse_ini_get_double,
               ini, "omega_l", "Cosmology");
    getFromIni(&(config->sigma8), parse_ini_get_double,
               ini, "sigma8", "Cosmology");
    getFromIni(&(config->Y), parse_ini_get_double,
               ini, "Y", "Cosmology");
    
    
    //Input
    getFromIni(&(config->grid_size), parse_ini_get_int32,
               ini, "gridsize", "Input");
    getFromIni(&(config->box_size), parse_ini_get_double,
               ini, "boxsize", "Input");
    
    getFromIni(&(config->input_doubleprecision), parse_ini_get_int32,
               ini, "inputFilesAreInDoublePrecision", "Input");
    getFromIni(&(config->inputfiles_comoving), parse_ini_get_int32,
               ini, "inputFilesAreComoving", "Input");

    getFromIni(&(config->igm_density_file), parse_ini_get_string,
               ini, "inputIgmDensityFile", "Input");
    getFromIni(&(config->dens_in_overdensity), parse_ini_get_int32,
               ini, "densityInOverdensity", "Input");
    getFromIni(&(config->mean_density), parse_ini_get_double,
               ini, "meanDensity", "Input");
    getFromIni(&(config->default_mean_density), parse_ini_get_int32,
               ini, "useDefaultMeanDensity", "Input");
    
    getFromIni(&(config->igm_clump_file), parse_ini_get_string,
               ini, "inputIgmClumpFile", "Input");
    
    
    getFromIni(&(config->sources_file), parse_ini_get_string,
               ini, "inputSourcesFile", "Input");
    getFromIni(&(config->nion_file), parse_ini_get_string,
               ini, "inputNionFile", "Input");

    getFromIni(&(config->padded_box), parse_ini_get_int32,
               ini, "paddedBox", "Input");
    
    //Bubble Model
    getFromIni(&(config->lin_scales), parse_ini_get_double,
               ini, "size_linear_scale", "BubbleModel");
    getFromIni(&(config->inc_log_scales), parse_ini_get_double,
               ini, "first_increment_in_logscale", "BubbleModel");
    getFromIni(&(config->max_scale), parse_ini_get_double,
               ini, "max_scale", "BubbleModel");
    getFromIni(&(config->ionize_sphere), parse_ini_get_int32,
               ini, "useIonizeSphereModel", "BubbleModel");
    
    //Photoionization
    getFromIni(&(config->use_web_model), parse_ini_get_int32,
               ini, "useWebModel", "PhotoionizationModel");
    getFromIni(&photHImodel, parse_ini_get_string,
               ini, "photHImodel", "PhotoionizationModel");
    getFromIni(&(config->calc_mfp), parse_ini_get_int32,
               ini, "calcMeanFreePath", "PhotoionizationModel");
    
    if(strcmp(photHImodel, "PHOTHI_CONST") == 0)
    {
        config->photHI_model = 0;
        getFromIni(&(config->photHI_bg), parse_ini_get_double,
                  ini, "photHI_bg", "PhotoionizationConst");
        config->photHI_bg_file = NULL;
        config->mfp = 0.;
        config->source_slope_index = 0.;
    }
    else if(strcmp(photHImodel, "PHOTHI_GIVEN") == 0)
    {
        config->photHI_model = 11;
        config->photHI_bg = 0.;
        getFromIni(&(config->photHI_bg_file), parse_ini_get_string,
                  ini, "photHI_bg_file", "PhotoionizationGiven");
        getFromIni(&(config->mfp), parse_ini_get_double,
                  ini, "meanFreePathInIonizedMedium", "PhotoionizationGiven");
        getFromIni(&(config->source_slope_index), parse_ini_get_double,
                  ini, "sourceSlopeIndex", "PhotoionizationMfp");
    }
    else if(strcmp(photHImodel, "PHOTHI_FLUX") == 0)
    {
        config->photHI_model = 1;
        config->photHI_bg = 0.;
        config->photHI_bg_file = NULL;
        getFromIni(&(config->mfp), parse_ini_get_double,
                  ini, "meanFreePathInIonizedMedium", "PhotoionizationFlux");
        getFromIni(&(config->source_slope_index), parse_ini_get_double,
                  ini, "sourceSlopeIndex", "PhotoionizationMfp");
    }
    else if(strcmp(photHImodel, "PHOTHI_MFP") == 0)
    {
        config->photHI_model = 2;
        config->photHI_bg = 0.;
        config->photHI_bg_file = NULL;
        config->mfp = 0.;
        getFromIni(&(config->source_slope_index), parse_ini_get_double,
                  ini, "sourceSlopeIndex", "PhotoionizationMfp");
    }
    else
    {
        printf("\n***WARNING***\nThis is not a supported photoionization model. Choose one of the following options:\n'PHOTHI_CONST': constant photoionization rate\n'PHOT_GIVEN': mean value of the photoionization rate is given by the values in the photHI_bg file, distribution follows the flux of sources\n'PHOTHI_FLUX': photoionization rate is calculated based on the flux from sources\n'PHOTHI_MFP': photoionization rate is given bz the mean free path of the inoized regions.\n***WARNING***\n");
        exit(EXIT_FAILURE);
    }
    
    //Recombinations
    getFromIni(&(config->calc_recomb), parse_ini_get_int32,
               ini, "calcRecombinations", "RecombinationModel");
    getFromIni(&recombModel, parse_ini_get_string,
              ini, "recombinationModel", "RecombinationModel");
            
    if(strcmp(recombModel, "RECOMB_DEFAULT") == 0)
    {
        config->const_recomb = 0;
        config->dnrec_dt = 0.;
        config->recomb_table = NULL;
        config->zmin = 0.;
        config->zmax = 0.;
        config->dz  = 0.;
        config->fmin = 0.;
        config->fmax  = 0.;
        config->df = 0.;
        config->dcellmin = 0.;
        config->dcellmax = 0.;
        config->ddcell = 0.;
    }
    else if(strcmp(recombModel, "RECOMB_CONST") == 0)
    {
        config->const_recomb = 1;
        getFromIni(&(config->dnrec_dt), parse_ini_get_double,
                  ini, "dnrec_dt", "RecombinationConst");
        config->recomb_table = NULL;
        config->zmin = 0.;
        config->zmax = 0.;
        config->dz  = 0.;
        config->fmin = 0.;
        config->fmax  = 0.;
        config->df = 0.;
        config->dcellmin = 0.;
        config->dcellmax = 0.;
        config->ddcell = 0.;
    }
    else if(strcmp(recombModel, "RECOMB_TABLE") == 0)
    {
        config->const_recomb = 0;
        config->dnrec_dt = 0.;
        getFromIni(&(config->recomb_table), parse_ini_get_string,
                  ini, "recombinationTable", "RecombinationTable");
        getFromIni(&(config->zmin), parse_ini_get_double,
                  ini, "zmin", "RecombinationTable");
        getFromIni(&(config->zmax), parse_ini_get_double,
                  ini, "zmax", "RecombinationTable");
        getFromIni(&(config->dz), parse_ini_get_double,
                  ini, "dz", "RecombinationTable");
        getFromIni(&(config->fmin), parse_ini_get_double,
                  ini, "fmin", "RecombinationTable");
        getFromIni(&(config->fmax), parse_ini_get_double,
                  ini, "fmax", "RecombinationTable");
        getFromIni(&(config->df), parse_ini_get_double,
                  ini, "df", "RecombinationTable");
        getFromIni(&(config->dcellmin), parse_ini_get_double,
                  ini, "dcellmin", "RecombinationTable");
        getFromIni(&(config->dcellmax), parse_ini_get_double,
                  ini, "dcellmax", "RecombinationTable");
        getFromIni(&(config->ddcell), parse_ini_get_double,
                  ini, "ddcell", "RecombinationTable");
    }
    else
    {
        printf("\n***WARNING***\nThis is not a supported recombination model. Choose one of the following options:\n'RECOMB_DEFAULT': recombination rate depends on density\n'RECOMB_CONST': recombination rate is given\n'RECOMB_TABLE': not supported currently.\n***WARNING***\n");
        exit(EXIT_FAILURE);
    }

    //Helium
    getFromIni(&(config->solve_He), parse_ini_get_int32,
               ini, "solveForHelium", "Helium");
    getFromIni(&(config->sources_HeI_file), parse_ini_get_string,
               ini, "inputSourcesHeIFile", "Helium");
    getFromIni(&(config->nion_HeI_file), parse_ini_get_string,
               ini, "inputNionHeIFile", "Helium");
    getFromIni(&(config->sources_HeII_file), parse_ini_get_string,
               ini, "inputSourcesHeIIFile", "Helium");
    getFromIni(&(config->nion_HeII_file), parse_ini_get_string,
               ini, "inputNionHeIIFile", "Helium");
    
    getFromIni(&(config->dnrec_HeI_dt), parse_ini_get_double,
               ini, "dnrec_HeI_dt", "Helium");
    getFromIni(&(config->dnrec_HeII_dt), parse_ini_get_double,
               ini, "dnrec_HeII_dt", "Helium");
    

    
    
    //Output
    getFromIni(&(config->out_XHII_file), parse_ini_get_string,
               ini, "output_XHII_file", "Output");
    getFromIni(&(config->write_photHI_file), parse_ini_get_int32,
               ini, "write_photHI_file", "Output");
    getFromIni(&(config->out_photHI_file), parse_ini_get_string,
               ini, "output_photHI_file", "Output");
    getFromIni(&(config->out_XHeII_file), parse_ini_get_string,
               ini, "output_XHeII_file", "Output");
    getFromIni(&(config->out_XHeIII_file), parse_ini_get_string,
               ini, "output_XHeIII_file", "Output");
    
    config->f = 0.69;
    config->factor = 0.69;
    
    free(photHImodel);
    free(recombModel);
    
    return config;
}

extern void
confObj_del(confObj_t *config)
{
    assert(config != NULL);
    assert(*config != NULL);
    
    xfree((*config)->sim_type);
    
    //General
    xfree((*config)->redshift_file);
    
    //Input
    xfree((*config)->igm_density_file);
    xfree((*config)->igm_clump_file);
    xfree((*config)->sources_file);
    xfree((*config)->nion_file);
    
    //Output
    xfree((*config)->out_XHII_file);
    xfree((*config)->out_photHI_file);
    
    //Photoionization
    xfree((*config)->photHI_bg_file);
    
    //Recombinations
    xfree((*config)->recomb_table);
 
    //Helium
    xfree((*config)->sources_HeI_file);
    xfree((*config)->nion_HeI_file);
    xfree((*config)->sources_HeII_file);
    xfree((*config)->nion_HeII_file);
    xfree((*config)->out_XHeII_file);
    xfree((*config)->out_XHeIII_file);

    xfree(*config);
    *config = NULL;
}

extern confObj_t
readConfObj(char *fileName) {
    confObj_t newConfObj;
    parse_ini_t ini;
    
    ini = parse_ini_open(fileName);
    if (ini == NULL) {
        fprintf(stderr, "FATAL:  Could not open %s for reading.\n",
                fileName);
        exit(EXIT_FAILURE);
    }
    
    newConfObj = confObj_new(ini);
    
    parse_ini_close(&ini);
    
    return newConfObj;
}

/*--- Implementations of local functions --------------------------------*/
