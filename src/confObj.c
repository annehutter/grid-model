/*
 *  confObj.c
 *  uvff
 *
 *  Created by Adrian Partl on 4/15/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

/*--- Includes ----------------------------------------------------------*/
#include "confObj.h"
#include <stdlib.h>
#include <stdio.h>
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
    
    //reading mandatory stuff
    
    //General
    getFromIni(&(config->num_snapshots), parse_ini_get_int32,
               ini, "numSnapshots", "General");    
    getFromIni(&(config->redshift_file), parse_ini_get_string,
               ini, "redshiftFile", "General");
    getFromIni(&(config->redshift_prev_snap), parse_ini_get_double,
               ini, "redshift_prevSnapshot", "General");
    getFromIni(&(config->redshift), parse_ini_get_double,
               ini, "finalRedshift", "General");
    getFromIni(&(config->evol_time), parse_ini_get_double,
               ini, "evolutionTime", "General");
    
    getFromIni(&(config->lin_scales), parse_ini_get_double,
               ini, "size_linear_scale", "General");
    getFromIni(&(config->inc_log_scales), parse_ini_get_double,
               ini, "first_increment_in_logscale", "General");
    
    getFromIni(&(config->default_mean_density), parse_ini_get_int32,
               ini, "useDefaultMeanDensity", "General");
    getFromIni(&(config->use_web_model), parse_ini_get_int32,
               ini, "useWebModel", "General");
    getFromIni(&(config->calc_ion_history), parse_ini_get_int32,
               ini, "calcIonHistory", "General");
    getFromIni(&(config->photHI_model), parse_ini_get_int32,
               ini, "photHImodel", "General");
    getFromIni(&(config->calc_mfp), parse_ini_get_int32,
               ini, "calcMeanFreePath", "General");
    getFromIni(&(config->const_recomb), parse_ini_get_int32,
               ini, "constantRecombinations", "General");
    getFromIni(&(config->calc_recomb), parse_ini_get_int32,
               ini, "calcRecombinations", "General");
    getFromIni(&(config->solve_He), parse_ini_get_int32,
               ini, "solveForHelium", "General");
    
    getFromIni(&(config->padded_box), parse_ini_get_int32,
               ini, "paddedBox", "General");
    
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
    
    getFromIni(&(config->igm_clump_file), parse_ini_get_string,
               ini, "inputIgmClumpFile", "Input");
    
    
    getFromIni(&(config->sources_file), parse_ini_get_string,
               ini, "inputSourcesFile", "Input");
    getFromIni(&(config->nion_file), parse_ini_get_string,
               ini, "inputNionFile", "Input");

    //Output
    getFromIni(&(config->out_XHII_file), parse_ini_get_string,
               ini, "output_XHII_file", "Output");
    getFromIni(&(config->write_photHI_file), parse_ini_get_int32,
               ini, "write_photHI_file", "Output");
    getFromIni(&(config->out_photHI_file), parse_ini_get_string,
               ini, "output_photHI_file", "Output");

    
    //Cosmology
    getFromIni(&(config->h), parse_ini_get_double,
               ini, "h", "Cosmology");
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

    //Photoionization
    getFromIni(&(config->photHI_bg_file), parse_ini_get_string,
               ini, "photHI_bg_file", "Photoionization")
    getFromIni(&(config->photHI_bg), parse_ini_get_double,
               ini, "photHI_bg", "Photoionization");
    getFromIni(&(config->mfp), parse_ini_get_double,
               ini, "meanFreePathInIonizedMedium", "Photoionization");
    getFromIni(&(config->source_slope_index), parse_ini_get_double,
               ini, "sourceSlopeIndex", "Photoionization");
    
    //Recombinations
    getFromIni(&(config->dnrec_dt), parse_ini_get_double,
               ini, "dnrec_dt", "Recombinations");
    getFromIni(&(config->recomb_table), parse_ini_get_string,
               ini, "recombinationTable", "Recombinations");
    getFromIni(&(config->zmin), parse_ini_get_double,
               ini, "zmin", "Recombinations");
    getFromIni(&(config->zmax), parse_ini_get_double,
               ini, "zmax", "Recombinations");
    getFromIni(&(config->dz), parse_ini_get_double,
               ini, "dz", "Recombinations");
    getFromIni(&(config->fmin), parse_ini_get_double,
               ini, "fmin", "Recombinations");
    getFromIni(&(config->fmax), parse_ini_get_double,
               ini, "fmax", "Recombinations");
    getFromIni(&(config->df), parse_ini_get_double,
               ini, "df", "Recombinations");
    getFromIni(&(config->dcellmin), parse_ini_get_double,
               ini, "dcellmin", "Recombinations");
    getFromIni(&(config->dcellmax), parse_ini_get_double,
               ini, "dcellmax", "Recombinations");
    getFromIni(&(config->ddcell), parse_ini_get_double,
               ini, "ddcell", "Recombinations");
    
    getFromIni(&(config->read_nrec_file), parse_ini_get_int32,
               ini, "readNrecFile", "Recombinations");
    getFromIni(&(config->nrec_file), parse_ini_get_string,
               ini, "inputRecombFile", "Recombinations");
    getFromIni(&(config->output_nrec_file), parse_ini_get_string,
               ini, "outputRecombFile", "Recombinations");

    //Helium
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
    
    getFromIni(&(config->out_XHeII_file), parse_ini_get_string,
               ini, "output_XHeII_file", "Helium");
    getFromIni(&(config->out_XHeIII_file), parse_ini_get_string,
               ini, "output_XHeIII_file", "Helium");
    
    config->f = 0.69;
    config->factor = 0.69;
    
    return config;
}

extern void
confObj_del(confObj_t *config)
{
    assert(config != NULL);
    assert(*config != NULL);
    
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
    xfree((*config)->nrec_file);
    xfree((*config)->output_nrec_file);
 
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
