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
	getFromIni(&(config->igm_density_file), parse_ini_get_string,
	           ini, "inputIgmDensityFile", "General");
	getFromIni(&(config->halo_density_file), parse_ini_get_string,
	           ini, "inputHaloDensityFile", "General");
	getFromIni(&(config->igm_clump_file), parse_ini_get_string,
	           ini, "inputIgmClumpFile", "General");
	
	getFromIni(&(config->grid_size), parse_ini_get_int32,
		   ini, "gridsize", "General");
	getFromIni(&(config->box_size), parse_ini_get_double,
	           ini, "boxsize", "General");
	
	getFromIni(&(config->sources_file), parse_ini_get_string,
	           ini, "inputSourcesFile", "General");
	getFromIni(&(config->redshift), parse_ini_get_double,
		   ini, "redshift", "General");
	getFromIni(&(config->evol_time), parse_ini_get_double,
		   ini, "evolutionTime", "General");
	getFromIni(&(config->dens_in_overdensity), parse_ini_get_int32,
		   ini, "densityInOverdensity", "General");
	getFromIni(&(config->mean_density), parse_ini_get_double,
		   ini, "meanDensity", "General");
	getFromIni(&(config->default_mean_density), parse_ini_get_int32,
		   ini, "useDefaultMeanDensity", "General");
	
	getFromIni(&(config->inputfiles_comoving), parse_ini_get_int32,
		   ini, "inputFilesAreComoving", "General");
	getFromIni(&(config->h), parse_ini_get_double,
	           ini, "h", "General");
	getFromIni(&(config->omega_b), parse_ini_get_double,
	           ini, "omega_b", "General");
	getFromIni(&(config->omega_m), parse_ini_get_double,
	           ini, "omega_m", "General");
	getFromIni(&(config->omega_l), parse_ini_get_double,
	           ini, "omega_l", "General");
	getFromIni(&(config->sigma8), parse_ini_get_double,
	           ini, "sigma8", "General");
	getFromIni(&(config->Y), parse_ini_get_double,
		   ini, "Y", "General");
	
	getFromIni(&(config->out_XHII_file), parse_ini_get_string,
		   ini, "output_XHII_file", "General");
	
	getFromIni(&(config->use_web_model), parse_ini_get_int32,
		   ini, "useWebModel", "General");
	getFromIni(&(config->photHI_bg), parse_ini_get_double,
		   ini, "photionRateHI_bg", "General");
	getFromIni(&(config->compute_photHIfield), parse_ini_get_int32,
		   ini, "computePhotionHIfield_opticallyThin", "General");
	getFromIni(&(config->mfp), parse_ini_get_double,
	           ini, "meanFreePathInIonizedMedium", "General");
	getFromIni(&(config->out_photHI_file), parse_ini_get_string,
		   ini, "output_photHI_file", "General");

	getFromIni(&(config->generate_recomb_tables), parse_ini_get_int32,
		   ini, "generateRecombinationTables", "General");
	return config;
}

extern void
confObj_del(confObj_t *config)
{
	assert(config != NULL);
	assert(*config != NULL);
	
	xfree((*config)->igm_density_file);
	xfree((*config)->halo_density_file);
	xfree((*config)->igm_clump_file);
	xfree((*config)->sources_file);
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
