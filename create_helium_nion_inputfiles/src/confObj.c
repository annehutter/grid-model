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
    
    //Input
    getFromIni(&(config->sources_file), parse_ini_get_string,
               ini, "inputSourcesFile", "Input");
    getFromIni(&(config->spectra_file), parse_ini_get_string,
               ini, "inputSpectraFile", "Input");
    getFromIni(&(config->dir_spectra), parse_ini_get_string,
               ini, "inputDirSpectra", "Input");
    getFromIni(&(config->factorNion), parse_ini_get_double,
               ini, "factorNion", "Input");
    getFromIni(&(config->factorFesc), parse_ini_get_double,
               ini, "factorFesc", "Input");
    
    //Cosmology
    getFromIni(&(config->Y), parse_ini_get_double,
               ini, "Y", "Cosmology");

    //Photoionization
    getFromIni(&(config->source_slope_index), parse_ini_get_double,
               ini, "sourceSlopeIndex", "Photoionization");
    
    //Temperature
    getFromIni(&(config->temperature), parse_ini_get_double,
               ini, "temperature", "Temperature");    
    
    //Output
    getFromIni(&(config->sources_file_new), parse_ini_get_string,
               ini, "sourcesFileNew", "Output");
    getFromIni(&(config->sources_file_HeI), parse_ini_get_string,
               ini, "sourcesFileHeI", "Output");
    getFromIni(&(config->sources_file_HeII), parse_ini_get_string,
               ini, "sourcesFileHeII", "Output");
    
    return config;
}

extern void
confObj_del(confObj_t *config)
{
    assert(config != NULL);
    assert(*config != NULL);
    
    xfree((*config)->sources_file);
    xfree((*config)->spectra_file);
    xfree((*config)->dir_spectra);

    xfree((*config)->sources_file_new);
    xfree((*config)->sources_file_HeI);
    xfree((*config)->sources_file_HeII);
    
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
