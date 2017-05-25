/*
 *  confObj.h
 *  uvff
 *
 *  Created by Adrian Partl on 4/15/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef CONFOBJ_H
#define CONFOBJ_H

/*--- Includes ----------------------------------------------------------*/
#include "parse_ini.h"
#include <stdint.h>
#include <stdbool.h>


/*--- ADT handle --------------------------------------------------------*/
typedef struct confObj_struct *confObj_t;


/*--- Implemention of main structure ------------------------------------*/
struct confObj_struct {
    //Input
    char           *sources_file;
    char           *spectra_file;
    char           *dir_spectra;
    double         factorNion;
    double         factorFesc;
 
    //Cosmology
    double         Y;
    
    //Photoionization
    double         source_slope_index;
    
    //Temperature
    double         temperature;
    
    //Output
    char           *sources_file_new;
    char           *sources_file_HeI;
    char           *sources_file_HeII;
};


/*--- Prototypes of exported functions ----------------------------------*/
extern confObj_t
readConfObj(char *fileName);

extern confObj_t
confObj_new(parse_ini_t ini);

extern void
confObj_del(confObj_t *config);


#endif
