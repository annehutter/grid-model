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
	//General
	char 			*igm_density_file;
	char 			*halo_density_file;
	char 			*igm_clump_file;
	
	int			grid_size;
	double			box_size;
	
	char 			*sources_file;
	double 			evol_time;
	double			redshift;
	
	int 			dens_in_overdensity;
	double			mean_density;
	
	int 			inputfiles_comoving;
	
	double			h;
	double			omega_b;
	double 			omega_m;
	double 			omega_l;
	double			sigma8;
	
	char			*out_XHII_file;
	
	int			use_web_model;
	double			photHI_bg;
	int			compute_photHIfield;
	double			mfp;
	char			*out_photHI_file;
};


/*--- Prototypes of exported functions ----------------------------------*/
extern confObj_t
readConfObj(char *fileName);

extern confObj_t
confObj_new(parse_ini_t ini);

extern void
confObj_del(confObj_t *config);


#endif
