/*
 *  confObj.h
 *  uvff
 *
 *  Created by 
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
  
    // Type
    char           *sim_type;
    int            calc_ion_history;
    int            num_snapshots;
    char           *redshift_file;
    double         redshift_prev_snap;
    double         redshift;
    double         evol_time;

    //Cosmology
    double         h;
    double         omega_b;
    double         omega_m;
    double         omega_l;
    double         sigma8;
    double         Y;
    
    //Input
    int            grid_size;
    double         box_size;
    int            input_doubleprecision;
    int            inputfiles_comoving;
    char           *igm_density_file;
    int            dens_in_overdensity;
    double         mean_density;
    int            default_mean_density;
    char           *igm_clump_file;
    char           *sources_file;
    char           *nion_file;
    int            padded_box;

    //BubbleModel
    double         lin_scales;
    double         inc_log_scales;
    double         max_scale;
    int            ionize_sphere;
        
    //Photoionization
    int            use_web_model;
    int            photHI_model;
    int            calc_mfp;
    double         photHI_bg;
    char           *photHI_bg_file;
    double         mfp;
    double         source_slope_index;

    //Recombinations
    int            calc_recomb;
    int            const_recomb;
    double         dnrec_dt;
    char           *recomb_table;
    double         zmin, zmax, dz;
    double         fmin, fmax, df;
    double         dcellmin, dcellmax, ddcell;
    
    //Helium
    int            solve_He;
    char           *sources_HeI_file;
    char           *nion_HeI_file;
    char           *sources_HeII_file;
    char           *nion_HeII_file;
    double         dnrec_HeI_dt;
    double         dnrec_HeII_dt;
  
    //Output
    char           *out_XHII_file;
    int            write_photHI_file;
    char           *out_photHI_file;
    char           *out_XHeII_file;
    char           *out_XHeIII_file;
    
    double         f;
    double         factor;
};


/*--- Prototypes of exported functions ----------------------------------*/
extern confObj_t
readConfObj(char *fileName);

extern confObj_t
confObj_new(parse_ini_t ini);

extern void
confObj_del(confObj_t *config);


#endif
