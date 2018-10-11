#ifndef MEAN_FREE_PATH_H
#define MEAN_FREE_PATH_H

double pdf_mfp(double x, void * p);
double calc_integral_mfp(pdf_params_t params, double upLim);
double calc_local_mfp(confObj_t simParam, double dens, double photHI, double temperature, double redshift);
void compute_mfp(grid_t *thisGrid, confObj_t simParam);
void compute_web_mfp(grid_t *thisGrid, confObj_t simParam);

#endif