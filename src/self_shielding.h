#ifndef SELF_SHIELDING_H
#define SELF_SHIELDING_H
#endif


double ss_calc_modPhotHI(double densH, double densSS);
double ss_calc_densSS(confObj_t simParam, double photHI, double temperature, double redshift);
double ss_calc_XHII(double dens, double photHI, double temp, double Y);
double calc_photHI_source(source_t *thisSource, double mfp_inv, double boxsize_Mpc, float x, float y, float z);
// void compute_photoionization_field(grid_t *thisGrid, sourcelist_t *thisSourcelist, confObj_t simParam);
// void construct_photHI_filter(fftw_complex *filter, grid_t *thisGrid, confObj_t simParam);
void convolve_fft_photHI(grid_t *thisGrid, fftw_complex *filter, fftw_complex *nion_smooth);
void compute_photHI(grid_t *thisGrid, confObj_t simParam);
void set_value_to_photoionization_field(grid_t *thisGrid, confObj_t simParam);
void set_value_to_photHI_bg(grid_t *thisGrid, confObj_t simParam, double value);
void compute_web_ionfraction(grid_t *thisGrid, confObj_t simParam);
