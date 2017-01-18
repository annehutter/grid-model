#ifndef FRACTION_Q_H
#define FRACTION_Q_H
#endif

void compute_cum_values(grid_t *thisGrid, confObj_t simParam, int specie);
void compute_Q(grid_t *thisGrid, fftw_complex *frac_Q, fftw_complex *nion, fftw_complex *nabs);
