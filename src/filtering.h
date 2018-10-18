#ifndef FILTERING_H
#define FILTERING_H

void determine_mfp_nion(fftw_complex *frac_Q_smooth, fftw_complex *nion_smooth, fftw_complex *mfp_nion, int nbins, ptrdiff_t local_n0, double scale);
void determine_mfp(fftw_complex *frac_Q_smooth, fftw_complex *mfp, int nbins, ptrdiff_t local_n0, double scale);
double determine_mean_mfp(fftw_complex *mfp, int nbins, ptrdiff_t local_n0);

void construct_tophat_filter(fftw_complex *filter, int nbins, ptrdiff_t local_0_start, ptrdiff_t local_n0, float smooth_scale);

void determine_ion_fractions(fftw_complex *nion_smooth, int nbins, ptrdiff_t local_n0, int smallest_scale);
void map_central_ionized_cell_to_sphere(fftw_complex *new_cum_nion_smooth, fftw_complex *cum_nion_smooth, fftw_complex *filter, grid_t *thisGrid);
void choose_ion_fraction(fftw_complex *nion_smooth, fftw_complex *XHII_tmp, grid_t *thisGrid);

void map_bubbles_to_nrec(fftw_complex *Xion_tmp, fftw_complex *nrec, grid_t *thisGrid);

void combine_bubble_and_web_model(fftw_complex *Xion_tmp, fftw_complex *Xion, grid_t *thisGrid);

void copy_grid_array(fftw_complex *Xion_tmp, fftw_complex *Xion, grid_t *thisGrid);

void update_web_model(grid_t *thisGrid, confObj_t simParam, photIonlist_t *photIonBgList, int thisRank);

void adapt_HeII_to_HeIII(grid_t *thisGrid);

void compute_ionization_field(confObj_t simParam, grid_t *thisGrid, photIonlist_t *photIonBgList, int specie, int thisRank);

#endif