#ifndef CONVOLUTION_FFT_H
#define CONVOLUTION_FFT_H

void convolve_fft(grid_t *thisGrid, fftw_complex *filter, fftw_complex *output, fftw_complex *input);
void convolve_fft_ktophat(grid_t *thisGrid, fftw_complex *kfilter, fftw_complex *output, fftw_complex *input);

#endif