#ifndef GRID_H
#define GRID_H
#endif

/* structure for grid */

typedef struct
{
	int nbins;
	float box_size;
	float lin_scales;
	float inc_log_scales;
	
	float xmin, ymin, zmin;
	
	fftw_complex *igm_density;
    fftw_complex *igm_clump;
    
    // hydrogen
	fftw_complex *nion;
	fftw_complex *cum_nion;
    fftw_complex *cum_nrec;
	fftw_complex *cum_nabs;
	fftw_complex *frac_Q;
	
	fftw_complex *XHII;
	fftw_complex *nrec;
	
	fftw_complex *photHI;
	double mean_photHI;
    
    //helium
    fftw_complex *nion_HeI;
    fftw_complex *nion_HeII;
    fftw_complex *cum_nion_HeI;
    fftw_complex *cum_nion_HeII;
    fftw_complex *cum_nrec_HeI;
    fftw_complex *cum_nrec_HeII;
    fftw_complex *cum_nabs_HeI;
    fftw_complex *cum_nabs_HeII;
    fftw_complex *frac_Q_HeI;
    fftw_complex *frac_Q_HeII;
    
    fftw_complex *XHeII;
    fftw_complex *XHeIII;
    fftw_complex *nrec_HeI;
    fftw_complex *nrec_HeII;
	
    //domain decomposition
	int local_n0;
	int local_0_start;
} grid_t;

/* functions */

grid_t *initGrid();
void read_files_to_grid(grid_t *thisGrid, confObj_t thisInput);
void read_array(fftw_complex *toThisArray, grid_t *thisGrid, char *filename, int double_precision);

#ifdef __MPI
void read_grid(fftw_complex *toThisArray, int nbins, int local_n0, int local_0_start, char *filename);
void read_grid_doubleprecision(fftw_complex *toThisArray, int nbins, int local_n0, int local_0_start, char *filename);
#else
void read_grid(fftw_complex *toThisArray, int nbins, int local_n0, char *filename);
void read_grid_doubleprecision(fftw_complex *toThisArray, int nbins, int local_n0, char *filename);
#endif

void initialize_grid(fftw_complex *thisArray, int nbins, int local_n0, double value);
void deallocate_grid(grid_t *thisGrid);

#ifdef __MPI
void write_grid_to_file_float(fftw_complex *thisArray, int nbins, int local_n0, int local_0_start, char *filename);
void write_grid_to_file_double(fftw_complex *thisArray, int nbins, int local_n0, int local_0_start, char *filename);
#else
void write_grid_to_file_float(fftw_complex *thisArray, int nbins, int local_n0, char *filename);
void write_grid_to_file_double(fftw_complex *thisArray, int nbins, int local_n0, char *filename);
#endif

void save_to_file(fftw_complex *thisArray, grid_t *thisGrid, char *filename);
