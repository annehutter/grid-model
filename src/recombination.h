#ifndef RECOMBINATION
#define RECOMBINATION
#endif

// typedef struct
// {
// 	double dens;
// 	double photHI;
// 	double temp;
// 	double redshift;
// 	confObj_t simParam;
// 	
// 	double amplitude;
// 	double constant;
// 	double beta;
// } nrec_params_t;

typedef struct
{
	double constant;
	double factor;
} dens_integrand_t;

typedef struct
{
	double z;
	double dens_cell;
	confObj_t simParam;
} redshift_integrand_t;

typedef struct
{
	double constant_min;
	double constant_max;
	double dconstant;
	
	double factor_min;
	double factor_max;
	double dfactor;
} dens_table_t;

typedef struct
{
	double dens_cell_min;
	double dens_cell_max;
	double ddens_cell;
	
	double zmin;
	double zmax;
	double dz;
} redshift_table_t;

typedef struct
{
	double zmin;
	double zmax;
	double dz;
	
	double fmin;
	double fmax;
	double df;
	
	double dcellmin;
	double dcellmax;
	double ddcell;
} integral_table_t;

// void compute_number_recombinations(grid_t *thisGrid, confObj_t simParam, char *filename, dens_table_t *thisDensTable, redshift_table_t *thisRedshiftTable);
// double get_nrec_history(confObj_t simParam, double *norm_pdf, dens_table_t *thisDensTable, double *dens_table, redshift_table_t *thisRedshiftTable, double *redshift_table, double dens, double photHI, double temp, double zstart, double redshift);
double get_Xe(double XHII, double XHeII, double XHeIII, double Y, int specie);

void compute_number_recombinations(grid_t *thisGrid, confObj_t simParam);
double get_nrec_history(confObj_t simParam, double dens, double clump, double photHI, double temp, double zstart, double redshift, double Xe);
double get_nrec_HeI_history(confObj_t simParam, double dens, double clump, double photHeI, double temp, double zstart, double redshift, double Xe);
double get_nrec_HeII_history(confObj_t simParam, double dens, double clump, double photHeII, double temp, double zstart, double redshift, double Xe);


void compute_number_recombinations_M2000(grid_t *thisGrid, confObj_t simParam, char *filename, const integral_table_t *thisIntegralTable);
double get_nrec_history_M2000(confObj_t simParam, const integral_table_t *thisIntegralTable, double *integral_table, double dens, double photHI, double temp, double zstart, double redshift);

void compute_number_recombinations_const(grid_t *thisGrid, confObj_t simParam, int specie);
double get_nrec_history_constantInTime(confObj_t simParam, double z, double zstart);
double get_nrec_HeI_history_constantInTime(confObj_t simParam, double z, double zstart);
double get_nrec_HeII_history_constantInTime(confObj_t simParam, double z, double zstart);

//------------------------------------------------------------------------------
// table for pdf
//------------------------------------------------------------------------------

double *read_table_norm_pdf(char *filename);
void compute_table_norm_pdf(double zmin, double zmax, double d, int rank, int size, char *filename);
double amplitude_norm_pdf(double z);
double constant_norm_pdf(double z);

//------------------------------------------------------------------------------
// table for integrals over density and redshift
//------------------------------------------------------------------------------

double *read_table_integrals(char *filename, const integral_table_t *thisIntegralTable);
integral_table_t * initIntegralTable(double zmin, double zmax, double dz, double fmin, double fmax, double df, double dcellmin, double dcellmax, double ddcell);

//------------------------------------------------------------------------------
// table to perform integral over density
//------------------------------------------------------------------------------

dens_table_t *initDensTable(double constant_min, double constant_max, double dconstant, double factor_min, double factor_max, double dfactor);
double dens_integrand(double x, void * p);
double calc_dens_integral(dens_integrand_t *params);
double *create_table_dens(dens_table_t *thisDensTable);
void compute_table_dens(double constant_min, double constant_max, double d1, double factor_min, double factor_max, double d2, int rank, int size, char *filename);

//------------------------------------------------------------------------------
// table to perform integral over redshift
//------------------------------------------------------------------------------

redshift_table_t *initRedshiftTable(double dens_cell_min, double dens_cell_max, double ddens_cell, double zmin, double zmax, double dz);
double redshift_integrand(double x, void * p);
double calc_redshift_integral(redshift_integrand_t *params, double d2);
double *create_table_redshift(redshift_table_t *thisRedshiftTable, confObj_t simParam);
void compute_table_redshift(double dens_cell_min, double dens_cell_max, double d1, double zmin, double zmax, double d2, confObj_t simParam, int rank, int size, char *filename);

//------------------------------------------------------------------------------
// write table
//------------------------------------------------------------------------------

#ifdef __MPI
void write_table(int num, int offset, double *array, char *filename);
#else
void write_table(int num, double *array, char *filename);
#endif
