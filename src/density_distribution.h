#ifndef MEAN_FREE_PATH
#define MEAN_FREE_PATH

typedef struct
{
	double redshift;
	double amplitude;
	double constant;
	double beta;
} pdf_params_t;

double Hubble(confObj_t simParam);
double dd_pdf(double x, void * p);
double dd_calc_integral(pdf_params_t params, double upLim);
double dd_pdf_mass(double x, void * p);
double dd_calc_integral_mass(pdf_params_t params, double upLim);
void dd_set_norm_pdf(pdf_params_t * params, double redshift);
double dd_frac_densSS(double densSS, confObj_t simParam);
double dd_calc_mfp(confObj_t simParam, double photHI, double temperature, double redshift);
void set_mfp_Miralda2000(confObj_t simParam);

#endif