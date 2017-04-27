#include <gsl/gsl_math.h>    //included because M_PI is not defined in <math.h>
#include "phys_const.h"

const double ev_to_erg = 1.6e-12;
const double planck_cgs = 6.626e-27; 		// Planck constant [erg/Hz]
const double boltzman_cgs = 1.38e-16;		// Boltzman constant [erg/K]
const double clight_cm = 3.0e10;
const double G = 6.67408e-8;		//gravitational constant in cm^3/(gramm s^2)

const double Mpc_cm = 3.0857e24;		//Mpc [cm]
const double Myr_s = 1.e6*365.*24.*60.*60.;	//Myr [s]

const double me_g = 9.11e-28;			//mass of eletron [gramms]
const double mp_g = 1.673e-24;		//mass of proton [gramms]

const double crossSecThom_cm = 6.6524e-25;	//Thomson cross sections [cm^2]

const double sigma_HI = 3.2e-18;//6.3e-18;	//in cm^2	if spectrally averaged it is lower...

const double recomb_HII = 4.18e-13;
const double recomb_HII_B = 2.59e-13;	//at T = 1.e4 K

const double recomb_HeII = 4.31e-13;
const double recomb_HeII_B = 2.72e-13; 
const double recomb_HeII_1 = 1.54e-13;

const double recomb_HeIII = 2.09e-12;

const double nu_HI = 3.284033e15;
const double nu_HeI = 5.937821e15;
const double nu_HeII = 13.140960e15;

const double gamma_gas = 5./3.;
const double rho_g_cm = 1.8791e-29;	// in gramms/cm^3

