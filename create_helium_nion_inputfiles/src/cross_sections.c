#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>    //included because M_PI is not defined in <math.h>
#include "phys_const.h"
#include "cross_sections.h"

// cross sections

double crossSecHI(double nu){
// 	double eps;

	if(nu < nu_HI){
		return 0.;
	}
	else{
// 		eps = sqrt(nu/nu_HI - 1.);
// 		return  6.3e-18*pow(nu/nu_HI,-4)*exp(4.-4.*atan(eps)/eps)/(1.-exp(-2.*M_PI/eps));
		return  6.3e-18*pow(nu/nu_HI,-3);
	}
}

double crossSecHeI(double nu){

	if(nu < nu_HeI){
		return 0.;
	}
	else{
		return 7.2e-18*(1.66*pow(nu/nu_HeI,-2.05)+0.66*pow(nu/nu_HeI,-3.05));
	}
}

double crossSecHeII(double nu){
	double eps;
	if(nu < nu_HeII){
		return 0.;
	}
	else{
		eps = sqrt(nu/nu_HeII - 1.);
		return 1.58e-18*pow(nu/nu_HeII,-4)*exp(4.-4.*atan(eps)/eps)/(1. - exp(-2.*M_PI/eps));
	}
}
