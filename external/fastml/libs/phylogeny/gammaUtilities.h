// $Id: gammaUtilities.h 4191 2008-06-12 19:03:36Z cohenofi $

 #ifndef ___GAMMA_UTILITIES
#define ___GAMMA_UTILITIES

#include "definitions.h"
#include "numRec.h" //fot the ITMAX	

/******************************************************************************
gamma utilities include calculating ln gamma and integral of gamma.
used mainly in building the gamma function and creating categories within it
******************************************************************************/

//gammln(xx): computes the ln of the Gamma function 
//the Gamma function is the integral from 0 to infinity of (t^(xx-1)*e^(-t)) dt.
MDOUBLE gammln(MDOUBLE xx);

//gammp(a, x): computes the incomplete Gamma function which is:
// 1/Gamma(a) * (the integral from 0 to x of (t^(a-1)*e^(-t)) dt)
//gammp can be computed in two different ways: by a series representation (gser(..)) 
//or by a continued fraction representation (gcf(..))
//gammp chooses to function will be used, according to the values of a and x 
MDOUBLE gammp(MDOUBLE a, MDOUBLE x);
void gser(MDOUBLE *gamser, MDOUBLE a, MDOUBLE x, MDOUBLE *gln);
void gcf(MDOUBLE *gammcf, MDOUBLE a, MDOUBLE x, MDOUBLE *gln);

MDOUBLE search_for_z_in_dis_with_any_beta(MDOUBLE alpha,MDOUBLE beta, MDOUBLE ahoson);
MDOUBLE search_for_z_in_dis_with_beta_1(MDOUBLE alpha, MDOUBLE ahoson);
MDOUBLE the_avarage_r_in_category_between_a_and_b(MDOUBLE a, MDOUBLE b, MDOUBLE alpha, MDOUBLE beta, int k);

//const int ITMAX = 100;
const MDOUBLE EPS = static_cast<MDOUBLE>(0.0000003);
const MDOUBLE FPMIN = static_cast<MDOUBLE>(1.0e-30);
const MDOUBLE ERR_FOR_GAMMA_CALC = static_cast<MDOUBLE>(0.00001);
const MDOUBLE MINIMUM_ALPHA_PARAM = static_cast<MDOUBLE>(0.05);	//was 0.05
const MDOUBLE MAXIMUM_ALPHA_PARAM = static_cast<MDOUBLE>(5.0);	
const MDOUBLE MINIMUM_BETA_PARAM = static_cast<MDOUBLE>(0.05);	//was 0.05
const MDOUBLE MAXIMUM_BETA_PARAM = static_cast<MDOUBLE>(5.0);	



//gammq(a, x) : computes 1 - the incomplete Gamma function (1-gammp(a,x)) which is:
//1/Gamma(a) * (the integral from infinite to x of (t^(a-1)*e^(-t)) dt).
//use for computing Chi-Square probability function (for the LRT):
//chiSquareProb(df,chiSquare) = gammq(df/2.0,chiSquare/2.0) 
MDOUBLE gammq(MDOUBLE a, MDOUBLE x);

#endif
