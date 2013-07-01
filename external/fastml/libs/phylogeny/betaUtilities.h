// 	$Id: betaUtilities.h 962 2006-11-07 15:13:34Z privmane $	
#ifndef ___BETA_UTILITIES
#define ___BETA_UTILITIES

#include "definitions.h"
#include "numRec.h"

/******************************************************************************
beta utilities include calculating inverse of the beta cdf and calculation of mean values
used mainly in building the gamma function and creating categories within it
******************************************************************************/

MDOUBLE inverseCDFBeta(MDOUBLE a, MDOUBLE b, MDOUBLE prob);
MDOUBLE computeAverage_r(MDOUBLE leftBound, MDOUBLE rightBound, MDOUBLE alpha, MDOUBLE beta, int k);
MDOUBLE incompleteBeta(MDOUBLE alpha, MDOUBLE beta, MDOUBLE x);
MDOUBLE betacf(MDOUBLE a, MDOUBLE b, MDOUBLE x);
MDOUBLE betaln(MDOUBLE alpha, MDOUBLE beta);



#endif
