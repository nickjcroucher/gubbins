// $Id: gammaDistribution.cpp 2862 2007-11-27 10:59:03Z itaymay $

 #include "definitions.h"
#include "gammaDistribution.h"
#include "gammaUtilities.h"
#include "logFile.h"
#include <cmath>


gammaDistribution::gammaDistribution(MDOUBLE alpha,int in_number_of_categories) :
	generalGammaDistribution(alpha,alpha,in_number_of_categories) {}

gammaDistribution::gammaDistribution(const gammaDistribution& other) :
	generalGammaDistribution(other) {}
	
void gammaDistribution::setAlpha(MDOUBLE in_alpha) {
	if (in_alpha == _alpha) return;
	setGammaParameters( categories(), in_alpha);
}

//this function builds the gamma distribution
void gammaDistribution::setGammaParameters(int in_number_of_categories, MDOUBLE in_alpha) {
	generalGammaDistribution::setGammaParameters(in_number_of_categories,in_alpha,in_alpha);
}

void gammaDistribution::change_number_of_categories(int in_number_of_categories) {
	if (in_number_of_categories == categories())
		return;
	setGammaParameters( in_number_of_categories, _alpha, _alpha);
}

void gammaDistribution::setGammaParameters(int numOfCategories ,MDOUBLE alpha, MDOUBLE beta) {
	if (alpha!=beta)
		errorMsg::reportError("gammaDistribution::setGammaParameters : can not set beta because alpha must be equal to beta");
	generalGammaDistribution::setGammaParameters(numOfCategories,alpha,beta);
}
