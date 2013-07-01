#include "gammaDistributionLaguerre.h"
#include "gammaUtilities.h"
#include "logFile.h"
#include <cmath>


gammaDistributionLaguerre::gammaDistributionLaguerre(MDOUBLE alpha,int in_number_of_categories) 
: generalGammaDistributionLaguerre(alpha,alpha,in_number_of_categories) 
{
}

gammaDistributionLaguerre::gammaDistributionLaguerre(const gammaDistributionLaguerre& other) 
: generalGammaDistributionLaguerre(other) 
{
}
	
void gammaDistributionLaguerre::setAlpha(MDOUBLE in_alpha) 
{
	if (in_alpha == _alpha) 
		return;
	setGammaParameters(categories(), in_alpha);
}

//this function builds the gamma distribution
void gammaDistributionLaguerre::setGammaParameters(int in_number_of_categories, MDOUBLE in_alpha) 
{
	generalGammaDistributionLaguerre::setGammaParameters(in_number_of_categories, in_alpha, in_alpha);
}

void gammaDistributionLaguerre::change_number_of_categories(int in_number_of_categories) 
{
	if (in_number_of_categories == categories())
		return;
	setGammaParameters(in_number_of_categories, _alpha, _alpha);
}

void gammaDistributionLaguerre::setGammaParameters(int numOfCategories ,MDOUBLE alpha, MDOUBLE beta) 
{
	if (alpha != beta)
		errorMsg::reportError("gammaDistributionLaguerre::setGammaParameters : can not set beta because alpha must be equal to beta");
	generalGammaDistributionLaguerre::setGammaParameters(numOfCategories, alpha, alpha);
}
