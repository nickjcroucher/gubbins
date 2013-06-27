#include "gammaDistributionFixedCategories.h"
#include "errorMsg.h"
#include "gammaUtilities.h"
#include "matrixUtils.h"

gammaDistributionFixedCategories::gammaDistributionFixedCategories(const Vdouble& fixedBoundaries, MDOUBLE alpha)
: generalGammaDistributionFixedCategories(fixedBoundaries,alpha,alpha) 
{

}

gammaDistributionFixedCategories::gammaDistributionFixedCategories(const gammaDistributionFixedCategories& other) 
: generalGammaDistributionFixedCategories(other) {
}

gammaDistributionFixedCategories::gammaDistributionFixedCategories(MDOUBLE alpha, int catNum)
: generalGammaDistributionFixedCategories(alpha, alpha,catNum) 
{
}

void gammaDistributionFixedCategories::setGammaParameters(int in_number_of_categories, MDOUBLE alpha)
{
	generalGammaDistributionFixedCategories::setGammaParameters(in_number_of_categories,alpha,alpha);
}


void gammaDistributionFixedCategories::setAlpha(MDOUBLE in_alpha) {
	if (in_alpha == _alpha) return;
	setGammaParameters( categories(), in_alpha);
}

void gammaDistributionFixedCategories::change_number_of_categories(int in_number_of_categories)
{
	generalGammaDistributionFixedCategories::change_number_of_categories(in_number_of_categories); 
}
