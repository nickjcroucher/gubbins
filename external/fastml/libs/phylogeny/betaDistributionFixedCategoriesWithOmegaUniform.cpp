#include "betaDistributionFixedCategoriesWithOmegaUniform.h"
#include "errorMsg.h"
#include "gammaUtilities.h"
#include "matrixUtils.h"


betaDistributionFixedCategoriesOmegaUniform::betaDistributionFixedCategoriesOmegaUniform(const betaDistributionFixedCategoriesOmegaUniform& other)
: _betaDistr(other._betaDistr),_omegaDistr(other._omegaDistr){
	
}

betaDistributionFixedCategoriesOmegaUniform::betaDistributionFixedCategoriesOmegaUniform(int betaDistrCatNum,MDOUBLE alpha,MDOUBLE beta,
																						 int omegaCatNum,MDOUBLE omegaLowerBound,MDOUBLE omegaUpperBound)
{
	_betaDistr.setBetaParameters(betaDistrCatNum,alpha,beta);
	_omegaDistr.setGlobalRate(1.0);
	_omegaDistr.setUniformParameters(omegaCatNum,omegaLowerBound,omegaUpperBound);
	
}

void betaDistributionFixedCategoriesOmegaUniform::setBetaParameters(int in_number_of_categories, MDOUBLE alpha,  MDOUBLE beta)
{
	_betaDistr.setBetaParameters(in_number_of_categories,alpha,beta);
}



void betaDistributionFixedCategoriesOmegaUniform::change_number_of_categories(int in_number_of_categories)
{
	_betaDistr.change_number_of_categories(in_number_of_categories); 
}


const MDOUBLE betaDistributionFixedCategoriesOmegaUniform::ratesProb(const int i_rate) const {
	int noBetaDistCat = _betaDistr.categories();
	if (i_rate < _betaDistr.categories())
		return _betaDistr.ratesProb(i_rate);
	else return _omegaDistr.ratesProb(i_rate - noBetaDistCat); //omega prob
}


const MDOUBLE betaDistributionFixedCategoriesOmegaUniform::rates(const int i) const {
	int noBetaDistCat = _betaDistr.categories();
	if (i < noBetaDistCat)
		return _betaDistr.rates(i);
	else return _omegaDistr.rates(i - noBetaDistCat); //omega

}

const MDOUBLE betaDistributionFixedCategoriesOmegaUniform::getCumulativeProb(const MDOUBLE x) const {
	return _betaDistr.getCumulativeProb(x);
}