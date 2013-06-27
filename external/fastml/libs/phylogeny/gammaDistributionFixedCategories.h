#ifndef ___GAMMA_DISTR_FIXED_CATEGORIES
#define ___GAMMA_DISTR_FIXED_CATEGORIES
/************************************************************
This class differ from the regular GammaDistribution in that 
the rateCategories are fixed according to the user's decision. 
Thus, only the probability of each category changes for each specific alpha value but 
the rate categories themselves are constant.
************************************************************/
#include "definitions.h"
#include "generalGammaDistributionFixedCategories.h"
#include "errorMsg.h"

class gammaDistributionFixedCategories : public generalGammaDistributionFixedCategories {

public:
	explicit gammaDistributionFixedCategories(const Vdouble& fixedBoundaries, MDOUBLE alpha);
	explicit gammaDistributionFixedCategories(const gammaDistributionFixedCategories& other);
	explicit gammaDistributionFixedCategories(MDOUBLE alpha, int catNum);
	virtual ~gammaDistributionFixedCategories() {}
	virtual distribution* clone() const { return new gammaDistributionFixedCategories(*this); }
	virtual void setGammaParameters(int in_number_of_categories, MDOUBLE alpha);
	virtual void setAlpha(MDOUBLE newAlpha);
	virtual void change_number_of_categories(int in_number_of_categories);
	// to prevent the user from using alpha!=beta
	virtual void setGammaParameters(int numOfCategories ,MDOUBLE alpha, MDOUBLE beta) {
		if (alpha!=beta)
			errorMsg::reportError("gammaDistributionFixedCategories::setGammaParameters : can not set beta because alpha must be equal to beta");
		generalGammaDistributionFixedCategories::setGammaParameters(numOfCategories,alpha,beta);
	}
	virtual void setBeta(MDOUBLE newBeta) {
		errorMsg::reportError("generalGammaDistributionFixedCategories::setBeta : can not set beta because alpha=beta");
	}
};



#endif

