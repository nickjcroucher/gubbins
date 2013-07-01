#ifndef ___GENERAL_GAMMA_DIST_LAGUERRE_FIXED_CATEGORIES
#define ___GENERAL_GAMMA_DIST_LAGUERRE_FIXED_CATEGORIES
/************************************************************
This class differ from the regular generalGammaDistribution in that 
the rateCategories are fixed according to the user's decision. 
Thus, only the probability of each category change for each specific alpha and beta values but 
the rate categories themselves are constant.
************************************************************/
#include "definitions.h"
#include "generalGammaDistribution.h"
#include "errorMsg.h"
class generalGammaDistributionFixedCategories : public generalGammaDistribution {

public:
	explicit generalGammaDistributionFixedCategories(const Vdouble& fixedBoundaries, MDOUBLE alpha, MDOUBLE beta);
	explicit generalGammaDistributionFixedCategories(const Vdouble& fixedRates, const Vdouble& boundaries, MDOUBLE alpha, MDOUBLE beta);
	explicit generalGammaDistributionFixedCategories(MDOUBLE alpha, MDOUBLE beta, int catNum);
	explicit generalGammaDistributionFixedCategories(const generalGammaDistributionFixedCategories& other);
	virtual ~generalGammaDistributionFixedCategories() {}
	virtual distribution* clone() const { return new generalGammaDistributionFixedCategories(*this); }
	virtual void change_number_of_categories(int in_number_of_categories); 
	virtual void setGammaParameters(int numOfCategories ,MDOUBLE alpha, MDOUBLE beta);
	virtual void setFixedCategories(const Vdouble& fixedBoundaries);

protected:
	virtual void setDefaultBoundaries(int catNum);
	virtual void setFixedCategories();
	virtual void fill_mean();
	virtual void computeRatesProbs();

};



#endif

