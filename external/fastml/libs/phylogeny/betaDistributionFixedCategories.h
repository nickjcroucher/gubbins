#ifndef ___BETA_FIXED_CATEGORIES_CATEGORIES
#define ___BETA_FIXED_CATEGORIES_CATEGORIES
/************************************************************
This class differ from the regular betaDistribution in that 
the rateCategories are fixed according to the user's decision. 
Thus, only the probability of each category change for each specific alpha and beta values but 
the rate categories themselves are constant.
************************************************************/
#include "definitions.h"
#include "betaDistribution.h"
#include "errorMsg.h"
class betaDistributionFixedCategories : public betaDistribution {

public:
	explicit betaDistributionFixedCategories(const Vdouble& fixedBoundaries, MDOUBLE alpha, MDOUBLE beta);
	explicit betaDistributionFixedCategories(const Vdouble& fixedRates, const Vdouble& boundaries, MDOUBLE alpha, MDOUBLE beta);
	explicit betaDistributionFixedCategories(MDOUBLE alpha, MDOUBLE beta, int catNum);
	explicit betaDistributionFixedCategories(const betaDistributionFixedCategories& other);
	explicit betaDistributionFixedCategories();
	virtual ~betaDistributionFixedCategories() {}
	virtual distribution* clone() const { return new betaDistributionFixedCategories(*this); }
	virtual void change_number_of_categories(int in_number_of_categories); 
	virtual void setBetaParameters(int numOfCategories ,MDOUBLE alpha, MDOUBLE beta);
	virtual void setFixedCategories(const Vdouble& fixedBoundaries);

protected:
	virtual void setDefaultBoundaries(int catNum);
	virtual void setFixedCategories();
	virtual void fill_mean();
	virtual void computeRatesProbs();

};



#endif

