// $Id: gammaDistribution.h 2768 2007-11-22 12:57:44Z osnatz $

#ifndef ___GAMMA_DIST_LAGUERRE
#define ___GAMMA_DIST_LAGUERRE
/************************************************************
This distribution can take several forms depending on its free parameter alpha 
(beta is assumed to be equal to alpha). For an extensive exlpanation of this distribution
see http://mathworld.wolfram.com/GammaDistribution.html.
please note that the borders of the categories are defined according to calculation of 
the gamma integral, according to numerical recipes in gammaUtilities
_globalRate represents the rate for two joint genes.
************************************************************/
#include "definitions.h"
#include "generalGammaDistributionLaguerre.h"
#include "errorMsg.h"

class gammaDistributionLaguerre : public generalGammaDistributionLaguerre {

public:
	explicit gammaDistributionLaguerre() {}
	explicit gammaDistributionLaguerre(MDOUBLE alpha,int in_number_of_categories);
	explicit gammaDistributionLaguerre(const gammaDistributionLaguerre& other);
	virtual ~gammaDistributionLaguerre() {}
	virtual distribution* clone() const { return new gammaDistributionLaguerre(*this); }

	virtual void setAlpha(MDOUBLE newAlpha);
	virtual void setGammaParameters(int numOfCategories=1 ,MDOUBLE alpha=1);
	virtual void change_number_of_categories(int in_number_of_categories);
	// to prevent the user from using alpha!=beta
	virtual void setGammaParameters(int numOfCategories ,MDOUBLE alpha, MDOUBLE beta); 
	virtual void setBeta(MDOUBLE newBeta) {errorMsg::reportError("gammaDistributionLaguerre::setBeta : can not set beta because alpha=beta");
	}
};
#endif
