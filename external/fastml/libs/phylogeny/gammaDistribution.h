// $Id: gammaDistribution.h 2862 2007-11-27 10:59:03Z itaymay $

#ifndef ___GAMMA_DIST
#define ___GAMMA_DIST
/************************************************************
This distribution can take several forms depending on its free parameter alpha 
(beta is assumed to be equal to alpha). For an extensive exlpanation of this distribution
see http://mathworld.wolfram.com/GammaDistribution.html.
please note that the borders of the categories are defined according to calculation of 
the gamma integral, according to numerical recipes in gammaUtilities
_globalRate represents the rate for two joint genes.
************************************************************/
#include "definitions.h"
#include "generalGammaDistribution.h"
#include "errorMsg.h"

class gammaDistribution : public generalGammaDistribution {

public:
	explicit gammaDistribution() {}
	explicit gammaDistribution(MDOUBLE alpha,int in_number_of_categories);
	explicit gammaDistribution(const gammaDistribution& other);
	virtual ~gammaDistribution() {}
	virtual distribution* clone() const { return new gammaDistribution(*this); }

	virtual void setAlpha(MDOUBLE newAlpha);
	virtual void setGammaParameters(int numOfCategories=1 ,MDOUBLE alpha=1);
	virtual void change_number_of_categories(int in_number_of_categories);
	// to prevent the user from using alpha!=beta
	virtual void setGammaParameters(int numOfCategories ,MDOUBLE alpha, MDOUBLE beta);
	virtual void setBeta(MDOUBLE newBeta) {errorMsg::reportError("gammaDistribution::setBeta : can not set beta because alpha=beta");}
};
#endif
