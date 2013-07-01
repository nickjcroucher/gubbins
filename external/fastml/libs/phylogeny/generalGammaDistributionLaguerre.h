// $Id: generalGammaDistributionLaguerre.h 2865 2007-11-27 11:00:26Z itaymay $
// version 1.00
// last modified Sep 2004

#ifndef ___GENERAL_GAMMA_DIST_LAGUERRE
#define ___GENERAL_GAMMA_DIST_LAGUERRE
/************************************************************
This class differ from the regular generalGammaDistribution in that 
the rateCategories and their probabilities are not constructed using Yang's quantile method.
Instead the general Guass-Laguerre quadrature method is used.
For example, if we want to compute the likelihood over the rate distribution, 
then we need to solve the integral

I[0_to_infinity]{P(data|r)*P(r)} 
	= I[0_to_infinity]{P(data|r)*b^a / Gamma(a)* exp(-b*r) * r^(a-1)dr}  //a = alpha, b = beta
	= b^(a)/Gamma(a) * I[0_to_infinity]{P(data|m/b) * exp(-m) * (m/b)^(a')/bdm}  ///substitute m=b*r, a'=a-1
	= 1/Gamma(a) * I[0_to_infinity]{P(data|m/b) * exp(-m) * m^a' dm}  //
Now - we can use the Guass-Laguerre formula, to get an approximation for the Integral.
The Xj and Wj are the absicassas and weights of the Laguerre polynoms
	= 1/Gamma(a) * sum[j = 0_to_catNum]{P(data|Xj/b) * Wj}  
  
The rates are the Xj/b and their priors is Wj/Gamma(a) 
The quadrature method is explained in Numerical Recipes (Press et al.; chapter 4.5) 
and is also mentioned in Felsenstein 2001 (JME 53: 447-455).
************************************************************/
#include "definitions.h"
#include "generalGammaDistribution.h"
class generalGammaDistributionLaguerre : public generalGammaDistribution {

public:
	explicit generalGammaDistributionLaguerre();
	explicit generalGammaDistributionLaguerre(MDOUBLE alpha, MDOUBLE beta, int in_number_of_categories);
	explicit generalGammaDistributionLaguerre(const generalGammaDistributionLaguerre& other);
	virtual ~generalGammaDistributionLaguerre();
	virtual void setGammaParameters(int numOfCategories ,MDOUBLE alpha, MDOUBLE beta);

	virtual distribution* clone() const { return new generalGammaDistributionLaguerre(*this); }
 	virtual MDOUBLE getBorder(const int i) const; 

protected:
	virtual void fillRatesAndProbs(int catNum);
};



#endif

