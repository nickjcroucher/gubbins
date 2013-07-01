// $Id: generalGammaDistributionLaguerre.cpp 2865 2007-11-27 11:00:26Z itaymay $
#include "generalGammaDistributionLaguerre.h"
#include "gammaUtilities.h"
#include "errorMsg.h"
#include "GLaguer.h"
#include <cmath>

generalGammaDistributionLaguerre::generalGammaDistributionLaguerre() 
: generalGammaDistribution()
{
}

generalGammaDistributionLaguerre::generalGammaDistributionLaguerre(const generalGammaDistributionLaguerre& other) :
	generalGammaDistribution(other)
{
}

generalGammaDistributionLaguerre::generalGammaDistributionLaguerre(MDOUBLE alpha,MDOUBLE beta,int in_number_of_categories)
: generalGammaDistribution()
{
	//The Laguerre function returns NULL values for very large numebr of categories (for example 700 categories with alpha = 1.5 and beta = 1.3)
//	if (in_number_of_categories > 200)
//		errorMsg::reportError("generalGammaDistributionLaguerre cannot work with more than 200 categories");
	_globalRate=1.0;
	setGammaParameters(in_number_of_categories,alpha,beta);
}

generalGammaDistributionLaguerre::~generalGammaDistributionLaguerre()
{
}


void generalGammaDistributionLaguerre::setGammaParameters(int in_number_of_categories, MDOUBLE in_alpha, MDOUBLE in_beta) {
	if ((in_alpha == _alpha) && (in_beta == _beta) && (in_number_of_categories == categories()))
		return;
	
	
	if (in_alpha < MINIMUM_ALPHA_PARAM)	
		in_alpha = MINIMUM_ALPHA_PARAM;// when alpha is very small there are underflaw problems
	if (in_beta < MINIMUM_ALPHA_PARAM)	
		in_beta = MINIMUM_ALPHA_PARAM;// when beta is very small there are underflaw problems

	_alpha = in_alpha;
	_beta = in_beta;
	_rates.clear();
	//_rates.resize(in_number_of_categories);
	_rates.resize(0);
	_ratesProb.clear();
	//_ratesProb.resize(in_number_of_categories);
	_ratesProb.resize(0);
	if (in_number_of_categories==1) {
		_rates.push_back(1.0);
		_ratesProb.push_back(1.0);
		return;
	}
	if (in_number_of_categories > 1) {	
		fillRatesAndProbs(in_number_of_categories);
		return ;
	}
	
}


MDOUBLE generalGammaDistributionLaguerre::getBorder(const int i) const
{
	errorMsg::reportError("With the Laguerre method the categories do not have a well defined border");
	return -1;
}	


void generalGammaDistributionLaguerre::fillRatesAndProbs(int catNum)
{
	Vdouble weights, abscissas;
	GLaguer lg(catNum, _alpha - 1, abscissas, weights);
	MDOUBLE sumP = 0.0;

	MDOUBLE gamAlpha = exp(gammln(_alpha));
	for (int i = 0; i < catNum; ++i)
	{
		//if (sumP > 0.99)
		//{
		//	_ratesProb.push_back(1-sumP);
		//	_rates.push_back(abscissas[i] / _beta); 
		//	break;
		//}

		_ratesProb.push_back(weights[i] / gamAlpha);
		_rates.push_back(abscissas[i] / _beta); 
		sumP += _ratesProb[i];
		//cerr<<i<<" rate = "<<_rates[i]<<" Pr = "<<_ratesProb[i]<<" sum = "<<sumP<<endl;
	}
	for (int j = 0; j < _ratesProb.size(); ++j)
	{
		_ratesProb[j] /= sumP;
	}
}


/*
void generalGammaDistributionLaguerre::fillRatesAndProbs(int catNum)
{
	Vdouble weights, abscissas;
	GLaguer lg(categories(), _alpha - 1, abscissas, weights);

	MDOUBLE gamAlpha = exp(gammln(_alpha));
	for (int i = 0; i < categories(); ++i)
	{
		_ratesProb[i] = weights[i] / gamAlpha;
		_rates[i] = abscissas[i] / _beta; 
	}
}
*/

