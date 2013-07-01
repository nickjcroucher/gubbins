#include "betaDistributionFixedCategories.h"
#include "errorMsg.h"
#include "gammaUtilities.h"


betaDistributionFixedCategories::betaDistributionFixedCategories(const Vdouble& fixedBoundaries, MDOUBLE alpha, MDOUBLE beta) :
betaDistribution()
{
	_alpha = alpha;
	_beta = beta;
	setFixedCategories(fixedBoundaries);
}


betaDistributionFixedCategories::betaDistributionFixedCategories(const Vdouble& fixedRates, const Vdouble& boundaries, MDOUBLE alpha, MDOUBLE beta) :
betaDistribution()
{
	if ((fixedRates.size() + 1) !=  boundaries.size())
		errorMsg::reportError("error in betaDistributionFixedCategories constructor");
	_alpha = alpha;
	_beta = beta;
	_rates = fixedRates;
	_boundary = boundaries;
	computeRatesProbs();
}



betaDistributionFixedCategories::betaDistributionFixedCategories(MDOUBLE alpha, MDOUBLE beta, int catNum)
: betaDistribution()
{
	_alpha = alpha;
	_beta = beta;
	setDefaultBoundaries(catNum);
}

betaDistributionFixedCategories::betaDistributionFixedCategories()
: betaDistribution()
{
	_alpha = 0.5;
	_beta = 0.5;
	setDefaultBoundaries(10);
}

betaDistributionFixedCategories::betaDistributionFixedCategories(const betaDistributionFixedCategories& other) 
: betaDistribution(other)
{}
void betaDistributionFixedCategories::change_number_of_categories(int in_number_of_categories) 
{
	setDefaultBoundaries(in_number_of_categories);
}


void betaDistributionFixedCategories::setFixedCategories(const Vdouble& fixedBoundaries){

	if (fixedBoundaries.size()<2)
		errorMsg::reportError("Error in generalGammaDistributionFixedCategories::setFixedCategories : at least two boundaries are required");
	if (fixedBoundaries[0] > 0.0)
		errorMsg::reportError("Error in generalGammaDistributionFixedCategories::setFixedCategories : first boundary should be zero");
	
	_boundary = fixedBoundaries;
	if (_boundary[_boundary.size()] > VERYBIG/10000.0)
		 _boundary[_boundary.size()] = VERYBIG/10000.0; // to avoid overflow 

	setFixedCategories();
}

void betaDistributionFixedCategories::setFixedCategories() {
	fill_mean();
	computeRatesProbs();
}

void betaDistributionFixedCategories::fill_mean()
{
	int numOfCategories = _boundary.size()-1;
	if (numOfCategories == 0)
		errorMsg::reportError("Error in gammaDistributionFixedCategories::fill_mean, fixed boundaries must be first initialized");
	_rates.clear();
	_rates.resize(numOfCategories,0.0);
	int cat;
	for (cat=0; cat<numOfCategories; ++cat) {
		_rates[cat] = (_boundary[cat]+_boundary[cat+1])/2.0; 
	}
	
}


// this function is here to override the inherited function
// note that the rates themselves and the boundaries do not change.
// the number of categories cannot be changed, since fixed categories must be given before
void betaDistributionFixedCategories::setBetaParameters (int in_number_of_categories, MDOUBLE in_alpha, MDOUBLE in_beta) {
	if (in_number_of_categories==1) {
		_rates[0] = 1.0;
		return;
	}
	if (in_number_of_categories != categories())
		errorMsg::reportError("betaDistributionFixedCategories::setGammaParameters: the number of categories cannot be changed, first call setFixedCategories");
	if ((in_alpha == _alpha) && (in_beta == _beta))
		return;

	if (in_alpha < MINIMUM_ALPHA_PARAM)	
		in_alpha = MINIMUM_ALPHA_PARAM;// when alpha is very small there are underflow problems
	if (in_beta < MINIMUM_ALPHA_PARAM)	
		in_beta = MINIMUM_ALPHA_PARAM;// when beta is very small there are underflaw problems

	_alpha = in_alpha;
	_beta = in_beta;
	computeRatesProbs();
}

void betaDistributionFixedCategories::computeRatesProbs(){
	MDOUBLE totalProb = 0.0;
	MDOUBLE catProb = 0.0;
	MDOUBLE lowerBoundaryProb = 0.0;
	MDOUBLE upperBoundaryProb = 0.0;
	int cat;
	_ratesProb.clear();
	_ratesProb.resize(categories());
	for (cat = 0; cat < categories()-1; ++cat) {
		upperBoundaryProb = getCumulativeProb(_boundary[cat+1]);
		catProb = upperBoundaryProb - lowerBoundaryProb;
		_ratesProb[cat] = catProb;
		totalProb += catProb;
		lowerBoundaryProb = upperBoundaryProb;
	}
	_ratesProb[cat] = 1.0 - totalProb;
}

void betaDistributionFixedCategories::setDefaultBoundaries(int catNum) 
{
	_boundary.clear();
	_boundary.resize(catNum+1,0.0);
	_boundary[0] = 0;
	_boundary[catNum] = 1.0; 
	switch (catNum)
	{
	case 1:
		break;
	case 2:
		_boundary[1] = 0.5;
		break;
	case 10:
		_boundary[1] = 0.1;
		_boundary[2] = 0.2;
		_boundary[3] = 0.3;
		_boundary[4] = 0.4;
		_boundary[5] = 0.5;
		_boundary[6] = 0.6;
		_boundary[7] = 0.7;
		_boundary[8] = 0.8;
		_boundary[9] = 0.9;
		break;
	default:
		errorMsg::reportError("error in betaDistributionFixedCategories::setDefaultBoundaries");
	}

	setFixedCategories();
}
