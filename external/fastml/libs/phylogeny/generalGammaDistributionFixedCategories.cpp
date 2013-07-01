#include "generalGammaDistributionFixedCategories.h"
#include "errorMsg.h"
#include "gammaUtilities.h"


generalGammaDistributionFixedCategories::generalGammaDistributionFixedCategories(const Vdouble& fixedBoundaries, MDOUBLE alpha, MDOUBLE beta) :
generalGammaDistribution()
{
	_alpha = alpha;
	_beta = beta;
	setFixedCategories(fixedBoundaries);
}

generalGammaDistributionFixedCategories::generalGammaDistributionFixedCategories(const Vdouble& fixedRates, const Vdouble& boundaries, MDOUBLE alpha, MDOUBLE beta) :
generalGammaDistribution()
{
	if ((fixedRates.size() + 1) !=  boundaries.size())
		errorMsg::reportError("error in generalGammaDistributionFixedCategories constructor");
	_alpha = alpha;
	_beta = beta;
	_rates = fixedRates;
	_bonderi = boundaries;
	computeRatesProbs();
}



generalGammaDistributionFixedCategories::generalGammaDistributionFixedCategories(MDOUBLE alpha, MDOUBLE beta, int catNum)
: generalGammaDistribution()
{
	_alpha = alpha;
	_beta = beta;
	setDefaultBoundaries(catNum);
}



generalGammaDistributionFixedCategories::generalGammaDistributionFixedCategories(const generalGammaDistributionFixedCategories& other) 
: generalGammaDistribution(other)
{}
void generalGammaDistributionFixedCategories::change_number_of_categories(int in_number_of_categories) 
{
	setDefaultBoundaries(in_number_of_categories);
}


void generalGammaDistributionFixedCategories::setFixedCategories(const Vdouble& fixedBoundaries){

	if (fixedBoundaries.size()<2)
		errorMsg::reportError("Error in generalGammaDistributionFixedCategories::setFixedCategories : at least two boundaries are required");
	if (fixedBoundaries[0] > 0.0)
		errorMsg::reportError("Error in generalGammaDistributionFixedCategories::setFixedCategories : first boundary should be zero");
	
	_bonderi = fixedBoundaries;
	if (_bonderi[_bonderi.size()] > VERYBIG/10000.0)
		 _bonderi[_bonderi.size()] = VERYBIG/10000.0; // to avoid overflow 

	setFixedCategories();
}

void generalGammaDistributionFixedCategories::setFixedCategories() {
	fill_mean();
	computeRatesProbs();
}

void generalGammaDistributionFixedCategories::fill_mean()
{
	int numOfCategories = _bonderi.size()-1;
	if (numOfCategories == 0)
		errorMsg::reportError("Error in gammaDistributionFixedCategories::fill_mean, fixed boundaries must be first initialized");
	_rates.clear();
	_rates.resize(numOfCategories,0.0);
	int cat;
	for (cat=0; cat<numOfCategories-1; ++cat) {
		_rates[cat] = (_bonderi[cat]+_bonderi[cat+1])/2.0;
	}
	if (numOfCategories>1) {
		//the rate of the last category cannot be the middle of its boundaries, since the upper bound is infinite
		MDOUBLE increment = _bonderi[cat] - _rates[cat-1];
		_rates[cat] = _bonderi[cat] + 2*increment;
	} else {
		_rates[0] = 1;
	}
}


// this function is here to override the inherited function
// note that the rates themselves and the boundaries do not change.
// the number of categories cannot be changed, since fixed categories must be given before
void generalGammaDistributionFixedCategories::setGammaParameters (int in_number_of_categories, MDOUBLE in_alpha, MDOUBLE in_beta) {
	if (in_number_of_categories==1) {
		_rates[0] = 1.0;
		return;
	}	
	if (in_number_of_categories != categories())
		errorMsg::reportError("generalGammaDistributionFixedCategories::setGammaParameters: the number of categories cannot be changed, first call setFixedCategories");
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

void generalGammaDistributionFixedCategories::computeRatesProbs(){
	MDOUBLE totalProb = 0.0;
	MDOUBLE catProb = 0.0;
	MDOUBLE lowerBoundaryProb = 0.0;
	MDOUBLE upperBoundaryProb = 0.0;
	int cat;
	_ratesProb.clear();
	_ratesProb.resize(categories());
	for (cat = 0; cat < categories()-1; ++cat) {
		upperBoundaryProb = getCumulativeProb(_bonderi[cat+1]);
		catProb = upperBoundaryProb - lowerBoundaryProb;
		_ratesProb[cat] = catProb;
		totalProb += catProb;
		lowerBoundaryProb = upperBoundaryProb;
	}
	_ratesProb[cat] = 1.0 - totalProb;
}

void generalGammaDistributionFixedCategories::setDefaultBoundaries(int catNum) 
{
	_bonderi.clear();
	_bonderi.resize(catNum+1,0.0);
	_bonderi[0] = 0;
	_bonderi[catNum] = VERYBIG/10000.0; //to avoid overflow
	switch (catNum)
	{
	case 1:
		break;
	case 2:
		_bonderi[1] = 1.0;
		break;
	case 3:
		_bonderi[1] = 0.5;
		_bonderi[2] = 1.0;
		break;
	case 4:
		_bonderi[1] = 0.5;
		_bonderi[2] = 1.0;
		_bonderi[3] = 1.5;
		break;
	case 5:
		_bonderi[1] = 0.4;
		_bonderi[2] = 0.8;
		_bonderi[3] = 1.2;
		_bonderi[4] = 1.6;
		break;
	case 10:
		_bonderi[1] = 0.01;
		_bonderi[2] = 0.1;
		_bonderi[3] = 0.25;
		_bonderi[4] = 0.55;
		_bonderi[5] = 0.95;
		_bonderi[6] = 1.5;
		_bonderi[7] = 3.0;
		_bonderi[8] = 5.0;
		_bonderi[9] = 7.0;
		break;
	case 16:
		_bonderi[1] = 0.001;
		_bonderi[2] = 0.01;
		_bonderi[3] = 0.1;
		_bonderi[4] = 0.15;
		_bonderi[5] = 0.35;
		_bonderi[6] = 0.55;
		_bonderi[7] = 0.75;
		_bonderi[8] = 0.95;
		_bonderi[9] = 1.5;
		_bonderi[10] = 3.0;
		_bonderi[11] = 4.5;
		_bonderi[12] = 6.0;
		_bonderi[13] = 7.5;
		_bonderi[14] = 9.0;
		_bonderi[15] = 12.0;
		break;
	default:
		errorMsg::reportError("error in generalGammaDistributionFixedCategories::setDefaultBoundaries");
	}

	setFixedCategories();
}

//void generalGammaDistributionFixedCategories::getDefaultRates(int catNum, Vdouble& fixedRates)
//{
//	fixedRates.resize(catNum, 0.0);
//	switch (catNum)
//	{
//	case 1:
//		fixedRates[0] = 1.0;
//		break;
//	case 2:
//		fixedRates[0] = 0.5;
//		fixedRates[1] = 1.5;
//		break;
//	case 3:
//		fixedRates[0] = 0.05;
//		fixedRates[1] = 0.5;
//		fixedRates[2] = 1.5;
//		break;
//	case 5:
//		fixedRates[0] = 0.05;
//		fixedRates[1] = 0.3;
//		fixedRates[2] = 0.6;
//		fixedRates[3] = 1.5;
//		fixedRates[4] = 5.0;
//		break;
//	case 8:
//		fixedRates[0] = 0.05;
//		fixedRates[1] = 0.15;
//		fixedRates[2] = 0.35;
//		fixedRates[3] = 0.6;
//		fixedRates[4] = 0.85;
//		fixedRates[5] = 1.5;
//		fixedRates[6] = 3.0;
//		fixedRates[7] = 5.0;
//		break;
//	case 12:
//		fixedRates[0] = 0.05;
//		fixedRates[1] = 0.15;
//		fixedRates[2] = 0.35;
//		fixedRates[3] = 0.55;
//		fixedRates[4] = 0.75;
//		fixedRates[5] = 0.95;
//		fixedRates[6] = 1.5;
//		fixedRates[7] = 3.0;
//		fixedRates[8] = 4.5;
//		fixedRates[9] = 6.0;
//		fixedRates[10] = 7.5;
//		fixedRates[11] = 9.0;
//		break;
//	case 16:
//		fixedRates[0] = 0.00000001;
//		fixedRates[1] = 0.001;
//		fixedRates[2] = 0.01;
//		fixedRates[3] = 0.1;
//		fixedRates[4] = 0.15;
//		fixedRates[5] = 0.35;
//		fixedRates[6] = 0.55;
//		fixedRates[7] = 0.75;
//		fixedRates[8] = 0.95;
//		fixedRates[9] = 1.5;
//		fixedRates[10] = 3.0;
//		fixedRates[11] = 4.5;
//		fixedRates[12] = 6.0;
//		fixedRates[13] = 7.5;
//		fixedRates[14] = 9.0;
//		fixedRates[15] = 12.0;
//		break;
//	case 24:
//		fixedRates[0] = 0.000000000000001;
//		fixedRates[1] = 1;
//		fixedRates[2] = 2;
//		fixedRates[3] = 3;
//		fixedRates[4] = 4;
//		fixedRates[5] = 5;
//		fixedRates[6] = 6;
//		fixedRates[7] = 7;
//		fixedRates[8] = 8;
//		fixedRates[9] = 9;
//		fixedRates[10] = 10;
//		fixedRates[11] = 11;
//		fixedRates[12] = 12;
//		fixedRates[13] = 13;
//		fixedRates[14] = 14;
//		fixedRates[15] = 15;
//		fixedRates[16] = 16;
//		fixedRates[17] = 17;
//		fixedRates[18] = 18;
//		fixedRates[19] = 19;
//		fixedRates[20] = 20;
//		fixedRates[21] = 21;
//		fixedRates[22] = 22;
//		fixedRates[23] = 23;
//		break;
//	case 32:
//		fixedRates[0] = 0.00000001;
//		fixedRates[1] = 0.0000001;
//		fixedRates[2] = 0.000001;
//		fixedRates[3] = 0.00001;
//		fixedRates[4] = 0.0001;
//		fixedRates[5] = 0.001;
//		fixedRates[6] = 0.01;
//		fixedRates[7] = 0.1;
//		fixedRates[8] = 0.15;
//		fixedRates[9] = 0.2;
//		fixedRates[10] = 0.25;
//		fixedRates[11] = 0.3;
//		fixedRates[12] = 0.35;
//		fixedRates[13] = 0.4;
//		fixedRates[14] = 0.45;
//		fixedRates[15] = 0.5;
//		fixedRates[16] = 0.6;
//		fixedRates[17] = 0.7;
//		fixedRates[18] = 0.8;
//		fixedRates[19] = 0.9;
//		fixedRates[20] = 1.0;
//		fixedRates[21] = 1.2;
//		fixedRates[22] = 1.4;
//		fixedRates[23] = 1.6;
//		fixedRates[24] = 1.8;
//		fixedRates[25] = 2.0;
//		fixedRates[26] = 2.5;
//		fixedRates[27] = 3.0;
//		fixedRates[28] = 4.0;
//		fixedRates[29] = 5.0;		
//		fixedRates[30] = 7.5;
//		fixedRates[31] = 15.0;
//		break;
//	case 36:
//		fixedRates[0] = 0.00000001;
//		fixedRates[1] = 0.0000001;
//		fixedRates[2] = 0.000001;
//		fixedRates[3] = 0.00001;
//		fixedRates[4] = 0.0001;
//		fixedRates[5] = 0.001;
//		fixedRates[6] = 0.01;
//		fixedRates[7] = 0.1;
//		fixedRates[8] = 0.15;
//		fixedRates[9] = 0.2;
//		fixedRates[10] = 0.25;
//		fixedRates[11] = 0.3;
//		fixedRates[12] = 0.35;
//		fixedRates[13] = 0.4;
//		fixedRates[14] = 0.45;
//		fixedRates[15] = 0.5;
//		fixedRates[16] = 0.6;
//		fixedRates[17] = 0.7;
//		fixedRates[18] = 0.8;
//		fixedRates[19] = 0.9;
//		fixedRates[20] = 1.0;
//		fixedRates[21] = 1.2;
//		fixedRates[22] = 1.4;
//		fixedRates[23] = 1.6;
//		fixedRates[24] = 1.8;
//		fixedRates[25] = 2.0;
//		fixedRates[26] = 2.5;
//		fixedRates[27] = 3.0;
//		fixedRates[28] = 4.0;
//		fixedRates[29] = 5.0;		
//		fixedRates[30] = 7.5;
//		fixedRates[31] = 10.0;
//		fixedRates[32] = 12.5;
//		fixedRates[33] = 15.0;
//		fixedRates[34] = 20.0;
//		fixedRates[35] = 30.0;
//		break;
//
//	default:
//		errorMsg::reportError("error in generalGammaDistributionFixedCategories::getFixedCategories");
//	}
//
//}
