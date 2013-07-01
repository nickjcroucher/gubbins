// $Id: betaDistribution.cpp 3985 2008-05-11 11:00:44Z adido $

#include "betaDistribution.h"
#include "gammaUtilities.h"
#include "betaUtilities.h"
#include "errorMsg.h"
#include "logFile.h"
#include <cmath>


betaDistribution::betaDistribution() 
{
	_alpha = 0.0;
	_beta = 0.0;
	_boundary.resize(0,0);
	_rates.resize(0,0);
	_ratesProb.resize(0,0);
	_globalRate = 1;//??? 0.5 or 1
	_discretizationType = MEDIAN;
}

// note that the order of initalization makes a diffrence.
betaDistribution::betaDistribution(const betaDistribution& other) : 
	_boundary(other._boundary),
	_alpha(other._alpha), 
	_beta(other._beta), 
	_rates(other._rates),
	_ratesProb(other._ratesProb),
	_globalRate(other._globalRate),
	_discretizationType(other._discretizationType){
}

betaDistribution::betaDistribution(MDOUBLE alpha,MDOUBLE beta,int in_number_of_categories,discretizationType in_discretizationType) :distribution(){
	_globalRate=1.0;
	_discretizationType = in_discretizationType;
	setBetaParameters(in_number_of_categories,alpha,beta);
}

betaDistribution::~betaDistribution() {
	_boundary.clear();
	_rates.clear();
	_ratesProb.clear();
}

void betaDistribution::setAlpha(MDOUBLE in_alpha) {
	if (in_alpha == _alpha) 
		return;
	setBetaParameters(categories(), in_alpha, _beta);
}

void betaDistribution::setBeta(MDOUBLE in_beta) {
	if (in_beta == _beta)
		return;
	setBetaParameters( categories(), _alpha, in_beta);
}

void betaDistribution::setDiscretizationType(discretizationType in_discretizationType) {
	if (in_discretizationType == _discretizationType)
		return;
	_discretizationType = in_discretizationType;
	if (categories() > 1) 	
		fill_rates();

}
void betaDistribution::change_number_of_categories(int in_number_of_categories) {
	if (in_number_of_categories == categories())
		return;
	setBetaParameters( in_number_of_categories, _alpha, _beta);
}

void betaDistribution::setBetaParameters(int in_number_of_categories, MDOUBLE in_alpha, MDOUBLE in_beta) {
	if ((in_alpha == _alpha) && (in_beta == _beta) && (in_number_of_categories == categories()))
		return;
	
	
	if (in_alpha < MINIMUM_ALPHA_PARAM)	
		in_alpha = MINIMUM_ALPHA_PARAM;// when alpha is very small there are underflaw problems
	if (in_beta < MINIMUM_ALPHA_PARAM)	
		in_beta = MINIMUM_ALPHA_PARAM;// when beta is very small there are underflaw problems

	_alpha = in_alpha;
	_beta = in_beta;
	_rates.clear();
	_rates.resize(in_number_of_categories);
	_ratesProb.clear();
	_ratesProb.resize(in_number_of_categories, 1.0/in_number_of_categories);
	_boundary.clear();
	_boundary.resize(in_number_of_categories+1);
	if (in_number_of_categories==1) {
		_rates[0] = 1.0;
		return;
	}
	if (categories() > 1) {	
		fill_rates();
		return ;
	}
	
}
int betaDistribution::fill_rates() {
	fill_boundaries();
	int i;
	//LOG(5,<<endl<<" alpha = "<<_alpha<<" beta = "<< _beta<<endl);
	//for (i=0; i<=categories(); ++i) cout<<endl<<_boundary[i];
		//LOG(5,<<"\n====== the r categories are =====\n");
	for (i=0; i<categories(); ++i) {
		if (_discretizationType == MEAN)
			_rates[i]=computeAverage_r(_boundary[i], _boundary[i+1], _alpha, _beta, categories());
		else //_discretizationType == MEDIAN
			_rates[i] =inverseCDFBeta(_alpha, _beta,static_cast<MDOUBLE>(i*2 +1)/(2*categories()));
		//LOG(5,<<_rates[i]<<endl);
	}
	//LOG(5,<<endl<<_alpha<<endl);
	return 0;
}

int betaDistribution::fill_boundaries() {
	int i;
	//LOG(5,<<endl<<"========BOUNDARY============="<<endl);
	for (i=1; i<categories(); ++i)
	{
		_boundary[i]=inverseCDFBeta(_alpha, _beta,static_cast<MDOUBLE>(i)/categories());
		//LOG(5,<<"_boundary[ "<<i<<"] ="<<_boundary[i]<<endl);
	}
	_boundary[0]=0;
	_boundary[i]=1;
	
	return 0;
}


const MDOUBLE betaDistribution::getCumulativeProb(const MDOUBLE x) const
{//	 
	//since r~gamma(alpha, beta) then beta*r~ gamma(alpha,1)=gammp
	//here we assume alpha=beta
	return incompleteBeta(_alpha,_beta,x);
}



