// $Id: generalGammaDistribution.cpp 2768 2007-11-22 12:57:44Z osnatz $

#include "generalGammaDistribution.h"
#include "gammaUtilities.h"
#include "errorMsg.h"
#include "logFile.h"
#include <cmath>


generalGammaDistribution::generalGammaDistribution() :
_alpha(0.0),
_beta(0.0),
_globalRate(1.0)
{
	_bonderi.resize(0,0);
	_rates.resize(0,0);
	_ratesProb.resize(0,0);
}

generalGammaDistribution::generalGammaDistribution(const generalGammaDistribution& other) : 
	
	_alpha(other._alpha),
	_beta(other._beta),
	_rates(other._rates),
	_ratesProb(other._ratesProb),
	_globalRate(other._globalRate),
	_bonderi(other._bonderi)
	{}


generalGammaDistribution::generalGammaDistribution(MDOUBLE alpha,MDOUBLE beta,int in_number_of_categories) :
	_globalRate(1.0) 
{
	setGammaParameters(in_number_of_categories,alpha,beta);
}

void generalGammaDistribution::setAlpha(MDOUBLE in_alpha) {
	if (in_alpha == _alpha) 
		return;
	setGammaParameters(categories(), in_alpha, _beta);
}

void generalGammaDistribution::setBeta(MDOUBLE in_beta) {
	if (in_beta == _beta)
		return;
	setGammaParameters( categories(), _alpha, in_beta);
}

void generalGammaDistribution::change_number_of_categories(int in_number_of_categories) {
	if (in_number_of_categories == categories())
		return;
	setGammaParameters( in_number_of_categories, _alpha, _beta);
}

void generalGammaDistribution::setGammaParameters(int in_number_of_categories, MDOUBLE in_alpha, MDOUBLE in_beta) {
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
	_bonderi.clear();
	_bonderi.resize(in_number_of_categories+1);
	if (in_number_of_categories==1) {
		_rates[0] = 1.0;
		return;
	}
	if (categories() > 1) {	
		fill_mean();
		return ;
	}
	
}
void generalGammaDistribution::fill_mean() {
	fill_bonderi();
	int i;
	//for (i=0; i<=categories(); ++i) cout<<endl<<bonderi[i];
	//LOG(5,<<"\n====== the r categories are =====\n");
	for (i=0; i<categories(); ++i) {
		_rates[i]=the_avarage_r_in_category_between_a_and_b(_bonderi[i], _bonderi[i+1], _alpha, _beta, categories());
		//LOG(5,<<meanG[i]<<endl);
	}
	//LOG(5,<<endl<<alpha<<endl);
	//return 0;
}

void generalGammaDistribution::fill_bonderi() {
	int i;
	for (i=1; i<categories(); ++i)
	{
		_bonderi[i]=search_for_z_in_dis_with_any_beta(_alpha, _beta,static_cast<MDOUBLE>(i)/categories());
	}
	_bonderi[0]=0;
	_bonderi[i]=VERYBIG/10000.0;// this is becuase we multiply bondei[i] by alpha or beta, and 
	// by this manipulation we avoid overflows...;
	
	//return 0;
}


const MDOUBLE generalGammaDistribution::getCumulativeProb(const MDOUBLE x) const
{//	 
	//since r~gamma(alpha, beta) then beta*r~ gamma(alpha,1)=gammp
	//here we assume alpha=beta
	return gammp(_alpha, x*_beta);
}
