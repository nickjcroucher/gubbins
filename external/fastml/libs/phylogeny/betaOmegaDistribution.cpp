// $Id: betaOmegaDistribution.cpp 962 2006-11-07 15:13:34Z privmane $

#include "betaOmegaDistribution.h"
#include "gammaUtilities.h"
#include "betaUtilities.h"
#include "errorMsg.h"
#include "logFile.h"
#include <cmath>


betaOmegaDistribution::betaOmegaDistribution() 
{
	_omega=1;
	_betaProb = 0.5;
}

// note that the order of initalization makes a diffrence.
betaOmegaDistribution::betaOmegaDistribution(const betaOmegaDistribution& other) : 
	_betaDistr(other._betaDistr),
	_omega(other._omega),
	_betaProb(other._betaProb){
}

betaOmegaDistribution::betaOmegaDistribution(MDOUBLE alpha,MDOUBLE beta,int in_number_of_categories,MDOUBLE betaProb,MDOUBLE omega) :distribution(){
	_omega = omega;
	_betaProb = betaProb;
	_betaDistr.setGlobalRate(1.0);
	_betaDistr.setBetaParameters(in_number_of_categories,alpha,beta);
}

betaOmegaDistribution::~betaOmegaDistribution() {}


void betaOmegaDistribution::setBetaOmegaParameters(int in_number_of_categories,MDOUBLE alpha, MDOUBLE beta,MDOUBLE betaProb,MDOUBLE omega){
	_omega = omega;
	_betaProb = betaProb;
	_betaDistr.setBetaParameters(in_number_of_categories, alpha,  beta);

}
const MDOUBLE betaOmegaDistribution::ratesProb(const int i) const {
	if (i < _betaDistr.categories())
		return _betaDistr.ratesProb(i)*_betaProb;
	else return (1-_betaProb); //omega prob
}


const MDOUBLE betaOmegaDistribution::rates(const int i) const {
	if (i < _betaDistr.categories())
		return _betaDistr.rates(i);
	else return _omega; //omega
}



const MDOUBLE betaOmegaDistribution::getCumulativeProb(const MDOUBLE x) const
{ return _betaDistr.getCumulativeProb(x);
}




