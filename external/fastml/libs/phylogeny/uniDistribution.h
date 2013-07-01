// $Id: uniDistribution.h 2812 2007-11-25 10:32:11Z itaymay $

 // version 2.00 
// last modified 21 Mar 2004
#ifndef ___UNIFORM_DIST
#define ___UNIFORM_DIST

#include "distribution.h"

/***********************************************************
 This represents a distribution of one line over the value 1:
		|
________|________
		1
_globalRate represents the rate for two joint genes.
************************************************************/

class uniDistribution : public distribution {

public:
	uniDistribution() {_globalRate=1;}
	virtual const int categories() const { return 1;}
	virtual void change_number_of_categories(int in_number_of_categories);
	virtual const MDOUBLE rates(const int i) const { return _globalRate;};
	virtual const MDOUBLE ratesProb(const int i) const { return 1.0;};
	virtual distribution* clone() const { return new uniDistribution(*this); }
 	virtual void setGlobalRate(const MDOUBLE x) {_globalRate = x;}
 	virtual MDOUBLE getGlobalRate() const{return _globalRate;}
	virtual const MDOUBLE getCumulativeProb(const MDOUBLE x) const 	{
		if (x<1.0) return 0.0;	else return 1.0;
	} 

	MDOUBLE _globalRate;
};

#endif

