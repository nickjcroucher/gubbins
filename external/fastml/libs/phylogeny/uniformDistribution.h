// $Id: uniformDistribution.h 5807 2009-01-20 09:23:51Z adido $

 // version 2.00 
// last modified 21 Mar 2004
#ifndef ___FLAT_DIST
#define ___FLAT_DIST

/************************************************************ 
This represents a uniform distribution of one column (rectangular distribution) between 
a (lower_bound) and b (upper_bound)

		|---|
________|___|_____
		a	b
the distribution (or rather (a,b)) is divided into categories (portions of the distribution)
, where _rates is a vector with the median value for each category. _ratesProb represents
the probability of each category.
_globalRate represents the rate for two joint genes.
************************************************************/


#include "definitions.h"
#include "distribution.h"

class uniformDistribution : public distribution {

public:
	explicit uniformDistribution(const int numOfCategories, MDOUBLE lowerBound, 
		MDOUBLE upperBound);
	explicit uniformDistribution(){_globalRate=1.0;};
	explicit uniformDistribution(const uniformDistribution& other);
        
	virtual ~uniformDistribution() {};

	const int categories() const {return _rates.size();}
	virtual void change_number_of_categories(int in_number_of_categories);
	virtual const MDOUBLE rates(const int i) const {return _rates[i]*_globalRate;}
	virtual const MDOUBLE ratesProb(const int i) const {return _ratesProb[i];}
	virtual distribution* clone() const { return new uniformDistribution(*this); }
	virtual void setGlobalRate(const MDOUBLE x) {_globalRate = x;}
 	virtual MDOUBLE getGlobalRate() const {return _globalRate;}

	virtual const MDOUBLE getCumulativeProb(const MDOUBLE x) const;
	MDOUBLE getBorder(const int i) const ; 	//return the ith border. Note:  _bonderi[0] = m_lowerLimit, _bondery[categories()] = m_upperLimit

	void setUniformParameters(const int numOfCategories, MDOUBLE lowerBound, MDOUBLE upperBound);



private:	
	Vdouble _rates;
	Vdouble _ratesProb;
	MDOUBLE _globalRate;

	MDOUBLE _interval;
	MDOUBLE _upperBound;
	MDOUBLE _lowerBound;
};


#endif

//TO DO:
//1. change categories() to numOfCategories()


