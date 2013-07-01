// $Id: betaDistribution.h 5803 2009-01-20 09:17:05Z adido $

#ifndef ___BETA_DIST
#define ___BETA_DIST
/************************************************************
This distribution can take several forms depending on its free parameters alpha,beta 
For an extensive exlpanation of this distribution
see http://mathworld.wolfram.com/BetaDistribution.html
************************************************************/
#include "definitions.h"
#include "distribution.h"

class betaDistribution : public distribution {

public:
	enum discretizationType{MEAN, MEDIAN};
	explicit betaDistribution(MDOUBLE alpha, MDOUBLE beta, int in_number_of_categories,discretizationType in_discretizationType = MEDIAN);
	explicit betaDistribution(const betaDistribution& other);
	explicit betaDistribution();
	virtual ~betaDistribution();
	virtual void setBetaParameters(int numOfCategories ,MDOUBLE alpha, MDOUBLE beta);

	virtual const int categories() const {return _rates.size();}
	virtual const MDOUBLE rates(const int i) const {return _rates[i]*_globalRate;}
	virtual const MDOUBLE ratesProb(const int i) const {return _ratesProb[i];}
	virtual distribution* clone() const { return new betaDistribution(*this); }
 	virtual void setGlobalRate(const MDOUBLE x) {_globalRate = x;}
 	virtual MDOUBLE getGlobalRate()const {return _globalRate;}
	virtual const MDOUBLE getCumulativeProb(const MDOUBLE x) const;
	virtual void setAlpha(MDOUBLE newAlpha);
	virtual MDOUBLE getAlpha() const {return _alpha;};
	virtual void setBeta(MDOUBLE newBeta);
	virtual MDOUBLE getBeta() const {return _beta;};
	virtual void setDiscretizationType(discretizationType in_discretizationType);
	virtual discretizationType getDiscretizationType() const {return _discretizationType;};

	virtual void change_number_of_categories(int in_number_of_categories);
	virtual MDOUBLE getBorder(const int i) const {return _boundary[i];}	//return the ith border. Note:  _bonderi[0] = 0, _bondery[categories()] = infinite


private:	
	int fill_rates();
	int fill_boundaries();
	
	
protected:
	MDOUBLE _alpha;
	MDOUBLE _beta;
	
	vector<MDOUBLE> _rates;
	vector<MDOUBLE> _ratesProb;
	MDOUBLE _globalRate;
	discretizationType _discretizationType;
	vector<MDOUBLE> _boundary;
	
};



#endif

