// $Id: generalGammaDistribution.h 3044 2007-12-18 15:54:50Z itaymay $

#ifndef ___GENERAL_GAMMA_DIST
#define ___GENERAL_GAMMA_DIST
/************************************************************
This distribution can take several forms depending on its free parameters alpha,beta 
(unalike gammaDist. alpha is not necessarily equal to beta). 
For an extensive exlpanation of this distribution
see http://mathworld.wolfram.com/GammaDistribution.html
************************************************************/
#include "definitions.h"
#include "distribution.h"

enum quadratureType {QUANTILE, LAGUERRE};

class generalGammaDistribution : public distribution {

public:
	explicit generalGammaDistribution();
	explicit generalGammaDistribution(MDOUBLE alpha, MDOUBLE beta, int in_number_of_categories);
	explicit generalGammaDistribution(const generalGammaDistribution& other);
	virtual ~generalGammaDistribution() {}
	virtual distribution* clone() const { return new generalGammaDistribution(*this); }
	
	virtual void setGammaParameters(int numOfCategories ,MDOUBLE alpha, MDOUBLE beta);
	virtual const int categories() const {return _rates.size();}
	virtual const MDOUBLE rates(const int i) const {return _rates[i]*_globalRate;}
	virtual const MDOUBLE ratesProb(const int i) const {return _ratesProb[i];}

 	virtual void setGlobalRate(const MDOUBLE x) {_globalRate = x;}
 	virtual MDOUBLE getGlobalRate()const {return _globalRate;}
	virtual const MDOUBLE getCumulativeProb(const MDOUBLE x) const;
	virtual void setAlpha(MDOUBLE newAlpha);
	virtual MDOUBLE getAlpha() const {return _alpha;}
	virtual void setBeta(MDOUBLE newBeta);
	virtual MDOUBLE getBeta() const {return _beta;}
	virtual void change_number_of_categories(int in_number_of_categories);
	virtual MDOUBLE getBorder(const int i) const {return _bonderi[i];}	//return the ith border. Note:  _bonderi[0] = 0, _bondery[categories()] = infinite

	virtual Vdouble getBorders() const {return _bonderi;}	
	virtual Vdouble getRates() const {return _rates;}	

protected:	
	virtual void fill_mean();
	virtual void fill_bonderi();
	
	
protected:
	MDOUBLE _alpha;
	MDOUBLE _beta;
	
	vector<MDOUBLE> _rates;
	vector<MDOUBLE> _ratesProb;
	MDOUBLE _globalRate;
	vector<MDOUBLE> _bonderi; //Note: _bonderi[0] = 0, _bondery[categories()] = infinite
};



#endif

