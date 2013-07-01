// $Id: betaOmegaDistribution.h 962 2006-11-07 15:13:34Z privmane $

#ifndef ___BETA_OMEGA_DIST
#define ___BETA_OMEGA_DIST
/************************************************************
This distribution can take several forms depending on its free parameters alpha,beta 
For an extensive exlpanation of this distribution
see http://mathworld.wolfram.com/BetaDistribution.html
************************************************************/
#include "definitions.h"
#include "distribution.h"
#include "betaDistribution.h"

#include "logFile.h"

using namespace std;


class betaOmegaDistribution : public distribution {

public:
	explicit betaOmegaDistribution(MDOUBLE alpha, MDOUBLE beta, int in_number_of_categories,MDOUBLE betaProb,MDOUBLE omega);
	explicit betaOmegaDistribution(const betaOmegaDistribution& other);
	explicit betaOmegaDistribution();
	virtual ~betaOmegaDistribution();
	virtual void setBetaOmegaParameters(int in_number_of_categories,MDOUBLE alpha, MDOUBLE beta,MDOUBLE betaProb,MDOUBLE omega);
	virtual void setBetaParameters(int numOfCategories ,MDOUBLE alpha, MDOUBLE beta){_betaDistr.setBetaParameters(numOfCategories,alpha,beta);}

	virtual const int categories() const {return _betaDistr.categories()+1;}
	virtual const MDOUBLE rates(const int i) const;
	virtual const MDOUBLE ratesProb(const int i) const;
	virtual distribution* clone() const { return new betaOmegaDistribution(*this); }
	virtual void setGlobalRate(const MDOUBLE x) {_betaDistr.setGlobalRate(x);}
	virtual MDOUBLE getGlobalRate()const {return  _betaDistr.getGlobalRate();}
	virtual const MDOUBLE getCumulativeProb(const MDOUBLE x) const;
	virtual void setAlpha(MDOUBLE newAlpha){ _betaDistr.setAlpha(newAlpha);}
	virtual MDOUBLE getAlpha() const {return _betaDistr.getAlpha();};
	virtual void setBeta(MDOUBLE newBeta){_betaDistr.setBeta(newBeta);}
	virtual MDOUBLE getBeta() const {return _betaDistr.getBeta();};
	virtual void change_number_of_categories(int in_number_of_categories){_betaDistr.change_number_of_categories(in_number_of_categories);}
	virtual MDOUBLE getBorder(const int i) const {return _betaDistr.getBorder(i);}	//return the ith border. Note:  _bonderi[0] = 0, _bondery[categories()] = infinite
	virtual MDOUBLE getOmega() const {return _omega;}
	virtual MDOUBLE getBetaProb() const {return _betaProb;};
	virtual void setOmega(MDOUBLE omega) { _omega = omega;};
	virtual void setBetaProb(MDOUBLE betaProb) { _betaProb = betaProb;};

private:	
	betaDistribution _betaDistr;
	MDOUBLE _omega;
	MDOUBLE _betaProb;
};



#endif

