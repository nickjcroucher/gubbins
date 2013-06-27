#ifndef __DISTPLUSINV
#define __DISTPLUSINV
/************************************************************
This class describes a combination of a predefined dsitrubtion ,
with an additional invariant category of probability _Pinv
This category is always the last rate category (i.e., rate(categories()) == 0)
************************************************************/
#include "definitions.h"
#include "distribution.h"

class distributionPlusInvariant : public distribution {
public:
	explicit distributionPlusInvariant(
		distribution* pDist, const MDOUBLE pInv, const MDOUBLE globalRate=1, MDOUBLE rateInvariantVal=1e-10);
	explicit distributionPlusInvariant();
	distributionPlusInvariant(const distributionPlusInvariant& other): _pBaseDist(NULL){(*this) = other;}	
	virtual distributionPlusInvariant& operator=(const distributionPlusInvariant& other);
	distributionPlusInvariant* clone() const {return new distributionPlusInvariant(*this);}

	virtual ~distributionPlusInvariant();

	distribution* getBaseDistribution(){return _pBaseDist;}
	//get/set the parameters of the mixture
	const int categories() const; 
 	void setGlobalRate(const MDOUBLE r) {_globalRate = r;}
 	MDOUBLE getGlobalRate() const {return _globalRate;}
	virtual void setInvProb(const MDOUBLE p) {_Pinv = p;}
	const MDOUBLE getInvProb() const {return _Pinv;}

	//get distribution statistics
	virtual const MDOUBLE getCumulativeProb(const MDOUBLE x) const;
	virtual const MDOUBLE rates(const int category) const;
	virtual const MDOUBLE ratesProb(const int i) const;

protected:
	MDOUBLE _globalRate;
	MDOUBLE _Pinv;
	MDOUBLE _rateInvariantVal;
	distribution* _pBaseDist;
};
#endif
