#ifndef __GENERAL_GAMMA_DIST_PLUSINV
#define __GENERAL_GAMMA_DIST_PLUSINV
/************************************************************
This class describes a combination of a predefined dsitrubtion ,
with an additional invariant category of probability _Pinv
This category is always the last rate category (i.e., rate(categories()) == 0)
************************************************************/
#include "definitions.h"
#include "distributionPlusInvariant.h"
#include "distribution.h"
#include "generalGammaDistribution.h"
#include "errorMsg.h"
#include "gammaUtilities.h"
#include "logFile.h"
#include <cmath>



class generalGammaDistributionPlusInvariant : public distributionPlusInvariant {
public:
	explicit generalGammaDistributionPlusInvariant(distribution* pDist, const MDOUBLE pInv, const MDOUBLE globalRate=1, MDOUBLE rateInvariantVal=1e-10): distributionPlusInvariant(pDist,pInv,globalRate,rateInvariantVal){}
	explicit generalGammaDistributionPlusInvariant();
	generalGammaDistributionPlusInvariant(const generalGammaDistributionPlusInvariant& other) {(*this) = other;}	
	//virtual generalGammaDistributionPlusInvariant& operator=(const generalGammaDistributionPlusInvariant& other);
	generalGammaDistributionPlusInvariant* clone() const {return new generalGammaDistributionPlusInvariant(*this);}
	virtual ~generalGammaDistributionPlusInvariant(){}

//	distribution* getBaseDistribution(){return _pBaseDist;}
////get/set the parameters of the mixture
//	const int categories() const; 
//	void setGlobalRate(const MDOUBLE r) {_globalRate = r;}
//	MDOUBLE getGlobalRate() const {return _globalRate;}
//	virtual void setInvProb(const MDOUBLE p) {_Pinv = p;}
//	const MDOUBLE getInvProb() const {return _Pinv;}
//
////get distribution statistics
//	virtual const MDOUBLE getCumulativeProb(const MDOUBLE x) const;
//	virtual const MDOUBLE rates(const int category) const;
//	virtual const MDOUBLE ratesProb(const int i) const;

// get generalGammaDistribution params
	virtual void setAlpha(MDOUBLE newAlpha) {return static_cast<generalGammaDistribution*>(_pBaseDist)->setAlpha(newAlpha);};
	virtual MDOUBLE getAlpha() const {return static_cast<generalGammaDistribution*>(_pBaseDist)->getAlpha();}
	virtual void setBeta(MDOUBLE newBeta) {return static_cast<generalGammaDistribution*>(_pBaseDist)->setBeta(newBeta);};
	virtual MDOUBLE getBeta() const {return static_cast<generalGammaDistribution*>(_pBaseDist)->getBeta();}
//protected:
	//MDOUBLE _globalRate;
	//MDOUBLE _Pinv;
	//distribution* _pBaseDist;
};
#endif
