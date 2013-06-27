#ifndef ___GAMMA_DIST_PLUSINV
#define ___GAMMA_DIST_PLUSINV
/************************************************************
This class describes a combination of a predefined dsitrubtion ,
with an additional invariant category of probability _Pinv
This category is always the last rate category (i.e., rate(categories()) == 0)
************************************************************/
#include "definitions.h"
#include "distributionPlusInvariant.h"
#include "distribution.h"
#include "gammaDistribution.h"
#include "errorMsg.h"
#include "gammaUtilities.h"
#include "logFile.h"
#include <cmath>



class gammaDistributionPlusInvariant : public distributionPlusInvariant {
public:
	explicit gammaDistributionPlusInvariant(distribution* pDist, const MDOUBLE pInv, const MDOUBLE globalRate=1, MDOUBLE rateInvariantVal=1e-10): distributionPlusInvariant(pDist,pInv,globalRate,rateInvariantVal){}
	explicit gammaDistributionPlusInvariant();
	gammaDistributionPlusInvariant(const gammaDistributionPlusInvariant& other) {(*this) = other;}	
	//virtual gammaDistributionPlusInvariant& operator=(const gammaDistributionPlusInvariant& other);
	gammaDistributionPlusInvariant* clone() const {return new gammaDistributionPlusInvariant(*this);}
	virtual ~gammaDistributionPlusInvariant(){}



// get GammaDistribution params
	virtual void setAlpha(MDOUBLE newAlpha) {return static_cast<gammaDistribution*>(_pBaseDist)->setAlpha(newAlpha);};
	virtual MDOUBLE getAlpha() const {return static_cast<gammaDistribution*>(_pBaseDist)->getAlpha();}

};
#endif
