// 	$Id: bestAlphaAndNu.h 1975 2007-04-22 13:47:28Z privmane $	
#ifndef ___BEST_ALPHA_AND_NU
#define ___BEST_ALPHA_AND_NU

#include "definitions.h"

#include "sequenceContainer.h"
#include "stochasticProcess.h"
#include "gammaDistribution.h"
#include "tree.h"
#include "replacementModelSSRV.h"
#include "tamura92.h"
#include "stochasticProcessSSRV.h"
#include "C_evalParamUSSRV.h"
#include "bestAlpha.h"
#include "numRec.h"
#include "bblEM.h"
#include "logFile.h"

// ******************
// *     USSRV      *
// ******************

// Nu is fixed. The tree is fixed
class bestAlphaFixedTreeUSSRV {
public:
	explicit bestAlphaFixedTreeUSSRV() {}
	MDOUBLE operator()(const tree& et,
					   const sequenceContainer& sc,
					   const sequenceContainer& baseSc,
					   ussrvModel& model,
					   const Vdouble * weights=NULL,
					   const MDOUBLE upperBoundOnAlpha = 15,
					   const MDOUBLE epsilonAlphaOptimization = 0.01);
	MDOUBLE getBestAlpha() {return _bestAlpha;}
	MDOUBLE getBestL() {return _bestL;}
	
	void setAlpha(MDOUBLE alpha, ussrvModel& model) const  
	{
		model.updateAlpha(alpha);
	}

	void setBestL(MDOUBLE bestL) { _bestL = bestL;} 

private:
	MDOUBLE _bestAlpha;
	MDOUBLE _bestL;
};

// Alpha is fixed
class bestNuFixedTreeUSSRV {
public:
	explicit bestNuFixedTreeUSSRV(){}
	MDOUBLE operator()(const tree& et,
					   const sequenceContainer& sc,
					   const sequenceContainer& baseSc,
					   ussrvModel& model,
					   const Vdouble * weights=NULL,
					   const MDOUBLE upperBoundOnNu = 15,
					   const MDOUBLE epsilonNuOptimization = 0.01);
	MDOUBLE getBestNu() {return _bestNu;}
	MDOUBLE getBestL() {return _bestL;}
	void setNu(MDOUBLE nu, ussrvModel& model) const
	{
		model.updateNu(nu);
	}
	void setBestL(MDOUBLE bestL) { _bestL = bestL;}

private:
	MDOUBLE _bestNu;
	MDOUBLE _bestL;
};

class bestFFixedTreeUSSRV {
public:
	explicit bestFFixedTreeUSSRV() {}
	MDOUBLE operator()(const tree& et,
					   const sequenceContainer& sc,
					   const sequenceContainer& baseSc,
					   ussrvModel& model,
					   const Vdouble * weights=NULL,
					   const MDOUBLE upperBoundOnF = 1,
					   const MDOUBLE epsilonFOptimization = 0.01);
	MDOUBLE getBestF() {return _bestF;}
	MDOUBLE getBestL() {return _bestL;}
	void setF(MDOUBLE f, ussrvModel& model) const
	{
		if ( (f>1) || (f < 0))
		{
			LOG(5,<<"bestFFixedTreeSSRV:setF, f must be between 0 to 1. f = " << f << endl);
			return;
		}
		model.updateF(f);
	}
	void setBestL(MDOUBLE bestL) { _bestL = bestL;}

private:
	MDOUBLE _bestF;
	MDOUBLE _bestL;
};


// ******************
// *     SSRV       *
// ******************

// Nu is fixed. The tree is fixed
class bestAlphaFixedTreeSSRV {
public:
	explicit bestAlphaFixedTreeSSRV() {}
	MDOUBLE operator()(const tree& et,
		const sequenceContainer& sc,
		stochasticProcessSSRV& ssrvSp,
		const Vdouble * weights=NULL,
		const MDOUBLE lowerBoundOnAlpha = 0,
		const MDOUBLE upperBoundOnAlpha = 10,
		const MDOUBLE epsilonAlphaOptimization = 0.01);
	MDOUBLE getBestAlpha() {return _bestAlpha;}
	MDOUBLE getBestL() {return _bestL;}

	void setAlpha(MDOUBLE alpha, stochasticProcessSSRV& ssrvSp) const  
	{
		if (alpha<0) 
			errorMsg::reportError("bestAlphaFixedTreeSSRV::setAlpha, alpha is < 0 ");
		
		replacementModelSSRV* pMulRM = static_cast<replacementModelSSRV*>(ssrvSp.getPijAccelerator()->getReplacementModel());
		gammaDistribution* gammaDist = static_cast<gammaDistribution*>(pMulRM->getDistribution()); 
		gammaDist->setAlpha(alpha);
		pMulRM->updateQ();
	}

	void setBestL(MDOUBLE bestL) { _bestL = bestL;} 

private:
	MDOUBLE _bestAlpha;
	MDOUBLE _bestL;
};

// Alpha is fixed
class bestNuFixedTreeSSRV {
public:
	explicit bestNuFixedTreeSSRV(){}
	MDOUBLE operator()(const tree& et,
		const sequenceContainer& sc,
		stochasticProcessSSRV& ssrvSp,
		const Vdouble * weights=NULL,
		const MDOUBLE lowerBoundOnNu = 0,
		const MDOUBLE upperBoundOnNu = 15,
		const MDOUBLE epsilonNuOptimization = 0.01);
	MDOUBLE getBestNu() {return _bestNu;}
	MDOUBLE getBestL() {return _bestL;}
	void setNu(MDOUBLE nu, stochasticProcessSSRV& ssrvSp) const
	{
		if (nu<0) 
			errorMsg::reportError("ussrvModel::updateNu , nu is < 0");
		
		static_cast<replacementModelSSRV*>(ssrvSp.getPijAccelerator()->getReplacementModel())->setRateOfRate(nu);
	}

	void setBestL(MDOUBLE bestL) { _bestL = bestL;}

private:
	MDOUBLE _bestNu;
	MDOUBLE _bestL;
};


class bestTamura92ParamFixedTreeSSRV {
public:
	explicit bestTamura92ParamFixedTreeSSRV(){}
	MDOUBLE operator()(const tree& et,
		const sequenceContainer& sc,
		stochasticProcessSSRV& ssrvSp,
		const Vdouble * weights=NULL,
		const int maxTotalIterations = 5,
		const MDOUBLE epsilonLikelihoodImprovment = 0.05,
		const MDOUBLE lowerBoundOnTrTv =  0.0,
		const MDOUBLE upperBoundOnTrTv = 10.0,
		const MDOUBLE lowerBoundOnTheta = 0.0,
		const MDOUBLE upperBoundOnTheta = 1.0,
		const MDOUBLE epsilonTrTvOptimization = 0.01,
		const MDOUBLE epsilonThetaOptimization = 0.01);
	MDOUBLE getBestTrTv() {return _bestTrTv;}
	MDOUBLE getBestTheta() {return _bestTheta;}
	MDOUBLE getBestL() {return _bestL;}
	void setTrTv(MDOUBLE TrTv, stochasticProcessSSRV& ssrvSp) const  {
		replacementModelSSRV* pMulRM = static_cast<replacementModelSSRV*>(ssrvSp.getPijAccelerator()->getReplacementModel());
		static_cast<tamura92*>(pMulRM->getBaseRM())->changeTrTv(TrTv);
		pMulRM->updateQ();
	}

	void setTheta(MDOUBLE theta, stochasticProcessSSRV& ssrvSp) const  {
		replacementModelSSRV* pMulRM = static_cast<replacementModelSSRV*>(ssrvSp.getPijAccelerator()->getReplacementModel());
		static_cast<tamura92*>(pMulRM->getBaseRM())->changeTheta(theta);
		pMulRM->updateFreq();
		pMulRM->updateQ();
	}

	void setTrTvAndTheta(MDOUBLE TrTv, MDOUBLE theta, stochasticProcessSSRV& ssrvSp) {
		replacementModelSSRV* pMulRM = static_cast<replacementModelSSRV*>(ssrvSp.getPijAccelerator()->getReplacementModel());
		tamura92* tamuraRM = static_cast<tamura92*>(pMulRM->getBaseRM());
		tamuraRM->changeTrTv(TrTv);
		tamuraRM->changeTheta(theta);
		pMulRM->updateFreq();
		pMulRM->updateQ();
	}

private:
	MDOUBLE _bestTrTv;
	MDOUBLE _bestTheta;
	MDOUBLE _bestL;
};


#endif // ___BEST_ALPHA_AND_NU
