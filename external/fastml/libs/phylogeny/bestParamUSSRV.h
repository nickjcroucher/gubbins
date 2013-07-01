// 	$Id: bestParamUSSRV.h 1975 2007-04-22 13:47:28Z privmane $	
#ifndef BEST_PARAM_USSRV
#define BEST_PARAM_USSRV

#include "definitions.h"
#include "sequenceContainer.h"
#include "stochasticProcess.h"
#include "gammaDistribution.h"
#include "tree.h"
#include "replacementModelSSRV.h"
#include "stochasticProcessSSRV.h"
#include "C_evalParamUSSRV.h"
#include "bestAlpha.h"
#include "numRec.h"
#include "bblEM.h"
#include "logFile.h"
#include "bestAlphaAndNu.h"
#include "bblEM2USSRV.h"
#include "someUtil.h"
#include <ctime> 

// ***************
// *    USSRV    *
// ***************

class bestParamUSSRV
{
public:
	explicit bestParamUSSRV(bool AlphaOptimization, bool NuOptimization,
							bool FOptimization, bool bblOptimization):
							_AlphaOptimizationFlag(AlphaOptimization),
							_NuOptimizationFlag(NuOptimization),
							_FOptimizationFlag(FOptimization),
							_bblOptimizationFlag(bblOptimization) {}
	
	MDOUBLE operator() (tree& et,
					   const sequenceContainer& sc,
					   const sequenceContainer& baseSc,
					   ussrvModel& model,
					   const Vdouble * weights=NULL,
					   const MDOUBLE AlphaUpperBound = 15, 
					   const MDOUBLE NuUpperBound = 15, 
					   const MDOUBLE FUpperBound = 1, 
					   const MDOUBLE epsilonParamOptimization = 0.01,
					   const MDOUBLE epsilonLikelihoodImprovment = 0.01,
					   const int maxIterations = 50,
					   const int maxOfParametersAndBblIterations = 40);

	MDOUBLE getBestAlpha() {return _bestAlpha;}
	MDOUBLE getBestNu() {return _bestNu;}
	MDOUBLE getBestF() {return _bestF;}
	MDOUBLE getBestL() {return _bestL;}

private:
	MDOUBLE _bestAlpha;
	MDOUBLE _bestNu;
	MDOUBLE _bestF;
	MDOUBLE _bestL;

	// flags
	bool _AlphaOptimizationFlag;
	bool _NuOptimizationFlag;
    bool _FOptimizationFlag;
	bool _bblOptimizationFlag;
};

// ***************
// *     SSRV    *
// ***************

class bestParamSSRV
{
public:
	explicit bestParamSSRV(bool AlphaOptimization, bool NuOptimization, bool tamura92Optimization,
						   bool bblOptimization):
		_AlphaOptimizationFlag(AlphaOptimization),
		_NuOptimizationFlag(NuOptimization),
		_tamura92OptimizationFlag(tamura92Optimization),
		_bblOptimizationFlag(bblOptimization) {}

	MDOUBLE operator() (tree& et,
						const sequenceContainer& sc,
						stochasticProcessSSRV& ssrvSp,
						const Vdouble * weights=NULL,
						const MDOUBLE AlphaUpperBound = 15, 
						const MDOUBLE NuUpperBound = 15,
						const MDOUBLE TrTvUpperBound = 10,
						const MDOUBLE epsilonParamOptimization = 0.01,
						const MDOUBLE epsilonLikelihoodImprovment = 0.01,
						const MDOUBLE epsilonBbl = 0.05,
						const int maxIterations = 50,
						const int maxOfParametersAndBblIterations = 40);

	// Variant that can work on a const tree - only if we're not doing BBL
	// WARNING: Running this with bblOptimization==true will give a fatal error
	MDOUBLE operator() (const tree& et,
						const sequenceContainer& sc,
						stochasticProcessSSRV& ssrvSp,
						const Vdouble * weights=NULL,
						const MDOUBLE AlphaUpperBound = 15, 
						const MDOUBLE NuUpperBound = 15,
						const MDOUBLE TrTvUpperBound = 10,
						const MDOUBLE epsilonParamOptimization = 0.01,
						const MDOUBLE epsilonLikelihoodImprovment = 0.01,
						const MDOUBLE epsilonBbl = 0.05,
						const int maxIterations = 50,
						const int maxOfParametersAndBblIterations = 40);

	MDOUBLE getBestAlpha() {return _bestAlpha;}
	MDOUBLE getBestNu() {return _bestNu;}
	MDOUBLE getBestTrTv() {return _bestTrTv;}
	MDOUBLE getBestTheta() {return _bestTheta;}
	MDOUBLE getBestL() {return _bestL;}

private:
	MDOUBLE _bestAlpha;
	MDOUBLE _bestNu;
	MDOUBLE _bestTrTv;
	MDOUBLE _bestTheta;
	MDOUBLE _bestL;

	// flags
	bool _AlphaOptimizationFlag;
	bool _NuOptimizationFlag;
	bool _tamura92OptimizationFlag;
	bool _bblOptimizationFlag;
};

#endif // BEST_PARAM_USSRV

