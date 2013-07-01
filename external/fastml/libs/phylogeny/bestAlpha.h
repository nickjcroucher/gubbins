// $Id: bestAlpha.h 5786 2009-01-19 22:22:48Z rubi $

#ifndef ___BEST_ALPHA
#define ___BEST_ALPHA

#include "definitions.h"

#include "likelihoodComputation.h"
#include "sequenceContainer.h"
#include "stochasticProcess.h"
#include "gammaDistribution.h"
#include "tree.h"
#include "logFile.h"

#ifndef VERBOS
#define VERBOS
#endif

class bestAlphaFixedTree {
public:
	explicit bestAlphaFixedTree(const tree& et,
					   const sequenceContainer& sc,
					   stochasticProcess& sp,
					   const Vdouble * weights=NULL,
					   const MDOUBLE upperBoundOnAlpha = 15,
					   const MDOUBLE epsilonAlphaOptimization = 0.01);
	MDOUBLE getBestAlpha() {return _bestAlpha;}
	MDOUBLE getBestL() {return _bestL;}
private:
	MDOUBLE _bestAlpha;
	MDOUBLE _bestL;
};

class bestAlphaAndBBL {
public:
	explicit bestAlphaAndBBL(tree& et, //find Best Alpha and best BBL
					   const sequenceContainer& sc,
					   stochasticProcess& sp,
					   const Vdouble * weights=NULL,
					   const MDOUBLE initAlpha = 1.5,
					   const MDOUBLE upperBoundOnAlpha = 5.0,
					   const MDOUBLE epsilonLoglikelihoodForAlphaOptimization= 0.01,
					   const MDOUBLE epsilonLoglikelihoodForBBL= 0.05,
					   const int maxBBLIterations=10,
					   const int maxTotalIterations=5);
	MDOUBLE getBestAlpha() {return _bestAlpha;}
	MDOUBLE getBestL() {return _bestL;}
private:
	MDOUBLE _bestAlpha;
	MDOUBLE _bestL;
};

class bestBetaAndBBL {
public:
	explicit bestBetaAndBBL(tree& et, //find Best Beta and best BBL
					   const sequenceContainer& sc,
					   stochasticProcess& sp,
					   const Vdouble * weights=NULL,
					   const MDOUBLE initBeta = 1.5,
					   const MDOUBLE upperBoundOnBeta = 5.0,
					   const MDOUBLE epsilonLoglikelihoodForBetaOptimization= 0.01,
					   const MDOUBLE epsilonLoglikelihoodForBBL= 0.05,
					   const int maxBBLIterations=10,
					   const int maxTotalIterations=5);
	MDOUBLE getBestBeta() {return _bestBeta;}
	MDOUBLE getBestL() {return _bestL;}
private:
	MDOUBLE _bestBeta;
	MDOUBLE _bestL;
};

class bestAlphaAndBetaAndBBL {
public:
	explicit bestAlphaAndBetaAndBBL(tree& et, //find Best Alpha and best BBL
					   const sequenceContainer& sc,
					   stochasticProcess& sp,
					   const Vdouble * weights=NULL,
					   const MDOUBLE initAlpha = 1.5,
					   const MDOUBLE initBeta = 1.5,
					   const MDOUBLE upperBoundOnAlpha = 5.0,
					   const MDOUBLE upperBoundOnBeta = 5.0,
					   const MDOUBLE epsilonLoglikelihoodForAlphaOptimization= 0.01,
					   const MDOUBLE epsilonLoglikelihoodForBetaOptimization = 0.01,
					   const MDOUBLE epsilonLoglikelihoodForBBL= 0.05,
					   const int maxBBLIterations=10,
					   const int maxTotalIterations=5);
	MDOUBLE getBestAlpha() {return _bestAlpha;}
	MDOUBLE getBestBeta() {return _bestBeta;}
	MDOUBLE getBestL() {return _bestL;}
private:
	MDOUBLE _bestAlpha;
	MDOUBLE _bestBeta;
	MDOUBLE _bestL;
};


class C_evalAlpha{
public:
  C_evalAlpha(	const tree& et,
				const sequenceContainer& sc,
				stochasticProcess& sp,
				const Vdouble * weights = NULL)
    : _et(et),_sc(sc),_weights(weights),_sp(sp){};
private:
	const tree& _et;
	const sequenceContainer& _sc;
	const Vdouble * _weights;
	stochasticProcess& _sp;
public:
	MDOUBLE operator() (MDOUBLE alpha) {
		if (_sp.categories() == 1) {
			errorMsg::reportError(" one category when trying to optimize alpha");
		}
		(static_cast<gammaDistribution*>(_sp.distr()))->setAlpha(alpha);
		
		MDOUBLE res = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(_et,_sc,_sp,_weights);
		//LOG(5,<<" with alpha = "<<alpha<<" logL = "<<res<<endl);
#ifdef VERBOS
		LOG(7,<<" while in brent: with alpha = "<<alpha<<" logL = "<<res<<endl);
#endif
		return -res;
	}
};

class C_evalBeta{
public:
  C_evalBeta(	const tree& et,
				const sequenceContainer& sc,
				stochasticProcess& sp,
				const Vdouble * weights = NULL)
    : _et(et),_sc(sc),_weights(weights),_sp(sp){};
private:
	const tree& _et;
	const sequenceContainer& _sc;
	const Vdouble * _weights;
	stochasticProcess& _sp;
public:
	MDOUBLE operator() (MDOUBLE beta) {
		if (_sp.categories() == 1) {
			errorMsg::reportError(" one category when trying to optimize beta");
		}
		(static_cast<generalGammaDistribution*>(_sp.distr()))->setBeta(beta);
		
		MDOUBLE res = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(_et,_sc,_sp,_weights);
		//LOG(5,<<" with alpha = "<<alpha<<" logL = "<<res<<endl);
#ifdef VERBOS
		LOG(7,<<" while in brent: with beta = "<<beta<<" logL = "<<res<<endl);
#endif
		return -res;
	}
};

#endif


