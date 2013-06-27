// $Id: bestTamura92param.h 962 2006-11-07 15:13:34Z privmane $

#ifndef ___BEST_TAMURA92_PARAM
#define ___BEST_TAMURA92_PARAM

#include "definitions.h"

#include "likelihoodComputation.h"
#include "sequenceContainer.h"
#include "stochasticProcess.h"
#include "gammaDistribution.h"
#include "tree.h"
#include "tamura92.h"


class bestTamura92ParamFixedTree {
public:
	explicit bestTamura92ParamFixedTree(const tree& et, // find best TrTv and theta
										const sequenceContainer& sc,
										stochasticProcess& sp,
										const Vdouble * weights,
										const int maxTotalIterations = 5,
										const MDOUBLE epsilonLikelihoodImprovment = 0.05,
										const MDOUBLE epsilonLoglikelihoodForTrTvOptimization = 0.01,
										const MDOUBLE epsilonLoglikelihoodForThetaOptimization = 0.01,
										const MDOUBLE upperBoundOnTrTv = 5.0);
	MDOUBLE getBestTrTv() {return _bestTrTv;}
	MDOUBLE getBestTheta() {return _bestTheta;}
	MDOUBLE getBestL() {return _bestL;}
private:
	MDOUBLE _bestTrTv;
	MDOUBLE _bestTheta;
	MDOUBLE _bestL;
};

class bestTamura92ParamAndBBL {
public:
	explicit bestTamura92ParamAndBBL(tree& et, //find best TrTv, theta and best BBL
									 const sequenceContainer& sc,
									 stochasticProcess& sp,
									 const Vdouble * weights=NULL,
									 const int maxTotalIterations=5,
									 const MDOUBLE epsilonLikelihoodImprovment = 0.05,
									 const MDOUBLE epsilonLoglikelihoodForTrTvOptimization = 0.01,
									 const MDOUBLE epsilonLoglikelihoodForThetaOptimization = 0.01,
									 const MDOUBLE epsilonLoglikelihoodForBBL = 0.01,
									 const MDOUBLE upperBoundOnTrTv = 5.0,
									 const int maxBBLIterations=10);
	MDOUBLE getBestTrTv() {return _bestTrTv;}
	MDOUBLE getBestTheta() {return _bestTheta;}
	MDOUBLE getBestL() {return _bestL;}
private:
	MDOUBLE _bestTrTv;
	MDOUBLE _bestTheta;
	MDOUBLE _bestL;
};

class bestTamura92ParamAlphaAndBBL {
public:
	explicit bestTamura92ParamAlphaAndBBL( //find best TrTv, theta, Alpha and best branch lengths
		tree& et,
		const sequenceContainer& sc,
		stochasticProcess& sp,
		const Vdouble * weights=NULL,
		const int maxTotalIterations=5,
		const MDOUBLE epsilonLikelihoodImprovment= 0.05,
		const MDOUBLE epsilonLoglikelihoodForTrTvOptimization= 0.01,
		const MDOUBLE epsilonLoglikelihoodForThetaOptimization= 0.01,
		const MDOUBLE epsilonLoglikelihoodForAlphaOptimization= 0.01,
		const MDOUBLE epsilonLoglikelihoodForBBL= 0.01,
		const MDOUBLE upperBoundOnTrTv = 5.0,
		const int maxBBLIterations=10,
		const MDOUBLE initAlpha = 1.5,
		const MDOUBLE upperBoundOnAlpha = 5.0);
	MDOUBLE getBestTrTv() {return _bestTrTv;}
	MDOUBLE getBestTheta() {return _bestTheta;}
	MDOUBLE getBestAlpha() {return _bestAlpha;}
	MDOUBLE getBestL() {return _bestL;}
private:
	MDOUBLE _bestTrTv;
	MDOUBLE _bestTheta;
	MDOUBLE _bestAlpha;
	MDOUBLE _bestL;
};



class C_evalTrTvParam{
public:
  C_evalTrTvParam( const tree& et,
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
	MDOUBLE operator() (MDOUBLE TrTv) {
		(static_cast<tamura92*>(_sp.getPijAccelerator()->getReplacementModel()))->changeTrTv(TrTv);
		
		MDOUBLE res = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(_et,_sc,_sp,_weights);
		LOG(5,<<" with TrTv = "<<TrTv<<" logL = "<<res<<endl);
		return -res;
	}
};

class C_evalTheta{
public:
  C_evalTheta(	const tree& et,
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
	MDOUBLE operator() (MDOUBLE theta) {
		(static_cast<tamura92*>(_sp.getPijAccelerator()->getReplacementModel()))->changeTheta(theta);
		
		MDOUBLE res = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(_et,_sc,_sp,_weights);
		LOG(5,<<" with theta = "<<theta<<" logL = "<<res<<endl);
		return -res;
	}
};




#endif


