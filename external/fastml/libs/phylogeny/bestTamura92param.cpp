// $Id: bestTamura92param.cpp 962 2006-11-07 15:13:34Z privmane $

#include "bestTamura92param.h"
#include <iostream>
using namespace std;

#include "bblEM.h"
#include "numRec.h"
#include "logFile.h"
#include "bestAlpha.h"

bestTamura92ParamFixedTree::bestTamura92ParamFixedTree(const tree& et, // find best TrTv and theta
													   const sequenceContainer& sc,
													   stochasticProcess& sp,
													   const Vdouble * weights,
													   const int maxTotalIterations,
													   const MDOUBLE epsilonLikelihoodImprovment,
													   const MDOUBLE epsilonLoglikelihoodForTrTvOptimization,
													   const MDOUBLE epsilonLoglikelihoodForThetaOptimization,
													   const MDOUBLE upperBoundOnTrTv) {
	LOG(5,<<"Starting bestTamura92ParamFixedTree: find Best TrTv and theta"<<endl);
	MDOUBLE oldL = VERYSMALL;
	MDOUBLE newL = VERYSMALL;

	// first guess for the parameters
	MDOUBLE prevTrTv = upperBoundOnTrTv*0.3;
	MDOUBLE prevTheta = 0.5;

	for (int i=0; i < maxTotalIterations; ++i) {
		// optimize TrTv
		newL = -brent(0.0, prevTrTv, upperBoundOnTrTv,
					  C_evalTrTvParam(et,sc,sp,weights),
					  epsilonLoglikelihoodForTrTvOptimization,
					  &_bestTrTv);

		// optimize Theta
		newL = -brent(0.0, prevTheta, 1.0,
					  C_evalTheta(et,sc,sp,weights),
					  epsilonLoglikelihoodForThetaOptimization,
					  &_bestTheta);

		// check for improvement in the likelihood
		if (newL > oldL+epsilonLikelihoodImprovment) {
			prevTrTv = _bestTrTv;
			prevTheta = _bestTheta;
			oldL = newL;
			_bestL = newL;
		} else {
			if (newL>oldL) {
				_bestL = newL;
			} else {
				_bestL = oldL;
				_bestTrTv = prevTrTv;
				_bestTheta = prevTheta;
			}
			break;
		}
	}
}

bestTamura92ParamAndBBL::bestTamura92ParamAndBBL(tree& et, //find best TrTv, theta and best BBL
												 const sequenceContainer& sc,
												 stochasticProcess& sp,
												 const Vdouble * weights,
												 const int maxTotalIterations,
												 const MDOUBLE epsilonLikelihoodImprovment,
												 const MDOUBLE epsilonLoglikelihoodForTrTvOptimization,
												 const MDOUBLE epsilonLoglikelihoodForThetaOptimization,
												 const MDOUBLE epsilonLoglikelihoodForBBL,
												 const MDOUBLE upperBoundOnTrTv,
												 const int maxBBLIterations){
	LOG(5,<<"Starting bestTamura92ParamAndBBL: find best TrTv, theta and BBL"<<endl);
	MDOUBLE oldL = VERYSMALL;
	MDOUBLE newL = VERYSMALL;

	// first guess for the parameters
	MDOUBLE prevTrTv = upperBoundOnTrTv*0.3;
	MDOUBLE prevTheta = 0.5;
	tree prevTree;

	for (int i=0; i < maxTotalIterations; ++i) {
		// optimize TrTv
		newL = -brent(0.0, prevTrTv, upperBoundOnTrTv,
					  C_evalTrTvParam(et,sc,sp,weights),
					  epsilonLoglikelihoodForTrTvOptimization,
					  &_bestTrTv);
		(static_cast<tamura92*>(sp.getPijAccelerator()->getReplacementModel()))->changeTrTv(_bestTrTv);

		// optimize Theta
		newL = -brent(0.0, prevTheta, 1.0,
					  C_evalTheta(et,sc,sp,weights),
					  epsilonLoglikelihoodForThetaOptimization,
					  &_bestTheta);
		(static_cast<tamura92*>(sp.getPijAccelerator()->getReplacementModel()))->changeTheta(_bestTheta);

		// optimize branch lengths
		bblEM bblEM1(et,sc,sp,NULL,maxBBLIterations,epsilonLoglikelihoodForBBL);//maxIterations=1000
		newL =bblEM1.getTreeLikelihood();

		// check for improvement in the likelihood
		if (newL > oldL+epsilonLikelihoodImprovment) {
			prevTrTv = _bestTrTv;
			prevTheta = _bestTheta;
			oldL = newL;
			_bestL = newL;
			prevTree = et;
		} else {
			if (newL>oldL) {
				_bestL = newL;
			} else {
				_bestL = oldL;
				_bestTrTv = prevTrTv;
				_bestTheta = prevTheta;
				et = prevTree;
			}
			break;
		}
	}
}
		
bestTamura92ParamAlphaAndBBL::bestTamura92ParamAlphaAndBBL( //find best TrTv, theta, Alpha and best branch lengths
	tree& et,
	const sequenceContainer& sc,
	stochasticProcess& sp,
	const Vdouble * weights,
	const int maxTotalIterations,
	const MDOUBLE epsilonLikelihoodImprovment,
	const MDOUBLE epsilonLoglikelihoodForTrTvOptimization,
	const MDOUBLE epsilonLoglikelihoodForThetaOptimization,
	const MDOUBLE epsilonLoglikelihoodForAlphaOptimization,
	const MDOUBLE epsilonLoglikelihoodForBBL,
	const MDOUBLE upperBoundOnTrTv,
	const int maxBBLIterations,
	const MDOUBLE initAlpha,
	const MDOUBLE upperBoundOnAlpha)

{
	MDOUBLE oldL = VERYSMALL;
	MDOUBLE newL = VERYSMALL;

	// first guess for the parameters
	MDOUBLE prevTrTv = static_cast<tamura92*>(sp.getPijAccelerator()->getReplacementModel())->getTrTv();
	MDOUBLE prevTheta = static_cast<tamura92*>(sp.getPijAccelerator()->getReplacementModel())->getTheta();
	MDOUBLE prevAlpha = initAlpha;
	tree prevTree;

	for (int i=0; i < maxTotalIterations; ++i) {

		// optimize TrTv
		newL = -brent(0.0, prevTrTv, upperBoundOnTrTv,
					  C_evalTrTvParam(et,sc,sp,weights),
					  epsilonLoglikelihoodForTrTvOptimization,
					  &_bestTrTv);
		(static_cast<tamura92*>(sp.getPijAccelerator()->getReplacementModel()))->changeTrTv(_bestTrTv);

		// optimize Theta
		newL = -brent(0.0, prevTheta, 1.0,
					  C_evalTheta(et,sc,sp,weights),
					  epsilonLoglikelihoodForThetaOptimization,
					  &_bestTheta);
		(static_cast<tamura92*>(sp.getPijAccelerator()->getReplacementModel()))->changeTheta(_bestTheta);

		// optimize Alpha
		newL = -brent(0.0, prevAlpha, upperBoundOnAlpha,
					  C_evalAlpha(et,sc,sp,weights),
					  epsilonLoglikelihoodForAlphaOptimization,
					  &_bestAlpha);
		(static_cast<gammaDistribution*>(sp.distr()))->setAlpha(_bestAlpha);
 
		LOG(5,<<"# bestTamura92ParamAlphaAndBBL::bestTamura92ParamAlphaAndBBL iteration " << i << ": after param optimization:" <<endl
		      <<"# old L = " << oldL << "\t"
		      <<"# new L = " << newL << endl
		      <<"# new Alpha = " << _bestAlpha << endl);

		// optimize branch lengths
		bblEM bblEM1(et,sc,sp,NULL,maxBBLIterations,epsilonLoglikelihoodForBBL);//maxIterations=1000
		newL =bblEM1.getTreeLikelihood();

		LOG(5,<<"# bestTamura92ParamAlphaAndBBL::bestTamura92ParamAlphaAndBBL iteration " << i << ": after branch lengths optimization:" <<endl 
		      <<"# After BBL new L = "<<newL<<" old L = "<<oldL<<endl
		      <<"# The tree:" );
		LOGDO(5,et.output(myLog::LogFile()));

		// check for improvement in the likelihood
		if (newL > oldL+epsilonLikelihoodImprovment) {
		    oldL = newL;
			_bestL = newL;
			prevTrTv = _bestTrTv;
			prevTheta = _bestTheta;
			prevAlpha = _bestAlpha;
			prevTree = et;
		} else {
			if (newL>oldL) {
				_bestL = newL;
			} else {
				_bestL = oldL;
				_bestTrTv = prevTrTv;
				_bestTheta = prevTheta;
				et = prevTree;
			}
		    break;
		}
	}
}

