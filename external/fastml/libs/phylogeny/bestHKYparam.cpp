// $Id: bestHKYparam.cpp 4314 2008-06-25 13:09:12Z itaymay $

#include "bestHKYparam.h"
#include <iostream>
using namespace std;

#include "bblEM.h"
#include "numRec.h"
#include "logFile.h"
#include "bestAlpha.h"

bestHkyParamAndBBL::bestHkyParamAndBBL(tree& et, //find Best HkyParam and best BBL
					   const sequenceContainer& sc,
					   stochasticProcess& sp,
					   const Vdouble * weights,
					   const MDOUBLE upperBoundOnHkyParam,
					   const MDOUBLE epsilonHkyParamOptimization,
					   const MDOUBLE epsilonLikelihoodImprovment,
					   const int maxBBLIterations,
					   const int maxTotalIterations){
	LOG(5,<<"find Best HkyParam and best BBL"<<endl);
//	LOG(5,<<" 1. bestHkyParam::findBestHkyParam"<<endl);
//	brLenOpt br1(*et,*pi,weights);
	MDOUBLE oldL = VERYSMALL;
	_bestL = VERYSMALL;
	const MDOUBLE bx=upperBoundOnHkyParam*0.3;
	const MDOUBLE ax=0.01;
	const MDOUBLE cx=upperBoundOnHkyParam;
	MDOUBLE bestA=0;
	for (int i=0; i < maxTotalIterations; ++i) {
		_bestL = -brent(ax,bx,cx,
		C_evalHkyParam(et,sc,sp,weights),
		epsilonHkyParamOptimization,
		&bestA);

		if (_bestL > oldL+epsilonLikelihoodImprovment) {
			oldL = _bestL;
		} 
		else {//LL converged
			if (_bestL > oldL)
				(static_cast<hky*>(sp.getPijAccelerator()->getReplacementModel()))->changeTrTv(bestA);
			else
                _bestL = oldL;
            break;
		}
		_bestHkyParam = bestA;
		(static_cast<hky*>(sp.getPijAccelerator()->getReplacementModel()))->changeTrTv(bestA);
		LOG(5,<<"bestHkyParamAndBBL: trtv = "<<_bestHkyParam<<endl);
		bblEM bblEM1(et,sc,sp,NULL,maxBBLIterations,epsilonLikelihoodImprovment);//maxIterations=1000
		_bestL =bblEM1.getTreeLikelihood();
		if (_bestL > oldL+epsilonLikelihoodImprovment) {
			oldL = _bestL;
		}
		else {
			_bestL = oldL;
			break;
		}
	}
}
		
bestHkyParamFixedTree::bestHkyParamFixedTree(const tree& et, //findBestHkyParamFixedTree
					   const sequenceContainer& sc,
					   stochasticProcess& sp,
					   const Vdouble * weights,
					   const MDOUBLE upperBoundOnHkyParam,
					   const MDOUBLE epsilonHkyParamOptimization){
	LOG(5,<<"findBestHkyParamFixedTree"<<endl);
	MDOUBLE bestA=0;
	const MDOUBLE cx=upperBoundOnHkyParam;// left, midle, right limit on HkyParam
	const MDOUBLE bx=cx*0.3;
	const MDOUBLE ax=0;

	
	_bestL = -brent(ax,bx,cx,
		C_evalHkyParam(et,sc,sp,weights),
		epsilonHkyParamOptimization,
		&bestA);
	_bestHkyParam= bestA;
	(static_cast<hky*>(sp.getPijAccelerator()->getReplacementModel()))->changeTrTv(bestA);
}



bestHkyParamAlphaAndBBL::bestHkyParamAlphaAndBBL( //find best TrTv (=HkyParam), Alpha and best branch lengths
	tree& et,
	const sequenceContainer& sc,
	stochasticProcess& sp,
	const Vdouble * weights,
	const int maxTotalIterations,
	const MDOUBLE epsilonLikelihoodImprovment,
	const MDOUBLE epsilonHkyParamOptimization,
	const MDOUBLE epsilonAlphaOptimization,
	const MDOUBLE epsilonBBL,
	const MDOUBLE upperBoundOnHkyParam,
	const int maxBBLIterations,
	const MDOUBLE initAlpha,
	const MDOUBLE upperBoundOnAlpha)

{
	MDOUBLE oldL = VERYSMALL;
	MDOUBLE newL = VERYSMALL;

	// first guess for the parameters
	MDOUBLE prevHkyParam = static_cast<hky*>(sp.getPijAccelerator()->getReplacementModel())->getTrTv();
	MDOUBLE prevAlpha = initAlpha;
	tree prevTree;

	for (int i=0; i < maxTotalIterations; ++i) {

		// optimize HkyParam
		newL = -brent(0.0, prevHkyParam, upperBoundOnHkyParam,
					  C_evalHkyParam(et,sc,sp,weights),
					  epsilonHkyParamOptimization,
					  &_bestHkyParam);
		(static_cast<hky*>(sp.getPijAccelerator()->getReplacementModel()))->changeTrTv(_bestHkyParam);
		LOG(5,<<"bestHkyParamAlphaAndBBL: trtv = "<<_bestHkyParam<<endl);
		// optimize Alpha
		newL = -brent(0.0, prevAlpha, upperBoundOnAlpha,
					  C_evalAlpha(et,sc,sp,weights),
					  epsilonAlphaOptimization,
					  &_bestAlpha);
		(static_cast<gammaDistribution*>(sp.distr()))->setAlpha(_bestAlpha);
		
		LOG(5,<<"# bestHkyParamAlphaAndBBL::bestHkyParamAlphaAndBBL iteration " << i << ": after param optimization:" <<endl
		      <<"# old L = " << oldL << "\t"
		      <<"# new L = " << newL << endl
			  <<"# new hkyParam = " << _bestHkyParam << endl
		      <<"# new Alpha = " << _bestAlpha << endl);

		// optimize branch lengths
		bblEM bblEM1(et,sc,sp,NULL,maxBBLIterations,epsilonBBL);
		newL =bblEM1.getTreeLikelihood();

		LOG(5,<<"# bestHkyParamAlphaAndBBL::bestHkyParamAlphaAndBBL iteration " << i << ": after branch lengths optimization:" <<endl 
		      <<"# After BBL new L = "<<newL<<" old L = "<<oldL<<endl
		      <<"# The tree:" );
		LOGDO(5,et.output(myLog::LogFile()));

		// check for improvement in the likelihood
		if (newL > oldL+epsilonLikelihoodImprovment) {
		    oldL = newL;
			_bestL = newL;
			prevHkyParam = _bestHkyParam;
			prevAlpha = _bestAlpha;
			prevTree = et;
		} else {
			if (newL>oldL) {
				_bestL = newL;
			} else {
				_bestL = oldL;
				_bestHkyParam = prevHkyParam;
				et = prevTree;
			}
		    break;
		}
	}
}

