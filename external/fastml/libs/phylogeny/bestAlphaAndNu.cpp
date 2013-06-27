// 	$Id: bestAlphaAndNu.cpp 1975 2007-04-22 13:47:28Z privmane $	
#include <iostream>
using namespace std;

#include "bestAlphaAndNu.h"

// ******************
// *     USSRV      *
// ******************

MDOUBLE bestFFixedTreeUSSRV::operator()(const tree& et, 
					   const sequenceContainer& sc,
					   const sequenceContainer& baseSc,
					   ussrvModel& model,
					   const Vdouble * weights,
					   const MDOUBLE upperBoundOnF,
					   const MDOUBLE epsilonFOptimization){
	
	MDOUBLE bestF=0;
	const MDOUBLE cx=upperBoundOnF;// left, middle, right limit on alpha
	const MDOUBLE bx=model.getF();
	const MDOUBLE ax=0.0;
	LOG(5,<<"****    Optimizing F    **** " << endl<< "bestFFixedTreeSSRV::operator() bx is :" << bx << endl);
	LOG(9,<<"ax is :" << ax << " cx is :" << cx << endl);
	_bestL = -brent(ax,bx,cx,
		C_evalFUSSRV(et,sc,baseSc,&model,weights),
		epsilonFOptimization,
		&bestF);
	setF(bestF,model);
	_bestF= bestF;
	return _bestL;
}

MDOUBLE bestAlphaFixedTreeUSSRV::operator()(const tree& et, //findBestAlphaFixedTree
					   const sequenceContainer& sc,
					   const sequenceContainer& baseSc,
					   ussrvModel& model,
					   const Vdouble * weights,
					   const MDOUBLE upperBoundOnAlpha,
					   const MDOUBLE epsilonAlphaOptimization){
	
	MDOUBLE bestA=0;
	const MDOUBLE cx=upperBoundOnAlpha;// left, middle, right limit on alpha
	const MDOUBLE bx=model.getAlpha();
	const MDOUBLE ax=0.0;
	LOG(5,<<"****    Optimizing Alpha    **** " << endl<< "bestAlphaFixedTreeSSRV::operator() bx is :" << bx << endl);
	_bestL = -brent(ax,bx,cx,
		C_evalAlphaUSSRV(et,sc,baseSc,&model,weights),
		epsilonAlphaOptimization,
		&bestA);
	setAlpha(bestA,model);
	_bestAlpha= bestA;
	return _bestL;
}

// Alpha is fixed
MDOUBLE bestNuFixedTreeUSSRV::operator()(const tree& et, 
					   const sequenceContainer& sc,
					   const sequenceContainer& baseSc,
					   ussrvModel& model,
					   const Vdouble * weights,
					   const MDOUBLE upperBoundOnNu,
					   const MDOUBLE epsilonNuOptimization){
		
	
	MDOUBLE bestN=0;
	// define the Nu bounds
	const MDOUBLE cx=upperBoundOnNu;// left, midle, right limit on alpha
	const MDOUBLE bx= model.getNu(); 
	const MDOUBLE ax=0.0;
	LOG(5,<<"****    Optimizing Nu    **** " << endl << "bestNuFixedTreeSSRV::operator() bx is : " << bx << endl);
	_bestL = -brent(ax,bx,cx, C_evalNuUSSRV(et,sc,baseSc,&model,weights), epsilonNuOptimization, &bestN);
	setNu(bestN,model);
	_bestNu= bestN;
	return _bestL;
}


// ******************
// *     SSRV       *
// ******************

MDOUBLE bestAlphaFixedTreeSSRV::operator()(const tree& et, //findBestAlphaFixedTree
	const sequenceContainer& sc, stochasticProcessSSRV& ssrvSp,	const Vdouble * weights,
	const MDOUBLE lowerBoundOnAlpha, const MDOUBLE upperBoundOnAlpha, const MDOUBLE epsilonAlphaOptimization){

	MDOUBLE bestA=0;
	const MDOUBLE cx=upperBoundOnAlpha;// left, midle, right limit on alpha
	replacementModelSSRV* pMulRM = static_cast<replacementModelSSRV*>(ssrvSp.getPijAccelerator()->getReplacementModel());
	gammaDistribution* gammaDist = static_cast<gammaDistribution*>(pMulRM->getDistribution()); 
	const MDOUBLE bx=gammaDist->getAlpha();
	const MDOUBLE ax=lowerBoundOnAlpha;
	LOG(5,<<"****    Optimizing Alpha    **** " << endl<< "bestAlphaFixedTreeSSRV::operator() bx is :" << bx << endl);
	_bestL = -brent(ax,bx,cx,
		C_evalAlphaSSRV(et,sc,ssrvSp,weights), epsilonAlphaOptimization, &bestA);
	
	setAlpha(bestA,ssrvSp);
	_bestAlpha= bestA;
	return _bestL;
}

// Alpha is fixed
MDOUBLE bestNuFixedTreeSSRV::operator()(const tree& et, const sequenceContainer& sc, 
	stochasticProcessSSRV& ssrvSp, const Vdouble * weights, const MDOUBLE lowerBoundOnNu, const MDOUBLE upperBoundOnNu,
	const MDOUBLE epsilonNuOptimization) {

	MDOUBLE bestN=0;
	// define the Nu bounds
	const MDOUBLE cx=upperBoundOnNu;// left, middle, right limit on alpha
	const MDOUBLE bx= static_cast<replacementModelSSRV*>(ssrvSp.getPijAccelerator()->getReplacementModel())->getRateOfRate();
	const MDOUBLE ax=lowerBoundOnNu;
	LOG(5,<<"****    Optimizing Nu    **** " << endl << "bestNuFixedTreeSSRV::operator() bx is : " << bx << endl);
	_bestL = -brent(ax,bx,cx, C_evalNuSSRV(et,sc,ssrvSp,weights), epsilonNuOptimization, &bestN);
	
	setNu(bestN,ssrvSp);
	_bestNu= bestN;
	return _bestL;
}


MDOUBLE bestTamura92ParamFixedTreeSSRV::operator()(const tree& et,
		const sequenceContainer& sc,
		stochasticProcessSSRV& ssrvSp,
		const Vdouble * weights/*= NULL */,
		const int maxTotalIterations /* = 5 */,
		const MDOUBLE epsilonLikelihoodImprovment /* = 0.05 */,
		const MDOUBLE lowerBoundOnTrTv /* = 0.0 */,
		const MDOUBLE upperBoundOnTrTv /* = 10.0 */,
		const MDOUBLE lowerBoundOnTheta /* = 0.0 */,
		const MDOUBLE upperBoundOnTheta /* = 1.0 */,
		const MDOUBLE epsilonTrTvOptimization /* = 0.01 */,
		const MDOUBLE epsilonThetaOptimization /* = 0.01 */){

	LOG(5,<<"Starting bestTamura92ParamFixedTreeSSRV::operator() :  find Best TrTv and theta"<<endl);
	MDOUBLE oldL = VERYSMALL;
	MDOUBLE newL = VERYSMALL;

	// first guess for the parameters
	MDOUBLE prevTrTv = static_cast<tamura92*>(static_cast<replacementModelSSRV*>(ssrvSp.getPijAccelerator()->getReplacementModel())->getBaseRM())->getTrTv();
	MDOUBLE prevTheta = static_cast<tamura92*>(static_cast<replacementModelSSRV*>(ssrvSp.getPijAccelerator()->getReplacementModel())->getBaseRM())->getTheta();

	for (int i=0; i < maxTotalIterations; ++i) {
		// optimize TrTv
		newL = -brent(lowerBoundOnTrTv, prevTrTv, upperBoundOnTrTv,
			C_evalTrTvSSRV(et,sc,ssrvSp,weights),
			epsilonTrTvOptimization,
			&_bestTrTv);
		setTrTv(_bestTrTv,ssrvSp);

		// optimize Theta
		newL = -brent(lowerBoundOnTheta, prevTheta, upperBoundOnTheta,
			C_evalThetaSSRV(et,sc,ssrvSp,weights),
			epsilonThetaOptimization,
			&_bestTheta);
		setTheta(_bestTheta,ssrvSp);

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
				LOG(5,<<"bestTamura92ParamFixedTreeSSRV::operator() likelihood went down!"<<endl<<"oldL = "<< oldL <<" newL= "<<newL<<endl);
				_bestL = oldL;
				_bestTrTv = prevTrTv;
				_bestTheta = prevTheta;
				setTrTvAndTheta(prevTrTv,prevTheta,ssrvSp);
			}
			break;
		}
	}
	return _bestL;
}
