// $Id: ssrvDistanceSeqs2Tree.cpp 962 2006-11-07 15:13:34Z privmane $

#include "ssrvDistanceSeqs2Tree.h"
//#include "bestAlphaAndNu.h"
#include "bestParamUSSRV.h"
#include "someUtil.h"
#include <float.h>

tree ssrvDistanceSeqs2Tree::seqs2TreeIterative(const sequenceContainer &sc, MDOUBLE initAlpha, MDOUBLE initNu, const Vdouble *weights, const tree* constraintTreePtr) {
	_constraintTreePtr=constraintTreePtr;
	_alpha = initAlpha;
	_newNu = _nu = initNu;
	_weights = weights;
	return seqs2TreeIterativeInternal(sc, true);
}

tree ssrvDistanceSeqs2Tree::seqs2TreeIterative(const sequenceContainer &sc, const Vdouble *weights, const tree* constraintTreePtr) {
	_constraintTreePtr=constraintTreePtr;
	_weights = weights;
	return seqs2TreeIterativeInternal(sc, false);
}

tree ssrvDistanceSeqs2Tree::seqs2TreeIterative(const sequenceContainer &sc, const tree &initTree, const Vdouble *weights, const tree* constraintTreePtr) {
	_constraintTreePtr=constraintTreePtr;
	_weights = weights;
	return seqs2TreeIterativeInternalInitTreeGiven(sc, initTree);
}

tree ssrvDistanceSeqs2Tree::seqs2TreeIterative(const sequenceContainer &sc, const tree &initTree, MDOUBLE initAlpha, const Vdouble *weights, const tree* constraintTreePtr) {
	_alpha = initAlpha;
	_weights = weights;

	_constraintTreePtr=constraintTreePtr;
	return seqs2TreeIterativeInternalInitTreeGiven(sc, false, initTree, initAlpha);
}

tree ssrvDistanceSeqs2Tree::seqs2TreeIterative(const sequenceContainer &sc, const tree &initTree, MDOUBLE initAlpha, MDOUBLE initNu, const Vdouble *weights, const tree* constraintTreePtr) {
	_alpha = initAlpha;
	_newNu = _nu = initNu;
	_weights = weights;

	_constraintTreePtr=constraintTreePtr;
	return seqs2TreeIterativeInternalInitTreeGiven(sc, true, initTree, initAlpha);
}

// NOTE! This version is a NON-ITERATIVE version that uses the side info supplied by the user
tree ssrvDistanceSeqs2Tree::seqs2Tree(const sequenceContainer &sc, MDOUBLE alpha, MDOUBLE nu, const Vdouble *weights, const tree* constraintTreePtr) {
	_weights = weights;
	_alpha = alpha;
	_newNu = _nu = nu;
	_constraintTreePtr=constraintTreePtr;
	seqs2TreeOneIterationInternal(sc, true);
	return _newTree;
}

tree ssrvDistanceSeqs2Tree::seqs2TreeBootstrap(const sequenceContainer &sc, const MDOUBLE alpha, MDOUBLE nu, const Vdouble *weights, const tree* constraintTreePtr) {
	_weights = weights;
	_alpha = alpha;
	_newNu = _nu = nu;
	return static_cast<iterativeDistanceSeqs2Tree *>(this)->seqs2TreeBootstrap(sc, weights, constraintTreePtr);
}

// NOTE! This version calls ITERATIVE seqs2Tree because side info is not given by the user, so we have to generate and optimize it
tree ssrvDistanceSeqs2Tree::seqs2Tree(const sequenceContainer &sc, const Vdouble *weights, const tree* constraintTreePtr) {
	return seqs2TreeIterative(sc,weights,constraintTreePtr);
}

MDOUBLE ssrvDistanceSeqs2Tree::optimizeSideInfo(const sequenceContainer &sc, tree &et)
{
	if (!dynamic_cast<tamura92*>(
			static_cast<replacementModelSSRV*>(_spPtr->getPijAccelerator()->getReplacementModel())
			->getBaseRM()
			)
		) {
		bestParamSSRV optimizer(true,true,false,true); // optimize alpha, nu, NOT tamura92 params, and bbl		
		optimizer(et,sc,*static_cast<stochasticProcessSSRV*>(_spPtr),_weights,
				  15,15,0.5,_epsilonLikelihoodImprovement4alphaOptimiz,_epsilonLikelihoodImprovement,
				  _epsilonLikelihoodImprovement4BBL,_maxIterationsBBL,5);
		_newAlpha=optimizer.getBestAlpha();
		_newNu=optimizer.getBestNu();
		return(optimizer.getBestL());
	} else {
		bestParamSSRV optimizer(true,true,true,true); // optimize alpha, nu, tamura92 params, and bbl		
		optimizer(et,sc,*static_cast<stochasticProcessSSRV*>(_spPtr),_weights,
				  15,15,0.5,_epsilonLikelihoodImprovement4alphaOptimiz,_epsilonLikelihoodImprovement,
				  _epsilonLikelihoodImprovement4BBL,_maxIterationsBBL,5);
		_newAlpha=optimizer.getBestAlpha();
		_newNu=optimizer.getBestNu();
		return(optimizer.getBestL());
	}
}

MDOUBLE ssrvDistanceSeqs2Tree::calcSideInfoGivenTreeAndAlpha(const sequenceContainer &sc, const tree &et, MDOUBLE alpha) 
{
	_newAlpha = alpha;
	(static_cast<gammaDistribution*>(_spPtr->distr()))->setAlpha(alpha);

	// optimize only nu (and tamura92 params, if relevant)
	if (!dynamic_cast<tamura92*>(
			static_cast<replacementModelSSRV*>(_spPtr->getPijAccelerator()->getReplacementModel())
			->getBaseRM()
			)
		) {
		bestParamSSRV optimizer(false,true,false,false);
		optimizer(et,sc,*(static_cast<stochasticProcessSSRV*>(_spPtr)),_weights,
				  15,15,_epsilonLikelihoodImprovement4alphaOptimiz,_epsilonLikelihoodImprovement,
				  _epsilonLikelihoodImprovement4BBL,_maxIterationsBBL,5);
		_newNu=optimizer.getBestNu();
		return(optimizer.getBestL());
	} else {
		bestParamSSRV optimizer(false,true,true,false);
		optimizer(et,sc,*(static_cast<stochasticProcessSSRV*>(_spPtr)),_weights,
				  15,15,_epsilonLikelihoodImprovement4alphaOptimiz,_epsilonLikelihoodImprovement,
				  _epsilonLikelihoodImprovement4BBL,_maxIterationsBBL,5);
		_newNu=optimizer.getBestNu();
		return(optimizer.getBestL());
	}
}

void ssrvDistanceSeqs2Tree::acceptSideInfo()
{
	_alpha = _newAlpha;
	_nu = _newNu;
}

void ssrvDistanceSeqs2Tree::utilizeSideInfo()
{
	// set new alpha value in the sp that is used in _distM
	LOG(10,<<"# utilizing alpha "<<_alpha<<" and nu "<<_nu<<endl);
	(static_cast<gammaDistribution*>(_spPtr->distr()))->setAlpha(_alpha);
	(static_cast<stochasticProcessSSRV*>(_spPtr))->setRateOfRate(_nu);
}

void ssrvDistanceSeqs2Tree::printSideInfo(ostream& out) const
{
	out<<"Alpha: "<< _alpha <<" Nu: "<< _nu <<endl;
}

// non virtual
void ssrvDistanceSeqs2Tree::setSideInfo(const MDOUBLE alpha, MDOUBLE nu)
{
	_alpha = alpha;
	_nu = nu;
}

ssrvDistanceSeqs2Tree::alphaAndNu ssrvDistanceSeqs2Tree::getSideInfo() const
{
	return alphaAndNu(_alpha, _nu);
}
