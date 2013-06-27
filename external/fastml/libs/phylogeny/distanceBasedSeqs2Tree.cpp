// $Id: distanceBasedSeqs2Tree.cpp 6002 2009-03-20 19:39:03Z privmane $

#include "distanceBasedSeqs2Tree.h"
#include "uniDistribution.h"
#include "distanceTable.h"
#include "bestAlpha.h"
#include "siteSpecificRate.h"
#include "someUtil.h"
#include "bblEM.h"
#include "tamura92.h"
#include "bestTamura92param.h"
#include "bestGtrModelParams.h"
#include <float.h>
#include "replacementModelSSRV.h"
#include "trivialAccelerator.h"

// **********************************************************************
// *** The basic non-iterative versions *********************************
// **********************************************************************

tree distanceBasedSeqs2Tree::seqs2Tree(const sequenceContainer &sc, const Vdouble *weights, const tree* constraintTreePtr) {
	_constraintTreePtr=constraintTreePtr;
	_weights = weights;

	// Calculate distance table
	tree et;
	VVdouble distTable;
	vector<string> vNames;
	giveDistanceTable(_distM,sc,distTable,vNames,_weights);

	// Build tree from the distance table
	et = _dist2et->computeTree(distTable, vNames, _constraintTreePtr);

	LOG(6,<<"# distanceBasedSeqs2Tree::seqs2Tree: The reconsructed tree:"<<endl);
	LOGDO(6,et.output(myLog::LogFile()));

	return et;
}

tree distanceBasedSeqs2Tree::seqs2TreeBootstrap(const sequenceContainer &sc, const Vdouble *weights, const tree* constraintTreePtr) {
	return seqs2Tree(sc, weights, constraintTreePtr);
}

// **********************************************************************
// *** iterativeDistanceSeqs2Tree ***************************************
// **********************************************************************

iterativeDistanceSeqs2Tree::iterativeDistanceSeqs2Tree(likeDist &distM, distances2Tree &dist2et, const Vdouble *weights,
													   const MDOUBLE epsilonLikelihoodImprovement, 
													   const MDOUBLE epsilonLikelihoodImprovement4alphaOptimiz, 
													   const MDOUBLE epsilonLikelihoodImprovement4BBL, 
													   const int maxIterationsBBL)
	: distanceBasedSeqs2Tree(distM, dist2et, weights),
	  _epsilonLikelihoodImprovement             ( epsilonLikelihoodImprovement             ),
	  _epsilonLikelihoodImprovement4alphaOptimiz( epsilonLikelihoodImprovement4alphaOptimiz),
	  _epsilonLikelihoodImprovement4BBL         ( epsilonLikelihoodImprovement4BBL         ),
	  _maxIterationsBBL                         ( maxIterationsBBL                         )
{
	// Check that the stochasticProcess in likeDist is not const
	if (distM.isTheInternalStochasticProcessConst()) {
		errorMsg::reportError("iterativeDistanceSeqs2Tree::iterativeDistanceSeqs2Tree: The stochasticProcess in the given likeDist object is const. A non-const stochasticProcess is required.");
	}

	// Keep a pointer to the stochasticProcess in distM, so that we will be able to change its alpha, etc.
	_spPtr = &(distM.getNonConstStochasticProcess());
	if (_spPtr->categories() >1)
		_alpha = (static_cast<gammaDistribution*>(_spPtr->distr()))->getAlpha();
	else
		_alpha=-99.9;				// this should never be used

}

// *** Iterative tree building ******************************************
tree iterativeDistanceSeqs2Tree::seqs2TreeIterativeInternal(const sequenceContainer &sc, bool initSideInfoGiven) {
	LOGDO(3,printTime(myLog::LogFile()));
	LOG(3,<<"# iterativeDistanceSeqs2Tree::seqs2TreeIterativeInternal:"<<endl<<"# Initial tree:"<<endl);
	seqs2TreeOneIterationInternal(sc, initSideInfoGiven);

	return seqs2TreeIterativeInternalInitTreeGiven(sc, true, _newTree, _newAlpha);
}

// *** Iterative tree building, given an initial tree and alpha *********
// *** Optimize branch lengths and sideInfo for the given tree topology
tree iterativeDistanceSeqs2Tree::seqs2TreeIterativeInternalInitTreeGiven(const sequenceContainer &sc, const tree &initTree) {
	LOG(7,<<"# iterativeDistanceSeqs2Tree::seqs2TreeIterativeInternalInitTreeGiven: Started optimizeSideInfo. ");
	LOGDO(7,printTime(myLog::LogFile()));
	_newTree=initTree;
	_newTreeLogLikelihood=optimizeSideInfo(sc, _newTree);
	LOG(7,<<"# iterativeDistanceSeqs2Tree::seqs2TreeIterativeInternalInitTreeGiven: Finished optimizeSideInfo. ");
	LOGDO(7,printTime(myLog::LogFile()));

	return seqs2TreeIterativeInternalInitTreeGiven(sc, true, _newTree, _newAlpha);
}

// *** Iterative tree building, given an initial tree and alpha *********
// *** If sideInfo is not given - calculate it for the fixed tree and alpha
tree iterativeDistanceSeqs2Tree::seqs2TreeIterativeInternalInitTreeGiven(const sequenceContainer &sc, bool initSideInfoGiven, const tree &initTree, MDOUBLE initAlpha) {
	_newTree=initTree;
	_newAlpha=initAlpha;

	LOGDO(3,printTime(myLog::LogFile()));
	LOG(3,<<"# iterativeDistanceSeqs2Tree::seqs2TreeIterativeInternalInitTreeGiven"<<endl);
	if (!initSideInfoGiven) {
		_newTreeLogLikelihood=calcSideInfoGivenTreeAndAlpha(sc, initTree, initAlpha);
	}
	int iterationNum = 0;
	LOGDO(3,printTime(myLog::LogFile()));
	LOG(3,<<"# iterativeDistanceSeqs2Tree::seqs2TreeIterativeInternalInitTreeGiven:"<<endl<<"# The given initial tree:"<<endl);
	LOGDO(3,_newTree.output(myLog::LogFile()));
  
	do {
		++iterationNum;
		LOGDO(5,printTime(myLog::LogFile()));
		LOG(3,<<"# Iteration "<<iterationNum<<":"<<endl);

		// save the best tree so far, and its likelihood and the sideInfo that was calculated for it
		_et=_newTree;				
		_treeLogLikelihood=_newTreeLogLikelihood;
		acceptSideInfo();
		LOG(7,<<"# Side info for the tree"<<endl);
		LOGDO(7,printSideInfo(myLog::LogFile()));

		seqs2TreeOneIterationInternal(sc, true);

	} while (_newTreeLogLikelihood > _treeLogLikelihood + _epsilonLikelihoodImprovement);

	LOGDO(3,printTime(myLog::LogFile()));
	LOG(3,<<"# iterativeDistanceSeqs2Tree::seqs2TreeIterativeInternalInitTreeGiven:"<<endl<<"# Finished iterative distance-based tree reconstruction, done "<<iterationNum<<" iterations"<<endl);
	return _et;
}

// *** Tree building procedure that is called iteratively **********************
void iterativeDistanceSeqs2Tree::seqs2TreeOneIterationInternal(const sequenceContainer &sc, const bool sideInfoSet) {

	// 1. Calculate distance table
	VVdouble distTable;
	vector<string> vNames;
	LOG(7,<<"# iterativeDistanceSeqs2Tree::seqs2TreeOneIterationInternal: Started giveDistanceTable. ");
	LOGDO(7,printTime(myLog::LogFile()));
	if (!sideInfoSet) { // Then use homogeneous rates

		// Create homogeneous likeDist
		_alpha = 1.5;		// Since no ASRV side info is known yet, we set an initial alpha for bestAlphaAndBBL optimizations
		uniDistribution distribution;
		stochasticProcess* uniDistSp = NULL;
		replacementModelSSRV* rmSSRV = 
			dynamic_cast<replacementModelSSRV*>(_spPtr->getPijAccelerator()->getReplacementModel());
		if (!rmSSRV) {
			uniDistSp = new stochasticProcess(&distribution, _spPtr->getPijAccelerator());
		} else {
			trivialAccelerator pijAcc(rmSSRV->getBaseRM());
			uniDistSp = new stochasticProcess(&distribution, &pijAcc);
		}
		likeDist homogeneousDist(*uniDistSp,static_cast<likeDist*>(_distM)->getToll());

		giveDistanceTable(&homogeneousDist,sc,distTable,vNames,_weights);
		delete uniDistSp;

	} else {			// use the side information
		utilizeSideInfo();
		giveDistanceTable(_distM,sc,distTable,vNames,_weights);
	}
	LOG(7,<<"# iterativeDistanceSeqs2Tree::seqs2TreeOneIterationInternal: Finished giveDistanceTable, started distances2Tree::computeTree. ");
	LOGDO(7,printTime(myLog::LogFile()));

	// 2. Build tree from the distance table
	_newTree = _dist2et->computeTree(distTable, vNames, _constraintTreePtr);
	LOG(7,<<"# iterativeDistanceSeqs2Tree::seqs2TreeOneIterationInternal: Finished distances2Tree::computeTree, started optimizeSideInfo. ");
	LOGDO(7,printTime(myLog::LogFile()));

	// 3. Optimize branch lengths and side info for the tree topology
	_newTreeLogLikelihood=optimizeSideInfo(sc, _newTree);
	LOG(7,<<"# iterativeDistanceSeqs2Tree::seqs2TreeOneIterationInternal: Finished distances2Tree::optimizeSideInfo. ");
	LOGDO(7,printTime(myLog::LogFile()));

	if (!sideInfoSet) {
		LOG(5,<<"# iterativeDistanceSeqs2Tree::seqs2TreeOneIterationInternal:"<<endl<<"# Homogeneous rates tree"<<endl);
	} else {
		LOG(5,<<"# iterativeDistanceSeqs2Tree::seqs2TreeOneIterationInternal:"<<endl<<"# Tree based on alpha"<<endl);
	}
	LOGDO(5,_newTree.output(myLog::LogFile()));
	LOG(5,<<"# Log likelihood:"<<endl<<_newTreeLogLikelihood<<endl);
}

// Perform one bootstrap iteration, assuming that side info has been set (as if acceptSideInfo has been called)
tree iterativeDistanceSeqs2Tree::seqs2TreeBootstrap(const sequenceContainer &sc, const Vdouble *weights, const tree* constraintTreePtr) {
	LOG(3,<<"# iterativeDistanceSeqs2Tree::seqs2TreeBootstrap: Started a single bootstrap iteration. ");
	LOGDO(3,printTime(myLog::LogFile()));
	_constraintTreePtr=constraintTreePtr;
	_weights = weights;

	// Calculate distance table
	tree localScopeEt;
	VVdouble distTable;
	vector<string> vNames;
	utilizeSideInfo();
	giveDistanceTable(_distM,sc,distTable,vNames,_weights);

	// Build tree from the distance table
	localScopeEt = _dist2et->computeTree(distTable,vNames, _constraintTreePtr);

	LOG(3,<<"# iterativeDistanceSeqs2Tree::seqs2TreeBootstrapInternal:"<<endl<<"# Bootstrap tree based on alpha, without optimizations"<<endl);
	LOGDO(3,localScopeEt.output(myLog::LogFile()));

	return localScopeEt;
}

/********************************
 * commonAlphaDistanceSeqs2Tree *
 ********************************/
tree commonAlphaDistanceSeqs2Tree::seqs2TreeIterative(const sequenceContainer &sc, MDOUBLE initAlpha, const Vdouble *weights, const tree* constraintTreePtr) {
	_constraintTreePtr=constraintTreePtr;
	_alpha = initAlpha;
	_weights = weights;
	return seqs2TreeIterativeInternal(sc, true);
}

tree commonAlphaDistanceSeqs2Tree::seqs2TreeIterative(const sequenceContainer &sc, const Vdouble *weights, const tree* constraintTreePtr) {
	_constraintTreePtr=constraintTreePtr;
	_weights = weights;
	return seqs2TreeIterativeInternal(sc, false);
}

tree commonAlphaDistanceSeqs2Tree::seqs2TreeIterative(const sequenceContainer &sc, const tree &initTree, const Vdouble *weights, const tree* constraintTreePtr) {
	_constraintTreePtr=constraintTreePtr;
	_weights = weights;
	return seqs2TreeIterativeInternalInitTreeGiven(sc, initTree);
}

tree commonAlphaDistanceSeqs2Tree::seqs2TreeIterative(const sequenceContainer &sc, const tree &initTree, MDOUBLE initAlpha, const Vdouble *weights, const tree* constraintTreePtr) {
	_alpha = initAlpha;
	_weights = weights;

	_constraintTreePtr=constraintTreePtr;
	return seqs2TreeIterativeInternalInitTreeGiven(sc, true, initTree, initAlpha);
}

// NOTE! This version is a NON-ITERATIVE version that uses the side info supplied by the user
tree commonAlphaDistanceSeqs2Tree::seqs2Tree(const sequenceContainer &sc, MDOUBLE alpha, const Vdouble *weights, const tree* constraintTreePtr) {
	_weights = weights;
	_alpha = alpha;
	_constraintTreePtr=constraintTreePtr;
	seqs2TreeOneIterationInternal(sc, true);
	return _newTree;
}

tree commonAlphaDistanceSeqs2Tree::seqs2TreeBootstrap(const sequenceContainer &sc, const MDOUBLE alpha, const Vdouble *weights, const tree* constraintTreePtr) {
	_weights = weights;
	_alpha = alpha;
	return static_cast<iterativeDistanceSeqs2Tree *>(this)->seqs2TreeBootstrap(sc, weights, constraintTreePtr);
}

// NOTE! This version calls ITERATIVE seqs2Tree because side info is not given by the user, so we have to generate and optimize it
tree commonAlphaDistanceSeqs2Tree::seqs2Tree(const sequenceContainer &sc, const Vdouble *weights, const tree* constraintTreePtr) {
	return seqs2TreeIterative(sc,weights,constraintTreePtr);
}

MDOUBLE commonAlphaDistanceSeqs2Tree::optimizeSideInfo(const sequenceContainer &sc, tree &et)
{
	if (dynamic_cast<tamura92*>(_spPtr->getPijAccelerator()->getReplacementModel())) {
		// Optimizing params of the tamura92 model
		bestTamura92ParamAlphaAndBBL optimizer(et, sc, *_spPtr, _weights, 5, _epsilonLikelihoodImprovement/*0.05*/,
											   _epsilonLikelihoodImprovement4alphaOptimiz/*0.01*/, 
											   _epsilonLikelihoodImprovement4alphaOptimiz/*0.01*/, 
											   _epsilonLikelihoodImprovement4alphaOptimiz/*0.01*/, 
											   _epsilonLikelihoodImprovement4BBL/*0.01*/,
											   5.0, _maxIterationsBBL, _alpha, 5.0 );
		_newAlpha=optimizer.getBestAlpha();
		return(optimizer.getBestL());

	} else if (dynamic_cast<gtrModel*>(_spPtr->getPijAccelerator()->getReplacementModel())) {
		// Optimizing params of the gtr model
		bestGtrModel optimizer(et, sc, *_spPtr, _weights, 5,
							   _epsilonLikelihoodImprovement,
							   _epsilonLikelihoodImprovement4alphaOptimiz,
							   true, true);
		_newAlpha=optimizer.getBestAlpha();
		return(optimizer.getBestL());

	} else {
		bestAlphaAndBBL optimizer(et, sc, *_spPtr, _weights, _alpha, 5.0,
								  _epsilonLikelihoodImprovement4BBL/*0.01*/, _epsilonLikelihoodImprovement4alphaOptimiz,
								  _maxIterationsBBL);
		_newAlpha=optimizer.getBestAlpha();
		return(optimizer.getBestL());
	}
}

MDOUBLE commonAlphaDistanceSeqs2Tree::calcSideInfoGivenTreeAndAlpha(const sequenceContainer &sc, const tree &et, MDOUBLE alpha) 
{
	_newAlpha = alpha;
	(static_cast<gammaDistribution*>(_spPtr->distr()))->setAlpha(alpha);
	return likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(et, sc, *_spPtr, _weights);
}

void commonAlphaDistanceSeqs2Tree::acceptSideInfo()
{
	_alpha = _newAlpha;
}

void commonAlphaDistanceSeqs2Tree::utilizeSideInfo()
{
	// set new alpha value in the sp that is used in _distM
	(static_cast<gammaDistribution*>(_spPtr->distr()))->setAlpha(_alpha);
	LOG(10,<<"# utilizing alpha"<<endl<<_alpha<<endl<<endl);

}

void commonAlphaDistanceSeqs2Tree::printSideInfo(ostream& out) const
{
	out<<"Alpha: "<<_alpha<<endl;
}

// non virtual
void commonAlphaDistanceSeqs2Tree::setSideInfo(const MDOUBLE alpha)
{
	_alpha=alpha;
}

MDOUBLE commonAlphaDistanceSeqs2Tree::getSideInfo() const
{
	return _alpha;
}

/******************************
 * rate4siteDistanceSeqs2Tree *
 ******************************/
tree rate4siteDistanceSeqs2Tree::seqs2TreeIterative(const sequenceContainer &sc, const Vdouble &initRates, const Vdouble *weights, const tree* constraintTreePtr) {
	_rates = initRates;
	_constraintTreePtr=constraintTreePtr;
	_weights = weights;
	return seqs2TreeIterativeInternal(sc, true);
}

tree rate4siteDistanceSeqs2Tree::seqs2TreeIterative(const sequenceContainer &sc, const Vdouble *weights, const tree* constraintTreePtr) {
	_constraintTreePtr=constraintTreePtr;
	_weights = weights;
	return seqs2TreeIterativeInternal(sc, false);
}

tree rate4siteDistanceSeqs2Tree::seqs2TreeIterative(const sequenceContainer &sc, const tree &initTree, const Vdouble *weights, const tree* constraintTreePtr) {
	_constraintTreePtr=constraintTreePtr;
	_weights = weights;
	return seqs2TreeIterativeInternalInitTreeGiven(sc, initTree);
}

tree rate4siteDistanceSeqs2Tree::seqs2TreeIterative(const sequenceContainer &sc, const tree &initTree, MDOUBLE initAlpha, const Vdouble *weights, const tree* constraintTreePtr) {
	_constraintTreePtr=constraintTreePtr;
	_weights = weights;
	return seqs2TreeIterativeInternalInitTreeGiven(sc, false, initTree, initAlpha);
}

// NOTE! This version is a NON-ITERATIVE version that uses the side info supplied by the user
tree rate4siteDistanceSeqs2Tree::seqs2Tree(const sequenceContainer &sc, const Vdouble &rates, const Vdouble *weights, const tree* constraintTreePtr) {
	_weights = weights;
	_rates = rates;
	_constraintTreePtr=constraintTreePtr;

	seqs2TreeOneIterationInternal(sc, true);
	return _newTree;
}

tree rate4siteDistanceSeqs2Tree::seqs2TreeBootstrap(const sequenceContainer &sc, const Vdouble &rates, const Vdouble *weights, const tree* constraintTreePtr) {
	_weights = weights;
	_rates = rates;
	return static_cast<iterativeDistanceSeqs2Tree *>(this)->seqs2TreeBootstrap(sc, weights, constraintTreePtr);
}

// NOTE! This version calls ITERATIVE seqs2Tree because side info is not given by the user, so we have to generate and optimize it
tree rate4siteDistanceSeqs2Tree::seqs2Tree(const sequenceContainer &sc, const Vdouble *weights, const tree* constraintTreePtr) {
	return seqs2TreeIterative(sc,weights,constraintTreePtr);
}

MDOUBLE rate4siteDistanceSeqs2Tree::optimizeSideInfo(const sequenceContainer &sc, tree &et)
{
	bblEM optimizer(et, sc, *_spPtr, _weights, _maxIterationsBBL, _epsilonLikelihoodImprovement4BBL);

	// Note: this verstion of ML rates computation can only use a uniDistribution stochasticProcess
	Vdouble likelihoods;
	MDOUBLE treeLogLikelihood = computeML_siteSpecificRate(_newRates, likelihoods, sc, *_spPtr, et,20,_epsilonLikelihoodImprovement);
	//computeEB_EXP_siteSpecificRate
	return(treeLogLikelihood);
}

MDOUBLE rate4siteDistanceSeqs2Tree::calcSideInfoGivenTreeAndAlpha(const sequenceContainer &sc, const tree &et, MDOUBLE alpha) 
{
	_newAlpha = alpha;
	Vdouble likelihoods;
	MDOUBLE treeLogLikelihood = computeML_siteSpecificRate(_newRates, likelihoods, sc, *_spPtr, et,20,_epsilonLikelihoodImprovement);
	//computeEB_EXP_siteSpecificRate
	return(treeLogLikelihood);
}

void rate4siteDistanceSeqs2Tree::acceptSideInfo()
{
	_alpha = _newAlpha;
	_rates = _newRates;
}

void rate4siteDistanceSeqs2Tree::utilizeSideInfo()
{
	(static_cast<givenRatesMLDistance*>(_distM))->setRates(_rates);
	LOG(10,<<"# utilizing rates"<<endl<<_rates<<endl<<endl);

	// set new alpha value in the sp that is used in _distM 
	//  (static_cast<gammaDistribution*>(_spPtr->distr()))->setAlpha(_alpha);
}

void rate4siteDistanceSeqs2Tree::printSideInfo(ostream& out) const
{
	if (_rates.size())
		out<<"ML rates: "<<_rates<<endl;
}

// non virtual
void rate4siteDistanceSeqs2Tree::setSideInfo(const Vdouble &rates)
{
	_rates = rates;
}

const Vdouble& rate4siteDistanceSeqs2Tree::getSideInfo() const
{
	return _rates;
}

/******************************
 * posteriorDistanceSeqs2Tree *
 ********************************/
tree posteriorDistanceSeqs2Tree::seqs2TreeIterative(const sequenceContainer &sc, MDOUBLE initAlpha, const VVdoubleRep &initPosterior, const Vdouble *weights, const tree* constraintTreePtr) {
	_alpha = initAlpha;
	_posterior = initPosterior;
	_weights = weights;
	_constraintTreePtr=constraintTreePtr;
	return seqs2TreeIterativeInternal(sc, true);
}

tree posteriorDistanceSeqs2Tree::seqs2TreeIterative(const sequenceContainer &sc, const Vdouble *weights, const tree* constraintTreePtr) {
	_constraintTreePtr=constraintTreePtr;
	_weights = weights;
	return seqs2TreeIterativeInternal(sc, false);
}

tree posteriorDistanceSeqs2Tree::seqs2TreeIterative(const sequenceContainer &sc, const tree &initTree, const Vdouble *weights, const tree* constraintTreePtr) {
	_constraintTreePtr=constraintTreePtr;
	_weights = weights;
	return seqs2TreeIterativeInternalInitTreeGiven(sc, initTree);
}

tree posteriorDistanceSeqs2Tree::seqs2TreeIterative(const sequenceContainer &sc, const tree &initTree, MDOUBLE initAlpha, const Vdouble *weights, const tree* constraintTreePtr) {
	_constraintTreePtr=constraintTreePtr;
	_weights = weights;
	return seqs2TreeIterativeInternalInitTreeGiven(sc, false, initTree, initAlpha);
}

tree posteriorDistanceSeqs2Tree::seqs2TreeIterative(const sequenceContainer &sc, const tree &initTree, MDOUBLE initAlpha, const VVdoubleRep &initPosterior, const Vdouble *weights, const tree* constraintTreePtr) {
	_alpha = initAlpha;
	_posterior = initPosterior;
	_weights = weights;
	_constraintTreePtr=constraintTreePtr;
	return seqs2TreeIterativeInternalInitTreeGiven(sc, true, initTree, initAlpha);
}

// NOTE! This version is a NON-ITERATIVE version that uses the side info supplied by the user
tree posteriorDistanceSeqs2Tree::seqs2Tree(const sequenceContainer &sc, const VVdoubleRep &posterior, const Vdouble *weights, const tree* constraintTreePtr) {
	_weights = weights;
	_posterior = posterior;
	_constraintTreePtr=constraintTreePtr;
	seqs2TreeOneIterationInternal(sc, true);
	return _newTree;
}

tree posteriorDistanceSeqs2Tree::seqs2TreeBootstrap(const sequenceContainer &sc, const VVdoubleRep &posterior, const Vdouble *weights, const tree* constraintTreePtr) {
	_weights = weights;
	_posterior = posterior;
	return static_cast<iterativeDistanceSeqs2Tree *>(this)->seqs2TreeBootstrap(sc, weights, constraintTreePtr);
}

// NOTE! This version calls ITERATIVE seqs2Tree because side info is not given by the user, so we have to generate and optimize it
tree posteriorDistanceSeqs2Tree::seqs2Tree(const sequenceContainer &sc, const Vdouble *weights, const tree* constraintTreePtr) {
  return seqs2TreeIterative(sc, weights, constraintTreePtr);
}

MDOUBLE posteriorDistanceSeqs2Tree::optimizeSideInfo(const sequenceContainer &sc, tree &et)
{
	if (dynamic_cast<tamura92*>(_spPtr->getPijAccelerator()->getReplacementModel())) {
		// Optimizing params of the tamura92 model
		bestTamura92ParamAlphaAndBBL optimizer(et, sc, *_spPtr, _weights, 5, _epsilonLikelihoodImprovement/*0.05*/,
											   _epsilonLikelihoodImprovement4alphaOptimiz/*0.01*/, 
											   _epsilonLikelihoodImprovement4alphaOptimiz/*0.01*/, 
											   _epsilonLikelihoodImprovement4alphaOptimiz/*0.01*/, 
											   _epsilonLikelihoodImprovement4BBL/*0.01*/,
											   5.0, _maxIterationsBBL, _alpha, 5.0 );
		_newAlpha=optimizer.getBestAlpha();
		return(optimizer.getBestL());

	} else if (dynamic_cast<gtrModel*>(_spPtr->getPijAccelerator()->getReplacementModel())) {
		// Optimizing params of the gtr model
		bestGtrModel optimizer(et, sc, *_spPtr, _weights, 5,
							   _epsilonLikelihoodImprovement,
							   _epsilonLikelihoodImprovement4alphaOptimiz,
							   true, true);
		_newAlpha=optimizer.getBestAlpha();
		return(optimizer.getBestL());

	} else {
		bestAlphaAndBBL optimizer(et, sc, *_spPtr, _weights, _alpha, 5.0,
								  _epsilonLikelihoodImprovement4BBL/*0.01*/, _epsilonLikelihoodImprovement4alphaOptimiz,
								  _maxIterationsBBL);
		_newAlpha=optimizer.getBestAlpha();	// cached only to make alpha optimization faster
	}

	// Compute posterior probabilities of rates per site
	return likelihoodComputation::getPosteriorOfRates(et, sc, *_spPtr, _newPosterior);
}

MDOUBLE posteriorDistanceSeqs2Tree::calcSideInfoGivenTreeAndAlpha(const sequenceContainer &sc, const tree &et, MDOUBLE alpha) 
{
	_newAlpha = alpha;
	(static_cast<gammaDistribution*>(_spPtr->distr()))->setAlpha(alpha);
	// Compute posterior probabilities of rates per site
	return likelihoodComputation::getPosteriorOfRates(et, sc, *_spPtr, _newPosterior);
}

void posteriorDistanceSeqs2Tree::acceptSideInfo()
{
	_alpha = _newAlpha;
	_posterior = _newPosterior;
}

void posteriorDistanceSeqs2Tree::utilizeSideInfo()
{
	(static_cast<posteriorDistance*>(_distM))->setPosterior(_posterior);
	LOG(10,<<"# utilizing posterior"<<endl<<_posterior<<endl<<endl);
	// set new alpha value in the sp that is used in _distM 
	//  (static_cast<gammaDistribution*>(_spPtr->distr()))->setAlpha(_alpha);
}

void posteriorDistanceSeqs2Tree::printSideInfo(ostream& out) const
{
	if (_posterior.size())
		out<<_posterior<<endl;
}

// non virtual
void posteriorDistanceSeqs2Tree::setSideInfo(const VVdoubleRep &posterior)
{
	_posterior = posterior;
}

const VVdoubleRep& posteriorDistanceSeqs2Tree::getSideInfo() const
{
	return _posterior;
}

