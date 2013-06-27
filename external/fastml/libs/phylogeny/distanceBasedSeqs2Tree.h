// $Id: distanceBasedSeqs2Tree.h 5989 2009-03-19 09:27:26Z privmane $

#ifndef ___DISTANCE_BASED_SEQS2TREE
#define ___DISTANCE_BASED_SEQS2TREE

#include "distanceMethod.h"
#include "sequenceContainer.h"
#include "stochasticProcess.h"
#include "likeDist.h"
#include "distances2Tree.h"
#include "givenRatesMLDistance.h"
#include "posteriorDistance.h"
#include "float.h"

// NOTE: These modules take sequenceContainer as argument, and do not
// manipulate it.  If you want to take care of gaps do it yourself!
class distanceBasedSeqs2Tree {
public:
    distanceBasedSeqs2Tree(distanceMethod &distM, distances2Tree &dist2et, const Vdouble *weights = NULL)
	  : _distM(distM.clone()), _dist2et(dist2et.clone()), _weights(weights), _treeLogLikelihood(VERYBIG) {}
  virtual ~distanceBasedSeqs2Tree() {delete (_distM);delete (_dist2et);}
    virtual tree seqs2Tree(const sequenceContainer &sc, const Vdouble *weights=NULL, const tree* constraintTreePtr=NULL);
    // Does one bootstrap iteration
    virtual tree seqs2TreeBootstrap(const sequenceContainer &sc, const Vdouble *weights=NULL, const tree* constraintTreePtr=NULL);
    virtual MDOUBLE getLogLikelihood() {return _treeLogLikelihood;}

protected:
    distanceMethod *_distM;
    distances2Tree *_dist2et;
    const Vdouble * _weights;
    MDOUBLE _treeLogLikelihood;
	const tree* _constraintTreePtr;
};

class iterativeDistanceSeqs2Tree : public distanceBasedSeqs2Tree {
public:
    iterativeDistanceSeqs2Tree(likeDist &distM, distances2Tree &dist2et, const Vdouble *weights = NULL,
			       const MDOUBLE epsilonLikelihoodImprovement = 0.001, 
			       const MDOUBLE epsilonLikelihoodImprovement4alphaOptimiz = 0.001, 
			       const MDOUBLE epsilonLikelihoodImprovement4BBL = 0.001, 
			       const int maxIterationsBBL = 10);
    virtual ~iterativeDistanceSeqs2Tree() {}
    virtual tree seqs2Tree(const sequenceContainer &sc, const Vdouble *weights=NULL, const tree* constraintTreePtr=NULL) = 0; // iterative
    virtual tree seqs2TreeIterative(const sequenceContainer &sc, const Vdouble *weights=NULL, const tree* constraintTreePtr=NULL) = 0;
    // Start from optimization of branch length and side info for a given initial topology
    virtual tree seqs2TreeIterative(const sequenceContainer &sc, const tree &initTree, const Vdouble *weights=NULL, const tree* constraintTreePtr=NULL) = 0;
    // Start from calculating side info for a given tree and alpha
    virtual tree seqs2TreeIterative(const sequenceContainer &sc, const tree &initTree, MDOUBLE initAlpha, const Vdouble *weights=NULL, const tree* constraintTreePtr=NULL) = 0;
    // Does one bootstrap iteration
  virtual tree seqs2TreeBootstrap(const sequenceContainer &sc, const Vdouble *weights=NULL, const tree* constraintTreePtr=NULL);
    tree getTree() {return _et;}

    // *** handling side info ***

    // Optimize nj tree (optimize alpha, branch lengths, etc.) and produce
    // side info based on the optimized tree
    virtual MDOUBLE optimizeSideInfo(const sequenceContainer &sc, tree &et) = 0;
    // Calculate side info without changing the given tree and alpha
    // (Optimization should be done in here for side info that includes other optimizable parameters
    //  e.g. ML rates, Nu...)
    virtual MDOUBLE calcSideInfoGivenTreeAndAlpha(const sequenceContainer &sc, const tree &et, MDOUBLE alpha) = 0;
    // Copy new side info (based on the new tree) to the "current" side info variable, before the next iteration
    virtual void acceptSideInfo() = 0;
    // Apply the optimized side info into _optimizedSp
    virtual void utilizeSideInfo() = 0;
    virtual void printSideInfo(ostream& out) const = 0;
    MDOUBLE getAlpha() const { return _alpha; }


protected:
    tree seqs2TreeIterativeInternal(const sequenceContainer &sc, bool initSideInfoGiven=false);
    tree seqs2TreeIterativeInternalInitTreeGiven(const sequenceContainer &sc, const tree &initTree);
    tree seqs2TreeIterativeInternalInitTreeGiven(const sequenceContainer &sc, bool initSideInfoGiven, const tree &initTree, MDOUBLE initAlpha);
    void seqs2TreeOneIterationInternal(const sequenceContainer &sc, const bool sideInfoSet);

    MDOUBLE _newTreeLogLikelihood;
    MDOUBLE _epsilonLikelihoodImprovement;
    MDOUBLE _epsilonLikelihoodImprovement4alphaOptimiz;
    MDOUBLE _epsilonLikelihoodImprovement4BBL;
    int _maxIterationsBBL;

    MDOUBLE _alpha;
    MDOUBLE _newAlpha;

    stochasticProcess *_spPtr;
    tree _et, _newTree;
};

class commonAlphaDistanceSeqs2Tree : public iterativeDistanceSeqs2Tree {
public:
    // Given likeDist is assumed to hold a gamma-distribution stochasticProcess
    commonAlphaDistanceSeqs2Tree(likeDist &distM, distances2Tree &dist2et, const Vdouble *weights = NULL,
				 const MDOUBLE epsilonLikelihoodImprovement = 0.001, 
				 const MDOUBLE epsilonLikelihoodImprovement4alphaOptimiz = 0.001, 
				 const MDOUBLE epsilonLikelihoodImprovement4BBL = 0.001, 
				 const int maxIterationsBBL = 50)
	: iterativeDistanceSeqs2Tree(distM, dist2et, weights, epsilonLikelihoodImprovement, epsilonLikelihoodImprovement4alphaOptimiz, epsilonLikelihoodImprovement4BBL, maxIterationsBBL) {}
    virtual ~commonAlphaDistanceSeqs2Tree() {}

    // NOTE! This version calls ITERATIVE seqs2Tree because side info is not given by the user, so we have to generate and optimize it
    virtual tree seqs2Tree(const sequenceContainer &sc, const Vdouble *weights=NULL, const tree* constraintTreePtr=NULL);
    // NOTE! This version is a NON-ITERATIVE version that uses the side info supplied by the user
            tree seqs2Tree(const sequenceContainer &sc, MDOUBLE alpha, const Vdouble *weights=NULL, const tree* constraintTreePtr=NULL);
    // Does one bootstrap iteration
            tree seqs2TreeBootstrap(const sequenceContainer &sc, const MDOUBLE alpha, const Vdouble *weights, const tree* constraintTreePtr=NULL);
    // Explicitly ask for iterations
    virtual tree seqs2TreeIterative(const sequenceContainer &sc, const Vdouble *weights=NULL, const tree* constraintTreePtr=NULL); // homogenous rates will be used for first iteration
            tree seqs2TreeIterative(const sequenceContainer &sc, MDOUBLE initAlpha, const Vdouble *weights=NULL, const tree* constraintTreePtr=NULL);
    virtual tree seqs2TreeIterative(const sequenceContainer &sc, const tree &initTree, const Vdouble *weights=NULL, const tree* constraintTreePtr=NULL);
    virtual tree seqs2TreeIterative(const sequenceContainer &sc, const tree &initTree, MDOUBLE initAlpha, const Vdouble *weights=NULL, const tree* constraintTreePtr=NULL);

    // handling side info
    virtual MDOUBLE optimizeSideInfo(const sequenceContainer &sc, tree &et);
    virtual MDOUBLE calcSideInfoGivenTreeAndAlpha(const sequenceContainer &sc, const tree &et, MDOUBLE alpha);
    virtual void acceptSideInfo();
    virtual void utilizeSideInfo();
    virtual void printSideInfo(ostream& out) const;
    void setSideInfo(const MDOUBLE alpha);
    MDOUBLE getSideInfo() const;
};

class rate4siteDistanceSeqs2Tree : public iterativeDistanceSeqs2Tree {
public:
    rate4siteDistanceSeqs2Tree(givenRatesMLDistance &distM, distances2Tree &dist2et, const Vdouble *weights = NULL,
			       const MDOUBLE epsilonLikelihoodImprovement = 0.001, 
			       const MDOUBLE epsilonLikelihoodImprovement4alphaOptimiz = 0.001, 
			       const MDOUBLE epsilonLikelihoodImprovement4BBL = 0.001, 
			       const int maxIterationsBBL = 50)
	: iterativeDistanceSeqs2Tree(distM, dist2et, weights, epsilonLikelihoodImprovement, epsilonLikelihoodImprovement4alphaOptimiz, epsilonLikelihoodImprovement4BBL, maxIterationsBBL) {}
    virtual ~rate4siteDistanceSeqs2Tree() {}

    // NOTE! This version calls ITERATIVE seqs2Tree because side info is not given by the user, so we have to generate and optimize it
    virtual tree seqs2Tree(const sequenceContainer &sc, const Vdouble *weights = NULL, const tree* constraintTreePtr=NULL);
    // NOTE! This version is a NON-ITERATIVE version that uses the side info supplied by the user
            tree seqs2Tree(const sequenceContainer &sc, const Vdouble &rates, const Vdouble *weights = NULL, const tree* constraintTreePtr=NULL);
    // Does one bootstrap iteration
            tree seqs2TreeBootstrap(const sequenceContainer &sc, const Vdouble &rates, const Vdouble *weights, const tree* constraintTreePtr=NULL);
    // Explicitly ask for iterations
    virtual tree seqs2TreeIterative(const sequenceContainer &sc, const Vdouble *weights=NULL, const tree* constraintTreePtr=NULL); // homogenous rates will be used for first iteration
            tree seqs2TreeIterative(const sequenceContainer &sc, const Vdouble &initRates, const Vdouble *weights=NULL, const tree* constraintTreePtr=NULL);
    virtual tree seqs2TreeIterative(const sequenceContainer &sc, const tree &initTree, const Vdouble *weights=NULL, const tree* constraintTreePtr=NULL);
    virtual tree seqs2TreeIterative(const sequenceContainer &sc, const tree &initTree, MDOUBLE initAlpha, const Vdouble *weights=NULL, const tree* constraintTreePtr=NULL);

    // handling side info
    virtual MDOUBLE optimizeSideInfo(const sequenceContainer &sc, tree &et);
    virtual MDOUBLE calcSideInfoGivenTreeAndAlpha(const sequenceContainer &sc, const tree &et, MDOUBLE alpha);
    virtual void acceptSideInfo();
    virtual void utilizeSideInfo();
    virtual void printSideInfo(ostream& out) const;
    void setSideInfo(const Vdouble &rates);
    const Vdouble& getSideInfo() const;

private:
    Vdouble _rates;
    Vdouble _newRates;
};

class posteriorDistanceSeqs2Tree : public iterativeDistanceSeqs2Tree {
public:
    posteriorDistanceSeqs2Tree(posteriorDistance &distM, distances2Tree &dist2et, const Vdouble *weights = NULL,
			       const MDOUBLE epsilonLikelihoodImprovement = 0.001, 
			       const MDOUBLE epsilonLikelihoodImprovement4alphaOptimiz = 0.001, 
			       const MDOUBLE epsilonLikelihoodImprovement4BBL = 0.001, 
			       const int maxIterationsBBL = 50)
	: iterativeDistanceSeqs2Tree(distM, dist2et, weights, epsilonLikelihoodImprovement, epsilonLikelihoodImprovement4alphaOptimiz, epsilonLikelihoodImprovement4BBL, maxIterationsBBL) {}
    virtual ~posteriorDistanceSeqs2Tree() {}

    // NOTE! This version calls ITERATIVE seqs2Tree because side info is not given by the user, so we have to generate and optimize it
    virtual tree seqs2Tree(const sequenceContainer &sc, const Vdouble *weights = NULL, const tree* constraintTreePtr=NULL);
    // NOTE! This version is a NON-ITERATIVE version that uses the side info supplied by the user
            tree seqs2Tree(const sequenceContainer &sc, const VVdoubleRep &posterior, const Vdouble *weights = NULL, const tree* constraintTreePtr=NULL);
    // Does one bootstrap iteration
            tree seqs2TreeBootstrap(const sequenceContainer &sc, const VVdoubleRep &posterior, const Vdouble *weights=NULL, const tree* constraintTreePtr=NULL);
    // Explicitly ask for iterations
    virtual tree seqs2TreeIterative(const sequenceContainer &sc, const Vdouble *weights=NULL, const tree* constraintTreePtr=NULL); // homogenous rates will be used for first iteration
            tree seqs2TreeIterative(const sequenceContainer &sc, MDOUBLE initAlpha, const VVdoubleRep &initPosterior, const Vdouble *weights=NULL, const tree* constraintTreePtr=NULL);
    virtual tree seqs2TreeIterative(const sequenceContainer &sc, const tree &initTree, const Vdouble *weights=NULL, const tree* constraintTreePtr=NULL);
    virtual tree seqs2TreeIterative(const sequenceContainer &sc, const tree &initTree, MDOUBLE initAlpha, const Vdouble *weights=NULL, const tree* constraintTreePtr=NULL);
            tree seqs2TreeIterative(const sequenceContainer &sc, const tree &initTree, MDOUBLE initAlpha, const VVdoubleRep &initPosterior, const Vdouble *weights=NULL, const tree* constraintTreePtr=NULL);

    // handling side info
    virtual MDOUBLE optimizeSideInfo(const sequenceContainer &sc, tree &et);
    virtual MDOUBLE calcSideInfoGivenTreeAndAlpha(const sequenceContainer &sc, const tree &et, MDOUBLE alpha);
    virtual void acceptSideInfo();
    virtual void utilizeSideInfo();
    virtual void printSideInfo(ostream& out) const;
    void setSideInfo(const VVdoubleRep &posterior);
    const VVdoubleRep& getSideInfo() const;

private:
    VVdoubleRep _posterior;
    VVdoubleRep _newPosterior;
};

#endif
