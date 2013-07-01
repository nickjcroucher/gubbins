// $Id: ssrvDistanceSeqs2Tree.h 962 2006-11-07 15:13:34Z privmane $

#ifndef ___SSRV_DISTANCE_SEQS2TREE
#define ___SSRV_DISTANCE_SEQS2TREE

#include "distanceBasedSeqs2Tree.h"
#include "tree.h"

/* class ssrvDistanceSeqs2Tree 
A type of distance-based tree reconstruction method like the iterative
method commonAlphaDistanceSeqs2Tree, but using a model with SSRV
(Site-Specific Rate Variation, AKA covarion model).  Compared to
commonAlphaDistanceSeqs2Tree, we change the distance method to use an
SSRV model, and in the optimizations we estimate ni in addition to
alpha.
*/
class ssrvDistanceSeqs2Tree : public iterativeDistanceSeqs2Tree {
public:
    // Given likeDist is assumed to hold a gamma-distribution, SSRV stochasticProcess
    ssrvDistanceSeqs2Tree(likeDist &distM, distances2Tree &dist2et, const Vdouble *weights = NULL,
			  const MDOUBLE epsilonLikelihoodImprovement = 0.001, 
			  const MDOUBLE epsilonLikelihoodImprovement4paramOptimiz = 0.001, 
			  const MDOUBLE epsilonLikelihoodImprovement4BBL = 0.001, 
			  const int maxIterationsBBL = 50)
	: iterativeDistanceSeqs2Tree(distM, dist2et, weights, epsilonLikelihoodImprovement, epsilonLikelihoodImprovement4paramOptimiz, epsilonLikelihoodImprovement4BBL, maxIterationsBBL) {}
    virtual ~ssrvDistanceSeqs2Tree () {}

    // Datastruct for handling side info for the SSRV model (used as return value)
    struct alphaAndNu {
	MDOUBLE alpha;
	MDOUBLE nu;
	alphaAndNu(){}
	alphaAndNu(MDOUBLE setAlpha, MDOUBLE setNu) : alpha(setAlpha), nu(setNu) {}
    };

    // NOTE! This version calls ITERATIVE seqs2Tree because side info is not given by the user, so we have to generate and optimize it
    virtual tree seqs2Tree(const sequenceContainer &sc, const Vdouble *weights=NULL, const tree* constraintTreePtr=NULL);
    // NOTE! This version is a NON-ITERATIVE version that uses the side info supplied by the user
            tree seqs2Tree(const sequenceContainer &sc, MDOUBLE alpha, MDOUBLE nu, const Vdouble *weights=NULL, const tree* constraintTreePtr=NULL);
    // Does one bootstrap iteration
            tree seqs2TreeBootstrap(const sequenceContainer &sc, const MDOUBLE alpha, MDOUBLE nu, const Vdouble *weights, const tree* constraintTreePtr=NULL);
    // Explicitly ask for iterations
    virtual tree seqs2TreeIterative(const sequenceContainer &sc, const Vdouble *weights=NULL, const tree* constraintTreePtr=NULL); // homogenous rates will be used for first iteration
            tree seqs2TreeIterative(const sequenceContainer &sc, MDOUBLE initAlpha, MDOUBLE initNu, const Vdouble *weights=NULL, const tree* constraintTreePtr=NULL);
    virtual tree seqs2TreeIterative(const sequenceContainer &sc, const tree &initTree, const Vdouble *weights=NULL, const tree* constraintTreePtr=NULL);
    virtual tree seqs2TreeIterative(const sequenceContainer &sc, const tree &initTree, MDOUBLE initAlpha, const Vdouble *weights=NULL, const tree* constraintTreePtr=NULL);
            tree seqs2TreeIterative(const sequenceContainer &sc, const tree &initTree, MDOUBLE initAlpha, MDOUBLE initNu, const Vdouble *weights=NULL, const tree* constraintTreePtr=NULL);

    // handling side info
    virtual MDOUBLE optimizeSideInfo(const sequenceContainer &sc, tree &et);
    virtual MDOUBLE calcSideInfoGivenTreeAndAlpha(const sequenceContainer &sc, const tree &et, MDOUBLE alpha);
    virtual void acceptSideInfo();
    virtual void utilizeSideInfo();
    virtual void printSideInfo(ostream& out) const;
    void setSideInfo(const MDOUBLE alpha, MDOUBLE nu);
    alphaAndNu getSideInfo() const;

protected:
    MDOUBLE _nu;
    MDOUBLE _newNu;
};

#endif
