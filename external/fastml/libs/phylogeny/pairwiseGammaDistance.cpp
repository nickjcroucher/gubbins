// $Id: pairwiseGammaDistance.cpp 962 2006-11-07 15:13:34Z privmane $

#include "pairwiseGammaDistance.h"
#include "numRec.h"
#include "countTableComponent.h"
#include "likeDist.h"
#include "uniDistribution.h"
#include <cmath>

// Local utility functions
MDOUBLE pairwiseGammaDistance::giveInitialGuessOfDistance(
    const sequence& s1,
    const sequence& s2,
    const vector<MDOUBLE>  * weights,
    MDOUBLE* score) const {
    uniDistribution ud;
    stochasticProcess uniSp(&ud,_sp.getPijAccelerator());
    likeDist ld(uniSp);
    return (ld.giveDistance(s1,s2,weights,score));
}

class C_eval_gammaMLAlpha{ 
private:
    const stochasticProcess& _sp;
    const sequence& _s1;
    const sequence& _s2;
    const MDOUBLE _distance;
    const Vdouble* _weights;
    //  const VVdouble& _posteriorProb; // pos, rate
public:
    C_eval_gammaMLAlpha(const stochasticProcess& sp,
			const sequence& s1,
			const sequence& s2,
			const MDOUBLE distance,
			//		      const VVdouble& posteriorProb,
			const Vdouble  * weights):  _sp(sp),
						    _s1(s1),
						    _s2(s2),
						    _distance(distance),
						    _weights(weights) 
	//						  _posteriorProb(posteriorProb)
	{};

    // this cast is required as the distribution within the
    // stochasticProcess is kept as the parent "distribution" class that
    // knows nothing of Alpha
    void setAlpha(MDOUBLE alpha) {
	(static_cast<gammaDistribution*>(_sp.distr()))->setAlpha(alpha);
    }


    MDOUBLE operator() (MDOUBLE alpha) {
	setAlpha(alpha);
	MDOUBLE likelihood = likeDist::evalLikelihoodForDistance(_sp,_s1,_s2,_distance,_weights);
	LOG(11,<<"check alpha="<<alpha<<", bl="<<_distance<<" gives "<<likelihood<<endl);
	return -likelihood;
    };
};

// returns the best alpha for a given distance
MDOUBLE pairwiseGammaDistance::optimizeAlphaFixedDist(const sequence& s1,
						      const sequence& s2,
						      stochasticProcess & sp,
						      const MDOUBLE branchL,
						      const vector<MDOUBLE>  * weights,
						      MDOUBLE* score) const { // changes sp.
    MDOUBLE bestA=0.0;
    MDOUBLE bestQ=0.0;
    const MDOUBLE upperBoundOnAlpha = 15.0;
    const MDOUBLE epsilonAlphaOptimization = 0.01;
    const MDOUBLE cx=upperBoundOnAlpha;// left, midle, right limit on alpha
    const MDOUBLE bx=cx*0.3;
    const MDOUBLE ax=0.0;

	
    bestQ = -brent(ax,bx,cx,
		   C_eval_gammaMLAlpha(sp,s1,s2,branchL,weights),
		   epsilonAlphaOptimization,
		   &bestA);
    (static_cast<gammaDistribution*>(sp.distr()))->setAlpha(bestA);
    if (score) *score = bestQ;
    return bestA;
}

class C_evalAlphaForPairOfSeq{
private:
    const countTableComponentGam& _ctc;
    stochasticProcess& _sp;
    const MDOUBLE _branchL;
public:
    C_evalAlphaForPairOfSeq(const countTableComponentGam& ctc,
			    const MDOUBLE branchL,
			    stochasticProcess& sp):_ctc(ctc), _sp(sp), _branchL(branchL) {};

    MDOUBLE operator() (MDOUBLE alpha) {
	(static_cast<gammaDistribution*>(_sp.distr()))->setAlpha(alpha);
	C_evalLikeDist cev(_ctc,_sp);
	MDOUBLE L=cev(_branchL);
	LOG(10,<<"check alpha="<<alpha<<", bl="<<_branchL<<" gives "<<L<<endl);
	return L;
    };
};

// returns the best alpha for a given distance
MDOUBLE pairwiseGammaDistance::optimizeAlphaFixedDist(stochasticProcess & sp,
						      const countTableComponentGam & ctc,
						      const MDOUBLE branchL,
						      const vector<MDOUBLE>  * weights,
						      MDOUBLE* score) const { // changes sp.
    MDOUBLE bestA=0.0;
    MDOUBLE bestQ=0.0;
    const MDOUBLE upperBoundOnAlpha = 15.0;
    const MDOUBLE epsilonAlphaOptimization = 0.01;
    const MDOUBLE cx=upperBoundOnAlpha;// left, midle, right limit on alpha
    const MDOUBLE bx=cx*0.3;
    const MDOUBLE ax=0.0;

	
    bestQ = -brent(ax,bx,cx,
		   C_evalAlphaForPairOfSeq(ctc,branchL,sp),
		   epsilonAlphaOptimization,
		   &bestA);
    (static_cast<gammaDistribution*>(sp.distr()))->setAlpha(bestA);
    if (score) *score = bestQ;
    return bestA;
}

const MDOUBLE pairwiseGammaDistance::giveDistance(const sequence& s1,
					    const sequence& s2,
					    const vector<MDOUBLE>  * weights,
					    MDOUBLE* score,
					    MDOUBLE* alpha) const {
  
    MDOUBLE resL = 0.0;
    MDOUBLE currentDistance = giveInitialGuessOfDistance(s1,s2,weights,&resL);
	
    countTableComponentGam ctc; // from technical reasons.
	
    stochasticProcess tmpSp(_sp);
	
    const int maxIter = 30;
    MDOUBLE newDist = 0.0;
    MDOUBLE lastBestAlpha = 0.0;
    for (int i=0; i < maxIter; ++i) {
	lastBestAlpha = optimizeAlphaFixedDist(s1, s2, tmpSp, currentDistance, weights, &resL); // changes sp.
	LOG(8,<<"lastBestAlpha="<<lastBestAlpha<<"("<<"\t L="<<resL<<"\t");
	likeDist tmpld(tmpSp); // we must create a new ld, that will include the stochastic process with the new alpha
	newDist = tmpld.giveDistance(s1, s2, weights, &resL);
	LOG(8,<<"dist="<<newDist<<"(L="<<resL<<")"<<endl);
	if (fabs(newDist-currentDistance)<_toll) break;
	currentDistance = newDist;
    }
    if (score) *score = resL;
    if (alpha) *alpha = lastBestAlpha;
    assert (newDist >=0);
    return newDist;
}

