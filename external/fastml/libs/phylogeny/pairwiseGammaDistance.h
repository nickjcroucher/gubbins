// $Id: pairwiseGammaDistance.h 962 2006-11-07 15:13:34Z privmane $

#ifndef PAIRWISE_GAMMA_DISTANCE_H
#define PAIRWISE_GAMMA_DISTANCE_H

#include "likeDist.h"
#include "stochasticProcess.h"
#include "definitions.h"
#include "sequence.h"
#include "gammaDistribution.h"
#include "logFile.h"

#include <cmath>
using namespace std;

// Finds ML distance with a gamma-ASRV stochasticProcess for a pair of
// sequences while optimizing the alpha parameter for the given pair of
// sequences.  
// Was called "njGamma::giveDistanceOptAlphaForPairOfSequences"
class pairwiseGammaDistance : public likeDist {
public:
    explicit pairwiseGammaDistance(const stochasticProcess & sp,
				   const MDOUBLE toll =0.0001,
				   const MDOUBLE maxPairwiseDistance = 5.0)
	: likeDist(sp,toll,maxPairwiseDistance) {}

    explicit pairwiseGammaDistance(stochasticProcess & sp,
				   const MDOUBLE toll =0.0001,
				   const MDOUBLE maxPairwiseDistance = 5.0)
	: likeDist(sp,toll,maxPairwiseDistance) {}

    const MDOUBLE giveDistance(const sequence& s1,
			       const sequence& s2,
			       const vector<MDOUBLE> * weights = NULL,
			       MDOUBLE* score=NULL,
			       MDOUBLE* alpha=NULL) const;
  
  virtual pairwiseGammaDistance* clone() const {return new pairwiseGammaDistance(*this);}

    void setAlpha(MDOUBLE alpha) {
	(static_cast<gammaDistribution*>(_sp.distr()))->setAlpha(alpha);
    }

  
protected:
    MDOUBLE giveInitialGuessOfDistance(const sequence& s1,
				       const sequence& s2,
				       const vector<MDOUBLE>  * weights,
				       MDOUBLE* score) const;
    MDOUBLE optimizeAlphaFixedDist(const sequence& s1,
				   const sequence& s2,
				   stochasticProcess & sp,
				   const MDOUBLE branchL,
				   const vector<MDOUBLE>  * weights,
				   MDOUBLE* score=NULL) const;
    MDOUBLE optimizeAlphaFixedDist(stochasticProcess & sp,
				   const countTableComponentGam & ctc,
				   const MDOUBLE branchL,
				   const vector<MDOUBLE>  * weights,
				   MDOUBLE* score=NULL) const;
};

#endif
