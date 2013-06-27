// $Id: posteriorDistance.h 1752 2007-02-26 14:01:09Z itaymay $


#ifndef POSTERIOR_DISTANCE_H
#define POSTERIOR_DISTANCE_H

#include "likeDist.h"
#include "stochasticProcess.h"
#include "definitions.h"
#include "sequence.h"
#include "gammaDistribution.h"
#include "logFile.h"

#include <cmath>
using namespace std;

class posteriorDistance : public likeDist {
public:
    explicit posteriorDistance(const stochasticProcess & sp,
			       const VVdoubleRep & posteriorProb, // pos * rate
			       const MDOUBLE toll =0.0001,
			       const MDOUBLE maxPairwiseDistance = 5.0);

    explicit posteriorDistance(stochasticProcess & sp,
			       const VVdoubleRep & posteriorProb, // pos * rate
			       const MDOUBLE toll =0.0001,
			       const MDOUBLE maxPairwiseDistance = 5.0);

    explicit posteriorDistance(const stochasticProcess & sp,
			       const MDOUBLE toll =0.0001,
			       const MDOUBLE maxPairwiseDistance = 5.0);

    explicit posteriorDistance(stochasticProcess & sp,
			       const MDOUBLE toll =0.0001,
			       const MDOUBLE maxPairwiseDistance = 5.0);
    posteriorDistance(const posteriorDistance& other);
  virtual posteriorDistance* clone() const {return new posteriorDistance(*this);}

    // distance is computed based on the posterior probability
    const MDOUBLE giveDistance(const sequence& s1,
			       const sequence& s2,
			       const vector<MDOUBLE>  * weights,
			       MDOUBLE* score=NULL) const;
  
    MDOUBLE giveDistanceOptAlphaForEachPairOfSequences(const sequence& s1,
						       const sequence& s2,
						       const vector<MDOUBLE>  * weights,
						       MDOUBLE* score=NULL,
						       MDOUBLE* alpha=NULL) const;

    MDOUBLE giveDistanceOptAlphaForPairOfSequences(const sequence& s1,
						   const sequence& s2,
						   const vector<MDOUBLE>  * weights,
						   MDOUBLE* score,
						   MDOUBLE* alpha) const;

    void setPosterior(VVdoubleRep posteriorProb) {_posteriorProb = posteriorProb;}
    void setAlpha(MDOUBLE alpha) {
	(static_cast<gammaDistribution*>(_sp.distr()))->setAlpha(alpha);
    }

private:
    VVdoubleRep _posteriorProb;
    MDOUBLE giveInitialGuessOfDistance(const sequence& s1,
				       const sequence& s2,
				       const vector<MDOUBLE>  * weights,
				       MDOUBLE* score) const;
};



#endif
