// $Id: givenRatesMLDistance.h 962 2006-11-07 15:13:34Z privmane $

#ifndef ___GIVEN_RATES_ML_DISTANCE_H
#define ___GIVEN_RATES_ML_DISTANCE_H

#include "definitions.h"
#include "countTableComponent.h"
#include "likeDist.h"
#include "stochasticProcess.h"
#include "logFile.h"
#include <cmath>
using namespace std;

class givenRatesMLDistance : public likeDist {
public:
  explicit givenRatesMLDistance(const stochasticProcess& sp,
				const Vdouble& rates,
				const MDOUBLE toll =0.0001,
				const MDOUBLE maxPairwiseDistance = 5.0
				)
    :  likeDist(sp, toll,maxPairwiseDistance),_rates(rates) {}

  explicit givenRatesMLDistance(stochasticProcess& sp,
				const Vdouble& rates,
				const MDOUBLE toll =0.0001,
				const MDOUBLE maxPairwiseDistance = 5.0
				)
    :  likeDist(sp, toll,maxPairwiseDistance),_rates(rates) {}

  explicit givenRatesMLDistance(const stochasticProcess& sp,
				const MDOUBLE toll =0.0001,
				const MDOUBLE maxPairwiseDistance = 5.0
				)
    :  likeDist(sp, toll,maxPairwiseDistance),_rates(0) {}

  explicit givenRatesMLDistance(stochasticProcess& sp,
				const MDOUBLE toll =0.0001,
				const MDOUBLE maxPairwiseDistance = 5.0
				)
    :  likeDist(sp, toll,maxPairwiseDistance),_rates(0) {}

  givenRatesMLDistance(const givenRatesMLDistance& other):
	likeDist(static_cast<likeDist>(other)), _rates(other._rates) {}

  virtual givenRatesMLDistance* clone() const {return new givenRatesMLDistance(*this);}

  void setRates(const Vdouble &rates) {_rates = rates;}

  // Returns the estimated ML distance between the 2 sequences.
  // if score is given, it will be assigned the log-likelihood.
  const MDOUBLE giveDistance(const sequence& s1,
			     const sequence& s2,
			     const vector<MDOUBLE>  * weights,
			     MDOUBLE* score=NULL) const;

private:
	Vdouble _rates;
};

#endif

