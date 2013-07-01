// $Id: getRandomWeights.h 962 2006-11-07 15:13:34Z privmane $

#ifndef __GET_RANDOM_WEIGHTS
#define __GET_RANDOM_WEIGHTS

#include "definitions.h"


class getRandomWeights {
public:
	// this function starts with a vector of weights like that (1,1,1,1,1,1,...1)
	// it then take two positions by random
	// add 1 to the first, and substract 1 from the second.
	// if it can not substract 1 from the second, it draw a new "second"
	static void randomWeights(Vdouble& weights,
				  const MDOUBLE expectedNumberOfSwapsPerPosition);

	// a position is chosen randomly and the weight of this position is
	// sampled from a gamma distribution with parameters alpha = 1/temperature
	// and beta = 1/temperature.
	static void randomWeightsGamma(Vdouble& weights,
				       const MDOUBLE temperature);

	// this function starts with a vector of weights like that (0,0,0,...,0)
	// a position is chosen randomly and the weight of this position
	// is increased by 1. This process is repeated weights.size() times.
	static void standardBPWeights(Vdouble& weights);
};

#endif 

