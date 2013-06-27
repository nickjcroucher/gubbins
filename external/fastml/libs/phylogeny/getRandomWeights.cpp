// $Id: getRandomWeights.cpp 962 2006-11-07 15:13:34Z privmane $

#include "getRandomWeights.h"
#include "talRandom.h"



void swapRand(Vdouble& weights) {
	int j;
    int i = talRandom::giveIntRandomNumberBetweenZeroAndEntry(weights.size());
	do {
        j = talRandom::giveIntRandomNumberBetweenZeroAndEntry(weights.size());
	} while ( weights[j] <= 0 );

	weights[i]++;
    weights[j]--;
}

void getRandomWeights::randomWeights(Vdouble& weights,
								const MDOUBLE expectedNumberOfSwapsPerPosition) {
	// note that some positions will change more than once, and some won't.
	// thus the second argument is an average of sites swaped
	int i;
	const double DefaultWeight = 1;
	for (i=0; i< weights.size(); ++i) weights[i] = DefaultWeight;

	for ( i = 0 ; i < expectedNumberOfSwapsPerPosition*weights.size() ; ++i ) {  
        swapRand(weights);							
	}
}

void getRandomWeights::standardBPWeights(Vdouble& weights) {
	int i;
	for (i=0; i< weights.size(); ++i) weights[i] = 0.0;
	for (i=0; i< weights.size(); ++i) {
	    int k = talRandom::giveIntRandomNumberBetweenZeroAndEntry(weights.size());
		weights[k]++;
	}
}

#define MIN_WEIGHT (0.00001)
void getRandomWeights::randomWeightsGamma(Vdouble& weights,
					  const MDOUBLE temperature) {
	int i;
	const double oneOverT = 1.0/temperature;
	for (i=0; i< weights.size(); ++i) {
		weights[i] = talRandom::SampleGamma(oneOverT,oneOverT);
		if (weights[i]<MIN_WEIGHT) {
			weights[i] = MIN_WEIGHT;
		}
	}
}

