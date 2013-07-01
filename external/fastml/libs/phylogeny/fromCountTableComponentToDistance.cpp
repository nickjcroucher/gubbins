// $Id: fromCountTableComponentToDistance.cpp 4742 2008-08-19 17:40:56Z cohenofi $

#include "fromCountTableComponentToDistance.h"
#include "likeDist.h"
#include <cassert>

fromCountTableComponentToDistance::fromCountTableComponentToDistance(
		const countTableComponentGam& ctc,
		const stochasticProcess &sp,
		const MDOUBLE toll,
		const MDOUBLE brLenIntialGuess,
		unObservableData* unObservableData_p) : _sp(sp), _ctc(ctc),_unObservableData_p(unObservableData_p) {
	_distance = brLenIntialGuess ;//0.03;
	_toll = toll;
}

void fromCountTableComponentToDistance::computeDistance() {
	MDOUBLE maxPairwiseDistance = 5.0; // The default
	likeDist likeDist1(_sp,_toll,maxPairwiseDistance,_unObservableData_p);
	MDOUBLE initGuess = _distance;
	_distance = likeDist1.giveDistance(_ctc,_likeDistance,initGuess);
	assert(_distance>=0);
}
