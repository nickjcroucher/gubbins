// $Id: fromCountTableComponentToDistance.cpp 4471 2008-07-17 15:38:50Z cohenofi $

#include "fromCountTableComponentToDistancefixRoot.h"
#include "likeDistfixRoot.h"
#include <cassert>

fromCountTableComponentToDistancefixRoot::fromCountTableComponentToDistancefixRoot(
		const vector<countTableComponentGam>& ctc,
		const stochasticProcess &sp,
		const MDOUBLE toll,
		const MDOUBLE brLenIntialGuess,
		unObservableData*  unObservableData_p) 
		: _sp(sp), _ctc(ctc) {
	_distance =brLenIntialGuess ;//0.03;
	_toll = toll;
	_unObservableData_p = unObservableData_p;

}

void fromCountTableComponentToDistancefixRoot::computeDistance() {
	MDOUBLE maxPairwiseDistance = 5.0; // The default
	likeDistfixRoot likeDist1(_sp,_toll,maxPairwiseDistance,_unObservableData_p);
	MDOUBLE initGuess = _distance;
	_distance = likeDist1.giveDistance(_ctc,_likeDistance,initGuess);
	assert(_distance>=0);
}
