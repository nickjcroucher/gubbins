// $Id: fromCountTableComponentToDistance2Codon.cpp 950 2006-10-19 12:12:34Z eyalprivman $

#include "fromCountTableComponentToDistance2Codon.h"
#include "likeDist2Codon.h"
#include "likeDist.h"
#include <cassert>

fromCountTableComponentToDistance2Codon::fromCountTableComponentToDistance2Codon(
		const countTableComponentGam& ctc,
		const vector<stochasticProcess> &spVec,
		const MDOUBLE toll,
		const MDOUBLE brLenIntialGuess ) : _spVec(spVec), _ctc(ctc) {
	_distance =brLenIntialGuess ;//0.03;
	_toll = toll;
}

void fromCountTableComponentToDistance2Codon::computeDistance() {
	likeDist2Codon likeDist1(_spVec,_toll);
	MDOUBLE initGuess = _distance;
	_distance = likeDist1.giveDistance(_ctc,_likeDistance,initGuess);
	assert(_distance>=0);
}
