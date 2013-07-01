// $Id: fromCountTableComponentToDistance2Codon.h 950 2006-10-19 12:12:34Z eyalprivman $

#ifndef ___FROM_COUNT_TABLE_COMPONENT_TO_DISTANCE_2_CODON
#define ___FROM_COUNT_TABLE_COMPONENT_TO_DISTANCE_2_CODON

#include "definitions.h"
#include "countTableComponent.h"
#include "stochasticProcess.h"

static const MDOUBLE startingGuessForTreeBrLen = 0.029;

class fromCountTableComponentToDistance2Codon {

public:
	explicit fromCountTableComponentToDistance2Codon(
		const countTableComponentGam& ctc,
		const vector<stochasticProcess> &spVec,
		const MDOUBLE toll,
		const MDOUBLE brLenIntialGuess);// =startingGuessForTreeBrLen

	void computeDistance();// return the likelihood
	MDOUBLE getDistance() { return _distance;} // return the distance.
	MDOUBLE getLikeDistance() { return _likeDistance;} // return the distance.
private:
	const vector<stochasticProcess> & _spVec;
	const countTableComponentGam& _ctc;
	MDOUBLE _toll;
	MDOUBLE _distance;
	MDOUBLE _likeDistance;
	int alphabetSize() {return _ctc.alphabetSize();}
};

#endif

