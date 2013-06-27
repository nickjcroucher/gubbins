// $Id: fromCountTableComponentToDistance.h 4742 2008-08-19 17:40:56Z cohenofi $

#ifndef ___FROM_COUNT_TABLE_COMPONENT_TO_DISTANCE
#define ___FROM_COUNT_TABLE_COMPONENT_TO_DISTANCE

#include "definitions.h"
#include "countTableComponent.h"
#include "stochasticProcess.h"
#include "unObservableData.h"

static const MDOUBLE startingGuessForTreeBrLen = 0.029;

class fromCountTableComponentToDistance {

public:
	explicit fromCountTableComponentToDistance(
		const countTableComponentGam& ctc,
		const stochasticProcess &sp,
		const MDOUBLE toll,
		const MDOUBLE brLenIntialGuess,				// =startingGuessForTreeBrLen
		unObservableData* unObservableData_p = NULL);	// a class used to for presence/absence

	void computeDistance();// return the likelihood
	MDOUBLE getDistance() { return _distance;} // return the distance.
	MDOUBLE getLikeDistance() { return _likeDistance;} // return the distance.
private:
	const stochasticProcess & _sp;
	const countTableComponentGam& _ctc;
	MDOUBLE _toll;
	MDOUBLE _distance;
	MDOUBLE _likeDistance;
	unObservableData*  _unObservableData_p;
	int alphabetSize() {return _ctc.alphabetSize();}
};

#endif

