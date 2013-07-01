// $Id: fromCountTableComponentToDistanceProp.h 962 2006-11-07 15:13:34Z privmane $

#ifndef ___FROM_COUNT_TABLE_COMPONENT_TO_DISTANCE_PROP
#define ___FROM_COUNT_TABLE_COMPONENT_TO_DISTANCE_PROP

#include "definitions.h"
#include "countTableComponent.h"
#include "stochasticProcess.h"


class fromCountTableComponentToDistanceProp {

public:
	explicit fromCountTableComponentToDistanceProp(
		const vector<countTableComponentGam>& ctc,
		const vector<stochasticProcess> &sp,
		const MDOUBLE toll,
		const MDOUBLE brLenIntialGuess = 0.029);// =startingGuessForTreeBrLen

	void computeDistance();// return the likelihood
	MDOUBLE getDistance() { return _distance;} // return the distance.
	MDOUBLE getLikeDistance() { return _likeDistance;} // return the distance.
private:
	const vector<stochasticProcess> & _sp;
	const vector<countTableComponentGam>& _ctc;
	MDOUBLE _toll;
	MDOUBLE _distance;
	MDOUBLE _likeDistance;
	int alphabetSize() {return (_ctc.empty()?0:_ctc[0].alphabetSize());}
};

#endif

