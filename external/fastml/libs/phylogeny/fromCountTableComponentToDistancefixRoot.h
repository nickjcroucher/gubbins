// $Id: fromCountTableComponentToDistance.h 4471 2008-07-17 15:38:50Z cohenofi $

#ifndef ___FROM_COUNT_TABLE_COMPONENT_TO_DISTANCE__FIX_ROOT
#define ___FROM_COUNT_TABLE_COMPONENT_TO_DISTANCE__FIX_ROOT

#include "definitions.h"
#include "countTableComponent.h"
#include "stochasticProcess.h"
#include "unObservableData.h"

static const MDOUBLE startingGuessForTreeBrLen = 0.029;

class fromCountTableComponentToDistancefixRoot {

public:
	explicit fromCountTableComponentToDistancefixRoot(
		const vector<countTableComponentGam>& ctc,
		const stochasticProcess &sp,
		const MDOUBLE toll,
		const MDOUBLE brLenIntialGuess,      // =startingGuessForTreeBrLen
		unObservableData*  unObservableData_p);  

	void computeDistance();// return the likelihood
	MDOUBLE getDistance() { return _distance;} // return the distance.
	MDOUBLE getLikeDistance() { return _likeDistance;} // return the distance.
private:
	const stochasticProcess & _sp;
	const vector<countTableComponentGam>& _ctc;			//_ctc[letterAtRoot][rate][alph][alph]
	MDOUBLE _toll;
	MDOUBLE _distance;
	MDOUBLE _likeDistance;
	unObservableData*  _unObservableData_p;

//	int alphabetSize() {return _ctc.alphabetSize();}
	int alphabetSize() {return _ctc[0].alphabetSize();}
};

#endif

