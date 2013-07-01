// 	$Id: fromCountTableComponentToDistance2USSRV.h 962 2006-11-07 15:13:34Z privmane $	

#ifndef ___FROM_COUNT_TABLE_COMPONENT_TO_DISTANCE_2_USSRV
#define ___FROM_COUNT_TABLE_COMPONENT_TO_DISTANCE_2_USSRV

#include "definitions.h"
#include "countTableComponent.h"
#include "stochasticProcess.h"
#include "ussrvModel.h"
#include "likeDist2USSRV.h"

static const MDOUBLE startingGuessForTreeBrLen = 0.029;

class fromCountTableComponentToDistance2USSRV {

public:
	explicit fromCountTableComponentToDistance2USSRV(
		const countTableComponentGam& ctcBase,
		const countTableComponentHom& ctcSSRV,
		const ussrvModel& model,
		MDOUBLE toll,
		MDOUBLE brLenIntialGuess);// =startingGuessForTreeBrLen

	void computeDistance();// return the likelihood
	MDOUBLE getDistance() { return _distance;} // return the distance.
	MDOUBLE getLikeDistance() { return _likeDistance;} // return the distance.

private:
	const ussrvModel & _model;
	const countTableComponentGam& _ctcBase;
	const countTableComponentHom& _ctcSSRV;
	MDOUBLE _toll;
	MDOUBLE _distance;
	MDOUBLE _likeDistance;
//	int alphabetSize() {return _ctc.alphabetSize();}
};

#endif //___FROM_COUNT_TABLE_COMPONENT_TO_DISTANCE_2_USSRV

