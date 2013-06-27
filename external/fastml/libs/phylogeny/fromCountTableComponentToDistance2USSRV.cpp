// 	$Id: fromCountTableComponentToDistance2USSRV.cpp 962 2006-11-07 15:13:34Z privmane $	

#include "fromCountTableComponentToDistance2USSRV.h"
#include "likeDist.h"
#include <cassert>

fromCountTableComponentToDistance2USSRV::fromCountTableComponentToDistance2USSRV(
		const countTableComponentGam& ctcBase,
		const countTableComponentHom& ctcSSRV,
		const ussrvModel &model,
		MDOUBLE toll,
		MDOUBLE brLenIntialGuess ) : _model(model), _ctcBase(ctcBase), _ctcSSRV(ctcSSRV) {
	_distance = brLenIntialGuess ;//0.03;
	_toll = toll;
}

void fromCountTableComponentToDistance2USSRV::computeDistance() {
	likeDist2USSRV likeDist1(_model,_toll); 
	MDOUBLE initGuess = _distance;
	_distance = likeDist1.giveDistance(_ctcBase,_ctcSSRV,_likeDistance,initGuess);
	assert(_distance>=0);
}
