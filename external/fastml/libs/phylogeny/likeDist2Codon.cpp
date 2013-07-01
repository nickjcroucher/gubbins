// $RCSfile$ $Revision: 4699 $ $Date: 2008-08-14 17:19:46 +0300 (ה, 14 אוגוסט 2008) $

#include "likeDist2Codon.h"
#include "numRec.h"


const MDOUBLE likeDist2Codon::giveDistance(	const countTableComponentGam& ctc,
										   MDOUBLE& resQ,
										   const MDOUBLE initialGuess) const {
	//return giveDistanceNR(ctc,resL,initialGuess);
	return giveDistanceBrent(ctc,resQ,initialGuess);
}

const MDOUBLE likeDist2Codon::giveDistanceBrent(	const countTableComponentGam& ctc,
										   MDOUBLE& resL,
										   const MDOUBLE initialGuess) const {
	const MDOUBLE ax=0,bx=initialGuess,cx=_maxPairwiseDistance,tol=_toll;
	MDOUBLE dist=-1.0;
	resL = -dbrent(ax,bx,cx,
		  C_evalLikeDist2Codon(ctc,_spVec),
		  C_evalLikeDist_d_2Codon(ctc,_spVec),
		  tol,
		  &dist);
	return dist;
}
