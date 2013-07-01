// 	$Id: likeDist2USSRV.cpp 962 2006-11-07 15:13:34Z privmane $	


#include "likeDist2USSRV.h"
#include "numRec.h"


const MDOUBLE likeDist2USSRV::giveDistance(	const countTableComponentGam& ctcBase,
										   const countTableComponentHom& ctcSSRV,
										   MDOUBLE& resQ,
										   const MDOUBLE initialGuess) const {
	return giveDistanceBrent(ctcBase,ctcSSRV,resQ,initialGuess);
}


const MDOUBLE likeDist2USSRV::giveDistanceBrent(const countTableComponentGam& ctcBase,
												const countTableComponentHom& ctcSSRV,
												MDOUBLE& resL,
												const MDOUBLE initialGuess) const {
	const MDOUBLE ax=0,bx=initialGuess,cx=_maxPairwiseDistance,tol=_toll;
	LOG(12,<<"ax: " << ax << " bx: " << bx << " cx: " << cx << endl);
	MDOUBLE dist=-1.0;
	resL = -brent(ax,bx,cx,
		  C_evalLikeDist2USSRV(ctcBase,ctcSSRV,_model),
		  tol,
		  &dist);
	
	
	LOG(9, <<"brent: resL = " << resL << " dist = " << dist << endl);
	
	return dist;
}	

// @@@@dbrent doesn't work. I should try fix this
//const MDOUBLE likeDist2USSRV::giveDistanceBrent(const countTableComponentGam& ctcBase,
//												const countTableComponentHom& ctcSSRV,
//												MDOUBLE& resL,
//												const MDOUBLE initialGuess) const {
//	const MDOUBLE ax=0,bx=initialGuess,cx=_maxPairwiseDistance,tol=_toll;
//	const MDOUBLE ax_debug=0,bx_debug=initialGuess,cx_debug=_maxPairwiseDistance,tol_debug=_toll;
//	MDOUBLE dist=-1.0;
//	// @@@@ debug OZ
//	MDOUBLE dist_debug=-1.0;
//	MDOUBLE resL_debug = -brent(ax_debug,bx_debug,cx_debug,
//		  C_evalLikeDist2USSRV(ctcBase,ctcSSRV,_model),
//		  tol_debug,
//		  &dist_debug);
//	
//	resL = -dbrent(ax,bx,cx,
//		  C_evalLikeDist2USSRV(ctcBase,ctcSSRV,_model),
//		  C_evalLikeDist_d_2USSRV(ctcBase,ctcSSRV,_model),
//		  tol,
//		  &dist);
//
//	MDOUBLE small = 0.001;
//	if ((resL < resL_debug - small) || (resL_debug < resL-small) ||
//		(dist < dist_debug - small) || (dist_debug < dist-small))
//	{
//		LOG(8,<<"likeDist2USSRV::giveDistanceBrent, different results when using brent and dbrent" << endl);
//		LOG(8,<<"dbrent resL = " << resL << " , brent resL = " << resL_debug << endl);
//		LOG(8,<<"dbrent dist = " << dist << " , brent dist = " << dist_debug << endl);
//	}
//	// end of debug OZ
//	return dist;
//}
