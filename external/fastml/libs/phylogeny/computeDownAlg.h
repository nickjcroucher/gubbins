// $Id: computeDownAlg.h 3107 2007-12-27 12:38:05Z adist $

#ifndef ___COMPUTE_DOWN_ALG
#define ___COMPUTE_DOWN_ALG

#include "definitions.h"
#include "tree.h"
#include "suffStatComponent.h"
#include "sequenceContainer.h"
#include "computePijComponent.h"


class computeDownAlg {
public: 
	void fillComputeDown(const tree& et,
					   const sequenceContainer& sc,
					   const int pos,
					   const computePijHom& pi,
					   suffStatGlobalHomPos& ssc,
					   const suffStatGlobalHomPos& cup);

	void fillComputeDown(const tree& et,
					   const sequenceContainer& sc,
					   const int pos,
					   const stochasticProcess& sp,
					   suffStatGlobalHomPos& ssc,
					   const suffStatGlobalHomPos& cup);

	void fillComputeDownSpecificRate(const tree& et,
					   const sequenceContainer& sc,
					   const int pos,
					   const stochasticProcess& sp,
					   suffStatGlobalHomPos& ssc,
					   const suffStatGlobalHomPos& cup,
					   const MDOUBLE gRate);

/** compute the down computation for a non-reversible model:
	each down computation is conditioned on the state at the root.
	This means that the vector field is of one additional dimension (the alphabet at the root)
	and hence the use of the suffStatGlobalGamPos (=vector<suffStatGlobalHomPos>)
**/
	void fillComputeDownNonReversible(const tree& et,
		const sequenceContainer& sc,
		const int pos,
		const computePijHom& pi,
		suffStatGlobalGamPos& sscGivenRoot,
		const suffStatGlobalHomPos& cup);
};
#endif
