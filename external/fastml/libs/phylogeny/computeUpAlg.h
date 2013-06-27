// $Id: computeUpAlg.h 962 2006-11-07 15:13:34Z privmane $

#ifndef ___COMPUTE_UP_ALG
#define ___COMPUTE_UP_ALG

#include "definitions.h"
#include "tree.h"
#include "suffStatComponent.h"
#include "sequenceContainer.h"
#include "computePijComponent.h"


class computeUpAlg {
public: 
	void fillComputeUp(const tree& et,
					   const sequenceContainer& sc,
					   const int pos,
					   const computePijHom& pi,
					   suffStatGlobalHomPos& ssc);

	void fillComputeUp(const tree& et,
					 const sequenceContainer & sc,
					 const computePijGam& pi,
					 suffStatGlobalGam& ssc);

	/*void fillComputeUp(const tree& et, // not to be used at all. problematic in case of a gamma function.
					   const sequenceContainer& sc,
					   const int pos,
					   const stochasticProcess& sp,
					   suffStatGlobalHomPos& ssc);*/

	/*void fillComputeUp(const tree& et, // not to be used, accept for debuging (very slow func.)
					   const sequenceContainer& sc,
					   const stochasticProcess& sp,
					   suffStatGlobalGam& ssc);*/

	void fillComputeUpSpecificGlobalRate(const tree& et,
				   const sequenceContainer& sc,
				   const int pos,
				   const stochasticProcess& sp,
				   suffStatGlobalHomPos& ssc,
				   const MDOUBLE gRate);

// my attemp to add factors
	void fillComputeUpWithFactors(const tree& et,
				   const sequenceContainer& sc,
				   const int pos,
				   const computePijHom& pi,
				   suffStatGlobalHomPos& ssc,
				   vector<MDOUBLE>& factors);
	void fillComputeUpWithFactors(const tree& et,
				   const sequenceContainer& sc,
				   const int pos,
				   const stochasticProcess& sp,
				   suffStatGlobalHomPos& ssc,
				   vector<MDOUBLE>& factors);
	void fillComputeUpSpecificGlobalRateFactors(const tree& et,
				   const sequenceContainer& sc,
				   const int pos,
				   const stochasticProcess& sp,
				   suffStatGlobalHomPos& ssc,
				   const MDOUBLE gRate,
				   vector<MDOUBLE>& factors);
};
#endif


