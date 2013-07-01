// $Id: computeCounts.h 4545 2008-07-30 18:37:25Z cohenofi $

// version 1.00
// last modified 3 Nov 2002

#ifndef ___COMPUTE_COUNTS
#define ___COMPUTE_COUNTS

#include "definitions.h"
#include "countTableComponent.h"
#include "sequenceContainer.h"
#include "computePijComponent.h"
#include "suffStatComponent.h"

// things included for the function "fillCountTableComponentGam"
#include "sequenceContainer.h"

class computeCounts {
public:
	explicit computeCounts() {};
	void computeCountsNodeFatherNodeSonHomPos(const sequenceContainer& sc,
										const computePijHom& pi,
										const stochasticProcess& sp,
										const suffStatGlobalHomPos& cup,
										const suffStatGlobalHomPos& cdown,
										const MDOUBLE weight,
										 const doubleRep posProb,
										const tree::nodeP nodeSon,
										countTableComponentHom& _ctc,
										const MDOUBLE rateCategorProb = 1.0); //CODE_RED
	void computeCountsNodeFatherNodeSonHomPos(const sequenceContainer& sc,
		const computePijHom& pi,
		const stochasticProcess& sp,
		const suffStatGlobalHomPos& cup,
		const suffStatGlobalHomPos& cdown,
		const MDOUBLE weight,
		const doubleRep posProb,
		const tree::nodeP nodeSon,
		countTableComponentHom& _ctc,
		const MDOUBLE rateCategorProb,
		const int letterInRoot); 



	void fillCountTableComponentGam(countTableComponentGam& ctcGam,
								const stochasticProcess& sp,
								const sequenceContainer& sc,
								const computePijGam& pij0,
								const suffStatGlobalGam& cup,
								const suffStatGlobalGam& cdown,
								const Vdouble * weights,
								tree::nodeP nodeSon,
								const VdoubleRep& posProbVec);

	void fillCountTableComponentGamSpecRateCategor(const int rateCategor,
											   countTableComponentHom& ctcHom,
											   const stochasticProcess& sp,
											   const sequenceContainer& sc,
											   const computePijHom& pi,
											   const suffStatGlobalGam& cup,
												const suffStatGlobalGam& cdown,
												const Vdouble * weights,
												const VdoubleRep& posProbVec, //prob of the position with gamma
												tree::nodeP nodeSon);
};


#endif
