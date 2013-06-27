// $Id: computeMarginalAlg.h 962 2006-11-07 15:13:34Z privmane $

#ifndef ___COMPUTE_MARGINAL_ALG
#define ___COMPUTE_MARGINAL_ALG

#include "definitions.h"
#include "suffStatComponent.h"
#include "sequenceContainer.h"
#include "computePijComponent.h"

// This function will give one (for DNA, for example)
// P(A | DATA), P (C | DATA), ... etc, for each node.
// This is the case in the homogenous model only.
// for the Gamma case, the marginal in a specific node, is in fact
// p(A | DATA, r), P( C | DATA, r), ... etc.

class computeMarginalAlg {
public: 
	void fillComputeMarginal(const tree& et,
					   const sequenceContainer& sc,
					   const stochasticProcess& sp,
					   const int pos,
					   const computePijHom& pi,
					   suffStatGlobalHomPos& ssc,
					   const suffStatGlobalHomPos& cup,
					   const suffStatGlobalHomPos& cdown,
					   doubleRep & posProb);
};
#endif
