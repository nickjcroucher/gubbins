// $Id: likelihoodComputation2Codon.h 4699 2008-08-14 14:19:46Z privmane $

#ifndef ___LIKELIHOOD_COMPUTATION_2_CODON
#define ___LIKELIHOOD_COMPUTATION_2_CODON

#include "definitions.h"
#include "computePijComponent.h"
#include "sequenceContainer.h"
#include "suffStatComponent.h"

namespace likelihoodComputation2Codon {

	MDOUBLE getTreeLikelihoodAllPosAlphTheSame(const tree& et,
							const sequenceContainer& sc,
							const vector<stochasticProcess>& spVec,  
							const distribution * distr);

	MDOUBLE getProbOfPosUpIsFilledSelectionGam(const int pos,const tree& et, //used for gamma model
						const sequenceContainer& sc,
						const stochasticProcess& sp,
						const suffStatGlobalGamPos& cup,
						const distribution * distr);

	MDOUBLE getTreeLikelihoodFromUp2(const tree& et,
						const sequenceContainer& sc,
						const stochasticProcess& sp,
						const suffStatGlobalGam& cup,
						Vdouble& posLike, // fill this vector with each position likelihood but without the weights.
						const distribution * distr,
						const Vdouble * weights=0);
};



#endif
