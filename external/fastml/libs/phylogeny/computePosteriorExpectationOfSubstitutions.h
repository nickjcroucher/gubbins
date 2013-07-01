
#ifndef ___COMPUTE_POSTERIOR_EXPECTATION_OF_SUBSTITUTIONS
#define ___COMPUTE_POSTERIOR_EXPECTATION_OF_SUBSTITUTIONS


/*
This is a father class where it implements the computePosteriorExpectationOfSubstitutions 
procedure for a reversible stochastic process. Its son, computePosteriorExpectationOfSubstitutions_nonReversibleSp
implements the computePosteriorExpectationOfSubstitutions for a non-reversible stochastic process. The implementation 
difference is in two functions: computePosteriorOfChangeGivenTerminals and computePosterioGivenTerminalsPerBranch
*/

#include "definitions.h"
#include "simulateJumps.h"
#include "tree.h"
#include "sequenceContainer.h"
#include "stochasticProcess.h"
#include "suffStatComponent.h"
#include "computePijComponent.h"
#include "simulateJumpsAbstract.h"

class computePosteriorExpectationOfSubstitutions {

public:
	explicit computePosteriorExpectationOfSubstitutions(const tree &tr, const sequenceContainer &sc, const stochasticProcess *sp);
	virtual ~computePosteriorExpectationOfSubstitutions(){};


	VVdouble computeExpectationAcrossTree(simulateJumpsAbstract &sim,  //input given from simulation studies
		const VVVdouble &posteriorProbs, VVVdouble &expForBranch);	
	VVdouble computePosteriorAcrossTree(simulateJumpsAbstract &sim, //input given from simulation studies
		const VVVdouble &posteriorProbsGivenTerminals,VVVdouble &probsForBranch);
	
	virtual void computePosteriorOfChangeGivenTerminals(VVVdouble &posteriorPerNodePer2States, int pos);

private:
	MDOUBLE computePosteriorOfChangePerBranch(
		simulateJumpsAbstract &sim, //input given from simulation studies
		const  VVVdouble &posteriorProbs,
		tree::nodeP node,
		int fromState, int toState);
	
	MDOUBLE computeExpectationOfChangePerBranch(
		simulateJumpsAbstract &sim, //input given from simulation studies
		const VVVdouble &posteriorProbsGivenTerminals,
		tree::nodeP node,
		int fromState, int toState);
	
	MDOUBLE computePosterioGivenTerminalsPerBranch	(int nodeId,int sonState, int fatherState,suffStatGlobalHomPos &sscUp,
		suffStatGlobalHomPos &sscDown,computePijHom &pi, MDOUBLE &LLData, const string nodeName);


protected:
	const tree &_tr;
	const sequenceContainer &_sc;
	const stochasticProcess *_sp;
};


#endif
