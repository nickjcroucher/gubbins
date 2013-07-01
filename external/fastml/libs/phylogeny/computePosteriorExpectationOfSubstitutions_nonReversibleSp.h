#ifndef ___COMPUTE_POSTERIOR_EXPECTATION_OF_SUBSTITUTIONS_NONREVERSIBLESP
#define ___COMPUTE_POSTERIOR_EXPECTATION_OF_SUBSTITUTIONS_NONREVERSIBLESP

#include "computePosteriorExpectationOfSubstitutions.h"

class computePosteriorExpectationOfSubstitutions_nonReversibleSp:public computePosteriorExpectationOfSubstitutions {
public:
	explicit computePosteriorExpectationOfSubstitutions_nonReversibleSp(const tree &tr, const sequenceContainer &sc, stochasticProcess *sp):computePosteriorExpectationOfSubstitutions(tr,sc,sp){}
	virtual ~computePosteriorExpectationOfSubstitutions_nonReversibleSp(){};
	
	void computePosteriorOfChangeGivenTerminals(VVVdouble &posteriorPerNodePer2States, int pos);

private:
		MDOUBLE computePosterioGivenTerminalsPerBranch	(int nodeId,int sonState, int fatherState,suffStatGlobalHomPos &sscUp,
		suffStatGlobalGamPos &sscDown,computePijHom &pi, MDOUBLE &LLData, const string nodeName);

};

#endif



