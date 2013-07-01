#ifndef ___JOINT_NO_GAMMA
#define ___JOINT_NO_GAMMA

#include "definitions.h"
#include "tree.h"
#include "stochasticProcess.h"
#include "sequenceContainer.h"
#include "computePijComponent.h"
#include "suffStatComponent.h"
#include "suffStatComponentJointNoGamma.h"

class jointNoGamma {
public:
	explicit jointNoGamma(
		const tree& et,
		const stochasticProcess& sp,
		const sequenceContainer& sc);

	void compute();
	void outputTheJointProbAtEachSite(const string & outputFileProbJoint);
	sequenceContainer getTheJointReconstruction() const {return _resultSec;}

private:
	void fillComputeUp(const int pos,
				   suffStatGlobalHomPos& ssc,
				   suffStatGlobalHomPosJointNoGamma& sscJointNoGam);
	vector<int> computeJointAncestralFromSSC(
				   const int pos,
				   const suffStatGlobalHomPos& ssc,
				   const suffStatGlobalHomPosJointNoGamma& sscFASTML,
				   doubleRep & likelihoodOfReconstruction);
	void fromJointReconstructionToSequenceContainer(const vector<string> & ancestralSequences);

	const tree& _et;
	const stochasticProcess& _sp;
	const sequenceContainer& _sc;
	sequenceContainer _resultSec;
	computePijHom _cpih;
	vector<doubleRep> _jointLikelihoodOfPositions;
};



#endif
