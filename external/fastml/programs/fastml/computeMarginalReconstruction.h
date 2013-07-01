#ifndef ___COMPUTE_MARGINAL_RECONSTRUCTION
#define ___COMPUTE_MARGINAL_RECONSTRUCTION

#include "definitions.h"
#include "tree.h"
#include "stochasticProcess.h"
#include "sequenceContainer.h"
#include "suffStatComponent.h"

class computeMarginalReconstruction {
public:
	explicit computeMarginalReconstruction(
		const tree& et,
		vector<stochasticProcess>& spVec,
		const sequenceContainer& sc);

	void compute(const distribution * forceDistr);
	void outputTheMarginalProbForEachCharForEachNode(const string& outputFileName);
	sequenceContainer getResultingMarginalReconstruction() const {return _resultSec;}
private:
	const tree& _et;
	vector<stochasticProcess>& _spVec;
	const sequenceContainer& _sc;
	sequenceContainer _resultSec;

	// this will be the marginal for each node, for each pos, for each letter
	VVVdouble _resultProb; //_resultProb[pos][node][letter]

	// this will be the marginal for each node, for each pos, of the best reconsturction.
	VVdouble _bestProb; //_resultProb[pos][node]

	void fillResultProb(const suffStatGlobalGamPos& ssc,const stochasticProcess & sp,const tree& et, const int pos);
	void fillMarginalReconstruction();
	void fillMarginalReconstructionSpecificNode(tree::nodeP mynode);
	void outputTheMarginalProbForEachCharForEachNodePos(ostream& out,const int pos);

};

#endif
